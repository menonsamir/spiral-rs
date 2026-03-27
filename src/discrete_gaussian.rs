use rand::distributions::WeightedIndex;
use rand::prelude::Distribution;
use rand::Rng;
use rand_chacha::ChaCha20Rng;
use subtle::ConditionallySelectable;
use subtle::ConstantTimeGreater;

use crate::poly::*;
use std::f64::consts::PI;

pub const NUM_WIDTHS: usize = 4;

/// Table of u64 values representing a Gaussian of width 6.4
/// (standard deviation = 6.4/sqrt(2*pi))
///
/// This is the cumulative distribution function of this distribution,
/// in the range [-26, 26], multiplied by 2^64. Values exactly equal to 2^64 have
/// been replaced with 2^64-1, for representation as u64's.
// const CDF_TABLE_GAUS_6_4: [u64; 53] = [
//     0,
//     0,
//     0,
//     7,
//     225,
//     6114,
//     142809,
//     2864512,
//     49349166,
//     730367088,
//     9288667698,
//     101545086850,
//     954617134063,
//     7720973857474,
//     53757667977838,
//     322436486442815,
//     1667499996257361,
//     7443566871362048,
//     28720140744863884,
//     95948302954529081,
//     278161926109627739,
//     701795634139702303,
//     1546646853635104741,
//     2991920295851131431,
//     5112721055115151939,
//     7782220156096217088,
//     10664523917613334528,
//     13334023018594399677,
//     15454823777858420185,
//     16900097220074446875,
//     17744948439569849313,
//     18168582147599923877,
//     18350795770755022535,
//     18418023932964687732,
//     18439300506838189568,
//     18445076573713294255,
//     18446421637223108801,
//     18446690316041573778,
//     18446736352735694142,
//     18446743119092417553,
//     18446743972164464766,
//     18446744064420883918,
//     18446744072979184528,
//     18446744073660202450,
//     18446744073706687104,
//     18446744073709408807,
//     18446744073709545502,
//     18446744073709551391,
//     18446744073709551609,
//     18446744073709551615,
//     18446744073709551615,
//     18446744073709551615,
//     18446744073709551615,
// ];

pub struct DiscreteGaussian {
    pub weighted_index: WeightedIndex<f64>,
    pub cdf_table: Vec<u64>,
    pub max_val: i64,
}

impl DiscreteGaussian {
    pub fn init(noise_width: f64) -> Self {
        let max_val = (noise_width * (NUM_WIDTHS as f64)).ceil() as i64;
        let mut table = Vec::new();
        let mut total = 0.0;

        // assign discrete probabilities to each possible integer output
        for i in -max_val..max_val + 1 {
            let p_val = f64::exp(-PI * f64::powi(i as f64, 2) / f64::powi(noise_width, 2));
            table.push(p_val);
            total += p_val;
        }

        // build a CDF table for possible outputs
        let mut cdf_table = Vec::new();
        let mut cum_prob = 0.0;

        for p_val in &table {
            cum_prob += p_val / total;
            let cum_prob_u64 = (cum_prob * (u64::MAX as f64)).round() as u64;
            cdf_table.push(cum_prob_u64);
        }

        Self {
            weighted_index: WeightedIndex::new(table).unwrap(),
            cdf_table,
            max_val,
        }
    }

    pub fn sample(&self, modulus: u64, rng: &mut ChaCha20Rng) -> u64 {
        let sampled_val = rng.gen::<u64>();
        let len = (2 * self.max_val + 1) as usize;
        let mut to_output = 0;

        for i in (0..len).rev() {
            let mut out_val = (i as i64) - self.max_val;
            // this branch is ok: not secret-dependent
            if out_val < 0 {
                out_val += modulus as i64;
            }
            let out_val_u64 = out_val as u64;

            // let point = CDF_TABLE_GAUS_6_4[i];
            let point = self.cdf_table[i];

            // if sampled_val <= point, set to_output := out_val
            // (in constant time)
            let cmp = !(sampled_val.ct_gt(&point));
            to_output.conditional_assign(&out_val_u64, cmp);
        }
        to_output
    }

    /// Sample from a discrete Gaussian distribution. THIS IS NOT CONSTANT TIME!
    pub fn fast_sample(&self, modulus: u64, rng: &mut ChaCha20Rng) -> u64 {
        let sampled_val = self.weighted_index.sample(rng);
        let mut val = (sampled_val as i64) - self.max_val;
        if val < 0 {
            val += modulus as i64;
        }
        val as u64
    }

    pub fn sample_matrix(&self, p: &mut PolyMatrixRaw, rng: &mut ChaCha20Rng) {
        let modulus = p.get_params().modulus;
        for r in 0..p.rows {
            for c in 0..p.cols {
                let poly = p.get_poly_mut(r, c);
                for z in 0..poly.len() {
                    let s = self.sample(modulus, rng);
                    poly[z] = s;
                }
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::util::*;
    use rand::SeedableRng;

    #[test]
    fn dg_seems_okay() {
        let params = get_test_params();
        let dg = DiscreteGaussian::init(params.noise_width);
        let mut rng = get_chacha_rng();
        let mut v = Vec::new();
        let trials = 10000;
        let mut sum = 0;
        for _ in 0..trials {
            let val = dg.sample(params.modulus, &mut rng);
            let mut val_i64 = val as i64;
            if val_i64 >= (params.modulus as i64) / 2 {
                val_i64 -= params.modulus as i64;
            }
            v.push(val_i64);
            sum += val_i64;
        }
        let expected_mean = 0;
        let computed_mean = sum as f64 / trials as f64;
        let expected_std_dev = params.noise_width / f64::sqrt(2f64 * std::f64::consts::PI);
        let std_dev_of_mean = expected_std_dev / f64::sqrt(trials as f64);
        println!("mean:: expected: {}, got: {}", expected_mean, computed_mean);
        assert!(f64::abs(computed_mean) < std_dev_of_mean * 5f64);

        let computed_variance: f64 = v
            .iter()
            .map(|x| (computed_mean - (*x as f64)).powi(2))
            .sum::<f64>()
            / (v.len() as f64);
        let computed_std_dev = computed_variance.sqrt();
        println!(
            "std_dev:: expected: {}, got: {}",
            expected_std_dev, computed_std_dev
        );
        assert!((computed_std_dev - expected_std_dev).abs() < (expected_std_dev * 0.1));
    }

    #[test]
    fn cdf_table_seems_okay() {
        let dg = DiscreteGaussian::init(6.4);
        println!("{:?}", dg.cdf_table);
    }

    // ----------------------------------------------------------------
    // CDF table structural invariants
    // ----------------------------------------------------------------

    #[test]
    fn cdf_table_is_monotonically_nondecreasing() {
        for width in [1.0, 3.2, 6.4, 10.0, 25.0] {
            let dg = DiscreteGaussian::init(width);
            for i in 1..dg.cdf_table.len() {
                assert!(
                    dg.cdf_table[i] >= dg.cdf_table[i - 1],
                    "CDF not monotonic at index {} for width {}: {} < {}",
                    i,
                    width,
                    dg.cdf_table[i],
                    dg.cdf_table[i - 1]
                );
            }
        }
    }

    #[test]
    fn cdf_table_last_entry_near_u64_max() {
        for width in [1.0, 3.2, 6.4, 10.0, 25.0] {
            let dg = DiscreteGaussian::init(width);
            let last = *dg.cdf_table.last().unwrap();
            let threshold = u64::MAX - (u64::MAX / 1000);
            assert!(
                last >= threshold,
                "CDF last entry {} not near u64::MAX for width {}",
                last,
                width
            );
        }
    }

    #[test]
    fn cdf_table_has_correct_length() {
        for width in [1.0, 3.2, 6.4, 10.0, 25.0] {
            let dg = DiscreteGaussian::init(width);
            let expected_len = (2 * dg.max_val + 1) as usize;
            assert_eq!(
                dg.cdf_table.len(),
                expected_len,
                "CDF table length mismatch for width {}: expected {}, got {}",
                width,
                expected_len,
                dg.cdf_table.len()
            );
        }
    }

    #[test]
    fn cdf_table_symmetric_around_center() {
        let dg = DiscreteGaussian::init(6.4);
        let len = dg.cdf_table.len();
        let center = len / 2;
        for i in 1..center {
            let left_mass = dg.cdf_table[center - i];
            let right_complement = u64::MAX - dg.cdf_table[center + i - 1];
            let diff = if left_mass > right_complement {
                left_mass - right_complement
            } else {
                right_complement - left_mass
            };
            let tolerance = u64::MAX / 1000;
            assert!(
                diff < tolerance,
                "CDF not symmetric at offset {}: left_mass={}, right_complement={}",
                i,
                left_mass,
                right_complement
            );
        }
    }

    // ----------------------------------------------------------------
    // max_val computation
    // ----------------------------------------------------------------

    #[test]
    fn max_val_scales_with_noise_width() {
        let dg_small = DiscreteGaussian::init(1.0);
        let dg_large = DiscreteGaussian::init(10.0);
        assert!(
            dg_large.max_val > dg_small.max_val,
            "larger noise width should produce larger max_val"
        );
        assert_eq!(
            dg_small.max_val,
            (1.0 * NUM_WIDTHS as f64).ceil() as i64
        );
        assert_eq!(
            dg_large.max_val,
            (10.0 * NUM_WIDTHS as f64).ceil() as i64
        );
    }

    // ----------------------------------------------------------------
    // Deterministic seeded output
    // ----------------------------------------------------------------

    #[test]
    fn sample_is_deterministic_with_same_seed() {
        let dg = DiscreteGaussian::init(6.4);
        let modulus = 268369921u64;
        let seed = get_chacha_static_seed();

        let mut rng1 = ChaCha20Rng::from_seed(seed);
        let mut rng2 = ChaCha20Rng::from_seed(seed);

        for i in 0..1000 {
            let v1 = dg.sample(modulus, &mut rng1);
            let v2 = dg.sample(modulus, &mut rng2);
            assert_eq!(v1, v2, "sample diverged at iteration {} with same seed", i);
        }
    }

    #[test]
    fn fast_sample_is_deterministic_with_same_seed() {
        let dg = DiscreteGaussian::init(6.4);
        let modulus = 268369921u64;
        let seed = get_chacha_static_seed();

        let mut rng1 = ChaCha20Rng::from_seed(seed);
        let mut rng2 = ChaCha20Rng::from_seed(seed);

        for i in 0..1000 {
            let v1 = dg.fast_sample(modulus, &mut rng1);
            let v2 = dg.fast_sample(modulus, &mut rng2);
            assert_eq!(
                v1, v2,
                "fast_sample diverged at iteration {} with same seed",
                i
            );
        }
    }

    #[test]
    fn different_seeds_produce_different_sequences() {
        let dg = DiscreteGaussian::init(6.4);
        let modulus = 268369921u64;

        let mut rng1 = ChaCha20Rng::from_seed([0u8; 32]);
        let mut rng2 = ChaCha20Rng::from_seed([1u8; 32]);

        let mut same_count = 0;
        let trials = 100;
        for _ in 0..trials {
            if dg.sample(modulus, &mut rng1) == dg.sample(modulus, &mut rng2) {
                same_count += 1;
            }
        }
        assert!(
            same_count < trials,
            "different seeds produced identical sequences"
        );
    }

    // ----------------------------------------------------------------
    // Output range validation
    // ----------------------------------------------------------------

    #[test]
    fn sample_output_always_in_modulus_range() {
        let dg = DiscreteGaussian::init(6.4);
        let mut rng = ChaCha20Rng::from_seed(get_chacha_static_seed());

        for modulus in [101u64, 1021, 65537, 268369921, 249561089] {
            for _ in 0..1000 {
                let val = dg.sample(modulus, &mut rng);
                assert!(
                    val < modulus,
                    "sample {} out of range [0, {}) for modulus {}",
                    val,
                    modulus,
                    modulus
                );
            }
        }
    }

    #[test]
    fn fast_sample_output_always_in_modulus_range() {
        let dg = DiscreteGaussian::init(6.4);
        let mut rng = ChaCha20Rng::from_seed(get_chacha_static_seed());

        for modulus in [101u64, 1021, 65537, 268369921, 249561089] {
            for _ in 0..1000 {
                let val = dg.fast_sample(modulus, &mut rng);
                assert!(
                    val < modulus,
                    "fast_sample {} out of range [0, {})",
                    val,
                    modulus
                );
            }
        }
    }

    // ----------------------------------------------------------------
    // Extreme noise widths
    // ----------------------------------------------------------------

    #[test]
    fn very_small_noise_width() {
        let dg = DiscreteGaussian::init(0.5);
        assert!(dg.max_val >= 1);
        assert!(dg.cdf_table.len() >= 3);

        let mut rng = ChaCha20Rng::from_seed(get_chacha_static_seed());
        let modulus = 268369921u64;
        for _ in 0..500 {
            let val = dg.sample(modulus, &mut rng);
            assert!(val < modulus);
        }
    }

    #[test]
    fn large_noise_width() {
        let dg = DiscreteGaussian::init(100.0);
        assert_eq!(dg.max_val, (100.0 * NUM_WIDTHS as f64).ceil() as i64);
        assert_eq!(dg.cdf_table.len(), (2 * dg.max_val + 1) as usize);

        let mut rng = ChaCha20Rng::from_seed(get_chacha_static_seed());
        let modulus = 268369921u64;
        for _ in 0..500 {
            let val = dg.sample(modulus, &mut rng);
            assert!(val < modulus);
        }
    }

    // ----------------------------------------------------------------
    // Boundary RNG values via constant-time sample path
    // ----------------------------------------------------------------

    #[test]
    fn sample_with_rng_returning_zero() {
        let dg = DiscreteGaussian::init(6.4);
        let modulus = 268369921u64;

        let mut rng = ChaCha20Rng::from_seed([0u8; 32]);
        let val = dg.sample(modulus, &mut rng);
        assert!(val < modulus);
    }

    #[test]
    fn sample_with_rng_returning_high_values() {
        let dg = DiscreteGaussian::init(6.4);
        let modulus = 268369921u64;

        let mut rng = ChaCha20Rng::from_seed([0xFFu8; 32]);
        let val = dg.sample(modulus, &mut rng);
        assert!(val < modulus);
    }

    // ----------------------------------------------------------------
    // sample_matrix fills all entries
    // ----------------------------------------------------------------

    #[test]
    fn sample_matrix_fills_all_entries() {
        let params = get_test_params();
        let dg = DiscreteGaussian::init(params.noise_width);
        let mut rng = ChaCha20Rng::from_seed(get_chacha_static_seed());

        let mut mat = PolyMatrixRaw::zero(&params, 2, 2);
        dg.sample_matrix(&mut mat, &mut rng);

        let total_entries = mat.data.len();
        let zero_count = mat.data.as_slice().iter().filter(|&&x| x == 0).count();
        assert!(
            zero_count < total_entries,
            "sample_matrix left all entries at zero"
        );
        let zero_ratio = zero_count as f64 / total_entries as f64;
        assert!(
            zero_ratio < 0.5,
            "sample_matrix left {:.1}% entries at zero, expected most to be nonzero",
            zero_ratio * 100.0
        );
    }

    #[test]
    fn sample_matrix_values_in_modulus_range() {
        let params = get_test_params();
        let dg = DiscreteGaussian::init(params.noise_width);
        let mut rng = ChaCha20Rng::from_seed(get_chacha_static_seed());

        let mut mat = PolyMatrixRaw::zero(&params, 2, 2);
        dg.sample_matrix(&mut mat, &mut rng);

        for &val in mat.data.as_slice() {
            assert!(
                val < params.modulus,
                "sample_matrix produced {} >= modulus {}",
                val,
                params.modulus
            );
        }
    }

    // ----------------------------------------------------------------
    // Distribution agreement between sample and fast_sample
    // ----------------------------------------------------------------

    #[test]
    fn sample_and_fast_sample_produce_similar_distributions() {
        let dg = DiscreteGaussian::init(6.4);
        let modulus = 268369921u64;
        let trials = 20000;

        let mut rng_ct = ChaCha20Rng::from_seed(get_chacha_static_seed());
        let mut rng_fast = ChaCha20Rng::from_seed([42u8; 32]);

        let to_signed = |val: u64| -> i64 {
            let v = val as i64;
            if v >= (modulus as i64) / 2 {
                v - modulus as i64
            } else {
                v
            }
        };

        let mut sum_ct: i64 = 0;
        let mut sum_sq_ct: f64 = 0.0;
        let mut sum_fast: i64 = 0;
        let mut sum_sq_fast: f64 = 0.0;

        for _ in 0..trials {
            let v_ct = to_signed(dg.sample(modulus, &mut rng_ct));
            sum_ct += v_ct;
            sum_sq_ct += (v_ct as f64).powi(2);

            let v_fast = to_signed(dg.fast_sample(modulus, &mut rng_fast));
            sum_fast += v_fast;
            sum_sq_fast += (v_fast as f64).powi(2);
        }

        let mean_ct = sum_ct as f64 / trials as f64;
        let mean_fast = sum_fast as f64 / trials as f64;
        let var_ct = sum_sq_ct / trials as f64 - mean_ct.powi(2);
        let var_fast = sum_sq_fast / trials as f64 - mean_fast.powi(2);

        let expected_std = 6.4 / f64::sqrt(2.0 * PI);
        let expected_var = expected_std.powi(2);

        assert!(
            mean_ct.abs() < 1.0,
            "constant-time sample mean {} too far from 0",
            mean_ct
        );
        assert!(
            mean_fast.abs() < 1.0,
            "fast_sample mean {} too far from 0",
            mean_fast
        );

        assert!(
            (var_ct - expected_var).abs() < expected_var * 0.2,
            "constant-time sample variance {} too far from expected {}",
            var_ct,
            expected_var
        );
        assert!(
            (var_fast - expected_var).abs() < expected_var * 0.2,
            "fast_sample variance {} too far from expected {}",
            var_fast,
            expected_var
        );
    }
}
