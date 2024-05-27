use crate::{arith::*, number_theory::*};

pub fn powers_of_primitive_root(root: u64, modulus: u64, poly_len_log2: usize) -> Vec<u64> {
    let poly_len = 1usize << poly_len_log2;
    let mut root_powers = vec![0u64; poly_len];
    let mut power = root;
    for i in 1..poly_len {
        let idx = reverse_bits(i as u64, poly_len_log2) as usize;
        root_powers[idx] = power;
        power = multiply_uint_mod(power, root, modulus);
    }
    root_powers[0] = 1;
    root_powers
}

pub fn scale_powers_u64(modulus: u64, poly_len: usize, inp: &[u64]) -> Vec<u64> {
    let mut scaled_powers = vec![0; poly_len];
    for i in 0..poly_len {
        let wide_val = (inp[i] as u128) << 64u128;
        let quotient = wide_val / (modulus as u128);
        scaled_powers[i] = quotient as u64;
    }
    scaled_powers
}

pub fn scale_powers_u32(modulus: u32, poly_len: usize, inp: &[u64]) -> Vec<u64> {
    let mut scaled_powers = vec![0; poly_len];
    for i in 0..poly_len {
        let wide_val = inp[i] << 32;
        let quotient = wide_val / (modulus as u64);
        scaled_powers[i] = (quotient as u32) as u64;
    }
    scaled_powers
}

pub fn build_ntt_tables_alt(
    poly_len: usize,
    moduli: &[u64],
    opt_roots: Option<&[u64]>,
) -> Vec<Vec<Vec<u64>>> {
    let poly_len_log2 = log2(poly_len as u64) as usize;
    let mut output: Vec<Vec<Vec<u64>>> = vec![Vec::new(); moduli.len()];
    for coeff_mod in 0..moduli.len() {
        let modulus = moduli[coeff_mod];
        let root = if let Some(roots) = opt_roots {
            roots[coeff_mod]
        } else {
            get_minimal_primitive_root(2 * poly_len as u64, modulus).unwrap()
        };
        let inv_root = invert_uint_mod(root, modulus).unwrap();

        let root_powers = powers_of_primitive_root(root, modulus, poly_len_log2);
        let scaled_root_powers = scale_powers_u64(modulus, poly_len, root_powers.as_slice());
        let mut inv_root_powers = powers_of_primitive_root(inv_root, modulus, poly_len_log2);
        for i in 0..poly_len {
            inv_root_powers[i] = div2_uint_mod(inv_root_powers[i], modulus);
        }
        let scaled_inv_root_powers =
            scale_powers_u64(modulus, poly_len, inv_root_powers.as_slice());

        output[coeff_mod] = vec![
            root_powers,
            scaled_root_powers,
            inv_root_powers,
            scaled_inv_root_powers,
        ];
    }
    output
}

pub fn build_ntt_tables(
    poly_len: usize,
    moduli: &[u64],
    opt_roots: Option<&[u64]>,
) -> Vec<Vec<Vec<u64>>> {
    let poly_len_log2 = log2(poly_len as u64) as usize;
    let mut output: Vec<Vec<Vec<u64>>> = vec![Vec::new(); moduli.len()];
    for coeff_mod in 0..moduli.len() {
        let modulus = moduli[coeff_mod];
        let modulus_as_u32 = modulus.try_into().unwrap();
        let root = if let Some(roots) = opt_roots {
            roots[coeff_mod]
        } else {
            get_minimal_primitive_root(2 * poly_len as u64, modulus).unwrap()
        };
        let inv_root = invert_uint_mod(root, modulus).unwrap();

        let root_powers = powers_of_primitive_root(root, modulus, poly_len_log2);
        let scaled_root_powers = scale_powers_u32(modulus_as_u32, poly_len, root_powers.as_slice());
        let mut inv_root_powers = powers_of_primitive_root(inv_root, modulus, poly_len_log2);
        for i in 0..poly_len {
            inv_root_powers[i] = div2_uint_mod(inv_root_powers[i], modulus);
        }
        let scaled_inv_root_powers =
            scale_powers_u32(modulus_as_u32, poly_len, inv_root_powers.as_slice());

        output[coeff_mod] = vec![
            root_powers,
            scaled_root_powers,
            inv_root_powers,
            scaled_inv_root_powers,
        ];
    }
    output
}
