use std::arch::x86_64::*;

use crate::params::*;

use super::{ntt_forward_alt, ntt_inverse_alt};

pub fn ntt_forward(params: &Params, operand_overall: &mut [u64]) {
    if params.crt_count == 1 {
        ntt_forward_alt(params, operand_overall);
        return;
    }
    let log_n = params.poly_len_log2;
    let n = 1 << log_n;

    for coeff_mod in 0..params.crt_count {
        let operand = unsafe {
            std::slice::from_raw_parts_mut(operand_overall.as_mut_ptr().add(coeff_mod * n), n)
        };

        let forward_table = params.get_ntt_forward_table(coeff_mod);
        let forward_table_prime = params.get_ntt_forward_prime_table(coeff_mod);
        let modulus_small = params.moduli[coeff_mod] as u32;
        let two_times_modulus_small: u32 = 2 * modulus_small;

        for mm in 0..log_n {
            let m = 1 << mm;
            let t = n >> (mm + 1);

            for i in 0..m {
                let w = unsafe { *forward_table.get_unchecked(m + i) };
                let w_prime = unsafe { *forward_table_prime.get_unchecked(m + i) };

                let op = unsafe {
                    std::slice::from_raw_parts_mut(operand.as_mut_ptr().add(2 * t * i), 2 * t)
                };

                if t < 4 || log_n <= 10 {
                    for j in 0..t {
                        let x: u32 = unsafe { *op.get_unchecked(j) as u32 };
                        let y: u32 = unsafe { *op.get_unchecked(t + j) as u32 };

                        let curr_x: u32 =
                            x - (two_times_modulus_small * ((x >= two_times_modulus_small) as u32));
                        let q_tmp: u64 = ((y as u64) * (w_prime as u64)) >> 32u64;
                        let q_new = w * (y as u64) - q_tmp * (modulus_small as u64);

                        unsafe {
                            *op.get_unchecked_mut(j) = curr_x as u64 + q_new;
                            *op.get_unchecked_mut(t + j) =
                                curr_x as u64 + ((two_times_modulus_small as u64) - q_new);
                        }
                    }
                } else {
                    // t >= 4
                    unsafe {
                        for j in (0..t).step_by(4) {
                            // Use AVX2 here
                            let p_x = op.get_unchecked_mut(j) as *mut u64;
                            let p_y = op.get_unchecked_mut(j + t) as *mut u64;
                            let x = _mm256_load_si256(p_x as *const __m256i);
                            let y = _mm256_load_si256(p_y as *const __m256i);

                            let cmp_val = _mm256_set1_epi64x(two_times_modulus_small as i64);
                            let gt_mask = _mm256_cmpgt_epi64(x, cmp_val);

                            let to_subtract = _mm256_and_si256(gt_mask, cmp_val);
                            let curr_x = _mm256_sub_epi64(x, to_subtract);

                            // uint32_t q_val = ((y) * (uint64_t)(Wprime)) >> 32;
                            let w_prime_vec = _mm256_set1_epi64x(w_prime as i64);
                            let product = _mm256_mul_epu32(y, w_prime_vec);
                            let q_val = _mm256_srli_epi64(product, 32);

                            // q_val = W * y - q_val * modulus_small;
                            let w_vec = _mm256_set1_epi64x(w as i64);
                            let w_times_y = _mm256_mul_epu32(y, w_vec);
                            let modulus_small_vec = _mm256_set1_epi64x(modulus_small as i64);
                            let q_scaled = _mm256_mul_epu32(q_val, modulus_small_vec);
                            let q_final = _mm256_sub_epi64(w_times_y, q_scaled);

                            let new_x = _mm256_add_epi64(curr_x, q_final);
                            let q_final_inverted = _mm256_sub_epi64(cmp_val, q_final);
                            let new_y = _mm256_add_epi64(curr_x, q_final_inverted);

                            _mm256_store_si256(p_x as *mut __m256i, new_x);
                            _mm256_store_si256(p_y as *mut __m256i, new_y);
                        }
                    }
                }
            }
        }

        if log_n <= 10 {
            for i in 0..n {
                operand[i] -= ((operand[i] >= two_times_modulus_small as u64) as u64)
                    * two_times_modulus_small as u64;
                operand[i] -= ((operand[i] >= modulus_small as u64) as u64) * modulus_small as u64;
            }
            continue;
        }

        for i in (0..n).step_by(4) {
            unsafe {
                let p_x = &mut operand[i] as *mut u64;

                let cmp_val1 = _mm256_set1_epi64x(two_times_modulus_small as i64);
                let mut x = _mm256_load_si256(p_x as *const __m256i);
                let mut gt_mask = _mm256_cmpgt_epi64(x, cmp_val1);
                let mut to_subtract = _mm256_and_si256(gt_mask, cmp_val1);
                x = _mm256_sub_epi64(x, to_subtract);

                let cmp_val2 = _mm256_set1_epi64x(modulus_small as i64);
                gt_mask = _mm256_cmpgt_epi64(x, cmp_val2);
                to_subtract = _mm256_and_si256(gt_mask, cmp_val2);
                x = _mm256_sub_epi64(x, to_subtract);
                _mm256_store_si256(p_x as *mut __m256i, x);
            }
        }
    }
}

pub fn ntt_inverse(params: &Params, operand_overall: &mut [u64]) {
    if params.crt_count == 1 {
        ntt_inverse_alt(params, operand_overall);
        return;
    }
    for coeff_mod in 0..params.crt_count {
        let n = params.poly_len;

        let operand = &mut operand_overall[coeff_mod * n..coeff_mod * n + n];

        let inverse_table = params.get_ntt_inverse_table(coeff_mod);
        let inverse_table_prime = params.get_ntt_inverse_prime_table(coeff_mod);
        let modulus = params.moduli[coeff_mod];
        let two_times_modulus: u64 = 2 * modulus;
        for mm in (0..params.poly_len_log2).rev() {
            let h = 1 << mm;
            let t = n >> (mm + 1);

            let mut it = operand.chunks_exact_mut(2 * t);

            for i in 0..h {
                let w = inverse_table[h + i];
                let w_prime = inverse_table_prime[h + i];

                let op = it.next().unwrap();

                if t < 4 {
                    for j in 0..t {
                        let x = op[j];
                        let y = op[t + j];

                        let t_tmp = two_times_modulus - y + x;
                        let curr_x = x + y - (two_times_modulus * (((x << 1) >= t_tmp) as u64));
                        let h_tmp = (t_tmp * w_prime) >> 32;

                        let res_x = (curr_x + (modulus * ((t_tmp & 1) as u64))) >> 1;
                        let res_y = w * t_tmp - h_tmp * modulus;

                        op[j] = res_x;
                        op[t + j] = res_y;
                    }
                } else {
                    unsafe {
                        for j in (0..t).step_by(4) {
                            // Use AVX2 here
                            let p_x = &mut op[j] as *mut u64;
                            let p_y = &mut op[j + t] as *mut u64;
                            let x = _mm256_load_si256(p_x as *const __m256i);
                            let y = _mm256_load_si256(p_y as *const __m256i);

                            let modulus_vec = _mm256_set1_epi64x(modulus as i64);
                            let two_times_modulus_vec =
                                _mm256_set1_epi64x(two_times_modulus as i64);
                            let mut t_tmp = _mm256_set1_epi64x(two_times_modulus as i64);
                            t_tmp = _mm256_sub_epi64(t_tmp, y);
                            t_tmp = _mm256_add_epi64(t_tmp, x);
                            let gt_mask = _mm256_cmpgt_epi64(_mm256_slli_epi64(x, 1), t_tmp);
                            let to_subtract = _mm256_and_si256(gt_mask, two_times_modulus_vec);
                            let mut curr_x = _mm256_add_epi64(x, y);
                            curr_x = _mm256_sub_epi64(curr_x, to_subtract);

                            let w_prime_vec = _mm256_set1_epi64x(w_prime as i64);
                            let mut h_tmp = _mm256_mul_epu32(t_tmp, w_prime_vec);
                            h_tmp = _mm256_srli_epi64(h_tmp, 32);

                            let and_mask = _mm256_set_epi64x(1, 1, 1, 1);
                            let eq_mask =
                                _mm256_cmpeq_epi64(_mm256_and_si256(t_tmp, and_mask), and_mask);
                            let to_add = _mm256_and_si256(eq_mask, modulus_vec);

                            let new_x = _mm256_srli_epi64(_mm256_add_epi64(curr_x, to_add), 1);

                            let w_vec = _mm256_set1_epi64x(w as i64);
                            let w_times_t_tmp = _mm256_mul_epu32(t_tmp, w_vec);
                            let h_tmp_times_modulus = _mm256_mul_epu32(h_tmp, modulus_vec);
                            let new_y = _mm256_sub_epi64(w_times_t_tmp, h_tmp_times_modulus);

                            _mm256_store_si256(p_x as *mut __m256i, new_x);
                            _mm256_store_si256(p_y as *mut __m256i, new_y);
                        }
                    }
                }
            }
        }

        // for i in 0..n {
        //     operand[i] -= ((operand[i] >= two_times_modulus) as u64) * two_times_modulus;
        //     operand[i] -= ((operand[i] >= modulus) as u64) * modulus;
        // }

        for i in (0..n).step_by(4) {
            unsafe {
                let p_x = &mut operand[i] as *mut u64;

                let cmp_val1 = _mm256_set1_epi64x(two_times_modulus as i64);
                let mut x = _mm256_load_si256(p_x as *const __m256i);
                let mut gt_mask = _mm256_or_si256(
                    _mm256_cmpgt_epi64(x, cmp_val1),
                    _mm256_cmpeq_epi64(x, cmp_val1),
                );
                let mut to_subtract = _mm256_and_si256(gt_mask, cmp_val1);
                x = _mm256_sub_epi64(x, to_subtract);

                let cmp_val2 = _mm256_set1_epi64x(modulus as i64);
                gt_mask = _mm256_or_si256(
                    _mm256_cmpgt_epi64(x, cmp_val2),
                    _mm256_cmpeq_epi64(x, cmp_val2),
                );
                to_subtract = _mm256_and_si256(gt_mask, cmp_val2);
                x = _mm256_sub_epi64(x, to_subtract);
                _mm256_store_si256(p_x as *mut __m256i, x);
            }
        }
    }
}

pub fn get_ntt_impl() -> String {
    "AVX2".to_string()
}
