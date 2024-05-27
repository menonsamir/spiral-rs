use super::{ntt_forward_alt, ntt_inverse_alt};
use crate::params::*;

pub fn ntt_forward(params: &Params, operand_overall: &mut [u64]) {
    if params.crt_count == 1 {
        ntt_forward_alt(params, operand_overall);
        return;
    }
    let log_n = params.poly_len_log2;
    let n = 1 << log_n;

    for coeff_mod in 0..params.crt_count {
        let forward_table = params.get_ntt_forward_table(coeff_mod).as_ptr();
        let forward_table_prime = params.get_ntt_forward_prime_table(coeff_mod).as_ptr();
        let modulus_small = params.moduli[coeff_mod] as u32;
        let two_times_modulus_small: u32 = 2 * modulus_small;
        let op_mut_ptr = unsafe { operand_overall.as_mut_ptr().add(coeff_mod * n) };

        for mm in 0..log_n {
            let m = 1 << mm;
            let t = n >> (mm + 1);

            for i in 0..m {
                let w = unsafe { *forward_table.add(m + i) };
                let w_prime = unsafe { *forward_table_prime.add(m + i) };

                let op = unsafe { op_mut_ptr.add(i * 2 * t) };

                for j in 0..t {
                    let x: u32 = unsafe { *op.add(j) as u32 };
                    let y: u32 = unsafe { *op.add(t + j) as u32 };

                    let mut curr_x = x;
                    if x >= two_times_modulus_small {
                        curr_x -= two_times_modulus_small;
                    }
                    let q_tmp = ((y as u64) * (w_prime as u64)) >> 32u64;
                    let q_new = w * (y as u64) - q_tmp * (modulus_small as u64);

                    let res_j = curr_x as u64 + q_new;
                    let res_t_j = curr_x as u64 + ((two_times_modulus_small as u64) - q_new);

                    unsafe {
                        *op.add(j) = res_j;
                        *op.add(t + j) = res_t_j;
                    }
                }
            }
        }

        for i in 0..n {
            unsafe {
                let operand_i = op_mut_ptr.add(i);
                *operand_i -= ((*operand_i >= two_times_modulus_small as u64) as u64)
                    * two_times_modulus_small as u64;
                *operand_i -= ((*operand_i >= modulus_small as u64) as u64) * modulus_small as u64;
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

        let inverse_table = params.get_ntt_inverse_table(coeff_mod).as_ptr();
        let inverse_table_prime = params.get_ntt_inverse_prime_table(coeff_mod).as_ptr();
        let modulus = params.moduli[coeff_mod];
        let two_times_modulus: u64 = 2 * modulus;
        let op_mut_ptr = unsafe { operand_overall.as_mut_ptr().add(coeff_mod * n) };

        for mm in (0..params.poly_len_log2).rev() {
            let h = 1 << mm;
            let t = n >> (mm + 1);

            for i in 0..h {
                let w = unsafe { *inverse_table.add(h + i) };
                let w_prime = unsafe { *inverse_table_prime.add(h + i) };

                let op = unsafe { op_mut_ptr.add(i * 2 * t) };

                for j in 0..t {
                    let x = unsafe { *op.add(j) };
                    let y = unsafe { *op.add(t + j) };

                    let t_tmp = two_times_modulus - y + x;
                    let curr_x = x + y - (two_times_modulus * (((x << 1) >= t_tmp) as u64));
                    let h_tmp = (t_tmp * w_prime) >> 32;

                    let res_x = (curr_x + (modulus * ((t_tmp & 1) as u64))) >> 1;
                    let res_y = w * t_tmp - h_tmp * modulus;

                    unsafe {
                        *op.add(j) = res_x;
                        *op.add(t + j) = res_y;
                    }
                }
            }
        }

        for i in 0..n {
            unsafe {
                let operand_i = op_mut_ptr.add(i);
                *operand_i -= ((*operand_i >= two_times_modulus) as u64) * two_times_modulus;
                *operand_i -= ((*operand_i >= modulus) as u64) * modulus;
            }
        }
    }
}

pub fn get_ntt_impl() -> String {
    "BASE".to_string()
}
