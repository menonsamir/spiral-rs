use crate::params::Params;

pub fn ntt_forward_alt(params: &Params, operand_overall: &mut [u64]) {
    let log_n = params.poly_len_log2;
    let n = 1 << log_n;

    for coeff_mod in 0..params.crt_count {
        let operand = &mut operand_overall[coeff_mod * n..coeff_mod * n + n];

        let forward_table = params.get_ntt_forward_table(coeff_mod);
        let forward_table_prime = params.get_ntt_forward_prime_table(coeff_mod);
        let modulus_small = params.moduli[coeff_mod];
        let two_times_modulus_small = 2 * modulus_small;

        for mm in 0..log_n {
            let m = 1 << mm;
            let t = n >> (mm + 1);

            let mut it = operand.chunks_exact_mut(2 * t);

            for i in 0..m {
                let w = forward_table[m + i];
                let w_prime = forward_table_prime[m + i];

                let op = it.next().unwrap();

                for j in 0..t {
                    let x: u64 = op[j] as u64;
                    let y: u64 = op[t + j] as u64;

                    let curr_x: u64 =
                        x - (two_times_modulus_small * ((x >= two_times_modulus_small) as u64));
                    let q_tmp = ((y as u128) * (w_prime as u128)) >> 64u64;
                    let q_new = (w as u128) * (y as u128) - q_tmp * (modulus_small as u128);
                    let q_new = (q_new % (modulus_small as u128)) as u64;

                    op[j] = curr_x as u64 + q_new;
                    op[t + j] = curr_x as u64 + ((two_times_modulus_small as u64) - q_new);
                }
            }
        }

        for i in 0..n {
            operand[i] -= ((operand[i] >= two_times_modulus_small as u64) as u64)
                * two_times_modulus_small as u64;
            operand[i] -= ((operand[i] >= modulus_small as u64) as u64) * modulus_small as u64;
        }
    }
}

pub fn ntt_inverse_alt(params: &Params, operand_overall: &mut [u64]) {
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

                for j in 0..t {
                    let x = op[j];
                    let y = op[t + j];

                    let t_tmp = two_times_modulus - y + x;
                    let curr_x = x + y - (two_times_modulus * (((x << 1) >= t_tmp) as u64));
                    let h_tmp = ((t_tmp as u128) * (w_prime as u128)) >> 64;

                    let res_x = (curr_x + (modulus * ((t_tmp & 1) as u64))) >> 1;
                    let res_y = ((w as u128) * (t_tmp as u128)) - (h_tmp * modulus as u128);

                    op[j] = res_x;
                    op[t + j] = (res_y % (modulus as u128)) as u64;
                }
            }
        }

        for i in 0..n {
            operand[i] -= ((operand[i] >= two_times_modulus) as u64) * two_times_modulus;
            operand[i] -= ((operand[i] >= modulus) as u64) * modulus;
        }
    }
}
