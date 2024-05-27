use std::time::Instant;

use super::*;
use crate::{aligned_memory::AlignedMemory64, ntt::build_ntt_tables, params::Params, util::*};
use rand::Rng;

fn get_params() -> Params {
    get_test_params()
}

const REF_VAL: u64 = 519370102;

#[test]
fn build_ntt_tables_correct() {
    println!("using AVX impl: {}", get_ntt_impl());
    let moduli = [268369921u64, 249561089u64];
    let poly_len = 2048usize;
    let res = build_ntt_tables(poly_len, moduli.as_slice(), None);
    assert_eq!(res.len(), 2);
    assert_eq!(res[0].len(), 4);
    assert_eq!(res[0][0].len(), poly_len);
    assert_eq!(res[0][2][0], 134184961u64);
    assert_eq!(res[0][2][1], 96647580u64);
    let mut x1 = 0u64;
    for i in 0..res.len() {
        for j in 0..res[0].len() {
            for k in 0..res[0][0].len() {
                x1 ^= res[i][j][k];
            }
        }
    }
    assert_eq!(x1, REF_VAL);
}

#[test]
fn ntt_forward_correct() {
    let params = get_params();
    let mut v1 = AlignedMemory64::new(2 * 2048);
    v1[0] = 100;
    v1[2048] = 100;
    ntt_forward(&params, v1.as_mut_slice());
    assert_eq!(v1[50], 100);
    assert_eq!(v1[2048 + 50], 100);
}

#[test]
fn ntt_inverse_correct() {
    let params = get_params();
    let mut v1 = AlignedMemory64::new(2 * 2048);
    for i in 0..v1.len() {
        v1[i] = 100;
    }
    ntt_inverse(&params, v1.as_mut_slice());
    assert_eq!(v1[0], 100);
    assert_eq!(v1[2048], 100);
    assert_eq!(v1[50], 0);
    assert_eq!(v1[2048 + 50], 0);
}

#[test]
fn ntt_correct() {
    let params = get_params();
    let trials = 10000;
    let mut v1 = AlignedMemory64::new(trials * params.crt_count * params.poly_len);
    let mut rng = rand::thread_rng();
    for trial in 0..trials {
        for i in 0..params.crt_count {
            for j in 0..params.poly_len {
                let idx = calc_index(&[trial, i, j], &[trials, params.crt_count, params.poly_len]);
                let val: u64 = rng.gen();
                v1[idx] = val % params.moduli[i];
            }
        }
    }
    let mut v2 = v1.clone();
    for chunk in v2
        .as_mut_slice()
        .chunks_exact_mut(params.crt_count * params.poly_len)
    {
        ntt_forward(&params, chunk);
    }

    let now = Instant::now();
    for chunk in v2
        .as_mut_slice()
        .chunks_exact_mut(params.crt_count * params.poly_len)
    {
        ntt_inverse(&params, chunk);
    }
    println!("{} trials of ntti took: {:?}", trials, now.elapsed());

    let now = Instant::now();
    for chunk in v2
        .as_mut_slice()
        .chunks_exact_mut(params.crt_count * params.poly_len)
    {
        ntt_forward(&params, chunk);
    }
    println!("{} trials of ntt took: {:?}", trials, now.elapsed());

    for chunk in v2
        .as_mut_slice()
        .chunks_exact_mut(params.crt_count * params.poly_len)
    {
        ntt_inverse(&params, chunk);
    }

    for i in 0..trials * params.crt_count * params.poly_len {
        assert_eq!(v1[i], v2[i]);
    }
}

fn get_alt_params() -> Params {
    Params::init(
        2048,
        &vec![180143985094819841u64],
        6.4,
        2,
        256,
        20,
        4,
        8,
        56,
        8,
        true,
        9,
        6,
        1,
        2048,
        0,
    )
}

#[test]
fn alt_ntt_correct() {
    let params = get_alt_params();
    let mut v1 = AlignedMemory64::new(params.crt_count * params.poly_len);
    let mut rng = rand::thread_rng();
    for i in 0..params.crt_count {
        for j in 0..params.poly_len {
            let idx = calc_index(&[i, j], &[params.crt_count, params.poly_len]);
            let val: u64 = rng.gen();
            v1[idx] = val % params.moduli[i];
        }
    }
    let mut v2 = v1.clone();
    ntt_forward(&params, v2.as_mut_slice());
    ntt_inverse(&params, v2.as_mut_slice());
    for i in 0..params.crt_count * params.poly_len {
        assert_eq!(v1[i], v2[i]);
    }
}
