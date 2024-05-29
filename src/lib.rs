#![cfg_attr(not(feature = "skip-avx512"), feature(avx512_target_feature))]
#![cfg_attr(not(feature = "skip-avx512"), feature(stdarch_x86_avx512))]

pub mod aligned_memory;
pub mod arith;
pub mod discrete_gaussian;
pub mod noise_estimate;
pub mod number_theory;
pub mod util;

pub mod gadget;
pub mod ntt;
pub mod params;
pub mod poly;

pub mod client;
pub mod key_value;

#[cfg(feature = "server")]
pub mod server;
