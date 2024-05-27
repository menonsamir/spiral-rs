mod alt;
mod tables;
pub use alt::*;
pub use tables::*;

// unusual, but we want to test whichever implementation is available
#[cfg(test)]
mod test;

// // if we have avx-512
#[cfg(target_feature = "avx512f")]
mod avx512;
#[cfg(target_feature = "avx512f")]
pub use avx512::*;

// if we have avx2 but not avx-512
#[cfg(all(target_feature = "avx2", not(target_feature = "avx512f")))]
mod avx2;
#[cfg(all(target_feature = "avx2", not(target_feature = "avx512f")))]
pub use avx2::*;

// otherwise
#[cfg(not(any(target_feature = "avx2", target_feature = "avx512f")))]
mod base;
#[cfg(not(any(target_feature = "avx2", target_feature = "avx512f")))]
pub use base::*;
