# spiral-rs

This is a Rust implementation of some functionality in the [Spiral PIR scheme](https://eprint.iacr.org/2022/368). This fork provides core routines for use in [YPIR](https://github.com/menonsamir/ypir).

For a complete, working version of spiral-rs, please see [this repository](https://github.com/blyssprivacy/sdk/tree/main/lib).

## Building

You must have AVX512 to build the main branch of this repository. For an implemnetation that does not require AVX512, switch to the `avoid-avx512` branch.