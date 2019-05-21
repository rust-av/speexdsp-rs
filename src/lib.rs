#[cfg(feature="sys")]
extern crate speexdsp_sys;
extern crate libc;

pub mod resampler;
mod speex_resample;
