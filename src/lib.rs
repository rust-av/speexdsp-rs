#![feature ( const_slice_as_ptr , extern_types , libc )]
#![feature(extern_prelude)]

#[cfg(feature="sys")]
extern crate speexdsp_sys;
extern crate libc;

pub mod resampler;
mod speex_resample;
