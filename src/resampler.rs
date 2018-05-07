use speexdsp_sys::resampler::*;

use std::ffi::CStr;
use std::fmt;

#[derive(Clone, Copy, Debug)]
pub enum Error {
    AllocFailed = 1,
    BadState = 2,
    InvalidArg = 3,
    PtrOverlap = 4,
    Overflow = 5,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let v = unsafe { CStr::from_ptr(speex_resampler_strerror(*self as i32)) };

        write!(f, "{}", v.to_string_lossy())
    }
}

impl Error {
    fn from_i32(v: i32) -> Self {
        match v as u32 {
            RESAMPLER_ERR_ALLOC_FAILED => Error::AllocFailed,
            RESAMPLER_ERR_BAD_STATE => Error::BadState,
            RESAMPLER_ERR_INVALID_ARG => Error::InvalidArg,
            RESAMPLER_ERR_PTR_OVERLAP => Error::PtrOverlap,
            RESAMPLER_ERR_OVERFLOW => Error::Overflow,
            _ => unreachable!(),
        }
    }
}

pub struct State {
    st: *mut SpeexResamplerState,
}

impl State {
    pub fn new(
        channels: usize,
        in_rate: usize,
        out_rate: usize,
        quality: u8,
    ) -> Result<Self, Error> {
        let mut err = 0;
        let st = unsafe {
            speex_resampler_init(
                channels as u32,
                in_rate as u32,
                out_rate as u32,
                quality as i32,
                &mut err,
            )
        };

        if st.is_null() {
            Err(Error::from_i32(err))
        } else {
            Ok(State { st })
        }
    }

    pub fn set_rate(&mut self, in_rate: usize, out_rate: usize) {
        unsafe { speex_resampler_set_rate(self.st, in_rate as u32, out_rate as u32) };
    }

    pub fn get_rate(&self) -> (usize, usize) {
        let mut in_rate = 0;
        let mut out_rate = 0;

        unsafe { speex_resampler_get_rate(self.st, &mut in_rate, &mut out_rate) };

        (in_rate as usize, out_rate as usize)
    }

    pub fn process_float(
        &mut self,
        index: usize,
        input: &[f32],
        output: &mut [f32],
    ) -> Result<(usize, usize), Error> {
        let mut in_len = input.len() as u32;
        let mut out_len = output.len() as u32;
        let ret = unsafe {
            speex_resampler_process_float(
                self.st,
                index as u32,
                input.as_ptr(),
                &mut in_len,
                output.as_mut_ptr(),
                &mut out_len,
            )
        };

        if ret != 0 {
            Err(Error::from_i32(ret))
        } else {
            Ok((in_len as usize, out_len as usize))
        }
    }

    pub fn skip_zeros(&mut self) {
        unsafe { speex_resampler_skip_zeros(self.st) };
    }

    pub fn reset(&mut self) {
        unsafe { speex_resampler_reset_mem(self.st) };
    }

    pub fn get_input_latency(&self) -> usize {
        unsafe { speex_resampler_get_input_latency(self.st) as usize }
    }

    pub fn get_output_latency(&self) -> usize {
        unsafe { speex_resampler_get_output_latency(self.st) as usize }
    }
}

impl Drop for State {
    fn drop(&mut self) {
        unsafe { speex_resampler_destroy(self.st) };
    }
}
