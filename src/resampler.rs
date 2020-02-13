#[derive(Clone, Copy, Debug)]
pub enum Error {
    AllocFailed = 1,
    BadState = 2,
    InvalidArg = 3,
    PtrOverlap = 4,
    Overflow = 5,
}

use std::fmt;

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let v = match self {
            Error::AllocFailed => "Memory allocation failed.",
            Error::BadState => "Bad resampler state.",
            Error::InvalidArg => "Invalid argument.",
            Error::PtrOverlap => "Input and output buffers overlap.",
            Error::Overflow => "Muldiv overflow.",
        };

        write!(f, "{}", v)
    }
}

pub trait Resampler: Sized {
    fn new(channels: usize, in_rate: usize, out_rate: usize, quality: usize) -> Result<Self, Error>;
    fn set_rate(&mut self, in_rate: usize, out_rate: usize) -> Result<(), Error>;
    fn get_rate(&self) -> (usize, usize);
    fn get_ratio(&self) -> (usize, usize);
    fn process_float(
        &mut self,
        index: usize,
        input: &[f32],
        output: &mut [f32],
    ) -> Result<(usize, usize), Error>;
    fn skip_zeros(&mut self);
    fn reset(&mut self);
    fn get_input_latency(&self) -> usize;
    fn get_output_latency(&self) -> usize;
    fn set_quality(&mut self, quality: usize) -> Result<(), Error>;
    fn get_quality(&self) -> usize;
}

#[cfg(feature = "sys")]
mod sys {
    use super::{Error, Resampler};
    use speexdsp_sys::resampler::*;

    impl From<i32> for Error {
        fn from(v: i32) -> Error {
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

    impl Resampler for State {
        fn new(
            channels: usize,
            in_rate: usize,
            out_rate: usize,
            quality: usize,
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
                Err(err.into())
            } else {
                Ok(State { st })
            }
        }

        fn set_rate(&mut self, in_rate: usize, out_rate: usize) -> Result<(), Error> {
            let ret = unsafe { speex_resampler_set_rate(self.st, in_rate as u32, out_rate as u32) };
            if ret != 0 {
                Err(ret.into())
            } else {
                Ok(())
            }
        }

        fn get_rate(&self) -> (usize, usize) {
            let mut in_rate = 0;
            let mut out_rate = 0;

            unsafe { speex_resampler_get_rate(self.st, &mut in_rate, &mut out_rate) };

            (in_rate as usize, out_rate as usize)
        }

        fn get_ratio(&self) -> (usize, usize) {
            let mut num = 0;
            let mut den = 0;

            unsafe { speex_resampler_get_ratio(self.st, &mut num, &mut den) };

            (num as usize, den as usize)
        }

        fn process_float(
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
                Err(ret.into())
            } else {
                Ok((in_len as usize, out_len as usize))
            }
        }

        fn skip_zeros(&mut self) {
            unsafe { speex_resampler_skip_zeros(self.st) };
        }

        fn reset(&mut self) {
            unsafe { speex_resampler_reset_mem(self.st) };
        }

        fn get_input_latency(&self) -> usize {
            unsafe { speex_resampler_get_input_latency(self.st) as usize }
        }

        fn get_output_latency(&self) -> usize {
            unsafe { speex_resampler_get_output_latency(self.st) as usize }
        }

        fn set_quality(&mut self, quality: usize) -> Result<(), Error> {
            let ret = unsafe { speex_resampler_set_quality(self.st, quality as i32) };
            if ret != 0 {
                Err(ret.into())
            } else {
                Ok(())
            }
        }

        fn get_quality(&self) -> usize {
            let mut c_get = 0;
            unsafe { speex_resampler_get_quality(self.st, &mut c_get) };
            c_get as usize
        }
    }

    impl Drop for State {
        fn drop(&mut self) {
            unsafe { speex_resampler_destroy(self.st) };
        }
    }
}

pub mod native {
    pub use speexdsp_resampler::State;
    use super::{Resampler, Error};

    impl From<speexdsp_resampler::Error> for Error {
        fn from(v: speexdsp_resampler::Error) -> Error {
            match v {
                speexdsp_resampler::Error::AllocFailed => Error::AllocFailed,
                speexdsp_resampler::Error::InvalidArg => Error::InvalidArg,
            }
        }
    }

    impl Resampler for State {
        fn new(
            channels: usize,
            in_rate: usize,
            out_rate: usize,
            quality: usize,
        ) -> Result<Self, Error> {
            State::new(channels, in_rate, out_rate, quality)
                .map_err(|e| e.into())
        }
        fn set_rate(&mut self, in_rate: usize, out_rate: usize) -> Result<(), Error> {
            State::set_rate(self, in_rate, out_rate).map_err(|e| e.into())
        }
        fn get_rate(&self) -> (usize, usize) {
            State::get_rate(self)
        }
        fn get_ratio(&self) -> (usize, usize) {
            State::get_ratio(self)
        }
        fn process_float(
            &mut self,
            index: usize,
            input: &[f32],
            output: &mut [f32],
        ) -> Result<(usize, usize), Error> {
            State::process_float(self, index, input, output)
                .map_err(|e| e.into())
        }
        fn skip_zeros(&mut self) {
            State::skip_zeros(self);
        }
        fn reset(&mut self) {
            State::reset(self);
        }
        fn get_input_latency(&self) -> usize {
            State::get_input_latency(self)
        }
        fn get_output_latency(&self) -> usize {
            State::get_output_latency(self)
        }
        fn set_quality(&mut self, quality: usize) -> Result<(), Error> {
            State::set_quality(self, quality).map_err(|e| e.into())
        }
        fn get_quality(&self) -> usize {
            State::get_quality(self)
        }
    }
}

#[cfg(feature = "sys")]
pub use self::sys::*;

#[cfg(not(feature = "sys"))]
pub use self::native::*;
