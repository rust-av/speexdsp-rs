#[cfg(feature = "sys")]
mod sys {
    use speexdsp_sys::resampler::*;
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
            unsafe {
                speex_resampler_set_rate(
                    self.st,
                    in_rate as u32,
                    out_rate as u32,
                )
            };
        }

        pub fn get_rate(&self) -> (usize, usize) {
            let mut in_rate = 0;
            let mut out_rate = 0;

            unsafe {
                speex_resampler_get_rate(self.st, &mut in_rate, &mut out_rate)
            };

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

        pub fn set_quality(&self, quality: usize) -> Result<(), Error> {
            let ret = unsafe {
                speex_resampler_set_quality(self.st, quality as i32)
            };
            if ret != 0 {
                Err(Error::from_i32(ret))
            } else {
                Ok(())
            }
        }

        pub fn get_quality(&self) -> usize {
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

#[cfg(feature = "sys")]
pub use self::sys::{Error, State};

pub mod native {
    use speex_resample::*;

    pub struct State {
        st: SpeexResamplerState,
    }

    impl State {
        pub fn new(
            channels: usize,
            in_rate: usize,
            out_rate: usize,
            quality: u8,
        ) -> Self {
            let st = speex_resampler_init(
                channels as u32,
                in_rate as u32,
                out_rate as u32,
                quality as i32,
            );
            Self { st }
        }

        pub fn set_rate(&mut self, in_rate: usize, out_rate: usize) {
            speex_resampler_set_rate(
                &mut self.st,
                in_rate as u32,
                out_rate as u32,
            );
        }

        pub fn get_rate(&mut self) -> (usize, usize) {
            let mut in_rate = 0;
            let mut out_rate = 0;

            speex_resampler_get_rate(
                &mut self.st,
                &mut in_rate,
                &mut out_rate,
            );

            (in_rate as usize, out_rate as usize)
        }

        pub fn process_float(
            &mut self,
            index: usize,
            input: &[f32],
            output: &mut [f32],
        ) -> (usize, usize) {
            let mut in_len = input.len() as u32;
            let mut out_len = output.len() as u32;
            speex_resampler_process_float(
                &mut self.st,
                index as u32,
                input,
                &mut in_len,
                output,
                &mut out_len,
            );

            (in_len as usize, out_len as usize)
        }

        pub fn skip_zeros(&mut self) {
            speex_resampler_skip_zeros(&mut self.st);
        }

        pub fn reset(&mut self) {
            speex_resampler_reset_mem(&mut self.st);
        }

        pub fn get_input_latency(&mut self) -> usize {
            speex_resampler_get_input_latency(&mut self.st) as usize
        }

        pub fn get_output_latency(&mut self) -> usize {
            speex_resampler_get_output_latency(&mut self.st) as usize
        }

        pub fn set_quality(&mut self, quality: usize) {
            speex_resampler_set_quality(&mut self.st, quality as i32);
        }

        pub fn get_quality(&mut self) -> usize {
            let mut c_get = 0;
            speex_resampler_get_quality(&mut self.st, &mut c_get);
            c_get as usize
        }
    }
}

#[cfg(not(feature = "sys"))]
pub use self::native::State;
