#[cfg(feature = "sys")]
mod sys {
    use echo::SpeexEcho;
    use speexdsp_sys::preprocess::*;
    use std::convert::From;
    use std::ffi::c_void;
    use std::fmt;

    #[derive(Clone, Copy, Debug)]
    pub enum SpeexPreprocessConst {
        SPEEX_PREPROCESS_SET_DENOISE = 0,
        SPEEX_PREPROCESS_GET_DENOISE = 1,
        SPEEX_PREPROCESS_SET_AGC = 2,
        SPEEX_PREPROCESS_GET_AGC = 3,
        SPEEX_PREPROCESS_SET_VAD = 4,
        SPEEX_PREPROCESS_GET_VAD = 5,
        SPEEX_PREPROCESS_SET_AGC_LEVEL = 6,
        SPEEX_PREPROCESS_GET_AGC_LEVEL = 7,
        SPEEX_PREPROCESS_SET_DEREVERB = 8,
        SPEEX_PREPROCESS_GET_DEREVERB = 9,
        SPEEX_PREPROCESS_SET_DEREVERB_LEVEL = 10,
        SPEEX_PREPROCESS_GET_DEREVERB_LEVEL = 11,
        SPEEX_PREPROCESS_SET_DEREVERB_DECAY = 12,
        SPEEX_PREPROCESS_GET_DEREVERB_DECAY = 13,
        SPEEX_PREPROCESS_SET_PROB_START = 14,
        SPEEX_PREPROCESS_GET_PROB_START = 15,
        SPEEX_PREPROCESS_SET_PROB_CONTINUE = 16,
        SPEEX_PREPROCESS_GET_PROB_CONTINUE = 17,
        SPEEX_PREPROCESS_SET_NOISE_SUPPRESS = 18,
        SPEEX_PREPROCESS_GET_NOISE_SUPPRESS = 19,
        SPEEX_PREPROCESS_SET_ECHO_SUPPRESS = 20,
        SPEEX_PREPROCESS_GET_ECHO_SUPPRESS = 21,
        SPEEX_PREPROCESS_SET_ECHO_SUPPRESS_ACTIVE = 22,
        SPEEX_PREPROCESS_GET_ECHO_SUPPRESS_ACTIVE = 23,
        SPEEX_PREPROCESS_SET_ECHO_STATE = 24,
        SPEEX_PREPROCESS_GET_ECHO_STATE = 25,
        SPEEX_PREPROCESS_SET_AGC_INCREMENT = 26,
        SPEEX_PREPROCESS_GET_AGC_INCREMENT = 27,
        SPEEX_PREPROCESS_SET_AGC_DECREMENT = 28,
        SPEEX_PREPROCESS_GET_AGC_DECREMENT = 29,
        SPEEX_PREPROCESS_SET_AGC_MAX_GAIN = 30,
        SPEEX_PREPROCESS_GET_AGC_MAX_GAIN = 31,
        SPEEX_PREPROCESS_GET_AGC_LOUDNESS = 33,
        SPEEX_PREPROCESS_GET_AGC_GAIN = 35,
        SPEEX_PREPROCESS_GET_PSD_SIZE = 37,
        SPEEX_PREPROCESS_GET_PSD = 39,
        SPEEX_PREPROCESS_GET_NOISE_PSD_SIZE = 41,
        SPEEX_PREPROCESS_GET_NOISE_PSD = 43,
        SPEEX_PREPROCESS_GET_PROB = 45,
        SPEEX_PREPROCESS_SET_AGC_TARGET = 46,
        SPEEX_PREPROCESS_GET_AGC_TARGET = 47,
    }

    #[derive(Clone, Copy, Debug)]
    pub enum Error {
        FailedInit,
        UnknownRequest,
    }

    impl fmt::Display for Error {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            let v = match self {
                Error::FailedInit => "Failed to initialize",
                Error::UnknownRequest => "The request is unknown",
            };

            write!(f, "{}", v)
        }
    }

    pub enum Variant {
        U32(u32),
        F32(f32),
        Echo(SpeexEcho),
    }

    impl From<u32> for Variant {
        fn from(item: u32) -> Self {
            Variant::U32(item)
        }
    }

    impl From<f32> for Variant {
        fn from(item: f32) -> Self {
            Variant::F32(item)
        }
    }

    impl From<&SpeexEcho> for Variant {
        fn from(item: &SpeexEcho) -> Self {
            Variant::Echo(item.clone())
        }
    }

    pub struct SpeexPreprocess {
        st: *mut SpeexPreprocessState,
    }

    impl SpeexPreprocess {
        pub fn new(
            frame_size: usize,
            sampling_rate: usize,
        ) -> Result<Self, Error> {
            let st = unsafe {
                speex_preprocess_state_init(
                    frame_size as i32,
                    sampling_rate as i32,
                )
            };

            if st.is_null() {
                Err(Error::FailedInit)
            } else {
                Ok(SpeexPreprocess { st })
            }
        }

        pub fn preprocess_run(&mut self, x: &mut [i16]) -> usize {
            unsafe { speex_preprocess_run(self.st, x.as_mut_ptr()) as usize }
        }

        pub fn preprocess(&mut self, x: &mut [i16], echo: usize) -> usize {
            unsafe {
                speex_preprocess(self.st, x.as_mut_ptr(), echo as *mut i32)
                    as usize
            }
        }

        pub fn preprocess_estimate_update(&mut self, x: &mut [i16]) {
            unsafe {
                speex_preprocess_estimate_update(self.st, x.as_mut_ptr())
            };
        }

        pub fn preprocess_ctl<T: Into<Variant>>(
            &mut self,
            request: SpeexPreprocessConst,
            value: T,
        ) -> Result<(), Error> {
            let ptr_v = match (request, value.into()) {
                (
                    SpeexPreprocessConst::SPEEX_PREPROCESS_SET_DEREVERB_DECAY,
                    Variant::F32(val),
                )
                | (
                    SpeexPreprocessConst::SPEEX_PREPROCESS_SET_DEREVERB_LEVEL,
                    Variant::F32(val),
                ) => &val as *const f32 as *mut c_void,
                (
                    SpeexPreprocessConst::SPEEX_PREPROCESS_SET_ECHO_STATE,
                    Variant::Echo(val),
                ) => val.get_ptr() as *mut c_void,
                (_, Variant::U32(val)) => val as *mut c_void,
                _ => panic!("This type is not accepted"),
            };
            let ret = unsafe {
                speex_preprocess_ctl(self.st, request as i32, ptr_v) as usize
            };
            if ret != 0 {
                Err(Error::UnknownRequest)
            } else {
                Ok(())
            }
        }
    }

    impl Drop for SpeexPreprocess {
        fn drop(&mut self) {
            unsafe { speex_preprocess_state_destroy(self.st) };
        }
    }
}

#[cfg(feature = "sys")]
pub use self::sys::{Error, SpeexPreprocess, SpeexPreprocessConst};
