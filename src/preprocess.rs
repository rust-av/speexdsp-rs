#[cfg(feature = "sys")]
mod sys {
    use speexdsp_sys::preprocess::*;
    use std::ffi::c_void;
    use std::fmt;

    #[derive(Clone, Copy, Debug)]
    pub enum Error {
        FailedInit,
    }

    impl fmt::Display for Error {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            let v = match self {
                Error::FailedInit => "Failed to initialize",
            };

            write!(f, "{}", v)
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
                    sampling_rate as i32
                )
            };

            if st.is_null() {
                Err(Error::FailedInit)
            } else {
                Ok(SpeexPreprocess { st })
            }
        }

        pub fn preprocess_run(&mut self, x: usize) -> usize {
            unsafe {
                speex_preprocess_run(self.st, x as *mut i16) as usize
            }
        }

        pub fn preprocess(&mut self, x: usize, echo: usize) -> usize {
            unsafe {
                speex_preprocess(
                    self.st,
                    x as *mut i16,
                    echo as *mut i32,
                ) as usize
            }
        }

        pub fn preprocess_estimate_update(&mut self, x: usize) {
            unsafe {
                speex_preprocess_estimate_update(
                    self.st,
                    x as *mut i16,
                )
            };
        }

        pub fn preprocess_ctl(&mut self, request: usize, ptr: &mut f32) -> usize {
            let ptr_v = ptr as *mut f32 as *mut c_void;
            unsafe {
                speex_preprocess_ctl(self.st, request as i32, ptr_v) as usize
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
pub use self::sys::{Error, SpeexPreprocess};
