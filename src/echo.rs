#![allow(non_camel_case_types)]

#[cfg(feature = "sys")]
mod sys {
    use speexdsp_sys::echo::*;
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

    #[derive(Clone)]
    pub struct SpeexEcho {
        pub(crate) st: *mut SpeexEchoState,
    }

    impl SpeexEcho {
        pub fn new(
            frame_size: usize,
            filter_length: usize,
        ) -> Result<Self, Error> {
            let st = unsafe {
                speex_echo_state_init(frame_size as i32, filter_length as i32)
            };

            if st.is_null() {
                Err(Error::FailedInit)
            } else {
                Ok(SpeexEcho { st })
            }
        }

        pub fn new_multi_channel(
            frame_size: usize,
            filter_length: usize,
            nb_mic: usize,
            nb_speakers: usize,
        ) -> Result<Self, Error> {
            let st = unsafe {
                speex_echo_state_init_mc(
                    frame_size as i32,
                    filter_length as i32,
                    nb_mic as i32,
                    nb_speakers as i32,
                )
            };

            if st.is_null() {
                Err(Error::FailedInit)
            } else {
                Ok(SpeexEcho { st })
            }
        }

        /// Performs echo cancellation a frame, based on the audio sent to the speaker (no delay is added to playback in this form).
        pub fn echo_cancellation(
            &mut self,
            rec: &[i16],
            play: &[i16],
            out: &mut [i16],
        ) {
            unsafe {
                speex_echo_cancellation(
                    self.st,
                    rec.as_ptr(),
                    play.as_ptr(),
                    out.as_mut_ptr(),
                )
            };
        }

        pub fn echo_capture(&mut self, rec: &[i16], out: &mut [i16]) {
            unsafe {
                speex_echo_capture(self.st, rec.as_ptr(), out.as_mut_ptr())
            };
        }

        pub fn echo_playback(&mut self, play: &[i16]) {
            unsafe { speex_echo_playback(self.st, play.as_ptr()) };
        }

        pub fn echo_reset(&mut self) {
            unsafe { speex_echo_state_reset(self.st) };
        }

        fn set_i32(&mut self, cmd: u32, val: i32) -> &mut Self {
            let mut v = val;
            unsafe {
                speex_echo_ctl(
                    self.st,
                    cmd as std::os::raw::c_int,
                    &mut v as *mut _ as *mut c_void,
                );
            }
            self
        }

        fn get_i32(&mut self, cmd: u32) -> i32 {
            let mut v: i32 = 0;
            unsafe {
                speex_echo_ctl(
                    self.st,
                    cmd as std::os::raw::c_int,
                    &mut v as *mut _ as *mut c_void,
                );
            }
            v
        }

        pub fn get_frame_size(&mut self) -> i32 {
            self.get_i32(SPEEX_ECHO_GET_FRAME_SIZE)
        }

        pub fn set_sampling_rate(&mut self, val: usize) -> &mut Self {
            self.set_i32(SPEEX_ECHO_SET_SAMPLING_RATE, val as i32)
        }

        pub fn get_sampling_rate(&mut self) -> i32 {
            self.get_i32(SPEEX_ECHO_GET_SAMPLING_RATE)
        }

        pub fn get_impulse_response_size(&mut self) -> i32 {
            self.get_i32(SPEEX_ECHO_GET_IMPULSE_RESPONSE_SIZE)
        }

        // todo: SPEEX_ECHO_GET_IMPULSE_RESPONSE
    }

    impl Drop for SpeexEcho {
        fn drop(&mut self) {
            unsafe { speex_echo_state_destroy(self.st) };
        }
    }

    pub struct SpeexDecorr {
        st: *mut SpeexDecorrState,
    }

    impl SpeexDecorr {
        pub fn new(
            rate: usize,
            channels: usize,
            frame_size: usize,
        ) -> Result<Self, Error> {
            let st = unsafe {
                speex_decorrelate_new(
                    rate as i32,
                    channels as i32,
                    frame_size as i32,
                )
            };

            if st.is_null() {
                Err(Error::FailedInit)
            } else {
                Ok(SpeexDecorr { st })
            }
        }

        pub fn decorrelate(
            &mut self,
            input: &[i16],
            out: &mut [i16],
            strength: usize,
        ) {
            unsafe {
                speex_decorrelate(
                    self.st,
                    input.as_ptr(),
                    out.as_mut_ptr(),
                    strength as i32,
                )
            };
        }
    }

    impl Drop for SpeexDecorr {
        fn drop(&mut self) {
            unsafe { speex_decorrelate_destroy(self.st) };
        }
    }
}

#[cfg(feature = "sys")]
pub use self::sys::{Error, SpeexDecorr, SpeexEcho};
