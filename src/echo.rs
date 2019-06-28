#[cfg(feature = "sys")]
mod sys {
    use speexdsp_sys::echo::*;
    use std::ffi::c_void;
    use std::fmt;

    #[derive(Clone, Copy, Debug)]
    pub enum SpeexEchoConst {
        SPEEX_ECHO_GET_FRAME_SIZE = 3,
        SPEEX_ECHO_SET_SAMPLING_RATE = 24,
        SPEEX_ECHO_GET_SAMPLING_RATE = 25,
        SPEEX_ECHO_GET_IMPULSE_RESPONSE_SIZE = 27,
        SPEEX_ECHO_GET_IMPULSE_RESPONSE = 29,
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

    #[derive(Clone)]
    pub struct SpeexEcho {
        st: *mut SpeexEchoState,
    }

    impl SpeexEcho {
        pub fn new(frame_size: usize, filter_length: usize) -> Result<Self, Error> {
            let st = unsafe { speex_echo_state_init(frame_size as i32, filter_length as i32) };

            if st.is_null() {
                Err(Error::FailedInit)
            } else {
                Ok(SpeexEcho { st })
            }
        }

        pub fn echo_init(
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

        pub fn echo_cancellation(&mut self, rec: &[i16], play: &[i16], out: &mut [i16]) {
            unsafe {
                speex_echo_cancellation(self.st, rec.as_ptr(), play.as_ptr(), out.as_mut_ptr())
            };
        }

        pub fn echo_cancel(
            &mut self,
            rec: &[i16],
            play: &[i16],
            out: &mut [i16],
            yout: &mut [i32],
        ) {
            unsafe {
                speex_echo_cancel(
                    self.st,
                    rec.as_ptr(),
                    play.as_ptr(),
                    out.as_mut_ptr(),
                    yout.as_mut_ptr(),
                )
            };
        }

        pub fn echo_capture(&mut self, rec: &[i16], out: &mut [i16]) {
            unsafe { speex_echo_capture(self.st, rec.as_ptr(), out.as_mut_ptr()) };
        }

        pub fn echo_playback(&mut self, play: &[i16]) {
            unsafe { speex_echo_playback(self.st, play.as_ptr()) };
        }

        pub fn echo_reset(&mut self) {
            unsafe { speex_echo_state_reset(self.st) };
        }

        pub fn echo_ctl(&mut self, request: SpeexEchoConst, ptr: usize) -> Result<(), Error> {
            let ret =
                unsafe { speex_echo_ctl(self.st, request as i32, ptr as *mut c_void) as usize };
            if ret != 0 {
                Err(Error::UnknownRequest)
            } else {
                Ok(())
            }
        }

        pub(crate) fn get_ptr(&self) -> *mut SpeexEchoState {
            self.st
        }
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
        pub fn new(rate: usize, channels: usize, frame_size: usize) -> Result<Self, Error> {
            let st =
                unsafe { speex_decorrelate_new(rate as i32, channels as i32, frame_size as i32) };

            if st.is_null() {
                Err(Error::FailedInit)
            } else {
                Ok(SpeexDecorr { st })
            }
        }

        pub fn decorrelate(&mut self, input: &[i16], out: &mut [i16], strength: usize) {
            unsafe {
                speex_decorrelate(self.st, input.as_ptr(), out.as_mut_ptr(), strength as i32)
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
pub use self::sys::{Error, SpeexDecorr, SpeexEcho, SpeexEchoConst};
