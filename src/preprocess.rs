#![allow(non_camel_case_types)]

#[cfg(feature = "sys")]
mod sys {
    use speexdsp_sys::echo::SpeexEchoState;
    use speexdsp_sys::preprocess::*;
    use std::ffi::c_void;
    use std::fmt;

    use crate::echo::SpeexEcho;

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

    /// State of the preprocessor (one per channel).
    pub struct SpeexPreprocess {
        st: *mut SpeexPreprocessState,
    }

    impl SpeexPreprocess {
        /// Creates a new preprocessing state. You MUST create one state per channel processed.
        ///
        /// frame_size is the number of samples to process at one time (should correspond to 10-20 ms). Must be
        /// the same value as that used for the echo canceller for residual echo cancellation to work.
        ///
        /// sampling_rate Sampling rate used for the input.
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

        /// Preprocess a frame.
        /// The buffer must be same size as specified in [`SpeexPreprocess::new`].
        ///
        /// Return a bool value for voice activity (true for speech, false for noise/silence), ONLY if VAD turned on.
        pub fn preprocess_run(&mut self, x: &mut [i16]) -> bool {
            unsafe { speex_preprocess_run(self.st, x.as_mut_ptr()) != 0 }
        }

        /// Update preprocessor state, but do not compute the output.
        /// The buffer must be same size as specified in [`SpeexPreprocess::new`].
        pub fn preprocess_estimate_update(&mut self, x: &mut [i16]) {
            unsafe {
                speex_preprocess_estimate_update(self.st, x.as_mut_ptr())
            };
        }
    }

    impl Drop for SpeexPreprocess {
        fn drop(&mut self) {
            unsafe { speex_preprocess_state_destroy(self.st) };
        }
    }

    impl SpeexPreprocess {
        fn set_bool(&mut self, cmd: u32, enable: bool) -> &mut Self {
            let mut v: i32 = if enable { 1 } else { 0 };
            unsafe {
                speex_preprocess_ctl(
                    self.st,
                    cmd as std::os::raw::c_int,
                    &mut v as *mut _ as *mut c_void,
                );
            }
            self
        }

        fn get_bool(&mut self, cmd: u32) -> bool {
            let mut v: i32 = 0;
            unsafe {
                speex_preprocess_ctl(
                    self.st,
                    cmd as std::os::raw::c_int,
                    &mut v as *mut _ as *mut c_void,
                );
            }
            v != 0
        }

        fn set_f32(&mut self, cmd: u32, mut val: f32) -> &mut Self {
            unsafe {
                speex_preprocess_ctl(
                    self.st,
                    cmd as std::os::raw::c_int,
                    &mut val as *mut _ as *mut c_void,
                );
            }
            self
        }

        fn get_f32(&mut self, cmd: u32) -> f32 {
            let mut v: f32 = 0.0;
            unsafe {
                speex_preprocess_ctl(
                    self.st,
                    cmd as std::os::raw::c_int,
                    &mut v as *mut _ as *mut c_void,
                );
            }
            v
        }

        fn set_i32(&mut self, cmd: u32, mut val: i32) -> &mut Self {
            unsafe {
                speex_preprocess_ctl(
                    self.st,
                    cmd as std::os::raw::c_int,
                    &mut val as *mut _ as *mut c_void,
                );
            }
            self
        }

        fn get_i32(&mut self, cmd: u32) -> i32 {
            let mut v: i32 = 0;
            unsafe {
                speex_preprocess_ctl(
                    self.st,
                    cmd as std::os::raw::c_int,
                    &mut v as *mut _ as *mut c_void,
                );
            }
            v
        }

        pub fn set_denoise(&mut self, enable: bool) -> &mut Self {
            self.set_bool(SPEEX_PREPROCESS_SET_DENOISE, enable)
        }
        pub fn get_denoise(&mut self) -> bool {
            self.get_bool(SPEEX_PREPROCESS_GET_DENOISE)
        }

        pub fn set_agc(&mut self, enable: bool) -> &mut Self {
            self.set_bool(SPEEX_PREPROCESS_SET_AGC, enable)
        }
        pub fn get_agc(&mut self) -> bool {
            self.get_bool(SPEEX_PREPROCESS_GET_AGC)
        }

        pub fn set_vad(&mut self, enable: bool) -> &mut Self {
            self.set_bool(SPEEX_PREPROCESS_SET_VAD, enable)
        }
        pub fn get_vad(&mut self) -> bool {
            self.get_bool(SPEEX_PREPROCESS_GET_VAD)
        }

        pub fn set_agc_level(&mut self, level: f32) -> &mut Self {
            self.set_f32(SPEEX_PREPROCESS_SET_AGC_LEVEL, level)
        }
        pub fn get_agc_level(&mut self) -> f32 {
            self.get_f32(SPEEX_PREPROCESS_GET_AGC_LEVEL)
        }

        pub fn set_dereverb(&mut self, enable: bool) -> &mut Self {
            self.set_bool(SPEEX_PREPROCESS_SET_DEREVERB, enable)
        }
        pub fn get_dereverb(&mut self) -> bool {
            self.get_bool(SPEEX_PREPROCESS_GET_DEREVERB)
        }

        pub fn set_dereverb_level(&mut self, level: f32) -> &mut Self {
            self.set_f32(SPEEX_PREPROCESS_SET_DEREVERB_LEVEL, level)
        }
        pub fn get_dereverb_level(&mut self) -> f32 {
            self.get_f32(SPEEX_PREPROCESS_GET_DEREVERB_LEVEL)
        }

        pub fn set_dereverb_decay(&mut self, decay: f32) -> &mut Self {
            self.set_f32(SPEEX_PREPROCESS_SET_DEREVERB_DECAY, decay)
        }
        pub fn get_dereverb_decay(&mut self) -> f32 {
            self.get_f32(SPEEX_PREPROCESS_GET_DEREVERB_DECAY)
        }

        pub fn set_prob_start(&mut self, val: i32) -> &mut Self {
            self.set_i32(SPEEX_PREPROCESS_SET_PROB_START, val)
        }
        pub fn get_prob_start(&mut self) -> i32 {
            self.get_i32(SPEEX_PREPROCESS_GET_PROB_START)
        }

        pub fn set_prob_continue(&mut self, val: i32) -> &mut Self {
            self.set_i32(SPEEX_PREPROCESS_SET_PROB_CONTINUE, val)
        }
        pub fn get_prob_continue(&mut self) -> i32 {
            self.get_i32(SPEEX_PREPROCESS_GET_PROB_CONTINUE)
        }

        pub fn set_noise_suppress(&mut self, val: i32) -> &mut Self {
            self.set_i32(SPEEX_PREPROCESS_SET_NOISE_SUPPRESS, val)
        }
        pub fn get_noise_suppress(&mut self) -> i32 {
            self.get_i32(SPEEX_PREPROCESS_GET_NOISE_SUPPRESS)
        }

        pub fn set_echo_suppress(&mut self, val: i32) -> &mut Self {
            self.set_i32(SPEEX_PREPROCESS_SET_ECHO_SUPPRESS, val)
        }
        pub fn get_echo_suppress(&mut self) -> i32 {
            self.get_i32(SPEEX_PREPROCESS_GET_ECHO_SUPPRESS)
        }

        pub fn set_echo_suppress_active(&mut self, val: i32) -> &mut Self {
            self.set_i32(SPEEX_PREPROCESS_SET_ECHO_SUPPRESS_ACTIVE, val)
        }
        pub fn get_echo_suppress_active(&mut self) -> i32 {
            self.get_i32(SPEEX_PREPROCESS_GET_ECHO_SUPPRESS_ACTIVE)
        }

        pub fn set_echo_state(
            &mut self,
            state: Option<&SpeexEcho>,
        ) -> &mut Self {
            let ptr = match state {
                Some(echo) => echo.st as *mut c_void,
                None => std::ptr::null_mut(),
            };

            unsafe {
                speex_preprocess_ctl(
                    self.st,
                    SPEEX_PREPROCESS_SET_ECHO_STATE as std::os::raw::c_int,
                    ptr,
                );
            }
            self
        }

        pub fn get_echo_state(&mut self) -> Option<SpeexEcho> {
            let mut ptr: *mut SpeexEchoState = std::ptr::null_mut();

            unsafe {
                speex_preprocess_ctl(
                    self.st,
                    SPEEX_PREPROCESS_GET_ECHO_STATE as std::os::raw::c_int,
                    &mut ptr as *mut _ as *mut c_void,
                );
            }

            if ptr.is_null() {
                None
            } else {
                Some(SpeexEcho { st: ptr })
            }
        }

        pub fn set_agc_increment(&mut self, val: i32) -> &mut Self {
            self.set_i32(SPEEX_PREPROCESS_SET_AGC_INCREMENT, val)
        }
        pub fn get_agc_increment(&mut self) -> i32 {
            self.get_i32(SPEEX_PREPROCESS_GET_AGC_INCREMENT)
        }

        pub fn set_agc_decrement(&mut self, val: i32) -> &mut Self {
            self.set_i32(SPEEX_PREPROCESS_SET_AGC_DECREMENT, val)
        }
        pub fn get_agc_decrement(&mut self) -> i32 {
            self.get_i32(SPEEX_PREPROCESS_GET_AGC_DECREMENT)
        }

        pub fn set_agc_max_gain(&mut self, val: i32) -> &mut Self {
            self.set_i32(SPEEX_PREPROCESS_SET_AGC_MAX_GAIN, val)
        }
        pub fn get_agc_max_gain(&mut self) -> i32 {
            self.get_i32(SPEEX_PREPROCESS_GET_AGC_MAX_GAIN)
        }

        pub fn get_agc_loudness(&mut self) -> i32 {
            self.get_i32(SPEEX_PREPROCESS_GET_AGC_LOUDNESS)
        }

        pub fn get_agc_gain(&mut self) -> i32 {
            self.get_i32(SPEEX_PREPROCESS_GET_AGC_GAIN)
        }

        pub fn get_psd_size(&mut self) -> i32 {
            self.get_i32(SPEEX_PREPROCESS_GET_PSD_SIZE)
        }

        // todo: SPEEX_PREPROCESS_GET_PSD

        pub fn get_noise_psd_size(&mut self) -> i32 {
            self.get_i32(SPEEX_PREPROCESS_GET_NOISE_PSD_SIZE)
        }

        // todo: SPEEX_PREPROCESS_GET_NOISE_PSD

        pub fn get_prob(&mut self) -> i32 {
            self.get_i32(SPEEX_PREPROCESS_GET_PROB)
        }

        pub fn set_agc_target(&mut self, val: i32) -> &mut Self {
            self.set_i32(SPEEX_PREPROCESS_SET_AGC_TARGET, val)
        }
        pub fn get_agc_target(&mut self) -> i32 {
            self.get_i32(SPEEX_PREPROCESS_GET_AGC_TARGET)
        }
    }
}

#[cfg(feature = "sys")]
pub use self::sys::{Error, SpeexPreprocess};
