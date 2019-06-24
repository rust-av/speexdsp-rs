#[cfg(feature = "sys")]
mod sys {
    use speexdsp_sys::jitter::*;
    use std::ffi::c_void;
    use std::fmt;

    #[derive(Clone, Copy, Debug)]
    pub enum SpeexJitterConst {
        JITTER_BUFFER_SET_MARGIN = 0,
        JITTER_BUFFER_GET_MARGIN = 1,
        JITTER_BUFFER_GET_AVAILABLE_COUNT = 3,
        JITTER_BUFFER_SET_DESTROY_CALLBACK = 4,
        JITTER_BUFFER_GET_DESTROY_CALLBACK = 5,
        JITTER_BUFFER_SET_DELAY_STEP = 6,
        JITTER_BUFFER_GET_DELAY_STEP = 7,
        JITTER_BUFFER_SET_CONCEALMENT_SIZE = 8,
        JITTER_BUFFER_GET_CONCEALMENT_SIZE = 9,
        JITTER_BUFFER_SET_MAX_LATE_RATE = 10,
        JITTER_BUFFER_GET_MAX_LATE_RATE = 11,
        JITTER_BUFFER_SET_LATE_COST = 12,
        JITTER_BUFFER_GET_LATE_COST = 13,
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

    pub struct SpeexBufferPacket {
        pt: *mut JitterBufferPacket,
    }

    impl SpeexBufferPacket {
        pub fn new(
            data: &mut [i8],
            len: usize,
            timestamp: usize,
            span: usize,
            sequence: usize,
            user_data: usize,
        ) -> Self {
            let mut pt = JitterBufferPacket {
                data: data.as_mut_ptr(),
                len: len as u32,
                timestamp: timestamp as u32,
                span: span as u32,
                sequence: sequence as u16,
                user_data: user_data as u32,
            };
            SpeexBufferPacket {
                pt: &mut pt as *mut JitterBufferPacket,
            }
        }
    }

    pub struct SpeexJitter {
        st: *mut JitterBuffer,
    }

    impl SpeexJitter {
        pub fn new(step_size: usize) -> Result<Self, Error> {
            let st = unsafe { jitter_buffer_init(step_size as i32) };

            if st.is_null() {
                Err(Error::FailedInit)
            } else {
                Ok(SpeexJitter { st })
            }
        }

        pub fn buffer_reset(&mut self) {
            unsafe { jitter_buffer_reset(self.st) };
        }

        pub fn buffer_put(&mut self, packet: SpeexBufferPacket) {
            unsafe { jitter_buffer_put(self.st, packet.pt) };
        }

        pub fn buffer_get(
            &mut self,
            packet: SpeexBufferPacket,
            desired_span: usize,
            offset: usize,
        ) -> usize {
            unsafe {
                jitter_buffer_get(self.st, packet.pt, desired_span as i32, offset as *mut i32)
                    as usize
            }
        }

        pub fn buffer_get_another(&mut self, packet: SpeexBufferPacket) -> usize {
            unsafe { jitter_buffer_get_another(self.st, packet.pt) as usize }
        }

        pub fn buffer_get_pointer_timestap(&mut self) -> usize {
            unsafe { jitter_buffer_get_pointer_timestamp(self.st) as usize }
        }

        pub fn buffer_tick(&mut self) {
            unsafe { jitter_buffer_tick(self.st) };
        }

        pub fn buffer_remaining_span(&mut self, rem: usize) {
            unsafe { jitter_buffer_remaining_span(self.st, rem as u32) };
        }

        pub fn buffer_ctl(&mut self, request: SpeexJitterConst, ptr: usize) -> Result<(), Error> {
            let ret =
                unsafe { jitter_buffer_ctl(self.st, request as i32, ptr as *mut c_void) as usize };
            if ret != 0 {
                Err(Error::UnknownRequest)
            } else {
                Ok(())
            }
        }

        pub fn buffer_update_delay(
            &mut self,
            packet: SpeexBufferPacket,
            start_offset: usize,
        ) -> usize {
            unsafe {
                jitter_buffer_update_delay(self.st, packet.pt, start_offset as *mut i32) as usize
            }
        }
    }

    impl Drop for SpeexJitter {
        fn drop(&mut self) {
            unsafe { jitter_buffer_destroy(self.st) };
        }
    }
}
#[cfg(feature = "sys")]
pub use self::sys::{Error, SpeexBufferPacket, SpeexJitter, SpeexJitterConst};
