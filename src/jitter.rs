#![allow(non_camel_case_types)]

#[cfg(feature = "sys")]
mod sys {
    use speexdsp_sys::jitter::*;
    use std::ffi::c_void;
    use std::fmt;

    const BUFFER_OK: i32 = JITTER_BUFFER_OK as i32;
    const BUFFER_MISSING: i32 = JITTER_BUFFER_MISSING as i32;
    const BUFFER_INSERTION: i32 = JITTER_BUFFER_INSERTION as i32;

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

    #[derive(Clone, Copy, Debug, PartialEq)]
    pub enum Error {
        BufferOk = 0,
        BufferMissing = 1,
        BufferInsertion = 2,
        BufferInternalError = -1,
        BufferBadArgument = -2,
    }

    impl fmt::Display for Error {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            let v = match self {
                Error::BufferOk => "The buffer is ok",
                Error::BufferMissing => "The buffer is missing",
                Error::BufferInsertion => "Bad insertion into the buffer",
                Error::BufferInternalError => "Internal Error",
                Error::BufferBadArgument => "Bad Argument...",
            };

            write!(f, "{}", v)
        }
    }

    impl Default for SpeexBufferPacket {
        fn default() -> Self {
            Self::new()
        }
    }

    impl Error {
        fn from_i32(v: i32) -> Self {
            match v {
                BUFFER_OK => Error::BufferOk,
                BUFFER_MISSING => Error::BufferMissing,
                BUFFER_INSERTION => Error::BufferInsertion,
                JITTER_BUFFER_INTERNAL_ERROR => Error::BufferInternalError,
                JITTER_BUFFER_BAD_ARGUMENT => Error::BufferBadArgument,
                _ => unreachable!(),
            }
        }
    }

    pub struct SpeexBufferPacket {
        pt: JitterBufferPacket,
    }

    impl SpeexBufferPacket {
        pub fn new() -> Self {
            let pt = JitterBufferPacket {
                data: std::ptr::null_mut(),
                len: 0,
                timestamp: 0,
                span: 0,
                sequence: 0,
                user_data: 0,
            };

            SpeexBufferPacket { pt }
        }

        pub fn create(
            &mut self,
            data: &mut [i8],
            len: usize,
            timestamp: usize,
            span: usize,
            sequence: usize,
            user_data: usize,
        ) {
            let pt = JitterBufferPacket {
                data: data.as_mut_ptr(),
                len: len as u32,
                timestamp: timestamp as u32,
                span: span as u32,
                sequence: sequence as u16,
                user_data: user_data as u32,
            };
            self.pt = pt;
        }

        pub fn len(&self) -> usize {
            self.pt.len as usize
        }

        pub fn is_empty(&self) -> bool {
            self.pt.len != 0
        }

        pub fn timestamp(&self) -> usize {
            self.pt.timestamp as usize
        }

        pub fn span(&self) -> usize {
            self.pt.span as usize
        }

        pub fn sequence(&self) -> usize {
            self.pt.sequence as usize
        }

        pub fn user_data(&self) -> usize {
            self.pt.user_data as usize
        }

        pub fn set_data(&mut self, data: &mut [i8]) {
            self.pt.data = data.as_mut_ptr();
        }

        pub fn set_len(&mut self, len: usize) {
            self.pt.len = len as u32;
        }
    }

    pub struct SpeexJitter {
        st: *mut JitterBuffer,
    }

    impl SpeexJitter {
        pub fn new(step_size: usize) -> Result<Self, Error> {
            let st = unsafe { jitter_buffer_init(step_size as i32) };

            if st.is_null() {
                Err(Error::BufferInternalError)
            } else {
                Ok(SpeexJitter { st })
            }
        }

        pub fn buffer_reset(&mut self) {
            unsafe { jitter_buffer_reset(self.st) };
        }

        pub fn buffer_put(&mut self, packet: &SpeexBufferPacket) {
            unsafe { jitter_buffer_put(self.st, &packet.pt) };
        }

        pub fn buffer_get(
            &mut self,
            packet: &mut SpeexBufferPacket,
            desired_span: usize,
            offset: usize,
        ) -> Error {
            let err_i32 = unsafe {
                jitter_buffer_get(
                    self.st,
                    &mut packet.pt,
                    desired_span as i32,
                    offset as *mut i32,
                )
            };
            Error::from_i32(err_i32)
        }

        pub fn buffer_get_another(
            &mut self,
            packet: &mut SpeexBufferPacket,
        ) -> Error {
            let err_i32 =
                unsafe { jitter_buffer_get_another(self.st, &mut packet.pt) };
            Error::from_i32(err_i32)
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

        pub fn buffer_ctl(
            &mut self,
            request: SpeexJitterConst,
            ptr: usize,
        ) -> Result<(), Error> {
            let ret = unsafe {
                jitter_buffer_ctl(self.st, request as i32, ptr as *mut c_void)
                    as usize
            };
            if ret != 0 {
                Err(Error::BufferBadArgument)
            } else {
                Ok(())
            }
        }

        pub fn buffer_update_delay(
            &mut self,
            packet: &mut SpeexBufferPacket,
            start_offset: usize,
        ) -> Error {
            let err_i32 = unsafe {
                jitter_buffer_update_delay(
                    self.st,
                    &mut packet.pt,
                    start_offset as *mut i32,
                )
            };
            Error::from_i32(err_i32)
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
