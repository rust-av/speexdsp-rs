extern crate byteorder;
extern crate speexdsp_sys;

use byteorder::{BigEndian, ByteOrder};
use speexdsp_sys::preprocess::*;
use std::ffi::c_void;
use std::io::Read;

macro_rules! preprocess_ctl {
    ($st:ident, $param:ident, $value:expr) => {
        let v = &mut ($value as f32) as *mut f32 as *mut c_void;
        unsafe { speex_preprocess_ctl($st, $param as i32, v) };
    };
}

const NN: usize = 160;

fn main() {
    let mut input: [i16; NN] = [0; NN];
    let mut buffer: [u8; NN * 2] = [0; NN * 2];

    let st = unsafe { speex_preprocess_state_init(NN as i32, 8000) };
    preprocess_ctl!(st, SPEEX_PREPROCESS_SET_DENOISE, 1);
    preprocess_ctl!(st, SPEEX_PREPROCESS_SET_AGC, 0);
    preprocess_ctl!(st, SPEEX_PREPROCESS_SET_AGC_LEVEL, 8000);
    preprocess_ctl!(st, SPEEX_PREPROCESS_SET_DEREVERB, 0);
    preprocess_ctl!(st, SPEEX_PREPROCESS_SET_DEREVERB_DECAY, 0f32);
    preprocess_ctl!(st, SPEEX_PREPROCESS_SET_DEREVERB_LEVEL, 0f32);

    while let Ok(n) = std::io::stdin().read(&mut buffer) {
        if n == 0 {
            break;
        }
        BigEndian::read_i16_into(&buffer, &mut input);
        unsafe { speex_preprocess_run(st, input.as_mut_ptr()) };
        println!("{:?}", &input[..]);
    }

    unsafe { speex_preprocess_state_destroy(st) };
}
