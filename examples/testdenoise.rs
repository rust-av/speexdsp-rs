#![allow(dead_code)]
#![allow(unused_imports)]
extern crate byteorder;
extern crate speexdsp;

use byteorder::{BigEndian, ByteOrder};
use std::io::Read;

#[cfg(feature = "sys")]
fn main() {
    use speexdsp::preprocess::*;

    const NN: usize = 160;
    let mut input: [i16; NN] = [0; NN];
    let mut buffer: [u8; NN * 2] = [0; NN * 2];

    let mut st = SpeexPreprocess::new(NN, 8000).unwrap();

    st.set_denoise(true)
        .set_agc(false)
        .set_agc_level(8000f32)
        .set_dereverb(false)
        .set_dereverb_decay(0f32)
        .set_dereverb_level(0f32);

    while let Ok(n) = std::io::stdin().read(&mut buffer) {
        if n == 0 {
            break;
        }
        BigEndian::read_i16_into(&buffer, &mut input);
        st.preprocess_run(&mut input);
        println!("{:?}", &input[..]);
    }
}

#[cfg(not(feature = "sys"))]
fn main() {
    unimplemented!();
}
