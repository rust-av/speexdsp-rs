#[cfg(feature = "sys")]
extern crate byteorder;
extern crate speexdsp;

#[cfg(feature = "sys")]
fn main() {
    use byteorder::{BigEndian, ByteOrder};
    use speexdsp::preprocess::SpeexPreprocessConst::*;
    use speexdsp::preprocess::*;
    use std::io::Read;

    const NN: usize = 160;
    let mut input: [i16; NN] = [0; NN];
    let mut buffer: [u8; NN * 2] = [0; NN * 2];

    let mut st = SpeexPreprocess::new(NN, 8000).unwrap();
    st.preprocess_ctl(SPEEX_PREPROCESS_SET_DENOISE, 1f32).unwrap();
    st.preprocess_ctl(SPEEX_PREPROCESS_SET_AGC, 0f32).unwrap();
    st.preprocess_ctl(SPEEX_PREPROCESS_SET_AGC_LEVEL, 8000f32).unwrap();
    st.preprocess_ctl(SPEEX_PREPROCESS_SET_DEREVERB, 0f32).unwrap();
    st.preprocess_ctl(SPEEX_PREPROCESS_SET_DEREVERB_DECAY, 0f32).unwrap();
    st.preprocess_ctl(SPEEX_PREPROCESS_SET_DEREVERB_LEVEL, 0f32).unwrap();

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
