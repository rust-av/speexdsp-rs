extern crate speexdsp;
#[macro_use]
extern crate proptest;

use proptest::collection::size_range;
use proptest::prelude::*;

use speexdsp::resampler::*;
use std::f32::consts::PI;

const PERIOD: f32 = 32f32;
const INBLOCK: usize = 1024;
const RATE: usize = 48000;

proptest! {
    #[test]
    fn process_resampler(quality in 0usize..10, ref fin in any_with::<Vec<f32>>(size_range(INBLOCK * 2).lift())) {
        resampler(quality, fin);
    }

}

fn resampler(quality: usize, fin: &Vec<f32>) {
    let mut rate = 1000;
    let mut off = 0;
    let mut avail = INBLOCK as isize;

    let fin: Vec<f32> = fin
        .iter()
        .map(|f| (f / PERIOD * 2.0 * PI).sin() * 0.9)
        .collect();
    let mut fout = vec![0f32; INBLOCK * 4];

    let mut st = State::new(1, RATE, RATE, 4).unwrap();

    st.set_rate(RATE, rate);
    st.skip_zeros();

    st.set_quality(quality).unwrap();

    loop {
        let in_len = avail as usize;
        let out_len = (in_len * rate + RATE - 1) / RATE;

        let (in_len, _out_len) = st
            .process_float(0, &fin[off..off + in_len], &mut fout[..out_len])
            .unwrap();

        off += in_len as usize;
        avail += INBLOCK as isize - in_len as isize;

        if off >= INBLOCK {
            off -= INBLOCK;
        }

        rate += 100;
        if rate > 128000 {
            break;
        }

        st.set_rate(RATE, rate);
    }
}
