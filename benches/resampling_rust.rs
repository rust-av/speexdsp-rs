use std::f32::consts::PI;

use criterion::{criterion_group, criterion_main, Criterion};
use speexdsp::resampler::*;

const PERIOD: f32 = 32f32;
const INBLOCK: usize = 1024;
const RATE: usize = 48000;

fn resample_rs() {
    let mut rate = 1000;
    let mut off = 0;
    let mut avail = INBLOCK as isize;

    let fin: Vec<f32> = (0..INBLOCK * 4)
        .map(|i| ((i as f32) / PERIOD * 2.0 * PI).sin() * 0.9)
        .collect();
    let mut fout = vec![0f32; INBLOCK * 8];

    let mut st = State::new(1, RATE, RATE, 8).unwrap();

    st.set_rate(RATE, rate).unwrap();
    st.skip_zeros();

    let mut data = Vec::new();

    loop {
        let in_len = avail as usize;
        let out_len = (in_len * rate + RATE - 1) / RATE;

        let (in_len, out_len) = st
            .process_float(0, &fin[off..off + in_len], &mut fout[..out_len])
            .unwrap();

        off += in_len as usize;
        avail += INBLOCK as isize - in_len as isize;

        if off >= INBLOCK {
            off -= INBLOCK;
        }

        data.push(fout[..out_len as usize].to_vec());

        rate += 5000;
        if rate > 128000 {
            break;
        }

        st.set_rate(RATE, rate).unwrap();
    }
}

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("resampler_rust", |b| b.iter(|| resample_rs()));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
