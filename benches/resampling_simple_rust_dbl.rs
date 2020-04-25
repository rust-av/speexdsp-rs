use std::f32::consts::PI;

use criterion::{criterion_group, criterion_main, Criterion};
use speexdsp::resampler::*;

const PERIOD: f32 = 32f32;
const RATE: usize = 48000;
const INBLOCK: usize = 1024 * 4;

fn resample_rs() {
    let fin: Vec<f32> = (0..INBLOCK * 4)
        .map(|i| ((i as f32) / PERIOD * 2.0 * PI).sin() * 0.9)
        .collect();
    let mut fout = vec![0f32; INBLOCK * 4 * 4];

    let mut st = State::new(1, RATE, RATE * 4, 10).unwrap();

    st.skip_zeros();

    let (_in_len, _out_len) = st.process_float(0, &fin, &mut fout).unwrap();
}

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("resampler_simple_rust_dbl", |b| {
        b.iter(|| resample_rs())
    });
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = criterion_benchmark
}
criterion_main!(benches);
