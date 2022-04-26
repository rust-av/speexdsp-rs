use criterion::{criterion_group, criterion_main, Criterion};

#[cfg(feature = "sys")]
use speexdsp_sys::resampler::*;

#[cfg(feature = "sys")]
fn resample_c() {
    use std::f32::consts::PI;
    use std::ptr;

    const PERIOD: f32 = 32f32;
    const RATE: usize = 48000;
    const INBLOCK: usize = 1024 * 4;

    let fin: Vec<f32> = (0..INBLOCK * 4)
        .map(|i| ((i as f32) / PERIOD * 2.0 * PI).sin() * 0.9)
        .collect();
    let mut fout = vec![0f32; INBLOCK * 4 * 4];

    let mut in_len = (INBLOCK * 4) as u32;
    let mut out_len = (INBLOCK * 4 * 4) as u32;

    let st = unsafe {
        speex_resampler_init(
            1,
            RATE as u32,
            (RATE * 4) as u32,
            10,
            ptr::null_mut(),
        )
    };
    unsafe { speex_resampler_skip_zeros(st) };

    unsafe {
        speex_resampler_process_float(
            st,
            0,
            fin.as_ptr(),
            &mut in_len,
            fout.as_mut_ptr(),
            &mut out_len,
        )
    };
}

#[cfg(not(feature = "sys"))]
fn resample_c() {}

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("resampler_simple_c_dbl", |b| b.iter(resample_c));
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = criterion_benchmark
}
criterion_main!(benches);
