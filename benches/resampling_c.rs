#[cfg(feature = "sys")]
extern crate speexdsp_sys;

use criterion::{criterion_group, criterion_main, Criterion};
#[cfg(feature = "sys")]
use speexdsp_sys::resampler::*;

#[cfg(feature = "sys")]
fn resample_c() {
    use std::f32::consts::PI;
    use std::ptr;

    const PERIOD: f32 = 32f32;
    const INBLOCK: usize = 1024;
    const RATE: u32 = 48000;

    let mut rate = 1000;
    let mut off = 0;
    let mut avail = INBLOCK as isize;

    let fin: Vec<f32> = (0..INBLOCK * 4)
        .map(|i| ((i as f32) / PERIOD * 2.0 * PI).sin() * 0.9)
        .collect();
    let mut fout = vec![0f32; INBLOCK * 4];

    let st =
        unsafe { speex_resampler_init(1, RATE, RATE, 8, ptr::null_mut()) };
    unsafe { speex_resampler_set_rate(st, RATE, rate) };
    unsafe { speex_resampler_skip_zeros(st) };

    loop {
        let mut in_len = avail as u32;
        let mut out_len = (in_len * rate + RATE - 1) / RATE;

        unsafe {
            speex_resampler_process_float(
                st,
                0,
                fin[off..].as_ptr(),
                &mut in_len,
                fout.as_mut_ptr(),
                &mut out_len,
            )
        };

        off += in_len as usize;
        avail += INBLOCK as isize - in_len as isize;

        if off >= INBLOCK {
            off -= INBLOCK;
        }

        rate += 5000;
        if rate > 128000 {
            break;
        }

        unsafe { speex_resampler_set_rate(st, RATE, rate) };
    }

    unsafe { speex_resampler_destroy(st) };
}

#[cfg(not(feature = "sys"))]
fn resample_c() {}

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("resampler_c", |b| b.iter(|| resample_c()));
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = criterion_benchmark
}
criterion_main!(benches);
