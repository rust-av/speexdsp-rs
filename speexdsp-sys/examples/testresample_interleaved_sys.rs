use speexdsp_sys::resampler::*;
use std::f32::consts::PI;
use std::ptr;

const PERIOD: f32 = 32f32;
const INBLOCK: usize = 1024;
const RATE: u32 = 48000;

fn main() {
    let mut rate = 1000;
    let mut off = 0;
    let mut avail = INBLOCK as isize;
    let channels = 2;

    let fin: Vec<f32> = (0..INBLOCK * 2)
        .map(|i| {
            let v = ((i as f32) / PERIOD * 2.0 * PI).sin() * 0.9;
            std::iter::repeat(v).take(channels)
        })
        .flatten()
        .collect();
    let mut fout = vec![0f32; INBLOCK * channels * 4];

    let st = unsafe {
        speex_resampler_init(channels as u32, RATE, RATE, 4, ptr::null_mut())
    };
    unsafe { speex_resampler_set_rate(st, RATE, rate) };
    unsafe { speex_resampler_skip_zeros(st) };

    loop {
        let mut in_len = avail as u32;
        let mut out_len = (in_len * rate + RATE - 1) / RATE;

        let prev_in_len = in_len;
        let prev_out_len = out_len;

        unsafe {
            speex_resampler_process_interleaved_float(
                st,
                fin[off * channels..].as_ptr(),
                &mut in_len,
                fout.as_mut_ptr(),
                &mut out_len,
            )
        };

        eprintln!(
            "{} {} {} {} -> {} {}",
            rate, off, prev_in_len, prev_out_len, in_len, out_len
        );

        off += in_len as usize;
        avail += INBLOCK as isize - in_len as isize;

        if off >= INBLOCK {
            off -= INBLOCK;
        }

        assert_eq!(fout[0], fout[1]);

        for v in &fout[..out_len as usize * channels] {
            println!("{:.05}", v);
        }

        rate += 5000;
        if rate > 128000 {
            break;
        }

        unsafe { speex_resampler_set_rate(st, RATE, rate) };
    }

    unsafe { speex_resampler_destroy(st) };
}
