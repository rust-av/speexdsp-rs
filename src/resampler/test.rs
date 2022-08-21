use crate::resampler::*;
use std::f32::consts::PI;

const PERIOD: f32 = 32f32;
const INBLOCK: usize = 1024;
const RATE: usize = 48000;

fn make_input<T>(block: usize, period: f32, channels: usize) -> Vec<T>
where
    T: From<f32> + Clone,
{
    (0..block * 2)
        .map(|i| {
            let v = (((i as f32) / period * 2.0 * PI).sin() * 0.9).into();
            std::iter::repeat(v).take(channels)
        })
        .flatten()
        .collect()
}

#[test]
fn interleaved() {
    let mut rate = 1000;
    let mut off = 0;
    let mut avail = INBLOCK as isize;
    let channels = 2;

    let fin = make_input(INBLOCK, PERIOD, channels);

    let mut fout = vec![0f32; INBLOCK * channels * 4];

    let mut st = State::new(channels, RATE, RATE, 4).unwrap();
    st.set_rate(RATE, rate);
    st.skip_zeros();

    loop {
        let mut in_len = avail as usize;
        let mut out_len = (in_len * rate + RATE - 1) / RATE;

        let prev_in_len = in_len;
        let prev_out_len = out_len;

        eprintln!("in_len {}, out_len {}", in_len, out_len);

        (in_len, out_len) = st
            .process_interleaved_float(
                &fin[off * channels..(off + in_len) * channels],
                &mut fout[..out_len * channels],
            )
            .unwrap();

        in_len /= channels;
        out_len /= channels;

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

        st.set_rate(RATE, rate);
    }
}
