use speexdsp_resampler::*;

use std::f32::consts::PI;

const PERIOD: f32 = 32f32;
const INBLOCK: usize = 1024;
const RATE: usize = 48000;

fn main() {
    let mut rate = 1000;
    let mut off = 0;
    let mut avail = INBLOCK as isize;

    let fin: Vec<f32> = (0..INBLOCK * 4)
        .map(|i| ((i as f32) / PERIOD * 2.0 * PI).sin() * 0.9)
        .collect();
    let mut fout = vec![0f32; INBLOCK * 8];

    let mut st = State::new(1, RATE, RATE, 4).unwrap();

    st.set_rate(RATE, rate).unwrap();
    st.skip_zeros();

    st.set_quality(10).unwrap();

    eprintln!("Quality: {}", st.get_quality());

    let mut data = Vec::new();

    loop {
        let in_len = avail as usize;
        let out_len = (in_len * rate + RATE - 1) / RATE;

        let prev_in_len = in_len;
        let prev_out_len = out_len;

        let (in_len, out_len) = st
            .process_float(0, &fin[off..off + in_len], &mut fout[..out_len])
            .unwrap();

        eprintln!(
            "{} {} {} {} -> {} {}",
            rate, off, prev_in_len, prev_out_len, in_len, out_len
        );

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

    println!("{:#?}", data);
}
