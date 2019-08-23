extern crate assert_approx_eq;
extern crate interpolate_name;
extern crate speexdsp;

#[cfg(feature = "sys")]
mod comparison {
    use assert_approx_eq::assert_approx_eq;
    use interpolate_name::interpolate_test;

    use speexdsp::resampler::*;
    use speexdsp::speex_resample::*;

    use std::f32::consts::PI;

    const PERIOD: f32 = 32f32;
    const INBLOCK: usize = 1024;
    const RATE: usize = 48000;

    #[inline(always)]
    fn process_float_native(
        st: &mut SpeexResamplerState,
        index: usize,
        input: &[f32],
        output: &mut [f32],
    ) -> (usize, usize) {
        let mut in_len = input.len() as u32;
        let mut out_len = output.len() as u32;
        st.process_float(
            index as u32,
            input,
            &mut in_len,
            output,
            &mut out_len,
        );

        (in_len as usize, out_len as usize)
    }

    #[interpolate_test(8, 8)]
    #[interpolate_test(10, 10)]
    fn resampling(quality: usize) {
        let start_rate = 1000;
        let mut off = 0;
        let mut avail = INBLOCK as isize;

        let fin: Vec<f32> = (0..INBLOCK * 4)
            .map(|i| ((i as f32) / PERIOD * 2.0 * PI).sin() * 0.9)
            .collect();
        let mut fout = vec![0f32; INBLOCK * 8];
        let mut fout_native = vec![0f32; INBLOCK * 8];

        let mut st = State::new(1, RATE, RATE, 4).unwrap();

        let mut st_native = SpeexResamplerState::new(1, RATE, RATE, 4);

        st.set_rate(RATE, start_rate);
        st.skip_zeros();

        st.set_quality(quality).unwrap();

        st_native.set_rate(RATE, start_rate);
        st_native.skip_zeros();

        st_native.set_quality(quality);

        for rate in (start_rate..128000).step_by(5000) {
            let in_len = avail as usize;
            let out_len = (in_len * rate + RATE - 1) / RATE;
            let prev_in_len = in_len;
            let prev_out_len = out_len;

            let (in_len_native, out_len_native) = process_float_native(
                &mut st_native,
                0,
                &fin[off..off + in_len],
                &mut fout_native[..out_len],
            );

            let (in_len, out_len) = st
                .process_float(
                    0,
                    &fin[off..off + in_len],
                    &mut fout[..out_len],
                )
                .unwrap();

            eprintln!(
                "{} {} {} {} -> {} {}",
                rate, off, prev_in_len, prev_out_len, in_len, out_len
            );

            assert_eq!(in_len, in_len_native);
            assert_eq!(out_len, out_len_native);

            off += in_len as usize;
            avail += INBLOCK as isize - in_len as isize;

            if off >= INBLOCK {
                off -= INBLOCK;
            }

            let fout_s = &fout[..out_len as usize];
            let fout_native_s = &fout_native[..out_len as usize];

            fout_s
                .iter()
                .zip(fout_native_s.iter())
                .for_each(|(&x, &y)| {
                    assert_approx_eq!(x, y, 1.0e-6);
                });

            st.set_rate(RATE, rate);
            st_native.set_rate(RATE, rate);
        }
    }
}
