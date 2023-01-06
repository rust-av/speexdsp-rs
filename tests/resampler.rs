#[cfg(feature = "sys")]
mod comparison {
    use assert_approx_eq::assert_approx_eq;
    use interpolate_name::interpolate_test;

    use speexdsp::resampler as sys;
    use speexdsp::resampler::native;
    use speexdsp::resampler::Resampler;

    use std::f32::consts::PI;

    const PERIOD: f32 = 32f32;
    const INBLOCK: usize = 1024;
    const RATE: usize = 48000;

    #[interpolate_test(quality8_num_gt_den, 8, true)]
    #[interpolate_test(quality10_num_gt_den, 10, true)]
    #[interpolate_test(quality8_num_lt_den, 8, false)]
    #[interpolate_test(quality10_num_lt_den, 10, false)]
    fn resampling(quality: usize, num_gt_den: bool) {
        let mut init_rate = RATE;
        let mut start_rate = 1000;
        let mut off = 0;
        let mut avail = INBLOCK as isize;

        let fin: Vec<f32> = (0..INBLOCK * 4)
            .map(|i| ((i as f32) / PERIOD * 2.0 * PI).sin() * 0.9)
            .collect();
        let mut fout = vec![0f32; INBLOCK * 8];
        let mut fout_native = vec![0f32; INBLOCK * 8];

        if !num_gt_den {
            init_rate = 96000;
            start_rate = 96000;
        }

        let mut st = sys::State::new(1, RATE, init_rate, 4).unwrap();

        let mut st_native = native::State::new(1, RATE, init_rate, 4).unwrap();

        st.set_rate(RATE, start_rate).unwrap();
        st.skip_zeros();

        st.set_quality(quality).unwrap();

        st_native.set_rate(RATE, start_rate).unwrap();
        st_native.skip_zeros();

        st_native.set_quality(quality).unwrap();

        for rate in (start_rate..128000).step_by(5000) {
            let in_len = avail as usize;
            let out_len = (in_len * rate + RATE - 1) / RATE;
            let prev_in_len = in_len;
            let prev_out_len = out_len;

            let (in_len_native, out_len_native) = st_native
                .process_float(
                    0,
                    &fin[off..off + in_len],
                    &mut fout_native[..out_len],
                )
                .unwrap();

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

            off += in_len;
            avail += INBLOCK as isize - in_len as isize;

            if off >= INBLOCK {
                off -= INBLOCK;
            }

            let fout_s = &fout[..out_len];
            let fout_native_s = &fout_native[..out_len];

            fout_s
                .iter()
                .zip(fout_native_s.iter())
                .for_each(|(&x, &y)| {
                    assert_approx_eq!(x, y, 1.0e-6);
                });

            if num_gt_den {
                st.set_rate(RATE, rate).unwrap();
                st_native.set_rate(RATE, rate).unwrap();
            } else {
                st.set_rate(rate, RATE).unwrap();
                st_native.set_rate(rate, RATE).unwrap();
            }
        }
    }
}
