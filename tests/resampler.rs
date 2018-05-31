extern crate speexdsp;

#[cfg(test)]
mod tests {
    use speexdsp::resampler::State;
    use std::f32::consts::PI;

    const RATE: usize = 48000;
    const INBLOCK: usize = 1024;
    const PERIOD: f32 = 32f32;

    #[test]
    fn test_rate() {
        let mut st = State::new(1, RATE, RATE, 4).unwrap();
        let expected_in_rate = 44100;
        let expected_out_rate = 44000;
        st.set_rate(expected_in_rate, expected_out_rate);
        assert_eq!((expected_in_rate, expected_out_rate), st.get_rate());
    }

    #[test]
    fn test_process_float() {
        let mut st = State::new(1, RATE, RATE, 4).unwrap();
        assert_eq!((1, 1), st.process_float(42, &[8.5], &mut [9.5]).unwrap());
    }

    #[test]
    fn test_skip_zeros() {
        let mut st = State::new(1, RATE, RATE, 4).unwrap();
        let mut avail = INBLOCK as isize;
        let mut rate = 1000;
        let off = 0;
        let fin: Vec<f32> = (0..INBLOCK * 2)
            .map(|i| ((i as f32) / PERIOD * 2.0 * PI).sin() * 0.9)
            .collect();
        let mut fout = vec![0f32; INBLOCK * 4];
        let in_len = avail as usize;
        let out_len = (in_len * rate + RATE - 1) / RATE;
        let f: f32 = 0.0;

        // Don't skip first zero value
        let (_, out_len) = st.process_float(0, &fin[off..off + in_len], &mut fout[..out_len])
            .unwrap();
        assert_eq!(&f, &fout[..out_len as usize][0]);

        // skip first zero value
        st.skip_zeros();
        let (_, out_len) = st.process_float(0, &fin[off..off + in_len], &mut fout[..out_len])
            .unwrap();
        assert_ne!(&f, &fout[..out_len as usize][0]);
    }

    #[test]
    fn test_latency() {
        let mut quality = 2;
        let mut st = State::new(1, RATE, RATE, quality).unwrap();
        assert_eq!((quality * 8) as usize, st.get_input_latency());
        assert_eq!((quality * 8) as usize, st.get_output_latency());

        quality = 4;
        st = State::new(1, RATE, RATE, quality).unwrap();
        assert_eq!((quality * 8) as usize, st.get_input_latency());
        assert_eq!((quality * 8) as usize, st.get_output_latency());
    }

    #[test]
    fn test_quality() {
        let st = State::new(1, RATE, RATE, 4).unwrap();
        st.set_quality(10).unwrap();
        assert_eq!(10, st.get_quality());
    }

}
