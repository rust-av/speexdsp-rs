extern crate speexdsp;

#[cfg(test)]
mod tests {
    use speexdsp::resampler::State;
    const RATE: usize = 48000;

    #[test]
    fn test_quality() {
        let st = State::new(1, RATE, RATE, 4).unwrap();
        st.set_quality(10).unwrap();
        assert_eq!(10, st.get_quality());
    }

    #[test]
    fn test_rate() {
        let mut st = State::new(1, RATE, RATE, 4).unwrap();
        let expected_in_rate = 44100;
        let expected_out_rate = 44000;
        st.set_rate(expected_in_rate, expected_out_rate);
        assert_eq!((expected_in_rate, expected_out_rate), st.get_rate());
    }

}
