extern crate speexdsp_fft;

mod orig;

use crate::orig::fftwrap as original;
use speexdsp_fft::*;

const INPUT: [f32; 64] = [
    1.0,
    1.999002,
    2.832922,
    3.3986773,
    3.6254587,
    3.484332,
    2.992106,
    2.2089396,
    1.230057,
    0.17277533,
    -0.8392913,
    -1.692595,
    -2.2982898,
    -2.6031098,
    -2.5950398,
    -2.303094,
    -1.7913916,
    -1.1485368,
    -0.47395423,
    0.13674068,
    0.60517603,
    0.8805469,
    0.94540364,
    0.8161983,
    0.53870463,
    0.17917383,
    -0.18727368,
    -0.4891553,
    -0.6700415,
    -0.69716847,
    -0.5659317,
    -0.29971862,
    0.05482897,
    0.4362133,
    0.77873474,
    1.0236322,
    1.1288913,
    1.0761548,
    0.8736177,
    0.55450416,
    0.1713717,
    -0.21276897,
    -0.53509253,
    -0.74354666,
    -0.8057217,
    -0.7144466,
    -0.48923856,
    -0.17328042,
    0.17359133,
    0.48484406,
    0.6984716,
    0.7677859,
    0.6698748,
    0.41039884,
    0.023795009,
    -0.43124664,
    -0.8804407,
    -1.2453988,
    -1.4562676,
    -1.4634898,
    -1.2468007,
    -0.8200231,
    -0.23070133,
    0.4454986,
];

const EPSILON: f32 = 1e-8;

macro_rules! test_fftwrap {
    ($func: ident) => {
        let mut output = [0.; 64];
        let mut table = spx_fft_init(INPUT.len());
        let mut input = INPUT.clone();

        $func(&mut table, &mut input, &mut output);

        let mut expected_output = [0.; 64];
        unsafe {
            let table = original::spx_fft_init(INPUT.len() as i32);
            let mut input = INPUT.clone();

            original::$func(
                table,
                input.as_mut_ptr(),
                expected_output.as_mut_ptr(),
            );
            original::spx_fft_destroy(table);
        };

        assert!(output
            .iter()
            .zip(expected_output.iter())
            .all(|(a, b)| (a - b).abs() < EPSILON));
    };
}

#[test]
fn fft() {
    test_fftwrap!(spx_fft);
}

#[test]
fn ifft() {
    test_fftwrap!(spx_ifft);
}
