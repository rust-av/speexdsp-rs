extern crate speexdsp_fft;

mod orig;

use std::ptr::null_mut;

use speexdsp_fft::*;

use crate::orig::fftwrap as original_fftwrap;
use crate::orig::smallft as original_smallft;

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
const N: usize = 2;
const SPLITCACHE_SIZE: usize = 32;

macro_rules! test_fftwrap {
    ($func: ident) => {
        let mut output = [0.; 64];
        let mut drft_lookup = DrftLookup::new(INPUT.len());
        let mut input = INPUT.clone();

        drft_lookup.$func(&mut input, &mut output);

        let mut expected_output = [0.; 64];
        unsafe {
            let drft_lookup =
                original_fftwrap::spx_fft_init(INPUT.len() as i32);
            let mut input = INPUT.clone();

            original_fftwrap::$func(
                drft_lookup,
                input.as_mut_ptr(),
                expected_output.as_mut_ptr(),
            );
            original_fftwrap::spx_fft_destroy(drft_lookup);
        };

        assert!(output
            .iter()
            .zip(expected_output.iter())
            .all(|(a, b)| (a - b).abs() < EPSILON));
    };
}

macro_rules! drft {
    ($func: ident) => {
        let mut drft_lookup = DrftLookup::new(N);
        let mut data = vec![0.; 32];

        drft_lookup.$func(&mut data);

        let mut original_drft_lookup = original_smallft::drft_lookup {
            n: 0,
            trigcache: null_mut(),
            splitcache: null_mut(),
        };
        let mut data_original = vec![0.; 32];
        unsafe {
            original_smallft::spx_drft_init(
                &mut original_drft_lookup,
                N as i32,
            );
            original_smallft::$func(
                &mut original_drft_lookup,
                data_original.as_mut_ptr(),
            );
        }

        let expected_trigcache = unsafe {
            Vec::from_raw_parts(
                original_drft_lookup.trigcache as *mut f32,
                3 * N,
                3 * N,
            )
        };

        let expected_splitcache = unsafe {
            Vec::from_raw_parts(
                original_drft_lookup.splitcache as *mut i32,
                SPLITCACHE_SIZE,
                SPLITCACHE_SIZE,
            )
        };

        assert!(drft_lookup
            .trigcache
            .iter()
            .zip(expected_trigcache.iter())
            .all(|(&a, &b)| (a - b).abs() < EPSILON));

        assert_eq!(&drft_lookup.splitcache, &expected_splitcache);
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

#[test]
fn fdrffti() {
    let drft_lookup = DrftLookup::new(N);

    let mut original_drft_lookup = original_smallft::drft_lookup {
        n: 0,
        trigcache: null_mut(),
        splitcache: null_mut(),
    };
    unsafe {
        original_smallft::spx_drft_init(&mut original_drft_lookup, N as i32);
    }

    let expected_trigcache = unsafe {
        Vec::from_raw_parts(
            original_drft_lookup.trigcache as *mut f32,
            3 * N,
            3 * N,
        )
    };

    let expected_splitcache = unsafe {
        Vec::from_raw_parts(
            original_drft_lookup.splitcache as *mut i32,
            SPLITCACHE_SIZE,
            SPLITCACHE_SIZE,
        )
    };

    assert!(drft_lookup
        .trigcache
        .iter()
        .zip(expected_trigcache.iter())
        .all(|(&a, &b)| (a - b).abs() < EPSILON));

    assert_eq!(&drft_lookup.splitcache, &expected_splitcache);
}

#[test]
fn drftf1() {
    drft!(spx_drft_forward);
}

#[test]
fn drftb1() {
    drft!(spx_drft_backward);
}
