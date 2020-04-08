extern crate speexdsp_fft;

mod orig;

use crate::orig::smallft as original;
use speexdsp_fft::*;

use std::os::raw::{c_float, c_int};
use std::ptr::null_mut;

const EPSILON: c_float = 1e-6;
const N: usize = 2;
const SPLITCACHE_SIZE: usize = 32;

macro_rules! drft {
    ($func: ident) => {
        let mut drft_lookup = drft_lookup::new(N);
        let mut data = vec![0.; 32];
        unsafe {
            $func(&mut drft_lookup, data.as_mut_ptr());
        }

        let mut original_drft_lookup = original::drft_lookup {
            n: 0,
            trigcache: null_mut(),
            splitcache: null_mut(),
        };
        let mut data_original = vec![0.; 32];
        unsafe {
            original::spx_drft_init(&mut original_drft_lookup, N as c_int);
            original::$func(
                &mut original_drft_lookup,
                data_original.as_mut_ptr(),
            );
        }

        let expected_trigcache = unsafe {
            Vec::from_raw_parts(
                original_drft_lookup.trigcache as *mut c_float,
                3 * N,
                3 * N,
            )
        };

        let expected_splitcache = unsafe {
            Vec::from_raw_parts(
                original_drft_lookup.splitcache as *mut c_int,
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
fn fdrffti() {
    let drft_lookup = drft_lookup::new(N);

    let mut original_drft_lookup = original::drft_lookup {
        n: 0,
        trigcache: null_mut(),
        splitcache: null_mut(),
    };
    unsafe {
        original::spx_drft_init(&mut original_drft_lookup, N as c_int);
    }

    let expected_trigcache = unsafe {
        Vec::from_raw_parts(
            original_drft_lookup.trigcache as *mut c_float,
            3 * N,
            3 * N,
        )
    };

    let expected_splitcache = unsafe {
        Vec::from_raw_parts(
            original_drft_lookup.splitcache as *mut c_int,
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
