#![allow(
    dead_code,
    mutable_transmutes,
    non_camel_case_types,
    non_snake_case,
    non_upper_case_globals,
    unused_assignments,
    unused_mut
)]

use std::{
    ffi::c_void,
    os::raw::{c_double, c_float, c_int},
};

use crate::smallft::*;

/* Copyright (C) 2005-2006 Jean-Marc Valin
   File: fftwrap.c

   Wrapper for various FFTs

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

   - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   - Neither the name of the Xiph.org Foundation nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
#[no_mangle]
pub unsafe extern "C" fn spx_fft_init(mut size: c_int) -> *mut c_void {
    let mut table = Box::new(drft_lookup::new(size as usize));

    return Box::into_raw(table) as *mut c_void;
}

#[no_mangle]
pub unsafe extern "C" fn spx_fft_destroy(mut table: *mut c_void) {
    let table = table as *mut drft_lookup;

    Box::from_raw(table);
}

#[no_mangle]
pub unsafe extern "C" fn spx_fft(
    mut table: *mut c_void,
    mut in_0: *mut c_float,
    mut out: *mut c_float,
) {
    if in_0 == out {
        let mut i: c_int = 0;
        let mut scale: c_float =
            (1.0f64 / (*(table as *mut drft_lookup)).n as c_double) as c_float;
        eprintln!("FFT should not be done in-place");
        i = 0 as c_int;
        while i < (*(table as *mut drft_lookup)).n {
            *out.offset(i as isize) = scale * *in_0.offset(i as isize);
            i += 1
        }
    } else {
        let mut i_0: c_int = 0;
        let mut scale_0: c_float =
            (1.0f64 / (*(table as *mut drft_lookup)).n as c_double) as c_float;
        i_0 = 0 as c_int;
        while i_0 < (*(table as *mut drft_lookup)).n {
            *out.offset(i_0 as isize) = scale_0 * *in_0.offset(i_0 as isize);
            i_0 += 1
        }
    }
    spx_drft_forward(table as *mut drft_lookup, out);
}
#[no_mangle]
pub unsafe extern "C" fn spx_ifft(
    mut table: *mut c_void,
    mut in_0: *mut c_float,
    mut out: *mut c_float,
) {
    if in_0 == out {
        eprintln!("FFT should not be done in-place");
    } else {
        let mut i: c_int = 0;
        i = 0 as c_int;
        while i < (*(table as *mut drft_lookup)).n {
            *out.offset(i as isize) = *in_0.offset(i as isize);
            i += 1
        }
    }
    spx_drft_backward(table as *mut drft_lookup, out);
}
#[no_mangle]
pub unsafe extern "C" fn spx_fft_float(
    mut table: *mut c_void,
    mut in_0: *mut c_float,
    mut out: *mut c_float,
) {
    spx_fft(table, in_0, out);
}
#[no_mangle]
pub unsafe extern "C" fn spx_ifft_float(
    mut table: *mut c_void,
    mut in_0: *mut c_float,
    mut out: *mut c_float,
) {
    spx_ifft(table, in_0, out);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::orig;

    const INPUT: [c_float; 64] = [
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

    #[test]
    fn fft() {
        let mut output = [0.; 64];
        unsafe {
            let table = spx_fft_init(INPUT.len() as i32);
            let mut input = INPUT.clone();

            spx_fft(table, input.as_mut_ptr(), output.as_mut_ptr());
            spx_fft_destroy(table);
        };

        let mut expected_output = [0.; 64];
        unsafe {
            let table = orig::spx_fft_init(INPUT.len() as i32);
            let mut input = INPUT.clone();

            orig::spx_fft(
                table,
                input.as_mut_ptr(),
                expected_output.as_mut_ptr(),
            );
            orig::spx_fft_destroy(table);
        };

        assert!(output
            .iter()
            .zip(expected_output.iter())
            .all(|(a, b)| (a - b).abs() < EPSILON));
    }

    #[test]
    fn ifft() {
        let mut output = [0.; 64];
        unsafe {
            let table = spx_fft_init(INPUT.len() as i32);
            let mut input = INPUT.clone();

            spx_ifft(table, input.as_mut_ptr(), output.as_mut_ptr());
            spx_fft_destroy(table);
        };

        let mut expected_output = [0.; 64];
        unsafe {
            let table = orig::spx_fft_init(INPUT.len() as i32);
            let mut input = INPUT.clone();

            orig::spx_ifft(
                table,
                input.as_mut_ptr(),
                expected_output.as_mut_ptr(),
            );
            orig::spx_fft_destroy(table);
        };

        assert!(output
            .iter()
            .zip(expected_output.iter())
            .all(|(a, b)| (a - b).abs() < EPSILON));
    }
}
