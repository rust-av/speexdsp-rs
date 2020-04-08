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
