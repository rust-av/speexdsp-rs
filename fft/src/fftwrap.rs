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
    os::raw::{c_double, c_float, c_int, c_ulong},
};

use crate::smallft::*;

extern "C" {
    #[no_mangle]
    fn calloc(_: c_ulong, _: c_ulong) -> *mut c_void;
    #[no_mangle]
    fn free(__ptr: *mut c_void);
}
/* Copyright (C) 2007 Jean-Marc Valin

   File: os_support.h
   This is the (tiny) OS abstraction layer. Aside from math.h, this is the
   only place where system headers are allowed.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:

   1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
   IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
   OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
   HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
   STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
/* * Speex wrapper for calloc. To do your own dynamic allocation, all you need to do is replace this function, speex_realloc and speex_free
NOTE: speex_alloc needs to CLEAR THE MEMORY */
#[inline]
unsafe extern "C" fn speex_alloc(mut size: c_int) -> *mut c_void {
    /* WARNING: this is not equivalent to malloc(). If you want to use malloc()
    or your own allocator, YOU NEED TO CLEAR THE MEMORY ALLOCATED. Otherwise
    you will experience strange bugs */
    return calloc(size as c_ulong, 1 as c_int as c_ulong);
}
/* * Speex wrapper for calloc. To do your own dynamic allocation, all you need to do is replace this function, speex_realloc and speex_alloc */
#[inline]
unsafe extern "C" fn speex_free(mut ptr: *mut c_void) {
    free(ptr);
}
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
    let mut table: *mut drft_lookup = 0 as *mut drft_lookup;
    table =
        speex_alloc(::std::mem::size_of::<drft_lookup>() as c_ulong as c_int) as *mut drft_lookup;
    spx_drft_init(table, size);
    return table as *mut c_void;
}
#[no_mangle]
pub unsafe extern "C" fn spx_fft_destroy(mut table: *mut c_void) {
    spx_drft_clear(table as *mut drft_lookup);
    speex_free(table);
}
#[no_mangle]
pub unsafe extern "C" fn spx_fft(
    mut table: *mut c_void,
    mut in_0: *mut c_float,
    mut out: *mut c_float,
) {
    if in_0 == out {
        let mut i: c_int = 0;
        let mut scale: c_float = (1.0f64 / (*(table as *mut drft_lookup)).n as c_double) as c_float;
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
        const OUTPUT: [c_float; 64] = [
            0.1049276,
            0.1307667,
            0.030051949,
            0.3844139,
            -0.2095685,
            0.21490711,
            -0.44597265,
            -0.089277595,
            -0.6891664,
            -0.04035188,
            0.019584754,
            -0.024931509,
            0.010252248,
            -0.017489789,
            0.006553323,
            -0.0132080875,
            0.004674957,
            -0.010483358,
            0.0035711597,
            -0.008629703,
            0.0028556336,
            -0.0073068417,
            0.0023576305,
            -0.0063280445,
            0.0019917116,
            -0.0055828393,
            0.0017111003,
            -0.0050026253,
            0.0014882982,
            -0.0045423717,
            0.0013063122,
            -0.004171742,
            0.0011540875,
            -0.0038695945,
            0.0010239538,
            -0.0036207065,
            0.0009108633,
            -0.0034139454,
            0.0008109063,
            -0.0032412987,
            0.00072141737,
            -0.0030964226,
            0.00064024096,
            -0.0029745854,
            0.000565944,
            -0.0028722635,
            0.0004970236,
            -0.0027864575,
            0.0004323877,
            -0.0027147215,
            0.00037147873,
            -0.0026555816,
            0.00031354558,
            -0.0026075542,
            0.00025797822,
            -0.0025695562,
            0.00020420551,
            -0.0025407895,
            0.00015196204,
            -0.0025205165,
            0.000100664794,
            -0.002508685,
            0.000050345436,
            -0.0025048554,
        ];

        unsafe {
            let table = spx_fft_init(INPUT.len() as i32);
            let mut input = INPUT.clone();
            let mut output = [0.; 64];

            spx_fft(table, input.as_mut_ptr(), output.as_mut_ptr());
            spx_fft_destroy(table);

            assert!(output
                .iter()
                .zip(OUTPUT.iter())
                .all(|(a, b)| (a - b).abs() < EPSILON));
        };
    }

    #[test]
    fn ifft() {
        const OUTPUT: [c_float; 64] = [
            7.4301796,
            19.377022,
            11.15199,
            36.036034,
            8.808279,
            7.028938,
            -22.96261,
            -18.215559,
            -65.566,
            -39.499527,
            0.25965786,
            -15.380432,
            0.38086653,
            -9.388211,
            0.50033426,
            -6.537902,
            0.58552694,
            -4.8504934,
            0.6434908,
            -3.7273202,
            0.6820408,
            -2.9216883,
            0.70659614,
            -2.31283,
            0.7206211,
            -1.8346872,
            0.7262783,
            -1.4483719,
            0.7247982,
            -1.1297483,
            0.7167704,
            -0.86352634,
            0.7022257,
            -0.640151,
            0.6806288,
            -0.4542675,
            0.650887,
            -0.30418348,
            0.61105156,
            -0.19215584,
            0.55802155,
            -0.12564278,
            0.48677158,
            -0.120605946,
            0.38896036,
            -0.20878506,
            0.25000048,
            -0.4552467,
            0.0424819,
            -1.0070467,
            -0.2909441,
            -2.259602,
            -0.89128554,
            -5.6861143,
            -2.2141702,
            -24.455938,
            22.67344,
            26.346508,
            34.141342,
            35.06725,
            35.645927,
            21.961308,
            7.3117857,
            15.947028,
        ];

        unsafe {
            let table = spx_fft_init(INPUT.len() as i32);
            let mut input = INPUT.clone();
            let mut output = [0.; 64];

            spx_ifft(table, input.as_mut_ptr(), output.as_mut_ptr());
            spx_fft_destroy(table);

            assert!(output
                .iter()
                .zip(OUTPUT.iter())
                .all(|(a, b)| (a - b).abs() < EPSILON));
        };
    }
}
