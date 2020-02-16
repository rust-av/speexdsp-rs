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
    os::raw::{c_char, c_double, c_float, c_int, c_long, c_schar, c_ulong, c_ushort},
};

extern "C" {
    pub type _IO_wide_data;
    pub type _IO_codecvt;
    pub type _IO_marker;
    #[no_mangle]
    static mut stderr: *mut FILE;
    #[no_mangle]
    fn fprintf(_: *mut FILE, _: *const c_char, _: ...) -> c_int;
    #[no_mangle]
    fn calloc(_: c_ulong, _: c_ulong) -> *mut c_void;
    #[no_mangle]
    fn free(__ptr: *mut c_void);
    #[no_mangle]
    fn spx_drft_forward(l: *mut drft_lookup, data: *mut c_float);
    #[no_mangle]
    fn spx_drft_backward(l: *mut drft_lookup, data: *mut c_float);
    #[no_mangle]
    fn spx_drft_init(l: *mut drft_lookup, n: c_int);
    #[no_mangle]
    fn spx_drft_clear(l: *mut drft_lookup);
}
pub type __off_t = c_long;
pub type __off64_t = c_long;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _IO_FILE {
    pub _flags: c_int,
    pub _IO_read_ptr: *mut c_char,
    pub _IO_read_end: *mut c_char,
    pub _IO_read_base: *mut c_char,
    pub _IO_write_base: *mut c_char,
    pub _IO_write_ptr: *mut c_char,
    pub _IO_write_end: *mut c_char,
    pub _IO_buf_base: *mut c_char,
    pub _IO_buf_end: *mut c_char,
    pub _IO_save_base: *mut c_char,
    pub _IO_backup_base: *mut c_char,
    pub _IO_save_end: *mut c_char,
    pub _markers: *mut _IO_marker,
    pub _chain: *mut _IO_FILE,
    pub _fileno: c_int,
    pub _flags2: c_int,
    pub _old_offset: __off_t,
    pub _cur_column: c_ushort,
    pub _vtable_offset: c_schar,
    pub _shortbuf: [c_char; 1],
    pub _lock: *mut c_void,
    pub _offset: __off64_t,
    pub _codecvt: *mut _IO_codecvt,
    pub _wide_data: *mut _IO_wide_data,
    pub _freeres_list: *mut _IO_FILE,
    pub _freeres_buf: *mut c_void,
    pub __pad5: c_int,
    pub _mode: c_int,
    pub _unused2: c_char,
}
pub type _IO_lock_t = ();
pub type FILE = _IO_FILE;
/* *******************************************************************
*                                                                  *
* THIS FILE IS PART OF THE OggVorbis SOFTWARE CODEC SOURCE CODE.   *
* USE, DISTRIBUTION AND REPRODUCTION OF THIS LIBRARY SOURCE IS     *
* GOVERNED BY A BSD-STYLE SOURCE LICENSE INCLUDED WITH THIS SOURCE *
* IN 'COPYING'. PLEASE READ THESE TERMS BEFORE DISTRIBUTING.       *
*                                                                  *
* THE OggVorbis SOURCE CODE IS (C) COPYRIGHT 1994-2001             *
* by the XIPHOPHORUS Company http://www.xiph.org/                  *
*                                                                  *
********************************************************************

function: fft transform
last mod: $Id: smallft.h,v 1.3 2003/09/16 18:35:45 jm Exp $

********************************************************************/
/* *
   @file smallft.h
   @brief Discrete Rotational Fourier Transform (DRFT)
*/
/* * Discrete Rotational Fourier Transform lookup */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct drft_lookup {
    pub n: c_int,
    pub trigcache: *mut c_float,
    pub splitcache: *mut c_int,
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
#[inline]
unsafe extern "C" fn speex_warning(mut str: *const c_char) {
    fprintf(
        stderr,
        b"warning: %s\n\x00" as *const u8 as *const c_char,
        str,
    );
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
        speex_warning(b"FFT should not be done in-place\x00" as *const u8 as *const c_char);
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
        speex_warning(b"FFT should not be done in-place\x00" as *const u8 as *const c_char);
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
