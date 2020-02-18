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

extern "C" {
    #[no_mangle]
    fn cos(_: c_double) -> c_double;
    #[no_mangle]
    fn sin(_: c_double) -> c_double;
    #[no_mangle]
    fn calloc(_: c_ulong, _: c_ulong) -> *mut c_void;
    #[no_mangle]
    fn free(__ptr: *mut c_void);
}
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
#[inline]
unsafe extern "C" fn speex_alloc(mut size: c_int) -> *mut c_void {
    return calloc(size as c_ulong, 1 as c_int as c_ulong);
}
#[inline]
unsafe extern "C" fn speex_free(mut ptr: *mut c_void) {
    free(ptr);
}
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

function: *unnormalized* fft transform
last mod: $Id: smallft.c,v 1.19 2003/10/08 05:12:37 jm Exp $

********************************************************************/
/* FFT implementation from OggSquish, minus cosine transforms,
 * minus all but radix 2/4 case.  In Vorbis we only need this
 * cut-down version.
 *
 * To do more than just power-of-two sized vectors, see the full
 * version I wrote for NetLib.
 *
 * Note that the packing is a little strange; rather than the FFT r/i
 * packing following R_0, I_n, R_1, I_1, R_2, I_2 ... R_n-1, I_n-1,
 * it follows R_0, R_1, I_1, R_2, I_2 ... R_n-1, I_n-1, I_n like the
 * FORTRAN version
 */
unsafe extern "C" fn drfti1(mut n: c_int, mut wa: *mut c_float, mut ifac: *mut c_int) {
    static mut ntryh: [c_int; 4] = [4 as c_int, 2 as c_int, 3 as c_int, 5 as c_int];
    static mut tpi: c_float = 6.28318530717958648f32;
    let mut arg: c_float = 0.;
    let mut argh: c_float = 0.;
    let mut argld: c_float = 0.;
    let mut fi: c_float = 0.;
    let mut ntry: c_int = 0 as c_int;
    let mut i: c_int = 0;
    let mut j: c_int = -(1 as c_int);
    let mut k1: c_int = 0;
    let mut l1: c_int = 0;
    let mut l2: c_int = 0;
    let mut ib: c_int = 0;
    let mut ld: c_int = 0;
    let mut ii: c_int = 0;
    let mut ip: c_int = 0;
    let mut is: c_int = 0;
    let mut nq: c_int = 0;
    let mut nr: c_int = 0;
    let mut ido: c_int = 0;
    let mut ipm: c_int = 0;
    let mut nfm1: c_int = 0;
    let mut nl: c_int = n;
    let mut nf: c_int = 0 as c_int;
    'c_10244: loop {
        j += 1;
        if j < 4 as c_int {
            ntry = ntryh[j as usize]
        } else {
            ntry += 2 as c_int
        }
        loop {
            nq = nl / ntry;
            nr = nl - ntry * nq;
            if nr != 0 as c_int {
                break;
            }
            nf += 1;
            *ifac.offset((nf + 1 as c_int) as isize) = ntry;
            nl = nq;
            if !(ntry != 2 as c_int) {
                if !(nf == 1 as c_int) {
                    i = 1 as c_int;
                    while i < nf {
                        ib = nf - i + 1 as c_int;
                        *ifac.offset((ib + 1 as c_int) as isize) = *ifac.offset(ib as isize);
                        i += 1
                    }
                    *ifac.offset(2 as c_int as isize) = 2 as c_int
                }
            }
            if !(nl != 1 as c_int) {
                break 'c_10244;
            }
        }
    }
    *ifac.offset(0 as c_int as isize) = n;
    *ifac.offset(1 as c_int as isize) = nf;
    argh = tpi / n as c_float;
    is = 0 as c_int;
    nfm1 = nf - 1 as c_int;
    l1 = 1 as c_int;
    if nfm1 == 0 as c_int {
        return;
    }
    k1 = 0 as c_int;
    while k1 < nfm1 {
        ip = *ifac.offset((k1 + 2 as c_int) as isize);
        ld = 0 as c_int;
        l2 = l1 * ip;
        ido = n / l2;
        ipm = ip - 1 as c_int;
        j = 0 as c_int;
        while j < ipm {
            ld += l1;
            i = is;
            argld = ld as c_float * argh;
            fi = 0.0f32;
            ii = 2 as c_int;
            while ii < ido {
                fi += 1.0f32;
                arg = fi * argld;
                let fresh0 = i;
                i = i + 1;
                *wa.offset(fresh0 as isize) = cos(arg as c_double) as c_float;
                let fresh1 = i;
                i = i + 1;
                *wa.offset(fresh1 as isize) = sin(arg as c_double) as c_float;
                ii += 2 as c_int
            }
            is += ido;
            j += 1
        }
        l1 = l2;
        k1 += 1
    }
}
unsafe extern "C" fn fdrffti(mut n: c_int, mut wsave: *mut c_float, mut ifac: *mut c_int) {
    if n == 1 as c_int {
        return;
    }
    drfti1(n, wsave.offset(n as isize), ifac);
}
unsafe extern "C" fn dradf2(
    mut ido: c_int,
    mut l1: c_int,
    mut cc: *mut c_float,
    mut ch: *mut c_float,
    mut wa1: *mut c_float,
) {
    let mut i: c_int = 0;
    let mut k: c_int = 0;
    let mut ti2: c_float = 0.;
    let mut tr2: c_float = 0.;
    let mut t0: c_int = 0;
    let mut t1: c_int = 0;
    let mut t2: c_int = 0;
    let mut t3: c_int = 0;
    let mut t4: c_int = 0;
    let mut t5: c_int = 0;
    let mut t6: c_int = 0;
    t1 = 0 as c_int;
    t2 = l1 * ido;
    t0 = t2;
    t3 = ido << 1 as c_int;
    k = 0 as c_int;
    while k < l1 {
        *ch.offset((t1 << 1 as c_int) as isize) = *cc.offset(t1 as isize) + *cc.offset(t2 as isize);
        *ch.offset(((t1 << 1 as c_int) + t3 - 1 as c_int) as isize) =
            *cc.offset(t1 as isize) - *cc.offset(t2 as isize);
        t1 += ido;
        t2 += ido;
        k += 1
    }
    if ido < 2 as c_int {
        return;
    }
    if !(ido == 2 as c_int) {
        t1 = 0 as c_int;
        t2 = t0;
        k = 0 as c_int;
        while k < l1 {
            t3 = t2;
            t4 = (t1 << 1 as c_int) + (ido << 1 as c_int);
            t5 = t1;
            t6 = t1 + t1;
            i = 2 as c_int;
            while i < ido {
                t3 += 2 as c_int;
                t4 -= 2 as c_int;
                t5 += 2 as c_int;
                t6 += 2 as c_int;
                tr2 = *wa1.offset((i - 2 as c_int) as isize)
                    * *cc.offset((t3 - 1 as c_int) as isize)
                    + *wa1.offset((i - 1 as c_int) as isize) * *cc.offset(t3 as isize);
                ti2 = *wa1.offset((i - 2 as c_int) as isize) * *cc.offset(t3 as isize)
                    - *wa1.offset((i - 1 as c_int) as isize)
                        * *cc.offset((t3 - 1 as c_int) as isize);
                *ch.offset(t6 as isize) = *cc.offset(t5 as isize) + ti2;
                *ch.offset(t4 as isize) = ti2 - *cc.offset(t5 as isize);
                *ch.offset((t6 - 1 as c_int) as isize) =
                    *cc.offset((t5 - 1 as c_int) as isize) + tr2;
                *ch.offset((t4 - 1 as c_int) as isize) =
                    *cc.offset((t5 - 1 as c_int) as isize) - tr2;
                i += 2 as c_int
            }
            t1 += ido;
            t2 += ido;
            k += 1
        }
        if ido % 2 as c_int == 1 as c_int {
            return;
        }
    }
    t1 = ido;
    t2 = t1 - 1 as c_int;
    t3 = t2;
    t2 += t0;
    k = 0 as c_int;
    while k < l1 {
        *ch.offset(t1 as isize) = -*cc.offset(t2 as isize);
        *ch.offset((t1 - 1 as c_int) as isize) = *cc.offset(t3 as isize);
        t1 += ido << 1 as c_int;
        t2 += ido;
        t3 += ido;
        k += 1
    }
}
unsafe extern "C" fn dradf4(
    mut ido: c_int,
    mut l1: c_int,
    mut cc: *mut c_float,
    mut ch: *mut c_float,
    mut wa1: *mut c_float,
    mut wa2: *mut c_float,
    mut wa3: *mut c_float,
) {
    static mut hsqt2: c_float = 0.70710678118654752f32;
    let mut i: c_int = 0;
    let mut k: c_int = 0;
    let mut t0: c_int = 0;
    let mut t1: c_int = 0;
    let mut t2: c_int = 0;
    let mut t3: c_int = 0;
    let mut t4: c_int = 0;
    let mut t5: c_int = 0;
    let mut t6: c_int = 0;
    let mut ci2: c_float = 0.;
    let mut ci3: c_float = 0.;
    let mut ci4: c_float = 0.;
    let mut cr2: c_float = 0.;
    let mut cr3: c_float = 0.;
    let mut cr4: c_float = 0.;
    let mut ti1: c_float = 0.;
    let mut ti2: c_float = 0.;
    let mut ti3: c_float = 0.;
    let mut ti4: c_float = 0.;
    let mut tr1: c_float = 0.;
    let mut tr2: c_float = 0.;
    let mut tr3: c_float = 0.;
    let mut tr4: c_float = 0.;
    t0 = l1 * ido;
    t1 = t0;
    t4 = t1 << 1 as c_int;
    t2 = t1 + (t1 << 1 as c_int);
    t3 = 0 as c_int;
    k = 0 as c_int;
    while k < l1 {
        tr1 = *cc.offset(t1 as isize) + *cc.offset(t2 as isize);
        tr2 = *cc.offset(t3 as isize) + *cc.offset(t4 as isize);
        t5 = t3 << 2 as c_int;
        *ch.offset(t5 as isize) = tr1 + tr2;
        *ch.offset(((ido << 2 as c_int) + t5 - 1 as c_int) as isize) = tr2 - tr1;
        t5 += ido << 1 as c_int;
        *ch.offset((t5 - 1 as c_int) as isize) = *cc.offset(t3 as isize) - *cc.offset(t4 as isize);
        *ch.offset(t5 as isize) = *cc.offset(t2 as isize) - *cc.offset(t1 as isize);
        t1 += ido;
        t2 += ido;
        t3 += ido;
        t4 += ido;
        k += 1
    }
    if ido < 2 as c_int {
        return;
    }
    if !(ido == 2 as c_int) {
        t1 = 0 as c_int;
        k = 0 as c_int;
        while k < l1 {
            t2 = t1;
            t4 = t1 << 2 as c_int;
            t6 = ido << 1 as c_int;
            t5 = t6 + t4;
            i = 2 as c_int;
            while i < ido {
                t2 += 2 as c_int;
                t3 = t2;
                t4 += 2 as c_int;
                t5 -= 2 as c_int;
                t3 += t0;
                cr2 = *wa1.offset((i - 2 as c_int) as isize)
                    * *cc.offset((t3 - 1 as c_int) as isize)
                    + *wa1.offset((i - 1 as c_int) as isize) * *cc.offset(t3 as isize);
                ci2 = *wa1.offset((i - 2 as c_int) as isize) * *cc.offset(t3 as isize)
                    - *wa1.offset((i - 1 as c_int) as isize)
                        * *cc.offset((t3 - 1 as c_int) as isize);
                t3 += t0;
                cr3 = *wa2.offset((i - 2 as c_int) as isize)
                    * *cc.offset((t3 - 1 as c_int) as isize)
                    + *wa2.offset((i - 1 as c_int) as isize) * *cc.offset(t3 as isize);
                ci3 = *wa2.offset((i - 2 as c_int) as isize) * *cc.offset(t3 as isize)
                    - *wa2.offset((i - 1 as c_int) as isize)
                        * *cc.offset((t3 - 1 as c_int) as isize);
                t3 += t0;
                cr4 = *wa3.offset((i - 2 as c_int) as isize)
                    * *cc.offset((t3 - 1 as c_int) as isize)
                    + *wa3.offset((i - 1 as c_int) as isize) * *cc.offset(t3 as isize);
                ci4 = *wa3.offset((i - 2 as c_int) as isize) * *cc.offset(t3 as isize)
                    - *wa3.offset((i - 1 as c_int) as isize)
                        * *cc.offset((t3 - 1 as c_int) as isize);
                tr1 = cr2 + cr4;
                tr4 = cr4 - cr2;
                ti1 = ci2 + ci4;
                ti4 = ci2 - ci4;
                ti2 = *cc.offset(t2 as isize) + ci3;
                ti3 = *cc.offset(t2 as isize) - ci3;
                tr2 = *cc.offset((t2 - 1 as c_int) as isize) + cr3;
                tr3 = *cc.offset((t2 - 1 as c_int) as isize) - cr3;
                *ch.offset((t4 - 1 as c_int) as isize) = tr1 + tr2;
                *ch.offset(t4 as isize) = ti1 + ti2;
                *ch.offset((t5 - 1 as c_int) as isize) = tr3 - ti4;
                *ch.offset(t5 as isize) = tr4 - ti3;
                *ch.offset((t4 + t6 - 1 as c_int) as isize) = ti4 + tr3;
                *ch.offset((t4 + t6) as isize) = tr4 + ti3;
                *ch.offset((t5 + t6 - 1 as c_int) as isize) = tr2 - tr1;
                *ch.offset((t5 + t6) as isize) = ti1 - ti2;
                i += 2 as c_int
            }
            t1 += ido;
            k += 1
        }
        if ido & 1 as c_int != 0 {
            return;
        }
    }
    t1 = t0 + ido - 1 as c_int;
    t2 = t1 + (t0 << 1 as c_int);
    t3 = ido << 2 as c_int;
    t4 = ido;
    t5 = ido << 1 as c_int;
    t6 = ido;
    k = 0 as c_int;
    while k < l1 {
        ti1 = -hsqt2 * (*cc.offset(t1 as isize) + *cc.offset(t2 as isize));
        tr1 = hsqt2 * (*cc.offset(t1 as isize) - *cc.offset(t2 as isize));
        *ch.offset((t4 - 1 as c_int) as isize) = tr1 + *cc.offset((t6 - 1 as c_int) as isize);
        *ch.offset((t4 + t5 - 1 as c_int) as isize) = *cc.offset((t6 - 1 as c_int) as isize) - tr1;
        *ch.offset(t4 as isize) = ti1 - *cc.offset((t1 + t0) as isize);
        *ch.offset((t4 + t5) as isize) = ti1 + *cc.offset((t1 + t0) as isize);
        t1 += ido;
        t2 += ido;
        t4 += t3;
        t6 += ido;
        k += 1
    }
}
unsafe extern "C" fn dradfg(
    mut ido: c_int,
    mut ip: c_int,
    mut l1: c_int,
    mut idl1: c_int,
    mut cc: *mut c_float,
    mut c1: *mut c_float,
    mut c2: *mut c_float,
    mut ch: *mut c_float,
    mut ch2: *mut c_float,
    mut wa: *mut c_float,
) {
    static mut tpi: c_float = 6.283185307179586f32;
    let mut idij: c_int = 0;
    let mut ipph: c_int = 0;
    let mut i: c_int = 0;
    let mut j: c_int = 0;
    let mut k: c_int = 0;
    let mut l: c_int = 0;
    let mut ic: c_int = 0;
    let mut ik: c_int = 0;
    let mut is: c_int = 0;
    let mut t0: c_int = 0;
    let mut t1: c_int = 0;
    let mut t2: c_int = 0;
    let mut t3: c_int = 0;
    let mut t4: c_int = 0;
    let mut t5: c_int = 0;
    let mut t6: c_int = 0;
    let mut t7: c_int = 0;
    let mut t8: c_int = 0;
    let mut t9: c_int = 0;
    let mut t10: c_int = 0;
    let mut dc2: c_float = 0.;
    let mut ai1: c_float = 0.;
    let mut ai2: c_float = 0.;
    let mut ar1: c_float = 0.;
    let mut ar2: c_float = 0.;
    let mut ds2: c_float = 0.;
    let mut nbd: c_int = 0;
    let mut dcp: c_float = 0.;
    let mut arg: c_float = 0.;
    let mut dsp: c_float = 0.;
    let mut ar1h: c_float = 0.;
    let mut ar2h: c_float = 0.;
    let mut idp2: c_int = 0;
    let mut ipp2: c_int = 0;
    arg = tpi / ip as c_float;
    dcp = cos(arg as c_double) as c_float;
    dsp = sin(arg as c_double) as c_float;
    ipph = ip + 1 as c_int >> 1 as c_int;
    ipp2 = ip;
    idp2 = ido;
    nbd = ido - 1 as c_int >> 1 as c_int;
    t0 = l1 * ido;
    t10 = ip * ido;
    if !(ido == 1 as c_int) {
        ik = 0 as c_int;
        while ik < idl1 {
            *ch2.offset(ik as isize) = *c2.offset(ik as isize);
            ik += 1
        }
        t1 = 0 as c_int;
        j = 1 as c_int;
        while j < ip {
            t1 += t0;
            t2 = t1;
            k = 0 as c_int;
            while k < l1 {
                *ch.offset(t2 as isize) = *c1.offset(t2 as isize);
                t2 += ido;
                k += 1
            }
            j += 1
        }
        is = -ido;
        t1 = 0 as c_int;
        if nbd > l1 {
            j = 1 as c_int;
            while j < ip {
                t1 += t0;
                is += ido;
                t2 = -ido + t1;
                k = 0 as c_int;
                while k < l1 {
                    idij = is - 1 as c_int;
                    t2 += ido;
                    t3 = t2;
                    i = 2 as c_int;
                    while i < ido {
                        idij += 2 as c_int;
                        t3 += 2 as c_int;
                        *ch.offset((t3 - 1 as c_int) as isize) = *wa
                            .offset((idij - 1 as c_int) as isize)
                            * *c1.offset((t3 - 1 as c_int) as isize)
                            + *wa.offset(idij as isize) * *c1.offset(t3 as isize);
                        *ch.offset(t3 as isize) = *wa.offset((idij - 1 as c_int) as isize)
                            * *c1.offset(t3 as isize)
                            - *wa.offset(idij as isize) * *c1.offset((t3 - 1 as c_int) as isize);
                        i += 2 as c_int
                    }
                    k += 1
                }
                j += 1
            }
        } else {
            j = 1 as c_int;
            while j < ip {
                is += ido;
                idij = is - 1 as c_int;
                t1 += t0;
                t2 = t1;
                i = 2 as c_int;
                while i < ido {
                    idij += 2 as c_int;
                    t2 += 2 as c_int;
                    t3 = t2;
                    k = 0 as c_int;
                    while k < l1 {
                        *ch.offset((t3 - 1 as c_int) as isize) = *wa
                            .offset((idij - 1 as c_int) as isize)
                            * *c1.offset((t3 - 1 as c_int) as isize)
                            + *wa.offset(idij as isize) * *c1.offset(t3 as isize);
                        *ch.offset(t3 as isize) = *wa.offset((idij - 1 as c_int) as isize)
                            * *c1.offset(t3 as isize)
                            - *wa.offset(idij as isize) * *c1.offset((t3 - 1 as c_int) as isize);
                        t3 += ido;
                        k += 1
                    }
                    i += 2 as c_int
                }
                j += 1
            }
        }
        t1 = 0 as c_int;
        t2 = ipp2 * t0;
        if nbd < l1 {
            j = 1 as c_int;
            while j < ipph {
                t1 += t0;
                t2 -= t0;
                t3 = t1;
                t4 = t2;
                i = 2 as c_int;
                while i < ido {
                    t3 += 2 as c_int;
                    t4 += 2 as c_int;
                    t5 = t3 - ido;
                    t6 = t4 - ido;
                    k = 0 as c_int;
                    while k < l1 {
                        t5 += ido;
                        t6 += ido;
                        *c1.offset((t5 - 1 as c_int) as isize) = *ch
                            .offset((t5 - 1 as c_int) as isize)
                            + *ch.offset((t6 - 1 as c_int) as isize);
                        *c1.offset((t6 - 1 as c_int) as isize) =
                            *ch.offset(t5 as isize) - *ch.offset(t6 as isize);
                        *c1.offset(t5 as isize) = *ch.offset(t5 as isize) + *ch.offset(t6 as isize);
                        *c1.offset(t6 as isize) = *ch.offset((t6 - 1 as c_int) as isize)
                            - *ch.offset((t5 - 1 as c_int) as isize);
                        k += 1
                    }
                    i += 2 as c_int
                }
                j += 1
            }
        } else {
            j = 1 as c_int;
            while j < ipph {
                t1 += t0;
                t2 -= t0;
                t3 = t1;
                t4 = t2;
                k = 0 as c_int;
                while k < l1 {
                    t5 = t3;
                    t6 = t4;
                    i = 2 as c_int;
                    while i < ido {
                        t5 += 2 as c_int;
                        t6 += 2 as c_int;
                        *c1.offset((t5 - 1 as c_int) as isize) = *ch
                            .offset((t5 - 1 as c_int) as isize)
                            + *ch.offset((t6 - 1 as c_int) as isize);
                        *c1.offset((t6 - 1 as c_int) as isize) =
                            *ch.offset(t5 as isize) - *ch.offset(t6 as isize);
                        *c1.offset(t5 as isize) = *ch.offset(t5 as isize) + *ch.offset(t6 as isize);
                        *c1.offset(t6 as isize) = *ch.offset((t6 - 1 as c_int) as isize)
                            - *ch.offset((t5 - 1 as c_int) as isize);
                        i += 2 as c_int
                    }
                    t3 += ido;
                    t4 += ido;
                    k += 1
                }
                j += 1
            }
        }
    }
    ik = 0 as c_int;
    while ik < idl1 {
        *c2.offset(ik as isize) = *ch2.offset(ik as isize);
        ik += 1
    }
    t1 = 0 as c_int;
    t2 = ipp2 * idl1;
    j = 1 as c_int;
    while j < ipph {
        t1 += t0;
        t2 -= t0;
        t3 = t1 - ido;
        t4 = t2 - ido;
        k = 0 as c_int;
        while k < l1 {
            t3 += ido;
            t4 += ido;
            *c1.offset(t3 as isize) = *ch.offset(t3 as isize) + *ch.offset(t4 as isize);
            *c1.offset(t4 as isize) = *ch.offset(t4 as isize) - *ch.offset(t3 as isize);
            k += 1
        }
        j += 1
    }
    ar1 = 1.0f32;
    ai1 = 0.0f32;
    t1 = 0 as c_int;
    t2 = ipp2 * idl1;
    t3 = (ip - 1 as c_int) * idl1;
    l = 1 as c_int;
    while l < ipph {
        t1 += idl1;
        t2 -= idl1;
        ar1h = dcp * ar1 - dsp * ai1;
        ai1 = dcp * ai1 + dsp * ar1;
        ar1 = ar1h;
        t4 = t1;
        t5 = t2;
        t6 = t3;
        t7 = idl1;
        ik = 0 as c_int;
        while ik < idl1 {
            let fresh2 = t7;
            t7 = t7 + 1;
            let fresh3 = t4;
            t4 = t4 + 1;
            *ch2.offset(fresh3 as isize) =
                *c2.offset(ik as isize) + ar1 * *c2.offset(fresh2 as isize);
            let fresh4 = t6;
            t6 = t6 + 1;
            let fresh5 = t5;
            t5 = t5 + 1;
            *ch2.offset(fresh5 as isize) = ai1 * *c2.offset(fresh4 as isize);
            ik += 1
        }
        dc2 = ar1;
        ds2 = ai1;
        ar2 = ar1;
        ai2 = ai1;
        t4 = idl1;
        t5 = (ipp2 - 1 as c_int) * idl1;
        j = 2 as c_int;
        while j < ipph {
            t4 += idl1;
            t5 -= idl1;
            ar2h = dc2 * ar2 - ds2 * ai2;
            ai2 = dc2 * ai2 + ds2 * ar2;
            ar2 = ar2h;
            t6 = t1;
            t7 = t2;
            t8 = t4;
            t9 = t5;
            ik = 0 as c_int;
            while ik < idl1 {
                let fresh6 = t8;
                t8 = t8 + 1;
                let fresh7 = t6;
                t6 = t6 + 1;
                *ch2.offset(fresh7 as isize) += ar2 * *c2.offset(fresh6 as isize);
                let fresh8 = t9;
                t9 = t9 + 1;
                let fresh9 = t7;
                t7 = t7 + 1;
                *ch2.offset(fresh9 as isize) += ai2 * *c2.offset(fresh8 as isize);
                ik += 1
            }
            j += 1
        }
        l += 1
    }
    t1 = 0 as c_int;
    j = 1 as c_int;
    while j < ipph {
        t1 += idl1;
        t2 = t1;
        ik = 0 as c_int;
        while ik < idl1 {
            let fresh10 = t2;
            t2 = t2 + 1;
            *ch2.offset(ik as isize) += *c2.offset(fresh10 as isize);
            ik += 1
        }
        j += 1
    }
    if ido < l1 {
        i = 0 as c_int;
        while i < ido {
            t1 = i;
            t2 = i;
            k = 0 as c_int;
            while k < l1 {
                *cc.offset(t2 as isize) = *ch.offset(t1 as isize);
                t1 += ido;
                t2 += t10;
                k += 1
            }
            i += 1
        }
    } else {
        t1 = 0 as c_int;
        t2 = 0 as c_int;
        k = 0 as c_int;
        while k < l1 {
            t3 = t1;
            t4 = t2;
            i = 0 as c_int;
            while i < ido {
                let fresh11 = t3;
                t3 = t3 + 1;
                let fresh12 = t4;
                t4 = t4 + 1;
                *cc.offset(fresh12 as isize) = *ch.offset(fresh11 as isize);
                i += 1
            }
            t1 += ido;
            t2 += t10;
            k += 1
        }
    }
    t1 = 0 as c_int;
    t2 = ido << 1 as c_int;
    t3 = 0 as c_int;
    t4 = ipp2 * t0;
    j = 1 as c_int;
    while j < ipph {
        t1 += t2;
        t3 += t0;
        t4 -= t0;
        t5 = t1;
        t6 = t3;
        t7 = t4;
        k = 0 as c_int;
        while k < l1 {
            *cc.offset((t5 - 1 as c_int) as isize) = *ch.offset(t6 as isize);
            *cc.offset(t5 as isize) = *ch.offset(t7 as isize);
            t5 += t10;
            t6 += ido;
            t7 += ido;
            k += 1
        }
        j += 1
    }
    if ido == 1 as c_int {
        return;
    }
    if nbd < l1 {
        t1 = -ido;
        t3 = 0 as c_int;
        t4 = 0 as c_int;
        t5 = ipp2 * t0;
        j = 1 as c_int;
        while j < ipph {
            t1 += t2;
            t3 += t2;
            t4 += t0;
            t5 -= t0;
            i = 2 as c_int;
            while i < ido {
                t6 = idp2 + t1 - i;
                t7 = i + t3;
                t8 = i + t4;
                t9 = i + t5;
                k = 0 as c_int;
                while k < l1 {
                    *cc.offset((t7 - 1 as c_int) as isize) = *ch.offset((t8 - 1 as c_int) as isize)
                        + *ch.offset((t9 - 1 as c_int) as isize);
                    *cc.offset((t6 - 1 as c_int) as isize) = *ch.offset((t8 - 1 as c_int) as isize)
                        - *ch.offset((t9 - 1 as c_int) as isize);
                    *cc.offset(t7 as isize) = *ch.offset(t8 as isize) + *ch.offset(t9 as isize);
                    *cc.offset(t6 as isize) = *ch.offset(t9 as isize) - *ch.offset(t8 as isize);
                    t6 += t10;
                    t7 += t10;
                    t8 += ido;
                    t9 += ido;
                    k += 1
                }
                i += 2 as c_int
            }
            j += 1
        }
        return;
    } else {
        t1 = -ido;
        t3 = 0 as c_int;
        t4 = 0 as c_int;
        t5 = ipp2 * t0;
        j = 1 as c_int;
        while j < ipph {
            t1 += t2;
            t3 += t2;
            t4 += t0;
            t5 -= t0;
            t6 = t1;
            t7 = t3;
            t8 = t4;
            t9 = t5;
            k = 0 as c_int;
            while k < l1 {
                i = 2 as c_int;
                while i < ido {
                    ic = idp2 - i;
                    *cc.offset((i + t7 - 1 as c_int) as isize) = *ch
                        .offset((i + t8 - 1 as c_int) as isize)
                        + *ch.offset((i + t9 - 1 as c_int) as isize);
                    *cc.offset((ic + t6 - 1 as c_int) as isize) = *ch
                        .offset((i + t8 - 1 as c_int) as isize)
                        - *ch.offset((i + t9 - 1 as c_int) as isize);
                    *cc.offset((i + t7) as isize) =
                        *ch.offset((i + t8) as isize) + *ch.offset((i + t9) as isize);
                    *cc.offset((ic + t6) as isize) =
                        *ch.offset((i + t9) as isize) - *ch.offset((i + t8) as isize);
                    i += 2 as c_int
                }
                t6 += t10;
                t7 += t10;
                t8 += ido;
                t9 += ido;
                k += 1
            }
            j += 1
        }
        return;
    };
}
unsafe extern "C" fn drftf1(
    mut n: c_int,
    mut c: *mut c_float,
    mut ch: *mut c_float,
    mut wa: *mut c_float,
    mut ifac: *mut c_int,
) {
    let mut i: c_int = 0;
    let mut k1: c_int = 0;
    let mut l1: c_int = 0;
    let mut l2: c_int = 0;
    let mut na: c_int = 0;
    let mut kh: c_int = 0;
    let mut nf: c_int = 0;
    let mut ip: c_int = 0;
    let mut iw: c_int = 0;
    let mut ido: c_int = 0;
    let mut idl1: c_int = 0;
    let mut ix2: c_int = 0;
    let mut ix3: c_int = 0;
    nf = *ifac.offset(1 as c_int as isize);
    na = 1 as c_int;
    l2 = n;
    iw = n;
    k1 = 0 as c_int;
    while k1 < nf {
        kh = nf - k1;
        ip = *ifac.offset((kh + 1 as c_int) as isize);
        l1 = l2 / ip;
        ido = n / l2;
        idl1 = ido * l1;
        iw -= (ip - 1 as c_int) * ido;
        na = 1 as c_int - na;
        if ip != 4 as c_int {
            if ip != 2 as c_int {
                if ido == 1 as c_int {
                    na = 1 as c_int - na
                }
                if na != 0 as c_int {
                    dradfg(
                        ido,
                        ip,
                        l1,
                        idl1,
                        ch,
                        ch,
                        ch,
                        c,
                        c,
                        wa.offset(iw as isize).offset(-(1 as c_int as isize)),
                    );
                    na = 0 as c_int
                } else {
                    dradfg(
                        ido,
                        ip,
                        l1,
                        idl1,
                        c,
                        c,
                        c,
                        ch,
                        ch,
                        wa.offset(iw as isize).offset(-(1 as c_int as isize)),
                    );
                    na = 1 as c_int
                }
            } else if na != 0 as c_int {
                dradf2(
                    ido,
                    l1,
                    ch,
                    c,
                    wa.offset(iw as isize).offset(-(1 as c_int as isize)),
                );
            } else {
                dradf2(
                    ido,
                    l1,
                    c,
                    ch,
                    wa.offset(iw as isize).offset(-(1 as c_int as isize)),
                );
            }
        } else {
            ix2 = iw + ido;
            ix3 = ix2 + ido;
            if na != 0 as c_int {
                dradf4(
                    ido,
                    l1,
                    ch,
                    c,
                    wa.offset(iw as isize).offset(-(1 as c_int as isize)),
                    wa.offset(ix2 as isize).offset(-(1 as c_int as isize)),
                    wa.offset(ix3 as isize).offset(-(1 as c_int as isize)),
                );
            } else {
                dradf4(
                    ido,
                    l1,
                    c,
                    ch,
                    wa.offset(iw as isize).offset(-(1 as c_int as isize)),
                    wa.offset(ix2 as isize).offset(-(1 as c_int as isize)),
                    wa.offset(ix3 as isize).offset(-(1 as c_int as isize)),
                );
            }
        }
        l2 = l1;
        k1 += 1
    }
    if na == 1 as c_int {
        return;
    }
    i = 0 as c_int;
    while i < n {
        *c.offset(i as isize) = *ch.offset(i as isize);
        i += 1
    }
}
unsafe extern "C" fn dradb2(
    mut ido: c_int,
    mut l1: c_int,
    mut cc: *mut c_float,
    mut ch: *mut c_float,
    mut wa1: *mut c_float,
) {
    let mut i: c_int = 0;
    let mut k: c_int = 0;
    let mut t0: c_int = 0;
    let mut t1: c_int = 0;
    let mut t2: c_int = 0;
    let mut t3: c_int = 0;
    let mut t4: c_int = 0;
    let mut t5: c_int = 0;
    let mut t6: c_int = 0;
    let mut ti2: c_float = 0.;
    let mut tr2: c_float = 0.;
    t0 = l1 * ido;
    t1 = 0 as c_int;
    t2 = 0 as c_int;
    t3 = (ido << 1 as c_int) - 1 as c_int;
    k = 0 as c_int;
    while k < l1 {
        *ch.offset(t1 as isize) = *cc.offset(t2 as isize) + *cc.offset((t3 + t2) as isize);
        *ch.offset((t1 + t0) as isize) = *cc.offset(t2 as isize) - *cc.offset((t3 + t2) as isize);
        t1 += ido;
        t2 = t1 << 1 as c_int;
        k += 1
    }
    if ido < 2 as c_int {
        return;
    }
    if !(ido == 2 as c_int) {
        t1 = 0 as c_int;
        t2 = 0 as c_int;
        k = 0 as c_int;
        while k < l1 {
            t3 = t1;
            t4 = t2;
            t5 = t4 + (ido << 1 as c_int);
            t6 = t0 + t1;
            i = 2 as c_int;
            while i < ido {
                t3 += 2 as c_int;
                t4 += 2 as c_int;
                t5 -= 2 as c_int;
                t6 += 2 as c_int;
                *ch.offset((t3 - 1 as c_int) as isize) =
                    *cc.offset((t4 - 1 as c_int) as isize) + *cc.offset((t5 - 1 as c_int) as isize);
                tr2 =
                    *cc.offset((t4 - 1 as c_int) as isize) - *cc.offset((t5 - 1 as c_int) as isize);
                *ch.offset(t3 as isize) = *cc.offset(t4 as isize) - *cc.offset(t5 as isize);
                ti2 = *cc.offset(t4 as isize) + *cc.offset(t5 as isize);
                *ch.offset((t6 - 1 as c_int) as isize) = *wa1.offset((i - 2 as c_int) as isize)
                    * tr2
                    - *wa1.offset((i - 1 as c_int) as isize) * ti2;
                *ch.offset(t6 as isize) = *wa1.offset((i - 2 as c_int) as isize) * ti2
                    + *wa1.offset((i - 1 as c_int) as isize) * tr2;
                i += 2 as c_int
            }
            t1 += ido;
            t2 = t1 << 1 as c_int;
            k += 1
        }
        if ido % 2 as c_int == 1 as c_int {
            return;
        }
    }
    t1 = ido - 1 as c_int;
    t2 = ido - 1 as c_int;
    k = 0 as c_int;
    while k < l1 {
        *ch.offset(t1 as isize) = *cc.offset(t2 as isize) + *cc.offset(t2 as isize);
        *ch.offset((t1 + t0) as isize) =
            -(*cc.offset((t2 + 1 as c_int) as isize) + *cc.offset((t2 + 1 as c_int) as isize));
        t1 += ido;
        t2 += ido << 1 as c_int;
        k += 1
    }
}
unsafe extern "C" fn dradb3(
    mut ido: c_int,
    mut l1: c_int,
    mut cc: *mut c_float,
    mut ch: *mut c_float,
    mut wa1: *mut c_float,
    mut wa2: *mut c_float,
) {
    static mut taur: c_float = -0.5f32;
    static mut taui: c_float = 0.8660254037844386f32;
    let mut i: c_int = 0;
    let mut k: c_int = 0;
    let mut t0: c_int = 0;
    let mut t1: c_int = 0;
    let mut t2: c_int = 0;
    let mut t3: c_int = 0;
    let mut t4: c_int = 0;
    let mut t5: c_int = 0;
    let mut t6: c_int = 0;
    let mut t7: c_int = 0;
    let mut t8: c_int = 0;
    let mut t9: c_int = 0;
    let mut t10: c_int = 0;
    let mut ci2: c_float = 0.;
    let mut ci3: c_float = 0.;
    let mut di2: c_float = 0.;
    let mut di3: c_float = 0.;
    let mut cr2: c_float = 0.;
    let mut cr3: c_float = 0.;
    let mut dr2: c_float = 0.;
    let mut dr3: c_float = 0.;
    let mut ti2: c_float = 0.;
    let mut tr2: c_float = 0.;
    t0 = l1 * ido;
    t1 = 0 as c_int;
    t2 = t0 << 1 as c_int;
    t3 = ido << 1 as c_int;
    t4 = ido + (ido << 1 as c_int);
    t5 = 0 as c_int;
    k = 0 as c_int;
    while k < l1 {
        tr2 = *cc.offset((t3 - 1 as c_int) as isize) + *cc.offset((t3 - 1 as c_int) as isize);
        cr2 = *cc.offset(t5 as isize) + taur * tr2;
        *ch.offset(t1 as isize) = *cc.offset(t5 as isize) + tr2;
        ci3 = taui * (*cc.offset(t3 as isize) + *cc.offset(t3 as isize));
        *ch.offset((t1 + t0) as isize) = cr2 - ci3;
        *ch.offset((t1 + t2) as isize) = cr2 + ci3;
        t1 += ido;
        t3 += t4;
        t5 += t4;
        k += 1
    }
    if ido == 1 as c_int {
        return;
    }
    t1 = 0 as c_int;
    t3 = ido << 1 as c_int;
    k = 0 as c_int;
    while k < l1 {
        t7 = t1 + (t1 << 1 as c_int);
        t5 = t7 + t3;
        t6 = t5;
        t8 = t1;
        t9 = t1 + t0;
        t10 = t9 + t0;
        i = 2 as c_int;
        while i < ido {
            t5 += 2 as c_int;
            t6 -= 2 as c_int;
            t7 += 2 as c_int;
            t8 += 2 as c_int;
            t9 += 2 as c_int;
            t10 += 2 as c_int;
            tr2 = *cc.offset((t5 - 1 as c_int) as isize) + *cc.offset((t6 - 1 as c_int) as isize);
            cr2 = *cc.offset((t7 - 1 as c_int) as isize) + taur * tr2;
            *ch.offset((t8 - 1 as c_int) as isize) = *cc.offset((t7 - 1 as c_int) as isize) + tr2;
            ti2 = *cc.offset(t5 as isize) - *cc.offset(t6 as isize);
            ci2 = *cc.offset(t7 as isize) + taur * ti2;
            *ch.offset(t8 as isize) = *cc.offset(t7 as isize) + ti2;
            cr3 = taui
                * (*cc.offset((t5 - 1 as c_int) as isize) - *cc.offset((t6 - 1 as c_int) as isize));
            ci3 = taui * (*cc.offset(t5 as isize) + *cc.offset(t6 as isize));
            dr2 = cr2 - ci3;
            dr3 = cr2 + ci3;
            di2 = ci2 + cr3;
            di3 = ci2 - cr3;
            *ch.offset((t9 - 1 as c_int) as isize) = *wa1.offset((i - 2 as c_int) as isize) * dr2
                - *wa1.offset((i - 1 as c_int) as isize) * di2;
            *ch.offset(t9 as isize) = *wa1.offset((i - 2 as c_int) as isize) * di2
                + *wa1.offset((i - 1 as c_int) as isize) * dr2;
            *ch.offset((t10 - 1 as c_int) as isize) = *wa2.offset((i - 2 as c_int) as isize) * dr3
                - *wa2.offset((i - 1 as c_int) as isize) * di3;
            *ch.offset(t10 as isize) = *wa2.offset((i - 2 as c_int) as isize) * di3
                + *wa2.offset((i - 1 as c_int) as isize) * dr3;
            i += 2 as c_int
        }
        t1 += ido;
        k += 1
    }
}
unsafe extern "C" fn dradb4(
    mut ido: c_int,
    mut l1: c_int,
    mut cc: *mut c_float,
    mut ch: *mut c_float,
    mut wa1: *mut c_float,
    mut wa2: *mut c_float,
    mut wa3: *mut c_float,
) {
    static mut sqrt2: c_float = 1.414213562373095f32;
    let mut i: c_int = 0;
    let mut k: c_int = 0;
    let mut t0: c_int = 0;
    let mut t1: c_int = 0;
    let mut t2: c_int = 0;
    let mut t3: c_int = 0;
    let mut t4: c_int = 0;
    let mut t5: c_int = 0;
    let mut t6: c_int = 0;
    let mut t7: c_int = 0;
    let mut t8: c_int = 0;
    let mut ci2: c_float = 0.;
    let mut ci3: c_float = 0.;
    let mut ci4: c_float = 0.;
    let mut cr2: c_float = 0.;
    let mut cr3: c_float = 0.;
    let mut cr4: c_float = 0.;
    let mut ti1: c_float = 0.;
    let mut ti2: c_float = 0.;
    let mut ti3: c_float = 0.;
    let mut ti4: c_float = 0.;
    let mut tr1: c_float = 0.;
    let mut tr2: c_float = 0.;
    let mut tr3: c_float = 0.;
    let mut tr4: c_float = 0.;
    t0 = l1 * ido;
    t1 = 0 as c_int;
    t2 = ido << 2 as c_int;
    t3 = 0 as c_int;
    t6 = ido << 1 as c_int;
    k = 0 as c_int;
    while k < l1 {
        t4 = t3 + t6;
        t5 = t1;
        tr3 = *cc.offset((t4 - 1 as c_int) as isize) + *cc.offset((t4 - 1 as c_int) as isize);
        tr4 = *cc.offset(t4 as isize) + *cc.offset(t4 as isize);
        t4 += t6;
        tr1 = *cc.offset(t3 as isize) - *cc.offset((t4 - 1 as c_int) as isize);
        tr2 = *cc.offset(t3 as isize) + *cc.offset((t4 - 1 as c_int) as isize);
        *ch.offset(t5 as isize) = tr2 + tr3;
        t5 += t0;
        *ch.offset(t5 as isize) = tr1 - tr4;
        t5 += t0;
        *ch.offset(t5 as isize) = tr2 - tr3;
        t5 += t0;
        *ch.offset(t5 as isize) = tr1 + tr4;
        t1 += ido;
        t3 += t2;
        k += 1
    }
    if ido < 2 as c_int {
        return;
    }
    if !(ido == 2 as c_int) {
        t1 = 0 as c_int;
        k = 0 as c_int;
        while k < l1 {
            t2 = t1 << 2 as c_int;
            t3 = t2 + t6;
            t4 = t3;
            t5 = t4 + t6;
            t7 = t1;
            i = 2 as c_int;
            while i < ido {
                t2 += 2 as c_int;
                t3 += 2 as c_int;
                t4 -= 2 as c_int;
                t5 -= 2 as c_int;
                t7 += 2 as c_int;
                ti1 = *cc.offset(t2 as isize) + *cc.offset(t5 as isize);
                ti2 = *cc.offset(t2 as isize) - *cc.offset(t5 as isize);
                ti3 = *cc.offset(t3 as isize) - *cc.offset(t4 as isize);
                tr4 = *cc.offset(t3 as isize) + *cc.offset(t4 as isize);
                tr1 =
                    *cc.offset((t2 - 1 as c_int) as isize) - *cc.offset((t5 - 1 as c_int) as isize);
                tr2 =
                    *cc.offset((t2 - 1 as c_int) as isize) + *cc.offset((t5 - 1 as c_int) as isize);
                ti4 =
                    *cc.offset((t3 - 1 as c_int) as isize) - *cc.offset((t4 - 1 as c_int) as isize);
                tr3 =
                    *cc.offset((t3 - 1 as c_int) as isize) + *cc.offset((t4 - 1 as c_int) as isize);
                *ch.offset((t7 - 1 as c_int) as isize) = tr2 + tr3;
                cr3 = tr2 - tr3;
                *ch.offset(t7 as isize) = ti2 + ti3;
                ci3 = ti2 - ti3;
                cr2 = tr1 - tr4;
                cr4 = tr1 + tr4;
                ci2 = ti1 + ti4;
                ci4 = ti1 - ti4;
                t8 = t7 + t0;
                *ch.offset((t8 - 1 as c_int) as isize) = *wa1.offset((i - 2 as c_int) as isize)
                    * cr2
                    - *wa1.offset((i - 1 as c_int) as isize) * ci2;
                *ch.offset(t8 as isize) = *wa1.offset((i - 2 as c_int) as isize) * ci2
                    + *wa1.offset((i - 1 as c_int) as isize) * cr2;
                t8 += t0;
                *ch.offset((t8 - 1 as c_int) as isize) = *wa2.offset((i - 2 as c_int) as isize)
                    * cr3
                    - *wa2.offset((i - 1 as c_int) as isize) * ci3;
                *ch.offset(t8 as isize) = *wa2.offset((i - 2 as c_int) as isize) * ci3
                    + *wa2.offset((i - 1 as c_int) as isize) * cr3;
                t8 += t0;
                *ch.offset((t8 - 1 as c_int) as isize) = *wa3.offset((i - 2 as c_int) as isize)
                    * cr4
                    - *wa3.offset((i - 1 as c_int) as isize) * ci4;
                *ch.offset(t8 as isize) = *wa3.offset((i - 2 as c_int) as isize) * ci4
                    + *wa3.offset((i - 1 as c_int) as isize) * cr4;
                i += 2 as c_int
            }
            t1 += ido;
            k += 1
        }
        if ido % 2 as c_int == 1 as c_int {
            return;
        }
    }
    t1 = ido;
    t2 = ido << 2 as c_int;
    t3 = ido - 1 as c_int;
    t4 = ido + (ido << 1 as c_int);
    k = 0 as c_int;
    while k < l1 {
        t5 = t3;
        ti1 = *cc.offset(t1 as isize) + *cc.offset(t4 as isize);
        ti2 = *cc.offset(t4 as isize) - *cc.offset(t1 as isize);
        tr1 = *cc.offset((t1 - 1 as c_int) as isize) - *cc.offset((t4 - 1 as c_int) as isize);
        tr2 = *cc.offset((t1 - 1 as c_int) as isize) + *cc.offset((t4 - 1 as c_int) as isize);
        *ch.offset(t5 as isize) = tr2 + tr2;
        t5 += t0;
        *ch.offset(t5 as isize) = sqrt2 * (tr1 - ti1);
        t5 += t0;
        *ch.offset(t5 as isize) = ti2 + ti2;
        t5 += t0;
        *ch.offset(t5 as isize) = -sqrt2 * (tr1 + ti1);
        t3 += ido;
        t1 += t2;
        t4 += t2;
        k += 1
    }
}
unsafe extern "C" fn dradbg(
    mut ido: c_int,
    mut ip: c_int,
    mut l1: c_int,
    mut idl1: c_int,
    mut cc: *mut c_float,
    mut c1: *mut c_float,
    mut c2: *mut c_float,
    mut ch: *mut c_float,
    mut ch2: *mut c_float,
    mut wa: *mut c_float,
) {
    static mut tpi: c_float = 6.283185307179586f32;
    let mut idij: c_int = 0;
    let mut ipph: c_int = 0;
    let mut i: c_int = 0;
    let mut j: c_int = 0;
    let mut k: c_int = 0;
    let mut l: c_int = 0;
    let mut ik: c_int = 0;
    let mut is: c_int = 0;
    let mut t0: c_int = 0;
    let mut t1: c_int = 0;
    let mut t2: c_int = 0;
    let mut t3: c_int = 0;
    let mut t4: c_int = 0;
    let mut t5: c_int = 0;
    let mut t6: c_int = 0;
    let mut t7: c_int = 0;
    let mut t8: c_int = 0;
    let mut t9: c_int = 0;
    let mut t10: c_int = 0;
    let mut t11: c_int = 0;
    let mut t12: c_int = 0;
    let mut dc2: c_float = 0.;
    let mut ai1: c_float = 0.;
    let mut ai2: c_float = 0.;
    let mut ar1: c_float = 0.;
    let mut ar2: c_float = 0.;
    let mut ds2: c_float = 0.;
    let mut nbd: c_int = 0;
    let mut dcp: c_float = 0.;
    let mut arg: c_float = 0.;
    let mut dsp: c_float = 0.;
    let mut ar1h: c_float = 0.;
    let mut ar2h: c_float = 0.;
    let mut ipp2: c_int = 0;
    t10 = ip * ido;
    t0 = l1 * ido;
    arg = tpi / ip as c_float;
    dcp = cos(arg as c_double) as c_float;
    dsp = sin(arg as c_double) as c_float;
    nbd = ido - 1 as c_int >> 1 as c_int;
    ipp2 = ip;
    ipph = ip + 1 as c_int >> 1 as c_int;
    if ido < l1 {
        t1 = 0 as c_int;
        i = 0 as c_int;
        while i < ido {
            t2 = t1;
            t3 = t1;
            k = 0 as c_int;
            while k < l1 {
                *ch.offset(t2 as isize) = *cc.offset(t3 as isize);
                t2 += ido;
                t3 += t10;
                k += 1
            }
            t1 += 1;
            i += 1
        }
    } else {
        t1 = 0 as c_int;
        t2 = 0 as c_int;
        k = 0 as c_int;
        while k < l1 {
            t3 = t1;
            t4 = t2;
            i = 0 as c_int;
            while i < ido {
                *ch.offset(t3 as isize) = *cc.offset(t4 as isize);
                t3 += 1;
                t4 += 1;
                i += 1
            }
            t1 += ido;
            t2 += t10;
            k += 1
        }
    }
    t1 = 0 as c_int;
    t2 = ipp2 * t0;
    t5 = ido << 1 as c_int;
    t7 = t5;
    j = 1 as c_int;
    while j < ipph {
        t1 += t0;
        t2 -= t0;
        t3 = t1;
        t4 = t2;
        t6 = t5;
        k = 0 as c_int;
        while k < l1 {
            *ch.offset(t3 as isize) =
                *cc.offset((t6 - 1 as c_int) as isize) + *cc.offset((t6 - 1 as c_int) as isize);
            *ch.offset(t4 as isize) = *cc.offset(t6 as isize) + *cc.offset(t6 as isize);
            t3 += ido;
            t4 += ido;
            t6 += t10;
            k += 1
        }
        t5 += t7;
        j += 1
    }
    if !(ido == 1 as c_int) {
        if nbd < l1 {
            t1 = 0 as c_int;
            t2 = ipp2 * t0;
            t7 = 0 as c_int;
            j = 1 as c_int;
            while j < ipph {
                t1 += t0;
                t2 -= t0;
                t3 = t1;
                t4 = t2;
                t7 += ido << 1 as c_int;
                t8 = t7;
                t9 = t7;
                i = 2 as c_int;
                while i < ido {
                    t3 += 2 as c_int;
                    t4 += 2 as c_int;
                    t8 += 2 as c_int;
                    t9 -= 2 as c_int;
                    t5 = t3;
                    t6 = t4;
                    t11 = t8;
                    t12 = t9;
                    k = 0 as c_int;
                    while k < l1 {
                        *ch.offset((t5 - 1 as c_int) as isize) = *cc
                            .offset((t11 - 1 as c_int) as isize)
                            + *cc.offset((t12 - 1 as c_int) as isize);
                        *ch.offset((t6 - 1 as c_int) as isize) = *cc
                            .offset((t11 - 1 as c_int) as isize)
                            - *cc.offset((t12 - 1 as c_int) as isize);
                        *ch.offset(t5 as isize) =
                            *cc.offset(t11 as isize) - *cc.offset(t12 as isize);
                        *ch.offset(t6 as isize) =
                            *cc.offset(t11 as isize) + *cc.offset(t12 as isize);
                        t5 += ido;
                        t6 += ido;
                        t11 += t10;
                        t12 += t10;
                        k += 1
                    }
                    i += 2 as c_int
                }
                j += 1
            }
        } else {
            t1 = 0 as c_int;
            t2 = ipp2 * t0;
            t7 = 0 as c_int;
            j = 1 as c_int;
            while j < ipph {
                t1 += t0;
                t2 -= t0;
                t3 = t1;
                t4 = t2;
                t7 += ido << 1 as c_int;
                t8 = t7;
                k = 0 as c_int;
                while k < l1 {
                    t5 = t3;
                    t6 = t4;
                    t9 = t8;
                    t11 = t8;
                    i = 2 as c_int;
                    while i < ido {
                        t5 += 2 as c_int;
                        t6 += 2 as c_int;
                        t9 += 2 as c_int;
                        t11 -= 2 as c_int;
                        *ch.offset((t5 - 1 as c_int) as isize) = *cc
                            .offset((t9 - 1 as c_int) as isize)
                            + *cc.offset((t11 - 1 as c_int) as isize);
                        *ch.offset((t6 - 1 as c_int) as isize) = *cc
                            .offset((t9 - 1 as c_int) as isize)
                            - *cc.offset((t11 - 1 as c_int) as isize);
                        *ch.offset(t5 as isize) =
                            *cc.offset(t9 as isize) - *cc.offset(t11 as isize);
                        *ch.offset(t6 as isize) =
                            *cc.offset(t9 as isize) + *cc.offset(t11 as isize);
                        i += 2 as c_int
                    }
                    t3 += ido;
                    t4 += ido;
                    t8 += t10;
                    k += 1
                }
                j += 1
            }
        }
    }
    ar1 = 1.0f32;
    ai1 = 0.0f32;
    t1 = 0 as c_int;
    t2 = ipp2 * idl1;
    t9 = t2;
    t3 = (ip - 1 as c_int) * idl1;
    l = 1 as c_int;
    while l < ipph {
        t1 += idl1;
        t2 -= idl1;
        ar1h = dcp * ar1 - dsp * ai1;
        ai1 = dcp * ai1 + dsp * ar1;
        ar1 = ar1h;
        t4 = t1;
        t5 = t2;
        t6 = 0 as c_int;
        t7 = idl1;
        t8 = t3;
        ik = 0 as c_int;
        while ik < idl1 {
            let fresh13 = t6;
            t6 = t6 + 1;
            let fresh14 = t7;
            t7 = t7 + 1;
            let fresh15 = t4;
            t4 = t4 + 1;
            *c2.offset(fresh15 as isize) =
                *ch2.offset(fresh13 as isize) + ar1 * *ch2.offset(fresh14 as isize);
            let fresh16 = t8;
            t8 = t8 + 1;
            let fresh17 = t5;
            t5 = t5 + 1;
            *c2.offset(fresh17 as isize) = ai1 * *ch2.offset(fresh16 as isize);
            ik += 1
        }
        dc2 = ar1;
        ds2 = ai1;
        ar2 = ar1;
        ai2 = ai1;
        t6 = idl1;
        t7 = t9 - idl1;
        j = 2 as c_int;
        while j < ipph {
            t6 += idl1;
            t7 -= idl1;
            ar2h = dc2 * ar2 - ds2 * ai2;
            ai2 = dc2 * ai2 + ds2 * ar2;
            ar2 = ar2h;
            t4 = t1;
            t5 = t2;
            t11 = t6;
            t12 = t7;
            ik = 0 as c_int;
            while ik < idl1 {
                let fresh18 = t11;
                t11 = t11 + 1;
                let fresh19 = t4;
                t4 = t4 + 1;
                *c2.offset(fresh19 as isize) += ar2 * *ch2.offset(fresh18 as isize);
                let fresh20 = t12;
                t12 = t12 + 1;
                let fresh21 = t5;
                t5 = t5 + 1;
                *c2.offset(fresh21 as isize) += ai2 * *ch2.offset(fresh20 as isize);
                ik += 1
            }
            j += 1
        }
        l += 1
    }
    t1 = 0 as c_int;
    j = 1 as c_int;
    while j < ipph {
        t1 += idl1;
        t2 = t1;
        ik = 0 as c_int;
        while ik < idl1 {
            let fresh22 = t2;
            t2 = t2 + 1;
            *ch2.offset(ik as isize) += *ch2.offset(fresh22 as isize);
            ik += 1
        }
        j += 1
    }
    t1 = 0 as c_int;
    t2 = ipp2 * t0;
    j = 1 as c_int;
    while j < ipph {
        t1 += t0;
        t2 -= t0;
        t3 = t1;
        t4 = t2;
        k = 0 as c_int;
        while k < l1 {
            *ch.offset(t3 as isize) = *c1.offset(t3 as isize) - *c1.offset(t4 as isize);
            *ch.offset(t4 as isize) = *c1.offset(t3 as isize) + *c1.offset(t4 as isize);
            t3 += ido;
            t4 += ido;
            k += 1
        }
        j += 1
    }
    if !(ido == 1 as c_int) {
        if nbd < l1 {
            t1 = 0 as c_int;
            t2 = ipp2 * t0;
            j = 1 as c_int;
            while j < ipph {
                t1 += t0;
                t2 -= t0;
                t3 = t1;
                t4 = t2;
                i = 2 as c_int;
                while i < ido {
                    t3 += 2 as c_int;
                    t4 += 2 as c_int;
                    t5 = t3;
                    t6 = t4;
                    k = 0 as c_int;
                    while k < l1 {
                        *ch.offset((t5 - 1 as c_int) as isize) =
                            *c1.offset((t5 - 1 as c_int) as isize) - *c1.offset(t6 as isize);
                        *ch.offset((t6 - 1 as c_int) as isize) =
                            *c1.offset((t5 - 1 as c_int) as isize) + *c1.offset(t6 as isize);
                        *ch.offset(t5 as isize) =
                            *c1.offset(t5 as isize) + *c1.offset((t6 - 1 as c_int) as isize);
                        *ch.offset(t6 as isize) =
                            *c1.offset(t5 as isize) - *c1.offset((t6 - 1 as c_int) as isize);
                        t5 += ido;
                        t6 += ido;
                        k += 1
                    }
                    i += 2 as c_int
                }
                j += 1
            }
        } else {
            t1 = 0 as c_int;
            t2 = ipp2 * t0;
            j = 1 as c_int;
            while j < ipph {
                t1 += t0;
                t2 -= t0;
                t3 = t1;
                t4 = t2;
                k = 0 as c_int;
                while k < l1 {
                    t5 = t3;
                    t6 = t4;
                    i = 2 as c_int;
                    while i < ido {
                        t5 += 2 as c_int;
                        t6 += 2 as c_int;
                        *ch.offset((t5 - 1 as c_int) as isize) =
                            *c1.offset((t5 - 1 as c_int) as isize) - *c1.offset(t6 as isize);
                        *ch.offset((t6 - 1 as c_int) as isize) =
                            *c1.offset((t5 - 1 as c_int) as isize) + *c1.offset(t6 as isize);
                        *ch.offset(t5 as isize) =
                            *c1.offset(t5 as isize) + *c1.offset((t6 - 1 as c_int) as isize);
                        *ch.offset(t6 as isize) =
                            *c1.offset(t5 as isize) - *c1.offset((t6 - 1 as c_int) as isize);
                        i += 2 as c_int
                    }
                    t3 += ido;
                    t4 += ido;
                    k += 1
                }
                j += 1
            }
        }
    }
    if ido == 1 as c_int {
        return;
    }
    ik = 0 as c_int;
    while ik < idl1 {
        *c2.offset(ik as isize) = *ch2.offset(ik as isize);
        ik += 1
    }
    t1 = 0 as c_int;
    j = 1 as c_int;
    while j < ip {
        t1 += t0;
        t2 = t1;
        k = 0 as c_int;
        while k < l1 {
            *c1.offset(t2 as isize) = *ch.offset(t2 as isize);
            t2 += ido;
            k += 1
        }
        j += 1
    }
    if nbd > l1 {
        is = -ido - 1 as c_int;
        t1 = 0 as c_int;
        j = 1 as c_int;
        while j < ip {
            is += ido;
            t1 += t0;
            t2 = t1;
            k = 0 as c_int;
            while k < l1 {
                idij = is;
                t3 = t2;
                i = 2 as c_int;
                while i < ido {
                    idij += 2 as c_int;
                    t3 += 2 as c_int;
                    *c1.offset((t3 - 1 as c_int) as isize) = *wa
                        .offset((idij - 1 as c_int) as isize)
                        * *ch.offset((t3 - 1 as c_int) as isize)
                        - *wa.offset(idij as isize) * *ch.offset(t3 as isize);
                    *c1.offset(t3 as isize) = *wa.offset((idij - 1 as c_int) as isize)
                        * *ch.offset(t3 as isize)
                        + *wa.offset(idij as isize) * *ch.offset((t3 - 1 as c_int) as isize);
                    i += 2 as c_int
                }
                t2 += ido;
                k += 1
            }
            j += 1
        }
        return;
    } else {
        is = -ido - 1 as c_int;
        t1 = 0 as c_int;
        j = 1 as c_int;
        while j < ip {
            is += ido;
            t1 += t0;
            idij = is;
            t2 = t1;
            i = 2 as c_int;
            while i < ido {
                t2 += 2 as c_int;
                idij += 2 as c_int;
                t3 = t2;
                k = 0 as c_int;
                while k < l1 {
                    *c1.offset((t3 - 1 as c_int) as isize) = *wa
                        .offset((idij - 1 as c_int) as isize)
                        * *ch.offset((t3 - 1 as c_int) as isize)
                        - *wa.offset(idij as isize) * *ch.offset(t3 as isize);
                    *c1.offset(t3 as isize) = *wa.offset((idij - 1 as c_int) as isize)
                        * *ch.offset(t3 as isize)
                        + *wa.offset(idij as isize) * *ch.offset((t3 - 1 as c_int) as isize);
                    t3 += ido;
                    k += 1
                }
                i += 2 as c_int
            }
            j += 1
        }
        return;
    };
}
unsafe extern "C" fn drftb1(
    mut n: c_int,
    mut c: *mut c_float,
    mut ch: *mut c_float,
    mut wa: *mut c_float,
    mut ifac: *mut c_int,
) {
    let mut i: c_int = 0;
    let mut k1: c_int = 0;
    let mut l1: c_int = 0;
    let mut l2: c_int = 0;
    let mut na: c_int = 0;
    let mut nf: c_int = 0;
    let mut ip: c_int = 0;
    let mut iw: c_int = 0;
    let mut ix2: c_int = 0;
    let mut ix3: c_int = 0;
    let mut ido: c_int = 0;
    let mut idl1: c_int = 0;
    nf = *ifac.offset(1 as c_int as isize);
    na = 0 as c_int;
    l1 = 1 as c_int;
    iw = 1 as c_int;
    k1 = 0 as c_int;
    while k1 < nf {
        ip = *ifac.offset((k1 + 2 as c_int) as isize);
        l2 = ip * l1;
        ido = n / l2;
        idl1 = ido * l1;
        if ip != 4 as c_int {
            if ip != 2 as c_int {
                if ip != 3 as c_int {
                    /*    The radix five case can be translated later..... */
                    /*    if(ip!=5)goto L112;

                      ix2=iw+ido;
                      ix3=ix2+ido;
                      ix4=ix3+ido;
                      if(na!=0)
                        dradb5(ido,l1,ch,c,wa+iw-1,wa+ix2-1,wa+ix3-1,wa+ix4-1);
                      else
                        dradb5(ido,l1,c,ch,wa+iw-1,wa+ix2-1,wa+ix3-1,wa+ix4-1);
                      na=1-na;
                      goto L115;

                    L112:*/
                    if na != 0 as c_int {
                        dradbg(
                            ido,
                            ip,
                            l1,
                            idl1,
                            ch,
                            ch,
                            ch,
                            c,
                            c,
                            wa.offset(iw as isize).offset(-(1 as c_int as isize)),
                        );
                    } else {
                        dradbg(
                            ido,
                            ip,
                            l1,
                            idl1,
                            c,
                            c,
                            c,
                            ch,
                            ch,
                            wa.offset(iw as isize).offset(-(1 as c_int as isize)),
                        );
                    }
                    if ido == 1 as c_int {
                        na = 1 as c_int - na
                    }
                } else {
                    ix2 = iw + ido;
                    if na != 0 as c_int {
                        dradb3(
                            ido,
                            l1,
                            ch,
                            c,
                            wa.offset(iw as isize).offset(-(1 as c_int as isize)),
                            wa.offset(ix2 as isize).offset(-(1 as c_int as isize)),
                        );
                    } else {
                        dradb3(
                            ido,
                            l1,
                            c,
                            ch,
                            wa.offset(iw as isize).offset(-(1 as c_int as isize)),
                            wa.offset(ix2 as isize).offset(-(1 as c_int as isize)),
                        );
                    }
                    na = 1 as c_int - na
                }
            } else {
                if na != 0 as c_int {
                    dradb2(
                        ido,
                        l1,
                        ch,
                        c,
                        wa.offset(iw as isize).offset(-(1 as c_int as isize)),
                    );
                } else {
                    dradb2(
                        ido,
                        l1,
                        c,
                        ch,
                        wa.offset(iw as isize).offset(-(1 as c_int as isize)),
                    );
                }
                na = 1 as c_int - na
            }
        } else {
            ix2 = iw + ido;
            ix3 = ix2 + ido;
            if na != 0 as c_int {
                dradb4(
                    ido,
                    l1,
                    ch,
                    c,
                    wa.offset(iw as isize).offset(-(1 as c_int as isize)),
                    wa.offset(ix2 as isize).offset(-(1 as c_int as isize)),
                    wa.offset(ix3 as isize).offset(-(1 as c_int as isize)),
                );
            } else {
                dradb4(
                    ido,
                    l1,
                    c,
                    ch,
                    wa.offset(iw as isize).offset(-(1 as c_int as isize)),
                    wa.offset(ix2 as isize).offset(-(1 as c_int as isize)),
                    wa.offset(ix3 as isize).offset(-(1 as c_int as isize)),
                );
            }
            na = 1 as c_int - na
        }
        l1 = l2;
        iw += (ip - 1 as c_int) * ido;
        k1 += 1
    }
    if na == 0 as c_int {
        return;
    }
    i = 0 as c_int;
    while i < n {
        *c.offset(i as isize) = *ch.offset(i as isize);
        i += 1
    }
}
#[no_mangle]
pub unsafe extern "C" fn spx_drft_forward(mut l: *mut drft_lookup, mut data: *mut c_float) {
    if (*l).n == 1 as c_int {
        return;
    }
    drftf1(
        (*l).n,
        data,
        (*l).trigcache,
        (*l).trigcache.offset((*l).n as isize),
        (*l).splitcache,
    );
}
#[no_mangle]
pub unsafe extern "C" fn spx_drft_backward(mut l: *mut drft_lookup, mut data: *mut c_float) {
    if (*l).n == 1 as c_int {
        return;
    }
    drftb1(
        (*l).n,
        data,
        (*l).trigcache,
        (*l).trigcache.offset((*l).n as isize),
        (*l).splitcache,
    );
}
#[no_mangle]
pub unsafe extern "C" fn spx_drft_init(mut l: *mut drft_lookup, mut n: c_int) {
    (*l).n = n;
    (*l).trigcache = speex_alloc(
        ((3 as c_int * n) as c_ulong).wrapping_mul(::std::mem::size_of::<c_float>() as c_ulong)
            as c_int,
    ) as *mut c_float;
    (*l).splitcache = speex_alloc(
        (32 as c_int as c_ulong).wrapping_mul(::std::mem::size_of::<c_int>() as c_ulong) as c_int,
    ) as *mut c_int;
    fdrffti(n, (*l).trigcache, (*l).splitcache);
}
#[no_mangle]
pub unsafe extern "C" fn spx_drft_clear(mut l: *mut drft_lookup) {
    if !l.is_null() {
        if !(*l).trigcache.is_null() {
            speex_free((*l).trigcache as *mut c_void);
        }
        if !(*l).splitcache.is_null() {
            speex_free((*l).splitcache as *mut c_void);
        }
    };
}

#[cfg(test)]
mod tests {
    use std::os::raw::{c_float, c_int};

    #[test]
    fn fdrffti_simple() {
        let mut trigcache = [42. as c_float; 3];
        let mut splitcache = [24 as c_int; 32];

        unsafe {
            super::fdrffti(1, trigcache.as_mut_ptr(), splitcache.as_mut_ptr());
        }
        assert!(trigcache.iter().all(|&x| x == 42.));
        assert!(splitcache.iter().all(|&x| x == 24));
    }

    #[test]
    fn fdrffti() {
        const SIZE: usize = 1024;
        #[rustfmt::skip]
        const TRIGCACHE: [c_float; 1018] = [
            0.99998116, 0.0061358847, 0.9999247, 0.012271538, 0.9998306, 0.01840673, 0.9996988,
            0.024541229, 0.9995294, 0.030674804, 0.99932235, 0.036807224, 0.99907774, 0.04293826,
            0.99879545, 0.049067676, 0.99847555, 0.055195246, 0.9981181, 0.06132074, 0.99772304,
            0.06744392, 0.99729043, 0.07356457, 0.9968203, 0.07968244, 0.9963126, 0.08579732,
            0.9957674, 0.091908954, 0.9951847, 0.09801714, 0.9945646, 0.10412164, 0.993907,
            0.110222206, 0.9932119, 0.116318636, 0.99247956, 0.12241068, 0.99170977, 0.12849812,
            0.99090266, 0.13458072, 0.9900582, 0.14065824, 0.9891765, 0.14673047, 0.9882576,
            0.15279719, 0.9873014, 0.15885815, 0.9863081, 0.16491313, 0.98527765, 0.1709619,
            0.9842101, 0.17700422, 0.9831055, 0.18303989, 0.9819639, 0.18906866, 0.98078525,
            0.19509032, 0.9795698, 0.20110464, 0.9783174, 0.20711139, 0.97702813, 0.21311033,
            0.9757021, 0.21910124, 0.97433937, 0.22508392, 0.97293997, 0.23105812, 0.9715039,
            0.2370236, 0.97003126, 0.2429802, 0.9685221, 0.24892761, 0.96697646, 0.25486568,
            0.96539444, 0.26079413, 0.96377605, 0.26671278, 0.9621214, 0.27262136, 0.9604305,
            0.2785197, 0.95870346, 0.28440756, 0.95694035, 0.29028466, 0.9551412, 0.2961509,
            0.953306, 0.30200595, 0.951435, 0.30784968, 0.94952816, 0.31368175, 0.9475856,
            0.31950203, 0.9456073, 0.32531032, 0.94359344, 0.3311063, 0.94154406, 0.33688986,
            0.9394592, 0.34266073, 0.937339, 0.34841868, 0.9351835, 0.35416353, 0.9329928,
            0.35989505, 0.93076694, 0.365613, 0.9285061, 0.3713172, 0.9262102, 0.37700742,
            0.9238795, 0.38268346, 0.92151403, 0.38834503, 0.9191139, 0.39399207, 0.9166791,
            0.3996242, 0.9142097, 0.40524134, 0.91170603, 0.4108432, 0.90916795, 0.41642958,
            0.9065957, 0.4220003, 0.9039893, 0.42755508, 0.9013488, 0.43309385, 0.8986745,
            0.43861625, 0.89596623, 0.44412217, 0.8932243, 0.44961134, 0.89044875, 0.45508358,
            0.88763964, 0.46053872, 0.8847971, 0.4659765, 0.88192123, 0.47139674, 0.8790122,
            0.47679925, 0.8760701, 0.48218375, 0.873095, 0.48755017, 0.87008697, 0.49289823,
            0.86704624, 0.49822766, 0.86397284, 0.50353837, 0.8608669, 0.5088302, 0.8577286,
            0.51410276, 0.854558, 0.519356, 0.8513552, 0.5245897, 0.84812033, 0.52980363,
            0.8448536, 0.53499764, 0.841555, 0.5401715, 0.8382247, 0.545325, 0.8348628, 0.550458,
            0.8314696, 0.55557024, 0.828045, 0.5606616, 0.82458925, 0.5657318, 0.8211025,
            0.57078075, 0.8175848, 0.5758082, 0.8140363, 0.580814, 0.81045717, 0.5857979,
            0.8068476, 0.5907597, 0.8032075, 0.5956993, 0.79953724, 0.6006165, 0.7958369,
            0.605511, 0.79210657, 0.6103828, 0.7883464, 0.61523163, 0.7845566, 0.6200572,
            0.7807372, 0.6248595, 0.77688843, 0.62963825, 0.77301043, 0.63439333, 0.76910335,
            0.63912445, 0.76516724, 0.64383155, 0.76120234, 0.64851445, 0.7572088, 0.65317285,
            0.7531868, 0.6578067, 0.7491364, 0.6624158, 0.74505776, 0.66699994, 0.7409511,
            0.671559, 0.7368166, 0.67609274, 0.7326543, 0.680601, 0.72846437, 0.6850837,
            0.7242471, 0.68954057, 0.7200025, 0.69397146, 0.7157308, 0.6983763, 0.7114322,
            0.70275474, 0.0, 0.0, 0.9999247, 0.012271538, 0.9996988, 0.024541229, 0.99932235,
            0.036807224, 0.99879545, 0.049067676, 0.9981181, 0.06132074, 0.99729043, 0.07356457,
            0.9963126, 0.08579732, 0.9951847, 0.09801714, 0.993907, 0.110222206, 0.99247956,
            0.12241068, 0.99090266, 0.13458072, 0.9891765, 0.14673047, 0.9873014, 0.15885815,
            0.98527765, 0.1709619, 0.9831055, 0.18303989, 0.98078525, 0.19509032, 0.9783174,
            0.20711139, 0.9757021, 0.21910124, 0.97293997, 0.23105812, 0.97003126, 0.2429802,
            0.96697646, 0.25486568, 0.96377605, 0.26671278, 0.9604305, 0.2785197, 0.95694035,
            0.29028466, 0.953306, 0.30200595, 0.94952816, 0.31368175, 0.9456073, 0.32531032,
            0.94154406, 0.33688986, 0.937339, 0.34841868, 0.9329928, 0.35989505, 0.9285061,
            0.3713172, 0.9238795, 0.38268346, 0.9191139, 0.39399207, 0.9142097, 0.40524134,
            0.90916795, 0.41642958, 0.9039893, 0.42755508, 0.8986745, 0.43861625, 0.8932243,
            0.44961134, 0.88763964, 0.46053872, 0.88192123, 0.47139674, 0.8760701, 0.48218375,
            0.87008697, 0.49289823, 0.86397284, 0.50353837, 0.8577286, 0.51410276, 0.8513552,
            0.5245897, 0.8448536, 0.53499764, 0.8382247, 0.545325, 0.8314696, 0.55557024,
            0.82458925, 0.5657318, 0.8175848, 0.5758082, 0.81045717, 0.5857979, 0.8032075,
            0.5956993, 0.7958369, 0.605511, 0.7883464, 0.61523163, 0.7807372, 0.6248595,
            0.77301043, 0.63439333, 0.76516724, 0.64383155, 0.7572088, 0.65317285, 0.7491364,
            0.6624158, 0.7409511, 0.671559, 0.7326543, 0.680601, 0.7242471, 0.68954057, 0.7157308,
            0.6983763, 0.70710677, 0.70710677, 0.69837624, 0.71573085, 0.6895405, 0.7242471,
            0.680601, 0.7326543, 0.6715589, 0.7409512, 0.66241574, 0.7491364, 0.6531728,
            0.7572089, 0.64383155, 0.7651673, 0.6343933, 0.77301043, 0.62485945, 0.7807373,
            0.6152316, 0.7883464, 0.605511, 0.7958369, 0.5956993, 0.8032075, 0.58579785,
            0.81045717, 0.57580817, 0.8175848, 0.5657318, 0.8245893, 0.5555702, 0.83146966,
            0.545325, 0.8382247, 0.53499764, 0.8448536, 0.52458966, 0.8513552, 0.5141027,
            0.85772866, 0.5035384, 0.86397284, 0.4928982, 0.87008697, 0.48218372, 0.87607014,
            0.47139665, 0.8819213, 0.46053872, 0.88763964, 0.4496113, 0.8932243, 0.4386162,
            0.8986745, 0.4275551, 0.9039893, 0.41642955, 0.909168, 0.40524128, 0.9142098,
            0.39399195, 0.9191139, 0.38268343, 0.9238795, 0.37131715, 0.9285061, 0.35989496,
            0.9329928, 0.34841868, 0.937339, 0.33688983, 0.94154406, 0.32531023, 0.94560736,
            0.31368166, 0.9495282, 0.30200595, 0.953306, 0.29028463, 0.95694035, 0.2785196,
            0.96043056, 0.26671275, 0.96377605, 0.25486565, 0.96697646, 0.24298012, 0.97003126,
            0.23105814, 0.97293997, 0.21910122, 0.9757021, 0.20711133, 0.9783174, 0.19509023,
            0.9807853, 0.18303989, 0.9831055, 0.17096186, 0.98527765, 0.15885808, 0.9873014,
            0.1467305, 0.9891765, 0.13458069, 0.99090266, 0.122410625, 0.99247956, 0.110222116,
            0.993907, 0.098017134, 0.9951847, 0.08579727, 0.9963126, 0.07356449, 0.99729043,
            0.06132075, 0.9981181, 0.04906765, 0.99879545, 0.036807165, 0.9993224, 0.024541136,
            0.9996988, 0.012271529, 0.9999247, 0.0, 0.0, 0.9998306, 0.01840673, 0.99932235,
            0.036807224, 0.99847555, 0.055195246, 0.99729043, 0.07356457, 0.9957674, 0.091908954,
            0.993907, 0.110222206, 0.99170977, 0.12849812, 0.9891765, 0.14673047, 0.9863081,
            0.16491312, 0.9831055, 0.18303989, 0.9795698, 0.20110464, 0.9757021, 0.21910124,
            0.9715039, 0.23702359, 0.96697646, 0.25486568, 0.9621214, 0.27262136, 0.95694035,
            0.29028466, 0.951435, 0.30784965, 0.9456073, 0.3253103, 0.9394592, 0.34266073,
            0.9329928, 0.35989505, 0.9262102, 0.3770074, 0.9191139, 0.39399207, 0.91170603,
            0.4108432, 0.9039893, 0.42755508, 0.89596623, 0.44412214, 0.88763964, 0.4605387,
            0.8790122, 0.47679925, 0.87008697, 0.49289823, 0.8608669, 0.5088302, 0.8513552,
            0.5245897, 0.841555, 0.5401715, 0.8314696, 0.55557024, 0.8211025, 0.57078075,
            0.8104572, 0.58579785, 0.7995373, 0.60061646, 0.78834647, 0.6152316, 0.77688843,
            0.62963825, 0.76516724, 0.64383155, 0.7531868, 0.6578067, 0.7409511, 0.671559,
            0.72846437, 0.6850837, 0.71573085, 0.69837624, 0.70275474, 0.71143216, 0.6895405,
            0.7242471, 0.6760927, 0.7368166, 0.66241574, 0.7491364, 0.6485144, 0.7612024,
            0.6343933, 0.77301043, 0.6200572, 0.78455657, 0.60551107, 0.7958369, 0.5907597,
            0.8068475, 0.5758082, 0.8175848, 0.56066155, 0.82804507, 0.545325, 0.8382247,
            0.52980363, 0.84812033, 0.5141027, 0.85772866, 0.4982277, 0.86704624, 0.48218372,
            0.87607014, 0.4659765, 0.8847971, 0.4496113, 0.8932243, 0.43309385, 0.9013488,
            0.41642955, 0.909168, 0.39962426, 0.916679, 0.38268343, 0.9238795, 0.36561295,
            0.930767, 0.34841868, 0.937339, 0.33110628, 0.94359344, 0.31368175, 0.94952816,
            0.29615086, 0.9551412, 0.27851972, 0.9604305, 0.2607941, 0.96539444, 0.24298024,
            0.97003126, 0.22508392, 0.97433937, 0.20711133, 0.9783174, 0.18906869, 0.9819639,
            0.17096186, 0.98527765, 0.15279722, 0.9882576, 0.13458069, 0.99090266, 0.11631868,
            0.9932119, 0.098017134, 0.9951847, 0.07968238, 0.99682033, 0.06132075, 0.9981181,
            0.042938218, 0.99907774, 0.024541255, 0.9996988, 0.006135858, 0.99998116,
            -0.012271497, 0.9999247, -0.030674815, 0.9995294, -0.04906774, 0.99879545,
            -0.067443915, 0.99772304, -0.08579736, 0.9963126, -0.10412162, 0.9945646, -0.12241071,
            0.9924795, -0.1406582, 0.99005824, -0.15885817, 0.9873014, -0.17700417, 0.9842101,
            -0.19509032, 0.98078525, -0.21311037, 0.97702813, -0.23105809, 0.97293997,
            -0.24892765, 0.9685221, -0.26671273, 0.96377605, -0.28440756, 0.95870346, -0.30200592,
            0.9533061, -0.31950206, 0.9475856, -0.3368898, 0.94154406, -0.35416353, 0.9351835,
            -0.37131724, 0.9285061, -0.38834503, 0.92151403, -0.40524134, 0.9142097, -0.42200023,
            0.9065957, -0.43861625, 0.8986744, -0.45508364, 0.8904487, -0.47139683, 0.88192123,
            -0.4875501, 0.87309504, -0.50353837, 0.86397284, -0.519356, 0.854558, -0.5349977,
            0.8448535, -0.5504579, 0.83486295, -0.56573176, 0.8245893, -0.58081394, 0.8140363,
            -0.59569937, 0.8032075, -0.6103829, 0.7921065, -0.62485945, 0.7807373, -0.63912445,
            0.76910335, -0.65317285, 0.7572088, -0.667, 0.7450577, -0.68060094, 0.73265433,
            -0.69397146, 0.72000253, 0.0, 0.0, 0.9996988, 0.024541229, 0.99879545, 0.049067676,
            0.99729043, 0.07356457, 0.9951847, 0.09801714, 0.99247956, 0.12241068, 0.9891765,
            0.14673047, 0.98527765, 0.1709619, 0.98078525, 0.19509032, 0.9757021, 0.21910124,
            0.97003126, 0.2429802, 0.96377605, 0.26671278, 0.95694035, 0.29028466, 0.94952816,
            0.31368175, 0.94154406, 0.33688986, 0.9329928, 0.35989505, 0.9238795, 0.38268346,
            0.9142097, 0.40524134, 0.9039893, 0.42755508, 0.8932243, 0.44961134, 0.88192123,
            0.47139674, 0.87008697, 0.49289823, 0.8577286, 0.51410276, 0.8448536, 0.53499764,
            0.8314696, 0.55557024, 0.8175848, 0.5758082, 0.8032075, 0.5956993, 0.7883464,
            0.61523163, 0.77301043, 0.63439333, 0.7572088, 0.65317285, 0.7409511, 0.671559,
            0.7242471, 0.68954057, 0.0, 0.0, 0.99879545, 0.049067676, 0.9951847, 0.09801714,
            0.9891765, 0.14673047, 0.98078525, 0.19509032, 0.97003126, 0.2429802, 0.95694035,
            0.29028466, 0.94154406, 0.33688986, 0.9238795, 0.38268346, 0.9039893, 0.42755508,
            0.88192123, 0.47139674, 0.8577286, 0.51410276, 0.8314696, 0.55557024, 0.8032075,
            0.5956993, 0.77301043, 0.63439333, 0.7409511, 0.671559, 0.70710677, 0.70710677,
            0.6715589, 0.7409512, 0.6343933, 0.77301043, 0.5956993, 0.8032075, 0.5555702,
            0.83146966, 0.5141027, 0.85772866, 0.47139665, 0.8819213, 0.4275551, 0.9039893,
            0.38268343, 0.9238795, 0.33688983, 0.94154406, 0.29028463, 0.95694035, 0.24298012,
            0.97003126, 0.19509023, 0.9807853, 0.1467305, 0.9891765, 0.098017134, 0.9951847,
            0.04906765, 0.99879545, 0.0, 0.0, 0.99729043, 0.07356457, 0.9891765, 0.14673047,
            0.9757021, 0.21910124, 0.95694035, 0.29028466, 0.9329928, 0.35989505, 0.9039893,
            0.42755508, 0.87008697, 0.49289823, 0.8314696, 0.55557024, 0.78834647, 0.6152316,
            0.7409511, 0.671559, 0.6895405, 0.7242471, 0.6343933, 0.77301043, 0.5758082,
            0.8175848, 0.5141027, 0.85772866, 0.4496113, 0.8932243, 0.38268343, 0.9238795,
            0.31368175, 0.94952816, 0.24298024, 0.97003126, 0.17096186, 0.98527765, 0.098017134,
            0.9951847, 0.024541255, 0.9996988, -0.04906774, 0.99879545, -0.12241071, 0.9924795,
            -0.19509032, 0.98078525, -0.26671273, 0.96377605, -0.3368898, 0.94154406, -0.40524134,
            0.9142097, -0.47139683, 0.88192123, -0.5349977, 0.8448535, -0.59569937, 0.8032075,
            -0.65317285, 0.7572088, 0.0, 0.0, 0.9951847, 0.09801714, 0.98078525, 0.19509032,
            0.95694035, 0.29028466, 0.9238795, 0.38268346, 0.88192123, 0.47139674, 0.8314696,
            0.55557024, 0.77301043, 0.63439333, 0.0, 0.0, 0.98078525, 0.19509032, 0.9238795,
            0.38268346, 0.8314696, 0.55557024, 0.70710677, 0.70710677, 0.5555702, 0.83146966,
            0.38268343, 0.9238795, 0.19509023, 0.9807853, 0.0, 0.0, 0.95694035, 0.29028466,
            0.8314696, 0.55557024, 0.6343933, 0.77301043, 0.38268343, 0.9238795, 0.098017134,
            0.9951847, -0.19509032, 0.98078525, -0.47139683, 0.88192123, 0.0, 0.0, 0.9238795,
            0.38268346, 0.0, 0.0, 0.70710677, 0.70710677, 0.0, 0.0, 0.38268343, 0.9238795,
        ];
        const SPLITCACHE: [c_int; 7] = [1024, 5, 4, 4, 4, 4, 4];
        const EPSILON: c_float = 1e-6;

        let mut trigcache = [0. as c_float; SIZE * 3];
        let mut splitcache = [0 as c_int; SIZE * 32];

        unsafe {
            super::fdrffti(
                SIZE as c_int,
                trigcache.as_mut_ptr(),
                splitcache.as_mut_ptr(),
            );
        }

        assert!(trigcache[..1024].iter().all(|&x| x == 0.));
        assert!(trigcache[1024..]
            .iter()
            .zip(TRIGCACHE.iter())
            .all(|(&a, &b)| (a - b).abs() < EPSILON));
        assert!(trigcache[1024 + 1018..].iter().all(|&x| x == 0.));

        assert!(splitcache[..7]
            .iter()
            .zip(SPLITCACHE.iter())
            .all(|(&a, &b)| a == b));
        assert!(splitcache[7..].iter().all(|&x| x == 0));
    }
}
