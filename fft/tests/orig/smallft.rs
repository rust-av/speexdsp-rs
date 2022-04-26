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

    fn calloc(_: c_ulong, _: c_ulong) -> *mut c_void;

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
pub unsafe fn speex_alloc(mut size: c_int) -> *mut c_void {
    calloc(size as c_ulong, 1_i32 as c_ulong)
}
#[inline]
pub unsafe fn speex_free(mut ptr: *mut c_void) {
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
pub unsafe fn drfti1(
    mut n: c_int,
    mut wa: *mut c_float,
    mut ifac: *mut c_int,
) {
    static mut ntryh: [c_int; 4] = [4_i32, 2_i32, 3_i32, 5_i32];
    static mut tpi: c_float = 6.283_185_5_f32;
    let mut arg: c_float = 0.;
    let mut argh: c_float = 0.;
    let mut argld: c_float = 0.;
    let mut fi: c_float = 0.;
    let mut ntry: c_int = 0_i32;
    let mut i: c_int = 0;
    let mut j: c_int = -1_i32;
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
    let mut nf: c_int = 0_i32;
    'c_10244: loop {
        j += 1;
        if j < 4_i32 {
            ntry = ntryh[j as usize]
        } else {
            ntry += 2_i32
        }
        loop {
            nq = nl / ntry;
            nr = nl - ntry * nq;
            if nr != 0_i32 {
                break;
            }
            nf += 1;
            *ifac.offset((nf + 1_i32) as isize) = ntry;
            nl = nq;
            if ntry == 2_i32 && nf != 1_i32 {
                i = 1_i32;
                while i < nf {
                    ib = nf - i + 1_i32;
                    *ifac.offset((ib + 1_i32) as isize) =
                        *ifac.offset(ib as isize);
                    i += 1
                }
                *ifac.offset(2_i32 as isize) = 2_i32
            }
            if nl == 1_i32 {
                break 'c_10244;
            }
        }
    }
    *ifac.offset(0_i32 as isize) = n;
    *ifac.offset(1_i32 as isize) = nf;
    argh = tpi / n as c_float;
    is = 0_i32;
    nfm1 = nf - 1_i32;
    l1 = 1_i32;
    if nfm1 == 0_i32 {
        return;
    }
    k1 = 0_i32;
    while k1 < nfm1 {
        ip = *ifac.offset((k1 + 2_i32) as isize);
        ld = 0_i32;
        l2 = l1 * ip;
        ido = n / l2;
        ipm = ip - 1_i32;
        j = 0_i32;
        while j < ipm {
            ld += l1;
            i = is;
            argld = ld as c_float * argh;
            fi = 0.0f32;
            ii = 2_i32;
            while ii < ido {
                fi += 1.0f32;
                arg = fi * argld;
                let fresh0 = i;
                i += 1;
                *wa.offset(fresh0 as isize) =
                    f64::cos(arg as c_double) as c_float;
                let fresh1 = i;
                i += 1;
                *wa.offset(fresh1 as isize) =
                    f64::sin(arg as c_double) as c_float;
                ii += 2_i32
            }
            is += ido;
            j += 1
        }
        l1 = l2;
        k1 += 1
    }
}
pub unsafe fn fdrffti(
    mut n: c_int,
    mut wsave: *mut c_float,
    mut ifac: *mut c_int,
) {
    if n == 1_i32 {
        return;
    }
    drfti1(n, wsave.offset(n as isize), ifac);
}
pub unsafe fn dradf2(
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
    t1 = 0_i32;
    t2 = l1 * ido;
    t0 = t2;
    t3 = ido << 1_i32;
    k = 0_i32;
    while k < l1 {
        *ch.offset((t1 << 1_i32) as isize) =
            *cc.offset(t1 as isize) + *cc.offset(t2 as isize);
        *ch.offset(((t1 << 1_i32) + t3 - 1_i32) as isize) =
            *cc.offset(t1 as isize) - *cc.offset(t2 as isize);
        t1 += ido;
        t2 += ido;
        k += 1
    }
    if ido < 2_i32 {
        return;
    }
    if ido != 2_i32 {
        t1 = 0_i32;
        t2 = t0;
        k = 0_i32;
        while k < l1 {
            t3 = t2;
            t4 = (t1 << 1_i32) + (ido << 1_i32);
            t5 = t1;
            t6 = t1 + t1;
            i = 2_i32;
            while i < ido {
                t3 += 2_i32;
                t4 -= 2_i32;
                t5 += 2_i32;
                t6 += 2_i32;
                tr2 = *wa1.offset((i - 2_i32) as isize)
                    * *cc.offset((t3 - 1_i32) as isize)
                    + *wa1.offset((i - 1_i32) as isize)
                        * *cc.offset(t3 as isize);
                ti2 = *wa1.offset((i - 2_i32) as isize)
                    * *cc.offset(t3 as isize)
                    - *wa1.offset((i - 1_i32) as isize)
                        * *cc.offset((t3 - 1_i32) as isize);
                *ch.offset(t6 as isize) = *cc.offset(t5 as isize) + ti2;
                *ch.offset(t4 as isize) = ti2 - *cc.offset(t5 as isize);
                *ch.offset((t6 - 1_i32) as isize) =
                    *cc.offset((t5 - 1_i32) as isize) + tr2;
                *ch.offset((t4 - 1_i32) as isize) =
                    *cc.offset((t5 - 1_i32) as isize) - tr2;
                i += 2_i32
            }
            t1 += ido;
            t2 += ido;
            k += 1
        }
        if ido % 2_i32 == 1_i32 {
            return;
        }
    }
    t1 = ido;
    t2 = t1 - 1_i32;
    t3 = t2;
    t2 += t0;
    k = 0_i32;
    while k < l1 {
        *ch.offset(t1 as isize) = -*cc.offset(t2 as isize);
        *ch.offset((t1 - 1_i32) as isize) = *cc.offset(t3 as isize);
        t1 += ido << 1_i32;
        t2 += ido;
        t3 += ido;
        k += 1
    }
}
pub unsafe fn dradf4(
    mut ido: c_int,
    mut l1: c_int,
    mut cc: *mut c_float,
    mut ch: *mut c_float,
    mut wa1: *mut c_float,
    mut wa2: *mut c_float,
    mut wa3: *mut c_float,
) {
    static mut hsqt2: c_float = 0.707_106_77_f32;
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
    t4 = t1 << 1_i32;
    t2 = t1 + (t1 << 1_i32);
    t3 = 0_i32;
    k = 0_i32;
    while k < l1 {
        tr1 = *cc.offset(t1 as isize) + *cc.offset(t2 as isize);
        tr2 = *cc.offset(t3 as isize) + *cc.offset(t4 as isize);
        t5 = t3 << 2_i32;
        *ch.offset(t5 as isize) = tr1 + tr2;
        *ch.offset(((ido << 2_i32) + t5 - 1_i32) as isize) = tr2 - tr1;
        t5 += ido << 1_i32;
        *ch.offset((t5 - 1_i32) as isize) =
            *cc.offset(t3 as isize) - *cc.offset(t4 as isize);
        *ch.offset(t5 as isize) =
            *cc.offset(t2 as isize) - *cc.offset(t1 as isize);
        t1 += ido;
        t2 += ido;
        t3 += ido;
        t4 += ido;
        k += 1
    }
    if ido < 2_i32 {
        return;
    }
    if ido != 2_i32 {
        t1 = 0_i32;
        k = 0_i32;
        while k < l1 {
            t2 = t1;
            t4 = t1 << 2_i32;
            t6 = ido << 1_i32;
            t5 = t6 + t4;
            i = 2_i32;
            while i < ido {
                t2 += 2_i32;
                t3 = t2;
                t4 += 2_i32;
                t5 -= 2_i32;
                t3 += t0;
                cr2 = *wa1.offset((i - 2_i32) as isize)
                    * *cc.offset((t3 - 1_i32) as isize)
                    + *wa1.offset((i - 1_i32) as isize)
                        * *cc.offset(t3 as isize);
                ci2 = *wa1.offset((i - 2_i32) as isize)
                    * *cc.offset(t3 as isize)
                    - *wa1.offset((i - 1_i32) as isize)
                        * *cc.offset((t3 - 1_i32) as isize);
                t3 += t0;
                cr3 = *wa2.offset((i - 2_i32) as isize)
                    * *cc.offset((t3 - 1_i32) as isize)
                    + *wa2.offset((i - 1_i32) as isize)
                        * *cc.offset(t3 as isize);
                ci3 = *wa2.offset((i - 2_i32) as isize)
                    * *cc.offset(t3 as isize)
                    - *wa2.offset((i - 1_i32) as isize)
                        * *cc.offset((t3 - 1_i32) as isize);
                t3 += t0;
                cr4 = *wa3.offset((i - 2_i32) as isize)
                    * *cc.offset((t3 - 1_i32) as isize)
                    + *wa3.offset((i - 1_i32) as isize)
                        * *cc.offset(t3 as isize);
                ci4 = *wa3.offset((i - 2_i32) as isize)
                    * *cc.offset(t3 as isize)
                    - *wa3.offset((i - 1_i32) as isize)
                        * *cc.offset((t3 - 1_i32) as isize);
                tr1 = cr2 + cr4;
                tr4 = cr4 - cr2;
                ti1 = ci2 + ci4;
                ti4 = ci2 - ci4;
                ti2 = *cc.offset(t2 as isize) + ci3;
                ti3 = *cc.offset(t2 as isize) - ci3;
                tr2 = *cc.offset((t2 - 1_i32) as isize) + cr3;
                tr3 = *cc.offset((t2 - 1_i32) as isize) - cr3;
                *ch.offset((t4 - 1_i32) as isize) = tr1 + tr2;
                *ch.offset(t4 as isize) = ti1 + ti2;
                *ch.offset((t5 - 1_i32) as isize) = tr3 - ti4;
                *ch.offset(t5 as isize) = tr4 - ti3;
                *ch.offset((t4 + t6 - 1_i32) as isize) = ti4 + tr3;
                *ch.offset((t4 + t6) as isize) = tr4 + ti3;
                *ch.offset((t5 + t6 - 1_i32) as isize) = tr2 - tr1;
                *ch.offset((t5 + t6) as isize) = ti1 - ti2;
                i += 2_i32
            }
            t1 += ido;
            k += 1
        }
        if ido & 1_i32 != 0 {
            return;
        }
    }
    t1 = t0 + ido - 1_i32;
    t2 = t1 + (t0 << 1_i32);
    t3 = ido << 2_i32;
    t4 = ido;
    t5 = ido << 1_i32;
    t6 = ido;
    k = 0_i32;
    while k < l1 {
        ti1 = -hsqt2 * (*cc.offset(t1 as isize) + *cc.offset(t2 as isize));
        tr1 = hsqt2 * (*cc.offset(t1 as isize) - *cc.offset(t2 as isize));
        *ch.offset((t4 - 1_i32) as isize) =
            tr1 + *cc.offset((t6 - 1_i32) as isize);
        *ch.offset((t4 + t5 - 1_i32) as isize) =
            *cc.offset((t6 - 1_i32) as isize) - tr1;
        *ch.offset(t4 as isize) = ti1 - *cc.offset((t1 + t0) as isize);
        *ch.offset((t4 + t5) as isize) = ti1 + *cc.offset((t1 + t0) as isize);
        t1 += ido;
        t2 += ido;
        t4 += t3;
        t6 += ido;
        k += 1
    }
}
pub unsafe fn dradfg(
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
    static mut tpi: c_float = 6.283_185_5_f32;
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
    dcp = f64::cos(arg as c_double) as c_float;
    dsp = f64::sin(arg as c_double) as c_float;
    ipph = (ip + 1_i32) >> 1_i32;
    ipp2 = ip;
    idp2 = ido;
    nbd = (ido - 1_i32) >> 1_i32;
    t0 = l1 * ido;
    t10 = ip * ido;
    if ido != 1_i32 {
        ik = 0_i32;
        while ik < idl1 {
            *ch2.offset(ik as isize) = *c2.offset(ik as isize);
            ik += 1
        }
        t1 = 0_i32;
        j = 1_i32;
        while j < ip {
            t1 += t0;
            t2 = t1;
            k = 0_i32;
            while k < l1 {
                *ch.offset(t2 as isize) = *c1.offset(t2 as isize);
                t2 += ido;
                k += 1
            }
            j += 1
        }
        is = -ido;
        t1 = 0_i32;
        if nbd > l1 {
            j = 1_i32;
            while j < ip {
                t1 += t0;
                is += ido;
                t2 = -ido + t1;
                k = 0_i32;
                while k < l1 {
                    idij = is - 1_i32;
                    t2 += ido;
                    t3 = t2;
                    i = 2_i32;
                    while i < ido {
                        idij += 2_i32;
                        t3 += 2_i32;
                        *ch.offset((t3 - 1_i32) as isize) = *wa
                            .offset((idij - 1_i32) as isize)
                            * *c1.offset((t3 - 1_i32) as isize)
                            + *wa.offset(idij as isize)
                                * *c1.offset(t3 as isize);
                        *ch.offset(t3 as isize) = *wa
                            .offset((idij - 1_i32) as isize)
                            * *c1.offset(t3 as isize)
                            - *wa.offset(idij as isize)
                                * *c1.offset((t3 - 1_i32) as isize);
                        i += 2_i32
                    }
                    k += 1
                }
                j += 1
            }
        } else {
            j = 1_i32;
            while j < ip {
                is += ido;
                idij = is - 1_i32;
                t1 += t0;
                t2 = t1;
                i = 2_i32;
                while i < ido {
                    idij += 2_i32;
                    t2 += 2_i32;
                    t3 = t2;
                    k = 0_i32;
                    while k < l1 {
                        *ch.offset((t3 - 1_i32) as isize) = *wa
                            .offset((idij - 1_i32) as isize)
                            * *c1.offset((t3 - 1_i32) as isize)
                            + *wa.offset(idij as isize)
                                * *c1.offset(t3 as isize);
                        *ch.offset(t3 as isize) = *wa
                            .offset((idij - 1_i32) as isize)
                            * *c1.offset(t3 as isize)
                            - *wa.offset(idij as isize)
                                * *c1.offset((t3 - 1_i32) as isize);
                        t3 += ido;
                        k += 1
                    }
                    i += 2_i32
                }
                j += 1
            }
        }
        t1 = 0_i32;
        t2 = ipp2 * t0;
        if nbd < l1 {
            j = 1_i32;
            while j < ipph {
                t1 += t0;
                t2 -= t0;
                t3 = t1;
                t4 = t2;
                i = 2_i32;
                while i < ido {
                    t3 += 2_i32;
                    t4 += 2_i32;
                    t5 = t3 - ido;
                    t6 = t4 - ido;
                    k = 0_i32;
                    while k < l1 {
                        t5 += ido;
                        t6 += ido;
                        *c1.offset((t5 - 1_i32) as isize) = *ch
                            .offset((t5 - 1_i32) as isize)
                            + *ch.offset((t6 - 1_i32) as isize);
                        *c1.offset((t6 - 1_i32) as isize) =
                            *ch.offset(t5 as isize) - *ch.offset(t6 as isize);
                        *c1.offset(t5 as isize) =
                            *ch.offset(t5 as isize) + *ch.offset(t6 as isize);
                        *c1.offset(t6 as isize) = *ch
                            .offset((t6 - 1_i32) as isize)
                            - *ch.offset((t5 - 1_i32) as isize);
                        k += 1
                    }
                    i += 2_i32
                }
                j += 1
            }
        } else {
            j = 1_i32;
            while j < ipph {
                t1 += t0;
                t2 -= t0;
                t3 = t1;
                t4 = t2;
                k = 0_i32;
                while k < l1 {
                    t5 = t3;
                    t6 = t4;
                    i = 2_i32;
                    while i < ido {
                        t5 += 2_i32;
                        t6 += 2_i32;
                        *c1.offset((t5 - 1_i32) as isize) = *ch
                            .offset((t5 - 1_i32) as isize)
                            + *ch.offset((t6 - 1_i32) as isize);
                        *c1.offset((t6 - 1_i32) as isize) =
                            *ch.offset(t5 as isize) - *ch.offset(t6 as isize);
                        *c1.offset(t5 as isize) =
                            *ch.offset(t5 as isize) + *ch.offset(t6 as isize);
                        *c1.offset(t6 as isize) = *ch
                            .offset((t6 - 1_i32) as isize)
                            - *ch.offset((t5 - 1_i32) as isize);
                        i += 2_i32
                    }
                    t3 += ido;
                    t4 += ido;
                    k += 1
                }
                j += 1
            }
        }
    }
    ik = 0_i32;
    while ik < idl1 {
        *c2.offset(ik as isize) = *ch2.offset(ik as isize);
        ik += 1
    }
    t1 = 0_i32;
    t2 = ipp2 * idl1;
    j = 1_i32;
    while j < ipph {
        t1 += t0;
        t2 -= t0;
        t3 = t1 - ido;
        t4 = t2 - ido;
        k = 0_i32;
        while k < l1 {
            t3 += ido;
            t4 += ido;
            *c1.offset(t3 as isize) =
                *ch.offset(t3 as isize) + *ch.offset(t4 as isize);
            *c1.offset(t4 as isize) =
                *ch.offset(t4 as isize) - *ch.offset(t3 as isize);
            k += 1
        }
        j += 1
    }
    ar1 = 1.0f32;
    ai1 = 0.0f32;
    t1 = 0_i32;
    t2 = ipp2 * idl1;
    t3 = (ip - 1_i32) * idl1;
    l = 1_i32;
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
        ik = 0_i32;
        while ik < idl1 {
            let fresh2 = t7;
            t7 += 1;
            let fresh3 = t4;
            t4 += 1;
            *ch2.offset(fresh3 as isize) =
                *c2.offset(ik as isize) + ar1 * *c2.offset(fresh2 as isize);
            let fresh4 = t6;
            t6 += 1;
            let fresh5 = t5;
            t5 += 1;
            *ch2.offset(fresh5 as isize) = ai1 * *c2.offset(fresh4 as isize);
            ik += 1
        }
        dc2 = ar1;
        ds2 = ai1;
        ar2 = ar1;
        ai2 = ai1;
        t4 = idl1;
        t5 = (ipp2 - 1_i32) * idl1;
        j = 2_i32;
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
            ik = 0_i32;
            while ik < idl1 {
                let fresh6 = t8;
                t8 += 1;
                let fresh7 = t6;
                t6 += 1;
                *ch2.offset(fresh7 as isize) +=
                    ar2 * *c2.offset(fresh6 as isize);
                let fresh8 = t9;
                t9 += 1;
                let fresh9 = t7;
                t7 += 1;
                *ch2.offset(fresh9 as isize) +=
                    ai2 * *c2.offset(fresh8 as isize);
                ik += 1
            }
            j += 1
        }
        l += 1
    }
    t1 = 0_i32;
    j = 1_i32;
    while j < ipph {
        t1 += idl1;
        t2 = t1;
        ik = 0_i32;
        while ik < idl1 {
            let fresh10 = t2;
            t2 += 1;
            *ch2.offset(ik as isize) += *c2.offset(fresh10 as isize);
            ik += 1
        }
        j += 1
    }
    if ido < l1 {
        i = 0_i32;
        while i < ido {
            t1 = i;
            t2 = i;
            k = 0_i32;
            while k < l1 {
                *cc.offset(t2 as isize) = *ch.offset(t1 as isize);
                t1 += ido;
                t2 += t10;
                k += 1
            }
            i += 1
        }
    } else {
        t1 = 0_i32;
        t2 = 0_i32;
        k = 0_i32;
        while k < l1 {
            t3 = t1;
            t4 = t2;
            i = 0_i32;
            while i < ido {
                let fresh11 = t3;
                t3 += 1;
                let fresh12 = t4;
                t4 += 1;
                *cc.offset(fresh12 as isize) = *ch.offset(fresh11 as isize);
                i += 1
            }
            t1 += ido;
            t2 += t10;
            k += 1
        }
    }
    t1 = 0_i32;
    t2 = ido << 1_i32;
    t3 = 0_i32;
    t4 = ipp2 * t0;
    j = 1_i32;
    while j < ipph {
        t1 += t2;
        t3 += t0;
        t4 -= t0;
        t5 = t1;
        t6 = t3;
        t7 = t4;
        k = 0_i32;
        while k < l1 {
            *cc.offset((t5 - 1_i32) as isize) = *ch.offset(t6 as isize);
            *cc.offset(t5 as isize) = *ch.offset(t7 as isize);
            t5 += t10;
            t6 += ido;
            t7 += ido;
            k += 1
        }
        j += 1
    }
    if ido == 1_i32 {
        return;
    }
    if nbd < l1 {
        t1 = -ido;
        t3 = 0_i32;
        t4 = 0_i32;
        t5 = ipp2 * t0;
        j = 1_i32;
        while j < ipph {
            t1 += t2;
            t3 += t2;
            t4 += t0;
            t5 -= t0;
            i = 2_i32;
            while i < ido {
                t6 = idp2 + t1 - i;
                t7 = i + t3;
                t8 = i + t4;
                t9 = i + t5;
                k = 0_i32;
                while k < l1 {
                    *cc.offset((t7 - 1_i32) as isize) = *ch
                        .offset((t8 - 1_i32) as isize)
                        + *ch.offset((t9 - 1_i32) as isize);
                    *cc.offset((t6 - 1_i32) as isize) = *ch
                        .offset((t8 - 1_i32) as isize)
                        - *ch.offset((t9 - 1_i32) as isize);
                    *cc.offset(t7 as isize) =
                        *ch.offset(t8 as isize) + *ch.offset(t9 as isize);
                    *cc.offset(t6 as isize) =
                        *ch.offset(t9 as isize) - *ch.offset(t8 as isize);
                    t6 += t10;
                    t7 += t10;
                    t8 += ido;
                    t9 += ido;
                    k += 1
                }
                i += 2_i32
            }
            j += 1
        }
    } else {
        t1 = -ido;
        t3 = 0_i32;
        t4 = 0_i32;
        t5 = ipp2 * t0;
        j = 1_i32;
        while j < ipph {
            t1 += t2;
            t3 += t2;
            t4 += t0;
            t5 -= t0;
            t6 = t1;
            t7 = t3;
            t8 = t4;
            t9 = t5;
            k = 0_i32;
            while k < l1 {
                i = 2_i32;
                while i < ido {
                    ic = idp2 - i;
                    *cc.offset((i + t7 - 1_i32) as isize) = *ch
                        .offset((i + t8 - 1_i32) as isize)
                        + *ch.offset((i + t9 - 1_i32) as isize);
                    *cc.offset((ic + t6 - 1_i32) as isize) = *ch
                        .offset((i + t8 - 1_i32) as isize)
                        - *ch.offset((i + t9 - 1_i32) as isize);
                    *cc.offset((i + t7) as isize) = *ch
                        .offset((i + t8) as isize)
                        + *ch.offset((i + t9) as isize);
                    *cc.offset((ic + t6) as isize) = *ch
                        .offset((i + t9) as isize)
                        - *ch.offset((i + t8) as isize);
                    i += 2_i32
                }
                t6 += t10;
                t7 += t10;
                t8 += ido;
                t9 += ido;
                k += 1
            }
            j += 1
        }
    };
}
pub unsafe fn drftf1(
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
    nf = *ifac.offset(1_i32 as isize);
    na = 1_i32;
    l2 = n;
    iw = n;
    k1 = 0_i32;
    while k1 < nf {
        kh = nf - k1;
        ip = *ifac.offset((kh + 1_i32) as isize);
        l1 = l2 / ip;
        ido = n / l2;
        idl1 = ido * l1;
        iw -= (ip - 1_i32) * ido;
        na = 1_i32 - na;
        if ip != 4_i32 {
            if ip != 2_i32 {
                if ido == 1_i32 {
                    na = 1_i32 - na
                }
                if na != 0_i32 {
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
                        wa.offset(iw as isize).offset(-(1_i32 as isize)),
                    );
                    na = 0_i32
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
                        wa.offset(iw as isize).offset(-(1_i32 as isize)),
                    );
                    na = 1_i32
                }
            } else if na != 0_i32 {
                dradf2(
                    ido,
                    l1,
                    ch,
                    c,
                    wa.offset(iw as isize).offset(-(1_i32 as isize)),
                );
            } else {
                dradf2(
                    ido,
                    l1,
                    c,
                    ch,
                    wa.offset(iw as isize).offset(-(1_i32 as isize)),
                );
            }
        } else {
            ix2 = iw + ido;
            ix3 = ix2 + ido;
            if na != 0_i32 {
                dradf4(
                    ido,
                    l1,
                    ch,
                    c,
                    wa.offset(iw as isize).offset(-(1_i32 as isize)),
                    wa.offset(ix2 as isize).offset(-(1_i32 as isize)),
                    wa.offset(ix3 as isize).offset(-(1_i32 as isize)),
                );
            } else {
                dradf4(
                    ido,
                    l1,
                    c,
                    ch,
                    wa.offset(iw as isize).offset(-(1_i32 as isize)),
                    wa.offset(ix2 as isize).offset(-(1_i32 as isize)),
                    wa.offset(ix3 as isize).offset(-(1_i32 as isize)),
                );
            }
        }
        l2 = l1;
        k1 += 1
    }
    if na == 1_i32 {
        return;
    }
    i = 0_i32;
    while i < n {
        *c.offset(i as isize) = *ch.offset(i as isize);
        i += 1
    }
}
pub unsafe fn dradb2(
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
    t1 = 0_i32;
    t2 = 0_i32;
    t3 = (ido << 1_i32) - 1_i32;
    k = 0_i32;
    while k < l1 {
        *ch.offset(t1 as isize) =
            *cc.offset(t2 as isize) + *cc.offset((t3 + t2) as isize);
        *ch.offset((t1 + t0) as isize) =
            *cc.offset(t2 as isize) - *cc.offset((t3 + t2) as isize);
        t1 += ido;
        t2 = t1 << 1_i32;
        k += 1
    }
    if ido < 2_i32 {
        return;
    }
    if ido != 2_i32 {
        t1 = 0_i32;
        t2 = 0_i32;
        k = 0_i32;
        while k < l1 {
            t3 = t1;
            t4 = t2;
            t5 = t4 + (ido << 1_i32);
            t6 = t0 + t1;
            i = 2_i32;
            while i < ido {
                t3 += 2_i32;
                t4 += 2_i32;
                t5 -= 2_i32;
                t6 += 2_i32;
                *ch.offset((t3 - 1_i32) as isize) = *cc
                    .offset((t4 - 1_i32) as isize)
                    + *cc.offset((t5 - 1_i32) as isize);
                tr2 = *cc.offset((t4 - 1_i32) as isize)
                    - *cc.offset((t5 - 1_i32) as isize);
                *ch.offset(t3 as isize) =
                    *cc.offset(t4 as isize) - *cc.offset(t5 as isize);
                ti2 = *cc.offset(t4 as isize) + *cc.offset(t5 as isize);
                *ch.offset((t6 - 1_i32) as isize) =
                    *wa1.offset((i - 2_i32) as isize) * tr2
                        - *wa1.offset((i - 1_i32) as isize) * ti2;
                *ch.offset(t6 as isize) = *wa1.offset((i - 2_i32) as isize)
                    * ti2
                    + *wa1.offset((i - 1_i32) as isize) * tr2;
                i += 2_i32
            }
            t1 += ido;
            t2 = t1 << 1_i32;
            k += 1
        }
        if ido % 2_i32 == 1_i32 {
            return;
        }
    }
    t1 = ido - 1_i32;
    t2 = ido - 1_i32;
    k = 0_i32;
    while k < l1 {
        *ch.offset(t1 as isize) =
            *cc.offset(t2 as isize) + *cc.offset(t2 as isize);
        *ch.offset((t1 + t0) as isize) = -(*cc.offset((t2 + 1_i32) as isize)
            + *cc.offset((t2 + 1_i32) as isize));
        t1 += ido;
        t2 += ido << 1_i32;
        k += 1
    }
}
pub unsafe fn dradb3(
    mut ido: c_int,
    mut l1: c_int,
    mut cc: *mut c_float,
    mut ch: *mut c_float,
    mut wa1: *mut c_float,
    mut wa2: *mut c_float,
) {
    static mut taur: c_float = -0.5f32;
    static mut taui: c_float = 0.866_025_4_f32;
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
    t1 = 0_i32;
    t2 = t0 << 1_i32;
    t3 = ido << 1_i32;
    t4 = ido + (ido << 1_i32);
    t5 = 0_i32;
    k = 0_i32;
    while k < l1 {
        tr2 = *cc.offset((t3 - 1_i32) as isize)
            + *cc.offset((t3 - 1_i32) as isize);
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
    if ido == 1_i32 {
        return;
    }
    t1 = 0_i32;
    t3 = ido << 1_i32;
    k = 0_i32;
    while k < l1 {
        t7 = t1 + (t1 << 1_i32);
        t5 = t7 + t3;
        t6 = t5;
        t8 = t1;
        t9 = t1 + t0;
        t10 = t9 + t0;
        i = 2_i32;
        while i < ido {
            t5 += 2_i32;
            t6 -= 2_i32;
            t7 += 2_i32;
            t8 += 2_i32;
            t9 += 2_i32;
            t10 += 2_i32;
            tr2 = *cc.offset((t5 - 1_i32) as isize)
                + *cc.offset((t6 - 1_i32) as isize);
            cr2 = *cc.offset((t7 - 1_i32) as isize) + taur * tr2;
            *ch.offset((t8 - 1_i32) as isize) =
                *cc.offset((t7 - 1_i32) as isize) + tr2;
            ti2 = *cc.offset(t5 as isize) - *cc.offset(t6 as isize);
            ci2 = *cc.offset(t7 as isize) + taur * ti2;
            *ch.offset(t8 as isize) = *cc.offset(t7 as isize) + ti2;
            cr3 = taui
                * (*cc.offset((t5 - 1_i32) as isize)
                    - *cc.offset((t6 - 1_i32) as isize));
            ci3 = taui * (*cc.offset(t5 as isize) + *cc.offset(t6 as isize));
            dr2 = cr2 - ci3;
            dr3 = cr2 + ci3;
            di2 = ci2 + cr3;
            di3 = ci2 - cr3;
            *ch.offset((t9 - 1_i32) as isize) =
                *wa1.offset((i - 2_i32) as isize) * dr2
                    - *wa1.offset((i - 1_i32) as isize) * di2;
            *ch.offset(t9 as isize) = *wa1.offset((i - 2_i32) as isize) * di2
                + *wa1.offset((i - 1_i32) as isize) * dr2;
            *ch.offset((t10 - 1_i32) as isize) =
                *wa2.offset((i - 2_i32) as isize) * dr3
                    - *wa2.offset((i - 1_i32) as isize) * di3;
            *ch.offset(t10 as isize) = *wa2.offset((i - 2_i32) as isize) * di3
                + *wa2.offset((i - 1_i32) as isize) * dr3;
            i += 2_i32
        }
        t1 += ido;
        k += 1
    }
}
pub unsafe fn dradb4(
    mut ido: c_int,
    mut l1: c_int,
    mut cc: *mut c_float,
    mut ch: *mut c_float,
    mut wa1: *mut c_float,
    mut wa2: *mut c_float,
    mut wa3: *mut c_float,
) {
    static mut sqrt2: c_float = std::f32::consts::SQRT_2;
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
    t1 = 0_i32;
    t2 = ido << 2_i32;
    t3 = 0_i32;
    t6 = ido << 1_i32;
    k = 0_i32;
    while k < l1 {
        t4 = t3 + t6;
        t5 = t1;
        tr3 = *cc.offset((t4 - 1_i32) as isize)
            + *cc.offset((t4 - 1_i32) as isize);
        tr4 = *cc.offset(t4 as isize) + *cc.offset(t4 as isize);
        t4 += t6;
        tr1 = *cc.offset(t3 as isize) - *cc.offset((t4 - 1_i32) as isize);
        tr2 = *cc.offset(t3 as isize) + *cc.offset((t4 - 1_i32) as isize);
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
    if ido < 2_i32 {
        return;
    }
    if ido != 2_i32 {
        t1 = 0_i32;
        k = 0_i32;
        while k < l1 {
            t2 = t1 << 2_i32;
            t3 = t2 + t6;
            t4 = t3;
            t5 = t4 + t6;
            t7 = t1;
            i = 2_i32;
            while i < ido {
                t2 += 2_i32;
                t3 += 2_i32;
                t4 -= 2_i32;
                t5 -= 2_i32;
                t7 += 2_i32;
                ti1 = *cc.offset(t2 as isize) + *cc.offset(t5 as isize);
                ti2 = *cc.offset(t2 as isize) - *cc.offset(t5 as isize);
                ti3 = *cc.offset(t3 as isize) - *cc.offset(t4 as isize);
                tr4 = *cc.offset(t3 as isize) + *cc.offset(t4 as isize);
                tr1 = *cc.offset((t2 - 1_i32) as isize)
                    - *cc.offset((t5 - 1_i32) as isize);
                tr2 = *cc.offset((t2 - 1_i32) as isize)
                    + *cc.offset((t5 - 1_i32) as isize);
                ti4 = *cc.offset((t3 - 1_i32) as isize)
                    - *cc.offset((t4 - 1_i32) as isize);
                tr3 = *cc.offset((t3 - 1_i32) as isize)
                    + *cc.offset((t4 - 1_i32) as isize);
                *ch.offset((t7 - 1_i32) as isize) = tr2 + tr3;
                cr3 = tr2 - tr3;
                *ch.offset(t7 as isize) = ti2 + ti3;
                ci3 = ti2 - ti3;
                cr2 = tr1 - tr4;
                cr4 = tr1 + tr4;
                ci2 = ti1 + ti4;
                ci4 = ti1 - ti4;
                t8 = t7 + t0;
                *ch.offset((t8 - 1_i32) as isize) =
                    *wa1.offset((i - 2_i32) as isize) * cr2
                        - *wa1.offset((i - 1_i32) as isize) * ci2;
                *ch.offset(t8 as isize) = *wa1.offset((i - 2_i32) as isize)
                    * ci2
                    + *wa1.offset((i - 1_i32) as isize) * cr2;
                t8 += t0;
                *ch.offset((t8 - 1_i32) as isize) =
                    *wa2.offset((i - 2_i32) as isize) * cr3
                        - *wa2.offset((i - 1_i32) as isize) * ci3;
                *ch.offset(t8 as isize) = *wa2.offset((i - 2_i32) as isize)
                    * ci3
                    + *wa2.offset((i - 1_i32) as isize) * cr3;
                t8 += t0;
                *ch.offset((t8 - 1_i32) as isize) =
                    *wa3.offset((i - 2_i32) as isize) * cr4
                        - *wa3.offset((i - 1_i32) as isize) * ci4;
                *ch.offset(t8 as isize) = *wa3.offset((i - 2_i32) as isize)
                    * ci4
                    + *wa3.offset((i - 1_i32) as isize) * cr4;
                i += 2_i32
            }
            t1 += ido;
            k += 1
        }
        if ido % 2_i32 == 1_i32 {
            return;
        }
    }
    t1 = ido;
    t2 = ido << 2_i32;
    t3 = ido - 1_i32;
    t4 = ido + (ido << 1_i32);
    k = 0_i32;
    while k < l1 {
        t5 = t3;
        ti1 = *cc.offset(t1 as isize) + *cc.offset(t4 as isize);
        ti2 = *cc.offset(t4 as isize) - *cc.offset(t1 as isize);
        tr1 = *cc.offset((t1 - 1_i32) as isize)
            - *cc.offset((t4 - 1_i32) as isize);
        tr2 = *cc.offset((t1 - 1_i32) as isize)
            + *cc.offset((t4 - 1_i32) as isize);
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
pub unsafe fn dradbg(
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
    static mut tpi: c_float = 6.283_185_5_f32;
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
    dcp = f64::cos(arg as c_double) as c_float;
    dsp = f64::sin(arg as c_double) as c_float;
    nbd = (ido - 1_i32) >> 1_i32;
    ipp2 = ip;
    ipph = (ip + 1_i32) >> 1_i32;
    if ido < l1 {
        t1 = 0_i32;
        i = 0_i32;
        while i < ido {
            t2 = t1;
            t3 = t1;
            k = 0_i32;
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
        t1 = 0_i32;
        t2 = 0_i32;
        k = 0_i32;
        while k < l1 {
            t3 = t1;
            t4 = t2;
            i = 0_i32;
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
    t1 = 0_i32;
    t2 = ipp2 * t0;
    t5 = ido << 1_i32;
    t7 = t5;
    j = 1_i32;
    while j < ipph {
        t1 += t0;
        t2 -= t0;
        t3 = t1;
        t4 = t2;
        t6 = t5;
        k = 0_i32;
        while k < l1 {
            *ch.offset(t3 as isize) = *cc.offset((t6 - 1_i32) as isize)
                + *cc.offset((t6 - 1_i32) as isize);
            *ch.offset(t4 as isize) =
                *cc.offset(t6 as isize) + *cc.offset(t6 as isize);
            t3 += ido;
            t4 += ido;
            t6 += t10;
            k += 1
        }
        t5 += t7;
        j += 1
    }
    if ido != 1_i32 {
        if nbd < l1 {
            t1 = 0_i32;
            t2 = ipp2 * t0;
            t7 = 0_i32;
            j = 1_i32;
            while j < ipph {
                t1 += t0;
                t2 -= t0;
                t3 = t1;
                t4 = t2;
                t7 += ido << 1_i32;
                t8 = t7;
                t9 = t7;
                i = 2_i32;
                while i < ido {
                    t3 += 2_i32;
                    t4 += 2_i32;
                    t8 += 2_i32;
                    t9 -= 2_i32;
                    t5 = t3;
                    t6 = t4;
                    t11 = t8;
                    t12 = t9;
                    k = 0_i32;
                    while k < l1 {
                        *ch.offset((t5 - 1_i32) as isize) = *cc
                            .offset((t11 - 1_i32) as isize)
                            + *cc.offset((t12 - 1_i32) as isize);
                        *ch.offset((t6 - 1_i32) as isize) = *cc
                            .offset((t11 - 1_i32) as isize)
                            - *cc.offset((t12 - 1_i32) as isize);
                        *ch.offset(t5 as isize) = *cc.offset(t11 as isize)
                            - *cc.offset(t12 as isize);
                        *ch.offset(t6 as isize) = *cc.offset(t11 as isize)
                            + *cc.offset(t12 as isize);
                        t5 += ido;
                        t6 += ido;
                        t11 += t10;
                        t12 += t10;
                        k += 1
                    }
                    i += 2_i32
                }
                j += 1
            }
        } else {
            t1 = 0_i32;
            t2 = ipp2 * t0;
            t7 = 0_i32;
            j = 1_i32;
            while j < ipph {
                t1 += t0;
                t2 -= t0;
                t3 = t1;
                t4 = t2;
                t7 += ido << 1_i32;
                t8 = t7;
                k = 0_i32;
                while k < l1 {
                    t5 = t3;
                    t6 = t4;
                    t9 = t8;
                    t11 = t8;
                    i = 2_i32;
                    while i < ido {
                        t5 += 2_i32;
                        t6 += 2_i32;
                        t9 += 2_i32;
                        t11 -= 2_i32;
                        *ch.offset((t5 - 1_i32) as isize) = *cc
                            .offset((t9 - 1_i32) as isize)
                            + *cc.offset((t11 - 1_i32) as isize);
                        *ch.offset((t6 - 1_i32) as isize) = *cc
                            .offset((t9 - 1_i32) as isize)
                            - *cc.offset((t11 - 1_i32) as isize);
                        *ch.offset(t5 as isize) =
                            *cc.offset(t9 as isize) - *cc.offset(t11 as isize);
                        *ch.offset(t6 as isize) =
                            *cc.offset(t9 as isize) + *cc.offset(t11 as isize);
                        i += 2_i32
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
    t1 = 0_i32;
    t2 = ipp2 * idl1;
    t9 = t2;
    t3 = (ip - 1_i32) * idl1;
    l = 1_i32;
    while l < ipph {
        t1 += idl1;
        t2 -= idl1;
        ar1h = dcp * ar1 - dsp * ai1;
        ai1 = dcp * ai1 + dsp * ar1;
        ar1 = ar1h;
        t4 = t1;
        t5 = t2;
        t6 = 0_i32;
        t7 = idl1;
        t8 = t3;
        ik = 0_i32;
        while ik < idl1 {
            let fresh13 = t6;
            t6 += 1;
            let fresh14 = t7;
            t7 += 1;
            let fresh15 = t4;
            t4 += 1;
            *c2.offset(fresh15 as isize) = *ch2.offset(fresh13 as isize)
                + ar1 * *ch2.offset(fresh14 as isize);
            let fresh16 = t8;
            t8 += 1;
            let fresh17 = t5;
            t5 += 1;
            *c2.offset(fresh17 as isize) = ai1 * *ch2.offset(fresh16 as isize);
            ik += 1
        }
        dc2 = ar1;
        ds2 = ai1;
        ar2 = ar1;
        ai2 = ai1;
        t6 = idl1;
        t7 = t9 - idl1;
        j = 2_i32;
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
            ik = 0_i32;
            while ik < idl1 {
                let fresh18 = t11;
                t11 += 1;
                let fresh19 = t4;
                t4 += 1;
                *c2.offset(fresh19 as isize) +=
                    ar2 * *ch2.offset(fresh18 as isize);
                let fresh20 = t12;
                t12 += 1;
                let fresh21 = t5;
                t5 += 1;
                *c2.offset(fresh21 as isize) +=
                    ai2 * *ch2.offset(fresh20 as isize);
                ik += 1
            }
            j += 1
        }
        l += 1
    }
    t1 = 0_i32;
    j = 1_i32;
    while j < ipph {
        t1 += idl1;
        t2 = t1;
        ik = 0_i32;
        while ik < idl1 {
            let fresh22 = t2;
            t2 += 1;
            *ch2.offset(ik as isize) += *ch2.offset(fresh22 as isize);
            ik += 1
        }
        j += 1
    }
    t1 = 0_i32;
    t2 = ipp2 * t0;
    j = 1_i32;
    while j < ipph {
        t1 += t0;
        t2 -= t0;
        t3 = t1;
        t4 = t2;
        k = 0_i32;
        while k < l1 {
            *ch.offset(t3 as isize) =
                *c1.offset(t3 as isize) - *c1.offset(t4 as isize);
            *ch.offset(t4 as isize) =
                *c1.offset(t3 as isize) + *c1.offset(t4 as isize);
            t3 += ido;
            t4 += ido;
            k += 1
        }
        j += 1
    }
    if ido != 1_i32 {
        if nbd < l1 {
            t1 = 0_i32;
            t2 = ipp2 * t0;
            j = 1_i32;
            while j < ipph {
                t1 += t0;
                t2 -= t0;
                t3 = t1;
                t4 = t2;
                i = 2_i32;
                while i < ido {
                    t3 += 2_i32;
                    t4 += 2_i32;
                    t5 = t3;
                    t6 = t4;
                    k = 0_i32;
                    while k < l1 {
                        *ch.offset((t5 - 1_i32) as isize) = *c1
                            .offset((t5 - 1_i32) as isize)
                            - *c1.offset(t6 as isize);
                        *ch.offset((t6 - 1_i32) as isize) = *c1
                            .offset((t5 - 1_i32) as isize)
                            + *c1.offset(t6 as isize);
                        *ch.offset(t5 as isize) = *c1.offset(t5 as isize)
                            + *c1.offset((t6 - 1_i32) as isize);
                        *ch.offset(t6 as isize) = *c1.offset(t5 as isize)
                            - *c1.offset((t6 - 1_i32) as isize);
                        t5 += ido;
                        t6 += ido;
                        k += 1
                    }
                    i += 2_i32
                }
                j += 1
            }
        } else {
            t1 = 0_i32;
            t2 = ipp2 * t0;
            j = 1_i32;
            while j < ipph {
                t1 += t0;
                t2 -= t0;
                t3 = t1;
                t4 = t2;
                k = 0_i32;
                while k < l1 {
                    t5 = t3;
                    t6 = t4;
                    i = 2_i32;
                    while i < ido {
                        t5 += 2_i32;
                        t6 += 2_i32;
                        *ch.offset((t5 - 1_i32) as isize) = *c1
                            .offset((t5 - 1_i32) as isize)
                            - *c1.offset(t6 as isize);
                        *ch.offset((t6 - 1_i32) as isize) = *c1
                            .offset((t5 - 1_i32) as isize)
                            + *c1.offset(t6 as isize);
                        *ch.offset(t5 as isize) = *c1.offset(t5 as isize)
                            + *c1.offset((t6 - 1_i32) as isize);
                        *ch.offset(t6 as isize) = *c1.offset(t5 as isize)
                            - *c1.offset((t6 - 1_i32) as isize);
                        i += 2_i32
                    }
                    t3 += ido;
                    t4 += ido;
                    k += 1
                }
                j += 1
            }
        }
    }
    if ido == 1_i32 {
        return;
    }
    ik = 0_i32;
    while ik < idl1 {
        *c2.offset(ik as isize) = *ch2.offset(ik as isize);
        ik += 1
    }
    t1 = 0_i32;
    j = 1_i32;
    while j < ip {
        t1 += t0;
        t2 = t1;
        k = 0_i32;
        while k < l1 {
            *c1.offset(t2 as isize) = *ch.offset(t2 as isize);
            t2 += ido;
            k += 1
        }
        j += 1
    }
    if nbd > l1 {
        is = -ido - 1_i32;
        t1 = 0_i32;
        j = 1_i32;
        while j < ip {
            is += ido;
            t1 += t0;
            t2 = t1;
            k = 0_i32;
            while k < l1 {
                idij = is;
                t3 = t2;
                i = 2_i32;
                while i < ido {
                    idij += 2_i32;
                    t3 += 2_i32;
                    *c1.offset((t3 - 1_i32) as isize) = *wa
                        .offset((idij - 1_i32) as isize)
                        * *ch.offset((t3 - 1_i32) as isize)
                        - *wa.offset(idij as isize) * *ch.offset(t3 as isize);
                    *c1.offset(t3 as isize) = *wa
                        .offset((idij - 1_i32) as isize)
                        * *ch.offset(t3 as isize)
                        + *wa.offset(idij as isize)
                            * *ch.offset((t3 - 1_i32) as isize);
                    i += 2_i32
                }
                t2 += ido;
                k += 1
            }
            j += 1
        }
    } else {
        is = -ido - 1_i32;
        t1 = 0_i32;
        j = 1_i32;
        while j < ip {
            is += ido;
            t1 += t0;
            idij = is;
            t2 = t1;
            i = 2_i32;
            while i < ido {
                t2 += 2_i32;
                idij += 2_i32;
                t3 = t2;
                k = 0_i32;
                while k < l1 {
                    *c1.offset((t3 - 1_i32) as isize) = *wa
                        .offset((idij - 1_i32) as isize)
                        * *ch.offset((t3 - 1_i32) as isize)
                        - *wa.offset(idij as isize) * *ch.offset(t3 as isize);
                    *c1.offset(t3 as isize) = *wa
                        .offset((idij - 1_i32) as isize)
                        * *ch.offset(t3 as isize)
                        + *wa.offset(idij as isize)
                            * *ch.offset((t3 - 1_i32) as isize);
                    t3 += ido;
                    k += 1
                }
                i += 2_i32
            }
            j += 1
        }
    };
}
pub unsafe fn drftb1(
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
    nf = *ifac.offset(1_i32 as isize);
    na = 0_i32;
    l1 = 1_i32;
    iw = 1_i32;
    k1 = 0_i32;
    while k1 < nf {
        ip = *ifac.offset((k1 + 2_i32) as isize);
        l2 = ip * l1;
        ido = n / l2;
        idl1 = ido * l1;
        if ip != 4_i32 {
            if ip != 2_i32 {
                if ip != 3_i32 {
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
                    if na != 0_i32 {
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
                            wa.offset(iw as isize).offset(-(1_i32 as isize)),
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
                            wa.offset(iw as isize).offset(-(1_i32 as isize)),
                        );
                    }
                    if ido == 1_i32 {
                        na = 1_i32 - na
                    }
                } else {
                    ix2 = iw + ido;
                    if na != 0_i32 {
                        dradb3(
                            ido,
                            l1,
                            ch,
                            c,
                            wa.offset(iw as isize).offset(-(1_i32 as isize)),
                            wa.offset(ix2 as isize).offset(-(1_i32 as isize)),
                        );
                    } else {
                        dradb3(
                            ido,
                            l1,
                            c,
                            ch,
                            wa.offset(iw as isize).offset(-(1_i32 as isize)),
                            wa.offset(ix2 as isize).offset(-(1_i32 as isize)),
                        );
                    }
                    na = 1_i32 - na
                }
            } else {
                if na != 0_i32 {
                    dradb2(
                        ido,
                        l1,
                        ch,
                        c,
                        wa.offset(iw as isize).offset(-(1_i32 as isize)),
                    );
                } else {
                    dradb2(
                        ido,
                        l1,
                        c,
                        ch,
                        wa.offset(iw as isize).offset(-(1_i32 as isize)),
                    );
                }
                na = 1_i32 - na
            }
        } else {
            ix2 = iw + ido;
            ix3 = ix2 + ido;
            if na != 0_i32 {
                dradb4(
                    ido,
                    l1,
                    ch,
                    c,
                    wa.offset(iw as isize).offset(-(1_i32 as isize)),
                    wa.offset(ix2 as isize).offset(-(1_i32 as isize)),
                    wa.offset(ix3 as isize).offset(-(1_i32 as isize)),
                );
            } else {
                dradb4(
                    ido,
                    l1,
                    c,
                    ch,
                    wa.offset(iw as isize).offset(-(1_i32 as isize)),
                    wa.offset(ix2 as isize).offset(-(1_i32 as isize)),
                    wa.offset(ix3 as isize).offset(-(1_i32 as isize)),
                );
            }
            na = 1_i32 - na
        }
        l1 = l2;
        iw += (ip - 1_i32) * ido;
        k1 += 1
    }
    if na == 0_i32 {
        return;
    }
    i = 0_i32;
    while i < n {
        *c.offset(i as isize) = *ch.offset(i as isize);
        i += 1
    }
}

pub unsafe fn spx_drft_forward(
    mut l: *mut drft_lookup,
    mut data: *mut c_float,
) {
    if (*l).n == 1_i32 {
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

pub unsafe fn spx_drft_backward(
    mut l: *mut drft_lookup,
    mut data: *mut c_float,
) {
    if (*l).n == 1_i32 {
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

pub unsafe fn spx_drft_init(mut l: *mut drft_lookup, mut n: c_int) {
    (*l).n = n;
    (*l).trigcache = speex_alloc(
        ((3_i32 * n) as c_ulong)
            .wrapping_mul(::std::mem::size_of::<c_float>() as c_ulong)
            as c_int,
    ) as *mut c_float;
    (*l).splitcache = speex_alloc(
        (32_i32 as c_ulong)
            .wrapping_mul(::std::mem::size_of::<c_int>() as c_ulong)
            as c_int,
    ) as *mut c_int;
    fdrffti(n, (*l).trigcache, (*l).splitcache);
}

pub unsafe fn spx_drft_clear(mut l: *mut drft_lookup) {
    if !l.is_null() {
        if !(*l).trigcache.is_null() {
            speex_free((*l).trigcache as *mut c_void);
        }
        if !(*l).splitcache.is_null() {
            speex_free((*l).splitcache as *mut c_void);
        }
    };
}
