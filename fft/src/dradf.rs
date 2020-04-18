pub(crate) fn dradf2(
    ido: i32,
    l1: i32,
    cc: &mut [f32],
    ch: &mut [f32],
    wa1: &mut [f32],
) {
    let mut t1 = 0;
    let mut t2 = l1 * ido;
    let t0 = t2;
    let mut t3 = ido << 1 as i32;

    for _ in 0..l1 {
        ch[(t1 as usize) << 1] = cc[t1 as usize] + cc[t2 as usize];
        ch[(t1 as usize) << 1 + t3 as usize - 1] =
            cc[t1 as usize] - cc[t2 as usize];
        t1 += ido;
        t2 += ido;
    }

    if ido < 2 {
        return;
    }

    if ido != 2 {
        t1 = 0;
        t2 = t0;
        for _ in 0..l1 {
            let mut t3 = t2;
            let mut t4 = (t1 << 1) + (ido << 1);
            let mut t5 = t1;
            let mut t6 = t1 + t1;
            for i in (2..ido as usize).step_by(2) {
                t3 += 2;
                t4 -= 2;
                t5 += 2;
                t6 += 2;
                let tr2 = wa1[i - 2] * cc[t3 as usize - 1]
                    + wa1[i - 1] * cc[t3 as usize];
                let ti2 = wa1[i - 2] * cc[t3 as usize]
                    - wa1[i - 1] * cc[t3 as usize - 1];
                ch[t6 as usize] = cc[t5 as usize] + ti2;
                ch[t4 as usize] = ti2 - cc[t5 as usize];
                ch[t6 as usize - 1] = cc[t5 as usize - 1] + tr2;
                ch[t4 as usize - 1] = cc[t5 as usize - 1] - tr2;
            }
            t1 += ido;
            t2 += ido;
        }

        if ido % 2 == 1 {
            return;
        }
    }

    t1 = ido;
    t2 = t1 - 1;
    t3 = t2;
    t2 += t0;

    for _ in 0..l1 {
        ch[t1 as usize] = -cc[t2 as usize];
        ch[t1 as usize - 1] = cc[t3 as usize];
        t1 += ido << 1;
        t2 += ido;
        t3 += ido;
    }
}

pub(crate) fn dradf4(
    ido: i32,
    l1: i32,
    cc: &mut [f32],
    ch: &mut [f32],
    wa1: &[f32],
    wa2: &[f32],
    wa3: &[f32],
) {
    const HSQT2: f32 = 0.707_106_781_186_547_52;

    let t0 = l1 * ido;
    let mut t1 = t0;
    let mut t4 = t1 << 1;
    let mut t2 = t1 + (t1 << 1);
    let mut t3 = 0;

    for _ in 0..l1 {
        let tr1 = cc[t1 as usize] + cc[t2 as usize];
        let tr2 = cc[t3 as usize] + cc[t4 as usize];
        let mut t5 = t3 << 2;
        ch[t5 as usize] = tr1 + tr2;
        ch[(ido << 2) as usize + t5 as usize - 1] = tr2 - tr1;
        t5 += ido << 1;
        ch[t5 as usize - 1] = cc[t3 as usize] - cc[t4 as usize];
        ch[t5 as usize] = cc[t2 as usize] - cc[t1 as usize];
        t1 += ido;
        t2 += ido;
        t3 += ido;
        t4 += ido;
    }

    if ido < 2 {
        return;
    }

    if ido != 2 {
        t1 = 0;
        for _ in 0..l1 {
            t2 = t1;
            t4 = t1 << 2;
            let t6 = ido << 1;
            let mut t5 = t6 + t4;
            for i in (2..ido as usize).step_by(2) {
                t2 += 2;
                t3 = t2;
                t4 += 2;
                t5 -= 2;
                t3 += t0;
                let cr2 = wa1[i - 2] * cc[t3 as usize - 1]
                    + wa1[i - 1] * cc[t3 as usize];
                let ci2 = wa1[i - 2] * cc[t3 as usize]
                    - wa1[i - 1] * cc[t3 as usize - 1];
                t3 += t0;
                let cr3 = wa2[i - 2] * cc[t3 as usize - 1]
                    + wa2[i - 1] * cc[t3 as usize];
                let ci3 = wa2[i - 2] * cc[t3 as usize]
                    - wa2[i - 1] * cc[t3 as usize - 1];
                t3 += t0;
                let cr4 = wa3[i - 2] * cc[t3 as usize - 1]
                    + wa3[i - 1] * cc[t3 as usize];
                let ci4 = wa3[i - 2] * cc[t3 as usize]
                    - wa3[i - 1] * cc[t3 as usize - 1];
                let tr1 = cr2 + cr4;
                let tr4 = cr4 - cr2;
                let ti1 = ci2 + ci4;
                let ti4 = ci2 - ci4;
                let ti2 = cc[t2 as usize] + ci3;
                let ti3 = cc[t2 as usize] - ci3;
                let tr2 = cc[t2 as usize - 1] + cr3;
                let tr3 = cc[t2 as usize - 1] - cr3;
                ch[t4 as usize - 1] = tr1 + tr2;
                ch[t4 as usize] = ti1 + ti2;
                ch[t5 as usize - 1] = tr3 - ti4;
                ch[t5 as usize] = tr4 - ti3;
                ch[t4 as usize + t6 as usize - 1] = ti4 + tr3;
                ch[t4 as usize + t6 as usize] = tr4 + ti3;
                ch[t5 as usize + t6 as usize - 1] = tr2 - tr1;
                ch[t5 as usize + t6 as usize] = ti1 - ti2;
            }
            t1 += ido;
        }

        if ido & 1 != 0 {
            return;
        }
    }

    t1 = t0 + ido - 1;
    t2 = t1 + (t0 << 1);
    t3 = ido << 2;
    t4 = ido;
    let t5 = ido << 1;
    let mut t6 = ido;

    for _ in 0..l1 {
        let ti1 = -HSQT2 * (cc[t1 as usize] + cc[t2 as usize]);
        let tr1 = HSQT2 * (cc[t1 as usize] - cc[t2 as usize]);
        ch[t4 as usize - 1] = tr1 + cc[t6 as usize - 1];
        ch[t4 as usize + t5 as usize - 1] = cc[t6 as usize - 1] - tr1;
        ch[t4 as usize] = ti1 - cc[t1 as usize + t0 as usize];
        ch[t4 as usize + t5 as usize] = ti1 + cc[t1 as usize + t0 as usize];
        t1 += ido;
        t2 += ido;
        t4 += t3;
        t6 += ido;
    }
}

#[inline(always)]
fn dradfg_l102(
    ido: i32,
    nbd: i32,
    ip: i32,
    l1: i32,
    t0: i32,
    cc: &[f32],
    wa: &[f32],
    ch: &mut [f32],
) {
    let mut is = -ido;
    let mut t1 = 0 as i32;
    if nbd > l1 {
        for _ in 1..ip {
            t1 += t0;
            is += ido;
            let mut t2 = -ido + t1;
            for _ in 0..l1 {
                let mut idij = is - 1 as i32;
                t2 += ido;
                let mut t3 = t2;
                for _ in (2..ido as usize).step_by(2) {
                    idij += 2 as i32;
                    t3 += 2 as i32;
                    ch[t3 as usize - 1] = wa[idij as usize - 1]
                        * cc[t3 as usize - 1]
                        + wa[idij as usize] * cc[t3 as usize];
                    ch[t3 as usize] = wa[idij as usize - 1] * cc[t3 as usize]
                        - wa[idij as usize] * cc[t3 as usize - 1];
                }
            }
        }
    } else {
        for _ in 0..ip {
            is += ido;
            let mut idij = is - 1 as i32;
            t1 += t0;
            let mut t2 = t1;
            for _ in (2..ido as usize).step_by(2) {
                idij += 2 as i32;
                t2 += 2 as i32;
                let mut t3 = t2;
                for _ in 0..l1 {
                    ch[t3 as usize - 1] = wa[idij as usize - 1]
                        * cc[t3 as usize - 1]
                        + wa[idij as usize] * cc[t3 as usize];
                    ch[t3 as usize] = wa[idij as usize - 1] * cc[t3 as usize]
                        - wa[idij as usize] * cc[t3 as usize - 1];
                    t3 += ido;
                }
            }
        }
    }
}

#[inline(always)]
fn dradfg_l103(
    ido: i32,
    nbd: i32,
    ipp2: i32,
    ipph: i32,
    l1: i32,
    t0: i32,
    ch: &[f32],
    cc: &mut [f32],
) {
    let mut t1 = 0 as i32;
    let mut t2 = ipp2 * t0;
    if nbd < l1 {
        for _ in 1..ipph {
            t1 += t0;
            t2 -= t0;
            let mut t3 = t1;
            let mut t4 = t2;
            for _ in (2..ido as usize).step_by(2) {
                t3 += 2 as i32;
                t4 += 2 as i32;
                let mut t5 = t3 - ido;
                let mut t6 = t4 - ido;
                for _ in 0..l1 {
                    t5 += ido;
                    t6 += ido;
                    cc[t5 as usize - 1] =
                        ch[t5 as usize - 1] + ch[t6 as usize - 1];
                    cc[t6 as usize - 1] = ch[t5 as usize] - ch[t6 as usize];
                    cc[t5 as usize] = ch[t5 as usize] + ch[t6 as usize];
                    cc[t6 as usize] =
                        ch[t6 as usize - 1] - ch[t5 as usize - 1];
                }
            }
        }
    } else {
        for _ in 1..ipph {
            t1 += t0;
            t2 -= t0;
            let mut t3 = t1;
            let mut t4 = t2;
            for _ in 0..l1 {
                let mut t5 = t3;
                let mut t6 = t4;
                for _ in (2..ido as usize).step_by(2) {
                    t5 += 2 as i32;
                    t6 += 2 as i32;
                    cc[t5 as usize - 1] =
                        ch[t5 as usize - 1] + ch[t6 as usize - 1];
                    cc[t6 as usize - 1] = ch[t5 as usize] - ch[t6 as usize];
                    cc[t5 as usize] = ch[t5 as usize] + ch[t6 as usize];
                    cc[t6 as usize] =
                        ch[t6 as usize - 1] - ch[t5 as usize - 1];
                }
                t3 += ido;
                t4 += ido;
            }
        }
    }
}

#[inline(always)]
fn dradfg_l104(
    ido: i32,
    nbd: i32,
    idp2: i32,
    ipp2: i32,
    ipph: i32,
    l1: i32,
    t0: i32,
    t2: i32,
    t10: i32,
    ch: &[f32],
    cc: &mut [f32],
) {
    if nbd < l1 {
        let mut t1 = -ido;
        let mut t3 = 0 as i32;
        let mut t4 = 0 as i32;
        let mut t5 = ipp2 * t0;
        for _ in 1..ipph {
            t1 += t2;
            t3 += t2;
            t4 += t0;
            t5 -= t0;
            for i in (2..ido as i32).step_by(2) {
                let mut t6 = idp2 + t1 - i;
                let mut t7 = i + t3;
                let mut t8 = i + t4;
                let mut t9 = i + t5;
                for _ in 0..l1 {
                    cc[t7 as usize - 1] =
                        ch[t8 as usize - 1] + ch[t9 as usize - 1];
                    cc[t6 as usize - 1] =
                        ch[t8 as usize - 1] - ch[t9 as usize - 1];
                    cc[t7 as usize] = ch[t8 as usize] + ch[t9 as usize];
                    cc[t6 as usize] = ch[t9 as usize] - ch[t8 as usize];
                    t6 += t10;
                    t7 += t10;
                    t8 += ido;
                    t9 += ido;
                }
            }
        }
    } else {
        let mut t1 = -ido;
        let mut t3 = 0 as i32;
        let mut t4 = 0 as i32;
        let mut t5 = ipp2 * t0;
        for _ in 1..ipph {
            t1 += t2;
            t3 += t2;
            t4 += t0;
            t5 -= t0;
            let mut t6 = t1;
            let mut t7 = t3;
            let mut t8 = t4;
            let mut t9 = t5;
            for _ in 0..l1 {
                for i in (2..ido as usize).step_by(2) {
                    let ic = idp2 - i as i32;
                    cc[i + t7 as usize - 1] =
                        ch[i + t8 as usize - 1] + ch[i + t9 as usize - 1];
                    cc[ic as usize + t6 as usize - 1] =
                        ch[i + t8 as usize - 1] - ch[i + t9 as usize - 1];
                    cc[i + t7 as usize] =
                        ch[i + t8 as usize] + ch[i + t9 as usize];
                    cc[ic as usize + t6 as usize] =
                        ch[i + t9 as usize] - ch[i + t8 as usize];
                }
                t6 += t10;
                t7 += t10;
                t8 += ido;
                t9 += ido;
            }
        }
    }
}

pub(crate) fn dradfg(
    ido: i32,
    ip: i32,
    l1: i32,
    idl1: i32,
    cc: &mut [f32],
    ch: &mut [f32],
    wa: &[f32],
) {
    const TPI: f32 = 6.283_185_307_179_586;

    let arg = TPI / ip as f32;
    let dcp = f64::cos(arg as f64) as f32;
    let dsp = f64::sin(arg as f64) as f32;
    let ipph = ip + 1 >> 1;
    let ipp2 = ip;
    let idp2 = ido;
    let nbd = ido - 1 >> 1;
    let t0 = l1 * ido;
    let t10 = ip * ido;

    if ido != 1 {
        for ik in 0..idl1 {
            ch[ik as usize] = cc[ik as usize];
        }

        let mut t1 = 0;
        for _ in 1..ip {
            t1 += t0;
            let mut t2 = t1;
            for _ in 0..l1 {
                ch[t2 as usize] = cc[t2 as usize];
                t2 += ido;
            }
        }

        dradfg_l102(ido, nbd, ip, l1, t0, cc, wa, ch);

        dradfg_l103(ido, nbd, ipp2, ipph, l1, t0, ch, cc);
    }

    for ik in 0..idl1 {
        cc[ik as usize] = ch[ik as usize];
    }

    let mut t1 = 0 as i32;
    let mut t2 = ipp2 * idl1;
    for _ in 1..ipph {
        t1 += t0;
        t2 -= t0;
        let mut t3 = t1 - ido;
        let mut t4 = t2 - ido;
        for _ in 0..l1 {
            t3 += ido;
            t4 += ido;
            cc[t3 as usize] = ch[t3 as usize] + ch[t4 as usize];
            cc[t4 as usize] = ch[t4 as usize] - ch[t3 as usize];
        }
    }

    let mut ar1: f32 = 1.0;
    let mut ai1: f32 = 0.0;
    t1 = 0 as i32;
    t2 = ipp2 * idl1;
    let mut t3 = (ip - 1 as i32) * idl1;

    for _ in 1..ipph {
        t1 += idl1;
        t2 -= idl1;
        let ar1h = dcp * ar1 - dsp * ai1;
        ai1 = dcp * ai1 + dsp * ar1;
        ar1 = ar1h;
        let mut t4 = t1;
        let mut t5 = t2;
        let mut t6 = t3;
        let mut t7 = idl1;
        for ik in 0..idl1 {
            let fresh2 = t7;
            t7 = t7 + 1;
            let fresh3 = t4;
            t4 = t4 + 1;
            ch[fresh3 as usize] = cc[ik as usize] + ar1 * cc[fresh2 as usize];
            let fresh4 = t6;
            t6 = t6 + 1;
            let fresh5 = t5;
            t5 = t5 + 1;
            ch[fresh5 as usize] = ai1 * cc[fresh4 as usize];
        }
        let dc2 = ar1;
        let ds2 = ai1;
        let mut ar2 = ar1;
        let mut ai2 = ai1;
        t4 = idl1;
        t5 = (ipp2 - 1 as i32) * idl1;
        for _ in 2..ipph {
            t4 += idl1;
            t5 -= idl1;
            let ar2h = dc2 * ar2 - ds2 * ai2;
            ai2 = dc2 * ai2 + ds2 * ar2;
            ar2 = ar2h;
            let mut t6 = t1;
            let mut t7 = t2;
            let mut t8 = t4;
            let mut t9 = t5;
            for _ in 0..idl1 {
                let fresh6 = t8;
                t8 = t8 + 1;
                let fresh7 = t6;
                t6 = t6 + 1;
                ch[fresh7 as usize] += ar2 * cc[fresh6 as usize];
                let fresh8 = t9;
                t9 = t9 + 1;
                let fresh9 = t7;
                t7 = t7 + 1;
                ch[fresh9 as usize] += ai2 * cc[fresh8 as usize];
            }
        }
    }

    t1 = 0 as i32;
    for _ in 1..ipph {
        t1 += idl1;
        t2 = t1;
        for ik in 0..idl1 {
            let fresh10 = t2;
            t2 = t2 + 1;
            ch[ik as usize] += cc[fresh10 as usize];
        }
    }

    if ido < l1 {
        for i in 0..ido {
            t1 = i;
            t2 = i;
            for _ in 0..l1 {
                cc[t2 as usize] = ch[t1 as usize];
                t1 += ido;
                t2 += t10;
            }
        }
    } else {
        t1 = 0 as i32;
        t2 = 0 as i32;
        for _ in 0..l1 {
            t3 = t1;
            let mut t4 = t2;
            for _ in 0..ido {
                let fresh11 = t3;
                t3 = t3 + 1;
                let fresh12 = t4;
                t4 = t4 + 1;
                cc[fresh12 as usize] = ch[fresh11 as usize];
            }
            t1 += ido;
            t2 += t10;
        }
    }

    t1 = 0 as i32;
    t2 = ido << 1 as i32;
    t3 = 0 as i32;
    let mut t4 = ipp2 * t0;
    for _ in 1..ipph {
        t1 += t2;
        t3 += t0;
        t4 -= t0;
        let mut t5 = t1;
        let mut t6 = t3;
        let mut t7 = t4;
        for _ in 0..l1 {
            cc[t5 as usize - 1] = ch[t6 as usize];
            cc[t5 as usize] = ch[t7 as usize];
            t5 += t10;
            t6 += ido;
            t7 += ido;
        }
    }

    if ido == 1 {
        return;
    }

    dradfg_l104(ido, nbd, idp2, ipp2, ipph, l1, t0, t2, t10, ch, cc);
}
