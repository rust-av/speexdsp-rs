pub(crate) fn dradb2(
    ido: i32,
    l1: i32,
    cc: &[f32],
    ch: &mut [f32],
    wa1: &[f32],
) {
    let t0 = l1 * ido;
    let mut t1 = 0;
    let mut t2 = 0;
    let mut t3 = (ido << 1) - 1_i32;

    for _ in 0..l1 {
        ch[t1 as usize] = cc[t2 as usize] + cc[t3 as usize + t2 as usize];
        ch[t1 as usize + t0 as usize] =
            cc[t2 as usize] - cc[t3 as usize + t2 as usize];
        t1 += ido;
        t2 = t1 << 1;
    }

    if ido < 2 {
        return;
    }

    if ido != 2 {
        t1 = 0;
        t2 = 0;
        for _ in 0..l1 {
            t3 = t1;
            let mut t4 = t2;
            let mut t5 = t4 + (ido << 1);
            let mut t6 = t0 + t1;
            for i in (2..ido as usize).step_by(2) {
                t3 += 2;
                t4 += 2;
                t5 -= 2;
                t6 += 2;
                ch[t3 as usize - 1] =
                    cc[t4 as usize - 1] + cc[t5 as usize - 1];
                let tr2 = cc[t4 as usize - 1] - cc[t5 as usize - 1];
                ch[t3 as usize] = cc[t4 as usize] - cc[t5 as usize];
                let ti2 = cc[t4 as usize] + cc[t5 as usize];
                ch[t6 as usize - 1] = wa1[i - 2] * tr2 - wa1[i - 1] * ti2;
                ch[t6 as usize] = wa1[i - 2] * ti2 + wa1[i - 1] * tr2;
            }
            t1 += ido;
            t2 = t1 << 1;
        }

        if ido % 2 == 1 {
            return;
        }
    }

    t1 = ido - 1;
    t2 = ido - 1;
    for _ in 0..l1 {
        ch[t1 as usize] = cc[t2 as usize] + cc[t2 as usize];
        ch[t1 as usize + t0 as usize] =
            -(cc[t2 as usize + 1] + cc[t2 as usize + 1]);
        t1 += ido;
        t2 += ido << 1;
    }
}

pub(crate) fn dradb3(
    ido: i32,
    l1: i32,
    cc: &[f32],
    ch: &mut [f32],
    wa1: &[f32],
    wa2: &[f32],
) {
    const TAUR: f32 = -0.5;
    const TAUI: f32 = 0.866_025_4;

    let t0 = l1 * ido;
    let mut t1 = 0;
    let t2 = t0 << 1;
    let mut t3 = ido << 1;
    let t4 = ido + (ido << 1);
    let mut t5 = 0;

    for _ in 0..l1 {
        let tr2 = cc[t3 as usize - 1] + cc[t3 as usize - 1];
        let cr2 = cc[t5 as usize] + TAUR * tr2;
        ch[t1 as usize] = cc[t5 as usize] + tr2;
        let ci3 = TAUI * (cc[t3 as usize] + cc[t3 as usize]);
        ch[t1 as usize + t0 as usize] = cr2 - ci3;
        ch[t1 as usize + t2 as usize] = cr2 + ci3;
        t1 += ido;
        t3 += t4;
        t5 += t4;
    }

    if ido == 1 {
        return;
    }

    t1 = 0;
    t3 = ido << 1;
    for _ in 0..l1 {
        let mut t7 = t1 + (t1 << 1);
        t5 = t7 + t3;
        let mut t6 = t5;
        let mut t8 = t1;
        let mut t9 = t1 + t0;
        let mut t10 = t9 + t0;
        for i in (2..ido as usize).step_by(2) {
            t5 += 2;
            t6 -= 2;
            t7 += 2;
            t8 += 2;
            t9 += 2;
            t10 += 2;
            let tr2 = cc[t5 as usize - 1] + cc[t6 as usize - 1];
            let cr2 = cc[t7 as usize - 1] + TAUR * tr2;
            ch[t8 as usize - 1] = cc[t7 as usize - 1] + tr2;
            let ti2 = cc[t5 as usize] - cc[t6 as usize];
            let ci2 = cc[t7 as usize] + TAUR * ti2;
            ch[t8 as usize] = cc[t7 as usize] + ti2;
            let cr3 = TAUI * (cc[t5 as usize - 1] - cc[t6 as usize - 1]);
            let ci3 = TAUI * (cc[t5 as usize] + cc[t6 as usize]);
            let dr2 = cr2 - ci3;
            let dr3 = cr2 + ci3;
            let di2 = ci2 + cr3;
            let di3 = ci2 - cr3;
            ch[t9 as usize - 1] = wa1[i - 2] * dr2 - wa1[i - 1] * di2;
            ch[t9 as usize] = wa1[i - 2] * di2 + wa1[i - 1] * dr2;
            ch[t10 as usize - 1] = wa2[i - 2] * dr3 - wa2[i - 1] * di3;
            ch[t10 as usize] = wa2[i - 2] * di3 + wa2[i - 1] * dr3;
        }
        t1 += ido;
    }
}

pub(crate) fn dradb4(
    ido: i32,
    l1: i32,
    cc: &[f32],
    ch: &mut [f32],
    wa1: &[f32],
    wa2: &[f32],
    wa3: &[f32],
) {
    const SQRT2: f32 = std::f32::consts::SQRT_2;

    let t0 = l1 * ido;
    let mut t1 = 0;
    let mut t2 = ido << 2;
    let mut t3 = 0;
    let t6 = ido << 1;

    for _ in 0..l1 {
        let mut t4 = t3 + t6;
        let mut t5 = t1;
        let tr3 = cc[t4 as usize - 1] + cc[t4 as usize - 1];
        let tr4 = cc[t4 as usize] + cc[t4 as usize];
        t4 += t6;
        let tr1 = cc[t3 as usize] - cc[t4 as usize - 1];
        let tr2 = cc[t3 as usize] + cc[t4 as usize - 1];
        ch[t5 as usize] = tr2 + tr3;
        t5 += t0;
        ch[t5 as usize] = tr1 - tr4;
        t5 += t0;
        ch[t5 as usize] = tr2 - tr3;
        t5 += t0;
        ch[t5 as usize] = tr1 + tr4;
        t1 += ido;
        t3 += t2;
    }

    if ido < 2 {
        return;
    }

    if ido != 2 {
        t1 = 0;
        for _ in 0..l1 {
            t2 = t1 << 2;
            t3 = t2 + t6;
            let mut t4 = t3;
            let mut t5 = t4 + t6;
            let mut t7 = t1;
            for i in (2..ido as usize).step_by(2) {
                t2 += 2;
                t3 += 2;
                t4 -= 2;
                t5 -= 2;
                t7 += 2;
                let ti1 = cc[t2 as usize] + cc[t5 as usize];
                let ti2 = cc[t2 as usize] - cc[t5 as usize];
                let ti3 = cc[t3 as usize] - cc[t4 as usize];
                let tr4 = cc[t3 as usize] + cc[t4 as usize];
                let tr1 = cc[t2 as usize - 1] - cc[t5 as usize - 1];
                let tr2 = cc[t2 as usize - 1] + cc[t5 as usize - 1];
                let ti4 = cc[t3 as usize - 1] - cc[t4 as usize - 1];
                let tr3 = cc[t3 as usize - 1] + cc[t4 as usize - 1];
                ch[t7 as usize - 1] = tr2 + tr3;
                let cr3 = tr2 - tr3;
                ch[t7 as usize] = ti2 + ti3;
                let ci3 = ti2 - ti3;
                let cr2 = tr1 - tr4;
                let cr4 = tr1 + tr4;
                let ci2 = ti1 + ti4;
                let ci4 = ti1 - ti4;
                let mut t8 = t7 + t0;
                ch[t8 as usize - 1] = wa1[i - 2] * cr2 - wa1[i - 1] * ci2;
                ch[t8 as usize] = wa1[i - 2] * ci2 + wa1[i - 1] * cr2;
                t8 += t0;
                ch[t8 as usize - 1] = wa2[i - 2] * cr3 - wa2[i - 1] * ci3;
                ch[t8 as usize] = wa2[i - 2] * ci3 + wa2[i - 1] * cr3;
                t8 += t0;
                ch[t8 as usize - 1] = wa3[i - 2] * cr4 - wa3[i - 1] * ci4;
                ch[t8 as usize] = wa3[i - 2] * ci4 + wa3[i - 1] * cr4;
            }
            t1 += ido;
        }

        if ido % 2 == 1 {
            return;
        }
    }

    t1 = ido;
    t2 = ido << 2;
    let mut t3 = ido - 1;
    let mut t4 = ido + (ido << 1);
    for _ in 0..l1 {
        let mut t5 = t3;
        let ti1 = cc[t1 as usize] + cc[t4 as usize];
        let ti2 = cc[t4 as usize] - cc[t1 as usize];
        let tr1 = cc[t1 as usize - 1] - cc[t4 as usize - 1];
        let tr2 = cc[t1 as usize - 1] + cc[t4 as usize - 1];
        ch[t5 as usize] = tr2 + tr2;
        t5 += t0;
        ch[t5 as usize] = SQRT2 * (tr1 - ti1);
        t5 += t0;
        ch[t5 as usize] = ti2 + ti2;
        t5 += t0;
        ch[t5 as usize] = -SQRT2 * (tr1 + ti1);
        t3 += ido;
        t1 += t2;
        t4 += t2;
    }
}

#[inline(always)]
fn dradbg_l102(
    ido: i32,
    nbd: i32,
    ipp2: i32,
    ipph: i32,
    l1: i32,
    t0: i32,
    t10: i32,
    cc: &[f32],
    ch: &mut [f32],
) {
    if ido != 1 {
        if nbd < l1 {
            let mut t1 = 0;
            let mut t2 = ipp2 * t0;
            let mut t7 = 0;
            for _ in 1..ipph {
                t1 += t0;
                t2 -= t0;
                let mut t3 = t1;
                let mut t4 = t2;
                t7 += ido << 1;
                let mut t8 = t7;
                let mut t9 = t7;
                for _ in (2..ido as usize).step_by(2) {
                    t3 += 2;
                    t4 += 2;
                    t8 += 2;
                    t9 -= 2;
                    let mut t5 = t3;
                    let mut t6 = t4;
                    let mut t11 = t8;
                    let mut t12 = t9;
                    for _ in 0..l1 {
                        ch[t5 as usize - 1] =
                            cc[t11 as usize - 1] + cc[t12 as usize - 1];
                        ch[t6 as usize - 1] =
                            cc[t11 as usize - 1] - cc[t12 as usize - 1];
                        ch[t5 as usize] = cc[t11 as usize] - cc[t12 as usize];
                        ch[t6 as usize] = cc[t11 as usize] + cc[t12 as usize];
                        t5 += ido;
                        t6 += ido;
                        t11 += t10;
                        t12 += t10;
                    }
                }
            }
        } else {
            let mut t1 = 0;
            let mut t2 = ipp2 * t0;
            let mut t7 = 0;
            for _ in 1..ipph {
                t1 += t0;
                t2 -= t0;
                let mut t3 = t1;
                let mut t4 = t2;
                t7 += ido << 1;
                let mut t8 = t7;
                for _ in 0..l1 {
                    let mut t5 = t3;
                    let mut t6 = t4;
                    let mut t9 = t8;
                    let mut t11 = t8;
                    for _ in (2..ido as usize).step_by(2) {
                        t5 += 2;
                        t6 += 2;
                        t9 += 2;
                        t11 -= 2;
                        ch[t5 as usize - 1] =
                            cc[t9 as usize - 1] + cc[t11 as usize - 1];
                        ch[t6 as usize - 1] =
                            cc[t9 as usize - 1] - cc[t11 as usize - 1];
                        ch[t5 as usize] = cc[t9 as usize] - cc[t11 as usize];
                        ch[t6 as usize] = cc[t9 as usize] + cc[t11 as usize];
                    }
                    t3 += ido;
                    t4 += ido;
                    t8 += t10;
                }
            }
        }
    }
}

#[inline(always)]
fn dradbg_l103(
    ido: i32,
    nbd: i32,
    ipp2: i32,
    ipph: i32,
    l1: i32,
    t0: i32,
    cc: &[f32],
    ch: &mut [f32],
) {
    if ido != 1 {
        if nbd < l1 {
            let mut t1 = 0;
            let mut t2 = ipp2 * t0;
            for _ in 1..ipph {
                t1 += t0;
                t2 -= t0;
                let mut t3 = t1;
                let mut t4 = t2;
                for _ in (2..ido as usize).step_by(2) {
                    t3 += 2;
                    t4 += 2;
                    let mut t5 = t3;
                    let mut t6 = t4;
                    for _ in 0..l1 {
                        ch[t5 as usize - 1] =
                            cc[t5 as usize - 1] - cc[t6 as usize];
                        ch[t6 as usize - 1] =
                            cc[t5 as usize - 1] + cc[t6 as usize];
                        ch[t5 as usize] =
                            cc[t5 as usize] + cc[t6 as usize - 1];
                        ch[t6 as usize] =
                            cc[t5 as usize] - cc[t6 as usize - 1];
                        t5 += ido;
                        t6 += ido;
                    }
                }
            }
        } else {
            let mut t1 = 0;
            let mut t2 = ipp2 * t0;
            for _ in 1..ipph {
                t1 += t0;
                t2 -= t0;
                let mut t3 = t1;
                let mut t4 = t2;
                for _ in 0..l1 {
                    let mut t5 = t3;
                    let mut t6 = t4;
                    for _ in (2..ido as usize).step_by(2) {
                        t5 += 2;
                        t6 += 2;
                        ch[t5 as usize - 1] =
                            cc[t5 as usize - 1] - cc[t6 as usize];
                        ch[t6 as usize - 1] =
                            cc[t5 as usize - 1] + cc[t6 as usize];
                        ch[t5 as usize] =
                            cc[t5 as usize] + cc[t6 as usize - 1];
                        ch[t6 as usize] =
                            cc[t5 as usize] - cc[t6 as usize - 1];
                    }
                    t3 += ido;
                    t4 += ido;
                }
            }
        }
    }
}

#[inline(always)]
fn dradbg_l104(
    ido: i32,
    nbd: i32,
    ip: i32,
    l1: i32,
    t0: i32,
    cc: &mut [f32],
    ch: &[f32],
    wa: &[f32],
) {
    if nbd > l1 {
        let mut is = -ido - 1;
        let mut t1 = 0;
        for _ in 1..ip {
            is += ido;
            t1 += t0;
            let mut t2 = t1;
            for _ in 0..l1 {
                let mut idij = is;
                let mut t3 = t2;
                for _ in (2..ido as usize).step_by(2) {
                    idij += 2;
                    t3 += 2;
                    cc[t3 as usize - 1] = wa[idij as usize - 1]
                        * ch[t3 as usize - 1]
                        - wa[idij as usize] * ch[t3 as usize];
                    cc[t3 as usize] = wa[idij as usize - 1] * ch[t3 as usize]
                        + wa[idij as usize] * ch[t3 as usize - 1];
                }
                t2 += ido;
            }
        }
    } else {
        let mut is = -ido - 1;
        let mut t1 = 0;
        for _ in 1..ip {
            is += ido;
            t1 += t0;
            let mut idij = is;
            let mut t2 = t1;
            for _ in (2..ido as usize).step_by(2) {
                t2 += 2;
                idij += 2;
                let mut t3 = t2;
                for _ in 0..l1 {
                    cc[t3 as usize - 1] = wa[idij as usize - 1]
                        * ch[t3 as usize - 1]
                        - wa[idij as usize] * ch[t3 as usize];
                    cc[t3 as usize] = wa[idij as usize - 1] * ch[t3 as usize]
                        + wa[idij as usize] * ch[t3 as usize - 1];
                    t3 += ido;
                }
            }
        }
    }
}

pub(crate) fn dradbg(
    ido: i32,
    ip: i32,
    l1: i32,
    idl1: i32,
    cc: &mut [f32],
    ch: &mut [f32],
    wa: &[f32],
) {
    const TPI: f32 = 6.283_185_5;

    let t10 = ip * ido;
    let t0 = l1 * ido;
    let arg = TPI / ip as f32;
    let dcp = f64::cos(arg as f64) as f32;
    let dsp = f64::sin(arg as f64) as f32;
    let nbd = (ido - 1) >> 1_i32;
    let ipp2 = ip;
    let ipph = (ip + 1) >> 1_i32;

    if ido < l1 {
        for t1 in 0..ido {
            let mut t2 = t1;
            let mut t3 = t1;
            for _ in 0..l1 {
                ch[t2 as usize] = cc[t3 as usize];
                t2 += ido;
                t3 += t10;
            }
        }
    } else {
        let mut t1 = 0;
        let mut t2 = 0;
        for _ in 0..l1 {
            let mut t3 = t1;
            let mut t4 = t2;
            for _ in 0..ido {
                ch[t3 as usize] = cc[t4 as usize];
                t3 += 1;
                t4 += 1;
            }
            t1 += ido;
            t2 += t10;
        }
    }

    let mut t1 = 0;
    let mut t2 = ipp2 * t0;
    let mut t5 = ido << 1;
    let t7 = t5;

    for _ in 1..ipph {
        t1 += t0;
        t2 -= t0;
        let mut t3 = t1;
        let mut t4 = t2;
        let mut t6 = t5;
        for _ in 0..l1 {
            ch[t3 as usize] = cc[t6 as usize - 1] + cc[t6 as usize - 1];
            ch[t4 as usize] = cc[t6 as usize] + cc[t6 as usize];
            t3 += ido;
            t4 += ido;
            t6 += t10;
        }
        t5 += t7;
    }

    dradbg_l102(ido, nbd, ipp2, ipph, l1, t0, t10, cc, ch);

    let ar1: f32 = 1.0;
    let ai1: f32 = 0.0;
    let mut t1 = 0;
    let mut t2 = ipp2 * idl1;
    let t9 = t2;
    let t3 = (ip - 1) * idl1;

    for _ in 1..ipph {
        t1 += idl1;
        t2 -= idl1;
        let ar1h = dcp * ar1 - dsp * ai1;
        let ai1 = dcp * ai1 + dsp * ar1;
        let ar1 = ar1h;
        let mut t4 = t1;
        let mut t5 = t2;
        let mut t7 = idl1;
        let mut t8 = t3;

        for t6 in 0..idl1 {
            let fresh13 = t6;
            let fresh14 = t7;
            t7 += 1;
            let fresh15 = t4;
            t4 += 1;
            cc[fresh15 as usize] =
                ch[fresh13 as usize] + ar1 * ch[fresh14 as usize];
            let fresh16 = t8;
            t8 += 1;
            let fresh17 = t5;
            t5 += 1;
            cc[fresh17 as usize] = ai1 * ch[fresh16 as usize];
        }
        let dcc = ar1;
        let ds2 = ai1;
        let ar2 = ar1;
        let ai2 = ai1;
        let mut t6 = idl1;
        let mut t7 = t9 - idl1;

        for _ in 2..ipph {
            t6 += idl1;
            t7 -= idl1;
            let ar2h = dcc * ar2 - ds2 * ai2;
            let ai2 = dcc * ai2 + ds2 * ar2;
            let ar2 = ar2h;
            let mut t4 = t1;
            let mut t5 = t2;
            let mut t11 = t6;
            let mut t12 = t7;

            for _ in 0..idl1 {
                let fresh18 = t11;
                t11 += 1;
                let fresh19 = t4;
                t4 += 1;
                cc[fresh19 as usize] += ar2 * ch[fresh18 as usize];
                let fresh20 = t12;
                t12 += 1;
                let fresh21 = t5;
                t5 += 1;
                cc[fresh21 as usize] += ai2 * ch[fresh20 as usize];
            }
        }
    }

    t1 = 0;
    for _ in 1..ipph {
        t1 += idl1;
        let mut t2 = t1;
        for ik in 0..idl1 {
            let fresh22 = t2;
            t2 += 1;
            ch[ik as usize] += ch[fresh22 as usize];
        }
    }

    t1 = 0;
    t2 = ipp2 * t0;
    for _ in 1..ipph {
        t1 += t0;
        t2 -= t0;
        let mut t3 = t1;
        let mut t4 = t2;
        for _ in 0..l1 {
            ch[t3 as usize] = cc[t3 as usize] - cc[t4 as usize];
            ch[t4 as usize] = cc[t3 as usize] + cc[t4 as usize];
            t3 += ido;
            t4 += ido;
        }
    }

    dradbg_l103(ido, nbd, ipp2, ipph, l1, t0, cc, ch);

    if ido == 1 {
        return;
    }

    for ik in 0..idl1 {
        cc[ik as usize] = ch[ik as usize];
    }

    t1 = 0;
    for _ in 1..ip {
        t1 += t0;
        t2 = t1;
        for _ in 0..l1 {
            cc[t2 as usize] = ch[t2 as usize];
            t2 += ido;
        }
    }

    dradbg_l104(ido, nbd, ip, l1, t0, cc, ch, wa);
}
