use crate::dradb::*;
use crate::dradf::*;

#[derive(Clone)]
pub struct DrftLookup {
    pub n: usize,
    pub trigcache: Vec<f32>,
    pub splitcache: Vec<i32>,
}

impl DrftLookup {
    pub fn new(n: usize) -> Self {
        let mut drft = Self {
            n: n,
            trigcache: vec![0.0; 3 * n],
            splitcache: vec![0; 32],
        };

        fdrffti(n, &mut drft.trigcache, &mut drft.splitcache);

        drft
    }
}

#[inline(always)]
fn drfti1_c_10244(ifac: &mut [i32], n: i32, nf: &mut i32) {
    const NTRYH: [i32; 4] = [4, 2, 3, 5];

    let mut j: i32 = -1;
    let mut ntry = 0;
    let mut nl = n;

    'c_10244: loop {
        j += 1;
        if j < 4 {
            ntry = NTRYH[j as usize]
        } else {
            ntry += 2
        }
        loop {
            let nq = nl / ntry;
            let nr = nl - ntry * nq;
            if nr != 0 {
                break;
            }
            *nf += 1;
            ifac[*nf as usize + 1] = ntry;
            nl = nq;
            if ntry == 2 {
                if *nf != 1 {
                    for i in 1..*nf {
                        let ib = *nf - i + 1 as i32;
                        ifac[ib as usize + 1] = ifac[ib as usize];
                    }
                    ifac[2] = 2;
                }
            }
            if nl == 1 {
                break 'c_10244;
            }
        }
    }
}

fn drfti1(wa: &mut [f32], ifac: &mut [i32]) {
    const TPI: f32 = 6.283_185_307_179_586_48;

    let n = wa.len() as i32;
    let mut nf = 0;

    drfti1_c_10244(ifac, n, &mut nf);

    ifac[0] = n;
    ifac[1] = nf;
    let nfm1 = nf - 1;
    let argh = TPI / n as f32;
    let mut is = 0;
    let mut l1 = 1;

    if nfm1 == 0 {
        return;
    }

    for k1 in 0..nfm1 {
        let ip = ifac[k1 as usize + 2];
        let l2 = l1 * ip;
        let ido = n / l2;
        let ipm = ip - 1;
        let mut ld = 0;
        let mut j = 0;
        while j < ipm {
            ld += l1;
            let argld = ld as f32 * argh;
            let mut i = is;
            let mut ii = 2;
            let mut fi = 0.0f32;
            while ii < ido {
                fi += 1.0f32;
                let arg = fi * argld;
                let fresh0 = i;
                i = i + 1;
                wa[fresh0] = f64::cos(arg as f64) as f32;
                let fresh1 = i;
                i = i + 1;
                wa[fresh1] = f64::sin(arg as f64) as f32;
                ii += 2 as i32
            }
            is += ido as usize;
            j += 1
        }
        l1 = l2;
    }
}

fn fdrffti(n: usize, wsave: &mut [f32], ifac: &mut [i32]) {
    if n == 1 {
        return;
    }
    drfti1(&mut wsave[n..n * 2], ifac);
}

#[inline(always)]
fn drftf1_l102(
    ip: i32,
    ido: i32,
    l1: i32,
    idl1: i32,
    iw: i32,
    na: &mut i32,
    c: &mut [f32],
    ch: &mut [f32],
    wa: &mut [f32],
) {
    let iw_temp = iw as usize - 1;
    if ip != 4 {
        if ip != 2 {
            if ido == 1 {
                *na = 1 - *na;
            }
            if *na != 0 {
                dradfg(ido, ip, l1, idl1, ch, c, &wa[iw_temp..]);
                *na = 0;
            } else {
                dradfg(ido, ip, l1, idl1, c, ch, &wa[iw_temp..]);
                *na = 1;
            }
        } else if *na != 0 {
            dradf2(ido, l1, ch, c, &mut wa[iw_temp..]);
        } else {
            dradf2(ido, l1, c, ch, &mut wa[iw_temp..]);
        }
    }
}

fn drftf1(
    n: i32,
    c: &mut [f32],
    ch: &mut [f32],
    wa: &mut [f32],
    ifac: &mut [i32],
) {
    let nf = ifac[1];
    let mut na = 1;
    let mut l2 = n;
    let mut iw = n;

    for k1 in 0..nf {
        let kh = nf - k1;
        let ip = ifac[kh as usize + 1];
        let l1 = l2 / ip;
        let ido = n / l2;
        let idl1 = ido * l1;
        iw -= (ip - 1) * ido;
        na = 1 - na;

        if ip != 4 {
            drftf1_l102(ip, ido, l1, idl1, iw, &mut na, c, ch, wa);
        } else {
            let ix2 = iw + ido;
            let ix3 = ix2 + ido;
            if na != 0 {
                dradf4(
                    ido,
                    l1,
                    ch,
                    c,
                    &wa[(iw as usize - 1)..],
                    &wa[(ix2 as usize - 1)..],
                    &wa[(ix3 as usize - 1)..],
                );
            } else {
                dradf4(
                    ido,
                    l1,
                    c,
                    ch,
                    &wa[(iw as usize - 1)..],
                    &wa[(ix2 as usize - 1)..],
                    &wa[(ix3 as usize - 1)..],
                );
            }
        }

        l2 = l1;
    }

    if na == 1 {
        return;
    }

    c[..n as usize].copy_from_slice(&ch[..n as usize]);
}

#[inline(always)]
fn drftb1_l102(
    ip: i32,
    ido: i32,
    l1: i32,
    idl1: i32,
    iw: i32,
    na: &mut i32,
    c: &mut [f32],
    ch: &mut [f32],
    wa: &mut [f32],
) {
    let iw_temp = iw as usize - 1;
    if ip != 2 {
        if ip != 3 {
            if *na != 0 {
                dradbg(ido, ip, l1, idl1, ch, c, &wa[iw_temp..]);
            } else {
                dradbg(ido, ip, l1, idl1, c, ch, &wa[iw_temp..]);
            }

            if ido == 1 {
                *na = 1 - *na;
            }
        } else {
            let ix2 = iw + ido - 1;
            if *na != 0 {
                dradb3(ido, l1, &ch, c, &wa[iw_temp..], &wa[ix2 as usize..]);
            } else {
                dradb3(ido, l1, &c, ch, &wa[iw_temp..], &wa[ix2 as usize..]);
            }
            *na = 1 - *na;
        }
    } else {
        if *na != 0 {
            dradb2(ido, l1, &ch, c, &wa[iw_temp..]);
        } else {
            dradb2(ido, l1, &c, ch, &wa[iw_temp..]);
        }
        *na = 1 - *na;
    }
}

fn drftb1(
    n: i32,
    c: &mut [f32],
    ch: &mut [f32],
    wa: &mut [f32],
    ifac: &mut [i32],
) {
    let nf = ifac[1];
    let mut l1 = 1;
    let mut na = 0;
    let mut iw = 1;

    for k1 in 0..nf {
        let ip = ifac[k1 as usize + 2];
        let l2 = ip * l1;
        let ido = n / l2;
        let idl1 = ido * l1;

        if ip != 4 {
            drftb1_l102(ip, ido, l1, idl1, iw, &mut na, c, ch, wa);
        } else {
            let ix2 = iw + ido;
            let ix3 = ix2 + ido;
            if na != 0 {
                dradb4(
                    ido,
                    l1,
                    &ch,
                    c,
                    &wa[(iw as usize - 1)..],
                    &wa[(ix2 as usize - 1)..],
                    &wa[(ix3 as usize - 1)..],
                );
            } else {
                dradb4(
                    ido,
                    l1,
                    &c,
                    ch,
                    &wa[(iw as usize - 1)..],
                    &wa[(ix2 as usize - 1)..],
                    &wa[(ix3 as usize - 1)..],
                );
            }
            na = 1 - na
        }

        l1 = l2;
        iw += (ip - 1) * ido;
    }

    if na == 0 {
        return;
    }

    c[..n as usize].copy_from_slice(&ch[..n as usize]);
}

pub fn spx_drft_forward(l: &mut DrftLookup, data: &mut [f32]) {
    if l.n == 1 {
        return;
    }

    let mut trigcache_temp = l.trigcache[l.n as usize..].to_vec();

    drftf1(
        l.n as i32,
        data,
        &mut l.trigcache,
        &mut trigcache_temp,
        &mut l.splitcache,
    );

    l.trigcache[l.n as usize..].copy_from_slice(&trigcache_temp);
}

pub fn spx_drft_backward(l: &mut DrftLookup, data: &mut [f32]) {
    if l.n == 1 {
        return;
    }

    let mut trigcache_temp = l.trigcache[l.n as usize..].to_vec();

    drftb1(
        l.n as i32,
        data,
        &mut l.trigcache,
        &mut trigcache_temp,
        &mut l.splitcache,
    );

    l.trigcache[l.n as usize..].copy_from_slice(&trigcache_temp);
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f32 = 1e-6;

    #[test]
    fn fdrffti_simple() {
        let mut trigcache = [42.; 3];
        let mut splitcache = [24; 32];

        fdrffti(1, &mut trigcache, &mut splitcache);

        assert!(trigcache.iter().all(|&x| (x - 42.).abs() < EPSILON));
        assert!(splitcache.iter().all(|&x| x == 24));
    }
}
