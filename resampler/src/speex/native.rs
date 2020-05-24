#[inline(always)]
pub fn cubic_coef(frac: f32, interp: &mut [f32]) {
    interp[0] = -0.166_67 * frac + 0.166_67 * frac * frac * frac;
    interp[1] = frac + 0.5 * frac * frac - 0.5f32 * frac * frac * frac;
    interp[3] =
        -0.333_33 * frac + 0.5 * frac * frac - 0.166_67 * frac * frac * frac;
    interp[2] =
        (1.0f64 - interp[0] as f64 - interp[1] as f64 - interp[3] as f64)
            as f32;
}

#[allow(clippy::too_many_arguments)]
#[inline(always)]
pub fn interpolate_step_single(
    in_slice: &[f32],
    out_slice: &mut [f32],
    out_stride: usize,
    out_sample: usize,
    oversample: usize,
    offset: usize,
    n: usize,
    sinc_table: &[f32],
    frac: f32,
) {
    let mut accum: [f32; 4] = [0.; 4];
    in_slice.iter().zip(0..n).for_each(|(&curr_in, j)| {
        let idx = (2 + (j + 1) * oversample as usize) - offset as usize;
        accum.iter_mut().zip(sinc_table.iter().skip(idx)).for_each(
            |(v, &s)| {
                *v += curr_in * s;
            },
        );
    });
    let mut interp: [f32; 4] = [0.; 4];
    cubic_coef(frac, &mut interp);
    out_slice[(out_stride * out_sample) as usize] = interp
        .iter()
        .zip(accum.iter())
        .map(|(&x, &y)| x * y)
        .fold(0., |acc, x| acc + x);
}

#[allow(clippy::too_many_arguments)]
#[inline(always)]
pub fn interpolate_step_double(
    in_slice: &[f32],
    out_slice: &mut [f32],
    out_stride: usize,
    out_sample: usize,
    oversample: usize,
    offset: usize,
    n: usize,
    sinc_table: &[f32],
    frac: f32,
) {
    let mut accum: [f64; 4] = [0.0; 4];
    in_slice.iter().zip(0..n).for_each(|(&curr_in, j)| {
        let idx = (2 + (j + 1) * oversample as usize) - offset as usize;
        accum.iter_mut().zip(sinc_table.iter().skip(idx)).for_each(
            |(v, &s)| {
                *v += (curr_in * s) as f64;
            },
        );
    });
    let mut interp: [f32; 4] = [0.; 4];
    cubic_coef(frac, &mut interp);
    out_slice[(out_stride * out_sample) as usize] = interp
        .iter()
        .zip(accum.iter())
        .map(|(&x, &y)| x * y as f32)
        .fold(0., |acc, x| acc + x);
}

#[inline(always)]
pub fn direct_step_single(
    in_slice: &[f32],
    out_slice: &mut [f32],
    out_stride: usize,
    out_sample: usize,
    n: usize,
    sinc_table: &[f32],
) {
    let mut sum: f32 = 0.0;
    let mut j = 0;
    while j < n {
        sum += sinc_table[j as usize] * in_slice[j as usize];
        j += 1
    }
    out_slice[(out_stride * out_sample) as usize] = sum;
}

#[inline(always)]
pub fn direct_step_double(
    in_slice: &[f32],
    out_slice: &mut [f32],
    out_stride: usize,
    out_sample: usize,
    n: usize,
    sinc_table: &[f32],
) {
    let mut accum: [f64; 4] = [0.0; 4];
    let mut j = 0;

    while j < n {
        accum[0usize] +=
            f64::from(sinc_table[j as usize] * in_slice[j as usize]);
        accum[1usize] += f64::from(
            sinc_table[(j + 1) as usize] * in_slice[(j + 1) as usize],
        );
        accum[2usize] += f64::from(
            sinc_table[(j + 2) as usize] * in_slice[(j + 2) as usize],
        );
        accum[3usize] += f64::from(
            sinc_table[(j + 3) as usize] * in_slice[(j + 3) as usize],
        );
        j += 4
    }
    let sum: f64 =
        accum[0usize] + accum[1usize] + accum[2usize] + accum[3usize];
    out_slice[(out_stride * out_sample) as usize] = sum as f32;
}
