#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

#[inline(always)]
pub unsafe fn swap_lo_and_hi(v: __m128) -> __m128 {
    _mm_shuffle_ps(v, v, 0b01001110)
}

#[inline(always)]
pub unsafe fn sum_m128_into_m128d(v: __m128) -> __m128d {
    let a = _mm_cvtps_pd(v);
    let b = _mm_cvtps_pd(swap_lo_and_hi(v));
    _mm_add_pd(a, b)
}

#[inline(always)]
pub unsafe fn split_m128_into_m128d(v: __m128) -> (__m128d, __m128d) {
    let a = _mm_cvtps_pd(v);
    let b = _mm_cvtps_pd(swap_lo_and_hi(v));
    (a, b)
}

pub unsafe fn hsum_m128d(v: __m128d) -> f64 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    _mm_cvtsd_f64(_mm_hadd_pd(v, v))
}

#[inline(always)]
pub unsafe fn cubic_coef(frac: f32, interp: &mut std::arch::x86_64::__m128) {
    let mut interp_v = [0.0, 0.0, 0.0, 0.0];

    interp_v[0] = -0.166_67 * frac + 0.166_67 * frac * frac * frac;
    interp_v[1] = frac + 0.5 * frac * frac - 0.5f32 * frac * frac * frac;
    interp_v[3] =
        -0.333_33 * frac + 0.5 * frac * frac - 0.166_67 * frac * frac * frac;
    interp_v[2] = (1.0f64
        - interp_v[0] as f64
        - interp_v[1] as f64
        - interp_v[3] as f64) as f32;

    *interp = _mm_loadu_ps(interp_v.as_ptr());
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
    unsafe {
        let mut accum = _mm_setzero_ps();
        in_slice.iter().zip(0..n).for_each(|(&curr_in, j)| {
            let idx = (2 + (j + 1) * oversample as usize) - offset as usize;
            let sinc_ptr: *const f32 = sinc_table[idx..].as_ptr();
            accum = _mm_add_ps(
                accum,
                _mm_mul_ps(_mm_loadu_ps(sinc_ptr), _mm_set1_ps(curr_in)),
            );
        });
        let mut interp = _mm_setzero_ps();
        cubic_coef(frac, &mut interp);
        let v = _mm_mul_ps(interp, accum);
        out_slice[(out_stride * out_sample) as usize] =
            _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(v, v), v));
    }
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
    unsafe {
        let mut accum_lo = _mm_setzero_pd();
        let mut accum_hi = _mm_setzero_pd();
        in_slice.iter().zip(0..n).for_each(|(&curr_in, j)| {
            let idx = (2 + (j + 1) * oversample as usize) - offset as usize;
            let sinct_ptr: *const f32 = sinc_table[idx..].as_ptr();
            let v = _mm_mul_ps(_mm_loadu_ps(sinct_ptr), _mm_set1_ps(curr_in));
            let (v64_lo, v64_hi) = split_m128_into_m128d(v);
            accum_lo = _mm_add_pd(accum_lo, v64_lo);
            accum_hi = _mm_add_pd(accum_hi, v64_hi);
        });

        let mut interp = _mm_setzero_ps();
        cubic_coef(frac, &mut interp);

        let accum32 = _mm_add_ps(
            _mm_cvtpd_ps(accum_lo),
            swap_lo_and_hi(_mm_cvtpd_ps(accum_hi)),
        );
        let v = _mm_mul_ps(accum32, interp);

        out_slice[(out_stride * out_sample) as usize] =
            _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(v, v), v));
    }
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
    unsafe {
        let accum = sinc_table
            .chunks_exact(8)
            .zip(in_slice.chunks_exact(8))
            .take((n as usize) / 8)
            .fold(_mm_setzero_ps(), |acc, (sinct_p, iptr_p)| {
                let sinct_v = _mm_loadu_ps(sinct_p.as_ptr());
                let iptr_v = _mm_loadu_ps(iptr_p.as_ptr());

                let acc = _mm_add_ps(acc, _mm_mul_ps(sinct_v, iptr_v));

                let sinct_v = _mm_loadu_ps(sinct_p[4..].as_ptr());
                let iptr_v = _mm_loadu_ps(iptr_p[4..].as_ptr());
                _mm_add_ps(acc, _mm_mul_ps(sinct_v, iptr_v))
            });

        let accum = _mm_add_ps(accum, _mm_movehl_ps(accum, accum));
        let accum =
            _mm_add_ss(accum, _mm_shuffle_ps(accum, accum, 0b01010101));

        _mm_store_ss(
            out_slice[((out_stride * out_sample) as usize)..].as_mut_ptr(),
            accum,
        )
    }
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
    unsafe {
        let mut accum = _mm_setzero_pd();
        let mut j = 0;

        while j < n {
            let sinct_v = _mm_loadu_ps(sinc_table[j as usize..].as_ptr());
            let iptr_v = _mm_loadu_ps(in_slice[j as usize..].as_ptr());
            let v = _mm_mul_ps(sinct_v, iptr_v);

            accum = _mm_add_pd(accum, _mm_cvtps_pd(v));
            j += 2;
        }

        out_slice[(out_stride * out_sample) as usize] =
            hsum_m128d(accum) as f32;
    }
}

#[cfg(test)]
mod tests {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    #[test]
    fn test_swap_lo_and_hi() {
        let input: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0];
        let mut output: Vec<f32> = vec![0.0; 4];

        unsafe {
            let v = _mm_loadu_ps(input.as_ptr());
            let result = super::swap_lo_and_hi(v);
            _mm_storeu_ps(output.as_mut_ptr(), result)
        }

        assert_eq!(output, vec![3.0, 4.0, 1.0, 2.0]);
    }

    #[test]
    fn test_sum_m128_into_m128d() {
        let input: Vec<f32> = vec![1.0, 2.0, 3.0, 4.0];
        let mut output: Vec<f64> = vec![0.0; 2];

        unsafe {
            let v = _mm_loadu_ps(input.as_ptr());
            let result = super::sum_m128_into_m128d(v);
            _mm_storeu_pd(output.as_mut_ptr(), result)
        }

        assert_eq!(output, vec![4.0, 6.0]);
    }

    #[test]
    fn test_hsum_m128d() {
        let input: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0];

        let result = unsafe {
            let v = _mm_loadu_pd(input.as_ptr());
            super::hsum_m128d(v)
        };

        assert_eq!(result as usize, 3 as usize);
    }
}
