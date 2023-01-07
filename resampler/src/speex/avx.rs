#[cfg(target_arch = "x86")]
use std::arch::x86::*;
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

#[inline(always)]
unsafe fn hsum_m256d(v: __m256d) -> f64 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    let v1 = _mm256_hadd_pd(v, v);
    let v2 = _mm256_add_pd(v1, _mm256_permute2f128_pd(v1, v1, 1));
    _mm256_cvtsd_f64(v2)
}

use super::sse3::cubic_coef;

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
            let idx = (2 + (j + 1) * oversample) - offset;
            let sinc_ptr: *const f32 = sinc_table[idx..].as_ptr();
            accum = _mm_add_ps(
                accum,
                _mm_mul_ps(_mm_loadu_ps(sinc_ptr), _mm_set1_ps(curr_in)),
            );
        });
        let mut interp = _mm_setzero_ps();
        cubic_coef(frac, &mut interp);
        let v = _mm_mul_ps(interp, accum);
        out_slice[(out_stride * out_sample)] =
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
        let mut accum = _mm256_setzero_pd();
        in_slice.iter().zip(0..n).for_each(|(&curr_in, j)| {
            let idx = (2 + (j + 1) * oversample) - offset;
            let sinct_ptr: *const f32 = sinc_table[idx..].as_ptr();
            let v = _mm_mul_ps(_mm_loadu_ps(sinct_ptr), _mm_set1_ps(curr_in));
            accum = _mm256_add_pd(accum, _mm256_cvtps_pd(v));
        });

        let mut interp = _mm_setzero_ps();
        cubic_coef(frac, &mut interp);

        let accum = _mm256_mul_pd(accum, _mm256_cvtps_pd(interp));

        out_slice[(out_stride * out_sample)] = hsum_m256d(accum) as f32;
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
            .take(n / 8)
            .fold(_mm256_setzero_ps(), |acc, (sinct_p, iptr_p)| {
                let sinct_v = _mm256_loadu_ps(sinct_p.as_ptr());
                let iptr_v = _mm256_loadu_ps(iptr_p.as_ptr());

                _mm256_add_ps(acc, _mm256_mul_ps(sinct_v, iptr_v))
            });

        let accum = _mm256_hadd_ps(accum, accum);
        let accum =
            _mm256_add_ps(accum, _mm256_permute2f128_ps(accum, accum, 1));
        let accum = _mm256_hadd_ps(accum, accum);

        _mm_store_ss(
            out_slice[(out_stride * out_sample)..].as_mut_ptr(),
            _mm256_castps256_ps128(accum),
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
        let mut accum = _mm256_setzero_pd();
        let mut j = 0;

        while j < n {
            let sinct_v = _mm_loadu_ps(sinc_table[j..].as_ptr());
            let iptr_v = _mm_loadu_ps(in_slice[j..].as_ptr());
            let v = _mm_mul_ps(sinct_v, iptr_v);

            accum = _mm256_add_pd(accum, _mm256_cvtps_pd(v));
            j += 4;
        }

        out_slice[(out_stride * out_sample)] = hsum_m256d(accum) as f32;
    }
}

#[cfg(test)]
mod tests {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    #[test]
    fn test_hsum_m256d() {
        let input: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0];

        let result = unsafe {
            let v = _mm256_loadu_pd(input.as_ptr());
            super::hsum_m256d(v)
        };

        assert_eq!(result as usize, 10_usize);
    }
}
