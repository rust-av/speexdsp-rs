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
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx") {
            unsafe {
                avx_wrapper::interpolate_step_single(
                    in_slice, out_slice, out_stride, out_sample, oversample,
                    offset, n, sinc_table, frac,
                );
            }
            return;
        }

        if is_x86_feature_detected!("sse3") {
            unsafe {
                sse3_wrapper::interpolate_step_single(
                    in_slice, out_slice, out_stride, out_sample, oversample,
                    offset, n, sinc_table, frac,
                );
            }
            return;
        }
    }

    super::native::interpolate_step_single(
        in_slice, out_slice, out_stride, out_sample, oversample, offset, n,
        sinc_table, frac,
    );
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
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx") {
            unsafe {
                avx_wrapper::interpolate_step_double(
                    in_slice, out_slice, out_stride, out_sample, oversample,
                    offset, n, sinc_table, frac,
                );
            }
            return;
        }

        if is_x86_feature_detected!("sse3") {
            unsafe {
                sse3_wrapper::interpolate_step_double(
                    in_slice, out_slice, out_stride, out_sample, oversample,
                    offset, n, sinc_table, frac,
                );
            }
            return;
        }
    }

    super::native::interpolate_step_double(
        in_slice, out_slice, out_stride, out_sample, oversample, offset, n,
        sinc_table, frac,
    );
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
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx") {
            unsafe {
                avx_wrapper::direct_step_single(
                    in_slice, out_slice, out_stride, out_sample, n, sinc_table,
                );
            }
            return;
        }

        if is_x86_feature_detected!("sse3") {
            unsafe {
                sse3_wrapper::direct_step_single(
                    in_slice, out_slice, out_stride, out_sample, n, sinc_table,
                );
            }
            return;
        }
    }

    super::native::direct_step_single(
        in_slice, out_slice, out_stride, out_sample, n, sinc_table,
    );
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
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx") {
            unsafe {
                avx_wrapper::direct_step_double(
                    in_slice, out_slice, out_stride, out_sample, n, sinc_table,
                );
            }
            return;
        }

        if is_x86_feature_detected!("sse3") {
            unsafe {
                sse3_wrapper::direct_step_double(
                    in_slice, out_slice, out_stride, out_sample, n, sinc_table,
                );
            }
            return;
        }
    }

    super::native::direct_step_double(
        in_slice, out_slice, out_stride, out_sample, n, sinc_table,
    );
}

mod avx_wrapper;
mod sse3_wrapper;
