#[allow(clippy::too_many_arguments)]
#[target_feature(enable = "avx")]
pub unsafe fn interpolate_step_single(
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
    crate::speex::avx::interpolate_step_single(
        in_slice, out_slice, out_stride, out_sample, oversample, offset, n,
        sinc_table, frac,
    );
}

#[allow(clippy::too_many_arguments)]
#[target_feature(enable = "avx")]
pub unsafe fn interpolate_step_double(
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
    crate::speex::avx::interpolate_step_double(
        in_slice, out_slice, out_stride, out_sample, oversample, offset, n,
        sinc_table, frac,
    );
}

#[target_feature(enable = "avx")]
pub unsafe fn direct_step_single(
    in_slice: &[f32],
    out_slice: &mut [f32],
    out_stride: usize,
    out_sample: usize,
    n: usize,
    sinc_table: &[f32],
) {
    crate::speex::avx::direct_step_single(
        in_slice, out_slice, out_stride, out_sample, n, sinc_table,
    );
}

#[target_feature(enable = "avx")]
pub unsafe fn direct_step_double(
    in_slice: &[f32],
    out_slice: &mut [f32],
    out_stride: usize,
    out_sample: usize,
    n: usize,
    sinc_table: &[f32],
) {
    crate::speex::avx::direct_step_double(
        in_slice, out_slice, out_stride, out_sample, n, sinc_table,
    );
}
