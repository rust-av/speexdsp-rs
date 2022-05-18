#![allow(
    dead_code,
    mutable_transmutes,
    non_camel_case_types,
    non_snake_case,
    non_upper_case_globals,
    unused_mut
)]

use libc::{calloc, exit, free, realloc};

pub const RESAMPLER_ERR_INVALID_ARG: unnamed = 3;
pub type __uint64_t = libc::c_ulong;
pub const RESAMPLER_ERR_OVERFLOW: unnamed = 5;
pub type uint32_t = libc::c_uint;
pub type _LIB_VERSION_TYPE = libc::c_int;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct __va_list_tag {
    pub gp_offset: libc::c_uint,
    pub fp_offset: libc::c_uint,
    pub overflow_arg_area: *mut libc::c_void,
    pub reg_save_area: *mut libc::c_void,
}
pub const RESAMPLER_ERR_PTR_OVERLAP: unnamed = 4;
pub type __off_t = libc::c_long;
pub const _SVID_: _LIB_VERSION_TYPE = 0;
pub type __uint16_t = libc::c_ushort;
pub type __off64_t = libc::c_long;
pub const _IEEE_: _LIB_VERSION_TYPE = -1;
pub type SpeexResamplerState = SpeexResamplerState_;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct SpeexResamplerState_ {
    pub in_rate: spx_uint32_t,
    pub out_rate: spx_uint32_t,
    pub num_rate: spx_uint32_t,
    pub den_rate: spx_uint32_t,
    pub quality: libc::c_int,
    pub nb_channels: spx_uint32_t,
    pub filt_len: spx_uint32_t,
    pub mem_alloc_size: spx_uint32_t,
    pub buffer_size: spx_uint32_t,
    pub int_advance: libc::c_int,
    pub frac_advance: libc::c_int,
    pub cutoff: libc::c_float,
    pub oversample: spx_uint32_t,
    pub initialised: libc::c_int,
    pub started: libc::c_int,
    pub last_sample: *mut spx_int32_t,
    pub samp_frac_num: *mut spx_uint32_t,
    pub magic_samples: *mut spx_uint32_t,
    pub mem: *mut spx_word16_t,
    pub sinc_table: *mut spx_word16_t,
    pub sinc_table_length: spx_uint32_t,
    pub resampler_ptr: resampler_basic_func,
    pub in_stride: libc::c_int,
    pub out_stride: libc::c_int,
}
pub type size_t = libc::c_ulong;
pub type spx_int32_t = int32_t;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _IO_marker {
    pub _next: *mut _IO_marker,
    pub _sbuf: *mut _IO_FILE_0,
    pub _pos: libc::c_int,
}
pub type __dev_t = libc::c_ulong;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct FuncDef {
    pub table: *const libc::c_double,
    pub oversample: libc::c_int,
}
pub type int32_t = libc::c_int;
pub type spx_int16_t = int16_t;
pub const _POSIX_: _LIB_VERSION_TYPE = 2;
pub const _ISOC_: _LIB_VERSION_TYPE = 3;
pub type spx_word16_t = libc::c_float;
pub type _IO_FILE = _IO_FILE_0;
pub type spx_word32_t = libc::c_float;
pub const _XOPEN_: _LIB_VERSION_TYPE = 1;
pub type int16_t = libc::c_short;
pub type __uint32_t = libc::c_uint;
pub type unnamed = libc::c_uint;
pub const RESAMPLER_ERR_SUCCESS: unnamed = 0;
pub type resampler_basic_func = Option<
    unsafe extern "C" fn(
        _: *mut SpeexResamplerState,
        _: spx_uint32_t,
        _: *const spx_word16_t,
        _: *mut spx_uint32_t,
        _: *mut spx_word16_t,
        _: *mut spx_uint32_t,
    ) -> libc::c_int,
>;
pub type _IO_lock_t = ();
pub const RESAMPLER_ERR_BAD_STATE: unnamed = 2;
pub const RESAMPLER_ERR_MAX_ERROR: unnamed = 6;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct QualityMapping {
    pub base_length: libc::c_int,
    pub oversample: libc::c_int,
    pub downsample_bandwidth: libc::c_float,
    pub upsample_bandwidth: libc::c_float,
    pub window_func: *const FuncDef,
}
pub type FILE = _IO_FILE_0;
#[derive(Copy, Clone)]
#[repr(C)]
pub union unnamed_0 {
    __l: libc::c_double,
    __i: [libc::c_int; 3],
}
pub type spx_uint32_t = uint32_t;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _IO_FILE_0 {
    pub _flags: libc::c_int,
    pub _IO_read_ptr: *mut libc::c_char,
    pub _IO_read_end: *mut libc::c_char,
    pub _IO_read_base: *mut libc::c_char,
    pub _IO_write_base: *mut libc::c_char,
    pub _IO_write_ptr: *mut libc::c_char,
    pub _IO_write_end: *mut libc::c_char,
    pub _IO_buf_base: *mut libc::c_char,
    pub _IO_buf_end: *mut libc::c_char,
    pub _IO_save_base: *mut libc::c_char,
    pub _IO_backup_base: *mut libc::c_char,
    pub _IO_save_end: *mut libc::c_char,
    pub _markers: *mut _IO_marker,
    pub _chain: *mut _IO_FILE_0,
    pub _fileno: libc::c_int,
    pub _flags2: libc::c_int,
    pub _old_offset: __off_t,
    pub _cur_column: libc::c_ushort,
    pub _vtable_offset: libc::c_schar,
    pub _shortbuf: [libc::c_char; 1],
    pub _lock: *mut libc::c_void,
    pub _offset: __off64_t,
    pub __pad1: *mut libc::c_void,
    pub __pad2: *mut libc::c_void,
    pub __pad3: *mut libc::c_void,
    pub __pad4: *mut libc::c_void,
    pub __pad5: size_t,
    pub _mode: libc::c_int,
    pub _unused2: [libc::c_char; 20],
}
pub const RESAMPLER_ERR_ALLOC_FAILED: unnamed = 1;
pub type __compar_fn_t = Option<
    unsafe extern "C" fn(
        _: *const libc::c_void,
        _: *const libc::c_void,
    ) -> libc::c_int,
>;
/* * Create a new resampler with integer input and output rates.
 * @param nb_channels Number of channels to be processed
 * @param in_rate Input sampling rate (integer number of Hz).
 * @param out_rate Output sampling rate (integer number of Hz).
 * @param quality Resampling quality between 0 and 10, where 0 has poor quality
 * and 10 has very high quality.
 * @return Newly created resampler state
 * @retval NULL Error: not enough memory
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_init(
    mut nb_channels: spx_uint32_t,
    mut in_rate: spx_uint32_t,
    mut out_rate: spx_uint32_t,
    mut quality: libc::c_int,
    mut err: *mut libc::c_int,
) -> *mut SpeexResamplerState {
    return speex_resampler_init_frac(
        nb_channels,
        in_rate,
        out_rate,
        in_rate,
        out_rate,
        quality,
        err,
    );
}
/* * Create a new resampler with fractional input/output rates. The sampling
 * rate ratio is an arbitrary rational number with both the numerator and
 * denominator being 32-bit integers.
 * @param nb_channels Number of channels to be processed
 * @param ratio_num Numerator of the sampling rate ratio
 * @param ratio_den Denominator of the sampling rate ratio
 * @param in_rate Input sampling rate rounded to the nearest integer (in Hz).
 * @param out_rate Output sampling rate rounded to the nearest integer (in Hz).
 * @param quality Resampling quality between 0 and 10, where 0 has poor quality
 * and 10 has very high quality.
 * @return Newly created resampler state
 * @retval NULL Error: not enough memory
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_init_frac(
    mut nb_channels: spx_uint32_t,
    mut ratio_num: spx_uint32_t,
    mut ratio_den: spx_uint32_t,
    mut in_rate: spx_uint32_t,
    mut out_rate: spx_uint32_t,
    mut quality: libc::c_int,
    mut err: *mut libc::c_int,
) -> *mut SpeexResamplerState {
    if nb_channels == 0i32 as libc::c_uint
        || ratio_num == 0i32 as libc::c_uint
        || ratio_den == 0i32 as libc::c_uint
        || quality > 10i32
        || quality < 0i32
    {
        if !err.is_null() {
            *err = RESAMPLER_ERR_INVALID_ARG as libc::c_int
        }
        return 0 as *mut SpeexResamplerState;
    } else {
        let mut st = speex_alloc(::std::mem::size_of::<SpeexResamplerState>()
            as libc::c_ulong as libc::c_int)
            as *mut SpeexResamplerState;
        if st.is_null() {
            if !err.is_null() {
                *err = RESAMPLER_ERR_ALLOC_FAILED as libc::c_int
            }
            return 0 as *mut SpeexResamplerState;
        } else {
            (*st).initialised = 0i32;
            (*st).started = 0i32;
            (*st).in_rate = 0i32 as spx_uint32_t;
            (*st).out_rate = 0i32 as spx_uint32_t;
            (*st).num_rate = 0i32 as spx_uint32_t;
            (*st).den_rate = 0i32 as spx_uint32_t;
            (*st).quality = -1i32;
            (*st).sinc_table_length = 0i32 as spx_uint32_t;
            (*st).mem_alloc_size = 0i32 as spx_uint32_t;
            (*st).filt_len = 0i32 as spx_uint32_t;
            (*st).mem = 0 as *mut spx_word16_t;
            (*st).resampler_ptr = None;
            (*st).cutoff = 1.0f32;
            (*st).nb_channels = nb_channels;
            (*st).in_stride = 1i32;
            (*st).out_stride = 1i32;
            (*st).buffer_size = 160i32 as spx_uint32_t;
            (*st).last_sample =
                speex_alloc((nb_channels as libc::c_ulong).wrapping_mul(
                    ::std::mem::size_of::<spx_int32_t>() as libc::c_ulong,
                ) as libc::c_int) as *mut spx_int32_t;
            if !(*st).last_sample.is_null() {
                (*st).magic_samples = speex_alloc(
                    (nb_channels as libc::c_ulong)
                        .wrapping_mul(::std::mem::size_of::<spx_uint32_t>()
                            as libc::c_ulong)
                        as libc::c_int,
                ) as *mut spx_uint32_t;
                if !(*st).magic_samples.is_null() {
                    (*st).samp_frac_num = speex_alloc(
                        (nb_channels as libc::c_ulong).wrapping_mul(
                            ::std::mem::size_of::<spx_uint32_t>()
                                as libc::c_ulong,
                        ) as libc::c_int,
                    )
                        as *mut spx_uint32_t;
                    if !(*st).samp_frac_num.is_null() {
                        speex_resampler_set_quality(st, quality);
                        speex_resampler_set_rate_frac(
                            st, ratio_num, ratio_den, in_rate, out_rate,
                        );
                        let filter_err = update_filter(st);
                        if filter_err == RESAMPLER_ERR_SUCCESS as libc::c_int {
                            (*st).initialised = 1i32
                        } else {
                            speex_resampler_destroy(st);
                            st = 0 as *mut SpeexResamplerState
                        }
                        if !err.is_null() {
                            *err = filter_err
                        }
                        return st;
                    }
                }
            }
            if !err.is_null() {
                *err = RESAMPLER_ERR_ALLOC_FAILED as libc::c_int
            }
            speex_resampler_destroy(st);
            return 0 as *mut SpeexResamplerState;
        }
    };
}
/* * Destroy a resampler state.
 * @param st Resampler state
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_destroy(
    mut st: *mut SpeexResamplerState,
) -> () {
    speex_free((*st).mem as *mut libc::c_void);
    speex_free((*st).sinc_table as *mut libc::c_void);
    speex_free((*st).last_sample as *mut libc::c_void);
    speex_free((*st).magic_samples as *mut libc::c_void);
    speex_free((*st).samp_frac_num as *mut libc::c_void);
    speex_free(st as *mut libc::c_void);
}
/* * Speex wrapper for calloc. To do your own dynamic allocation, all you need to do is replace this function, speex_realloc and speex_alloc */
unsafe extern "C" fn speex_free(mut ptr: *mut libc::c_void) -> () {
    free(ptr);
}
unsafe extern "C" fn update_filter(
    mut st: *mut SpeexResamplerState,
) -> libc::c_int {
    let mut current_block: u64;
    let mut old_length: spx_uint32_t = (*st).filt_len;
    let mut old_alloc_size: spx_uint32_t = (*st).mem_alloc_size;
    let mut use_direct: libc::c_int = 0;
    let mut min_sinc_table_length: spx_uint32_t = 0;
    let mut min_alloc_size: spx_uint32_t = 0;
    (*st).int_advance =
        (*st).num_rate.wrapping_div((*st).den_rate) as libc::c_int;
    (*st).frac_advance =
        (*st).num_rate.wrapping_rem((*st).den_rate) as libc::c_int;
    (*st).oversample =
        quality_map[(*st).quality as usize].oversample as spx_uint32_t;
    (*st).filt_len =
        quality_map[(*st).quality as usize].base_length as spx_uint32_t;
    if (*st).num_rate > (*st).den_rate {
        (*st).cutoff = quality_map[(*st).quality as usize]
            .downsample_bandwidth
            * (*st).den_rate as libc::c_float
            / (*st).num_rate as libc::c_float;
        if _muldiv(
            &mut (*st).filt_len as *mut spx_uint32_t,
            (*st).filt_len,
            (*st).num_rate,
            (*st).den_rate,
        ) != RESAMPLER_ERR_SUCCESS as libc::c_int
        {
            current_block = 13465332693182510667;
        } else {
            (*st).filt_len =
                ((*st).filt_len.wrapping_sub(1i32 as libc::c_uint)
                    & !7i32 as libc::c_uint)
                    .wrapping_add(8i32 as libc::c_uint);
            if (2i32 as libc::c_uint).wrapping_mul((*st).den_rate)
                < (*st).num_rate
            {
                (*st).oversample >>= 1i32
            }
            if (4i32 as libc::c_uint).wrapping_mul((*st).den_rate)
                < (*st).num_rate
            {
                (*st).oversample >>= 1i32
            }
            if (8i32 as libc::c_uint).wrapping_mul((*st).den_rate)
                < (*st).num_rate
            {
                (*st).oversample >>= 1i32
            }
            if (16i32 as libc::c_uint).wrapping_mul((*st).den_rate)
                < (*st).num_rate
            {
                (*st).oversample >>= 1i32
            }
            if (*st).oversample < 1i32 as libc::c_uint {
                (*st).oversample = 1i32 as spx_uint32_t;
                current_block = 8258075665625361029;
            } else {
                current_block = 8258075665625361029;
            }
        }
    } else {
        (*st).cutoff = quality_map[(*st).quality as usize].upsample_bandwidth;
        current_block = 8258075665625361029;
    }
    match current_block {
        8258075665625361029 => {
            use_direct = ((*st).filt_len.wrapping_mul((*st).den_rate)
                <= (*st)
                    .filt_len
                    .wrapping_mul((*st).oversample)
                    .wrapping_add(8i32 as libc::c_uint)
                && (2147483647i32 as libc::c_ulong)
                    .wrapping_div(
                        ::std::mem::size_of::<spx_word16_t>() as libc::c_ulong
                    )
                    .wrapping_div((*st).den_rate as libc::c_ulong)
                    >= (*st).filt_len as libc::c_ulong)
                as libc::c_int;
            if 0 != use_direct {
                min_sinc_table_length =
                    (*st).filt_len.wrapping_mul((*st).den_rate);
                current_block = 14523784380283086299;
            } else if (2147483647i32 as libc::c_ulong)
                .wrapping_div(
                    ::std::mem::size_of::<spx_word16_t>() as libc::c_ulong
                )
                .wrapping_sub(8i32 as libc::c_ulong)
                .wrapping_div((*st).oversample as libc::c_ulong)
                < (*st).filt_len as libc::c_ulong
            {
                current_block = 13465332693182510667;
            } else {
                min_sinc_table_length = (*st)
                    .filt_len
                    .wrapping_mul((*st).oversample)
                    .wrapping_add(8i32 as libc::c_uint);
                current_block = 14523784380283086299;
            }
            match current_block {
                13465332693182510667 => {}
                _ => {
                    if (*st).sinc_table_length < min_sinc_table_length {
                        let mut sinc_table: *mut spx_word16_t = speex_realloc(
                            (*st).sinc_table as *mut libc::c_void,
                            (min_sinc_table_length as libc::c_ulong)
                                .wrapping_mul(
                                    ::std::mem::size_of::<spx_word16_t>()
                                        as libc::c_ulong,
                                ) as libc::c_int,
                        )
                            as *mut spx_word16_t;
                        if sinc_table.is_null() {
                            current_block = 13465332693182510667;
                        } else {
                            (*st).sinc_table = sinc_table;
                            (*st).sinc_table_length = min_sinc_table_length;
                            current_block = 11650488183268122163;
                        }
                    } else {
                        current_block = 11650488183268122163;
                    }
                    match current_block {
                        13465332693182510667 => {}
                        _ => {
                            if 0 != use_direct {
                                let mut i: spx_uint32_t = 0;
                                i = 0i32 as spx_uint32_t;
                                while i < (*st).den_rate {
                                    let mut j: spx_int32_t = 0;
                                    j = 0i32;
                                    while (j as libc::c_uint) < (*st).filt_len
                                    {
                                        *(*st).sinc_table.offset(
                                            i.wrapping_mul((*st).filt_len)
                                                .wrapping_add(
                                                    j as libc::c_uint,
                                                )
                                                as isize,
                                        ) = sinc(
                                            (*st).cutoff,
                                            (j - (*st).filt_len as spx_int32_t
                                                / 2i32
                                                + 1i32)
                                                as libc::c_float
                                                - i as libc::c_float
                                                    / (*st).den_rate
                                                        as libc::c_float,
                                            (*st).filt_len as libc::c_int,
                                            quality_map
                                                [(*st).quality as usize]
                                                .window_func,
                                        );
                                        j += 1
                                    }
                                    i = i.wrapping_add(1)
                                }
                                if (*st).quality > 8i32 {
                                    (*st).resampler_ptr =
                                        Some(resampler_basic_direct_double)
                                } else {
                                    (*st).resampler_ptr =
                                        Some(resampler_basic_direct_single)
                                }
                            } else {
                                let mut i_0: spx_int32_t = 0;
                                i_0 = -4i32;
                                while i_0
                                    < (*st)
                                        .oversample
                                        .wrapping_mul((*st).filt_len)
                                        .wrapping_add(4i32 as libc::c_uint)
                                        as spx_int32_t
                                {
                                    *(*st)
                                        .sinc_table
                                        .offset((i_0 + 4i32) as isize) = sinc(
                                        (*st).cutoff,
                                        i_0 as libc::c_float
                                            / (*st).oversample
                                                as libc::c_float
                                            - (*st).filt_len.wrapping_div(
                                                2i32 as libc::c_uint,
                                            )
                                                as libc::c_float,
                                        (*st).filt_len as libc::c_int,
                                        quality_map[(*st).quality as usize]
                                            .window_func,
                                    );
                                    i_0 += 1
                                }
                                if (*st).quality > 8i32 {
                                    (*st).resampler_ptr = Some(
                                        resampler_basic_interpolate_double,
                                    )
                                } else {
                                    (*st).resampler_ptr = Some(
                                        resampler_basic_interpolate_single,
                                    )
                                }
                            }
                            min_alloc_size = (*st)
                                .filt_len
                                .wrapping_sub(1i32 as libc::c_uint)
                                .wrapping_add((*st).buffer_size);
                            if min_alloc_size > (*st).mem_alloc_size {
                                let mut mem: *mut spx_word16_t =
                                    0 as *mut spx_word16_t;
                                if (2147483647i32 as libc::c_ulong)
                                    .wrapping_div(::std::mem::size_of::<
                                        spx_word16_t,
                                    >(
                                    )
                                        as libc::c_ulong)
                                    .wrapping_div(
                                        (*st).nb_channels as libc::c_ulong,
                                    )
                                    < min_alloc_size as libc::c_ulong
                                {
                                    current_block = 13465332693182510667;
                                } else {
                                    mem = speex_realloc(
                                        (*st).mem as *mut libc::c_void,
                                        ((*st)
                                            .nb_channels
                                            .wrapping_mul(min_alloc_size)
                                            as libc::c_ulong)
                                            .wrapping_mul(
                                                ::std::mem::size_of::<
                                                    spx_word16_t,
                                                >(
                                                )
                                                    as libc::c_ulong,
                                            )
                                            as libc::c_int,
                                    )
                                        as *mut spx_word16_t;
                                    if mem.is_null() {
                                        current_block = 13465332693182510667;
                                    } else {
                                        (*st).mem = mem;
                                        (*st).mem_alloc_size = min_alloc_size;
                                        current_block = 15089075282327824602;
                                    }
                                }
                            } else {
                                current_block = 15089075282327824602;
                            }
                            match current_block {
                                13465332693182510667 => {}
                                _ => {
                                    if 0 == (*st).started {
                                        let mut i_1: spx_uint32_t = 0;
                                        i_1 = 0i32 as spx_uint32_t;
                                        while i_1
                                            < (*st).nb_channels.wrapping_mul(
                                                (*st).mem_alloc_size,
                                            )
                                        {
                                            *(*st).mem.offset(i_1 as isize) =
                                                0i32 as spx_word16_t;
                                            i_1 = i_1.wrapping_add(1)
                                        }
                                    } else if (*st).filt_len > old_length {
                                        let mut i_2: spx_uint32_t = 0;
                                        i_2 = (*st).nb_channels;
                                        loop {
                                            let fresh0 = i_2;
                                            i_2 = i_2.wrapping_sub(1);
                                            if !(0 != fresh0) {
                                                break;
                                            }
                                            let mut j_0: spx_uint32_t = 0;
                                            let mut olen: spx_uint32_t =
                                                old_length;
                                            olen = old_length.wrapping_add(
                                                (2i32 as libc::c_uint)
                                                    .wrapping_mul(
                                                        *(*st)
                                                            .magic_samples
                                                            .offset(
                                                                i_2 as isize,
                                                            ),
                                                    ),
                                            );
                                            j_0 = old_length
                                                .wrapping_sub(
                                                    1i32 as libc::c_uint,
                                                )
                                                .wrapping_add(
                                                    *(*st)
                                                        .magic_samples
                                                        .offset(i_2 as isize),
                                                );
                                            loop {
                                                let fresh1 = j_0;
                                                j_0 = j_0.wrapping_sub(1);
                                                if !(0 != fresh1) {
                                                    break;
                                                }
                                                *(*st).mem.offset(
                                                    i_2.wrapping_mul(
                                                        (*st).mem_alloc_size,
                                                    )
                                                    .wrapping_add(j_0)
                                                    .wrapping_add(
                                                        *(*st)
                                                            .magic_samples
                                                            .offset(
                                                                i_2 as isize,
                                                            ),
                                                    )
                                                        as isize,
                                                ) = *(*st).mem.offset(
                                                    i_2.wrapping_mul(
                                                        old_alloc_size,
                                                    )
                                                    .wrapping_add(j_0)
                                                        as isize,
                                                )
                                            }
                                            j_0 = 0i32 as spx_uint32_t;
                                            while j_0
                                                < *(*st)
                                                    .magic_samples
                                                    .offset(i_2 as isize)
                                            {
                                                *(*st).mem.offset(
                                                    i_2.wrapping_mul(
                                                        (*st).mem_alloc_size,
                                                    )
                                                    .wrapping_add(j_0)
                                                        as isize,
                                                ) = 0i32 as spx_word16_t;
                                                j_0 = j_0.wrapping_add(1)
                                            }
                                            *(*st)
                                                .magic_samples
                                                .offset(i_2 as isize) =
                                                0i32 as spx_uint32_t;
                                            if (*st).filt_len > olen {
                                                j_0 = 0i32 as spx_uint32_t;
                                                while j_0
                                                    < olen.wrapping_sub(
                                                        1i32 as libc::c_uint,
                                                    )
                                                {
                                                    *(*st).mem.offset(i_2.wrapping_mul((*st).mem_alloc_size).wrapping_add((*st).filt_len.wrapping_sub(2i32
                                                                                                                                                          as
                                                                                                                                                          libc::c_uint).wrapping_sub(j_0))
                                                                          as
                                                                          isize)
                                                        =
                                                        *(*st).mem.offset(i_2.wrapping_mul((*st).mem_alloc_size).wrapping_add(olen.wrapping_sub(2i32
                                                                                                                                                    as
                                                                                                                                                    libc::c_uint).wrapping_sub(j_0))
                                                                              as
                                                                              isize);
                                                    j_0 = j_0.wrapping_add(1)
                                                }
                                                while j_0
                                                    < (*st)
                                                        .filt_len
                                                        .wrapping_sub(
                                                        1i32 as libc::c_uint,
                                                    )
                                                {
                                                    *(*st).mem.offset(i_2.wrapping_mul((*st).mem_alloc_size).wrapping_add((*st).filt_len.wrapping_sub(2i32
                                                                                                                                                          as
                                                                                                                                                          libc::c_uint).wrapping_sub(j_0))
                                                                          as
                                                                          isize)
                                                        =
                                                        0i32 as spx_word16_t;
                                                    j_0 = j_0.wrapping_add(1)
                                                }
                                                let ref mut fresh2 = *(*st)
                                                    .last_sample
                                                    .offset(i_2 as isize);
                                                *fresh2 =
                                                    (*fresh2 as
                                                         libc::c_uint).wrapping_add((*st).filt_len.wrapping_sub(olen).wrapping_div(2i32
                                                                                                                                       as
                                                                                                                                       libc::c_uint))
                                                        as spx_int32_t as
                                                        spx_int32_t
                                            } else {
                                                *(*st)
                                                    .magic_samples
                                                    .offset(i_2 as isize) =
                                                    olen.wrapping_sub(
                                                        (*st).filt_len,
                                                    )
                                                    .wrapping_div(
                                                        2i32 as libc::c_uint,
                                                    );
                                                j_0 = 0i32 as spx_uint32_t;
                                                while j_0 < (*st)
                                                    .filt_len
                                                    .wrapping_sub(
                                                        1i32 as libc::c_uint,
                                                    )
                                                    .wrapping_add(
                                                        *(*st)
                                                            .magic_samples
                                                            .offset(
                                                                i_2 as isize,
                                                            ),
                                                    )
                                                {
                                                    *(*st).mem.offset(i_2.wrapping_mul((*st).mem_alloc_size).wrapping_add(j_0)
                                                                          as
                                                                          isize)
                                                        =
                                                        *(*st).mem.offset(i_2.wrapping_mul((*st).mem_alloc_size).wrapping_add(j_0).wrapping_add(*(*st).magic_samples.offset(i_2
                                                                                                                                                                                as
                                                                                                                                                                                isize))
                                                                              as
                                                                              isize);
                                                    j_0 = j_0.wrapping_add(1)
                                                }
                                            }
                                        }
                                    } else if (*st).filt_len < old_length {
                                        let mut i_3: spx_uint32_t = 0;
                                        i_3 = 0i32 as spx_uint32_t;
                                        while i_3 < (*st).nb_channels {
                                            let mut j_1: spx_uint32_t = 0;
                                            let mut old_magic: spx_uint32_t =
                                                *(*st)
                                                    .magic_samples
                                                    .offset(i_3 as isize);
                                            *(*st)
                                                .magic_samples
                                                .offset(i_3 as isize) =
                                                old_length
                                                    .wrapping_sub(
                                                        (*st).filt_len,
                                                    )
                                                    .wrapping_div(
                                                        2i32 as libc::c_uint,
                                                    );
                                            j_1 = 0i32 as spx_uint32_t;
                                            while j_1
                                                < (*st)
                                                    .filt_len
                                                    .wrapping_sub(
                                                        1i32 as libc::c_uint,
                                                    )
                                                    .wrapping_add(
                                                        *(*st)
                                                            .magic_samples
                                                            .offset(
                                                                i_3 as isize,
                                                            ),
                                                    )
                                                    .wrapping_add(old_magic)
                                            {
                                                *(*st).mem.offset(
                                                    i_3.wrapping_mul(
                                                        (*st).mem_alloc_size,
                                                    )
                                                    .wrapping_add(j_1)
                                                        as isize,
                                                ) = *(*st).mem.offset(
                                                    i_3.wrapping_mul(
                                                        (*st).mem_alloc_size,
                                                    )
                                                    .wrapping_add(j_1)
                                                    .wrapping_add(
                                                        *(*st)
                                                            .magic_samples
                                                            .offset(
                                                                i_3 as isize,
                                                            ),
                                                    )
                                                        as isize,
                                                );
                                                j_1 = j_1.wrapping_add(1)
                                            }
                                            let ref mut fresh3 = *(*st)
                                                .magic_samples
                                                .offset(i_3 as isize);
                                            *fresh3 = (*fresh3 as libc::c_uint)
                                                .wrapping_add(old_magic)
                                                as spx_uint32_t
                                                as spx_uint32_t;
                                            i_3 = i_3.wrapping_add(1)
                                        }
                                    }
                                    return RESAMPLER_ERR_SUCCESS
                                        as libc::c_int;
                                }
                            }
                        }
                    }
                }
            }
        }
        _ => {}
    }
    (*st).resampler_ptr = Some(resampler_basic_zero);
    (*st).filt_len = old_length;
    return RESAMPLER_ERR_ALLOC_FAILED as libc::c_int;
}
unsafe extern "C" fn resampler_basic_zero(
    mut st: *mut SpeexResamplerState,
    mut channel_index: spx_uint32_t,
    mut in_0: *const spx_word16_t,
    mut in_len: *mut spx_uint32_t,
    mut out: *mut spx_word16_t,
    mut out_len: *mut spx_uint32_t,
) -> libc::c_int {
    let mut out_sample: libc::c_int = 0i32;
    let mut last_sample: libc::c_int =
        *(*st).last_sample.offset(channel_index as isize);
    let mut samp_frac_num: spx_uint32_t =
        *(*st).samp_frac_num.offset(channel_index as isize);
    let out_stride: libc::c_int = (*st).out_stride;
    let int_advance: libc::c_int = (*st).int_advance;
    let frac_advance: libc::c_int = (*st).frac_advance;
    let den_rate: spx_uint32_t = (*st).den_rate;
    while !(last_sample >= *in_len as spx_int32_t
        || out_sample >= *out_len as spx_int32_t)
    {
        let fresh4 = out_sample;
        out_sample = out_sample + 1;
        *out.offset((out_stride * fresh4) as isize) = 0i32 as spx_word16_t;
        last_sample += int_advance;
        samp_frac_num = (samp_frac_num as libc::c_uint)
            .wrapping_add(frac_advance as libc::c_uint)
            as spx_uint32_t as spx_uint32_t;
        if !(samp_frac_num >= den_rate) {
            continue;
        }
        samp_frac_num = (samp_frac_num as libc::c_uint).wrapping_sub(den_rate)
            as spx_uint32_t as spx_uint32_t;
        last_sample += 1
    }
    *(*st).last_sample.offset(channel_index as isize) = last_sample;
    *(*st).samp_frac_num.offset(channel_index as isize) = samp_frac_num;
    return out_sample;
}
/* * Speex wrapper for realloc. To do your own dynamic allocation, all you need to do is replace this function, speex_alloc and speex_free */
unsafe extern "C" fn speex_realloc(
    mut ptr: *mut libc::c_void,
    mut size: libc::c_int,
) -> *mut libc::c_void {
    return realloc(ptr, size as usize);
}
unsafe extern "C" fn resampler_basic_interpolate_single(
    mut st: *mut SpeexResamplerState,
    mut channel_index: spx_uint32_t,
    mut in_0: *const spx_word16_t,
    mut in_len: *mut spx_uint32_t,
    mut out: *mut spx_word16_t,
    mut out_len: *mut spx_uint32_t,
) -> libc::c_int {
    let N: libc::c_int = (*st).filt_len as libc::c_int;
    let mut out_sample: libc::c_int = 0i32;
    let mut last_sample: libc::c_int =
        *(*st).last_sample.offset(channel_index as isize);
    let mut samp_frac_num: spx_uint32_t =
        *(*st).samp_frac_num.offset(channel_index as isize);
    let out_stride: libc::c_int = (*st).out_stride;
    let int_advance: libc::c_int = (*st).int_advance;
    let frac_advance: libc::c_int = (*st).frac_advance;
    let den_rate: spx_uint32_t = (*st).den_rate;
    let mut sum: spx_word32_t = 0.;
    while !(last_sample >= *in_len as spx_int32_t
        || out_sample >= *out_len as spx_int32_t)
    {
        let mut iptr: *const spx_word16_t =
            &*in_0.offset(last_sample as isize) as *const spx_word16_t;
        let offset: libc::c_int = samp_frac_num
            .wrapping_mul((*st).oversample)
            .wrapping_div((*st).den_rate)
            as libc::c_int;
        let frac: spx_word16_t = samp_frac_num
            .wrapping_mul((*st).oversample)
            .wrapping_rem((*st).den_rate)
            as libc::c_float
            / (*st).den_rate as libc::c_float;
        let mut interp: [spx_word16_t; 4] = [0.; 4];
        let mut j: libc::c_int = 0;
        let mut accum: [spx_word32_t; 4] = [
            0i32 as spx_word32_t,
            0i32 as spx_word32_t,
            0i32 as spx_word32_t,
            0i32 as spx_word32_t,
        ];
        j = 0i32;
        while j < N {
            let curr_in: spx_word16_t = *iptr.offset(j as isize);
            accum[0usize] += curr_in
                * *(*st).sinc_table.offset(
                    (4i32 as libc::c_uint)
                        .wrapping_add(
                            ((j + 1i32) as libc::c_uint)
                                .wrapping_mul((*st).oversample),
                        )
                        .wrapping_sub(offset as libc::c_uint)
                        .wrapping_sub(2i32 as libc::c_uint)
                        as isize,
                );
            accum[1usize] += curr_in
                * *(*st).sinc_table.offset(
                    (4i32 as libc::c_uint)
                        .wrapping_add(
                            ((j + 1i32) as libc::c_uint)
                                .wrapping_mul((*st).oversample),
                        )
                        .wrapping_sub(offset as libc::c_uint)
                        .wrapping_sub(1i32 as libc::c_uint)
                        as isize,
                );
            accum[2usize] += curr_in
                * *(*st).sinc_table.offset(
                    (4i32 as libc::c_uint)
                        .wrapping_add(
                            ((j + 1i32) as libc::c_uint)
                                .wrapping_mul((*st).oversample),
                        )
                        .wrapping_sub(offset as libc::c_uint)
                        as isize,
                );
            accum[3usize] += curr_in
                * *(*st).sinc_table.offset(
                    (4i32 as libc::c_uint)
                        .wrapping_add(
                            ((j + 1i32) as libc::c_uint)
                                .wrapping_mul((*st).oversample),
                        )
                        .wrapping_sub(offset as libc::c_uint)
                        .wrapping_add(1i32 as libc::c_uint)
                        as isize,
                );
            j += 1
        }
        cubic_coef(frac, interp.as_mut_ptr());
        sum = interp[0usize] * accum[0usize]
            + interp[1usize] * accum[1usize]
            + interp[2usize] * accum[2usize]
            + interp[3usize] * accum[3usize];
        sum = sum;
        let fresh5 = out_sample;
        out_sample = out_sample + 1;
        *out.offset((out_stride * fresh5) as isize) = sum;
        last_sample += int_advance;
        samp_frac_num = (samp_frac_num as libc::c_uint)
            .wrapping_add(frac_advance as libc::c_uint)
            as spx_uint32_t as spx_uint32_t;
        if !(samp_frac_num >= den_rate) {
            continue;
        }
        samp_frac_num = (samp_frac_num as libc::c_uint).wrapping_sub(den_rate)
            as spx_uint32_t as spx_uint32_t;
        last_sample += 1
    }
    *(*st).last_sample.offset(channel_index as isize) = last_sample;
    *(*st).samp_frac_num.offset(channel_index as isize) = samp_frac_num;
    return out_sample;
}
unsafe extern "C" fn cubic_coef(
    mut frac: spx_word16_t,
    mut interp: *mut spx_word16_t,
) -> () {
    *interp.offset(0isize) = -0.16666999459266663f32 * frac
        + 0.16666999459266663f32 * frac * frac * frac;
    *interp.offset(1isize) =
        frac + 0.5f32 * frac * frac - 0.5f32 * frac * frac * frac;
    *interp.offset(3isize) = -0.3333300054073334f32 * frac
        + 0.5f32 * frac * frac
        - 0.16666999459266663f32 * frac * frac * frac;
    *interp.offset(2isize) = (1.0f64
        - *interp.offset(0isize) as libc::c_double
        - *interp.offset(1isize) as libc::c_double
        - *interp.offset(3isize) as libc::c_double)
        as spx_word16_t;
}
unsafe extern "C" fn resampler_basic_interpolate_double(
    mut st: *mut SpeexResamplerState,
    mut channel_index: spx_uint32_t,
    mut in_0: *const spx_word16_t,
    mut in_len: *mut spx_uint32_t,
    mut out: *mut spx_word16_t,
    mut out_len: *mut spx_uint32_t,
) -> libc::c_int {
    let N: libc::c_int = (*st).filt_len as libc::c_int;
    let mut out_sample: libc::c_int = 0i32;
    let mut last_sample: libc::c_int =
        *(*st).last_sample.offset(channel_index as isize);
    let mut samp_frac_num: spx_uint32_t =
        *(*st).samp_frac_num.offset(channel_index as isize);
    let out_stride: libc::c_int = (*st).out_stride;
    let int_advance: libc::c_int = (*st).int_advance;
    let frac_advance: libc::c_int = (*st).frac_advance;
    let den_rate: spx_uint32_t = (*st).den_rate;
    let mut sum: spx_word32_t = 0.;
    while !(last_sample >= *in_len as spx_int32_t
        || out_sample >= *out_len as spx_int32_t)
    {
        let mut iptr: *const spx_word16_t =
            &*in_0.offset(last_sample as isize) as *const spx_word16_t;
        let offset: libc::c_int = samp_frac_num
            .wrapping_mul((*st).oversample)
            .wrapping_div((*st).den_rate)
            as libc::c_int;
        let frac: spx_word16_t = samp_frac_num
            .wrapping_mul((*st).oversample)
            .wrapping_rem((*st).den_rate)
            as libc::c_float
            / (*st).den_rate as libc::c_float;
        let mut interp: [spx_word16_t; 4] = [0.; 4];
        let mut j: libc::c_int = 0;
        let mut accum: [libc::c_double; 4] = [
            0i32 as libc::c_double,
            0i32 as libc::c_double,
            0i32 as libc::c_double,
            0i32 as libc::c_double,
        ];
        j = 0i32;
        while j < N {
            let curr_in: libc::c_double =
                *iptr.offset(j as isize) as libc::c_double;
            accum[0usize] += (curr_in as spx_word32_t
                * *(*st).sinc_table.offset(
                    (4i32 as libc::c_uint)
                        .wrapping_add(
                            ((j + 1i32) as libc::c_uint)
                                .wrapping_mul((*st).oversample),
                        )
                        .wrapping_sub(offset as libc::c_uint)
                        .wrapping_sub(2i32 as libc::c_uint)
                        as isize,
                )) as libc::c_double;
            accum[1usize] += (curr_in as spx_word32_t
                * *(*st).sinc_table.offset(
                    (4i32 as libc::c_uint)
                        .wrapping_add(
                            ((j + 1i32) as libc::c_uint)
                                .wrapping_mul((*st).oversample),
                        )
                        .wrapping_sub(offset as libc::c_uint)
                        .wrapping_sub(1i32 as libc::c_uint)
                        as isize,
                )) as libc::c_double;
            accum[2usize] += (curr_in as spx_word32_t
                * *(*st).sinc_table.offset(
                    (4i32 as libc::c_uint)
                        .wrapping_add(
                            ((j + 1i32) as libc::c_uint)
                                .wrapping_mul((*st).oversample),
                        )
                        .wrapping_sub(offset as libc::c_uint)
                        as isize,
                )) as libc::c_double;
            accum[3usize] += (curr_in as spx_word32_t
                * *(*st).sinc_table.offset(
                    (4i32 as libc::c_uint)
                        .wrapping_add(
                            ((j + 1i32) as libc::c_uint)
                                .wrapping_mul((*st).oversample),
                        )
                        .wrapping_sub(offset as libc::c_uint)
                        .wrapping_add(1i32 as libc::c_uint)
                        as isize,
                )) as libc::c_double;
            j += 1
        }
        cubic_coef(frac, interp.as_mut_ptr());
        sum = (interp[0usize] as libc::c_double * accum[0usize]
            + interp[1usize] as libc::c_double * accum[1usize]
            + interp[2usize] as libc::c_double * accum[2usize]
            + interp[3usize] as libc::c_double * accum[3usize])
            as spx_word32_t;
        let fresh6 = out_sample;
        out_sample = out_sample + 1;
        *out.offset((out_stride * fresh6) as isize) = sum;
        last_sample += int_advance;
        samp_frac_num = (samp_frac_num as libc::c_uint)
            .wrapping_add(frac_advance as libc::c_uint)
            as spx_uint32_t as spx_uint32_t;
        if !(samp_frac_num >= den_rate) {
            continue;
        }
        samp_frac_num = (samp_frac_num as libc::c_uint).wrapping_sub(den_rate)
            as spx_uint32_t as spx_uint32_t;
        last_sample += 1
    }
    *(*st).last_sample.offset(channel_index as isize) = last_sample;
    *(*st).samp_frac_num.offset(channel_index as isize) = samp_frac_num;
    return out_sample;
}
static mut quality_map: [QualityMapping; 11] = unsafe {
    [
        QualityMapping {
            base_length: 8i32,
            oversample: 4i32,
            downsample_bandwidth: 0.8299999833106995f32,
            upsample_bandwidth: 0.8600000143051148f32,
            window_func: &_KAISER6 as *const FuncDef,
        },
        QualityMapping {
            base_length: 16i32,
            oversample: 4i32,
            downsample_bandwidth: 0.8500000238418579f32,
            upsample_bandwidth: 0.8799999952316284f32,
            window_func: &_KAISER6 as *const FuncDef,
        },
        QualityMapping {
            base_length: 32i32,
            oversample: 4i32,
            downsample_bandwidth: 0.8820000290870667f32,
            upsample_bandwidth: 0.9100000262260437f32,
            window_func: &_KAISER6 as *const FuncDef,
        },
        QualityMapping {
            base_length: 48i32,
            oversample: 8i32,
            downsample_bandwidth: 0.8949999809265137f32,
            upsample_bandwidth: 0.9169999957084656f32,
            window_func: &_KAISER8 as *const FuncDef,
        },
        QualityMapping {
            base_length: 64i32,
            oversample: 8i32,
            downsample_bandwidth: 0.9210000038146973f32,
            upsample_bandwidth: 0.9399999976158142f32,
            window_func: &_KAISER8 as *const FuncDef,
        },
        QualityMapping {
            base_length: 80i32,
            oversample: 16i32,
            downsample_bandwidth: 0.921999990940094f32,
            upsample_bandwidth: 0.9399999976158142f32,
            window_func: &_KAISER10 as *const FuncDef,
        },
        QualityMapping {
            base_length: 96i32,
            oversample: 16i32,
            downsample_bandwidth: 0.9399999976158142f32,
            upsample_bandwidth: 0.9449999928474426f32,
            window_func: &_KAISER10 as *const FuncDef,
        },
        QualityMapping {
            base_length: 128i32,
            oversample: 16i32,
            downsample_bandwidth: 0.949999988079071f32,
            upsample_bandwidth: 0.949999988079071f32,
            window_func: &_KAISER10 as *const FuncDef,
        },
        QualityMapping {
            base_length: 160i32,
            oversample: 16i32,
            downsample_bandwidth: 0.9599999785423279f32,
            upsample_bandwidth: 0.9599999785423279f32,
            window_func: &_KAISER10 as *const FuncDef,
        },
        QualityMapping {
            base_length: 192i32,
            oversample: 32i32,
            downsample_bandwidth: 0.9679999947547913f32,
            upsample_bandwidth: 0.9679999947547913f32,
            window_func: &_KAISER12 as *const FuncDef,
        },
        QualityMapping {
            base_length: 256i32,
            oversample: 32i32,
            downsample_bandwidth: 0.9750000238418579f32,
            upsample_bandwidth: 0.9750000238418579f32,
            window_func: &_KAISER12 as *const FuncDef,
        },
    ]
};
static mut _KAISER12: FuncDef = unsafe {
    FuncDef {
        table: kaiser12_table.as_ptr(),
        oversample: 64i32,
    }
};
static mut kaiser12_table: [libc::c_double; 68] = {
    [
        0.99859849f64,
        1.0f64,
        0.99859849f64,
        0.99440475f64,
        0.98745105f64,
        0.97779076f64,
        0.9654977f64,
        0.95066529f64,
        0.93340547f64,
        0.91384741f64,
        0.89213598f64,
        0.86843014f64,
        0.84290116f64,
        0.81573067f64,
        0.78710866f64,
        0.75723148f64,
        0.7262997f64,
        0.69451601f64,
        0.66208321f64,
        0.62920216f64,
        0.59606986f64,
        0.56287762f64,
        0.52980938f64,
        0.49704014f64,
        0.46473455f64,
        0.43304576f64,
        0.40211431f64,
        0.37206735f64,
        0.343018f64,
        0.3150649f64,
        0.28829195f64,
        0.26276832f64,
        0.23854851f64,
        0.21567274f64,
        0.19416736f64,
        0.17404546f64,
        0.15530766f64,
        0.13794293999999999f64,
        0.12192957f64,
        0.10723616f64,
        0.09382272f64,
        0.08164178f64,
        0.0706395f64,
        0.06075685f64,
        0.05193064f64,
        0.04409466f64,
        0.03718069f64,
        0.03111947f64,
        0.02584161f64,
        0.02127838f64,
        0.0173625f64,
        0.01402878f64,
        0.01121463f64,
        0.00886058f64,
        0.00691064f64,
        0.00531256f64,
        0.00401805f64,
        0.00298291f64,
        0.00216702f64,
        0.00153438f64,
        0.00105297f64,
        0.00069463f64,
        0.00043489f64,
        0.00025272f64,
        0.00013031f64,
        0.0000527734f64,
        0.00001f64,
        0.0f64,
    ]
};
static mut _KAISER10: FuncDef = unsafe {
    FuncDef {
        table: kaiser10_table.as_ptr(),
        oversample: 32i32,
    }
};
static mut kaiser10_table: [libc::c_double; 36] = {
    [
        0.99537781f64,
        1.0f64,
        0.99537781f64,
        0.98162644f64,
        0.95908712f64,
        0.92831446f64,
        0.89005583f64,
        0.84522401f64,
        0.79486424f64,
        0.74011713f64,
        0.68217934f64,
        0.62226347f64,
        0.56155915f64,
        0.5011968f64,
        0.44221549f64,
        0.38553619f64,
        0.33194107f64,
        0.28205962f64,
        0.23636152f64,
        0.19515633f64,
        0.15859932f64,
        0.1267028f64,
        0.09935205f64,
        0.07632451f64,
        0.05731132f64,
        0.0419398f64,
        0.02979584f64,
        0.0204451f64,
        0.01345224f64,
        0.00839739f64,
        0.00488951f64,
        0.00257636f64,
        0.00115101f64,
        0.00035515f64,
        0.0f64,
        0.0f64,
    ]
};
static mut _KAISER8: FuncDef = unsafe {
    FuncDef {
        table: kaiser8_table.as_ptr(),
        oversample: 32i32,
    }
};
static mut kaiser8_table: [libc::c_double; 36] = {
    [
        0.99635258f64,
        1.0f64,
        0.99635258f64,
        0.98548012f64,
        0.96759014f64,
        0.943022f64,
        0.91223751f64,
        0.87580811f64,
        0.83439927f64,
        0.78875245f64,
        0.73966538f64,
        0.68797126f64,
        0.6345175f64,
        0.58014482f64,
        0.52566725f64,
        0.47185369f64,
        0.4194115f64,
        0.36897272f64,
        0.32108304f64,
        0.27619388f64,
        0.23465776f64,
        0.1967267f64,
        0.1625538f64,
        0.13219758f64,
        0.10562887f64,
        0.08273982f64,
        0.06335451f64,
        0.04724088f64,
        0.03412321f64,
        0.0236949f64,
        0.01563093f64,
        0.00959968f64,
        0.00527363f64,
        0.00233883f64,
        0.0005f64,
        0.0f64,
    ]
};
static mut _KAISER6: FuncDef = unsafe {
    FuncDef {
        table: kaiser6_table.as_ptr(),
        oversample: 32i32,
    }
};
static mut kaiser6_table: [libc::c_double; 36] = {
    [
        0.99733006f64,
        1.0f64,
        0.99733006f64,
        0.98935595f64,
        0.97618418f64,
        0.95799003f64,
        0.93501423f64,
        0.90755855f64,
        0.87598009f64,
        0.84068475f64,
        0.80211977f64,
        0.76076565f64,
        0.71712752f64,
        0.67172623f64,
        0.62508937f64,
        0.57774224f64,
        0.53019925f64,
        0.48295561f64,
        0.43647969f64,
        0.39120616f64,
        0.34752997f64,
        0.30580127f64,
        0.26632152f64,
        0.22934058f64,
        0.19505503f64,
        0.16360756f64,
        0.13508755f64,
        0.10953262f64,
        0.0869312f64,
        0.067226f64,
        0.0503182f64,
        0.03607231f64,
        0.02432151f64,
        0.01487334f64,
        0.00752f64,
        0.0f64,
    ]
};
unsafe extern "C" fn sinc(
    mut cutoff: libc::c_float,
    mut x: libc::c_float,
    mut N: libc::c_int,
    mut window_func: *const FuncDef,
) -> spx_word16_t {
    let mut xx: libc::c_float = x * cutoff;
    if (x as libc::c_double).abs() < 0.000001f64 {
        return cutoff;
    } else if (x as libc::c_double).abs() > 0.5f64 * N as libc::c_double {
        return 0i32 as spx_word16_t;
    } else {
        return (cutoff as libc::c_double
            * (3.141592653589793f64 * xx as libc::c_double).sin()
            / (3.141592653589793f64 * xx as libc::c_double)
            * compute_func(
                (2.0f64 * x as libc::c_double / N as libc::c_double).abs()
                    as libc::c_float,
                window_func,
            )) as spx_word16_t;
    };
}
unsafe extern "C" fn compute_func(
    mut x: libc::c_float,
    mut func: *const FuncDef,
) -> libc::c_double {
    let mut y: libc::c_float = 0.;
    let mut frac: libc::c_float = 0.;
    let mut interp: [libc::c_double; 4] = [0.; 4];
    let mut ind: libc::c_int = 0;
    y = x * (*func).oversample as libc::c_float;
    ind = (y as libc::c_double).floor() as libc::c_int;
    frac = y - ind as libc::c_float;
    interp[3usize] = -0.1666666667f64 * frac as libc::c_double
        + 0.1666666667f64 * (frac * frac * frac) as libc::c_double;
    interp[2usize] = frac as libc::c_double
        + 0.5f64 * (frac * frac) as libc::c_double
        - 0.5f64 * (frac * frac * frac) as libc::c_double;
    interp[0usize] = -0.3333333333f64 * frac as libc::c_double
        + 0.5f64 * (frac * frac) as libc::c_double
        - 0.1666666667f64 * (frac * frac * frac) as libc::c_double;
    interp[1usize] = 1.0f32 as libc::c_double
        - interp[3usize]
        - interp[2usize]
        - interp[0usize];
    return interp[0usize] * *(*func).table.offset(ind as isize)
        + interp[1usize] * *(*func).table.offset((ind + 1i32) as isize)
        + interp[2usize] * *(*func).table.offset((ind + 2i32) as isize)
        + interp[3usize] * *(*func).table.offset((ind + 3i32) as isize);
}
unsafe extern "C" fn resampler_basic_direct_single(
    mut st: *mut SpeexResamplerState,
    mut channel_index: spx_uint32_t,
    mut in_0: *const spx_word16_t,
    mut in_len: *mut spx_uint32_t,
    mut out: *mut spx_word16_t,
    mut out_len: *mut spx_uint32_t,
) -> libc::c_int {
    let N: libc::c_int = (*st).filt_len as libc::c_int;
    let mut out_sample: libc::c_int = 0i32;
    let mut last_sample: libc::c_int =
        *(*st).last_sample.offset(channel_index as isize);
    let mut samp_frac_num: spx_uint32_t =
        *(*st).samp_frac_num.offset(channel_index as isize);
    let mut sinc_table: *const spx_word16_t = (*st).sinc_table;
    let out_stride: libc::c_int = (*st).out_stride;
    let int_advance: libc::c_int = (*st).int_advance;
    let frac_advance: libc::c_int = (*st).frac_advance;
    let den_rate: spx_uint32_t = (*st).den_rate;
    let mut sum: spx_word32_t = 0.;
    while !(last_sample >= *in_len as spx_int32_t
        || out_sample >= *out_len as spx_int32_t)
    {
        let mut sinct: *const spx_word16_t = &*sinc_table
            .offset(samp_frac_num.wrapping_mul(N as libc::c_uint) as isize)
            as *const spx_word16_t;
        let mut iptr: *const spx_word16_t =
            &*in_0.offset(last_sample as isize) as *const spx_word16_t;
        let mut j: libc::c_int = 0;
        sum = 0i32 as spx_word32_t;
        j = 0i32;
        while j < N {
            sum += *sinct.offset(j as isize) * *iptr.offset(j as isize);
            j += 1
        }
        sum = sum;
        let fresh7 = out_sample;
        out_sample = out_sample + 1;
        *out.offset((out_stride * fresh7) as isize) = sum;
        last_sample += int_advance;
        samp_frac_num = (samp_frac_num as libc::c_uint)
            .wrapping_add(frac_advance as libc::c_uint)
            as spx_uint32_t as spx_uint32_t;
        if !(samp_frac_num >= den_rate) {
            continue;
        }
        samp_frac_num = (samp_frac_num as libc::c_uint).wrapping_sub(den_rate)
            as spx_uint32_t as spx_uint32_t;
        last_sample += 1
    }
    *(*st).last_sample.offset(channel_index as isize) = last_sample;
    *(*st).samp_frac_num.offset(channel_index as isize) = samp_frac_num;
    return out_sample;
}
unsafe extern "C" fn resampler_basic_direct_double(
    mut st: *mut SpeexResamplerState,
    mut channel_index: spx_uint32_t,
    mut in_0: *const spx_word16_t,
    mut in_len: *mut spx_uint32_t,
    mut out: *mut spx_word16_t,
    mut out_len: *mut spx_uint32_t,
) -> libc::c_int {
    let N: libc::c_int = (*st).filt_len as libc::c_int;
    let mut out_sample: libc::c_int = 0i32;
    let mut last_sample: libc::c_int =
        *(*st).last_sample.offset(channel_index as isize);
    let mut samp_frac_num: spx_uint32_t =
        *(*st).samp_frac_num.offset(channel_index as isize);
    let mut sinc_table: *const spx_word16_t = (*st).sinc_table;
    let out_stride: libc::c_int = (*st).out_stride;
    let int_advance: libc::c_int = (*st).int_advance;
    let frac_advance: libc::c_int = (*st).frac_advance;
    let den_rate: spx_uint32_t = (*st).den_rate;
    let mut sum: libc::c_double = 0.;
    while !(last_sample >= *in_len as spx_int32_t
        || out_sample >= *out_len as spx_int32_t)
    {
        let mut sinct: *const spx_word16_t = &*sinc_table
            .offset(samp_frac_num.wrapping_mul(N as libc::c_uint) as isize)
            as *const spx_word16_t;
        let mut iptr: *const spx_word16_t =
            &*in_0.offset(last_sample as isize) as *const spx_word16_t;
        let mut j: libc::c_int = 0;
        let mut accum: [libc::c_double; 4] = [
            0i32 as libc::c_double,
            0i32 as libc::c_double,
            0i32 as libc::c_double,
            0i32 as libc::c_double,
        ];
        j = 0i32;
        while j < N {
            accum[0usize] += (*sinct.offset(j as isize)
                * *iptr.offset(j as isize))
                as libc::c_double;
            accum[1usize] += (*sinct.offset((j + 1i32) as isize)
                * *iptr.offset((j + 1i32) as isize))
                as libc::c_double;
            accum[2usize] += (*sinct.offset((j + 2i32) as isize)
                * *iptr.offset((j + 2i32) as isize))
                as libc::c_double;
            accum[3usize] += (*sinct.offset((j + 3i32) as isize)
                * *iptr.offset((j + 3i32) as isize))
                as libc::c_double;
            j += 4i32
        }
        sum = accum[0usize] + accum[1usize] + accum[2usize] + accum[3usize];
        let fresh8 = out_sample;
        out_sample = out_sample + 1;
        *out.offset((out_stride * fresh8) as isize) = sum as spx_word16_t;
        last_sample += int_advance;
        samp_frac_num = (samp_frac_num as libc::c_uint)
            .wrapping_add(frac_advance as libc::c_uint)
            as spx_uint32_t as spx_uint32_t;
        if !(samp_frac_num >= den_rate) {
            continue;
        }
        samp_frac_num = (samp_frac_num as libc::c_uint).wrapping_sub(den_rate)
            as spx_uint32_t as spx_uint32_t;
        last_sample += 1
    }
    *(*st).last_sample.offset(channel_index as isize) = last_sample;
    *(*st).samp_frac_num.offset(channel_index as isize) = samp_frac_num;
    return out_sample;
}
unsafe extern "C" fn _muldiv(
    mut result: *mut spx_uint32_t,
    mut value: spx_uint32_t,
    mut mul: spx_uint32_t,
    mut div: spx_uint32_t,
) -> libc::c_int {
    if result.is_null() {
        _speex_fatal(
            (*::std::mem::transmute::<&[u8; 25], &mut [libc::c_char; 25]>(
                b"assertion failed: result\x00",
            ))
            .as_mut_ptr(),
            (*::std::mem::transmute::<&[u8; 11], &mut [libc::c_char; 11]>(
                b"resample.c\x00",
            ))
            .as_mut_ptr(),
            594i32,
        );
    }
    let mut major: spx_uint32_t = value.wrapping_div(div);
    let mut remainder: spx_uint32_t = value.wrapping_rem(div);
    if remainder > 4294967295u32.wrapping_div(mul)
        || major > 4294967295u32.wrapping_div(mul)
        || major.wrapping_mul(mul)
            > 4294967295u32
                .wrapping_sub(remainder.wrapping_mul(mul).wrapping_div(div))
    {
        return RESAMPLER_ERR_OVERFLOW as libc::c_int;
    } else {
        *result = remainder
            .wrapping_mul(mul)
            .wrapping_div(div)
            .wrapping_add(major.wrapping_mul(mul));
        return RESAMPLER_ERR_SUCCESS as libc::c_int;
    };
}
/* * For n elements worth of memory, set every byte to the value of c, starting at address dst */
/* * Copy n elements from src to dst, allowing overlapping regions. The 0* term
provides compile-time type checking */
/* * Copy n elements from src to dst. The 0* term provides compile-time type checking  */
unsafe extern "C" fn _speex_fatal(
    mut st: *const libc::c_char,
    mut file: *const libc::c_char,
    mut line: libc::c_int,
) -> () {
    use std::ffi::CStr;

    eprintln!(
        "Fatal (internal) error in {}, line {}: {}",
        CStr::from_ptr(file).to_string_lossy(),
        line,
        CStr::from_ptr(st).to_string_lossy()
    );
    exit(1i32);
}
/* * Set (change) the input/output sampling rates and resampling ratio
 * (fractional values in Hz supported).
 * @param st Resampler state
 * @param ratio_num Numerator of the sampling rate ratio
 * @param ratio_den Denominator of the sampling rate ratio
 * @param in_rate Input sampling rate rounded to the nearest integer (in Hz).
 * @param out_rate Output sampling rate rounded to the nearest integer (in Hz).
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_set_rate_frac(
    mut st: *mut SpeexResamplerState,
    mut ratio_num: spx_uint32_t,
    mut ratio_den: spx_uint32_t,
    mut in_rate: spx_uint32_t,
    mut out_rate: spx_uint32_t,
) -> libc::c_int {
    let mut fact: spx_uint32_t = 0;
    let mut old_den: spx_uint32_t = 0;
    let mut i: spx_uint32_t = 0;
    if ratio_num == 0i32 as libc::c_uint || ratio_den == 0i32 as libc::c_uint {
        return RESAMPLER_ERR_INVALID_ARG as libc::c_int;
    } else if (*st).in_rate == in_rate
        && (*st).out_rate == out_rate
        && (*st).num_rate == ratio_num
        && (*st).den_rate == ratio_den
    {
        return RESAMPLER_ERR_SUCCESS as libc::c_int;
    } else {
        old_den = (*st).den_rate;
        (*st).in_rate = in_rate;
        (*st).out_rate = out_rate;
        (*st).num_rate = ratio_num;
        (*st).den_rate = ratio_den;
        fact = _gcd((*st).num_rate, (*st).den_rate);
        (*st).num_rate = ((*st).num_rate as libc::c_uint).wrapping_div(fact)
            as spx_uint32_t as spx_uint32_t;
        (*st).den_rate = ((*st).den_rate as libc::c_uint).wrapping_div(fact)
            as spx_uint32_t as spx_uint32_t;
        if old_den > 0i32 as libc::c_uint {
            i = 0i32 as spx_uint32_t;
            while i < (*st).nb_channels {
                if _muldiv(
                    &mut *(*st).samp_frac_num.offset(i as isize)
                        as *mut spx_uint32_t,
                    *(*st).samp_frac_num.offset(i as isize),
                    (*st).den_rate,
                    old_den,
                ) != RESAMPLER_ERR_SUCCESS as libc::c_int
                {
                    return RESAMPLER_ERR_OVERFLOW as libc::c_int;
                } else {
                    if *(*st).samp_frac_num.offset(i as isize)
                        >= (*st).den_rate
                    {
                        *(*st).samp_frac_num.offset(i as isize) =
                            (*st).den_rate.wrapping_sub(1i32 as libc::c_uint)
                    }
                    i = i.wrapping_add(1)
                }
            }
        }
        if 0 != (*st).initialised {
            return update_filter(st);
        } else {
            return RESAMPLER_ERR_SUCCESS as libc::c_int;
        }
    };
}
unsafe extern "C" fn _gcd(
    mut a: spx_uint32_t,
    mut b: spx_uint32_t,
) -> spx_uint32_t {
    let mut temp: spx_uint32_t = 0;
    while b != 0i32 as libc::c_uint {
        temp = a;
        a = b;
        b = temp.wrapping_rem(b)
    }
    return a;
}
/* * Set (change) the conversion quality.
 * @param st Resampler state
 * @param quality Resampling quality between 0 and 10, where 0 has poor
 * quality and 10 has very high quality.
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_set_quality(
    mut st: *mut SpeexResamplerState,
    mut quality: libc::c_int,
) -> libc::c_int {
    if quality > 10i32 || quality < 0i32 {
        return RESAMPLER_ERR_INVALID_ARG as libc::c_int;
    } else if (*st).quality == quality {
        return RESAMPLER_ERR_SUCCESS as libc::c_int;
    } else {
        (*st).quality = quality;
        if 0 != (*st).initialised {
            return update_filter(st);
        } else {
            return RESAMPLER_ERR_SUCCESS as libc::c_int;
        }
    };
}
/* * Speex wrapper for calloc. To do your own dynamic allocation, all you need to do is replace this function, speex_realloc and speex_free
NOTE: speex_alloc needs to CLEAR THE MEMORY */
unsafe extern "C" fn speex_alloc(mut size: libc::c_int) -> *mut libc::c_void {
    return calloc(size as usize, 1);
}
/* * Resample a float array. The input and output buffers must *not* overlap.
 * @param st Resampler state
 * @param channel_index Index of the channel to process for the multi-channel
 * base (0 otherwise)
 * @param in Input buffer
 * @param in_len Number of input samples in the input buffer. Returns the
 * number of samples processed
 * @param out Output buffer
 * @param out_len Size of the output buffer. Returns the number of samples written
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_process_float(
    mut st: *mut SpeexResamplerState,
    mut channel_index: spx_uint32_t,
    mut in_0: *const libc::c_float,
    mut in_len: *mut spx_uint32_t,
    mut out: *mut libc::c_float,
    mut out_len: *mut spx_uint32_t,
) -> libc::c_int {
    let mut j: libc::c_int = 0;
    let mut ilen: spx_uint32_t = *in_len;
    let mut olen: spx_uint32_t = *out_len;
    let mut x: *mut spx_word16_t = (*st)
        .mem
        .offset(channel_index.wrapping_mul((*st).mem_alloc_size) as isize);
    let filt_offs: libc::c_int =
        (*st).filt_len.wrapping_sub(1i32 as libc::c_uint) as libc::c_int;
    let xlen: spx_uint32_t =
        (*st).mem_alloc_size.wrapping_sub(filt_offs as libc::c_uint);
    let istride: libc::c_int = (*st).in_stride;
    if 0 != *(*st).magic_samples.offset(channel_index as isize) {
        olen = (olen as libc::c_uint).wrapping_sub(speex_resampler_magic(
            st,
            channel_index,
            &mut out as *mut *mut libc::c_float,
            olen,
        ) as libc::c_uint) as spx_uint32_t as spx_uint32_t
    }
    if 0 == *(*st).magic_samples.offset(channel_index as isize) {
        while 0 != ilen && 0 != olen {
            let mut ichunk: spx_uint32_t =
                if ilen > xlen { xlen } else { ilen };
            let mut ochunk: spx_uint32_t = olen;
            if !in_0.is_null() {
                j = 0i32;
                while (j as libc::c_uint) < ichunk {
                    *x.offset((j + filt_offs) as isize) =
                        *in_0.offset((j * istride) as isize);
                    j += 1
                }
            } else {
                j = 0i32;
                while (j as libc::c_uint) < ichunk {
                    *x.offset((j + filt_offs) as isize) = 0i32 as spx_word16_t;
                    j += 1
                }
            }
            speex_resampler_process_native(
                st,
                channel_index,
                &mut ichunk as *mut spx_uint32_t,
                out,
                &mut ochunk as *mut spx_uint32_t,
            );
            ilen = (ilen as libc::c_uint).wrapping_sub(ichunk) as spx_uint32_t
                as spx_uint32_t;
            olen = (olen as libc::c_uint).wrapping_sub(ochunk) as spx_uint32_t
                as spx_uint32_t;
            out = out
                .offset(ochunk.wrapping_mul((*st).out_stride as libc::c_uint)
                    as isize);
            if in_0.is_null() {
                continue;
            }
            in_0 = in_0
                .offset(ichunk.wrapping_mul(istride as libc::c_uint) as isize)
        }
    }
    *in_len = (*in_len as libc::c_uint).wrapping_sub(ilen) as spx_uint32_t
        as spx_uint32_t;
    *out_len = (*out_len as libc::c_uint).wrapping_sub(olen) as spx_uint32_t
        as spx_uint32_t;
    return if (*st).resampler_ptr == Some(resampler_basic_zero) {
        RESAMPLER_ERR_ALLOC_FAILED as libc::c_int
    } else {
        RESAMPLER_ERR_SUCCESS as libc::c_int
    };
}
unsafe extern "C" fn speex_resampler_process_native(
    mut st: *mut SpeexResamplerState,
    mut channel_index: spx_uint32_t,
    mut in_len: *mut spx_uint32_t,
    mut out: *mut spx_word16_t,
    mut out_len: *mut spx_uint32_t,
) -> libc::c_int {
    let mut j: libc::c_int = 0i32;
    let N: libc::c_int = (*st).filt_len as libc::c_int;
    let mut out_sample: libc::c_int = 0i32;
    let mut mem: *mut spx_word16_t = (*st)
        .mem
        .offset(channel_index.wrapping_mul((*st).mem_alloc_size) as isize);
    let mut ilen: spx_uint32_t = 0;
    (*st).started = 1i32;
    out_sample = (*st).resampler_ptr.expect("non-null function pointer")(
        st,
        channel_index,
        mem,
        in_len,
        out,
        out_len,
    );
    if *(*st).last_sample.offset(channel_index as isize)
        < *in_len as spx_int32_t
    {
        *in_len =
            *(*st).last_sample.offset(channel_index as isize) as spx_uint32_t
    }
    *out_len = out_sample as spx_uint32_t;
    let ref mut fresh9 = *(*st).last_sample.offset(channel_index as isize);
    *fresh9 = (*fresh9 as libc::c_uint).wrapping_sub(*in_len) as spx_int32_t
        as spx_int32_t;
    ilen = *in_len;
    j = 0i32;
    while j < N - 1i32 {
        *mem.offset(j as isize) =
            *mem.offset((j as libc::c_uint).wrapping_add(ilen) as isize);
        j += 1
    }
    return RESAMPLER_ERR_SUCCESS as libc::c_int;
}
unsafe extern "C" fn speex_resampler_magic(
    mut st: *mut SpeexResamplerState,
    mut channel_index: spx_uint32_t,
    mut out: *mut *mut spx_word16_t,
    mut out_len: spx_uint32_t,
) -> libc::c_int {
    let mut tmp_in_len: spx_uint32_t =
        *(*st).magic_samples.offset(channel_index as isize);
    let mut mem: *mut spx_word16_t = (*st)
        .mem
        .offset(channel_index.wrapping_mul((*st).mem_alloc_size) as isize);
    let N: libc::c_int = (*st).filt_len as libc::c_int;
    speex_resampler_process_native(
        st,
        channel_index,
        &mut tmp_in_len as *mut spx_uint32_t,
        *out,
        &mut out_len as *mut spx_uint32_t,
    );
    let ref mut fresh10 = *(*st).magic_samples.offset(channel_index as isize);
    *fresh10 = (*fresh10 as libc::c_uint).wrapping_sub(tmp_in_len)
        as spx_uint32_t as spx_uint32_t;
    if 0 != *(*st).magic_samples.offset(channel_index as isize) {
        let mut i: spx_uint32_t = 0;
        i = 0i32 as spx_uint32_t;
        while i < *(*st).magic_samples.offset(channel_index as isize) {
            *mem.offset(
                ((N - 1i32) as libc::c_uint).wrapping_add(i) as isize
            ) = *mem.offset(
                ((N - 1i32) as libc::c_uint)
                    .wrapping_add(i)
                    .wrapping_add(tmp_in_len) as isize,
            );
            i = i.wrapping_add(1)
        }
    }
    *out =
        (*out)
            .offset(out_len.wrapping_mul((*st).out_stride as libc::c_uint)
                as isize);
    return out_len as libc::c_int;
}
/* * Resample an int array. The input and output buffers must *not* overlap.
 * @param st Resampler state
 * @param channel_index Index of the channel to process for the multi-channel
 * base (0 otherwise)
 * @param in Input buffer
 * @param in_len Number of input samples in the input buffer. Returns the number
 * of samples processed
 * @param out Output buffer
 * @param out_len Size of the output buffer. Returns the number of samples written
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_process_int(
    mut st: *mut SpeexResamplerState,
    mut channel_index: spx_uint32_t,
    mut in_0: *const spx_int16_t,
    mut in_len: *mut spx_uint32_t,
    mut out: *mut spx_int16_t,
    mut out_len: *mut spx_uint32_t,
) -> libc::c_int {
    let mut j: libc::c_int = 0;
    let istride_save: libc::c_int = (*st).in_stride;
    let ostride_save: libc::c_int = (*st).out_stride;
    let mut ilen: spx_uint32_t = *in_len;
    let mut olen: spx_uint32_t = *out_len;
    let mut x: *mut spx_word16_t = (*st)
        .mem
        .offset(channel_index.wrapping_mul((*st).mem_alloc_size) as isize);
    let xlen: spx_uint32_t = (*st)
        .mem_alloc_size
        .wrapping_sub((*st).filt_len.wrapping_sub(1i32 as libc::c_uint));
    let ylen: libc::c_uint = if olen < 8192i32 as libc::c_uint {
        olen
    } else {
        8192i32 as libc::c_uint
    };
    let vla = ylen as usize;
    let mut ystack: Vec<spx_word16_t> = ::std::vec::from_elem(0., vla);
    (*st).out_stride = 1i32;
    while 0 != ilen && 0 != olen {
        let mut y: *mut spx_word16_t = ystack.as_mut_ptr();
        let mut ichunk: spx_uint32_t = if ilen > xlen { xlen } else { ilen };
        let mut ochunk: spx_uint32_t = if olen > ylen { ylen } else { olen };
        let mut omagic: spx_uint32_t = 0i32 as spx_uint32_t;
        if 0 != *(*st).magic_samples.offset(channel_index as isize) {
            omagic = speex_resampler_magic(
                st,
                channel_index,
                &mut y as *mut *mut spx_word16_t,
                ochunk,
            ) as spx_uint32_t;
            ochunk = (ochunk as libc::c_uint).wrapping_sub(omagic)
                as spx_uint32_t as spx_uint32_t;
            olen = (olen as libc::c_uint).wrapping_sub(omagic) as spx_uint32_t
                as spx_uint32_t
        }
        if 0 == *(*st).magic_samples.offset(channel_index as isize) {
            if !in_0.is_null() {
                j = 0i32;
                while (j as libc::c_uint) < ichunk {
                    *x.offset(
                        (j as libc::c_uint)
                            .wrapping_add((*st).filt_len)
                            .wrapping_sub(1i32 as libc::c_uint)
                            as isize,
                    ) = *in_0.offset((j * istride_save) as isize)
                        as spx_word16_t;
                    j += 1
                }
            } else {
                j = 0i32;
                while (j as libc::c_uint) < ichunk {
                    *x.offset(
                        (j as libc::c_uint)
                            .wrapping_add((*st).filt_len)
                            .wrapping_sub(1i32 as libc::c_uint)
                            as isize,
                    ) = 0i32 as spx_word16_t;
                    j += 1
                }
            }
            speex_resampler_process_native(
                st,
                channel_index,
                &mut ichunk as *mut spx_uint32_t,
                y,
                &mut ochunk as *mut spx_uint32_t,
            );
        } else {
            ichunk = 0i32 as spx_uint32_t;
            ochunk = 0i32 as spx_uint32_t
        }
        j = 0i32;
        while (j as libc::c_uint) < ochunk.wrapping_add(omagic) {
            *out.offset((j * ostride_save) as isize) =
                (if *ystack.as_mut_ptr().offset(j as isize) < -32767.5f32 {
                    -32768i32
                } else if *ystack.as_mut_ptr().offset(j as isize) > 32766.5f32
                {
                    32767i32
                } else {
                    (0.5f64
                        + *ystack.as_mut_ptr().offset(j as isize)
                            as libc::c_double)
                        .floor() as spx_int16_t
                        as libc::c_int
                }) as spx_int16_t;
            j += 1
        }
        ilen = (ilen as libc::c_uint).wrapping_sub(ichunk) as spx_uint32_t
            as spx_uint32_t;
        olen = (olen as libc::c_uint).wrapping_sub(ochunk) as spx_uint32_t
            as spx_uint32_t;
        out = out.offset(
            ochunk
                .wrapping_add(omagic)
                .wrapping_mul(ostride_save as libc::c_uint)
                as isize,
        );
        if in_0.is_null() {
            continue;
        }
        in_0 = in_0
            .offset(ichunk.wrapping_mul(istride_save as libc::c_uint) as isize)
    }
    (*st).out_stride = ostride_save;
    *in_len = (*in_len as libc::c_uint).wrapping_sub(ilen) as spx_uint32_t
        as spx_uint32_t;
    *out_len = (*out_len as libc::c_uint).wrapping_sub(olen) as spx_uint32_t
        as spx_uint32_t;
    return if (*st).resampler_ptr == Some(resampler_basic_zero) {
        RESAMPLER_ERR_ALLOC_FAILED as libc::c_int
    } else {
        RESAMPLER_ERR_SUCCESS as libc::c_int
    };
}
/* * Resample an interleaved float array. The input and output buffers must *not* overlap.
 * @param st Resampler state
 * @param in Input buffer
 * @param in_len Number of input samples in the input buffer. Returns the number
 * of samples processed. This is all per-channel.
 * @param out Output buffer
 * @param out_len Size of the output buffer. Returns the number of samples written.
 * This is all per-channel.
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_process_interleaved_float(
    mut st: *mut SpeexResamplerState,
    mut in_0: *const libc::c_float,
    mut in_len: *mut spx_uint32_t,
    mut out: *mut libc::c_float,
    mut out_len: *mut spx_uint32_t,
) -> libc::c_int {
    let mut i: spx_uint32_t = 0;
    let mut bak_out_len: spx_uint32_t = *out_len;
    let mut bak_in_len: spx_uint32_t = *in_len;
    let mut istride_save = (*st).in_stride;
    let mut ostride_save = (*st).out_stride;
    (*st).out_stride = (*st).nb_channels as libc::c_int;
    (*st).in_stride = (*st).out_stride;
    while i < (*st).nb_channels {
        *out_len = bak_out_len;
        *in_len = bak_in_len;
        if !in_0.is_null() {
            speex_resampler_process_float(
                st,
                i,
                in_0.offset(i as isize),
                in_len,
                out.offset(i as isize),
                out_len,
            );
        } else {
            speex_resampler_process_float(
                st,
                i,
                0 as *const libc::c_float,
                in_len,
                out.offset(i as isize),
                out_len,
            );
        }
        i = i.wrapping_add(1)
    }
    (*st).in_stride = istride_save;
    (*st).out_stride = ostride_save;
    return if (*st).resampler_ptr == Some(resampler_basic_zero) {
        RESAMPLER_ERR_ALLOC_FAILED as libc::c_int
    } else {
        RESAMPLER_ERR_SUCCESS as libc::c_int
    };
}
/* * Resample an interleaved int array. The input and output buffers must *not* overlap.
 * @param st Resampler state
 * @param in Input buffer
 * @param in_len Number of input samples in the input buffer. Returns the number
 * of samples processed. This is all per-channel.
 * @param out Output buffer
 * @param out_len Size of the output buffer. Returns the number of samples written.
 * This is all per-channel.
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_process_interleaved_int(
    mut st: *mut SpeexResamplerState,
    mut in_0: *const spx_int16_t,
    mut in_len: *mut spx_uint32_t,
    mut out: *mut spx_int16_t,
    mut out_len: *mut spx_uint32_t,
) -> libc::c_int {
    let mut i: spx_uint32_t = 0;
    let mut bak_out_len: spx_uint32_t = *out_len;
    let mut bak_in_len: spx_uint32_t = *in_len;
    let mut istride_save = (*st).in_stride;
    let mut ostride_save = (*st).out_stride;
    (*st).out_stride = (*st).nb_channels as libc::c_int;
    (*st).in_stride = (*st).out_stride;
    while i < (*st).nb_channels {
        *out_len = bak_out_len;
        *in_len = bak_in_len;
        if !in_0.is_null() {
            speex_resampler_process_int(
                st,
                i,
                in_0.offset(i as isize),
                in_len,
                out.offset(i as isize),
                out_len,
            );
        } else {
            speex_resampler_process_int(
                st,
                i,
                0 as *const spx_int16_t,
                in_len,
                out.offset(i as isize),
                out_len,
            );
        }
        i = i.wrapping_add(1)
    }
    (*st).in_stride = istride_save;
    (*st).out_stride = ostride_save;
    return if (*st).resampler_ptr == Some(resampler_basic_zero) {
        RESAMPLER_ERR_ALLOC_FAILED as libc::c_int
    } else {
        RESAMPLER_ERR_SUCCESS as libc::c_int
    };
}
/* * Set (change) the input/output sampling rates (integer value).
 * @param st Resampler state
 * @param in_rate Input sampling rate (integer number of Hz).
 * @param out_rate Output sampling rate (integer number of Hz).
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_set_rate(
    mut st: *mut SpeexResamplerState,
    mut in_rate: spx_uint32_t,
    mut out_rate: spx_uint32_t,
) -> libc::c_int {
    return speex_resampler_set_rate_frac(
        st, in_rate, out_rate, in_rate, out_rate,
    );
}
/* * Get the current input/output sampling rates (integer value).
 * @param st Resampler state
 * @param in_rate Input sampling rate (integer number of Hz) copied.
 * @param out_rate Output sampling rate (integer number of Hz) copied.
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_get_rate(
    mut st: *mut SpeexResamplerState,
    mut in_rate: *mut spx_uint32_t,
    mut out_rate: *mut spx_uint32_t,
) -> () {
    *in_rate = (*st).in_rate;
    *out_rate = (*st).out_rate;
}
/* * Get the current resampling ratio. This will be reduced to the least
 * common denominator.
 * @param st Resampler state
 * @param ratio_num Numerator of the sampling rate ratio copied
 * @param ratio_den Denominator of the sampling rate ratio copied
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_get_ratio(
    mut st: *mut SpeexResamplerState,
    mut ratio_num: *mut spx_uint32_t,
    mut ratio_den: *mut spx_uint32_t,
) -> () {
    *ratio_num = (*st).num_rate;
    *ratio_den = (*st).den_rate;
}
/* * Get the conversion quality.
 * @param st Resampler state
 * @param quality Resampling quality between 0 and 10, where 0 has poor
 * quality and 10 has very high quality.
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_get_quality(
    mut st: *mut SpeexResamplerState,
    mut quality: *mut libc::c_int,
) -> () {
    *quality = (*st).quality;
}
/* * Set (change) the input stride.
 * @param st Resampler state
 * @param stride Input stride
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_set_input_stride(
    mut st: *mut SpeexResamplerState,
    mut stride: spx_uint32_t,
) -> () {
    (*st).in_stride = stride as libc::c_int;
}
/* * Get the input stride.
 * @param st Resampler state
 * @param stride Input stride copied
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_get_input_stride(
    mut st: *mut SpeexResamplerState,
    mut stride: *mut spx_uint32_t,
) -> () {
    *stride = (*st).in_stride as spx_uint32_t;
}
/* * Set (change) the output stride.
 * @param st Resampler state
 * @param stride Output stride
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_set_output_stride(
    mut st: *mut SpeexResamplerState,
    mut stride: spx_uint32_t,
) -> () {
    (*st).out_stride = stride as libc::c_int;
}
/* * Get the output stride.
 * @param st Resampler state copied
 * @param stride Output stride
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_get_output_stride(
    mut st: *mut SpeexResamplerState,
    mut stride: *mut spx_uint32_t,
) -> () {
    *stride = (*st).out_stride as spx_uint32_t;
}
/* * Get the latency introduced by the resampler measured in input samples.
 * @param st Resampler state
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_get_input_latency(
    mut st: *mut SpeexResamplerState,
) -> libc::c_int {
    return (*st).filt_len.wrapping_div(2i32 as libc::c_uint) as libc::c_int;
}
/* * Get the latency introduced by the resampler measured in output samples.
 * @param st Resampler state
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_get_output_latency(
    mut st: *mut SpeexResamplerState,
) -> libc::c_int {
    return (*st)
        .filt_len
        .wrapping_div(2i32 as libc::c_uint)
        .wrapping_mul((*st).den_rate)
        .wrapping_add((*st).num_rate >> 1i32)
        .wrapping_div((*st).num_rate) as libc::c_int;
}
/* * Make sure that the first samples to go out of the resamplers don't have
 * leading zeros. This is only useful before starting to use a newly created
 * resampler. It is recommended to use that when resampling an audio file, as
 * it will generate a file with the same length. For real-time processing,
 * it is probably easier not to use this call (so that the output duration
 * is the same for the first frame).
 * @param st Resampler state
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_skip_zeros(
    mut st: *mut SpeexResamplerState,
) -> libc::c_int {
    let mut i: spx_uint32_t = 0;
    while i < (*st).nb_channels {
        *(*st).last_sample.offset(i as isize) =
            (*st).filt_len.wrapping_div(2i32 as libc::c_uint) as spx_int32_t;
        i = i.wrapping_add(1)
    }
    return RESAMPLER_ERR_SUCCESS as libc::c_int;
}
/* * Reset a resampler so a new (unrelated) stream can be processed.
 * @param st Resampler state
 */
#[no_mangle]
pub unsafe extern "C" fn speex_resampler_reset_mem(
    mut st: *mut SpeexResamplerState,
) -> libc::c_int {
    let mut i: spx_uint32_t = 0;
    while i < (*st).nb_channels {
        *(*st).last_sample.offset(i as isize) = 0i32;
        *(*st).magic_samples.offset(i as isize) = 0i32 as spx_uint32_t;
        *(*st).samp_frac_num.offset(i as isize) = 0i32 as spx_uint32_t;
        i = i.wrapping_add(1)
    }
    i = 0i32 as spx_uint32_t;
    while i
        < (*st)
            .nb_channels
            .wrapping_mul((*st).filt_len.wrapping_sub(1i32 as libc::c_uint))
    {
        *(*st).mem.offset(i as isize) = 0i32 as spx_word16_t;
        i = i.wrapping_add(1)
    }
    return RESAMPLER_ERR_SUCCESS as libc::c_int;
}
