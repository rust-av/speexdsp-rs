#![allow(dead_code)]

/// [speexdsp][1]-derived audio resampler
///
/// [1]: https://github.com/xiph/speexdsp
use std::mem;

use std::f64::consts::PI as PI_64;

pub const RESAMPLER_ERR_SUCCESS: usize = 0;
pub const RESAMPLER_ERR_ALLOC_FAILED: usize = 1;
pub const RESAMPLER_ERR_INVALID_ARG: usize = 3;
pub const RESAMPLER_ERR_OVERFLOW: usize = 5;

#[derive(Clone)]
pub struct SpeexResamplerState {
    in_rate: u32,
    out_rate: u32,
    num_rate: u32,
    den_rate: u32,
    quality: u32,
    nb_channels: u32,
    filt_len: u32,
    mem_alloc_size: u32,
    buffer_size: u32,
    int_advance: u32,
    frac_advance: u32,
    cutoff: f32,
    oversample: u32,
    initialised: u32,
    started: u32,
    last_sample: Vec<u32>,
    samp_frac_num: Vec<u32>,
    magic_samples: Vec<u32>,
    mem: Vec<f32>,
    sinc_table: Vec<f32>,
    sinc_table_length: u32,
    resampler_ptr: ResamplerBasicFunc,
    in_stride: u32,
    out_stride: u32,
}

#[derive(Copy, Clone)]
pub struct FuncDef {
    table: &'static [f64],
    oversample: usize,
}

impl FuncDef {
    pub const fn new(table: &'static [f64], oversample: usize) -> Self {
        Self { table, oversample }
    }
}

#[derive(Copy, Clone)]
pub struct QualityMapping {
    base_length: usize,
    oversample: usize,
    downsample_bandwidth: f32,
    upsample_bandwidth: f32,
    window_func: &'static FuncDef,
}

impl QualityMapping {
    pub const fn new(
        base_length: usize,
        oversample: usize,
        downsample_bandwidth: f32,
        upsample_bandwidth: f32,
        window_func: &'static FuncDef,
    ) -> Self {
        Self {
            base_length,
            oversample,
            downsample_bandwidth,
            upsample_bandwidth,
            window_func,
        }
    }
}

pub type ResamplerBasicFunc = Option<
    fn(
        _: &mut SpeexResamplerState,
        _: u32,
        _: &[f32],
        _: &mut u32,
        _: &mut [f32],
        _: &mut u32,
    ) -> i32,
>;

macro_rules! chunk_assign {
    ($ch_mut:ident, $lbound_mut:expr, $ubound_mut:expr, $val:expr) => {
        $ch_mut[$lbound_mut as usize..$ubound_mut as usize]
            .iter_mut()
            .for_each(|x| *x = $val);
    };
}

macro_rules! chunk_copy {
    ($ch_mut:ident, $lbound_mut:expr, $ubound_mut:expr,
     $ch:ident, $lbound:expr, $ubound:expr) => {
        $ch_mut[$lbound_mut as usize..$ubound_mut as usize]
            .iter_mut()
            .zip($ch[$lbound as usize..$ubound as usize].iter())
            .for_each(|(x, y)| *x = *y);
    };
}

macro_rules! algo {
    ($self:ident, $ch_mut:ident, $ch:ident,
     $old_length:ident, $magic:expr) => {
        let olen = $old_length + 2 * $magic;
        let filt_len = $self.filt_len - 1;
        if $self.filt_len > olen {
            let new_filt_len = $self.filt_len - olen;
            for new_last_sample in &mut $self.last_sample {
                chunk_copy!($ch_mut, new_filt_len, filt_len, $ch, 0, olen - 1);
                chunk_assign!($ch_mut, 0, new_filt_len, 0.0);
                $magic = 0;
                *new_last_sample += new_filt_len / 2;
            }
        } else {
            $magic = (olen - $self.filt_len) / 2;
            let ubound_mut = filt_len + $magic;
            let ubound = ubound_mut + $magic;
            chunk_copy!($ch_mut, 0, ubound_mut, $ch, $magic, ubound);
        }
    };
}

impl SpeexResamplerState {
    /* * Create a new resampler with integer input and output rates.
     * @param nb_channels number of channels to be processed
     * @param in_rate Input sampling rate (integer number of Hz).
     * @param out_rate Output sampling rate (integer number of Hz).
     * @param quality Resampling quality between 0 and 10, where 0 has poor quality
     * and 10 has very high quality.
     * @return newly created resampler state
     */
    pub fn new(
        nb_channels: usize,
        in_rate: usize,
        out_rate: usize,
        quality: usize,
    ) -> Self {
        Self::init_frac(
            nb_channels,
            in_rate,
            out_rate,
            in_rate,
            out_rate,
            quality,
        )
    }

    /* * Create a new resampler with fractional input/output rates. The sampling
     * rate ratio is an arbitrary rational number with both the numerator and
     * denominator being 32-bit integers.
     * @param nb_channels number of channels to be processed
     * @param ratio_num numerator of the sampling rate ratio
     * @param ratio_den Denominator of the sampling rate ratio
     * @param in_rate Input sampling rate rounded to the nearest integer (in Hz).
     * @param out_rate Output sampling rate rounded to the nearest integer (in Hz).
     * @param quality Resampling quality between 0 and 10, where 0 has poor quality
     * and 10 has very high quality.
     * @return newly created resampler state
     * @retval nULL Error: not enough memory
     */
    pub fn init_frac(
        nb_channels: usize,
        ratio_num: usize,
        ratio_den: usize,
        in_rate: usize,
        out_rate: usize,
        quality: usize,
    ) -> Self {
        if nb_channels == 0 || ratio_num == 0 || ratio_den == 0 || quality > 10
        {
            panic!("Set the correct parameters!");
        }
        let mut st = Self {
            initialised: 0,
            started: 0,
            in_rate: 0,
            out_rate: 0,
            num_rate: 0,
            den_rate: 0,
            quality: 0,
            sinc_table: Vec::new(),
            sinc_table_length: 0,
            mem: Vec::new(),
            frac_advance: 0,
            int_advance: 0,
            mem_alloc_size: 0,
            filt_len: 0,
            resampler_ptr: None,
            cutoff: 1.0,
            nb_channels: nb_channels as u32,
            in_stride: 1,
            out_stride: 1,
            buffer_size: 160,
            oversample: 0,
            last_sample: vec![0; nb_channels as usize],
            magic_samples: vec![0; nb_channels as usize],
            samp_frac_num: vec![0; nb_channels as usize],
        };
        st.set_quality(quality);
        st.set_rate_frac(ratio_num, ratio_den, in_rate, out_rate);
        let filter_err = st.update_filter();
        if filter_err == RESAMPLER_ERR_SUCCESS {
            st.initialised = 1;
        } else {
            panic!("Error");
        }
        st
    }

    /* * Resample a float array. The input and output buffers must *not* overlap.
     * @param st Resampler state
     * @param channel_index Index of the channel to process for the multi-channel
     * base (0 otherwise)
     * @param in Input buffer
     * @param in_len number of input samples in the input buffer. Returns the
     * number of samples processed
     * @param out Output buffer
     * @param out_len Size of the output buffer. Returns the number of samples written
     */
    pub fn process_float(
        &mut self,
        channel_index: u32,
        mut in_0: &[f32],
        in_len: &mut u32,
        mut out: &mut [f32],
        out_len: &mut u32,
    ) -> usize {
        if in_0.is_empty() {
            panic!("Empty input slice is not allowed");
        }
        let mut ilen = *in_len;
        let mut olen = *out_len;
        let channel_idx = channel_index as usize;
        let filt_offs = (self.filt_len - 1) as usize;
        let mem_idx = filt_offs + channel_idx * self.mem_alloc_size as usize;
        let xlen = self.mem_alloc_size - self.filt_len - 1;
        let istride = self.in_stride as usize;
        if self.magic_samples[channel_idx] != 0 {
            olen -= speex_resampler_magic(self, channel_index, &mut out, olen);
        }
        if self.magic_samples[channel_idx] == 0 {
            while 0 != ilen && 0 != olen {
                let mut ichunk: u32 = if ilen > xlen { xlen } else { ilen };
                let mut ochunk: u32 = olen;
                let mem_slice = &mut self.mem[mem_idx..];
                let mem_iter = mem_slice.iter_mut();
                let in_iter = in_0.iter().step_by(istride);
                mem_iter.zip(in_iter).take(ichunk as usize).for_each(
                    |(x, &y)| {
                        *x = y;
                    },
                );
                speex_resampler_process_native(
                    self,
                    channel_index,
                    &mut ichunk,
                    out,
                    &mut ochunk,
                );
                ilen -= ichunk;
                olen -= ochunk;
                out = &mut out[(ochunk * self.out_stride) as usize..][..];
                in_0 = &in_0[(ichunk * self.in_stride) as usize..][..];
            }
        }
        *in_len -= ilen;
        *out_len -= olen;
        let resampler = self.resampler_ptr.unwrap();
        if resampler as usize == resampler_basic_zero as usize {
            RESAMPLER_ERR_ALLOC_FAILED
        } else {
            RESAMPLER_ERR_SUCCESS
        }
    }

    /* * Resample an int array. The input and output buffers muself *not* overlap.
     * @param self Resampler selfate
     * @param channel_index Index of the channel to process for the multi-channel
     * base (0 otherwise)
     * @param in Input buffer
     * @param in_len number of input samples in the input buffer. Returns the number
     * of samples processed
     * @param out Output buffer
     * @param out_len Size of the output buffer. Returns the number of samples written
     */
    pub fn process_int(
        &mut self,
        channel_index: u32,
        mut in_0: &[i16],
        in_len: &mut u32,
        mut out: &mut [i16],
        out_len: &mut u32,
    ) -> usize {
        if in_0.is_empty() {
            panic!("Empty input slice is not allowed");
        }
        let mut j: u32;
        let istride_save = self.in_stride;
        let ostride_save = self.out_stride;
        let mut ilen = *in_len;
        let mut olen = *out_len;
        let mem_idx = (channel_index * self.mem_alloc_size) as usize;
        let xlen: u32 = self.mem_alloc_size - self.filt_len - 1;
        let ylen: u32 = if olen < 8192 { olen } else { 8192 };
        let mut yselfack: Vec<f32> = vec![0.; ylen as usize];
        self.out_stride = 1;
        while 0 != ilen && 0 != olen {
            let mut ichunk: u32 = if ilen > xlen { xlen } else { ilen };
            let mut ochunk: u32 = if olen > ylen { ylen } else { olen };
            let mut omagic: u32 = 0;
            if self.magic_samples[channel_index as usize] != 0 {
                omagic = speex_resampler_magic(
                    self,
                    channel_index,
                    &mut yselfack.as_mut_slice(),
                    ochunk,
                ) as u32;
                ochunk -= omagic;
                olen -= omagic
            }
            if 0 == self.magic_samples[channel_index as usize] {
                j = 0u32;
                while j < ichunk {
                    self.mem[mem_idx
                        + j as usize
                        + (self.filt_len - 1) as usize] =
                        in_0[(j * istride_save) as usize] as f32;
                    j += 1
                }
                speex_resampler_process_native(
                    self,
                    channel_index,
                    &mut ichunk,
                    yselfack.as_mut_slice(),
                    &mut ochunk,
                );
            } else {
                ichunk = 0i32 as u32;
                ochunk = 0i32 as u32
            }
            j = 0u32;
            while j < ochunk + omagic {
                out[(j * ostride_save) as usize] =
                    (if yselfack[j as usize] < -32767.5 {
                        -32768
                    } else if yselfack[j as usize] > 32766.5 {
                        32767
                    } else {
                        (0.5 + yselfack[j as usize]).floor() as i32
                    }) as i16;
                j += 1
            }
            ilen -= ichunk;
            olen -= ochunk;
            out = &mut out[(ochunk + omagic * ostride_save as u32) as usize..];
            in_0 = &in_0[(ichunk * istride_save as u32) as usize..];
        }
        self.out_stride = ostride_save;
        *in_len -= ilen;
        *out_len -= olen;
        let resampler = self.resampler_ptr.unwrap();
        if resampler as usize == resampler_basic_zero as usize {
            RESAMPLER_ERR_ALLOC_FAILED
        } else {
            RESAMPLER_ERR_SUCCESS
        }
    }

    /* * Resample an interleaved float array. The input and output buffers muself *not* overlap.
     * @param self Resampler selfate
     * @param in Input buffer
     * @param in_len number of input samples in the input buffer. Returns the number
     * of samples processed. This is all per-channel.
     * @param out Output buffer
     * @param out_len Size of the output buffer. Returns the number of samples written.
     * This is all per-channel.
     */
    pub fn process_interleaved_float(
        &mut self,
        in_0: &[f32],
        in_len: &mut u32,
        out: &mut [f32],
        out_len: &mut u32,
    ) -> usize {
        if in_0.is_empty() {
            panic!("Empty input slice is not allowed");
        }
        let istride_save = self.in_stride;
        let ostride_save = self.out_stride;
        let bak_out_len = *out_len;
        let bak_in_len = *in_len;

        self.out_stride = self.nb_channels;
        self.in_stride = self.out_stride;
        for i in 0..self.nb_channels as usize {
            *out_len = bak_out_len;
            *in_len = bak_in_len;
            self.process_float(
                i as u32,
                &in_0[i..],
                in_len,
                &mut out[i..],
                out_len,
            );
        }
        self.in_stride = istride_save;
        self.out_stride = ostride_save;
        let resampler = self.resampler_ptr.unwrap();
        if resampler as usize == resampler_basic_zero as usize {
            RESAMPLER_ERR_ALLOC_FAILED
        } else {
            RESAMPLER_ERR_SUCCESS
        }
    }

    /* * Resample an interleaved int array. The input and output buffers muself *not* overlap.
     * @param self Resampler selfate
     * @param in Input buffer
     * @param in_len number of input samples in the input buffer. Returns the number
     * of samples processed. This is all per-channel.
     * @param out Output buffer
     * @param out_len Size of the output buffer. Returns the number of samples written.
     * This is all per-channel.
     */
    pub fn process_interleaved_int(
        &mut self,
        in_0: &[i16],
        in_len: &mut u32,
        out: &mut [i16],
        out_len: &mut u32,
    ) -> usize {
        if in_0.is_empty() {
            panic!("Empty input slice is not allowed");
        }
        let istride_save = self.in_stride;
        let ostride_save = self.out_stride;
        let bak_out_len = *out_len;
        let bak_in_len = *in_len;

        self.out_stride = self.nb_channels;
        self.in_stride = self.out_stride;
        for i in 0..self.nb_channels as usize {
            *out_len = bak_out_len;
            *in_len = bak_in_len;
            self.process_int(
                i as u32,
                &in_0[i..],
                in_len,
                &mut out[i..],
                out_len,
            );
        }
        self.in_stride = istride_save;
        self.out_stride = ostride_save;
        let resampler = self.resampler_ptr.unwrap();
        if resampler as usize == resampler_basic_zero as usize {
            RESAMPLER_ERR_ALLOC_FAILED
        } else {
            RESAMPLER_ERR_SUCCESS
        }
    }

    /* * Set (change) the conversion quality.
     * @param st Resampler state
     * @param quality Resampling quality between 0 and 10, where 0 has poor
     * quality and 10 has very high quality.
     */
    pub fn set_quality(&mut self, quality: usize) -> usize {
        if quality > 10 {
            RESAMPLER_ERR_INVALID_ARG
        } else if self.quality as usize == quality {
            RESAMPLER_ERR_SUCCESS
        } else {
            self.quality = quality as u32;
            if self.initialised != 0 {
                self.update_filter()
            } else {
                RESAMPLER_ERR_SUCCESS
            }
        }
    }

    /* * Set (change) the input/output sampling rates (integer value).
     * @param st Resampler state
     * @param in_rate Input sampling rate (integer number of Hz).
     * @param out_rate Output sampling rate (integer number of Hz).
     */
    pub fn set_rate(&mut self, in_rate: usize, out_rate: usize) -> usize {
        self.set_rate_frac(in_rate, out_rate, in_rate, out_rate)
    }

    /* * Set (change) the input/output sampling rates and resampling ratio
     * (fractional values in Hz supported).
     * @param st Resampler state
     * @param ratio_num numerator of the sampling rate ratio
     * @param ratio_den Denominator of the sampling rate ratio
     * @param in_rate Input sampling rate rounded to the nearest integer (in Hz).
     * @param out_rate Output sampling rate rounded to the nearest integer (in Hz).
     */
    pub fn set_rate_frac(
        &mut self,
        ratio_num: usize,
        ratio_den: usize,
        in_rate: usize,
        out_rate: usize,
    ) -> usize {
        if ratio_num == 0 || ratio_den == 0 {
            RESAMPLER_ERR_INVALID_ARG
        } else if self.in_rate == in_rate as u32
            && self.out_rate == out_rate as u32
            && self.num_rate == ratio_num as u32
            && self.den_rate == ratio_den as u32
        {
            RESAMPLER_ERR_SUCCESS
        } else {
            let old_den = self.den_rate;
            self.in_rate = in_rate as u32;
            self.out_rate = out_rate as u32;
            self.num_rate = ratio_num as u32;
            self.den_rate = ratio_den as u32;
            let fact = _gcd(self.num_rate, self.den_rate);
            self.num_rate = self.num_rate / fact;
            self.den_rate = self.den_rate / fact;
            if old_den > 0 {
                for val in &mut self.samp_frac_num {
                    let res = _muldiv(val, *val, self.den_rate, old_den);
                    if res != RESAMPLER_ERR_SUCCESS {
                        return RESAMPLER_ERR_OVERFLOW;
                    } else {
                        if *val >= self.den_rate {
                            *val = self.den_rate - 1;
                        }
                    }
                }
            }
            if self.initialised != 0 {
                self.update_filter()
            } else {
                RESAMPLER_ERR_SUCCESS
            }
        }
    }

    /* * Get the current input/output sampling rates (integer value).
     * @param st Resampler state
     */
    pub fn get_rate(&self) -> (usize, usize) {
        (self.in_rate as usize, self.out_rate as usize)
    }

    /* * Get the current resampling ratio. This will be reduced to the least
     * common denominator.
     * @param st Resampler state
     */
    pub fn get_ratio(&self) -> (usize, usize) {
        (self.num_rate as usize, self.den_rate as usize)
    }

    /* * Get the conversion quality.
     * @param st Resampler state
     * @return Resampling quality between 0 and 10, where 0 has poor
     * quality and 10 has very high quality.
     */
    pub fn get_quality(&self) -> usize {
        self.quality as usize
    }

    /* * Set (change) the input stride.
     * @param st Resampler state
     * @param stride Input stride
     */
    pub fn set_input_stride(&mut self, stride: usize) {
        self.in_stride = stride as u32;
    }

    /* * Get the input stride.
     * @param st Resampler state
     * @return Input stride copied
     */
    pub fn get_input_stride(&mut self) -> usize {
        self.in_stride as usize
    }

    /* * Set (change) the output stride.
     * @param st Resampler state
     * @param stride Output stride
     */
    pub fn set_output_stride(&mut self, stride: usize) {
        self.out_stride = stride as u32;
    }

    /* * Get the output stride.
     * @param st Resampler state copied
     * @return Output stride
     */
    pub fn get_output_stride(&mut self) -> usize {
        self.out_stride as usize
    }

    /* * Get the latency introduced by the resampler measured in input samples.
     * @param st Resampler state
     */
    pub fn get_input_latency(&self) -> usize {
        (self.filt_len / 2) as usize
    }

    /* * Get the latency introduced by the resampler measured in output samples.
     * @param st Resampler state
     */
    pub fn get_output_latency(&self) -> usize {
        (((self.filt_len / 2) * self.den_rate + (self.num_rate >> 1))
            / self.num_rate) as usize
    }

    /* * Make sure that the first samples to go out of the resamplers don't have
     * leading zeros. This is only useful before starting to use a newly created
     * resampler. It is recommended to use that when resampling an audio file, as
     * it will generate a file with the same length. For real-time processing,
     * it is probably easier not to use this call (so that the output duration
     * is the same for the first frame).
     * @param st Resampler state
     */
    pub fn skip_zeros(&mut self) {
        let filt_len = self.filt_len / 2;
        self.last_sample.iter_mut().for_each(|v: &mut u32| {
            *v = filt_len;
        });
    }

    /* * Reset a resampler so a new (unrelated) stream can be processed.
     * @param st Resampler state
     */
    pub fn reset_mem(&mut self) {
        self.last_sample.iter_mut().for_each(|elem| *elem = 0);
        self.magic_samples.iter_mut().for_each(|elem| *elem = 0);
        self.samp_frac_num.iter_mut().for_each(|elem| *elem = 0);

        self.mem.iter_mut().for_each(|elem| *elem = 0.);
    }

    #[inline]
    fn num_den(&mut self) {
        self.cutoff = QUALITY_MAP[self.quality as usize].downsample_bandwidth
            * self.den_rate as f32
            / self.num_rate as f32;
        let pass = self.filt_len;
        _muldiv(&mut self.filt_len, pass, self.num_rate, self.den_rate);
        self.filt_len = ((self.filt_len - 1) & (!7)) + 8;
        self.oversample = (1..5)
            .filter(|x| 2u32.pow(*x) * self.den_rate < self.num_rate)
            .fold(self.oversample, |acc, _| acc >> 1);

        if self.oversample < 1 {
            self.oversample = 1;
        }
    }

    #[inline]
    fn use_direct(&mut self) {
        let iter_chunk = self.sinc_table.chunks_mut(self.filt_len as usize);
        for (i, chunk) in iter_chunk.enumerate() {
            for (j, elem) in chunk.iter_mut().enumerate() {
                *elem = sinc(
                    self.cutoff,
                    (j as f32 - self.filt_len as f32 / 2.0 + 1.0)
                        - (i as f32) / self.den_rate as f32,
                    self.filt_len as i32,
                    QUALITY_MAP[self.quality as usize].window_func,
                );
            }
        }
        if self.quality > 8 {
            self.resampler_ptr = Some(resampler_basic_direct_double);
        } else {
            self.resampler_ptr = Some(resampler_basic_direct_single);
        }
    }

    #[inline]
    fn not_use_direct(&mut self) {
        let quality = self.quality as usize;
        let cutoff = self.cutoff;
        let oversample = self.oversample;
        let filt_len = self.filt_len;
        self.sinc_table
            .iter_mut()
            .enumerate()
            .take((oversample * filt_len + 8) as usize)
            .for_each(|(i, x)| {
                *x = sinc(
                    cutoff,
                    (i as i32 - 4) as f32 / oversample as f32
                        - filt_len as f32 / 2.0,
                    filt_len as i32,
                    QUALITY_MAP[quality].window_func,
                )
            });
        if self.quality > 8 {
            self.resampler_ptr = Some(resampler_basic_interpolate_double);
        } else {
            self.resampler_ptr = Some(resampler_basic_interpolate_single);
        }
    }

    #[inline(always)]
    fn chunks_iterator(
        &mut self,
        old_length: u32,
        alloc_size: usize,
        algo: usize,
    ) {
        let mem_copy = self.mem.clone();

        let mut_mem = self.mem.chunks_mut(self.mem_alloc_size as usize);
        let mem = mem_copy.chunks(alloc_size);

        for (ch_mut, ch) in mut_mem.zip(mem) {
            for magic in &mut self.magic_samples {
                if algo == 0 {
                    let range = old_length - 1 + *magic;
                    chunk_copy!(ch_mut, *magic, range, ch, 0, range);
                    chunk_assign!(ch_mut, 0, *magic, 0.0);
                } else if algo == 1 {
                    algo!(self, ch_mut, ch, old_length, *magic);
                } else {
                    let skip = (old_length - self.filt_len) / 2;
                    let ubound = self.filt_len - 1 + skip + *magic;
                    chunk_copy!(ch_mut, 0, ubound, ch, skip, ubound + skip);
                    *magic += skip;
                }
            }
        }
    }

    fn update_filter(&mut self) -> usize {
        let old_length = self.filt_len;
        let quality = self.quality as usize;
        let old_alloc_size = self.mem_alloc_size as usize;
        self.int_advance = self.num_rate / self.den_rate;
        self.frac_advance = self.num_rate % self.den_rate;
        self.oversample = QUALITY_MAP[quality].oversample as u32;
        self.filt_len = QUALITY_MAP[quality].base_length as u32;
        if self.num_rate > self.den_rate {
            self.num_den();
        } else {
            self.cutoff = QUALITY_MAP[quality].upsample_bandwidth;
        }

        let use_direct = self.filt_len * self.den_rate
            <= self.filt_len * self.oversample + 8
            && 2147483647 as u64
                / ::std::mem::size_of::<f32>() as u64
                / self.den_rate as u64
                >= self.filt_len as u64;

        let mut min_sinc_table_length = self.filt_len * self.den_rate;
        if !use_direct {
            min_sinc_table_length = self.filt_len * self.oversample + 8;
        }

        if self.sinc_table_length < min_sinc_table_length {
            self.sinc_table = vec![0.0; min_sinc_table_length as usize];
            self.sinc_table_length = min_sinc_table_length;
        }

        if use_direct {
            self.use_direct();
        } else {
            self.not_use_direct();
        }

        let min_alloc_size = self.filt_len - 1 + self.buffer_size;
        if min_alloc_size > self.mem_alloc_size {
            let mem = self.mem.clone();
            self.mem = vec![0.0; (self.nb_channels * min_alloc_size) as usize];
            self.mem[0..mem.len()].copy_from_slice(&mem);
            self.mem_alloc_size = min_alloc_size;
        }

        if self.started == 0 {
            let dim = (self.nb_channels * self.mem_alloc_size) as usize;
            self.mem = vec![0.0; dim];
        } else if self.filt_len > old_length {
            self.chunks_iterator(old_length, old_alloc_size, 0);
            self.chunks_iterator(old_length, self.mem_alloc_size as usize, 1);
        } else if self.filt_len < old_length {
            self.chunks_iterator(old_length, self.mem_alloc_size as usize, 2);
        }
        return RESAMPLER_ERR_SUCCESS;
    }
}

fn resampler_basic_zero(
    st: &mut SpeexResamplerState,
    channel_index: u32,
    _in_0: &[f32],
    in_len: &mut u32,
    out: &mut [f32],
    out_len: &mut u32,
) -> i32 {
    let mut out_sample: u32 = 0;
    let mut last_sample = st.last_sample[channel_index as usize];
    let mut samp_frac_num = st.samp_frac_num[channel_index as usize];
    let out_stride = st.out_stride;
    let int_advance = st.int_advance;
    let frac_advance = st.frac_advance;
    let den_rate: u32 = st.den_rate;
    while !(last_sample >= *in_len || out_sample >= *out_len) {
        out[(out_stride * out_sample) as usize] = 0.0;
        out_sample = out_sample + 1;
        last_sample += int_advance;
        samp_frac_num += frac_advance as u32;
        if samp_frac_num >= den_rate {
            samp_frac_num -= den_rate as u32;
            last_sample += 1
        }
    }
    st.last_sample[channel_index as usize] = last_sample;
    st.samp_frac_num[channel_index as usize] = samp_frac_num;
    out_sample as i32
}

fn resampler_basic_interpolate_single(
    st: &mut SpeexResamplerState,
    channel_index: u32,
    in_0: &[f32],
    in_len: &mut u32,
    out: &mut [f32],
    out_len: &mut u32,
) -> i32 {
    let n = st.filt_len as usize;
    let channel_idx = channel_index as usize;
    let mut last_sample = st.last_sample[channel_idx];
    let mut samp_frac_num = st.samp_frac_num[channel_idx];
    let out_stride = st.out_stride;
    let int_advance = st.int_advance;
    let frac_advance = st.frac_advance;
    let den_rate = st.den_rate;
    let oversample = st.oversample;
    let sinc_table = &st.sinc_table;

    let mut out_sample: u32 = 0;
    while !(last_sample >= *in_len || out_sample >= *out_len) {
        let iptr = &in_0[last_sample as usize..];
        let offset = samp_frac_num * oversample / den_rate;
        let frac =
            ((samp_frac_num * oversample) % den_rate) as f32 / den_rate as f32;
        let mut accum: [f32; 4] = [0.; 4];
        iptr.iter().zip(0..n).for_each(|(&curr_in, j)| {
            let idx = (2 + (j + 1) * oversample as usize) - offset as usize;
            accum.iter_mut().zip(sinc_table.iter().skip(idx)).for_each(
                |(v, &s)| {
                    *v += curr_in * s;
                },
            );
        });
        let mut interp: [f32; 4] = [0.; 4];
        cubic_coef(frac, &mut interp);
        out[(out_stride * out_sample) as usize] = interp
            .iter()
            .zip(accum.iter())
            .map(|(&x, &y)| x * y)
            .fold(0., |acc, x| acc + x);
        out_sample += 1;
        last_sample += int_advance;
        samp_frac_num += frac_advance as u32;
        if samp_frac_num >= den_rate {
            samp_frac_num -= den_rate;
            last_sample += 1;
        }
    }
    st.last_sample[channel_idx] = last_sample;
    st.samp_frac_num[channel_idx] = samp_frac_num;
    out_sample as i32
}

fn cubic_coef(frac: f32, interp: &mut [f32]) {
    interp[0] =
        -0.16666999459266663 * frac + 0.16666999459266663 * frac * frac * frac;
    interp[1] = frac + 0.5 * frac * frac - 0.5f32 * frac * frac * frac;
    interp[3] = -0.3333300054073334 * frac + 0.5 * frac * frac
        - 0.16666999459266663 * frac * frac * frac;
    interp[2] =
        (1.0f64 - interp[0] as f64 - interp[1] as f64 - interp[3] as f64)
            as f32;
}

fn resampler_basic_interpolate_double(
    st: &mut SpeexResamplerState,
    channel_index: u32,
    in_0: &[f32],
    in_len: &mut u32,
    out: &mut [f32],
    out_len: &mut u32,
) -> i32 {
    let n = st.filt_len as usize;
    let channel_idx = channel_index as usize;
    let mut last_sample = st.last_sample[channel_idx];
    let mut samp_frac_num = st.samp_frac_num[channel_idx];
    let out_stride = st.out_stride;
    let int_advance = st.int_advance;
    let oversample = st.oversample;
    let frac_advance = st.frac_advance;
    let den_rate = st.den_rate;
    let sinc_table = &st.sinc_table;

    let mut out_sample: u32 = 0;

    while !(last_sample >= *in_len || out_sample >= *out_len) {
        let iptr: &[f32] = &in_0[last_sample as usize..];
        let offset = samp_frac_num * st.oversample / st.den_rate;
        let frac =
            (samp_frac_num * oversample % den_rate) as f32 / den_rate as f32;
        let mut accum: [f64; 4] = [0.0; 4];
        iptr.iter().zip(0..n).for_each(|(&curr_in, j)| {
            let idx = (2 + (j + 1) * oversample as usize) - offset as usize;
            accum.iter_mut().zip(sinc_table.iter().skip(idx)).for_each(
                |(v, &s)| {
                    *v += (curr_in * s) as f64;
                },
            );
        });
        let mut interp: [f32; 4] = [0.; 4];
        cubic_coef(frac, &mut interp);
        out[(out_stride * out_sample) as usize] = interp
            .iter()
            .zip(accum.iter())
            .map(|(&x, &y)| x * y as f32)
            .fold(0., |acc, x| acc + x);
        out_sample = out_sample + 1;
        last_sample += int_advance;
        samp_frac_num += frac_advance;
        if samp_frac_num >= den_rate {
            samp_frac_num -= den_rate;
            last_sample += 1
        }
    }
    st.last_sample[channel_index as usize] = last_sample;
    st.samp_frac_num[channel_index as usize] = samp_frac_num;
    out_sample as i32
}

static QUALITY_MAP: [QualityMapping; 11] = [
    QualityMapping::new(
        8,
        4,
        0.8299999833106995,
        0.8600000143051148,
        &_KAISER6,
    ),
    QualityMapping::new(
        16,
        4,
        0.8500000238418579,
        0.8799999952316284,
        &_KAISER6,
    ),
    QualityMapping::new(
        32,
        4,
        0.8820000290870667,
        0.9100000262260437,
        &_KAISER6,
    ),
    QualityMapping::new(
        48,
        8,
        0.8949999809265137,
        0.9169999957084656,
        &_KAISER8,
    ),
    QualityMapping::new(
        64,
        8,
        0.9210000038146973,
        0.9399999976158142,
        &_KAISER8,
    ),
    QualityMapping::new(
        80,
        16,
        0.921999990940094,
        0.9399999976158142,
        &_KAISER10,
    ),
    QualityMapping::new(
        96,
        16,
        0.9399999976158142,
        0.9449999928474426,
        &_KAISER10,
    ),
    QualityMapping::new(
        128,
        16,
        0.949999988079071,
        0.949999988079071,
        &_KAISER10,
    ),
    QualityMapping::new(
        160,
        16,
        0.9599999785423279,
        0.9599999785423279,
        &_KAISER10,
    ),
    QualityMapping::new(
        192,
        32,
        0.9679999947547913,
        0.9679999947547913,
        &_KAISER12,
    ),
    QualityMapping::new(
        256,
        32,
        0.9750000238418579,
        0.9750000238418579,
        &_KAISER12,
    ),
];

static _KAISER12: FuncDef = FuncDef::new(&KAISER12_TABLE, 64);

static KAISER12_TABLE: [f64; 68] = {
    [
        0.99859849,
        1.0,
        0.99859849,
        0.99440475,
        0.98745105,
        0.97779076,
        0.9654977,
        0.95066529,
        0.93340547,
        0.91384741,
        0.89213598,
        0.86843014,
        0.84290116,
        0.81573067,
        0.78710866,
        0.75723148,
        0.7262997,
        0.69451601,
        0.66208321,
        0.62920216,
        0.59606986,
        0.56287762,
        0.52980938,
        0.49704014,
        0.46473455,
        0.43304576,
        0.40211431,
        0.37206735,
        0.343018,
        0.3150649,
        0.28829195,
        0.26276832,
        0.23854851,
        0.21567274,
        0.19416736,
        0.17404546,
        0.15530766,
        0.13794293999999999,
        0.12192957,
        0.10723616,
        0.09382272,
        0.08164178,
        0.0706395,
        0.06075685,
        0.05193064,
        0.04409466,
        0.03718069,
        0.03111947,
        0.02584161,
        0.02127838,
        0.0173625,
        0.01402878,
        0.01121463,
        0.00886058,
        0.00691064,
        0.00531256,
        0.00401805,
        0.00298291,
        0.00216702,
        0.00153438,
        0.00105297,
        0.00069463,
        0.00043489,
        0.00025272,
        0.00013031,
        0.0000527734,
        0.00001,
        0.0,
    ]
};

static _KAISER10: FuncDef = FuncDef::new(&KAISER10_TABLE, 32);

static KAISER10_TABLE: [f64; 36] = {
    [
        0.99537781, 1.0, 0.99537781, 0.98162644, 0.95908712, 0.92831446,
        0.89005583, 0.84522401, 0.79486424, 0.74011713, 0.68217934,
        0.62226347, 0.56155915, 0.5011968, 0.44221549, 0.38553619, 0.33194107,
        0.28205962, 0.23636152, 0.19515633, 0.15859932, 0.1267028, 0.09935205,
        0.07632451, 0.05731132, 0.0419398, 0.02979584, 0.0204451, 0.01345224,
        0.00839739, 0.00488951, 0.00257636, 0.00115101, 0.00035515, 0.0, 0.0,
    ]
};

static _KAISER8: FuncDef = FuncDef::new(&KAISER8_TABLE, 32);

static KAISER8_TABLE: [f64; 36] = {
    [
        0.99635258, 1.0, 0.99635258, 0.98548012, 0.96759014, 0.943022,
        0.91223751, 0.87580811, 0.83439927, 0.78875245, 0.73966538,
        0.68797126, 0.6345175, 0.58014482, 0.52566725, 0.47185369, 0.4194115,
        0.36897272, 0.32108304, 0.27619388, 0.23465776, 0.1967267, 0.1625538,
        0.13219758, 0.10562887, 0.08273982, 0.06335451, 0.04724088,
        0.03412321, 0.0236949, 0.01563093, 0.00959968, 0.00527363, 0.00233883,
        0.0005, 0.0,
    ]
};

static _KAISER6: FuncDef = FuncDef::new(&KAISER6_TABLE, 32);

static KAISER6_TABLE: [f64; 36] = {
    [
        0.99733006, 1.0, 0.99733006, 0.98935595, 0.97618418, 0.95799003,
        0.93501423, 0.90755855, 0.87598009, 0.84068475, 0.80211977,
        0.76076565, 0.71712752, 0.67172623, 0.62508937, 0.57774224,
        0.53019925, 0.48295561, 0.43647969, 0.39120616, 0.34752997,
        0.30580127, 0.26632152, 0.22934058, 0.19505503, 0.16360756,
        0.13508755, 0.10953262, 0.0869312, 0.067226, 0.0503182, 0.03607231,
        0.02432151, 0.01487334, 0.00752, 0.0,
    ]
};

fn sinc(cutoff: f32, x: f32, n: i32, window_func: &FuncDef) -> f32 {
    let xx = f64::from(x * cutoff);
    let x_abs = f64::from(x).abs();
    let n_64 = f64::from(n);
    let cutoff_64 = f64::from(cutoff);
    if x_abs < 0.000001 {
        cutoff
    } else if x_abs > 0.5 * n_64 {
        0.0
    } else {
        let first_factor = cutoff_64 * (PI_64 * xx).sin() / (PI_64 * xx);
        let second_factor = compute_func(
            (2.0 * f64::from(x) / n_64).abs() as f32,
            window_func,
        );
        (first_factor * second_factor) as f32
    }
}

fn compute_func(x: f32, func: &FuncDef) -> f64 {
    let mut interp: [f64; 4] = [0.0; 4];
    let y = x * func.oversample as f32;
    let ind = y.floor() as usize;
    let frac = f64::from(y - ind as f32);
    interp[3] = -0.1666666667 * frac + 0.1666666667 * frac.powi(3);
    interp[2] = frac + 0.5 * frac.powi(2) - 0.5 * frac.powi(3);
    interp[0] = -0.3333333333 * frac + 0.5 * frac.powi(2)
        - 0.1666666667 * frac.powi(3);
    interp[1] = 1.0 - interp[3] - interp[2] - interp[0];

    interp
        .iter()
        .zip(func.table.iter().skip(ind))
        .map(|(&x, &y)| x * y)
        .sum()
}

fn resampler_basic_direct_single(
    st: &mut SpeexResamplerState,
    channel_index: u32,
    in_0: &[f32],
    in_len: &mut u32,
    out: &mut [f32],
    out_len: &mut u32,
) -> i32 {
    let n: i32 = st.filt_len as i32;
    let mut out_sample: u32 = 0;
    let mut last_sample = st.last_sample[channel_index as usize];
    let mut samp_frac_num = st.samp_frac_num[channel_index as usize];
    let out_stride = st.out_stride;
    let int_advance = st.int_advance;
    let frac_advance = st.frac_advance;
    let den_rate: u32 = st.den_rate;
    while !(last_sample >= *in_len || out_sample >= *out_len) {
        let sinct: &[f32] =
            &st.sinc_table[(samp_frac_num * n as u32) as usize..];
        let iptr: &[f32] = &in_0[last_sample as usize..];
        let mut sum: f32 = 0.0;
        let mut j: i32 = 0;
        while j < n {
            sum += sinct[j as usize] * iptr[j as usize];
            j += 1
        }
        out[(out_stride * out_sample) as usize] = sum;
        out_sample += 1;
        last_sample += int_advance;
        samp_frac_num += frac_advance as u32;
        if samp_frac_num >= den_rate {
            samp_frac_num -= den_rate as u32;
            last_sample += 1
        }
    }
    st.last_sample[channel_index as usize] = last_sample;
    st.samp_frac_num[channel_index as usize] = samp_frac_num;
    out_sample as i32
}

fn resampler_basic_direct_double(
    st: &mut SpeexResamplerState,
    channel_index: u32,
    in_0: &[f32],
    in_len: &mut u32,
    out: &mut [f32],
    out_len: &mut u32,
) -> i32 {
    let n: i32 = st.filt_len as i32;
    let mut out_sample: u32 = 0;
    let mut last_sample = st.last_sample[channel_index as usize];
    let mut samp_frac_num = st.samp_frac_num[channel_index as usize];
    let out_stride = st.out_stride;
    let int_advance = st.int_advance;
    let frac_advance = st.frac_advance;
    let den_rate: u32 = st.den_rate;
    while !(last_sample >= *in_len || out_sample >= *out_len) {
        let sinct: &[f32] =
            &st.sinc_table[(samp_frac_num * n as u32) as usize..];
        let iptr: &[f32] = &in_0[last_sample as usize..];
        let mut accum: [f64; 4] = [0.0; 4];
        let mut j: i32 = 0;
        while j < n {
            accum[0usize] += f64::from(sinct[j as usize] * iptr[j as usize]);
            accum[1usize] +=
                f64::from(sinct[(j + 1) as usize] * iptr[(j + 1) as usize]);
            accum[2usize] +=
                f64::from(sinct[(j + 2) as usize] * iptr[(j + 2) as usize]);
            accum[3usize] +=
                f64::from(sinct[(j + 3) as usize] * iptr[(j + 3) as usize]);
            j += 4
        }
        let sum: f64 =
            accum[0usize] + accum[1usize] + accum[2usize] + accum[3usize];
        out[(out_stride * out_sample) as usize] = sum as f32;
        out_sample += 1;
        last_sample += int_advance;
        samp_frac_num += frac_advance as u32;
        if samp_frac_num >= den_rate {
            samp_frac_num -= den_rate as u32;
            last_sample += 1;
        }
    }
    st.last_sample[channel_index as usize] = last_sample as u32;
    st.samp_frac_num[channel_index as usize] = samp_frac_num;
    out_sample as i32
}

fn _muldiv(result: &mut u32, value: u32, mul: u32, div: u32) -> usize {
    let major: u32 = value / div;
    let remainder: u32 = value % div;
    if remainder > 4294967295 / mul
        || major > 4294967295 / mul
        || major * mul > 4294967295 - remainder * mul / div
    {
        RESAMPLER_ERR_OVERFLOW
    } else {
        *result = remainder * mul / div + major * mul;
        RESAMPLER_ERR_SUCCESS
    }
}

fn _gcd(mut a: u32, mut b: u32) -> u32 {
    while b != 0 {
        let temp = a;
        a = b;
        b = temp % b;
    }
    a
}

fn speex_resampler_process_native(
    st: &mut SpeexResamplerState,
    channel_index: u32,
    in_len: &mut u32,
    out: &mut [f32],
    out_len: &mut u32,
) -> usize {
    let n: i32 = st.filt_len as i32;
    let mem_idx = (channel_index * st.mem_alloc_size) as usize;
    st.started = 1;
    let mem = &st.mem.clone();
    let out_sample: i32 = st.resampler_ptr.expect("non-null function pointer")(
        st,
        channel_index,
        mem,
        in_len,
        out,
        out_len,
    );
    if st.last_sample[channel_index as usize] < *in_len {
        *in_len = st.last_sample[channel_index as usize] as u32;
    }
    *out_len = out_sample as u32;
    st.last_sample[channel_index as usize] -= *in_len;
    let ilen: u32 = *in_len;
    let mut j: i32 = 0;
    while j < n - 1 {
        st.mem[mem_idx + j as usize] =
            st.mem[mem_idx + (j as u32 + ilen) as usize];
        j += 1
    }
    RESAMPLER_ERR_SUCCESS
}

fn speex_resampler_magic<'a, 'b>(
    st: &mut SpeexResamplerState,
    channel_index: u32,
    out: &'a mut &'b mut [f32],
    mut out_len: u32,
) -> u32 {
    let channel_idx = channel_index as usize;
    let mut tmp_in_len = st.magic_samples[channel_idx];
    let mem_idx = (st.filt_len + channel_index * st.mem_alloc_size) as usize;
    speex_resampler_process_native(
        st,
        channel_index,
        &mut tmp_in_len,
        *out,
        &mut out_len,
    );
    st.magic_samples[channel_idx] -= tmp_in_len;
    if st.magic_samples[channel_idx] != 0 {
        let mem = &st.mem[mem_idx - 1 + tmp_in_len as usize..].to_vec();
        st.mem
            .iter_mut()
            .skip(mem_idx - 1)
            .zip(mem.iter())
            .take(st.magic_samples[channel_idx] as usize)
            .for_each(|(x, &y)| *x = y);
    }
    let value: &'b mut [f32] = mem::replace(out, &mut []);
    *out = &mut value[(out_len * st.out_stride as u32) as usize..];
    out_len
}
