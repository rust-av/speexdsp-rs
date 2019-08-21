use std::mem;

pub const RESAMPLER_ERR_SUCCESS: u32 = 0;
pub const RESAMPLER_ERR_ALLOC_FAILED: u32 = 1;
pub const RESAMPLER_ERR_INVALID_ARG: u32 = 3;
pub const RESAMPLER_ERR_OVERFLOW: u32 = 5;

#[derive(Clone)]
pub struct SpeexResamplerState {
    in_rate: u32,
    out_rate: u32,
    num_rate: u32,
    den_rate: u32,
    quality: i32,
    nb_channels: u32,
    filt_len: u32,
    mem_alloc_size: u32,
    buffer_size: u32,
    int_advance: i32,
    frac_advance: i32,
    cutoff: f32,
    oversample: u32,
    initialised: i32,
    started: i32,
    last_sample: Vec<i32>,
    samp_frac_num: Vec<u32>,
    magic_samples: Vec<u32>,
    mem: Vec<f32>,
    sinc_table: Vec<f32>,
    sinc_table_length: u32,
    resampler_ptr: ResamplerBasicFunc,
    in_stride: i32,
    out_stride: i32,
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
            quality: -1,
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
        if filter_err == RESAMPLER_ERR_SUCCESS as i32 {
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
    ) -> i32 {
        let mut j: i32;
        let mut ilen: u32 = *in_len;
        let mut olen: u32 = *out_len;
        let mem_idx = (channel_index * self.mem_alloc_size) as usize;
        let filt_offs: i32 = (self.filt_len - 1) as i32;
        let xlen: u32 = self.mem_alloc_size - filt_offs as u32;
        let istride: i32 = self.in_stride;
        if 0 != self.magic_samples[channel_index as usize] {
            olen -= speex_resampler_magic(self, channel_index, &mut out, olen)
                as u32;
        }
        if 0 == self.magic_samples[channel_index as usize] {
            while 0 != ilen && 0 != olen {
                let mut ichunk: u32 = if ilen > xlen { xlen } else { ilen };
                let mut ochunk: u32 = olen;
                if !in_0.is_empty() {
                    j = 0i32;
                    while (j as u32) < ichunk {
                        self.mem[mem_idx + (j + filt_offs) as usize] =
                            in_0[((j * istride) as usize)];
                        j += 1
                    }
                } else {
                    j = 0i32;
                    while (j as u32) < ichunk {
                        self.mem[mem_idx + (j + filt_offs) as usize] = 0.0;
                        j += 1
                    }
                }
                speex_resampler_process_native(
                    self,
                    channel_index,
                    &mut ichunk,
                    out,
                    &mut ochunk,
                );
                ilen -= ichunk;
                olen -= ochunk;
                out = &mut out[(ochunk * self.out_stride as u32) as usize..];
                if !in_0.is_empty() {
                    in_0 = &in_0[(ichunk * istride as u32) as usize..];
                }
            }
        }
        *in_len -= ilen;
        *out_len -= olen;
        let resampler = self.resampler_ptr.unwrap();
        if resampler as usize == resampler_basic_zero as usize {
            RESAMPLER_ERR_ALLOC_FAILED as i32
        } else {
            RESAMPLER_ERR_SUCCESS as i32
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
    ) -> i32 {
        let mut j: i32;
        let istride_save: i32 = self.in_stride;
        let ostride_save: i32 = self.out_stride;
        let mut ilen: u32 = *in_len;
        let mut olen: u32 = *out_len;
        let mem_idx = (channel_index * self.mem_alloc_size) as usize;
        let xlen: u32 = self.mem_alloc_size - self.filt_len - 1;
        let ylen: u32 = if olen < 8192 { olen } else { 8192 };
        let mut yselfack: Vec<f32> = vec![0.; ylen as usize];
        self.out_stride = 1;
        while 0 != ilen && 0 != olen {
            let mut ichunk: u32 = if ilen > xlen { xlen } else { ilen };
            let mut ochunk: u32 = if olen > ylen { ylen } else { olen };
            let mut omagic: u32 = 0;
            if 0 != self.magic_samples[channel_index as usize] {
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
                if !in_0.is_empty() {
                    j = 0i32;
                    while (j as u32) < ichunk {
                        self.mem[mem_idx
                            + j as usize
                            + (self.filt_len - 1) as usize] =
                            in_0[(j * istride_save) as usize] as f32;
                        j += 1
                    }
                } else {
                    j = 0i32;
                    while (j as u32) < ichunk {
                        self.mem[mem_idx
                            + j as usize
                            + (self.filt_len - 1) as usize] = 0.0;
                        j += 1
                    }
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
            j = 0i32;
            while (j as u32) < ochunk + omagic {
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
            if !in_0.is_empty() {
                in_0 = &in_0[(ichunk * istride_save as u32) as usize..];
            }
        }
        self.out_stride = ostride_save;
        *in_len -= ilen;
        *out_len -= olen;
        let resampler = self.resampler_ptr.unwrap();
        if resampler as usize == resampler_basic_zero as usize {
            RESAMPLER_ERR_ALLOC_FAILED as i32
        } else {
            RESAMPLER_ERR_SUCCESS as i32
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
    ) -> i32 {
        let istride_save: i32 = self.in_stride;
        let ostride_save: i32 = self.out_stride;
        let bak_out_len: u32 = *out_len;
        let bak_in_len: u32 = *in_len;

        self.out_stride = self.nb_channels as i32;
        self.in_stride = self.out_stride;
        let mut i = 0;
        while i < self.nb_channels as usize {
            *out_len = bak_out_len;
            *in_len = bak_in_len;
            if !in_0.is_empty() {
                self.process_float(
                    i as u32,
                    &in_0[i..],
                    in_len,
                    &mut out[i..],
                    out_len,
                );
            } else {
                self.process_float(
                    i as u32,
                    &[],
                    in_len,
                    &mut out[i..],
                    out_len,
                );
            }
            i += 1;
        }
        self.in_stride = istride_save;
        self.out_stride = ostride_save;
        let resampler = self.resampler_ptr.unwrap();
        if resampler as usize == resampler_basic_zero as usize {
            RESAMPLER_ERR_ALLOC_FAILED as i32
        } else {
            RESAMPLER_ERR_SUCCESS as i32
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
    ) -> i32 {
        let istride_save = self.in_stride;
        let ostride_save = self.out_stride;
        let bak_out_len = *out_len;
        let bak_in_len = *in_len;

        self.out_stride = self.nb_channels as i32;
        self.in_stride = self.out_stride;
        let mut i = 0;
        while i < self.nb_channels as usize {
            *out_len = bak_out_len;
            *in_len = bak_in_len;
            if in_0.is_empty() {
                self.process_int(
                    i as u32,
                    &in_0[i..],
                    in_len,
                    &mut out[i..],
                    out_len,
                );
            } else {
                self.process_int(
                    i as u32,
                    &[],
                    in_len,
                    &mut out[i..],
                    out_len,
                );
            }
            i += 1;
        }
        self.in_stride = istride_save;
        self.out_stride = ostride_save;
        let resampler = self.resampler_ptr.unwrap();
        if resampler as usize == resampler_basic_zero as usize {
            RESAMPLER_ERR_ALLOC_FAILED as i32
        } else {
            RESAMPLER_ERR_SUCCESS as i32
        }
    }

    /* * Set (change) the conversion quality.
     * @param st Resampler state
     * @param quality Resampling quality between 0 and 10, where 0 has poor
     * quality and 10 has very high quality.
     */
    pub fn set_quality(&mut self, quality: usize) -> i32 {
        if quality > 10 {
            RESAMPLER_ERR_INVALID_ARG as i32
        } else if self.quality as usize == quality {
            RESAMPLER_ERR_SUCCESS as i32
        } else {
            self.quality = quality as i32;
            if 0 != self.initialised {
                self.update_filter()
            } else {
                RESAMPLER_ERR_SUCCESS as i32
            }
        }
    }

    /* * Set (change) the input/output sampling rates (integer value).
     * @param st Resampler state
     * @param in_rate Input sampling rate (integer number of Hz).
     * @param out_rate Output sampling rate (integer number of Hz).
     */
    pub fn set_rate(&mut self, in_rate: usize, out_rate: usize) -> i32 {
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
    ) -> i32 {
        if ratio_num == 0 || ratio_den == 0 {
            RESAMPLER_ERR_INVALID_ARG as i32
        } else if self.in_rate == in_rate as u32
            && self.out_rate == out_rate as u32
            && self.num_rate == ratio_num as u32
            && self.den_rate == ratio_den as u32
        {
            RESAMPLER_ERR_SUCCESS as i32
        } else {
            let old_den = self.den_rate;
            self.in_rate = in_rate as u32;
            self.out_rate = out_rate as u32;
            self.num_rate = ratio_num as u32;
            self.den_rate = ratio_den as u32;
            let fact: u32 = _gcd(self.num_rate, self.den_rate);
            self.num_rate = self.num_rate / fact;
            self.den_rate = self.den_rate / fact;
            if old_den > 0 {
                let mut i = 0;
                while i < self.nb_channels {
                    let pass = self.samp_frac_num[i as usize];
                    let res = _muldiv(
                        &mut self.samp_frac_num[i as usize],
                        pass,
                        self.den_rate,
                        old_den,
                    );
                    if res != RESAMPLER_ERR_SUCCESS as i32 {
                        return RESAMPLER_ERR_OVERFLOW as i32;
                    } else {
                        if self.samp_frac_num[i as usize] >= self.den_rate {
                            self.samp_frac_num[i as usize] = self.den_rate - 1;
                        }
                        i += 1;
                    }
                }
            }
            if 0 != self.initialised {
                self.update_filter()
            } else {
                RESAMPLER_ERR_SUCCESS as i32
            }
        }
    }

    /* * Get the current input/output sampling rates (integer value).
     * @param st Resampler state
     */
    pub fn get_rate(&mut self) -> (usize, usize) {
        (self.in_rate as usize, self.out_rate as usize)
    }

    /* * Get the current resampling ratio. This will be reduced to the least
     * common denominator.
     * @param st Resampler state
     */
    pub fn get_ratio(&mut self) -> (usize, usize) {
        (self.num_rate as usize, self.den_rate as usize)
    }

    /* * Get the conversion quality.
     * @param st Resampler state
     * @return Resampling quality between 0 and 10, where 0 has poor
     * quality and 10 has very high quality.
     */
    pub fn get_quality(&mut self) -> usize {
        self.quality as usize
    }

    /* * Set (change) the input stride.
     * @param st Resampler state
     * @param stride Input stride
     */
    pub fn set_input_stride(&mut self, stride: u32) {
        self.in_stride = stride as i32;
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
    pub fn set_output_stride(&mut self, stride: u32) {
        self.out_stride = stride as i32;
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
    pub fn get_input_latency(&mut self) -> i32 {
        (self.filt_len / 2) as i32
    }

    /* * Get the latency introduced by the resampler measured in output samples.
     * @param st Resampler state
     */
    pub fn get_output_latency(&mut self) -> i32 {
        (((self.filt_len / 2) * self.den_rate + (self.num_rate >> 1))
            / self.num_rate) as i32
    }

    /* * Make sure that the first samples to go out of the resamplers don't have
     * leading zeros. This is only useful before starting to use a newly created
     * resampler. It is recommended to use that when resampling an audio file, as
     * it will generate a file with the same length. For real-time processing,
     * it is probably easier not to use this call (so that the output duration
     * is the same for the first frame).
     * @param st Resampler state
     */
    pub fn skip_zeros(&mut self) -> i32 {
        let filt_len = (self.filt_len / 2) as i32;
        self.last_sample.iter_mut().for_each(|v| {
            *v = filt_len;
        });

        RESAMPLER_ERR_SUCCESS as i32
    }

    /* * Reset a resampler so a new (unrelated) stream can be processed.
     * @param st Resampler state
     */
    pub fn reset_mem(&mut self) -> i32 {
        let nb_channels = self.nb_channels as usize;
        self.last_sample = vec![0; nb_channels];
        self.magic_samples = vec![0; nb_channels];
        self.samp_frac_num = vec![0; nb_channels];

        self.mem = vec![0.0; nb_channels * (self.filt_len - 1) as usize];
        RESAMPLER_ERR_SUCCESS as i32
    }

    fn update_filter(&mut self) -> i32 {
        let mut current_block: u64;
        let old_length = self.filt_len;
        let old_alloc_size = self.mem_alloc_size;
        let mut min_sinc_table_length = 0;
        self.int_advance = (self.num_rate / self.den_rate) as i32;
        self.frac_advance = (self.num_rate % self.den_rate) as i32;
        self.oversample = QUALITY_MAP[self.quality as usize].oversample as u32;
        self.filt_len = QUALITY_MAP[self.quality as usize].base_length as u32;
        if self.num_rate > self.den_rate {
            self.cutoff = QUALITY_MAP[self.quality as usize]
                .downsample_bandwidth
                * self.den_rate as f32
                / self.num_rate as f32;
            let pass = self.filt_len;
            let ret = _muldiv(
                &mut self.filt_len,
                pass,
                self.num_rate,
                self.den_rate,
            );
            if ret != RESAMPLER_ERR_SUCCESS as i32 {
                current_block = 13465332693182510667;
            } else {
                self.filt_len = (self.filt_len.wrapping_sub(1i32 as u32)
                    & !7i32 as u32)
                    .wrapping_add(8i32 as u32);
                if (2i32 as u32).wrapping_mul(self.den_rate) < self.num_rate {
                    self.oversample >>= 1i32
                }
                if (4i32 as u32).wrapping_mul(self.den_rate) < self.num_rate {
                    self.oversample >>= 1i32
                }
                if (8i32 as u32).wrapping_mul(self.den_rate) < self.num_rate {
                    self.oversample >>= 1i32
                }
                if (16i32 as u32).wrapping_mul(self.den_rate) < self.num_rate {
                    self.oversample >>= 1i32
                }
                if self.oversample < 1i32 as u32 {
                    self.oversample = 1i32 as u32;
                    current_block = 8258075665625361029;
                } else {
                    current_block = 8258075665625361029;
                }
            }
        } else {
            self.cutoff =
                QUALITY_MAP[self.quality as usize].upsample_bandwidth;
            current_block = 8258075665625361029;
        }
        match current_block {
            8258075665625361029 => {
                let use_direct: i32 =
                    (self.filt_len.wrapping_mul(self.den_rate)
                        <= self
                            .filt_len
                            .wrapping_mul(self.oversample)
                            .wrapping_add(8i32 as u32)
                        && (2147483647i32 as u64)
                            .wrapping_div(::std::mem::size_of::<f32>() as u64)
                            .wrapping_div(self.den_rate as u64)
                            >= self.filt_len as u64)
                        as i32;
                if 0 != use_direct {
                    min_sinc_table_length =
                        self.filt_len.wrapping_mul(self.den_rate);
                    current_block = 14523784380283086299;
                } else if (2147483647i32 as u64)
                    .wrapping_div(::std::mem::size_of::<f32>() as u64)
                    .wrapping_sub(8i32 as u64)
                    .wrapping_div(self.oversample as u64)
                    < self.filt_len as u64
                {
                    current_block = 13465332693182510667;
                } else {
                    min_sinc_table_length = self
                        .filt_len
                        .wrapping_mul(self.oversample)
                        .wrapping_add(8i32 as u32);
                    current_block = 14523784380283086299;
                }
                match current_block {
                    13465332693182510667 => {}
                    _ => {
                        if self.sinc_table_length < min_sinc_table_length {
                            self.sinc_table =
                                vec![0.0; min_sinc_table_length as usize];
                            self.sinc_table_length = min_sinc_table_length;
                        }
                        current_block = 11650488183268122163;
                        match current_block {
                            13465332693182510667 => {}
                            _ => {
                                if 0 != use_direct {
                                    let mut i: u32 = 0;
                                    while i < self.den_rate {
                                        let mut j: i32 = 0;
                                        while (j as u32) < self.filt_len {
                                            self.sinc_table[i
                                                .wrapping_mul(self.filt_len)
                                                .wrapping_add(j as u32)
                                                as usize] = sinc(
                                                self.cutoff,
                                                (j - self.filt_len as i32
                                                    / 2i32
                                                    + 1i32)
                                                    as f32
                                                    - i as f32
                                                        / self.den_rate as f32,
                                                self.filt_len as i32,
                                                QUALITY_MAP
                                                    [self.quality as usize]
                                                    .window_func,
                                            );
                                            j += 1
                                        }
                                        i = i.wrapping_add(1)
                                    }
                                    if self.quality > 8i32 {
                                        self.resampler_ptr =
                                            Some(resampler_basic_direct_double)
                                    } else {
                                        self.resampler_ptr =
                                            Some(resampler_basic_direct_single)
                                    }
                                } else {
                                    let mut i_0: i32 = -4;
                                    while i_0
                                        < self
                                            .oversample
                                            .wrapping_mul(self.filt_len)
                                            .wrapping_add(4i32 as u32)
                                            as i32
                                    {
                                        self.sinc_table
                                            [(i_0 + 4i32) as usize] = sinc(
                                            self.cutoff,
                                            i_0 as f32
                                                / self.oversample as f32
                                                - self
                                                    .filt_len
                                                    .wrapping_div(2i32 as u32)
                                                    as f32,
                                            self.filt_len as i32,
                                            QUALITY_MAP[self.quality as usize]
                                                .window_func,
                                        );
                                        i_0 += 1
                                    }
                                    if self.quality > 8i32 {
                                        self.resampler_ptr = Some(
                                            resampler_basic_interpolate_double,
                                        )
                                    } else {
                                        self.resampler_ptr = Some(
                                            resampler_basic_interpolate_single,
                                        )
                                    }
                                }
                                let min_alloc_size: u32 = self
                                    .filt_len
                                    .wrapping_sub(1i32 as u32)
                                    .wrapping_add(self.buffer_size);
                                if min_alloc_size > self.mem_alloc_size {
                                    if (2147483647i32 as u64)
                                        .wrapping_div(
                                            ::std::mem::size_of::<f32>()
                                                as u64,
                                        )
                                        .wrapping_div(self.nb_channels as u64)
                                        < min_alloc_size as u64
                                    {
                                        current_block = 13465332693182510667;
                                    } else {
                                        self.mem = Vec::with_capacity(
                                            self.nb_channels
                                                .wrapping_mul(min_alloc_size)
                                                as usize,
                                        );
                                        self.mem_alloc_size = min_alloc_size;
                                        current_block = 15089075282327824602;
                                    }
                                } else {
                                    current_block = 15089075282327824602;
                                }
                                match current_block {
                                    13465332693182510667 => {}
                                    _ => {
                                        if 0 == self.started {
                                            let dim: usize =
                                                self.nb_channels.wrapping_mul(
                                                    self.mem_alloc_size,
                                                )
                                                    as usize;
                                            for _ in 0..dim {
                                                self.mem.push(0.0);
                                            }
                                        } else if self.filt_len > old_length {
                                            let mut i_2: u32 =
                                                self.nb_channels;
                                            loop {
                                                let fresh0 = i_2;
                                                i_2 = i_2.wrapping_sub(1);
                                                if !(0 != fresh0) {
                                                    break;
                                                }
                                                let olen: u32 = old_length
                                                    .wrapping_add(
                                                        (2i32 as u32)
                                                            .wrapping_mul(
                                                            self.magic_samples
                                                                [i_2 as usize],
                                                        ),
                                                    );
                                                let mut j_0: u32 = old_length
                                                    .wrapping_sub(1i32 as u32)
                                                    .wrapping_add(
                                                        self.magic_samples
                                                            [i_2 as usize],
                                                    );
                                                loop {
                                                    let fresh1 = j_0;
                                                    j_0 = j_0.wrapping_sub(1);
                                                    if !(0 != fresh1) {
                                                        break;
                                                    }
                                                    println!(
                                                    "{}",
                                                    i_2.wrapping_mul(
                                                        self.mem_alloc_size
                                                    )
                                                    .wrapping_add(j_0)
                                                    .wrapping_add(
                                                        self.magic_samples
                                                            [i_2 as usize]
                                                    )
                                                );
                                                    self.mem[i_2
                                                    .wrapping_mul(
                                                        self.mem_alloc_size,
                                                    )
                                                    .wrapping_add(j_0)
                                                    .wrapping_add(
                                                        self.magic_samples
                                                            [i_2 as usize],
                                                    )
                                                    as usize] = self.mem[i_2
                                                    .wrapping_mul(
                                                        old_alloc_size,
                                                    )
                                                    .wrapping_add(j_0)
                                                    as usize];
                                                    println!("here");
                                                }
                                                j_0 = 0i32 as u32;
                                                while j_0
                                                    < self.magic_samples
                                                        [i_2 as usize]
                                                {
                                                    self.mem[i_2
                                                    .wrapping_mul(
                                                        self.mem_alloc_size,
                                                    )
                                                    .wrapping_add(j_0)
                                                    as usize] = 0i32 as f32;
                                                    j_0 = j_0.wrapping_add(1)
                                                }
                                                self.magic_samples
                                                    [i_2 as usize] =
                                                    0i32 as u32;
                                                if self.filt_len > olen {
                                                    j_0 = 0i32 as u32;
                                                    while j_0
                                                        < olen.wrapping_sub(
                                                            1i32 as u32,
                                                        )
                                                    {
                                                        self.mem[i_2
                                                        .wrapping_mul(
                                                            self.mem_alloc_size,
                                                        )
                                                        .wrapping_add(
                                                            self.filt_len
                                                                .wrapping_sub(
                                                                    2i32
                                                                        as u32,
                                                                )
                                                                .wrapping_sub(
                                                                    j_0,
                                                                ),
                                                        )
                                                        as usize] = self.mem[i_2
                                                        .wrapping_mul(
                                                            self.mem_alloc_size,
                                                        )
                                                        .wrapping_add(
                                                            olen.wrapping_sub(
                                                                2i32 as u32,
                                                            )
                                                            .wrapping_sub(j_0),
                                                        )
                                                        as usize];
                                                        j_0 =
                                                            j_0.wrapping_add(1)
                                                    }
                                                    while j_0
                                                        < self
                                                            .filt_len
                                                            .wrapping_sub(
                                                                1i32 as u32,
                                                            )
                                                    {
                                                        self.mem[i_2
                                                        .wrapping_mul(
                                                            self.mem_alloc_size,
                                                        )
                                                        .wrapping_add(
                                                            self.filt_len
                                                                .wrapping_sub(
                                                                    2i32
                                                                        as u32,
                                                                )
                                                                .wrapping_sub(
                                                                    j_0,
                                                                ),
                                                        )
                                                        as usize] =
                                                        0i32 as f32;
                                                        j_0 =
                                                            j_0.wrapping_add(1)
                                                    }
                                                    let ref mut fresh2 = self
                                                        .last_sample
                                                        [i_2 as usize];
                                                    *fresh2 = (*fresh2 as u32)
                                                        .wrapping_add(
                                                            self.filt_len
                                                                .wrapping_sub(
                                                                    olen,
                                                                )
                                                                .wrapping_div(
                                                                    2i32
                                                                        as u32,
                                                                ),
                                                        )
                                                        as i32
                                                        as i32
                                                } else {
                                                    self.magic_samples
                                                        [i_2 as usize] = olen
                                                        .wrapping_sub(
                                                            self.filt_len,
                                                        )
                                                        .wrapping_div(
                                                            2i32 as u32,
                                                        );
                                                    j_0 = 0i32 as u32;
                                                    while j_0 < self
                                                        .filt_len
                                                        .wrapping_sub(
                                                            1i32 as u32,
                                                        )
                                                        .wrapping_add(
                                                            self.magic_samples
                                                                [i_2 as usize],
                                                        )
                                                    {
                                                        self.mem[i_2
                                                        .wrapping_mul(
                                                            self.mem_alloc_size,
                                                        )
                                                        .wrapping_add(j_0)
                                                        as usize] = self.mem[i_2
                                                        .wrapping_mul(
                                                            self.mem_alloc_size,
                                                        )
                                                        .wrapping_add(j_0)
                                                        .wrapping_add(
                                                            self.magic_samples
                                                                [i_2 as usize],
                                                        )
                                                        as usize];
                                                        j_0 =
                                                            j_0.wrapping_add(1)
                                                    }
                                                }
                                            }
                                        } else if self.filt_len < old_length {
                                            let mut i_3: u32 = 0;
                                            while i_3 < self.nb_channels {
                                                let old_magic: u32 = self
                                                    .magic_samples
                                                    [i_3 as usize];
                                                self.magic_samples
                                                    [i_3 as usize] =
                                                    old_length
                                                        .wrapping_sub(
                                                            self.filt_len,
                                                        )
                                                        .wrapping_div(
                                                            2i32 as u32,
                                                        );
                                                let mut j_1: u32 = 0;
                                                while j_1
                                                    < self
                                                        .filt_len
                                                        .wrapping_sub(
                                                            1i32 as u32,
                                                        )
                                                        .wrapping_add(
                                                            self.magic_samples
                                                                [i_3 as usize],
                                                        )
                                                        .wrapping_add(
                                                            old_magic,
                                                        )
                                                {
                                                    self.mem[i_3
                                                    .wrapping_mul(
                                                        self.mem_alloc_size,
                                                    )
                                                    .wrapping_add(j_1)
                                                    as usize] = self.mem[i_3
                                                    .wrapping_mul(
                                                        self.mem_alloc_size,
                                                    )
                                                    .wrapping_add(j_1)
                                                    .wrapping_add(
                                                        self.magic_samples
                                                            [i_3 as usize],
                                                    )
                                                    as usize];
                                                    j_1 = j_1.wrapping_add(1)
                                                }
                                                let ref mut fresh3 = self
                                                    .magic_samples
                                                    [i_3 as usize];
                                                *fresh3 = (*fresh3 as u32)
                                                    .wrapping_add(old_magic)
                                                    as u32
                                                    as u32;
                                                i_3 = i_3.wrapping_add(1)
                                            }
                                        }
                                        return RESAMPLER_ERR_SUCCESS as i32;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            _ => {}
        }
        self.resampler_ptr = Some(resampler_basic_zero);
        self.filt_len = old_length;
        RESAMPLER_ERR_ALLOC_FAILED as i32
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
    let mut out_sample: i32 = 0;
    let mut last_sample: i32 = st.last_sample[channel_index as usize];
    let mut samp_frac_num: u32 = st.samp_frac_num[channel_index as usize];
    let out_stride: i32 = st.out_stride;
    let int_advance: i32 = st.int_advance;
    let frac_advance: i32 = st.frac_advance;
    let den_rate: u32 = st.den_rate;
    while !(last_sample >= *in_len as i32 || out_sample >= *out_len as i32) {
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
    out_sample
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
    let mut last_sample = st.last_sample[channel_index as usize];
    let mut samp_frac_num = st.samp_frac_num[channel_index as usize];
    let out_stride = st.out_stride;
    let int_advance = st.int_advance;
    let frac_advance = st.frac_advance;
    let den_rate = st.den_rate;
    let oversample = st.oversample;
    let sinc_table = &st.sinc_table;

    let mut out_sample = 0;
    while !(last_sample >= *in_len as i32 || out_sample >= *out_len as i32) {
        let iptr = &in_0[last_sample as usize..];
        let offset = samp_frac_num * oversample / den_rate;
        let frac: f32 =
            ((samp_frac_num * oversample) % den_rate) as f32 / den_rate as f32;
        let mut accum: [f32; 4] = [0.; 4];
        iptr.iter().zip(0..n).for_each(|(&curr_in, j)| {
            let idx = (2 + ((j as u32 + 1) * oversample) - offset) as usize;
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
    st.last_sample[channel_index as usize] = last_sample;
    st.samp_frac_num[channel_index as usize] = samp_frac_num;
    out_sample
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
    let n: i32 = st.filt_len as i32;
    let mut out_sample: i32 = 0;
    let mut last_sample: i32 = st.last_sample[channel_index as usize];
    let mut samp_frac_num: u32 = st.samp_frac_num[channel_index as usize];
    let out_stride: i32 = st.out_stride;
    let int_advance: i32 = st.int_advance;
    let frac_advance: i32 = st.frac_advance;
    let den_rate: u32 = st.den_rate;
    while !(last_sample >= *in_len as i32 || out_sample >= *out_len as i32) {
        let iptr: &[f32] = &in_0[last_sample as usize..];
        let offset: i32 = (samp_frac_num * st.oversample / st.den_rate) as i32;
        let frac: f32 = (samp_frac_num * st.oversample % st.den_rate) as f32
            / st.den_rate as f32;
        let mut interp: [f32; 4] = [0.; 4];
        let mut accum: [f64; 4] = [0.0; 4];
        let mut j: i32 = 0;
        while j < n {
            let curr_in: f64 = iptr[j as usize] as f64;
            accum[0usize] += curr_in
                * st.sinc_table[4
                    + ((j + 1) as usize * st.oversample as usize)
                    - offset as usize
                    - 2 as usize] as f64;
            accum[1usize] += curr_in
                * st.sinc_table[4
                    + ((j + 1) as usize * st.oversample as usize)
                    - offset as usize
                    - 1 as usize] as f64;
            accum[2usize] += curr_in
                * st.sinc_table[4
                    + ((j + 1) as usize * st.oversample as usize)
                    - offset as usize as usize] as f64;
            accum[3usize] += curr_in
                * st.sinc_table[4
                    + ((j + 1) as usize * st.oversample as usize)
                    - offset as usize
                    + 1 as usize] as f64;
            j += 1
        }
        cubic_coef(frac, &mut interp);
        let sum: f32 = (interp[0usize] as f64 * accum[0usize]
            + interp[1usize] as f64 * accum[1usize]
            + interp[2usize] as f64 * accum[2usize]
            + interp[3usize] as f64 * accum[3usize])
            as f32;
        out[(out_stride * out_sample) as usize] = sum;
        out_sample = out_sample + 1;
        last_sample += int_advance;
        samp_frac_num += frac_advance as u32;
        if samp_frac_num >= den_rate {
            samp_frac_num -= den_rate;
            last_sample += 1
        }
    }
    st.last_sample[channel_index as usize] = last_sample;
    st.samp_frac_num[channel_index as usize] = samp_frac_num;
    out_sample
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

static _KAISER10: FuncDef = FuncDef::new(&KAISER10_TABLE, 32);

static KAISER10_TABLE: [f64; 36] = {
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

static _KAISER8: FuncDef = FuncDef::new(&KAISER8_TABLE, 32);

static KAISER8_TABLE: [f64; 36] = {
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

static _KAISER6: FuncDef = FuncDef::new(&KAISER6_TABLE, 32);

static KAISER6_TABLE: [f64; 36] = {
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

fn sinc(cutoff: f32, x: f32, n: i32, window_func: &FuncDef) -> f32 {
    let xx: f32 = x * cutoff;
    if (f64::from(x)).abs() < 0.000001f64 {
        cutoff
    } else if (f64::from(x)).abs() > 0.5f64 * f64::from(n) {
        0.0
    } else {
        (f64::from(cutoff) * (std::f64::consts::PI * f64::from(xx)).sin()
            / (std::f64::consts::PI * f64::from(xx))
            * compute_func(
                (2.0f64 * f64::from(x) / f64::from(n)).abs() as f32,
                window_func,
            )) as f32
    }
}

fn compute_func(x: f32, func: &FuncDef) -> f64 {
    let mut interp: [f64; 4] = [0.; 4];
    let y: f32 = x * func.oversample as f32;
    let ind: i32 = (f64::from(y)).floor() as i32;
    let frac: f32 = y - ind as f32;
    interp[3usize] = -0.1666666667f64 * f64::from(frac)
        + 0.1666666667f64 * f64::from(frac * frac * frac);
    interp[2usize] = f64::from(frac) + 0.5f64 * f64::from(frac * frac)
        - 0.5f64 * f64::from(frac * frac * frac);
    interp[0usize] = -0.3333333333f64 * f64::from(frac)
        + 0.5f64 * f64::from(frac * frac)
        - 0.1666666667f64 * f64::from(frac * frac * frac);
    interp[1usize] =
        f64::from(1.0f32) - interp[3usize] - interp[2usize] - interp[0usize];

    interp[0usize] * func.table[ind as usize]
        + interp[1usize] * func.table[(ind + 1) as usize]
        + interp[2usize] * func.table[(ind + 2) as usize]
        + interp[3usize] * func.table[(ind + 3) as usize]
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
    let mut out_sample: i32 = 0;
    let mut last_sample: i32 = st.last_sample[channel_index as usize];
    let mut samp_frac_num: u32 = st.samp_frac_num[channel_index as usize];
    let out_stride: i32 = st.out_stride;
    let int_advance: i32 = st.int_advance;
    let frac_advance: i32 = st.frac_advance;
    let den_rate: u32 = st.den_rate;
    while !(last_sample >= *in_len as i32 || out_sample >= *out_len as i32) {
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
    out_sample
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
    let mut out_sample: i32 = 0;
    let mut last_sample: i32 = st.last_sample[channel_index as usize];
    let mut samp_frac_num: u32 = st.samp_frac_num[channel_index as usize];
    let out_stride: i32 = st.out_stride;
    let int_advance: i32 = st.int_advance;
    let frac_advance: i32 = st.frac_advance;
    let den_rate: u32 = st.den_rate;
    while !(last_sample >= *in_len as i32 || out_sample >= *out_len as i32) {
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
    st.last_sample[channel_index as usize] = last_sample;
    st.samp_frac_num[channel_index as usize] = samp_frac_num;
    out_sample
}

fn _muldiv(result: &mut u32, value: u32, mul: u32, div: u32) -> i32 {
    let major: u32 = value / div;
    let remainder: u32 = value % div;
    if remainder > 4294967295 / mul
        || major > 4294967295 / mul
        || major * mul > 4294967295 - remainder * mul / div
    {
        RESAMPLER_ERR_OVERFLOW as i32
    } else {
        *result = remainder * mul / div + major * mul;
        RESAMPLER_ERR_SUCCESS as i32
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
) -> i32 {
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
    if st.last_sample[channel_index as usize] < *in_len as i32 {
        *in_len = st.last_sample[channel_index as usize] as u32;
    }
    *out_len = out_sample as u32;
    st.last_sample[channel_index as usize] -= *in_len as i32;
    let ilen: u32 = *in_len;
    let mut j: i32 = 0;
    while j < n - 1 {
        st.mem[mem_idx + j as usize] =
            st.mem[mem_idx + (j as u32 + ilen) as usize];
        j += 1
    }
    RESAMPLER_ERR_SUCCESS as i32
}

fn speex_resampler_magic<'a, 'b>(
    st: &mut SpeexResamplerState,
    channel_index: u32,
    out: &'a mut &'b mut [f32],
    mut out_len: u32,
) -> i32 {
    let mut tmp_in_len = st.magic_samples[channel_index as usize];
    let mem_idx = (channel_index * st.mem_alloc_size) as usize;
    let n: i32 = st.filt_len as i32;
    speex_resampler_process_native(
        st,
        channel_index,
        &mut tmp_in_len,
        *out,
        &mut out_len,
    );
    st.magic_samples[channel_index as usize] -= tmp_in_len;
    if 0 != st.magic_samples[channel_index as usize] {
        let mut i: u32 = 0;
        while i < st.magic_samples[channel_index as usize] {
            st.mem[mem_idx + (n - 1) as usize + i as usize] = st.mem[mem_idx
                + (n - 1) as usize
                + i as usize
                + tmp_in_len as usize];
            i += 1;
        }
    }
    let value: &'b mut [f32] = mem::replace(out, &mut []);
    *out = &mut value[(out_len * st.out_stride as u32) as usize..];
    out_len as i32
}
