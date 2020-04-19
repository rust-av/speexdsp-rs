/* Copyright (C) 2005-2006 Jean-Marc Valin
   File: fftwrap.c

   Wrapper for various FFTs

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

   - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   - Neither the name of the Xiph.org Foundation nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

use crate::smallft::*;

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

    pub fn spx_fft(&mut self, in_0: &[f32], out: &mut [f32]) {
        let scale = (1.0f64 / self.n as f64) as f32;
        if in_0 == out {
            eprintln!("FFT should not be done in-place");
        }

        out.iter_mut()
            .zip(in_0.iter())
            .take(self.n as usize)
            .for_each(|(o, i)| *o = scale * *i);

        self.spx_drft_forward(out);
    }

    pub fn spx_ifft(&mut self, in_0: &[f32], out: &mut [f32]) {
        if in_0 == out {
            eprintln!("FFT should not be done in-place");
        } else {
            out.copy_from_slice(&in_0[..self.n as usize]);
        }

        self.spx_drft_backward(out);
    }

    pub fn spx_fft_float(&mut self, in_0: &[f32], out: &mut [f32]) {
        self.spx_fft(in_0, out);
    }

    pub fn spx_ifft_float(&mut self, in_0: &[f32], out: &mut [f32]) {
        self.spx_ifft(in_0, out);
    }

    pub fn spx_drft_forward(&mut self, data: &mut [f32]) {
        if self.n == 1 {
            return;
        }

        let mut trigcache_temp = self.trigcache[self.n as usize..].to_vec();

        drftf1(
            self.n as i32,
            data,
            &mut self.trigcache,
            &mut trigcache_temp,
            &mut self.splitcache,
        );

        self.trigcache[self.n as usize..].copy_from_slice(&trigcache_temp);
    }

    pub fn spx_drft_backward(&mut self, data: &mut [f32]) {
        if self.n == 1 {
            return;
        }

        let mut trigcache_temp = self.trigcache[self.n as usize..].to_vec();

        drftb1(
            self.n as i32,
            data,
            &mut self.trigcache,
            &mut trigcache_temp,
            &mut self.splitcache,
        );

        self.trigcache[self.n as usize..].copy_from_slice(&trigcache_temp);
    }
}
