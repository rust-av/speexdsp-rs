use crate::smallft::*;

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
pub fn spx_fft_init(size: usize) -> DrftLookup {
    DrftLookup::new(size)
}

pub fn spx_fft(table: &mut DrftLookup, in_0: &mut [f32], out: &mut [f32]) {
    let scale = (1.0f64 / table.n as f64) as f32;
    if in_0 == out {
        eprintln!("FFT should not be done in-place");
    }

    out.iter_mut()
        .zip(in_0.iter())
        .take(table.n as usize)
        .for_each(|(o, i)| *o = scale * *i);

    spx_drft_forward(table, out);
}

pub fn spx_ifft(table: &mut DrftLookup, in_0: &mut [f32], out: &mut [f32]) {
    if in_0 == out {
        eprintln!("FFT should not be done in-place");
    } else {
        out.copy_from_slice(&in_0[..table.n as usize]);
    }

    spx_drft_backward(table, out);
}

pub fn spx_fft_float(
    table: &mut DrftLookup,
    in_0: &mut [f32],
    out: &mut [f32],
) {
    spx_fft(table, in_0, out);
}

pub fn spx_ifft_float(
    table: &mut DrftLookup,
    in_0: &mut [f32],
    out: &mut [f32],
) {
    spx_ifft(table, in_0, out);
}
