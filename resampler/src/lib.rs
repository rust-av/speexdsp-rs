mod speex;

pub struct State {
    st: speex::SpeexResamplerState,
}

#[derive(Debug, Clone, Copy)]
pub enum Error {
    Fail,
    InvalidArg,
    Overflow,
}

impl State {
    pub fn new(
        channels: usize,
        in_rate: usize,
        out_rate: usize,
        quality: usize,
    ) -> Result<Self, Error> {
        let st = speex::SpeexResamplerState::new(channels, in_rate, out_rate, quality);

        Ok(State { st })
    }

    pub fn set_rate(&mut self, in_rate: usize, out_rate: usize) -> Result<(), Error> {
        if self.st.set_rate(in_rate, out_rate) != 0 {
            Err(Error::InvalidArg)
        } else {
            Ok(())
        }
    }

    pub fn get_rate(&self) -> (usize, usize) {
        self.st.get_rate()
    }

    pub fn process_float(
        &mut self,
        index: usize,
        input: &[f32],
        output: &mut [f32],
    ) -> Result<(usize, usize), Error> {
        let mut in_len = input.len() as u32;
        let mut out_len = output.len() as u32;
        let ret = self
            .st
            .process_float(index as u32, input, &mut in_len, output, &mut out_len);

        if ret != 0 {
            Err(Error::Fail)
        } else {
            Ok((in_len as usize, out_len as usize))
        }
    }

    pub fn skip_zeros(&mut self) {
        self.st.skip_zeros();
    }

    pub fn reset(&mut self) {
        self.st.reset_mem();
    }

    pub fn get_input_latency(&self) -> usize {
        self.st.get_input_latency()
    }

    pub fn get_output_latency(&self) -> usize {
        self.st.get_output_latency()
    }

    pub fn set_quality(&mut self, quality: usize) -> Result<(), Error> {
        if self.st.set_quality(quality) != 0 {
            Err(Error::InvalidArg)
        } else {
            Ok(())
        }
    }

    pub fn get_quality(&self) -> usize {
        self.st.get_quality()
    }
}
