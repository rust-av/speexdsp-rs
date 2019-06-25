#![allow(dead_code)]
#![allow(unused_imports)]
extern crate byteorder;
extern crate speexdsp;
extern crate structopt;

use byteorder::{BigEndian, ByteOrder};
use std::fs::File;
use std::io::{Read, Write};
use std::path::PathBuf;
use structopt::StructOpt;

#[derive(Debug, StructOpt)]
struct CliArgs {
    #[structopt(parse(from_os_str))]
    echo_fd_path: PathBuf,
    #[structopt(parse(from_os_str))]
    ref_fd_path: PathBuf,
    #[structopt(parse(from_os_str))]
    e_fd_path: PathBuf,
}

#[cfg(feature = "sys")]
fn main() -> std::io::Result<()> {
    use speexdsp::echo::SpeexEchoConst::*;
    use speexdsp::echo::*;
    use speexdsp::preprocess::SpeexPreprocessConst::*;
    use speexdsp::preprocess::*;

    const NN: usize = 160;
    const TAIL: usize = 1024;

    let sample_rate: usize = 8000;

    let mut echo_buf: [i16; NN] = [0; NN];
    let mut echo_read_buf: [u8; NN * 2] = [0; NN * 2];
    let mut ref_buf: [i16; NN] = [0; NN];
    let mut ref_read_buf: [u8; NN * 2] = [0; NN * 2];
    let mut e_buf: [i16; NN] = [0; NN];

    let opts = CliArgs::from_args();

    let mut ref_fd = File::open(opts.ref_fd_path)?;
    let mut echo_fd = File::open(opts.echo_fd_path)?;
    let mut e_fd = File::create(opts.e_fd_path)?;

    let mut st = SpeexEcho::new(NN, TAIL).unwrap();
    let mut den = SpeexPreprocess::new(NN, sample_rate).unwrap();
    st.echo_ctl(SPEEX_ECHO_SET_SAMPLING_RATE, sample_rate)
        .unwrap();
    den.preprocess_ctl(SPEEX_PREPROCESS_SET_ECHO_STATE, &st)
        .unwrap();

    loop {
        let n = echo_fd.read(&mut echo_read_buf)?;
        let nn = ref_fd.read(&mut ref_read_buf)?;
        if n == 0 && nn == 0 {
            break;
        }
        BigEndian::read_i16_into(&echo_read_buf, &mut echo_buf);
        BigEndian::read_i16_into(&ref_read_buf, &mut ref_buf);
        st.echo_cancellation(&ref_buf, &echo_buf, &mut e_buf);
        den.preprocess_run(&mut e_buf);
        BigEndian::write_i16_into(&e_buf, &mut ref_read_buf);
        e_fd.write_all(&ref_read_buf)?;
    }

    Ok(())
}

#[cfg(not(feature = "sys"))]
fn main() {
    unimplemented!();
}
