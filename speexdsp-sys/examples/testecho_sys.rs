extern crate byteorder;
extern crate speexdsp_sys;
extern crate structopt;

use byteorder::{BigEndian, ByteOrder};
use speexdsp_sys::echo::*;
use speexdsp_sys::preprocess::*;
use std::ffi::c_void;
use std::fs::File;
use std::io::{Read, Write};
use std::path::PathBuf;
use structopt::StructOpt;

macro_rules! preprocess_ctl {
    ($st:ident, $param:ident, $value:expr, $type:tt) => {
        let v = $value as *mut $type as *mut c_void;
        unsafe { speex_preprocess_ctl($st, $param as i32, v) };
    };
}

#[derive(Debug, StructOpt)]
struct CliArgs {
    #[structopt(parse(from_os_str))]
    echo_fd_path: PathBuf,
    #[structopt(parse(from_os_str))]
    ref_fd_path: PathBuf,
    #[structopt(parse(from_os_str))]
    e_fd_path: PathBuf,
}

fn main() -> std::io::Result<()> {
    const NN: usize = 160;
    const TAIL: usize = 1024;

    let sample_rate: i32 = 8000;

    let mut echo_buf: [i16; NN] = [0; NN];
    let mut echo_read_buf: [u8; NN * 2] = [0; NN * 2];
    let mut ref_buf: [i16; NN] = [0; NN];
    let mut ref_read_buf: [u8; NN * 2] = [0; NN * 2];
    let mut e_buf: [i16; NN] = [0; NN];

    let opts = CliArgs::from_args();

    let mut ref_fd = File::open(opts.ref_fd_path)?;
    let mut echo_fd = File::open(opts.echo_fd_path)?;
    let mut e_fd = File::create(opts.e_fd_path)?;

    let st = unsafe { speex_echo_state_init(NN as i32, TAIL as i32) };
    let den = unsafe { speex_preprocess_state_init(NN as i32, sample_rate) };
    unsafe {
        speex_echo_ctl(
            st,
            SPEEX_ECHO_SET_SAMPLING_RATE as i32,
            sample_rate as *mut c_void,
        )
    };
    preprocess_ctl!(
        den,
        SPEEX_PREPROCESS_SET_ECHO_STATE,
        st,
        SpeexPreprocessState
    );

    loop {
        let n = echo_fd.read(&mut echo_read_buf)?;
        let nn = ref_fd.read(&mut ref_read_buf)?;
        if n == 0 && nn == 0 {
            break;
        }
        BigEndian::read_i16_into(&echo_read_buf, &mut echo_buf);
        BigEndian::read_i16_into(&ref_read_buf, &mut ref_buf);
        unsafe {
            speex_echo_cancellation(
                st,
                ref_buf.as_ptr(),
                echo_buf.as_ptr(),
                e_buf.as_mut_ptr(),
            )
        };
        unsafe { speex_preprocess_run(den, e_buf.as_mut_ptr()) };
        BigEndian::write_i16_into(&e_buf, &mut ref_read_buf);
        e_fd.write_all(&ref_read_buf)?;
    }

    unsafe { speex_echo_state_destroy(st) };
    unsafe { speex_preprocess_state_destroy(den) };

    Ok(())
}
