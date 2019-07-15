#![allow(unused_imports)]
extern crate byteorder;
#[cfg(feature = "sys")]
extern crate speexdsp;

use byteorder::{BigEndian, ByteOrder};

#[cfg(feature = "sys")]
use speexdsp::jitter::*;

#[cfg(feature = "sys")]
fn synth_in(input: &mut SpeexBufferPacket, idx: usize, span: usize) {
    let mut buf = [0; 4];
    BigEndian::write_u32(&mut buf, idx as u32);
    let mut buf_i8: [i8; 4] = [0; 4];
    for i in 0..3 {
        buf_i8[i] = buf[i] as i8;
    }
    input.create(&mut buf_i8, 32, idx * 10, span * 10, idx, 0);
}

#[cfg(feature = "sys")]
fn jitter_fill(jb: &mut SpeexJitter) {
    let mut buffer: [i8; 65536] = [0; 65536];

    let mut input = SpeexBufferPacket::new();
    let mut output = SpeexBufferPacket::new();

    output.set_data(&mut buffer);

    jb.buffer_reset();

    for i in 0..100 {
        synth_in(&mut input, i, 1);
        jb.buffer_put(&input);

        output.set_len(65536);
        if jb.buffer_get(&mut output, 10, 0) != Error::BufferOk {
            eprintln!("Fill test failed iteration {}", i);
        }
        if output.timestamp() != i * 10 {
            println!(
                "Fill test expected {} got {}",
                i * 10,
                output.timestamp()
            );
        }
        jb.buffer_tick();
    }
}

#[cfg(feature = "sys")]
fn main() {
    let mut buffer: [i8; 65536] = [0; 65536];
    let mut jb = SpeexJitter::new(10).unwrap();

    let mut input = SpeexBufferPacket::new();
    let mut output = SpeexBufferPacket::new();

    output.set_data(&mut buffer);

    jitter_fill(&mut jb);
    for _ in 0..100 {
        output.set_len(65536);
        jb.buffer_get(&mut output, 10, 0);
        jb.buffer_tick();
    }

    synth_in(&mut input, 100, 1);
    jb.buffer_put(&input);
    output.set_len(65536);
    if jb.buffer_get(&mut output, 10, 0) != Error::BufferOk {
        eprintln!("Failed frozen sender resynchronize");
    } else {
        println!("Frozen sender: Jitter {}", output.timestamp() - 100 * 10);
    }
}

#[cfg(not(feature = "sys"))]
fn main() {
    unimplemented!();
}
