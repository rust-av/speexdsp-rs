use byteorder::{BigEndian, ByteOrder};
use speexdsp_sys::jitter::*;
use std::ptr;

macro_rules! null_struct {
    ($v:ident) => {
        let mut $v = JitterBufferPacket {
            data: std::ptr::null_mut::<i8>(),
            len: 0,
            timestamp: 0,
            span: 0,
            sequence: 0,
            user_data: 0,
        };
    };
}

fn synth_in(input: &mut JitterBufferPacket, idx: usize, span: usize) {
    let mut buf = [0; 4];
    BigEndian::write_u32(&mut buf, idx as u32);
    input.data = buf.as_mut_ptr() as *mut i8;
    input.len = 32;
    input.timestamp = (idx * 10) as u32;
    input.span = (span * 10) as u32;
    input.sequence = idx as u16;
    input.user_data = 0;
}

fn jitter_fill(jb: *mut JitterBuffer) {
    let mut buffer: [i8; 65536] = [0; 65536];

    null_struct!(input);
    null_struct!(output);

    output.data = buffer.as_mut_ptr();

    unsafe {
        jitter_buffer_reset(jb);
    }

    for i in 0..100 {
        synth_in(&mut input, i, 1);
        unsafe {
            jitter_buffer_put(jb, &input);
        }

        output.len = 65536;
        let err = unsafe {
            jitter_buffer_get(jb, &mut output, 10, ptr::null_mut::<i32>())
        };
        if err != (JITTER_BUFFER_OK as i32) {
            eprintln!("Fill test failed iteration {}", i);
        }
        if output.timestamp != (i * 10) as u32 {
            println!("Fill test expected {} got {}", i * 10, output.timestamp);
        }
        unsafe {
            jitter_buffer_tick(jb);
        }
    }
}

fn main() {
    let mut buffer: [i8; 65536] = [0; 65536];
    let jb = unsafe { jitter_buffer_init(10) };

    null_struct!(input);
    null_struct!(output);

    output.data = buffer.as_mut_ptr();

    jitter_fill(jb);

    for _ in 0..100 {
        output.len = 65536;
        unsafe {
            jitter_buffer_get(jb, &mut output, 10, ptr::null_mut::<i32>());
            jitter_buffer_tick(jb);
        }
    }

    synth_in(&mut input, 100, 1);
    unsafe {
        jitter_buffer_put(jb, &input);
    }
    output.len = 65536;
    let err = unsafe {
        jitter_buffer_get(jb, &mut output, 10, ptr::null_mut::<i32>())
    };
    if err != (JITTER_BUFFER_OK as i32) {
        eprintln!("Failed frozen sender resynchronize");
    } else {
        println!("Frozen sender: Jitter {}", output.timestamp - 100 * 10);
    }
}
