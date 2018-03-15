extern crate autotools;
extern crate bindgen;
extern crate metadeps;

use std::env;
use std::fs::OpenOptions;
use std::io::Write;

fn format_write(builder: bindgen::Builder, output: &str) {
    let s = builder
        .generate()
        .unwrap()
        .to_string()
        .replace("/**", "/*")
        .replace("/*!", "/*");

    let mut file = OpenOptions::new()
        .write(true)
        .truncate(true)
        .create(true)
        .open(output)
        .unwrap();

    let _ = file.write(s.as_bytes());
}

fn common_builder() -> bindgen::Builder {
    bindgen::builder()
        .raw_line("#![allow(dead_code)]")
        .raw_line("#![allow(non_camel_case_types)]")
        .raw_line("#![allow(non_snake_case)]")
        .raw_line("#![allow(non_upper_case_globals)]")
}

fn main() {
    if cfg!(feature="build") {
        // TODO: decide how to fetch the source
        let dst = autotools::build("speexdsp");

        env::set_var("PKG_CONFIG_PATH", dst.join("lib/pkgconfig"));
    }

    let libs = metadeps::probe().unwrap();

    let headers = libs.get("speexdsp").unwrap().include_paths.clone();

    for e in ["echo", "jitter", "preprocess", "resampler"].iter() {
        let mut builder = common_builder().header(format!("data/{}.h", e));

        for header in headers.iter() {
            builder = builder.clang_arg("-I").clang_arg(header.to_str().unwrap());
        }

        // Manually fix the comment so rustdoc won't try to pick them
        format_write(builder, &format!("src/{}.rs", e));
    }
}
