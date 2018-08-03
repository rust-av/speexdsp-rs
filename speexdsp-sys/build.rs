extern crate autotools;
extern crate bindgen;
extern crate metadeps;

use std::path::PathBuf;
use std::io::Write;
use std::env;
use std::fs::File;

fn format_write(builder: bindgen::Builder) -> String {
    builder
        .generate()
        .unwrap()
        .to_string()
        .replace("/**", "/*")
        .replace("/*!", "/*")
}

fn main() {
    if cfg!(feature="build") {
        // TODO: decide how to fetch the source
        let dst = autotools::build("speexdsp");

        env::set_var("PKG_CONFIG_PATH", dst.join("lib/pkgconfig"));
    }

    let libs = metadeps::probe().unwrap();

    let headers = libs.get("speexdsp").unwrap().include_paths.clone();

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());

    for e in ["echo", "jitter", "preprocess", "resampler"].iter() {
        let mut builder = bindgen::builder().header(format!("data/{}.h", e));

        for header in headers.iter() {
            builder = builder.clang_arg("-I").clang_arg(header.to_str().unwrap());
        }

        // Manually fix the comment so rustdoc won't try to pick them
        let s = format_write(builder);

        let lib = format!("{}.rs", e);

        let mut file = File::create(out_path.join(&lib)).unwrap();

        let _ = file.write(s.as_bytes());
    }
}
