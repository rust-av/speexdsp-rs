extern crate autotools;
extern crate bindgen;
extern crate system_deps;

use std::env;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

fn format_write(builder: bindgen::Builder) -> String {
    builder
        .generate()
        .unwrap()
        .to_string()
        .replace("/**", "/*")
        .replace("/*!", "/*")
}

fn main() {
    let libs = system_deps::Config::new()
        .add_build_internal("speexdsp", |lib, version| {
            // TODO: decide how to fetch the source
            let dst = autotools::build("speexdsp");
            system_deps::Library::from_internal_pkg_config(&dst, lib, version)
        })
        .probe()
        .unwrap();

    let headers = libs.get_by_name("speexdsp").unwrap().include_paths.clone();

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());

    for e in ["echo", "jitter", "preprocess", "resampler"].iter() {
        let mut builder = bindgen::builder()
            .size_t_is_usize(true)
            .layout_tests(false)
            .header(format!("data/{}.h", e));

        for header in headers.iter() {
            builder =
                builder.clang_arg("-I").clang_arg(header.to_str().unwrap());
        }

        // Manually fix the comment so rustdoc won't try to pick them
        let s = format_write(builder);

        let lib = format!("{}.rs", e);

        let mut file = File::create(out_path.join(&lib)).unwrap();

        let _ = file.write(s.as_bytes());
    }
}
