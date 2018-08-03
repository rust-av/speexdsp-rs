# speexdsp bindings and c2rust version

[![LICENSE](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

It is a simple safe abstraction based on [speexdsp][2].

It is available as [binding][1] or as pure-rust implementation.

## Building

By default the pure-rust implementation is used, optionally the simd-accelerated original
C version can be used instead using the feature `sys`.

The bindings are generated using the headers and libraries that ought to be present in the system.

- Make sure you have `clang` and `libclang` installed.
- Make sure the `speexdsp` C headers and pkg-config files are installed.

## TODO
- [ ] Source build speexdsp
- [x] Simple bindings
- [x] Safe abstraction
- [x] Examples
- [ ] Clean pure-rust reimplementation

[1]: https://github.com/servo/rust-bindgen
[2]: https://github.com/xiph/speexdsp
