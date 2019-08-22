# speexdsp bindings and c2rust version

[![LICENSE](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Actions Status](https://github.com/rust-av/speexdsp-rs/workflows/speexdsp/badge.svg)](https://github.com/rust-av/speexdsp-rs/actions)
[![dependency status](https://deps.rs/repo/github/rust-av/speexdsp-rs/status.svg)](https://deps.rs/repo/github/rust-av/speexdsp-rs)
[![IRC](https://img.shields.io/badge/irc-%23rust--av-blue.svg)](http://webchat.freenode.net?channels=%23rust-av&uio=d4)

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

## Testing

Currently we have only an integration test to compare the C and the Rust implementation.
To run it issue:

``` sh
$ cargo test --features=sys
```

[1]: https://github.com/servo/rust-bindgen
[2]: https://github.com/xiph/speexdsp
