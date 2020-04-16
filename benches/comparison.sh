#!/usr/bin/env sh

set -e

echo "Cloning speexdsp"
git clone --depth 1  https://github.com/xiph/speexdsp target/speexdsp_repo || true
cd target/speexdsp
echo "Running ./autogen.sh from speexdsp repo"
./autogen.sh
echo "\nRunning ./configure with disabled sse from speexdsp repo\n"
./configure --disable-sse --quiet
echo "\nCompiling speexdsp\n"
make -s
export DSP_LIB=$(pwd)/libspeexdsp/.libs/libspeexdsp.so.1
cd ../../
echo "\nRunning resampler_c with LD_PRELOAD=${DSP_LIB}\n"
LD_PRELOAD=$DSP_LIB cargo bench -q --features sys -- resampler_c
echo "\nRunning resampler_rust\n"
cargo bench -q -- resampler_rust
