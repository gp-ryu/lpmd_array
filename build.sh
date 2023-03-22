#!/usr/bin/env bash

set -o xtrace -o nounset -o pipefail -o errexit

PREFIX=$(realpath $(dirname ${BASH_SOURCE[0]}))

cmake -S . -G Ninja -B build \
    -DCMAKE_INSTALL_PREFIX=${PREFIX} \
    -DCMAKE_INSTALL_LIBDIR=${PREFIX}/lib \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS="-I/home/gryu/micromamba/envs/lpmd_dev/include" 
    #-DCMAKE_LD_FLAGS="-L/home/gryu/micromamba/envs/lpmd_dev/lib -Wl,-rpath,/home/gryu/micromamba/envs/lpmd_dev/lib"
cmake --build build
cmake --install build
