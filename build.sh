#!/bin/bash
#usage: . ./build.sh

if [ $# -eq 0 ]; then
    echo "Usage: . ./build.sh arg"
    return
else
    platform=$1
fi

source env.sh

build_dir="cmake-build-release-intel"

cmake -B ${build_dir}/
cmake --build ${build_dir} -- VERBOSE=1
