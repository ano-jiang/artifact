#!/bin/bash

source env.sh

build_dir="cmake-build-release-intel"

cmake -B ${build_dir}/
cmake --build ${build_dir} -- VERBOSE=1
