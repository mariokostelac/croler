#!/usr/bin/env bash
set -ex
cd pipeline/qpid
make
cd ../brahle_assembly
if [[ ! -d bin ]]; then
    mkdir bin
fi
make
