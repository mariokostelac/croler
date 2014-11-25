#!/usr/bin/env bash
set -ex

# build overlap phase
cd pipeline/qpid
make

# build layout phase
cd ../brahle_assembly
if [[ ! -d bin ]]; then
    mkdir bin
fi
make

# build consensus phase
cd ../msa
if [[ ! -d bin ]]; then
    mkdir bin
fi
make
