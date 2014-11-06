#!/usr/bin/env bash
set -ex
cd pipeline/qpid
make
cd ../brahle_assembly
make
