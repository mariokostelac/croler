#!/usr/bin/env bash

if [[ $# -ne 2 ]]; then
    echo "usage: $0 <reads_file.fasta> <reads_file.afg>"
    exit 1
fi

TMP_OVERLAPS=overlaps.afg

./pipeline/qpid/bin/overlap $1 -o $TMP_OVERLAPS
./pipeline/brahle_assembly/bin/main_layout $2 $TMP_OVERLAPS
