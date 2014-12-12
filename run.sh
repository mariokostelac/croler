#!/usr/bin/env bash

if [[ $# -ne 1 ]]; then
    echo "usage: $0 <reads_file.fasta>"
    exit 1
fi

TMP_OVERLAPS=overlaps.afg
TMP_LAYOUT=layout.afg

./pipeline/qpid/bin/overlap $1 -o $TMP_OVERLAPS
./pipeline/brahle_assembly/bin/main_layout $1 $TMP_OVERLAPS
grbin-msa $1 $TMP_LAYOUT
