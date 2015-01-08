#!/usr/bin/env bash

if [[ $# -ne 1 ]]; then
    echo "usage: $0 <reads_file.afg>"
    exit 1
fi

line="---------------------------------------------------------------"

BASE=$(basename $1 .afg)

TMP_OVERLAPS=${BASE}_overlaps.afg
TMP_LAYOUT=${BASE}_layout.afg
TMP_CONSENSUS=${BASE}_consensus.fasta

echo "files"
echo $TMP_OVERLAPS
echo $TMP_LAYOUT
echo $TMP_CONSENSUS

set -e

echo "OVERLAP phase"
echo $line
./bin/overlap $1 -o $TMP_OVERLAPS
echo $line
echo

echo "LAYOUT phase"
echo $line
./bin/layout $1 $TMP_OVERLAPS
mv layout.afg $TMP_LAYOUT
echo $line
echo

echo "CONSENSUS phase"
echo $line
./bin/consensus $1 $TMP_LAYOUT
mv consensus.fasta $TMP_CONSENSUS
echo $line
