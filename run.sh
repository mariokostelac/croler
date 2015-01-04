#!/usr/bin/env bash

if [[ $# -ne 1 ]]; then
    echo "usage: $0 <reads_file.afg>"
    exit 1
fi

line="---------------------------------------------------------------"

TMP_OVERLAPS=overlaps.afg
TMP_LAYOUT=layout.afg

set -e

echo "OVERLAP phase"
echo $line
./bin/overlap $1 -o $TMP_OVERLAPS
echo $line
echo

echo "LAYOUT phase"
echo $line
./bin/layout $1 $TMP_OVERLAPS
echo $line
echo

echo "CONSENSUS phase"
echo $line
./bin/consensus $1 $TMP_LAYOUT
echo $line
