#!/usr/bin/env bash

if [[ $# -ne 1 ]]; then
    echo "usage: $0 <reads_file.afg>"
    exit 1
fi

line="---------------------------------------------------------------"

BASE=$(basename $1 .afg)

LOG_FILE="${BASE}.log"
TMP_OVERLAPS=${BASE}_overlaps.afg
TMP_LAYOUT=${BASE}_layout.afg
TMP_CONSENSUS=${BASE}_consensus.fasta

suff=$(date +%s)

# prepare log file
if [[ -f $LOG_FILE ]]; then
  echo "Found existing ${LOG_FILE}"
  echo "Moving ${LOG_FILE} to ${LOG_FILE}_${suff}"
  mv ${LOG_FILE} ${LOG_FILE}_${suff}
fi
touch $LOG_FILE

# move old results if exist
if [[ -f $TMP_CONSENSUS ]]; then
  echo "Found existing ${TMP_CONSENSUS}"
  echo "Moving ${TMP_CONSENSUS} to ${TMP_CONSENSUS}_${suff}"
  mv ${TMP_CONSENSUS} ${TMP_CONSENSUS}_${suff}
fi

set -e

echo
echo "Feel free to put this in background and start 'tail -f $LOG_FILE'"
echo

echo "running OVERLAP phase..."
echo "OVERLAP" >> $LOG_FILE
echo $line >> $LOG_FILE
time ./bin/overlap $1 -o $TMP_OVERLAPS &>> $LOG_FILE
echo $line >> $LOG_FILE

echo "running LAYOUT phase..."
echo "LAYOUT" >> $LOG_FILE
echo $line >> $LOG_FILE
time ./bin/layout $1 $TMP_OVERLAPS &>> $LOG_FILE
mv layout.afg $TMP_LAYOUT
echo $line >> $LOG_FILE

echo "running CONSENSUS phase..."
echo "CONSENSUS" >> $LOG_FILE
echo $line >> $LOG_FILE
time ./bin/consensus $1 $TMP_LAYOUT &>> $LOG_FILE
mv consensus.fasta $TMP_CONSENSUS
echo $line >> $LOG_FILE

echo $line 
echo "RESULTS  $TMP_CONSENSUS"
echo "LOGS     $LOG_FILE"
echo "OVERLAPS $TMP_OVERLAPS"
echo "LAYOUT   $TMP_LAYOUT"

