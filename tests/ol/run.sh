#!/usr/bin/env bash
set -e

input=$1
shift

bank_we="bank.we"
overlaps_we="overlaps.we"
layout_we="layout.we"
contigs_we="contigs.we"

bank_minimus="bank.minimus"
overlaps_minimus="overlaps.minimus"
layout_minimus="layout.minimus"
contigs_minimus="contigs.minimus"

# move our
touch blank.we
for i in *.we; do
  if [[ -e $i.old ]]; then
    rm -r $i.old
  fi
  mv $i $i.old
done
rm blank.we.old

# move minimus
touch blank.minimus
for i in *.minimus; do
  if [[ -e $i.old ]]; then
    rm -r $i.old
  fi
  mv $i $i.old
done
rm blank.minimus.old

###### OUR PIPELINE ###############################
bank-transact -c -m $input -b $bank_we

# our overlap
croler-overlap $input -o $overlaps_we
bank-transact -m $overlaps_we -b $bank_we

# our layout
croler-layout $input $overlaps_we
mv layout.afg $layout_we
bank-transact -m $layout_we -b $bank_we

# minimus consensus
make-consensus -B -b $bank_we &> consensus.log.we
bank2fasta -b $bank_we > contigs.we

###### MINIMUS PIPELINE ##########################
bank-transact -c -m $input -b $bank_minimus

# minimus overlap
hash-overlap $bank_minimus -B $bank_minimus

# minimus layout
tigger -b $bank_minimus

# minimus consensus
make-consensus -B -b $bank_minimus &> consensus.log.minimus
bank2fasta -b $bank_minimus > contigs.minimus
