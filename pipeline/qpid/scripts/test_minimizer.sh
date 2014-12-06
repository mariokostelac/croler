#!/bin/bash

for test in test_data/minimizer/test*.fasta
do
    echo -n $test"..."
    ./bin/test_minimizer $test > tmp.out
    res=$(diff -w $test".out" tmp.out && echo "PASS" || echo "FAIL")
    echo $res
done
rm tmp.out
