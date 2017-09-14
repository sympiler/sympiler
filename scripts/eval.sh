#!/bin/bash
binFile=$1
matrixPath=$2
bigMat=$3
echo "Running $binFile for the set in $matrixPath"
for i in {1..5}
do
    for f in ${matrixPath}*.mtx
    do
        $binFile $f $bigMat
        echo ""
    done
    echo ""
done