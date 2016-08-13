#!/bin/bash

for x in getrf lu_gemm trsm_left trsm_right
do
    infile=${x}profile.avg.csv
    outfile=${x}.csv
    tmpfile=${x}.tmp
    tail -n +2 ../gemm.csv > $tmpfile
    echo "nthreads,energy,time,power,speedup" > $outfile
    paste -d ',' $tmpfile $infile | awk -F "," 'BEGIN {OFS=","} {print $1,$2,$9,$8,$5}' >> $outfile
    rm $tmpfile
done
