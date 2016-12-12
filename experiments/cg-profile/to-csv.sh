#!/bin/bash

for x in assign comp_alpha comp_beta copy dot spmv start subtract sum_contribs
do
    infile=${x}.avg.csv
    outfile=${x}.csv
    tmpfile=${x}.tmp
    tail -n +2 ../gemm.csv > $tmpfile
    echo "nthreads,energy,time,power,speedup" > $outfile
    paste -d ',' $tmpfile $infile | awk -F "," 'BEGIN {OFS=","} {print $1,$2,$9,$8,$5}' >> $outfile
    rm $tmpfile
done
