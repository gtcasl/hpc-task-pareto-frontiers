#!/bin/bash

for nx in 100 120 140 160 180 200
do
    for x in `seq 1 228`
    do
        OMP_NUM_THREADS=$x KMP_AFFINITY=scatter ./spmv 100 100 $nx
    done
    mv spmv.csv spmv.$nx.csv
done
