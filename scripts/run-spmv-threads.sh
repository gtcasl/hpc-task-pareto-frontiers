#!/bin/bash

for x in `seq 1 228`
do
    LD_LIBRARY_PATH=/home/eric OMP_NUM_THREADS=$x KMP_AFFINITY=scatter ./spmv 100 100 100
done
