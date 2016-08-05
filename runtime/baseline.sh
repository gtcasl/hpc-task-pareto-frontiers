#!/bin/bash

#export NUMTHREADS=112
#export POWERLIMIT=110

# Warmup
#for y in `seq 1 2`
#do
#    #./cg.sh
#    ./cholesky.sh
#done

for x in `seq 15 5 30`
do
    #POWERLIMIT=$x ./cg.sh
    #mv scheduler.cg.log cg.$x.log
    NUMTHREADS=$x ./run -s baseline cholesky 6 1400
    mv scheduler.log cholesky.$x.baseline.log
done
