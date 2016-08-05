#!/bin/bash

export NUMTHREADS=112
#export POWERLIMIT=110

# Warmup
#for y in `seq 1 2`
#do
#    #./cg.sh
#    ./cholesky.sh
#done

for x in 40 50 60 1000
do
    #POWERLIMIT=$x ./cg.sh
    #mv scheduler.cg.log cg.$x.log
    POWERLIMIT=$x ./cholesky.sh
    mv scheduler.cholesky.log cholesky.$x.log
done
