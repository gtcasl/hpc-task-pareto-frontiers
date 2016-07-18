#!/bin/bash

#export NUMTHREADS=228
#export POWERLIMIT=110

# Warmup
for y in `seq 1 2`
do
    #./cg.profile.sh
    ./cholesky.profile.sh
done

for x in `seq 1 228`
do
    echo "Running with $x threads"
    #NUMTHREADS=$x ./cg.profile.sh
    #mv profile.cg.log cg.$x.log
    NUMTHREADS=$x ./cholesky.profile.sh
    mv profile.cholesky.log cholesky.$x.log
done
