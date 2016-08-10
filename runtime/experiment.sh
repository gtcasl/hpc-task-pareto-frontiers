#!/bin/bash

#export NUMTHREADS=112
#export POWERLIMIT=110

# for testing the baseline
for x in `seq 112 -1 1`
do
    echo "---------------------- CONTROL $x"
    NUMTHREADS=$x ./cholesky.control.sh
    mv cholesky.log cholesky.control.$x.log
done

# for testing the experiment
for x in `seq 30 5 220`
do
    echo "---------------------- EXPERIMENTAL $x"
    #POWERLIMIT=$x ./cg.sh
    #mv scheduler.cg.log cg.$x.log
    NUMTHREADS=112 POWERLIMIT=$x ./cholesky.sh
    mv cholesky.log cholesky.experimental.$x.log
done
