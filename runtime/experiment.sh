#!/bin/bash

#export NUMTHREADS=228
#export POWERLIMIT=110

# Warmup
for y in `seq 1 2`
do
    ./cg.sh
    ./cholesky.sh
done

for x in 1 5 10 20 30 40 50 60
do
    POWERLIMIT=$x ./cg.sh
    mv scheduler.cg.log cg.$x.log
    POWERLIMIT=$x ./cholesky.sh
    mv scheduler.cholesky.log cholesky.$x.log
done
./cg.sh
mv scheduler.cg.log cg.none.log
./cholesky.sh
mv scheduler.cholesky.log cholesky.none.log
