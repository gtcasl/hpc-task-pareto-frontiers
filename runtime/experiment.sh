#!/bin/bash

#export NUMTHREADS=228
#export POWERLIMIT=110

for y in `seq 1 10`
do
    #for x in 1 5 10 20 30 40 50 60
    #do
    #    POWERLIMIT=$x ./cg.sh
    #    POWERLIMIT=$x ./cholesky.sh
    #    tail -n 3 scheduler.cg.log >> cg.$x.log
    #    tail -n 3 scheduler.cholesky.log >> cholesky.$x.log
    #done
    #./cg.sh
    ./cholesky.sh
    #tail -n 3 scheduler.cg.log >> cg.none.log
    tail -n 3 scheduler.cholesky.log >> cholesky.none.log
done

