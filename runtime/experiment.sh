#!/bin/bash

#export NUMTHREADS=228
#export POWERLIMIT=110

./cg.sh
./cholesky.sh
tail -n 3 scheduler.cg.log > cg.power.log
tail -n 3 scheduler.cholesky.log > cholesky.power.log

for x in 1 5 10 20 30 40 50 60
do
    POWERLIMIT=$x ./cg.sh
    POWERLIMIT=$x ./cholesky.sh
    tail -n 3 scheduler.cg.log >> cg.power.log
    tail -n 3 scheduler.cholesky.log >> cholesky.power.log
done
