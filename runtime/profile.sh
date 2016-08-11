#!/bin/bash

#export NUMTHREADS=228
#export POWERLIMIT=110

for x in `seq 112 -1 1`
do
    echo "---------------------- PROFILE $x"
    export NUMTHREADS=$x
    for i in 1 2 3 4 5
    do
        if timeout 300s ./run -s profiling choleskyprofiling 6 1400
        then
            break
        fi
    done
    mv scheduler.log cholesky.profile.$x.log
done
