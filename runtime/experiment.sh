#!/bin/bash

#export NUMTHREADS=112
#export POWERLIMIT=110

for x in `seq 112 -1 1`
do
    echo "---------------------- FAIR $x"
    export NUMTHREADS=$x
    for i in 1 2 3 4 5
    do
        if timeout 120s ./run -s fair cholesky 6 1400
        then
            break
        fi
    done
    mv scheduler.log cholesky.fair.$x.log
done

#for x in 60 84
#do
#    echo "---------------------- CORE CONSTRAINED $x"
#    export NUMTHREADS=$x
#    for i in `seq 1 10`
#    do
#        if timeout 300s ./run -s coreconstrained cholesky 6 1400
#        then
#            break
#        fi
#    done
#    mv scheduler.log cholesky.coreconstrained.$x.log
#done

#for x in 295
#do
#    echo "---------------------- POWER AWARE $x"
#    export NUMTHREADS=112
#    export POWERLIMIT=$x
#    for i in `seq 1 10`
#    do
#        if timeout 120s ./run -s poweraware cholesky 6 1400
#        then
#            break
#        fi
#    done
#    mv scheduler.log cholesky.poweraware.$x.log
#done

#for x in 320
#do
#    echo "---------------------- SLACK AWARE $x"
#    export NUMTHREADS=112
#    export POWERLIMIT=$x
#    for i in `seq 1 10`
#    do
#        if timeout 120s ./run -s slackaware cholesky 6 1400
#        then
#            break
#        fi
#    done
#    mv scheduler.log cholesky.slackaware.$x.log
#done
