#!/bin/bash

#export NUMTHREADS=228
#export POWERLIMIT=110

export MIC_ENV_PREFIX=PHI
export PHI_KMP_AFFINITY=granularity=fine,balanced
export PHI_KMP_PLACE_THREADS=2t
for freq in 600000 650000 700000 800000 900000 1000000 1100000
do
    ssh mic0 'bash -s' < ../scripts/setfreq.sh $freq
    echo "FREQ: $freq"
    for x in `seq 1 1 112`
    do
        echo "---------------------- $x cg"
        export NUMTHREADS=$x
        for i in 1 2 3 4 5
        do
            if timeout 300s ./run -s profiling cgprofiling 170 170 170 8
            then
                break
            fi
        done
        mv scheduler.log freq-experiment/cgprofile.$freq.$x.log
    done
    for x in `seq 1 1 112`
    do
        echo "---------------------- $x lu"
        export NUMTHREADS=$x
        for i in 1 2 3 4 5
        do
            if timeout 300s ./run -s profiling luprofiling 6 1400
            then
                break
            fi
        done
        mv scheduler.log freq-experiment/luprofile.$freq.$x.log
    done
    for x in `seq 1 1 112`
    do
        echo "---------------------- $x cholesky"
        export NUMTHREADS=$x
        for i in 1 2 3 4 5
        do
            if timeout 300s ./run -s profiling choleskyprofiling 6 1400
            then
                break
            fi
        done
        mv scheduler.log freq-experiment/choleskyprofile.$freq.$x.log
    done
done
