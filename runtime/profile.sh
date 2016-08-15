#!/bin/bash

#export NUMTHREADS=228
#export POWERLIMIT=110

export MIC_ENV_PREFIX=PHI
export PHI_KMP_AFFINITY=granularity=fine,balanced
export PHI_KMP_PLACE_THREADS=2t
for x in `seq 1 1 112`
do
    echo "---------------------- balanced $x"
    export NUMTHREADS=$x
    for i in 1 2 3 4 5
    do
        if timeout 300s ./run -s profiling cgprofiling 170 170 170 8
        then
            break
        fi
    done
    mv scheduler.log balanced.$x.log
done
