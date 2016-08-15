#!/bin/bash

echo "LU"
mkdir lu-results
for x in `seq 112 -1 1`
do
    echo "---------------------- FAIR $x"
    export NUMTHREADS=$x
    for i in 1 2 3 4 5
    do
        if timeout 120s ./run -s fair lu 6 1400
        then
            break
        fi
    done
    mv scheduler.log lu.fair.$x.log
done

for x in `seq 112 -1 1`
do
    echo "---------------------- CORE CONSTRAINED $x"
    export NUMTHREADS=$x
    for i in `seq 1 10`
    do
        if timeout 300s ./run -s coreconstrained lu 6 1400
        then
            break
        fi
    done
    mv scheduler.log lu.coreconstrained.$x.log
done

for x in `seq 30 5 330`
do
    echo "---------------------- POWER AWARE $x"
    export NUMTHREADS=112
    export POWERLIMIT=$x
    for i in `seq 1 10`
    do
        if timeout 120s ./run -s poweraware lu 6 1400
        then
            break
        fi
    done
    mv scheduler.log lu.poweraware.$x.log
done

for x in `seq 30 5 330`
do
    echo "---------------------- SLACK AWARE $x"
    export NUMTHREADS=112
    export POWERLIMIT=$x
    for i in `seq 1 10`
    do
        if timeout 120s ./run -s slackaware lu 6 1400
        then
            break
        fi
    done
    mv scheduler.log lu.slackaware.$x.log
done
mv *.log lu-results

echo "CG"
mkdir cg-results
for x in `seq 112 -1 1`
do
    echo "---------------------- FAIR $x"
    export NUMTHREADS=$x
    for i in 1 2 3 4 5
    do
        if timeout 120s ./run -s fair cg 170 170 170 8
        then
            break
        fi
    done
    mv scheduler.log cg.fair.$x.log
done

for x in `seq 112 -1 1`
do
    echo "---------------------- CORE CONSTRAINED $x"
    export NUMTHREADS=$x
    for i in `seq 1 10`
    do
        if timeout 300s ./run -s coreconstrained cg 170 170 170 8
        then
            break
        fi
    done
    mv scheduler.log cg.coreconstrained.$x.log
done

for x in `seq 30 5 330`
do
    echo "---------------------- POWER AWARE $x"
    export NUMTHREADS=112
    export POWERLIMIT=$x
    for i in `seq 1 10`
    do
        if timeout 120s ./run -s poweraware cg 170 170 170 8
        then
            break
        fi
    done
    mv scheduler.log cg.poweraware.$x.log
done

for x in `seq 30 5 330`
do
    echo "---------------------- SLACK AWARE $x"
    export NUMTHREADS=112
    export POWERLIMIT=$x
    for i in `seq 1 10`
    do
        if timeout 120s ./run -s slackaware cg 170 170 170 8
        then
            break
        fi
    done
    mv scheduler.log cg.slackaware.$x.log
done
mv *.log cg-results
