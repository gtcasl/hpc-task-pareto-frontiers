#!/bin/bash

#echo "LU"
#mkdir lu-results
#for x in `seq 112 -1 1`
#do
#    echo "---------------------- FAIR $x"
#    export NUMTHREADS=$x
#    for i in 1 2 3 4 5
#    do
#        if timeout 120s ./run -s fair lu 6 1400
#        then
#            break
#        fi
#    done
#    mv scheduler.log lu.fair.$x.log
#done
#
#for x in `seq 112 -1 1`
#do
#    echo "---------------------- CORE CONSTRAINED $x"
#    export NUMTHREADS=$x
#    for i in `seq 1 10`
#    do
#        if timeout 300s ./run -s coreconstrained lu 6 1400
#        then
#            break
#        fi
#    done
#    mv scheduler.log lu.coreconstrained.$x.log
#done
#
#for x in `seq 30 5 330`
#do
#    echo "---------------------- POWER AWARE $x"
#    export NUMTHREADS=112
#    export POWERLIMIT=$x
#    for i in `seq 1 10`
#    do
#        if timeout 120s ./run -s poweraware lu 6 1400
#        then
#            break
#        fi
#    done
#    mv scheduler.log lu.poweraware.$x.log
#done
#
#for x in `seq 30 5 330`
#do
#    echo "---------------------- SLACK AWARE $x"
#    export NUMTHREADS=112
#    export POWERLIMIT=$x
#    for i in `seq 1 10`
#    do
#        if timeout 120s ./run -s slackaware lu 6 1400
#        then
#            break
#        fi
#    done
#    mv scheduler.log lu.slackaware.$x.log
#done
#mv *.log lu-results

echo "CG"
#mkdir cg-results
for x in `seq 32 1 112`
do
    mv cg.poweraware.$x.log cg-results/cg.fair.$x.log
done
for x in `seq 31 -1 1`
do
    echo "---------------------- FAIR $x"
    export NUMTHREADS=$x
    tmax=0
    for i in 1 2 3 4 5
    do
        if timeout 120s ./run -s fair cg 192 192 192 6
        then
            newt=`tail -n 2 scheduler.log | head -n 1 | awk -F "[=,]" '{print $2}'`
            echo "tmax: $tmax new: $newt"
            if (( $(echo $newt'>'$tmax | bc -l) )); then
                tmax=$newt
                mv scheduler.log cg.fair.$x.log
            fi
        fi
    done
    #tail -n 5 scheduler.log
done

#for x in `seq 112 -1 1`
#do
#    echo "---------------------- CORE CONSTRAINED $x"
#    export NUMTHREADS=$x
#    for i in `seq 1 10`
#    do
#        if timeout 300s ./run -s coreconstrained cg 192 192 192 6
#        then
#            break
#        fi
#    done
#    mv scheduler.log cg.coreconstrained.$x.log
#done

# stopped at 55
for x in `seq 55 5 160`
do
    echo "---------------------- POWER AWARE $x"
    export NUMTHREADS=112
    export POWERLIMIT=$x
    tmax=1000
    for i in `seq 1 10`
    do
        if timeout 120s ./run -s poweraware cg 192 192 192 6
        then
            newt=`tail -n 2 scheduler.log | head -n 1 | awk -F "[=,]" '{print $2}'`
            echo "tmax: $tmax new: $newt"
            if (( $(echo $newt'<'$tmax | bc -l) )); then
                tmax=$newt
                mv scheduler.log cg.poweraware.$x.log
            fi
        fi
    done
    #mv scheduler.log cg.poweraware.$x.log
done

for x in `seq 30 5 330`
do
    echo "---------------------- SLACK AWARE $x"
    export NUMTHREADS=112
    export POWERLIMIT=$x
    tmax=1000
    for i in `seq 1 10`
    do
        if timeout 120s ./run -s slackaware cg 192 192 192 6
        then
            newt=`tail -n 2 scheduler.log | head -n 1 | awk -F "[=,]" '{print $2}'`
            echo "tmax: $tmax new: $newt"
            if (( $(echo $newt'<'$tmax | bc -l) )); then
                tmax=$newt
                mv scheduler.log cg.slackaware.$x.log
            fi
        fi
    done
    #mv scheduler.log cg.slackaware.$x.log
done
rm scheduler.log
mv *.log cg-results
cd cg-results
./aggregate.sh
cd -
