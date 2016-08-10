#!/bin/bash

#cholesky
cmd(){
    ret=1
    for i in 1 2 3 4 5
    do
        if timeout 120s ./run -s advanced cholesky 6 1400 ; then
            ret=0
            break
        fi
    done
    return $ret
}

min_time=10000000
for i in `seq 1 10`
do
    if cmd ; then
        # succeeded
        grep time_s scheduler.log
        new_time=`grep time_s scheduler.log | awk -F "[=,]" '{print $2}'`
        best=`echo $new_time'<'$min_time | bc -l`
        if [ $best -eq 1 ] ; then
            cp scheduler.log cholesky.log
            min_time=$new_time
        fi
    fi
done
