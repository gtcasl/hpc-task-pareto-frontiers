#!/bin/bash

app=cg

for x in coreconstrained fair 
do
    for y in `seq 1 1 112`
    do
        ../scripts/get-num-threads.py $app.$x.$y.log >> $app.$x.threads
    done
done
for x in poweraware slackaware
do
    for y in `seq 30 5 330`
    do
        ../scripts/get-num-threads.py $app.$x.$y.log >> $app.$x.threads
    done
done
