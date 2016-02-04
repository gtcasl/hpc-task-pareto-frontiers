#!/bin/bash

# available frequencies: 600000 650000 700000 800000 900000 1000000 1100000

for x in `seq 0 227`
do 
    echo $1 > /sys/devices/system/cpu/cpu$x/cpufreq/scaling_setspeed
done
