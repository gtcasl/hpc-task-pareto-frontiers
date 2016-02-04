#!/bin/bash

for freq in 600000 650000 700000 800000 900000 1000000 1100000
do
    ssh mic0 'bash -s' < ./setfreq.sh $freq
    ssh mic0 'cd /home/eric && bash -s' < ./run-spmv-threads.sh
    scp mic0:/home/eric/spmv.csv spmv.$freq.csv
done
