#!/bin/bash

nthreads=228

#cg
NUMTHREADS=$nthreads ./go.sh 6 cg 150 150 150 2
mv schedule.log schedule.cg.log

#cholesky
NUMTHREADS=$nthreads ./go.sh 17 cholesky 4 64
mv schedule.log schedule.cholesky.log
