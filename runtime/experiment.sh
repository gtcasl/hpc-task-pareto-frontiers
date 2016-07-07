#!/bin/bash

nthreads=228
powerlimit=109

#cg
POWERLIMIT=$powerlimit NUMTHREADS=$nthreads ./go.sh 6 cg 150 150 150 2
mv scheduler.log scheduler.cg.log

#cholesky
POWERLIMIT=$powerlimit NUMTHREADS=$nthreads ./go.sh 17 cholesky 4 64
mv scheduler.log scheduler.cholesky.log
