#!/bin/bash

# Note: to set the maximum usable threads, set environment variable 'NUMTHREADS'

if [[ $1 -le 1 ]]; then
    echo "Error: Must launch more than one process"
    exit 1
fi

export I_MPI_MIC=enable
export DAPL_DBG_TYPE=0
np=`expr $1 - 1`
shift 1
app=run

scp run.mic mic0:

mpirun -iface mic0 -np 1 -host `hostname` ./run.host $@ : -np $np -host mic0 /home/eric/run.mic $@
