#!/bin/bash

if [[ $1 -le 1 ]]; then
    echo "Error: Must launch more than one process"
    exit 1
fi

export I_MPI_MIC=enable
np=`expr $1 - 1`
shift 1
app=run

scp run.mic mic0:

mpirun -iface mic0 -np 1 -host `hostname` ./run.host $@ : -np $np -host mic0 /home/eric/run.mic $@
