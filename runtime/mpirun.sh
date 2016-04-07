#!/bin/bash

if [[ $1 -le 1 ]]; then
    echo "Error: Must launch more than one process"
    exit 1
fi

export I_MPI_MIC=enable
np=`expr $1 - 1`
shift 1
app=$1
shift 1
mpirun -np 1 -host `hostname` $app.host $@ : -np $np -host mic0 $app.mic $@
