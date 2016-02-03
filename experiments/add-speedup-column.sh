#!/bin/bash

paste -d "," $1 <( echo "speedup" ; tail -n +2 <$1 | awk -F "," 'NR == 1 {t=$3}; {print t/$3}' ) > tmpfile
mv tmpfile $1
