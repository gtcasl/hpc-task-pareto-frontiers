#!/bin/bash

paste -d "," $1 <( echo "power" ; tail -n +2 <$1 | awk -F "," '{print $2/$3}' ) > tmpfile
mv tmpfile $1
