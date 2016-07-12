#!/usr/bin/env python

for i in range(1,229):
    fname = "cholesky." + str(i) + ".log"
    with open(fname, 'r') as f:
        for line in f:
            if "elapsed_seconds" in line:
                d = {}
                for e in line.split(',')[:-1]:
                    [key,val] = e.split('=')
                    d[key] = val
                elapsed = d['elapsed_seconds']
                name = d['name']
                avgpow = d['avg_power_W']
                maxpow = d['max_power_W']
                with open(name + ".csv", 'a') as res:
                    res.write(str(i) + "," + elapsed + "," + avgpow + "," + maxpow + "\n")
