#!/usr/bin/env python

import numpy as np

for i in range(1,229):
    fname = "cg." + str(i) + ".log"
    with open(fname, 'r') as f:
        avg_pows = {}
        max_pows = {}
        for line in f:
            if "elapsed_seconds" in line:
                d = {}
                for e in line.split(',')[:-1]:
                    [key,val] = e.split('=')
                    d[key] = val
                elapsed = d['elapsed_seconds']
                name = d['name']
                avg_pows[name] = avg_pows.get(name, []) + [float(d['avg_power_W'])]
                max_pows[name] = max_pows.get(name, []) + [float(d['max_power_W'])]

        for name,_ in avg_pows.items():
            std_avg_pow = np.std(avg_pows[name])
            std_max_pow = np.std(max_pows[name])
            with open(name + ".std.csv", 'a') as res:
                res.write(str(i) + "," + str(std_avg_pow) + "," + str(std_max_pow) + "\n")
