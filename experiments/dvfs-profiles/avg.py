#!/usr/bin/env python

import sys

for freq in [600000,650000,700000,800000,900000,1000000,1100000]:
    for app in ['cg', 'cholesky', 'lu']:
        # results = {name : {num, string} }
        results = {}

        for i in range(1,113):
            arg = app + 'profile.' + str(freq) + '.' + str(i) + '.log'
            with open(arg, 'r') as f:
                sum_avg_pows = {}
                sum_max_pows = {}
                sum_avg_times = {}
                count = {}
                for line in f:
                    if "elapsed_seconds" in line:
                        d = {}
                        for e in line.split(',')[:-1]:
                            [key,val] = e.split('=')
                            d[key] = val
                        elapsed = d['elapsed_seconds']
                        name = d['name']
                        sum_avg_pows[name] = sum_avg_pows.get(name, 0) + float(d['avg_power_W'])
                        sum_max_pows[name] = sum_max_pows.get(name, 0) + float(d['max_power_W'])
                        sum_avg_times[name] = sum_avg_times.get(name, 0) + float(d['avg_time_s'])
                        count[name] = count.get(name, 0) + 1

                for name, n in count.items():
                    avg_pow = sum_avg_pows[name] / n
                    max_pow = sum_max_pows[name] / n
                    avg_time = sum_avg_times[name] / n
                    line = str(i) + "," + str(avg_pow) + "," + str(max_pow) + "," + str(avg_time) + "\n"
                    results[name] = results.get(name,{})
                    results[name][i] = line

        for name,res in results.items():
            with open(name + '.' + str(freq) + ".avg.csv", 'a') as resfile:
                for s in sorted(res):
                    resfile.write(res[s])
                        
