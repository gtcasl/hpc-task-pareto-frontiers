#!/usr/bin/env python

for i in range(1,229):
    fname = "cg." + str(i) + ".log"
    with open(fname, 'r') as f:
        sum_avg_pows = {}
        sum_max_pows = {}
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
                count[name] = count.get(name, 0) + 1

        for name, n in count.items():
            avg_pow = sum_avg_pows[name] / n
            max_pow = sum_max_pows[name] / n
            with open(name + ".avg.csv", 'a') as res:
                res.write(str(i) + "," + str(avg_pow) + "," + str(max_pow) + "\n")
