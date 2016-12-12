#!/usr/bin/env python
import sys

fname = sys.argv[1]

with open(fname, 'r') as fin:
    data = fin.read().splitlines(True)

lines = []
for line in data:
    lines.append(dict(u.split('=') for u in line.rstrip().split(",")[:-1]))

state = 'start'
start_time = -1
time = 0
tot_time = 0
sum_energy = 0
for line in lines:
    if state == 'start' and 'start_task' in line:
        time = float(line['start_time'])
        state = 'stats'
        if start_time == -1:
            start_time = time
    elif state == 'stats' and 'message' in line and line['message'] == 'tick stats':
        power = float(line['estimated_current_power'])
        sum_energy = sum_energy + (time - start_time) * power
        state = 'start'
        start_time = time
    if 'time_s' in line:
        tot_time = float(line['time_s'])


print sum_energy / tot_time
