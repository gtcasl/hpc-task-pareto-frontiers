#!/usr/bin/env python
import sys

fname = sys.argv[1]

with open(fname, 'r') as fin:
    data = fin.read().splitlines(True)

lines = []
for line in data:
    lines.append(dict(u.split('=') for u in line.rstrip().split(",")[:-1]))

max_threads = 0
cur_threads = 0
for line in lines:
    if "start_task" in line:
        cur_threads = cur_threads + int(line['nthreads'])
        if cur_threads > max_threads:
            max_threads = cur_threads
    if "end_task" in line:
        cur_threads = cur_threads - int(line['nthreads'])

print max_threads
