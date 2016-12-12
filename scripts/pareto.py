#!/usr/bin/env python
'''
Method to take two equally-sized lists and return just the elements which lie 
on the Pareto frontier, sorted into order.
Default behaviour is to find the maximum for both X and Y, but the option is
available to specify maxX = False or maxY = False to find the minimum for either
or both of the parameters.
'''
def pareto_frontier(Xs, Ys, indexes, maxX = True, maxY = True):
    # Sort the list in either ascending or descending order of X
    myList = sorted(enumerate([[Xs[i], Ys[i]] for i in range(len(Xs))]),
                    reverse=maxX,
                    key=lambda x:x[1])
    myList = [(int(indexes[i]),j) for (i,j) in myList]
# Start the Pareto frontier with the first value in the sorted list
    p_front = [myList[0]]    
# Loop through the sorted list
    for pair in myList[1:]:
        if maxY: 
            if pair[1][1] >= p_front[-1][1][1]: # Look for higher values of Y
                p_front.append(pair) #  and add them to the Pareto frontier
        else:
            if pair[1][1] <= p_front[-1][1][1]: # Look for lower values of Y
                p_front.append(pair) #  and add them to the Pareto frontier
# Turn resulting pairs back into a list of Xs and Ys
    return [(p[1][0],p[1][1],p[0]) for p in p_front]

if __name__ == "__main__":
    import sys
    import csv

    time = [] # third column
    power = [] # fourth column
    nthreads = [] # first column
    with open(sys.argv[1], "r") as csvfile:
        csvreader = csv.reader(csvfile)
        first_line = csvreader.next()
        pow_col = first_line.index("power")
        time_col = first_line.index("time")
        nthreads_col = first_line.index("nthreads")
        for row in csvreader:
            time.append(float(row[time_col]))
            power.append(float(row[pow_col]) - 106.9)
            nthreads.append(float(row[nthreads_col]))
    speedups = [time[0] / t for t in time]
    frontier = pareto_frontier(power, speedups, nthreads, maxX=False, maxY=True)
    print "power,speedup,nthreads"
    for x in frontier:
        print ",".join([str(i) for i in x])
