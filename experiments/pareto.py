#!/usr/bin/env python
'''
Method to take two equally-sized lists and return just the elements which lie 
on the Pareto frontier, sorted into order.
Default behaviour is to find the maximum for both X and Y, but the option is
available to specify maxX = False or maxY = False to find the minimum for either
or both of the parameters.
'''
def pareto_frontier(Xs, Ys, maxX = True, maxY = True):
    # Sort the list in either ascending or descending order of X
    myList = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse=maxX)
# Start the Pareto frontier with the first value in the sorted list
    p_front = [myList[0]]    
# Loop through the sorted list
    for pair in myList[1:]:
        if maxY: 
            if pair[1] >= p_front[-1][1]: # Look for higher values of Y
                p_front.append(pair) #  and add them to the Pareto frontier
        else:
            if pair[1] <= p_front[-1][1]: # Look for lower values of Y
                p_front.append(pair) #  and add them to the Pareto frontier
# Turn resulting pairs back into a list of Xs and Ys
    p_frontX = [pair[0] for pair in p_front]
    p_frontY = [pair[1] for pair in p_front]
    return p_frontX, p_frontY


if __name__ == "__main__":
    import sys
    import csv

    speedup = [] # third column
    power = [] # fourth column
    with open(sys.argv[1], "r") as csvfile:
        csvreader = csv.reader(csvfile)
        first_line = csvreader.next()
        pow_col = first_line.index("power")
        speedup_col = first_line.index("speedup")
        for row in csvreader:
            speedup.append(float(row[speedup_col]))
            power.append(float(row[pow_col]))
    frontier = pareto_frontier(power, speedup, maxX=False, maxY=True)
    print "power,speedup"
    for p,t in zip(frontier[0],frontier[1]):
        print "%s,%s" % (p,t)
