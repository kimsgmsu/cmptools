#!/usr/bin/env python
"""
derivative-lineplot: Differentiate line plot data
"""

import os
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

def read_line_data_from_file(fname, args):
    """
    Read line data from a file
    """
    # open input file
    try:
        f = open(fname, 'r')
    except IOError:
        util.Error('Failed to open input file "{}"'.format(fname))

    xcol = args.xcol
    ycol = args.ycol
    xlist = []
    ylist = []
    for line in f:
        line = line.strip()

        # skip empty line
        if len(line) < 1:
            continue
    
        # remove comment
        if line[0:1] == "#":
            continue

        words = line.split()
        x = float(words[xcol])
        y = float(words[ycol])
        xlist.append(x)
        ylist.append(y)

    f.close()
    return np.array(xlist), np.array(ylist)

def write_line_data_to_file(fname, xdata, ydata, dydx, append):
    """
    Write line data to a file
    """
    # open output file
    try:
        f = open(fname, 'w')
    except IOError:
        util.Error('Failed to open output file "{}"'.format(fname))

    nx = len(xdata)
    if append:
        f.write('# X Y dYdX \n')
    else:
        f.write('# X dYdX \n')

    for ix in range(nx):
        if append:
            f.write('{} {} {}\n'.format(xdata[ix], ydata[ix], dydx[ix]))
        else:
            f.write('{} {}\n'.format(xdata[ix], dydx[ix]))

    f.close()

def derivative_line_data(xdata, ydata, args):

    nx = len(xdata)
    print('*** Differentiating line data on {} grid ***'.format(nx))

    dydx = np.zeros(nx)
    for i in range(nx):
        if i > 0:
            h1 = xdata[i] - xdata[i-1]
            y1 = ydata[i-1]
        else:
            h1 = xdata[i+1] - xdata[i]
            y1 = ydata[i+1] - 2*(ydata[i+1]-ydata[i])
        if i < (nx-1):
            h2 = h1
            y2 = ydata[i+1]
        else:
            h2 = h1
            y2 = ydata[i-1] - 2*(ydata[i-1]-ydata[i])
        dydx[i] = (y2 - y1)/(h1 + h2)

    return dydx

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('=================================')
        print('    DERIVATIVE LINE PLOT DATA    ')
        print('=================================')

        # read input data
        xdata, ydata = read_line_data_from_file(args.infile, args)

        # generate differentiated data
        dydx = derivative_line_data(xdata, ydata, args)

        # plot
        plt.plot(xdata, ydata, 'b-')
        plt.plot(xdata, dydx, 'r-')
        plt.show()

        # write output data
        write_line_data_to_file(args.outfile, xdata, ydata, dydx, args.append)

if __name__ == "__main__":
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='derivative-lineplot: Differentiate line plot data')
    parser.add_argument('infile', help='Input file for line data')
    parser.add_argument('outfile', help='Output file for differentiated data')
    parser.add_argument('--xcol', type=int, default=0,
                        help='column for X-data')
    parser.add_argument('--ycol', type=int, default=1,
                        help='column for Y-data')
    parser.add_argument("-a", "--append", help="append the result",
                    action="store_true")
    args = parser.parse_args()

    app = myapp()
    app.run(args)
