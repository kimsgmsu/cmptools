#!/usr/bin/env python
"""
integrate-lineplot: Integrate line plot data
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

def write_line_data_to_file(fname, xdata, ydata, y_int, append):
    """
    Write line data to a file
    """
    # open output file
    try:
        f = open(fname, 'w')
    except IOError:
        util.Error('Failed to open output file "{}"'.format(fname))

    nx = len(xdata)
    f.write('# X Y INT \n')
    if append:
        f.write('# X Y INT \n')
    else:
        f.write('# X INT \n')

    for ix in range(nx):
        if append:
            f.write('{} {} {}\n'.format(xdata[ix], ydata[ix], y_int[ix]))
        else:
            f.write('{} {}\n'.format(xdata[ix], y_int[ix]))

    f.close()

def lower_limit_index(xdata, x1):
    if x1 is None:
        return 0
    else:
        for ix, x in enumerate(xdata):
            ix1 = ix
            if x > x1:
                break
    return ix1
    
def upper_limit_index(xdata, x2):
    if x2 is None:
        return len(xdata)
    else:
        for ix, x in enumerate(xdata):
            ix2 = ix+1
            if x > x2:
                break
    return ix2
    
def integrate_line_data(xdata, ydata, args):

    nx = len(xdata)
    ix1 = lower_limit_index(xdata, args.x1)
    ix2 = upper_limit_index(xdata, args.x2)
    print('*** Integrating line data on grid ({},{}) ***'.format(ix1, ix2))

    y_int1 = np.zeros(ix1)
    y_int2 = integrate.cumtrapz(ydata[ix1:ix2], xdata[ix1:ix2], initial=0)
    y_int3 = np.full(nx-ix2, y_int2[len(y_int2)-1])
    return np.concatenate((y_int1, y_int2, y_int3))
                        
class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('================================')
        print('    INTEGRATE LINE PLOT DATA    ')
        print('================================')

        # read input data
        xdata, ydata = read_line_data_from_file(args.infile, args)

        # generate integrated data
        y_int = integrate_line_data(xdata, ydata, args)

        # plot
        plt.plot(xdata, ydata, 'b-')
        plt.plot(xdata, y_int, 'r-')
        plt.show()

        # write output data
        write_line_data_to_file(args.outfile, xdata, ydata, y_int, args.append)

if __name__ == "__main__":
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='integrate-lineplot: Integrate line plot data')
    parser.add_argument('infile', help='Input data file for DOS')
    parser.add_argument('outfile', help='Output data file for integrated-DOS')
    parser.add_argument('--x1', type=float, default=None,
                        help='lower limit for integration')
    parser.add_argument('--x2', type=float, default=None,
                        help='upper limit for integration')
    parser.add_argument('--xcol', type=int, default=0,
                        help='column for X-data')
    parser.add_argument('--ycol', type=int, default=1,
                        help='column for Y-data')
    parser.add_argument("-a", "--append", help="append the result",
                    action="store_true")

    args = parser.parse_args()

    app = myapp()
    app.run(args)
