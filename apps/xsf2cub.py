#!/usr/bin/env python

"""
xsf2cub: Convert xsf format to cube format
"""

import sys
from ase.io.xsf import read_xsf
from ase.io.cube import write_cube
from ase        import Atoms

class xsf2cub:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('====================')
        print('   XSF to CUB   ')
        print('====================')

        # open input file
        try:
            infile = open(args.infile, 'r')
            print('Reading input file "{}"'.format(args.infile))
        except IOError:
            print('Failed to open the input file "{}"'.format(args.infile))
            raise

        # read input file
        xsf = read_xsf(infile, read_data=True)

        # open output file
        try:
            outfile = open(args.outfile, 'w')
            print('Opened file "{}" for writing'.format(args.outfile))
        except IOError:
            print('Failed to open the output file "{}"'.format(args.outfile))
            raise

        # write output file
        write_cube(outfile, xsf[1], xsf[0])

if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='xsf2cub: Convert xsf format to cube format')
    parser.add_argument('infile', help='input file')
    parser.add_argument('outfile', help='output file')
    args = parser.parse_args()
    app = xsf2cub()
    app.run(args)
