#!/usr/bin/env python

"""
xsf2chgcar: Convert xsf format to CHGCAR format
"""

import sys
from ase.io.xsf import read_xsf
from ase import Atoms
import vasp

class xsf2chgcar:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('====================')
        print('   XSF to CHGCAR   ')
        print('====================')

        # open input file
        try:
            infile = open(args.infile, 'r')
            print('Reading input file "{}"'.format(args.infile))
        except IOError:
            print('Failed to open the input file "{}"'.format(args.infile))
            raise

        # read input file
        array, images = read_xsf(infile, read_data=True)

        # open output file
        try:
            outfile = open(args.outfile, 'w')
            print('Opened file "{}" for writing'.format(args.outfile))
        except IOError:
            print('Failed to open the output file "{}"'.format(args.outfile))
            raise

        # write output file
        vasp.write_chgcar(outfile, images[0], array)

if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='xsf2cub: Convert xsf format to CHGCAR format')
    parser.add_argument('infile', help='input file')
    parser.add_argument('outfile', help='output file')
    args = parser.parse_args()
    app = xsf2chgcar()
    app.run(args)
