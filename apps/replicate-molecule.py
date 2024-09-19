#!/usr/bin/env python

"""
replicate-molecule: Replicate molecular structure
"""

import numpy as np
import cmptlib

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('========================')
        print('   REPLICATE MOLECULE   ')
        print('========================')

        # read input file
        mol = cmptlib.readMol(args.infile, fmt=args.infmt,
                              verbose=args.verbose)

        # convert nrep parameters
        try:
            nrep = eval(args.nrep)
        except IOError:
            print('Invalid NREP descriptor')
            raise

        # replicate molecule
        result = cmptlib.replicateMol(mol, nrep)
        result.report()

        # write output file
        cmptlib.writeMol(result, args.outfile, fmt=args.outfmt,
                         coord=args.coord, verbose=args.verbose)

if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='replicate-molecule: Replicate molecular structure')
    parser.add_argument('infile', help='input file')
    parser.add_argument('outfile', help='output file')
    parser.add_argument('--nrep', type=str, default='[2,2,2]',
                        help='No. of cells to replicate in each direction. Given as a string of a valid python int list.  Ex) "[3,1,2]"')
    parser.add_argument('-i', '--infmt',  default=None,
                        help='format for input file')
    parser.add_argument('-o', '--outfmt',  default=None,
                        help='format for output file')
    parser.add_argument('--coord', default=None, choices=['frac', 'cart'],
                        help='coord system for output file')
    parser.add_argument('--verbose', type=int, default=4,
                        help='Verbose level index')

    args = parser.parse_args()

    app = myapp()
    app.run(args)
