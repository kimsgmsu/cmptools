#!/usr/bin/env python

"""
centerMol: Move the molecule so that the center-of-mass is 
           at the center of the bounding box
"""

import cmptlib

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('=====================')
        print('   CENTER MOLECULE   ')
        print('=====================')

        # read input file
        mol = cmptlib.readMol(args.infile, fmt=args.infmt,
                              verbose=args.verbose)

        # center the molecule
        cmptlib.centerMol(mol)

        # write output file
        cmptlib.writeMol(mol, args.outfile, fmt=args.outfmt,
                         coord=args.coord, verbose=args.verbose)

if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='centerMol: Move a molecule to the center of the bounding box')
    parser.add_argument('infile', help='input file')
    parser.add_argument('outfile', help='output file')
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
