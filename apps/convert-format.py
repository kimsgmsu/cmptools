#!/usr/bin/env python

"""
convert-format: Convert a molecular structure file to a different format
"""

import cmptlib

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('====================')
        print('   CONVERT FORMAT   ')
        print('====================')

        # read input file
        mol = cmptlib.readMol(args.infile, fmt=args.infmt,
                              verbose=args.verbose)

        # write output file
        cmptlib.writeMol(mol, args.outfile, fmt=args.outfmt,
                         coord=args.coord, verbose=args.verbose)

if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='convert-format: Convert a molecular structure file format')
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
