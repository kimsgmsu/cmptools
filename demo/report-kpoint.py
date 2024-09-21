#!/usr/bin/env python

"""
report-kpoint: 
"""

import cmptlib

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('======================')
        print('    REPORT KPOINT     ')
        print('======================')

        # read input molecule
        mol = cmptlib.readMol(args.infile, fmt=args.infmt,
                              verbose=args.verbose)

        mol.lattice.report()

if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='report-kpoint: Report on kpoints')
    parser.add_argument('infile', help='input file')
    parser.add_argument('-i', '--infmt',  default=None,
                        help='format for input file')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbose level index')

    args = parser.parse_args()

    app = myapp()
    app.run(args)
