#!/usr/bin/env python

"""
report-molecule: Report on a molecule
"""

import cmptlib

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('======================')
        print('    REPORT MOLECULE   ')
        print('======================')

        # read input file
        mol = cmptlib.readMol(args.infile, fmt=args.infmt, verbose=args.verbose)

        # report the structure
        mol.report('==== MOLECULE REPORT ====')

if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='report-molecule: Report on a molecule')
    parser.add_argument('infile', help='input file')
    parser.add_argument('-i', '--infmt',  default=None,
                        help='format for input file')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbose level index')

    args = parser.parse_args()

    app = myapp()
    app.run(args)
