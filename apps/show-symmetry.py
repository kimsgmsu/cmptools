#!/usr/bin/env python
"""
show-symmetry: Show the space group of a crystal structure
"""

import math
import numpy as np
import cmptlib
import Molecule

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('===================')
        print('   SHOW SYMMETRY   ')
        print('===================')

        # read input file
        mol = cmptlib.readMol(args.infile, fmt=args.infmt,
                               verbose=args.verbose)

        # report symmetry
        pcmol = cmptlib.reportSymmetry(mol, args.prec)

        # write primitive cell file
        if args.pcfile is not None:
            cmptlib.writeMol(pcmol, args.pcfile)

if __name__ == "__main__":
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='show-symmetry: Show the space group of a crystal structure')
    parser.add_argument('infile', help='input file')
    parser.add_argument('-i', '--infmt',  default=None,
                        help='format for input file')
    parser.add_argument('--verbose', type=int, default=4,
                        help='Verbose level index')
    parser.add_argument('--prec', type=float, default=1.0e-4,
                        help='Precision for symmetry checking (1e-4)')
    parser.add_argument('--pcfile',
                        help='Output file for primitive cell')

    args = parser.parse_args()

    app = myapp()
    app.run(args)
