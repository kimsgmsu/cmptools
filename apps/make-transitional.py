#!/usr/bin/env python
"""
make-transitional: Make transitional molecules
"""

import os
import numpy as np
import cmptlib

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('=================================')
        print('   MAKE TRANSITIONAL MOLECULES   ')
        print('=================================')

        # read molecules
        mol1 = cmptlib.readMol(args.mol1, verbose=args.verbose)
        mol2 = cmptlib.readMol(args.mol2, verbose=args.verbose)

        # make transitional molecules
        mol3 = cmptlib.makeTransitionalMol(mol1, mol2, nmol=args.nmol,
                                           prefix=args.prefix, fmt=args.fmt)

if __name__ == "__main__":
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='make-transitional: Make trasitional molecules')
    parser.add_argument('mol1', help='molecule 1')
    parser.add_argument('mol2', help='molecule 2')
    parser.add_argument('nmol', type=int, help='No. of transitional molecules to make')
    parser.add_argument('--prefix', default='tmp', 
                        help='Prefix for file names for transitional molecules')
    parser.add_argument('--fmt', default='vasp', help='format of the output files')
    parser.add_argument('--verbose', type=int, default=4,
                        help='Verbose level index')

    args = parser.parse_args()

    app = myapp()
    app.run(args)
