#!/usr/bin/env python

"""
print-displacement: Print displacement of atoms between two molecules
"""

import numpy as np
import cmptlib

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('========================')
        print('   PRINT DISPLACEMENT   ')
        print('========================')

        # read molecules
        mol1 = cmptlib.readMol(args.mol1, verbose=args.verbose)
        mol2 = cmptlib.readMol(args.mol2, verbose=args.verbose)

        # report displacement
        cmptlib.printDisplacement(mol1, mol2, args.cutoff)

if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='print-displacement: Print displacement of atoms between two molecules.  Ex) print-displacement.py mol1.msf mol2.msf --cutoff=2.0 > displacement.out')
    parser.add_argument('mol1', help='molecule 1')
    parser.add_argument('mol2', help='molecule 2')
    parser.add_argument('--cutoff', type=float, default=None,
                        help='cutoff displacement. Displacement less than this will not be printed')
    parser.add_argument('--verbose', type=int, default=4,
                        help='Verbose level index')

    args = parser.parse_args()

    app = myapp()
    app.run(args)
