#!/usr/bin/env python

"""
print-nbrmap: Print neighbor map of a molecule
"""

import numpy as np
import cmptlib

class Neighbor:
    """
    class for a neighboring atom.
    """
    def __init__(self, atom, cell, dist, cpos):
        self.atom = atom
        self.cell = cell
        self.dist = dist
        self.cpos = cpos

class print_nbrmap:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('==================')
        print('   PRINT NBRMAP   ')
        print('==================')

        # read input file
        mol = cmptlib.readMol(args.infile, fmt=args.infmt,
                              verbose=args.verbose)

        # convert padding
        try:
            padding = eval(args.padding)
        except IOError:
            print('Invalid padding descriptor')
            raise

        # convert onlyatom list
        try:
            onlyatoms = eval(args.onlyatoms)
        except IOError:
            print('Invalid ONLYATOMS descriptor')
            raise

        # print a neighbor map
        cmptlib.printNeighborMap(mol, args.cutoff, padding, onlyatoms,
                                 args.tol, only_closer_dist=args.only_closer_dist)

if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='print-nbrmap: Print a neighbor map.  Ex) print-nbrmap.py GdN.msf 10.0 --padding="[3,2,2]" --onlyatoms="[1,18]" > GdN.nbrmap.out')
    parser.add_argument('infile', help='input file')
    parser.add_argument('cutoff', type=float, default=5.0,
                        help='cutoff radius')
    parser.add_argument('-i', '--infmt',  default=None,
                        help='format for input file')
    parser.add_argument('--padding', type=str, default='[1, 1, 1]',
                        help='No. of cells to pad in each direction. Given as a string of a valid python int list.  Ex) "[3,3,3]"')
    parser.add_argument('--onlyatoms', type=str, default='[]',
                        help='List of only atom indices to check neighbors. Given as a string of a valid python int list.  Ex) "[1,2,8,10]".  "[]"=ALL')
    parser.add_argument('--dtol', type=float, default=1.0e-4,
                        help='tolerance for distance to be considered as zero')
    parser.add_argument('--only_closer_dist', type=float, default=None,
                        help='print atoms only if neighbors are closer than this value')
    parser.add_argument('--verbose', type=int, default=4,
                        help='Verbose level index')

    args = parser.parse_args()

    app = print_nbrmap()
    app.run(args)
