#!/usr/bin/env python

"""
canonizeMol: Canonize a molecular structure file
"""

from __future__ import print_function
import os
import numpy as np
import fileUtil

def std_lattice(mol):
    """
    Standardize the lattice vectors
    """
    pass
    
class canonizeMol:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('=======================')
        print('   CANONIZE MOLECULE   ')
        print('=======================')

        # read input file
        mol = fileUtil.read(args.infile)
        mol.report()

        # lattice vectors
        if args.lattice:
            std_lattice(mol)

        # write output file
        fileUtil.write(mol, args.outfile)

if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='canonizeMol: Canonize a molecular structure file')
    parser.add_argument('infile', help='input file')
    parser.add_argument('outfile', help='output file')
    parser.add_argument('--lattice',
                        help='standardize lattice vectors',
                        action="store_true")
    parser.add_argument('--position',
                        help='standardize stom positions',
                        action="store_true")

    args = parser.parse_args()
    print(args)
    app = canonizeMol()
    app.run(args)
