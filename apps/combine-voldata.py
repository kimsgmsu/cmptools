#!/usr/bin/env python
"""
combine-voldata: Combine volumetric data of two molecules into one
"""

import os
import numpy as np
import cmptlib

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('=============================')
        print('   COMBINE TWO VOLUME DATA   ')
        print('=============================')

        # read molecules
        mol1 = cmptlib.readMol(args.mol1, verbose=args.verbose)
        mol2 = cmptlib.readMol(args.mol2, verbose=args.verbose)

        # combine two molecule by combining volumetric data
        mol3 = cmptlib.combineVolumeData(mol1, mol2, expr=args.expr)

        # write output molecule
        print('Writing output molecule "{}"'.format(args.outfile))
        cmptlib.writeMol(mol3, args.outfile, coord=args.coord,
                         verbose=args.verbose)

if __name__ == "__main__":
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='combine-voldata: Combine volumetric data of two molecule files')
    parser.add_argument('mol1', help='molecule 1')
    parser.add_argument('mol2', help='molecule 2')
    parser.add_argument('outfile', help='output file')
    parser.add_argument('expr',
                        help='Expression to combine two volumetric data. Must be enclosed in quotation marks. Ex) "a-b", "(a+b)/2"')
    parser.add_argument('--coord', default=None, choices=['frac', 'cart'],
                        help='coord system for output file')
    parser.add_argument('--verbose', type=int, default=4,
                        help='Verbose level index')

    args = parser.parse_args()

    app = myapp()
    app.run(args)
