#!/usr/bin/env python
"""
split-voldata: Split volumetric data into separate molecules
"""

import os
import cmptlib

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('======================')
        print('    SPLIT VOL DATA    ')
        print('======================')

        # read input file
        mol0 = cmptlib.readMol(args.infile, fmt=args.infmt,
                               verbose=args.verbose)

        # split volumetric data
        mol_list = mol0.split_voldata()

        # write output file
        if args.rootname is None:
            rootname = os.path.splitext(args.infile)[0]
        else:
            rootname = args.rootname

        outfmt = args.outfmt
        if outfmt is None:
            (base,ext) = os.path.splitext(args.infile)
            outfmt = ext[1:]
            if outfmt == '':
                outfmt = base.lower()

        for i, mol in enumerate(mol_list):
            fname = rootname + '.voldata.{}.{}'.format(i,outfmt)
            cmptlib.writeMol(mol, fname, fmt=args.outfmt,
                             coord=args.coord, verbose=args.verbose)

if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='split-voldata: Split a msf file to separate volumetric data files')
    parser.add_argument('infile', help='input file')
    parser.add_argument('-i', '--infmt',  default=None,
                        help='format for input file')
    parser.add_argument('--rootname',  default=None,
                        help='root name for output files')
    parser.add_argument('-o', '--outfmt',  default=None,
                        help='format for output files')
    parser.add_argument('--coord', default=None, choices=['frac', 'cart'],
                        help='coord system for output files')
    parser.add_argument('--verbose', type=int, default=4,
                        help='Verbose level index')

    args = parser.parse_args()

    app = myapp()
    app.run(args)
