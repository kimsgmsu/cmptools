#!/usr/bin/env python
"""
generate-slab: Generate a slab of atoms
"""

import math
import numpy as np
import cmptlib
import Molecule, Atom, Lattice
from pymatgen.core.surface import Slab, SlabGenerator

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('===================')
        print('   GENERATE SLAB   ')
        print('===================')

        # read input file
        mol0 = cmptlib.readMol(args.infile, verbose=args.verbose)

        # convert miller index
        try:
            miller = eval(args.miller_index)
        except IOError:
            print('Invalid Miller index')
            raise

        mol1 = cmptlib.generateSlab(mol=mol0, miller=miller,
                                    min_slab_size=args.min_slab_size,
                                    min_vacuum_size=args.min_vacuum_size,
                                    lll_reduce=args.lll_reduce,
                                    center_slab=args.center_slab,
                                    in_unit_planes=args.in_unit_planes,
                                    primitive=args.primitive,
                                    max_normal_search=args.max_normal_search,
                                    reorient_lattice=args.reorient_lattice)

        # output slab
        cmptlib.writeMol(mol1, args.outfile, verbose=args.verbose)

if __name__ == "__main__":
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='generate-slab: Generate a slab of atoms')
    parser.add_argument('infile', help='initial structure file')
    parser.add_argument('outfile', help='output file')
    parser.add_argument('--miller_index', type=str, default='[1, 0, 0]',
                        help='Miller index of plane parallel to surface. A string of a valid python int list.  Ex) "[1,0,0]"')
    parser.add_argument('--min_slab_size', type=float, default=10.0,
                        help='Minimum size of layers containing atoms. (default=10.0)')
    parser.add_argument('--min_vacuum_size', type=float, default=10.0,
                        help='Minimum size of layers containing vacuum. (default=10.0)')
    parser.add_argument('--lll_reduce', type=bool, default=False,
                        help='Flag to perform LLL reduction to orthogonalize the slab')
    parser.add_argument('--center_slab', type=bool, default=False,
                        help='Flag to center the slab')
    parser.add_argument('--in_unit_planes', type=bool, default=False,
                        help='Lengths in units of hkl plane')
    parser.add_argument('--primitive', type=bool, default=False,
                        help='Flag to reduce the slab to a primitive cell')
    parser.add_argument('--max_normal_search', type=int, default=0,
                        help='Number of searches for a normal lattice vector')
    parser.add_argument('--reorient_lattice', type=bool, default=False,
                        help='Flag to reorient the lattice vectors such that the c-direction is the third vector of the lattice matrix')
    parser.add_argument('--verbose', type=int, default=4,
                        help='Verbose level index')

    args = parser.parse_args()

    app = myapp()
    app.run(args)
