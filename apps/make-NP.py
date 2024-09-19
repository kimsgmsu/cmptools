#!/usr/bin/env python

"""
make-NP: Draw spheres on atoms.
"""

import numpy as np
import math
import cmptlib, Lattice, Plane

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('=============')
        print('   MAKE NP   ')
        print('=============')

        # read the unit cell
        unitcell = cmptlib.readMol(args.infile, fmt=args.infmt,
                                   verbose=args.verbose)

        # make a NP
        molNP = make_NP(unitcell, args.radius, args.dshell,
                        args.dfix, args.center)

        # report the molecule
        molNP.report()

        # write output file
        cmptlib.writeMol(molNP, args.outfile, fmt=args.outfmt,
                         coord=args.coord, verbose=args.verbose)

def make_NP(mol0, radius=10.0, dshell=6.0, dfix=3.0, center=True):

    # set radii
    r3 = max(radius, 0.0)
    r1 = min(max(r3 - dshell, 0.0), r3)
    r2 = min(max(r1 + dfix, r1), r3)
    print('r1={}, r2={}, r3={}'.format(r1, r2, r3))

    if center:
        cmptlib.centerMol(mol0)

    # determine the cell range
    cell_range = set_cell_range(mol0, radius)

    # set molecule
    mol1 = mol0.copy()
    # lattice
    cfactor = 2.0
    cmin = 30.0
    cmax = np.amax(np.array(cell_range))
    matrix = np.eye(3)
    for idir in range(3):
        matrix[idir,idir] = cmax*2*cfactor+cmin
    mol1.lattice = Lattice.Lattice(matrix)

    cobox = np.zeros(3)
    for idir in range(3):
        cobox += mol1.lattice.matrix[idir]*0.5
    print('cobox =', cobox)

    mol1.atomList = []
    for c0 in range(-cell_range[0], cell_range[0]+1):
        for c1 in range(-cell_range[1], cell_range[1]+1):
            for c2 in range(-cell_range[2], cell_range[2]+1):
                dpos = mol0.to_cartesian(np.array([c0, c1, c2]))
                for atom0 in mol0.atomList:
                    atom1 = atom0.copy()
                    atom1.pos = cobox + atom0.pos + dpos
                    dist = atom1.distance(cobox, coord='cart', latt=mol1.lattice)
                    if dist < r1 or dist > r3:
                        continue
                    if dist < r2:
                        atom1.symbol = 'Fx'
                    mol1.atomList.append(atom1)

    return mol1

def set_cell_range(mol, radius):
    cell_range = [0, 0, 0]
    for idir in range(3):
        miller = np.zeros(3)
        miller[idir] = 1.0
        point = mol.to_cartesian(miller)
        plane = Plane.Plane.from_miller_and_point(miller, mol.lattice.matrix, point)
        cell_range[idir] = math.ceil(radius/math.fabs(plane.d0))
        print('idir={}: miller={}, radius={}, d0={}, cell_range={}'.format(
            idir, miller, radius, plane.d0, cell_range[idir]))
    return cell_range

if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='make-NP: Make a nano-particle.  Ex) make-NP.py Cu-111.bulk.msf Cu-NP.vasp --radius=12.0 --dshell=6.0 --dfix=3.0 --center')
    parser.add_argument('infile', help='input file')
    parser.add_argument('outfile', help='output file')
    parser.add_argument('-i', '--infmt',  default=None,
                        help='format for input file')
    parser.add_argument('-o', '--outfmt',  default=None,
                        help='format for output file')
    parser.add_argument('--radius', type=float, default=10.0,
                        help='Radius of the particle, default=10.0')
    parser.add_argument('--dshell', type=float, default=6.0,
                        help='Thickness of the shell, default=6.0')
    parser.add_argument('--dfix', type=float, default=3.0,
                        help='Thickness of the fixed shell, default=3.0')
    parser.add_argument('--center', action='store_true',
                        help='Flag to center the molecule')
    parser.add_argument('--coord', default=None, choices=['frac', 'cart'],
                        help='coord system for output file')
    parser.add_argument('--verbose', type=int, default=3,
                        help='Verbose level index')

    args = parser.parse_args()

    app = myapp()
    app.run(args)
