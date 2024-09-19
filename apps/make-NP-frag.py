#!/usr/bin/env python

"""
make-NP: Draw spheres on atoms.
"""

import numpy as np
import math
from pymatgen.core.operations import SymmOp
import cmptlib, util, Lattice, Plane

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('==================')
        print('   MAKE NP-part   ')
        print('==================')

        # convert slice range
        if args.xslice == 'None':
            xslice = None
        else:
            try:
                xslice = eval(args.xslice)
            except IOError:
                print('Invalid xslice descriptor')
                raise
        print('xslice = ', xslice)

        # read the unit cell
        unitcell = cmptlib.readMol(args.infile, fmt=args.infmt,
                                   verbose=args.verbose)

        # make a nano-particle
        mol = make_NP(unitcell, args.radius, args.dshell,
                      args.dfix, args.nocenter, args.theta0, args.phi0,
                      args.dalpha, args.dalpha_fix, xslice,
                      args.reorient, args.boundingbox, args.padding)

        # write intermediate file
        cmptlib.writeMol(mol, 'tmp-NP.vasp')

        """
        """

        # report the molecule
        mol.report()

        # write output file
        cmptlib.writeMol(mol, args.outfile, fmt=args.outfmt,
                         coord=args.coord, verbose=args.verbose)

def make_NP(mol0, radius=10.0, dshell=-1.0, dfix=-1.0, nocenter=False,
            theta0=0.0, phi0=0.0, dalpha=math.pi, dalpha_fix=0.0, xslice=None,
            reorient=True, boundingbox=True, padding=10.0):

    # set radii
    r3 = max(radius, 0.0)
    if dshell < 0:
        r1 = 0.0
    else:
        r1 = max(r3 - dshell, 0.0)
    if dfix < 0:
        r2 = r1
    else:
        r2 = min(r1 + dfix, r3)
    print('r1={}, r2={}, r3={}'.format(r1, r2, r3))

    # set angles
    if dalpha < 0:
        alpha1 = 200
    else:
        alpha1 = min(max(dalpha, 0.0), 180)
    if dalpha_fix < 0:
        alpha2 = alpha1 + 10
    else:
        alpha2 = alpha1 + dalpha_fix
    print('alpha1={}, alpha2={}'.format(alpha1, alpha2))

    if not nocenter:
        cmptlib.centerMol(mol0)

    # determine the cell range
    cell_range = set_cell_range(mol0, radius)
    print('cell_range =', cell_range)

    # set molecule
    mol1 = mol0.copy()
    # lattice
    matrix = np.eye(3)
    for idir in range(3):
        matrix[idir,idir] = radius*2 + padding
    mol1.lattice = Lattice.Lattice(matrix)

    cobox = np.zeros(3)
    for idir in range(3):
        cobox += mol1.lattice.matrix[idir]*0.5
    print('cobox =', cobox)
    dpos0 = util.asCartesian([r3, theta0, phi0])
    print('dpos0 = ', dpos0)

    # generate atoms
    mol1.atomList = []
    for c0 in range(-cell_range[0], cell_range[0]+1):
        for c1 in range(-cell_range[1], cell_range[1]+1):
            for c2 in range(-cell_range[2], cell_range[2]+1):
                dpos = mol0.to_cartesian(np.array([c0, c1, c2]))
                for atom0 in mol0.atomList:
                    atom1 = atom0.copy()
                    atom1.pos = cobox + atom0.pos + dpos
                    dpos1 = atom1.pos - cobox
                    r = np.linalg.norm(dpos1)
                    # r
                    if r < r1 or r > r3:
                        continue
                    if r < r2:
                        atom1.symbol = 'Fx'
                    # alpha
                    alpha = util.angle_between(dpos1, dpos0)*180/math.pi
                    if alpha < alpha1:
                        pass
                    elif alpha > alpha2:
                        continue
                    else:
                        atom1.symbol = 'Fx'
                    # xslice
                    if xslice is None:
                        pass
                    else:
                        if dpos1[0] < xslice[0] or dpos1[0] > xslice[1]:
                            continue
                    mol1.atomList.append(atom1)

    # reorient along c-axis
    if reorient:
        # Symm op. for rotating the molecule
        origin = cobox
        axis = [1, 0, 0]
        angle = theta0
        op = SymmOp.from_origin_axis_angle(origin, axis, angle)
        mol1.apply_operation(op)
        
    # make a bounding box
    if boundingbox:
        mol1.make_bounding_box(np.full((3,2),padding))
        
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
    parser.add_argument('--dshell', type=float, default=-1.0,
                        help='Thickness of the shell (Neg = full radius), default=-1.0')
    parser.add_argument('--dfix', type=float, default=-1.0,
                        help='Thickness of the fixed shell (Neg = no fixed shell), default=-1.0')
    parser.add_argument('--nocenter', action='store_true',
                        help='Flag not to center the molecule')
    parser.add_argument('--theta0', type=float, default=0.0,
                        help='Theta0 (deg), default=0.0')
    parser.add_argument('--phi0', type=float, default=0.0,
                        help='Phi0 (deg), default=0.0')
    parser.add_argument('--dalpha', type=float, default=-1.0,
                        help='dalpha (deg) (Neg = no slice), default=-1.0')
    parser.add_argument('--dalpha_fix', type=float, default=-1.0,
                        help='dalpha_fix (deg) (Neg = no fixed slice), default=0.0')
    parser.add_argument('--reorient', action='store_true',
                        help='Flag to reorient the molecule along c-axis')
    parser.add_argument('--boundingbox', action='store_true',
                        help='Flag to put the final molecule in a bounding box')
    parser.add_argument('--padding', type=float, default=10.0,
                        help='Thickness of padding around the bounding box, default=10.0')
    parser.add_argument('--xslice', type=str, default='None',
                        help='Range for X slice. Ex) "[-1.5, 2.3]"')
    parser.add_argument('--coord', default=None, choices=['frac', 'cart'],
                        help='coord system for output file')
    parser.add_argument('--verbose', type=int, default=3,
                        help='Verbose level index')

    args = parser.parse_args()

    app = myapp()
    app.run(args)
