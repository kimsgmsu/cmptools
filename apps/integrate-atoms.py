#!/usr/bin/env python

"""
integrate-atoms: Integrate volumetric data on atoms.
"""

import numpy as np
import util, cmptlib, Molecule
import Plane
import math

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('=====================')
        print('   INTEGRATE ATOMS   ')
        print('=====================')

        # read volumetric data file
        mol_vdata = cmptlib.readMol(args.vdata)

        # read atoms
        mol_atoms = cmptlib.readMol(args.atoms)
        mol_atoms.report('Atoms to integrate:')

        # integrate data on atoms
        for kvd, vdata in enumerate(mol_vdata.volumeDataList):
            integrate_atoms(vdata, mol_vdata, mol_atoms)

def integrate_atoms(vdata, mol_vdata, mol_atoms):

    from scipy.interpolate import RegularGridInterpolator

    ndim = vdata.data.shape
    nx = ndim[0]
    ny = ndim[1]
    nz = ndim[2]
    xlist = np.linspace(0.0, 1.0, nx)
    ylist = np.linspace(0.0, 1.0, ny)
    zlist = np.linspace(0.0, 1.0, nz)
#    my_func = RegularGridInterpolator((xlist, ylist, zlist), vdata.data)
    hx = 1.0/np.float_(nx-1)
    hy = 1.0/np.float_(ny-1)
    hz = 1.0/np.float_(nz-1)
    sfac = hx*hy*hz

    print('*** Integrating volume data on atoms: {} on ({},{},{}) grid ***'.format(
        vdata.name, nx, ny, nz))
    lbar = 80
    print('-'*lbar)
    print('  ID TYPE (             POSITION            )', end='')
    print('  INTEGRAL  ')
    print('-'*lbar)
    Isum = 0
    for ka, atom in enumerate(mol_atoms.atomList):
        radius = mol_atoms.get_Element(atom.symbol).radius
        center = atom.pos
        bbox = get_bounding_box(center, radius, mol_vdata, ndim)
        Ival = 0
        for ix in range(bbox[0][0], bbox[0][1]):
            x = (ix+0.5)*hx
            for iy in range(bbox[1][0], bbox[1][1]):
                y = (iy+0.5)*hy
                for iz in range(bbox[2][0], bbox[2][1]):
                    z = (iz+0.5)*hz
                    fpos = np.array([x, y, z])
                    cpos = mol_vdata.to_cartesian(fpos)
                    dist = np.linalg.norm(cpos - center)
                    if dist < radius:
#                        rpos = util.reduced_vector2(fpos)
#                        fval = my_func(rpos)[0]
                        idx = [ix, iy, iz]
                        ridx = util.reduced_index(idx, ndim)
                        fval = vdata.data[ridx[0],ridx[1],ridx[2]]
                        Ival += fval
        Ival *= sfac
        Isum += Ival
        print('{:>4}'.format(ka+1), end='')
        print(' {:>4}'.format(atom.symbol), end='')
        print(' (', end='')
        for c in atom.pos:
            print(' {:>10.6f}'.format(c), end='')
        print(') ', end='')
        print(' {}'.format(Ival))
    print('-'*lbar)
    print('  SUM:', ' '*38, ' {}'.format(Isum))
    print('-'*lbar)
        
def get_bounding_box(center, radius, mol, ndim):
    bbox = []
    for idir in range(3):
        rng1 = get_bounding_range(idir, center, radius, mol, ndim)
        bbox.append(rng1)
    return bbox

def get_bounding_range(kdir, center, radius, mol, ndim):
    miller = [ 0, 0, 0]
    miller[kdir] = 1
    plane = Plane.Plane.from_miller_and_pos(miller, mol.lattice.matrix, center)
    pos1 = center - plane.normal*radius
    pos2 = center + plane.normal*radius
    fpos0 = mol.to_fractional(center)
    fpos1 = mol.to_fractional(pos1)
    fpos2 = mol.to_fractional(pos2)
    ip0 = math.floor(fpos0[kdir]*ndim[kdir])
    ip1 = math.floor(fpos1[kdir]*ndim[kdir])
    ip2 = math.ceil(fpos2[kdir]*ndim[kdir])
    kp1 = min(ip1, ip2) - 1
    kp2 = max(ip1, ip2) + 1
    return (kp1, kp2)

if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='integrate-atoms: Integrate volumetric data on atoms.  Ex) integrate-atoms.py Gd2C.chgden.vasp Gd2C.spheres.msf > out.dat')
    parser.add_argument('vdata', help='volumetric data file')
    parser.add_argument('atoms', help='atoms file')
    parser.add_argument('--verbose', type=int, default=3,
                        help='Verbose level index')

    args = parser.parse_args()

    app = myapp()
    app.run(args)
