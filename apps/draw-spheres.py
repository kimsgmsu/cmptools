#!/usr/bin/env python

"""
draw-spheres: Draw spheres on atoms.
"""

import numpy as np
import util, cmptlib, Molecule
import Plane
import math

class draw_spheres:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('==================')
        print('   DRAW SPHERES   ')
        print('==================')

        # read input file
        mol = cmptlib.readMol(args.infile)
        mol.report('Input molecule:')

        # convert the DIM parameters
        try:
            ndim = util.str_to_list(args.ndim)
        except:
            print('Invalid DIM descriptor')
            raise

        # count spheres
        print('inverse=', args.inverse)
        count_spheres(mol, ndim, args)
        mol.report('Output molecule:')

        # write output file
        cmptlib.writeMol(mol, args.outfile)

def count_spheres(mol, ndim, args):

    if args.inverse:
        data = np.ones( ndim, dtype=np.float_ )
    else:
        data = np.zeros( ndim, dtype=np.float_ )

    nx = ndim[0]
    ny = ndim[1]
    nz = ndim[2]
    hx = 1.0/np.float_(nx-1)
    hy = 1.0/np.float_(ny-1)
    hz = 1.0/np.float_(nz-1)
    for katom, atom in enumerate(mol.atomList):
        util.progress(katom, mol.natom - 1)
        radius = mol.get_Element(atom.symbol).radius
        center = atom.pos
        bbox = get_bounding_box(center, radius, mol, ndim)
        print('Atom #{}: {} at {} r={}'.format(katom, atom.symbol,
                                          mol.to_fractional(center),
                                          radius))
        print('bbox = ', bbox)
        for ix in range(bbox[0][0], bbox[0][1]):
            ix0 = ix % nx
            x = (ix+0.5)*hx
            for iy in range(bbox[1][0], bbox[1][1]):
                iy0 = iy % ny
                y = (iy+0.5)*hy
                for iz in range(bbox[2][0], bbox[2][1]):
                    iz0 = iz % nz
                    z = (iz+0.5)*hz
                    fpos = np.array([x, y, z])
                    cpos = mol.to_cartesian(fpos)
                    dist = np.linalg.norm(cpos - center)
                    if dist < radius:
                        if args.accumulate:
                            data[ix0,iy0,iz0] += 1
                        elif args.inverse:
                            data[ix0,iy0,iz0] = 0
                        else:
                            data[ix0,iy0,iz0] = 1

    vd = Molecule.VolumetricData(name='Overlap index', data=data, extra=None)
    mol.volumeDataList = [vd]
    mol.report('In draw-spheres')

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
    print('bbox: center={}, radius={}, normal={}, pos1={}, pos2={}'.format(
        center, radius, plane.normal, pos1, pos2))
    print('bbox: fpos0={}, fpos1={}, fpos2={}'.format(fpos0, fpos1, fpos2))
    print('bbox: ip0={}, ip1={}, ip2={}, (kp1,kp2)=({},{})'.format(ip0, ip1, ip2, kp1, kp2))
    return (kp1, kp2)

if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='draw-spheres: Draw spheres on each site.  Ex) draw-spheres.py Gd2C.csfx Gd2C-spheres.vasp --ndim="10 10 10" > Gd2C.spheres.out')
    parser.add_argument('infile', help='input file')
    parser.add_argument('outfile', help='output file')
    parser.add_argument('--ndim', type=str, default='10 10 10',
                        help='Dimension of the 3D grid. default="10 10 10"')
    parser.add_argument('--accumulate', default=False, action="store_true",
                        help='Draw spheres accumulatively')
    parser.add_argument('--inverse', default=False, action="store_true",
                        help='Draw inverse of spheres (where no spheres)')
    parser.add_argument('--verbose', type=int, default=3,
                        help='Verbose level index')

    args = parser.parse_args()

    app = draw_spheres()
    app.run(args)
