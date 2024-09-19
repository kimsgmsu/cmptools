#!/usr/bin/env python

"""
chgcar-R2H: Convert a chgcar file from rhombohedral to hexagonal coord system
"""

from __future__ import print_function
import os
import numpy as np
import Molecule
import vasp
import Cmptlib

def Rhom2Hex(mol0):
    """
    Convert a molecule from rhombohedral to hexagonal lattice vectors
    """

    mol = mol0.copy()
    # new lattice vectors
    mol.unitcell.type = "Rhombohedral"
    mol.unitcell.avec[:,0] = -mol0.unitcell.avec[:,0] +mol0.unitcell.avec[:,1]
    mol.unitcell.avec[:,1] =  mol0.unitcell.avec[:,0] -mol0.unitcell.avec[:,2]
    mol.unitcell.avec[:,2] =  mol0.unitcell.avec[:,0] +mol0.unitcell.avec[:,1] +mol0.unitcell.avec[:,2] 
    mol.unitcell.construct()
    print('New rhombohedral unit cell:')
    mol.unitcell.report()

    # position generators
    pgen = np.zeros((3,3))
    pgen[:,0] = [0, 0, 0]
    pgen[:,1] = [0, 1, 0]
    pgen[:,2] = [1, 1, 0]

    # generate atom positions
    mol.unitcell.atomcoord = "fractional"
    mol.atomlist = []
    for atom0 in mol0.atomlist:
        for ig in range(3):
            gv = pgen[:,ig]
            atom = atom0.copy()
            fpos = np.zeros(3)
            for i in range(3):
                fpos[i] = atom0.fpos[i] + gv[i]
            atom.cpos = np.dot(mol0.unitcell.f2c, fpos)
            atom.fpos = np.dot(mol.unitcell.c2f, atom.cpos)
            mol.atomlist.append(atom)

    mol.construct()

    from scipy.interpolate import RegularGridInterpolator
    xlist0 = np.linspace(0.0, 1.0, mol0.volumedata.npnt[0])
    ylist0 = np.linspace(0.0, 1.0, mol0.volumedata.npnt[1])
    zlist0 = np.linspace(0.0, 1.0, mol0.volumedata.npnt[2])
    my_func = RegularGridInterpolator((xlist0, ylist0, zlist0), mol0.volumedata.data)

    nx = mol0.volumedata.npnt[0]
    ny = mol0.volumedata.npnt[1]
    nz = mol0.volumedata.npnt[2]
    
    xfactor=1.0/np.float_(nx-1)
    yfactor=1.0/np.float_(ny-1)
    zfactor=1.0/np.float_(nz-1)

    mol.volumedata = Molecule.VolumeData(nx,ny,nz)
    for iz in xrange(nz):
        z = iz*zfactor
        for iy in xrange(ny):
            y = iy*yfactor
            for ix in xrange(nx):
                x = ix*xfactor
                fpos = np.array([x,y,z])
                cpos = np.dot(mol.unitcell.f2c, fpos)
                fpos0 = np.dot(mol0.unitcell.c2f, cpos)
                fpos0 = normalize_fraction(fpos0)
                mol.volumedata.data[ix,iy,iz] = my_func(fpos0)
    
    mol.report()
    return mol

class chgcar_R2H:

    def __init__(self, reportlevel):

        self.reportlevel = reportlevel
        self.header = []

    def run(self, infile, outfile):

        # open input file
        try:
            self.infile = open(infile, 'r')
            print('Reading input file "{}"'.format(infile))
        except IOError:
            print('Failed to open the input file "{}"'.format(infile))
            raise

        # read vasp file
        mol0 = vasp.read(self.infile, self.reportlevel)

        # convert to hexagonal lattice
        mol1 = Cmptlib.Rhom2Hex(mol0)

        # open output file
        try:
            self.outfile = open(outfile, 'w')
            print('Opened file "{}" for writing chgden'.format(outfile))
        except IOError:
            print('Failed to open the output file "{}"'.format(outfile))
            raise

        # write output file
        vasp.write(mol1, self.outfile)

        # close files
        self.infile.close()
        self.outfile.close()

if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='chgcar-R2H: '+
        'Convert a chgcar file from rhombohedral to hexagonal coord system')
    parser.add_argument('--reportlevel', type=int, default=4,
                        help='Report level index (0-5)')
    parser.add_argument('infile', help='input chgcar file')
    parser.add_argument('outfile', help='output file')

    args = parser.parse_args()

    app = chgcar_R2H(reportlevel=args.reportlevel)
    app.run(args.infile, args.outfile)
