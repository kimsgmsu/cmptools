#!/usr/bin/env python
"""
Module for XYZ file I/O
"""

import sys
import os
import numpy as np
import Molecule, Atom, Lattice
import util

def read(f, verbose=0):

    mol = Molecule.Molecule()

    # lattice
    scale = np.ones(3)
    matrix = np.eye(3)
    pbc = [False, False, False]
    mol.lattice = Lattice.Lattice(matrix, pbc=pbc, output_scale=scale)
    
    # no of atoms
    try:
        words = f.readline().strip().split()
        natom = int(words[0])
    except Exception:
        raise RuntimeError('Invalid or missing no of atoms')

    # header
    mol.name = f.readline().strip()
    mol.info = ''

    # atom list
    mol.atomList = []
    for iatom in range(natom):
        line = f.readline().strip()
        words = line.split()
        atom = Atom.Atom()
        atom.symbol = words[0]
        s1 = words[1]
        s2 = words[2]
        s3 = words[3]
        atom.pos = np.array([float(s1),float(s2),float(s3)])
        if len(words) > 4:
            atom.info = ' '.join(words[4:])
        else:
            atom.info = ''
        mol.atomList.append(atom)

    # ensure the lattice is larger than the molecule
    maxpos = np.zeros(3)
    for atom in mol.atomList:
        maxpos = np.maximum(maxpos, atom.pos)

    sfac = 2.0
    minleng = 10.0
    lengths = np.zeros(3)
    for idir in range(3):
        lengths[idir] = (maxpos[idir]+minleng)*sfac

    angles = np.full((3), 90.0)
    pbc = [False, False, False]
    mol.lattice = Lattice.Lattice.from_parameters(lengths, angles, pbc=pbc)
    return mol

def write(mol, f, verbose=0):
    # no of atoms
    f.write(' {}'.format(len(mol.atomList)))
    # Molecule name
    f.write('{}\n'.format(mol.name))
    # atomlist
    for atom in mol.atomList:
        f.write(' {}'.format(str(atom.pos)[1:-1]))
        f.write(' {}\n'.format(atom.info))

if __name__ == "__main__":
    # stand-alone mode
    if len(sys.argv) < 3:
        print('Usage: vasp.py infile outfile')
    else:
        infile = sys.argv[1]
        f = open(infile, 'r')
        print('Reading input XYZ file "{}"'.format(infile))
        mol = read(f, verbose=5)
        f.close()
        outfile = sys.argv[2]
        print('Writing output XYZ file "{}"'.format(outfile))
        f = open(outfile, "w")
        write(mol, f, verbose=5)
        f.close()
