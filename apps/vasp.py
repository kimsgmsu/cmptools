#!/usr/bin/env python
"""
Module for vasp file I/O
"""

import sys
import os
import numpy as np
import Molecule, Atom, Lattice
import util

def next_block_is_volumetric_data(f):
    """
    Test if the next block is a volumetric data
    """
    pos = f.tell()
    nbline = util.get_next_nb_line(f)
    try:
        words = nbline.split()
        ndim = [int(words[i]) for i in range(3)]
        ans = True
    except:
        ans = False
    f.seek(pos)
    return ans
    
def read(f, verbose=0):

    try:
        # poscar portion
        mol = read_poscar(f, verbose=verbose)

        # volumetric data
        mol.volumeDataList = read_volumetric_data_list(f)
    except:
        raise Exception('%vasp-error: Invalid formatted vasp file')
    return mol

def read_poscar(f, verbose=0):
    """
    Read the "poscar" portion
    """

    mol = Molecule.Molecule()

    # header
    line = f.readline()
    mol.name = line.strip()
    mol.info = ''

    # lattice
    mol.latt = read_lattice(f)

    # element symbol list
    symbol_list, natom_of_symbol = read_symbols(f)

    # atom list
    mol.atomList, frac = read_atoms(f, mol, symbol_list, natom_of_symbol)

    mol.output_coord = 'frac' if frac else 'frac'

    # complete element table
    mol.make_ElementTable()

    # velocities
    natom = len(mol.atomList)
    velocityList, pcHeader, pcList = read_velocities(f, natom)

    mol.velocityList = velocityList
    mol.pcHeader = pcHeader
    mol.pcList = pcList

    return mol

def read_lattice(f):
    # scale
    line = f.readline().strip()
    scale = float(line.split()[0])

    # lattice
    pbc = True
    matrix = np.eye(3)
    for aid in range(3):
        line = f.readline().strip()
        slist = line.split()[0:3]
        matrix[aid] = [float(s) for s in slist]

    if scale < 0:
        # In vasp, a negative scale factor is treated as a volume.
        # Need to translate this to a proper lattice vector scaling.
        vol = abs(np.linalg.det(matrix))
        try:
            scale = -scale
            matrix *= (scale/vol)**(1/3)
        except:
            raise ValueError('%vasp: Invalid lattice vectors')
    else:
        matrix *= scale

    latt_scale = np.full((3), scale)
    return Lattice.Lattice(matrix, pbc=pbc)

def read_symbols(f):

    symbol_list = []

    # test if this is the line for number of atoms
    line = f.readline().strip()
    words = line.split()
    try:
        natom = int(words[0])
        vasp5 = False
    except:
        vasp5 = True

    for w1 in words:
        if vasp5:
            symbol = w1
        else:
            symbol = 'XX' + str(len(symbol_list)+1)
        symbol_list.append(symbol)

    # number of atoms for each atom type
    if vasp5:
        line = f.readline().strip()
        words = line.split()

    natom_of_symbol = []
    for w1 in words:
        natom_of_symbol.append(int(w1))

    return symbol_list, natom_of_symbol

def read_atoms(f, mol, symbol_list, natom_of_symbol):

    # selective dynamics
    line = f.readline()
    sdyn_used = line[0] in 'sS'
    if sdyn_used:
        line = f.readline()

    # Is coords fractional?
    wd = line.split()[0].lower()
    if wd[0] == 'd':
        frac = True
    else:
        frac = False

    # atom list
    atomList = []
    for symbol, natom in zip(symbol_list, natom_of_symbol):
        for i in range(natom):

            atom = Atom.Atom()

            atom.symbol = symbol

            line = f.readline().strip()
            words = line.split()
            assert len(words) >= 3, \
                'Invalid position vector for atom #{}'.format(len(atomList)+1)
            pos = np.array([float(s) for s in words[0:3]])
            if frac:
                atom.pos = mol.to_cartesian(pos)
            else:
                atom.pos = pos

            if sdyn_used:
                words = words[3:]
                assert len(words) >= 3, \
                    'Invalid sdyn vector for atom #{}'.format(len(atomList)+1)
                atom.sdyn = ' '.join(words)
            else:
                atom.sdyn = None

            words = words[3:]
            atom.info = ' '.join(words)

            atomList.append(atom)

    return atomList, frac

def read_velocities(f, natom):

    velocityList = None
    pcHeader = None
    pcList = None

    if util.end_of_file(f) or next_block_is_volumetric_data(f):
        return velocityList, pcHeader, pcList

    # velocity data
    util.skip_blank_lines(f)
    velocityList = []
    for iatom in range(natom):
        line = f.readline().strip()
        words = line.split()
        assert len(words) >= 3, 'Invalid velocity vector for site #{}'.format(iatom+1)
        vel = np.array([float(s) for s in words[0:3]])
        velocityList.append(vel)

    if util.end_of_file(f) or next_block_is_volumetric_data(f):
        return velocityList, pcHeader, pcList

    # Predictor-corrector header
    # First line in chunk is a key in CONTCAR
    # Second line is POTIM
    # Third line is the thermostat parameters
    try:
        pcHeader = [f.readline() for i in range(3)]
    except:
        raise SyntaxError('Incorrect predictor-corrector header')

    # Parse the predictor-corrector data
    # There are 3 sets of 3xN Predictor corrector parameters
    # Rest is three sets of parameters, each set contains
    # x, y, z predictor-corrector parameters for every atom in orde

    # x-components
    pc_x = []
    for iatom in range(natom):
        line = f.readline().strip()
        pc1 = [float(w) for w in line.split()]
        pc_x.append(pc1)

    # y-components
    pc_y = []
    for iatom in range(natom):
        line = f.readline().strip()
        pc1 = [float(w) for w in line.split()]
        pc_y.append(pc1)

    # z-components
    pc_z = []
    for iatom in range(natom):
        line = f.readline().strip()
        pc1 = [float(w) for w in line.split()]
        pc_z.append(pc1)

    pcList = []
    for iatom in range(natom):
        pcList.append([pc_x[iatom], pc_y[iatom], pc_z[iatom]])
    return velocityList, pcHeader, pcList

def read_volumetric_data_list(f):
    """
    Read volumetric data list
    """
    vdList = []

    # get next non-blank line
    dimline = util.get_next_nb_line(f)

    # process lines until EOF
    while True:

        # if EOF, we are done
        if len(dimline) < 1:
            return vdList

        kvd = len(vdList)
        # read one volumetric data
        vd, dimline = read_volumetric_data(f, kvd, dimline)

        vdList.append(vd)

    return vdList

def read_volumetric_data(f, kvd, dimline):

    # grid dimensions
    ndim = []
    try:
        words = dimline.split()
        for i in range(3):
            ndim.append(int(words[i]))
    except Exception:
        raise RuntimeError('Invalid or missing dimension for volumetric data')

    nx = ndim[0]
    ny = ndim[1]
    nz = ndim[2]
    print('*** Reading volume data set {} on ({},{},{}) grid ***'.format(kvd, nx, ny, nz))
    #read data on grid points
    data = np.zeros( (nx, ny, nz), dtype=np.float_ )
    words = []
    for i2 in range(nz):
        for i1 in range(ny):
            for i0 in range(nx):
                if len(words) < 1:
                    line = f.readline()
                    words = line.split()
                data[i0,i1,i2] = np.float_(words.pop(0))

    # extra data
    lines = []
    while True:
        line = util.get_next_nb_line(f)
        if len(line) < 1 or line == dimline:
            break
        lines.append(line)
    extra = '\n'.join(lines)
    name = 'set {}'.format(kvd)
    return Molecule.VolumetricData(name, data, extra), line

def write(mol, f, coord=None, verbose=3):

    # poscar portion
    write_poscar(mol, f, coord, verbose)

    # volumetric data list
    write_volumetric_data_list(f, mol.volumeDataList)

def write_poscar(mol, f, coord=None, verbose=3):

    if coord == None:
        coord = mol.output_coord

    # format for float numericals
    ffmt = '{:>12.6f}'

    # Header
    f.write('{}\n'.format(mol.name))

    # Lattice
    # NOTE: vesta has a bug in reading vasp file in cartesian coord.
    # As a workaround, use scale=1.0 for cartesian coord.
    latt = mol.latt

    f.write('{}\n'.format(1.0))
    for v in latt.matrix:
        f.write(' '.join([ffmt.format(c) for c in v]))
        f.write('\n')

    # Element list
    s = ' '.join([symbol for symbol in mol.distinct_atomic_symbols])
    f.write(s+'\n')
    s = ' '.join([str(mol.natom_of_symbol(symbol)) for symbol in mol.distinct_atomic_symbols])
    f.write(s+'\n')

    # selective dynamics
    sdyn_used = mol.sdyn_used
    if sdyn_used:
        f.write('Selective dynamics\n')
    
    # Atom position
    s = 'Direct' if coord == 'frac' else 'Cartesian'
    f.write(s+'\n')

    # atom list grouped by element symbols
    for symbol in mol.distinct_atomic_symbols:
        for atom in mol.atomList_of_symbol(symbol):
            if coord == 'cart':
                v = atom.pos
            else:
                v = mol.to_fractional(atom.pos)
            line = ' '.join([ffmt.format(c) for c in v])
            if sdyn_used:
                if atom.sdyn is not None:
                    line += ' ' + str(atom.sdyn)
                else:
                    line += ' ' + 'T T T'
            if atom.info is not None:
                line += ' ' + atom.info
            f.write(line +'\n')

def write_volumetric_data_list(f, vdList):
    if vdList is None:
        return

    f.write('\n')
    for vd in vdList:
        write_volumetric_data(f, vd)

def write_volumetric_data(f, vd):

    dim = vd.data.shape
    line = ' '.join([str(n) for n in dim])
    f.write(line +'\n')

    ffmt = ' {:17.10E}'
    ndata_per_line = 5
    ndata = 0
    for i2 in range(dim[2]):
        for i1 in range(dim[1]):
            for i0 in range(dim[0]):
                f.write(ffmt.format(vd.data[i0,i1,i2]))
                ndata += 1
                if ndata >= ndata_per_line:
                    f.write('\n')
                    ndata = 0
    if ndata > 0:
        f.write('\n')

    # extra data
    if vd.extra is not None:
        f.write(' {}\n'.format(vd.extra))

if __name__ == "__main__":
    # stand-alone mode
    if len(sys.argv) < 3:
        print('Usage: vasp.py infile outfile')
    else:
        infile = sys.argv[1]
        f = open(infile, 'r')
        print('Reading input VASP file "{}"'.format(infile))
        mol = read(f, reportlevel=5)
        f.close()
        outfile = sys.argv[2]
        print('Writing output VASP file "{}"'.format(outfile))
        f = open(outfile, "w")
        write(mol, f, reportlevel=5)
        f.close()
