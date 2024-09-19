#!/usr/bin/env python
"""
Module for fdf (siesta) file I/O
"""

from __future__ import print_function
import sys
import os
import numpy as np
import Molecule

# characters to ignore
IGNORECHARS = '.-_'

def to_Ang(vstr):
    value, unit = vstr.split(' ', maxsplit=1)
    try:
        v = float(value)
    except:
        raise RuntimeError('Invalid value "{}" for length'.format(value))
    
    if unit == '':
        pass
    elif unit == 'ang':
        pass
    elif unit == 'bohr':
        v *= BOHR
    elif unit == 'm':
        v *= METER
    elif unit == 'cm':
        v *= CENTIMETER
    elif unit == 'nm':
        v *= NANOMETER
    else:
        raise RuntimeError('Invalid unit "{}" for length'.format(unit))
    return v

def readline_fdf(f):
    while True:
        try:
            s = f.readline()
        except:
            s = None
            return s
        # ignore comment
        s = s.split('#', maxsplit=1)[0]
        # ignore case
        s = s.strip().lower()
        # ignore some chars
        for c in IGNORECHARS:
            s = s.replace(c, '')
        # remove multiple whitespaces
        s = ' '.join(s.split())
        # done if not blank
        if len(s) > 0:
            return s

def read(f, reportlevel=0):

    mol = Molecule.Molecule()
    mol.unitcell = Molecule.UnitCell()

    # process until done
    while True:
        line = readline_fdf(f)
        if line is None:
            break
        # keyword and its value
        keyword, value = line.split(' ', maxsplit=1)

        if keyword == '%block':
            _read_block(value, f, mol)
        elif keyword == 'systemname':
            mol.info = value
        elif keyword == 'systemlabel':
            mol.name = value
        elif keyword == 'numberofatoms':
            mol.natom = int(value)
        elif keyword == 'latticeconstant':
            mol.unitcell.lattParam = to_Ang(value)
            
    mol.finalize()
    return mol

def _read_block(name, f, mol):
    if name == 'latticevectors':
        for aid in range(3):
            line = readline_fdf(f)
            slist = line.split()
            for i,s in enumerate(slist):
                mol.unitcell.avec[aid][i] = float(s)*mol.unitcell.lattParam
        keyword = None
    elif name == 'chemicalspecieslabel':
        mol.typelist = []
        for aid in range(3):
            line = readline_fdf(f)
            slist = line.split()
            for i,s in enumerate(slist):
                mol.unitcell.avec[aid][i] = float(s)*mol.unitcell.lattParam
        keyword = '%endblock'

    # skip to the end of the block
    while True:
        if keyword == '%endblock':
            break
        line = readline_fdf(f)
        keyword, value = line.split(' ', maxsplit=1)

def _read_unitcell(f):
    uc = Molecule.Unitcell()
    # lattice parameter
    line = f.readline().strip()
    uc.lattParam = float(line)
    # primitive axes
    for aid in range(3):
        line = f.readline().strip()
        slist = line.split()[0:3]
        for i,s in enumerate(slist):
            uc.avec[aid][i] = float(s)*uc.lattParam
    uc.finalize()
    return uc

def _read_typelist(f):
    typelist = []
    # test if this is the line for number of atoms
    line = f.readline().strip()
    words = line.split()
    try:
        natom = int(words[0])
        vasp5 = False
    except:
        vasp5 = True

    for word1 in words:
        if vasp5:
            name = word1
        else:
            name = 'XX' + str(len(typelist)+1)
        typelist.append(name)

    # number of atoms for each element
    if vasp5:
        line = f.readline().strip()
        words = line.split()

    nalist = []
    for sn in words:
        nalist.append(int(sn))
    print("nalist=", nalist)
    return typelist, nalist

def _read_atomlist(f, typelist, nalist):
    # coord system
    line = f.readline().strip()
    coord = line.split()[0].lower()
    if coord[0] == 'd':
        coord = 'frac'
    else:
        coord = 'cart'
    atomlist = []
    for index, type in enumerate(typelist):
        for iatom in range(nalist[index]):
            line = f.readline().strip()
            words = line.split()
            atom = Molecule.Atom()
            atom.id = len(atomlist) + 1
            atom.type = type
            s1 = words[0]
            s2 = words[1]
            s3 = words[2]
            if coord == 'frac':
                atom.fpos = np.array([float(s1),float(s2),float(s3)])
            else:
                atom.cpos = np.array([float(s1),float(s2),float(s3)])

            if len(words) > 3:
                atom.info = ' '.join(words[3:])
            else:
                atom.info = ''
            atomlist.append(atom)
    return atomlist, coord

def _read_volume_data(f):
    """
    Read volumetric data
    """
    vd = Molecule.VolumeData()
    # skip one blank line
    if f.readline() == '':
        return vd
    # grid dimensions
    line = f.readline()
    if line == '':
        return vd

    vd.npnt = []
    try:
        words = line.split()
        for i in range(3):
            vd.npnt.append(int(words[i]))
    except Exception:
        raise RuntimeError('Invalid or missing dimension for volume data')

    #read data on grid points
    vd.data = np.zeros( (vd.npnt[0], vd.npnt[1], vd.npnt[2]), dtype=np.float_ )
    words = []
    for i2 in xrange(vd.npnt[2]):
        for i1 in xrange(vd.npnt[1]):
            for i0 in xrange(vd.npnt[0]):
                if len(words) < 1:
                    line = f.readline()
                    words = line.split()
                vd.data[i0,i1,i2] = np.float_(words.pop(0))
    return vd

def write(mol, f, coord=None, reportlevel=0):
    # Molecule name
    f.write('{}\n'.format(mol.name))
    # unit cell
    f.write("%s\n" % mol.unitcell.lattParam)
    for aid in range(3):
        avec = mol.unitcell.avec[aid]/mol.unitcell.lattParam
        f.write("%s\n" % str(avec)[1:-1])

    # set of different atom types
    typelist = []
    for atom in mol.atomlist:
        if not (atom.type in typelist):
            typelist.append(atom.type)

    # sort atoms
    alist = []
    nalist = []
    for type in typelist:
        na = 0
        for atom in mol.atomlist:
            if atom.type == type:
                na += 1
                alist.append(atom)
        nalist.append(na)

    # type list
    s = ' '.join([str(t) for t in typelist])
    f.write(' {}\n'.format(s))
    s = ' '.join([str(na) for na in nalist])
    f.write(' {}\n'.format(s))
    
    # atomlist
    if coord == None:
        coord = mol.atomcoord
    if coord == "frac":
        s = "Direct"
    else:
        s = "Cartesian"
    f.write('{}\n'.format(s))
    for atom in alist:
        if coord == "frac":
            pos = atom.fpos
        else:
            pos = atom.cpos
        f.write(' {}'.format(str(pos)[1:-1]))
        f.write(' {}\n'.format(atom.info))

    # volume data
    vd = mol.volumedata
    if len(vd.npnt) < 1:
        pass
    else:
        f.write('\n')
        for npt in vd.npnt:
            f.write(' {}'.format(npt))
        f.write('\n')
        max_count_per_line = 5
        line = ''
        count = 0
        for i2 in xrange(vd.npnt[2]):
            for i1 in xrange(vd.npnt[1]):
                for i0 in xrange(vd.npnt[0]):
                    word = '{:17.10E}'.format(vd.data[i0, i1, i2])
                    line = line + ' ' + word
                    count = count + 1
                    if count >= max_count_per_line:
                        f.write(line+'\n')
                        line = ''
                        count = 0

        if len(line) > 0:
            f.write(line+'\n')
                        

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
