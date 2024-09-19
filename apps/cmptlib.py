#!/usr/bin/env python

"""
cmptlib: library module for cmptools
"""

import os, math
import numpy as np
from pymatgen.core.surface import Slab, SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import util, Transform3d
import Lattice, Molecule, Atom
import msf, vasp, pdb, xsf, xyz

class Neighbor:
    """
    class for a neighboring atom.
    """
    def __init__(self, atom, aid, cell, dist, pos):
        self.atom = atom
        self.aid = aid
        self.cell = cell
        self.dist = dist
        self.pos = pos

def addAtom(symbol, mol, posList, overlap, dtol):

    # set molecule
    result = mol.copy()

    # add new atomic symbol
    if symbol not in result.distinct_atomic_symbols:
        result.add_Element(symbol)

    # add new atoms
    for pos in posList:
        minaid, mindist = result.closest_atom(pos, coord='cart')
        if mindist < dtol:
            if overlap == 'discard':
                continue
            elif overlap == 'replace':
                result.delete_atom(minaid)
            elif overlap == 'keep':
                pass
            else:
                util.Error('%addAtom: Invalid overlap option, "{}"'.format(overlap))
        atom = Atom.Atom()
        atom.symbol = symbol
        atom.pos = pos
        result.atomList.append(atom)

    return result

def attachMol(mol1, anchor1, mol2, anchor2, overlap, dtol):

    # set molecule
    result = mol1.copy()

    # copy elements
    for symbol in mol2.distinct_atomic_symbols:
        if symbol not in result.distinct_atomic_symbols:
            result.add_Element(symbol)

    # find the affine transformation
    trans = Transform3d.find_affine_transformation(anchor2, anchor1)

    print('testing transformation')
    for pos in anchor2.pointList:
        pos2 = trans(pos.reshape((1,3))).reshape(3)
        print('{} => {}'.format(pos, pos2))
    
    # attach atoms
    for ka2, atom2 in enumerate(mol2.atomList):
        print('%attachMol: attaching atom #{} of the mol2:'.format(ka2+1))
        pos2 = atom2.pos
        pos3 = trans(pos2.reshape((1,3))).reshape(3)
        minaid, mindist = result.closest_atom(pos3, coord='cart')

        print('{:>4}'.format(ka2+1), end='')
        print(' {:>4}'.format(atom2.symbol), end='')
        print(' (', end='')
        for c in pos2:
            print(' {:>8.4f}'.format(c), end='')
        print(') => (', end='')
        for c in pos3:
            print(' {:>8.4f}'.format(c), end='')
        print(') ', end='')
        print('min.dist.atom={} min.dist={:>8.4f}'.format(minaid, mindist), end='')
        atom1 = result.atom_of_index(minaid)
        print(' {:>4} ('.format(atom1.symbol), end='')
        pos1 = atom1.pos
        for c in pos1:
            print(' {:>8.4f}'.format(c), end='')
        print(')')

        if mindist < dtol:
            if overlap == 'discard':
                print('%attachMol: Discarding atom #{} of the mol2 molecule'.format(ka2+1))
                continue
            elif overlap == 'replace':
                print('%attachMol: Replacing atom #{} of the mol1 molecule'.format(minaid))
                result.delete_atom(minaid)
            elif overlap == 'keep':
                print('%attachMol: Keeping overlapping atom #{} of the mol2'.format(ka2+1))
                pass
            else:
                util.Error('%attachMol: Invalid overlap option, "{}"'.format(overlap))
        atom3 = atom2.copy()
        atom3.pos = pos3
        result.atomList.append(atom3)

    return result

def boxMol(mol, padding):
    """
    Put a bounding box around a molecule
    """
    import collections
    Range = collections.namedtuple('Range', 'min max')
    mol_range = []
    for idir in range(3):
        rng = Range(min=math.inf, max=-math.inf)
        mol_range.append(rng)
    for atom in mol.atomList:
        for idir in range(3):
            if atom.pos[idir] < mol_range[idir].min:
                mol_range[idir].min = atom.pos[idir]
            if atom.pos[idir] > mol_range[idir].max:
                mol_range[idir].max = atom.pos[idir]

    # set bounding box
    lengths = []
    angles = []
    for idir in range(3):
        leng = (mol_range[idir].max - mol_range[idir].min)
        leng += padding[idir]
        ang = 90
        lengths.append(leng)
        angles.append(ang)
    
    mol.lattice = Lattice.from_parameters(lengths, angles, pbc=False)

    # center of unit cell
    cobox = np.zeros(3)
    for idir in range(3):
        cobox += mol.lattice.matrix[idir]*0.5

    # Move the center of mass to the origin
    dvec = cobox - mol.center_of_mass()
    for atom in mol.atomList:
        atom.pos = atom.pos + dvec

def centerMol(mol):

    # move atoms into the unit cell
    # center of unit cell
    cobox = np.zeros(3)
    for idir in range(3):
        cobox += mol.lattice.matrix[idir]*0.5
    print('centerMol: cobox =', cobox)

    # displacement vector
    dvec = cobox - mol.center_of_mass
    for atom in mol.atomList:
        atom.pos = atom.pos + dvec

def combineVolumeData(mol1, mol2, expr='a+b'):

    print('combination expression: ', expr)
    def comb_func(a,b):
        return eval(expr)

    # Use the structure of mol1
    mol3 = mol1.copy()

    vdlist1 = mol1.volumeDataList
    vdlist2 = mol2.volumeDataList
    vdlist3 = mol3.volumeDataList
    for (vd1, vd2, vd3) in zip(vdlist1, vdlist2, vdlist3):
        ndim = vd1.data.shape
        for i2 in range(ndim[2]):
            for i1 in range(ndim[1]):
                for i0 in range(ndim[0]):
                    d1 = vd1.data[i0,i1,i2]
                    d2 = vd2.data[i0,i1,i2]
                    d3 = comb_func(d1, d2)
                    vd3.data[i0,i1,i2] = d3
    return mol3

def copyMol(mol, sel):

    result = mol.copy()

    # copy atoms 
    result.atomList = []
    for atom in mol.selected_atoms(sel):
        result.atomList.append(atom.copy())

    return result

def deleteAtoms(mol, sel):

    result = mol.copy()

    # copy atoms except the deleted ones
    result.atomList = []
    for atom in mol.unselected_atoms(sel):
        result.atomList.append(atom.copy())

    return result

def drawSpheres(mol, ndim, accumulate, inverse, fname):

    if inverse:
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
        radius = mol.get_Element(atom.symbol).radius
        center = atom.pos
        bbox = util.get_bounding_box(center, radius, mol, ndim)
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
                        if accumulate:
                            data[ix0,iy0,iz0] += 1
                        elif inverse:
                            data[ix0,iy0,iz0] = 0
                        else:
                            data[ix0,iy0,iz0] = 1

    vd = Molecule.VolumetricData(name='Overlap index', data=data, extra=None)
    mol.volumeDataList = [vd]
    mol.report('In draw-spheres')
    writeMol(mol, fname)

def generateSlab(mol, miller, min_slab_size, min_vacuum_size,
                 lll_reduce, center_slab, in_unit_planes,
                 primitive, max_normal_search, reorient_lattice):

    print('%generateSlab: miller_index={}'.format(miller))
    print('center_slab={}'.format(center_slab))

    # pymatgen.Structure from a molecule
    mgStr = mol.to_pymatgen_structure()

    slabgen = SlabGenerator(initial_structure=mgStr,
                            miller_index=miller,
                            min_slab_size=min_slab_size,
                            min_vacuum_size=min_vacuum_size,
                            lll_reduce=lll_reduce,
                            center_slab=center_slab,
                            in_unit_planes=in_unit_planes,
                            primitive=primitive,
                            max_normal_search=max_normal_search,
                            reorient_lattice=reorient_lattice)
    slab = slabgen.get_slab()
    return Molecule.from_pymatgen_structure(slab)

def interpolate_height(mol, sel, free, fixed, coord='cart', hdir=2, rcut=5.0, hcut=1.0):
    """
    Interpolate the height of "free" atoms in regerence to "fixed" atoms
    """
    crange = [1, 1, 1]
    dtol = 0.001
    for atom in mol.selected_atoms(sel):
        if atom.symbol != free:
            continue
        pos1 = atom.pos
        fpos1 = mol.to_fractional(pos1)
        if coord == 'frac':
            h1 = fpos1[hdir]
        else:
            h1 = pos1[hdir]

        hsum = 0.0
        rsum = 0.0
        for atom2 in mol.atomList:
            if atom2.symbol != fixed:
                continue
            for c0 in range(-crange[0], crange[0]+1):
                for c1 in range(-crange[1], crange[1]+1):
                    for c2 in range(-crange[2], crange[2]+1):
                        cell = np.array([c0, c1, c2])
                        fpos2 = mol.to_fractional(atom2.pos) + cell
                        pos2 = mol.to_cartesian(fpos2)
                        dr = np.linalg.norm(pos2 - pos1)
                        if coord == 'frac':
                            h2 = mol.to_fractional(pos2)[hdir]
                        else:
                            h2 = pos2[hdir]
                        dh = h2 - h1
                        if dr < dtol:
                            continue
                        elif dr > rcut:
                            continue
                        elif abs(dh) > hcut:
                            continue
                        hsum += h2/dr
                        rsum += 1/dr
        h = hsum/rsum
        if coord == 'frac':
            fpos1[hdir] = h
            pos1 = mol.to_cartesian(fpos1)
        else:
            pos1[hdir] = h
        atom.pos = pos1

def makeSupercell(mol0, supercell):

    print('Making a supercell')
    # set molecule
    mol1 = mol0.copy()
    latt0 = mol0.lattice
    mol0.lattice.report('orig lattice:')
    # new lattice vector
    avec = np.matmul(supercell, latt0.matrix)
    latt1 = Lattice.Lattice(avec)
    mol1.lattice = latt1
    mol1.lattice.report('new lattice:')

    # get cell ranges
    crange = util.get_cell_range_from_box(supercell)
    print('crange=', crange)

    # Populate atoms
    mol1.atomList = []
    for c0 in range(crange[0][0], crange[0][1]):
        for c1 in range(crange[1][0], crange[1][1]):
            for c2 in range(crange[2][0], crange[2][1]):
                dpos = mol0.to_cartesian(np.array([c0, c1, c2]))
                for atom0 in mol0.atomList:
                    atom1 = atom0.copy()
                    pos = atom0.pos + dpos
                    fpos = latt1.to_fractional(pos)
                    if util.within_lattice(fpos):
                        atom1.pos = pos
                        mol1.atomList.append(atom1)

    return mol1

def makeTransitionalMol(mol1, mol2, nmol, prefix, fmt):

    # check if the two molecules have the same lattice structure
    if not Lattice.equivalent(mol1.lattice, mol2.lattice):
        util.Error('%makeTransitionalMol: Inequivalent lattices')

    # check if the two molecules have the same no. of atoms
    if len(mol1.atomList) != len(mol2.atomList):
        util.Error('%makeTransitionalMol: Different no. of atoms')

    # transitional structures
    for imol in range(nmol):
        kmol = imol + 1

        mol = mol1.copy()
        # displace atom positions
        for (atom, atom1, atom2) in \
            zip(mol.atomList, mol1.atomList, mol2.atomList):

            dpos = atom2.pos - atom1.pos
            dfpos = mol1.to_fractional(dpos)
            rdfpos = util.reduced_vector(dfpos, mol1.lattice.pbc)
            dpos = np.array(mol1.to_cartesian(rdfpos))*(kmol/(nmol+1))
            atom.pos = atom1.pos + dpos

        fname = '{}{}.{}'.format(prefix, kmol, fmt)
        mol.name = fname
        print('Saving transitional structure #{} to "{}"'.format(kmol, fname))
        writeMol(mol, fname)

def moveAtoms(mol, sel=None, vec=np.zeros(3),
              coord='cart', overlap='keep', dtol=0.001):

    if coord == 'frac':
        vec = mol.to_cartesian(vec)

    atomList = mol.atomList
    mol.atomList = []

    # add atoms
    for ka, atom in enumerate(atomList):

        if sel is None or ka+1 in sel:
            atom.pos += vec

        minaid, mindist = mol.closest_atom(atom.pos, coord='cart')
        if minaid is not None and mindist < dtol:
            if overlap == 'discard':
                print('%moveMol: Discarding atom #{}'.format(ka+1))
                continue
            elif overlap == 'replace':
                print('%moveMol: Replacing atom #{}'.format(minaid))
                mol.delete_atom(minaid)
            elif overlap == 'keep':
                print('%moveMol: Keeping overlapping atom #{}'.format(ka+1))
                pass
            else:
                util.Error('%moveMol: Invalid overlap option, "{}"'.format(overlap))

        mol.atomList.append(atom)

def mergeMol(mol, host, sel=None, vec=np.zeros(3),
            overlap='discard', dtol=0.001, verbose=1):

    print('%mergeMol: vec = {}'.format(vec))
    # set molecule
    result = host.copy()

    # add atoms
    for ka, atom in enumerate(mol.atomList):
        if sel is None:
            pass
        elif ka+1 in sel:
            pass
        else:
            continue

        pos = atom.pos + vec

        minaid, mindist = result.closest_atom(pos, coord='cart')
        if mindist < dtol:
            if overlap == 'discard':
                print('%mergeMol: Discarding atom #{}'.format(ka+1))
                continue
            elif overlap == 'replace':
                print('%mergeMol: Replacing atom #{}'.format(minaid))
                result.delete_atom(minaid)
            elif overlap == 'keep':
                print('%mergeMol: Keeping overlapping atom #{}'.format(ka+1))
                pass
            else:
                util.Error('%mergeMol: Invalid overlap option, "{}"'.format(overlap))
        atom2 = atom.copy()
        atom2.pos = pos
        result.atomList.append(atom2)

    return result

def moveMol(mol, sel=None, vec=np.zeros(3),
            overlap='discard', dtol=0.001, verbose=1):

    # set molecule
    result = mol.copy()
    if sel is None:
        sel = set({})
        for ka in range(mol.natom):
            sel.update({ka+1})

    # move atoms
    result.atomList = []
    for ka, atom in enumerate(mol.atomList):
        if (ka+1) in sel:
            pos = atom.pos + vec
        else:
            pos = atom.pos

        minaid, mindist = result.closest_atom(pos, coord='cart')
        if mindist < dtol:
            if overlap == 'discard':
                print('%moveMol: Discarding atom #{}'.format(ka+1))
                continue
            elif overlap == 'replace':
                print('%moveMol: Replacing atom #{}'.format(minaid))
                result.delete_atom(minaid)
            elif overlap == 'keep':
                print('%moveMol: Keeping overlapping atom #{}'.format(ka+1))
                pass
            else:
                util.Error('%moveMol: Invalid overlap option, "{}"'.format(overlap))
        atom3 = atom.copy()
        atom3.pos = pos
        result.atomList.append(atom3)

    return result

def printDisplacement(mol1, mol2, cutoff=None):

    # check if two molecules have the same lattice structure
    if not Lattice.equivalent(mol1.lattice, mol2.lattice):
        util.Error('%printDisplacement: Inequivalent lattices')

    # check if the two molecules have the same no. of atoms
    if len(mol1.atomList) != len(mol2.atomList):
        util.Error('%printDisplacement: Different no. of atoms')

    lhbar = 24
    print('='*lhbar)
    print('   PRINT DISPLACEMENT   ')
    print('='*lhbar)

    print('Name1: ', mol1.name)
    print('Name2: ', mol2.name)

    latt = mol1.lattice
    print('*** Lattice:')
    ffmt = '{:>10.6f}'
    lhbar = 59
    print('-'*lhbar)
    print('indx (   avec[x]   avec[y]   avec[z])     length      angle')
    print('-'*lhbar)
    for idir in range(3):
        print('{:>3}: ('.format(idir+1), end='')
        for v in latt.matrix[idir]:
            print(ffmt.format(v), end='')
        print(') ', end='')
        print(ffmt.format(latt.abc[idir]), end=' ')
        print(ffmt.format(latt.angles[idir]))
    print('-'*lhbar)
    print('  Cell volume = {:>10.6f}'.format(latt.volume))

    print('*** Displacement List ({})'.format(len(mol1.atomList)), end='')
    print(' cutoff = {}'.format(cutoff))
    sdyn_used = mol1.sdyn_used
    lhbar = 84
    if sdyn_used:
        lhbar += 8
    print('-'*lhbar)
    print('  ID TYPE (          CART-COORDS         )', end='')
    print(' (           FRAC-COORD         )', end='')
    print(' ( SDYN)' if sdyn_used else '', end='')
    print(' DISTANCE')
    print('-'*lhbar)
    
    # print displacement
    for ka1, (atom1, atom2) in enumerate(zip(mol1.atomList, mol2.atomList)):

        dpos = atom2.pos - atom1.pos
        dfpos = mol1.to_fractional(dpos)
        rdfpos = util.reduced_vector(dfpos, mol1.lattice.pbc)
        dpos = mol1.to_cartesian(rdfpos)
        dist = np.linalg.norm(dpos)
        if cutoff is not None:
            if dist < cutoff:
                continue

        # atom1
        print('{:>4}'.format(ka1+1), end='')
        print(' {:>4}'.format(atom1.symbol), end='')
        print(' (', end='')
        for c in atom1.pos:
            print(' {:>9.5f}'.format(c), end='')
        print(') (', end='')
        fpos = mol1.to_fractional(atom1.pos)
        for c in fpos:
            print(' {:>9.5f}'.format(c), end='')
        if sdyn_used:
            if atom1.sdyn is None:
                s = 'T T T'
            else:
                s = str(atom1.sdyn)
            print(') ({})'.format(s), end='')
        else:
            print(')', end='')
        print()

        # atom2
        print('{:>4}'.format(''), end='')
        print(' {:>4}'.format(atom2.symbol), end='')
        print(' (', end='')
        for c in atom2.pos:
            print(' {:>9.5f}'.format(c), end='')
        print(') (', end='')
        fpos = mol2.to_fractional(atom2.pos)
        for c in fpos:
            print(' {:>9.5f}'.format(c), end='')
        if sdyn_used:
            if atom2.sdyn is None:
                s = 'T T T'
            else:
                s = str(atom2.sdyn)
            print(') ({})'.format(s), end='')
        else:
            print(')', end='')
        print()

        # diff
        print('{:>4}'.format(''), end='')
        print(' {:>4}'.format('--'), end='')
        print(' (', end='')
        for c in dpos:
            print(' {:>9.5f}'.format(c), end='')
        print(') (', end='')
        for c in rdfpos:
            print(' {:>9.5f}'.format(c), end='')
        if sdyn_used:
            print(')       ', end='')
        else:
            print(')', end='')
        print(' {:>9.5f}'.format(dist))

    print('='*lhbar)

def printNeighborMap(mol, cutoff, padding, onlyatoms, tol, only_closer_dist=None):

    latt = mol.lattice
    crange = padding

    lhbar = 90
    print('='*lhbar)
    print('   PRINT NBRMAP   ')
    
    for ka1, atom1 in enumerate(mol.atomList):

        #check if we need to do this atom
        if len(onlyatoms) < 1:
            pass
        elif ka1+1 not in onlyatoms:
            continue

        print('-'*lhbar)
        pos1 = atom1.pos
        nbrlist1 = []
        for ka2, atom2 in enumerate(mol.atomList):
            for c0 in range(-crange[0], crange[0]+1):
                for c1 in range(-crange[1], crange[1]+1):
                    for c2 in range(-crange[2], crange[2]+1):
                        cell = np.array([c0, c1, c2])
                        fpos2 = mol.to_fractional(atom2.pos) + cell
                        pos2 = mol.to_cartesian(fpos2)
                        dist = np.linalg.norm(pos2 - pos1)
                        if dist < tol:
                            pass
                        elif dist > cutoff:
                            pass
                        else:
                            nbr = Neighbor(atom2, ka2+1, cell, dist, pos2)
                            nbrlist1.append(nbr)
        nbrlist1.sort(key=lambda nbr: nbr.dist)
        # skip this atom if distance are larger than minimum
        if only_closer_dist is not None:
            if nbrlist1[0].dist > only_closer_dist:
                continue
        print('({} {})'.format(ka1+1, atom1.symbol), end=' ')
        print('({:>12.8f} {:>12.8f} {:>12.8f})'.format(pos1[0], pos1[1], pos1[2]))
        for inbr, nbr in enumerate(nbrlist1):
            print('    {:>3} {:>7.4f}'.format(inbr+1, nbr.dist), end=' ')
            print('({:>4} {:>2})'.format(nbr.aid, nbr.atom.symbol), end=' ')
            print('[{:>3}'.format(nbr.cell[0]), end=' ')
            print('{:>3}'.format(nbr.cell[1]), end=' ')
            print('{:>3}]'.format(nbr.cell[2]), end=' ')
            print('({:>12.8f}'.format(nbr.pos[0]), end=' ')
            print('{:>12.8f}'.format(nbr.pos[1]), end=' ')
            print('{:>12.8f})'.format(nbr.pos[2]))

    print('='*lhbar)

def randomizePos(mol, sel=None, dmax=1.0):

    import random

    for atom in mol.selected_atoms(sel):
        for kdir in range(3):
            atom.pos += random.uniform(-dmax, dmax)

def readMol(fname, fmt=None, option=None, verbose=0):
    """
    Read a molecule from an input file
    """
    # file format
    if fmt == None:
        (head,tail) = os.path.split(fname)
        (base,ext) = os.path.splitext(fname)
        fmt = ext[1:]
        if fmt == '':
            fmt = base.lower()

    # open input file
    try:
        f = open(fname, 'r')
    except IOError:
        util.Error('%readMol: Failed to open input file "{}"'.format(fname))

    if fmt == 'msf':
        mol = msf.read(f, verbose)
    elif fmt == "pdb":
        mol = pdb.read(f, verbose)
    elif fmt == "fdf":
        mol = fdf.read(f, verbose)
    elif fmt in ['vasp', 'poscar', 'chgcar']:
        mol = vasp.read(f, verbose)
    elif fmt == "xsf":
        mol = xsf.read(f, verbose)
    elif fmt == "xyz":
        mol = xyz.read(f, verbose)
    else:
        util.Error('Unrecognized format "{}"'.format(fmt))

    return mol

def reduceAtomPos(mol):
    mol.reduceAtomPos()

    for atom in mol.atomList:
        fpos = util.reduced_vector2(mol.to_fractional(atom.pos))
        atom.pos = mol.to_cartesian(fpos)

def remove_overlap_atoms(mol, cutoff=1.5, method='average', sel=None):

    # set molecule
    result = mol.copy()

    atomList = mol.atomList

    # remove overlapped atoms
    for ka1, atom1 in enumerate(mol.atomList):

        if atomList[ka1] is None:
            continue

        for ka2, atom2 in enumerate(mol.atomList):

            if ka2 <= ka1:
                continue

            if atomList[ka2] is None:
                continue

            dpos = atom2.pos - atom1.pos
            dfpos = mol.to_fractional(dpos)
            rdfpos = util.reduced_vector(dfpos, mol.lattice.pbc)
            dpos = mol.to_cartesian(rdfpos)
            dist = np.linalg.norm(dpos)

            if dist > cutoff:
                continue
            
            if method == 'average':
                atom1.pos = (atom1.pos + atom2.pos)/2
                atomList[ka2] = None
            elif method == 'select':
                if ka1+1 in sel:
                    atomList[ka2] = None
                else:
                    atomList[ka1] = None


    result.atomList = []
    for atom in atomList:
        if atom is not None:
            result.atomList.append(atom)

    return result

def replicateMol(mol0, nrep):

    print('Replicating molecule with nrep = {}'.format(nrep))
    # set molecule
    mol1 = mol0.copy()
    # lattice
    lattice0 = mol0.lattice.copy()
    matrix0 = lattice0.matrix
    matrix = matrix0.copy()
    for aid in range(3):
        matrix[aid] = matrix0[aid]*nrep[aid]
    mol1.lattice = Lattice.Lattice(matrix)

    # Generate atoms
    mol1.atomList = []
    for c0 in range(nrep[0]):
        for c1 in range(nrep[1]):
            for c2 in range(nrep[2]):
                dpos = mol0.to_cartesian(np.array([c0, c1, c2]))
                for atom0 in mol0.atomList:
                    atom1 = atom0.copy()
                    atom1.pos = atom0.pos + dpos
                    mol1.atomList.append(atom1)

    return mol1

def reproduceMol(mol0, latt):

    print('Reproducing molecule')
    # set molecule
    mol1 = mol0.copy()
    latt0 = mol0.lattice
    # new lattice vector
    avec = np.matmul(latt, latt0.matrix)
    latt1 = Lattice.Lattice(avec)
    mol1.lattice = latt1
    mol1.lattice.report('new lattice:')

    # get cell ranges
    origin=np.zeros(3)
    bbox=[[origin,avec[0]],[origin,avec[1]],[origin,avec[2]]]
    crange = util.get_cell_range(latt0, bbox)

    # Reproduce atoms
    mol1.atomList = []
    for c0 in range(crange[0][0], crange[0][1]):
        for c1 in range(crange[1][0], crange[1][1]):
            for c2 in range(crange[2][0], crange[2][1]):
                dpos = mol0.to_cartesian(np.array([c0, c1, c2]))
                for atom0 in mol0.atomList:
                    atom1 = atom0.copy()
                    pos = atom0.pos + dpos
                    fpos = latt1.to_fractional(pos)
                    if util.within_lattice(fpos):
                        atom1.pos = pos
                        mol1.atomList.append(atom1)

    return mol1

def standardizeLattice(mol0, method=None):

    mol1 = mol0.copy()
    # standard lattice
    a, b, c = mol0.lattice.abc
    angles = mol0.lattice.angles
    angles_r = np.radians(angles)
    cos_alpha, cos_beta, cos_gamma = np.cos(angles_r)
    sin_alpha, sin_beta, sin_gamma = np.sin(angles_r)

    if method == None or method == 'a-along-x':
        c1 = c * cos_beta
        c2 = (c * (cos_alpha - (cos_beta * cos_gamma))) / sin_gamma

        vector_a = [float(a), 0.0, 0.0]
        vector_b = [b * cos_gamma, b * sin_gamma, 0]
        vector_c = [c1, c2, math.sqrt(c ** 2 - c1 ** 2 - c2 ** 2)]

    elif method == 'c-along-z':
        val = (cos_alpha * cos_beta - cos_gamma) / (sin_alpha * sin_beta)
        # Sometimes rounding errors result in values slightly > 1.
        val = abs_cap(val)
        gamma_star = np.arccos(val)

        vector_a = [a * sin_beta, 0.0, a * cos_beta]
        vector_b = [
            -b * sin_alpha * np.cos(gamma_star),
            b * sin_alpha * np.sin(gamma_star),
            b * cos_alpha,
        ]
        vector_c = [0.0, 0.0, float(c)]

    else:
        raise Exception('%standardizeLattice: Invalid method "{}"'.format(method))

    # New lattice
    matrix = [vector_a, vector_b, vector_c]
    mol1.lattice = Lattice.Lattice(matrix, pbc=mol0.pbc)
    for (atom1, atom0) in zip(mol1.atomList, mol0.atomList):
        fpos = mol0.to_fractional(atom0.pos)
        atom1.pos = mol1.to_cartesian(fpos)

    return mol1

def symmetry_map(poslist, numlist, rot, trans, pbc, prec=1.0e-4):

    natom = len(poslist)
    nop = len(rot)
    sym_map = np.full((natom, nop), -1, dtype=int)
    for i0 in range(natom):
        p0 = poslist[i0]
        n0 = numlist[i0]
        for iop in range(nop):
            p1 = rot[iop].dot(p0) + trans[iop]
            for i2 in range(natom):
                n2 = numlist[i2]
                if n2 == n0:
                    p2 = poslist[i2]
                    dp = p1 - p2
                    dp = util.reduced_vector(dp, pbc)
                    dist = np.linalg.norm(dp)
                    if dist < prec:
                        sym_map[i0,iop] = i2+1
                        break

    return sym_map

def reportSymmetry(mol, prec=1.0e-5):

    # banner
    print('=====================')
    print('   REPORT SYMMETRY   ')
    print('=====================')

    # pymatgen.Structure from a molecule
    mgStr = mol.to_pymatgen_structure()

    analyzer = SpacegroupAnalyzer(mgStr,
                                  symprec=prec,
                                  angle_tolerance=5)

    print('==========================')
    print('   SPACE GROUP ANALYSIS   ')
    print('==========================')
    dataset = analyzer.get_symmetry_dataset()
    print('\nSpace Group: #{} {}'.format(
        dataset['number'],
        dataset['international']
    ))
    print('Precision =', prec)
    print('Transformation matrix to the standardized basis vectors:')
    for vec in dataset['transformation_matrix']:
        print('  {:>10.6f} {:>10.6f} {:>10.6f} '
              .format(vec[0], vec[1], vec[2]))
    shift = dataset['origin_shift']
    print('Origin shift: {:>10.6f} {:>10.6f} {:>10.6f} '
          .format(shift[0], shift[1], shift[2]))
    equivat = dataset['equivalent_atoms']
    # symmetry operations
    rot = dataset['rotations']
    trans = dataset['translations']
    poslist = []
    for atom in mol.atomList:
        poslist.append(mol.to_fractional(atom.pos))
    symblist = mol.atomic_symbols
    smap = symmetry_map(poslist, symblist, rot, trans, mol.pbc, prec=prec)
    print('Symmetry operations: ', len(rot))
    print('----------------------------------------------------------')
    print('  OP    ROT[1]     ROT[2]     ROT[3]          TRANS       ')
    print('----- ---------- ---------- ---------- -------------------')
    for iop, (r, t) in enumerate(zip(rot, trans)):
        print('{:>4}:'.format(iop+1), end=' ')
        print('({:>2} {:>2} {:>2})'
              .format(r[0][0], r[0][1], r[0][2]), end=' ')
        print('({:>2} {:>2} {:>2})'
              .format(r[1][0], r[1][1], r[1][2]), end=' ')
        print('({:>2} {:>2} {:>2})'
              .format(r[2][0], r[2][1], r[2][2]), end=' ')
        print('({:>7.4f} {:>7.4f} {:>7.4f})'
              .format(t[0], t[1], t[2]))

    # set multiplicity
    multipl = np.zeros(len(equivat), dtype=np.int)
    for ia in range(len(equivat)):
        multipl[equivat[ia]] += 1
    for ia in range(len(equivat)):
        multipl[ia] = multipl[equivat[ia]]
    wyckoff = dataset['wyckoffs']
    print('---------------------------------------------------', end='')
    print('----'*len(rot))
    print('INDX SYMB    x        y        z     EQAT MULT WYKO', end='')
    for iop in range(len(rot)):
        print(' {:>3}'.format(iop+1), end='')
    print()
    print('---- ---- -------- -------- -------- ---- ---- ----', end='')
    print('----'*len(rot))
    for ia, atom in enumerate(mol.atomList):
        fpos = mol.to_fractional(atom.pos)
        print('{:>4}'.format(ia+1), end=' ')
        print('{:>4}'.format(atom.symbol), end=' ')
        print('{:>8.4f}'.format(fpos[0]), end=' ')
        print('{:>8.4f}'.format(fpos[1]), end=' ')
        print('{:>8.4f}'.format(fpos[2]), end=' ')
        print('{:>4}'.format(equivat[ia]+1), end=' ')
        print('{:>4}'.format(multipl[ia]), end=' ')
        print('{:>4}'.format(wyckoff[ia]), end='')
        for iop in range(len(rot)):
            print(' {:>3}'.format(smap[ia,iop]), end='')
        print()
    print('---------------------------------------------------', end='')
    print('----'*len(rot))

    # primitive cell
    pc = analyzer.find_primitive()
    pcmol = Molecule.from_pymatgen_structure(pc)
    pcmol.report('\n=== Primitive cell ===')

    return pcmol

def wrapMol(mol, pads=None):
    """
    Wrap a molecule in a box
    """
    if pads is None:
        pads = np.full((3,2), -1)
    else:
        pads = np.array(pads)

    # original lattice
    latt0 = mol.lattice
    abc0 = latt0.abc

    # fractional position range
    frange = mol.frac_range()

    # flag for wrapping
    ifwrap = [False for i in range(3)]
    for idir in range(3):
        pad1 = pads[idir]
        if pad1[0] < 0 or pad1[1] < 0:
            ifwrap[idir] = False
        else:
            ifwrap[idir] = True

    # new lengths
    abc = np.zeros(3)
    for idir in range(3):
        pad1 = pads[idir]
        if ifwrap[idir]:
            width = (frange[idir,1]-frange[idir,0])*abc0[idir]
            abc[idir] = width+pads[idir,1] - pads[idir,0]
        else:
            abc[idir] = abc0[idir]
            frange[idir] = [0, 1]

    # new lattice vectors
    matrix = latt0.matrix
    for idir in range(3):
        if ifwrap[idir]:
            matrix[idir] *= abc[idir]/abc0[idir]

    # new lattice
    latt = Lattice.Lattice(matrix, pbc=latt0.pbc)

    # lower-left corner of padding
    pad00 = np.zeros(3)
    for idir in range(3):
        if ifwrap[idir]:
            pad00[idir] = pads[idir,0]/abc0[idir]

    # move atoms
    cpad00 = latt0.to_cartesian(pad00)
    crange0 = latt0.to_cartesian(frange[:,0])
    dpos = cpad00 - crange0
    for atom in mol.atomList:
        atom.pos += dpos

    mol.lattice = latt

def updateMol(mol, fname, fmt=None, coord=None, option=None, verbose=0):
    """
    Update a molecule to an output file
    """
    writeMol(mol, fname, fmt=fmt, coord=coord, option=option, verbose=verbose)
    return readMol(fname, fmt=fmt, option=option, verbose=verbose)

def writeMol(mol, fname, fmt=None, coord=None, option=None, verbose=0):
    """
    Write a molecule to an output file
    """
    if verbose > 3:
        print('Output file =', fname)

    # file format
    if fmt == None:
        (head,tail) = os.path.split(fname)
        (base,ext) = os.path.splitext(fname)
        fmt = ext[1:]
        if fmt == '':
            fmt = base.lower()

    if verbose > 3:
        print('Output format =', fmt)

    # open output file
    try:
        f = open(fname, 'w')
    except IOError:
        Error('Failed to open output file "{}"'.format(fname))

    if fmt == 'msf':
        msf.write(mol, f, coord, verbose)
    elif fmt == "fdf":
        fdf.write(mol, f, coord, verbose)
    elif fmt == "pdb":
        pdb.write(mol, f, verbose)
    elif fmt in ['vasp', 'poscar', 'chgcar']:
        vasp.write(mol, f, coord, verbose)
    elif fmt == "xsf":
        xsf.write(mol, f, verbose)
    elif fmt == "xyz":
        xyz.write(mol, f, verbose)
    else:
        util.Error('%writeMol: Unrecognized format "{}"'.format(fmt))

