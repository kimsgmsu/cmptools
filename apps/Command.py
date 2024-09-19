#!/usr/bin/env python
"""
Command processing module
"""

import sys
import numpy as np
import util, cmptlib
import Plane, Selection, Anchor, Point, Vector
import Lattice, Molecule, Kpoint
import msf

def anchor(system, comm):
    Anchor.process(system, comm)

def addAtom(system, comm):
    """
    Add atoms to a molecule
    """
    symbol = util.get_string("symbol", comm)
    mol = util.get_molecule("mol", comm, system)
    overlap = util.get_string("overlap", comm)
    dtol = util.get_float("dtol", comm, req=False, default=0.01)
    method = util.get_string('method', comm)
    coord = util.get_string("coord", comm, req=False, default='cart')
    if method == 'point':
        pidList = util.get_name_list('point', comm)
        posList = []
        for pid in pidList:
            pos = system.get_point(pid).pos
            posList.append(pos)
    else:
        raise Exception('%addAtoms: Invalid method "{}"'.format(method))

    result = cmptlib.addAtom(symbol, mol, posList, overlap=overlap, dtol=dtol)
    mid = util.get_string("result", comm)
    system.moleculeList.update({mid : result})

def attachMol(system, comm):
    mol1 = util.get_molecule("mol1", comm, system)
    anchor1 = util.get_anchor("anchor1", comm, system)
    mol2 = util.get_molecule("mol2", comm, system)
    anchor2 = util.get_anchor("anchor2", comm, system)
    overlap = util.get_string("overlap", comm)
    dtol = util.get_float("dtol", comm, req=False, default=0.01)
    result = cmptlib.attachMol(mol1, anchor1, mol2, anchor2, overlap, dtol)
    mid = util.get_string("result", comm)
    system.moleculeList.update({mid : result})

def boxMol(system, comm):
    mol = util.get_molecule("mol", comm, system)
    padding = util.get_value("padding", comm, system)
    mol.make_bounding_box(padding)

def centerMol(system, comm):
    mol = util.get_molecule("mol", comm, system)
    cmptlib.centerMol(mol)

def copyMol(system, comm):
    mol = util.get_molecule("mol", comm, system)
    sel = util.get_selection("sel", comm, system, req=False, default=None)
    result = cmptlib.copyMol(mol, sel)
    mid = util.get_string("result", comm)
    system.moleculeList.update({mid:result})

def deleteAtoms(system, comm):
    mol = util.get_molecule("mol", comm, system)
    sel = util.get_selection("sel", comm, system)
    result = cmptlib.deleteAtoms(mol, sel)
    mid = util.get_string("result", comm)
    system.moleculeList.update({mid:result})

def drawSpheres(system, comm):
    mol = util.get_molecule("mol", comm, system)
    ndim = util.get_value_list("ndim", comm, req=False, default="[1,1,1]")
    accumulate = util.get_bool("accumulate", comm, req=False, default=None)
    inverse = util.get_bool("inverse", comm, req=False, default=None)
    fname = util.get_string("file", comm)
    cmptlib.drawSpheres(mol, ndim, accumulate, inverse, fname)

def editAtoms(system, comm):
    mol = util.get_molecule("mol", comm, system)
    sel = util.get_selection("sel", comm, system)
    item = util.get_string("item", comm)
    if item == 'symbol':
        symbol = util.get_string("symbol", comm)
        for atom in mol.selected_atoms(sel):
            atom.symbol = symbol
    elif item == 'sdyn':
        sdyn = util.get_string("sdyn", comm)
        for atom in mol.selected_atoms(sel):
            atom.sdyn = sdyn
    else:
        util.Error('%makemol: Invalid item "{}"'.format(item))


def generateSlab(system, comm):
    mol0 = util.get_molecule('mol', comm, system)
    miller = util.get_value_list('miller_index', comm)
    min_slab_size = util.get_float('min_slab_size', comm)
    min_vacuum_size = util.get_float('min_vacuum_size', comm)
    lll_reduce = util.get_bool('lll_reduce', comm, req=False, default='False')
    center_slab = util.get_bool('center_slab', comm, req=False, default='False')
    in_unit_planes = util.get_bool('in_unit_planes', comm, req=False, default='False')
    primitive = util.get_bool('primitive', comm, req=False, default='False')
    max_normal_search = util.get_int('max_normal_search', comm,
                                     req=False, default='0')
    reorient_lattice = util.get_bool('reorient_lattice', comm,
                                     req=False, default='False')
    mol1 = cmptlib.generateSlab(mol0, miller, min_slab_size, min_vacuum_size,
                                lll_reduce, center_slab, in_unit_planes,
                                primitive, max_normal_search, reorient_lattice)
    mol1.elementTable = mol0.elementTable
    if system.verbose >= 4:
        mol1.report("%generateSlab: molecule")
    mid = util.get_string("result", comm)
    system.moleculeList.update({mid : mol1})

def interpolate_height(system, comm):
    mol = util.get_molecule("mol", comm, system)
    sel = util.get_selection("sel", comm, system)
    free = util.get_string("free", comm)
    fixed = util.get_string("fixed", comm)
    coord = util.get_string("coord", comm)
    hdir = util.get_int("hdir", comm)
    rcut = util.get_float('rcut', comm)
    hcut = util.get_float('hcut', comm)
    cmptlib.interpolate_height(mol, sel, free, fixed, coord, hdir, rcut, hcut)

def kpoint(system, comm):
    Kpoint.process(system, comm)

def lattice(system, comm):
    Lattice.process(system, comm)

def makeSupercell(system, comm):
    mol = util.get_molecule('mol', comm, system)
    supercell = util.get_value('supercell', comm, system)
    result = cmptlib.makeSupercell(mol, supercell)
    mid = util.get_string('result', comm)
    system.moleculeList.update({mid : result})

def makeTransitionalMol(system, comm):
    mol1 = util.get_molecule("mol1", comm, system)
    mol2 = util.get_molecule("mol2", comm, system)
    nmol = util.get_int("nmol", comm)
    prefix = util.get_string("prefix", comm)
    fmt = util.get_string("format", comm)
    cmptlib.makeTransitionalMol(mol1, mol2, nmol, prefix, fmt)

def mergeMol(system, comm):
    """
    merge a molecule to a "host" molecule
    """
    mol = util.get_molecule("mol", comm, system)
    host = util.get_molecule("host", comm, system)
    sel = util.get_selection("sel", comm, system, req=False, default=None)
    overlap = util.get_string("overlap", comm)
    dtol = util.get_float("dtol", comm, req=False, default=0.01)
    method = util.get_string('method', comm)
    coord = util.get_string("coord", comm, req=False, default='cart')
    if method == 'none':
        vec = np.zeros(3)
    elif method == 'vector':
        vec = util.get_vector("vec", comm, system).value
        if coord == 'frac':
            vec = host.to_cartesian(vec)
    elif method == 'atom-to-point':
        atom = util.get_atom('atom', comm, mol)
        point = util.get_point("point", comm, system)
        vec = point.pos - atom.pos
    elif method == 'atom-to-atom':
        atom1 = util.get_atom('atom1', comm, mol)
        atom2 = util.get_atom('atom2', comm, host)
        vec = atom2.pos - atom1.pos
    else:
        raise Exception('%mergeMol: Invalid method "{}"'.format(method))

    result = cmptlib.mergeMol(mol, host, sel, vec=vec,
                              overlap=overlap, dtol=dtol, verbose=system.verbose)
    mid = util.get_string("result", comm)
    system.moleculeList.update({mid : result})

def molecule(system, comm):
    Molecule.process(system, comm)

def move_atoms_to_unitcell(system, comm):
    mol = util.get_molecule("mol", comm, system)
    overlap = util.get_string("overlap", comm, req=False, default='keep')
    dtol = util.get_value("dtol", comm, system, req=False, default=0.01)
    mol.move_atoms_to_unitcell(overlap=overlap, dtol=dtol)

def moveAtoms(system, comm):
    mol = util.get_molecule("mol", comm, system)
    sel = util.get_selection("sel", comm, system, req=False, default=None)
    vec = util.get_vector("vec", comm, system)
    coord = util.get_string("coord", comm, req=False, default='cart')
    overlap = util.get_string("overlap", comm, req=False, default='keep')
    dtol = util.get_value("dtol", comm, system, req=False, default='0.01')
    cmptlib.moveAtoms(mol, sel, vec.value,
                      coord=coord, overlap=overlap, dtol=dtol)

def moveMol(system, comm):
    host = util.get_molecule("host", comm, system)
    mol = util.get_molecule("mol", comm, system)
    sel = util.get_selection("sel", comm, system, req=False)
    method = util.get_string('method', comm)
    coord = util.get_string("coord", comm, req=False, default='cart')
    overlap = util.get_string("overlap", comm, req=False, default='ignore')
    dtol = util.get_value("dtol", comm, system, req=False, default=0.01)
    if method == 'vector':
        vec = util.get_vector("vector", comm, system).value
        if coord == 'frac':
            vec = mol.to_cartesian(value)
    elif method == 'atom-to-pos':
        atom = util.get_atom('atom', comm, mol)
        pos = np.array(util.get_value_list("pos", comm))
        if coord == 'frac':
            pos = mol.to_cartesian(pos)
        vec = pos - atom.pos
    elif method == 'atom-to-atom':
        atom1 = util.get_atom('atom1', comm, mol)
        atom2 = util.get_atom('atom2', comm, mol)
        vec = atom2.pos - atom1.pos
    else:
        raise Exception('%moveMol: undefined method "{}"'.format(method))

    result = cmptlib.mergeMol(host, mol, sel, vec=vec,
                             overlap=overlap, dtol=dtol, verbose=system.verbose)
    mid = util.get_string("result", comm)
    system.moleculeList.update({mid : result})

def newMol(system, comm):
    result = msf.read_molecule(comm)
    mid = util.get_string("result", comm)
    system.moleculeList.update({mid : result})

def plane(system, comm):
    Plane.process(system, comm)

def point(system, comm):
    Point.process(system, comm)

def printDisplacement(system, comm):
    mol1 = util.get_molecule("mol1", comm, system)
    mol2 = util.get_molecule("mol2", comm, system)
    cutoff = util.get_float("cutoff", comm, req=False, default=None)
    cmptlib.printDisplacement(mol1, mol2, cutoff)

def printNeighborMap(system, comm):
    mol = util.get_molecule("mol", comm, system)
    cutoff = util.get_float("cutoff", comm, req=False, default="5.0")
    padding = util.get_value_list("padding", comm, req=False, default="[1,1,1]")
    onlyatoms = util.get_value_list("onlyatoms", comm, req=False, default=None)
    tol = util.get_float("dtol", comm, req=False, default="1.0e-4")
    only_closer_dist = util.get_float("only_closer_dist", comm, req=False, default=None)
    cmptlib.printNeighborMap(mol, cutoff, padding, onlyatoms, tol)

def randomizePos(system, comm):
    mol = util.get_molecule("mol", comm, system)
    sel = util.get_selection("sel", comm, system, req=False, default=None)
    dmax = util.get_value("dmax", comm, system, req=False, default='1.0')
    cmptlib.randomizePos(mol, sel, dmax)

def readMol(system, comm):
    fname = util.get_string("file", comm)
    fmt = util.get_string("format", comm, req=False)
    option = util.get_string("option", comm, req=False)
    mol = cmptlib.readMol(fname, fmt=fmt, option=option, verbose=system.verbose)
    mid = util.get_string("result", comm)
    system.moleculeList.update({mid : mol})
    if system.verbose >= 1:
        print('%readMol: "{}" {}'.format(fname, mol))

def remove_overlap_atoms(system, comm):
    mol = util.get_molecule("mol", comm, system)
    cutoff = util.get_value("cutoff", comm, system, req=False, default='1.5')
    method = util.get_string("method", comm, req=False, default='average')
    sel = util.get_selection("sel", comm, system, req=False, default=None)
    result = cmptlib.remove_overlap_atoms(mol, cutoff=cutoff, method=method, sel=sel)
    mid = util.get_string("result", comm)
    system.moleculeList.update({mid : result})

def replicateMol(system, comm):
    mol = util.get_molecule("mol", comm, system)
    nrep = util.get_value("nrep", comm, system)
    result = cmptlib.replicateMol(mol, nrep)
    print('%makemol: Replicated molecule: ', result)
    mid = util.get_string("result", comm)
    system.moleculeList.update({mid : result})

def reproduceMol(system, comm):
    mol = util.get_molecule("mol", comm, system)
    latt = util.get_value("latt", comm, system)
    result = cmptlib.reproduceMol(mol, latt)
    mid = util.get_string("result", comm)
    system.moleculeList.update({mid : result})

def report(system, comm):
    item = util.get_string("item", comm)
    header = util.get_string('header', comm, req=False)
    if item == 'selection':
        sel = util.get_selection("sel", comm, system)
        print('Selection {}: {}'.format(header, sel))
    elif item == 'anchor':
        anchor = util.get_anchor("anchor", comm, system)
        print('Anchor {}: {}'.format(header, anchor))
    elif item == 'plane':
        plane = util.get_plane("plane", comm, system)
        print('Plane {}: {}'.format(header, plane))
    elif item == 'molecule':
        mol = util.get_molecule("mol", comm, system)
        mol.report(header)
    elif item == 'variable':
        (var, value) = util.get_variable("var", comm, system)
        print('Variable {}={}'.format(var, value))
    else:
        util.Error('%makemol: Unrecognized item "{}"'.format(item))

def reportSymmetry(system, comm):
    mol = util.get_molecule("mol", comm, system)
    prec = util.get_float("prec", comm, req=False, default="1.0e-4")
    pcmol = cmptlib.reportSymmetry(mol, prec=prec)
    mid = util.get_string("pcmol", comm, req=False, default=None)
    if mid is not None:
        system.moleculeList.update({mid : pcmol})

def selection(system, comm):
    Selection.process(system, comm)

def setMol(system, comm):
    mol = util.get_molecule("mol", comm, system)
    item = util.get_string("item", comm)
    if item == 'name':
        value = util.get_string('value', comm)
        mol.name = value
    elif item == 'info':
        value = util.get_string('value', comm)
        mol.info = value
    elif item == 'lattice':
        rescale = util.get_bool('rescale', comm)
        print('setMol: rescale = ', rescale)
        latt0 = mol.latt
        latt1 = msf.read_lattice(comm, VARS=system.variableList)
        if rescale:
            print('setMol: rescaling atom positions')
            for atom in mol.atomList:
                fpos = latt0.to_fractional(atom.pos)
                atom.pos = latt1.to_cartesian(fpos)
        mol.latt = latt1
    else:
        raise Exception('%setMol: invalid item to set "{}"'.format(item))

def set_param(system, comm):
    for name, value in comm.items():
        print('%makemol.set_param: {0} = {1}'.format(name, value))
        if name == 'title':
            system.title = value
        elif name == 'verbose':
            system.verbose = int(value)
        else:
            util.Error('%set: Undefined parameter "{}"'.format(name))

def setVariable(system, comm):
    value = util.get_value("value", comm, system)
    vid = util.get_string("var", comm)
    system.variableList.update({vid : value})

def standardizeLattice(system, comm):
    mol = util.get_molecule("mol", comm, system)
    mid0 = util.get_string("mol", comm)
    method = util.get_string('method', comm)
    mid = util.get_string("result", comm, req=False, default=mid0)
    mol2 = cmptlib.standardizeLattice(mol, method=method)
    mid = util.get_string("result", comm, req=False, default=mid0)
    system.moleculeList.update({mid : mol2})

def updateMol(system, comm):
    mol = util.get_molecule("mol", comm, system)
    mid0 = util.get_string("mol", comm)
    fname = util.get_string("file", comm)
    fmt = util.get_string("format", comm, req=False)
    coord = util.get_string("coord", comm, req=False)
    option = util.get_string("option", comm, req=False)
    verbose = util.get_int("verbose", comm, req=False, default=system.verbose)
    mol2 = cmptlib.updateMol(mol, fname, fmt=fmt, coord=coord,
                             option=option, verbose=verbose)
    mid = util.get_string("result", comm, req=False, default=mid0)
    system.moleculeList.update({mid : mol2})
    if system.verbose >= 1:
        print('%updateMol: "{}" {}'.format(fname, mol2))

def vector(system, comm):
    Vector.process(system, comm)

def wrapMol(system, comm):
    mol = util.get_molecule("mol", comm, system)
    pads = util.get_value_array("pads", comm)
    cmptlib.wrapMol(mol, np.array(pads))

def writeMol(system, comm):
    mol = util.get_molecule("mol", comm, system)
    fname = util.get_string("file", comm)
    fmt = util.get_string("format", comm, req=False)
    coord = util.get_string("coord", comm, req=False)
    option = util.get_string("option", comm, req=False)
    cmptlib.writeMol(mol, fname, fmt=fmt, coord=coord,
                     option=option, verbose=system.verbose)
    if system.verbose >= 1:
        print('%writeMol: "{}" {}'.format(fname, mol))

