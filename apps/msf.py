#!/usr/bin/env python
"""
Module for msf file I/O
"""

from __future__ import print_function
import xml.etree.ElementTree as ET
import numpy as np
import Molecule, Atom, Lattice
import util

def read(f, verbose=0):

    # get the element tree
    tree = ET.parse(f)
    if verbose > 8:
        print("%msf.read: printing nodes")
        for node in tree.iter():
            print(node.tag, node.attrib)

    # root
    root = tree.getroot()
    if root.tag != "molecule":
        raise Exception('<molecule> tag missing. Not an msf file.')
        
    # read molecule
    return read_molecule(tree)

def read_molecule(xmltree, VARS=None):

    print('before a mol ctreated')
    mol = Molecule.Molecule()
    print('mol created')

    # Header
    name = ''
    node = xmltree.find("name")
    if node != None:
        if node.text != None:
            name = node.text.strip()
    mol.name = name

    info = ''
    node = xmltree.find("info")
    if node != None:
        if node.text != None:
            info = node.text.strip()
    mol.info = info

    # Lattice
    mol.latt = read_lattice(xmltree, VARS)

    # Element list
    read_elementList(xmltree, mol)

    # Atom list
    mol.atomList, mol.output_coord = read_atomList(xmltree, mol)

    # Volumetric data
    mol.vdList = read_volumetricDataList(xmltree)

    # complete element table
    mol.make_ElementTable()

    return mol

def read_lattice(xmltree, VARS=None):

    xmlnode = xmltree.find('lattice')
    if xmlnode is None:
        raise Exception("<lattice> tag missing")
    method = util.get_string('method', xmlnode, req=False, default='vectors')
    pbc = util.get_bool('pbc', xmlnode, req=False, default=True)

    latt_scale = read_lattice_scale(xmlnode)

    if method == "vectors":
        matrix = np.eye(3)
        for xmlitem in xmlnode:
            if xmlitem.tag == "axis":
                aname = xmlitem.get("name")
                try:
                    aid = ["a1","a2","a3"].index(aname)
                except ValueError:
                    raise Exception(
                        '%msf-error: Invalid axis name "{}"'.format(aname))
                svec = xmlitem.get("vector")
                if svec == None:
                    raise Exception('%msf-error: "vector" is missing in <axis> tag')
                avec = np.array(util.str_to_list(svec, VARS), dtype=float)
                matrix[aid] = avec*latt_scale[aid]
        lattice = Lattice.Lattice(matrix, pbc)

    elif method == 'parameters':
        abc = np.array([1.0, 1.0, 1.0])
        angles = np.array([90.0, 90.0, 90.0])
        for xmlitem in xmlnode:
            if xmlitem.tag == "lengths":
                svec = xmlitem.get("vector")
                if svec == None:
                    raise Exception('%msf-error: "vector" is missing in <lengths> tag')
                lengths = np.array(util.str_to_list(svec))*latt_scale
            elif xmlitem.tag == "angles":
                svec = xmlitem.get("vector")
                if svec == None:
                    raise Exception('%msf-error: "vector" is missing in <angles> tag')
                angles = np.array(util.str_to_list(svec), dtype=float)

        lattice = Lattice.Lattice.from_parameters(lengths, angles, pbc)

    return lattice

def read_lattice_scale(xmltree):

    scale = np.ones(3)
    xmlnode = xmltree.find('scale')
    if xmlnode is None:
        return scale

    s = xmlnode.get("value", default="1.0 1.0 1.0")
    try:
        v = util.str_to_list(s)
    except IOError:
        Error('Invalid list value "{}"'.format(s))

    for idir in range(3):
        try:
            scale[idir] = v[idir]
        except:
            scale[idir] = v[0]

    return scale

def read_elementList(tree, mol):

    xmlnode = tree.find("elementList")
    if xmlnode == None:
        return

    for elem in xmlnode.iter():
        if elem.tag == 'element':
            symbol = util.get_string('symbol', elem)
            radius = util.get_float('radius', elem, req=False, default=None)
            mass = util.get_float('mass', elem, req=False, default=None)
            info = util.get_string('info', elem, req=False, default='')
            element = Molecule.Element(symbol, info=info, mass=mass, radius=radius)
            mol.elementTable.update({symbol : element})

def read_atomList(tree, mol):

    xmlnode = tree.find("atomList")
    if xmlnode == None:
        raise Exception("<atomList> tag missing")
    coord = xmlnode.get("coord")
    if coord == None:
        raise Exception('"<atomList coord" attribute missing')
    elif coord == 'frac':
        cart = False
    elif coord == 'cart':
        cart = True
    else:
        raise Exception('%msf-error: invalid coord value "{}"'.format(coord))

    atomList = []
    for elem in xmlnode.iter():
        if elem.tag == 'atom':
            atom = Atom.Atom()
            atom.symbol = util.get_string('symbol', elem)
            pos = np.array(util.get_value_list('pos', elem))
            if coord == 'cart':
                atom.pos = pos
            else:
                atom.pos = mol.to_cartesian(pos)
            atom.sdyn = util.get_string('sdyn', elem, req=False)
            atom.info = util.get_string('info', elem, req=False)
            atomList.append(atom)

    return atomList, coord

def read_volumetricDataList(tree):

    vdList = []
    xmlnode = tree.find('volumetricDataList')
    if xmlnode is None:
        return vdList

    for elem in xmlnode.iter():
        if elem.tag == 'volumetricData':
            vd = read_volumetricData(elem)
            vdList.append(vd)

    return vdList

def read_volumetricData(xmltree):

    name = xmltree.get('name')
    if name is None:
        raise Exception('%msf-error: "name" is missing in <volumetricData> tag')

    # main data
    elem = xmltree.find('data')
    if elem == None:
        raise Exception('<data> tag missing in <volumetricData> xmltree')
    ndim = np.array(util.get_value_list('dim', elem), dtype=int)
    
    #read data on grid points
    data = np.zeros( (ndim[0], ndim[1], ndim[2]), dtype=np.float_ )
    words = elem.text.strip().split()
    for i2 in range(ndim[2]):
        for i1 in range(ndim[1]):
            for i0 in range(ndim[0]):
                try:
                    data[i0,i1,i2] = np.float_(words.pop(0))
                except:
                    raise Exception(
                        '%msf: Invalid volumetric data at "({},{},{})"'.format(i0,i1,i2))

    # extra data
    elem = xmltree.find('extra')
    if elem == None:
        extra = None
    else:
        extra = elem.text.strip()

    return Molecule.VolumetricData(name=name, data=data, extra=extra)

def write(mol, f, coord=None, verbose=3):

    # Header
    f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    f.write('<molecule>\n')

    # structure
    write_structure(mol, f, coord, verbose)

    # volumetric data
    write_volumetricDataList(mol, f)

    # trailer
    f.write('</molecule>\n')

def write_structure(mol, f, coord=None, verbose=3):

    # Header
    f.write('  <name>{}</name>\n'.format(mol.name))
    f.write('  <info>{}</info>\n'.format(mol.info))

    # Lattice
    ffmt = ' {:>10.6f}'
    latt = mol.latt
    f.write('  <lattice method="vectors" pbc="{}">\n'.format(latt.pbc))
    f.write('    <scale value="[1.0, 1.0, 1.0]" />\n')
    for i in range(3):
        v = latt.matrix[i]
        f.write('    <axis name="a{}" vector="{}" />\n'.format(
            i+1, util.list_to_str(v, ffmt)))
    f.write('  </lattice>\n')

    # Element list
    f.write('  <elementList>\n')
    for symbol in mol.distinct_atomic_symbols:
        element = mol.get_Element(symbol)
        s = '    <element symbol="{}"'.format(symbol)
        s += ' radius="{}"'.format(element.radius)
        s += ' mass="{}"'.format(element.mass)
        s += ' info="{}" />\n'.format(element.info)
        f.write(s)
    f.write('  </elementList>\n')
        
    # Atom list
    if coord == None:
        coord = mol.output_coord
    f.write('  <atomList coord="{}" >\n'.format(coord))
    for ka, atom in enumerate(mol.atomList):
        s = '    <atom id="{}"'.format(ka+1)
        s += ' symbol="{}"'.format((atom.symbol))
        if coord == 'cart':
            pos = atom.pos
        else:
            pos = mol.to_fractional(atom.pos)
        s += ' pos="{}"'.format(util.list_to_str(pos, ffmt))
        if atom.sdyn is not None:
            s += ' sdyn="{}"'.format(str(atom.sdyn))
        s += ' info="{}" />\n'.format(atom.info)
        f.write(s)
    f.write('  </atomList>\n')

def write_volumetricDataList(mol, f):

    if mol.volumeDataList is None or len(mol.volumeDataList) < 1:
        return

    # volumetric data
    f.write('<volumetricDataList>\n')
    for vd in mol.volumeDataList:
        write_volumetricData(f, vd)
    f.write('</volumetricDataList>\n')

def write_volumetricData(f, vd):

    if vd is None:
        return

    f.write('<volumetricData name="{}">\n'.format(vd.name))
    dim = vd.data.shape
    f.write('  <data dim="{}">\n'.format(util.list_to_str(dim)))
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
    f.write('  </data>\n')

    # extra data
    if vd.extra is not None:
        f.write('  <extra>\n')
        f.write(vd.extra+'\n')
        f.write('  </extra>\n')
    f.write('</volumetricData>\n')
    
if __name__ == "__main__":
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='msf.py: Read and write a msf format file')
    parser.add_argument('infile', help='input file')
    parser.add_argument('outfile', help='output file')
    parser.add_argument('--coord', default=None, choices=['frac', 'cart'],
                        help='coord system for output file')

    args = parser.parse_args()

    f = open(args.infile, 'r')
    mol = read(f, verbose=5)
    f.close()
    f = open(args.outfile, 'w')
    write(mol, f, verbose=5)
    f.close()
