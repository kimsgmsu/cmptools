#!/usr/bin/env python3

import sys
from xml.etree import ElementTree as ET
import os
import Molecule

def read(f, verbose=0):
  molecule = Molecule.Molecule()
  molecule.unitlength = "angstrom"
  molecule.unitcell = Molecule.Unitcell()
  molecule.atomcoord = "cartesian"
  molecule.atomlist = []
  for line in f.readlines():
    if verbose > 4:
      print('line = "'+line.rstrip()+'"')
    key = line[:6].strip()
    if verbose > 4:
      print('key   = "'+key+'"')
    if key == "COMPND":
      molecule.name = line[6:].strip()
    elif key == "CRYST1":
      molecule.unitcell = _read_unitcell(line)
    elif key == "ATOM":
      atom = _read_atom(line, molecule.unitcell.c2f)
      molecule.atomlist.append(atom)

  # report the molecule molecule read
  molecule.report()
  return molecule

def _read_unitcell(line):
  uc = Molecule.Unitcell()
  # primitive axes
  a = float(line[9:16])
  b = float(line[16:25])
  c = float(line[25:34])
  alpha = float(line[34:41])
  beta = float(line[41:48])
  gamma = float(line[48:55])
  uc.axis = cmptlib.box_vectors(numpy.array([a,b,c,alpha,beta,gamma]))
  uc.setup()
  return uc

def _read_atom(line, c2f):
  atom = Molecule.Atom()
  s = line[7:12].strip()
  atom.id = int(s)
  atom.name = line[13:17].strip()
  s1 = line[31:39]
  s2 = line[39:47]
  s3 = line[47:54]
  pos = numpy.array([float(s1),float(s2),float(s3)])
  atom.pos = numpy.dot(c2f, pos)
  return atom

def write(molecule, f, verbose=0):
  # COMPND
  s = 'COMPND ' + molecule.name + '\n'
  f.write(s)
  # CRYST1
  s = _write_cryst1(molecule.unitcell.axis)
  f.write(s)
  # atomlist
  for atom in molecule.atomlist:
    pos = np.dot(molecule.unitcell.f2c, atom.fpos)
    s = _write_atom(atom, pos)
    f.write(s)
  # TER
  s = 'TER   {:>5}'.format(len(molecule.atomlist)+1)
  s += '      UNK\n'
  f.write(s)
  # END
  f.write("END\n")

def _write_cryst1(axis):
  dims = cmptlib.box_dimensions(axis)
  print("dims=", dims[:])
  # Record name
  s = 'CRYST1'
  for d in dims:
    s += '{:9.3f}'.format(d)
  return s+'\n'

def _write_atom(atom, pos):
  # Record name
  s = 'ATOM  '
  # serial number
  s += '{:>5}'.format(atom.id)
  # space
  s += ' '
  # atom name
  s += '{:>4.4}'.format(atom.name)
  # alt. location
  s += ' '
  # residue name
  s += 'UNK'
  # space
  s += ' '
  # chain id
  s += ' '
  # residue seq number
  s += '{:>5}'.format(0)
  # 3 spaces
  s += '   '
  # coordinates
  for aid in range(3):
    s += "{:8.3f}".format(pos)
  # occupancy
  s += '{:6.2f}'.format(1.0)
  # temp factor
  s += '{:6.2f}'.format(0.0)
  # 6 spaces
  s += '      '
  # seg. id
  s += '    '
  # element
  s += '{:>2}'.format(atom.element)
  # charge
  s += '  '
  return s+'\n'

if __name__ == "__main__":
  # stand-alone mode
  if len(sys.argv) < 2:
    print('pdb: Usage: pdb.py infile [outfile]')
  else:
    infile = sys.argv[1]
    f = open(infile, 'r')
    print('pdb: Reading input PDB file "' + infile + '"')
    molecule = read(f, verbose=5)
    f.close()

  if len(sys.argv) > 1:
    outfile = sys.argv[2]
    print('pdb: Writing output PDB file "' + outfile + '"')
    f = open(outfile, "w")
    write(molecule, f, verbose=5)
    f.close()
