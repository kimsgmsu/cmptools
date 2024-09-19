#!/usr/bin/env python

"""
rhom2hex: Convert a molecule from a hexagonal unit cell to a rhombohedral unit cell 
"""

import os
import numpy as np
import Molecule, Atom, Element, Lattice
import util

dTOL = 1.0e-3

def hex2rhom(mol0):
    """
    Convert a molecule from a hexagonal to a rhombohedral unit cell 
    """

    # new lattice vectors
    matt0 = mol0.lattice.matrix
    matt1 = np.eye(3)
    matt1[0] = (2*matt0[0] +matt0[1] +matt0[2])/3.0
    matt1[1] = (-matt0[0] +matt0[1] +matt0[2])/3.0
    matt1[2] = (-matt0[0] -2*matt0[1] +matt0[2])/3.0
    latt1 = Lattice.Lattice(matt1)
    print('New rhombohedral unit cell:')
    latt1.report()

    # Generate atom positions
    atomList = []
    for atom0 in mol0.atomList:
        atom1 = atom0.copy()
        atom1.indx = len(atomList)+1
        atom1.lattice = latt1
        atom1.cpos = atom0.cpos
        atom1.move_to_unit_cell()
        if len(atomList) < 1:
            atomList.append(atom1)
        else:
            min_dist = Atom.min_distance_to_atom_list(atomList, atom1, pbc=True)
            if min_dist > dTOL:
                atomList.append(atom1)

    mol1 = Molecule.Molecule(name=mol0.name, lattice=latt1,
                             elementList=mol0.elementList,
                             atomList=atomList,
                             to_unit_cell=True)

    mol1.report()
    return mol1

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('====================')
        print('      Hex2Rhom      ')
        print('====================')

        # read input file
        mol0 = util.read_from_file(args.infile, verbose=args.verbose)

        # convert to hexagonal lattice
        mol1 = hex2rhom(mol0)

        # write output file
        util.write_to_file(mol1, args.outfile, verbose=args.verbose)

if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='hex2rhom: Convert a molecule from hexagonal to rhombohedral unit cell')
    parser.add_argument('infile', help='input file')
    parser.add_argument('outfile', help='output file')
    parser.add_argument('--verbose', type=int, default=4,
                        help='Verbose level index')

    args = parser.parse_args()

    app = myapp()
    app.run(args)
