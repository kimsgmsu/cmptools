#!/usr/bin/env python

"""
rhom2hex: Convert a molecule from rhombohedral unit cell to hexagonal unit cell 
"""

import os
import numpy as np
import Molecule, Atom, Element, Lattice
import util

def rhom2hex(mol0):
    """
    Convert a molecule from rhombohedral to hexagonal unit cell 
    """

    # new lattice vectors
    matt0 = mol0.lattice.matrix
    matt1 = np.eye(3)
    matt1[0] = matt0[0] - matt0[1]
    matt1[1] = matt0[1] - matt0[2]
    matt1[2] = matt0[0] + matt0[1] + matt0[2] 
    latt1 = Lattice.Lattice(matt1, output_scale=mol0.lattice.output_scale)
    print('New hexagonal unit cell:')
    latt1.report()

    # position generators
    dfposList = []
    dfpos = np.array([0, 0, 0])
    dfposList.append(dfpos)
    dfpos = np.array([0, 1, 0])
    dfposList.append(dfpos)
    dfpos = np.array([1, 1, 0])
    dfposList.append(dfpos)

    # generate atom positions
    atomList = []
    for dfpos in dfposList:
        for atom0 in mol0.atomList:
            atom1 = atom0.copy()
            atom1.lattice = latt1
            fpos = atom0.fpos + dfpos
            atom1.cpos = mol0.lattice.to_cpos(fpos)
            atom1.indx = len(atomList)+1
            atomList.append(atom1)

    mol1 = Molecule.Molecule(name=mol0.name, lattice=latt1,
                             elementList=mol0.elementList,
                             atomList=atomList,
                             output_coord=mol0.output_coord,
                             to_unit_cell=True)

    mol1.report()
    return mol1

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('====================')
        print('      RHOM2HEX      ')
        print('====================')

        # read input file
        mol0 = util.read_from_file(args.infile, verbose=args.verbose)

        # convert to hexagonal lattice
        mol1 = rhom2hex(mol0)

        # write output file
        util.write_to_file(mol1, args.outfile, verbose=args.verbose)

if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='rhom2hex: Convert a molecule from rhombohedral to hexagonal unit cell')
    parser.add_argument('infile', help='input file')
    parser.add_argument('outfile', help='output file')
    parser.add_argument('--verbose', type=int, default=4,
                        help='Verbose level index')

    args = parser.parse_args()

    app = myapp()
    app.run(args)
