#!/usr/bin/env python

"""
replicate-molecule: Replicate molecular structure
"""

import numpy as np
import util
import Lattice

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('========================')
        print('   REPLICATE MOLECULE   ')
        print('========================')

        # read input file
        mol0 = util.read_from_file(args.infile, fmt=args.infmt,
                                   verbose=args.verbose)

        # convert nrep parameters
        try:
            nrep = eval(args.nrep)
        except IOError:
            print('Invalid NREP descriptor')
            raise

        # replicate molecule
        mol1 = self.replicate(mol0, nrep)
        mol1.report()

        # write output file
        util.write_to_file(mol1, args.outfile, fmt=args.outfmt,
                           coord=args.coord, verbose=args.verbose)

    def replicate(self, mol0, nrep):

        print('nrep = {}'.format(nrep))
        # copy molecule
        mol1 = mol0.copy()
        # lattice
        matrix = np.eye(3)
        for aid in range(3):
            matrix[aid] = mol0.lattice.matrix[aid]*nrep[aid]

        mol1.lattice = Lattice.Lattice(matrix,
                                       output_scale=mol0.lattice.output_scale)

        # Generate atoms
        mol1.atomList = []
        for c0 in range(nrep[0]):
            for c1 in range(nrep[1]):
                for c2 in range(nrep[2]):
                    dfpos = np.array([c0, c1, c2])
                    for atom0 in mol0.atomList:
                        atom1 = atom0.copy()
                        atom1.lattice = mol1.lattice
                        atom1.indx = len(mol1.atomList)+1
                        fpos = atom0.fpos + dfpos
                        cpos = atom0.lattice.to_cpos(fpos)
                        atom1.cpos = cpos
                        mol1.atomList.append(atom1)

        return mol1

if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='replicate-molecule: Replicate molecular structure')
    parser.add_argument('infile', help='input file')
    parser.add_argument('outfile', help='output file')
    parser.add_argument('--nrep', type=str, default='[2,2,2]',
                        help='No. of cells to replicate in each direction. Given as a string of a valid python int list.  Ex) "[3,1,2]"')
    parser.add_argument('-i', '--infmt',  default=None,
                        help='format for input file')
    parser.add_argument('-o', '--outfmt',  default=None,
                        help='format for output file')
    parser.add_argument('--coord', default=None, choices=['frac', 'cart'],
                        help='coord system for output file')
    parser.add_argument('--verbose', type=int, default=4,
                        help='Verbose level index')

    args = parser.parse_args()

    app = myapp()
    app.run(args)
