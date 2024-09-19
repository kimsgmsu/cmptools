#!/usr/bin/env python

"""
transplant: Transplant an atom using three base atoms
"""

import numpy as np
import util, cmptlib

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('================')
        print('   TRANSPLANT   ')
        print('================')

        # Original positions of base atoms
        list0 = []
        print('Enter the original positions of base atoms:')
        for iatom in range(3):
            line = input('pos({})?'.format(iatom))
            pos = np.array(cmptlib.str_to_list(line), dtype=float)
            list0.append(pos)

        # position of atom to transplant
        line = input('pos of the atom to transplant?')
        pos0 = np.array(cmptlib.str_to_list(line), dtype=float)
        
        # New positions of base atoms
        list1 = []
        print('Enter the new positions of base atoms:')
        for iatom in range(3):
            line = input('pos({})?'.format(iatom))
            pos = np.array(cmptlib.str_to_list(line), dtype=float)
            list1.append(pos)

        # Pad the data with ones, so that our transformation can do
        # translations too
        array0 = np.array(list0)
        print('array0=',array0)
        print('pos0=',pos0)
        array1 = np.array(list1)
        print('array1=',array1)

        pad = lambda x: np.hstack([x, np.ones((x.shape[0], 1))])
        unpad = lambda x: x[:,:-1]
        X = pad(array0)
        Y = pad(array1)

        # Solve the least squares problem X * A = Y
        # to find our transformation matrix A
        A, res, rank, s = np.linalg.lstsq(X, Y, rcond=None)

        transform = lambda x: unpad(np.dot(pad(x), A))

        print('Target:')
        print(array1)
        print('Result:')
        print(transform(array0))
        print('Max error:', np.abs(array1 - transform(array0)).max())

        apos0 = np.array(pos0).reshape(1,3)
        print('Transplanted pos=',transform(apos0))

if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='transplant: Transplant an atom using three base atoms')
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
