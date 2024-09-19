#!/usr/bin/env python

import numpy as np
import util

position_atol = 1.0e-6

class Atom:
    """
    module for an Atom object.
    """
    def __init__(self):
        """
        Create an Atom.
        """
        self.symbol = ''
        self.molecule = None
        self.pos = None
        self.sdyn = None
        self.info = ''

    def copy(self):
        cp = self.__class__()
        cp.symbol = self.symbol
        cp.molecule = self.molecule
        cp.pos = self.pos
        cp.sdyn = self.sdyn
        cp.info = self.info
        return cp

    def move_to_unit_cell(self):
        """
        Move the atom to inside the unit cell
        """
        self.fpos = util.reduced_vector(self.fpos, pbc=True)

    def distance(self, pos, coord, latt):
        """
        Distance from a given position
        """
        if coord == 'cart':
            pos0 = pos
        else:
            pos0 = latt.to_cartesian(pos)

        dpos = self.pos - pos0
        dfpos = latt.to_fractional(dpos)
        rdfpos = util.reduced_vector(dfpos, latt.pbc)
        dpos = latt.to_cartesian(rdfpos)
        return np.linalg.norm(dpos)

    def __repr__(self):
        ffmt= '{:>10.6f}'
        return '<atom symbol="{}" pos="{}" info="{}" />'.format(
            self.symbol, util.list_to_str(self.pos, ffmt), self.info)

    def __str__(self):
        return self.__repr__()

if __name__ == "__main__":
    # stand-alone mode
    fpos = [1, 0, 0]
    atom1 = Atom('H', latt, coord=fpos, info='atom1')
    print('atom1=', atom1)
    atom2 = Atom('C', latt, coord=[2.1, 0.5, -1.8], cart=True, info='atom2')
    print('atom2=', atom2)

    atom1.fpos = [0, 1, 0]
    print('after reassignment, atom1={}, cpos={}'.format(atom1, atom1.cpos))
    atom1.cpos = [0, 1.0, 0]
    print('after reassignment, atom1={}, cpos={}'.format(atom1, atom1.cpos))
    
