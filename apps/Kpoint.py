#!/usr/bin/env python

import numpy as np
import util

position_tol = 1.0e-6

class Kpoint:
    """
    module for a K-point object.
    """
    def __init__(self, rpos, name=''):
        """
        Create a Kpoint.
        """
        self.name = name
        self.rpos = rpos

    def cpos(self, mol):
        return mol.recip_to_cart(self.rpos)

    def leng(self, mol):
        return np.linalg.norm(np.array(self.cpos(mol)))

    def unit_vector(self, mol):
        cpos = self.cpos(mol)
        try:
            uvec = cpos / self.leng(mol)
        except:
            uvec = cpos
        return uvec

    def copy(self):
        cp = self.__class__()
        cp.name = self.name
        cp.rpos = self.rpos
        return cp

    def __repr__(self):
        ffmt= '{:>10.6f}'
        return '<kpoint name="{}" rpos="{}" />'.format(
            self.name, util.list_to_str(self.rpos, ffmt))

    def __str__(self):
        return self.__repr__()

    def report(self, mol=None, header=None):
        if header is not None:
            print(header)
        ffmt = '{:>10.6f}'
        print('<kpoint name="{}" rpos="{}"'.format(
            self.name, util.list_to_str(self.rpos, ffmt)))
        if mol is not None:
            leng = self.leng(mol)
            cpos = self.cpos(mol)
            uvec = self.unit_vector(mol)
            print('     leng="{}" cpos="{}" uvec="{}" />'.format(
                leng, util.list_to_str(cpos, ffmt),
                util.list_to_str(uvec, ffmt)))

def process(system, node):
    """
    Process a "kpoint" command
    """
    action = node.get('action', default=None)

    if action == 'new':
        new(system, node)
    elif action == 'report':
        report(system, node)
    else:
        raise Exception('%Kpoint: Undefined action "{}"'.format(action))
        
def new(system, node):

    method = util.get_string('method', node)
    if method == 'pos':
        pos = np.array(util.get_value('pos', node, system))
        coord = util.get_string('coord', node, req=False, default='frac')
        if coord == 'cart':
            mol = util.get_molecule('mol', node, system)
            rpos = mol.cart_to_recip(pos)
        elif coord == 'frac':
            rpos = pos
        else:
            raise Exception('%Kpoint: Undefined coord "{}"'.format(coord))
    else:
        raise Exception('%Kpoint: Undefined method "{}"'.format(method))

    name = util.get_string('name', node, req=False)
    kpoint = Kpoint(rpos=rpos, name=name)
    kid = util.get_string('result', node, system)
    system.kpointList.update({kid : kpoint})

def report(system, node):

    kpoint = util.get_kpoint('kpoint', node, system)
    mol = util.get_molecule('mol', node, system, req=False)
    header = util.get_string('header', node, req=False)
    kpoint.report(mol,header)

if __name__ == "__main__":
    # stand-alone mode
    import Lattice, Molecule
    import numpy as np
    latt = Lattice.Lattice([[1.1, 0.12, 0.13], [0.21, 2.2, 0.23], [0.31, 0.32, 3.3]])
    latt.report()
    mol = Molecule.Molecule(latt=latt)
    mol.report()
    kpoint1 = Kpoint([0.5, 0, 0], name='M')
    kpoint1.report(mol=mol, header='kpoint1:')
    kpoint2 = Kpoint([1, 0, 0], name='2M')
    kpoint2.report(mol=mol, header='kpoint2:')
    kpoint3 = Kpoint([0.3333, 0.3333, 0], name='K')
    kpoint3.report(mol=mol, header='kpoint3:')
