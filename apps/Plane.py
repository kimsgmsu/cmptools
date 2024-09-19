#!/usr/bin/env python

"""
3D plane objects
"""

import math
import numpy as np
import util

# tolerance for zero
ZTOL = 1.0e-10

class Plane:
    """
    module for a plane object.
    """
    def __init__(self, normal, d0):
        """
        Create a 3D plane

        Args:
            normal: (np.array) normal vector to the plane
            d0: (float) distance from the origin
        """
        self._normal = util.normalized_vector(normal)
        self._d0 = d0

    @classmethod
    def from_plane_paramters(cls, a, b, c, d):
        """
        Create a plane from the eqn of a plane:
            a*x + b*y + c*z + d = 0
        """
        norm = math.sqrt(a*a + b*b + c*c)
        if norm < ZTOL:
            raise Exception('%Geom.Plane: invalid plane parameters')
        sfactor = 1.0/norm
        normal = np.array([a, b, c])*sfactor
        d0 = d*sfactor
        return cls(normal, d0)

    @classmethod
    def from_two_vectors_and_dist0(cls, v1, v2, d0):
        """
        Create a plane from two vectors and a point on the plane:
        """
        normal = util.normalized_vector(np.cross(v1, v2))
        return cls(normal, d0)

    @classmethod
    def from_normal_and_pos(cls, normal, pos):
        """
        Create a plane from a normal vector and a position on the plane:
        """
        nv = util.normalized_vector(normal)
        d0 = -np.dot(nv, pos)
        return cls(nv, d0)

    @classmethod
    def from_two_vectors_and_pos(cls, v1, v2, pos):
        """
        Create a plane from three positions on the plane:
        """
        normal = np.cross(v1, v2)
        return cls.from_normal_and_pos(normal, pos)

    @classmethod
    def from_three_pos(cls, posList):
        """
        Create a plane from three positions on the plane:
        """
        p0 = posList[0]
        p1 = posList[1]
        p2 = posList[2]
        v1 = np.array(p1) - np.array(p0)
        v2 = np.array(p2) - np.array(p0)
        return cls.from_two_vectors_and_pos(v1, v2, p0)

    @classmethod
    def from_miller_and_pos(cls, miller, avec, pos):
        """
        Create a plane from miller index and a position
        """
        vol = np.dot(avec[0], np.cross(avec[1], avec[2]))
        bvec = np.zeros((3,3))
        for i1 in range(3):
            i2 = (i1+1) % 3
            i3 = (i1+2) % 3
            bvec[i1] = np.cross(avec[i2], avec[i3]) / vol

        gvec = np.matmul(np.transpose(bvec), miller)
        return cls.from_normal_and_pos(gvec, pos)

    def copy(self):
        return self.__class__(normal=self.normal, d0=self.d0)

    @property
    def normal(self):
        return self._normal

    @normal.setter
    def normal(self, val):
        self._normal = val

    @property
    def d0(self):
        return self._d0

    @d0.setter
    def d0(self, val):
        self._d0 = val

    def move(self, factor=1.0, shift=0.0):
        self.d0 = factor*self.d0 - shift

    def distance_from_point(self, pos):
        """
        Returns distance from a point to the plane.

        Args:
            pos: Cartesian coordinates of the point.

        Returns:
            Signed distance (float)
        """
        return np.dot(self.normal, pos) + self.d0

    def projected_point(self, pos):
        """
        Returns a point that is a projection to the plane.

        Args:
            pos: Cartesian coordinates of the point.

        Returns:
            position of the projected point
        """
        dist = self.distance_from_point(pos)
        return pos - dist*self.normal

    def __repr__(self):
        return '<plane normal="{}" d0="{:>8.3f}" />'.format(self.normal, self.d0)

    def __str__(self):
        return self.__repr__()

def process(system, node):
    """
    Process a "plane" command
    """
    action = node.get('action', default=None)

    if action == 'new':
        new(system, node)
    elif action == 'check':
        check(system, node)
    elif action == 'move':
        move(system, node)
    elif action == 'report':
        report(system, node)
    else:
        raise Exception('%Plane: Undefined action "{}"'.format(action))
        
def new(system, node):

    method = util.get_string('method', node)
    if method == '3-atoms':
        mol = util.get_molecule('mol', node, system)
        aidList = util.get_value_list('atom', node)
        if len(aidList) != 3:
            raise Exception('%Plane: 3-atoms method requires exactly three atoms')
        posList = []
        for aid in aidList:
            atom = mol.atom_of_index(aid)
            posList.append(atom.pos)
        result = Plane.from_three_pos(posList)

    elif method == '3-points':
        pidList = util.get_name_list('points', node)
        if len(pidList) != 3:
            raise Exception('%Plane: 3-points method requires exactly three points')
        posList = []
        for pid in pidList:
            point = system.get_point(pid)
            posList.append(point.pos)
        result = Plane.from_three_pos(posList)

    elif method == 'normal-and-atom':
        normal = util.get_value('normal', node, system)
        mol = util.get_molecule("mol", node, system)
        coord = util.get_string("coord", node, req=False, default='cart')
        if coord == 'frac':
            normal = mol.to_cartesian(normal)
        atom = util.get_atom('atom', node, mol)
        result = Plane.from_normal_and_pos(normal, atom.pos)

    elif method == 'normal-and-point':
        normal = util.get_value('normal', node, system)
        coord = util.get_string("coord", node, req=False, default='cart')
        if coord == 'frac':
            mol = util.get_molecule("mol", node, system)
            normal = mol.to_cartesian(normal)
        point = util.get_point('point', node, system)
        result = Plane.from_normal_and_pos(normal, point.pos)

    elif method == 'miller-and-atom':
        miller = np.array(util.get_value('miller', node, system))
        mol = util.get_molecule('mol', node, system)
        atom = util.get_atom('atom', node, mol)
        result = Plane.from_miller_and_pos(miller, mol.lattice.matrix, atom.pos)

    elif method == 'miller-and-point':
        miller = np.array(util.get_value('miller', node, system))
        mol = util.get_molecule('mol', node, system)
        point = util.get_point('point', node, system)
        result = Plane.from_miller_and_pos(miller, mol.lattice.matrix, point.pos)

    else:
        raise Exception('%Plane: Undefined method "{}"'.format(method))

    pid = util.get_string('result', node, system)
    shift = util.get_float('shift', node, req=False, default=0.0)
    result.move(shift=shift)
    system.planeList.update({pid : result})

def check(system, node):

    pid = util.get_string('plane', node, system)
    plane = system.get_plane(pid)
    item = util.get_string('item', node)
    if item == 'distance-to-atom':
        mol = util.get_molecule('mol', node, system)
        aid = util.get_int('atom', node)
        atom = mol.atom_of_index(aid)
        cpos = mol.to_cpos(atom.fpos)
        dist = plane.distance_from_point(cpos)
        print('Plane.check: distance to atom from a plane:')
        print('Plane {} = {}'.format(pid, plane))
        print('Atom {} = {}'.format(aid, atom))
        print('distance = {}'.format(dist))
    else:
        raise Exception('%Plane: invalid item to check "{}"'.format(item))

def report(system, node):

    pid = util.get_string('plane', node, system)
    plane = system.get_plane(pid)
    header = util.get_string('header', node, req=False)
    if header is not None:
        print(header)
    print('Plane {} = {}'.format(pid, plane))

def move(system, node):
    """
    Move a plane
    """
    plane = util.get_plane('plane', node, system)
    factor = util.get_float('factor', node)
    shift = util.get_float('shift', node)
    plane.move(factor=factor, shift=shift)

if __name__ == "__main__":
    # stand-alone mode
    plane1 = Plane([1.0, 2.0, 3.0], 2.0)
    print('plane1 = ', plane1)

    p1 = np.array([1, 0, 0])
    p2 = np.array([1, 1, 0])
    p3 = np.array([1, 0, 1])
    plane2 = Plane.from_two_vectors_and_pos(p2-p1, p3-p1, p1)
    print('plane2 = ', plane2)

    p1 = np.array([1, 0, 0])
    p2 = np.array([0, 1, 0])
    p3 = np.array([1, 0, 1])
    plane2 = Plane.from_two_vectors_and_pos(p2-p1, p3-p1, p1)
    print('plane2 = ', plane2)

    plane3 = Plane.from_three_pos((p1, p2, p3))
    print('plane3 = ', plane3)
    
    a = 1.0
    c = 2.0
    a1 = np.array([1, 0, 0])*a
    a2 = np.array([-1.0, math.sqrt(3), 0])*(a/2)
    a3 = np.array([0, 0, 1])*c
    avec = np.array([a1, a2, a3])
    print('avec=',avec)
    plane4 = Plane.from_miller_and_pos([1, 0, 0], avec, [0.3, 0, 0])
    print('plane4 = ', plane4)
    plane5 = Plane.from_miller_and_pos([0, 1, 0], avec, [0, 0.5, 0])
    print('plane5 = ', plane5)
    plane6 = Plane.from_miller_and_pos([0, 0, 1], avec, [0, 0, 0.7])
    print('plane6 = ', plane6)
    plane7 = Plane.from_miller_and_pos([1, 0, 1], avec, [1, 0, 0])
    print('plane7 = ', plane7)
