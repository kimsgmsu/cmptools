#!/usr/bin/env python

"""
3D vector object
"""

import math
import numpy as np
import util, Geom
import Molecule

# tolerance for zero
ZTOL = 1.0e-14

class Vector:
    """
    module for a vector object.
    """
    def __init__(self, value):
        """
        Create a 3D vector object

        """
        self._value = np.array(value)

    def copy(self):
        return self.__class__(self.value)

    @property
    def value(self):
        return self._value

    def add(self, other):
        self._value = np.add(self.value, other.value)

    def move(self, delta):
        self._value += delta

    def scale(self, sfac):
        self._value = np.multiply(self.value, sfac)

    def to_fractional(self, mol):
        return mol.to_fractional(self._value)

    def __repr__(self):
        return '<vector value="{}" />'.format(self.value)

    def __str__(self):
        return self.__repr__()

def process(system, node):
    """
    Process a "vector" command
    """
    action = node.get('action', default=None)

    if action == 'new':
        new(system, node)
    elif action == 'add':
        add(system, node)
    elif action == 'copy':
        copy(system, node)
    elif action == 'scale':
        scale(system, node)
    elif action == 'translate':
        translate(system, node)
    elif action == 'report':
        report(system, node)
    else:
        raise Exception('%Vector: Undefined action "{}"'.format(action))
        
def new(system, node):

    method = util.get_string('method', node)
    if method == 'atom':
        mol = util.get_molecule('mol', node, system)
        cell = np.array(util.get_value_array('cell', node, req=False, default=[0,0,0]))
        aid = util.get_int('atom', node)
        value = mol.atom_position_in_cell(aid, cell)

    elif method == 'center-of-vectors':
        vidList = util.get_name_list('vector', node)
        posList = []
        for vid in vidList:
            vec = system.get_vector(vid)
            posList.append(vec.value)
        center, normal = center_of_many_positions(posList)
        shift = util.get_value('shift', node, system, req=False, default=0.0)
        value = center + shift*normal

    elif method == 'pos':
        pos = np.array(util.get_value("pos", comm))
        if coord == 'frac':
            mol = util.get_molecule('mol', node, system)
            pos = mol.to_cartesian(pos)

    elif method == 'value':
        value = np.array(util.get_value('value', node, system))

    elif method == 'point-to-point':
        point1 = util.get_point('point1', node, system)
        point2 = util.get_point('point2', node, system)
        value = point2.pos - point1.pos

    elif method == 'vector-to-vector':
        vec1 = util.get_vector('vec1', node, system)
        vec2 = util.get_vector('vec2', node, system)
        value = vec2.value - vec1.value

    else:
        raise Exception('%Vector: Undefined method "{}"'.format(method))

    vid = util.get_string('result', node, system)
    vec = Vector(value)
    system.vectorList.update({vid : vec})

def add(system, node):
    """
    Add two vectors
    """
    vec1 = util.get_vector('vec1', node, system)
    vec2 = util.get_vector('vec2', node, system)
    vec = vec1.add(vec2)
    vid = util.get_string('result', node, system)
    system.vectorList.update({vid : vec})

def scale(system, node):
    """
    Scale a vector
    """
    vec = util.get_vector('vec', node, system).copy()
    scale = np.array(util.get_value('scale', node, system))
    vec.scale(scale)
    vid = util.get_string('result', node, system)
    system.vectorList.update({vid : vec})

def translate(system, node):
    """
    Translate a vector
    """
    vec = util.get_vector('vec', node, system).copy()
    dvec = np.array(util.get_value('dvec', node, system))
    coord = util.get_string('coord', node, req=False, default='cart')
    if coord == 'frac':
        mol = util.get_molecule('mol', node, system)
        dvec = mol.to_cartesian(dvec)
    vec.move(dvec)
    vid = util.get_string('result', node, system)
    system.vectorList.update({vid : vec})

def copy(system, node):
    """
    Copy a vector
    """
    vec = util.get_vector('vec', node, system)
    vid = util.get_string('result', node, system)
    system.vectorList.update({vid : vec.copy()})

def report(system, node):

    vid = util.get_string('vec', node, system)
    vector = system.get_vector(vid)
    mol = util.get_molecule("mol", node, system, req=False, default=None)
    print('Vector {} = {}'.format(vid, vector), end='')
    if mol is not None:
        fpos = vector.to_fractional(mol)
        print(' ( frac: {})'.format(fpos), end='')
    print('')

def center_of_many_positions(posList):
    """
    Center of several positions
    """
    center = np.zeros(3)
    normal = np.array([0, 0, 1])
    npos = len(posList)
    if npos < 1:
        raise Exception('%Point.center_of_many_positions: empty list of positions')
    elif npos < 2:
        center = posList[0]
        normal = np.array([0, 0, 1])
    elif npos < 3:
        p0 = np.array(posList[0])
        p1 = np.array(posList[1])
        center = (p0 + p1)/2
        normal = util.normalized_vector(p1 - p0)
    else:
        sum_area = 0
        sum_vec = np.zeros(3)
        sum_center = np.zeros(3)
        for kp, pos in enumerate(posList):
            if kp == 0:
                p0 = pos
                continue
            elif kp == 1:
                p2 = pos
                continue
            p1 = p2
            p2 = pos
            t = Geom.Triangle.from_three_points(p0, p1, p2)
            sum_area += t.area
            sum_vec += t.vector
            sum_center += t.center*t.area
        center = sum_center/sum_area
        normal = util.normalized_vector(sum_vec)
    return center, normal

if __name__ == "__main__":
    # stand-alone mode
    vec1 = Vector([1.0, 2.0, 3.0])
    print('vec1 = {}'.format(vec1))
