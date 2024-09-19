#!/usr/bin/env python

"""
3D point objects
"""

import math
import numpy as np
import util, Geom
import Plane

# tolerance for zero
ZTOL = 1.0e-14

class Point:
    """
    module for a point object.
    """
    def __init__(self, pos):
        """
        Create a 3D point object

        """
        self._pos = np.array(pos)

    def copy(self):
        return self.__class__(pos=self.pos)

    @property
    def pos(self):
        return self._pos

    @pos.setter
    def pos(self, pos):
        self._pos = np.array(pos)

    def move(self, dpos):
        self._pos += dpos

    def distance_from_pos(self, pos):
        return np.linalg.norm(self.pos - pos)

    def distance_from_point(self, pt):
        return self.distance_from_pos(pt.pos)

    def __repr__(self):
        return '<point pos="{}" />'.format(self.pos)

    def __str__(self):
        return self.__repr__()

def process(system, node):
    """
    Process a "point" command
    """
    action = node.get('action', default=None)

    if action == 'new':
        new(system, node)
    elif action == 'copy':
        copy(system, node)
    elif action == 'move':
        move(system, node)
    elif action == 'project':
        project(system, node)
    elif action == 'report':
        report(system, node)
    else:
        raise Exception('%Point: Undefined action "{}"'.format(action))
        
def new(system, node):

    method = util.get_string('method', node)
    if method == 'atom':
        mol = util.get_molecule('mol', node, system)
        cell = np.array(util.get_value_array('cell', node, req=False, default=[0,0,0]))
        aid = util.get_int('atom', node)
        pos = mol.atom_position_in_cell(aid, cell)

    elif method == 'vector':
        vec = util.get_vector('vec', node, system)
        pos = vec.value

    elif method == 'center-of-points':
        pidList = util.get_name_list('point', node)
        posList = []
        for pid in pidList:
            point = system.get_point(pid)
            posList.append(point.pos)
        pos, normal = center_of_many_positions(posList)
        shift = util.get_value('shift', node, system, req=False, default=0.0)
        pos += shift*normal

    elif method == 'center-of-atoms':
        mol = util.get_molecule('mol', node, system)
        aidList = util.get_value_list('atom', node)
        posList = []
        for aid in aidList:
            atom = mol.atom_of_index(aid)
            posList.append(atom.pos)
        pos, normal = center_of_many_positions(posList)
        shift = util.get_value('shift', node, system, req=False, default=0.0)
        pos += shift*normal

    elif method == 'pos':
        pos = np.array(util.get_value('pos', node, system))
        coord = util.get_string('coord', node, req=False, default='cart')
        if coord == 'frac':
            mol = util.get_molecule('mol', node, system)
            pos = mol.to_cartesian(pos)

    else:
        raise Exception('%Point: Undefined method "{}"'.format(method))

    point = Point(pos)
    pid = util.get_string('result', node, system)
    system.pointList.update({pid : point})

def copy(system, node):
    """
    Copy a point
    """
    point = util.get_point('point', node, system).copy()
    pid = util.get_string('result', node, system)
    system.pointList.update({pid : point})

def translate(system, node):
    """
    Translate a point
    """
    point = util.get_point('point', node, system)
    dpos = np.array(util.get_value('dpos', node, system))
    coord = util.get_string('coord', node, req=False, default='cart')
    if coord == 'frac':
        mol = util.get_molecule('mol', node, system)
        dpos = mol.to_cartesian(dpos)
    point.move(dpos)
    pid = util.get_string('result', node, system)
    system.pointList.update({pid : point})

def move(system, node):
    """
    Move a point
    """
    point = util.get_point('point', node, system)
    vec = util.get_vector('vec', node, system)
    point.move(vec.value)

def project(system, node):
    """
    Project a point
    """
    point = util.get_point('point', node, system)
    method = util.get_string('method', node)
    if method == 'to-plane':
        plane = util.get_plane('plane', node, system)
        pos = plane.projected_point(point.pos)

    else:
        raise Exception('%Point.project: Undefined method "{}"'.format(method))

    point.pos = pos
    pid = util.get_string('point', node, system)
    system.pointList.update({pid : point})

def report(system, node):

    pid = util.get_string('point', node, system)
    point = system.get_point(pid)
    mol = util.get_molecule("mol", node, system, req=False, default=None)
    print('Point {} = {}'.format(pid, point), end='')
    if mol is not None:
        fpos = mol.to_fractional(point.pos)
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
    pos1 = Point([1.0, 2.0, 3.0])
