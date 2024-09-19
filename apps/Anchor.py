#!/usr/bin/env python

"""
Anchor  objects
"""

import math
import numpy as np
import util
from Geom import Tetrahedron
import Atom

# tolerance for zero
ZTOL = 1.0e-10

class Anchor(Tetrahedron):
    """
    module for a anchor object.
    """
    def __init__(self, vertex):
        """
        Create a anchor
        """
        super().__init__(vertex)

    def __repr__(self):
        return '<anchor vertex="{}" volume="{}" pos="{}" dir="{}" />'.format(
            self.vertex, self.volume, self.position, self.direction)

    @classmethod
    def from_points(cls, pList, direction=np.array([0,0,1]), scale=1.0, shift=0.0):
        """
        Create a anchor from list of points.
        """
        if len(pList) == 4:
            return cls.from_four_points(pList[0], pList[1], pList[2], pList[3],
                                          scale=scale, shift=shift)
        elif len(pList) == 3:
            return cls.from_three_points(pList[0], pList[1], pList[2], 
                                          scale=scale, shift=shift)
        elif len(pList) == 2:
            return cls.from_two_points(pList[0], pList[1],
                                       scale=scale, shift=shift)
        elif len(pList) == 1:
            return cls.from_one_point(pList[0], direction=direction,
                                      scale=scale, shift=shift)
        else:
            raise Exception('%anchor: atoms method requires one to four atoms')

    @classmethod
    def from_four_points(cls, p0, p1, p2, p3, scale=1.0, shift=0.0):
        """
        Create a anchor from three points.
        """
        tet = super().from_four_points(p0, p1, p2, p3)
        apex = tet.center + (tet.apex - tet.center)*scale
        tet = super().from_four_points(p0, p1, p2, apex)
        return tet.move_by(tet.direction*shift)

    @classmethod
    def from_three_points(cls, p0, p1, p2, scale=1.0, shift=0.0):
        """
        Create a anchor from three points.
        """
        tet = super().from_three_points(p0, p1, p2, scale=scale)
        return tet.move_by(tet.direction*shift)

    @classmethod
    def from_two_points(cls, p0, p1, scale=1.0, shift=0.0):
        """
        Create a anchor from two points.
        """
        tet = super().from_two_points(p0, p1, scale=scale)
        return tet.move_by(tet.direction*shift)

    @classmethod
    def from_one_point(cls, point, direction=None, scale=1.0, shift=0.0):
        """
        Create a anchor from a point.
        """
        if direction is None:
            direction = np.array([0, 0, 1])
        tet = super().from_one_point(point, direction, scale=scale)
        return tet.move_by(tet.direction*shift)

def process(system, node):
    """
    Process a "anchor" command
    """
    action = node.get('action', default=None)
    print('action = "{}"'.format(action))

    if action == 'new':
        new(system, node)
    elif action == 'check':
        check(system, node)
    elif action == 'scale':
        scale(system, node)
    else:
        raise Exception('%makemol.anchor: Undefined action "{}"'.format(action))

def new(system, node):

    method = util.get_string('method', node)
    if method == 'atom':
        mol = util.get_molecule('mol', node, system)
        aidList = util.get_value_list('atom', node, system)
        posList = []
        for aid in aidList:
            atom = mol.atom_of_index(aid)
            if atom is None:
                util.Error('%anchor.new: Invalid atom index, "{}"'.format(aid))
            posList.append(atom.pos)

    elif method == 'vector':
        vidList = util.get_name_list('vector', node, system)
        posList = []
        for vid in vidList:
            vec = system.get_vector(vid)
            if vec is None:
                util.Error('%anchor.new: Invalid vector, "{}"'.format(vid))
            posList.append(vec.value)

    else:
        raise Exception('%anchor.new: Undefined method "{}"'.format(method))

    # generate an anchor
    direction = util.get_value_list('direction', node, system, req=False, default=None)
    scale = util.get_float('scale', node, req=False, default=1.0)
    shift = util.get_float('shift', node, req=False, default=0.0)
    anchor = Anchor.from_points(posList, direction=direction,
                                scale=scale, shift=shift)
    aid = util.get_string('result', node, system)
    system.anchorList.update({aid : anchor})

if __name__ == "__main__":
    # stand-alone mode
    p0 = np.array([0, 0, 0])
    p1 = np.array([1, 0, 0])
    p2 = np.array([0, -1, 0])
    p3 = np.array([0, 0, 1])
    vertex = np.stack((p0, p1, p2, p3))
    anchor1 = Anchor(vertex)
    print('anchor1 = ', anchor1)

    anchor2 = Anchor.from_four_points(p0, p1, p2, p3)
    print('anchor2 = ', anchor2)

    p0 = np.array([0, 0, 1])
    p1 = np.array([1, 0, 0])
    p2 = np.array([1, -1, 0])
    anchor3 = Anchor.from_three_points(p0, p1, p2)
    print('anchor3 = ', anchor3)
