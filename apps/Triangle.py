#!/usr/bin/env python

"""
3D triangle objects
"""

import math
import numpy as np
import util

# tolerance for zero
ZTOL = 1.0e-14

class Triangle:
    """
    module for a 3D triangle object.
    """
    def __init__(self, vertex):
        """
        Create a triangle
        """
        self._vertex = np.array(vertex)
        v1 = vertex[1] - vertex[0]
        v2 = vertex[2] - vertex[0]
        self._vector = np.cross(v1, v2)
        norm = np.linalg.norm(self._vector)
        try:
            self._normal = self._vector/norm
        except:
            raise Exception('%Triangle: Invalid vertices for a triangle "{}"'.format(vertex))
        self._area = 0.5*norm
        self._centroid = vertex[0] + (v1+v2)*(1/3)

    @property
    def vertex(self):
        return self._vertex

    @property
    def area(self):
        return self._area

    @property
    def normal(self):
        return self._normal

    @property
    def centroid(self):
        return self._centroid

    @property
    def center(self):
        return self.centroid

    def copy(self):
        return self.__class__(vertex=self.vertex)

    def __repr__(self):
        s = '<triangle>\n'
        s += ' <vertex="{}"/>\n'.format(self.vertex)
        s += ' area="{}" normal="{}" centroid="{}"\n'.format(
            self.area, self.normal, self.centroid)
        s += '</triangle>'
        return s

    def __str__(self):
        return self.__repr__()

    @classmethod
    def from_two_vectors(cls, v1, v2):
        """
        Create a triangle from two vectors
        """
        return cls(np.stack((np.zeros(3), v1, v2)))

    @classmethod
    def from_three_points(cls, p1, p2, p3):
        return cls(np.stack((p1, p2, p3)))

    @classmethod
    def from_point_and_two_vectors(cls, pt, v1, v2):
        """
        Create a triangle from a point and two vectors
        """
        return cls.from_three_points(np.stack((pt,
                                               np.array(pt)+np.array(v1),
                                               np.array(pt)+np.array(v2))))

def process(system, node):
    """
    Process a "triangle" command
    """
    action = node.get('action', default=None)
    print('action = "{}"'.format(action))

    if action == 'new':
        new(system, node)
    elif action == 'check':
        check(system, node)
    elif action == 'move':
        move(system, node)
    elif action == 'report':
        report(system, node)
    else:
        raise Exception('%Triangle: Undefined action "{}"'.format(action))
        
def new(system, node):

    method = util.get_string('method', node)
    if method == 'three-atoms' or method == '3-atoms':
        mol = util.get_molecule('mol', node, system)
        aidList = util.get_list('atom', node)
        if len(aidList) != 3:
            raise Exception('%Triangle: 3-atoms method requires exactly three atoms')
        pts = np.zeros((3,3))
        for k,aid in enumerate(aidList):
            atom = mol.atom_of_index(aid)
            pts[k] = atom.pos
        result = Triangle.from_three_points(pts[0], pts[1], pts[2])

    else:
        raise Exception('%Triangle: Undefined method "{}"'.format(method))

    tid = util.get_string('result', node, system)
    system.triangleList.update({tid : result})

def report(system, node):

    tid = util.get_string('triangle', node, system)
    triangle = system.get_triangle(tid)
    header = util.get_string('header', node, req=False)
    if header is not None:
        print(header)
    print('Triangle {} = {}'.format(tid, triangle))

def move(system, node):
    """
    Move a triangle
    """
    triangle = util.get_triangle('triangle', node, system)
    factor = util.get_float('factor', node)
    shift = util.get_float('shift', node)
    triangle.move(factor=factor, shift=shift)

if __name__ == "__main__":
    # stand-alone mode
    triangle1 = Triangle([1.0, 2.0, 3.0], 2.0)
    print('triangle1 = ', triangle1)

    p1 = np.array([1, 0, 0])
    p2 = np.array([1, 1, 0])
    p3 = np.array([1, 0, 1])
    triangle2 = Triangle.from_two_vectors_and_point(p2-p1, p3-p1, p1)
    print('triangle2 = ', triangle2)

    p1 = np.array([1, 0, 0])
    p2 = np.array([0, 1, 0])
    p3 = np.array([1, 0, 1])
    triangle2 = Triangle.from_two_vectors_and_point(p2-p1, p3-p1, p1)
    print('triangle2 = ', triangle2)

    triangle3 = Triangle.from_three_points(p1, p2, p3)
    print('triangle3 = ', triangle3)
    
    a = 1.0
    c = 2.0
    a1 = np.array([1, 0, 0])*a
    a2 = np.array([-1.0, math.sqrt(3), 0])*(a/2)
    a3 = np.array([0, 0, 1])*c
    latt = np.array([a1, a2, a3])
    print('latt=',latt)
    triangle4 = Triangle.from_miller_index([1, 0, 0], latt)
    print('triangle4 = ', triangle4)
    triangle5 = Triangle.from_miller_index([0, 1, 0], latt)
    print('triangle5 = ', triangle5)
    triangle6 = Triangle.from_miller_index([0, 0, 1], latt)
    print('triangle6 = ', triangle6)
    triangle6 = Triangle.from_miller_index([1, 0, 1], latt)
    print('triangle6 = ', triangle6)
