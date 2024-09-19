#!/usr/bin/env python

"""
3D tetrahedron  objects

    A 3D object is manipulated by a tertrahedron associated with it.
    (e.g. a cube --> a tetrahedron made of three vertices of the base 
                     + one vertex closest to the first vertex
          a sphere --> a terahedron for the circumscribed cube)

"""

import math
import numpy as np
import util

# small number to be considered as zero
ZTOL = 1.0e-14

class Arrow:
    """
    module for an Arrow object.
    """
    def __init__(self, vertex):
        """
        Create an Arrow
        """
        self._vertex = np.array(vertex)
        self._position = vertex[0]
        self._vector = vertex[1] - vertex[0]
        self._length = np.linalg.norm(self._vector)
        try:
            self._direction = self._vector/self._length
        except:
            raise ValueError('Invalid vertex for Arrow: {}'.format(vertex))

    @property
    def vertex(self):
        return self._vertex

    @property
    def length(self):
        return self._length

    @property
    def position(self):
        return self._position

    @property
    def vector(self):
        return self._vector

    @property
    def direction(self):
        return self._direction

    @classmethod
    def from_two_vectors(cls, v1, v2):
        """
        Create a triangle from two vectors
        """
        p0 = np.zeros(3)
        return cls(np.stack((p0, v1, v2)))

    @classmethod
    def from_three_points(cls, p1, p2, p3):
        return cls(np.stack((p1, p2, p3)))

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

class Triangle:
    """
    module for a triangle object.
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
        self._area = 0.5*norm
        self._centroid = vertex[0] + (v1+v2)*(1/3)

    @property
    def vertex(self):
        return self._vertex

    @property
    def area(self):
        return self._area

    @property
    def vector(self):
        return self._vector

    @property
    def normal(self):
        return util.normalized_vector(self.vector)

    @property
    def centroid(self):
        return self._centroid

    @property
    def center(self):
        return self.centroid

    @classmethod
    def from_two_vectors(cls, v1, v2):
        """
        Create a triangle from two vectors
        """
        p0 = np.zeros(3)
        return cls(np.stack((p0, v1, v2)))

    @classmethod
    def from_three_points(cls, p1, p2, p3):
        return cls(np.stack((p1, p2, p3)))

    def copy(self):
        return self.__class__(vertex=self.vertex)

    def __repr__(self):
        s = '<triangle>\n'
        s += ' <vertex="{}"/>\n'.format(self.vertex)
        s += ' area="{}" vector="{}" centroid="{}"\n'.format(
            self.area, self.vector, self.centroid)
        s += '</triangle>'
        return s

    def __str__(self):
        return self.__repr__()

class Tetrahedron:
    """
    module for a tetrahedron object.
    """
    def __init__(self, vertex):
        """
        Create a tetrahedron
        """
        self._vertex = np.array(vertex)
        self._base = Triangle(vertex[0:3])
        self._apex = vertex[3]
        v1 = vertex[1] - vertex[0]
        v2 = vertex[2] - vertex[0]
        v3 = vertex[3] - vertex[0]
        v1v2 = np.cross(v1, v2)
        self._volume = np.dot(v1v2, v3)/6
        self._centroid = vertex[0] + (v1+v2+v3)/4
        dvec = self._apex - self._centroid
        self._direction = dvec/np.linalg.norm(dvec)

    @property
    def vertex(self):
        return self._vertex

    @property
    def base(self):
        return self._base

    @property
    def apex(self):
        return self._apex

    @property
    def centroid(self):
        return self._centroid

    @property
    def center(self):
        return self.centroid

    @property
    def position(self):
        return self.center

    @property
    def vector(self):
        return self.apex - self.center

    @property
    def direction(self):
        return util.normalized_vector(self.vector)

    @classmethod
    def from_three_vectors(cls, v1, v2, v3):
        """
        Create a tetrahedron from three vectors
        """
        p0 = np.zeros(3)
        return cls(np.stack((p0, v1, v2, v3)))

    @classmethod
    def from_base_and_apex(cls, base, apex):
        """
        Create a tetrahedron from a base (triangle) and an apex
        """
        return cls(np.stack((base.vertex, apex)))

    @classmethod
    def from_four_points(cls, p0, p1, p2, p3):
        """
        Create a tetrahedron from 4 points.
        """
        return cls(np.stack((p0, p1, p2, p3)))

    @classmethod
    def from_three_points(cls, p0, p1, p2, scale=1.0):
        """
        Create a tetrahedron from three points.

        The three points form the base and
        the apex is above the centroid of the base 
        at the distance of close-packing

        """
        base = Triangle.from_three_points(p0, p1, p2)
        hhat = base.normal
        a = math.sqrt(4*base.area/math.sqrt(3))
        h = (math.sqrt(6)/3)*a*scale
        apex = base.centroid + h*hhat
        return cls.from_four_points(p0, p1, p2, apex)

    @classmethod
    def from_two_points(cls, p0, p1, scale=1.0):
        """
        Create a tetrahedron from two points.

        The first point is at the centroid of the base and
        the second point at the apex of the tetrahedron.

        """
        orig = p0
        hvec = p1 - p0
        h = np.linalg.norm(hvec)*scale
        a = (3/math.sqrt(6))*h
        hhat = util.normalized_vector(hvec)
        ghat = util.perp_direction(hhat)
        gvec = (a/math.sqrt(3))*ghat
        dvec = (1/2)*gvec
        fvec = (a/2)*np.cross(ghat, hhat)
        v0 = orig + gvec
        v1 = orig - fvec - dvec
        v2 = orig + fvec - dvec
        v3 = orig + hvec
        return cls(np.array([v0, v1, v2, v3]))

    @classmethod
    def from_one_point(cls, point, direction=np.array([0, 0, 1]), scale=1.0):
        """
        Create a tetrahedron associated with a point
        The given point is at the centroid of the base and
        the vector is from the point to the apex.

        """
        hhat = util.normalized_vector(direction)
        a = scale
        h = (math.sqrt(6)/3)*a
        hvec = h*hhat
        ghat = util.perp_direction(hhat)
        gvec = (a/math.sqrt(3))*ghat
        dvec = (1/2)*gvec
        fvec = (a/2)*np.cross(ghat, hhat)
        orig = point
        v0 = orig + gvec
        v1 = orig - fvec - dvec
        v2 = orig + fvec - dvec
        v3 = orig + hvec
        return cls(np.array([v0, v1, v2, v3]))

    def copy(self):
        return self.__class__(vertex=self.vertex)

    def scale(self, factor):
        center = self.center
        vertex = []
        for v0 in self.vertex:
            v1 = center + (v0 - center)*factor
            vertex.append(v1)
        return self.__class__(vertex=np.array(vertex))

    def move_by(self, disp):
        vertex = []
        for v0 in self.vertex:
            vertex.append(v0 + disp)
        return self.__class__(vertex=np.array(vertex))

    @property
    def volume(self):
        return self._volume

    @property
    def pointList(self):
        plist=[]
        for i in range(4):
            plist.append(self.vertex[i])
        return plist

    def __repr__(self):
        s = '<tetrahedron>\n'
        s += ' <vertex="{}"/>\n'.format(self.vertex)
        s += ' volume="{}" centroid="{}" />\n'.format(self.volume, self.centroid)
        s += ' base.area="{}" base.centeroid="{}" \n'.format(self.base.area, self.base.centroid)
        s += ' apex="{}" direction="{}"\n'.format(self.apex, self.direction)
        s += '</tetrahedron>'
        return s

    def __str__(self):
        return self.__repr__()

if __name__ == "__main__":

    # stand-alone mode
    orig = np.zeros(3)
    hvec = np.array([0, 0, 1])
    h = np.linalg.norm(hvec)
    a = (3/math.sqrt(6))*h
    hhat = util.normalized_vector(hvec)
    ghat = util.perp_direction(hhat)
    gvec = (a/math.sqrt(3))*ghat
    dvec = (1/2)*gvec
    fvec = (a/2)*np.cross(ghat, hhat)
    p0 = orig + gvec
    p1 = orig - fvec - dvec
    p2 = orig + fvec - dvec
    p3 = orig + hvec

    base = Triangle.from_three_points(p0, p1, p2)
    print('base = ', base)

    tet1 = Tetrahedron.from_four_points(p0, p1, p2, p3)
    print('tet1 = ', tet1)
    print('Expected:')
    base_area = (math.sqrt(3)/4)*(a**2)
    vol = (math.sqrt(2)/12)*(a**3)
    centroid = base.centroid + hvec*(1/4)
    print('base area = {}, volume={}'.format(base.area, vol))
    print('centroid = ', centroid)

    print('')
    tet12 = tet1.move_by(-tet1.center)
    print('tet12 after move = ', tet12)

    print('')
    tet13 = tet12.scale(2.0)
    print('tet13 after scale = ', tet13)

    tet2 = Tetrahedron.from_three_points(p0, p1, p2)
    print('')
    print('tet2 from 3 points = ', tet2)
    print('')
    print('Expected:')
    print('volume=', vol)
    print('centroid = ', centroid)

    tet3 = Tetrahedron.from_two_points(base.centroid, p3)
    print('')
    print('tet3 from two points = ', tet3)

    tet4 = Tetrahedron.from_one_point(base.centroid, direction=hvec, scale=1.0)
    print('')
    print('\ntet4 from one point = ', tet4)
