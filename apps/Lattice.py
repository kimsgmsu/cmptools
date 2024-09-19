#!/usr/bin/env python

from __future__ import print_function

import copy
import numpy as np
import pymatgen
import util

class Lattice:
    """
    A lattice object.
    """
    def __init__(self, matrix, pbc=None):
        """
        Create a lattice from any sequence of 9 numbers. Note that the sequence
        is assumed to be read one row at a time. Each row represents one
        lattice vector.

        Args:
            matrix: Sequence of numbers in any form. Examples of acceptable
                input.
                i) An actual numpy array.
                ii) [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
                iii) [1, 0, 0 , 0, 1, 0, 0, 0, 1]
                iv) (1, 0, 0, 0, 1, 0, 0, 0, 1)
                Each row should correspond to a lattice vector.
                E.g., [[10, 0, 0], [20, 10, 0], [0, 0, 30]] specifies a lattice
                with lattice vectors [10, 0, 0], [20, 10, 0] and [0, 0, 30].
            scale: (float) scale factor for output.
        """
        mat = np.array(matrix, dtype=np.float64).reshape((3, 3))
        lengths = np.sqrt(np.sum(mat ** 2, axis=1))
        angles = np.zeros(3)
        for i in range(3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            angles[i] = util.abs_cap(np.dot(mat[j], mat[k]) / (lengths[j] * lengths[k]))

        self._angles = np.arccos(angles) * 180. / np.pi
        self._lengths = lengths
        self._matrix = mat
        self._f2c = np.transpose(self._matrix)
        self._c2f = np.linalg.inv(self._f2c)
        self._metric_tensor = np.dot(self._f2c.T, self._f2c)
        if pbc is None:
            self._pbc = True
        else:
            self._pbc = pbc

    def __format__(self, fmt_spec=''):
        """
        Support format printing. Supported formats are:

        1. "l" for a list format that can be easily copied and pasted, e.g.,
           ".3fl" prints something like
           "[[10.000, 0.000, 0.000], [0.000, 10.000, 0.000], [0.000, 0.000, 10.000]]"
        2. "p" for lattice parameters ".1fp" prints something like
           "{10.0, 10.0, 10.0, 90.0, 90.0, 90.0}"
        3. Default will simply print a 3x3 matrix form. E.g.,
           10.000 0.000 0.000
           0.000 10.000 0.000
           0.000 0.000 10.000
        """
        m = self.matrix.tolist()
        if fmt_spec.endswith("l"):
            fmt = "[[{}, {}, {}], [{}, {}, {}], [{}, {}, {}]]"
            fmt_spec = fmt_spec[:-1]
        elif fmt_spec.endswith("p"):
            fmt = "{{{}, {}, {}, {}, {}, {}}}"
            fmt_spec = fmt_spec[:-1]
            m = self.lengths_and_angles
        else:
            fmt = "{} {} {}\n{} {} {}\n{} {} {}"
        return fmt.format(*[format(c, fmt_spec) for row in m
                            for c in row])

    @classmethod
    def from_parameters(cls, lengths, angles, pbc=None):
        """
        Create a Lattice using unit cell lengths and angles (in degrees).
        """
        mgLatt = pymatgen.core.lattice.Lattice.from_parameters(
            lengths[0], lengths[1], lengths[2],
            angles[0], angles[1], angles[2])
        return cls(mgLatt.matrix, pbc=pbc)

    def copy(self):
        """Deep copy of self."""
        return self.__class__(self.matrix, pbc=self.pbc)

    @property
    def matrix(self):
        """Copy of matrix representing the Lattice"""
        return np.copy(self._matrix)

    @property
    def pbc(self):
        """
        Get the periodic boundary condition flag
        """
        return self._pbc

    @pbc.setter
    def pbc(self, val):
        """
        Set the periodic boundary condition flag
        """
        self._pbc = val

    @property
    def f2c_matrix(self):
        """
        Matrix for fractional to cartesian transformation
        """
        return self._f2c

    @property
    def c2f_matrix(self):
        """
        Matrix for cartesian to fractional transformation
        """
        return self._c2f

    @property
    def metric_tensor(self):
        """
        The metric tensor of the lattice.
        """
        return self._metric_tensor

    def to_cartesian(self, fpos):
        """
        Convert the fractional coords to cartesian coords.

        Args:
            fpos (3x1 array): Fractional coords.

        Returns:
            Cartesian coords
        """
        return np.dot(self._f2c, fpos)

    def to_fractional(self, cpos):
        """
        Convert the cartesian coords to fractional coords.

        Args:
            cpos (3x1 array): Cartesian coords.

        Returns:
            Fractional coordinates.
        """
        return np.dot(self._c2f, cpos)

    def dot_fractional(self, fpos1, fpos2):
        """
        Returns the dot product of two fractional coords via the metric tensor.

        Args:
            fpos1, fpos2 (3x1 array): fractional coords.

        Returns:
            dot product
        """
        return np.dot(fpos1, np.dot(self._metric_tensor, fpos2))

    def distance_between_frac_points(self, fpos1, fpos2, pbc):
        """
        Returns distance between two fractional coordinates
        """
        dfpos = util.reduced_vector(fpos2 - fpos1, pbc)
        return self.dot_fractional(dfpos, dfpos)

    @property
    def angles(self):
        """
        Returns the angles (alpha, beta, gamma) of the lattice.
        """
        return tuple(self._angles)

    @property
    def a(self):
        """
        *a* lattice parameter.
        """
        return self._lengths[0]

    @property
    def b(self):
        """
        *b* lattice parameter.
        """
        return self._lengths[1]

    @property
    def c(self):
        """
        *c* lattice parameter.
        """
        return self._lengths[2]

    @property
    def abc(self):
        """
        Lengths of the lattice vectors, i.e. (a, b, c)
        """
        return tuple(self._lengths)

    @property
    def alpha(self):
        """
        Angle alpha of lattice in degrees.
        """
        return self._angles[0]

    @property
    def beta(self):
        """
        Angle beta of lattice in degrees.
        """
        return self._angles[1]

    @property
    def gamma(self):
        """
        Angle gamma of lattice in degrees.
        """
        return self._angles[2]

    @property
    def volume(self):
        """
        Volume of the unit cell.
        """
        m = self._matrix
        return abs(np.dot(np.cross(m[0], m[1]), m[2]))

    @property
    def lengths_and_angles(self):
        """
        Returns (lattice lengths, lattice angles).
        """
        return tuple(self._lengths), tuple(self._angles)

    @property
    def reciprocal_lattice(self):
        """
        Return the reciprocal lattice. Note that this is the standard
        reciprocal lattice used for solid state physics with a factor of 2 *
        pi. If you are looking for the crystallographic reciprocal lattice,
        use the reciprocal_lattice_crystallographic property.
        The property is lazily generated for efficiency.
        """
        try:
            return self._reciprocal_lattice
        except AttributeError:
            v = np.linalg.inv(self._matrix).T
            self._reciprocal_lattice = Lattice(v * 2 * np.pi)
            return self._reciprocal_lattice

    @property
    def reciprocal_lattice_crystallographic(self):
        """
        Returns the *crystallographic* reciprocal lattice, i.e., no factor of
        2 * pi.
        """
        return Lattice(self.reciprocal_lattice.matrix / (2 * np.pi))

    def report(self, header=None):

        if header is not None:
            print(header)
        print('*** Lattice:')
        ffmt = '{:>10.6f}'
        lhbar = 59
        print('PBC = ', self.pbc)
        print('-'*lhbar)
        print('indx (   avec[x]   avec[y]   avec[z])     length      angle')
        print('-'*lhbar)
        for idir in range(3):
            print('{:>3}: ('.format(idir+1), end='')
            for v in self.matrix[idir]:
                print(ffmt.format(v), end='')
            print(') ', end='')
            print(ffmt.format(self.abc[idir]), end=' ')
            print(ffmt.format(self.angles[idir]))
        print('-'*lhbar)
        print('  Cell volume = {:>10.6f}'.format(self.volume))
        reclatt = self.reciprocal_lattice
        print('*** Reciprocal lattice:')
        print('-'*lhbar)
        print('indx (   bvec[x]   bvec[y]   bvec[z])     length      angle')
        print('-'*lhbar)
        for idir in range(3):
            print('{:>3}: ('.format(idir+1), end='')
            for v in reclatt.matrix[idir]:
                print(ffmt.format(v), end='')
            print(') ', end='')
            print(ffmt.format(reclatt.abc[idir]), end=' ')
            print(ffmt.format(reclatt.angles[idir]))
        print('-'*lhbar)
        print('  Reciprocal cell volume = {:>10.6f}'.format(reclatt.volume))

def equivalent(latt1, latt2):
    """
    Test if two lattices are equivalent
    """
    return np.allclose(latt1.matrix, latt2.matrix, atol=1.0e-6)
        
def process(system, node):
    """
    Process a "lattice" command
    """
    action = node.get('action', default=None)

    if action == 'new':
        new(system, node)
    elif action == 'report':
        lid = util.get_string('lattice', node, system)
        latt = system.get_lattice(lid)
        header = util.get_string('header', node, req=False)
        latt.report(header)
    else:
        raise Exception('%Lattice: Undefined action "{}"'.format(action))
        
def new(system, node):

    method = util.get_string('method', node)
    if method == 'molecule':
        mol = util.get_molecule('mol', node, system)
        latt = mol.latt

    else:
        raise Exception('%Lattice: Undefined method "{}"'.format(method))

    lid = util.get_string('result', node, system)
    system.latticeList.update({lid : latt})

if __name__ == "__main__":
    # stand-alone mode
    latt = Lattice([[1.1, 0.12, 0.13], [0.21, 2.2, 0.23], [0.31, 0.32, 3.3]])
    latt.report()
    print('f2c matrix=',latt.f2c_matrix)
    print('c2f matrix=',latt.c2f_matrix)
    print('metric tensor=',latt.metric_tensor)

    fpos1 = np.array([1, 0, 0])
    print('fpos1=',fpos1)
    cpos1 = latt.to_cpos(fpos1)
    print('cpos1=',cpos1)
    fpos2 = latt.to_fpos(cpos1)
    print('fpos2=',fpos2)

    fpos1 = np.array([1, 0.5, 0])
    print('fpos1=',fpos1)
    cpos1 = latt.to_cpos(fpos1)
    print('cpos1=',cpos1)
    fpos2 = latt.to_fpos(cpos1)
    print('fpos2=',fpos2)

    print('leng using cpos1 = ',np.sum(cpos1*cpos1))
    print('leng using fpos1 = ',latt.dot_fpos(fpos1, fpos1))
    
