#!/usr/bin/env python

import math
import numpy as np
import copy
import util, cmptlib
import Lattice, Atom
import pymatgen
from pymatgen.core.periodic_table import Element as pgElement

class VolumetricData:
    """
    class for a volumetric data
    """
    def __init__(self, name, data, extra=''):
        self.name = name
        self.data = data
        self.extra = extra

    def report(self):
        print('VolumeData: {} ({})'.format(self.name, self.data.shape))

class Element:
    """
    module for an element.
    """
    def __init__(self, symbol, info='', mass=None, radius=None):
        """
        Create an element.
        """
        # Element in the database
        try:
            pge = pgElement(symbol)
        except:
            pge = pgElement('H')
        self.symbol = symbol
        self.info = info
        if mass is None:
            s = pge.atomic_mass.__str__()
            self.mass = float(s.split()[0])
        else:
            self.mass = mass
        if radius is None:
            s = pge.atomic_radius.__str__()
            self.radius = float(s.split()[0])
        else:
            self.radius = radius

    def __repr__(self):
        return '<element symbol="{}" radius="{}" mass="{}" info="{}" />'.format(
            self.symbol, self.radius, self.mass, self.info)

    def __str__(self):
        return self.__repr__()

class Molecule:
    """
    module for a molecule/crystal structure.
    """
    def __init__(self, name='', info='', latt=None,
                 elementTable={}, atomList=[], output_coord='frac', volumeDataList=[]):
        """
        Create a molecule.
        """
        self.name = name
        self.info = info
        self.latt = latt
        self.elementTable = elementTable
        self.atomList = atomList
        self.output_coord = output_coord
        self.volumeDataList = volumeDataList

    def apply_operation(self, symmop):
        """
        Apply a symmetry operation to the molecule
        """
        for atom in self.atomList:
            atom.pos = symmop.operate(atom.pos)

    @property
    def natom(self):
        return len(self.atomList)

    @property
    def pbc(self):
        return self.latt.pbc

    def natom_of_symbol(self, symbol):
        """
        No. of atoms of given element symbol
        """
        natom = 0
        for atom in self.atomList:
            if atom.symbol == symbol:
                natom += 1
        return natom

    def make_ElementTable(self):
        """
        Fill in missing elements in the table
        """
        for symbol in self.distinct_atomic_symbols:
            elem = self.elementTable.get(symbol)
            if elem is None:
                self.add_Element(symbol)

    def add_Element(self, symbol):
        """
        Add a chemical element to the table
        """
        elem = self.elementTable.get(symbol)
        if elem is None:
            elem = Element(symbol)
            self.elementTable.update({symbol : elem})

    def get_Element(self, symbol):
        """
        Get the chemical element of the given atomic symbol
        """
        elem = self.elementTable.get(symbol)
        if elem is None:
            elem = Element(symbol)
        return elem

    def get_Element_id(self, symbol):
        """
        Get the id number of the given atomic symbol
        """
        try:
            return tuple(self.elementTable).index(symbol)
        except:
            return -10

    @property
    def distinct_atomic_symbols(self):
        """
        List of distinct atomic symbols
        """
        alist = []
        for atom in self.atomList:
            if atom.symbol not in alist:
                alist.append(atom.symbol)
        return alist

    @property
    def atomic_symbols(self):
        """
        List of atom symbols
        """
        alist = []
        for atom in self.atomList:
            alist.append(atom.symbol)
        return alist

    def atomList_of_symbol(self, symbol):
        """
        Iterate over atoms of the given symbol
        """
        for atom in self.atomList:
            if atom.symbol == symbol:
                yield atom

    def atom_of_index(self, aid):
        """
        Get the atom of the given index
        """
        try:
            return self.atomList[aid-1]
        except:
            raise Exception('invalid atom index: {}'.format(aid))

    def atom_position_in_cell(self, aid, cell):
        """
        position of atom in the given periodic cell
        """
        cpos = self.atom_of_index(aid).pos
        fpos = self.to_fractional(cpos) + np.array(cell)
        return self.to_cartesian(fpos)

    def delete_atom(self, aid):
        """
        Delete an atom from the list
        """
        try:
            print('deleting atom #{} {}'.format(aid, self.atom_of_index(aid)))
            del self.atomList[aid-1]
        except:
            raise Exception('Invalid atom id #{}'.format(aid))

    def selected_atoms(self, selection):
        """
        Iterate over atoms in the selection
        """
        for ka, atom in enumerate(self.atomList):
            if selection is None:
                yield atom
            elif ka+1 in selection:
                yield atom

    def unselected_atoms(self, selection):
        """
        Iterate over atoms not in the selection
        """
        for ka, atom in enumerate(self.atomList):
            if selection is None:
                continue
            elif ka+1 not in selection:
                yield atom

    @property
    def grouped_by_symbol(self):
        """
        Iterate over atoms grouped by symbols
        """
        for symbol in self.distinct_atomic_symbols:
            for atom in self.atomList_of_symbol(symbol):
                yield atom

    def closest_atom(self, pos, coord='cart'):
        """
        Find the closes atom from a given point
        """
        if coord == 'cart':
            pos0 = pos
        else:
            pos0 = self.to_cartesian(pos)

        minaid = None
        mindist = math.inf
        for ka, atom in enumerate(self.atomList):
            dpos = atom.pos - pos0
            dfpos = self.to_fractional(dpos)
            rdfpos = util.reduced_vector(dfpos, self.latt.pbc)
            dpos = self.to_cartesian(rdfpos)
            dist = np.linalg.norm(dpos)
            if dist < mindist:
                minaid = ka+1
                mindist = dist
        return minaid, mindist

    def overlapped_atoms(self, cutoff):
        """
        Find overlapped atoms that are closer than the cutoff distance
        """
        aidset = set()
        for ka1, atom1 in enumerate(self.atomList):
            for ka2, atom2 in enumerate(self.atomList):
                dpos = atom2.pos - atom1.pos
                dfpos = self.to_fractional(dpos)
                rdfpos = util.reduced_vector(dfpos, self.latt.pbc)
                dpos = self.to_cartesian(rdfpos)
                dist = np.linalg.norm(dpos)
                if dist < cutoff:
                    aidset.add(ka2+1)
        return aidset

    @property
    def center_position(self):
        """
        Center position of the molecule
        """
        psum = np.zeros(3)
        for atom in self.atomList:
            psum += atom.pos
        try:
            return psum/self.natom
        except:
            raise Exception('Empty atom list')

    @property
    def center_of_mass(self):
        """
        Center of mass of the molecule
        """
        mpsum = np.zeros(3)
        for atom in self.atomList:
            mass = self.get_Element(atom.symbol).mass
            mpsum += mass*atom.pos
        if self.total_mass > 0:
            return mpsum/self.total_mass
        else:
            return mpsum

    @property
    def total_mass(self):
        """
        Total mass of the structure
        """
        tm = 0
        for atom in self.atomList:
            elem = self.get_Element(atom.symbol)
            tm += elem.mass
        return tm

    @property
    def sdyn_used(self):
        """
        check if sdyn is used for any atom
        """
        return any(atom.sdyn is not None for atom in self.atomList)

    def to_cartesian(self, fpos):
        """
        Get cartesian position from the given fractional position
        """
        return self.latt.to_cartesian(fpos)

    def to_fractional(self, cpos):
        """
        Get fractional coordinates from the given cartesian coordinates
        """
        return self.latt.to_fractional(cpos)

    def cart_to_recip(self, cpos):
        """
        Get reciprocal fractional coordinates from cartesian coordinates
        """
        return self.latt.reciprocal_lattice.to_fractional(cpos)

    def recip_to_cart(self, rpos):
        """
        Get cartesian coord from the reciprocal fractional coord
        """
        return self.latt.reciprocal_lattice.to_cartesian(rpos)

    def move_atoms_to_unitcell(self, overlap='keep', dtol=1.0e-6):
        """
        Move atoms to unit cell
        """
        mol = self.copy()
        mol.atomList = []
        for atom in self.atomList:
            fpos = self.to_fractional(atom.pos)
            pos2 = self.to_cartesian(util.reduced_vector2(fpos))

            minaid, mindist = mol.closest_atom(pos2, coord='cart')

            if mindist < dtol:
                if overlap == 'discard':
                    continue
                elif overlap == 'keep':
                    pass
                else:
                    util.Error('%move_atom_to_unitcell: Invalid overlap option, "{}"'.format(overlap))
            atom2 = atom.copy()
            atom2.pos = pos2
            mol.atomList.append(atom2)

    def frac_range(self):
        """
        ranges of fractional coordinates
        """
        frange = np.zeros((3,2))
        if len(self.atomList) > 0:
            for idir in range(3):
                frange[idir] = [math.inf, -math.inf]
            for atom in self.atomList:
                fpos = self.to_fractional(atom.pos)
                for idir in range(3):
                    if fpos[idir] < frange[idir,0]:
                        frange[idir,0] = fpos[idir]
                    if fpos[idir] > frange[idir,1]:
                        frange[idir,1] = fpos[idir]
        return frange

    def pos_range(self):
        """
        ranges of atom position of the molecule
        """
        prange = np.zeros((3,2))
        if len(self.atomList) > 0:
            for idir in range(3):
                prange[idir] = [math.inf, -math.inf]
            for atom in self.atomList:
                for idir in range(3):
                    if atom.pos[idir] < prange[idir,0]:
                        prange[idir,0] = atom.pos[idir]
                    if atom.pos[idir] > prange[idir,1]:
                        prange[idir,1] = atom.pos[idir]
        return prange

    def make_bounding_box(self, padding=None):
        """
        Make a bounding box for the molecule
        """
        if padding is None:
            padding = np.zeros(3,2)
        else:
            padding = np.array(padding)

        # position range
        prange = self.pos_range()

        # move atoms
        dvec = padding[:,0] - prange[:,0]
        for atom in self.atomList:
            atom.pos = atom.pos + dvec

        # size of the bounding box
        boxsz = np.zeros(3)
        for idir in range(3):
            boxsz[idir] = (prange[idir,1]-prange[idir,0])
            boxsz[idir] += (padding[idir,0]+padding[idir,1])

        # new lattice vectors
        matrix = np.eye(3)
        for idir in range(3):
            matrix[idir,idir] = boxsz[idir]
        self.latt = Lattice.Lattice(matrix)

    def split_voldata(self):
        """
        split the molecule for each volumetric data
        """
        mlist = []
        m0 = self.copy()
        vdlist = self.volumeDataList
        m0.volumeDataList = []
        for vd in vdlist:
            print('splitting {} {}'.format(vd.name, vd.data.shape))
            m = m0.copy()
            m.volumeDataList = [vd]
            mlist.append(m)

        if len(mlist) < 1:
            mlist.append(m0)

        return mlist

    def __repr__(self):
        return '<molecule name="{}" natom={} />'.format(
            self.name, self.natom)

    def report(self, header=None):

        if header is not None:
            print(header)
        print('Name: ', self.name)
        print('Info: ', self.info)

        self.latt.report()

        print('*** Element List ({})'.format(len(self.distinct_atomic_symbols)))
        lbar = 60
        print('-'*lbar)
        print('SYMB (MULT   MASS   RADIUS  MASS%   VOL%) INFO ')
        print('-'*lbar)
        total_mass = self.total_mass
        pop_sum = 0
        vol_sum = 0
        vf_sum = 0
        mass_sum = 0
        mf_sum = 0
        for symbol in self.distinct_atomic_symbols:
            elem = self.get_Element(symbol)
            pop = self.natom_of_symbol(symbol)
            pop_sum += pop
            v1 = (4/3)*np.pi*(elem.radius)**3
            vs = v1*pop
            vol_sum += vs
            vf = (vs/self.latt.volume)*100
            vf_sum += vf
            m1 = elem.mass
            ms = m1*pop
            mass_sum += ms
            mf = (ms/total_mass)*100
            mf_sum += mf
            print('{:>4} ({:>4} {:>6.2f} {:>8.3f} {:>6.2f} {:>6.2f})'.format(
                symbol, pop, elem.mass, elem.radius, mf, vf))
        print('-'*lbar)
        print('Total({:>4} {:>6.2f} {:>8} {:>6.2f} {:>6.2f})'.format(
                pop_sum, mass_sum, ' ', mf_sum, vf_sum))
        print()
        print('*** Atom List ({})'.format(len(self.atomList)))
        lbar = 95
        print('-'*lbar)
        sdyn_used = self.sdyn_used
        print('  ID TYPE (           CART-COORDS           )', end='')
        print(' (            FRAC-COORD           )', end='')
        print(' ( SDYN)' if sdyn_used else '', end='')
        print(' INFO  ')
        print('-'*lbar)
        for i, atom in enumerate(self.atomList):
            print('{:>4}'.format(i+1), end='')
            print(' {:>4}'.format(atom.symbol), end='')
            print(' (', end='')
            for c in atom.pos:
                print(' {:>10.6f}'.format(c), end='')
            print(') (', end='')
            fpos = self.to_fractional(atom.pos)
            for c in fpos:
                print(' {:>10.6f}'.format(c), end='')
            if sdyn_used:
                if atom.sdyn is None:
                    s = 'T T T'
                else:
                    s = util.bool_list_to_str(atom.sdyn)
                print(') ({})'.format(s), end='')
            else:
                print(')', end='')
            print(' {}'.format(atom.info))
        print('-'*lbar)
        # volumetric data
        if len(self.volumeDataList) > 0:
            print('*** Volumetric data: ({})'.format(len(self.volumeDataList)))
            lbar = 40
            print('-'*lbar)
            print('      KEY        ( DIMENSION )')
            print('-'*lbar)
            for vd in self.volumeDataList:
                print('{} {}'.format(vd.name, vd.data.shape))
            print('-'*lbar)

    def copy(self):
        cp = self.__class__()
        cp.name = self.name
        cp.info = self.info
        cp.latt = self.latt
        cp.elementTable = copy.deepcopy(self.elementTable)
        cp.atomList = copy.deepcopy(self.atomList)
        cp.output_coord = self.output_coord
        cp.volumeDataList = copy.deepcopy(self.volumeDataList)
        return cp

    def to_pymatgen_structure(self):
        """
        Convert to pymatgen Structure object
        """
        latt = self.latt.matrix
        slist = self.atomic_symbols
        plist = []
        for atom in self.atomList:
            plist.append(atom.pos)
        return pymatgen.core.structure.Structure(latt, slist, plist,
                                                 coords_are_cartesian=True)

def from_pymatgen_structure(mgStr):
    """
    Return a molecule from pymatgen Structure object
    """
    import itertools

    mol = Molecule()
    mol.latt = Lattice.Lattice(mgStr.latt.matrix)
    mol.atomList = []
    for site in mgStr.sites:
        atom = Atom.Atom()
        atom.symbol = ''.join(itertools.takewhile(str.isalpha, str(site.species)))
        atom.pos = site.coords
        mol.atomList.append(atom)

    mol.make_ElementTable()

    return mol

def process(system, node):
    """
    Process a "molecule" command
    """
    action = node.get('action', default=None)

    if action == 'new':
        new(system, node)
    elif action == 'report':
        report(system, node)
    else:
        raise Exception('%Molecule: Undefined action "{}"'.format(action))
        
def new(system, node):

    method = util.get_string('method', node)
    if method == 'read':
        fname = util.get_string("file", node)
        fmt = util.get_string("format", node, req=False)
        option = util.get_string("option", node, req=False)
        mol = cmptlib.readMol(fname, fmt=fmt, option=option,
                              verbose=system.verbose)

    else:
        raise Exception('%Molecule: Undefined method "{}"'.format(method))

    mid = util.get_string("result", node)
    system.moleculeList.update({mid : mol})
    if system.verbose >= 1:
        print('%Molecule: Read "{}"'.format(fname))

def report(system, node):

    mol = util.get_molecule('mol', node, system)
    header = util.get_string('header', node, req=False)
    mol.report(header)

if __name__ == "__main__":

    # stand-alone mode

    latt = Lattice.Lattice([[10.1, 1.2, 1.3], [2.1, 20.2, 2.3], [3.1, 3.2, 30.3]])
    latt.report()

    elemtable = {}
    elem = Element.Element('H', info='Hydrogen', mass=1.008, radius=1.20)
    print(elem)
    elemtable.update({elem.symbol: elem})
    elem = Element.Element('C', info='Carbon', mass=12.011, radius=1.70)
    print(elem)
    elemtable.update({elem.symbol: elem})
    print('elementList=',elemtable)

    atomList = []
    atom = Atom.Atom('H', latt, [1, 0, 0], info='atom1')
    atomList.append(atom)
    atom = Atom.Atom('C', latt, [2.1, 0.5, -1.8], cart=True, info='atom2')
    atomList.append(atom)
    atom = Atom.Atom('C', latt, [0, 1, 0.5], info='atom3')
    atomList.append(atom)
    atom = Atom.Atom('H', latt, [2.1, 0.5, -1.8], cart=False, info='atom4')
    atomList.append(atom)
    print('atomList=',atomList)

    mol = Molecule(name='my crystal', info='test crystal', latt=latt,
                   elementTable=elemtable, atomList=atomList)
    mol.report('== Molecule ==')

    st = mol.get_pymatgen_structure()
    print('pymatgen.Structure=',st)
