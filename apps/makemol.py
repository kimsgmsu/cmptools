#!/usr/bin/env python
"""
makemol: Make/edit molecular/crystal structure file
"""

import sys
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET
import os
import math
import numpy as np
import spglib
import util, cmptlib, Command

class myapp:

    def __init__(self):

        self.moleculeList = {}
        self.latticeList = {}
        self.selectionList = {}
        self.planeList = {}
        self.pointList = {}
        self.vectorList = {}
        self.anchorList = {}
        self.kpointList = {}
        self.variableList = {}
        self.verbose = 1

    def run(self, args):

        # banner
        print('========================')
        print('         MAKEMOL        ')
        print('========================')

        # set parameters
        self.ctrlfile = args.ctrlfile
        self.verbose = args.verbose
    
        # open control file
        try:
            self.ctrl = open(self.ctrlfile, 'r')
            if self.verbose > 0:
                print('%makemol: Reading control file "{}"'.format(
                    self.ctrlfile))
        except IOError:
            util.Error('%makemol: Failed to open control file "{}"'.format(
                self.ctrlfile))

        # get the element tree
        self.tree = ET.parse(self.ctrl)
        self.ctrl.close()
        if self.verbose > 4:
            for node in tree.iter():
                print(node.tag, node.attrib)

        # root
        root = self.tree.getroot()
        print("root=", root.tag)
        if root.tag != "makemol":
            print("%makemol: Invalid control file for makemol.")
            util.Error('%makemol: The root tag must be "makemol"')

        # process each command
        for comm in root:
            stat = self.proc_cmd(comm)
            if not stat:
                break

        # done
        print("%makemol: END")

    def proc_cmd(self, comm):
        """
        Process one command
        """
        command = comm.tag
        if self.verbose > 3:
            print('command = {}'.format(command))

        if command == "makemol":
            pass

        elif command == "stop":
            return False

        elif command == "addAtom":
            Command.addAtom(self, comm)

        elif command == "anchor":
            Command.anchor(self, comm)

        elif command == "applySYMOP":
            mol = comm.get("mol")
            op = comm.get("op")
            dest = comm.get("dest")
            self.applySYMOP(mol, op, dest)

        elif command == "attachMol":
            Command.attachMol(self, comm)

        elif command == "boxMol":
            Command.boxMol(self, comm)

        elif command == "centerMol":
            Command.centerMol(self, comm)

        elif command == "checkSym":
            symGrp = comm.get("group")
            self.checkSym(molecule, symGrp)

        elif command == "collect_atoms_to_unitcell":
            Command.collect_atoms_to_unitcell(self, comm)

        elif command == "comment":
            pass

        elif command == "copyMol":
            Command.copyMol(self, comm)

        elif command == "deleteAtoms":
            Command.deleteAtoms(self, comm)

        elif command == "diffMol":
            mol0 = comm.get("mol0")
            mol1 = comm.get("mol1")
            self.diffMol(mol0, mol1)

        elif command == "drawSpheres":
            Command.drawSpheres(self, comm)

        elif command == "editAtoms":
            Command.editAtoms(self, comm)

        elif command == "gen_equiv_atoms":
            self.gen_equiv_atoms(comm)

        elif command == "generateSlab":
            Command.generateSlab(self, comm)

        elif command == "interpolate_height":
            Command.interpolate_height(self, comm)

        elif command == "kpoint":
            Command.kpoint(self, comm)

        elif command == "lattice":
            Command.lattice(self, comm)

        elif command == "makeSupercell":
            Command.makeSupercell(self, comm)

        elif command == "makeTransitionalMol":
            Command.makeTransitionalMol(self, comm)

        elif command == "mergeMol":
            Command.mergeMol(self, comm)

        elif command == "molecule":
            Command.molecule(self, comm)

        elif command == "move_atoms_to_unitcell":
            Command.move_atoms_to_unitcell(self, comm)

        elif command == "moveAtoms":
            Command.moveAtoms(self, comm)

        elif command == "moveMol":
            Command.moveMol(self, comm)

        elif command == "newMol":
            Command.newMol(self, comm)

        elif command == "padMol":
            Command.padMol(self, comm)

        elif command == "plane":
            Command.plane(self, comm)

        elif command == "point":
            Command.point(self, comm)

        elif command == "printDisplacement":
            Command.printDisplacement(self, comm)

        elif command == "printNeighborMap":
            Command.printNeighborMap(self, comm)

        elif command == "randomizePos":
            Command.randomizePos(self, comm)

        elif command == "readMol":
            Command.readMol(self, comm)

        elif command == "report":
            Command.report(self, comm)

        elif command == "reportSymmetry":
            Command.reportSymmetry(self, comm)

        elif command == "reduceAtomPos":
            self.reduceAtomPos(comm)

        elif command == "remove_overlap_atoms":
            Command.remove_overlap_atoms(self, comm)

        elif command == "replicateMol":
            Command.replicateMol(self, comm)

        elif command == "reproduceMol":
            Command.reproduceMol(self, comm)

        elif command == "rotateMol":
            self.rotateMol(comm)

        elif command == "select":
            Command.selection(self, comm)

        elif command == "set":
            Command.set_param(self, comm)

        elif command == "setMol":
            Command.setMol(self, comm)

        elif command == "setVariable":
            Command.setVariable(self, comm)

        elif command == "standardizeLattice":
            Command.standardizeLattice(self, comm)

        elif command == "updateMol":
            Command.updateMol(self, comm)

        elif command == "vector":
            Command.vector(self, comm)

        elif command == "wrapMol":
            Command.wrapMol(self, comm)

        elif command == "writeMol":
            Command.writeMol(self, comm)

        else:
            util.Error('%makemol: Unrecognized command "{}"'.format(command))

        return True

    def diffMol(self, mid0, mid1):
        # difference between two molecules
        mol0 = self.getMol(mid0)
        if mol0 is None:
            return False
        mol1 = self.getMol(mid1)
        if mol1 is None:
            return False
        return cmptlib.diffMol(mol0, mol1)

    def gen_equiv_atoms(self, comm):
        mid = comm.get("mol")
        if mid is None:
            raise EX.inputError("Missing input molecule name.")
        mol = self.moleculeList.get(mid)
        if mol is None:
            raise EX.inputError('Invalid molecule "{}"'.format(mid))
        spcgrp = comm.get("spacegroup")
        if spcgrp is None:
            raise EX.inputError("Missing space group name.")
        cmptlib.gen_equiv_atoms(mol, spcgrp)

    def reduceAtomPos(self, comm):
        mol = util.get_molecule("mol", comm, self)
        cmptlib.reduceAtomPos(mol)

    def sanitizeMol(self, comm):
        mol = util.get_molecule("mol", comm, self)
        mol.sanitize()

    def get_kpoint(self, key):
        kpoint = self.kpointList.get(key)
        if kpoint is None:
            raise Exception('Invalid kpoint name "{}"'.format(key))
        return kpoint

    def get_lattice(self, key):
        latt = self.latticeList.get(key)
        if latt is None:
            raise Exception('Invalid lattice name "{}"'.format(key))
        return latt

    def get_molecule(self, key):
        mol = self.moleculeList.get(key)
        if mol is None:
            raise Exception('Invalid molecule name "{}"'.format(key))
        return mol

    def get_plane(self, key):
        plane = self.planeList.get(key)
        if plane is None:
            raise Exception('Invalid plane name "{}"'.format(key))
        return plane

    def get_point(self, key):
        pos = self.pointList.get(key)
        if pos is None:
            raise Exception('Invalid point name "{}"'.format(key))
        return pos

    def get_selection(self, key):
        if key is None:
            return None
        sel = self.selectionList.get(key)
        if sel is None:
            raise Exception('Invalid selection name "{}"'.format(key))
        return sel

    def get_anchor(self, key):
        anchor = self.anchorList.get(key)
        if anchor is None:
            raise Exception('Invalid anchor name "{}"'.format(key))
        return anchor

    def get_vector(self, key):
        vec = self.vectorList.get(key)
        if vec is None:
            raise Exception('Invalid vector name "{}"'.format(key))
        return vec

    def get_variable(self, key):
        val = self.variableList.get(key)
        if val is None:
            raise Exception('Invalid variable name "{}"'.format(key))
        return (key,val)

if __name__ == "__main__":
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='Edit a molecular structure file')
    parser.add_argument('ctrlfile', help='control file')
    parser.add_argument('--verbose', type=int, default=1,
                        help='Verbose level index (0-5)')

    args = parser.parse_args()

    app = myapp()
    app.run(args)
