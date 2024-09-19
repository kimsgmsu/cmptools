#!/usr/bin/env python

import util
import cmptlib
import numpy as np

def process(system, node):
    """
    Process a select command
    """
    action = node.get('action')

    if action == 'new':
        new(system, node)
    elif action == 'add':
        add(system, node)
    elif action == 'difference':
        difference(system, node)
    elif action == 'intersection':
        intersection(system, node)
    elif action == 'report':
        report(system, node)
    elif action == 'union':
        union(system, node)
    else:
        raise Exception('%Selection: Undefined action "{}"'.format(action))
        
def add(system, node):

    atoms = get_atomIndex_list(node, system)
    sel = util.get_selection('sel', node, system)
    sel.update(atoms)

def difference(system, node):

    first = util.get_selection('first', node, system)
    second = util.get_selection('second', node, system)
    sel = first.difference(second)
    sid = util.get_string('result', node)
    system.selectionList.update({sid : sel})

def intersection(system, node):

    setnames = util.get_name_list('sets', node, system)
    sel = None
    for sid in setnames:
        s1 = system.get_selection(sid)
        if sel is None:
            sel = s1
        else:
            sel = sel.intersection(s1)

    sid = util.get_string('result', node)
    system.selectionList.update({sid:sel})

def new(system, node):

    mol = util.get_molecule('mol', node, system)
    method = util.get_string('method', node)
    if method == 'atom':
        slist = util.get_string('atom', node)
        if slist == 'All' or slist == 'ALL':
            aidset = set(range(1,len(mol.atomList)+1))
        else:
            try:
                aidset = set(util.str_to_list(slist))
            except IOError:
                raise Exception('%Selection: Invalid atoms list "{}"'.format(slist))
    elif method == 'distance':
        cutoff = util.get_value('cutoff', node, system)
        aidset = mol.overlapped_atoms(cutoff)
    elif method == 'plane':
        plane = util.get_plane('plane', node, system)
        side = util.get_choice('side', ['up', 'down'], node)
        aidset = set()
        for ka, atom in enumerate(mol.atomList):
            dist = plane.distance_from_point(atom.pos)
            if side == "up":
                if dist > 0:
                    aidset.add(ka+1)
            elif side == "down":
                if dist < 0:
                    aidset.add(ka+1)
    elif method == 'symbol':
        symbol = util.get_string('symbol', node)
        aidset = set()
        for ka, atom in enumerate(mol.atomList):
            aid = ka+1
            if atom.symbol == symbol:
                aidset.add(aid)
    else:
        raise Exception('%Set: Undefined method "{}"'.format(method))

    sid = util.get_string('result', node)
    sel = set(aidset)
    system.selectionList.update({sid : sel})

def report(system, node):

    sel, sid = util.get_selection('sel', node, system, return_sid=True)
    header = util.get_string('header', node, req=False, default=sid)
    print('Selection {}: {}'.format(header, sel))

def union(system, node):

    setnames = util.get_name_list('sets', node, system)
    sel = set()
    for sid in setnames:
        sel = sel.union(system.get_selection(sid))

    sid = util.get_string('result', node)
    system.selectionList.update({sid:sel})

if __name__ == '__main__':
    # stand-alone mode
    sel1 = Selection('sel1')
    sel1.add({1, 3, 5, 14})
    sel1.report()
    sel1.add({2, 3, 4, 5, 14})
    sel1.report()
    sel1.add({1, 2, 3, 6})
    sel1.report()

    sel2 = Selection('sel2')
    sel2.add({2, 3, 7, 8, 12, 14, 15})
    sel2.report()
    sel2.add({3, 5, 6, 8})
    sel2.report()

    sel3 = sel1.intersection(sel2)
    sel3.sid = 'sel1.intersection.sel2'
    sel3.report()

    sel4 = sel1.union(sel2)
    sel4.sid = 'sel1.union.sel2'
    sel4.report()

    sel5 = sel1.difference(sel2)
    sel5.sid = 'sel1.difference.sel2'
    sel5.report()
