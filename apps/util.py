#!/usr/bin/env python

"""
util: convenience utility module
"""

import sys, os
import math, re
import numpy as np
import numpy.linalg as la
import Plane

# tolerance for zero
ZTOL = 1.0e-9

def Error(msg):
    print(msg)
    sys.exit(0)

def bool_to_str(v):
    """
    Convert a boolean to a string
    """
    if v:
        return 'T'
    else:
        return 'F'

def bool_list_to_str(vlist):
    """
    Convert a list of booleans to a string
    """
    try:
        return ', '.join([bool_to_str(v) for v in vlist])
    except:
        Error('Invalid boolean list "{}"'.format(vlist))

def list_to_str(vlist, fmt=None):
    """
    Convert a list to a string (eg, "[1.0, 3, -5.1]")
    """
    try:
        if fmt == None:
            s = ', '.join([str(v) for v in vlist])
        else:
            s = ', '.join([fmt.format(v) for v in vlist])
        return '[' + re.sub('\s+', ' ', s).strip() + ']'
    except:
        Error('Invalid list "{}"'.format(vlist))

def sanitize_list_string(st):
    """
    Add missing "[", "]", and "," to the given string
    """
    print('string="{}"'.format(st))
    s2 = re.sub('\s+', ' ', st).strip()
    if s2[0] != '[':
        s2 = '[' + s2
    if s2[-1] != ']':
        s2 = s2 + ']'
    if ',' not in s2:
        s2 = s2.replace(' ', ', ')
    return s2

def cleanup_list_string(s):
    """
    Replace "[", "]", and "," to " " 
    """
    return s.replace('[', '').replace(']', '').replace(',', ' ').strip()

def str_to_bool(s):
    """
    Convert a string to boolean value
    """
    if str(s).lower() in ("yes", "true", "t", "1"):
        v = True
    elif str(s).lower() in ("no", "false", "f", "0"):
        v = False
    else:
        raise Exception('Invalid boolean value "{}"'.format(s))
    return v

def str_to_bool_list(s):
    """
    Convert a string to a boolean list
    """
    s2 = cleanup_list_string(s)
    try:
        return [str_to_bool(w) for w in s2.split()]
    except:
        Error('Invalid boolean list "{}"'.format(s))
        
def str_to_list(s, VARS=None):
    """
    Convert a string (eg, "[1.0, 3, -5.1]") to a list
    """
    s2 = cleanup_list_string(s)
    try:
        return [eval(w, {}, VARS) for w in s2.split()]
    except:
        Error('Invalid number in "{}"'.format(s))
        
def str_to_name_list(s):
    """
    Convert a string (eg, "[mol1, mol23, point3]") to a list of names
    """
    s2 = cleanup_list_string(s)
    try:
        return [w for w in s2.split()]
    except:
        Error('Invalid list of names "{}"'.format(s))
        
def abs_cap(val, max_abs_val=1):
    """
    Returns the value with its absolute value capped at max_abs_val.
    Particularly useful in passing values to trignometric functions where
    numerical errors may result in an argument > 1 being passed in.

    Args:
        val (float): Input value.
        max_abs_val (float): The maximum absolute value for val. Defaults to 1.

    Returns:
        val if abs(val) < 1 else sign of val * max_abs_val.
    """
    return max(min(val, max_abs_val), -max_abs_val)

def reduced_vector(vec, pbc=True):
    """
    Reduced vector with components in the range (-0.5, 0.5]
    """
    return np.array([reduced_number(v1, pbc) for v1 in np.nditer(vec)]).reshape(vec.shape)

def reduced_number(v, pbc=True):
    """
    Reduced number in the range (-0.5, 0.5]
    """
    if pbc:
        v0 = v - math.ceil(v-0.5)
        if v0 < -0.5:
            v0 = 0.5
        elif math.isclose(v0, -0.5, rel_tol=1.0e-6):
            v0 = 0.5
        return v0
    else:
        return v1

def reduced_vector2(vec):
    """
    Reduced vector with components in the range [0, 1)
    """
    return np.array([reduced_number2(v1) for v1 in np.nditer(vec)]).reshape(vec.shape)

def reduced_number2(v):
    """
    Reduced number in the range [0.0, 1.0)
    """
    v2 = v - math.floor(v)
    if v2 > 1.0:
        v2 = 0.0
    elif math.isclose(v2, 1.0, rel_tol=1.0e-6):
        v2 = 0.0
    return v2

def reduced_index(idx, ndim):
    """
    Reduced index with components in the range [0, ndim[i])
    """
    return [i % ndim[k] for (k,i) in enumerate(idx)]

def normalized_vector(v):
    vn = np.linalg.norm(np.array(v))
    try:
        return v / vn
    except:
        return v
    
def perp_direction(v):
    """
    Get a vector most perpendicular to a given vector
    """
    v_cross_u = np.zeros((3,3))
    vcu_norm = np.zeros(3)
    u = np.eye(3)
    for idir in range(3):
        v_cross_u[idir] = np.cross(v, u[idir])
        vcu_norm[idir] = np.linalg.norm(v_cross_u[idir])

    maxdir = np.argmax(vcu_norm)
    return normalized_vector(v_cross_u[maxdir])

def angle_between(v1, v2):
    """ 
    Returns the angle in radians between vectors 'v1' and 'v2'
    """
    arg1 = la.norm(np.cross(v1, v2))
    arg2 = np.dot(v1, v2)
    angle = np.arctan2(arg1, arg2)
    return angle

def cell_volume(avec):
    """
    Volume of a cell defined by three vectors
    """
    v01 = np.cross(avec[0], avec[1])
    v012 = np.dot(v01, avec[2])
    return v012

def get_cell_range_from_box(box):
    abox = np.array(box)
    frange = [(1000,-1000), (1000,-1000), (1000,-1000)]
    # corners
    frange = get_fractional_range(abox[0], frange)
    frange = get_fractional_range(abox[1], frange)
    frange = get_fractional_range(abox[2], frange)
    # diagonals
    frange = get_fractional_range(abox[0]+abox[1], frange)
    frange = get_fractional_range(abox[1]+abox[2], frange)
    frange = get_fractional_range(abox[2]+abox[0], frange)
    frange = get_fractional_range(abox[0]+abox[1]+abox[2], frange)
    crange = []
    for idir in range(3):
        cr1 = (math.floor(frange[idir][0]), math.ceil(frange[idir][1]))
        crange.append(cr1)
    return crange

def get_fractional_range(vec, frange):
    origin = np.zeros(3)
    for idir in range(3):
        fr1 = frange[idir]
        p0 = min(origin[idir], vec[idir], fr1[0])
        p1 = max(origin[idir], vec[idir], fr1[1])
        frange[idir] = (p0, p1)
    return frange

def get_bounding_box(center, radius, mol, ndim):
    bbox = []
    for idir in range(3):
        rng1 = get_bounding_range(idir, center, radius, mol, ndim)
        bbox.append(rng1)
    return bbox

def get_bounding_range(kdir, center, radius, mol, ndim):
    miller = [ 0, 0, 0]
    miller[kdir] = 1
    plane = Plane.Plane.from_miller_and_pos(miller, mol.lattice.matrix, center)
    pos1 = center - plane.normal*radius
    pos2 = center + plane.normal*radius
    fpos0 = mol.to_fractional(center)
    fpos1 = mol.to_fractional(pos1)
    fpos2 = mol.to_fractional(pos2)
    ip0 = math.floor(fpos0[kdir]*ndim[kdir])
    ip1 = math.floor(fpos1[kdir]*ndim[kdir])
    ip2 = math.ceil(fpos2[kdir]*ndim[kdir])
    kp1 = min(ip1, ip2) - 1
    kp2 = max(ip1, ip2) + 1
    print('bbox: center={}, radius={}, normal={}, pos1={}, pos2={}'.format(
        center, radius, plane.normal, pos1, pos2))
    print('bbox: fpos0={}, fpos1={}, fpos2={}'.format(fpos0, fpos1, fpos2))
    print('bbox: ip0={}, ip1={}, ip2={}, (kp1,kp2)=({},{})'.format(ip0, ip1, ip2, kp1, kp2))
    return (kp1, kp2)

def degree(rad):
    return rad*180.0/math.pi

def radian(deg):
    return deg*math.pi/180.0

def asSpherical(xyz):
    #takes list xyz (single coord)
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    xy2 = x*x + y*y
    r = math.sqrt(xy2 + z*z)
    theta = math.atan2(z, math.sqrt(xy2))*180/math.pi #to degrees
    phi = math.atan2(y,x)*180/math.pi
    return [r,theta,phi]

def asCartesian(rthetaphi):
    #takes list rthetaphi (single coord)
    r = rthetaphi[0]
    theta = rthetaphi[1]*math.pi/180 # to radian
    phi = rthetaphi[2]*math.pi/180
    x = r * math.sin( theta ) * math.cos( phi )
    y = r * math.sin( theta ) * math.sin( phi )
    z = r * math.cos( theta )
    return [x,y,z]

def latt_leng_angle(avec):
    vlen = np.array([0.0, 0.0, 0.0])
    vang = np.array([0.0, 0.0, 0.0])
    for idir in range(3):
        jdir = (idir+1) % 3
        kdir = (jdir+1) % 3
        vlen[idir] = la.norm(avec[idir])
        vang[idir] = angle_between(avec[jdir], avec[kdir])
    return vlen, vang

def progress(count, total, status=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    sys.stdout.flush()

def get_next_nb_line(f):
    while True:
        line = f.readline()
        if len(line) < 1:
            return line
        line = line.strip()
        if len(line) > 0:
            return line

def end_of_file(f):
    pos = f.tell()
    line = get_next_nb_line(f)
    if len(line) < 1:
        ans = True
    else:
        ans = False
    f.seek(pos)
    return ans

def skip_blank_lines(f):
    while True:
        pos = f.tell()
        line = f.readline()
        if len(line) < 1:
            break
        line = line.strip()
        if len(line) > 0:
            break
    f.seek(pos)
    return

def get_string(key, node, req=True, default=None):

    s = node.get(key, default=default)
    if req and s is None:
        raise Exception('Missing required key "{}"'.format(key))
    return s

def get_choice(key, choiceList, node, req=True, default=None):

    s = get_string(key, node, req=req, default=default)
    if s in choiceList:
        return s
    else:
        raise Exception('Invalid choice "{}" is not in {}'.format(s, choiceList))

def get_bool(key, node, req=True, default=None):

    s = get_string(key, node, req=req, default=default)
    return str_to_bool(s)

def get_int(key, node, req=True, default=None):

    s = get_string(key, node, req=req, default=default)
    try:
        v = int(s)
    except IOError:
        raise Exception('Invalid int value "{}"'.format(s))
    return v

def get_float(key, node, req=True, default=None):

    s = get_string(key, node, req=req, default=default)
    if s is None:
        return s
    try:
        v = float(s)
    except IOError:
        raise Exception('Invalid float value "{}"'.format(s))
    return v

def get_value(key, node, system, req=True, default=None):

    s = get_string(key, node, req=req, default=default)
    if s is None:
        if req:
            raise Exception('Missing required key "{}"'.format(key))
        else:
            return default
    try:
        return eval(s, {}, system.variableList)
    except:
        Error('Invalid value in "{}"'.format(s))

def get_value_array(key, node, req=True, default=None):

    s = get_string(key, node, req=req, default=None)
    if s is None:
        return default
    try:
        return np.array(eval(s))
    except:
        Error('Invalid value in "{}"'.format(s))

def get_value_list(key, node, system=None, req=True, default=None):

    s = get_string(key, node, req=req, default=None)
    if s is None:
        return default
    s2 = cleanup_list_string(s)
    try:
        return [eval(w) for w in s2.split()]
    except:
        Error('Invalid value in "{}"'.format(s))

def get_name_list(key, node, req=True, default=None):

    s = get_string(key, node, req=req, default=None)
    if s is None:
        return default
    s2 = cleanup_list_string(s)
    try:
        return [w for w in s2.split()]
    except:
        Error('Invalid list of names "{}"'.format(s))

def get_bool_list(key, node, req=True, default=None):

    s = get_string(key, node, req=req, default=default)
    if s is None:
        return None
    try:
        v = str_to_bool_list(s)
    except IOError:
        Error('Invalid boolean list "{}"'.format(s))
    return v

def get_anchor(key, node, system):

    oid = get_string(key, node, req=True)
    return system.get_anchor(oid)

def get_kpoint(key, node, system):

    kid = get_string(key, node, req=True)
    return system.get_kpoint(kid)

def get_lattice(key, node, system, req=True, default=None):

    lid = get_string(key, node, req=req, default=default)
    if lid is None:
        return default
    return system.get_lattice(lid)

def get_molecule(key, node, system, req=True, default=None):

    mid = get_string(key, node, req=req, default=default)
    if mid is None:
        return default
    return system.get_molecule(mid)

def get_atom(key, node, mol):

    aid = get_int(key, node, req=True)
    try:
        return mol.atom_of_index(aid)
    except:
        raise Exception('Invalid atom index #{}'.format(aid))

def get_plane(key, node, system):

    pid = get_string(key, node, req=True)
    return system.get_plane(pid)

def get_point(key, node, system):

    pid = get_string(key, node, req=True)
    return system.get_point(pid)

def get_selection(key, node, system,
                  req=True, default=None, return_sid=False):

    sid = get_string(key, node, req=req, default=None)
    if sid is None:
        sel = default
    else:
        sel = system.get_selection(sid)

    if return_sid:
        return sel, sid
    else:
        return sel

def get_variable(key, node, system):

    vid = get_string(key, node, req=True)
    return system.get_variable(vid)

def get_vector(key, node, system):

    vid = get_string(key, node, req=True)
    return system.get_vector(vid)

def within_lattice(fpos, tol=1.0e-4):
    """
    Determine if the given atom position is inside a lattice
    """
    inside = True
    for idir in range(3):
        if fpos[idir] < -tol:
            inside = False
            break
        if fpos[idir] > 1.0-tol:
            inside = False
            break
        
    return inside

def read_next_line(f):
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue  # Discard comments and empty lines
        return line

def peek_line(f):
    """
    Read ahead one line without advancing the file pointer
    """
    pos = f.tell()
    line = f.readline()
    f.seek(pos)
    return line

if __name__ == "__main__":
    # stand-alone mode
    a = [1.0, 3, -5.1]
    s = list_to_str(a)
    print('{} ==> {}'.format(a,list_to_str(a)))
    print('{} ==> {}'.format(a,list_to_str(a, fmt='{:>8.2f}')))

    s = '[1.0,2, -3.0]'
    print('cleanup {} ==> {}'.format(s,cleanup_list_string(s)))

    b = np.array([[-0.5, 0.2, 0.7, 2.0, 0.9999999],
                  [-0.7, 2.3, -1.6, -1.0, 1.0000001]])
    print('original vector:', b)
    print('reduced vector:', reduced_vector(b))
    print('reduced vector2:', reduced_vector2(b))

    hvec = np.array([1, 0, 0])
    ghat = perp_direction(hvec)
    print('hvec = ', hvec)
    print('perp_dir = ', ghat)

    hvec = np.array([0, 0, 1])
    ghat = perp_direction(hvec)
    print('hvec = ', hvec)
    print('perp_dir = ', ghat)

    hvec = np.array([0, 1, 1])
    ghat = perp_direction(hvec)
    print('hvec = ', hvec)
    print('perp_dir = ', ghat)
    
    hvec = np.array([0, 1, 0.0001])
    ghat = perp_direction(hvec)
    print('hvec = ', hvec)
    print('perp_dir = ', ghat)
    
    v1 = [1.0, 0.0, 0.0]
    v2 = [0.1, -2.0, 0.0]
    print('v1 = {}'.format(v1))
    print('v2 = {}'.format(v2))
    print('angle(v1,v2) = {} rad'.format(angle_between(v1,v2)))
    print('angle(v1,v2) = {} deg'.format(degree(angle_between(v1,v2))))

    v1 = [1.0, 0.0, 0.0]
    v2 = [1.2, 0.0, -0.01]
    print('v1 = {}'.format(v1))
    print('v2 = {}'.format(v2))
    m1 = la.norm(v1)
    m2 = la.norm(v2)
    print('norm(v1), norm(v2) = {}, {}'.format(m1, m2))
    ang12 = angle_between(v1,v2)
    print('angle(v1,v2) = {} rad'.format(ang12))
    print('angle(v1,v2) = {} deg'.format(degree(ang12)))
    print('dot(v1,v2) = {}'.format(np.dot(v1,v2)))
    print('|v1|.|v2|.cos(ang(v1,v2)) = {}'.format(m1*m2*np.cos(ang12)))
