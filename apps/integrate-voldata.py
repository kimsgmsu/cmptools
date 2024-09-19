#!/usr/bin/env python
"""
integrate-voldata: Integrate volumetric data
"""

import os
import numpy as np
import cmptlib, util

def integrate_volume_data_list(vdList):
    """
    Integrate a volumetric data list
    """
    if vdList is None:
        return

    for k, vd in enumerate(vdList):
        Ip, Im = integrate_volume_data(vd)
        print('Volume data {}: I_plus = {}, I_minus = {}, I = {}'.format(k, Ip, Im, Ip+Im))

def integrate_volume_data(vd):

    from scipy.interpolate import RegularGridInterpolator

    ndim = vd.data.shape
    nx = ndim[0]
    ny = ndim[1]
    nz = ndim[2]
    xlist = np.linspace(0.0, 1.0, nx)
    ylist = np.linspace(0.0, 1.0, ny)
    zlist = np.linspace(0.0, 1.0, nz)
    my_func = RegularGridInterpolator((xlist, ylist, zlist), vd.data)

    hx = 1.0/np.float_(nx-1)
    hy = 1.0/np.float_(ny-1)
    hz = 1.0/np.float_(nz-1)

    print('*** Integrating volume data on ({},{},{}) grid ***'.format(nx, ny, nz))
    Ip = 0
    Im = 0
    for iz in range(nz):
        z = (iz+0.5)*hz
        for iy in range(ny):
            y = (iy+0.5)*hy
            for ix in range(nx):
                x = (ix+0.5)*hx
                fval = vd.data[ix, iy, iz]
                if fval < 0:
                    Im += fval
                else:
                    Ip += fval
        util.progress(iz, nz-1)

    print()
    sfac = hx*hy*hz
    Ip *= sfac
    Im *= sfac
    
    return Ip, Im

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('=============================')
        print('    INTEGRATE VOLUME DATA    ')
        print('=============================')

        # read molecule
        mol = cmptlib.readMol(args.mol)

        # integrate the volumetric data
        integrate_volume_data_list(mol.volumeDataList)

if __name__ == "__main__":
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='integrate-voldata: Integrate volumetric data')
    parser.add_argument('mol', help='Input molecule')

    args = parser.parse_args()

    app = myapp()
    app.run(args)
