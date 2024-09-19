#!/usr/bin/env python
"""
regrid-voldata: regrid volumetric data
"""

import os
import numpy as np
import util, cmptlib, Molecule

def regrid_volumetric_data_list(vdList0, ndim):
    """
    Regrid a volumetric data list
    """
    if vdList0 is None:
        return None

    vdList = []
    for vd0 in vdList0:
        vd = regrid_volumetric_data(vd0, ndim)
        vdList.append(vd)

    return vdList

def regrid_volumetric_data(vd0, ndim):

    from scipy.interpolate import RegularGridInterpolator

    ndim0 = vd0.data.shape
    xlist0 = np.linspace(0.0, 1.0, ndim0[0])
    ylist0 = np.linspace(0.0, 1.0, ndim0[1])
    zlist0 = np.linspace(0.0, 1.0, ndim0[2])
    my_func = RegularGridInterpolator((xlist0, ylist0, zlist0), vd0.data)

    data = np.zeros( (ndim[0], ndim[1], ndim[2]), dtype=np.float_ )

    xfactor=1.0/np.float_(ndim[0]-1)
    yfactor=1.0/np.float_(ndim[1]-1)
    zfactor=1.0/np.float_(ndim[2]-1)

    print('*** Generating data on new grid {} ***'.format(ndim))
    for iz in range(ndim[2]):
        util.progress(iz, ndim[2], 'Done')
        z = iz*zfactor
        for iy in range(ndim[1]):
            y = iy*yfactor
            for ix in range(ndim[0]):
                x = ix*xfactor
                pos = np.array([x,y,z])
                data[ix,iy,iz] = my_func(pos)

    print('.')
    return Molecule.VolumetricData(vd0.name, data, vd0.extra)
                        
class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('=======================')
        print('    REGRID VOL DATA    ')
        print('=======================')

        # read input file
        mol0 = cmptlib.readMol(args.infile, verbose=args.verbose)

        # regrid volumetric data
        ndim = [args.nx, args.ny, args.nz]
        mol1 = mol0.copy()
        mol1.volumeDataList = regrid_volumetric_data_list(mol0.volumeDataList, ndim)

        # write output file
        cmptlib.writeMol(mol1, args.outfile, verbose=args.verbose)

if __name__ == "__main__":
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='regrid-voldata: regrid volumetric data. Ex) python regrid-voldata.py large.msf small.msf 50 40 30')
    parser.add_argument('infile', help='input file')
    parser.add_argument('outfile', help='output file')
    parser.add_argument('nx', type=int, help='nx')
    parser.add_argument('ny', type=int, help='ny')
    parser.add_argument('nz', type=int, help='nz')
    parser.add_argument('--verbose', type=int, default=4,
                        help='Verbose level index')

    args = parser.parse_args()

    app = myapp()
    app.run(args)
