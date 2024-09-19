#!/usr/bin/env python
"""
plot-voldata: Plot volumetric data along a line
"""

import os
import numpy as np
import cmptlib, util

def integrate_volume_data_list(vdList, args):
    """
    Plot a volumetric data list
    """
    plotdata = None
    if vdList is None:
        return plotdata

    for k, vd in enumerate(vdList):
        if k == args.dataset:
            plotdata = plot_volume_data(vd, args)

    return plotdata

def plot_volume_data(vd, args):

    from scipy.interpolate import RegularGridInterpolator

    ndim = vd.data.shape
    nx = ndim[0]
    ny = ndim[1]
    nz = ndim[2]
    xlist = np.linspace(0.0, 1.0, nx)
    ylist = np.linspace(0.0, 1.0, ny)
    zlist = np.linspace(0.0, 1.0, nz)
    my_func = RegularGridInterpolator((xlist, ylist, zlist), vd.data)

    data = np.zeros( (nx, ny, nz), dtype=np.float_ )

    hx = 1.0/np.float_(nx-1)
    hy = 1.0/np.float_(ny-1)
    hz = 1.0/np.float_(nz-1)

    print('*** Plotting volume data on {}x{}x{} grid ***'.format(nx, ny, nz))
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
        print('========================')
        print('    PLOT VOLUME DATA    ')
        print('========================')

        # read molecule
        mol = cmptlib.readMol(args.mol)

        # plot the volumetric data
        plot_volume_data_list(mol.volumeDataList, args)

if __name__ == "__main__":
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='plot-voldata: Plot volumetric data')
    parser.add_argument('mol', help='Input molecule')
    parser.add_argument('outfile', help='Output data file')
    parser.add_argument('--start', default="[0,0,0]",
                        help='start point')
    parser.add_argument('--end', default="[0,0,1]",
                        help='end point')
    parser.add_argument('--coord', default='frac', choices=['frac', 'cart'],
                        help='coord system for start and end points')
    parser.add_argument('--npoint', type=int, default=100,
                        help='number of points to plot')
    parser.add_argument('--dataset', type=int, default=0,
                        help='volume data set to plot: 0 = first data set')

    args = parser.parse_args()

    app = myapp()
    app.run(args)
