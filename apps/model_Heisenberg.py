#!/usr/bin/env python

"""
 Heisenberg spin model
"""

import numpy as np
import matplotlib.pyplot as plt

from lmfit import Minimizer, Parameters, report_fit

def read_datafile(fname):
    data = np.loadtxt(fname)
    xdata = data[:, 0]
    ydata = data[:, 1]
    ixdata = np.indices([len(xdata)])
    return dict({'ixdata':ixdata, 'xdata':xdata, 'ydata':ydata})

# Model function
def func(params, x):
    """
    model function
    """
    amp = params['amp']
    cen = params['cen']
    wid = params['wid']
    y = (amp / (np.sqrt(2*np.pi) * wid)) * np.exp(-(x-cen)**2 / (2*wid**2))
    return y

# define objective function: returns the array to be minimized
def fcn2min(params, xdata, ydata):
    """
    Objective function to minimize
    """
    model = func(params, xdata)
    return model - ydata

class Heisenberg:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('======================')
        print('   Heisenberg Model   ')
        print('======================')

        # read input data file
        data = read_datafile(args.datafile)
        ixdata = data['ixdata']
        xdata = data['xdata']
        ydata = data['ydata']

        # create a set of Parameters
        params = Parameters()
        params.add('amp', value=5.0)
        params.add('cen', value=5.0)
        params.add('wid', value=1.0, min=0.0)

        # do fit, here with leastsq model
        minner = Minimizer(fcn2min, params, fcn_args=(xdata, ydata))
        result = minner.minimize()

        # calculate final result
        final = ydata + result.residual

        # write error report
        report_fit(result)

        # try to plot results
        try:
            import matplotlib.pyplot as plt
            plt.plot(xdata, ydata, 'bo')
#            plt.plot(xdata, final, 'r')
            plt.plot(xdata, final, 'ro')
            plt.show()
        except ImportError:
            print('matplot not available')


if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='Heisenberg spin model')
    parser.add_argument('--reportlevel', type=int, default=4,
                        help='Report level index (0-5)')
    parser.add_argument('datafile', help='input data file')

    args = parser.parse_args()

    app = Heisenberg()
    app.run(args)
        
