#!/usr/bin/env python

"""
report-wavecar: Report on a wavecar file
"""

import cmptlib
import wavecar

class myapp:

    def __init__(self):

        pass

    def run(self, args):

        # banner
        print('======================')
        print('    REPORT EAVECAR    ')
        print('======================')

        # read input file
        wfc = wavecar.vaspwfc(args.infile, lsorbit=args.lsorbit)

        # report the WAVECAR
        print('==== WAVECAR REPORT ====')
        print('Filename: ', wfc._fname)
        print('Record length: ', wfc._recl)
        print('LSORBIT: ', wfc._lsoc)
        print('NSPIN: ', wfc._nspin)
        print('NKPTS: ', wfc._nkpts)
        print('NBANDS: ', wfc._nbands)
        print('ENCUT: ', wfc._encut)
        print('ACELL: ', wfc._Acell)
        print('Volume: ', wfc._Omega)
        print('BCELL: ', wfc._Bcell)
        print('NGRID: ', wfc._ngrid)

        ispin = args.ispin
        ikpt = args.ikpt
        iband = args.iband
        
        print('(ispin, ikpt, iband) = ({}, {}, {})'.format(ispin, ikpt, iband))
        wfc.checkIndex(ispin, ikpt, iband)

        if wfc._lsoc:
            spinor = []
            dump = wfc.readBandCoeff(ispin, ikpt, iband)
            nplw = dump.shape[0] // 2
            Cg1 = dump[:nplw]
            Cg2 = dump[nplw:]
            spinor = [Cg1,Cg2]
            print('NPLW: ', nplw

if __name__ == '__main__':
    # stand-alone mode
    import argparse
    parser = argparse.ArgumentParser(
        description='report-wavecar: Report on a WAVECAR file')
    parser.add_argument('infile', help='input file')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbose level index')
    parser.add_argument('--lsorbit', help='LSORBIT or NONCOLLINEAR',
                        action="store_false")
    parser.add_argument('--ispin', type=int, default=0, help='ispin (default=0)')
    parser.add_argument('--ikpt', type=int, default=0, help='ikpt (default=0)')
    parser.add_argument('--iband', type=int, default=0, help='iband (default=0)')

    args = parser.parse_args()

    app = myapp()
    app.run(args)
