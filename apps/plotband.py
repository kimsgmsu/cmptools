#!/usr/bin/env python
# -*- coding=utf-8 -*-

import sys
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec

from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import Spin, OrbitalType


def rgbline(ax, k, e, red, green, blue, alpha=1.):
    # creation of segments based on
    # http://nbviewer.ipython.org/urls/raw.github.com/dpsanders/matplotlib-examples/master/colorline.ipynb
    pts = np.array([k, e]).T.reshape(-1, 1, 2)
    seg = np.concatenate([pts[:-1], pts[1:]], axis=1)

    nseg = len(k) - 1
    r = [0.5 * (red[i] + red[i + 1]) for i in range(nseg)]
    g = [0.5 * (green[i] + green[i + 1]) for i in range(nseg)]
    b = [0.5 * (blue[i] + blue[i + 1]) for i in range(nseg)]
    a = np.ones(nseg, np.float) * alpha
    lc = LineCollection(seg, colors=list(zip(r, g, b, a)), linewidth=2)
    ax.add_collection(lc)

class myKpoint:
    def _init__(self, pos, label):
        self.pos = pos
        self.label = label
        
class myKpointPath:
    def _init__(self, kpt1, kpt2):
        self.kpt1 = kpt1
        self.kpt2 = kpt2
        
class plotband:

    def __init__(self):

        pass

    def _get_kpoint_path_labels(self, kpointfile):

        # kpoints labels
        kls = Kpoints.from_file(args.kpoints).labels
        p2 = kls[0]
        labels = [p2]
        ifp1 = True
        for k in kls:
            if ifp1:
                if not k == p2:
                    labels[-1] += '/{}'.format(k)
            else:
                p2 = k
                labels.append(p2)
            ifp1 = not ifp1

        print('lables = ', labels)
        return labels

    def run(self, args):

        # banner
        print('=========================')
        print('   PLOT BAND STRUCTURE   ')
        print('=========================')

        # make band diagram
        self.make_band_diagram(args)

        # plot the band diagram
        if args.format == "gnuplot":
            self.make_gnuplot(args)
        else:
            self.make_matplotlib(args)

    def make_band_diagram(self, args):

        # dosrun
        self.dosrun = Vasprun(args.dosrun, parse_potcar_file=False)

        # vasprun
        bandrun = Vasprun(args.bandrun, parse_projected_eigen=True, parse_potcar_file=False)
        self.bands = bandrun.get_band_structure(args.kpoints,
                                                line_mode=True)
        self.efermi = self.bands.efermi
        # energy shift
        if args.shift_energy:
            self.eshift = self.efermi
        else:
            self.eshift = 0.0
    
        # kpoints labels
        sklabel_list = Kpoints.from_file(args.kpoints).labels
        assert len(sklabel_list) % 2 == 0, "No. of special kpoints must be even in line-mode"
        num_paths = len(sklabel_list) // 2
        self.kpath_list = []
        for ip in range(num_paths):
            self.kpath_list.append(
                [sklabel_list[ip*2], sklabel_list[ip*2+1]])
        self.label = []
        lab2 = self.kpath_list[0][0]
        for path in self.kpath_list:
            lab1 = path[0]
            if path[0]==lab2:
                lab1 = path[0]
            else:
                lab1 = '{}|{}'.format(lab2,path[0])
            self.label.append(lab1)
            lab2 = path[1]
        self.label.append(lab2)

        # full kpoints
        self.kpoint_list = bandrun.actual_kpoints
        self.num_kpts_per_path = len(self.kpoint_list) // num_paths
        assert num_paths * self.num_kpts_per_path == len(self.kpoint_list), "Total number of kpts must be a multiple of number of paths"
        
        # set k-coordinates
        reclatt = self.bands.lattice_rec
        self.xval=np.zeros(len(self.kpoint_list))
        self.label_kx = [0.0]
        kx0 = 0.0
        for ip in range(len(self.kpath_list)):
            ik1 = ip * self.num_kpts_per_path
            ik2 = ik1 + self.num_kpts_per_path -1
            if args.fixed_pathsize:
                dk = 1.0 / self.num_kpts_per_path
            else:
                kpt1 = np.array(self.kpoint_list[ik1])
                kpt2 = np.array(self.kpoint_list[ik2])
                dkvec = (kpt2 - kpt1) / self.num_kpts_per_path
                dk = reclatt.dot(dkvec, dkvec, frac_coords=True)

            for jkpt in range(self.num_kpts_per_path):
                kx = kx0 + jkpt*dk
                self.xval[ik1+jkpt] = kx
            kx0 = self.xval[ik2]
            self.label_kx.append(kx0)

        # Set x-y range
        kmin = 1.0e+99
        kmax = -1.0e+99
        emin =  1.0e+99
        emax = -1.0e+99
        for b in range(self.bands.nb_bands):
            for k in range(len(self.bands.kpoints)):
                x = self.xval[k]
                y = self.bands.bands[Spin.up][b][k]
                kmin = min(kmin, x)
                kmax = max(kmax, x)
                emin = min(emin, y)
                emax = max(emax, y)
                if self.bands.is_spin_polarized:
                    y = self.bands.bands[Spin.down][b][k]
                    emin = min(emin, y)
                    emax = max(emax, y)

        # Add margins
        margin_frac = 0.1
        self.kmin = kmin
        self.kmax = kmax
        self.emin = emin - (emax-emin)*margin_frac
        self.emax = emax + (emax-emin)*margin_frac

    def make_matplotlib(self, args):

        # Make matplotlib plot

        # kpoints labels
        kls = Kpoints.from_file(args.kpoints).labels
        p2 = kls[0]
        labels = [p2]
        ifp1 = True
        for k in kls:
            if ifp1:
                if not k == p2:
                    labels[-1] += '/{}'.format(k)
            else:
                p2 = k
                labels.append(p2)
            ifp1 = not ifp1

        print('lables = ', labels)
        #    labels = [r"$L$", r"$\Gamma$", r"$X$", r"$U,K$", r"$\Gamma$"]

        # density of states
        dosrun = Vasprun(args.dosrun)
        spd_dos = dosrun.complete_dos.get_spd_dos()

        # bands
        run = Vasprun(args.bandrun, parse_projected_eigen=True)
        bands = run.get_band_structure(args.kpoints,
                                       line_mode=True,
                                       efermi=dosrun.efermi)

        # set up matplotlib plot
        # ----------------------

        # general options for plot
        font = {'family': 'serif', 'size': 24}
        plt.rc('font', **font)

        # set up 2 graph with aspec ration 2/1
        # plot 1: bands diagram
        # plot 2: Density of States
        gs = GridSpec(1, 2, width_ratios=[2, 1])
        fig = plt.figure(figsize=(11.69, 8.27))
        fig.suptitle("Bands diagram of silicon")
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])  # , sharey=ax1)

        # set ylim for the plot
        # ---------------------
        emin = 1e100
        emax = -1e100
        for spin in bands.bands.keys():
            for b in range(bands.nb_bands):
                emin = min(emin, min(bands.bands[spin][b]))
                emax = max(emax, max(bands.bands[spin][b]))

        emin -= bands.efermi + 1
        emax -= bands.efermi - 1
        ax1.set_ylim(emin, emax)
        ax2.set_ylim(emin, emax)

        # Band Diagram
        # ------------
        name = "Si"
        pbands = bands.get_projections_on_elements_and_orbitals({name: ["s", "p", "d"]})

        # compute s, p, d normalized contributions
        contrib = np.zeros((bands.nb_bands, len(bands.kpoints), 3))
        for b in range(bands.nb_bands):
            for k in range(len(bands.kpoints)):
                sc = pbands[Spin.up][b][k][name]["s"]**2
                pc = pbands[Spin.up][b][k][name]["p"]**2
                dc = pbands[Spin.up][b][k][name]["d"]**2
                tot = sc + pc + dc
                if tot != 0.0:
                    contrib[b, k, 0] = sc / tot
                    contrib[b, k, 1] = pc / tot
                    contrib[b, k, 2] = dc / tot

        # plot bands using rgb mapping
        for b in range(bands.nb_bands):
            rgbline(ax1,
                    range(len(bands.kpoints)),
                    [e - bands.efermi for e in bands.bands[Spin.up][b]],
                    contrib[b, :, 0],
                    contrib[b, :, 1],
                    contrib[b, :, 2])

        # style
        ax1.set_xlabel("k-points")
        ax1.set_ylabel(r"$E - E_f$   /   eV")
        ax1.grid()

        # fermi level at 0
        ax1.hlines(y=0, xmin=0, xmax=len(bands.kpoints), color="k", lw=2)

        # labels
        nlabs = len(labels)
        step = len(bands.kpoints) / (nlabs - 1)
        for i, lab in enumerate(labels):
            ax1.vlines(i * step, emin, emax, "k")
        ax1.set_xticks([i * step for i in range(nlabs)])
        ax1.set_xticklabels(labels)

        ax1.set_xlim(0, len(bands.kpoints))
        #ax1.set_title("Bands diagram")

        # Density of states
        # -----------------

        ax2.set_yticklabels([])
        ax2.grid()
        ax2.set_xlim(1e-6, 3)
        ax2.set_xticklabels([])
        ax2.hlines(y=0, xmin=0, xmax=8, color="k", lw=2)
        ax2.set_xlabel("Density of States", labelpad=28)
        #ax2.set_title("Density of States")

        # spd contribution
        ax2.plot(spd_dos[OrbitalType.s].densities[Spin.up],
                 dosrun.tdos.energies - dosrun.efermi,
                 "r-", label="3s", lw=2)
        ax2.plot(spd_dos[OrbitalType.p].densities[Spin.up],
                 dosrun.tdos.energies - dosrun.efermi,
                 "g-", label="3p", lw=2)
        ax2.plot(spd_dos[OrbitalType.d].densities[Spin.up],
                 dosrun.tdos.energies - dosrun.efermi,
                 "b-", label="3d", lw=2)

        # total dos
        ax2.fill_between(dosrun.tdos.densities[Spin.up],
                         0,
                         dosrun.tdos.energies - dosrun.efermi,
                         color=(0.7, 0.7, 0.7),
                         facecolor=(0.7, 0.7, 0.7))

        ax2.plot(dosrun.tdos.densities[Spin.up],
                 dosrun.tdos.energies - dosrun.efermi,
                 color=(0.3, 0.3, 0.3),
                 label="total DOS")

        # plot format style
        # -----------------
        ax2.legend(fancybox=True, shadow=True, prop={'size': 18})
        plt.subplots_adjust(wspace=0)

        plt.show()
        plt.savefig(sys.argv[0].strip(".py") + ".pdf", format="pdf")

    def make_gnuplot(self, args):

        # Make gnuplot plot

        # gnuplot data file
        with open('{}-band.dat'.format(args.job), 'w') as datafile:
            for b in range(self.bands.nb_bands):
                for k in range(len(self.bands.kpoints)):
                    x = self.xval[k]
                    y = self.bands.bands[Spin.up][b][k]
                    if self.bands.is_spin_polarized:
                        y2 = self.bands.bands[Spin.down][b][k]
                        datafile.write('{:12.5g} {:12.5g} {:12.5g} \n'.format(x,y,y2))
                    else:
                        datafile.write('{:12.5g} {:12.5g}\n'.format(x,y))
                datafile.write('\n')

        # gnuplot control file
        try:
            fg = open('{}-band.gpl'.format(args.job), 'w')
            print('%plotband: Writing gnuplot script file "{}"'.format(fg.name))
        except IOError:
            util.Error('%plotband: Failed to open output file "{}"'.format('{}-band.gpl'.format(args.job)))

        # preamble
        print('set style data dots', file=fg)
        print('set xlabel "Special k-points"', file=fg)
        if args.shift_energy:
            print('set ylabel "Energy - Ef (eV)"', file=fg)
        else:
            print('set ylabel "Energy (eV)"', file=fg)
        print('set key top right', file=fg)
        print('set style arrow 1 nohead lw 1 lt 1 lc 0', file=fg)
        print('set style arrow 2 nohead lw 2 lt 0 lc 2 dt ".. "', file=fg)
        print('yshift = {}'.format(self.eshift), file=fg)
        # data range
        print('xmin = {}'.format(self.kmin), file=fg)
        print('xmax = {}'.format(self.kmax), file=fg)
        print('ymin = {} - yshift'.format(self.emin), file=fg)
        print('ymax = {} - yshift'.format(self.emax), file=fg)
        print('set xrange [xmin:xmax]', file=fg)
        print('set yrange [ymin:ymax]', file=fg)
        # x ticks
        print('set xtics (', file=fg, end='')
        for i in range(len(self.label)):
            if i > 0:
                print(', ', file=fg, end='')
            print('"{}" {:12.5g}'.format(
                self.label[i], self.label_kx[i]), file=fg, end='')
        print(')', file=fg)

        # vertical lines
        fmt='set arrow from {0:12.5g},ymin to {0:12.5g},ymax as 1'
        for i in range(len(self.label)-1):
            print(fmt.format(self.label_kx[i]), file=fg)

        # Fermi energy
        fmt='set arrow from xmin,({0:12.5g}-yshift) to xmax,({0:12.5g}-yshift) as 2'
        print(fmt.format(self.efermi), file=fg)

        # plot data
        if args.plottype == 'points':
            plot_com = 'points pt 6'
        elif args.plottype == 'lines':
            plot_com = 'lines'
        else:
            plot_com = 'linespoints pt 6 '
        s = 'plot'
        s += ' "{}-band.dat"'.format(args.job)
        s += ' using 1:($2-yshift)'
        if self.bands.is_spin_polarized:
            s += ' with {} lc rgb "red" title "spin.up"'.format(plot_com)
        else:
            s += ' with {} lc rgb "blue" title "spin.up"'.format(plot_com)
        print(s, end='', file=fg)
        if self.bands.is_spin_polarized:
            print(', \\', file=fg)
            s = ' "{}-band.dat"'.format(args.job)
            s += ' using 1:($3-yshift) '
            s += ' with {} lc rgb "blue" title "spin.dn"'.format(plot_com)
            print(s, end='', file=fg)
        print('', file=fg)

        # close gnuplot script file
        fg.close()

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description='plotband: Plot band structure.  Ex) plotband dosrun bandrun kpoints')
    parser.add_argument('dosrun', help='vasprun.xml for DOS')
    parser.add_argument('bandrun', help='vasprun.xml for BAND')
    parser.add_argument('kpoints', type=str,
                        help='KPOINTS file for BAND')
    parser.add_argument('--job', type=str, default="PLOTBAND",
                        help='Job name. Used for making output file names')
    parser.add_argument('--format',  default="gnuplot",
                        help='format of the output plot')
    parser.add_argument('--fixed-pathsize', dest='fixed_pathsize', action='store_true',
                        help='make the pathsize fixed')
    parser.add_argument('--no-fixed-pathsize', dest='fixed_pathsize', action='store_false',
                        help='make pathsize vary proportional to its length')
    parser.set_defaults(fixed_pathsize=False)
    parser.add_argument('--shift-energy-to-Ef', dest='shift_energy', action='store_true',
                        help='Shift energies to Fermi energy')
    parser.add_argument('--no-shift-energy-to-Ef', dest='shift_energy', action='store_false',
                        help='Do not shift energies to Fermi energy')
    parser.set_defaults(shift_energy=True)
    parser.add_argument('--plottype', type=str, default="lines",
                        choices=['points', 'lines', 'linespoints'],
                        help='plot type')

    args = parser.parse_args()

    app = plotband()
    app.run(args)
