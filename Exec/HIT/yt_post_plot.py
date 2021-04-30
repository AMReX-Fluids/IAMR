################################################################
################################################################
####
#### Make Plots from Post Processing Scrapes
####
#### Authors:       Emmanuel Motheau and John Wakefield
#### Last Edited:   May 23, 2018
####
####
####
################################################################
################################################################


#### Imports
from functools import partial
from functools import reduce
import itertools

import argparse

import pickle

import numpy as np
import scipy.linalg as la

import matplotlib.pyplot as plt


#### Hard Coded Arguments

plot_kolm = False


#### Plot Configuration
font = {
        #'family':   'sans-serif',
        #'weight':   'normal',
        #'size':     '10'
        }
legendfont = {
        #'family':   'sans-serif',
        #'weight':   'normal',
        'size':     '6'
        }
plot_dpi = 400
figargs = {
        'figsize':      (7, 5)
        }
plotargs = {
        'linewidth':    0.6,
        'markersize':   0.9
        }
plt.rc('text', usetex = True)
plt.rc('text.latex', preamble = r'\usepackage{amsfonts}')
plt.rc('font', **font)


#### Main Script
if __name__ == '__main__':
    ## Parse Input Arguments
    parser = argparse.ArgumentParser(
            description = 'Plots slices results from scraping one or more ``Combustion'' type simulations.'
            )
    parser.add_argument(
            'action',
            type = str,
            help = 'Which plots to make, currently ``slices'' or ``tseries''.'
            )
    parser.add_argument(
            'outname',
            type = str,
            help = 'Output file name prefix.'
            )
    parser.add_argument(
            'files',
            type = str,
            nargs = '+',
            help = 'Pickle files to plot from.'
            )
    parser.add_argument(
            '-n', '--no-name',
            help = 'Do not look for name information, do not make legends.',
            action = 'store_true'
            )
    parser.add_argument(
            '-d', '--pickle-dir',
            help = 'Name of pickle directory.',
            default = '.'
            )
    parser.add_argument(
            '-p', '--plots-dir',
            help = 'Name of subfolder in which to put plots.',
            default = 'plots'
            )
    parser.add_argument(
            '-r', '--ref-soln',
            help = 'If provided, the plotfile containing the reference solution.',
            default = None
            )
    parser.add_argument(
            '--use-set-formats',
            help = 'If present, plot line formats are found in a table; they are not assigned blindly.',
            action = 'store_true',
            default = False
            )
    parser.add_argument(
            '--allow-dil-sol',
            help = 'If present, enables the plotting of the dilatational and solenoidal parts of spectra seperately, if they are present in the input file.',
            action = 'store_true'
            )
    args = parser.parse_args()

    ## Prepare filename prefix
    fnameprefix = '/'.join([args.plots_dir, args.outname]) + '-'

    ## Colors
    if args.use_set_formats:
        colors = {
                'MDCD':     'g',
                'PeleC':    'c',
                'WENO':     'm',
                'WENOJS':   'm',
                'WENOM':    'C0',
                'WENOT':    'C4',
                'TENO':     'C4',
                'WENOZ':    'C3',
                'WENO3':    'C5'
                }
        dots = {
                '032':      '',
                '064':      '-',
                '128':      ':',
                '256':      '--'
                }
        symbs = {
                'HLL':      '>',
                'HLLC':     'x',
                'JBB':      '.',
                'CGF':      's'
                }
        altsymbs = {
                'HLL':      '<',
                'HLLC':     '+',
                'JBB':      'o',
                'CGF':      'd'
                }
    else:
        colors = [
            'b', 'r', 'c', 'm', 'y'
            ]
        dots = [
            ':', '--', '-'
            ]
        symbs = [
            '.', 'x', 'd', '^', 'v', '<', '>', 's'
            ]
    fmts = [colors, dots, symbs]
    ref_fmt = 'k-'

    ## Print Task
    print("Making {} plots with name prefix {}.".format(args.action, args.outname))

    ## Make Filename List
    fnames = list(args.files)
    if args.ref_soln is not None:
        fnames += [args.ref_soln]

    ## Load Pickles
    data = {}
    for fname in fnames:
        with open('/'.join([args.pickle_dir, fname]), 'rb') as f:
            data[fname.replace('.', '')] = pickle.load(f)
    print("Loaded data from {} methods.".format(len(data.keys())))
    method0 = data[list(data.keys())[0]]

    ## Parse filenames
    def getparts(name):
        return (name
                .replace('.', '')
                .replace('_', '')
                .replace('-final_slicepickle', '')
                .replace('pickle', '')
                .split('-')
                )
    parts = []
    for f in fnames:
        for k, g in enumerate(getparts(f)):
            if k == len(parts):
                parts.append(set([]))
            parts[k].add(g)
    parts = [dict((k, v) for v, k in enumerate(p)) for p in parts]
    def getfmt(name, alt = False):
        p = getparts(name)
        if args.ref_soln is not None and p == getparts(args.ref_soln):
            return ref_fmt
        if args.use_set_formats:
            if alt:
                return colors[p[1]] + dots[p[2]] + altsymbs[p[0]]
            else:
                return colors[p[1]] + dots[p[2]] + symbs[p[0]]
        else:
            unique_parts = 0
            for dct in parts:
                if len(dct.keys()) > 1:
                    unique_parts += 1
            fmt = ''
            j = 0   # format index
            k = 0   # part index
            while k < len(p) and j < len(fmts):
                if len(parts[k]) == 1:
                    k += 1
                else:
                    fmt += fmts[j][parts[k][p[k]] % len(fmts[j])]
                    j += 1
                    if j < unique_parts:    # don't go past last unique part
                        k += 1
            return fmt

    ## Perform requested action
    if args.action.lower() == 'slices':
        # Plotting of Timeslice Statistics and Graphics
        for method, res in data.items():
            print("{} reported final time {}.".format(
                method,
                res['stats']['time_adim']
                ))

        # Fix Naming from PeleC Files
        for f, res in data.items():
            if 'lines_Temp' in res['sliceplots'].keys():
                res['sliceplots']['lines_temperature'] = res['sliceplots']['lines_Temp']
                print('Corrected temp to temperature in line plot.')

        # Line plots
        lines = [
                'lines_divu',
                'lines_density',
                'lines_temperature',
                'lines_magvort'
                ]
        lines_titles = [    # Done to correct from an error in scraping
                'Divergence',
                'Density',
                'Temperature',
                'Vorticity'
                ]
        lines_fig, lines_axs = plt.subplots(2, 2, **figargs)
        lines_axs = lines_axs.flatten()
        methodnames = []
        for i, line in enumerate(lines):
            for f, res in data.items():
                lines_axs[i].plot(
                        res['sliceplots'][line]['x'],
                        res['sliceplots'][line]['y'],
                        getfmt(f),
                        **plotargs
                        )
                methodnames += [res['methodname'].replace('\'', '')]
            # Labels come from whichever method was last; this is fine.
            lines_axs[i].set_xlabel(res['sliceplots'][line]['xlabel'])
            lines_axs[i].set_ylabel(res['sliceplots'][line]['ylabel'])
            #lines_axs[i].set_title(res['sliceplots'][line]['title'])
            lines_axs[i].set_title(lines_titles[i])
        methodnames = methodnames[0:len(args.files)]
        if not args.no_name:
            lines_axs[0].legend(methodnames, loc = 'upper left', prop = legendfont)
        lines_fig.tight_layout()
        lines_fig.savefig(fnameprefix + 'spatial_quad' + '.pdf', dpi = plot_dpi)

        # Wavenumber plots
        spects = [
                ['spect_kinetic_energy', 'spect_velocity'],
                ['spect_vorticity'],
                ['spect_divergence', 'spect_dilatation'],
                ['spect_density']
                ]
        spect_fig, spect_axs = plt.subplots(2, 2, **figargs)
        spect_axs = spect_axs.flatten()
        methodnames = []
        legend_extra = 0
        first = True
        for i, spect_opts in enumerate(spects):
            xmax = []
            ymin = []
            ymax = []
            for f, res in data.items():
                spect = set(spect_opts).intersection(set(res['sliceplots'].keys())).pop()
                if  args.allow_dil_sol                          \
                and 'dil_' + spect in res['sliceplots'].keys()  \
                and 'sol_' + spect in res['sliceplots'].keys():
                    spect_axs[i].loglog(
                            res['sliceplots']['dil_' + spect]['x'][:-1],
                            res['sliceplots']['dil_' + spect]['y'][:-1],
                            getfmt(f, False),
                            **plotargs
                            )
                    methodnames += [res['methodname'].replace('\'', '') + ' (Comp)']
                    spect_axs[i].loglog(
                            res['sliceplots']['sol_' + spect]['x'][:-1],
                            res['sliceplots']['sol_' + spect]['y'][:-1],
                            getfmt(f, True),
                            **plotargs
                            )
                    methodnames += [res['methodname'].replace('\'', '') + ' (Incomp)']
                    legend_extra += 1
                    xmax += [res['sliceplots']['dil_' + spect]['x'][-1]]
                    xmax += [res['sliceplots']['sol_' + spect]['x'][-1]]
                    ymin += [res['sliceplots']['dil_' + spect]['y'][-2]]
                    ymin += [res['sliceplots']['sol_' + spect]['y'][-2]]
                    ymax += [max(res['sliceplots']['dil_' + spect]['y'])]
                    ymax += [max(res['sliceplots']['dil_' + spect]['y'])]
                else:
                    spect_axs[i].loglog(
                            res['sliceplots'][spect]['x'][:-1],
                            res['sliceplots'][spect]['y'][:-1],
                            getfmt(f),
                            **plotargs
                            )
                    methodnames += [res['methodname'].replace('\'', '')]
                    xmax += [res['sliceplots'][spect]['x'][-1]]
                    ymin += [res['sliceplots'][spect]['y'][-2]]
                    ymax += [max(res['sliceplots'][spect]['y'])]
                # Plot Kolmogorov
                if first and not i:
                    kolm_k = res['sliceplots'][spect]['x']
                    kolm_Ek = res['sliceplots'][spect]['y']
                    beg_pt = len(kolm_k) // 3
                    ref_pt = 2 * len(kolm_k) // 3
                    exp_pt = -3.0
                    first = False
            if plot_kolm:
                spect_axs[0].loglog(
                        kolm_k[beg_pt:],
                        kolm_Ek[ref_pt] * kolm_k[beg_pt:]**exp_pt / kolm_k[ref_pt]**exp_pt,
                        'k:'
                        )
            # Labels come from whichever method was last; this is fine.
            spect_axs[i].set_xlabel(res['sliceplots'][spect]['xlabel'])
            spect_axs[i].set_ylabel(res['sliceplots'][spect]['ylabel'])
            spect_axs[i].set_title( res['sliceplots'][spect]['title' ])
            xmax.sort()
            ymin.sort()
            ymax.sort()
            xmax = 1.2 * xmax[-2]
            ymin = 0.2 * ymin[ 1]
            ymax = 1.2 * ymax[-2]
            spect_axs[i].set_xlim([10**(-15) , xmax])
            spect_axs[i].set_ylim([ymin, ymax])
            plt.xticks(
                    [2, 3, 4, 6],
                    ['2', '3', '4', '6']
                    )
        methodnames = methodnames[0:len(args.files) + legend_extra]
        if not args.no_name:
            spect_axs[0].legend(methodnames, loc = 'lower left', prop = legendfont)
        spect_fig.tight_layout()
        spect_fig.savefig(fnameprefix + 'wavenumber_quad' + '.pdf', dpi = plot_dpi)
    elif args.action.lower() == 'tseries':
        ## Make time series plot
        tvars = [
                'kin_energy_avg_adim',
                'magvort_sq_avg_adim',
                'temp_var_sq_avg_adim',
                'divu_sq_avg_adim'
                ]
        tseries_fig, tseries_axs = plt.subplots(2, 2, **figargs)
        tseries_axs = tseries_axs.flatten()
        methodnames = []
        for i, tvar in enumerate(tvars):
            for f, res in data.items():
                tseries_axs[i].plot(
                        res['tseries'][tvar]['x'],
                        res['tseries'][tvar]['y'],
                        getfmt(f),
                        **plotargs
                        )
                methodnames += [res['methodname']]
            # Labels come from whichever method was last; this is fine.
            tseries_axs[i].set_xlabel(res['tseries'][tvar]['xlabel'])
            tseries_axs[i].set_ylabel(res['tseries'][tvar]['ylabel'])
            tseries_axs[i].set_title(res['tseries'][tvar]['title'])
        methodnames = methodnames[0:len(args.files)]
        if not args.no_name:
            tseries_axs[0].legend(methodnames, loc = 'upper right', prop = legendfont)
        tseries_fig.tight_layout()
        tseries_fig.savefig(fnameprefix + 'tseries_quad' + '.pdf', dpi = plot_dpi)
    else:
        raise Exception('Plotting method (action) not implemented.')




