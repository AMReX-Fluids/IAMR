################################################################
################################################################
####
#### Post Processing for HIT Simulations
####
#### Authors:       Emmanuel Motheau and John Wakefield
#### Last Edited:   May 10, 2018
####
#### Written for Python 3.
####
####
################################################################
################################################################


#### Imports
from functools import partial
from functools import reduce

import argparse

import os

import pickle

import yt
import numpy as np

from yt_post_scrape_common import *

import spectra


#### Configuration
yt.enable_parallelism()


#### Main Script
if __name__ == '__main__':
    ## Parse Input Arguments
    parser = argparse.ArgumentParser(
            description = 'Post Processing script for Boxlib slices.'
            )
    parser.add_argument(
            'plotfile',
            type = str,
            help = 'First argument must be the name of the plotfile to get slices from.'
            )
    parser.add_argument(
            '-d', '--root-dir',
            help = 'Directory to look in for simulation results.',
            default = '.'
            )
    args = parser.parse_args()
    # Remove Trailing Slashes
    args.plotfile = '/'.join(list(filter(lambda x: len(x), args.plotfile.split('/'))))

    ## Read ic.txt for nondimensionalization constants
    dim_consts = load_dim_consts(args.root_dir)
    dim_consts.update({
        'ke_factor':            3 * dim_consts['urms0']**2,
        'vort_mag_factor':          dim_consts['urms0'] / dim_consts['lambda0'],
        'dilatation_factor':    dim_consts['urms0'] / dim_consts['lambda0']
        })

    ## Load data, initialize storage
    results = {
            'stats':        {},
            'sliceplots':   {}
            }
    ds = yt.load(args.root_dir + '/' + args.plotfile)

    # Add computed fields
    computed_fields = [  # field name, units, displayname
            ('kin_energy', 'cm**2 / s**2', 'Kinetic Energy')
            ]
    for name, unitstr, dispname in computed_fields:
        ds.add_field(
                ('boxlib', name),
                units = unitstr,
                function = eval(name + '_func'),
                sampling_type = 'cell',
                display_name = dispname
                )

    # Nondimensionalize appropriate results
    results['stats'].update({
        'time_adim':            ds.current_time / dim_consts['tau']
        })

    # Line plots
    ray = ds.ortho_ray(
            0,          # x-axis cut
            (1, 1)      # position of cut
            )
    srt = np.argsort(ray['x']) # order by x
    lines = [   # field, adim factor, ylabel
            (
                'divu',
                dim_consts['dilatation_factor'],
                r'\( \theta / (u_\mathrm{rms, 0} / \lambda_0) \)'
                ),
            (
                'mag_vort',
                dim_consts['vort_mag_factor'],
                r'\( \|\nabla \times \mathbf{u}\|_2 / (u_\mathrm{rms, 0} / \lambda_0) \)'
                )
            ]
    for field, adimconst, ylab in lines:
        results['sliceplots'].update({
            'lines_' + field:   {
                'x':    np.array(ray['x'][srt]),
                'y':    np.array(ray[field][srt]) / adimconst,
                'xlabel':   r'\( x \)',
                'ylabel':   ylab,
                'title':    ds.field_info['boxlib', field].get_latex_display_name()
                }
            })

    spect_vars = [  # field list, ylabel, title, precomp filename, dil/sol decomp
            (
                ['x_velocity', 'y_velocity', 'z_velocity'],
                r'\( E_\mathbf{u} \)',
                'Kinetic Energy',
                'vel_spectrum.dat',
                True
                ),
            (
                ['x_vort', 'y_vort', 'z_vort'],
                r'\( E_{\omega} \)',
                'Vorticity',
                'vort_spectrum.dat',
                False
                ),
            (
                ['divu'],
                r'\( E_{\theta} \)',
                'Divergence',
                'divu_spectrum.dat',
                False
                ),
            ]
    precomp_spectra = True
    for _, _, _, specFileName, _ in spect_vars:
        specFile = '/'.join([args.root_dir, args.plotfile, specFileName])
        if not os.path.isfile(specFile):
            precomp_spectra = False
            print("Precomputed spectra file {} not found.".format(specFileName))
    dims = ds.domain_dimensions * int(np.product(ds.ref_factors[:ds.index.max_level]))
    nx, ny, nz = dims
    L = (ds.domain_right_edge - ds.domain_left_edge).d
    if precomp_spectra:
        print("Using precomputed spectra.")
        for i, var_info in enumerate(spect_vars):
            fields, ylabel, title, specFileName, _ = var_info
            specFile = '/'.join([args.root_dir, args.plotfile, specFileName])
            specData = spectra.getData(specFile)
            k = specData[1:, 0]
            E_spect = np.sum(specData[1:, 1:(2 * len(fields)):2], 1)
            #E_spect /= 1e8          # units
            # Plot the stuff
            results['sliceplots'].update({
                'spect_' + title.replace(' ', '_').lower():  {
                    'x':    k,
                    'y':    E_spect,
                    'xlabel':   r'\( k \)',
                    'ylabel':   ylabel,
                    'title':    title
                    }
                })
    else:
        print("Computing spectra with yt and python's (serial) FFT. This may be slow.")
        # Wavenumber plots
        # http://yt-project.org/doc/cookbook/calculating_information.html#making-a-turbulent-kinetic-energy-power-spectrum
        # Make Covering Cube
        cube = ds.covering_grid(
            ds.index.max_level,
            ds.domain_left_edge,
            dims,
            fields = [('boxlib', y) for x in spect_vars for y in x[0]]
            #num_ghost_zones = 1
            )
        # Do FFTs
        for i, var_info in enumerate(spect_vars):
            fields, ylabel, title, _, decomp = var_info
            # Scale to domain
            kx = np.fft.rfftfreq(nx) * nx / L[0]
            ky = np.fft.rfftfreq(ny) * ny / L[1]
            kz = np.fft.rfftfreq(nz) * nz / L[2]
            # Get FFTs
            Kk = np.zeros((nx // 2 + 1, ny // 2 + 1, nz // 2 + 1, 3))
            for fi, field in enumerate(fields):
                dat = cube[('boxlib', field)].d
                Kk[:, :, :, fi] = np.abs(
                        1.0 / nx / ny / nz *
                        np.fft.fftn(dat)[:(nx // 2 + 1), :(ny // 2 + 1), :(nz // 2 + 1)]
                        )
            # Compute Desired Quantities
            mg = np.meshgrid(kx, ky, kz, indexing = 'ij')
            k = np.sqrt(reduce(
                lambda a, b: a + b,
                map(
                    lambda c: c**2,
                    mg
                    )
                ))
            vals = 0.5 * sum([Kk[:, :, :, ik]**2 for ik in range(len(fields))])
            if decomp:
                dil_part_coeffs = np.divide(
                        sum([
                            np.multiply(
                                Kk[:, :, :, ik],
                                mg[ik]
                                )
                            for ik in range(3)
                            ]),
                        k**2)
                dil_part_x = np.multiply(dil_part_coeffs, mg[0])
                dil_part_y = np.multiply(dil_part_coeffs, mg[1])
                dil_part_z = np.multiply(dil_part_coeffs, mg[2])
                sol_part_x = Kk[:, :, :, 0] - dil_part_x
                sol_part_y = Kk[:, :, :, 0] - dil_part_y
                sol_part_z = Kk[:, :, :, 0] - dil_part_z
                dil_vals = 0.5 * (
                        dil_part_x**2 +
                        dil_part_y**2 +
                        dil_part_z**2
                        )
                sol_vals = 0.5 * (
                        sol_part_x**2 +
                        sol_part_y**2 +
                        sol_part_z**2
                        )
            # Bin radially
            kmin = np.min(1.0 / L)
            kmax = np.min(0.5 * dims / L)
            kbins = np.arange(kmin, kmax, kmin)
            N = len(kbins)
            which = np.digitize(k.flat, kbins)
            count = np.bincount(which)
            # Add powers accross wavenumbers
            E_spect = np.zeros(len(count) - 1)
            for n in range(len(count) - 1):
                E_spect[n] = np.sum(vals.flat[which == n + 1])
            if decomp:
                dil_spect = np.zeros(len(count) - 1)
                sol_spect = np.zeros(len(count) - 1)
                for n in range(len(count) - 1):
                    dil_spect[n] = np.sum(dil_vals.flat[which == n + 1])
                    sol_spect[n] = np.sum(sol_vals.flat[which == n + 1])
            # Get middle of bins
            k = 0.5 * (kbins[:N - 1] + kbins[1:])
            # Toss bin 0
            E_spect = E_spect[1:]
            if decomp:
                dil_spect = dil_spect[1:]
                sol_spect = sol_spect[1:]
            # Plot the stuff
            results['sliceplots'].update({
                'spect_' + title.replace(' ', '_').lower():  {
                    'x':    k,
                    'y':    E_spect,
                    'xlabel':   r'\( k \)',
                    'ylabel':   ylabel,
                    'title':    title
                    }
                })
            if decomp:
                results['sliceplots'].update({
                    'dil_spect_' + title.replace(' ', '_').lower():  {
                        'x':    k,
                        'y':    dil_spect,
                        'xlabel':   r'\( k \)',
                        'ylabel':   ylabel,
                        'title':    title
                        }
                    })
                results['sliceplots'].update({
                    'sol_spect_' + title.replace(' ', '_').lower():  {
                        'x':    k,
                        'y':    sol_spect,
                        'xlabel':   r'\( k \)',
                        'ylabel':   ylabel,
                        'title':    title
                        }
                    })

    ## Save Pickle
    results.update({
        'dim_consts':   dim_consts
        })
    with open(args.root_dir + '/' +  args.plotfile + '-slice-scrape.pickle', 'wb') as f:
        pickle.dump(results, f, pickle.HIGHEST_PROTOCOL)

    ## Save CSV
    spectra = list(filter(
            lambda s: "spect_" in s,
            results['sliceplots'].keys()
            ))
    spectra = {k: results['sliceplots'][k] for k in spectra}
    lines = list(filter(
            lambda s: "lines_" in s,
            results['sliceplots'].keys()
            ))
    lines = {k: results['sliceplots'][k] for k in lines}
    for d, xh, cn in [
            (spectra, 'k', 'spectra'),
            (lines  , 'x', 'slice'  )
            ]:
        csv = [[xh] + list(str(x) for x in next(iter(d.values()))['x'])]
        for k, v in d.items():
            csv += [[k] + list(str(x) for x in v['y'])]
        with open(args.root_dir + '/' +  args.plotfile + '-slice-' + cn + '.csv', 'w') as f:
            f.write('\n'.join(','.join(l) for l in transpose(csv)))
    






