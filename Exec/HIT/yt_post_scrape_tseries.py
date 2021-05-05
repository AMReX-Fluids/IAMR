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

import pickle

import yt
import numpy as np

from yt_post_scrape_common import *


#### Configuration
yt.enable_parallelism()


#### Main Script
if __name__ == '__main__':
    ## Parse Input Arguments
    parser = argparse.ArgumentParser(
            description = 'Post Processing script for ``combustion'' type simulations.'
            )
    parser.add_argument(
            'root_dir',
            type = str,
            help = 'First argument must be the name of the plotfile to analyze.'
            )
    parser.add_argument(
            '-m', '--method-name',
            type = str,
            help = 'Name of the method in use.',
            default = 'Unknown Method'
            )
    args = parser.parse_args()


    ## Read ic.txt for nondimensionalization constants
    dim_consts = load_dim_consts(args.root_dir)
    dim_consts.update({
        'ke_factor':            3 * dim_consts['urms0']**2,
        'vort_mag_factor':          dim_consts['urms0'] / dim_consts['lambda0'],
        })

    ## Load data, initialize storage
    ts = yt.load(args.root_dir + '/plt_[0-9]*[0-9]')
    storage = {}

    ## Iterate through timeseries
    for store, ds in ts.piter(storage=storage):
        results = {}

        # Add computed fields
        computed_fields = [  # field name, units, displayname
                ('kin_energy', 'cm**2 / s**2', 'Kinetic Energy'),
                ('magvort_sq', '1', 'Square of Magnitude of Vorticity'),
                ('divu_sq', '1', 'Square of Divergence')
                ]
        for name, unitstr, dispname in computed_fields:
            ds.add_field(
                    ('boxlib', name),
                    units = unitstr,
                    function = eval(name + '_func'),
                    sampling_type = 'cell',
                    display_name = dispname
                    )

        # Compute averages
        ad = ds.all_data()
        fields_to_average = [
            'kin_energy',
            'magvort_sq',
            'divu_sq'
            ]
        avgs = ad.quantities.weighted_average_quantity(
                fields_to_average, 'cell_volume'
                )
        for k, name in enumerate(fields_to_average):
            results.update({
                name + '_avg': avgs[k]
                })

        # Nondimensionalize appropriate results
        results.update({
            'time_adim':            ds.current_time / dim_consts['tau'],
            'kin_energy_avg_adim':  results['kin_energy_avg'] / (3 * dim_consts['urms0']**2),
            'magvort_sq_avg_adim':    results['magvort_sq_avg'] /
                (dim_consts['urms0'] / dim_consts['lambda0'])**2,
            'divu_sq_avg_adim': results['divu_sq_avg'] /
                (dim_consts['urms0'] / dim_consts['lambda0'])**2
            })

        # Store frame of aggregated data
        store.result = results

    ## Shuffle timeseries data into useable format
    keys = dict(storage[0]).keys()
    stats = {}
    for k in keys:
        stats[k] = np.array([storage[ti][k] for ti in range(len(storage))])

    ## Prep Timeseries Plots
    tseries_plots = {
            'kin_energy_avg_adim':  {
                'x' : stats['time_adim'],
                'y' : stats['kin_energy_avg_adim'],
                'xlabel' : r'Time (\( t / \tau \))' ,
                'ylabel' : r'\( \mathbb{E}\left[\|\mathbf{u}\|_2^2\right] / (3 u^2_\mathrm{rms}) \)',
                'title' : 'Kinetic Energy'
                },
            'magvort_sq_avg_adim':  {
                'x' : stats['time_adim'],
                'y' : stats['magvort_sq_avg_adim'],
                'xlabel' : r'Time (\( t / \tau \))' ,
                'ylabel' : r'\( \mathbb{E}\left[\|\nabla \times \mathbf{u}\|_2^2\right] / (u_\mathrm{rms} / \lambda_0)^2 \)',
                'title' : 'Enstrophy'
                },
            'divu_sq_avg_adim' : {
                'x' : stats['time_adim'],
                'y' : stats['divu_sq_avg_adim'],
                'xlabel' : r'Time (\( t / \tau \))' ,
                'ylabel' : r'\( \mathbb{E}\left[\left(\sum_i (\nabla \mathbf{u})_i\right)^2\right] / (u_\mathrm{rms} / \lambda_0)^2 \)',
                'title' : 'Dilatation'
                }
            }

    ## Save Pickle
    with open(args.root_dir + '/tseries' + '.pickle', 'wb') as f:
        pickle.dump(
                {
                    'dim_consts':   dim_consts,
                    'tseries':      tseries_plots,
                    'methodname':   args.method_name,
                    'stats':        stats
                    },
                f,
                pickle.HIGHEST_PROTOCOL
                )

    ## Save CSV
    csv = [['tau'] + list(str(x) for x in next(iter(tseries_plots.values()))['x'])]
    for k, v in tseries_plots.items():
        csv += [[k] + list(str(x) for x in v['y'])]
    with open(args.root_dir + '/tseries' + '.csv', 'w') as f:
        f.write('\n'.join(','.join(l) for l in transpose(csv)))




