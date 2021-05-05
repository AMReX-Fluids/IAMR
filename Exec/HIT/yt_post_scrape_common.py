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

import yt
import numpy as np


#TODO: migrate to ic.txt or elsewhere
gamma = 1.4



#### Helper functions

# Read ic.txt, put values in dictionary
def load_dim_consts(rootdir):
    text_file = open(rootdir + '/ic.txt', 'r')
    lines = text_file.read().replace(' ', '').split('\n')
    headings = lines[0].split(',')
    vals = map(float, lines[1].split(','))
    text_file.close()
    return dict(zip(headings, vals))

# Compute magnitude of field
def mag_sq_func(data_list):
    return reduce(
            lambda a, b: a + b,
            map(lambda c: c**2, data_list)
            )
def mag_func(data_list):
    return np.sqrt(mag_sq_func(data_list))

# Compute kinetic energy
def kin_energy_func(field, data):
    return mag_sq_func([
            data['boxlib', 'x_velocity'],
            data['boxlib', 'y_velocity'],
            data['boxlib', 'z_velocity']
            ])

# Compute vorticity
def magvort_sq_func(field, data):
    return data['boxlib', 'mag_vort'] * data['boxlib', 'mag_vort']

def vort_x_func(field, data):
    return (
        + data['boxlib', 'z_velocity_gradient_y']
        - data['boxlib', 'y_velocity_gradient_z']
        )
def vort_y_func(field, data):
    return (
        - data['boxlib', 'z_velocity_gradient_x']
        + data['boxlib', 'x_velocity_gradient_z']
        )
def vort_z_func(field, data):
    return (
        + data['boxlib', 'y_velocity_gradient_x']
        - data['boxlib', 'x_velocity_gradient_y']
        )
def vort_mag_func(field, data):
    return mag_func([
        data['boxlib', 'z_velocity_gradient_y']
            - data['boxlib', 'y_velocity_gradient_z'],
        - data['boxlib', 'z_velocity_gradient_x']
            + data['boxlib', 'x_velocity_gradient_z'],
        data['boxlib', 'y_velocity_gradient_x']
            - data['boxlib', 'x_velocity_gradient_y']
        ])
def vort_mag_sq_func(field, data):
    return mag_sq_func([
        data['boxlib', 'z_velocity_gradient_y']
            - data['boxlib', 'y_velocity_gradient_z'],
        - data['boxlib', 'z_velocity_gradient_x']
            + data['boxlib', 'x_velocity_gradient_z'],
        data['boxlib', 'y_velocity_gradient_x']
            - data['boxlib', 'x_velocity_gradient_y']
        ])

# Compute dilatation
def divu_sq_func(field, data):
    return data['boxlib', 'divu'] * data['boxlib', 'divu']

def dilatation_func(field, data):
    return (
        data['boxlib', 'x_velocity_gradient_x'] +
        data['boxlib', 'y_velocity_gradient_y'] +
        data['boxlib', 'z_velocity_gradient_z']
        )
def dilatation_sq_func(field, data):
    return (
        data['boxlib', 'x_velocity_gradient_x'] +
        data['boxlib', 'y_velocity_gradient_y'] +
        data['boxlib', 'z_velocity_gradient_z']
        )**2


# Transpose Lists
def transpose(l):
    m = [[x] for x in l[0]]
    for s in l[1:]:
        for k, x in enumerate(s):
            m[k] += [x]
    return m



