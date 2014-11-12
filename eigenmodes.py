#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')   # headless
import matplotlib.pyplot as plt
import os
import numpy
import math
from sympy import mpmath
import numpy as np

import fields 
import data_manipulation 
import signals 
import shotput
import geometry

def resistivity_matrix(count):
    rho_stainless_steel = 6.90e-7
    return rho_stainless_steel*numpy.identity(count)


def self_inductance(x, z, a):
    return x*(math.log(8.*x/a)-1.75)


def green(x, z, xc, zc):
    k2 = 4*xc*x/((xc + x)**2 + (z - zc)**2)
    return -math.sqrt(x*xc/k2)*((2-k2)*mpmath.ellipk(k2)-2*mpmath.ellipe(k2))


def inductance((x1, z1), (x2, z2), a):
    return self_inductance(x1, z1, a) if (x1,z1) == (x2,z2) \
           else -green(x1, z1, x2, z2)


def inductance_matrix(filaments):
    count = len(filaments)
    inductances = np.zeros((count, count))
    for i in range(count):
        for j in range(count):
            inductances[i][j] = inductance(filaments[i], filaments[j], geometry.asize)
    return inductances
    

geometry.plot_geometry()

a = inductance_matrix(geometry.top_shell_filaments(10)
                        +geometry.bottom_shell_filaments_mirror(10))
print a
print a.shape


#for i, sensor in enumerate(sensors):
#    startTime = time.time()
#    (sensor_time, sensor_signal) = get_sensor_time_signal(shot_num, sensor.name)
#   # Fix length of integrated data
#    (sensor_time, sensor_signal) = clip(sensor_time, sensor_signal) 
#
#   # Radial sensors point inwards
#    if 'TA' in sensor.name: 
#        if 'R' in sensor.name:
#            sensor.n_r = -1.0
#
#    # Calculate field.
#    field_vals = B_signal(vf_signal, VFR, VFZ, sensor.r, sensor.z, sensor.n_r, sensor.n_z)
#
#    # Text output.
#    print '------------------------------------------' + sensor.name + ':'
#    # Image output.
#    endTime = time.time()
#    elapsed = endTime-startTime
#    print '%.2f seconds elapsed.' % elapsed

