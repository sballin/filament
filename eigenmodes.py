#!/usr/bin/env python

import time
import matplotlib
matplotlib.use('agg')   # headless
import matplotlib.pyplot as plt
from pylab import *
import MDSplus
import pickle
import csv
import os
import numpy
import math
from fields import *
from data_manipulation import *
from signals import *
from shotput import *

def top_shell_filaments(count):
    major_radius = 0.92 
    nominal_minor_radius = 0.16
    shell_thickness = 0.004763
    minor_radius = nominal_minor_radius + shell_thickness
    shell_length = math.pi/float(2)*minor_radius
    top_part = 0.0635
    gap = (top_part + shell_length)/float(count)
    top_part_filaments = [(i, minor_radius) for i in 
                          numpy.arange(major_radius-top_part, major_radius, gap)]
    angle_increment = gap/shell_length*math.pi/float(2)
    curve_filaments = [(major_radius+minor_radius*math.cos(theta),
                      minor_radius*math.sin(theta)) for theta in
                      numpy.arange(0, math.pi/float(2), angle_increment)]
    return curve_filaments + top_part_filaments

def bottom_shell_filaments(count):
    major_radius = 0.92 
    nominal_minor_radius = 0.16
    shell_thickness = 0.004763
    minor_radius = nominal_minor_radius + shell_thickness
    shell_length = math.pi/float(2)*minor_radius
    gap = shell_length/float(count)
    curve_filaments = [(major_radius+minor_radius*math.cos(theta),
                       -minor_radius*math.sin(theta)) for theta in
                       numpy.arange(0, math.pi/float(2), math.pi/(float(2)*count))]
    return curve_filaments

coords = top_shell_filaments(40)+bottom_shell_filaments(40)
x,y = zip(*coords)
plt.plot(x, y, '.')
plt.axis('equal')
savefig(os.getcwd() + '/output/shell_filaments.png')

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

