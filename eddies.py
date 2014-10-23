#!/usr/bin/env python

import time
import matplotlib
matplotlib.use('Agg')   # headless
from pylab import *
import MDSplus
import pickle
import csv
import os
import numpy
from fields import *
from data_manipulation import *
from signals import *


def plot_comparison(sensor, sensor_time, sensor_signal, vf_time, field_vals, scaling_threshold, avg_scaling, shot):
    difference = [field_vals[i]-sensor_signal[i] for i in range(len(sensor_time))]
    plt.figure()
    plt.plot(vf_time, field_vals, label='calculated field', color='blue')
    plt.plot(sensor_time, sensor_signal, label='sensor data', color='red')
    plt.plot(sensor_time, difference, label='difference', color='green')
    plt.ylabel('Tesla')
    plt.xlabel('s')
    title = sensor + ", shot no. " + str(shot)
    plt.suptitle(title)
    plt.legend()
    plt.xlim([.0015, .03])
    #plt.ylim([-0.002,.015])
    current_dir = os.getcwd()
    savefig(current_dir + '/output/eddy_images/' + sensor + '.png')
    plt.clf()


# Get coil location data.
with open('./coil_R_Z','r') as f:
    ((OHR,OHZ),(VFR,VFZ),(SHR,SHZ)) = pickle.load(f)

sensor_file = open('./sensors_unique.csv') # was sensors_fb_p.csv
sensor_specs = csv.reader(sensor_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
for i in range(1): # skips first line
    next(sensor_specs)

with open('./sensor_blacklistQian.txt') as f:
    bad_sensors = f.read()

# [Sensor_ID]  [loc_x]  [loc_y]  [loc_z]  [n_x]  [n_y]  [n_z]
sensors = [Sensor(row[0], row[1], row[2], row[3], row[4], row[5], row[6]) for row in sensor_specs if row[0] not in bad_sensors];

shot = 81077
scaling_threshold = 0.0005
print 'Calculated only using values above %s' % scaling_threshold
(vf_time, vf_signal) = get_coil_time_signal(shot, 'VF')
(vf_time, vf_signal) = clip(vf_time, vf_signal) # fix length of integrated data

scaling_vals = []
avg_scaling_vals = []
bad_calculation_sensors = []

for i,sensor in enumerate(sensors):
    startTime = time.time()
    (sensor_time, sensor_signal) = get_sensor_time_signal(shot, sensor.name)
    (sensor_time, sensor_signal) = clip(sensor_time, sensor_signal) # fix length of integrated data
    if 'TA' in sensor.name: # radial TA sensors point inwards
        if 'R' in sensor.name:
            sensor.n_r = -1.0
    # Calculate field.
    field_vals = B_signal(vf_signal, VFR, VFZ, sensor.r, sensor.z, sensor.n_r, sensor.n_z)
    # Text output.
    avg_scaling = get_avg_scaling(field_vals, sensor_signal, scaling_threshold)
    print '------------------------------------------' + sensor.name + ':'
    print 'Measured field = %.6s * calculated field' % avg_scaling
    if avg_scaling < 100:
        avg_scaling_vals.append(avg_scaling)
    if avg_scaling > 3:
        bad_calculation_sensors.append(sensor.name)
    # Image output.
    plot_comparison(sensor.name, sensor_time, sensor_signal, vf_time, field_vals, scaling_threshold, avg_scaling, shot)
    # Report time elapsed per sensor.
    endTime = time.time()
    elapsed = endTime-startTime
    print '%.2f seconds elapsed.' % elapsed

avg_avg = sum(avg_scaling_vals)/float(len(avg_scaling_vals))
print 'Average of average scaling values: %.6s' % avg_avg
print 'bad calculations for these sensors:', bad_calculation_sensors

