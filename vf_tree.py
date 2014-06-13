#!/usr/bin/env python

import time
import matplotlib
matplotlib.use('Agg')   # headless
from pylab import *
import MDSplus
from flux import *
from data_manipulation import *
import pickle
import flux
import csv
import os
from scipy import special, constants, integrate


def get_vf_time_signal(shot):
    tree = MDSplus.Tree('hbtep2', shot)
    node_vf = tree.getNode('sensors.VF_CURRENT')
    (vf_time, vf_signal) = (node_vf.dim_of().data(), node_vf.data())
    return (vf_time, vf_signal)

def get_sensor_time_signal(shot, sensor):
    tree = MDSplus.Tree('hbtep2', shot)
    node_magnetic = tree.getNode('sensors.magnetic.' + sensor)
    (sensor_time, sensor_signal) = (node_magnetic.dim_of().data(), node_magnetic.data())
    return (sensor_time, sensor_signal)

def get_avg_scaling(field_vals, sensor_signal, scaling_threshold):
    i = 0
    ratios_sum = 0
    ratios = []
    for f in field_vals:
        if sensor_signal[i] > scaling_threshold:
            ratio = sensor_signal[i]/f
            ratios.append(ratio)
            ratios_sum += ratio
        i += 1
    try:
        avg_ratio = ratios_sum/len(ratios)
    except ZeroDivisionError:
        return 999999
    return avg_ratio

def plot_comparison(sensor, sensor_time, sensor_signal, vf_time, field_vals, scaling_threshold, avg_scaling, shot):
    plt.figure()
    plt.plot(vf_time, field_vals, label='calculated field', color='blue')
    plt.plot(sensor_time, sensor_signal, label='sensor data', color='red')
    plt.plot([0,.03], [scaling_threshold, scaling_threshold], label='scaling threshold', color='gray')
    plt.xlim([0,.03])
    #plt.ylim([-0.001,.012])
    plt.ylabel('Tesla')
    plt.xlabel('s')
    title = sensor + ", scaling value " + str(avg_scaling) + ", shot no. " + str(shot)
    plt.suptitle(title)
    plt.legend()
    current_dir = os.getcwd()
    savefig(current_dir + '/FB_images_integrate/' + sensor + '.png')
    plt.clf()

# Get coil location data.
with open('./coil_R_Z','r') as f:
    ((OHR,OHZ),(VFR,VFZ),(SHR,SHZ)) = pickle.load(f)

sensor_file = open('./sensors_fb_p.csv')
sensor_specs = csv.reader(sensor_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
for i in range(1): # skips first line
    next(sensor_specs)

sensor_names = []
pos_x = []
pos_y = []
pos_z = []
n_x = []
n_y = []
n_z = []

with open('./sensor_blacklistQian.txt') as f:
    bad_sensors = f.read()

# [Sensor_ID]  [loc_x]  [loc_y]  [loc_z]  [n_x]  [n_y]  [n_z]
for row in sensor_specs:
    if row[0] not in bad_sensors:
        sensor_names.append(row[0])
        pos_x.append(float(row[1]))
        pos_y.append(float(row[2]))
        pos_z.append(float(row[3]))
        n_x.append(float(row[4]))
        n_y.append(float(row[5]))
        n_z.append(float(row[6]))

shot = 81077
scaling_threshold = 0.0005
print 'Calculated only using values above %s' % scaling_threshold
(vf_time, vf_signal) = get_vf_time_signal(shot)
(vf_time, vf_signal) = clip(vf_time, vf_signal) # fix length of integrated data

scaling_vals = [0, ] #
avg_scaling_vals = []
bad_calculation_sensors = []

for i in range(len(sensor_names)):
    startTime = time.time()
    sensor = sensor_names[i]
    (sensor_time, sensor_signal) = get_sensor_time_signal(shot, sensor)
    (sensor_time, sensor_signal) = clip(sensor_time, sensor_signal) # fix length of integrated data
    # Prep specs.
    pos_r = sqrt(pos_x[i]**2+pos_y[i]**2)
    n_r = sqrt(n_x[i]**2+n_y[i]**2)
    # Calculate field.
    (field_vals, sensor_signal) = B_signal(vf_signal, sensor_signal, VFR, VFZ, pos_r, pos_z[i], n_r, n_z[i])
    # Text output.
    avg_scaling = get_avg_scaling(field_vals, sensor_signal, scaling_threshold)
    print '------------------------------------------' + sensor + ':'
    print 'Measured field = %.6s * calculated field' % avg_scaling
    if avg_scaling < 100:
        avg_scaling_vals.append(avg_scaling)
    if avg_scaling > 3:
        bad_calculation_sensors.append(sensor)

    # Image output.
    plot_comparison(sensor, sensor_time, sensor_signal, vf_time, field_vals, scaling_threshold, avg_scaling, shot)
    # Report time elapsed per sensor.
    endTime = time.time()
    elapsed = endTime-startTime
    print '%.2f seconds elapsed.' % elapsed

avg_avg = sum(avg_scaling_vals)/float(len(avg_scaling_vals))
print 'Average of average scaling values: %.6s' % avg_avg
print 'bad calculations for these sensors:', bad_calculation_sensors

def plot_all():
    plt.figure(1)

    plt.subplot(211)
    plt.plot(vf_time, field_vals, label='calculated field')
    plt.plot(sensor_time, sensor_signal, label='sensor data')
    # plt.plot(sensor_time, scaling_vals, label='scaling threshold')
    plt.xlim([0,.06])
    plt.ylabel('Tesla')
    plt.legend()

    plt.subplot(212)
    plt.plot(vf_time, vf_signal, label='vf current')
    plt.xlim([0,.06])
    plt.ylabel('A')
    plt.xlabel('s')
    plt.legend()

    plt.show()

#plot_flux()
