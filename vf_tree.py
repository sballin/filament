#!/usr/bin/env python

import math
from scipy import special, constants
from pylab import *
import MDSplus
from flux import *
from data_manipulation import *
import pickle

f = open('./coil_R_Z','r')
((OHR,OHZ),(VFR,VFZ),(SHR,SHZ)) = pickle.load(f)
f.close()

tree = MDSplus.Tree('hbtep2', 69523)

node_magnetic = tree.getNode('sensors.magnetic.FB01_S2P')
sensor_signal = node_magnetic.data()
sensor_time = node_magnetic.dim_of().data()

node_vf = tree.getNode('sensors.VF_CURRENT')
vf_signal = node_vf.data()
vf_time = node_vf.dim_of().data()

(sensor_time, sensor_signal) = clip(sensor_time, sensor_signal)
(vf_time, vf_signal) = clip(vf_time, vf_signal)

flux_vals = []
for i in range(len(vf_signal)):
	flux_vals.append(0)

"""FB01_S2P data."""
pos_x = -0.51261392
pos_y = -0.92477993
sensor_r = sqrt(pos_x*pos_x+pos_y*pos_y)
sensor_z = -0.077077734
sensor_z_component = 0.87206924

"""Calculate flux."""
for coil_index in range(len(VFR)-1):
	sig_index = 0
	for current  in vf_signal:
		flux_vals[sig_index] += flux(current, sensor_r, sensor_z, VFR[coil_index], VFZ[coil_index])	
		sig_index += 1

def plot_magnetic():
	plt.figure(1)
	plt.plot(sensor_time, sensor_signal)
	plt.show()

def plot_vf():
	plt.figure(1)
	plt.plot(vf_time, vf_signal)
	plt.show()

def plot_flux():
	plt.figure(1)
	plt.plot(vf_time, flux_vals, label='calculated flux')
	plt.plot(sensor_time, sensor_signal, label='sensor data')
	plt.legend()
	plt.show()

print '-------------------------'
print len(vf_time)
print len(flux_vals)
print len(VFR)
print len(VFZ)
print '-------------------------'
plot_flux()

