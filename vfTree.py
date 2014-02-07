#!/usr/bin/env python

import math
from scipy import special, constants
from pylab import *
import MDSplus
from hbtFlux import *
from dataManipulation import *
import pickle

f = open('./coil_R_Z','r')
((OHR,OHZ),(VFR,VFZ),(SHR,SHZ)) = pickle.load(f)
f.close()

tree = MDSplus.Tree('hbtep2',69523)

node_magnetic = tree.getNode('sensors.magnetic.FB01_S1P')
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

# DUMMY VALUES
SENSR = .5
SENSZ = .5 
for coil_index in range(len(VFR)-1):
	sig_index = 0
	for current  in vf_signal:
		flux_vals[sig_index] += flux(current, SENSR, SENSZ, VFR[coil_index], VFZ[coil_index])	
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
	plt.plot(vf_time, flux_vals)
	plt.plot(sensor_time, sensor_signal)
	plt.show()

print '-------------------------'
print len(vf_time)
print len(flux_vals)
print len(VFR)
print len(VFZ)
print '-------------------------'
plot_flux()

