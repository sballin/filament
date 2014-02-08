#!/usr/bin/env python

from pylab import *
import MDSplus
from flux import *
from data_manipulation import *
import pickle
import flux

"""Get coil location data."""
f = open('./coil_R_Z','r')
((OHR,OHZ),(VFR,VFZ),(SHR,SHZ)) = pickle.load(f)
f.close()

"""FB01_S2P data."""
pos_x = -0.51261392
pos_y = -0.92477993
pos_z = -0.077077734
sensor_r = sqrt(pos_x*pos_x+pos_y*pos_y)
sens_z_norm = 0.87206924
coils = 20.0
sens_len = .1
sens_width = .004

"""Get shot data from tree."""
tree = MDSplus.Tree('hbtep2', 69523)
node_magnetic = tree.getNode('sensors.magnetic.FB01_S2P')
sensor_signal = node_magnetic.data()
sensor_time = node_magnetic.dim_of().data()
node_vf = tree.getNode('sensors.VF_CURRENT')
vf_signal = node_vf.data()
vf_time = node_vf.dim_of().data()

"""Fix length of integrated data."""
(sensor_time, sensor_signal) = clip(sensor_time, sensor_signal)
(vf_time, vf_signal) = clip(vf_time, vf_signal)

"""Create empty array for flux values."""
flux_vals = []
for i in range(len(vf_signal)):
	flux_vals.append(0)

"""Calculate flux."""
for coil_index in range(len(VFR)-1):
	sig_index = 0
	for current in vf_signal:
		flux_vals[sig_index] += coils*sens_z_norm*(frac_flux(sens_len, current, sensor_r, pos_z, VFR[coil_index], VFZ[coil_index]) - frac_flux(sens_len, current, sensor_r-sens_width, pos_z, VFR[coil_index], VFZ[coil_index]))
		sig_index += 1

"""Calculate average ratio between sensor signal and calculated flux."""
i = 0
sum_ratios = 0
ratios = []
for f in flux_vals:
	if sensor_signal[i] > 0.001:
		ratio = sensor_signal[i]/f
		ratios.append(ratio)
		sum_ratios += ratio
	i += 1
avg_ratio = sum_ratios/len(ratios)
print 'Inductance: %s' % avg_ratio

"""Multiply flux values by inductance."""
for i in range(len(flux_vals)):
	flux_vals[i] *= avg_ratio

def plot_flux():
	plt.figure(1)
	plt.plot(vf_time, flux_vals, label='calculated flux')
	plt.plot(sensor_time, sensor_signal, label='sensor data')
	plt.legend()
	plt.show()

plot_flux()
