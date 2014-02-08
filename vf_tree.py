#!/usr/bin/env python

import time
from pylab import *
import MDSplus
from flux import *
from data_manipulation import *
import pickle
import flux

startTime = time.time()

"""Get coil location data."""
f = open('./coil_R_Z','r')
((OHR,OHZ),(VFR,VFZ),(SHR,SHZ)) = pickle.load(f)
f.close()

"""FB01_S2P data."""
pos_x = -0.51261392
pos_y = -0.92477993
pos_z = -0.077077734
sens_r = sqrt(pos_x*pos_x+pos_y*pos_y)
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

# """Create empty array for flux values."""
# flux_vals = []
# for i in range(len(vf_signal)):
# 	flux_vals.append(0)

"""Calculate flux."""
flux_vals = flux_val_calc(vf_signal, VFR, VFZ, sens_r, pos_z, sens_len, sens_width, sens_z_norm, coils)
# for i_coil in range(len(VFR)-1):
#     i_sig = 0
#     for current in vf_signal:
#         flux_vals[i_sig] += coils*sens_z_norm*(frac_flux(sens_len, current, sens_r, pos_z, VFR[i_coil], VFZ[i_coil]) - frac_flux(sens_len, current, sens_r-sens_width, pos_z, VFR[i_coil], VFZ[i_coil]))
#         i_sig += 1

"""Calculate average ratio between sensor signal and calculated flux."""
i = 0
ratios_sum = 0
scaling_threshold = 0.001
ratios = []
for f in flux_vals:
	if sensor_signal[i] > scaling_threshold:
		ratio = sensor_signal[i]/f
		ratios.append(ratio)
		ratios_sum += ratio
	i += 1
avg_ratio = ratios_sum/len(ratios)

scaling_vals = []
for i in sensor_time:
	scaling_vals.append(scaling_threshold)

print '------------------------'
print 'Average scaling factor between calculated and measured flux: %s' % avg_ratio
print 'Calculated only using values above %s' % scaling_threshold
print '------------------------'

"""Multiply flux values by inductance."""
for i in range(len(flux_vals)):
	flux_vals[i] *= avg_ratio

def plot_flux():
	plt.figure(1)
	plt.plot(vf_time, flux_vals, label='calculated flux')
	plt.plot(sensor_time, sensor_signal, label='sensor data')
	plt.plot(sensor_time, scaling_vals, label='scaling threshold')
	plt.legend()
	plt.show()

endTime = time.time()
elapsed = endTime-startTime
print '%.2f seconds elapsed.' % elapsed

plot_flux()
