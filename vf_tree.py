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

"""FB01_S2P specs."""
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

"""Calculate flux."""
# field_vals = sensor_flux_calc(vf_signal, VFR, VFZ, sens_r, pos_z, sens_len, sens_width, sens_z_norm, coils)
theta = acos(0.87206924)
field_vals = green_times_current(vf_signal, VFR, VFZ, sens_r, pos_z, theta)

"""Calculate average ratio between sensor signal and calculated flux."""
i = 0
ratios_sum = 0
scaling_threshold = 0.001
ratios = []
for f in field_vals:
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
print 'Measured field = %.6s * calculated field' % avg_ratio
print 'Calculated only using values above %s' % scaling_threshold
print '------------------------'

"""Multiply field values by inductance."""
for i in range(len(field_vals)):
	field_vals[i] *= -1

def plot_flux():
	plt.figure(1)

	plt.subplot(211)
	plt.plot(vf_time, field_vals, label='calculated field')
	plt.plot(sensor_time, sensor_signal, label='sensor data')
	plt.plot(sensor_time, scaling_vals, label='scaling threshold')
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

endTime = time.time()
elapsed = endTime-startTime
print '%.2f seconds elapsed.' % elapsed

print '%f' % green_base()

plot_flux()
