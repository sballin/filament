#!/usr/bin/env python

from scipy import linalg
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy
import MDSplus
import pickle
import csv
import os.path
import sys
from fields import *
from data_manipulation import *
import tokamak


sensors = tokamak.sensors_unique()

shot = 85140 # was 81077 (vf-only)

# Tree, VF, and OH coil data.
tree = MDSplus.Tree('hbtep2', shot)
(vf_time, vf_signal) = tokamak.get_coil_time_signal(shot, 'VF')
(oh_time, oh_signal) = tokamak.get_coil_time_signal(shot, 'OH')
with open('./resources/coil_R_Z','r') as f:
    ((OHR, OHZ), (VFR, VFZ), (SHR, SHZ)) = pickle.load(f)

if len(sys.argv) > 1:
    edge_gridpoints = int(sys.argv[1])
else:
    sys.exit('Please enter grid side dimension as an argument.')
rows = len(sensors)
columns = edge_gridpoints**2
side_in_m = 0.10 # 0.15 for full view of cross-section
grid_spacing = side_in_m/float(edge_gridpoints)
r_start = 0.92
z_start = -0.5*side_in_m
r_end = r_start+grid_spacing*edge_gridpoints
z_end = z_start+grid_spacing*edge_gridpoints

print 'Grid points:', columns
print 'Grid spacing:', grid_spacing
print 'R bounds:', r_start, r_end
print 'Z bounds:', z_start, z_end
print '--------------------------'

execfile('g_matrix.py')

G_inv = linalg.pinv(G_array, sys.float_info.min)
print 'G_inv dimensions:', G_inv.shape
print 'Matrix rank:', linalg.matrix_rank(G_inv)
print 'Condition number:', linalg.cond(G_inv)
diff=G_array-linalg.pinv(G_inv, sys.float_info.min)
print 'Min, max differences:', numpy.min(diff), numpy.max(diff)
print 'Average value:', numpy.average(G_inv)
if numpy.allclose(linalg.pinv(G_inv, sys.float_info.min), G_array):
    print 'Double inversion returns original matrix.'
else:
    print 'This G matrix might not be amenable to inversion.'
print '--------------------------'

plt.figure() # do not remove with rest of block
plt.subplot(121)
plt.title('G%d' % edge_gridpoints)
plt.imshow(G_array, interpolation='nearest')
plt.ylabel('Sensors (0-%d)' % len(sensors))
plt.xlabel('Current grid points')
plt.colorbar()
plt.subplot(122)
plt.title('pinv(G%d)' % edge_gridpoints)
plt.imshow(G_inv, interpolation='nearest')
plt.colorbar()

savefig(os.getcwd() + '/output/G_matrices/G%d.png' % edge_gridpoints)
print 'Created output/G_matrices/G%d.png' % edge_gridpoints
plt.clf()

signal_dict = tokamak.sensor_signal_dict(shot, sensors, vf_signal, oh_signal, True)

print 'Outputting current images.'
B = [0.0 for x in sensors]
(arbitrary_time, sensor_signal) = tokamak.get_sensor_time_signal(shot, sensors[0].name)
(arbitrary_time, sensor_signal) = clip(arbitrary_time, sensor_signal) # fix length
end_time = 0.010 # when things become boring
(times, ip) = read_data(tree, '.sensors.rogowskis:ip', t_start=0, t_end=end_time)
i = 0
image_index = 0
gs = gridspec.GridSpec(1, 2, height_ratios=[2,1]) # was 2,1
while arbitrary_time[i] < end_time:
    if i % 20 == 0: # frame spacing
        B = [signal_dict[sensor.name][i] for sensor in sensors]
        I_array = numpy.dot(G_inv, B)
        I_array = I_array.reshape(edge_gridpoints, edge_gridpoints)
        I_array *= 1e7 # because (cA)^+ = (A^+)/c and we want just A^+
        if i == 0:
            print 'I:', I_array.shape
#        I_sum = 0
#        for x in numpy.nditer(I_array):
#            I_sum += x
#        plt.text(-3, 1.5, 'Sum of currents: %.2f'%I_sum)

        plt.subplot(211)#211 gs[0]
        plt.suptitle('Current profile for shot %d at %f s' % (shot, arbitrary_time[i]))
        plt.imshow(I_array)
        plt.colorbar()

        plt.subplot(212)#212 gs[1]
        plt.plot(oh_time, oh_signal, label='OH')
        plt.plot(vf_time, vf_signal, label='VF')
        plt.plot(times, ip, label='Ip')
        plt.legend(prop={'size':6})
        plt.xlabel('s')
        plt.ylabel('A')
        plt.xlim(-.001, end_time)
        plt.axvline(arbitrary_time[i], color='r')

        savefig(os.getcwd() + '/output/currents/%05d.png' % image_index)
        image_index += 1
        plt.clf()
        sys.stdout.write("\r%05d.png (%f s)" % (image_index, arbitrary_time[i]))
        sys.stdout.flush()
    i += 1
print

# Make video
os.system('/usr/bin/ffmpeg -i output/currents/%%05d.png output/currents/shot%d_%dgridpoints.mp4' % (shot, edge_gridpoints**2))

print '--------------------------'
print 'Output output/currents/shot%d_%dgridpoints.mp4' % (shot, edge_gridpoints**2)
