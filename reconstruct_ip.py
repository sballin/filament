#!/usr/bin/env python

#from scipy import linalg
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy
import numpy.linalg as linalg
import MDSplus
import pickle
import csv
import os.path
import sys
import fields 
from data_manipulation import *
import tokamak
import eigenmodes


sensors = tokamak.sensors_unique()

shot = 85140 # was 81077 (vf-only)

# Tree, VF, and OH coil data.
tree = MDSplus.Tree('hbtep2', shot)
(vf_time, vf_signal) = tokamak.get_coil_time_signal(shot, 'VF')
(oh_time, oh_signal) = tokamak.get_coil_time_signal(shot, 'OH')

if len(sys.argv) > 1:
    edge_gridpoints = int(sys.argv[1])
else:
    sys.exit('Please enter grid side dimension as an argument.')
rows = len(sensors)
columns = edge_gridpoints**2
side_in_m = 0.19 # 0.15 for full view of cross-section
grid_spacing = side_in_m/float(edge_gridpoints)
r_start = 0.85
z_start = -0.40*side_in_m
r_end = r_start+grid_spacing*edge_gridpoints
z_end = z_start+grid_spacing*edge_gridpoints
print 'Grid points:', columns
print 'Grid spacing:', grid_spacing
print 'R bounds:', r_start, r_end
print 'Z bounds:', z_start, z_end
print '--------------------------'

# Load/create G matrix
if os.path.isfile('output/G_matrices/G%d.p' % edge_gridpoints):
    G_array = pickle.load(open('output/G_matrices/G%d.p' % edge_gridpoints, 'rb'))
    print 'DANGER: loaded G%d.p from disk.' % edge_gridpoints
else:
    print 'Making G%d.p' % edge_gridpoints
    G = [[0.0 for x in xrange(columns)] for x in xrange(rows)]
    for i,sensor in enumerate(sensors):
        for w in xrange(edge_gridpoints):
            for h in xrange(edge_gridpoints):
                r = r_start+w*grid_spacing
                z = z_start+h*grid_spacing
                G[i][h*edge_gridpoints+w] = 1e7*fields.greens_function(r, z-sensor.z, sensor.r, sensor.n_r, sensor.n_z) # 1e7 multiplication for matrix stuff
        progress = i/float(rows)*100
        sys.stdout.write('\r%.2f%%' %progress)
        sys.stdout.flush()
    G_array = numpy.array(G)
    pickle.dump(G_array, open('output/G_matrices/G%d.p'%edge_gridpoints, 'wb'))
    print
print 'Average value:', numpy.average(G_array)
print 'Matrix rank:', linalg.matrix_rank(G_array)

# Verify G matrix inversion
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

# Create comparison image for G and G_inv 
plt.figure() # do not remove with rest of block
plt.subplot(223)
plt.title('G%d' % edge_gridpoints)
plt.imshow(G_array, interpolation='nearest')
plt.ylabel('Sensors (0-%d)' % len(sensors))
plt.xlabel('Current grid points')
plt.colorbar()
plt.subplot(224)
plt.title('pinv(G%d)' % edge_gridpoints)
plt.imshow(G_inv, interpolation='nearest')
plt.colorbar()
plt.subplot(221)
plt.title('Gridpoints')
rs = []
zs = []
filaments = eigenmodes.ss_filaments(60)
for w in xrange(edge_gridpoints):
    for h in xrange(edge_gridpoints):
        rs.append(r_start+w*grid_spacing)
        zs.append(z_start+h*grid_spacing)
plt.plot(rs, zs, 'ro')
plt.plot([f.r for f in filaments], [f.z for f in filaments], 'o')
plt.gcf().gca().add_artist(plt.Circle((0.92, 0), .15, color='turquoise'))
plt.xlim([0.67, 1.17])
plt.ylim([-0.19, 0.19])
ax = plt.subplot(222, projection='3d')
plt.title('Sensors used')
ax.scatter([s.x for s in sensors], [s.y for s in sensors], [s.z for s in sensors], c='red', marker='o')
plt.savefig(os.getcwd() + '/output/G_matrices/G%d.png' % edge_gridpoints)
plt.clf()
print 'Created output/G_matrices/G%d.png' % edge_gridpoints

signal_dict = tokamak.sensor_signal_dict(shot, sensors, vf_signal, oh_signal, True)

# Write image sequence
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

        plt.suptitle('Current profile for shot %d at %f s' % (shot, arbitrary_time[i]))
        plt.subplot(222)
        plt.title('Plasma current (A)')
        plt.imshow(I_array)
        plt.clim(-5e7, 5e7)
        plt.colorbar()

        plt.subplot(221)
        plt.title('Coil currents and total Ip')
        plt.plot(oh_time, oh_signal, label='OH')
        plt.plot(vf_time, vf_signal, label='VF')
        plt.plot(times, ip, label='Ip')
        plt.legend(prop={'size':6})
        plt.xlabel('s')
        plt.ylabel('A')
        plt.xlim(-.001, end_time)
        plt.axvline(arbitrary_time[i], color='r')

        ax = plt.subplot(223, projection='3d')
        plt.title('Sensors used')
        ax.scatter([s.x for s in sensors], [s.y for s in sensors], [s.z for s in sensors], c='red', marker='o')
        #ax.view_init(30, angle)
        #ax.axis('off')
        plt.axis('equal')

        plt.subplot(224)
        plt.title('Gridpoints')
        plt.plot(rs, zs, 'ro')
        plt.plot([f.r for f in filaments], [f.z for f in filaments], 'o')
        plt.gcf().gca().add_artist(plt.Circle((0.92, 0), .15, color='turquoise'))
        plt.xlim([0.67, 1.17])
        plt.ylim([-0.19, 0.19])
        plt.xlabel('m')
        plt.savefig(os.getcwd() + '/output/currents/%05d.png' % image_index)
        image_index += 1
        plt.clf()
        sys.stdout.write("\r%05d.png (%f s)" % (image_index, arbitrary_time[i]))
        sys.stdout.flush()
    i += 1
print
print '--------------------------'

# Make video
os.system('/usr/bin/ffmpeg -i output/currents/%%05d.png output/currents/shot%d_%dgridpoints.mp4' % (shot, edge_gridpoints**2))
print 'Created output/currents/shot%d_%dgridpoints.mp4' % (shot, edge_gridpoints**2)
