#!/usr/bin/env python

from scipy import linalg # or scipy?
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
from pylab import *
import numpy
import MDSplus
import pickle
import csv
import os.path
import sys
from flux import *
from data_manipulation import *

def get_coil_time_signal(shot, coil):
    tree = MDSplus.Tree('hbtep2', shot)
    node = tree.getNode('sensors.'+coil+'_CURRENT')
    (time, signal) = (node.dim_of().data(), node.data())
    return (time, signal)

def get_sensor_time_signal(shot, sensor):
    tree = MDSplus.Tree('hbtep2', shot)
    node_magnetic = tree.getNode('sensors.magnetic.' + sensor)
    (sensor_time, sensor_signal) = (node_magnetic.dim_of().data(), node_magnetic.data())
    return (sensor_time, sensor_signal)

current_dir = os.getcwd()

sensor_file = open('./sensors_fb_p.csv')
blacklist = ['TA01_S2R', 'TA02_S2R', 'TA03_S2R', 'TA04_S2R', 'TA05_S2R', 'TA06_S2R', 'TA07_S2R', 'TA08_S2R', 'TA09_S2R', 'TA10_S2R', 'FB02_S1P', 'FB05_S1P', 'FB06_S2P', 'FB03_S4R', 'FB08_S3R', 'PA1_S09P', 'PA1_S25P', 'PA2_S08P', 'PA2_S09P', 'PA2_S25P', 'PA1_S29R']
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

# [Sensor_ID]  [loc_x]  [loc_y]  [loc_z]  [n_x]  [n_y]  [n_z]
for row in sensor_specs:
    if row[0] not in blacklist:
        sensor_names.append(row[0])
        pos_x.append(float(row[1]))
        pos_y.append(float(row[2]))
        pos_z.append(float(row[3]))
        n_x.append(float(row[4]))
        n_y.append(float(row[5]))
        n_z.append(float(row[6]))

shot = 85140 # 85140 # used 81077 vf-only

tree = MDSplus.Tree('hbtep2', shot)
# VF and OH coil data.
(vf_time, vf_signal) = get_coil_time_signal(shot, 'VF')
(oh_time, oh_signal) = get_coil_time_signal(shot, 'OH')
with open('./coil_R_Z','r') as f:
    ((OHR,OHZ),(VFR,VFZ),(SHR,SHZ)) = pickle.load(f)

rows = len(sensor_names)
edge_gridpoints = int(sys.argv[1])
columns = edge_gridpoints**2
side_in_m = 0.10
grid_spacing = side_in_m/float(edge_gridpoints) # 0.15 for full view of cross-section
r_start = 0.92
z_start = -0.5*side_in_m
r_end = r_start+grid_spacing*edge_gridpoints
z_end = z_start+grid_spacing*edge_gridpoints

G = [[0.0 for x in xrange(columns)] for x in xrange(rows)]

print 'Grid points:', columns
print 'Grid spacing:', grid_spacing
print 'R bounds:', r_start, r_end
print 'Z bounds:', z_start, z_end
print '--------------------------'

# Make/get G matrix
if not os.path.isfile('G%d.p'%edge_gridpoints):
    print 'Making G%d.p'%edge_gridpoints
    for i in range(rows):
        pos_r = sqrt(pos_x[i]**2+pos_y[i]**2)
        n_r = sqrt(n_x[i]**2+n_y[i]**2)
        for w in range(edge_gridpoints):
            for h in range(edge_gridpoints):
                r = r_start+w*grid_spacing
                z = z_start+h*grid_spacing
                G[i][h*edge_gridpoints+w] = 1e7*greens_function(r, z-pos_z[i], pos_r, n_r, n_z[i])
        progress = i/float(rows)*100
        sys.stdout.write('\r%.2f%%' %progress)
        sys.stdout.flush()
    G_array = numpy.array(G)
    pickle.dump(G_array, open('G%d.p'%edge_gridpoints, 'wb'))
    print
else:
    G_array = pickle.load(open('G%d.p'%edge_gridpoints, 'rb'))
    print 'WARNING: loaded G%d.p'%edge_gridpoints
#print 'Min value:', numpy.min(G_array)
print 'Average value:', numpy.average(G_array)
print 'Matrix rank:', linalg.matrix_rank(G_array)

G_inv = linalg.pinv(G_array, sys.float_info.min)
print 'G_inv dimensions:', G_inv.shape
#print 'Min value:', numpy.min(G_inv)
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
plt.title('G%d'%edge_gridpoints)
plt.imshow(G_array, interpolation='nearest')
plt.ylabel('Sensors (0-%d)'%len(sensor_names))
plt.xlabel('Current grid points')
plt.colorbar()
plt.subplot(122)
plt.title('pinv(G%d)'%edge_gridpoints)
plt.imshow(G_inv, interpolation='nearest')
plt.colorbar()
savefig(current_dir + '/Gs/G%d.png'%edge_gridpoints)
plt.clf()

if not os.path.isfile('signals%d.p'%shot):
    print 'Indexing and cleaning up sensor signals.'
    signal_dict = dict()
    for i,name in enumerate(sensor_names):
        (sensor_time, sensor_signal) = get_sensor_time_signal(shot, name)
        (sensor_time, sensor_signal) = clip(sensor_time, sensor_signal) # fix length
        sens_r = sqrt(pos_x[i]**2+pos_y[i]**2)
        n_r = sqrt(n_x[i]**2+n_y[i]**2)
        coils_signal = B_signal(vf_signal, VFR, VFZ, sens_r, pos_z[i], n_r, n_z[i])+OH_field(oh_signal, OHR, OHZ, sens_r, pos_z[i], n_r, n_z[i])
        for j in range(len(sensor_signal)):
            sensor_signal[j] -= coils_signal[j]
        signal_dict[sensor_names[i]] = sensor_signal
    print 'Sensors used:', len(sensor_names)
    pickle.dump(signal_dict, open('signals%d.p'%shot, 'wb'))
else:
    signal_dict = pickle.load(open('signals%d.p'%shot, 'rb'))
    print 'WARNING: loaded signals%d.p'%shot
print 'Outputting current images.'
B = [0.0 for x in xrange(len(sensor_names))]
(arbitrary_time, sensor_signal) = get_sensor_time_signal(shot, sensor_names[0])
(arbitrary_time, sensor_signal) = clip(arbitrary_time, sensor_signal) # fix length
end_time = 0.010 # when things become boring
(times, ip) = read_data(tree, '.sensors.rogowskis:ip', t_start=0, t_end=end_time)
i = 0
image_index = 0
gs = gridspec.GridSpec(2, 1,height_ratios=[5,1])
while arbitrary_time[i] < end_time:
    if i%20 == 0: # frame spacing
        for j in range(len(sensor_names)):
            signal = signal_dict[sensor_names[j]][i]
            B[j] = signal_dict[sensor_names[j]][i]
        I_array = numpy.dot(G_inv, B)
        I_array = I_array.reshape(edge_gridpoints, edge_gridpoints)
        I_array *= 1e7 # because (cA)^+ = (A^+)/c and we want just A^+
        if i == 0:
            print 'I:', I_array.shape
        plt.subplot(gs[0])#211
        plt.suptitle('Current profile for shot %d at %f s'%(shot, arbitrary_time[i]))
        I_sum = 0
        for x in numpy.nditer(I_array):
            I_sum += x
#        plt.text(-3, 1.5, 'Sum of currents: %.2f'%I_sum)
        plt.imshow(I_array)
        plt.colorbar()
        plt.subplot(gs[1])#212
        plt.plot(times, ip)
        plt.xlabel('s')
        plt.ylabel('Ip (A)')
        plt.xlim(-.001, end_time)
        plt.axvline(arbitrary_time[i], color='r')
        savefig(current_dir + '/currents/%05d.png'%image_index)
        image_index += 1
        plt.clf()
        sys.stdout.write("\r%05d.png (%f s)"%(image_index, arbitrary_time[i]))
        sys.stdout.flush()
    i += 1
print

# Make video
os.system('/usr/bin/ffmpeg -i currents/%%05d.png currents/shot%d_%dgridpoints.mp4'%(shot, edge_gridpoints**2))

print '--------------------------'
print 'Output currents/shot%d_%dgridpoints.mp4'%(shot, edge_gridpoints**2)
