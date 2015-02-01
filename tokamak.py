import matplotlib.pyplot as plt
import math
import MDSplus
import pickle
import csv
import sys
import os
import numpy
import data_manipulation 
import fields


#shot_num = 81077

conductSS = 1/72.0e-8
mu_0 = 4.0e-7*math.pi
asize = 0.005

# Get coil location data
with open('./resources/coil_R_Z', 'r') as f:
    ((OHR, OHZ), (VFR, VFZ), (SHR, SHZ)) = pickle.load(f)


class Sensor:
    def __init__(self, name, x, y, z, n_x, n_y, n_z):
        self.name = name
        self.x    = float(x)
        self.y    = float(y)
        self.z    = float(z)
        self.r    = math.sqrt(float(x)**2+float(y)**2)
        self.n_r  = math.sqrt(float(n_x)**2+float(n_y)**2)
        self.n_z  = float(n_z)


def get_coil_time_signal(shot, coil):
    tree = MDSplus.Tree('hbtep2', shot)
    node = tree.getNode('sensors.'+coil+'_CURRENT')
    (time, signal) = (node.dim_of().data(), node.data())
    return (time, signal)


def get_sensor_time_signal(shot, sensor):
    tree = MDSplus.Tree('hbtep2', shot)
    node_magnetic = tree.getNode('sensors.magnetic.' + sensor)
    (sensor_time, sensor_signal) = (node_magnetic.dim_of().data(), node_magnetic.data())
    # Fix length of integrated data
    (sensor_time, sensor_signal) = data_manipulation.clip(sensor_time, sensor_signal) 
    return (sensor_time, sensor_signal)


def sensors_unique():
    return read_sensor_data('./resources/sensors_unique.csv')


def sensors_PA(number):
    return [sensor for sensor 
            in read_sensor_data('./resources/sensors.csv')
            if 'PA' + str(number) in sensor.name]


def sensors_PA():
    return [sensor for sensor 
            in read_sensor_data('./resources/sensors.csv') if 'PA' in sensor.name]


def sensors_TA():
    return [sensor for sensor in read_sensor_data('./resources/sensors.csv')
            if 'TA' in sensor.name]
    

def sensors_FB():
    return [sensor for sensor in read_sensor_data('./resources/sensors.csv')
            if 'FB' in sensor.name]


def sensors_all():
    return read_sensor_data('./resources/sensors.csv')


def read_sensor_data(filename):
    sensor_file = open(filename) # was sensors_fb_p.csv
    sensor_specs = csv.reader(sensor_file, delimiter=',', quotechar='"', 
                              quoting=csv.QUOTE_MINIMAL)
    for i in range(1): # skips first line
        next(sensor_specs)
    
    bad_sensors = blacklist_sean()
    
    # [Sensor_ID]  [loc_x]  [loc_y]  [loc_z]  [n_x]  [n_y]  [n_z]
    return [Sensor(row[0], row[1], row[2], row[3], row[4], row[5], row[6]) 
            for row in sensor_specs if row[0] not in bad_sensors]


def blacklist_quian():
    with open('./resources/sensor_blacklistQian.txt') as f:
        return f.read()


def blacklist_sean():
    '''Based on data from certain VF-only shots, including 85140 and 81077'''
    return {'TA01_S2R', 'TA02_S2R', 'TA03_S2R', 'TA04_S2R', 'TA05_S2R', 'TA06_S2R', 
            'TA07_S2R', 'TA08_S2R', 'TA09_S2R', 'TA10_S2R', 'FB02_S1P', 'FB05_S1P',
            'FB06_S2P', 'FB03_S4R', 'FB08_S3R', 'PA1_S09P', 'PA1_S25P', 'PA2_S08P', 
            'PA2_S09P', 'PA2_S25P', 'PA1_S29R'}


def vf_data(shot_num):
    (vf_time, vf_signal) = get_coil_time_signal(shot_num, 'VF')
    # Fix length of integrated data
    (vf_time, vf_signal) = data_manipulation.clip(vf_time, vf_signal) 
    return vf_time, vf_signal


def sensor_signal_dict(shot, sensors, vf_signal, oh_signal, subtract):
    if subtract:
        outfile  = 'output/signals%dsubtracted.p' % shot
    else:
        outfile  = 'output/signals%d.p' % shot

    if os.path.isfile(outfile):
        signal_dict = pickle.load(open(outfile, 'rb'))
        print 'DANGER: loaded %s from disk.' % outfile 
    else:
        print 'Indexing and cleaning up sensor signals.'
        signal_dict = dict()
        for i, sensor in enumerate(sensors):
            (sensor_time, sensor_signal) = get_sensor_time_signal(shot, sensor.name)
            # Fix length of integrated data
            (sensor_time, sensor_signal) = data_manipulation.clip(sensor_time, sensor_signal) 
            if subtract:
                coils_signal = fields.B_VF(vf_signal, VFR, VFZ, sensor)[:-1] \
                               + fields.OH_field(oh_signal, OHR, OHZ, sensor)
                sensor_signal -= coils_signal[:len(sensor_signal)]
            signal_dict[sensor.name] = sensor_signal
            progress = i/float(len(sensors))*100
            sys.stdout.write('\r%.2f%%' % progress)
            sys.stdout.flush()
        print '\nSensors used:', len(sensors)
        pickle.dump(signal_dict, open(outfile, 'wb'))
        print 'Created ' + outfile
    return signal_dict
    

def top_shell_filaments(count):
    '''Return coordinates of evenly spaced filaments along top shell'''
    major_radius = 0.92 
    nominal_minor_radius = 0.16
    shell_thickness = 0.004763
    minor_radius = nominal_minor_radius + shell_thickness
    shell_length = math.pi/2.*minor_radius
    top_part = 0.0635
    gap = (top_part + shell_length)/float(count)
    angle_increment = gap/shell_length*math.pi/2.

    curve_filaments = [(major_radius+minor_radius*math.cos(theta),
                       minor_radius*math.sin(theta)) for theta in
                       numpy.arange(0, math.pi/2., angle_increment)]
    top_part_filaments = [(i, minor_radius) for i in 
                          numpy.arange(major_radius-0.5*gap, major_radius-top_part, -gap)]
    all_filaments = (curve_filaments + top_part_filaments)#[:count]
    all_filaments.reverse()
    return all_filaments 


def bot_shell_filaments_mirror(count):
    return [(x,-z) for (x,z) in top_shell_filaments(count)]


def bot_shell_filaments(count):
    '''Return coordinates of evenly spaced filaments along bottom shell'''
    major_radius = 0.92 
    nominal_minor_radius = 0.16
    shell_thickness = 0.004763
    minor_radius = nominal_minor_radius + shell_thickness
    shell_length = math.pi/2.*minor_radius
    gap = shell_length/float(count)

    curve_filaments = [(major_radius+minor_radius*math.cos(theta),
                       minor_radius*math.sin(theta)) for theta in
                       numpy.arange(-math.pi/2., 0, math.pi/2.*count)]
    return curve_filaments[:count]


def plot_geometry(filaments):
    (x,y) = zip(*filaments)
    plt.plot(x, y, '.')
    #plt.plot(VFR, VFZ, 'x')
    plt.axis('equal')
    
