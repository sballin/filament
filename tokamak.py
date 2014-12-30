from pylab import *
import MDSplus
import pickle
import csv
import os
import numpy
import data_manipulation 


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
        self.r    = sqrt(float(x)**2+float(y)**2)
        self.z    = float(z)
        self.n_r  = sqrt(float(n_x)**2+float(n_y)**2)
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
    return (sensor_time, sensor_signal)


def sensors_unique():
    return read_sensor_data('./resources/sensors_unique.csv')


def sensors_PA(number):
    return [sensor for sensor 
            in read_sensor_data('./resources/sensors.csv')
            if 'PA' + str(number) in sensor.name]


def read_sensor_data(filename):
    sensor_file = open(filename) # was sensors_fb_p.csv
    sensor_specs = csv.reader(sensor_file, delimiter=',', quotechar='"', 
                              quoting=csv.QUOTE_MINIMAL)
    for i in range(1): # skips first line
        next(sensor_specs)
    
    with open('./resources/sensor_blacklistQian.txt') as f:
        bad_sensors = f.read()
    
    # [Sensor_ID]  [loc_x]  [loc_y]  [loc_z]  [n_x]  [n_y]  [n_z]
    return [Sensor(row[0], row[1], row[2], row[3], row[4], row[5], row[6]) 
            for row in sensor_specs if row[0] not in bad_sensors]
    

def vf_shot_data(shot_num):
    (vf_time, vf_signal) = get_coil_time_signal(shot_num, 'VF')
    # Fix length of integrated data
    (vf_time, vf_signal) = data_manipulation.clip(vf_time, vf_signal) 
    return vf_time, vf_signal


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
    
