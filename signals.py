from math import *
import MDSplus


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

