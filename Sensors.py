#!/usr/bin/python

import numpy as np #numerical package
import pylab #a matlab-similar code, includes plotting support
import cPickle as pickle # a C-based I/O allowing compression and read/write
from MDSplus import * #standard MDSplus package
import Rogowski_calibration
import os
import sys #adding locations to download packages from
sys.path.append('/home/byrne/byrne_edited_code')
import Integration #code of my own to calculate the Integral of a function
import eyeguides
from collections import OrderedDict
'''This package will output the pickup from any of the HBT-EP magnetic sensor
arrays.  Created on 11/21/13 as an extension of TA.py.'''


def TA(shotnum):
    #allowing multiple time bases only because we sometimes drop samples.
    #if we ever reliably fix this, we can only need one time array, not one for
    #each trace.

    tree = Tree('hbtep2',shotnum)
    sensor_p = OrderedDict()
    sensor_r = OrderedDict()
    time_r = OrderedDict()
    time_p = OrderedDict()

    for i in range(30):
        node = tree.getNode('sensors.magnetic:TA{0:02d}_S{1:d}P'
                            .format((i/3+1),(i%3+1)))
        sensor_p['TA{0:02d}_S{1:d}P'.format((i/3+1),(i%3+1))] = node.data()
        #time_p['TA{0:02d}_S{1:d}P'.format((i/3+1),(i%3+1))] = node.dim_of().data()
        if (i-1)%3 ==0:
           node = tree.getNode('sensors.magnetic:TA{0:02d}_S{1:d}R'
                            .format((i/3+1),(i%3+1)))
           sensor_r['TA{0:02d}_S{1:d}R'.format((i/3+1),(i%3+1))] = node.data()
           #time_r['TA{0:02d}_S{1:d}R'.format((i/3+1),(i%3+1))] = node.dim_of().data() 
    time = node.dim_of().data() 
    #return(time_p,time_r, sensor_p, sensor_r)
    return(time, sensor_p, sensor_r)

def PA(shotnum):
    #allowing multiple time bases only because we sometimes drop samples.
    #if we ever reliably fix this, we can only need one time array, not one for
    #each trace.

    tree = Tree('hbtep2',shotnum)
    PA1_sensor_p = OrderedDict()
    PA1_sensor_r = OrderedDict()
    PA1_time_r = OrderedDict()
    PA1_time_p = OrderedDict()

    PA2_sensor_p = OrderedDict()
    PA2_sensor_r = OrderedDict()
    PA2_time_r = OrderedDict()
    PA2_time_p = OrderedDict()

    for i in range(32):
        key1 = 'PA1_S{0:02d}'.format(i+1)
        key2 = 'PA2_S{0:02d}'.format(i+1)
        name1 = 'sensors.magnetic:'+key1
        name2 = 'sensors.magnetic:'+key2

        node1 = tree.getNode(name1+'P')
        PA1_sensor_p[key1+'P'] = node1.data()
        PA1_time_p[key1+'P'] = node1.dim_of().data()

        node2 = tree.getNode(name2+'P')
        PA2_sensor_p[key2+'P'] = node2.data()
        PA2_time_p[key2+'P'] = node2.dim_of().data()


        if (i)%2 ==0:
           node1 = tree.getNode(name1+'R')
           PA1_sensor_r[key1+'R'] = node1.data()
           PA1_time_r[key1+'R'] = node1.dim_of().data() 

           node2 = tree.getNode(name2+'R')
           PA2_sensor_r[key2+'R'] = node2.data()
           PA2_time_r[key2+'R'] = node2.dim_of().data() 
    time = node2.dim_of().data() 
#    return(PA1_time_p,PA1_time_r, PA1_sensor_p, PA1_sensor_r,
#           PA2_time_p,PA2_time_r, PA2_sensor_p, PA2_sensor_r)
    return(time, PA1_sensor_p, PA1_sensor_r,PA2_sensor_p, PA2_sensor_r)
