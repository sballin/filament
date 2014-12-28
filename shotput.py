#!/usr/bin/env python

from pylab import *
import MDSplus
import pickle
import csv
import os
import numpy

import data_manipulation 
import signals 

#shot_num = 81077

def shot_data(shot_num):
    sensor_file = open('./resources/sensors_unique.csv') # was sensors_fb_p.csv
    sensor_specs = csv.reader(sensor_file, delimiter=',', quotechar='"', 
                              quoting=csv.QUOTE_MINIMAL)
    for i in range(1): # skips first line
        next(sensor_specs)
    
    with open('./resources/sensor_blacklistQian.txt') as f:
        bad_sensors = f.read()
    
    # [Sensor_ID]  [loc_x]  [loc_y]  [loc_z]  [n_x]  [n_y]  [n_z]
    sensors = [signals.Sensor(row[0], row[1], row[2], row[3], row[4], row[5], row[6]) 
               for row in sensor_specs if row[0] not in bad_sensors]
    
    (vf_time, vf_signal) = signals.get_coil_time_signal(shot_num, 'VF')
    # Fix length of integrated data
    (vf_time, vf_signal) = data_manipulation.clip(vf_time, vf_signal) 
    return sensors, vf_time, vf_signal

