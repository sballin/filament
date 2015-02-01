#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')   # headless
import matplotlib.pyplot as plt
import os
import fields
import data_manipulation 
import eigenmodes
import tokamak
import numpy as np


shot = 81077
vf_time, vf_signal = tokamak.vf_data(shot)
sensor_groups = []
TA_suffixes = ['S1P', 'S2P', 'S3P', 'S4P']
sensor_groups += [[s for s in tokamak.sensors_TA() if n in s.name] for n in TA_suffixes]
FB_suffixes = ['S1P', 'S2P', 'S3P', 'S4P', 'S1R', 'S2R', 'S3R', 'S4R', 'S4R']
sensor_groups += [[s for s in tokamak.sensors_FB() if n in s.name] for n in FB_suffixes]
PA_suffixes = ['S01R', 'S02R', 'S03R', 'S04R', 'S05R', 'S06R', 
'S07R', 'S08R', 'S09R', 'S10R', 'S11R', 'S12R', 'S13R', 'S14R', 'S15R', 'S16R', 'S17R', 'S18R',
'S19R', 'S20R', 'S21R', 'S01P', 'S02P', 'S03P', 'S04P', 'S05P', 'S06P', 'S07P', 'S08P', 'S09P', 
'S10P', 'S11P', 'S12P', 'S13P', 'S14P', 'S15P', 'S16P', 'S17P', 'S18P', 'S19P', 'S20P', 'S21P']
sensor_groups += [[s for s in tokamak.sensors_PA() if n in s.name] for n in PA_suffixes]

filaments = eigenmodes.ss_filaments(60)

(sensor_time, sensor_signal) = tokamak.get_sensor_time_signal(shot, sensor_groups[0][0].name)
chi_squares = np.zeros(len(sensor_time))
for sensors in sensor_groups:
    if len(sensors) == 0: continue
 
    signals = tokamak.sensor_signal_dict(shot, sensors, None, None, False)
    (sensor_time, sensor_signal) = tokamak.get_sensor_time_signal(shot, sensors[0].name)
 
    for i, sensor in enumerate(sensors):
         print sensor.name
         sensor_signal = signals[sensor.name]
     
         # Radial TA sensors point inwards
         if 'TA' in sensor.name: 
             if 'R' in sensor.name:
                 sensor.n_r = -1.0
     
         # Calculate field

         field_vf = fields.B_VF(vf_signal, tokamak.VFR, tokamak.VFZ, sensor)
         sensor.vf = field_vf
     
         # Calculate magnitude eddy current needed to create B_diff
         B_diff = field_vf - sensor_signal
         B_eddy1 = 0
         for f in filaments:
             greens = fields.greens_function(f.r, f.z-sensor.z, sensor.r,
                                             sensor.n_r, sensor.n_z)
             B_eddy1 += f.current_1*greens
         I_mag = [B/B_eddy1 for B in B_diff]
         sensor.eddy = I_mag

    group_signals = [signals[s.name] for s in sensors] 
    for i, time in enumerate(sensor_time):
        measurement_sum = 0
        for signal in group_signals: 
            measurement_sum += signal[time]
        measure_avg = measurement_sum/float(len(group_signals))
        measure_uncertainty = np.std([s[time] for s in group_signals])
        vf = sensors[0].vf
        eddy_sum = 0
        for sensor in sensors: 
            eddy_sum += sensor.eddy[i]
        eddy_avg = eddy_sum/float(len(sensors))

        if measure_uncertainty != 0:
            chi_squares[i] += (vf[i]+eddy_avg-measure_avg)**2/measure_uncertainty**2

plt.figure()
plt.xlabel('s')
plt.ylabel('Chi^2')
plt.plot(sensor_time, chi_squares)
plt.savefig(os.getcwd() + '/output/chi_squares.png')

