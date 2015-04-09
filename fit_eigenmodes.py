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
import math
import pickle
import scipy
import pylab


shot = 81077
coil = 'VF'
coil_time, coil_signal = tokamak.coil_current(shot, coil)
sensor_groups = []

#TA_suffixes = ['S1P', 'S2P', 'S3P', 'S4P']
#sensor_groups += [[s for s in tokamak.sensors('TA') if n in s.name] for n in TA_suffixes]

FB_suffixes = ['S1P', 'S2P', 'S3P', 'S4P', 'S1R', 'S2R', 'S3R', 'S4R', 'S4R']
FB_bad = ['S1P', 'S4P']
sensor_groups += [[s for s in tokamak.sensors('FB') if n in s.name] for n in FB_suffixes
                  if n not in FB_bad]

#PA_suffixes = ['S01R', 'S02R', 'S03R', 'S04R', 'S05R', 'S06R', 
#'S07R', 'S08R', 'S09R', 'S10R', 'S11R', 'S12R', 'S13R', 'S14R', 'S15R', 'S16R', 'S17R', 'S18R',
#'S19R', 'S20R', 'S21R', 'S01P', 'S02P', 'S03P', 'S04P', 'S05P', 'S06P', 'S07P', 'S08P', 'S09P', 
#'S10P', 'S11P', 'S12P', 'S13P', 'S14P', 'S15P', 'S16P', 'S17P', 'S18P', 'S19P', 'S20P', 'S21P']
#PA_bad = ['S10P', 'S11P', 'S13R', 'S14P', 'S08P']
#sensor_groups += [[s for s in tokamak.sensors('PA') if n in s.name] for n in PA_suffixes
#                  if n not in PA_bad]

no_filaments = 20 
filaments = eigenmodes.ss_filaments(no_filaments)

(sensor_time, sensor_signal) = tokamak.get_sensor_time_signal(shot, sensor_groups[0][0].name)

all_sensors = []
for sensor_group in sensor_groups:
	all_sensors += sensor_group
signals = tokamak.sensor_signal_dict(shot, all_sensors, None, None, False)

outfile = 'output/G_eigs_1.p'
loaded_G_eigs_1 = True
try: G_eigs_1 = pickle.load(open(outfile, 'rb'))
except IOError: 
    G_eigs_1 = {}
    loaded_G_eigs_1 = False

dChi_sq_dI = np.zeros((len(sensor_time), 2))

for sensors in sensor_groups:
    if len(sensors) == 0: continue

    print sensors[0].name[0:2] + sensors[0].name[4:8]
    (sensor_time, sensor_signal) = tokamak.get_sensor_time_signal(shot, sensors[0].name)
    coil_field = fields.coil_field(coil_signal, coil, sensors[0])

    if loaded_G_eigs_1:
        G_eig_1 = G_eigs_1[sensors[0].name]
    else:
        G_eig_1 = 0
        for f in filaments:
            G_eig_1 += f.current_1*fields.greens_function(f.r, f.z-sensors[0].z, sensors[0].r,
                                                          sensors[0].n_r, sensors[0].n_z)
        G_eigs_1[sensors[0].name] = G_eig_1

    group_signals = [signals[s.name] for s in sensors]
    for i, time in enumerate(sensor_time):
        measurement_sum = 0
        for signal in group_signals: 
            measurement_sum += signal[i]
        measure_avg = measurement_sum/float(len(group_signals))
        measure_stddev = np.std([s[i] for s in group_signals])
        dChi_sq_dI[i][0] += -1*G_eig_1*(coil_field[i]-measure_avg)/measure_stddev**2
        dChi_sq_dI[i][1] += G_eig_1**2/measure_stddev**2

I_eddy_time = np.zeros(len(sensor_time))
for i in range(len(sensor_time)):
    I_eddy_time[i] = dChi_sq_dI[i][0]/dChi_sq_dI[i][1]

if not loaded_G_eigs_1: pickle.dump(G_eigs_1, open(outfile, 'wb'))

print 'Caculating chi-squared values'

chi_squares = np.zeros(len(I_eddy_time))
for sensors in sensor_groups:
    if len(sensors) == 0: continue

    print sensors[0].name[0:2] + sensors[0].name[4:8]
    (sensor_time, sensor_signal) = tokamak.get_sensor_time_signal(shot, sensors[0].name)
    coil_field = fields.coil_field(coil_signal, coil, sensors[0])

    if loaded_G_eigs_1:
        G_eig_1 = G_eigs_1[sensors[0].name]
    else:
        G_eig_1 = 0
        for f in filaments:
            G_eig_1 += f.current_1*fields.greens_function(f.r, f.z-sensors[0].z, sensors[0].r,
                                                          sensors[0].n_r, sensors[0].n_z)
        G_eigs_1[sensors[0].name] = G_eig_1

    group_signals = [signals[s.name] for s in sensors]
    for i, time in enumerate(sensor_time):
        measurement_sum = 0
        for signal in group_signals: 
            measurement_sum += signal[i]
        measure_avg = measurement_sum/float(len(group_signals))
        measure_stddev = np.std([s[i] for s in group_signals])
        chi_squares[i] += (coil_field[i] + I_eddy_time[i] - measure_avg)**2/measure_stddev**2

print 'Making plots'

plt.figure()
plt.plot(sensor_time, chi_squares)
plt.ylabel('Chi-squared value')
plt.xlabel('Seconds')
plt.xlim([.0015, .027])
plt.title('Chi-squared values over shot ' + shot)
plt.savefig('output/chi_squared_time.png')
plt.clf()

plt.figure()
plt.plot(sensor_time, I_eddy_time)
plt.ylabel('Eddy current magnitude (A)')
plt.xlabel('Seconds')
plt.xlim([.0015, .027])
plt.title('Eddy current magnitude over shot ' + shot)
plt.savefig('output/I_eddy_time.png')
plt.clf()

for sensors in sensor_groups:
    if len(sensors) == 0: continue
    
    suffix = sensors[0].name[0:2] + sensors[0].name[4:8]
    print suffix

    B_diffs = np.zeros((len(sensors), len(sensor_time)))
    B_coil = fields.coil_field(coil_signal, coil, sensors[0])
    for i, sensor in enumerate(sensors):
        (sensor_time, sensor_signal) = tokamak.get_sensor_time_signal(shot, sensor.name)
        B_diffs[i] = sensor_signal

    B_diff_avg = np.mean(B_diffs, axis=0)
    B_diff_sem = scipy.stats.sem(B_diffs, axis=0)

    B_eddy = fields.filament_field(I_eddy_time, filaments, sensors[0])
    B_composite = B_eddy + B_coil

    fig = plt.figure()
    plt.errorbar(sensor_time, B_diff_avg, yerr=B_diff_sem, fmt='-', linewidth=2, ecolor='#006BB2', alpha=0.25, label='<B_sensor>')
    plt.plot(sensor_time, B_composite, 'r-', label='B_' + coil + '+B_eddy')
    plt.plot(sensor_time, B_coil, '--', label='B_' + coil)
    plt.xlim([.0015, .027])
    plt.xlabel('Seconds')
    plt.ylabel('Magnetic field strength (Tesla)')
    plt.legend()
    plt.title(suffix)
    plt.savefig('output/eddy_images/' + suffix + '.png')
    fig.tight_layout()
    plt.clf()
 
os.system('rm -f output/eddy_images/montage.jpg')
os.system('montage -geometry +0+0 output/eddy_images/* output/eddy_images/montage.jpg')

