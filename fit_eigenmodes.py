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
import pickle
import scipy


shot = 81077
vf_time, vf_signal = tokamak.vf_data(shot)
sensor_groups = []

TA_suffixes = ['S1P', 'S2P', 'S3P', 'S4P']
sensor_groups += [[s for s in tokamak.sensors('TA') if n in s.name] for n in TA_suffixes]

FB_suffixes = ['S1P', 'S2P', 'S3P', 'S4P', 'S1R', 'S2R', 'S3R', 'S4R', 'S4R']
FB_bad = ['S1P', 'S4P']
sensor_groups += [[s for s in tokamak.sensors('FB') if n in s.name] for n in FB_suffixes
                  if n not in FB_bad]

PA_suffixes = ['S01R', 'S02R', 'S03R', 'S04R', 'S05R', 'S06R', 
'S07R', 'S08R', 'S09R', 'S10R', 'S11R', 'S12R', 'S13R', 'S14R', 'S15R', 'S16R', 'S17R', 'S18R',
'S19R', 'S20R', 'S21R', 'S01P', 'S02P', 'S03P', 'S04P', 'S05P', 'S06P', 'S07P', 'S08P', 'S09P', 
'S10P', 'S11P', 'S12P', 'S13P', 'S14P', 'S15P', 'S16P', 'S17P', 'S18P', 'S19P', 'S20P', 'S21P']
PA_bad = ['S10P', 'S11P', 'S13R', 'S14P', 'S08P']
sensor_groups += [[s for s in tokamak.sensors('PA') if n in s.name] for n in PA_suffixes
                  if n not in PA_bad]

STEP = 1
filaments = eigenmodes.ss_filaments(60)

(sensor_time, sensor_signal) = tokamak.get_sensor_time_signal(shot, sensor_groups[0][0].name)
I_eddy_mags = np.arange(-2000, 2000, 5)
chi_squares = np.zeros((len(sensor_time[0::STEP]),(len(I_eddy_mags))))

signals = tokamak.sensor_signal_dict(shot, None, None, None, False)

outfile = 'output/G_eigs_1.p'
loaded_G_eigs_1 = True
try: G_eigs_1 = pickle.load(open(outfile, 'rb'))
except IOError: 
    G_eigs_1 = {}
    loaded_G_eigs_1 = False

for sensors in sensor_groups:
    if len(sensors) == 0: continue

    (sensor_time, sensor_signal) = tokamak.get_sensor_time_signal(shot, sensors[0].name)

    print sensors[0].name[0:2] + sensors[0].name[4:8]

    vf = fields.B_VF(vf_signal[0::STEP], tokamak.VFR, tokamak.VFZ, sensors[0])

    if not loaded_G_eigs_1:
        G_eig_1 = 0
        for f in filaments:
            G_eig_1 += f.current_1*fields.greens_function(f.r, f.z-sensors[0].z, sensors[0].r,
                                                      sensors[0].n_r, sensors[0].n_z)
        G_eigs_1[sensors[0].name] = G_eig_1
    else: G_eig_1 = G_eigs_1[sensors[0].name]

    group_signals = [signals[s.name] for s in sensors]
    for i, time in enumerate(sensor_time[0::STEP]):
        measurement_sum = 0
        for signal in group_signals: 
            measurement_sum += signal[i]
        measure_avg = measurement_sum/float(len(group_signals))
        measure_stddev = np.std([s[i] for s in group_signals])
        for j, I_eddy_mag in enumerate(I_eddy_mags):
            if measure_stddev != 0:
                chi_squares[i][j] += (vf[i]+I_eddy_mag*G_eig_1-measure_avg)**2/measure_stddev**2

print 'Making plots'

if not loaded_G_eigs_1: pickle.dump(G_eigs_1, open(outfile, 'wb'))

I_eddy_time = np.array([I_eddy_mags[np.argmin(chi_squares[i])] for i in range(chi_squares.shape[0])])

plt.figure()
plt.plot(sensor_time[0::STEP], I_eddy_time)
plt.ylabel('Eddy current magnitude (A)')
plt.xlabel('Seconds')
plt.title('Eddy current magnitude over shot 81077')
plt.savefig(os.getcwd() + '/output/I_eddy_time.png')

for sensors in sensor_groups:
    if len(sensors) == 0: continue
    
    suffix = sensors[0].name[0:2] + sensors[0].name[4:8]
    print suffix

    B_diffs = np.zeros((len(sensors), len(sensor_time)))
    B_VF = fields.B_VF(vf_signal, tokamak.VFR, tokamak.VFZ, sensors[0])
    for i, sensor in enumerate(sensors):
        (sensor_time, sensor_signal) = tokamak.get_sensor_time_signal(shot, sensor.name)
        B_diffs[i] = [B_VF[j] - sensor_signal[j] for j in xrange(len(sensor_time))]

    B_diff_avg = np.mean(B_diffs, axis=0)
    B_diff_sem = scipy.stats.sem(B_diffs, axis=0)
    B_eddy = fields.B_filaments(I_eddy_time, filaments, sensors[0])

    plt.figure()
    plt.errorbar(sensor_time, B_diff_avg, yerr=B_diff_sem, fmt='-o', label='B_sensor-B_VF')
    plt.plot(sensor_time[0::STEP], B_eddy, color='r', label='B_eddy calculated')
    plt.xlim([.0015, .03])
    plt.xlabel('Seconds')
    plt.ylabel('Magnetic field strength (Tesla)')
    plt.legend()
    plt.title(suffix)
    plt.savefig('output/eddy_images/' + suffix + '.png')
    plt.clf()
 
os.system('rm -f output/eddy_images/montage.jpg')
os.system('montage -geometry +0+0 output/eddy_images/* output/eddy_images/montage.jpg')

