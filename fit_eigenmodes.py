#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')   # headless
import matplotlib.pyplot as plt
import os
import fields
import data_manipulation 
import eigenmodes
import tokamak
import plotly.plotly as plotly
import time


def save_plot(sensor, sensor_time, sensor_signal, vf_time, field_vf, shot):
    difference = [field_vf[i]-sensor_signal[i] for i in xrange(len(sensor_time))]
    plt.figure()
    plt.plot(vf_time, field_vf, label='calculated field', color='blue')
    plt.plot(sensor_time, sensor_signal, label='sensor data', color='red')
    plt.plot(sensor_time, difference, label='difference', color='green')
    plt.ylabel('Tesla')
    plt.xlabel('s')
    title = sensor + ", shot no. " + str(shot)
    plt.suptitle(title)
    plt.legend()
    if 'PA' in sensor:
        plt.legend(frameon = 1).get_frame().set_facecolor('purple')
    elif 'TA' in sensor:
        plt.legend(frameon = 1).get_frame().set_facecolor('turquoise')
    elif 'FB' in sensor:
        plt.legend(frameon = 1).get_frame().set_facecolor('forestgreen')
    plt.xlim([.0015, .03])
    #plt.ylim([-0.002,.015])
    plt.savefig(os.getcwd() + '/output/eddy_images/' + sensor + '.png')
    plt.clf()


os.system('rm -f output/eddy_images/*')

shot = 81077
vf_time, vf_signal = tokamak.vf_data(shot)
sensors = tokamak.sensors_unique()
filaments = eigenmodes.ss_filaments(60)

I_mags = []
B_diffs = []

signals = tokamak.sensor_signal_dict(shot, sensors, None, None, False)
(sensor_time, sensor_signal) = tokamak.get_sensor_time_signal(shot, sensors[0].name)
# Fix length of integrated data
(sensor_time, sensor_signal) = data_manipulation.clip(sensor_time, sensor_signal) 

for i, sensor in enumerate(sensors):
    print sensor.name
    sensor_signal = signals[sensor.name]

    # Radial TA sensors point inwards
    if 'TA' in sensor.name: 
        if 'R' in sensor.name:
            sensor.n_r = -1.0

    # Calculate field
    field_vf = fields.B_VF(vf_signal, tokamak.VFR, tokamak.VFZ, sensor)

    # Calculate magnitude eddy current needed to create B_diff
    #B_diff = [field_vf[i]-sensor_signal[i] for i in xrange(len(sensor_time))]
    B_diff = field_vf - sensor_signal
    B_diffs.append(B_diff)
    B_eddy = 0
    for f in filaments:
        B_eddy += f.current_1*fields.greens_function(f.r, f.z-sensor.z, sensor.r,
                                          sensor.n_r, sensor.n_z)
    I_mag = [B/B_eddy for B in B_diff]
    I_mags.append(I_mag)

    save_plot(sensor.name, sensor_time, sensor_signal, vf_time, field_vf, shot)

# Plot current magnitudes
mpl_fig = plt.figure()
for i in xrange(len(I_mags)):
    name = sensors[i].name
    #if 'PA1_S32P' != name and 'PA1_S31P' != name and 'PA1_S30P' != name:
    if 'PA1' in name:
        pltcolor = 'purple'
    elif 'PA2' in name:
        pltcolor = 'pink'
    elif 'TA' in name:
        pltcolor = 'turquoise'
    elif 'FB' in name:
        pltcolor = 'forestgreen'
    plt.plot(sensor_time[0::20], I_mags[i][0::20], label=sensors[i].name, color=pltcolor)
    #plt.plot(sensor_time, B_diffs[i], color=pltcolor)
plt.title('Eddy current magnitude')
plt.xlabel('s')
plt.ylabel('A')
plt.savefig(os.getcwd() + '/output/I_mags.png')
#print plotly.plot_mpl(mpl_fig, filename="plotly version of an mpl figure")

# Delete/recreate overview image
os.system('rm -f output/eddy_images/montage.jpg')
os.system('montage -geometry 300x230+0+0 output/eddy_images/* output/eddy_images/montage.jpg')

