#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')   # headless
import matplotlib.pyplot as plt
import os
import fields
import data_manipulation 
import eigenmodes
import tokamak


def save_plot(sensor, sensor_time, sensor_signal, vf_time, field_vals, shot):
    difference = [field_vals[i]-sensor_signal[i] for i in range(len(sensor_time))]
    plt.figure()
    plt.plot(vf_time, field_vals, label='calculated field', color='blue')
    plt.plot(sensor_time, sensor_signal, label='sensor data', color='red')
    plt.plot(sensor_time, difference, label='difference', color='green')
    plt.ylabel('Tesla')
    plt.xlabel('s')
    title = sensor + ", shot no. " + str(shot)
    plt.suptitle(title)
    plt.legend()
    plt.xlim([.0015, .03])
    #plt.ylim([-0.002,.015])
    plt.savefig(os.getcwd() + '/output/eddy_images/' + sensor + '.png')
    plt.clf()


shot = 81077
vf_time, vf_signal = tokamak.vf_shot_data(shot)
sensors = tokamak.sensors_PA(2)
filaments = eigenmodes.ss_filaments(60)

I_mags = []

for i, sensor in enumerate(sensors):
    (sensor_time, sensor_signal) = tokamak.get_sensor_time_signal(shot, sensor.name)
    # Fix length of integrated data
    (sensor_time, sensor_signal) = data_manipulation.clip(sensor_time, sensor_signal) 
    # Radial TA sensors point inwards
    if 'TA' in sensor.name: 
        if 'R' in sensor.name:
            sensor.n_r = -1.0

    # Calculate field
    field_vals = fields.B_VF(vf_signal, tokamak.VFR, tokamak.VFZ, 
                             sensor.r, sensor.z, sensor.n_r, sensor.n_z)

    # Calculate magnitude eddy current needed to create B_diff
    B_diff = [field_vals[i]-sensor_signal[i] 
              for i in range(len(sensor_time))]
    B_eddy = 0
    for f in filaments:
        B_eddy += (f.current_1 + f.current_2) \
                  *fields.greens_function(f.r, f.z-sensor.z, sensor.r,
                                          sensor.n_r, sensor.n_z)
    I_mag = [B/B_eddy for B in B_diff]
    I_mags.append(I_mag)

    # Print sensor name and save plot
    print sensor.name
    save_plot(sensor.name, sensor_time, sensor_signal, vf_time, field_vals, shot)

# Plot current magnitudes
plt.figure()
for i in range(len(I_mags)):
    name = sensors[i].name
    #if 'PA1_S32P' != name and 'PA1_S31P' != name and 'PA1_S30P' != name:
    plt.plot(sensor_time, I_mags[i], label=sensors[i].name)

plt.savefig(os.getcwd() + '/output/I_mags.png')


