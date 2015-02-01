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
sensor_groups = [[s for s in tokamak.sensors_TA() if 'S1P' in s.name],
                 [s for s in tokamak.sensors_TA() if 'S2P' in s.name],
                 [s for s in tokamak.sensors_TA() if 'S3P' in s.name],
                 [s for s in tokamak.sensors_TA() if 'S4P' in s.name],
                 [s for s in tokamak.sensors_FB() if 'S1P' in s.name],
                 [s for s in tokamak.sensors_FB() if 'S2P' in s.name],
                 [s for s in tokamak.sensors_FB() if 'S3P' in s.name],
                 [s for s in tokamak.sensors_FB() if 'S4P' in s.name],
                 [s for s in tokamak.sensors_FB() if 'S1R' in s.name],
                 [s for s in tokamak.sensors_FB() if 'S2R' in s.name],
                 [s for s in tokamak.sensors_FB() if 'S3R' in s.name],
                 [s for s in tokamak.sensors_FB() if 'S4R' in s.name],
                 [s for s in tokamak.sensors_FB() if 'S4R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S01R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S02R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S03R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S04R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S05R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S06R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S07R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S08R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S09R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S10R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S11R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S12R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S13R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S14R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S15R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S16R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S17R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S18R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S19R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S20R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S21R' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S01P' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S02P' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S03P' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S04P' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S05P' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S06P' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S07P' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S08P' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S09P' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S10P' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S11P' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S12P' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S13P' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S14P' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S15P' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S16P' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S17P' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S18P' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S19P' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S20P' in s.name],
                 [s for s in tokamak.sensors_PA() if 'S21P' in s.name]]

filaments = eigenmodes.ss_filaments(60)

for sensors in sensor_groups:

    if len(sensors) == 0:
        continue
 
    I_mags = []
    B_diffs = []
    
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
     
         # Calculate magnitude eddy current needed to create B_diff
         #B_diff = [field_vf[i]-sensor_signal[i] for i in xrange(len(sensor_time))]
         B_diff = field_vf - sensor_signal
         B_diffs.append(B_diff)
         B_eddy1 = 0
         #B_eddy2 = 0
         for f in filaments:
             greens = fields.greens_function(f.r, f.z-sensor.z, sensor.r,
                                             sensor.n_r, sensor.n_z)
             B_eddy1 += f.current_1*greens
             #B_eddy2 += f.current_1*greens
         I_mag = [B/B_eddy1 for B in B_diff]
         I_mags.append(I_mag)
     
         #save_plot(sensor.name, sensor_time, sensor_signal, vf_time, field_vf, shot)
    
    # Plot current magnitudes
    mpl_fig = plt.figure()

    for i in xrange(len(I_mags)):
        plt.plot(sensor_time, I_mags[i], label=sensors[i].name)

    plt.legend()
    plt.title('Eddy current magnitudes')
    plt.xlabel('s')
    plt.ylabel('A')
    plt.savefig(os.getcwd() + '/output/I_mags' + sensors[0].name[0:2] 
                + sensors[0].name[5:]+ '.png')
    plt.close()
#print plotly.plot_mpl(mpl_fig, filename="plotly version of an mpl figure")

# Delete/recreate overview image
os.system('rm -f output/eddy_images/montage.jpg')
os.system('montage -geometry 300x230+0+0 output/eddy_images/* output/eddy_images/montage.jpg')

