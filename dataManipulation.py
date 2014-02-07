#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Copyright (C) Nikolaus Rath <Nikolaus@rath.org>

This program can be distributed under the terms of the GNU LGPL.

Utility functions for HBT-EP related programming.
'''
from __future__ import print_function, absolute_import, division
import sys
sys.path.append('/home/byrne/repos/hbt/python/hbtep')
from mdsplus import read_data, signals_from_shot
from fnmatch import fnmatch
from scipy.integrate import cumtrapz
from scipy.interpolate.interpolate import interp1d
from scipy.linalg import LinAlgError, cholesky
from scipy.stats.stats import nanmean, nanstd
import MDSplus
import argparse
import idlsave
import logging
import math
import numpy as np
import os
import re
import pylab as plt
import scipy.optimize as optimization

#def main():
tstart = -.0015
t0 =  -0.001
t1 =  0.06

def plot_subtraction(shot):
    ((C1_VF_std,C1_OH_std,C1_SH_std),(C1_VF_fft,C1_OH_fft,C1_SH_fft)) = responses()

    tree = MDSplus.Tree('hbtep2',shot)
    #tree = MDSplus.Tree('hbtep2',82646)
    #tree = MDSplus.Tree('hbtep2',81077)
    #tree = MDSplus.Tree('hbtep2',73056)
 
    (IVF, time_VF) = pull_data(tree, 'sensors.VF_Current')
    (IOH, time_OH) = pull_data(tree, 'sensors.OH_Current')
    (ISH, time_SH) = pull_data(tree, 'sensors.SH_Current')

    IVF = zero_out(time_VF,IVF)   
    IOH = zero_out(time_OH,IOH)
    ISH = zero_out(time_SH,ISH)

    (time_VF,IVF) = clip(time_VF,IVF)
    (time_OH,IOH) = clip(time_OH,IOH)
    (time_SH,ISH) = clip(time_SH,ISH)

    (C1, time_C1) = pull_data(tree,'sensors.rogowskis.cos_1')
    (IP, time_IP) = pull_data(tree,'sensors.rogowskis.ip')
    
    C1 = zero_out(time_C1,C1)
    IP = zero_out(time_IP,IP)

    (time_C1, C1) = clip(time_C1,C1)
    (time_IP, IP) = clip(time_C1,IP)

    plt.figure()
    plt.title('Standard Subtraction vs FFt for shot #' + str(shot))
    plt. xlabel('time (s)')
    plt.ylabel('error (pickup - subtraction)')

    plt.plot(time_C1,C1)

    C1_fft = fft_subtract(C1,C1_OH_fft,IOH)
    plt.plot(time_C1,C1_fft,'c')

    C1_fft = fft_subtract(C1_fft,C1_VF_fft,IVF)
    plt.plot(time_C1,C1_fft,'k')

    C1_fft = fft_subtract(C1_fft,C1_SH_fft,ISH)
    plt.plot(time_C1,C1_fft,'m')

    C1 = C1 - IOH*C1_OH_std
    plt.plot(time_C1,C1,'g--')

    C1 = C1 - IVF*C1_VF_std
    plt.plot(time_C1,C1,'r--')

    C1 = C1 - ISH*C1_SH_std
    plt.plot(time_C1,C1,'y--')

    plt.figure(9)
    plt.title('VF traces')
    plt.xlabel('time (s)')
    plt.ylabel('current (A)')
    plt.plot(time_VF,IVF,label = str(shot))
    plt.legend()

def fft_subtract(signal,response,drive):
    pickup_signal = response*np.fft.fft(drive)
    subtracted = np.fft.ifft(np.fft.fft(signal) - pickup_signal)
    return(subtracted)

def pull_data(tree,nodename):
    signal = tree.getNode(nodename).data()
    time = tree.getNode(nodename).dim_of().data()
    return(signal, time)

def zero_out(time,signal):
    signal -= np.mean(signal[np.where(time < -.0005)[0]])
    return(signal)

def clip(time, signal):
    index = np.where((time > t0) * (time < t1))[0]
    signal = signal[index]
    time = time[index]
    return (time,signal)

def responses():
    tree_VF = MDSplus.Tree('hbtep2',73056)
    tree_OH = MDSplus.Tree('hbtep2',73062)
    ##need to tailor values of time window to individul shot
    ##These numbers ensure ISH is within 95% of max during this window
    tree_SH = MDSplus.Tree('hbtep2',81759)
    #tSH0 = .0032
    #tSH1 = .0065

    (IVF, time_VF) = pull_data(tree_VF, 'sensors.VF_Current')
    (IOH, time_OH) = pull_data(tree_OH, 'sensors.OH_Current')
    (ISH, time_SH) = pull_data(tree_SH, 'sensors.SH_Current')

    IVF = zero_out(time_VF,IVF)   
    IOH = zero_out(time_OH,IOH)
    ISH = zero_out(time_SH,ISH)

    (time_VF,IVF) = clip(time_VF,IVF)
    (time_OH,IOH) = clip(time_OH,IOH)
    (time_SH,ISH) = clip(time_SH,ISH)

    (C1_VF, time_C1_VF) = pull_data(tree_VF,'sensors.rogowskis.cos_1')
    (IP_VF, time_IP_VF) = pull_data(tree_VF,'sensors.rogowskis.ip')
    
    C1_VF = zero_out(time_C1_VF,C1_VF)
    IP_VF = zero_out(time_IP_VF,IP_VF)

    (time_C1_VF, C1_VF) = clip(time_C1_VF,C1_VF)
    (time_IP_VF, IP_VF) = clip(time_C1_VF,IP_VF)

    (C1_OH, time_C1_OH) = pull_data(tree_OH,'sensors.rogowskis.cos_1')
    (IP_OH, time_IP_OH) = pull_data(tree_OH,'sensors.rogowskis.ip')
    
    C1_OH = zero_out(time_C1_OH,C1_OH)
    IP_OH = zero_out(time_IP_OH,IP_OH)

    (time_C1_OH, C1_OH) = clip(time_C1_OH,C1_OH)
    (time_IP_OH, IP_OH) = clip(time_C1_OH,IP_OH)

    (C1_SH, time_C1_SH) = pull_data(tree_SH,'sensors.rogowskis.cos_1')
    (IP_SH, time_IP_SH) = pull_data(tree_SH,'sensors.rogowskis.ip')
    
    C1_SH = zero_out(time_C1_SH,C1_SH)
    IP_SH = zero_out(time_IP_SH,IP_SH)

    (time_C1_SH, C1_SH) = clip(time_C1_SH,C1_SH)
    (time_IP_SH, IP_SH) = clip(time_C1_SH,IP_SH)

    (C1_VF_std) =  optimization.curve_fit(func, IVF, C1_VF)[0][0]

    (C1_OH_std) =  optimization.curve_fit(func, IOH, C1_OH)[0][0]

    (C1_SH_std) =  optimization.curve_fit(func, ISH, C1_SH)[0][0]
    
    (IP_VF_std) =  optimization.curve_fit(func, IVF, IP_VF)[0][0]

    (IP_OH_std) =  optimization.curve_fit(func, IOH, IP_OH)[0][0]

    (IP_SH_std) =  optimization.curve_fit(func, ISH, IP_SH)[0][0]

    C1_VF_fft = ( np.fft.fft(C1_VF) / np.fft.fft(IVF) )
                        
    C1_OH_fft = ( np.fft.fft(C1_OH) / np.fft.fft(IOH) )
    
    C1_SH_fft = ( np.fft.fft(C1_SH) / np.fft.fft(ISH) )
 
    IP_VF_fft = ( np.fft.fft(IP_VF) / np.fft.fft(IVF) )
                        
    IP_OH_fft = ( np.fft.fft(IP_OH) / np.fft.fft(IOH) )  

    IP_SH_fft = ( np.fft.fft(IP_SH) / np.fft.fft(ISH) )
    #return(time,IOH,IVF,C1_VF,C1_OH,C1_VF_std,C1_OH_std)
    return((C1_VF_std,C1_OH_std,C1_SH_std),(C1_VF_fft,C1_OH_fft,C1_SH_fft))

def func(x,a):
    return a*x

def plot_subtraction_qian():
    (C1_VF_std) = responses_std_qian()
    (C1_VF_fft) = responses_fft_qian()

    tree = MDSplus.Tree('hbtep2',82646)
    #tree = MDSplus.Tree('hbtep2',81077)
    #tree = MDSplus.Tree('hbtep2',73056)
 
    IVF = tree.getNode('sensors.VF_Current').data()
    time_VF = tree.getNode('sensors.VF_Current').dim_of().data()
 
    indexVF = np.where((time_VF>t0) * (time_VF <t1))[0]
    IVF -= np.mean(IVF[np.where(time_VF < -.0005)[0]])
 
    IVF = IVF[indexVF]
    time_VF = time_VF[indexVF]
 
    C1_raw = tree.getNode('sensors.rogowskis.cos_1.raw').data()
    C1_time_raw = tree.getNode('sensors.rogowskis.cos_1.raw').dim_of().data()
  
    C1 = tree.getNode('sensors.rogowskis.cos_1').data()
    C1_time = tree.getNode('sensors.rogowskis.cos_1').dim_of().data()
  
    C1_raw -= np.mean(C1_raw[np.where(C1_time_raw<-.0005)])
    C1 -= np.mean(C1[np.where(C1_time<-.0005)])

    C1_raw_index = np.where((C1_time_raw>t0) * (C1_time_raw<t1))[0]
    C1_raw = C1_raw[C1_raw_index]
    C1_time_raw = C1_time_raw[C1_raw_index]

    C1_index = np.where((C1_time>t0) * (C1_time <t1))[0]
    C1 = C1[C1_index]
    C1_time = C1_time[C1_index]

    plt.figure(1)
    plt.plot(C1_time_raw,C1_raw)

    C1_raw = np.fft.ifft( np.fft.fft(C1_raw)-C1_VF_fft*np.fft.fft(IVF))
    
    plt.plot(C1_time_raw,C1_raw)

    plt.figure(2)
    plt.plot(C1_time,C1)
    plt.plot(C1_time,C1-IVF*C1_VF_std)

    C1_raw_int = (cumtrapz(C1_raw, C1_time_raw) + C1_raw[:-1] * 
                  tree.getNode('.sensors.rogowskis:cos_1:rc_time').data()) 

    plt.plot(C1_time[1:],C1_raw_int)

def responses_std_qian():
    tree_VF = MDSplus.Tree('hbtep2',73056)
 
    IVF = tree_VF.getNode('sensors.VF_Current').data()

    time_IVF = tree_VF.getNode('sensors.VF_Current').dim_of().data()

    indexIVF = np.where((time_IVF>t0) * (time_IVF <t1))[0]
 
    IVF -= np.mean(IVF[np.where(time_IVF < -.0005)[0]])

    C1_VF = tree_VF.getNode('sensors.rogowskis.cos_1').data()

    time_C1_VF = tree_VF.getNode('sensors.rogowskis.cos_1').dim_of().data()

    indexC1VF = np.where((time_C1_VF>t0) * (time_C1_VF <t1))[0]
 
    C1_VF -= np.mean(C1_VF[np.where(time_C1_VF<-.0005)])
  
    (C1_VF_std) =  optimization.curve_fit(func, IVF[indexIVF], 
                                          C1_VF[indexC1VF], 0)[0][0]

    return(C1_VF_std)


def responses_fft_qian():
    tree_VF = MDSplus.Tree('hbtep2',73056)
 
    IVF = tree_VF.getNode('sensors.VF_Current').data()
 
    timeVF = tree_VF.getNode('sensors.VF_Current').dim_of().data()

    indexVF = np.where((timeVF>t0) * (timeVF<t1))[0]

    IVF -= np.mean(IVF[np.where(timeVF < -.0005)[0]])

    C1_VF = tree_VF.getNode('sensors.rogowskis.cos_1.raw').data()
    time_C1_VF = tree_VF.getNode('sensors.rogowskis.cos_1.raw').dim_of().data()
    indexc1vf = np.where((time_C1_VF>t0) * (time_C1_VF <t1))[0]

    C1_VF -= np.mean(C1_VF[np.where(timeVF<-.0005)])
 
    C1_VF_fft = ( np.fft.fft(C1_VF[indexc1vf]) / 
                  np.fft.fft(IVF[indexVF]) )
    
    return(C1_VF_fft)
