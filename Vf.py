#!/usr/bin/python

import numpy as np
from numpy import linalg
import pylab
from MDSplus import *
import sys
sys.path.append('/home/byrne/APS_Work/2011_Shaping_Coil/flux_surface_plots')
from hbtep import mdsplus

def Vf(shotnum):
    tree = Tree('hbtep2',shotnum)
    node = tree.getNode('sensors.Vf_Current')
    Vf = node.data()
    time = node.dim_of().data()
    timeindex1 = np.where(time >=0.0008)[0]
    timeindex2 = np.where(time <=.008)[0]
    Vfmax = max(Vf)
    #print len(time)
    t = time[min(timeindex1):max(timeindex2)]
    #print len(time)
    t = t-t[0]
    #print len(t)
    #print len(Vf[min(timeindex1):max(timeindex2)])
    #pylab.ylabel('Plasma Current (kA)')
    #pylab.plot(t*1000,Ip[min(timeindex1):max(timeindex2)]/1000,'g')
    #pylab.ylim([0,16])

    return(time, Vf)
    #pylab.legend(('Plasma Current'))#pylab.ylim([-V0,V0])
    #pylab.plot(time*1000,time*0.0+13, '--')

def niko(shotnum,time):
    '''Reads off the values at a series of given times, not merely the
    default digitization times in the tree'''
    tree = Tree('hbtep2',shotnum)
    node = 'sensors.Vf_Current'
    Vf = mdsplus.read_data(tree,node,times = time)
    return(Vf)
