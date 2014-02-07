#!/usr/bin/env python

import math
from scipy import special, constants
from pylab import *
import MDSplus

tree = MDSplus.Tree('hbtep2',69523)
node = tree.getNode('sensors.magnetic.FB01_S1P')
sensor_signal = tree.getNode('sensors.Magnetic:FB01_S1P').data()
sensor_time = tree.getNode('sensors.Magnetic:FB01_S1P').dim_of().data()

def green(x, z, xc, zc):
    zDiff = math.fabs(z-zc)
    denom = xc*xc + x*x + zDiff*zDiff + 2*x*xc
    k2 = 4.0 * x * xc / denom
    k = math.sqrt(k2)
    ellipK = special.ellipk(k)
    ellipE = special.ellipe(k)
    try:
        fGreen = 4.0*xc*((2-k2)*ellipK-2*ellipE)/(constants.c*math.sqrt(denom)*k2)
    except ZeroDivisionError:
        fGreen = 0    
    return fGreen

def flux(I, x, z, xc, zc):
    psi = 2*math.pi*x*I*green(x, z, xc, zc)
    return psi

"""Specify step size of range."""
def drange(start, stop, step):
    r = []
    rMax = start
    while rMax < stop:
        r.append(rMax)
        rMax += step
    return r

def plotVf():
	plt.figure(1)
	plt.plot(sensor_time, sensor_signal)
	plt.show()

print 'STUFF BEGIN'
sensor_signal.append(0)
print len(sensor_signal)
print len(sensor_time)
print 'STUFF END'
plotVf()

