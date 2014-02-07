#!/usr/bin/env python

import math
from scipy import special, constants
from pylab import *
import MDSplus
#from hbtFlux import green, flux, drange
from dataManipulation import *

tree = MDSplus.Tree('hbtep2',69523)
node = tree.getNode('sensors.magnetic.FB01_S1P')
sensor_signal = tree.getNode('sensors.Magnetic:FB01_S1P').data()
sensor_time = tree.getNode('sensors.Magnetic:FB01_S1P').dim_of().data()

(sensor_time, sensor_signal) = clip(sensor_time, sensor_signal)

def plotVf():
	plt.figure(1)
	plt.plot(sensor_time, sensor_signal)
	plt.show()

print '-------------------------'
print len(sensor_signal)
print len(sensor_time)
print '-------------------------'
plotVf()

