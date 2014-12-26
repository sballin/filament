#!/usr/bin/env python

import matplotlib.pyplot as plt
from sympy import mpmath
import os
import numpy
import math
import pickle

import fields 
import data_manipulation 
import signals
import shotput

conductSS = 1/72.0e-8
mu_0 = 4.0e-7*math.pi
asize = 0.005

# Get coil location data
with open('./coil_R_Z', 'r') as f:
    ((OHR, OHZ), (VFR, VFZ), (SHR, SHZ)) = pickle.load(f)


def top_shell_filaments(count):
    '''Return coordinates of evenly spaced filaments along top shell'''
    major_radius = 0.92 
    nominal_minor_radius = 0.16
    shell_thickness = 0.004763
    minor_radius = nominal_minor_radius + shell_thickness
    shell_length = math.pi/2.*minor_radius
    top_part = 0.0635
    gap = (top_part + shell_length)/float(count)
    angle_increment = gap/shell_length*math.pi/2.

    curve_filaments = [(major_radius+minor_radius*math.cos(theta),
                       minor_radius*math.sin(theta)) for theta in
                       numpy.arange(0, math.pi/2., angle_increment)]
    top_part_filaments = [(i, minor_radius) for i in 
                          numpy.arange(major_radius-0.5*gap, major_radius-top_part, -gap)]
    all_filaments = (curve_filaments + top_part_filaments)#[:count]
    all_filaments.reverse()
    return all_filaments 


def bottom_shell_filaments_mirror(count):
    return [(x,-z) for (x,z) in top_shell_filaments(count)]


def bottom_shell_filaments(count):
    '''Return coordinates of evenly spaced filaments along bottom shell'''
    major_radius = 0.92 
    nominal_minor_radius = 0.16
    shell_thickness = 0.004763
    minor_radius = nominal_minor_radius + shell_thickness
    shell_length = math.pi/2.*minor_radius
    gap = shell_length/float(count)

    curve_filaments = [(major_radius+minor_radius*math.cos(theta),
                       minor_radius*math.sin(theta)) for theta in
                       numpy.arange(-math.pi/2., 0, math.pi/2.*count)]
    return curve_filaments[:count]


def plot_geometry(filaments):
    (x,y) = zip(*filaments)
    plt.plot(x, y, '.')
    #plt.plot(VFR, VFZ, 'x')
    plt.axis('equal')
    
