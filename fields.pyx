from math import *
from pylab import *
from scipy import special, constants, integrate
from sympy import mpmath


def greens_function_integrate(R, z, y, n_y, n_z):
    f = lambda phi: 1e-7*R*(n_y*z*sin(phi)+n_z*(R-y*sin(phi)))/pow(R**2+y**2+z**2-2*y*R*sin(phi), 3/2.0)
    return integrate.quad(f, 0, 2*pi)[0]


def R_integral(R, y, z):
    return (4*R*sqrt(R**2+y**2+z**2)*(float(mpmath.ellipe(pi/4, -((4*R*y)/(R**2-2*R*y+y**2+z**2))))/sqrt((R**2+y**2+z**2)/(R**2-2*R*y+y**2+z**2))+float(mpmath.ellipe(pi/4, (4*R*y)/(R**2+2*R*y+y**2+z**2)))/sqrt((R**2+y**2+z**2)/(R**2+2*R*y+y**2+z**2))))/((R**2-2*R*y+y**2+z**2)*(R**2+2*R*y+y**2+z**2))


def sin_integral(R, y, z):
    return (2*(sqrt((R**2+y**2+z**2)/(R**2-2*R*y+y**2+z**2))*(R**4-2*R**3*y+2*R**2*(y**2+z**2)-2*R*y*(y**2+z**2)+(y**2+z**2)**2)*float(mpmath.ellipe(pi/4, -((4*R*y)/(R**2-2*R*y+y**2+z**2))))+(R**2+2*R*y+y**2+z**2)*((R**2+y**2+z**2)*sqrt((R**2+y**2+z**2)/(R**2+2*R*y+y**2+z**2))*float(mpmath.ellipe(pi/4, (4*R*y)/(R**2+2*R*y+y**2+z**2)))-(R**2-2*R*y+y**2+z**2)*(sqrt((R**2+y**2+z**2)/(R**2-2*R*y+y**2+z**2))*float(mpmath.ellipf(pi/4, -((4*R*y)/(R**2-2*R*y+y**2+z**2))))+sqrt((R**2+y**2+z**2)/(R**2+2*R*y+y**2+z**2))*float(mpmath.ellipf(pi/4, (4*R*y)/(R**2+2*R*y+y**2+z**2)))))))/(R*y*sqrt(R**2+y**2+z**2)*(R**2-2*R*y+y**2+z**2)*(R**2+2*R*y+y**2+z**2))


def greens_function(R, z, y, n_y, n_z):
    return 1e-7*R*(n_z*R_integral(R, y, z)+(n_y*z-n_z*y)*sin_integral(R, y, z))


def B_signal(vf_signal, VFR, VFZ, sens_r, pos_z, n_r, n_z):
    field_vals = [0 for i in vf_signal]
    for coil in range(len(VFR)):
        if VFR[coil] < 1.0:
            current_dir = -1.0   # current reversed in inner coils
        else:
            current_dir = 1.0
        time = 0
        factor = greens_function(VFR[coil], VFZ[coil]-pos_z, sens_r, n_r, n_z)
        for current in vf_signal:
            field = current*current_dir*factor
            field_vals[time] += field
            time += 1
    return field_vals


def OH_field(oh_signal, OHR, OHZ, sensor):
    field_vals = [0 for i in oh_signal]
    for coil in range(len(OHR)):
        current_dir = 1.0
        time = 0
        factor = greens_function(OHR[coil], OHZ[coil]-sensor.z, sensor.r, sensor.n_r, sensor.n_z)
        for current in oh_signal:
            field = current*current_dir*factor
            field_vals[time] += field
            time += 1
    return field_vals

