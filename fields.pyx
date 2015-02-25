import math
import numpy as np
from scipy import integrate
from sympy import mpmath


def greens_function_integrate(R, z, y, n_y, n_z):
    f = lambda phi: 1e-7*R*(n_y*z*math.sin(phi)+n_z*(R-y*math.sin(phi)))/pow(R**2+y**2+z**2-2*y*R*math.sin(phi), 3/2.0)
    return integrate.quad(f, 0, 2*math.pi)[0]


def R_integral(R, y, z):
    return (4*R*math.sqrt(R**2+y**2+z**2)*(float(mpmath.ellipe(math.pi/4, -((4*R*y)/(R**2-2*R*y+y**2+z**2))))/math.sqrt((R**2+y**2+z**2)/(R**2-2*R*y+y**2+z**2))+float(mpmath.ellipe(math.pi/4, (4*R*y)/(R**2+2*R*y+y**2+z**2)))/math.sqrt((R**2+y**2+z**2)/(R**2+2*R*y+y**2+z**2))))/((R**2-2*R*y+y**2+z**2)*(R**2+2*R*y+y**2+z**2))


def sin_integral(R, y, z):
    return (2*(math.sqrt((R**2+y**2+z**2)/(R**2-2*R*y+y**2+z**2))*(R**4-2*R**3*y+2*R**2*(y**2+z**2)-2*R*y*(y**2+z**2)+(y**2+z**2)**2)*float(mpmath.ellipe(math.pi/4, -((4*R*y)/(R**2-2*R*y+y**2+z**2))))+(R**2+2*R*y+y**2+z**2)*((R**2+y**2+z**2)*math.sqrt((R**2+y**2+z**2)/(R**2+2*R*y+y**2+z**2))*float(mpmath.ellipe(math.pi/4, (4*R*y)/(R**2+2*R*y+y**2+z**2)))-(R**2-2*R*y+y**2+z**2)*(math.sqrt((R**2+y**2+z**2)/(R**2-2*R*y+y**2+z**2))*float(mpmath.ellipf(math.pi/4, -((4*R*y)/(R**2-2*R*y+y**2+z**2))))+math.sqrt((R**2+y**2+z**2)/(R**2+2*R*y+y**2+z**2))*float(mpmath.ellipf(math.pi/4, (4*R*y)/(R**2+2*R*y+y**2+z**2)))))))/(R*y*math.sqrt(R**2+y**2+z**2)*(R**2-2*R*y+y**2+z**2)*(R**2+2*R*y+y**2+z**2))


def greens_function(R, z, y, n_y, n_z):
    return 1e-7*R*(n_z*R_integral(R, y, z)+(n_y*z-n_z*y)*sin_integral(R, y, z))


def B_VF(vf_signal, VFR, VFZ, sensor):
    G = 0
    for coil in xrange(len(VFR)):
        if VFR[coil] < 1.0: current_dir = -1.0   # current reversed in inner coils
        else: current_dir = 1.0
        G += current_dir*greens_function(VFR[coil], VFZ[coil]-sensor.z, sensor.r, sensor.n_r, sensor.n_z)
    return G*np.array(vf_signal)


def B_filaments(I_magnitude_timeseries, filaments, sensor):
    G = 0
    for f in filaments:
        G += f.current_1*greens_function(f.r, f.z-sensor.z, sensor.r, sensor.n_r, sensor.n_z)
    return G*np.array(I_magnitude_timeseries) 


def OH_field(oh_signal, OHR, OHZ, sensor):
    G = 0
    for coil in xrange(len(OHR)):
        G += greens_function(OHR[coil], OHZ[coil]-sensor.z, sensor.r, sensor.n_r, sensor.n_z)
    return G*np.array(oh_signal)

