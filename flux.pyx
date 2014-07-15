from math import *
from scipy import special, constants, integrate
from pylab import *
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
    field_vals = []
    for i in range(len(vf_signal)):
        field_vals.append(0)
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


def OH_field(oh_signal, OHR, OHZ, sens_r, pos_z, n_r, n_z):
    field_vals = []
    for i in range(len(oh_signal)):
        field_vals.append(0)
    for coil in range(len(OHR)):
        current_dir = 1.0
        time = 0
        factor = greens_function(OHR[coil], OHZ[coil]-pos_z, sens_r, n_r, n_z)
        for current in oh_signal:
            field = current*current_dir*factor
            field_vals[time] += field
            time += 1
    return field_vals


def grad_green(x, z, xc, zc):
    cdef double k2 = 4.0 * x * xc / ((x + xc) * (x + xc) + (z - zc) * (z - zc))
    cdef double m1 = 1.0 - k2
    cdef double k = sqrt(k2)
    cdef double EllipK = special.ellipk(k)
    cdef double EllipE = special.ellipe(k)

    cdef double GG = -sqrt(x * xc / k2) * ((2.0 - k2) * EllipK - 2.0 * EllipE)
    cdef double dKE = EllipE / m1 - EllipK
    cdef double dGx = 0.25 * (GG * k2 * (x + xc) / (x * xc) - sqrt(k2 / (x * xc)) * (2.0 * xc - k2 * (x + xc)) * dKE)
    cdef double dGz = 0.25 * k2 * ((z - zc) / (x * xc)) * (GG + sqrt(k2 * x * xc) * dKE)
    return (dGx, dGz)


def green_times_current(vf_signal, sensor_signal, VFR, VFZ, sens_r, pos_z, n_r, n_z):
    field_vals = []
    cdef double field
    cdef double current_dir
    cdef int time
    for i in range(len(vf_signal)):
        field_vals.append(0)
    for coil in range(len(VFR)):
        if VFR[coil] < 1.0:
            current_dir = -1.0   # current reversed in inner coils
        else:
            current_dir = 1.0
        (dGx, dGz) = grad_green(sens_r, pos_z, VFR[coil], VFZ[coil])
        time = 0
        for current in vf_signal:
            field = 4*pi*1e-7*current_dir*current*n_z/(2*pi*VFR[coil])*(n_z*dGx-n_r*dGz)
            field_vals[time] += field
            time += 1
        # neg_field = 0
        # neg_signal = 0
        # for i in range(len(field_vals)):
        #     if field_vals[i] < 0:
        #         neg_field += 1
        #     if sensor_signal[i] < 0:
        #         neg_signal += 1
        # if neg_field > 100:
        #     for i in range(len(field_vals)):
        #         field_vals[i] *= -1
        # if neg_signal > 100:
        #     for i in range(len(sensor_signal)):
        #         sensor_signal[i] *= -1
    for i in range(len(field_vals)):
        field_vals[i] = abs(field_vals[i])
        sensor_signal[i] = abs(sensor_signal[i])
    return (field_vals, sensor_signal)


def green_base():
    (dGx, dGz) = grad_green(0.00000001, 0.00000001, 7, 0.00000001)
    field_val = 2*10e-7/(0.00000001)*dGz
    return field_val


"""Green's function for magnetic vector potential."""
def green(x, z, xc, zc):
    zDiff = math.fabs(z-zc)
    denom = xc*xc + x*x + zDiff*zDiff + 2*x*xc
    k2 = 4.0*x*xc/denom
    k = math.sqrt(k2)
    ellipK = special.ellipk(k)
    ellipE = special.ellipe(k)
    try:
        return 4.0*xc*((2-k2)*ellipK-2*ellipE)/(constants.c*math.sqrt(denom)*k2)
    except ZeroDivisionError:
        return 0

"""Flux Phi=2*pi*A."""
def flux(I, x, z, xc, zc):
    return 2*math.pi*x*I*green(x, z, xc, zc)

"""Get flux through a slice of a loop by specifying length."""
def frac_flux(length, I, x, z, xc, zc):
    return length*I*green(x, z, xc, zc)

def drange(start, stop, step):
    r = []
    rMax = start
    while rMax < stop:
        r.append(rMax)
        rMax += step
    return r

"""Flux through sensor due to VF coil current."""
def sensor_flux_calc(vf_signal, VFR, VFZ, sens_r, pos_z, sens_len, sens_width, sens_z_norm, coils):
    flux_vals = []
    for i in range(len(vf_signal)):
        flux_vals.append(0)
    for i_coil in range(len(VFR)-1):
        i_sig = 0
        for current in vf_signal:
            flux_vals[i_sig] += coils*sens_z_norm*(frac_flux(sens_len, current, sens_r, pos_z, VFR[i_coil], VFZ[i_coil]) - frac_flux(sens_len, current, sens_r-sens_width, pos_z, VFR[i_coil], VFZ[i_coil]))
            i_sig += 1
    return flux_vals

def full_flux_calc(vf_signal, VFR, VFZ, sens_r, pos_z):
    flux_vals = []
    for i in range(len(vf_signal)):
        flux_vals.append(0)
    for i_coil in range(len(VFR)-1):
        i_sig = 0
        for current in vf_signal:
            flux_vals[i_sig] += flux(current, sens_r, pos_z, VFR[i_coil], VFZ[i_coil])
            i_sig += 1
    return flux_vals

"""Show behavior of Green's function at edge cases."""
def samplePlots():
    # Flux values for radius of flux loop expanding through xc
    xValsThrough = drange(0, 20, 0.01)
    fluxXThrough = []
    for x in xValsThrough:
        # Current loop: xc = 10,   zc = 0
        # Flux loop:    x  = 0-20, z  = 0
        fluxXThrough.append(flux(1, x, 0, 10, 0))

    # Flux values for position of flux loop going through zc
    zValsThrough = drange(-200, 200, 0.01)
    fluxZThrough = []
    for z in zValsThrough:
        # Current loop: xc = 10, zc = 10
        # Flux loop:    x  = 10, z  = -200-200
        fluxZThrough.append(flux(1, 10, z, 10, 0))

    # Flux values for radius of flux loop expanding outward from xc
    xValsAway = drange(0, 50, 0.1)
    fluxXAway = []
    for x in xValsAway:
        # Current loop: xc = 10,   zc = 10
        # Flux loop:    x  = 0-50, z  = 0
        fluxXAway.append(flux(1, x, 0, 10, 10))

    # Flux values for position of flux loop going away from zc
    zValsAway = drange(0, 100, 0.01)
    fluxZAway = []
    for z in zValsAway:
        # Current loop: xc = 10,   zc = 0
        # Flux loop:    x  = 10,   z  = 0-100
        fluxZAway.append(flux(1, 10, z, 10, 0))

    plt.figure(1)

    plt.subplot(221)
    plt.plot(xValsThrough, fluxXThrough)
    plt.xlabel('flux loop radius (x)')
    plt.ylabel('flux(1, x, 0, 10, 0)')
    plt.title("radius goes through xc")

    plt.subplot(222)
    plt.plot(zValsThrough, fluxZThrough)
    plt.xlabel('flux loop position (z)')
    plt.ylabel('flux(1, 10, z, 10, 0)')
    plt.title('z goes through zc')

    plt.subplot(223)
    plt.plot(xValsAway, fluxXAway)
    plt.xlabel('flux loop radius (x)')
    plt.ylabel('flux(1, x, 0, 10, 10)')
    plt.title("radius goes away from xc")

    plt.subplot(224)
    plt.plot(zValsAway, fluxZAway)
    plt.xlabel('flux loop position (z)')
    plt.ylabel('flux(1, 10, z, 10, 0)')
    plt.title('z goes away from zc')

    plt.tight_layout()
    plt.show()
