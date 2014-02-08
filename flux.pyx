import math
from scipy import special, constants
from pylab import *

"""Green's function for magnetic vector potential."""
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

"""Flux Phi=2*pi*A."""
def flux(I, x, z, xc, zc):
    psi = 2*math.pi*x*I*green(x, z, xc, zc)
    return psi

"""Get flux through a slice of a loop by specifying length."""
def frac_flux(length, I, x, z, xc, zc):
    psi = length*I*green(x, z, xc, zc)
    return psi

"""Specify step size of range."""
def drange(start, stop, step):
    r = []
    rMax = start
    while rMax < stop:
        r.append(rMax)
        rMax += step
    return r

"""Flux through sensor due to VF coil current."""
def flux_val_calc(vf_signal, VFR, VFZ, sens_r, pos_z, sens_len, sens_width, sens_z_norm, coils):
    flux_vals = []
    for i in range(len(vf_signal)):
        flux_vals.append(0)
    for i_coil in range(len(VFR)-1):
        i_sig = 0
        for current in vf_signal:
            flux_vals[i_sig] += coils*sens_z_norm*(frac_flux(sens_len, current, sens_r, pos_z, VFR[i_coil], VFZ[i_coil]) - frac_flux(sens_len, current, sens_r-sens_width, pos_z, VFR[i_coil], VFZ[i_coil]))
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