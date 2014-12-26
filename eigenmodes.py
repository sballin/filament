#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')   # headless
import matplotlib.pyplot as plt
import os
import numpy
import math
from sympy import mpmath
import numpy as np
import pickle

import fields 
import data_manipulation 
import signals 
import shotput
import geometry


def resistivity_matrix(count):
    rho_stainless_steel = 6.90e-7
    length_shell = 0.622655
    tauCoverageReduction = (360 - 10*28)/360.0
    conductSS = 1.38889e6
    thickSS = 0.0032
    shell_count = 8 # GUESS
    conductance = (length_shell/float(shell_count)) * geometry.mu_0 * tauCoverageReduction * conductSS * thickSS
    resistance = 1/float(conductance)
    return resistance*numpy.identity(count)


def self_inductance(x, z, a):
    return x*(math.log(8.*x/a)-1.75)


def green(x, z, xc, zc):
    k2 = 4*xc*x/((xc + x)**2 + (z - zc)**2)
    return -math.sqrt(x*xc/k2)*((2-k2)*mpmath.ellipk(k2)-2*mpmath.ellipe(k2))


def inductance((x1, z1), (x2, z2), a):
    if (x1, z1) == (x2, z2):
        return self_inductance(x1, z1, a) 
    else:
        return -green(x1, z1, x2, z2)


def inductance_matrix(filaments):
    count = len(filaments)
    if not os.path.isfile('output/inductances%d.p' % count):
        inductances = np.zeros((count, count))
        for i in range(count):
            for j in range(count):
                inductances[i][j] = inductance(filaments[i], filaments[j], geometry.asize)
        pickle.dump(inductances, open('output/inductances%d.p' % count, 'wb'))
        return inductances
    else:
        print 'DANGER: loaded inductances%d.p instead of calculating' % count
        return pickle.load(open('output/inductances%d.p' % count, 'rb'))


def get_currents(inductances, resistances):
    print "AVG(ORIG-PINV(PINV(ORIG))):" + str(np.average(resistances-np.linalg.pinv(np.linalg.pinv(resistances))))

    m = inductances.dot(np.linalg.pinv(resistances))

    m = np.append(m, np.ones((m.shape[0], 1)), axis=1)
    m = np.append(m, np.ones((1, m.shape[1])), axis=0)
    m[m.shape[0]-1][m.shape[1]-1] = 0

    eigenvals, eigenvecs = np.linalg.eig(m)
    return eigenvals, eigenvecs


filaments = geometry.top_shell_filaments(63)
geometry.plot_geometry(filaments)
inductances = inductance_matrix(filaments)
                      
print inductances.shape

plt.figure()

plt.subplot(221)
plt.title('Filaments')
geometry.plot_geometry(filaments)

plt.subplot(222)
plt.title('Inductances')
plt.ylim(0, inductances.shape[0])
plt.xlim(0, inductances.shape[1])
plt.pcolor(inductances)
plt.colorbar()

eigenvals, eigenvecs = get_currents(inductances, resistivity_matrix(inductances.shape[0]))

plt.subplot(223)
plt.title('L/R times')
plt.plot(range(len(eigenvals)), eigenvals, '.')
plt.yscale('log')

plt.subplot(224)
plt.title('Eigenmodes for top 3 L/R')
plt.xlabel('Toward outer midplane ->')
modes_1 = eigenvecs[:, 2] 
modes_2 = eigenvecs[:, 3] 
modes_3 = eigenvecs[:, 4] 
plt.plot(range(len(modes_1)), modes_1, '.')
plt.plot(range(len(modes_2)), modes_2, '.')
plt.plot(range(len(modes_3)), modes_3, '.')
plt.savefig(os.getcwd() + '/output/eigenmodes.png')


