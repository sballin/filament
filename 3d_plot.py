#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, PathPatch
# register Axes3D class with matplotlib by importing Axes3D
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.text import TextPath
from matplotlib.transforms import Affine2D
import matplotlib.cm as cm
import colorsys
import tokamak
import eigenmodes
import plotly.plotly as py
import numpy as np


def text3d(ax, xyz, s, zdir="z", size=None, angle=0, usetex=False, **kwargs):
    x, y, z = xyz
    if zdir == "y":
        xy1, z1 = (x, z), y
    elif zdir == "y":
        xy1, z1 = (y, z), x
    else:
        xy1, z1 = (x, y), z
    text_path = TextPath((0, 0), s, size=size, usetex=usetex)
    trans = Affine2D().rotate(angle).translate(xy1[0], xy1[1])
    p1 = PathPatch(trans.transform_path(text_path), **kwargs)
    ax.add_patch(p1)
    art3d.pathpatch_2d_to_3d(p1, z=z1, zdir=zdir)


filaments = eigenmodes.ss_filaments(10)
sensors = tokamak.sensors_including_blacklist()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d


class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


for s in sensors:
    scale = 7.0
    if s.name[-1] == 'R': c = '#1E90FF'
    else: c = 'orange'
    if s.n_r == -1.0:
        a = Arrow3D([s.x,s.n_r/2.**.5/scale+s.x],[s.y,s.n_r/2.**.5/scale+s.y],[s.z,s.n_z/scale+s.z], mutation_scale=10, lw=1, arrowstyle="-|>", color=c)
    else:
        a = Arrow3D([s.x,s.n_x/scale+s.x],[s.y,s.n_y/scale+s.y],[s.z,s.n_z/scale+s.z], mutation_scale=10, lw=1, arrowstyle="-|>", color=c)
    ax.add_artist(a)

ax.scatter([s.x for s in sensors], [s.y for s in sensors], [s.z for s in sensors], c='red', marker='o', s=10)

norm = matplotlib.colors.Normalize(vmin=-0.5, vmax=1)
cmap = cm.hot
m = cm.ScalarMappable(norm=norm, cmap=cmap)


def pseudocolor(val, minval, maxval):
    # convert val in range minval..maxval to the range 0..120 degrees which
    # correspond to the colors red..green in the HSV colorspace
    h = (-1*(float(val-minval) / (maxval-minval)) * 120) % 360
    # convert hsv color (h,1,1) to its rgb equivalent
    # note: the hsv_to_rgb() function expects h to be in the range 0..1 not 0..360
    r, g, b = colorsys.hsv_to_rgb(h/360, 1., 1.)
    return r, g, b 


for f in filaments:
    c = pseudocolor(-f.current_1, -.30, .30)
    print str(c) + " " + str(f.current_1)
    p = Circle((0, 0), radius=f.r, color=c, linewidth=1.5, fill=False, alpha=0.5)
    ax.add_patch(p)
    art3d.pathpatch_2d_to_3d(p, z=f.z, zdir="z")

#for i in range(len(tokamak.VFR)):
#    p = Circle((0, 0), radius=tokamak.VFR[i], color='#8ffe09', linewidth=1.5, fill=False, alpha=0.5)
#    ax.add_patch(p)
#    art3d.pathpatch_2d_to_3d(p, z=tokamak.VFZ[i], zdir="z")


#for s in sensors:
#    name = s.name[0:2] + s.name[4:8]
#    text3d(ax, (s.x, s.y, s.z),
#           name,
#           zdir="y", size=0.07, usetex=False,
#           ec="none", fc="k")

lim = 0.8 
fig.tight_layout()
ax.set_xlim3d(-lim, lim)
ax.set_ylim3d(-lim, lim)
ax.set_zlim3d(-0.6, 0.6)
#ax.axis('off')

ax.view_init(30, 30)

plt.savefig('output/3d.png', dpi=150)

plt.show()
