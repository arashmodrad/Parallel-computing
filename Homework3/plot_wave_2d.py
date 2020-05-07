import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.interpolate
from mpl_toolkits.mplot3d.axes3d import Axes3D


filename = str(sys.argv[1])
with open(filename, 'r') as f:
    lines = f.readlines()
    x = [float(line.split()[0]) for line in lines]
    y = [float(line.split()[1]) for line in lines]
    h = [float(line.split()[2]) for line in lines]
    u = [float(line.split()[3]) for line in lines]
    v = [float(line.split()[4]) for line in lines]

nx = int(math.sqrt(len(x)))

# Set up a regular grid of interpolation points
xi, yi = np.linspace(min(x), max(x), nx), np.linspace(min(y), max(y), nx)
xi, yi = np.meshgrid(xi, yi)

# Interpolate; there's also method='cubic' for 2-D data such as here
hi = scipy.interpolate.griddata((x, y), h, (xi, yi), method='linear')


fig = plt.figure()

# `ax` is a 3D-aware axis instance, because of the projection='3d' keyword argument to add_subplot
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.plot_surface(xi, yi, hi, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, alpha=0.75 )
ax.set_zlim(0.75, 1.25);

ax.set_xlim(min(x), max(x));
ax.set_ylim(min(y), max(y));

ax.set_xlabel('x [m]');
ax.set_ylabel('y [m]');
#plt.contour(zi, vmin=rho.min(), vmax=rho.max(), origin='lower',
#           extent=[x.min(), x.max(), y.min(), y.max()])

plt.show()
