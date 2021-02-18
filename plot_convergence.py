#!/usr/bin/env python

# Generates a convergence plot for Spiner.
# Requires scientific python stack

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import rc
import os
rc('font',size=14)
# mpl.rcParams['xtick.minor.size'] = 0
# mpl.rcParams['xtick.minor.width'] = 0

DFILE = "convergence.dat"
PLOTPATH = os.path.join("figs","convergence")

KX,KY,KZ = 2,3,4
xmin,xmax = 0,1
def test_function(y,x):
    return np.sin(2*np.pi*KX*x)*np.sin(2*np.pi*KY*y)

x = np.linspace(0,1,100)
X,Y = np.meshgrid(x,x)

data = np.loadtxt(DFILE)

fig,axarr = plt.subplots(1,2,figsize=(17,8))

axarr[0].plot(np.log2(data[:,0]),data[:,1],'ko-',lw=2)
axarr[0].set_yscale('log')
axarr[0].set_xticks(np.log2(data[:,0]))
axarr[0].set_xticklabels([str(d) for d in data[:,0]])
axarr[0].set_xlabel('Points per axis')
axarr[0].set_ylabel('Error')
axarr[0].grid(True)

mesh = axarr[1].pcolormesh(X,Y,test_function(Y,X))
mesh.set_edgecolor('face')
cbar = plt.colorbar(mesh)
cbar.set_label(r'$\sin(2\pi k_x x)\sin(2\pi k_y y)$')
axarr[1].set_xlabel(r'$x$')
axarr[1].set_ylabel(r'$y$')

for postfix in ['.pdf','.png']:
    plt.savefig(PLOTPATH+postfix,
                bbox_inches='tight')
                
