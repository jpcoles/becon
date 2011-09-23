import numpy as np
import pylab as pl

pl.ion()

fname = '/tmp/becon.%05i.rho.txt'

x,y,z = np.loadtxt(fname % 1, unpack=True)
line, = pl.semilogx(x,y, 'k-')
for i in range(1,468,1):
    print i
    x,y,z = np.loadtxt(fname % i, unpack=True)
    line.set_ydata(y)
    pl.draw()

