import sys
import numpy as np
import pylab as pl

pl.ion()

l0 = None
l1 = None
l2 = None
for i,line in enumerate(sys.stdin):
    if i % 100 != 0: continue
    d = eval("np.array(%s)" % line)
    if l2 is None:
        d0 = np.empty_like(d)
        d0[:len(d)/2] = d[len(d)/2:]
        d0[len(d)/2:] = d[:len(d)/2]
        l2, = pl.plot(np.abs(d0)**2, 'b,')
        pl.ylim(0,1)
    else:
        d0[:len(d)/2] = d[len(d)/2:]
        d0[len(d)/2:] = d[:len(d)/2]
        l2.set_ydata(np.abs(d0)**2)

    pl.draw()

pl.show()
