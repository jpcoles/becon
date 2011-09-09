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
    d = np.abs(d)**2
    if l2 is None:
        #l0, = pl.plot(d.real, 'r-')
        #l1, = pl.plot(d.imag, 'g-')

        l2, = pl.plot(d, 'b-')
	pl.ylim(0,np.amax(d)*3)
    else:
        #l0.set_ydata(d.real)
        #l1.set_ydata(d.imag)
        l2.set_ydata(d)

    #pl.clf()
    #pl.plot(d, 'k-', alpha=0.4)
    pl.draw()

pl.show()
