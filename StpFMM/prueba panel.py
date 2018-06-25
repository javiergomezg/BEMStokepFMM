import numpy as np
import cuadratura
from matplotlib import pyplot, cm
from matplotlib import rcParams


nq = 5

X1 = 1.0
X2 = 3.0
Y1 = 1.0
Y2 = 2.0

X = (X1, X2)
Y = (Y1, Y2) 

t = nodos_gauss_chebyshev(nq)
Xt, Yt = nodos_gauss_panel(t, X1, Y1, X2, Y2)


print t
print 
print Xt
print 
print Yt


fig2 = pyplot.figure(num = None, figsize = (8, 6), dpi = 80, 
                 facecolor = 'w', edgecolor = 'k')
pyplot.plot(X,Y, 'bo-', linewidth = 3)
pyplot.plot(Xt,Yt, 'ro', linewidth = 3)
#    pyplot.plot(XX,YY, 'bo-', linewidth = 3)
pyplot.grid(True)
pyplot.axis([X1-0.5,X2+0.5,Y1-0.5,Y2+0.5])
#    pyplot.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
#    pyplot.savefig('circulito.png')
pyplot.show()

