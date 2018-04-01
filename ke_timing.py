#! /usr/bin/ipython
import numpy as np
import ke
from gplot import *
from pause import *

ke._Ecs(0.0001, 1., n=15)
ke.E(0.0001, 1., n=15, typ='atr')

x=np.arange(0,10,0.001); y=1*x; ke._ke_newton.v_cos(x,y, x.size); gplot(x,y-np.cos(x)) 
x=np.arange(0,10,0.001); y=1*x; ke._ke_newton.v_sin(x,y, x.size); gplot(x,y-np.sin(x)) 

e = 1.0
E = np.arange(-10,6,0.001)
E = np.arange(0,3,0.001); M = ke.M(E, e)

if 0:
   gplot(M, E, ke._E_newton(M, e), ke.E(M, e, n=55, typ='atr'), ' us 1:3, "" us 1:4')
   gplot('"" us 1:($3-$2), "" us 1:($4-$2)')
   gplot(M, np.log10(np.abs(ke._E(M, e, n=55)-E)))

   pause()


M = np.arange(0, np.pi, np.pi/1000);

timeit ke._E(M, e, n=29)
#1000 loops, best of 3: 1.51 ms per loop

timeit ke.E(M, e, n=29)
#1000 loops, best of 3: 1.51 ms per loop

timeit ke._Ecs(M, e, n=29)
#1000 loops, best of 3: 1.54 ms per loop

timeit ke.E(M, e, n=29, typ='pn')
#1000 loops, best of 3: 1.7 ms per loop

timeit ke.E(M, e, n=29, typ='atr')
#100 loops, best of 3: 2.1 ms per loop

timeit ke._E1N(M, e, n=19)
#1000 loops, best of 3: 966 µs per loop

timeit ke._E_newton(M, e)
#100 loops, best of 3: 2.95 ms per loop

timeit ke._E_newton(M, e)
#100 loops, best of 3: 3.22 ms per loop

timeit ke._E_newton(M, e, typ="my")
#100 loops, best of 3: 2.43 ms per loop

timeit ke.E(M, e, n=29)
timeit ke.E(M, e, n=2, typ="N")
# ~ same speed =>
# 29/2


