#! /usr/bin/python

__author__ = 'Mathias Zechmeister'
__version__ = '2018-06-28'

__all__ = ['Ecs', 'Hcs']


# initialise cos-sin look-up table
from math import pi, cos, sin, floor
divtwopi = 1 / (2*pi)
acs = [(a, cos(a), sin(a)) for a in \
        [pi/2**i for i in range(1,60)]]

def Ecs(M, e, n=29):
   '''
   Compute the eccentric anomaly and its cos and sin-value CORDIC-like.
   
   Parameters
   ----------
   M : float
      Mean anomaly.
   e : float
      Eccentricity.
   n : integer
      Number of iterations.
   
   Example
   -------
   >>> Ecs(2-sin(2), 1)
   (1.99999999538762, -0.4161468323531165, 0.9092974287451092)
   '''
   E = 2 * pi * floor(M*divtwopi+0.5)
   cosE, sinE = 1., 0.
   for a,cosa,sina in acs[:n]:
      if E-e*sinE > M:
         a, sina = -a, -sina
      E += a
      cosE, sinE = cosE*cosa - sinE*sina,\
                   cosE*sina + sinE*cosa
   return E, cosE, sinE


# initialise cosh-sinh look-up table
from math import log, cosh, sinh, frexp, ldexp
ln2 = log(2)
acsh = [(a, cosh(a), sinh(a)) for a in \
        [4*ln2/2**i for i in range(1,60)]]

def Hcs(M, e, n=29):
   '''
   Compute the hyperbolic anomaly and its cosh and sinh-value CORDIC-like.
   
   Parameters
   ----------
   M : float
      Mean anomaly.
   e : float
      Eccentricity.
   n : integer
      Number of iterations.
   
   Example
   -------
   >>> Hcs(sinh(2)-2, 1)
   (1.9999999991222275, 3.7621956879000753, 3.626860404544669)
   >>> Hcs(-(sinh(5)-5), 1)
   (-5.000000000387744, 74.20994855355976, -74.20321060656327)
   '''
   m = max(0, frexp(M/e)[1])
   if M < 0: m = -m
   H = m * ln2
   coshH = ldexp(1, m-1) + ldexp(1, -m-1)
   sinhH = ldexp(1, m-1) - ldexp(1, -m-1)
   for a,cosha,sinha in acsh[:n]:
      if e*sinhH-H > M:
         a, sinha = -a, -sinha
      H += a
      coshH, sinhH = coshH*cosha + sinhH*sinha,\
                     coshH*sinha + sinhH*cosha
   return H, coshH, sinhH


if __name__ == "__main__":
   '''
   Running from command line.

   Without arguments the examples are checked with doctest:
   ~> ./ke_cordic.py

   With two arguments (mean anomaly and eccentricity) the eccentric:
   ~> ./ke_cordic.py 1.07 1.
   (1.9853116307811327, -0.4027463298095706, 0.915311637544788)
   (1.7649039106648186, 3.006107154374343, 2.834903917874572)

   With one arguments a graphic representation is plotted.
   '''
   import sys
   if len(sys.argv)>2:
      args = map(float, sys.argv[1:])
      if args[1] <= 1:
         print(Ecs(*args))
      if args[1] >= 1:
         print(Hcs(*args))
   elif len(sys.argv)==1:
      import doctest
      doctest.testmod()
      exec(doctest.script_from_examples(Ecs.__doc__))
      exec(doctest.script_from_examples(Hcs.__doc__))
   else:
      from gplot import *
      E = [i/1000. for i in range(-10000, 10000)]
      M = [t - sin(t)  for t in E]
      En = [Ecs(m, 1.)[0] for m in M]
      #gplot(M, E, En, 'w l, "" us 1:($3)')
      H = [i/1000. for i in range(-10000, 10000)]
      Mh = [sinh(h) - h for h in H]
      H0 = [Hcs(m, 1., n=0)[0] for m in Mh]
      Hn = [Hcs(m, 1.)[0] for m in Mh]
      gplot(M, E, En, 'w l, "" us 1:3,', Mh, H, Hn, H0, ' w l, "" us 1:3, "" us 1:4')
      #example()

