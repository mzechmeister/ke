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
   No output if successful.

   With two arguments (mean anomaly and eccentricity) the eccentric and hyperbolic
   anomaly are returned:
   ~> ./ke_cordic.py 1.07 1.
   (1.9853116307811327, -0.4027463298095706, 0.915311637544788)
   (1.7649039106648186, 3.006107154374343, 2.834903917874572)

   With one argument a graphic representation is plotted.
   '''
   import sys
   args = map(float, sys.argv[1:])
   if not args:
      import doctest
      doctest.testmod()
   elif len(args) == 1:
      from gplot import Gplot
      gplot = Gplot("-p")
      e = 1.0
      E = [i/1000. for i in range(-10000, 10000)]
      M = [t - e*sin(t)  for t in E]
      En = [Ecs(m, e)[0] for m in M]
      H = E[5000:-5000]
      Mh = [e*sinh(t) - t for t in H]
      H0 = [Hcs(m, e, n=0)[0] for m in Mh]
      Hn = [Hcs(m, e)[0] for m in Mh]
      gplot(M, E, En, 'w l, "" us 1:3,', Mh,
            H, Hn, H0, ' w l t "H", "" us 1:3 t "Hn", "" us 1:4 t "H0', tmp="$")
   else:
      if args[1] <= 1:
         print(Ecs(*args))
      if args[1] >= 1:
         print(Hcs(*args))
