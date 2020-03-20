#! /usr/bin/python

#
# This a demo for the computation of elementary functions with CORDIC.
#

from __future__ import division, print_function
from math import atan, atanh

__version__ = '2020-03-15'
__author__ = 'Mathias Zechmeister'

# shift sequence for linear, circular, and hyperbolic mode
# repeated in hyperbolic: i = 4, 13, 40, 121, ..., k, 3k+1, ...
Sl = list(range(60))
Sc = Sl
Sdbl = sorted(Sc + Sc[1:])     # double rotation
Sh = sorted(list(range(1,60)) + [4,13,40])
S = Sl, Sc, Sdbl, Sh

# rotation angles
al = [2**-i for i in Sl]
ac = [atan(2**-i) for i in Sc]
acd = [atan(2**-i) for i in Sdbl]
ah = [atanh(2**-i) for i in Sh]
a = al, ac, acd, ah

def cumprod(x):
   p = [x[0]]
   for xi in x[1:]:
      p.append(p[-1]*xi)
   return p

# scale corrections
Kl = [1] * len(Sl)
Kc = [1] + cumprod([(1+4**-i)**-0.5 for i in Sc[:-1]])
Kcd = [1] + cumprod([(1+4**-i)**-0.5 for i in Sdbl[:-1]])
Kh = [1] + cumprod([(1-4**-i)**-0.5 for i in Sh[:-1]])

K = Kl, Kc, Kcd, Kh


def rot(x, y, z, m=1, n=29, func=lambda x,y,z: z>0, verbose=False):
   '''
   Pure CORDIC core (floating point).

   Parameters
   ----------

   x, y, z : coordinate triple
   m : int
       Mode (0 linear, 1 circular, 3 (-1) hyperbolic).
   n : int
       Number of rotations.
   func : boolean callback function
       Function to decide for the rotation direction. If true is returned, then z is increased.

   Notes
   -----
    m Mode      Initial values     Functions
                xin   yin   zin    xR      yR      zR
    1 rotation  1     0     x      cos x   sin x   0
   -1 rotation  1     0     x      cosh x  sinh x  0
   -1 rotation  a     a     x      ae^x    ae^x    0
    1 vectoring 1     a     pi/2   ?a^2+1  0       acot(a)
   -1 vectoring a     1     0      ?a^2-1  0       acoth(a)
   -1 vectoring a+1   a-1   0      2?a     0       0.5 ln(a)
   -1 vectoring a+1/4 a-1/4 0      ?a      0       ln(1/4a)
   -1 vectoring a+b   a-b   0      2?ab    0       0.5 ln(a/b)

   Example
   -------
   # cos(pi/4)
   >>> rot(1., 0., math.pi/4), math.cos(math.pi/4)
   ((0.7071067796084353, 0.7071067827646602, -2.231788178302477e-09), 0.7071067811865476)

   '''
   mm = 1 if m==2 else m
   for i,(j,ai) in enumerate(list(zip(S[m], a[m]))[:n]):
      sgn = 1 if func(x,y,z) else -1
      if verbose:
         print("%2d %2d %9d x=%11.8f y=%11.8f  K=%10.8f ai=%11.8f zi=%11.8f" %
                   (i, j, 1<<j, x, y, K[m][i], sgn*ai, z), K[m][i]*x, K[m][i]*y)
      z -= sgn * ai
      x, y = x - mm*sgn*y/(1<<j),\
             y +    sgn*x/(1<<j)
   return K[m][n]*x, K[m][n]*y, z


def rect(r, x, **kwargs):
   '''
   Rectangular coordinates = rotation mode.

   Example
   -------
   >>> rect(2**0.5, math.pi/4)
   (0.9999999977682124, 1.0000000022317888)

   '''
   x, y, _ = rot(r, 0., x, func=lambda x,y,z: z>0, **kwargs)
   return x, y

def recth(r, x, **kwargs):
   '''
   Example
   -------
   '''
   x, y, _ = rot(r, 0., x, m=-1, func=lambda x,y,z: z>0, **kwargs)
   return x, y

def polar(x, y, **kwargs):
   '''
   Example
   -------
   # cos(pi/4)
   >>> polar(*cossin(math.pi/4)), math.pi/4
   ((1.0000000000000004, 0.7853981686162409), 0.7853981633974483)

   '''
   r, _, z = rot(x, y, 0., func=lambda x,y,z: y<0, **kwargs)
   return r, z

def polarh(x, y, **kwargs):
   '''
   Example
   -------
   >>> polarh(*coshsinh(1))
   (1.000000000000001, 0.9999999925533603)

   '''
   r, _, z = rot(x, y, 0., m=-1, func=lambda x,y,z: y<0, **kwargs)
   return r, z


def cossin(x, **kwargs):
   '''
   Example
   -------
   # cos(pi/4)
   >>> cossin(math.pi/4), math.cos(math.pi/4)
   ((0.7071067796084353, 0.7071067827646602), 0.7071067811865476)

   # cos(pi/8)
   >>> cossin(-math.pi/8), math.cos(-math.pi/8)
   ((0.9238795333214502, -0.38268343040918296), 0.9238795325112867)

   '''
   cosx, sinx = rect(1., x, **kwargs)
   return cosx, sinx

def coshsinh(x, **kwargs):
   '''
   Example
   -------
   >>> coshsinh(0)
   (1.0000000000000004, 4.121908780714506e-09)

   >>> coshsinh(1), (math.cosh(1), math.sinh(1))
   ((1.543080626063944, 1.1752011821530355), (1.5430806348152437, 1.1752011936438014))

   '''
   coshx, sinhx = recth(1., x, **kwargs)
   return coshx, sinhx


def cos(x, **kwargs):
   return cossin(x, **kwargs)[0]

def cosh(x, **kwargs):
   return coshsinh(x, **kwargs)[0]

def sin(x, **kwargs):
   return cossin(x, **kwargs)[1]

def sinh(x, **kwargs):
   return coshsinh(x, **kwargs)[1]

def aexp(x, a=1, **kwargs):
   ''' a*exp(x) '''
   expx, _, _ = rot(a, a, x, m=-1, func=lambda x,y,z: z>0, **kwargs)
   return expx

def exp(x, **kwargs):
   '''
   Example
   -------
   >>> exp(math.log(2))
   2.0000000088359626

   '''
   # expx = cosh(x) + sinh(x)
   # sum(coshsinh(x, **kwargs))
   return aexp(x, **kwargs)

def atan2(x, y, **kwargs):
   _, z = polar(x, y, **kwargs)
   return z

def atanh2(x, y, **kwargs):
   _, z = polarh(x, y, **kwargs)
   return z

def ln(x, **kwargs):
   z = 2 * atanh2(x+1, x-1, **kwargs)
   return z

def atan(y, **kwargs):
   '''
   Example
   -------
   >>> atan(1), math.pi/4
   (0.7853981611656603, 0.7853981633974483)
   '''
   _, z = polar(1., y, **kwargs)
   return z

def atanh(y, **kwargs):
   '''
   Example
   -------
   >>> atanh((2-1/2)/(2+1/2)), math.log(2)
   (0.693147184977926, 0.6931471805599453)
   '''
   _, z = polarh(1., y, **kwargs)
   return z

def asin(arg, n=29, **kwargs):
   '''
   Example
   -------
   >>> asin(0.5**0.5), math.pi/4
   (0.7853981611656603, 0.7853981633974483)

   >>> asin(0.5), math.pi/6
   (0.5235987769420017, 0.5235987755982988)

   '''
   _, _, z = rot(Kc[n], 0., 0., n=n, func=lambda x,y,z: y<arg, **kwargs)
   return -z

def acos(arg, n=29, **kwargs):
   _, _, z = rot(Kc[n], 0., 0., n=n, func=lambda x,y,z: x>arg, **kwargs)
   return -z

def asinh(arg, n=29):
   '''
   Example
   -------
   >>> asin(0.5**0.5), math.pi/4
   (0.7853981611656603, 0.7853981633974483)

   >>> asin(0.5), math.pi/6
   (0.5235987769420017, 0.5235987755982988)

   '''
   _, _, z = rot(Kh[n], 0., 0., m=-1, func=lambda x,y,z: y<arg, n=n)
   return -z


def sqrt(x, **kwargs):
   '''
   Compute sqrt with CORDIC.

   Example
   -------
   >>> sqrt(0.25)
   0.5000000000000002

   '''
   r, z = polarh(x+1/4, x-1/4, **kwargs)
   return r

def sqrt_short(r, n=42):
   '''
   Compute sqrt with CORDIC.

   Does not need z (and therefore no look-up table).

   With start vector x0 = r + 1/4 and y0 = r - 1/4 and hyperbolic mode,
   the sqrt of r will be in x when y -> 0.
   sqrt(r) = sqrt(x**2 - y**2) * Kn^-1
           = sqrt(x_0**2 - y_0**2)
           = sqrt(2*r/4 + 2*r/4) = sqrt(r)
   0.03 < r < 2.34

   Example
   -------
   >>> sqrt_short(0.25)
   0.5000000000000002

   '''
   if not 0.03 <= r <= 2.33:
      print('DomainError')
   x = r + 1/4.
   y = r - 1/4.
   for i in Sh[:n]:
      sgn = 1 if y<0 else -1
      x, y = x + sgn*y*2**-i,\
             y + sgn*x*2**-i
   return x * Kh[n]


# Double iteration

def dasin(arg, n=32, **kwargs):
   '''
   Double iteration.

   Example
   -------
   >>> asin(0.5**0.5), math.pi/4
   (0.7853981611656603, 0.7853981633974483)

   >>> asin(0.5), math.pi/6
   (0.5235987769420017, 0.5235987755982988)
   '''
   _, _, z = rot(Kcd[n], 0., 0., m=2, func=lambda x,y,z: y<arg, n=n, **kwargs)
   return -z

def dacos(arg, n=29, **kwargs):
   _, _, z = rot(Kcd[n], 0., 0., m=2, func=lambda x,y,z: x>arg, n=n, **kwargs)
   return -z


if __name__ == "__main__":
    '''
    Example
    -------
    ./cordic.py

    '''
    import doctest
    import math
    doctest.testmod()
