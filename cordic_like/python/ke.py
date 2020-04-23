import numpy as np

from ctypes import c_int, c_long, c_double, CFUNCTYPE
import ctypes

ptr_double = np.ctypeslib.ndpointer(dtype=np.float)

import os
libdir = os.path.dirname(os.path.abspath(__file__)) + os.sep + 'lib' + os.sep

_ke = np.ctypeslib.load_library(libdir+'ke', '.')
"%x"%ctypes.addressof(_ke._E)
_ke_newton = np.ctypeslib.load_library(libdir+'ke_newton', '.')

# https://hakantiftikci.wordpress.com/2009/11/15/an-exampe-for-callbacks-in-ctypes/
# M = [3.1, 1, 3, 4,5]; E = ke._E(M, 1.)
# _ke.v_Efunc(ctypes.addressof(_ke._E), np.asarray(M), En, e, n, nM)
# CB_FUNC_TYPE = ctypes.CFUNCTYPE(None, ctypes.c_int)(ke._E)

_ke.v_E.argtypes = [ ptr_double,   # M
                     ptr_double,   # En
                     c_double,     # e
                     c_int,        # n
                     c_long]       # N

corefunc = CFUNCTYPE(c_double, c_double, c_long)
corefunc = CFUNCTYPE(c_double)

# vectorising function
_ke.v_Efunc.argtypes = [corefunc,  # kernel function to be vectorised
                     ptr_double,   # M
                     ptr_double,   # En
                     c_double,     # e
                     c_int,        # n
                     c_long]       # N

#_ke.v_Efunc(ctypes.cast(_ke._E,corefunc), np.asarray(M), E, 1., 7, 5)

_ke_newton.v_cos.argtypes = \
_ke_newton.v_sin.argtypes = [ptr_double, # M
                     ptr_double,   # En
                     c_long]       # N

_ke_newton.vectorise.argtypes = [corefunc, # kernel function to be vectorised
                     ptr_double,   # M
                     ptr_double,   # En
                     c_double,     # e
                     c_double,     # eps
                     c_long]       # N

_ke.v_Ecs.argtypes = [ptr_double,  # M
                     ptr_double,   # En
                     ptr_double,   # cosE
                     ptr_double,   # sinE
                     c_double,     # e
                     c_int,        # n
                     c_long]       # N
'''
   >>> ke._Eldexp(M, 1.)
   array([3.12079558, 1.93456321, 3.07076673, 3.57764002, 4.15262143])
'''


def M(E, e):
   '''Kepler's equations.'''
   return E - e*np.sin(E)

def _E(M, e, n=29):
   '''
   Wrapper to c function of oneside CORDIC with half angle.
   
   Example
   -------
   >>> M = [1.0907, 2, 3.1, 4]
   >>> ke._E(M, 1.)
   array([1.99999818, 2.55419595, 3.12079558, 3.57764002])

   '''
   En = np.empty_like(M)
   _ke.v_E(np.asarray(M), En, e, n, np.size(M))
   return En

def E(M, e, n=29, typ='_E'):
   '''
   Generic wrapper to c function.
   
   Parameter
   ---------
   M : scalar or array_like
       Mean anomaly [rad].
   e : scalar
       Eccentricity [rad].
   n : int, optional
       Number of iterations.
   typ : str, optional
       Available methods:
       _E - onesided CORDIC with half angle base
       'atr' - twosided CORDIC with arctan radix
       'atrone': onesided CORDIC with atrctan radix
       'pn' - twosided CORDIC with half angle base
       'ldexp' - use ldexp for exponent manipulation (not faster)
       'E1N' - as _E with one Newton step
       'N' - Newton
   
   Example
   -------
   >>> M = ke.M(range(5), 1.)
   >>> M
   array([0.        , 0.15852902, 1.09070257, 2.85887999, 4.7568025 ])
   >>> ke.E(M, 1.)
   array([0.        , 0.99999999, 2.        , 3.        , 4.        ])
   >>> ke.E(M, 1., typ='ldexp')
   array([-1.21478072e-06,  9.99999999e-01,  2.00000000e+00,  3.00000000e+00,
           4.00000000e+00])

   '''
   En = np.empty_like(M)
   func = {'_E':    _ke._E,
           'atr':   _ke._Eatr,
           'atrone': _ke._Eatr,
           'pn':    _ke._Epn,
           'ldexp': _ke._Eldexp,
           'E1N': _ke._E1N,
           'N': _ke_newton._E_N,
           'myN': _ke_newton._E_myN
          }[typ]
   #print "%x"%ctypes.addressof(_ke._E)
   #_ke.v_Efunc(ctypes.addressof(_ke._E), np.asarray(M), En, e, n, nM)
   _ke.v_Efunc(ctypes.cast(func,corefunc), np.asarray(M), En, e, n, np.size(M))
   return En

def _E1N(M, e, n=29):
   '''
   Wrapper to c function v_E1N.
   
   Vectorised computation of E using onesided CORDIC with half angle base
   and one additional Newton iteration.
   
   Example
   -------
   >>> M = [3.1,1,3,4,5]
   >>> ke._E1N(M, 1.)
   
   '''
   En = np.empty_like(M)
   _ke.v_E1N(np.asarray(M), En, e, n, np.size(M))
   return En

def _E_newton(M, e, eps=1e-7, typ='_E'):
   '''
   Wrapper to c function.
   
   Newton methods.
   
   Example
   -------
   >>> M = [3.1,1,3,4,5]
   >>> ke._E_newton(M, 1.)
   
   '''
   func = {'_E': _ke_newton._E,
           'my': _ke_newton._E_my}[typ]
   En = np.empty_like(M)
   _ke_newton.vectorise(ctypes.cast(func,corefunc), np.asarray(M), En, e, eps, np.size(M))
   return En

def _Ecs(M, e, n=29):
   '''
   Wrapper to c function which returns E, cos E, and sin E.
   
   Example
   -------
   >>> M = [3.1,1,3,4,5]
   >>> ke._Ecs(M, 1.)
   
   '''
   En = np.empty_like(M)
   cosE = np.empty_like(M)
   sinE = np.empty_like(M)
   _ke.v_Ecs(np.asarray(M), En, cosE, sinE, e, n, np.size(M))
   return En, cosE, sinE



