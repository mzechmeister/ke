#! /usr/bin/python

from ctypes import byref, c_int, c_double, POINTER, CFUNCTYPE
import ctypes
import os
import numpy as np

c_double_p = POINTER(c_double)
ptr_double = np.ctypeslib.ndpointer(dtype=np.float)


libdir = os.path.dirname(os.path.abspath(__file__)) + os.sep + 'lib' + os.sep
_ke_dbl = ctypes.CDLL(libdir+'ke_dbl_i.so')

args = [ptr_double,   # M
        c_double,     # e
        ptr_double,   # En
        ptr_double,   # cosE
        ptr_double,   # sinE
        c_int         # N
       ]

corefunc = CFUNCTYPE(c_double)

# generic vectorising function
_ke_dbl.v_Efunc.argtypes = [corefunc] + args # kernel function to be vectorised
#_ke.v_Efunc(ctypes.cast(_ke.Ecs,corefunc), np.asarray(M), E, 1., 7, 5)

dummy_arg = byref(c_double(0.))
_ke_dbl.i_Ecs.argtypes = [c_double, c_double, c_double_p, c_double_p, c_double_p, c_int]
_ke_dbl.i_Ecs.restype = c_double

def _M(E, e):
   '''Kepler's equations.'''
   return E - e*np.sin(E)

def _i_E(M, e, N=56):
   '''
   Wrapper to scalar c function i_Ecs.
   
   Example
   -------
   >>> _i_E(1.09, 1., 56)
   1.9995038046090956

   '''
   return _ke_dbl.i_Ecs(M, e, dummy_arg, dummy_arg, dummy_arg, N)


def Ecs(M, e, N=56, typ='i_Ecs'):
   '''
   Generic wrapper to vectorised c function, which returns E, cos E, and sin E.
   
   Parameters
   ----------
   M : scalar or array_like
       Mean anomaly [rad].
   e : scalar
       Eccentricity.
   typ : str, optional
       Available methods:
       i64_Ecs - fix-point with 64 bits integer
   
   Example
   -------
   >>> M = _M([0, 1], 1.)
   >>> M
   array([0.        , 0.15852902])
   >>> Ecs(M, 1.).T
   array([[1.88127160e-06, 1.00000000e+00, 1.88127160e-06],
          [9.99999999e-01, 5.40302308e-01, 8.41470984e-01]])
   >>> Ecs(M, 1., typ='i_Ecs').T
   array([[1.88127160e-06, 1.00000000e+00, 1.88127160e-06],
          [9.99999999e-01, 5.40302308e-01, 8.41470984e-01]])

   '''
   func = getattr(_ke_dbl, typ)
   n = np.size(M)
   En, ecosE, esinE = ECS = np.empty((3,n), dtype='float64')
   _ke_dbl.v_Efunc(ctypes.cast(func,corefunc), np.asarray(M, dtype='float64'), e, En, ecosE, esinE, n, N)
   return ECS


def i_E(*args, **kwargs):
   '''
   Wrapper to c function which returns E, cos E, and sin E.
   
   Example
   -------
   >>> M = [3.1, 1, 3, 4, 5]
   >>> i_E(M, 1.)
   array([3.12079558, 1.93456321, 3.07076673, 4.33817839, 5.33817839])
   
   '''
   return Ecs(*args, **kwargs)[0]


if __name__ == "__main__":
    '''
    Example
    -------
    ./_ke_dbl.py

    '''
    import doctest
    #from math import sin, sinh
    doctest.testmod()


