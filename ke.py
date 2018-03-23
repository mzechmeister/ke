import numpy as np

from ctypes import c_int, c_long, c_double, CFUNCTYPE
import ctypes

ptr = np.ctypeslib.ndpointer

_ke = np.ctypeslib.load_library('lib/ke', '.')
"%x"%ctypes.addressof(_ke._E)
_ke_newton = np.ctypeslib.load_library('lib/ke_newton', '.')

# https://hakantiftikci.wordpress.com/2009/11/15/an-exampe-for-callbacks-in-ctypes/
# M = [3.1, 1, 3, 4,5]; E = ke._E(M, 1.)
# _ke.v_Efunc(ctypes.addressof(_ke._E), np.asarray(M), En, e, n, nM)
# CB_FUNC_TYPE = ctypes.CFUNCTYPE(None, ctypes.c_int)(ke._E)

#_cbspline.polyfit.restype = ctypes.c_int
# double v_E_ap(double *M, double *En, double e, int N, long nM) {
_ke.v_E.argtypes = [ ptr(dtype=np.float), # M
                     ptr(dtype=np.float), # En
                     c_double,            # e
                     c_int,               # n
                     c_long]              # N

corefunc = CFUNCTYPE(c_double, c_double, c_long)
corefunc = CFUNCTYPE(c_double)

# vectorising function
_ke.v_Efunc.argtypes = [corefunc,         # kernel function to be vectorised
                     ptr(dtype=np.float), # M
                     ptr(dtype=np.float), # En
                     c_double,            # e
                     c_int,               # n
                     c_long]              # N

#_ke.v_Efunc(ctypes.cast(_ke._E,corefunc), np.asarray(M), E, 1., 7, 5)

_ke_newton.v_cos.argtypes = _ke_newton.v_sin.argtypes = [ ptr(dtype=np.float), # M
                     ptr(dtype=np.float), # En
                     c_long]              # N

_ke_newton.vectorise.argtypes = [corefunc, # kernel function to be vectorised
                     ptr(dtype=np.float), # M
                     ptr(dtype=np.float), # En
                     c_double,            # e
                     c_double,            # n
                     c_long]              # N

_ke.v_Ecs.argtypes = [ ptr(dtype=np.float), # M
                     ptr(dtype=np.float), # En
                     ptr(dtype=np.float), # cosE
                     ptr(dtype=np.float), # sinE
                     c_double,            # e
                     c_int,               # n
                     c_long]              # N

def M(E, e):
   return E - e*np.sin(E)

def _E(M, e, n=29):
   '''Wrapper to c function.
   
   Example:
   --------
   >>> M = [3.1,1,3,4,5]
   >>> ke._E(M, 1.)
   
   '''
   En = np.empty_like(M)
   _ke.v_E(np.asarray(M), En, e, n, np.size(M))
   return En

def E(M, e, n=29, typ='_E'):
   '''Generic wrapper to c function.
   
   Example:
   --------
   >>> M = [3.1, 1, 3, 4,5]
   >>> ke._Eldexp(M, 1.)
   array([3.12079558, 1.93456321, 3.07076673, 3.57764002, 4.15262143])

   '''
   En = np.empty_like(M)
   func = {'_E':    _ke._E,
           'atr':   _ke._Eatr,
           'atrone': _ke._Eatr,
           'pn':    _ke._Epn,
           'ldexp': _ke._Eldexp,
           'N': _ke_newton._E_N,
           'myN': _ke_newton._E_myN}[typ]
   #print "%x"%ctypes.addressof(_ke._E)
   #_ke.v_Efunc(ctypes.addressof(_ke._E), np.asarray(M), En, e, n, nM)
   _ke.v_Efunc(ctypes.cast(func,corefunc), np.asarray(M), En, e, n, np.size(M))
   return En

def _E1N(M, e, n=29):
   ''' Wrapper to c function.
   
   Example:
   --------
   >>> M = [3.1,1,3,4,5]
   >>> ke._E(M, 1.)
   
   '''
   nM = np.size(M)
   En = np.empty_like(M)
   _ke.v_E1N(np.asarray(M), En, e, n, nM)
   return En

def _E_newton(M, e, eps=1e-7, typ='_E'):
   ''' Wrapper to c function.
   
   Example:
   --------
   >>> M = [3.1,1,3,4,5]
   >>> ke._E_newton(M, 1.)
   
   '''
   func = {'_E':    _ke_newton._E,
           'my':  _ke_newton._E_my}[typ]
   En = np.empty_like(M)
   _ke_newton.vectorise(ctypes.cast(func,corefunc), np.asarray(M), En, e, eps, np.size(M))
   return En

def _Ecs(M, e, n=29):
   ''' Wrapper to c function.
   
   Example:
   --------
   >>> M = [3.1,1,3,4,5]
   >>> ke._E(M, 1.)
   
   '''
   nM = np.size(M)
   En = np.empty_like(M)
   cosE = np.empty_like(M)
   sinE = np.empty_like(M)
   _ke.v_Ecs(np.asarray(M), En, cosE, sinE, e, n, nM)
   return En, cosE, sinE



