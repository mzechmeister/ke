'''
A wrapper to vectorise the pure python code of ke_dbl for numpy.
Only some builtin functions need their numpy analogs.
Note these changes are global in the python session!

This wrapper is only to demonstrate that the code can be vectorised and parallelised.
It is faster that the pure python version, but slower than the c wrapper,
which is also vectorised and meant for normal use.

A numpy version could be obtain by modifying the math import to
from numpy import copysign, ldexp, round

Example
-------
>>> from ke_dbl_np import i64_Ecs
>>> i64_Ecs([1,2,3], 1.)
(array([1.93456321, 2.55419595, 3.07076673]), array([-0.35579714, -0.83238624, -0.99749289]), array([0.93456321, 0.55419595, 0.07076673]))

'''

import numpy as np
import builtins
import math

# patch the functions before loading ke_cordic_dbl
#org = builtins.round, math.copysign, math.ldexp   # attempt to keep the changes local
#int_org = builtins.int 
builtins.round = lambda x: np.round(x).astype(int)
#builtins.int = lambda x: np.array(x).astype(int)   # frozen importlib._bootstrap_external
math.copysign = np.copysign
math.ldexp = np.ldexp
math.tau = np.array(2*np.pi)   
# the first division with tau will convert a list for M to an array, i.e.
# i64_Ecs([1,2,3], 1.) vs i64_Ecs(np.array([1,2,3]), 1.)

from ke_cordic_dbl import *

#builtins.round, math.copysign, math.ldexp = org    # does not work, also the 





