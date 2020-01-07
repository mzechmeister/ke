from __future__ import print_function
import ctypes

class ieee(float):
   """
   Analyse a float according to IEEE 754.

   Information about bits, sign, mantissa, and exponents are added as attributes.

   x : float

   Examples
   --------
   >>> x = ieee(18.4, bits=32)
   >>> x
   0 10000011 00100110011001100110011
   >>> x.d
   18.4
   >>> x.bin
   '01000001100100110011001100110011'
   >>> x.i
   1100165939L
   >>> x * 2
   0 10000000100 0010011001100110011001100110011001100110011001100110
   >>> x.d
   18.4
   >>> x.hex()
   '0x1.2666666666666p+4'
   >>> from math import frexp
   >>> frexp(x)
   (0.575, 5)
   >>> x.frexp
   (0.574999988079071, 5L)

   # (Not exact because 32 bit vs. 64 bit.)

   Notes
   -----
   x = SEEEEEEEEMMMMMMMMMMMMMMMMMMMMMMM
     = s * m * 2**e

   """
   def __new__(self, x, **kwargs):
       # required because float class is immutable
       return super(ieee, self).__new__(ieee, x)
   def __init__(self, x, bits=64):
       self.d = x
       self.bits = bits
       if bits == 64:
          self.r = 11   # number of exponent bits
          self.p = 52   # number of mantissa bits
          # get binary presentation
          self.i = ctypes.c_uint64.from_buffer(ctypes.c_double(x)).value
       else:
          self.r = 8
          self.p = 23
          self.i = ctypes.c_uint32.from_buffer(ctypes.c_float(x)).value
       self.bin = "{:0{p.bits}b}".format(self.i, p=self)
       self.B = (1 << self.r-1) - 1   # exponent bias
       self.E = (self.i >> self.p) & ((1 << self.r) - 1)
       self.e = self.E - self.B
       self.S = self.i >> 63                   # signbit [-1|0]
       self.s = (-1)**self.S                   # sign [-1|1]
       self.M = self.i & ((1 << self.p) - 1)   # mantisse bits
       self.m = 1 + self.M/float(1<<self.p)    # mantisse
       self.frexp = self.m/2, self.e+1
       self.sme = self.s, self.m, self.e       # sign, mantissa, exponent for s*m*2**e
       self.SEM = self.S, self.E, self.M
   def __repr__(self):
       fmtbin = "{:b} {:0{p.r}b} {:0{p.p}b}".format(*self.SEM, p=self)
       return fmtbin

def as_ieee(func):
    def function_wrapper(*x):
        return ieee(func(*x))
    return function_wrapper

# ieee.__add__ = as_ieee(float.__add__)
# '__divmod__', '__float__'  '__floordiv__',
# monkeypatch decorator
for m in ['__abs__', '__add__', '__div__', '__mod__', '__mul__',
          '__neg__', '__pos__', '__pow__', '__radd__', '__rdiv__',
          '__rmod__', '__rmul__', '__rpow__', '__rsub__', '__rtruediv__',
          '__sub__', '__truediv__']:
   setattr(ieee, m, as_ieee(getattr(float, m)))


class fx():
   """Fix point representation.

   Takes the mantissa bits, prepends the hidden bit, scales with the exponent
   by e left shifts, and puts the fix point by m left shifts.
   fx(x) is a bitwise version of int(x * 2**(52+m)).

   m : int
       Position of virtual binary point relative to matissa bits.

   Examples
   --------
   >>> x = fx(18.4, m=3)
   >>> x.d
   18.4
   >>> x.bin
   '0000100100110011001100110011001100110011001100110011001100000000'
   >>> "{:064b}".format(int(x.d * 2**(52+3)))
   '0000100100110011001100110011001100110011001100110011001100000000'

   """
   def __init__(self, x, m=0, bits=64):
       self.x = x                  # input value
       self.d = x                  # float value
       self.bits = bits
       self.m = m
       self.ix = ieee(x, bits=bits)
       if bits == 64:
          self.r = 11              # number of exponent bits
          self.p = 52              # number of mantissa bits
       else:
          self.r = 8
          self.p = 23
       i = self.ix.M + (1<<self.p) # add hidden bit
       m += self.ix.e              # normalise with exponent
       if m > 0: i <<= m           # scale with m
       elif m < 0: i >>= -m
       s = -self.ix.S
       i = s ^ s + i               # apply sign
       self.si = i                 # signed integer
       self.i = i & ((1<<self.bits)-1)    # make unsigned
       self.bin = "{:0{p.bits}b}".format(self.i, p=self)
   def __repr__(self):
       fmtbin = "{}: {}".format(self.x, self.bin)
       return fmtbin


def fx2float(x, m=0):
    """
    Floating point representation.

    Bitwise version of x * 2**-(52+m).

    x : int
        Fix point value. Must be python's signed integer.

    Example
    -------

    >>> x = fx(-.000000000007)
    >>> x; x.si/2.**52; fx2float(x.si)
    -7e-12: 1111111111111111111111111111111111111111111111111000010011011011
    -6.999956170261612e-12
    -6.999956170261612e-12
    >>> i = fx(-0.1, m=3)
    >>> i
    -0.1: 1111111111110011001100110011001100110011001100110011001100110011
    >>> fx2float(i.si, m=3)
    -0.1

    """
    r = 11
    p = 52
    s = x < 0               # m >> 63
    M = (x^-s) + s          # positive mantissa
    e = M.bit_length() - p - 1
    if e < 0: M <<= -e
    elif e > 0: M >>= e
    B = (1 << r-1) - 1      # exponent bias
    E = e - m + B
    i = (s<<63) | (E<<52) | M^(1<<52)   # removes hidden bit in M
    x = ctypes.c_double.from_buffer(ctypes.c_uint64(i)).value
    return x

