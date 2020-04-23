# CORDIC floating point with gnuplot
#
# author: Mathias Zechmeister
# version: 2020-04-09
#
# This code compares various attemps to compute the arcsine.

load "cordic.gnu"   # provides already _asin(x) from subsequent CORDIC operations!


# cordic_sgl is an enhanced cordic variant with the arcsine mode.
# It has an additional keyword "arg".

cordic_sgl(x, y, t, m, v, arg) = (\
    X = x,\
    Y = y,\
    T = t,\
    sum[i=1+_N*m+_N:N+_N*m+_N](\
       k = k_[i],\
       s = (v?arg-Y:T) < 0 ? -1. : 1.,\
       u = X,\
       X = X - m*s*Y/2.**k,\
       Y = Y +   s*u/2.**k,\
       T = T -   s*a_[i]),\
    k\
)

_asin_sgl(x) = (cordic_sgl(Kc, 0, 0, 1, 1, x), -T)   # single iterations are not accurate!
_cath_sgl(x) = (cordic_sgl(Kc, 0, 0, 1, 1, x), X)    # single iterations are not accurate!
#pl _asin_sgl(x)
#pl (_asin_sgl(x), x-T)

# plot [-1.2:1.2] asin(x), _asin_sgl(x), sqrt(1-x**2), _cath_sgl(x)



# CORDIC with double rotation and on the fly scale correction.
#
#  Reference:
#
#    Jean-Michel Muller,
#    Elementary Functions: Algorithms and Implementation,
#    Second Edition,
#    Birkhaeuser, 2006,
#    ISBN13: 978-0-8176-4372-0,
#    LC: QA331.M866.

cordic_dblw(x, y, t, m, v, arg) = (\
    X = x,\
    Y = y,\
    T = t,\
    w = arg,\
    sum[i=1+N*m+N:N+N*m+N](\
       k = k_[i],\
       s = (v?w-Y:T) < 0 ? -1. : 1.,\
       u = X,\
       X = X - m*s*Y/2.**k,\
       Y = Y +   s*u/2.**k,\
       u = X,\
       X = X - m*s*Y/2.**k,\
       Y = Y +   s*u/2.**k,\
       w = w + w/2.**(2*k),\
       T = T - s*2*a_[i]),\
    k\
)

_asin_dblw(x) = (cordic_dblw(1, 0, 0, 1, 1, x), -T)   # double iterations are accurate
_cath_dblw(x) = (cordic_dblw(Kc**2, 0, 0, 1, 1, Kc**2*x), X)    # double iterations are accurate
_cath_dblw(x) = (cordic_dblw(1, 0, 0, 1, 1, x), X*Kc**2)    # double iterations are accurate

# plot [-1.2:1.2] asin(x), _asin_dblw(x), sqrt(1-x**2), _cath_dblw(x)



# CORDIC with double rotation, but without scale correction.

Ndbl = 2 * _N

array kdbl[_N*3+Ndbl*4]   # shift sequence
array adbl[_N*3+Ndbl*4]   # base angles

N1 = _N*3
N2 = N1 + Ndbl
N3 = N2 + Ndbl
N4 = N3 + Ndbl

# scale correction
Kdbl1 = 1   #  => Kc^2
Kdbl2 = 1   #  => sqrt(2) * Kc^2
Kdbl3 = 1   #  => 2 * Kc^2
Kdbl4 = 1


do for [i=1:Ndbl] {
   kdbl[i+N1] = k = (i-1) / 2       # 0,0,1,1,2,2,3,3,4,4,5,5
   adbl[i+N1] = atan(2**-k)
   Kdbl1      = Kdbl1 / sqrt(1+4**-k)
   kdbl[i+N2] = k = (i  ) / 2       # 0,1,1,2,2,3,3,4,4,5,5
   adbl[i+N2] = atan(2**-k)
   Kdbl2        = Kdbl2 / sqrt(1+4**-k)
   kdbl[i+N3] = k = (i+1) / 2       # 1,1,2,2,3,3,4,4,5,5
   adbl[i+N3] = atan(2**-k)
   Kdbl3      = Kdbl3 / sqrt(1+4**-k)
   kdbl[i+N4] = k = (i - 1) - floor(2*(i-3)/5)   # [0,1,2,3,4,4,5,5,6,7,7,8,8,9,10,10,11,11,12,13,13,14,14,15,16,16,17,17,18,19,19,20]
   adbl[i+N4] = atan(2**-k)
   Kdbl4      = Kdbl4 / sqrt(1+4**-k)
}

cordic_dbl(x, y, t, m, v, arg) = (\
    X = x,\
    Y = y,\
    T = t,\
    sum[i=1+N0:N+N0](\
       k = kdbl[i],\
       s = (v?arg-Y:T) < 0 ? -1. : 1.,\
       s = v!=2? s : X<0?-1:s,\
       u = X,\
       X = X - m*s*Y/2.**k,\
       Y = Y +   s*u/2.**k,\
       T = T - s*adbl[i]),\
    k\
)

_asin_dbl1(x) = (N0=N1, cordic_dbl(Kdbl1, 0, 0, 1, 1, x), -T)   # convergence range: |x| < 2*Kc**2 = 0.737512
_asin_dbl2(x) = (N0=N2, cordic_dbl(Kdbl2, 0, 0, 1, 1, x), -T)   # convergence range: |x| < 0.91262 # < 2*(1+4**-1)*Kc**2 = 0.92189
_asin_dbl3(x) = (N0=N3, cordic_dbl(Kdbl3, 0, 0, 1, 1, x), -T)   # convergence range: |x| < 0.9867,0.9903
_asin_dbl4(x) = (N0=N4, cordic_dbl(Kdbl4, 0, 0, 1, 1, x), -T)   # convergence range: |x| < 0.9938

# plot asin(x), _asin_dbl1(x), _asin_dbl2(x), _asin_dbl3(x), _asin_dbl4(x)


# extend convergence range with check to prevent X<0 (mode: v=2)
_Asin_dbl1(x) = (N0=N1, cordic_dbl(Kdbl1, 0, 0, 1, 2, x), -T)   # convergence range: |x| < 2*Kc**2 = 0.737512
_Asin_dbl2(x) = (N0=N2, cordic_dbl(Kdbl2, 0, 0, 1, 2, x), -T)   # convergence range: |x| < 0.91262 # < 2*(1+4**-1)*Kc**2 = 0.92189
_Asin_dbl3(x) = (N0=N3, cordic_dbl(Kdbl3, 0, 0, 1, 2, x), -T)   # convergence range: |x| < 0.9867,0.9903
_Asin_dbl4(x) = (N0=N4, cordic_dbl(Kdbl4, 0, 0, 1, 2, x), -T)   # convergence range: |x| < 0.9938

# plot asin(x), _Asin_dbl1(x), _Asin_dbl2(x), _Asin_dbl3(x), _Asin_dbl4(x)


# special sqrt versions

# simplified hyperbolic version without t channel
_sqrt_xy(x) = (\
    X = x + Kh**2*0.25,\
    Y = x - Kh**2*0.25,\
    sum[i=1:N](\
       k = k_[i],\
       s = -Y<0 ? -1. : 1.,\
       u = X,\
       X = X + s*Y/2.**k,\
       Y = Y + s*u/2.**k),\
    X\
)

# Burkardt version (not really CORDIC)
_sqrt_bk(x) = (\
    X = x,\
    sum[i=1:N](\
       p = 1/2.**i,\
       s = (X+p)*(X+p)<x,\
       X = X + s*p),\
    X\
)

# plot [-0.5:5] sqrt(x), _sqrt_xy(x), _sqrt_bk(x)

