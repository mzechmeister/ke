# CORDIC floating point with gnuplot
#
# author: Mathias Zechmeister
# version: 2020-04-09
#
# This code emulates the functionality of CORDIC.
# It defines several elementary function without and with
# range extension.

if (!exists("_N")) _N = 32
N = _N

array k_[N*3]   # shift sequence
array a_[N*3]   # base angles

# scale correction
Kc = 1
Kh = 1

do for [i=1:N] {
   k_[i]     = k = i - (i>4) - (i>14)
   a_[i]     = atanh(2**-k)
   Kh        = Kh / sqrt(1-4**-k)
   k_[i+N]   = k = i - 1
   a_[i+N]   = 2**-k
   k_[i+N+N] = k = i - 1
   a_[i+N+N] = atan(2**-k)
   Kc        = Kc / sqrt(1+4**-k)
}

# Parameter: v - vectoring
cordic(x, y, t, m, v) = (\
    X = x,\
    Y = y,\
    T = t,\
    sum[i=1+_N*m+_N:N+_N*m+_N](\
       k = k_[i],\
       s = (v?-Y:T) < 0 ? -1. : 1.,\
       u = X,\
       X = X - m*s*Y/2.**k,\
       Y = Y +   s*u/2.**k,\
       T = T -   s*a_[i]),\
    k\
)

# rotating
rot_l(x, y, t) = cordic(x, y, t, 0, 0)    # x,  y+xt,  0
rot_c(x, y, t) = cordic(x, y, t, 1, 0)    # x cos t - y sin t,  x sin z + y cos t,  0
rot_h(x, y, t) = cordic(x, y, t, -1, 0)   # x cosh t + y sinh t,  x sinh z + y cosh t,  0

# vectoring
vec_l(x, y, t) = cordic(x, y, t, 0, 1)    # x,  0,  t + y/x
vec_c(x, y, t) = cordic(x, y, t, 1, 1)    # sqrt(x^2 + y^2),  0,  t + atan(y,x)
vec_h(x, y, t) = cordic(x, y, t, -1, 1)   # sqrt(x^2 - y^2),  0,  t + atanh(y,x)


# operators from linear mode
_mul(x, t)    = (rot_l(x, 0, t), Y)
_sqr(x)       = (rot_l(x, 0, x), Y)
_div(y, x)    = (vec_l(x, y, 0), T)
_inv(x)       = (vec_l(x, 1, 0), T)

# operators from circular mode
_cos(x)       = (rot_c(Kc, 0, x), X)
_sin(x)       = (rot_c(Kc, 0, x), Y)
_hypot(x, y)  = (vec_c(x, y, 0), Kc*X)   # hypotenuse sqrt(x^2+y^2)
_atan2(y, x)  = (vec_c(x, y, 0), T)
_atan(x)      = (vec_c(1, x, 0), T)
_acot(x)      = (vec_c(x, 1, 0), T)

# operators from hyperbolic mode
_cosh(x)      = (rot_h(Kh, 0, x), X)
_sinh(x)      = (rot_h(Kh, 0, x), Y)
_exp(x)       = (rot_h(Kh, Kh, x), X)
_cath(x,y)    = (vec_h(x, y, 0), Kh*X)   # cathetus  sqrt(x^2-y^2)
_atanh2(y,x)  = (vec_h(x, y, 0), T)
_atanh(x)     = (vec_h(1, x, 0), T)
_acoth(x)     = (vec_h(x, 1, 0), T)
_ln(x)        = (vec_h(x+1, x-1, 0), 2*T)
_sqrt(x)      = (vec_h(x+Kh**2/4, x-Kh**2/4, 0), X)


# With range extension
ROT_l(x, y, t) = (E=floor(abs(t)/log(2)), rot_l(x*2**E, y, t/2.**E))

ROT_c(x, y, t) = (m=floor(t/pi*2), Q=(m%4+4)%4, rot_c(Q==0?x:Q==1?-y:Q==2?-x:y, Q==0?y:Q==1?x:Q==2?-y:-x, t-pi/2*m))  # Walther (1971) four quadrants
ROT_c(x, y, t) = (m=floor(t/pi+0.5), Q=m&1, rot_c(Q?-x:x, Q?-y:y, t-pi*m))   # two branches
ROT_c(x, y, t) = (T=t-2*pi*floor(t/(2*pi))-pi, s=T<0?-1:1, rot_c(s*y, -s*x, T-s*pi/2))   # two branches
 
ROT_h(x, y, t) = (E=floor(t/log(2)), D=t-E*log(2), rot_h(x*(2**E+2**-E)/2 + y*(2**E-2**-E)/2, x*(2**E-2**-E)/2 + y*(2**E+2**-E)/2, D))
ROT_h(x, y, t) = (E=floor(t/log(2)), D=t-E*log(2), rot_h(((x+y)*2**E + (x-y)*2**-E)/2, ((x+y)*2**E - (x-y)*2**-E)/2, D))
ROT_h(x, y, t) = (E=floor(t/log(2)), D=t-E*log(2), u=(x+y)*2**E/2, v=(x-y)*2**-E/2, rot_h(u+v, u-v, D))


VEC_l(x, y, t) = (S=x<0?-1:1, Ex=floor(log(abs(x))/log(2)), Ey=floor(log(y)/log(2)), vec_l(S*x/2.**Ex, S*y/2**Ey, 0), T=t+T/2**(Ex-Ey))
VEC_l(x, y, t) = (S=x<0?1.:-1., Ex=floor(log(abs(x))/log(2)), Ey=floor(log(y)/log(2)), vec_l(-S*x, -S*y-x/2**(Ex-Ey), t-S/2**(Ex-Ey)))
VEC_l(x, y, t) = (S=x<0?-1:1, E=floor(log(abs(x))/log(2))-floor(log(y)/log(2)), vec_l(S*x/2.**E, S*y, 0), T=t+T/2**E)

VEC_c(x, y, t) = (Q=x<0, vec_c(Q?-x:x, Q?-y:y, t+pi*Q*(y<0?-1:1)))

VEC_h(x, y, t) = (S=x<0?-1:1, u=1-S*x, E=floor(log(u)/log(2)), M=u/2**E, vec_h(2+M-u, 2-M-u, 0), S*(T-log(2)*E/2), X*2**(E/2.)*Kh/2)
VEC_h(x, y, t) = (Ex=floor(log(x)/log(2)), vec_h((x+y)*2**E + (x-y)*2**-E, (x+y)*2**E - (x-y)*2**-E, t))
VEC_h(x, y, t) = (u=x-y, E=floor(log(u)/log(2)), M=u/2**E, vec_h(2*x+M-u, 2*x-M-u, 0))   # works (but not left)
VEC_h(x, y, t) = (E=floor(log(x-y)/log(2)), u=(x+y)*2**(E/2.)/2, v=(x-y)*2**(-E/2.)/2, vec_h(u+v, u-v, t-E/2.*log(2)))   # works
VEC_h(x, y, t) = (S=y<0?-1:1, E=S*floor(log(x-S*y)/log(2)), u=(x+y)*2**(E/2.)/2, v=(x-y)*2**(-E/2.)/2, vec_h(u+v, u-v, t-E/2.*log(2)))  # works both y but not for x<0
VEC_h(x, y, t) = (S=y<0 ^ x<0 ?-1:1, E=floor(log(abs(x)-abs(y))/log(4)), u=(abs(x)+abs(y))*2**E/2, v=(abs(x)-abs(y))*2**-E/2, vec_h(u+v, S*(u-v), t-S*E*log(2)))
VEC_h(x, y, t) = (S=y<0 ^ x<0 ?-1:1, E=floor((log(abs(x)-abs(y))-log(abs(x)))/log(4)), u=(abs(x)+abs(y))*2**E/2, v=(abs(x)-abs(y))*2**-E/2, vec_h(u+v, S*(u-v), t-S*E*log(2)))


MUL(x, t)    = (ROT_l(x, 0, t), Y)
SQR(x)       = (ROT_l(x, 0, x), Y)
DIV(y, x)    = (VEC_l(x, y, 0), T)
INV(x)       = (VEC_l(x, 1, 0), T)

COS(x)       = (ROT_c(Kc, 0, x), X)
SIN(x)       = (ROT_c(Kc, 0, x), Y)
HYPOT(x, y)  = (VEC_c(x, y, 0), Kc*X)
ATAN2(y, x)  = (VEC_c(x, y, 0), T)
ATAN(x)      = (VEC_c(1, x, 0), T)
ACOT(x)      = (VEC_c(x, 1, 0), T)

COSH(x)      = (ROT_h(Kh, 0, x), X)
SINH(x)      = (ROT_h(Kh, 0, x), Y)
EXP(x)       = (ROT_h(Kh, Kh, x), X)
CATH(x, y)   = (VEC_h(x, y, 0), Kh*X)   # cathetus  sqrt(x^2-y^2)
ATANH2(y, x) = (VEC_h(x, y, 0), T)
ATANH(x)     = (VEC_h(1, x, 0), T)
ACOTH(x)     = (VEC_h(x, 1, 0), T)
LN(x)        = (E=floor(log(x)/log(2)), _ln(x/2**E)+E*log(2))   # Walther (1971)
LN(x)        = (VEC_h(x+1, x-1, 0), 2*T)
SQRT(x)      = (E=floor(log(x)/log(2)), _sqrt(x/2**E)*2**(E/2.))   # Walther (1971)
SQRT(x)      = (VEC_h(x+Kh**2/4, x-Kh**2/4, 0), X)

# Operations derived from two CORDIC operations

_asin(x)   = _atan2(x, _cath(1, x))     # atan(x/sqrt(1-x^2))
_acos(x)   = _atan2(_cath(1, x), x)     # atan(sqrt(1-x^2)/x)
_tan(x)    = _div((rot_c(1, 0, x), Y), X)
_cot(x)    = _div((rot_c(1, 0, x), X), Y)
_asinh(x)  = _atanh2(x, _hypot(1, x))   # atanh(x/sqrt(x^2+1)) = ln(x+sqrt(x^2+1))
_acosh(x)  = _atanh2(_cath(x,1), x)     # atanh(sqrt(x^2-1)/x) = ln(x+sqrt(x^2-1))
_tanh(x)   = _div((rot_h(1, 0, x), Y), X)
_coth(x)   = _div((rot_h(1, 0, x), X), Y)

ASIN(x)    = ATAN2(x, CATH(1, x))
ACOS(x)    = ATAN2(CATH(1, x), x)
TAN(x)     = DIV((ROT_c(1, 0, x), Y), X)
COT(x)     = DIV((ROT_c(1, 0, x), X), Y)
ASINH(x)   = ATANH2(x, HYPOT(1,x))
ACOSH(x)   = ATANH2(CATH(x, 1), x)
TANH(x)    = DIV((ROT_h(1, 0, x), Y), X)
COTH(x)    = DIV((ROT_h(1, 0, x), X), Y)


# Operations derived from three CORDIC operations

POW(x, y)  = EXP(MUL(LN(x),y))
ROOT(x, y) = EXP(DIV(LN(x),y))



