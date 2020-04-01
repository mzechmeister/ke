# CORDIC floating point with gnuplot
#
# author: Mathias Zechmeister
# version: 2020-04-01
#
# This code emulates the functionality of CORDIC.
# It defines several elementary function without and with
# range extension.

if (!exists("N")) N = 32

# shift sequence
array kl[N]
array kc[N]
array kh[N]

# base angles
array al[N]
array ac[N]
array ah[N]

# scale correction
Kl = 1
Kc = 1
Kh = 1

do for [i=1:N] {
   kl[i] = i - 1
   kc[i] = i - 1
   kh[i] = i - (i>4) - (i>14)
   al[i] = 2**-kc[i]
   ac[i] = atan(2**-kc[i])
   ah[i] = atanh(2**-kh[i])
   Kc = Kc * (1+4**-kc[i])**-0.5
   Kh = Kh * (1-4**-kh[i])**-0.5
}

# Parameter: v - vectoring
cordic(x, y, t, m, v) = (\
    X = x,\
    Y = y,\
    T = t,\
    sum[i=1:N](\
       k = m==1? kc[i] : m? kh[i] : kl[i],\
       a = m==1? ac[i] : m? ah[i] : al[i],\
       u = v? -Y: T,\
       s = u<0 ? -1. : 1.,\
       u = X,\
       X = X - m*s*Y/2.**k,\
       Y = Y +   s*u/2.**k,\
       T = T - s*a),\
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
_mul(x, t)   = (rot_l(x, 0, t), Y)
_sqr(x)      = (rot_l(x, 0, x), Y)
_div(y, x)   = (vec_l(x, y, 0), T)
_inv(x)      = (vec_l(x, 1, 0), T)

# operators from circular mode
_cos(x)      = (rot_c(Kc, 0, x), X)
_sin(x)      = (rot_c(Kc, 0, x), Y)
_hypot(x,y)  = (vec_c(Kc*x, Kc*y, 0), X)   # hypotenuse sqrt(x^2+y^2)
_atan2(y,x)  = (vec_c(Kc*x, Kc*y, 0), T)
_atan(x)     = (vec_c(Kc, Kc*x, 0), T)
_acot(x)     = (vec_c(Kc*x, Kc, 0), T)

# operators from hyperbolic mode
_cosh(x)     = (rot_h(Kh, 0, x), X)
_sinh(x)     = (rot_h(Kh, 0, x), Y)
_exp(x)      = (rot_h(Kh, Kh, x), X)
_cath(x,y)   = (vec_h(x, y, 0), Kh*X)   # cathetus  sqrt(x^2-y^2)
_atanh2(y,x) = (vec_h(x, y, 0), T)
_atanh(x)    = (vec_h(1, x, 0), T)
_acoth(x)    = (vec_h(x, 1, 0), T)
_ln(x)       = (vec_h(x+1, x-1, 0), 2*T)
_sqrt(x)     = (vec_h(x+0.25, x-0.25, 0), Kh*X)


# With range extension
ROT_l(x, y, t) = (E=floor(abs(t)/log(2)), rot_l(x*2**E, y, t/2**E))

ROT_c(x, y, t) = (m=floor(t/pi*2), Q=(m%4+4)%4, rot_c(Q==0?x:Q==1?-y:Q==2?-x:y, Q==0?y:Q==1?x:Q==2?-y:-x, t-pi/2*m))  # Walther (1971) four quadrants
ROT_c(x, y, t) = (m=floor(t/pi+0.5), Q=m&1, rot_c(Q?-x:x, Q?-y:y, t-pi*m))   # two branches
ROT_c(x, y, t) = (m=floor(t/(2*pi)+0.5), Q=t-2*pi*m<0, rot_c(Q?y:-y, Q?-x:x, t-2*pi*m+(Q?pi/2:-pi/2)))   # two branches

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


MUL(x, t)   = (ROT_l(x, 0, t), Y)
SQR(x)      = (ROT_l(x, 0, x), Y)
DIV(y, x)   = (VEC_l(x, y, 0), T)
INV(x)      = (VEC_l(x, 1, 0), T)

COS(x)      = (ROT_c(Kc, 0, x), X)
SIN(x)      = (ROT_c(Kc, 0, x), Y)
HYPOT(x,y)  = (VEC_c(Kc*x, Kc*y, 0), X)
ATAN2(y,x)  = (VEC_c(Kc*x, Kc*y, 0), T)
ATAN(x)     = (VEC_c(Kc, Kc*x, 0), T)
ACOT(x)     = (VEC_c(Kc*x, Kc, 0), T)

COSH(x)     = (ROT_h(Kh, 0, x), X)
SINH(x)     = (ROT_h(Kh, 0, x), Y)
EXP(x)      = (ROT_h(Kh, Kh, x), X)
CATH(x,y)   = (VEC_h(x, y, 0), Kh*X)   # cathetus  sqrt(x^2-y^2)
ATANH2(y,x) = (VEC_h(x, y, 0), T)
ATANH(x)    = (VEC_h(1, x, 0), T)
ACOTH(x)    = (VEC_h(x, 1, 0), T)
LN(x)       = (E=floor(log(x)/log(2)), M=x/2**E, vec_h(M+1, M-1, 0), 2*T+E*log(2))
LN(x)       = (E=floor(log(x)/log(2)), _ln(x/2**E)+E*log(2))   # Walther (1971)
LN(x)       = (VEC_h(x+1, x-1, 0), 2*T)
SQRT(x)     = (E=floor(log(x)/log(2)), _sqrt(x/2**E)*2**(E/2.))   # Walther (1971)
SQRT(x)     = (VEC_h(x+0.25, x-0.25, 0), Kh*X)


