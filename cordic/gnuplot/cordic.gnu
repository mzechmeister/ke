# CORDIC floating point for gnuplot
# author: Mathias Zechmeister
# version: 2020-03-24

if (!exists("N")) N = 32

N = 32

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

# v: vectoring
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
cordic_l(x, y, t) = cordic(x, y, t, 0, 0)    # x, y+xz, 0
cordic_c(x, y, t) = cordic(x, y, t, 1, 0)    # x cos z - y sin z, x sin z + y cos z, 0
cordic_h(x, y, t) = cordic(x, y, t, -1, 0)   # x cosh z + y sinh z, x sinh z + y cosh z, 0

# vectoring
cordic_vl(x, y, t) = cordic(x, y, t, 0, 1)   # x, y+xz, 0
cordic_vc(x, y, t) = cordic(x, y, t, 1, 1)   # x cos z - y sin z, x sin z + y cos z, 0
cordic_vh(x, y, t) = cordic(x, y, t, -1, 1)  # x cosh z + y sinh z, x sinh z + y cosh z, 0


# operators from linear mode
_mul(x, t) = (cordic_l(x, 0, t), Y)
_sqr(x) = (cordic_l(x, 0, x), Y)
_div(y, x) = (cordic_vl(x, y, 0), T)
_inv(x) = (cordic_vl(x, 1, 0), T)

# operators from circular mode
_cos(x) = (cordic_c(Kc, 0, x), X)
_sin(x) = (cordic_c(Kc, 0, x), Y)
_atan(x) = (cordic_vc(Kc, Kc*x, 0), T)
_atan2(y,x) = (cordic_vc(Kc*x, Kc*y, 0), T)
_r(x,y) = (cordic_vc(Kc*x, Kc*y, 0), X)

# operators from hyperbolic mode
_cosh(x) = (cordic_h(Kh, 0, x), X)
_sinh(x) = (cordic_h(Kh, 0, x), Y)
_exp(x) = (cordic_h(Kh, Kh, x), X)
_atanh(x) = (cordic_vh(Kh, Kh*x, 0), T)
_ln(x) = (cordic_vh(Kh*(x+1), Kh*(x-1), 0), 2*T)
_sqrt(x) = (cordic_vh(Kh*(x+0.25), Kh*(x-0.25), x), X)


