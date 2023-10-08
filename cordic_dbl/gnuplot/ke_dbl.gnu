# Solving Kepler's equation with CORDIC double iteration
# by Mathias Zechmeister (2020-04-23)
# for gnuplot v5.2+ (provides array)

if (GPVAL_VERSION >= 5.4) {
   unset overflow
}

N = 50
N32 = 30

# set up lookup tables for circular mode

array kseq[N*1.6]
array kseq32[N32*1.6]
array i32_a[N32]
array f64_a[N]

# scale correction factor
K = 1.
K_sgl = 1.
K32 = 1.

n = 0
do for [k=0:N-1]{
   n = n + 1
   kseq[n] = k
   f64_a[k+1] = atan(2**-k)
   sec2a = 1 + 4**-k
   K_sgl = K_sgl / sec2a**.5
   if (sec2a>1) {
      n = n + 1
      K = K / sec2a
      kseq[n] = k
   }
}
Nn = n


# overflow for |x| > 4, so pi is included
fx(x) = int(x*(1 << 29))  # fix point
fp(x) = x*1./(1 << 29)    # floating point


n = 0
do for [k=0:N32-1]{
   n = n + 1
   kseq32[n] = k
   i32_a[k+1] = fx(atan(2**-k))
   if (2*k<23) {
      n = n + 1
      K32 = K32 / (1 + 4**-k)
      kseq32[n] = k
   }
}
Nn32 = n

# phase folding (rem2pi)
E0(M) = 2 * pi * floor(M/(2*pi)+0.5)

# function to emulate arithmetic right shift, since gnuplot has logical right shift
sar(x, n) = (_s=-(x<0), (_s^x) >> n ^ _s)


# floating point version

f64_E(M, e) = (En=E0(M), X=K*e, Y=0.,\
   sum[ii=1:Nn] (\
       k = kseq[ii],\
       sgn = M>En-Y? 1:-1,\
       En = En + sgn*f64_a[k+1],\
       u = X,\
       X = X - sgn*Y/2**k,\
       Y = Y + sgn*u/2**k),\
   En)


# single rotation version (just to demo non-converging!)

f64_E_sgl(M, e) = (En=E0(M), X=K_sgl*e, Y=0.,\
   sum[ii=1:N] (\
       k = ii-1,\
       sgn = M>En-Y? 1:-1,\
       En = En + sgn*f64_a[k+1],\
       u = X,\
       X = X - sgn*Y/2**k,\
       Y = Y + sgn*u/2**k),\
   En)


# emulated fix point version

i32_E(M, e) = (\
   T = fx(M-2*pi*floor(M/(2*pi)+0.5)),\
   X = fx(K32*e),\
   Y = 0,\
   sum[ii=1:Nn32] (\
       k = kseq32[ii],\
       s = sar(Y+T, 31),\
       T = T - (s ^ s+i32_a[k+1]),\
       u = X,\
       X = X - (s ^ s+sar(Y,k)),\
       Y = Y + (s ^ s+sar(u,k)) ),\
   M + fp(Y))

