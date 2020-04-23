# Solving Kepler's equation with CORDIC double iteration
# by Mathias Zechmeister (2020-04-23)
# for gnuplot v5.2+

N = 50
N32 = 30

# set up lookup tables for circular mode

array nseq[N*1.6]
array nseq32[N32*1.6]
array i32_a[N32]
array f64_a[N]

# scale correction factor
K = 1.
K_nodbl = 1.
K32 = 1.

i = 0
do for [n=0:N-1]{
   i = i + 1
   nseq[i] = n
   f64_a[n+1] = atan(2**-n)
   sec2a = 1 + 4**-n
   K_nodbl = K_nodbl / sec2a**.5
   if (sec2a>1) {
      i = i + 1
      K = K / sec2a
      nseq[i] = n
   }
}
Nn = i


# overflow for |x| > 4, so pi is included
fx(x) = int(x*(1 << 29))  # fix point
fp(x) = x*1./(1 << 29)    # floating point


i = 0
do for [n=0:N32-1]{
   i = i + 1
   nseq32[i] = n
   i32_a[n+1] = fx(atan(2**-n))
   if (2*n<23) {
      i = i + 1
      K32 = K32 / (1 + 4**-n)
      nseq32[i] = n
   }
}
Nn32 = i

# phase folding (rem2pi)
E0(M) = 2 * pi * floor(M/(2*pi)+0.5)

# function to emulate arithmetic right shift, since gnuplot has logical right shift
sar(x, n) = (_s=-(x>>31), (_s^x) >> n ^ _s)


# floating point version

f64_E(M, e) = (En=E0(M), ecn=K*e, esn=0.,\
   sum[ii=1:Nn] (\
       n = nseq[ii],\
       sgn = M>En-esn? 1:-1,\
       En = En + sgn*f64_a[n+1],\
       t = ecn,\
       ecn = ecn - sgn*esn/2**n,\
       esn = esn + sgn*t/2**n),\
   En)

f64_E_sgl(M, e) = (En=E0(M), ecn=K_nodbl*e, esn=0.,\
   sum[ii=1:N] (\
       n = ii-1,\
       sgn = M>En-esn? 1:-1,\
       En = En + sgn*f64_a[n+1],\
       t = ecn,\
       ecn = ecn - sgn*esn/2**n,\
       esn = esn + sgn*t/2**n),\
   En)

# emulated fix point version

i32_E(M, e) = (z=fx(2*pi*floor(M/(2*pi)+0.5)-M), ecn=fx(K32*e), esn=0,\
   sum[ii=1:Nn32] (\
       n = nseq32[ii],\
       s = sar(esn-z, 31),\
       z = z + (s ^ s+i32_a[n+1]),\
       t = ecn,\
       ecn = ecn - (s ^ s+sar(esn,n)),\
       esn = esn + (s ^ s+sar(t,n)) ),\
   M+fp(esn))


plot f64_E(x,1.), '+' us ($1-sin($1)):1, f64_E_sgl(x,1.), i32_E(x,1.)

print K, nseq, i32_a

