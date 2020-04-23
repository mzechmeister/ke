# Solving Kepler's equation with CORDIC double iteration
# by Mathias Zechmeister (2020-04-23)
# for gnuplot v5.2+

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
sar(x, n) = (_s=-(x>>31), (_s^x) >> n ^ _s)


# floating point version

f64_E(M, e) = (En=E0(M), ecn=K*e, esn=0.,\
   sum[ii=1:Nn] (\
       k = kseq[ii],\
       sgn = M>En-esn? 1:-1,\
       En = En + sgn*f64_a[k+1],\
       t = ecn,\
       ecn = ecn - sgn*esn/2**k,\
       esn = esn + sgn*t/2**k),\
   En)

f64_E_sgl(M, e) = (En=E0(M), ecn=K_sgl*e, esn=0.,\
   sum[ii=1:N] (\
       k = ii-1,\
       sgn = M>En-esn? 1:-1,\
       En = En + sgn*f64_a[k+1],\
       t = ecn,\
       ecn = ecn - sgn*esn/2**k,\
       esn = esn + sgn*t/2**k),\
   En)

# emulated fix point version

i32_E(M, e) = (z=fx(2*pi*floor(M/(2*pi)+0.5)-M), ecn=fx(K32*e), esn=0,\
   sum[ii=1:Nn32] (\
       k = kseq32[ii],\
       s = sar(esn-z, 31),\
       z = z + (s ^ s+i32_a[k+1]),\
       t = ecn,\
       ecn = ecn - (s ^ s+sar(esn,k)),\
       esn = esn + (s ^ s+sar(t,k)) ),\
   M+fp(esn))


plot f64_E(x,1.), '+' us ($1-sin($1)):1, f64_E_sgl(x,1.), i32_E(x,1.)

print K, kseq, i32_a

