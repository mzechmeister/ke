#! /usr/bin/gnuplot
# This generate Fig. A.2 in the cordic-dbl paper.

load 'cordic_asin.gnu'

set term postscript enh por col 10 size 10cm,6cm
set out "demo_asin.ps"

set key Left left rev title "CORDIC arcsine                        " spacing 1.2 samplen 2
set mxtics
set mytics
set rmar 0.5
set tmar 0.3
set lmar 4


# pl [0:1.1][-1.9:1.9] asin(x), _asin(x) # needs range extension
pl [0:1.1][0:1.9] asin(x) lc 7 lw 2,\
     _asin_sgl(x)   lc 1 lw 1.8 t "k_n = 0,1,2,3,... (single rotation)",\
     keyentry   w l lc 2 lw 1.2 t "k_n = 0,0,1,1,2,2,...",\
     keyentry   w l lc 5 lw 1.4 t "k_n = 0,1,1,2,2,...",\
     _asin_dbl3(x)  lc 3 lw 1.6 t "k_n = 1,1,2,2,...",\
     _asin_dbl2(x)  lc 5 lw 1.4 t "",\
     _asin_dbl1(x)  lc 2 lw 1.2 t ""

set out



