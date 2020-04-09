# Demo of CORDIC for gnuplot

load "cordic.gnu"

# output last iteration
print "k=",rot_h(Kh,0,1), "   X=", X, "   Y=",Y, "  T=", T, "  cosh(1)=", cosh(1)

# output the iteration sequence
do for [N=1:12] {print sprintf("k=%2d X=%- 16.12f Y=%- 16.12f T=%- 16.12f ", rot_h(1,0,1), X,Y,T), cosh(1)/Kh}

N = 32
plot [-1.5:1.5] _mul(1.3,x)-1.3*x
plot [-1.5:1.5] _cos(x)-cos(x), _sin(x)-sin(x)
plot [-1.1:1.1] _cosh(x)-cosh(x), _sinh(x)-sinh(x), _exp(x)-exp(x)
pl _atan(x)-atan(x)
print atan2(0.5,0.7), _atan2(0.5,0.7)
plot _hypot(1,x)-sqrt(1+x**2)
plot [0.027:2] _sqrt(x)-sqrt(x)
plot [0.12:2] _ln(x)-log(x)
plot [2.5:] _div(1.5, x), 1.5/x
plot [-2:2] _sqr(x) - x**2

plot _asin(x), ASIN(x), asin(x), _acos(x), ACOS(x), acos(x)
plot _tan(x), TAN(x), tan(x), _cot(x), COT(x), 1/tan(x)

plot _asinh(x), ASINH(x), asinh(x), _acosh(x), ACOSH(x), acosh(x)

set par
plot for [N=1:6] COS(t),SIN(t)


