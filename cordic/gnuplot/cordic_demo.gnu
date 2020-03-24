# demo
print cordic_h(1,0,1), X,Y,T, cosh(1)
do for [N=1:2] {print cordic_h(1,0,1), ' x:',X,Y,T, cosh(1)}

plot [-1.5:1.5] _mul(1.3,x)-1.3*x
plot [-1.5:1.5] _cos(x)-cos(x), _sin(x)-sin(x)
plot [-1.1:1.1] _cosh(x)-cosh(x), _sinh(x)-sinh(x), _exp(x)-exp(x)
pl _atan(x)-atan(x)
print atan2(0.5,0.7), _atan2(0.5,0.7)
plot _r(1,x)- sqrt(1+x**2)
plot _sqrt(x)- sqrt(x)
plot _ln(x) - log(x)
plot [2.5:] _div(5, x) - 5/x
plot [-2:2]_sqr(x) - x**2

