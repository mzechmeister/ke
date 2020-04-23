To plot Kepler's equation in gnuplot run

```gnuplot
load "ke.gnu"
plot [0:pi] _E(x, 1., 29)
```

And to compare with E-e*sinE, E:
```gnuplot
plot [0:pi] _E(x, 1., 29), "+" us (M($1)):1 w l t "E-e*sinE,E"
```
