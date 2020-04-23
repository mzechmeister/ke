load "ke.gnu"

set parametric

N = 15
e = 1.0

# define the algorithm to solve for E
E(M,e) = _Epn(M,e,N)   # twosided, half-angle

print Ms=pi/2, E(Ms,1.0), E(Ms,1.0)-sin(E(Ms,1.0))

# plot [-pi:pi] M(t),t, M(t),E(M(t),e),  M(t),(E(M(t),e)-t)*1e11

c0 = '#a2142f' # red
c1 = '#d95319' # orange
c2 = '#edb120' # yellow
c3 = '#77ac30' # green
c4 = '#4dbeee' # light-blue
c5 = '#0072bd' # blue
c6 = '#7e2f8e' # purple

set palette defined (0 c0, 1 c1, 2 c2, 3 c3, 4 c4, 5 c5) maxcolors 6

set macro

min(a,b) = a<b? a:b
max(a,b) = a>b? a:b


#E(M,N) = (M0=floor(M/(2*pi)+0.5)*2*pi, _E(M,e,N,1,M0,{1,0}))
E(M,N) = _E(M,e,N,1,M0,{1,0})
E(M,N) = _Epn(M,e,int(N))

# define the point to be shown
Eg = 2.0
g = Eg - sin(Eg) #1.0
#g=1; Eg=E(g,40)
print  Eg,g

set term post enh por sol col 12 size 10cm, 8cm
set out "ke_cordic.ps"

set xrange [0:pi]
set yrange [0:pi]
set cbrange [-0.5:5.5]
set cbtics scale 0 off -2.2,0
set xlabel "mean anomaly M"
set ylabel "eccentric anomaly E" off 1.5,0
set cblabel "iteration n" off -1,0
set tmar 0.3
#set rmar 0.95 does not work
set rmar at screen 0.88

set xtics (0, pi/8 1, "@^{/Symbol p}@\261_4" pi/4, 3*pi/8 1, "@^{/Symbol p}@\261_2" pi/2, 5*pi/8 1, "@^{3{/Symbol p}}{&{.}@\261_4}" 3*pi/4, 7*pi/8 1, "{/Symbol p}" pi)
set ytics (0, pi/8 1, "@^{/Symbol p}@\261_4" pi/4, 3*pi/8 1, "@^{/Symbol p}@\261_2" pi/2, 5*pi/8 1, "@^{3{/Symbol p}}{&{.}@\261_4}" 3*pi/4, 7*pi/8 1, "{/Symbol p}" pi)

set key bottom spacing 2 samplen 3

plot [0:pi] M(t),t lt 7 lw 2 t "M(E), E",\
     M(t),E(M(t),3) lt 1 lc rgb c3 t "E_3(M)",\
     M(t),E(M(t),5) lt 1 lc rgb c5 t "E_5(M)",\
     "<seq 6" us (E(g,$0)-sin(E(g,$0))):(E(g,$0)):0 pt 7 palette t "E_n(M_n)",\
     "<echo M E" us (g):(Eg) pt 7 lc 7 t "M=2-sin 2, E=2", '' us (g):(pi) w i lc 7 t""


#pl "<seq 5" us (g):(1):(g-0.05):(g+0.05):(E(g,$1-1)):(E(g,$1)):1 w boxxy lc var fill t "{/Symbol a}_n" boxxy +palette+4.6 not recognised
set out
print "ke_cordic.ps"

reset
set par
x1 = -1.3; x2 = 1.3
y1 = -0.1; y2 = 1.4
ratio = (y2-y1) / (x2-x1)
set term post enh por dash col 10 size 10cm, (10*ratio)cm
set out "cordic_path.ps"

unset border
unset xtics
unset ytics
#set rmar at screen 0.91
#set lmar at screen 0.03
set tmar at screen 0.95
set bmar 0.0
set tmar 0.0
set rmar 0
set lmar 0.
#set size squ
set size 1.,1.; set ori 0,0

set xrange [x1:x2]
set yrange [y1:y2]

set trange [0:1]
#set cbrange [:7] # let black for E_inf

unset colorbox
unset key

set palette defined (0 c0, 1 c1, 2 c2, 3 c3, 4 c4, 5 c5)

r = 1.19
rM = 1.09
set style fill solid .8 #border 0
rc(n) = 0.1 + n*0.05

set macro
arcrng = "min(E(g,n),E(g,n-1))/pi*180 : max(E(g,n),E(g,n-1))/pi*180"
n=1; set obj 8-n circle size rc(n) arc [ @arcrng ] fc rgb c1
n=2; set obj 8-n circle size rc(n) arc [ @arcrng ] fc rgb c2
n=3; set obj 8-n circle size rc(n) arc [ @arcrng ] fc rgb c3
n=4; set obj 8-n circle size rc(n) arc [ @arcrng ] fc rgb c4
n=5; set obj 8-n circle size rc(n) arc [ @arcrng ] fc rgb c5

pl cos(2*pi*t),sin(2*pi*t) lt 1 lc 9,\
   "<seq 6" us (0):(0):(cos(E(g,$0))):(sin(E(g,$0))):0 w vec filled palette lt 1,\
   "" us (r*cos(E(g,$0))):(r*sin(E(g,$0))):(sprintf("E_%d",$0)):0 w labels tc palette,\
   "" us (0):(0):(cos(M(E(g,$0)))):(sin(M(E(g,$0)))):0 w vec nohead palett dt 2,\
   "" us (cos(M(E(g,$0)))):(sin(M(E(g,$0)))):0 pt 6 palett,\
   "" us (rM*(cos(M(E(g,$0))))):(rM*(sin(M(E(g,$0))))):(sprintf("M_%d",$0)):0 w labels tc palette,\
   "<echo M" us (0):(0):(1.3*cos(g)):(1.3*sin(g)) w vec nohead lt 1 lc 7,\
   "" us (1.4*cos(g)):(1.4*sin(g)):1 w labels,\
   "<echo E" us (0):(0):(1.3*cos(Eg)):(1.3*sin(Eg)) w vec nohead lt 1 lc 7,\
   "" us (1.4*cos(Eg)):(1.4*sin(Eg)):1 w labels

#   for [n=2:6] (t0=E(g,n),t1=E(g,n-1),0.1+n*0.05)*cos((t1-t0)*t+t0),(0.1+n*0.05)*sin((t1-t0)*t+t0) lt 1 lc 2

set out
print "cordic_path.ps"
