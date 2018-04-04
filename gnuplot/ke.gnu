# Algorithms to solve Kepler's equation CORDIC-like
# by M. Zechmeister

# Define the look-up tables
# Arrays are emulated with the word function (one based!)

# half angle table
a(n) = n ? a(n-1).sprintf("%.16g ",pi/2**n) : ""
ak = a(60)
a(k) = word(ak,k)

# cos table
cosa(n) = n ? cosa(n-1).sprintf("%.16g ",cos(pi/2**n)) : ""
ck = cosa(60)
cosa(k) = word(ck,k)

# sin table
sina(n) = n ? sina(n-1).sprintf("%.16g ",sin(pi/2**n)) : ""
sk = sina(60)
sina(k) = word(sk,k)

# arctangent table
atr(n) = n ? atr(n-1).sprintf("%.16g ",atan(2**-(n-1))) : ""
atrk = atr(60)
atr(k) = word(atrk,k)

# scale factor table
K(n) = (K=1, s="", sum[k=0:n](s=s.sprintf("%.16g ", K), K=K/sqrt(1+4**-k)), s)
Kk = K(60)
K(k) = word(Kk,k)

# twosided, half angle
# recursive implementation
# we need to pass all temporary variables
# _Epn(M,e,N,n,E,cosE,sinE) = n>N? E : (sgn=(E-e*sinE>M?-1:1), _Epn(M,e,N,n=n+1,E+sgn*a(n),cosE*cosa(n)-sgn*sinE*sina(n),sgn*cosE*sina(n)+sinE*cosa(n)))
# E(M) = (E0=floor(M/(2*pi)+0.5)*2*pi, _Epn(M,e,N,1,E0,1,0))
#
_Epn(M,e,N) = (E=floor(M/(2*pi)+0.5)*2*pi, cosE=1, sinE=0, sum [_n=1:N] (E=E+(sgn=(E-e*sinE>M?-1:1))*a(_n), cosE=(cn=cosE)*(ca=cosa(_n))-sinE*(sa=sgn*sina(_n)), sinE=cn*sa+sinE*ca), E)

# twosided, half angle (semi-complex version)
_E(M,e,N) = (E=floor(M/(2*pi)+0.5)*2*pi, z={1,0}, sum [_n=1:N] (sgn=(E-e*imag(z)>M?-1:1), z=z*cosa(_n)+ {0,1}*sgn*sina(_n), E=E+sgn*a(_n)), E)

# onesided, half angle (only positive shift)
_E(M,e,N) = (E = floor(M/(2*pi)+0.5)*2*pi,\
             cosE = 1,\
             sinE = 0,\
             sgn = M>E?1:-1,\
             E = E*sgn,\
             m = M*sgn,\
             sum [_n=1:N] (\
                st = cosE*sina(_n)+sinE*cosa(_n),\
                Et = E+a(_n),\
                (Et-e*st<m) ? (cosE=(cosE*cosa(_n)-sinE*sina(_n)),sinE=st,E=Et):0\
             ),\
             sinE = sgn*sinE,\
             E*sgn\
            )

# twosided, arctangent radix (atr version)
_Eatr(M,e,N) = (E=floor(M/(2*pi)+0.5)*2*pi, cosE=0, sinE=M>E?1:-1, E=E+pi/2*sinE, sum [_n=0:N-1] (sgn=(E-e*K(_n+1)*sinE>M)?-1:1, cosE=(ct=cosE)-sinE/(sgn<<_n), sinE=sinE+ct/(sgn<<_n), E=E+sgn*atr(_n)), cosE=K(n)*cosE, sinE=K(_n)*sinE, E)

# 
_Ec(M,e,N) = (E=floor(M/(2*pi)+0.5)*2*pi, cosE=1, scE=0, sgn=M>E?1:-1, E=E*sgn, m=M*sgn, sum [_n=1:N] (sct=cosE*sina(_n)+scE*cosa(_n)+E*(cosa(_n)-1)-a(_n),Et=E+a(_n), (Et*(1-e)-e*sct<m)?(cosE=(cosE*cosa(_n)-(E+scE)*sina(_n)),scE=sct,E=Et):0), sinE=sgn*sinE, E*sgn)

# CORDIC-like coupled with one Newton iteration
E29N(M) = (E=_E(M,e,29), E+(M-(E-e*sinE))/(1-e*cosE))

# CORCIC-like coupled with one Halley iteration
E19H(M) = (E=_E(M,e,19), dM=M-E+e*sinE, g=1-e*cosE, E=E+dM*g/(g**2+0.5*dM*e*sinE))

# Kepler's equation with sin
_M(E) = E - e*sin(E)

# accurate version with Taylor expansion including range reduction
#Macc(E) = E-e*E+ e*(E**3/3!-E**5/5!+E**7/7!-E**9/9!+E**11/11!-E**11/11!)
#Macc(E) = E-e*E+ e*sum[i=0:7](n=3+2*i, E**n/n!*(i&1?-1:1))
_Macc(E) = (x=E-(E0=floor(E/(2*pi)+0.5)*2*pi), x=(x>pi/2)?(E0=E0-pi+2*x, pi-x):(x<-pi/2)?(E0=E0+pi+2*x, -pi-x):x, x - e*x -e*(sum[_i=0:8](x=-x, _n=3+2*_i, x**_n/_n!)) +E0)


_EN(M,E) = (En=E+(M-(E-e*sin(E)))/(1-e*cos(E)), (abs(En-E)<1e-4)||((_n=_n+1)>20)?E:_EN(M,En))
EN(M) = (E=M+0.85*e*sgn(sin(M)), _n=0, _EN(M,E))

M(E) = _M(E)
E(M) = _Epn(M,e,N)
E(M) = _E(M,e,N)
E55(M) = _E(M,e,55)

# if N is not integer: "range specifiers of sum must have integer values"

N = 55
e = 1.

