#! /usr/bin/python3
from math import atan, tau
from math import atanh, frexp, log
ln2 = log(2)

kmax = 53     # largest shift value
R = 1 << 61   # binary point

ak = []       # atan lookup table
K = 1         # scale correction factor
for k in range(kmax+1):
    ak += [(atan(2**-k), k)]
    if 2*k <= kmax:
        K /= 1 + 4**-k
        ak += ak[-1:]

KR = K * R
i64_ak = [(round(a*R), k) for a,k in ak]

def i64_Ecs(M, e):
    """
    Solve Kepler's equation with CORDIC double rotation (fix point version).

    M = E - e*sin(E)

    Args:
        M (float): Mean anomaly [rad].
        e (float): Eccentricity.

    Returns:
        The triple (E, e*cos(E), e*sin(E)).

    Example:
        >>> i64_Ecs(2-sin(2), 1)
        (2.0, -0.41614683654714246, 0.9092974268256817)

    """
    t = M - tau*round(M/tau)
    t = round(t * R)
    x, y = round(KR*e), 0
    for a,k in i64_ak:
        s = t+y >> 63
        t -= s ^ s + a
        x, y = x - (s^s+(y>>k)),\
               y + (s^s+(x>>k))
    return M+y/R, x/R, y/R


from math import copysign, ldexp

def Ecs(M, e):
    """
    Solve Kepler's equation with CORDIC double rotation.

    M = E - e*sin(E)

    Args:
        M (float): Mean anomaly [rad].
        e (float): Eccentricity.

    Returns:
        The triple (E, e*cos(E), e*sin(E)).

    Example
    -------
    >>> Ecs(2-sin(2), 1)
    (2.0, -0.4161468365471421, 0.9092974268256817)

    """
    M = M - tau*round(M/tau)
    x, y = K*e, 0
    E = 0
    for a,k in ak:
        sgn = copysign(1, M-(E-y))
        E += sgn * a
        x, y = x - sgn*ldexp(y, -k),\
               y + sgn*ldexp(x, -k)
    return E, x, y


def Ecs_z(M, e):
    """
    Solve Kepler's equation with CORDIC double rotation.

    Uses cumulator z (saves a subtraction, lower precision at zero).

    M = E - e*sin(E)

    Args:
        M (float): Mean anomaly [rad].
        e (float): Eccentricity.

    Returns:
        The triple (E, e*cos(E), e*sin(E)).

    Example:
        >>> Ecs_z(2-sin(2), 1)
        (1.9999999999999998, -0.4161468365471422, 0.9092974268256815)

    """
    t = M - tau*round(M/tau)
    x, y = K*e, 0
    for a,k in ak:
        sgn = copysign(1, t+y)
        t -= sgn * a
        x, y = x - sgn*ldexp(y, -k),\
               y + sgn*ldexp(x, -k)
    return M+y, x, y


ahk = []      # atanh lookup table
Kh = 1        # scale correction factor
for k in range(1, kmax+1):
    ahk += [(atanh(2**-k), k)]
    if 2*k <= kmax:
        Kh /= 1 - 4**-k
        ahk += ahk[-1:]

KhR = Kh * R
i64_ahk = [(round(a*R), k) for a,k in ahk]

def i64_Hcs(M, e):
    """
    Solve hyperbolic Kepler's equation with CORDIC double rotation (fix point version).

    M = e*sinh(H) - H

    Args:
        M (float): Mean anomaly.
        e (float): Eccentricity.

    Returns:
        The triple (E, e*cos(E), e*sin(E)).

    Examples:
        >>> i64_Hcs(sinh(2)-2, 1)
        (2.0, 3.762195691083632, 3.626860407847019)

    """
    m = max(0, frexp(M/e)[1])
    eK = round(e * KhR) >> 1
    x = (eK<<m) + (eK>>m)
    y = (eK<<m) - (eK>>m)
    if M < 0: m, y = -m, -y
    t = M + m*ln2
    t = -round(t * R)
    for a,k in i64_ahk:
        s = -t-y >> 63
        t -= s ^ s + a
        x, y = x + (s^s+(y>>k)),\
               y + (s^s+(x>>k))
    return y/R-M, x/R, y/R

#gplot(M, H, [Hcs(Mi, e)[0] for Mi in M], [i64_Hcs(Mi, e)[0] for Mi in M], [Hcs_z(Mi, e)[0] for Mi in M], 'w l, "" us 1:3 w l t "Hcs", "" us 1:4 w l t "Hcs_z", "" us 1:($3-$2) axis x1y2, "" us 1:($4-$2) axis x1y2')

def Hcs(M, e):
    """
    Solve hyperbolic Kepler's equation with CORDIC double rotation.

    M = e*sinh(H) - H

    Args:
        M (float): Mean anomaly.
        e (float): Eccentricity.

    Returns:
        The triple (H, e*cos(H), e*sin(H)).

    Examples:
        >>> Hcs(sinh(2)-2, 1)
        (2.0, 3.762195691083633, 3.626860407847019)

    """
    m = max(0, frexp(M/e)[1])
    if M < 0: m = -m
    H = m * ln2
    eK = e * Kh
    x = ldexp(eK, m-1) + ldexp(eK, -m-1)
    y = ldexp(eK, m-1) - ldexp(eK, -m-1)
    for a,k in ahk:
        sgn = copysign(1, M-(y-H))
        H += sgn * a
        x, y = x + sgn*ldexp(y, -k),\
               y + sgn*ldexp(x, -k)
    return H, x, y


def Hcs_z(M, e):
    """
    Solve hyperbolic Kepler's equation with CORDIC double rotation.

    M = e*sinh(H) - H

    Args:
        M (float): Mean anomaly.
        e (float): Eccentricity.

    Returns:
        The triple (H, e*cos(H), e*sin(H)).

    Examples:
        >>> Hcs_z(sinh(2)-2, 1)
        (1.9999999999999973, 3.7621956910836305, 3.6268604078470164)

    """
    m = max(0, frexp(M/e)[1])
    if M < 0: m = -m
    eK = e * Kh
    x = ldexp(eK, m-1) + ldexp(eK, -m-1)
    y = ldexp(eK, m-1) - ldexp(eK, -m-1)
    t = M + m*ln2
    for a,k in ahk:
        sgn = copysign(1, t-y)
        t += sgn * a
        x, y = x + sgn*ldexp(y, -k),\
               y + sgn*ldexp(x, -k)
    return y-M, x, y



if __name__ == "__main__":
    '''
    Example
    -------
    ./ke_cordic_dbl.py

    '''
    import doctest
    from math import sin, sinh
    doctest.testmod()







