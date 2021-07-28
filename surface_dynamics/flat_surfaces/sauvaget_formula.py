r""""
Sauvaget's formulas for the volume of the minimal strata.

The article of reference is [Sau18]_ see also [CMSZ20]_.
"""

from sage.rings.all import ZZ, QQ
from sage.misc.misc import newton_method_sizes

def minimal_strata_CMSZ(gmax):
    n = 2*gmax-1
    R = QQ['u']
    u = R.gen()
    # B(u) = (u/2) / sinh(u/2)
    # NOTE: coincide with formula (15) in CMSZ
    B = (2 * (u/2)._sinh_series(n+1).shift(-1)).inverse_series_trunc(n+1)
    print('B = {}'.format(B))
    Q = u * (sum(factorial(j-1) * B[j] * u**(j) for j in range(1,n)))._exp_series(n+1)
    # A = formula (14) in [CSMZ20] (up to a shift by t)
    tA = Q.revert_series(n+1).shift(-1).inverse_series_trunc(n)
    print('t*A = {}'.format(tA))

    # normalized values of volumes in CMSZ are
    # v(m_1, ..., m_n) = (2m_1+1) (2m_2+1) ... (2m_n+1) Vol(m_1, m_2, ..., m_n)
    return [2 * (2*pi)**(2*g) * (-1)**g / (2*g-2+1) / factorial(2*g-2+1) * tA[2*g] for g in range(1,gmax)]

def check():
    R = QQ['t']
    S = (2 * (t/2)._sin_series(11).shift(-1)).inverse_series_trunc(11)
    F = 1 + QQ((1,24))*t**2  + QQ((1,640))*3*t**4 + QQ((305,580608))*5*t**6
    for g in range(1,4):
        print(F._power_trunc(2*g,2*g+1)[2*g], (ZZ(2)*g).factorial() * S[2*g])

def compose_mod(f, g, prec):
    r"""
    return f(g(x)) mod x^prec

    EXAMPLES::

        sage: x = polygen(QQ, 'x')
        sage: f = x + x**2
        sage: g = x - 2*x**3 + 3*x**4

        sage: assert compose_mod(f, g, 3) == f(g).truncate(3)
        sage: assert compose_mod(f, g, 4) == f(g).truncate(4)
        sage: assert compose_mod(f, g, 5) == f(g).truncate(5)

        sage: f = x/2 + x**2 + x**5
        sage: g = x - 1/3*x**3 + 3*x**4 + x**9

        sage: assert compose_mod(f, g, 3) == f(g).truncate(3)
        sage: assert compose_mod(f, g, 4) == f(g).truncate(4)
        sage: assert compose_mod(f, g, 5) == f(g).truncate(5)
    """
    if not f[0].is_zero() or not g[0].is_zero():
        raise ValueError('invalid input')
    if f.parent() is not g.parent():
        raise ValueError('invalid input')

    R = f.parent()
    # Hörner scheme for g mod x^prec
    # f(x) * (a_1 + f(x) * (a_2 + f(x) * ... (a_n)))
    f = f.truncate(prec)
    d = f.degree()
    res = R(f[d])
    for i in range(d-1, -1, -1):
        res = res._mul_trunc_(g, prec)
        res += f[i]
    return res

def compositional_inverse(f, prec):
    r"""
    Return the compositional of ``f`` up to precision ``prec``.
    """
    if not f[0].is_zero():
        raise ValueError('no compositional inverse')
    R = f.parent()
    x = R.gen()
    g = R([0, f[1].inverse_of_unit()])  # approximation up to O(x^2)
    for prec in newton_method_sizes(prec)[1:]:
        g -= g.derivative() * (compose_mod(f, g, prec) - x)
    return g

def volumes(n):
    t = QQ['t'].gen()

    Shat = sum((ZZ(2)**(2*g-1) - 1)/ZZ(2)**(2*g-1) * (2*g+1) * t**(2*g+1) for g in range(n))
    print("Shat=", Shat)
    y = Shat.revert_series(n)
    print("y=", y)
    return (y//x).inverse_series_trunc(n)

