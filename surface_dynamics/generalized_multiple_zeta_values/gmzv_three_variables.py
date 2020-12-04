r"""
Three variables generalized multiple zeta values

This module is about the sums

.. MATH::

    Z3(a, b, c, d, e, f, g)
    =
    \sum_{x,y,z \geq 1} \frac{1}{x^a y^b z^c (x+y)^d (x+z)^e (y+z)^f (x+y+z)^g}
"""
#*****************************************************************************
#       Copyright (C) 2019-2023 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.rational_field import QQ
from sage.arith.misc import binomial
from sage.modular.multiple_zeta import Multizetas

from .options import VERBOSE, DIVERGENT_MZV

def Z3_sort_abc(a,b,c,d,e,f,g):
    if b > a:
        if c > b:
            # c > b > a
            a,b,c,d,e,f,g = c,b,a,f,e,d,g
        elif c > a:
            # b >= c > a
            a,b,c,d,e,f,g = b,c,a,f,d,e,g
        else:
            # b > a >= c
            a,b,c,d,e,f,g = b,a,c,d,f,e,g
    elif c > a:
        # c > a >= b
        a,b,c,d,e,f,g = c,a,b,e,f,d,g
    elif c > b:
        # a >= c >= b
        a,b,c,d,e,f,g = a,c,b,e,d,f,g
    else:
        pass
    assert a >= b >= c
    return a,b,c,d,e,f,g


def is_Z3_convergent(a, b, c, d, e, f, g):
    r"""
    Return whether `Z3(a,b,c,d,e,f,g)` is convergent.

    TESTS::

        sage: from surface_dynamics.misc.generalized_multiple_zeta_values import is_Z3_convergent

    Convergent examples::

        sage: assert is_Z3_convergent(2,0,2,2,0,0,0)
        sage: assert is_Z3_convergent(0,0,0,0,0,0,4)
        sage: assert is_Z3_convergent(1,0,0,1,0,0,2)
        sage: assert is_Z3_convergent(0,1,0,1,0,0,2)
        sage: assert is_Z3_convergent(0,1,0,0,0,1,2)

    Divergent examples::

        sage: assert not is_Z3_convergent(0,0,0,1,1,1,0)
        sage: assert not is_Z3_convergent(0,0,0,0,0,0,3)

    """
    from sage.geometry.polyhedron.constructor import Polyhedron
    x, y, z = ZZ['x,y,z'].gens()
    poly = x**a * y**b * z**c * (x+y)**d * (x+z)**e * (y+z)**f * (x+y+z)**g
    newton_polytope = Polyhedron(vertices=poly.exponents(), rays=[(-1,0,0),(0,-1,0),(0,0,-1)])
    V = newton_polytope.intersection(Polyhedron(rays=[(1,1,1)])).vertices()
    r = max(max(v.vector()) for v in V)
    return r > 1


def Z3(a, b, c, d, e, f, g, check_convergence=True):
    r"""
    The function ``Z3(a,b,c,d,e,f,g)``.

    The reduction algorithm was designed by Bill Allombert.

    EXAMPLES::

        sage: from surface_dynamics.generalized_multiple_zeta_values import Z3

        sage: M = Multizetas(QQ) # optional: mzv

        sage: Z3(1,1,1,1,1,1,1) # optional: mzv
        21/2*ζ(1,1,5) + 9/2*ζ(1,2,4) - 3/2*ζ(1,3,3) - 3/2*ζ(1,4,2) + 9/2*ζ(1,6)
        sage: Z3(3,0,0,0,0,3,0) # optional: mzv
        6*ζ(1,4) - 12*ζ(1,5) + 3*ζ(2,3) - 6*ζ(2,4) + ζ(3,2) - 2*ζ(3,3)

        sage: assert Z3(2,3,4,0,0,0,0) == M((2,)) * M((3,)) * M((4,)) # optional: mzv
        sage: assert Z3(1,0,0,2,0,0,3) == M((1,2,3)) # optional: mzv

        sage: assert Z3(0,0,0,2,0,1,1) == M((4,)) / 2 # optional: mzv
        sage: assert Z3(1,0,1,1,0,0,1) == 3 * M((1,1,2)) # optional: mzv

        sage: assert Z3(0,0,0,0,0,0,4) == 1/2 * M((2,)) - 3/2 * M((3,)) + M((4,)) # optional: mzv

        sage: assert Z3(0,0,0,2,0,1,1) == 2 * M((1,1,2)) - M((2,2)) - 3 *M((1,3)) # optional: mzv
    """
    M = Multizetas(QQ)
    CHECK_CONVERGENCE = False

    # x^a y^b z^c (x+y)^d (x+z)^e (y+z)^f (x+y+z)^g
    if VERBOSE:
        print("Z3({},{},{},{},{},{},{})".format(a,b,c,d,e,f,g))
    if a < 0 or b < 0 or c < 0 or d < 0 or e < 0 or f < 0 or g < 0:
        raise ValueError("invalid exponents for Z3: a={} b={} c={} d={} e={} f={} g={}".format(a,b,c,d,e,f,g))

    if not DIVERGENT_MZV and CHECK_CONVERGENCE and not is_Z3_convergent(a,b,c,d,e,f,g):
        raise DivergentZetaError("divergent Z3({},{},{},{},{},{},{})".format(a,b,c,d,e,f,g))

    # step 1: try to get rid of the terms (x+y), (x+z), (y+z)
    if d and e and f:
        if VERBOSE:
            print("reduction (x+y+z) = ((x+y) + (x+z) + (y+z)) / 2")
        return (Z3(a,b,c,d-1,e,f,g+1,CHECK_CONVERGENCE) + Z3(a,b,c,d,e-1,f,g+1,CHECK_CONVERGENCE) + Z3(a,b,c,d,e,f-1,g+1,CHECK_CONVERGENCE)) / 2
    if a and f:
        # x^a (y+z)^f -> (x+y+z)^g
        # additive version Z3(a-1,b,c,d,e,f,g+1) + Z3(a,b,c,d,e,f-1,g+1)
        if VERBOSE:
            print("reduction (x+y+z) = (x) + (y+z)")
        return sum(binomial(a+k-1, k) * Z3(0,b,c,d,e,f-k,g+a+k,CHECK_CONVERGENCE) for k in range(f)) + \
               sum(binomial(f+k-1, k) * Z3(a-k,b,c,d,e,0,g+f+k,CHECK_CONVERGENCE) for k in range(a))
    if b and e:
        # y^b (x+z)^e -> (x+y+z)^g
        # additive version Z3(a,b-1,c,d,e,f,g+1) + Z3(a,b,c,d,e-1,f,g+1)
        if VERBOSE:
            print("reduction (x+y+z) = (y) + (x+z)")
        return sum(binomial(b+k-1, k) * Z3(a,0,c,d,e-k,f,g+b+k,CHECK_CONVERGENCE) for k in range(e)) + \
               sum(binomial(e+k-1, k) * Z3(a,b-k,c,d,0,f,g+e+k,CHECK_CONVERGENCE) for k in range(b))
    if c and d:
        # z^c (x+y)^d -> (x+y+z)^g
        # additive version Z3(a,b,c-1,d,e,f,g+1) + Z3(a,b,c,d-1,e,f,g+1)
        if VERBOSE:
            print("reduction (x+y+z) = (z) + (x+y)")
        return sum(binomial(c+k-1, k) * Z3(a,b,0,d-k,e,f,g+c+k,CHECK_CONVERGENCE) for k in range(d)) + \
               sum(binomial(d+k-1, k) * Z3(a,b,c-k,0,e,f,g+d+k,CHECK_CONVERGENCE) for k in range(c))

    assert d*e*f == a*f == b*e == c*d == 0

    a,b,c,d,e,f,g = Z3_sort_abc(a,b,c,d,e,f,g)

    # step 2: kill c
    if c:
        if VERBOSE:
            print("reduction (x+y+z) = (x) + (y) + (z)")
        return Z3(a-1,b,c,d,e,f,g+1,CHECK_CONVERGENCE) + Z3(a,b-1,c,d,e,f,g+1,CHECK_CONVERGENCE) + Z3(a,b,c-1,d,e,f,g+1,CHECK_CONVERGENCE)

    assert c == 0

    if b:
        # x^a y^b -> (x+y)^d
        # additive version Z3(a-1,b,0,d+1,e,f,g) + Z3(a,b-1,0,d+1,e,f,g)
        if VERBOSE:
            print("reduction (x+y) = (x) + (y)")
        return sum(binomial(a+k-1, k) * Z3(0,b-k,c,d+a+k,e,f,g,CHECK_CONVERGENCE) for k in range(b)) + \
               sum(binomial(b+k-1, k) * Z3(a-k,0,c,d+b+k,e,f,g,CHECK_CONVERGENCE) for k in range(a))

    assert b == c == 0

    assert b == c == d*e*f == 0

    if a and f:
        raise RuntimeError

    assert b == c == d*e*f == a*f == 0

    if a == 0:
        d,e,f = sorted([d,e,f], reverse=True)

    if a and d and e:
        if VERBOSE:
            print("reduction (x+y+z) = (x+y) + (x+z) - (x)")
        return Z3(a,b,c,d-1,e,f,g+1,CHECK_CONVERGENCE) + Z3(a,b,c,d,e-1,f,g+1,CHECK_CONVERGENCE) - Z3(a-1,b,c,d,e,f,g+1,CHECK_CONVERGENCE)

    assert b == c == d*e*f == a*f == a*d*e == 0

    # x^a (x+y)^d (x+y+z)^g
    if e == f == 0:
        from .generalized_multiple_zeta_values import convergent_multizeta
        return convergent_multizeta((a,d,g))
    # x^a (x+z)^e (x+y+z)^g
    if d == f == 0:
        from .generalized_multiple_zeta_values import convergent_multizeta
        return convergent_multizeta((a,e,g))
    if a == f == 0:
        # (x+y)^d (x+z)^e (x+y+z)^g
        if g == 0:
            # (x+y)^d (x+z)^e
            assert d and e
            return M((d-1,e)) - M((d,e)) + M((e-1,d)) - M((e,d)) + M((d+e-1,)) - M((d+e,))
        elif d == e == 1:
            if VERBOSE:
                print("Bill formulas (x+y) (x+z) (x+y+z)^n")
            return sum(M((1,g+1-i,i)) for i in range(2,g+1)) - 3*M((1,g+1))
        elif g == 1:
            if VERBOSE:
                print("[0,0,0,d,0,f,1] = [f,0,d-1,1,0,0,1] - [f,d,1] - zeta([d,1,f]) - zeta([1+d,f])")
            from .gmzv_two_variables import Z2
            f,d = sorted([d,e])
            return Z3(f,0,d-1,1,0,0,1,CHECK_CONVERGENCE) - Z2(f,d,1,CHECK_CONVERGENCE) - M((f,1,d)) - M((f,d+1))
        else:
            if VERBOSE:
                print("[0,0,0,0,e,f,g] = [e,0,f,g,0,0,0] - [e,0,f,0,0,0,g] - Z(e,g,f) - Z(e,f+g) - Z(e+f,g)")
            from .gmzv_two_variables import Z2
            f,d = sorted([d,e])
            return Z3(f,0,d,g,0,0,0,CHECK_CONVERGENCE) - Z3(f,d,0,0,0,0,g,CHECK_CONVERGENCE) - Z2(f,d,g,CHECK_CONVERGENCE) - M((f,g,d)) - M((f,g+d))

    raise RuntimeError
