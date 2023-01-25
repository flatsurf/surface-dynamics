r"""
Two variables generalized multiple zeta values

This module is about the sums

.. MATH::

    Z2(a, b, c) = sum_{x, y} 1 / (x^a y^b (x + y)^c)
"""
#*****************************************************************************
#       Copyright (C) 2019-2023 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************


from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.arith.misc import binomial

import cypari2.handle_error

try:
    from sage.modular.multiple_zeta import Multizetas
except (ImportError, cypari2.handle_error.PariError):
    def Multizetas(*args, **kwds):
        raise ValueError('your sage version does not support multiple zeta values')

from .options import VERBOSE, DIVERGENT_MZV

def is_Z2_convergent(a, b, c):
    r"""
    Return whether Z2(a, b, c) is convergent.

    TESTS::

        sage: from surface_dynamics.generalized_multiple_zeta_values import is_Z2_convergent

    Convergent examples::

        sage: assert is_Z2_convergent(1,1,1)
        sage: assert is_Z2_convergent(2,2,0)
        sage: assert is_Z2_convergent(0,0,3)

    Divergent examples::

        sage: assert not is_Z2_convergent(0,0,2)
        sage: assert not is_Z2_convergent(1,2,0)
    """
    from sage.geometry.polyhedron.constructor import Polyhedron
    x, y = ZZ['x,y'].gens()
    poly = x**a * y**b * (x+y)**c
    newton_polytope = Polyhedron(vertices=poly.exponents(), rays=[(-1,0),(0,-1)])
    V = newton_polytope.intersection(Polyhedron(rays=[(1,1)])).vertices()
    r = max(max(v.vector()) for v in V)
    return r > 1


def Z2(a, b, c, check_convergence=True):
    M = Multizetas(QQ)
    if a == 0 and b == 0:
        return M((c-1,)) - M((c,))
    else:
        return sum(binomial(a+i-1,i) * M((b-i,c+a+i)) for i in range(b)) + sum(binomial(b+i-1,i) * M((a-i,c+b+i)) for i in range(a))
