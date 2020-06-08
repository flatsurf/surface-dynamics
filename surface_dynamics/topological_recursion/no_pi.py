#*****************************************************************************
#       Copyright (C) 2020 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************
from sage.misc.cachefunc import cached_function
from sage.rings.all import ZZ, QQ
from sage.arith.misc import bernoulli, factorial
ZZ_0 = ZZ.zero()
ZZ_1 = ZZ.one()
ZZ_2 = ZZ(2)

@cached_function
def zeta_no_pi(k):
    r"""
    zeta(2m+2) = (-1)^m B_{2m+2} (2pi)^(2m+2) / (2 (2m+2)!)

    EXAMPLES::

        sage: from surface_dynamics.topological_recursion.no_pi import zeta_no_pi
        sage: zeta_no_pi(2)
        1/6
        sage: zeta_no_pi(4)
        1/90
        sage: zeta_no_pi(6)
        1/945

        sage: zeta(2)
        1/6*pi^2
        sage: zeta(4)
        1/90*pi^4
        sage: zeta(6)
        1/945*pi^6
    """
    k = ZZ(k)
    if k % 2 == 1:
        raise ValueError
    m = (k - 2) // 2
    return ZZ(-1)**m * bernoulli(k) * ZZ_2**(k-1) / k.factorial()

