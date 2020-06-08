r"""
Kontsevich topological recursion that computes intersection of psi classes.

EXAMPLES::

    sage: from surface_dynamics import KontsevichTR
    sage: K = KontsevichTR()

    sage: for p in [(4,),(5,0),(4,1),(3,2)]:
    ....:     n = len(p)
    ....:     assert sum(p) == n + 3
    ....:     c = K.F(2, n, p) / prod((2*i + 1).multifactorial(2) for i in p)
    ....:     print(p, c)
    (4,) 1/1152
    (5, 0) 1/1152
    (4, 1) 1/384
    (3, 2) 29/5760
"""
#*****************************************************************************
#       Copyright (C) 2020 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************
from sage.misc.misc_c import prod
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.all import ZZ, QQ
from .topological_recursion import TopologicalRecursion

ZZ_1 = ZZ.one()

def psi_correlator(*args):
    r"""
    Return the integral of psi classes

    EXAMPLES:

    Examples in genus 0::

        sage: from surface_dynamics.topological_recursion.kontsevich import psi_correlator
        sage: psi_correlator(0,0,0)
        1
        sage: psi_correlator(1,0,0,0)
        1
        sage: psi_correlator(2,0,0,0,0)
        1
        sage: psi_correlator(1,1,0,0,0)
        2

    Examples in genus 1::

        sage: psi_correlator(1)
        1/24
        sage: psi_correlator(2, 0)
        1/24
        sage: psi_correlator(1, 1)
        1/24
        sage: psi_correlator(3, 0, 0)
        1/24
        sage: psi_correlator(2, 1, 0)
        1/12
        sage: psi_correlator(1, 1, 1)
        1/12

    genus 2::

        sage: psi_correlator(7)
        1/82944
        sage: psi_correlator(7, 1)
        5/82944
        sage: psi_correlator(6, 2)
        77/414720
        sage: psi_correlator(5, 3)
        503/1451520
        sage: psi_correlator(4, 4)
        607/1451520
    """
    assert all(x in ZZ and x >= 0 for x in args)
    n = len(args)
    s = sum(args)
    if (s - n) % 3:
        raise ValueError("the composition should sum up to 3*g - 3 + n")
    g = (s -n)//3 + 1
    return KontsevichTR().F(g, n, args) / prod((2*i + 1).multifactorial(2) for i in args)

class KontsevichTR(UniqueRepresentation, TopologicalRecursion):
    r"""
    Topological recursion for intersection of psi classes (Witten's conjecture)

    `\int_{Mgn} exp(L_i^2/2 psi_i)`

    EXAMPLES::

        sage: from surface_dynamics import KontsevichTR
        sage: K = KontsevichTR()

        sage: [K.F(2, 1, (i,)) for i in range(5)]
        [0, 0, 0, 0, 105/128]

        sage: K.polynomial(0,3)
        1
        sage: K.polynomial(0,4)
        3*b0 + 3*b1 + 3*b2 + 3*b3
        sage: K.polynomial(0,5)
        15*b0^2 + 18*b0*b1 + ...

        sage: K.polynomial(1,1)
        1/8*b0
        sage: K.polynomial(1,2)
        5/8*b0^2 + 3/8*b0*b1 + 5/8*b1^2

        sage: K.write(0, 4)
        [1, 0, 0, 0] 1/2
        sage: K.write(0, 5)
        [2, 0, 0, 0, 0] 1/8
        [1, 1, 0, 0, 0] 1/2
        sage: K.write(1, 2)
        [2, 0] 1/192
        [1, 1] 1/96
        sage: K.write(2, 1)
        [4] 1/442368
    """
    def __init__(self):
        TopologicalRecursion.__init__(self)

    def A(self, i, j, k):
        r"""
        A-data for Kontsevich recursion
        """
        if i == 0 and j == 0 and k == 0:
            return ZZ_1
        else:
            raise ValueError

    def B(self, g, n, i, j):
        r"""
        B-data for Kontsevich topological recursion
        """
        # before we had s=I[0]+I[m]=i+j
        # i + j = k + 1
        k = i + j - 1
        if k >= 0:
            yield (k, ZZ(2*j + 1))

    def C(self, i, jmax, kmax, smax):
        r"""
        C-data for Konstevich topological recursion
        """
        if i < 2 or jmax > i-2 or kmax > i-2:
            # this is a waste of time...
            return

        # iterate through indices (j,k) with j + k = i-2
        jmin = max(0, i - 2 - kmax)
        jmax = min(i - 2, jmax)
        for j in range(jmin, jmax+1):
            yield (j, i-j-2, ZZ_1)

    def D(self, i):
        r"""
        D-data for Kontsevich topological recursion
        """
        return QQ(((i == 1), 8))

