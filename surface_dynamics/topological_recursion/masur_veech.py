r"""
Topological recursion for Masur-Veech volumes.

This is the topological recursion as developped in Andersen et al.

EXAMPLES:

Masur-Veech volumes (without the pi power)::

    sage: from surface_dynamics import MasurVeechTR
    sage: MV = MasurVeechTR()
    sage: for g,n in [(0,4),(0,5),(1,1),(1,2),(1,3),(2,1),(2,2)]:
    ....:     coeff = 2**(4*g-2+n) * (4*g-4+n).factorial() / (6*g-7+2*n).factorial()
    ....:     v = coeff * MV.F(g, n, (0,)*n)
    ....:     print(g, n, v)
    0 4 2
    0 5 1
    1 1 2/3
    1 2 1/3
    1 3 11/60
    2 1 29/840
    2 2 337/18144

REFERENCES:

- J. E. Andersen, G. Borot, S. Charbonnier, V. Delecroix, A. Giacchetto,
  D. Lewanski and C. Wheeler
  "Topological recursion for Masur-Veech volumes"
  arXiv:1905.10352


.. TODO::

    - double check: write the formula that gives MasurVeech TR as a sum over stable
      graphs of the product of weighted Kontsevich
      (needs the list of stable graphs)
"""
#*****************************************************************************
#       Copyright (C) 2020 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************
from sage.misc.cachefunc import cached_method
from sage.rings.all import ZZ, QQ
from sage.structure.element import parent, get_coercion_model

from .topological_recursion import TopologicalRecursion
from .no_pi import zeta_no_pi

ZZ_1 = ZZ.one()
cm = get_coercion_model()

class MasurVeechTR(TopologicalRecursion):
    r"""
    Topological recursion for Masur-Veech volumes (Anderssen et al.)

    EXAMPLES:

    Below is the distribution of cylinders in genus zero and `n` equal 4, 5, 6, 7::

        sage: from surface_dynamics.topological_recursion import MasurVeechTR
        sage: MV = MasurVeechTR(polygen(QQ, 't'))
        sage: for n in range(4, 8):
        ....:     p = MV.F(0, n, (0,)*n)
        ....:     print(p / p(1))
        t
        5/9*t^2 + 4/9*t
        7/27*t^3 + 4/9*t^2 + 8/27*t
        1/9*t^4 + 8/27*t^3 + 256/675*t^2 + 16/75*t
    """
    def __init__(self, edge_weight=ZZ_1, vertex_weight=ZZ_1, cache_all=True):
        r"""
        INPUT:

        - ``edge_weight`` - a edge weight to be put on each cylinder
        """
        if edge_weight != ZZ_1 or vertex_weight != ZZ_1:
            P = cm.common_parent(vertex_weight / ZZ_1, edge_weight / ZZ_1)
            TopologicalRecursion.__init__(self, cache_all, base_ring=P)
        else:
            TopologicalRecursion.__init__(self, cache_all)
        self._edge_weight = edge_weight
        self._vertex_weight = vertex_weight

        if not self._vertex_weight.is_one():
            raise NotImplementedError

    def A(self, i, j, k):
        if i == 0 and j == 0 and k == 0:
            return ZZ_1

    def B(self, g, n, i, j):
        if i == 0 and j == 0:
            # twist
            kbound = 3*g - 3 + n
            for k in range(kbound + 1):
                yield (k, self._edge_weight * zeta_no_pi(2*k + 2))
        else:
            # Kontsevich initial data
            k = i + j - 1
            if k >= 0:
                yield (k, ZZ(2*j + 1))

    def C(self, i, jmax, kmax, smax):
        r"""
        Iterate through the non-zero ``(j, k, C(i, j, k))`` given the ``i``.

        INPUT:

        - ``jmax`` - max value for ``j``

        - ``kmax`` - max value for ``k``

        - ``smax`` - max value for ``j+k``

        TESTS::

            sage: from surface_dynamics.topological_recursion import MasurVeechTR
            sage: MV = MasurVeechTR()
            sage: for i in range(5):
            ....:     for s in range(10):
            ....:         for j,k,_ in MV.C(0, 5, 5, s):
            ....:             assert j+k <= s
        """
        cew = self._edge_weight

        if i == 0:
            for j in range(jmax+1):
                for k in range(min(smax-j, kmax) + 1):
                    yield (j, k, cew**2 * zeta_no_pi(2*j+2) * zeta_no_pi(2*k+2))

        if i >= 2 and i-2 <= smax:
            # j + k = i - 2
            jjmin = max(0, i - 2 - kmax)
            jjmax = min(i - 2, jmax)
            for j in range(jjmin, jjmax + 1):
                yield (j, i - j - 2, ZZ_1)

        for j in range(max(i-1, 0), jmax+1):
            # j + 1 >= i
            a = j + 1 - i
            for k in range(min(smax-j, kmax) + 1):
                yield (j, k, ZZ(2*k + 2*a + 1).factorial() * cew * zeta_no_pi(2*k + 2*a + 2) / ZZ(2*k + 1).factorial() / ZZ(2*a).factorial())

        for k in range(max(i-1, 0), kmax+1):
            # k + 1 >= i
            a = k + 1 - i
            for j in range(min(smax-k, jmax) + 1):
                yield (j, k, ZZ(2*j + 2*a + 1).factorial() * cew * zeta_no_pi(2*j + 2*a + 2) / ZZ(2*j + 1).factorial() / ZZ(2*a).factorial())

    def D(self, i):
        if i == 0:
            return self._edge_weight * zeta_no_pi(2) / 2
        elif i == 1:
            return QQ((1,8))
        else:
            return ZZ_0
