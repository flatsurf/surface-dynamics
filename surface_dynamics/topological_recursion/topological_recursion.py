r"""
Topological recursion

This module implements topological recursion. In particular it makes
available the computation of various integral of cohomology classes
on M_{g,n} (psi, lambda and kappa classes).

REFS:

- G. Borot
  "Topological recursion and geometry"
  arXiv:1705.09986

- J. E. Andersen, G. Borot, L. O. Chekhov, N. Orantin
  "The ABCD of topological recursion"
  arXiv:1703.03307

- R. Pandharipande
  "Cohomological field theory calculations"
  arXiv:1712.02528

TODO:

- other recursions:
  - ELSV, Eynard Bouchard-Marino conjectures and simple Hurwitz numbers (lambda classes)
  - bundle of quadratic differentials and co
- generating series
- simplification of the recursion via the generalized string and dilaton
"""
#*****************************************************************************
#       Copyright (C) 2020 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function
from six.moves import range

from collections import defaultdict
from itertools import combinations

from sage.misc.misc_c import prod
from sage.misc.cachefunc import cached_method

from sage.rings.all import ZZ, QQ
from sage.arith.misc import factorial, binomial
from sage.combinat.integer_lists.invlex import IntegerListsLex
from sage.combinat.composition import Compositions
from sage.combinat.partition import Partitions

ZZ_0 = ZZ.zero()
ZZ_2 = ZZ(2)

def compositions(g, n):
    r"""
    Return compositions of sum at most 3g-3+n and length n
    """
    for k in range(0, 3*g-3+n+1):
        for c in Compositions(k+n, length=n, min_part=1):
            yield [i-1 for i in c]

class TopologicalRecursion:
    r"""
    A generic topological recursion.

    For concrete examples see:

    - :class:`~surface_dynamics.topological_recursion.kontsevich.KontsevichTR`
    - :class:`~surface_dynamics.topological_recursion.masur_veech.MasurVeechTR`
    - :class:`~surface_dynamics.topological_recursion.symbolic.SymbolicTR`
    """
    def __init__(self, cache_all=False, base_ring=None):
        self._cache = {}
        self._cache_all = cache_all
        self._base_ring = QQ if base_ring is None else base_ring
        self._one = None

    # The functions A, B, C, D of the topological recursion
    # They have to be implemented in subclasses
    def A(self, i, j, k): raise NotImplementedError

    def B(self, g, n, i, j): raise NotImplementedError
    # iterator (k, B(i, j, k))

    def C(self, g, n, i): raise NotImplementedError
    # iterator (j, k, C(i, j, k))

    def D(self, i): raise NotImplementedError

    def F(self, g, n, I):
        r"""
        Return the value of the partition function F_{g,n}(I) where I = (i_0, i_1, ..., i_{n-1}) is
        a tuple of n non-negative integers whose sum is <= 3g-3+n.

        We assume here that the topological recursion given by A, B, C and D produces a partition
        function symmetric in the variables i_0, i_1, ..., i_{n-1}.
        """
        # we restrict to TR with sum(I) <= 3g-3+n
        assert len(I) == n
        if len(I) != n or sum(I) > 3*g-3+n:
            raise ValueError("length must be n = %d and sum at most 3g-3+n = %d, got I = %s" % (n, 3*g-3+n, I))

        I = tuple(sorted(I))

        if g == 0 and n == 3:
            # M_{0,3} initial condition
            if self._one is None:
                self._one = self._base_ring(self.A(I[0], I[1], I[2]))
            return self._one
        elif g == 1 and n == 1:
            # M_{1,1} initial condition
            return self.D(I[0])

        return self._virasoro(g, n, I)

    def _virasoro(self, g, n, I):
        # TODO: we should use exponential notations for I
        if (g, n, I) in self._cache:
            return self._cache[(g, n, I)]

        Idict = defaultdict(int)
        for i in I[1:]:
            Idict[i] += 1
        Ituple = sorted(Idict.items())

        # set to True for debugging information
        verbose = False

        if verbose:
            print("Frec({}, {}, {}".format(g, n, I))

        B = self.B
        C = self.C
        F = self.F

        if verbose:
            print("compute S1")
        J = I[1:]
        S1 = ZZ_0
        for im,_ in Ituple:
            fac = Idict[im]
            J = []
            for j,mult in Ituple:
                if j == im:
                    J.extend([j] * (mult - 1))
                else:
                    J.extend([j] * mult)
            J = tuple(J)

            s = sum(J) # = i + j
            for (a, value) in B(g, n-s-1, I[0], im):
                if verbose:
                    print("[S1] B({}, {}, {}) F({}, {}, {})".format(I[0], im, a, g, n-1, (a,) + J))
                S1 += fac * value * F(g, n-1, (a,) + J)
        if verbose:
            print("S1 = {}".format(S1))

        if verbose:
            print("S2")
        S2 = ZZ_0
        if g >= 1:
            # TODO: here the bound is actually for the sum since the new tuple
            # is (a, b) + I[1:].
            I1 = I[1:]
            abbound = 3 * (g-1) - 3 + (n+1) - sum(I1)
            for (a, b, value) in C(I[0], abbound, abbound, abbound):
                if verbose:
                    print("[S2] C({}, {}, {}) F({}, {}, {})".format(I[0], a, b, g-1, n+1, (a,b) + I1))
                S2 += value * F(g-1, n+1, (a, b) + I1)
        if verbose:
            print("S2 = {}".format(S2))

        if verbose:
            print("S3")
        S3 = ZZ_0
        for n1 in range(n):
            n2 = n - n1 - 1

            g1min = int(n1 < 2)
            g2min = int(n2 < 2)

            for M1 in IntegerListsLex(min_sum=n1, max_sum=n1, length=len(Ituple), ceiling=[v for k,v in Ituple]):
                I1 = []
                I2 = []
                fac = 1
                for m,(k,v) in zip(M1, Ituple):
                    I1.extend([k]*m)
                    I2.extend([k]*(v-m))
                    fac *= binomial(v, m)
                I1 = tuple(I1)
                I2 = tuple(I2)

                s1 = sum(I1)
                s2 = sum(I2)

                # the above choices of combinations already forces the genus
                # ie, 3g1-3+n1 >= s1 and 3g2-3+n2 >= s2
                # g1 >= (s1 + 3 - n1)/3
                gg1min = max(g1min, (s1-n1+4)//3)
                gg2min = max(g2min, (s2-n2+4)//3)
                gg1max = g - gg2min
                for g1 in range(gg1min, gg1max+1):
                    g2 = g - g1
                    a1bound = 3*g1 - 3 + (n1+1) - s1
                    a2bound = 3*g2 - 3 + (n2+1) - s2
                    for (a1, a2, value) in C(I[0], a1bound, a2bound, a1bound + a2bound):
                        f1 = F(g1, n1 + 1, (a1,) + I1)
                        f2 = F(g2, n2 + 1, (a2,) + I2)
                        if verbose:
                            print("[S3] C({}, {}, {}) F({}, {}, {}) F({}, {}, {})".format(
                            I[0], a1, a2, g1, n1+1, (a1,) + I1, g2, n2+1, (a2,) + I2))
                        S3 += fac * value * f1 * f2

        if verbose:
            print("S3 = {}".format(S3))

        value = S1 + (S2 + S3) / ZZ_2
        if self._cache_all or I[0] >= 2:
            # We do not cache string/dilaton if _cache_all = False
            self._cache[(g, n, I)] = value
        return value

    def polynomial(self, g, n, ring=None):
        if ring is None:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            if n == 1:
                # Otherwise, the b is not numbered!
                ring = PolynomialRing(self._base_ring, 'b0', n)
            else:
                ring = PolynomialRing(self._base_ring, 'b', n)
            b = ring.gens()
        p = ring.zero()
        p += sum(self.F(g,n,I) * prod(b[i]**j for i,j in enumerate(I)) for I in compositions(g,n))
        return p

    def write(self, g, n, s=None):
        r"""
        Print normalized values of the polynomials F_{g,n}
        """
        g = ZZ(g)
        n = ZZ(n)
        if n < 0:
            raise ValueError
        if s is None:
            s = range(3*g-3+n+1)
        elif isinstance(s, numbers.Integral):
            s = [s]

        for t in s:
            for p in Partitions(t+n, length=n):
                p = [i-1 for i in p]
                c = prod(factorial(2*i+1) for i in p)
                f = self.F(g, n, p)
                if f:
                    print("{} {}".format([i for i in p], self.F(g, n, p) / c))

