# ****************************************************************************
#       Copyright (C) 2011-2019 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from .origami_dense import PillowcaseCover_dense_pyx

try:
    # Trac #28652: Rework the constructor of PermutationGroupElement
    from sage.groups.perm_gps.constructor import PermutationGroupElement
except ImportError:
    from sage.groups.perm_gps.permgroup import PermutationGroupElement


def PillowcaseCover(g0, g1, g2, g3=None,
        sparse=False,
        check=True,
        as_tuple=False,
        positions=None, name=None):
    r"""
    Pillowcase cover constructor.

    The chosen flat structure is as follows::

        3-----2-----3
        |     .     |
        |     .     |
        |     .     |
        0-----1-----0
    """
    if not as_tuple:
        g0 = PermutationGroupElement(g0, check=check)
        g1 = PermutationGroupElement(g1, check=check)
        g2 = PermutationGroupElement(g2, check=check)
        if g3 is None:
            g3 = (~g2) * (~g1) * (~g0)
        else:
            g3 = PermutationGroupElement(g3, check=check)

        g0 = [i-1 for i in g0.domain()]
        g1 = [i-1 for i in g1.domain()]
        g2 = [i-1 for i in g2.domain()]
        g3 = [i-1 for i in g3.domain()]

        N = max([len(g0),len(g1),len(g2),len(g3)])
        g0.extend(range(len(g0),N))
        g1.extend(range(len(g1),N))
        g2.extend(range(len(g2),N))
        g3.extend(range(len(g3),N))

    elif check:
        s0 = set(g0)
        s1 = set(g1)
        s2 = set(g2)
        s3 = set(g3)
        N = len(g0)
        if len(g1) != N or len(g2) != N or len(g3) != N:
            raise ValueError("the four tuples must be of the same length")
        for i in range(N):
            if i not in g0:
                raise ValueError("%d is not in g0=%s" % (i, str(g0)))
            if i not in g1:
                raise ValueError("%d is not in g1=%s" % (i, str(g1)))
            if i not in g2:
                raise ValueError("%d is not in g2=%s" % (i, str(g2)))
            if i not in g3:
                raise ValueError("%d is not in g3=%s" % (i, str(g3)))

    pcc = PillowcaseCover_dense(tuple(g0),tuple(g1),tuple(g2),tuple(g3))

    if name is not None:
        pcc.rename(name)
#    if positions is not None:
#        o.set_positions(positions)
    if check:
        pcc._check()
    return pcc

class PillowcaseCover_dense(PillowcaseCover_dense_pyx):
    r"""
    Generic class for pillowcase cover.
    """
    def __repr__(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: PillowcaseCover('(1,2)(3,4)', '(1,3)', '()')   # indirect doctest
            g0 = (1,2)(3,4)
            g1 = (1,3)(2)(4)
            g2 = (1)(2)(3)(4)
            g3 = (1,4,3,2)
        """
        return '\n'.join("g%d = %s"%(i,self.g(i).cycle_string(True)) for i in range(4))

    def g(self, i=None):
        r"""
        Return the ``i``-th permutation that defines this pillowcase cover.
        """
        if i is None:
            return self.g(0), self.g(1), self.g(2), self.g(3)

        i = int(i)
        if i < 0 or i > 3:
            raise IndexError("the index i (={}) must be in {{0,1,2,3}}".format(i))

        return PermutationGroupElement([j+1 for j in self.g_tuple(i)], check=False)

    def _check(self):
        x = self.g(0) * self.g(1) * self.g(2) * self.g(3)
        if not x.is_one():
            raise ValueError

    def monodromy(self):
        r"""
        Return the monodromy group of the pillowcase cover.

        The monodromy group of an origami is the group generated by the
        permutations `g_i` for `i` in 0,1,2,3.
        """
        from sage.groups.perm_gps.permgroup import PermutationGroup
        return PermutationGroup(self.g())

    def as_graph(self):
        r"""
        Return the graph associated to self
        """
        from sage.graphs.digraph import DiGraph

        G = DiGraph(multiedges=True,loops=True)
        d = self.degree()
        g = [self.g_tuple(i) for i in range(4)]
        for i in range(d):
            for j in range(4):
                G.add_edge(i,g[j][i],j)
        return G

    def is_connected(self):
        r"""
        Check whether the origami is connected or not

        It is equivalent to ask whether the group generated by `r` and `u` acts
        transitively on the `\{1,\dots,n\}`.
        """
        return self.as_graph().is_connected()

    def connected_components(self):
        r"""
        Return the list of connected origami that composes this origami.
        """
        cc = self.as_graph().connected_components()
        g = [self.g_tuple(i) for i in range(4)]
        if len(cc) == 1:
            return [self]
        l = []
        for c in cc:
            gg = [[None] * len(c) for _ in range(4)]
            d = dict((c[i],i) for i in range(len(c)))
            for i in c:
                for j in range(4):
                    gg[j][d[i]] = d[g[j][i]]
            l.append(Pillowcase_cover(g[0],g[1],g[2],g[3],check=True,as_tuple=True))
        return l

    def is_orientable(self):
        r"""
        Test whether the foliation is orientable.
        """
        return self.as_graph().to_undirected().is_bipartite()

    def profile(self,i=None):
        r"""
        Return the profile (= ramification type above each pole).
        """
        if i is None:
            return [self.profile(i) for i in range(4)]
        return sorted((len(c) for c in self.g(i).cycle_tuples(singletons=True)),reverse=True)

    def stratum(self,fake_zeros=False):
        r"""
        Return the stratum of self. It may be either a stratum of Abelian or
        quadratic differentials.

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: PillowcaseCover('(1,2)(3,4)', '(1,3)', '()').stratum()
            Q_0(2, -1^6)
        """
        p = sum(self.profile(),[])
        if self.is_orientable():
            from surface_dynamics.flat_surfaces.abelian_strata import Stratum
            if fake_zeros:
                zeros = [(i-2)//2 for i in p]
            else:
                zeros = [(i-2)//2 for i in p if i != 2]
            if not zeros:
                return Stratum([0], k=1)
            return Stratum(zeros, k=1)

        else:
            from surface_dynamics.flat_surfaces.quadratic_strata import Stratum
            if fake_zeros:
                zeros = [i-2 for i in p]
            else:
                zeros = [i-2 for i in p if i != 2]
            return Stratum(zeros, k=2)

    def is_primitive(self, return_base=False):
        r"""
        A pillowcase cover is primitive if it does not cover an other pillowcase
        cover.
        """
        from sage.arith.all import is_prime
        if is_prime(self.degree()):
            return True

        return bool(gap.IsPrimitive(self.monodromy()))

    def orientation_cover(self):
        r"""
        Return the orientation cover as an origami.

        EXAMPLES:

        The pillowcase itself has cover a torus (made from 4 squares)::

            sage: from surface_dynamics import *

            sage: p0 = p1 = p2 = p3 = [0]
            sage: pc = PillowcaseCover(p0, p1, p2, p3, as_tuple=True)
            sage: pc.stratum()
            Q_0(-1^4)
            sage: o = pc.orientation_cover()
            sage: o
            (1,2)(3,4)
            (1,4)(2,3)
            sage: o.stratum()
            H_1(0)

        An example in Q(1,-1^5) whose cover belongs to H(2)::

            sage: p0 = [2,1,0]
            sage: p1 = [2,0,1]
            sage: p2 = [1,0,2]
            sage: p3 = [0,1,2]
            sage: pc = PillowcaseCover(p0, p1, p2, p3,as_tuple=True)
            sage: pc.stratum()
            Q_0(1, -1^5)
            sage: o = pc.orientation_cover()
            sage: o
            (1,2,3,4)(5,6)(7,10,9,8)(11,12)
            (1,10,5,12)(2,9)(3,8)(4,7,6,11)
            sage: o.stratum()
            H_2(2)

        A last example in Q(2^2)::

            sage: q = QuadraticCylinderDiagram('(0,1)-(2,3) (0,3)-(1,2)')
            sage: pc = q.cylcoord_to_pillowcase_cover([1,1,1,1], [2,2], [0,1])
            sage: pc.orientation_cover().stratum()
            H_3(1^4)
        """
        from surface_dynamics.misc.permutation import perm_invert
        from surface_dynamics.flat_surfaces.origamis.origami import Origami

        g0 = self.g_tuple(0)
        g1 = self.g_tuple(1)
        g2 = self.g_tuple(2)
        g3 = self.g_tuple(3)

        n = len(g0)

        r = [None] * (4*n)
        u = [None] * (4*n)

        h1 = perm_invert(g1)
        h2 = perm_invert(g2)

        for i in range(n):
            r[2*i] = 2*i+1
            r[2*i+1] = 2*g3[g2[i]]
            r[2*n+2*i+1] = 2*n + 2*i
            r[2*n+2*i] = 2*n + 2*g1[g0[i]] + 1

            u[2*i] = 2*n + 2*h2[i] + 1
            u[2*i+1] = 2*n + 2*g2[i]
            u[2*n+2*i] = 2*g1[i] + 1
            u[2*n+2*i+1] = 2*h1[i]

        return Origami(r, u, as_tuple=True)
