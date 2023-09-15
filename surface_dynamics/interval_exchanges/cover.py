# coding: utf8
r"""
Permutation cover

This module deals with combinatorial data for covering of connected
components of strata of Abelian and quadratic differentials. The main
feature is to be able to compute Lyapunov exponents.

.. TODO::

    It should be possible to compute the KZ action directly on isotypical
    components. That would dramatically reduce the dimension of the space!
"""
# *************************************************************************
# Copyright (C) 2015-20l6 Charles Fougeron <charlesfougeron@gmail.com>
#               2015-2021 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *************************************************************************
from __future__ import print_function, absolute_import

import numpy as np

from six.moves import range, map, filter, zip

from sage.misc.misc_c import prod

from sage.categories.additive_groups import AdditiveGroups
from sage.categories.groups import Groups
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.groups.libgap_group import GroupLibGAP

from sage.misc.cachefunc import cached_method
from sage.libs.gap.libgap import libgap
from sage.libs.gap.element import GapElement
from sage.rings.integer import Integer

from surface_dynamics.misc.permutation import perm_invert
from surface_dynamics.interval_exchanges.template import PermutationIET

_MulGroups = Groups()
_AddGroups = AdditiveGroups()

try:
    libgap_fail = libgap.fail
except AttributeError:
    # broken in old SageMath
    libgap_fail = libgap.eval("fail")

def to_gap(g):
    try:
        return libgap(g)
    except ValueError:
        # NOTE: the matrix interface is weird in SageMath < 9.0
        import sage.libs.gap.element
        g = g._gap_()
        if not isinstance(g, sage.libs.gap.element.GapElement):
            raise
        return g

class PermutationCover(object):
    r"""
    An interval exchange permutation together with covering data.

    Let `\pi` be the combinatorial data of an interval exchange transformation
    (or linear involution) on the alphabet `{1, 2, \ldots, m\}`. A cover of
    degree `d` is given by a list of permutations `\sigma_i \in S_d` for each `i
    \in \{1, 2, \ldots, m\}`.

    In order to do so, each interval on the base surface should come with an
    orientation. This orientation is automatically chosen by convention
    according to a clockwise orientation of the surface. The two copies of any
    interval have to be oriented alternatively in this chosen orientation and
    in the opposite to it.

    This convention is made such that the permutation associated to each interval
    is the action of the path going into the interval oriented according to
    the clockwise orientation and going out of the other one.

    This class store three attributes

    - ``_base`` -- combinatorial data of an i.e.t. or a l.i.

    - ``_degree_cover`` -- (integer) degree of the cover

    - ``_permut_cover`` -- list of permutations describing the cover

    - ``_inv_permut_cover`` -- list of inverses of ``_permut_cover``
    """
    def __init__(self, base, degree, perms):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: from surface_dynamics.interval_exchanges.cover import PermutationCover
            sage: p1 = iet.Permutation('a b c', 'c b a')
            sage: PermutationCover(p1, 2, [[0,1],[1,0],[1,0]])
            Covering of degree 2 of the permutation:
            a b c
            c b a

            sage: PermutationCover(p1, 2, [[0,1],[0,1],[0,1]])
            Traceback (most recent call last):
            ...
            ValueError: the cover is not connected

            sage: p1 = iet.Permutation('a b c d', 'b a d c')
            sage: PermutationCover(p1, 2, [[0,1],[1,0],[1,0]])
            Traceback (most recent call last):
            ...
            ValueError: the base must be irreducible

            sage: p2 = iet.Permutation('a b', 'b a')
            sage: PermutationCover(p2, 0, [[], []])
            Traceback (most recent call last):
            ...
            ValueError: the degree of the cover must be positive
        """
        if not base.is_irreducible():
            raise ValueError("the base must be irreducible")
        if degree == 0:
            raise ValueError("the degree of the cover must be positive")

        from surface_dynamics.misc.permutation import perms_are_transitive
        if not perms_are_transitive(perms):
            raise ValueError("the cover is not connected")

        self._base = base.__copy__()
        self._degree_cover = degree
        self._permut_cover = perms[:]
        self._inv_permut_cover = [perm_invert(p) for p in perms]

    def degree(self):
        return self._degree_cover

    def __repr__(self):
        r"""
        A representation of the generalized permutation cover.

        INPUT:

        - ``sep`` - (default: '\n') a separator for the two intervals

        OUTPUT:

        string -- the string that represents the permutation


        EXAMPLES::

            sage: from surface_dynamics import *
            sage: p1 = iet.Permutation('a b c', 'c b a')
            sage: p1.cover(['(1,2)', '(1,3)', '(2,3)'])
            Covering of degree 3 of the permutation:
            a b c
            c b a
        """
        s = 'Covering of degree %i of the permutation:\n' % (self.degree())
        s += str(self._base)
        return s

    def __eq__(self, other):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: p1 = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: p2 = iet.GeneralizedPermutation('a a b',' b c c')
            sage: p3 = iet.GeneralizedPermutation('a a b b', 'c c')
            sage: p1.cover(['(1)', '(1)', '(1)']) == p2.cover(['(1)', '(1)', '(1)'])
            True
            sage: p1.cover(['(1,2)', '', '']) == p1.cover(['(1,2)', '(1,2)', ''])
            False
            sage: p1.cover(['(1)', '(1)', '(1)']) == p3.cover(['(1)', '(1)', '(1)'])
            False

            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: G = Zmod(5)
            sage: p.regular_cover(G, [0, 1, 2]) == p.regular_cover(G, [0, 1, 2])
            True
            sage: p.regular_cover(G, [0, 1, 2]) == p.regular_cover(G, [1, 2, 0])
            False
        """
        return type(self) == type(other) and \
               self._base == other._base and \
               self.degree() == other.degree() and \
               self._permut_cover == other._permut_cover

    def __ne__(self, other):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: p1 = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: p2 = iet.GeneralizedPermutation('a a b',' b c c')
            sage: p3 = iet.GeneralizedPermutation('a a b b', 'c c')
            sage: p1.cover(['(1)', '(1)', '(1)']) != p2.cover(['(1)', '(1)', '(1)'])
            False
            sage: p1.cover(['(1,2)', '', '']) != p1.cover(['(1,2)', '(1,2)', ''])
            True
            sage: p1.cover(['(1)', '(1)', '(1)']) != p3.cover(['(1)', '(1)', '(1)'])
            True

            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: G = Zmod(5)
            sage: p.regular_cover(G, [0, 1, 2]) != p.regular_cover(G, [0, 1, 2])
            False
            sage: p.regular_cover(G, [0, 1, 2]) != p.regular_cover(G, [1, 2, 0])
            True
        """
        return type(self) != type(other) or \
               self._base != other._base or \
               self.degree() != other.degree() or \
               self._permut_cover != other._permut_cover

    def __len__(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: p1 = iet.Permutation('a b c', 'c b a')
            sage: p2 = p1.cover(['(1,2)', '(1,3)', '(2,3)'])
            sage: len(p2)
            3
        """
        return len(self._base)

    def __getitem__(self,i):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: p1 = iet.Permutation('a b c', 'c b a')
            sage: p2 = p1.cover(['(1,2)', '(1,3)', '(2,3)'])
            sage: p2[0]
            ['a', 'b', 'c']
            sage: p2[1]
            ['c', 'b', 'a']
        """
        return self._base[i]

    def base(self):
        r"""
        Return the combinatorial data corresponding to the base of this cover

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: p = iet.GeneralizedPermutation('a b b','c c a')
            sage: c = p.cover(['(1,2,3)','(1,3)','(1,2)'])
            sage: q = c.base()
            sage: q
            a b b
            c c a
            sage: q == p
            True
        """
        return self._base

    def __copy__(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: p1 = iet.Permutation('a b c', 'c b a')
            sage: p2 = p1.cover(['(1,2)', '(1,3)', '(2,3)'])
            sage: p2 == p2.__copy__()
            True
        """
        q = PermutationCover(self._base.__copy__(), self.degree(), self._permut_cover[:])
        return q

    def covering_data(self, label):
        r"""
        Returns the permutation associated to the given ``label``.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = QuadraticStratum([1,1,-1,-1]).components()[0].permutation_representative()
            sage: pc = p.orientation_cover()
            sage: pc
            Covering of degree 2 of the permutation:
            0 1 2 3 3
            2 1 4 4 0

            sage: pc.covering_data(1)
            ()
            sage: pc.covering_data(3)
            (1,2)
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        S = SymmetricGroup(self.degree())
        return S([i+1 for i in self.covering_data_tuple(label)])

    def covering_data_tuple(self, label):
        r"""
        Returns the permutation associated to the given ``label`` as a tuple on
        `\{0, 1, \ldots, d-1\}` where `d` is the degree of the cover.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = QuadraticStratum([1,1,-1,-1]).components()[0].permutation_representative()
            sage: pc = p.orientation_cover()
            sage: pc
            Covering of degree 2 of the permutation:
            0 1 2 3 3
            2 1 4 4 0

            sage: pc.covering_data_tuple(1)
            [0, 1]
            sage: pc.covering_data_tuple(3)
            [1, 0]
        """
        return self._permut_cover[self._base.alphabet().rank(label)]

    def interval_diagram(self, sign=False):
        r"""
        Return the interval diagram.

        This is the permutation induced on the subintervals of this cover while
        we turn around the singularities. This is mainly used to compute the
        stratum of this permutation.

        OUTPUT:

        A list of lists of pairs ``(label, index of interval)`` if
        ``sign=False`` or a list of triples ``(label, index of interval,
        sign)``.

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: c = p.cover(['(1,2)','(1,3)','(1,4)'])
            sage: c.interval_diagram()
            [[('a', 1), ('a', 0)],
             [('a', 2)],
             [('a', 3)],
             [('a', 0), ('b', 1), ('a', 1), ('b', 0), ('a', 2), ('b', 2)],
             [('a', 3), ('b', 3)],
             [('b', 2), ('c', 2), ('b', 0), ('c', 0), ('b', 3), ('c', 3)],
             [('b', 1), ('c', 1)],
             [('c', 3), ('c', 0)],
             [('c', 1)],
             [('c', 2)]]

            sage: c.interval_diagram(sign=True)
                [[('a', 1, -1), ('a', 0, -1)],
                 [('a', 2, -1)],
                 [('a', 3, -1)],
                 [('a', 0, 1), ('b', 1, 1), ..., ('b', 2, 1)],
                 [('a', 3, 1), ('b', 3, 1)],
                 [('b', 2, -1), ('c', 2, 1), ..., ('c', 3, 1)],
                 [('b', 1, -1), ('c', 1, 1)],
                 [('c', 3, -1), ('c', 0, -1)],
                 [('c', 1, -1)],
                 [('c', 2, -1)]]
        """
        base_diagram = self._base.interval_diagram(glue_ends=False, sign=True)
        singularities = []

        alphabet = self._base.alphabet()
        rank = alphabet.rank

        def perm(ss, label):
            return self._permut_cover[rank(label)] if ss == 1 else \
                self._inv_permut_cover[rank(label)]

        for orbit in base_diagram:
            cover_copies = set(range(self.degree()))
            while cover_copies:
                d = d_init = cover_copies.pop()
                singularity = []
                while True:
                    # lift a loop from downstair
                    for base_singularity in orbit:
                        label, s = base_singularity
                        if s == -1:
                            dd = perm(s,label)[d]
                        else:
                            dd = d
                        singularity.append((label,dd,s) if sign else (label,dd))
                        d = perm(s,label)[d]

                    # if it closes, break the loop
                    if d == d_init:
                        break
                    else:
                        cover_copies.remove(d)

                singularities.append(singularity)

        return singularities

    def _delta2(self):
        r"""
        Matrix `\delta_2: C_2 -> C_1` from the `2`-cycles to the `1`-cycles.

        The matrix acts on the left. The basis of `C1` (corresponding to
        columns) is ordered first by index in the cover and then by index in the
        base.

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: p = iet.GeneralizedPermutation('a b b', 'c c a')
            sage: c = p.cover(['(1)', '(1)', '(1)'])
            sage: c._delta2()
            [0 0 0]
            sage: c = p.cover(['(1,2)', '(1,3)', '(1,4)'])
            sage: c._delta2()
            [ 1  1  1 -1  0  0  0 -1  0  0  0 -1]
            [-1  0  0  1  0  0  0  0  0  0  0  0]
            [ 0 -1  0  0  0  0  0  1  0  0  0  0]
            [ 0  0 -1  0  0  0  0  0  0  0  0  1]
        """
        gens = [(i,a) for i in range(self.degree()) for a in self._base.letters()]
        gen_indices = {g:i for i,g in enumerate(gens)}
        B = []
        p = self._base
        signs = p._canonical_signs()[1]
        for k in range(self.degree()):
            border = [0] * len(gens)
            for i in range(2):
                for j in range(len(p[i])):
                    if signs[i][j] == -1:
                        perm_cover = self.covering_data_tuple(p[i][j])
                        border[gen_indices[(perm_cover.index(k),p[i][j])]] += -1
                    elif signs[i][j] == 1:
                        border[gen_indices[(k,p[i][j])]] += 1
                    else:
                        RuntimeError
            B.append(border)

        from sage.matrix.constructor import matrix
        from sage.rings.integer_ring import ZZ
        return matrix(ZZ, B)

    def _delta1(self):
        r"""
        Matrix `\delta_1: C_1 -> C_0` from the `1`-chains to the `0`-chains

        The matrix acts on the left. The basis of `C_1` (corresponding to
        rows) is ordered first by index in the cover and then by index in the
        base. The basis of `C_0` is ordered as given by the method
        :meth:`interval_diagram`.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('a b b','c c a')
            sage: c = p.cover(['(1)', '(1)', '(1)'])
            sage: m = c._delta1()
            sage: m
            [-1  1  0  0]
            [ 1  0 -1  0]
            [ 0  1  0 -1]
            sage: m.ncols() == len(c.profile())
            True
            sage: m.nrows() == len(c) * c.degree()
            True

            sage: c = p.cover(['(1,2)', '(1,3)', '(1,4)'])
            sage: (c._delta2() * c._delta1()).is_zero()
            True

            sage: a,b = QuaternionGroup().gens()
            sage: p = iet.Permutation('a b', 'b a')
            sage: c = p.cover([a, b])
            sage: assert c._delta1().nrows() == 2 * 8  # number of edges
            sage: assert c._delta1().ncols() == 4      # number of vertices
            sage: assert c._delta2().nrows() == 8      # number of faces
            sage: assert c._delta2().ncols() == 2 * 8  # number of edges
            sage: (c._delta2() * c._delta1()).is_zero()
            True

            sage: p = iet.GeneralizedPermutation('a a b c d b e', 'e f d c f')
            sage: c = p.cover(['(1,2,3,4)', '(1,3)', '(1,2)', '(2,3)', '(1,3,4)', '()'])
            sage: (c._delta2() * c._delta1()).is_zero()
            True
        """
        singularities = self.interval_diagram(sign=True)
        sing_to_index = {s:i for i,sing in enumerate(singularities) for s in sing}
        nb_sing = len(singularities)
        borders = []
        for d in range(self.degree()):
            for a in self._base.alphabet():
                border_side = [0] * nb_sing
                border_side[sing_to_index[(a,d,+1)]] += 1
                border_side[sing_to_index[(a,d,-1)]] += -1
                borders.append(border_side)

        from sage.matrix.constructor import matrix
        from sage.rings.integer_ring import ZZ
        return matrix(ZZ, borders)

    def profile(self):
        r"""
        Return the profile of the surface.

        The *profile* of a translation surface is the list of angles of
        singularities in the surface divided by pi.

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: p = iet.Permutation('a b', 'b a')
            sage: p.cover(['(1,2)', '(1,3)']).profile()
            [6]
            sage: p.cover(['(1,2,3)','(1,4)']).profile()
            [6, 2]
            sage: p.cover(['(1,2,3)(4,5,6)','(1,4,7)(2,5)(3,6)']).profile()
            [6, 2, 2, 2, 2]

            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: p.cover(['(1,2)', '()', '(1,2)']).profile()
            [2, 2, 2, 2]
        """
        from sage.combinat.partition import Partition
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup

        s = []

        base_diagram = self._base.interval_diagram(sign=True, glue_ends=True)
        p_id = SymmetricGroup(self.degree()).identity()
        for orbit in base_diagram:
            flat_orbit = []
            for x in orbit:
                if isinstance(x[0], tuple):
                    flat_orbit.extend(x)
                else:
                    flat_orbit.append(x)
            p = p_id
            for lab, sign in flat_orbit:
                q = self.covering_data(lab)
                if sign == -1: q = q.inverse()
                p = p*q
            for c in p.cycle_type():
                s.append(len(orbit)*c)

        return Partition(sorted(s,reverse=True))

    def is_orientable(self):
        r"""
        Test whether this permutation cover has an orientable foliation.

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: from itertools import product
            sage: it = iter(product(('()', '(1,2)'), repeat=3))
            sage: next(it)
            ('()', '()', '()')

            sage: for cov in it:
            ....:     c = p.cover(cov)
            ....:     print("%28s %s" % (cov, c.is_orientable()))
            ('()', '()', '(1,2)') False
            ('()', '(1,2)', '()') False
            ('()', '(1,2)', '(1,2)') False
            ('(1,2)', '()', '()') False
            ('(1,2)', '()', '(1,2)') True
            ('(1,2)', '(1,2)', '()') False
            ('(1,2)', '(1,2)', '(1,2)') False
        """
        from surface_dynamics.interval_exchanges.template import PermutationIET
        if isinstance(self._base, PermutationIET):
            return True
        elif any(x%2 for x in self.profile()):
            return False
        else:
            # here we unfold the holonomy, i.e. try to orient each copy with +1
            # or -1
            signs = [+1] + [None] * (self.degree() - 1)
            todo = [0]
            A = self._base.alphabet()
            p0 = set(map(A.rank, self._base[0]))
            p1 = set(map(A.rank, self._base[1]))
            inv_letters = p0.symmetric_difference(p1)
            while todo:
                i = todo.pop()
                s = signs[i]
                assert s is not None
                for j in inv_letters:
                    ii = self._permut_cover[j][i]
                    if signs[ii] is None:
                        signs[ii] = -s
                        todo.append(ii)
                    elif signs[ii] != -s:
                        return False
            return True

    @cached_method
    def monodromy(self):
        r"""
        Return the monodromy of this covering.

        That it to say the permutation group generated by the action of the
        fundamental group.

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: p.cover(['(1,2,3)', '(1,3,2)', '']).monodromy()
            Permutation Group with generators [(1,2,3), (1,3,2), ()]
            sage: p.cover(['(1,2)', '(1,3)', '']).monodromy()
            Permutation Group with generators [(1,2), (1,3), ()]
        """
        return PermutationGroup([self.covering_data(a) for a in self._base.letters()], canonicalize=False)

    @cached_method
    def automorphism_group(self):
        r"""
        Return the Deck group of the cover.

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: p.cover(['(1,2,3)', '(1,3,2)', '']).automorphism_group()
            Permutation Group with generators [(1,2,3), (1,3,2)]
            sage: p.cover(['(1,2)', '(1,3)', '']).automorphism_group()
            Permutation Group with generators [()]
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup

        Sd = SymmetricGroup(self.degree())
        G = libgap.Subgroup(Sd, [self.covering_data(a) for a in self._base.letters()])
        C = libgap.Centralizer(Sd, G)

        return PermutationGroup(libgap.GeneratorsOfGroup(C).sage())

    group = automorphism_group

    def stratum(self, fake_zeros=True):
        r"""
        Return the stratum of the covering translation surface.

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')

            sage: p.cover(['(1,2)', '()', '(1,2)']).stratum()
            H_1(0^4)
            sage: p.cover(['(1,2)', '()', '(1,2)']).stratum(fake_zeros=False)
            H_1(0)

            sage: p.cover(['(1,2)', '(1,2)', '(1,2)']).stratum()
            Q_0(0^2, -1^4)
            sage: p.cover(['(1,2)', '(1,2)', '(1,2)']).stratum(fake_zeros=False)
            Q_0(-1^4)

        TESTS::

            sage: from surface_dynamics import *
            sage: p1 = Permutation('(1,2,3)(4,5,6,7,8,9)(10,11)')
            sage: p2 = Permutation('(1,3,5,7,9)(6,10)')
            sage: Origami(p1, p2).stratum()
            H_6(10)
            sage: iet.Permutation('a b', 'b a').cover([p1, p2]).stratum()
            H_6(10)

            sage: p1 = Permutation('(1,2)(3,4)(5,6)(7,8,9,10)')
            sage: p2 = Permutation('(2,3)(4,5)(6,7)')
            sage: Origami(p1, p2).stratum()
            H_4(3^2)
            sage: iet.Permutation('a b', 'b a').cover([p1, p2]).stratum(fake_zeros=False)
            H_4(3^2)

            sage: p1 = Permutation('(1,2,3,4)(5,6,7,8,9,10)')
            sage: p2 = Permutation('(4,5)(8,10)')
            sage: Origami(p1, p2).stratum()
            H_3(2, 1^2)
            sage: iet.Permutation('a b', 'b a').cover([p1, p2]).stratum(fake_zeros=False)
            H_3(2, 1^2)
        """
        Z = [x-2 for x in self.profile() if fake_zeros or x != 2]
        if not Z:
            Z = [0]
        if self.is_orientable():
            from surface_dynamics.flat_surfaces.abelian_strata import AbelianStratum
            return AbelianStratum([z//2 for z in Z])
        else:
            from surface_dynamics.flat_surfaces.quadratic_strata import QuadraticStratum
            return QuadraticStratum(Z)

    def genus(self):
        r"""
        Genus of the covering translation surface

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: p.cover(['(1,2)', '()', '(1,2)']).genus()
            1
            sage: p.cover(['(1,2)', '(1,2)', '(1,2)']).genus()
            0

        TESTS::

            sage: from surface_dynamics import *
            sage: o = AbelianStratum([1,2,3,4]).one_component().one_origami()
            sage: assert(o.genus() == AbelianStratum([1,2,3,4]).genus())
            sage: qc = QuadraticStratum([1,2,3,4,-1,-1]).one_component()
            sage: p = qc.permutation_representative()
            sage: assert(p.orientation_cover().genus() == qc.orientation_cover_component().genus())
        """
        p = self.profile()
        return Integer((sum(p)-2*len(p))/4+1)

    def isotypic_projection_matrix(self, i, floating_point=False):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: p = iet.GeneralizedPermutation('a b b','c c a')
            sage: c = p.cover(['(1,2,3)','(1,3,2)','()'])
            sage: c.isotypic_projection_matrix(0)
            [1/3   0   0 1/3   0   0 1/3   0   0]
            [  0 1/3   0   0 1/3   0   0 1/3   0]
            [  0   0 1/3   0   0 1/3   0   0 1/3]
            [1/3   0   0 1/3   0   0 1/3   0   0]
            [  0 1/3   0   0 1/3   0   0 1/3   0]
            [  0   0 1/3   0   0 1/3   0   0 1/3]
            [1/3   0   0 1/3   0   0 1/3   0   0]
            [  0 1/3   0   0 1/3   0   0 1/3   0]
            [  0   0 1/3   0   0 1/3   0   0 1/3]
            sage: c.isotypic_projection_matrix(0, True)
            array([[0.33333333, 0.        , 0.        , 0.33333333, 0.        ,
                    0.        , 0.33333333, 0.        , 0.        ],
                   [0.        , 0.33333333, 0.        , 0.        , 0.33333333,
                    0.        , 0.        , 0.33333333, 0.        ],
                   [0.        , 0.        , 0.33333333, 0.        , 0.        ,
                    0.33333333, 0.        , 0.        , 0.33333333],
                   [0.33333333, 0.        , 0.        , 0.33333333, 0.        ,
                    0.        , 0.33333333, 0.        , 0.        ],
                   [0.        , 0.33333333, 0.        , 0.        , 0.33333333,
                    0.        , 0.        , 0.33333333, 0.        ],
                   [0.        , 0.        , 0.33333333, 0.        , 0.        ,
                    0.33333333, 0.        , 0.        , 0.33333333],
                   [0.33333333, 0.        , 0.        , 0.33333333, 0.        ,
                    0.        , 0.33333333, 0.        , 0.        ],
                   [0.        , 0.33333333, 0.        , 0.        , 0.33333333,
                    0.        , 0.        , 0.33333333, 0.        ],
                   [0.        , 0.        , 0.33333333, 0.        , 0.        ,
                    0.33333333, 0.        , 0.        , 0.33333333]])

        TESTS::

            sage: from surface_dynamics import *
            sage: import numpy as np
            sage: p = iet.Permutation('a b', 'b a')

            sage: Q = QuaternionGroup()
            sage: a,b = Q.gens()
            sage: pp = p.regular_cover(Q, [a, b])
            sage: for i in range(5):
            ....:     m1 = pp.isotypic_projection_matrix(i, True)
            ....:     m2 = pp.isotypic_projection_matrix(i, False).n().numpy()
            ....:     assert m1.shape == m2.shape and np.allclose(m1, m2)

            sage: G = SL(2, 4)
            sage: a, b = G.gens()
            sage: pp = p.regular_cover(G, [a, b])
            sage: for i in range(5):
            ....:     m1 = pp.isotypic_projection_matrix(i, True)
            ....:     m2 = pp.isotypic_projection_matrix(i, False).n().numpy()
            ....:     assert m1.shape == m2.shape and np.allclose(m1, m2)
        """
        from sage.matrix.special import identity_matrix
        from surface_dynamics.misc.group_representation import isotypic_projection_matrix
        if floating_point:
            m = isotypic_projection_matrix(
                   self.automorphism_group(),
                   self.degree(),
                   self._real_characters()[0][i],
                   self._real_characters()[1][i],
                   self._cc_mats(),
                   True)
            k = len(self._base) ** 2
            return np.hstack(np.hstack(np.tensordot(m, np.identity(len(self._base)), 0)))
        else:
            m = isotypic_projection_matrix(
                   self.automorphism_group(),
                   self.degree(),
                   self._real_characters()[0][i],
                   self._real_characters()[1][i],
                   self._cc_mats(),
                   False)
            return m.tensor_product(identity_matrix(len(self._base)), subdivide=False)

    @cached_method
    def _cc_mats(self):
        r"""
        The projection given by the conjugacy class.

        This is cached to speed up the computation of the projection matrices.
        See :meth:`~flatsurf.misc.group_representation.conjugacy_class_matrix`
        for more information.

        TESTS::

            sage: from surface_dynamics import *
            sage: p = iet.GeneralizedPermutation('a b b','c c a')
            sage: pp = p.cover(['(1,2,3)','(1,3,2)','()'])
            sage: pp._cc_mats()
            (
            [1 0 0]  [0 1 0]  [0 0 1]
            [0 1 0]  [0 0 1]  [1 0 0]
            [0 0 1], [1 0 0], [0 1 0]
            )

            sage: G = SL(2, 2)
            sage: a = G([[1,1],[0,1]])
            sage: b = G([[1,0],[1,1]])
            sage: c = G([[0,1],[1,0]])
            sage: pp = p.regular_cover(G, [a,b,c])
            sage: pp._cc_mats()
            (
            [1 0 0 0 0 0]  [0 1 0 1 0 1]  [0 0 1 0 1 0]
            [0 1 0 0 0 0]  [1 0 1 0 1 0]  [0 0 0 1 0 1]
            [0 0 1 0 0 0]  [0 1 0 1 0 1]  [1 0 0 0 1 0]
            [0 0 0 1 0 0]  [1 0 1 0 1 0]  [0 1 0 0 0 1]
            [0 0 0 0 1 0]  [0 1 0 1 0 1]  [1 0 1 0 0 0]
            [0 0 0 0 0 1], [1 0 1 0 1 0], [0 1 0 1 0 0]
            )
        """
        from surface_dynamics.misc.group_representation import conjugacy_class_matrix, \
                regular_conjugacy_class_matrix

        G = self.group()
        if G is self.automorphism_group():
            mats = []
            for cl in libgap(G).ConjugacyClasses():
                m = conjugacy_class_matrix(cl, self.degree())
                m.set_immutable()
                mats.append(m)
        else:
            G = self._gap_grp
            Glist = G.AsList()
            mats = []
            for cl in G.ConjugacyClasses():
                m = regular_conjugacy_class_matrix(cl, G)
                m.set_immutable()
                mats.append(m)
        return tuple(mats)

    @cached_method
    def _real_characters(self):
        r"""
        The real characters of the automorphism group

        OUTPUT:

        - table of characters

        - degrees

        For more information see
        :func:`~flatsurf.misc.group_representation.real_characters`.

        TESTS::

            sage: from surface_dynamics import *
            sage: p = iet.GeneralizedPermutation('a b b','c c a')
            sage: c = p.cover(['(1,2,3)','(1,3,2)','()'])
            sage: c._real_characters()
            ([(1, 1, 1), (2, -1, -1)], [1, 1])
        """
        from surface_dynamics.misc.group_representation import real_characters
        return real_characters(self.group())

    def isotypic_projectors(self, characters=None, floating_point=True):
        r"""
        Return the list of projectors on isotypical components in the canonical basis
        as well as the half ranks of the projection in absolute homology.

        INPUT:

        - ``characters`` - (default ``None``), if set to ``None`` compute projectors for
          all characters, otherwise if set to a list of characters return only the projectors
          for these characters.

        - ``floating_point`` - (default ``True``) if set to ``True`` then all computations are
          performed over floating point numbers. Otherwise, exact cyclotomic elements are kept
          until the very last steps of the computation.

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: p = iet.Permutation('a b', 'b a')
            sage: Q = QuaternionGroup()
            sage: a,b = Q.gens()
            sage: pp = p.regular_cover(Q, [a, b])
            sage: X = pp.isotypic_projectors()
            sage: X[0][0]
            array([[0.125, 0.   , 0.125, 0.   , 0.125, 0.   , 0.125, 0.   , 0.125,
                    0.   , 0.125, 0.   , 0.125, 0.   , 0.125, 0.   ],
                   [0.   , 0.125, 0.   , 0.125, 0.   , 0.125, 0.   , 0.125, 0.   ,
                    0.125, 0.   , 0.125, 0.   , 0.125, 0.   , 0.125],
            ...
                    [0.   , 0.125, 0.   , 0.125, 0.   , 0.125, 0.   , 0.125, 0.   ,
                     0.125, 0.   , 0.125, 0.   , 0.125, 0.   , 0.125]])
            sage: import numpy as np
            sage: Y = pp.isotypic_projectors(floating_point=False)
            sage: np.allclose(X[0], Y[0])
            True

            sage: X[1]
            [1, 0, 0, 0, 2]
            sage: sum(X[1]) == pp.genus()
            True

        TESTS::

            sage: from surface_dynamics import iet
            sage: p = iet.GeneralizedPermutation('c a a', 'b b c', alphabet='abc')
            sage: def cyclic(n,a):
            ....:     return [(i+a)%n + 1 for i in range(n)]
            sage: def cyclic_cover(n, a, b, c):
            ....:     return p.cover([cyclic(n,c), cyclic(n,a), cyclic(n, b)])
            sage: def cyclic_cover_regular(n, a, b, c):
            ....:     return p.regular_cover(Zmod(n), [c, a, b])

            sage: for (dat, ans) in [((7,1,1,2), [0,2,2,2]),
            ....:                    ((7,1,3,3), [0,1,1,1]),
            ....:                    ((8,1,2,4), [0,0,1,2,2])]:
            ....:     c1 = cyclic_cover(*dat)
            ....:     c2 = cyclic_cover_regular(*dat)
            ....:     assert c1.isotypic_projectors(floating_point=True)[1] == ans
            ....:     assert c1.isotypic_projectors(floating_point=False)[1] == ans
            ....:     assert c2.isotypic_projectors(floating_point=True)[1] == ans
            ....:     assert c2.isotypic_projectors(floating_point=False)[1] == ans
            ....:     assert c1.genus() == sum(ans)
        """
        if characters is None or characters is True:
            characters = self._real_characters()[0]
        elif isinstance(characters, (tuple,list)):
            if isinstance(characters[0], (tuple,list)):
                characters = map(tuple, characters)
            else:
                characters = [tuple(characters)]
                flatten = True

        all_characters = self._real_characters()[0]
        for c in characters:
            if c not in all_characters:
                raise ValueError("{} is an invalid character".format(c))

        nc = len(characters)
        na = len(all_characters)
        dim = len(self._base) * self.degree()
        projections = np.zeros((nc, dim, dim))
        dimensions = list(range(nc))
        ii = 0

        H = self._delta1().left_kernel().matrix()
        B = self._delta2()
        for i_char in range(na):
            if all_characters[i_char] not in characters:
                continue

            if floating_point:
                M = self.isotypic_projection_matrix(i_char, floating_point=True)
                dimension = np.linalg.matrix_rank(np.matmul(H.numpy(), M)) - \
                            np.linalg.matrix_rank(np.matmul(B.numpy(), M))
            else:
                M = self.isotypic_projection_matrix(i_char, floating_point=False)
                dimension = (H*M).rank() - (B*M).rank()

            if dimension%2:
                raise RuntimeError("get a subspace of odd dimension\n" +
                        "  i_char = {}\n".format(i_char) +
                        "  dim    = {}".format(dimension))
            dimensions[ii] = dimension//2

            if floating_point:
                projections[ii,:,:] = np.real_if_close(np.matrix.transpose(M))
            else:
                for i in range(dim):
                    for j in range(dim):
                        projections[ii,i,j] = float(M[j][i])

            ii += 1

        return projections, dimensions

    def lyapunov_exponents_H_plus(self, nb_vectors=None, nb_experiments=10,
                                  nb_iterations=65536, output_file=None,
                                  return_speed=False, isotypic_decomposition=False,
                                  return_char=False, verbose=False, floating_point=True):
        r"""
        Compute the H^+ Lyapunov exponents in  the covering locus.

        This method calls a C-library that performs renormalization (i.e. Rauzy
        induction and orthogonalization of vectors under the Kontsevich-Zorich
        cocycle). The computation might be significantly faster if
        ``nb_vectors=1`` (or if it is not provided but genus is 1) as in this
        case no orthogonalization is needed.

        INPUT:

         - ``nb_vectors`` -- the number of exponents to compute. The number of
           vectors must not exceed the dimension of the space!

         - ``nb_experiments`` -- the number of experiments to perform. It might
           be around 100 (default value) in order that the estimation of
           confidence interval is accurate enough.

         - ``nb_iterations`` -- the number of iteration of the Rauzy-Zorich
           algorithm to perform for each experiments. The default is 2^15=32768
           which is rather small but provide a good compromise between speed and
           quality of approximation.

         - ``output_file`` -- if provided (as a file object or a string) output
           the additional information in the given file rather than on the
           standard output.

         - ``return_speed`` -- whether or not return the lyapunov exponents list
           in a pair with the speed of the geodesic.

         - ``isotypic_decomposition`` -- either a boolean or a character or a
           list of characters.

         - ``return_char`` -- whether or not return the character corresponding to
           the isotypic component.

         - ``verbose`` -- if ``True`` provide additional information rather than
           returning only the Lyapunov exponents (i.e. elapsed time, confidence
           intervals, ...)

         - ``float`` -- whether the isotypical decomposition and projectors are computed
           over exact or floating point numbers

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: q = QuadraticStratum([1,1,-1,-1]).one_component()
            sage: q.lyapunov_exponents_H_plus(nb_iterations=2**19)  # abs tol 0.1
            [0.666666]
            sage: p = q.permutation_representative(reduced=False).orientation_cover()
            sage: c = p.lyapunov_exponents_H_plus(isotypic_decomposition=True, nb_iterations=2**19)[0]
            sage: c[0]  # abs tol 0.1
            1.000000

            sage: p = iet.GeneralizedPermutation('e a a', 'b b c c d d e')
            sage: p.stratum()
            Q_0(1, -1^5)
            sage: p.alphabet()
            {'e', 'a', 'b', 'c', 'd'}
            sage: c = p.cover(['()', '(1,2)', '()', '(1,2)', '(1,2)'])
            sage: c.stratum(fake_zeros=True)
            Q_1(1^2, 0^4, -1^2)
            sage: c.lyapunov_exponents_H_plus(nb_iterations=2**19)  # abs tol 0.1
            [0.666666]

        Some cyclic covers (see [EskKonZor11]_ for the formulas)::

            sage: p = iet.GeneralizedPermutation('c a a', 'b b c', alphabet='abc')
            sage: def cyclic(n,a):
            ....:     return [(i+a)%n + 1 for i in range(n)]
            sage: def cyclic_cover(n, a, b, c):
            ....:     return p.cover([cyclic(n,c), cyclic(n,a), cyclic(n, b)])
            sage: def cyclic_cover_regular(n, a, b, c):
            ....:     G = groups.misc.MultiplicativeAbelian([n])
            ....:     x, = G.gens()
            ....:     return p.regular_cover(G, [x**c, x**a, x**b])

            sage: c = cyclic_cover(7,1,1,2)
            sage: lexp = c.lyapunov_exponents_H_plus(isotypic_decomposition=True, nb_iterations=2**19)
            sage: lexp  # abs tol 0.1
            [[],
             [0.2857, 0.2857],
             [0.5714, 0.5714],
             [0.2857, 0.2857]]

            sage: c = cyclic_cover(7, 1, 2, 3)
            sage: lexp = c.lyapunov_exponents_H_plus(isotypic_decomposition=True, return_char=True, nb_iterations=2**19)
            sage: lexp[0]
            ([], (1, 1, 1, 1, 1, 1, 1))
            sage: lexp[1][0]  # abs tol 0.1
            [0.2857, 0.2857]
            sage: lexp[1][1]
            (2, E(7) + E(7)^6, ..., E(7) + E(7)^6)
            sage: lexp[2][0]  # abs tol 0.1
            [0.2857, 0.2857]
            sage: lexp[2][1]
            (2, E(7)^2 + E(7)^5, ...,  E(7)^2 + E(7)^5)
            sage: lexp[3][0]  # abs tol 0.1
            [0.5714, 0.5714]
            sage: lexp[3][1]
            (2, E(7)^3 + E(7)^4, ..., E(7)^3 + E(7)^4)

        The Eierlegendewollmilchsau as a quaternionic cover of the once
        punctured torus::

            sage: p = iet.Permutation('a b', 'b a')
            sage: Q = QuaternionGroup()
            sage: a,b = Q.gens()
            sage: c = p.cover([a, b])
            sage: c.lyapunov_exponents_H_plus(nb_iterations=2**19)  # abs tol 0.05
            [1.0, 0.0, 0.0]
        """
        import surface_dynamics.interval_exchanges.lyapunov_exponents as lyapunov_exponents
        import time
        from surface_dynamics.misc.statistics import mean_and_std_dev
        from sage.matrix.special import zero_matrix

        n = len(self)

        if nb_vectors is None:
            nb_vectors = self.genus()

        if output_file is None:
            from sys import stdout
            output_file = stdout
        elif isinstance(output_file, str):
            output_file = open(output_file, "w")

        nb_vectors = int(nb_vectors)
        nb_experiments = int(nb_experiments)
        nb_iterations = int(nb_iterations)

        if nb_vectors < 0:
            raise ValueError("the number of vectors must be positive")
        if nb_vectors == 0:
            return []
        if nb_experiments <= 0:
            raise ValueError("the number of experiments must be positive")
        if nb_iterations <= 0:
            raise ValueError("the number of iterations must be positive")

        if verbose:
            output_file.write("Stratum: {}\n".format(self.stratum()))

        k = int(len(self[0]))
        gp = list(range(2*n))
        twin = list(range(2*n))
        base_twin = self._base.twin_list()
        rank = self._base.alphabet().rank
        for i in range(2):
            for j,a in enumerate(self._base[i]):
                gp[i*k+j] = rank(a)
                ii,jj = base_twin[i][j]
                twin[i*k+j] = ii*k+jj

        flatten = False
        if isotypic_decomposition:
            projections, dimensions = self.isotypic_projectors(characters=isotypic_decomposition,
                                                               floating_point=True)
        else:
            projections = None
            dimensions = [nb_vectors]

        t0 = time.time()
        res = lyapunov_exponents.lyapunov_exponents_H_plus_cover(
            gp, k, twin, self._permut_cover, nb_experiments, nb_iterations,
            dimensions, projections, None, verbose)
        t1 = time.time()

        res_final = []

        m,d = mean_and_std_dev(res[0])
        lexp = m

        if verbose:
            from math import log, floor, sqrt
            output_file.write("sample of %d experiments\n"%nb_experiments)
            output_file.write("%d iterations (~2^%d)\n"%(
                    nb_iterations,
                    floor(log(nb_iterations) / log(2))))
            output_file.write("elapsed time %s\n"%time.strftime("%H:%M:%S",time.gmtime(t1-t0)))
            output_file.write("Lexp Rauzy-Zorich: %f (std. dev. = %f, conf. rad. 0.01 = %f)\n"%(
                    m,d, 2.576*d/sqrt(nb_experiments)))

        if isotypic_decomposition:
            i_0 = 1
            nc = len(dimensions)
            for i_char in range(nc):
                res_int = []
                if verbose:
                    output_file.write("##### char_%d #####\n" % (i_char))
                    output_file.write("chi = {}\n".format(self._real_characters()[0][i_char]))
                    output_file.write("dim = {}\n".format(dimensions[i_char]))
                for i in range(i_0, i_0 + dimensions[i_char]):
                    m,d = mean_and_std_dev(res[i])
                    res_int.append(m)
                    if verbose:
                        output_file.write("theta_%d          : %f (std. dev. = %f, conf. rad. 0.01 = %f)\n"%(
                            i,m,d, 2.576*d/sqrt(nb_experiments)))
                res_final.append((res_int, self._real_characters()[0][i_char]) if return_char else res_int)
                i_0 += dimensions[i_char]

        else:
            for i in range(1,nb_vectors+1):
                m,d = mean_and_std_dev(res[i])
                if verbose:
                    output_file.write("theta%d           : %f (std. dev. = %f, conf. rad. 0.01 = %f)\n"%(i,m,d, 2.576*d/sqrt(nb_experiments)))
                res_final.append(m)

        res_final = res_final[0] if flatten else res_final
        if return_speed:
            return (lexp, res_final)
        else:
            return res_final

    def masur_polygon(self, lengths, heights):
        r"""
        Return the Masur polygon for the given ``lengths`` and ``heights``.

        EXAMPLES::

            sage: from surface_dynamics import iet

            sage: p = iet.Permutation('a b c', 'c b a').cover(['(1,2)','(1,3)','(1,4)'])
            sage: S = p.masur_polygon([1,4,2], [2,0,-1])  # optional: sage_flatsurf
            sage: S.stratum()                             # optional: sage_flatsurf
            H_4(3^2)
            sage: p.stratum()                             # optional: sage_flatsurf
            H_4(3^2)
        """
        base_ring, triangles, tops, bots, mids = self._base._masur_polygon_helper(lengths, heights)

        from flatsurf import MutableOrientedSimilaritySurface
        S = MutableOrientedSimilaritySurface(base_ring)
        n = len(self._base)
        nt = len(triangles)
        d = self._degree_cover
        for _ in range(d):
            for t in triangles:
                S.add_polygon(t)
        for i in range(n):
            p1, e1 = tops[i]
            p2, e2 = bots[self._base._twin[0][i]]
            perm = self._permut_cover[self._base._labels[0][i]]
            for j in range(d):
                jj = perm[j]
                S.glue((p1 + nt*j, e1), (p2 + nt*jj, e2))
        for i in range(0,len(mids),2):
            p1, e1 = mids[i]
            p2, e2 = mids[i+1]
            for j in range(d):
                S.glue((p1 + nt*j, e1), (p2 + nt*j, e2))
        S.set_immutable()
        return S

class RegularCover(PermutationCover):
    r"""
    An interval exchange permutation together with a regular covering data.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: p = iet.Permutation('a b c', 'c b a')
        sage: G = Zmod(3)
        sage: p.regular_cover(G, [1, 0, 2])
        Regular cover of degree 3 with group Multiplicative Abelian group isomorphic to C3 of the permutation:
        a b c
        c b a
        sage: G = AbelianGroup([2,4])
        sage: p.regular_cover(G, [(1,0), (0,0), (0,1)])
        Regular cover of degree 8 with group Multiplicative Abelian group isomorphic to C2 x C4 of the permutation:
        a b c
        c b a

        sage: G = GL(2, 3)
        sage: g1, g2 = G.gens()
        sage: p.regular_cover(G, [g1, g2, g1 * g2])
        Regular cover of degree 48 with group General Linear Group of degree 2 over Finite Field of size 3 of the permutation:
        a b c
        c b a

    An example using a semi-direct product built with GAP::

        sage: C1 = libgap.CyclicGroup(2^5)
        sage: C2 = libgap.CyclicGroup(2)
        sage: gens1 = libgap.GeneratorsOfGroup(C1)
        sage: gens2 = libgap.GeneratorsOfGroup(C2)
        sage: alpha = libgap.GroupHomomorphismByImages(C1, C1, [gens1[0]], [gens1[0]^(-1)])
        sage: phi = libgap.GroupHomomorphismByImages(C2, libgap.AutomorphismGroup(C1), [gens2[0]], [alpha])
        sage: G = libgap.SemidirectProduct(C2, phi, C1)
        sage: x = libgap.Image(libgap.eval("Embedding")(G, 2), gens1[0])
        sage: y = libgap.Image(libgap.eval("Embedding")(G, 1), gens2[0])
        sage: p = iet.Permutation('a b', 'b a')
        sage: pp = p.regular_cover(G, [x, y])
        sage: pp
        Regular cover of degree 64 with group <pc group of size 64 with 6 generators> of the permutation:
        a b
        b a
        sage: pp.isotypic_projectors()[1]
        [1, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        sage: pp.lyapunov_exponents_H_plus(isotypic_decomposition=True)  # not tested (too long)
        [[1.000],
         [],
         [],
         [],
         [0.5065137173173205, 0.5064242515256905],
         [0.25364063157823696, 0.25344903617958725],
         [0.7562155727190676, 0.7562259157608896],
         [0.1267784864027805, 0.12676433404844129],
         [0.3806034038221321, 0.38056305175722355],
         [0.6287868199859215, 0.6287551837770511],
         [0.8797056039630107, 0.8797049392731162],
         [0.06522712186764967, 0.06491210232533576],
         [0.1903390673409358, 0.19013594965668126],
         [0.3173095347090416, 0.3171889234380455],
         [0.44497853721724967, 0.444804933849061],
         [0.5662613367551405, 0.5662491535341443],
         [0.6906239438432811, 0.6905716259455005],
         [0.8165503109407333, 0.8164624067719244],
         [0.9385523471932377, 0.9383981565833286]]
    """
    def __init__(self, base, grp, elts, check=True):
        gap_elts = elts
        gap_grp = None
        if isinstance(grp, GapElement):
            if not grp.IsGroup() or not grp.IsFinite():
                raise ValueError("grp must be a finite multiplicative group")
            gap_grp = grp
            grp = GroupLibGAP(grp)
        elif grp in _AddGroups:
            # dummy conversion of additive groups to GAP multiplicative groups
            from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
            if isinstance(grp, IntegerModRing_generic):
                from sage.groups.abelian_gps.abelian_group import AbelianGroup
                n = grp.order()
                gap_grp = libgap.AbelianGroup([n])
                grp = AbelianGroup([n])
                g, = grp.gens()
                pows = elts
                elts = [g**int(i) for i in pows]
                g, = gap_grp.GeneratorsOfGroup()
                gap_elts = [g**int(i) for i in pows]
            else:
                raise NotImplementedError("additive group unsupported")
        elif grp in _MulGroups:
            # the gap conversion of multiplicative group is a bit broken
            from sage.groups.abelian_gps.abelian_group import AbelianGroup_class
            if isinstance(grp, AbelianGroup_class):
                gap_grp = libgap.AbelianGroup(grp.gens_orders())
                gens = gap_grp.GeneratorsOfGroup()
                elts = [grp(x) for x in elts]
                gap_elts = [prod(g**int(i) for g,i in zip(gens, x.exponents())) for x in elts]
        else:
            raise RuntimeError("only multiplicative groups are supported")

        self._grp = grp
        self._elts = [grp(x) for x in elts]
        self._invs = [~x for x in self._elts]

        self._gap_grp = libgap(grp) if gap_grp is None else gap_grp
        self._gap_elts = list(map(to_gap, gap_elts))

        # monodromy as permutations
        try:
            d = self._grp.cardinality()
        except AttributeError:
            # cardinality not implemented yet on GroupLibGAP
            d = self._grp._libgap.Size().sage()

        plist = []
        Glist = self._gap_grp.AsList()
        for g in self._gap_elts:
            p = libgap.ListPerm(libgap.Permutation(g, Glist, libgap.OnRight), d)
            plist.append([i.sage() - 1 for i in p])
        PermutationCover.__init__(self, base, d, plist)

    def group(self):
        return self._grp

    def __repr__(self):
        return 'Regular cover of degree {} with group {} of the permutation:\n{}'.format(self.degree(), self.group(), str(self._base))

    def monodromy(self):
        return self._elts

    def degree(self):
        try:
            d = self._grp.cardinality()
        except AttributeError:
            # cardinality not implemented yet on GroupLibGAP
            d = self._grp._libgap.Size().sage()
        return d

    def profile(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: from surface_dynamics.interval_exchanges.cover import PermutationCover
            sage: p = iet.GeneralizedPermutation('a a b c d', 'b e d c f f e')
            sage: G = Zmod(2)
            sage: pp = p.regular_cover(G, [1, 1, 1, 1, 1, 1])
            sage: pp.profile()
            [4, 4, 3, 3, 2, 2, 1, 1]
            sage: PermutationCover.profile(pp)
            [4, 4, 3, 3, 2, 2, 1, 1]

            sage: pp = p.regular_cover(G, [0, 1, 1, 1, 0, 0])
            sage: pp.profile()
            [6, 4, 4, 2, 1, 1, 1, 1]
            sage: PermutationCover.profile(pp)
            [6, 4, 4, 2, 1, 1, 1, 1]

            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: G = GL(2, GF(3))
            sage: g1 = G([[2,0],[0,1]])
            sage: g2 = G([[2,1],[2,0]])
            sage: pp = p.regular_cover(G, [g1, g2, g1])
            sage: pp.profile() == PermutationCover.profile(pp)
            True
            sage: pp = p.regular_cover(G, [g1, g1*~g2, g2])
            sage: pp.profile() == PermutationCover.profile(pp)
            True
        """
        from sage.combinat.partition import Partition

        d = self.degree()
        alphabet = self._base.alphabet()
        rank = alphabet.rank

        p = []
        for orbit in self._base.interval_diagram(sign=True, glue_ends=True):
            flat_orbit = []
            for x in orbit:
                if isinstance(x[0], tuple):
                    flat_orbit.extend(x)
                else:
                    flat_orbit.append(x)

            g = self._grp.one()
            for label, s in flat_orbit:
                if s == 1:
                    g = g * self._elts[rank(label)]
                elif s == -1:
                    g = g * self._invs[rank(label)]
                else:
                    raise RuntimeError

            a = len(orbit)
            o = g.order()
            assert d % o == 0
            p.extend([o * a] * (d // o))

        return Partition(sorted(p, reverse=True))

    def is_orientable(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics import *
            sage: p = iet.GeneralizedPermutation('a a b c d', 'b e d c f f e')
            sage: G = Zmod(2)
            sage: p.regular_cover(G, [1, 1, 1, 1, 1, 1]).is_orientable()
            False
            sage: p.regular_cover(G, [0, 1, 1, 1, 0, 0]).is_orientable()
            False
            sage: p.regular_cover(G, [1, 0, 0, 0, 1, 1]).is_orientable()
            True

            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: G = GL(2, GF(3))
            sage: g1 = G([[2,0],[0,1]])
            sage: g2 = G([[2,1],[2,0]])
            sage: p.regular_cover(G, [g1, g2, g1]).is_orientable()
            True
        """
        if isinstance(self._base, PermutationIET):
            return True
        if self.degree() % 2:
            return False

        G = self._gap_grp
        C = libgap.CyclicGroup(2)
        g, = C.GeneratorsOfGroup()
        o = C.Identity()
        p = self._base
        images = [None] * len(p)
        for i in range(2):
            for j in range(len(p._labels[i])):
                lab = p._labels[i][j]
                if images[lab] is None:
                    if p._twin[i][j][0] == i:
                        images[lab] = g
                    else:
                        images[lab] = o
        return libgap.GroupHomomorphismByImages(G, C, self._gap_elts, images) != libgap_fail
