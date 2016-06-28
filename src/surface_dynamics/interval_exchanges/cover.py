# coding: utf8
r"""
Permutation cover

This module deals with combinatorial data for covering of connected components of
strata of Abelian and quadratic differentials. The main feature is to be able to
compute Lyapunov exponents.
"""
# *************************************************************************
# Copyright (C) 2015-2016 Charles Fougeron <charlesfougeron@gmail.com>
#                         Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at http://www.gnu.org/licenses/
# *************************************************************************

from sage.structure.sage_object import SageObject

from sage.misc.cachefunc import cached_method
from sage.libs.gap.libgap import libgap
from sage.rings.integer import Integer


class PermutationCover(SageObject):
    r"""
    A permutation with covering data

    Let `\pi` be the combinatorial data of an interval exchange transformation
    (or linear involution) on the alphabet `{1, 2, \ldots, m\}`. A cover of
    degree `d` is given by a list of permutations `\sigma_i \in S_d` for each `i
    \in \{1, 2, \ldots, m\}`.

    In order to do so, each interval on the base surface should come with an
    orientation. This orientation is automaticaly choosen by convention
    according to a clockwise orientation of the surface. The two copies of any
    interval have to be oriented alternatively in this choosen orientation and
    in the oposite to it.

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

            sage: from surface_dynamics.all import *
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
        """
        if not base.is_irreducible():
            raise ValueError("the base must be irreducible")
        from surface_dynamics.misc.permutation import perms_are_transitive
        if not perms_are_transitive(perms):
            raise ValueError("the cover is not connected")

        from surface_dynamics.misc.permutation import perm_invert

        self._base = base.__copy__()
        self._degree_cover = degree
        self._permut_cover = perms[:]
        self._inv_permut_cover = [perm_invert(p) for p in perms]

    def __len__(self):
        r"""
        TESTS::

            sage: from surface_dynamics.all import *
            sage: p1 = iet.Permutation('a b c', 'c b a')
            sage: p2 = p1.cover(['(1,2)', '(1,3)', '(2,3)'])
            sage: len(p2)
            3
        """
        return len(self._base)

    def __repr__(self):
        r"""
        A representation of the generalized permutation cover.

        INPUT:

        - ``sep`` - (default: '\n') a separator for the two intervals

        OUTPUT:

        string -- the string that represents the permutation


        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p1 = iet.Permutation('a b c', 'c b a')
            sage: p1.cover(['(1,2)', '(1,3)', '(2,3)'])
            Covering of degree 3 of the permutation:
            a b c
            c b a
        """
        s = 'Covering of degree %i of the permutation:\n'%(self._degree_cover)
        s += str(self._base)
        return s

    def __getitem__(self,i):
        r"""
        TESTS::

            sage: from surface_dynamics.all import *
            sage: p1 = iet.Permutation('a b c', 'c b a')
            sage: p2 = p1.cover(['(1,2)', '(1,3)', '(2,3)'])
            sage: p2[0]
            ['a', 'b', 'c']
            sage: p2[1]
            ['c', 'b', 'a']
        """
        return self._base[i]

    def __eq__(self, other):
        r"""
        TESTS::

            sage: from surface_dynamics.all import *
            sage: p1 = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: p2 = iet.GeneralizedPermutation('a a b',' b c c')
            sage: p3 = iet.GeneralizedPermutation('a a b b', 'c c')
            sage: p1.cover(['', '', '']) == p2.cover(['', '', ''])
            True
            sage: p1.cover(['(1,2)', '', '']) == p1.cover(['(1,2)', '(1,2)', ''])
            False
            sage: p1.cover(['', '', '']) == p3.cover(['', '', ''])
            False
        """
        return type(self) == type(other) and \
               self._base == other._base and \
               self._degree_cover == other._degree_cover and \
               self._permut_cover == other._permut_cover

    def __ne__(self, other):
        r"""
        TESTS::

            sage: from surface_dynamics.all import *
            sage: p1 = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: p2 = iet.GeneralizedPermutation('a a b',' b c c')
            sage: p3 = iet.GeneralizedPermutation('a a b b', 'c c')
            sage: p1.cover(['', '', '']) != p2.cover(['', '', ''])
            False
            sage: p1.cover(['(1,2)', '', '']) != p1.cover(['(1,2)', '(1,2)', ''])
            True
            sage: p1.cover(['', '', '']) != p3.cover(['', '', ''])
            True
        """
        return type(self) != type(other) or \
               self._base != other._base or \
               self._degree_cover != other._degree_cover or \
               self._permut_cover != other._permut_cover

    def __copy__(self):
        r"""
        TESTS::

            sage: from surface_dynamics.all import *
            sage: p1 = iet.Permutation('a b c', 'c b a')
            sage: p2 = p1.cover(['(1,2)', '(1,3)', '(2,3)'])
            sage: p2 == p2.__copy__()
            True
        """
        q = PermutationCover(self._base.__copy__(), self._degree_cover, self._permut_cover[:])
        return q

    def base(self):
        r"""
        Return the combinatorial data corresponding to the base of this cover

        EXAMPLES::

            sage: from surface_dynamics.all import *
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

    def covering_data(self, label):
        r"""
        Returns the permutation associated to the given ``label``.

        EXAMPLES::

            sage: from surface_dynamics.all import *

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
        S = SymmetricGroup(self._degree_cover)
        return S([i+1 for i in self.covering_data_tuple(label)])

    def covering_data_tuple(self, label):
        r"""
        Returns the permutation associated to the given ``label`` as a tuple on
        `\{0, 1, \ldots, d-1\}` where `d` is the degree of the cover.

        EXAMPLES::

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *
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
            [[('a', 1, 1), ('a', 0, 1)],
             [('a', 2, 1)],
             [('a', 3, 1)],
             [('a', 0, 0), ..., ('b', 2, 0)],
             [('a', 3, 0), ('b', 3, 0)],
             [('b', 2, 1), ..., ('c', 3, 0)],
             [('b', 1, 1), ('c', 1, 0)],
             [('c', 3, 1), ('c', 0, 1)],
             [('c', 1, 1)],
             [('c', 2, 1)]]
        """
        twins, orientation = self._base._canonical_signs()
        base_diagram = self._base.interval_diagram(glue_ends=False, sign=True)
        singularities = []

        alphabet = self._base.alphabet()
        rank = alphabet.rank
        perm = lambda sign,label: self._inv_permut_cover[rank(label)] if sign else \
                                  self._permut_cover[rank(label)]

        for orbit in base_diagram:
            init_label = orbit[0][0]
            cover_copies = set(range(self._degree_cover))
            while cover_copies:
                d = d_init = cover_copies.pop()
                singularity = []
                while True:
                    # lift a loop from downstair
                    for base_singularity in orbit:
                        label,s = base_singularity
                        if s:
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

            sage: from surface_dynamics.all import *
            sage: p = iet.GeneralizedPermutation('a b b', 'c c a')
            sage: c = p.cover(['()', '()', '()'])
            sage: c._delta2()
            [0 0 0]
            sage: c = p.cover(['(1,2)', '(1,3)', '(1,4)'])
            sage: c._delta2()
            [ 1  1  1 -1  0  0  0 -1  0  0  0 -1]
            [-1  0  0  1  0  0  0  0  0  0  0  0]
            [ 0 -1  0  0  0  0  0  1  0  0  0  0]
            [ 0  0 -1  0  0  0  0  0  0  0  0  1]
        """
        gens = [(i,a) for i in range(self._degree_cover) for a in self._base.letters()]
        gen_indices = {g:i for i,g in enumerate(gens)}
        B = []
        p = self._base
        signs = p._canonical_signs()[1]
        for k in xrange(self._degree_cover):
            border = [0] * len(gens)
            for i in xrange(2):
                for j in xrange(len(p[i])):
                    if signs[i][j]:  # -1 sign
                        perm_cover = self.covering_data_tuple(p[i][j])
                        border[gen_indices[(perm_cover.index(k),p[i][j])]] += -1
                    else:            # +1 sign
                        border[gen_indices[(k,p[i][j])]] += 1
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

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('a b b','c c a')
            sage: c = p.cover(['', '', ''])
            sage: m = c._delta1()
            sage: m
            [-1  1  0  0]
            [ 1  0 -1  0]
            [ 0  1  0 -1]
            sage: m.ncols() == len(c.profile())
            True
            sage: m.nrows() == len(c) * c._degree_cover
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
        for d in range(self._degree_cover):
            for a in self._base.alphabet():
                border_side = [0] * nb_sing
                border_side[sing_to_index[(a,d,0)]] += 1
                border_side[sing_to_index[(a,d,1)]] += -1
                borders.append(border_side)

        from sage.matrix.constructor import matrix
        from sage.rings.integer_ring import ZZ
        return matrix(ZZ, borders)

    def profile(self):
        r"""
        Return the profile of the surface.

        Return the list of angles of singularities in the surface divided by pi.

        EXAMPLES::

            sage: from surface_dynamics.all import *
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
        p_id = SymmetricGroup(self._degree_cover).identity()
        for orbit in base_diagram:
            flat_orbit = []
            for x in orbit:
                if isinstance(x[0], tuple):
                    flat_orbit.extend(x)
                else:
                    flat_orbit.append(x)
            p = p_id
            for lab,sign in flat_orbit:
                q = self.covering_data(lab)
                if sign: q = q.inverse()
                p = p*q
            for c in p.cycle_type():
               s.append(len(orbit)*c)

        return Partition(sorted(s,reverse=True))

    def is_orientable(self):
        r"""
        Test whether this permutation cover has an orientable foliation.

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: from itertools import product
            sage: it = iter(product(('()', '(1,2)'), repeat=3))
            sage: it.next()
            ('()', '()', '()')

            sage: for cov in it:
            ....:     c = p.cover(cov)
            ....:     print cov, c.is_orientable()
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
            signs = [+1] + [None] * (self._degree_cover - 1)
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

    def stratum(self):
        r"""
        Stratum of the covering translation surface

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: p.cover(['(1,2)', '()', '(1,2)']).stratum()
            H_1(0^4)
            sage: p.cover(['(1,2)', '(1,2)', '(1,2)']).stratum()
            Q_0(0^2, -1^4)
        """
        if self.is_orientable():
            from surface_dynamics.flat_surfaces.abelian_strata import AbelianStratum
            return AbelianStratum([(x-2)/2 for x in self.profile()])
        else:
            from surface_dynamics.flat_surfaces.quadratic_strata import QuadraticStratum
            return QuadraticStratum([x-2 for x in self.profile()])

    def genus(self):
        r"""
        Genus of the covering translation surface

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: p.cover(['(1,2)', '()', '(1,2)']).genus()
            1
            sage: p.cover(['(1,2)', '(1,2)', '(1,2)']).genus()
            0

        TESTS::

            sage: from surface_dynamics.all import *
            sage: o = AbelianStratum([1,2,3,4]).one_component().one_origami()
            sage: assert(o.genus() == AbelianStratum([1,2,3,4]).genus())
            sage: qc = QuadraticStratum([1,2,3,4,-1,-1]).one_component()
            sage: p = qc.permutation_representative()
            sage: assert(p.orientation_cover().genus() == qc.orientation_cover_component().genus())
        """
        p = self.profile()
        return Integer((sum(p)-2*len(p))/4+1)

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

            sage: from surface_dynamics.all import *
            sage: p = iet.GeneralizedPermutation('a b b','c c a')
            sage: c = p.cover(['(1,2,3)','(1,3,2)','()'])
            sage: c._real_characters()
            ([(1, 1, 1), (2, -1, -1)], [1, 1])
        """
        from surface_dynamics.misc.group_representation import real_characters
        return real_characters(self.automorphism_group())

    @cached_method
    def _cc_mats(self):
        r"""
        The projection given by the conjugacy class.

        This is cached to speed up the computation of the projection matrices.
        See :meth:`~flatsurf.misc.group_representation.conjugacy_class_matrix`
        for more informations.

        TESTS::

            sage: from surface_dynamics.all import *
            sage: p = iet.GeneralizedPermutation('a b b','c c a')
            sage: c = p.cover(['(1,2,3)','(1,3,2)','()'])
            sage: c._cc_mats()
            (
            [1 0 0]  [0 1 0]  [0 0 1]
            [0 1 0]  [0 0 1]  [1 0 0]
            [0 0 1], [1 0 0], [0 1 0]
            )
        """
        from surface_dynamics.misc.group_representation import conjugacy_class_matrix
        G = self.automorphism_group()
        mats = []
        for cl in libgap(G).ConjugacyClasses():
            m = conjugacy_class_matrix(cl, self._degree_cover)
            m.set_immutable()
            mats.append(m)
        return tuple(mats)

    def isotypic_projection_matrix(self, i):
        r"""
        TESTS::

            sage: from surface_dynamics.all import *
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
        """
        from sage.matrix.special import identity_matrix
        from surface_dynamics.misc.group_representation import isotypic_projection_matrix
        m = isotypic_projection_matrix(
               self.automorphism_group(),
               self._degree_cover,
               self._real_characters()[0][i],
               self._real_characters()[1][i],
               self._cc_mats())
        return m.tensor_product(identity_matrix(len(self._base)), subdivide=False)

    def lyapunov_exponents_H_plus(self, nb_vectors=None, nb_experiments=10,
                                  nb_iterations=65536, output_file=None,
                                  return_speed=False, isotypic_decomposition=False,
                                  return_char=False, verbose=False):
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

         - ``return_speed`` -- wether or not return the lyapunov exponents list
         in a pair with the speed of the geodesic.

         - ``isotypic_decomposition`` -- either a boolean or a character or a list of characters.

         - ``return_char`` -- whether or not return the character corresponding to
         the isotypic component.

         - ``verbose`` -- if ``True`` provide additional informations rather than
         returning only the Lyapunov exponents (i.e. ellapsed time, confidence
         intervals, ...)


        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: q = QuadraticStratum([1,1,-1,-1]).one_component()
            sage: q.lyapunov_exponents_H_plus() # abs tol 0.01
            [0.666666]
            sage: p = q.permutation_representative(reduced=False).orientation_cover()
            sage: c = p.lyapunov_exponents_H_plus(isotypic_decomposition=True)[0]
            sage: print c[0] # abs tol 0.05
            1.000000

            sage: p = iet.GeneralizedPermutation('e a a', 'b b c c d d e')
            sage: p.stratum()
            Q_0(1, -1^5)
            sage: p.alphabet()
            {'e', 'a', 'b', 'c', 'd'}
            sage: c = p.cover(['()', '(1,2)', '()', '(1,2)', '(1,2)'])
            sage: c.stratum()
            Q_1(1^2, 0^4, -1^2)
            sage: c.lyapunov_exponents_H_plus() # abs tol 0.01
            [0.666666]

        Some cyclic covers (see [EKZ]_ for the formulas)::

            sage: p = iet.GeneralizedPermutation('c a a', 'b b c', alphabet='abc')
            sage: def cyclic(n,a):
            ....:     return [(i+a)%n + 1 for i in range(n)]
            sage: def cyclic_cover(n, a, b, c):
            ....:     return p.cover([cyclic(n,c), cyclic(n,a), cyclic(n, b)])

            sage: c = cyclic_cover(7,1,1,2)
            sage: c.lyapunov_exponents_H_plus(isotypic_decomposition=True) # abs tol 0.05
            [[],
             [0.2857, 0.2857],
             [0.5714, 0.5714],
             [0.2857, 0.2857]]

            sage: c = cyclic_cover(7, 1, 2, 3)
            sage: lexp = c.lyapunov_exponents_H_plus(isotypic_decomposition=True, return_char=True)
            sage: lexp[0]
            ([], (1, 1, 1, 1, 1, 1, 1))
            sage: lexp[1][0] # abs tol 0.05
            [0.2857, 0.2857]
            sage: lexp[1][1]
            (2, E(7) + E(7)^6, ..., E(7) + E(7)^6)
            sage: lexp[2][0] # abs tol 0.05
            [0.2857, 0.2857]
            sage: lexp[2][1]
            (2, E(7)^2 + E(7)^5, ...,  E(7)^2 + E(7)^5)
            sage: lexp[3][0] # abs tol 0.05
            [0.5714, 0.5714]
            sage: lexp[3][1]
            (2, E(7)^3 + E(7)^4, ..., E(7)^3 + E(7)^4)

        The Eierlegendewollmilchsau as a quaternionic cover of the once
        punctured torus::

            sage: p = iet.Permutation('a b', 'b a')
            sage: Q = QuaternionGroup()
            sage: a,b = Q.gens()
            sage: c = p.cover([a, b])
            sage: c.lyapunov_exponents_H_plus()  # abs tol 0.05
            [1.0, 0.0, 0.0]

        REFERENCES:

        .. [EKZ] A. Eskin, M. Kontsevich, A. Zorich "Lyapunov spectrum of
                 square-tiled cyclic cover", J. Mod. Dyn. 5 (2011), no. 2, 319-353.
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

        if nb_vectors < 0 :     raise ValueError("the number of vectors must be positive")
        if nb_vectors == 0:     return []
        if nb_experiments <= 0: raise ValueError("the number of experiments must be positive")
        if nb_iterations <= 0 : raise ValueError("the number of iterations must be positive")

        if verbose:
            output_file.write("Stratum: {}\n".format(self.stratum()))

        k = int(len(self[0]))
        gp, twin = range(2*n), range(2*n)
        base_twin = self._base.twin_list()
        rank = self._base.alphabet().rank
        for i in range(2):
            for j,a in enumerate(self._base[i]):
                gp[i*k+j] = rank(a)
                ii,jj = base_twin[i][j]
                twin[i*k+j] = ii*k+jj

        flatten = False
        if isotypic_decomposition:
            if isotypic_decomposition is True:
                characters = self._real_characters()[0]
            elif isinstance(isotypic_decomposition, (tuple,list)):
                if isinstance(isotypic_decomposition[0], (tuple,list)):
                    characters = map(tuple, isotypic_decomposition)
                else:
                    characters = [tuple(isotypic_decomposition)]
                    flatten = True
            else:
                raise ValueError("isotypic_decomposition must be a boolean or a list of characters")

            all_characters = self._real_characters()[0]
            for c in characters:
                if c not in all_characters:
                    raise ValueError("{} is an invalid character".format(c))

            nc = len(characters)
            na = len(all_characters)
            dim = len(self._base) * self._degree_cover
            projections = range(dim**2 * nc)
            dimensions = range(nc)
            ii = 0

            H = self._delta1().left_kernel().matrix()
            B = self._delta2()
            for i_char in xrange(na):
                if all_characters[i_char] not in characters:
                    continue
                M = self.isotypic_projection_matrix(i_char)
                dimension = (H*M).rank() - (B*M).rank()

                if dimension%2:
                    raise RuntimeError("get a subspace of odd dimension\n" +
                            "  i_char = {}\n".format(i_char) +
                            "  dim    = {}".format(2*dimension))
                dimensions[ii] = dimension//2
                for i in xrange(dim):
                    for j in xrange(dim):
                        projections[ii * (dim**2) + i * dim + j] = float(M[j][i])
                ii += 1
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
            output_file.write("ellapsed time %s\n"%time.strftime("%H:%M:%S",time.gmtime(t1-t0)))
            output_file.write("Lexp Rauzy-Zorich: %f (std. dev. = %f, conf. rad. 0.01 = %f)\n"%(
                    m,d, 2.576*d/sqrt(nb_experiments)))

        if isotypic_decomposition:
            i_0 = 1
            for i_char in xrange(nc):
                res_int = []
                if verbose:
                    output_file.write("##### char_%d #####\n"%(i_char))
                    output_file.write("chi = {}\n".format(self._real_characters()[0][i_char]))
                    output_file.write("dim = {}\n".format(dimensions[i_char]))
                for i in xrange(i_0, i_0 + dimensions[i_char]):
                    m,d = mean_and_std_dev(res[i])
                    res_int.append(m)
                    if verbose:
                        output_file.write("theta_%d          : %f (std. dev. = %f, conf. rad. 0.01 = %f)\n"%(
                            i,m,d, 2.576*d/sqrt(nb_experiments)))
                res_final.append((res_int, self._real_characters()[0][i_char]) if return_char else res_int)
                i_0 += dimensions[i_char]

        else:
            for i in xrange(1,nb_vectors+1):
                m,d = mean_and_std_dev(res[i])
                if verbose:
                    output_file.write("theta%d           : %f (std. dev. = %f, conf. rad. 0.01 = %f)\n"%(i,m,d, 2.576*d/sqrt(nb_experiments)))
                res_final.append(m)

        res_final = res_final[0] if flatten else res_final
        if return_speed:
            return (lexp, res_final)
        else:
            return res_final

    @cached_method
    def monodromy(self):
        r"""
        Return the monodromy of this covering.

        That it to say the permutation group generated by the action of the
        fundamental group.

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: p.cover(['(1,2,3)', '(1,3,2)', '']).monodromy()
            Permutation Group with generators [(1,2,3), (1,3,2), ()]
            sage: p.cover(['(1,2)', '(1,3)', '']).monodromy()
            Permutation Group with generators [(1,2), (1,3), ()]
        """
        from sage.groups.perm_gps.permgroup import PermutationGroup
        return PermutationGroup([self.covering_data(a) for a in self._base.letters()], canonicalize=False)

    @cached_method
    def automorphism_group(self):
        r"""
        Return the Deck group of the cover.

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: p.cover(['(1,2,3)', '(1,3,2)', '']).automorphism_group()
            Permutation Group with generators [(1,2,3), (1,3,2)]
            sage: p.cover(['(1,2)', '(1,3)', '']).automorphism_group()
            Permutation Group with generators [()]
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        from sage.groups.perm_gps.permgroup import PermutationGroup

        Sd = SymmetricGroup(self._degree_cover)
        G = libgap.Subgroup(Sd, [self.covering_data(a) for a in self._base.letters()])
        C = libgap.Centralizer(Sd, G)

        return PermutationGroup(libgap.GeneratorsOfGroup(C).sage())

