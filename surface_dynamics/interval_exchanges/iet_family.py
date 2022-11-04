r"""
Linear families of interval exchange transformations

EXAMPLES:

A family is built from a permutation and a cone of admissible lengths::

    sage: from surface_dynamics import iet

    sage: p = iet.Permutation('a b c d e f g h', 'd c f e h g b a')
    sage: C = Polyhedron(rays=[(1, 1, 1, 1, 1, 0, 0, 1), (0, 0, 0, 1, 1, 1, 1, 1)])
    sage: F = iet.IETFamily(p, C)
    sage: F
    Linear iet family of dimension 2 in RR^8
    top a b c d e f g h
    bot d c f e h g b a
    0 0 0 1 1 1 1 1
    1 1 1 1 1 0 0 1

To have an overview of what a typical iet in your family looks like (ie minimal
or not) you can use
:meth:`~surface_dynamics.interval_exchanges.iet_family_pyx.IETFamily_pyx.random_element_statistics`::

    sage: (n1, n2, n3) = F.random_element_statistics(NumberField(x^3 - 2, 'a', embedding=AA(2)**(1/3)))
    sage: (n1, n2, n3)  # random
    (100, 0, 0)
    sage: assert n1 > 90 and n2 < 5 and n3 < 5, (n1, n2, n3)

The output in the above call is a triple of numbers ``(num_minimal_iets,
num_iets_with_saddle_connections, num_iets_with_unknown_behaviour)`` on a
random sample of 100 elements. Here, the family ``F`` seems to be mostly
generically made of minimal iet, however Boshernitzan criterion is not able to
certify it directly::

    sage: F.is_boshernitzan()
    False

The method :meth:`~IETFamily.rauzy_induction` applies Rauzy induction to get
rid of candidate saddle connections. For the family ``F`` it succeeds after a
single iteration (which is a top Rauzy induction)::

    sage: list(F.rauzy_induction(2))
    [(Linear iet family of dimension 2 in RR^8
      top a b c d e f g h
      bot d c f e h a g b
      0 0 0 1 1 1 1 1
      1 1 1 1 1 0 0 0,
      'boshernitzan',
      't')]

We now illustrate a more advanced feature. Using the iterator
:func:`~surface_dynamics.misc.linalg.isotropic_subspaces`. One can explore the
linear subspaces with vanishing Sah-Arnoux-Fathi invariant as follows::

    sage: from surface_dynamics.misc.linalg import isotropic_subspaces

    sage: p = iet.Permutation([1,2,0,3], [2,1,3,0], alphabet=[0,1,2,3])
    sage: x = polygen(QQ)
    sage: K.<cbrt3> = NumberField(x^3 - 3, embedding=AA(3)**(1/3))
    sage: for vectors in isotropic_subspaces(p.intersection_matrix(), 2, bound=1, contains_positive_vector=True):
    ....:     P = Polyhedron(lines=vectors).intersection(Polyhedron(rays=(ZZ**4).basis()))
    ....:     assert P.dimension() == 2
    ....:     F = iet.IETFamily(p, P)
    ....:     T = F.random_element(K)
    ....:     assert T.sah_arnoux_fathi_invariant().is_zero()

For each member of such family, one can look for iet with non-trivial dynamics (here none)::

    sage: x = polygen(QQ)
    sage: K.<cbrt3> = NumberField(x^3 - 3, embedding=AA(3)**(1/3))
    sage: p = iet.Permutation('a b c d e f', 'f e d c b a')
    sage: for vectors in isotropic_subspaces(p.intersection_matrix(), 3, bound=1, contains_positive_vector=True):
    ....:     P = Polyhedron(lines=vectors).intersection(Polyhedron(rays=(ZZ**6).basis()))
    ....:     assert P.dimension() == 3
    ....:     F = iet.IETFamily(p, P)
    ....:     T = F.random_element(K)
    ....:     assert T.sah_arnoux_fathi_invariant().is_zero()
    ....:     n_minimals, n_saddles, n_unknowns = F.random_element_statistics(K, num_exp=10, num_iterations=4096)
    ....:     if n_saddles != 10:
    ....:         print(vectors, n_minimals, n_saddles)
"""
#*****************************************************************************
#       Copyright (C) 2022 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************



from sage.misc.cachefunc import cached_method
from .iet_family_pyx import IETFamily_pyx

class IETFamily(IETFamily_pyx):
    def __init__(self, p, C):
        IETFamily_pyx.__init__(self, p, C)
        self._alphabet = p.alphabet()

    def __repr__(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: F = iet.IETFamily(p, [(2,3,0,0), (0,1,1,1)])
            sage: repr(F)  # indirect doctest
            'Linear iet family of dimension 2 in RR^4\ntop a b c d\nbot d c b a\n0 1 1 1\n2 3 0 0'
        """
        s = []

        s.append('Linear iet family of dimension {} in RR^{}'.format(self.dimension(), self.ambient_dimension()))
        perm = self.permutation()
        s.append('top ' + ' '.join(map(str, perm[0])))
        s.append('bot ' + ' '.join(map(str, perm[1])))
        for r in self.rays():
            s.append(' '.join(map(str,r)))
        return '\n'.join(s)

    def _new(self):
        F = IETFamily.__new__(IETFamily)
        F._alphabet = self._alphabet
        return F

    def permutation(self):
        r"""
        Return the permutation of this family.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: F = iet.IETFamily(p, Polyhedron(rays=(ZZ**4).basis()))
            sage: F.permutation()
            a b c d
            d c b a
        """
        from .constructors import Permutation

        top, bot = IETFamily_pyx.permutation(self)
        p = Permutation(top, bot, alphabet=range(self.ambient_dimension()))
        p.alphabet(self._alphabet)
        return p

    def rays(self):
        r"""
        Return the rays as vectors

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: F = iet.IETFamily(p, [(2,3,0,0), (0,1,1,1)])
            sage: F.rays()
            [(0, 1, 1, 1), (2, 3, 0, 0)]
        """
        F = self.free_module()
        return [F([self.ray_coefficient(i, j) for j in range(self.ambient_dimension())]) for i in range(self.n_rays())]

    def has_zero_saf(self):
        r"""
        Return whether this family only consists of interval exchanges with vanishing SAF invariant.

        EXAMPLES::

            sage: from surface_dynamics import iet

            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: F = iet.IETFamily(p, [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
            sage: F.has_zero_saf()
            False

            sage: path = 'bbbbbtbttbttbbttbtttbbbbttbttbtbbtbbtbbtt'
            sage: p = iet.Permutation('a b c d e f g h i', 'i h g f e d c b a')
            sage: R = p.rauzy_diagram()
            sage: T = R.path(p, *path).self_similar_iet()
            sage: F = iet.IETFamily(T)
            sage: F.has_zero_saf()
            True
        """
        p = self.permutation()
        F = self.free_module()
        rays = [F(r) for r in self.rays()]
        omega = p.intersection_matrix()
        return all((rays[i] * omega * rays[j]).is_zero() for i in range(len(rays)) for j in range(i))

    def is_periodic_boshernitzan(self, solver="PPL"):
        r"""
        Return whether Boshernitzan criterion holds for periodic trajectories.

        The output is a boolean. If ``True`` then the family is generically
        without periodic trajectories. If ``False`` then there are candidate
        homology classes that could be stable periodic trajectories in the
        family. Let us emphasize that a ``False`` output does not imply the
        presence of (stable) periodic trajectories.

        .. SEEALSO::

            :meth:`IETFamily.is_boshernitzan` tests the absence of saddle connections
            (aka Keane property i.d.o.c. condition) in the family.
        """
        from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException

        p = self.permutation()

        F = self.free_module()
        E = F.subspace(self.rays()).basis_matrix()
        R = E.right_kernel_matrix()
        omega = p.intersection_matrix()

        M = MixedIntegerLinearProgram(solver=solver)
        # a1 tau_1 + ... a_k tau_k = b_1 r_1 + ... b_s _s
        a = M.new_variable(nonnegative=True)
        b = M.new_variable()

        M.add_constraint(M.sum(a[j] for j in range(len(p))) >= 1)
        for i in range(len(p)):
            t_sum = M.sum(a[j] * omega[j,i] for j in range(len(p)))
            r_sum = M.sum(b[j] * R[j,i] for j in range(R.nrows()))
            M.add_constraint(t_sum == r_sum)
        try:
            M.solve(objective_only=True)
            return False
        except MIPSolverException:
            pass

        return True

    def _intersection_polytope(self, v):
        from sage.geometry.polyhedron.constructor import Polyhedron

        p = self.permutation()

        F = self.free_module()
        E = F.subspace(self.rays()).basis_matrix()
        R = E.right_kernel_matrix()
        omega = p.intersection_matrix()

        # a1 tau_1 + ... a_k tau_k = b_1 r_1 + ... b_s _s
        poscone = Polyhedron(rays=omega.rows())
        linspace = Polyhedron(vertices=[v], lines=R.rows())
        return poscone.intersection(linspace)

    def cylinder_candidates(self):
        return self._intersection_polytope(self.free_module().zero())

    def saddle_connection_candidates(self, itop=None, ibot=None):
        r"""
        Iterate through the list of saddle connection candidates.

        The list is presented as a list of triples ``(ibot, itop, polytope)``
        where ``ibot`` and ``itop`` are respectively indices of singularities
        in the bottom (range) and top (domain) intervals and polytope is
        the polytope whose integer points are the candidates for Abelianization
        of trajectories from ``ibot`` to ``itop``.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: p = iet.Permutation('a e c d b', 'e d c b a')
            sage: pos = Polyhedron(rays=(QQ**5).basis())
            sage: C = Polyhedron(eqns=[(0, 1, 0, 0, 0, -1)]).intersection(pos)
            sage: F = iet.IETFamily(p, C)
            sage: list(F.saddle_connection_candidates())
            [(4,
              4,
              A 0-dimensional polyhedron in QQ^5 defined as the convex hull of 1 vertex)]
        """
        F = self.free_module()
        p = self.permutation()

        if itop is None or ibot is None:
            if itop is not None or ibot is not None:
                raise ValueError('itop and ibot must be simultaneously zero')
            itops = range(1, len(p))
            ibots = range(1, len(p))
        else:
            itops = [itop]
            ibots = [ibot]

        for ibot in ibots:
            bot = sum(F.gen(p._labels[1][i]) for i in range(ibot))
            for itop in itops:
                top = sum(F.gen(p._labels[0][i]) for i in range(itop))
                P = self._intersection_polytope(top - bot)
                if not P.is_empty():
                    yield (ibot, itop, P)

    def is_boshernitzan(self, itop=None, ibot=None, certificate=False, solver="PPL"):
        r"""
        Return whether this slice satisfies Boshernitzan conditions.

        The output is a boolean. If ``True`` the the family is generically without
        saddle connections. If ``False`` then there are candidate relative homology
        classes that could be stable saddle connections in the family. Let us
        emphasize that a ``False`` output does not imply the presence of (stable)
        saddle connections.

        .. SEEALSO::

            :meth:`~surface_dynamics.interval_exchanges.iet_family_pyx.IETFamily_pyx.has_zero_connection`
            test whether there is a connection of zero length (ie a bottom
            singularity that coincides with a top one)

            :meth:`IETFamily.is_periodic_boshernitzan` test for absence of periodic
            trajectories (which is a weaker than absence of saddle connections)

        EXAMPLES::

            sage: from surface_dynamics import iet

        Two examples with connections::

            sage: p = iet.Permutation('a e c d b', 'e d c b a')
            sage: pos = Polyhedron(rays=(QQ**5).basis())
            sage: C = Polyhedron(eqns=[(0, 1, 0, 0, 0, -1)]).intersection(pos)
            sage: F = iet.IETFamily(p, C)
            sage: F.is_boshernitzan()
            False
            sage: F.has_zero_connection()
            True

            sage: C = Polyhedron(eqns=[(0, 1, -1, 0, 0, 0)]).intersection(pos)
            sage: F = iet.IETFamily(p, C)
            sage: F.is_boshernitzan()
            False
            sage: F.has_zero_connection()
            True

        An example that satisfies Boshernitzan condition::

            sage: C = Polyhedron(eqns=[(0, 1, -2, 0, 0, 0)]).intersection(pos)
            sage: F = iet.IETFamily(p, C)
            sage: F.is_boshernitzan()
            True
            sage: F.has_zero_connection()
            False

        TESTS::

            sage: from surface_dynamics import iet
            sage: p = iet.Permutation('a b c', 'c b a')
            sage: for C in [[[1, 0, 2], [0, 1, 1]],
            ....:           [[1, 0, 3], [0, 1, 1]],
            ....:           [[2, 0, 1], [0, 1, 0]]]:
            ....:     F = iet.IETFamily(p, C)
            ....:     assert not F.has_zero_connection() and not F.is_boshernitzan()
            sage: for C in [[[3, 0, 1], [0, 1, 0]],
            ....:           [[3, 1, 0], [0, 0, 1]],
            ....:           [[3, 0, 1], [0, 1, 1]],
            ....:           [[2, 0, 1], [0, 1, 2]]]:
            ....:     F = iet.IETFamily(p, C)
            ....:     assert not F.has_zero_connection() and F.is_boshernitzan()
            sage: for C in [[[1, 0, 1], [0, 1, 0]],
            ....:           [[1, 0, 1], [0, 1, 1]]]:
            ....:     F = iet.IETFamily(p, C)
            ....:     assert F.has_zero_connection() and not F.is_boshernitzan()

            sage: p = iet.Permutation('a b c d', 'c a d b')
            sage: for C in [[[0, 1, 2, 1], [1, 0, 1, 2]],
            ....:           [[0, 1, 2, 1], [1, 1, 0, 2]],
            ....:           [[0, 1, 2, 2], [2, 0, 1, 1]],
            ....:           [[0, 2, 1, 1], [2, 0, 2, 1]]]:
            ....:     F = iet.IETFamily(p, C)
            ....:     assert not F.has_zero_connection() and not F.is_boshernitzan()
            sage: for C in [[[0, 1, 1, 1], [1, 0, 1, 1]],
            ....:           [[0, 1, 1, 1], [1, 1, 0, 1]],
            ....:           [[0, 1, 1, 1], [1, 1, 1, 0]],
            ....:           [[1, 0, 1, 1], [1, 1, 1, 0]]]:
            ....:     F = iet.IETFamily(p, C)
            ....:     assert F.has_zero_connection() and not F.is_boshernitzan()
            sage: for C in [[[1, 0, 1, 1], [1, 1, 0, 1]]]:
            ....:     F = iet.IETFamily(p, C)
            ....:     assert not F.has_zero_connection() and F.is_boshernitzan()
        """
        from sage.modules.free_module import FreeModule
        from sage.rings.rational_field import QQ
        from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException

        p = self.permutation()

        if itop is None or ibot is None:
            if itop is not None or ibot is not None:
                raise ValueError('itop and ibot must be simultaneously None')
            itops = range(1, len(p))
            ibots = range(1, len(p))
        else:
            itops = [itop]
            ibots = [ibot]

        F = FreeModule(QQ, len(p))
        E = F.subspace(self.rays()).basis_matrix()
        R = E.right_kernel_matrix()
        assert all(r.dot_product(e).is_zero() for e in E.rows() for r in R.rows())
        omega = p.intersection_matrix()

        for itop in itops:
            for ibot in ibots:
                top = sum(F.gen(p._labels[0][i]) for i in range(itop))
                bot = sum(F.gen(p._labels[1][i]) for i in range(ibot))

                M = MixedIntegerLinearProgram(solver=solver)
                # a1 tau_1 + ... a_k tau_k = phi + b_1 r_1 + ... b_s _s
                a = M.new_variable(nonnegative=True)
                b = M.new_variable()
                for i in range(len(p)):
                    t_sum = M.sum(a[j] * omega[j,i] for j in range(len(p)))
                    r_sum = M.sum(b[j] * R[j,i] for j in range(R.nrows()))
                    M.add_constraint(t_sum == top[i] - bot[i] + r_sum)
                try:
                    M.solve(objective_only=True)
                    return False
                except MIPSolverException:
                    pass

        return True

    def rauzy_induction(self, max_depth=5, verbose=False):
        r"""
        Iterate through the new families obtained by performing Rauzy
        induction.

        Performing Rauzy induction cut the simplex of lengths into two
        subsimplices. Performing Rauzy induction allows to zoom in different
        part of the families and make the search for dynamical behaviour more
        accurate (eg saddle connections, periodic components, minimalit
        components). Each step of the Rauzy induction might split the family
        into two familes giving rise to a tree of possibilities. The search is
        cut when

        - either reaching a decidable state (ie family satisfying Boshernitzan
          condition or with the presence of a zero connection).
        - when the depth reaches ``max_depth`` (default to ``5``)

        OUTPUT: Each element of the output is a triple ``(family, state,
        induction_path)`` where

        - ``family`` is the family obtained after performing Rauzy induction
        - ``state`` is a string indicating what kind of family was obtained. It
          is either ``"unknown"``, ``"autosim"`` (when ``family`` was already
          encountered somewhere else in the search), ``"saddle"`` (when all
          elements of ``family`` share has a zero connection) or
          ``"boshernitzan"`` (when ``family`` satisfies Boshernitzan criterion
          (in particular, the iets in ``family`` generically have no saddle
          connections)).
        - ``induction_path`` a string made of letters ``"t"`` and ``"b"`` which
          is the sequence of top and bottom Rauzy induction done to arrive to
          ``family`` from this family

        EXAMPLES::

            sage: from surface_dynamics import iet

        A zero SAF example (which is very likely to mostly have stable
        connections)::

            sage: p = iet.Permutation([0,1,2,3,4,5],[5,4,3,2,1,0])
            sage: rays = [[5, 1, 0, 0, 3, 8], [2, 1, 0, 3, 0, 5], [1, 0, 1, 2, 0, 3], [3, 0, 1, 0, 2, 5]]
            sage: F = iet.IETFamily(p, rays)
            sage: for family, state, path in F.rauzy_induction(10):
            ....:     print(state, path)
            unknown tttbbttbb
            unknown tttbtbttb
            unknown tttbttbtt
            unknown tttbtttbt
            unknown tttbttttb
            unknown tttbttttt
            unknown ttttbbbbb
            unknown ttttbbbbt
            unknown ttttbbbtb
            saddle ttttbbtbb
            saddle ttttbtbb

        An example where part of the family has a stable connection and part of
        it satisfies Boshernitzan condition::

            sage: p = iet.Permutation("a b c d e f g h", "d c f e h g b a")
            sage: rays = [[1, 2, 0, 1, 1, 0, 0, 2],
            ....:         [2, 0, 2, 0, 1, 1, 1, 1],
            ....:         [2, 2, 0, 0, 0, 2, 2, 1]]
            sage: F = iet.IETFamily(p, rays)
            sage: [state for family, state, path in F.rauzy_induction(3)]
            ['boshernitzan', 'saddle']
            sage: x = polygen(QQ)
            sage: K = NumberField(x^3 - 2, 'a', embedding=AA(2)**(1/3))
            sage: n1, n2, n3 = F.random_element_statistics(K)
            sage: (n1, n2, n3)  # random
            (19, 81, 0)
            sage: assert n1 > 0 and n2 > 0 and n3 == 0, (n1, n2, n3)
        """
        saf_zero = self.has_zero_saf()
        s = ''
        f = self
        branch = [[(s, f)]]
        seen = set([f])
        while True:
            if verbose:
                print("branch:")
                for ss,ff in branch[-1]:
                    print(ss)
                    print(ff)
                    print()
            while len(branch) < max_depth:
                branch.append([])
                for ss, ff in f.children():
                    if verbose:
                        print("looking children", ss)
                        print("cone", ff)
                    # check for saddles among the children
                    if ff.has_zero_connection():
                        yield ff, 'saddle', s+ss
                    elif not saf_zero and ff.is_boshernitzan():
                        yield ff, 'boshernitzan', s+ss
                    # check for auto-simlarity
                    elif ff in seen:
                        yield ff, 'autosim', s+ss
                    else:
                        branch[-1].append((s+ss, ff))

                if not branch[-1]:
                    branch.pop(-1)
                    break
                else:
                    s, f = branch[-1][-1]
                    assert f not in seen
                    seen.add(f)

            if len(branch) == max_depth:
                s, f = branch[-1][-1]
                yield f, 'unknown', s

            # backtrack
            while branch and len(branch[-1]) == 1:
                s, f = branch.pop()[0]
                if f not in seen:
                    raise RuntimeError("s = {}\nf =\n{}".format(s, f))
                seen.remove(f)

            if not branch:
                return

            s, f = branch[-1].pop(-1)
            seen.remove(f)
            s, f = branch[-1][-1]
            seen.add(f)
