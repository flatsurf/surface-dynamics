r"""
Strata of differentials on Riemann surfaces

The space of Abelian differentials on Riemann surfaces of a given genus is
stratified by degrees of zeros. Each stratum has one, two or three connected
components each of which is associated to an extended Rauzy class. The
:meth:`~surface_dynamics.flat_surfaces.strata.Stratum.components` method
lists connected components of a stratum.

The work for Abelian differentials was done by Maxim Kontsevich and Anton
Zorich in [KonZor03]_ and for quadratic differentials by Erwan Lanneau in
[Lan08]_. Zorich gave an algorithm to pass from a connected component of a
stratum to the associated Rauzy class (for both interval exchange
transformations and linear involutions) in [Zor08]_ and is implemented for
Abelian stratum at different level (approximately one for each component):

- for connected stratum :meth:`~AbelianStratumComponent.permutation_representative`

- for hyperellitic component
  :meth:`~HypAbelianStratumComponent.permutation_representative`

- for non hyperelliptic component, the algorithm is the same as for connected
  component

- for odd component :meth:`~OddAbelianStratumComponent.permutation_representative`

- for even component :meth:`~EvenAbelianStratumComponent.permutation_representative`

The inverse operation (pass from an interval exchange transformation to
the connected component) is partially written in [KonZor03]_ and
simply named here
?

Some of the code here was first available on Mathematica [ZS]_.

A refinement of Zorich representatives was worked out by L. Jefreys in
[Jef19]_. Namely, for each connected component of Abelian differential
his construction provides a square-tiled surface with both in horizontal
and vertical direction a decomposition with single cylinder of height one.
The implementation is available as

- for connected stratum :meth:`~AbelianStratumComponent.single_cylinder_representative`

- for hyperelliptic component
  :meth:`~HypAbelianStratumComponent.single_cylinder_representative`

- for odd component :meth:`~OddAbelianStratumComponent.single_cylinder_representative`

- for even component :meth:`~EvenAbelianStratumComponent.single_cylinder_representative`

AUTHORS:

- Vincent Delecroix (2009-09-29): initial version

EXAMPLES:

    sage: from surface_dynamics import *

Construction of a stratum from a list of singularity degrees::

    sage: a = Stratum([1,1], k=1)
    sage: a
    H_2(1^2)
    sage: a.surface_genus()
    2
    sage: a.dimension()
    5

::

    sage: a = Stratum([4,3,2,1], k=1)
    sage: a
    H_6(4, 3, 2, 1)
    sage: a.surface_genus()
    6
    sage: a.dimension()
    15

By convention, the degrees are always written in decreasing order::

    sage: a1 = Stratum([4,3,2,1], k=1)
    sage: a1
    H_6(4, 3, 2, 1)
    sage: a2 = Stratum([2,3,1,4], k=1)
    sage: a2
    H_6(4, 3, 2, 1)
    sage: a1 == a2
    True

It is possible to lis strata and their connected components::

    sage: Stratum([10], k=1).components()
    (H_6(10)^hyp, H_6(10)^odd, H_6(10)^even)

Get a list of strata with constraints on genus or on the number of intervals
of a representative::

    sage: AbelianStrata(genus=3).list()
    [H_3(4), H_3(3, 1), H_3(2^2), H_3(2, 1^2), H_3(1^4)]

Obtains the connected components of a stratum::

    sage: a = Stratum([0], k=1)
    sage: a.components()
    (H_1(0)^hyp,)

::

    sage: @cached_function
    ....: def nb_irred_perm(n):
    ....:     if n == 0 or n == 1: return 1
    ....:     return factorial(n) - sum(nb_irred_perm(k) * factorial(n - k) for k in range(1,n))
    sage: [nb_irred_perm(i) for i in range(10)]
    [1, 1, 1, 3, 13, 71, 461, 3447, 29093, 273343]

::

    sage: A = AbelianStrata(dimension=5, fake_zeros=True)
    sage: N = 0
    sage: for a in A:
    ....:    for cc in a.components():
    ....:       for z in set(a.signature()):
    ....:           p = cc.permutation_representative(left_degree=z)
    ....:           n = p.rauzy_diagram().cardinality()
    ....:           print("%13s, %d  :  %d"%(cc, z, n))
    ....:           print(p)
    ....:           N += n
    H_2(2, 0)^hyp, 0  :  11
    0 1 2 3 4
    4 2 1 3 0
    H_2(2, 0)^hyp, 2  :  35
    0 1 2 3 4
    4 1 3 2 0
     H_2(1^2)^hyp, 1  :  15
    0 1 2 3 4
    4 3 2 1 0
     H_1(0^4)^hyp, 0  :  10
    0 1 2 3 4
    4 0 1 2 3
    sage: N
    71
    sage: nb_irred_perm(5)
    71

"""
#*****************************************************************************
#       Copyright (C) 2009-2019 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function, absolute_import, division
from six.moves import range, map, filter, zip
from six import iteritems

import numbers

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets

from sage.combinat.partition import Partitions, Partition
from sage.rings.integer import Integer
from sage.rings.rational import Rational
from sage.rings.infinity import Infinity

from surface_dynamics.flat_surfaces.strata import Stratum, StratumComponent, Strata

def _cylinder_diagrams_with_symmetric(iterator):
    r"""
    Iterate through all possible reflections of the cylinder diagrams in ``iterator``.

    TESTS::

        sage: from surface_dynamics import CylinderDiagram
        sage: from surface_dynamics.flat_surfaces.abelian_strata import _cylinder_diagrams_with_symmetric

        sage: cd_000 = CylinderDiagram("(0,1)-(0,3,4) (2,3)-(1) (4)-(2)")
        sage: cd_100 = CylinderDiagram("(0,2,3,1)-(0,2,1,4) (4)-(3)")
        sage: cd_010 = CylinderDiagram("(0,1)-(0,2,4,5) (2,5)-(3) (3,4)-(1)")
        sage: cd_001 = CylinderDiagram("(0,3)-(0,4,5) (1,4,2)-(1,3) (5)-(2)")
        sage: cd_111 = CylinderDiagram("(0,3,1,2)-(0,4,1,5) (4)-(2) (5)-(3)")

        sage: sum(1 for _ in _cylinder_diagrams_with_symmetric([cd_000]))
        4
        sage: sum(1 for _ in _cylinder_diagrams_with_symmetric([cd_100]))
        2
        sage: sum(1 for _ in _cylinder_diagrams_with_symmetric([cd_010]))
        2
        sage: sum(1 for _ in _cylinder_diagrams_with_symmetric([cd_001]))
        2
        sage: sum(1 for _ in _cylinder_diagrams_with_symmetric([cd_111]))
        1
    """
    for cd in iterator:
        yield cd
        sym = cd.symmetries()
        if sym == (True, True, True):
            pass
        elif sym == (False, False, False):
            yield cd.horizontal_symmetry()
            yield cd.vertical_symmetry()
            yield cd.inverse()
        elif sym == (True, False, False):
            yield cd.vertical_symmetry()
        elif sym == (False, True, False):
            yield cd.horizontal_symmetry()
        elif sym == (False, False, True):
            yield cd.horizontal_symmetry()
        else:
            raise RuntimeError


def DeprecatedAbelianStratumConstructor(*l):
    """
    TESTS::

        sage: from surface_dynamics import AbelianStratum, Stratum

        sage: s = AbelianStratum(0)
        doctest:warning
        ...
        UserWarning: AbelianStratum has changed its arguments in order to handle meromorphic and higher-order differentials; use Stratum instead
        sage: s is loads(dumps(s))
        True

        sage: AbelianStratum(1,1,1,1) is AbelianStratum(1,1,1,1)
        True
        sage: AbelianStratum(1,1,1,1) is Stratum((1, 1, 1, 1), 1)
        True
    """
    import warnings

    warnings.warn('AbelianStratum has changed its arguments in order to handle meromorphic and higher-order differentials; use Stratum instead')

    if len(l) == 1:
        try:
            Integer(l[0])
        except TypeError:
            l = l[0]
    elif len(l) == 2 and isinstance(l[0], (tuple,list)) and isinstance(1, numbers.Integral):
        l = tuple(l[0]) + (0,) * l[1]
    if isinstance(l, dict):
        l = sum(([i] * mult for i, mult in l.items()), [])

    l = tuple(sorted(map(Integer, l), reverse=True))
    return Stratum(l, 1)


class AbelianStratum(Stratum):
    """
    Stratum of Abelian differentials.

    A stratum with a marked outgoing separatrix corresponds to Rauzy diagram
    with left induction, a stratum with marked incoming separatrix correspond
    to Rauzy diagram with right induction.
    If there is no marked separatrix, the associated Rauzy diagram is the
    extended Rauzy diagram (consideration of the
    :meth:`surface_dynamics.interval_exchanges.template.Permutation.symmetric`
    operation of Boissy-Lanneau).

    When you want to specify a marked separatrix, the degree on which it is is
    the first term of your degrees list.

    INPUT:

    - ``marked_separatrix`` - ``None`` (default) or 'in' (for incoming
      separatrix) or 'out' (for outgoing separatrix).

    EXAMPLES::

        sage: from surface_dynamics import *

    Creation of an Abelian stratum and get its connected components::

        sage: a = Stratum((2, 2), k=1)
        sage: a
        H_3(2^2)
        sage: a.components()
        (H_3(2^2)^hyp, H_3(2^2)^odd)

    Get a permutation representative of a connected component::

        sage: a = Stratum((2,2), k=1)
        sage: a_hyp, a_odd = a.components()
        sage: a_hyp.permutation_representative()
        0 1 2 3 4 5 6
        6 5 4 3 2 1 0
        sage: a_odd.permutation_representative()
        0 1 2 3 4 5 6
        3 2 4 6 5 1 0

    You can specify the alphabet::

        sage: a_odd.permutation_representative(alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        A B C D E F G
        D C E G F B A
    """
    _name = 'H'
    _latex_name = '\\mathcal{H}'

    def has_odd_component(self):
        r"""
        Test whether this stratum has an odd spin component.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: Stratum((2,), k=1).has_odd_component()
            False
            sage: Stratum((4,), k=1).has_odd_component()
            True
            sage: Stratum((4,), k=1).odd_component()
            H_3(4)^odd

            sage: Stratum((0,), k=1).has_odd_component()
            False
        """
        if any(z < 0 for z in self.signature()):
            raise NotImplementedError('meromorphic stratum')
        return all(z % 2 == 0 for z in self.signature()) and self.surface_genus() >= 3

    def has_even_component(self):
        r"""
        Test whether this stratum has an even spin component.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: Stratum((2,2), k=1).has_even_component()
            False
            sage: Stratum((6,), k=1).has_even_component()
            True
            sage: Stratum((6,), k=1).even_component()
            H_4(6)^even

            sage: Stratum((0,), k=1).has_even_component()
            False
        """
        if any(z < 0 for z in self.signature()):
            raise NotImplementedError('meromorphic stratum')
        return all(z % 2 == 0 for z in self.signature()) and self.surface_genus() >= 4

    def has_hyperelliptic_component(self):
        r"""
        Test whether this stratum has an hyperelliptic component.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: Stratum((2,1,1), k=1).has_hyperelliptic_component()
            False
            sage: Stratum((2,2), k=1).has_hyperelliptic_component()
            True
            sage: Stratum((0,0,0), k=1).has_hyperelliptic_component()
            True
            sage: Stratum((2,0), k=1).has_hyperelliptic_component()
            True
            sage: Stratum((2,2), k=1).hyperelliptic_component()
            H_3(2^2)^hyp
        """
        if any(z < 0 for z in self.signature()):
            raise NotImplementedError('meromorphic stratum')
        z = [z for z in self.signature() if z > 0]
        return len(z) <= 1 or (len(z) == 2 and z[0] == z[1])

    def has_non_hyperelliptic_component(self):
        r"""
        Test whether this stratum has a non-hyperelliptic component.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: Stratum((1,1), k=1).has_non_hyperelliptic_component()
            False
            sage: Stratum((3,3), k=1).has_non_hyperelliptic_component()
            True
            sage: Stratum((3,3,0), k=1).has_non_hyperelliptic_component()
            True

            sage: Stratum((3,3), k=1).non_hyperelliptic_component()
            H_4(3^2)^nonhyp
        """
        if any(z < 0 for z in self.signature()):
            raise NotImplementedError('meromorphic stratum')
        z = [z for z in self.signature() if z > 0]
        return len(z) == 2 and z[0] == z[1] and z[0] % 2 == 1 and z[0] > 1

    # TODO: redesign this!!!
    def odd_component(self):
        r"""
        Return the odd component of self (if any).

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: a = Stratum([2,2], k=1)
            sage: a
            H_3(2^2)
            sage: a.odd_component()
            H_3(2^2)^odd
        """
        if not self.has_odd_component():
            raise ValueError("No odd spin component in this stratum")
        return OddASC(self)

    def even_component(self):
        r"""
        Return the even component of self (if any)

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: a = Stratum({2:4}, k=1)
            sage: a
            H_5(2^4)
            sage: a.even_component()
            H_5(2^4)^even
        """
        if not self.has_even_component():
            raise ValueError("No even spin component in this stratum")
        return EvenASC(self)

    def hyperelliptic_component(self):
        r"""
        Return the hyperelliptic component of self (if any)

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: a = Stratum([10], k=1)
            sage: a
            H_6(10)
            sage: a.hyperelliptic_component()
            H_6(10)^hyp
        """
        if not self.has_hyperelliptic_component():
            raise ValueError("No hyperelliptic component in this stratum")
        return HypASC(self)

    def non_hyperelliptic_component(self):
        r"""
        Return the non hyperelliptic component of self (if any)

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: a = Stratum((3,3), k=1)
            sage: a
            H_4(3^2)
            sage: a.non_hyperelliptic_component()
            H_4(3^2)^nonhyp
        """
        if not self.has_non_hyperelliptic_component():
            raise ValueError("No non hyperelliptic component in this stratum")
        return NonHypASC(self)

    #
    # Quadratic cover
    #

    # TODO: this should also be done for higher order differentials
    def orientation_quotients(self, fake_zeros=False):
        r"""
        Return the list of quadratic strata such that their orientation cover
        are contained in this stratum.

        If ``fake_zeros`` (default: False) is True we do care about poles which
        becomes a marked zero.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

        The stratum H(2g-2) has one conic singularities of angle `2(2g-1)pi`. The
        only way a surface in H(2g-2) covers a quadratic differential is that
        the quadratic differential has as unique zeros a conical singularity of
        angle `(2g-1) \pi`. The number of poles may vary and give a collection
        of possibilities::

            sage: Stratum([2], k=1).orientation_quotients()
            [Q_0(1, -1^5)]
            sage: Stratum([4], k=1).orientation_quotients()
            [Q_1(3, -1^3), Q_0(3, -1^7)]
            sage: Stratum([6], k=1).orientation_quotients()
            [Q_2(5, -1), Q_1(5, -1^5), Q_0(5, -1^9)]

        A stratum with two zeros may or may not have orientation quotients::

            sage: Stratum([1,1], k=1).orientation_quotients()
            [Q_1(2, -1^2), Q_0(2, -1^6)]
            sage: Stratum([2,2], k=1).orientation_quotients()
            [Q_1(1^2, -1^2), Q_0(1^2, -1^6), Q_1(4, -1^4), Q_0(4, -1^8)]
            sage: Stratum([3,1], k=1).orientation_quotients()
            []

        To impose that covering of poles are fake zeros, switch option
        ``fake_zeros`` to ``True``::

            sage: Stratum([2,2,0,0], k=1).orientation_quotients(fake_zeros=True)
            [Q_1(1^2, -1^2)]
        """
        e = {}
        z = tuple(m for m in self.signature() if m)
        for i in z:
            if i not in e: e[i] = 0
            e[i] += 1

        # the odd degrees (corresponding to angles 2((2m+1)+1) times pi should
        # be non ramified and hence come by pair.
        if any(e[i]%2 for i in e if i%2):
            return []

        pairings = []
        for d, m in iteritems(e):
            if d % 2: # if the degree is odd it is necessarily non ramified
                pairings.append([(d, m//2)])
            else: # if the degree is even ramified and non ramified are possible
                pairings.append([(d, k) for k in range(m//2 + 1)])

        import itertools
        res = []

        for p in itertools.product(*pairings):
            ee = dict((d - 1, 0) for d in e)
            ee.update((2 * d, 0) for d in e)
            for d,m in p:
                ee[d - 1] += e[d] - 2 * m
                ee[2 * d] += m

            degrees = []
            for d in ee: degrees.extend([d] * ee[d])

            s = sum(degrees)
            self_nb_fake_zeros = self.signature().count(0)
            for nb_poles in range(s % 4, s + 5, 4):
                q = Stratum(degrees + [-1] * nb_poles, k=2)
                q_nb_poles = q.signature().count(-1)
                if not q.is_empty() and (not fake_zeros or q_nb_poles <= self_nb_fake_zeros):
                    res.append(q)

        return res

    #
    # Separatrix and cylinder diagrams
    #

    # TODO: to be removed. Make separatrix_diagrams an iterator.
    def separatrix_diagram_iterator(self, ncyls=None):
        r"""
        Return an iterator over the separatrix diagrams of this stratum.

        For strata of small dimension, it could be faster to use the method
        separatrix_diagrams.

        INPUT:

        - ``ncyls`` -- an optional number of cylinders
        """
        from .separatrix_diagram import separatrix_diagram_iterator
        return separatrix_diagram_iterator([m+1 for m in self.signature()],ncyls)

    def separatrix_diagrams(self, ncyls=None):
        r"""
        Returns the list of separatrix diagrams that appears in this stratum.

        INPUT:

        - ``database`` - boolean (default: True) - if True, use the
          FlatSurfacesDatabase

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: a = Stratum([2], k=1)
            sage: a
            H_2(2)
            sage: for s in a.separatrix_diagrams(): print(s)
            (0,1,2)-(0,1,2)
            (0)(1,2)-(0,1)(2)

        TESTS::

            sage: from surface_dynamics import Stratum

            sage: for (zeros, ncyl) in [((4,), 3), ((2,2), 4)]:
            ....:     S = Stratum(zeros, k=1).separatrix_diagrams(3)
            ....:     for i in range(len(S)):
            ....:         for j in range(i):
            ....:              assert not S[i].is_isomorphic(S[j])
        """
        return sorted(self.separatrix_diagram_iterator(ncyls))

    def cylinder_diagram_iterator(self, ncyls=None, up_to_symmetry=True, force_computation=False):
        r"""
        Iterator over all cylinder diagram of this stratum.

        The generation is up to isomorphism and horizontal/vertical symmetry
        (and they are in standard form).

        INPUT:

        - ``ncyls`` -- an optional number of cylinders

        - ``up_to_symmetry`` - (boolean, default ``True``) to return only
          cylinder diagrams up to horizontal and vertical symmetry.

        - ``force_computation`` -- if ``True`` do no use the database of
          cylinder diagrams (default is ``False``)

        EXAMPLES::

            sage: from surface_dynamics import Stratum, CylinderDiagram

            sage: A = Stratum([4], k=1)
            sage: C1 = [CylinderDiagram('(0,2,1)-(0,3,4) (3)-(2) (4)-(1)'),
            ....:       CylinderDiagram('(0,2,1)-(0,3,4) (3)-(1) (4)-(2)'),
            ....:       CylinderDiagram('(0,1)-(0,3,4) (2,3)-(1) (4)-(2)'),
            ....:       CylinderDiagram('(0,2)-(4) (1,4)-(2,3) (3)-(0,1)'),
            ....:       CylinderDiagram('(0,2)-(0,3) (1,3)-(1,4) (4)-(2)'),
            ....:       CylinderDiagram('(0,1)-(0,3) (2,3)-(1,4) (4)-(2)')]
            sage: C2 = list(A.cylinder_diagram_iterator(3, force_computation=True))
            sage: assert len(C1) == len(C2)
            sage: for (c1, c2) in zip(C1, C2):
            ....:     assert c1.is_isomorphic(c2) or \
            ....:            c1.is_isomorphic(c2.horizontal_symmetry()) or \
            ....:            c1.is_isomorphic(c2.vertical_symmetry()) or \
            ....:            c1.is_isomorphic(c2.inverse())

            sage: sum(1 for _ in A.cylinder_diagram_iterator(3, True, True))
            6
            sage: sum(1 for _ in A.cylinder_diagram_iterator(3, True, False))
            6
            sage: sum(1 for _ in A.cylinder_diagram_iterator(3, False, True))
            9
            sage: sum(1 for _ in A.cylinder_diagram_iterator(3, False, False))
            9
        """
        if ncyls is not None:
            if not isinstance(ncyls, numbers.Integral):
                raise TypeError("ncyls should be None or an integer")
            if ncyls < 0 or ncyls > self.surface_genus() + len(self.signature()) - 1:
                raise ValueError("ncyls is not valid")

        if not force_computation:
            for cc in self.components():
                yield from cc.cylinder_diagram_iterator(ncyls, up_to_symmetry, False)
        else:
            for sd in self.separatrix_diagram_iterator(ncyls):
                iterator = sd.cylinder_diagram_iterator(up_to_symmetry=True)
                if not up_to_symmetry:
                    iterator = _cylinder_diagrams_with_symmetric(iterator)
                yield from iterator

    def cylinder_diagrams(self, ncyls=None, up_to_symmetry=True, force_computation=False):
        r"""
        Return a list of cylinder diagram of this stratum.

        INPUT::

        - ``ncyls`` -- an optional number of cylinders

        - ``up_to_symmetry`` - (boolean, default ``True``) to return only
          cylinder diagrams up to horizontal and vertical symmetry.

        - ``force_computation`` -- If ``True`` then do not use the database of
          cylinder diagrams (default is ``False``).

        EXAMPLES::

            sage: from surface_dynamics import Stratum, CylinderDiagram

            sage: A = Stratum([2,2], k=1)
            sage: C1 = [CylinderDiagram('(0,1)-(0,5) (2)-(4) (3,4)-(1) (5)-(2,3)'),
            ....:       CylinderDiagram('(0,2,1)-(3,4,5) (3)-(1) (4)-(2) (5)-(0)'),
            ....:       CylinderDiagram('(0,2,1)-(3,5,4) (3)-(1) (4)-(2) (5)-(0)'),
            ....:       CylinderDiagram('(0,3)-(5) (1)-(0) (2,5)-(3,4) (4)-(1,2)'),
            ....:       CylinderDiagram('(0,3)-(0,5) (1,2)-(1,4) (4)-(3) (5)-(2)'),
            ....:       CylinderDiagram('(0,5)-(3,4) (1,4)-(2,5) (2)-(0) (3)-(1)'),
            ....:       CylinderDiagram('(0,5)-(3,4) (1,4)-(2,5) (2)-(1) (3)-(0)')]
            sage: C2 = A.cylinder_diagrams(4)
            sage: assert len(C1) == len(C2)
            sage: isoms = []
            sage: for c in A.cylinder_diagrams(4):
            ....:     isom = []
            ....:     for i,cc in enumerate(C1):
            ....:         if c.is_isomorphic(cc) or \
            ....:            c.is_isomorphic(cc.horizontal_symmetry()) or \
            ....:            c.is_isomorphic(cc.vertical_symmetry()) or \
            ....:            c.is_isomorphic(cc.inverse()):
            ....:              isom.append(i)
            ....:     assert len(isom) == 1, isom
            ....:     isoms.extend(isom)
            sage: assert sorted(isoms) == [0, 1, 2, 3, 4, 5, 6]

            sage: len(A.cylinder_diagrams(4, up_to_symmetry=False))
            7
            sage: sum(4 / (1 + sum(cd.symmetries())) for cd in A.cylinder_diagrams(4, up_to_symmetry=True))
            7

        Recovering the multiplicity of the symmetric versions::

            sage: total = 0
            sage: for c in Stratum([2,1,1], k=1).cylinder_diagrams(2):
            ....:     total += 4 // (1 + sum(c.symmetries()))
            sage: total
            61
            sage: len(Stratum([2, 1, 1], k=1).cylinder_diagrams(2, up_to_symmetry=False))
            61

        You obtain the same number directly::

            sage: Stratum([2, 1, 1], k=1).cylinder_diagrams_number(2, up_to_symmetry=False)
            61
        """
        return sorted(self.cylinder_diagram_iterator(ncyls, up_to_symmetry, force_computation))

    def cylinder_diagrams_by_component(self, ncyls=None, up_to_symmetry=True, force_computation=False):
        r"""
        Return a dictionary component -> list of cylinder diagrams.

        INPUT:

        - ``ncyls`` - None or integer (default: None) - the number of cylinders

        - ``up_to_symmetry`` - (boolean, default ``True``) to return only
          cylinder diagrams up to horizontal and vertical symmetry.

        - ``force_computation`` - boolean (default: ``False``) - if ``False``,
          then try to use the database.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: A = Stratum([4], k=1)
            sage: cyls = A.cylinder_diagrams_by_component(ncyls=2, force_computation=True)
            sage: A_hyp = A.hyperelliptic_component()
            sage: A_odd = A.odd_component()
            sage: len(cyls[A_odd])
            4
            sage: len(cyls[A_hyp])
            2

            sage: all(c.ncyls() == 2 for c in cyls[A_hyp])
            True
            sage: all(c.stratum_component() == A_hyp for c in cyls[A_hyp])
            True

            sage: all(c.ncyls() == 2 for c in cyls[A_odd])
            True
            sage: all(c.stratum_component() == A_odd for c in cyls[A_odd])
            True

            sage: for ncyls in range(1, 4):
            ....:     for up_to_symmetry in [True, False]:
            ....:         cd1 = A.cylinder_diagrams_by_component(ncyls, up_to_symmetry, True)
            ....:         cd2 = A.cylinder_diagrams_by_component(ncyls, up_to_symmetry, False)
            ....:         assert len(cd1[A_hyp]) == len(cd2[A_hyp])
            ....:         assert len(cd1[A_odd]) == len(cd2[A_odd])
        """
        if ncyls is not None:
            if not isinstance(ncyls, (int,Integer)):
                raise TypeError("ncyls should be None or an integer")
            if ncyls < 0 or ncyls > self.surface_genus() + len(self.signature()) - 1:
                raise ValueError

        if not force_computation:
            return dict((c,c.cylinder_diagrams(ncyls, up_to_symmetry, False)) for c in self.components())

        d = dict((c,[]) for c in self.components())
        for cyl in self.cylinder_diagrams(ncyls, True, True):
            sc = cyl.stratum_component()
            if not up_to_symmetry:
                d[sc].extend(_cylinder_diagrams_with_symmetric([cyl]))
            else:
                d[sc].append(cyl)

        return d

    def one_cylinder_diagram(self):
        r"""
        Return a diagram with one cylinder in this connected component.

        The diagram returned is the one deduced from the method representative.

        INPUT:

        - ``ncyls`` - the number of cylinders

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: a = Stratum([3,2,1], k=1)
            sage: a
            H_4(3, 2, 1)
            sage: c = a.one_cylinder_diagram();c
            (0,8,3,2,1,6,5,4,7)-(0,8,7,6,5,4,3,2,1)
            sage: c.stratum()
            H_4(3, 2, 1)
        """
        return self.one_component().one_cylinder_diagram()

    def separatrix_diagrams_number(self, ncyls=None):
        r"""
        Return the number of separatrix diagram that belongs to this stratum.
        """
        return sum(1 for _ in self.separatrix_diagram_iterator(ncyls))

    def cylinder_diagrams_number(self, ncyls=None, up_to_symmetry=True, force_computation=False):
        r"""
        Return the number of cylinder diagram that belongs to this stratum.

        INPUT:

        - ``ncyls`` -- an optional number of cylinders

        - ``up_to_symmetry`` - (boolean, default ``True``) to return only
          cylinder diagrams up to horizontal and vertical symmetry.

        - ``force_computation`` -- if ``True`` do no use the database of
          cylinder diagrams (default is ``False``)

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: H22 = Stratum([2,2], k=1)
            sage: H22.cylinder_diagrams_number(3)
            18
            sage: H22.cylinder_diagrams_number(4)
            7

        If ``force_computation`` is set to ``True`` then the database is not
        used. It might be slower for large strata::

            sage: H22.cylinder_diagrams_number(3, force_computation=True)
            18
            sage: H22.cylinder_diagrams_number(4, force_computation=True)
            7

            sage: H31 = Stratum([3,1], k=1)
            sage: for d in range(1,5):
            ....:     print("%d %d" %(H31.cylinder_diagrams_number(d, True, False),
            ....:                     H31.cylinder_diagrams_number(d, True, True)))
            2 2
            12 12
            16 16
            4 4

            sage: H211 = Stratum([2,1,1], k=1)
            sage: for d in range(1,6):
            ....:     print("%d %d" % (H211.cylinder_diagrams_number(d, True, False),
            ....:               H211.cylinder_diagrams_number(d, True, True)))
            5 5
            29 29
            53 53
            27 27
            8 8
        """
        if ncyls is not None and ncyls > self.surface_genus() + len(self.signature()) - 1:
            return 0
        if not force_computation:
            from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
            CDB = CylinderDiagrams()
            if all(CDB.has_component(cc) for cc in self.components()):
                if up_to_symmetry:
                    return CDB.count(self, ncyls)

            return sum(cc.cylinder_diagrams_number(ncyls, up_to_symmetry, False) for cc in self.components())

        if up_to_symmetry:
            return sum(1 for _ in self.cylinder_diagram_iterator(ncyls, True, True))
        else:
            return sum(4 // (1 + sum(cd.symmetries())) for cd in self.cylinder_diagram_iterator(ncyls, True, True))

    def single_cylinder_representative(self, alphabet=None, reduced=True):
        r"""
        Returns a single cylinder permutation representative.

        Returns a permutation representative of a square-tiled surface in this
        component having a single vertical cylinder and a single horizontal cylinder.

        Such representatives were constructed for every stratum of Abelian
        differentials by Jeffreys [Jef19]_.

        INPUT:

        - ``alphabet`` -- an optional alphabet for the permutation representative

        - ``reduced`` (boolean, default ``True``) -- whether to return a reduced
          permutation (ie without labels)

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: C = Stratum([2,0], k=1)
            sage: p = C.single_cylinder_representative()
            sage: p
            0 1 2 3 4
            4 3 1 2 0
            sage: p.stratum() == C
            True

            sage: C = Stratum([3,1], k=1)
            sage: p = C.single_cylinder_representative(alphabet=Alphabet(name='lower'))
            sage: p
            a b c d e f g
            c f b g e d a
            sage: p.stratum() == C
            True

            sage: C = Stratum([2], k=1)
            sage: C.single_cylinder_representative()
            Traceback (most recent call last):
            ...
            ValueError: no 1,1-square-tiled surfaces in this stratum try again with H_2(2, 0)
            sage: C = Stratum([1,1], k=1)
            sage: C.single_cylinder_representative()
            Traceback (most recent call last):
            ...
            ValueError: no 1,1-square-tiled surfaces in this stratum try again with H_2(1^2, 0^2)
        """
        genus = self.surface_genus()
        nb_real_zeros = sum(map(bool, self.signature()), 0)
        nb_fake_zeros = self.signature().count(0)

        if genus == 2 and nb_real_zeros == 1 and nb_fake_zeros < 1:
            raise ValueError("no 1,1-square-tiled surfaces in this stratum try again with H_2(2, 0)")
        elif genus == 2 and nb_real_zeros == 2 and nb_fake_zeros < 2:
            raise ValueError("no 1,1-square-tiled surfaces in this stratum try again with H_2(1^2, 0^2)")

        return self.one_component().single_cylinder_representative(alphabet, reduced)

    def single_cylinder_origami(self):
        r"""
        Returns an origami associated to a single cylinder permutation representative.

        Returns an origami in this connected component having a single vertical
        cylinder and a single horizontal cylinder.

        Examples::

            sage: from surface_dynamics import Stratum

            sage: C = Stratum([4], k=1)
            sage: O = C.single_cylinder_origami()
            sage: O
            (1,2,3,4,5)
            (1,4,3,5,2)
            sage: O.stratum() == Stratum([4], k=1)
            True
            sage: C = Stratum([2,0], k=1)
            sage: O = C.single_cylinder_origami()
            sage: O
            (1,2,3,4)
            (1,3,2,4)
            sage: O.stratum() == Stratum([2], k=1)
            True
        """
        return self.single_cylinder_representative(reduced=False).to_origami()


class AbelianStratumComponent(StratumComponent):
    r"""
    Connected component of Abelian stratum.

    .. warning::

        Internal class! Do not use directly!
    """
    _name = 'c'

    def spin(self):
        r"""
        Return ``None`` since surfaces in this component have no spin.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: c = Stratum([1,1,1,1], k=1).unique_component(); c
            H_3(1^4)^c
            sage: c.spin() is None
            True
        """
        return None

    def permutation_representative(self, left_degree=None, reduced=True, alphabet=None, relabel=True):
        r"""
        Returns the Zorich representative of this connected component.

        Zorich constructs explicitly interval exchange
        transformations for each stratum in [Zor08]_.

        INPUT:

        - ``reduced`` - boolean (default: ``True``): whether you
          obtain a reduced or labelled permutation

        - ``alphabet`` - an alphabet or ``None``: whether you want to
          specify an alphabet for your permutation

        - ``left_degree`` - the degree of the singularity on the left of the
          interval.

        OUTPUT:

        permutation -- a permutation which lives in this component

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: c = Stratum([1,1,1,1], k=1).unique_component()
            sage: p = c.permutation_representative(alphabet="abcdefghi")
            sage: p
            a b c d e f g h i
            e d c f i h g b a
            sage: p.stratum_component()
            H_3(1^4)^c

            sage: cc = Stratum([3,2,1,0], k=1).unique_component()
            sage: p = cc.permutation_representative(left_degree=3); p
            0 1 2 3 4 5 6 7 8 9 10
            4 3 7 6 5 10 9 8 2 0 1
            sage: p.stratum_component()
            H_4(3, 2, 1, 0)^c
            sage: p.marking().left()
            4
            sage: p.rauzy_diagram()  # not tested
            Rauzy diagram with 1060774 permutations

            sage: p = cc.permutation_representative(left_degree=2); p
            0 1 2 3 4 5 6 7 8 9 10
            4 3 5 7 6 10 9 8 2 0 1
            sage: p.stratum_component()
            H_4(3, 2, 1, 0)^c
            sage: p.marking().left()
            3
            sage: p.rauzy_diagram()  # not tested
            Rauzy diagram with 792066 permutations

            sage: p = cc.permutation_representative(left_degree=1); p
            0 1 2 3 4 5 6 7 8 9 10
            5 4 3 7 6 8 10 9 2 0 1
            sage: p.stratum_component()
            H_4(3, 2, 1, 0)^c
            sage: p.marking().left()
            2
            sage: p.rauzy_diagram()  # not tested
            Rauzy diagram with 538494 permutations

            sage: p = cc.permutation_representative(left_degree=0); p
            0 1 2 3 4 5 6 7 8 9 10
            4 2 7 6 5 10 9 8 1 3 0
            sage: p.stratum_component()
            H_4(3, 2, 1, 0)^c
            sage: p.marking().left()
            1
            sage: p.rauzy_diagram()  # not tested
            Rauzy diagram with 246914 permutations
        """
        stratum = self.stratum()

        g = stratum.surface_genus()
        zeros = list(m for m in stratum.signature() if m)
        n = stratum.signature().count(0)

        if left_degree is not None and left_degree != 0:
            i = zeros.index(left_degree)
            zeros.insert(0, zeros.pop(i))

        if alphabet is None:
            alphabet = range(stratum.dimension())

        l0 = list(range(0, 4*g-3))
        l1 = [4, 3, 2]
        for k in range(5, 4*g-6, 4):
            l1 += [k, k+3, k+2, k+1]
        l1 += [1, 0]
        k = 3
        for d in zeros:
            for i in range(d-1):
                del l0[l0.index(k)]
                del l1[l1.index(k)]
                k += 2
            k += 2

        if n != 0:
            interval = range(4*g-3, 4*g-3+n)

            if left_degree == 0:
                k = l0.index(4)
                l0[k:k] = interval
                l1[-1:-1] = interval
            else:
                l0[1:1] = interval
                l1.extend(interval)

        if reduced:
            from surface_dynamics.interval_exchanges.reduced import ReducedPermutationIET
            p = ReducedPermutationIET([l0, l1])

        else:
            from surface_dynamics.interval_exchanges.labelled import LabelledPermutationIET
            p = LabelledPermutationIET([l0, l1])

        p.alphabet(alphabet)
        return p

    def lyapunov_exponents_approx(self, **kargs):
        r"""
        Return the approximate Lyapunov exponents of the KZ-cocycle.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: Stratum([2], k=1).unique_component().lyapunov_exponents_approx(nb_iterations=2**21)  # abs tol .05
            [1.000, 0.333]

            sage: H4hyp, H4odd = Stratum([4], k=1).components()
            sage: H4hyp.lyapunov_exponents_approx(nb_iterations=2**21) # abs tol .05
            [1.000, 0.616, 0.184]
            sage: H4odd.lyapunov_exponents_approx(nb_iterations=2**21) # abs tol .05
            [1.000, 0.418, 0.182]
        """
        perm = self.permutation_representative(reduced=False)
        return perm.lyapunov_exponents_approx(**kargs)

    lyapunov_exponents = lyapunov_exponents_approx

    # TODO
    # def sum_of_lyapunov_exponents(self)
    # TODO
    # def volume(self)
    # TODO
    # def carea(self)

    def random_standard_permutation(self, nsteps=64):
        r"""
        Perform a random walk on rauzy diagram stopped on a standard permutation.

        INPUT:

        - ``nsteps`` - integer or None - perform nsteps and then stops as soon
          as a Strebel differential is found.

        At each step, with probability 1/3 we perform one of the following
        moves:

        - exchange top,bottom and left,right (proba 1/10)

        - top rauzy move (proba 9/20)

        - bot rauzy move (proba 9/20)

        EXAMPLES:

            sage: from surface_dynamics import Stratum

            sage: C = Stratum([10], k=1).hyperelliptic_component()
            sage: p = C.random_standard_permutation(); p   # random
            0 1 2 3 4 5 6 7 8 9 10 11
            11 10 9 8 7 6 5 4 3 2 1 0
            sage: p.stratum_component()
            H_6(10)^hyp

            sage: C = Stratum([6,4,2], k=1).odd_component(); C
            H_7(6, 4, 2)^odd
            sage: p = C.random_standard_permutation(); p  # random
            0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
            15 2 14 12 3 11 6 10 8 5 9 13 7 4 1 0
            sage: p.stratum_component()
            H_7(6, 4, 2)^odd

            sage: C = Stratum([2,2,2,2], k=1).even_component(); C
            H_5(2^4)^even
            sage: p = C.random_standard_permutation(); p  # random
            0 1 2 3 4 5 6 7 8 9 10 11 12
            12 4 9 11 8 3 7 6 1 10 2 5 0
            sage: p.stratum_component()
            H_5(2^4)^even

            sage: C = Stratum([32], k=1).odd_component(); C
            H_17(32)^odd
            sage: p = C.random_standard_permutation(); p  # random
            0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33
            33 30 10 3 32 19 11 28 4 14 24 15 21 20 9 12 25 6 2 29 26 23 27 13 8 1 18 17 16 31 7 22 5 0
            sage: p.stratum_component()
            H_17(32)^odd
        """
        import sage.misc.prandom as prandom

        p = self.permutation_representative()
        if nsteps is None:
            nsteps = 64 * self.stratum().dimension()
        d = len(p)-1

        for _ in range(nsteps):
            rd = prandom.random()
            if rd < 0.1:   # (inplace) symmetric with proba 1/10
                p._inversed_twin()
                p._reversed_twin()
            elif rd < .55: # (inplace) rauzy move top with proba 9/20
                p._move(1, 0, 1, p._twin[0][0])
            else:          # (inplace) rauzy move bot with proba 9/20
                p._move(0, 0, 0, p._twin[1][0])

        while not p.is_standard():
            rd = prandom.random()
            if rd < 0.1:   # (inplace) symmetric with proba 1/10
                p._inversed_twin()
                p._reversed_twin()
            elif rd < .55: # (inplace) rauzy move top with proba 9/20
                p._move(1, 0, 1, p._twin[0][0])
            else:          # (inplace) rauzy move bot with proba 9/20
                p._move(0, 0, 0, p._twin[1][0])

        return p

    def rauzy_diagram(self, *args, **kwds):
        r"""
        Returns the extended Rauzy diagram associated to this connected component.

        OUTPUT:

        rauzy diagram -- the Rauzy diagram associated to this stratum

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: c = Stratum([0], k=1).components()[0]
            sage: r = c.rauzy_diagram()
        """
        if kwds.get("left_degree",None) is not None:
            return self.permutation_representative(*args, **kwds).rauzy_diagram()
        return self.permutation_representative(*args,**kwds).rauzy_diagram(extended=True)

    def rauzy_class_cardinality(self, left_degree=None, reduced=True):
        r"""
        Rauzy diagram cardinality for connected components.

        Returns the cardinality of the extended Rauzy diagram associated to this
        connected component.

        If left_degree is provided then it returns the cardinality of the Rauzy
        diagram with a singularity of that degree attached on the left.
        Otherwise it returns the cardinality of the extended Rauzy diagram.

        INPUT:

        - ``left_degree`` - the degree to be attached to the singularity on the
          left

        - ``reduced`` - boolean (default: True) - consider the cardinality of
          reduced or extended Rauzy diagram

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: a = Stratum({1:4}, k=1).unique_component()
            sage: a
            H_3(1^4)^c
            sage: a.rauzy_diagram()
            Rauzy diagram with 1255 permutations
            sage: a.rauzy_class_cardinality()
            1255

            sage: cc = Stratum([3,2,1], k=1).unique_component()
            sage: cc.rauzy_diagram(left_degree=3)   # long time
            Rauzy diagram with 96434 permutations
            sage: cc.rauzy_class_cardinality(left_degree=3)
            96434

            sage: cc.rauzy_diagram(left_degree=2)   # long time
            Rauzy diagram with 72006 permutations
            sage: cc.rauzy_class_cardinality(left_degree=2)
            72006

            sage: cc.rauzy_diagram(left_degree=1)   # long time
            Rauzy diagram with 48954 permutations
            sage: cc.rauzy_class_cardinality(left_degree=1)
            48954

            sage: a = Stratum({1:8}, k=1).unique_component()
            sage: a
            H_5(1^8)^c
            sage: a.rauzy_class_cardinality()
            55184875

        Cardinalities for labeled Rauzy classes instead of reduced::

            sage: cc = Stratum([2,1,1], k=1).unique_component()
            sage: cc.rauzy_diagram(left_degree=2, reduced=False)
            Rauzy diagram with 3676 permutations
            sage: cc.rauzy_class_cardinality(left_degree=2, reduced=False)
            3676

            sage: cc.rauzy_diagram(left_degree=1, reduced=False)
            Rauzy diagram with 3774 permutations
            sage: cc.rauzy_class_cardinality(left_degree=1,reduced=False)
            3774

            sage: cc = Stratum([2,1,1,0], k=1).unique_component()
            sage: cc.rauzy_diagram(left_degree=2, reduced=False) # long time
            Rauzy diagram with 33084 permutations
            sage: cc.rauzy_diagram(left_degree=1, reduced=False) # long time
            Rauzy diagram with 33966 permutations
            sage: cc.rauzy_diagram(left_degree=0, reduced=False) # long time
            Rauzy diagram with 30828 permutations

            sage: cc.rauzy_class_cardinality(left_degree=2, reduced=False)
            33084
            sage: cc.rauzy_class_cardinality(left_degree=1, reduced=False)
            33966
            sage: cc.rauzy_class_cardinality(left_degree=0, reduced=False)
            30828
        """
        import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rdc

        profile = list(map(lambda x: x+1, self.stratum().signature()))
        s = len(self.stratum().signature())

        if left_degree is not None:
            assert isinstance(left_degree, (int,Integer)), "if not None, left_degree should be an integer"
            left_degree = int(left_degree) + 1
            assert left_degree in profile, "if not None, the degree should be one of the degree of the stratum"

            if reduced:
                return Integer(rdc.gamma_irr(profile,left_degree))
            return Rational((1,2)) * (Partition(profile).centralizer_size() /
                    (left_degree * profile.count(left_degree)) *
                    Integer(rdc.gamma_irr(profile,left_degree)))

        if reduced:
            return Integer(rdc.gamma_irr(profile))

        raise NotImplementedError("not known formula for labeled extended Rauzy classes")

    def standard_permutations_number(self, left_degree=None):
        r"""
        Return the number of standard permutations in the Rauzy class associated
        to this connected component.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: cc = Stratum([3,1], k=1).unique_component()
            sage: sum(1 for p in cc.rauzy_diagram() if p.is_standard())
            24
            sage: cc.standard_permutations_number()
            24

            sage: sum(1 for p in cc.rauzy_diagram(left_degree=3) if p.is_standard())
            16
            sage: cc.standard_permutations_number(left_degree=3)
            16

            sage: sum(1 for p in cc.rauzy_diagram(left_degree=1) if p.is_standard())
            8
            sage: cc.standard_permutations_number(left_degree=1)
            8

            sage: cc = Stratum({1:10}, k=1).unique_component()
            sage: cc
            H_6(1^10)^c
            sage: cc.standard_permutations_number()
            59520825
        """
        import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rdc

        profile = [x+1 for x in self.stratum().signature()]

        if left_degree is not None:
            assert isinstance(left_degree, (int,Integer)), "if not None, left_degree should be an integer"
            left_degree = int(left_degree) + 1
            assert left_degree in profile, "if not None, the degree should be one of the degree of the stratum"
            return Integer(rdc.number_of_standard_permutations(profile,left_degree))

        return Integer(rdc.number_of_standard_permutations(profile))

    def standard_permutations(self):
        r"""
        Return the set of standard permutations.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: C = Stratum([4], k=1).odd_component()
            sage: C
            H_3(4)^odd
            sage: for p in C.standard_permutations(): print("%s\n***********" % p)
            0 1 2 3 4 5
            5 2 1 4 3 0
            ***********
            0 1 2 3 4 5
            5 3 1 4 2 0
            ***********
            0 1 2 3 4 5
            5 4 1 3 2 0
            ***********
            0 1 2 3 4 5
            5 2 4 1 3 0
            ***********
            0 1 2 3 4 5
            5 4 2 1 3 0
            ***********
            0 1 2 3 4 5
            5 2 4 3 1 0
            ***********
            0 1 2 3 4 5
            5 3 2 4 1 0
            ***********
        """
        p = self.permutation_representative(reduced=True)
        return sorted(q for q in p.rauzy_diagram(symmetric=True) if q.is_standard())

    def one_cylinder_diagram(self):
        r"""
        Return a diagram with one cylinder in this connected component.

        The diagram returned is the one deduced from the method
        permutation_representative.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: A = Stratum([2,2], k=1).odd_component()
            sage: c = A.one_cylinder_diagram()
            sage: c
            (0,5,1,3,2,4)-(0,5,4,3,2,1)
            sage: c.stratum_component()
            H_3(2^2)^odd

            sage: A = Stratum([3,3], k=1).non_hyperelliptic_component()
            sage: c = A.one_cylinder_diagram()
            sage: c
            (0,7,3,2,1,5,4,6)-(0,7,6,5,4,3,2,1)
            sage: c.stratum_component()
            H_4(3^2)^nonhyp
        """
        from .separatrix_diagram import CylinderDiagram
        t = self.permutation_representative(reduced=True).to_standard()
        return CylinderDiagram([(t[1][1:],t[0][-2::-1])],check=True)

    def cylinder_diagram_iterator(self, ncyls=None, up_to_symmetry=True, force_computation=False):
        r"""
        An iterator over the cylinder diagrams.

        INPUT::

        - ``ncyls`` -- (optional) a fixed number of cylinders

        - ``up_to_symmetry`` - (boolean, default ``True``) to return only
          cylinder diagrams up to horizontal and vertical symmetry.

        - ``force_computation`` -- (default ``False``) whether the database
          should be used or not

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: A = Stratum([1,1,1,1], k=1)
            sage: cc = A.unique_component()
            sage: it = cc.cylinder_diagram_iterator(3)
            sage: cyl = next(it); cyl
            (0,7)-(0,5) (1,6)-(1,4) (2,4,3,5)-(2,7,3,6)
            sage: cyl.stratum_component()
            H_3(1^4)^c
            sage: cyl.ncyls()
            3

        Note that if you set ``force_computation`` to ``True`` the order of the
        iteration might be different and you might obtain cylinder diagram with some
        symmetries applied::

            sage: # long time
            sage: C1 = list(cc.cylinder_diagram_iterator(3, force_computation=False))
            sage: C2 = list(cc.cylinder_diagram_iterator(3, force_computation=True))
            sage: assert len(C1) == len(C2)
            sage: isoms = []
            sage: for c in C1:
            ....:     isom = []
            ....:     for i,cc in enumerate(C2):
            ....:         if c.is_isomorphic(cc) or \
            ....:            c.is_isomorphic(cc.horizontal_symmetry()) or \
            ....:            c.is_isomorphic(cc.vertical_symmetry()) or \
            ....:            c.is_isomorphic(cc.inverse()):
            ....:              isom.append(i)
            ....:     assert len(isom) == 1, isom
            ....:     isoms.extend(isom)
            sage: assert sorted(isoms) == list(range(len(C1)))
        """
        if ncyls is not None:
            if not isinstance(ncyls, (int,Integer)):
                raise TypeError("ncyls should be None or an integer")
            if ncyls < 0 or ncyls > self.stratum().surface_genus() + len(self.stratum().signature()) - 1:
                raise ValueError("ncyls is not valid")

        if not force_computation:
            from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
            CDB = CylinderDiagrams()
            if CDB.has_component(self):
                iterator = CDB.get_iterator(self, ncyls)
                if up_to_symmetry:
                    return iterator
                else:
                    return _cylinder_diagrams_with_symmetric(iterator)

        iterator = self._cylinder_diagram_iterator(ncyls)
        if up_to_symmetry:
            return iterator
        else:
            return _cylinder_diagrams_with_symmetric(iterator)

    def _cylinder_diagram_iterator(self, ncyls):
        r"""
        Default implementation for cylinder diagrams for connected stratum.
        """
        return self.stratum().cylinder_diagram_iterator(ncyls, True, True)

    def cylinder_diagrams(self, ncyls=None, up_to_symmetry=True, force_computation=False):
        r"""
        Return the list of cylinder diagrams associated to this component.

        INPUT:

        - ``ncyls`` - integer or list of integers (default: None) - consider
          only the cylinder diagrams with a given number of cylinders.

        - ``up_to_symmetry`` - (boolean, default ``True``) to return only
          cylinder diagrams up to horizontal and vertical symmetry.

        - ``force_computation`` -- boolean (default ``False``). If ``True``, the
          database of cylinder diagrams is not used.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: C = Stratum([1,1,1,1], k=1).unique_component(); C
            H_3(1^4)^c
            sage: for c in C.cylinder_diagrams(6): print(c)
            (0,1)-(7) (2)-(0) (3)-(4) (4,7)-(5,6) (5)-(1) (6)-(2,3)
            (0,1)-(7) (2)-(1) (3)-(0) (4,7)-(5,6) (5)-(4) (6)-(2,3)
            (0,3)-(6,7) (1,2)-(4,5) (4)-(1) (5)-(3) (6)-(2) (7)-(0)
            (0,3)-(6,7) (1,2)-(4,5) (4)-(3) (5)-(0) (6)-(2) (7)-(1)
        """
        return sorted(self.cylinder_diagram_iterator(ncyls, up_to_symmetry, force_computation))

    def cylinder_diagrams_number(self, ncyls=None, up_to_symmetry=True, force_computation=False):
        r"""
        Return the number of cylinder diagrams.

        INPUT:

        - ``ncyls`` - integer or list of integers (default: None) - restrict the
          counting to a given number of cylinders.

        - ``up_to_symmetry`` - (boolean, default ``True``) to count only
          cylinder diagrams up to horizontal and vertical symmetry.

        - ``force_computation`` - (default: ``False``) whether we use the
          database or compute explicitly using the generation algorithm.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: C = Stratum([3,1], k=1).unique_component()
            sage: C.cylinder_diagrams_number(1)
            2
            sage: C.cylinder_diagrams_number(2)
            12
            sage: C.cylinder_diagrams_number(3)
            16
            sage: C.cylinder_diagrams_number(4)
            4

        Note that when setting ``force_computation`` to ``True`` we got the
        same numbers::

            sage: for i in range(1,5):
            ....:     print(C.cylinder_diagrams_number(i, force_computation=True))
            2
            12
            16
            4

            sage: C = Stratum([6], k=1)
            sage: C_hyp = C.hyperelliptic_component()
            sage: C_odd = C.odd_component()
            sage: C_even = C.even_component()
            sage: for i in range(1,5): print(C.cylinder_diagrams_number(i))
            16
            76
            130
            67
            sage: for i in range(1,5): print(C_hyp.cylinder_diagrams_number(i))
            1
            3
            8
            4
            sage: for i in range(1,5): print(C_odd.cylinder_diagrams_number(i))
            11
            49
            80
            42

            sage: for i in range(1,5):
            ....:     print(C_even.cylinder_diagrams_number(i, True, False))
            4
            24
            42
            21
            sage: for i in range(1,5):                                      # long time
            ....:     print(C_even.cylinder_diagrams_number(i, True, True)) # long time
            4
            24
            42
            21
        """
        if ncyls is not None and ncyls > self.stratum().surface_genus() + len(self.stratum().signature()) - 1:
            return 0

        if not force_computation:
            from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
            CDB = CylinderDiagrams()
            if CDB.has_component(self):
                if up_to_symmetry:
                    return CDB.count(self, ncyls)
                else:
                    return sum(4 // (1 + sum(cd.symmetries())) for cd in CDB.get_iterator(self, ncyls))

        if up_to_symmetry:
            return sum(1 for _ in self.cylinder_diagram_iterator(ncyls, True, True))
        else:
            return sum(4 // (1 + sum(cd.symmetries())) for cd in self.cylinder_diagram_iterator(ncyls, True, True))

    def one_origami(self):
        r"""
        Returns an origami in this component

        The origami returned has the minimal number of squares and one
        cylinder. It is obtained from the permutation representative of the
        stratum.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: a = Stratum([2,2], k=1).one_component()
            sage: a.one_origami().stratum()
            H_3(2^2)

            sage: Stratum([3,2,1], k=1).unique_component().one_origami().stratum()
            H_4(3, 2, 1)
        """
        from surface_dynamics.flat_surfaces.origamis.origami_dense import Origami_dense_pyx
        t = self.permutation_representative(reduced=True).to_standard()
        t.alphabet(range(len(t)))
        h_perm = list(range(1,len(t)-1)) + [0]
        v_perm = t[1][1:]
        return Origami_dense_pyx(tuple(h_perm), tuple(v_perm))

    def origami_iterator(self, n, reduced=True, primitive=False):
        r"""
        Iterator through the set of origamis with ``n`` squares in this stratum.

        The output origamis are in normal form. But be careful as there may be
        repetition in the output!

        INPUT:

        - ``n`` - integer - the number of squares

        - ``reduced`` - boolean (default: ``True``)

        - ``primitive`` - boolean (default: ``False``)

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: cc = Stratum([6], k=1).even_component()
            sage: it = cc.origami_iterator(13)
            sage: o = next(it)
            sage: o
            (1,2,3,4,5,6,7,8,9,10,11,12,13)
            (1,8,2)(3,12,6,11,5,13,7,9)(4,10)
            sage: o.stratum_component()
            H_4(6)^even
        """
        reduced = not reduced
        primitive = not primitive
        for c in self.cylinder_diagram_iterator():
            do_h, do_v, do_i = c.symmetries()
            do_h = not do_h
            do_v = not do_v
            do_i = not do_i
            for o in c.origami_iterator(n):
                if ((reduced or o.is_reduced()) and (primitive or o.is_primitive())):
                    o.relabel(inplace=True)
                    yield o

                    if do_h:
                        oo = o.horizontal_symmetry()
                        oo.relabel(inplace=True)
                        yield oo

                    if do_v:
                        oo = o.vertical_symmetry()
                        oo.relabel(inplace=True)
                        yield oo

                    if do_i:
                        oo = o.inverse()
                        oo.relabel(inplace=True)
                        yield oo

    def origamis(self, n, reduced=True, primitive=False):
        r"""
        Return the set of origamis with n squares in this stratum.

        INPUT:

        - ``n`` - integer - the number of squares

        - ``reduced`` - boolean (default: True)

        - ``primitive`` - boolean (default: False)

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: H11_hyp = Stratum([1,1], k=1).hyperelliptic_component()
            sage: len(H11_hyp.origamis(6))
            88

            sage: T6 = H11_hyp.arithmetic_teichmueller_curves(6)
            sage: len(T6)
            5
            sage: sum(t.veech_group().index() for t in T6)
            88

            sage: H4_odd = Stratum([4], k=1).odd_component()
            sage: len(H4_odd.origamis(6))
            155
            sage: T6 = H4_odd.arithmetic_teichmueller_curves(6)
            sage: sum(t.veech_group().index() for t in T6)
            155
        """
        return set(self.origami_iterator(n, reduced, primitive))

    def arithmetic_teichmueller_curves(self, n, primitive=False):
        r"""
        Return the arithmetic Teichmueller curves in that component of stratum.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: A = Stratum([2], k=1).hyperelliptic_component(); A
            H_2(2)^hyp
            sage: for i in range(3,10):
            ....:     print("%d %d" % (i,len(A.arithmetic_teichmueller_curves(i))))
            3 1
            4 1
            5 2
            6 1
            7 2
            8 1
            9 2

            sage: A = Stratum([1,1], k=1).hyperelliptic_component(); A
            H_2(1^2)^hyp
            sage: for i in range(4,10):
            ....:    T = A.arithmetic_teichmueller_curves(i)
            ....:    T_prim = list(filter(lambda t:t.origami().is_primitive(), T))
            ....:    print("%d %d %d" % (i,len(T),len(T_prim)))
            4 2 1
            5 1 1
            6 5 2
            7 2 2
            8 4 2
            9 4 2

            sage: A = Stratum([4], k=1).hyperelliptic_component(); A
            H_3(4)^hyp
            sage: for i in range(5,10):
            ....:    print("%d %d" % (i,len(A.arithmetic_teichmueller_curves(i))))
            5 2
            6 4
            7 3
            8 3
            9 4
        """
        from .origamis.teichmueller_curve import TeichmuellerCurvesOfOrigamis
        tcurves = TeichmuellerCurvesOfOrigamis(self.origamis(n),assume_normal_form=True)
        if primitive:
            return [tcurve for tcurve in tcurves if tcurve.origami().is_primitive()]
        return tcurves

    def single_cylinder_representative(self, alphabet=None, reduced=True):
        r"""
        Returns a single cylinder permutation representative.

        Returns a cylindric permutation representative of this connected
        stratum (or non-hyperelliptic component) such that the associated
        square-tiled surface is made of a single cylinder of height one in both
        horizontal and vertical direction.

        Such representatives were constructed for every stratum of Abelian
        differentials by Jeffreys [Jef19]_.

        INPUT:

        - ``alphabet`` -- an optional alphabet for the permutation representative

        - ``reduced`` (boolean, default ``True``) -- whether to return a reduced
          permutation (ie without labels)

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: cc = Stratum([1,1,1,1], k=1).unique_component()
            sage: p = cc.single_cylinder_representative()
            sage: p
            0 1 2 3 4 5 6 7 8
            2 6 5 3 1 8 4 7 0
            sage: p.stratum_component() == cc
            True

            sage: cc = Stratum([2,1,1], k=1).unique_component()
            sage: p = cc.single_cylinder_representative()
            sage: p
            0 1 2 3 4 5 6 7
            2 6 4 1 7 5 3 0
            sage: p.stratum_component() == cc
            True

            sage: cc = Stratum([3,3], k=1).non_hyperelliptic_component()
            sage: p = cc.single_cylinder_representative(alphabet=Alphabet(name='lower'))
            sage: p
            a b c d e f g h i
            c i g f h e b d a
            sage: p.stratum_component() == cc
            True
        """
        from surface_dynamics.flat_surfaces.single_cylinder import (cylinder_concatenation,
                only_even_2, only_odds_11, odd_zeros_one_one)
        from surface_dynamics.interval_exchanges.constructors import GeneralizedPermutation

        zeros = self.stratum().signature()
        real_zeros = [z for z in zeros if z != 0]
        nb_fake_zeros = len(zeros) - len(real_zeros)
        odd_zeros = [z for z in real_zeros if z % 2 == 1]
        even_zeros = [z for z in real_zeros if z % 2 == 0]

        fk_zeros_perm = GeneralizedPermutation([0],[0])
        mk_pt_perm = GeneralizedPermutation([0,1],[1,0])
        for i in range(nb_fake_zeros):
            fk_zeros_perm = cylinder_concatenation(fk_zeros_perm, mk_pt_perm)

        if even_zeros == [2]:
            perm = only_even_2(odd_zeros)
        elif odd_zeros == [1,1]:
            perm = only_odds_11(even_zeros)
        else:
            if even_zeros:
                even_perm = Stratum(even_zeros, k=1).odd_component().single_cylinder_representative()
            else:
                even_perm = GeneralizedPermutation([0],[0])
            odd_perm = odd_zeros_one_one(odd_zeros)
            perm = cylinder_concatenation(even_perm,odd_perm)

        perm = cylinder_concatenation(fk_zeros_perm, perm)

        if alphabet is not None:
            perm.alphabet(alphabet)

        return perm

    def single_cylinder_origami(self):
        r"""
        Returns an origami associated to a single cylinder permutation representative.

        Returns an origami in this connected (or non-hyperelliptic) component
        having a single vertical cylinder and a single horizontal cylinder.

        Such representatives were constructed for every stratum of Abelian
        differentials by Jeffreys [Jef19]_.

        Examples::

            sage: from surface_dynamics import Stratum

            sage: cc = Stratum([4], k=1).odd_component()
            sage: O = cc.single_cylinder_origami()
            sage: O
            (1,2,3,4,5)
            (1,4,3,5,2)
            sage: O.stratum_component() == cc
            True
            sage: cc = Stratum([5,3], k=1).unique_component()
            sage: O = cc.single_cylinder_origami()
            sage: O
            (1,2,3,4,5,6,7,8,9,10)
            (1,9,8,10,6,7,4,3,5,2)
            sage: O.stratum_component() == cc
            True
            sage: cc = Stratum([4,2], k=1).even_component()
            sage: O = cc.single_cylinder_origami()
            sage: O
            (1,2,3,4,5,6,7,8)
            (1,3,7,5,6,8,4,2)
            sage: O.stratum_component() == cc
            True
        """
        return self.single_cylinder_representative(reduced=False).to_origami()

ASC = AbelianStratumComponent

class HypAbelianStratumComponent(ASC):
    """
    Hyperelliptic component of Abelian stratum.
    """
    _name = 'hyp'

    def spin(self):
        r"""
        Return the spin parity of hyperelliptic stratum.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

        For the strata `H(2g-2)`::

            sage: c = Stratum([0], k=1).hyperelliptic_component()
            sage: c.spin()
            1
            sage: p = c.permutation_representative()
            sage: p.arf_invariant()
            1

            sage: c = Stratum([2], k=1).hyperelliptic_component()
            sage: c.spin()
            1
            sage: p = c.permutation_representative()
            sage: p.arf_invariant()
            1

            sage: c = Stratum([4], k=1).hyperelliptic_component()
            sage: c.spin()
            0
            sage: p = c.permutation_representative()
            sage: p.arf_invariant()
            0

        For the strata `H(g-1,g-1)`::

            sage: c = Stratum([2,2], k=1).hyperelliptic_component()
            sage: c.spin()
            0
            sage: p = c.permutation_representative()
            sage: p.arf_invariant()
            0

            sage: c = Stratum([4,4], k=1).hyperelliptic_component()
            sage: c.spin()
            1
            sage: p = c.permutation_representative()
            sage: p.arf_invariant()
            1
        """
        z = tuple(m for m in self.stratum().signature() if m > 0)
        if not z:
            return Integer(1)
        elif len(z) == 1:
            return Integer(((self.stratum().surface_genus()+1)//2) % 2)
        elif len(z) == 2:
            if z[0] % 2:
                return None
            return Integer(((self.stratum().surface_genus()+1)//2) %2)

    def permutation_representative(self, left_degree=None, reduced=True, alphabet=None, relabel=True):
        r"""
        Returns the Zorich representative of this connected component.

        Zorich constructs explicitly interval exchange
        transformations for each stratum in [Zor08]_.

        INPUT:

        - ``reduced`` - boolean (default: ``True``): whether you obtain
          a reduced or labelled permutation

        - ``alphabet`` - alphabet or ``None`` (default: ``None``):
          whether you want to specify an alphabet for your
          representative

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: c = Stratum([0], k=1).hyperelliptic_component()
            sage: p = c.permutation_representative()
            sage: p
            0 1
            1 0
            sage: p.stratum_component()
            H_1(0)^hyp

            sage: c = Stratum([0,0], k=1).hyperelliptic_component()
            sage: p = c.permutation_representative(alphabet="abc")
            sage: p
            a b c
            c b a
            sage: p.stratum_component()
            H_1(0^2)^hyp

            sage: c = Stratum([2,2], k=1).hyperelliptic_component()
            sage: p = c.permutation_representative(alphabet="ABCDEFGHIJKL")
            sage: p
            A B C D E F G
            G F E D C B A
            sage: c = Stratum([1,1,0], k=1).hyperelliptic_component()
            sage: p = c.permutation_representative(left_degree=1); p
            0 1 2 3 4 5
            5 1 4 3 2 0
            sage: p.marking().left()
            2
            sage: p.rauzy_diagram()
            Rauzy diagram with 90 permutations

            sage: p = c.permutation_representative(left_degree=0); p
            0 1 2 3 4 5
            5 3 2 1 4 0
            sage: p.marking().left()
            1
            sage: p.rauzy_diagram()
            Rauzy diagram with 20 permutations
        """
        g = self._stratum.surface_genus()
        n = self._stratum.signature().count(0)
        m = sum(x != 0 for x in self._stratum.signature())

        if left_degree is not None:
            if not isinstance(left_degree, (int,Integer)) or left_degree not in self.stratum().signature():
                raise ValueError("left_degree (=%d) should be one of the degree"%left_degree)

        if m == 0:  # on the torus
            if n == 1:
                l0 = [0, 1]
                l1 = [1, 0]
            elif n == 2:
                l0 = [0, 1, 2]
                l1 = [2, 1, 0]
            else:
                l0 = list(range(1, n+2))
                l1 = [n+1] + list(range(1, n+1))

        elif m == 1:  # H(2g-2,0^n) or H(0,2g-2,0^(n-1))
            l0 = list(range(1, 2*g+1))
            l1 = list(range(2*g, 0, -1))
            interval = list(range(2*g+1, 2*g+n+1))

            if left_degree == 0:
                l0[-1:-1] = interval
                l1[-1:-1] = interval
            else:
                l0[1:1] = interval
                l1[1:1] = interval

        else:  # H(g-1,g-1,0^n) or H(0,g-1,g-1,0^(n-1))
            l0 = list(range(1, 2*g+2))
            l1 = list(range(2*g+1, 0, -1))
            interval = list(range(2*g+2, 2*g+n+2))

            if left_degree == 0:
                l0[-1:-1] = interval
                l1[-1:-1] = interval
            else:
                l0[1:1] = interval
                l1[1:1] = interval

        if reduced:
            from surface_dynamics.interval_exchanges.reduced import ReducedPermutationIET
            p = ReducedPermutationIET([l0, l1])

        else:
            from surface_dynamics.interval_exchanges.labelled import LabelledPermutationIET
            p = LabelledPermutationIET([l0, l1])
        if alphabet is not None:
            p.alphabet(alphabet)
        elif relabel:
            p.alphabet(range(len(p)))
        return p

    def rauzy_class_cardinality(self, left_degree=None, reduced=True):
        r"""
        Return the cardinality of the extended Rauzy diagram associated to the
        hyperelliptic component

        The cardinality of the Rauzy diagram or extended Rauzy diagram
        associated to `H_{hyp}(2g-2,0^k)` or `H_{hyp}(g-1,g-1,0^k)` depends only
        on the dimension `d` of the initial stratum `\mathcal{H}_{hyp}(2g-2)`
        for which `d=2g` or `\mathcal{H}_{hyp}(g-1,g-1)` for which
        `d=2g+1` and the number of fake zeros `k`. The formula is

        .. MATH::

            \binom{d+k+1}{k} (2^{d-1}-1) + d \binom{d+k}{k-1}

        INPUT:

        - ``left_degree`` - integer - the degree of the singularity attached at
          the left of the interval.

        - ``reduced`` - boolean (default: True) - if False, consider labeled
          Rauzy diagrams instead of reduced.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

        The case of the torus is a little bit different::

            sage: c = Stratum([0], k=1).hyperelliptic_component()
            sage: c.rauzy_diagram()
            Rauzy diagram with 1 permutation
            sage: c.rauzy_class_cardinality()
            1
            sage: c = Stratum([0,0], k=1).hyperelliptic_component()
            sage: c.rauzy_diagram()
            Rauzy diagram with 3 permutations
            sage: c.rauzy_class_cardinality()
            3

        Examples in genus 2::

            sage: c = Stratum([2,0], k=1).hyperelliptic_component()
            sage: c.rauzy_diagram()
            Rauzy diagram with 46 permutations
            sage: c.rauzy_class_cardinality()
            46

            sage: c.rauzy_diagram(left_degree=2)
            Rauzy diagram with 35 permutations
            sage: c.rauzy_class_cardinality(left_degree=2)
            35

            sage: c.rauzy_diagram(left_degree=0)
            Rauzy diagram with 11 permutations
            sage: c.rauzy_class_cardinality(left_degree=0)
            11
            sage: c.rauzy_diagram(left_degree=0, reduced=False)
            Rauzy diagram with 33 permutations
            sage: c.rauzy_class_cardinality(left_degree=0, reduced=False)
            33

            sage: c = Stratum([1,1,0,0], k=1).hyperelliptic_component()
            sage: c.rauzy_diagram()
            Rauzy diagram with 455 permutations
            sage: c.rauzy_class_cardinality()
            455

            sage: c.rauzy_diagram(left_degree=1)
            Rauzy diagram with 315 permutations
            sage: c.rauzy_class_cardinality(left_degree=1)
            315
            sage: c.rauzy_diagram(left_degree=1, reduced=False)
            Rauzy diagram with 630 permutations
            sage: c.rauzy_class_cardinality(left_degree=1, reduced=False)
            630

            sage: c.rauzy_diagram(left_degree=0)
            Rauzy diagram with 140 permutations
            sage: c.rauzy_class_cardinality(left_degree=0)
            140
            sage: c.rauzy_diagram(left_degree=0, reduced=False)
            Rauzy diagram with 560 permutations
            sage: c.rauzy_class_cardinality(left_degree=0, reduced=False)
            560

        Other examples in higher genus::

            sage: c = Stratum([12,0,0], k=1).hyperelliptic_component()
            sage: c.rauzy_class_cardinality()
            1114200
            sage: c.rauzy_class_cardinality(left_degree=12, reduced=False)
            1965840

            sage: c = Stratum([14], k=1).hyperelliptic_component()
            sage: c.rauzy_class_cardinality()
            32767
        """
        from sage.arith.all import binomial

        if left_degree is not None:
            assert isinstance(left_degree, (int,Integer)), "if not None, left_degree should be an integer"
            assert left_degree in self.stratum().signature(), "if not None, the degree should be one of the degree of the stratum"

        if reduced is False:
            if left_degree is None:
                raise NotImplementedError("no formula known for cardinality of labeled extended Rauzy classes")
            zeros = self.stratum().signature()
            profile = Partition([x+1 for x in zeros])
            if sum(bool(z) for z in self.stratum().signature()) == 1:
                epsilon = 1
            else:
                epsilon = Rational((1,self.stratum().surface_genus()))
            return epsilon * (profile.centralizer_size() /
                    ((left_degree+1) * zeros.count(left_degree)) *
                    self.rauzy_class_cardinality(left_degree=left_degree,reduced=True))

        k = self.stratum().signature().count(0)
        dd = self.stratum().dimension()  # it is d+k
        d = dd-k

        if self.stratum().surface_genus() == 1:
            if k == 0: return 1
            return binomial(dd,2)

        if left_degree is None:
            return binomial(dd+1,k) * (2**(d-1)-1) + d * binomial(dd,k-1)

        if left_degree == 0:
            return binomial(dd,k-1) * (2**(d-1)-1 + d)
        else:
            return binomial(dd,k) * (2**(d-1)-1)

    def random_standard_permutation(self, nsteps=None):
        r"""
        In hyperelliptic component there is only one standard permutation.
        """
        if 0 not in self.stratum().signature():
            return self.permutation_representative()

        raise NotImplementedError("not implemented when there are fake zeros")

    def standard_permutations(self, reduced=True):
        r"""
        Return the standard permutations in this hyperelliptic component.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: Stratum([6], k=1).hyperelliptic_component().standard_permutations()
            [0 1 2 3 4 5 6 7
             7 6 5 4 3 2 1 0]
        """
        if 0 not in self.stratum().signature():
            d = self.stratum().dimension()
            l0 = list(range(d))
            l1 = list(range(d-1,-1,-1))

            if reduced:
                from surface_dynamics.interval_exchanges.reduced import ReducedPermutationIET
                p = ReducedPermutationIET([l0, l1])
            else:
                from surface_dynamics.interval_exchanges.labelled import LabelledPermutationIET
                p = LabelledPermutationIET([l0, l1])

            return [p]

        raise NotImplementedError("not implemented when there are fake zeros")

    def standard_permutations_number(self):
        r"""
        Return the number of standard permutations in this hyperelliptic
        component.
        """
        if 0 not in self.stratum().signature():
            return Integer(1)

        raise NotImplementedError("not implemented when there are fake zeros")

    def _cylinder_diagram_iterator(self, ncyls=None):
        r"""
        Returns the list of cylinder diagrams associated to this hyperelliptic
        component.

        INPUT:

        - ``ncyls`` -- an optional number of cylinders

        EXAMPLES::

            sage: from surface_dynamics import Stratum, CylinderDiagram

            sage: C = Stratum([2,2], k=1).hyperelliptic_component()
            sage: [sum(1 for c in C.cylinder_diagram_iterator(n)) for n in range(1,5)]
            [1, 3, 5, 2]

        When ``ncyls`` is set to ``None``, the iterator can reasonably be used
        with very large data::

            sage: C = Stratum([10,10], k=1).hyperelliptic_component()
            sage: it = C.cylinder_diagram_iterator()
            sage: c = next(it)
            sage: c.is_isomorphic(CylinderDiagram('(0,2,5,1)-(0,2,21,1) (3,4)-(3,6) (6,19)-(4,20) (7,9)-(8,10) (8,12)-(7,11) (10,14)-(9,13) (11,15)-(12,16) (13,17)-(14,18) (16,20)-(15,19) (18,21)-(5,17)'))
            True
            sage: c.stratum_component()
            H_11(10^2)^hyp
            sage: c.ncyls()
            10
        """
        from .separatrix_diagram import hyperelliptic_cylinder_diagram_iterator

        if ncyls is not None:
            ncyls = int(ncyls)
            return filter(lambda c: c.ncyls() == ncyls, self.cylinder_diagram_iterator())

        stratum = self.stratum()

        if 0 in stratum.signature():
            raise ValueError("the stratum has fake zeros")

        z = stratum.signature()

        return hyperelliptic_cylinder_diagram_iterator(len(z)+sum(z))

    def single_cylinder_representative(self, alphabet=None, reduced=True):
        r"""
        Returns a single cylinder permutation representative.

        Returns a permutation representative of a square-tiled surface in this
        component having a single vertical cylinder and a single horizontal cylinder.

        Such representatives were constructed for every stratum of Abelian
        differentials by Jeffreys [Jef19]_.

        INPUT:

        - ``alphabet`` -- an optional alphabet for the permutation representative

        - ``reduced`` (boolean, default ``True``) -- whether to return a reduced
          permutation (ie without labels)

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: cc = Stratum([2,0], k=1).hyperelliptic_component()
            sage: p = cc.single_cylinder_representative(alphabet=Alphabet(name='upper'))
            sage: p
            A B C D E
            E D B C A
            sage: p.stratum_component() == cc
            True

            sage: cc = Stratum({3:2,0:6}, k=1).hyperelliptic_component()
            sage: p = cc.single_cylinder_representative()
            sage: p
            0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
            14 12 13 10 11 8 9 7 5 6 3 4 1 2 0
            sage: p.stratum_component() == cc
            True

            sage: cc = Stratum([2], k=1).hyperelliptic_component()
            sage: cc.single_cylinder_representative()
            Traceback (most recent call last):
            ...
            ValueError: no 1,1-square-tiled surfaces in this connected component try again with H_2(2, 0)^hyp
            sage: cc = Stratum({3:2,0:5}, k=1).hyperelliptic_component()
            sage: cc.single_cylinder_representative()
            Traceback (most recent call last):
            ...
            ValueError: no 1,1-square-tiled surfaces in this connected component try again with H_4(3^2, 0^6)^hyp
        """
        stratum = self.stratum()
        genus = stratum.surface_genus()
        nb_fk_zeros = sum(m == 0 for m in stratum.signature())
        nb_real_zeros = sum(m != 0 for m in stratum.signature())
        add_fk_zeros = nb_fk_zeros - 2 * genus + 4 - nb_real_zeros

        from surface_dynamics.interval_exchanges.constructors import GeneralizedPermutation

        if nb_real_zeros == 1 and add_fk_zeros < 0:
            raise ValueError("no 1,1-square-tiled surfaces in this connected component try again with %s^hyp" %(str(Stratum({2*genus-2:1,0:2*genus-3}, k=1))))
        elif nb_real_zeros == 2 and add_fk_zeros < 0:
            raise ValueError("no 1,1-square-tiled surfaces in this connected component try again with %s^hyp" %(str(Stratum({genus-1:2,0:2*genus-2}, k=1))))
        elif not nb_real_zeros:
            from surface_dynamics.flat_surfaces.single_cylinder import cylinder_concatenation
            fk_zeros_perm = GeneralizedPermutation([0],[0])
            mk_pt_perm = GeneralizedPermutation([0,1],[1,0])
            for i in range(nb_fk_zeros):
                fk_zeros_perm = cylinder_concatenation(fk_zeros_perm, mk_pt_perm)
            if alphabet is not None:
                fk_zeros_perm.alphabet(alphabet)

            return fk_zeros_perm.reduced() if reduced else fk_zeros_perm

        else:
            top_row = list(range(0, 4*genus-3+2*(nb_real_zeros-1)+add_fk_zeros))
            bot_row = [4*genus-4+2*(nb_real_zeros-1)+add_fk_zeros]
            for i in range(4*genus-6+2*(nb_real_zeros-1)+add_fk_zeros,2*genus-2+add_fk_zeros,-2):
                bot_row.append(i)
                bot_row.append(i+1)
            bot_row.extend(2*genus-1+i for i in range(add_fk_zeros+1))
            for i in range(2*genus-3,-1,-2):
                bot_row.append(i)
                bot_row.append(i+1)
            bot_row.append(0)

            perm = GeneralizedPermutation(top_row, bot_row, reduced=reduced)
            if alphabet is not None:
                perm.alphabet(alphabet)
            return perm


HypASC = HypAbelianStratumComponent

class NonHypAbelianStratumComponent(ASC):
    """
    Non hyperelliptic component of Abelian stratum.
    """
    _name = 'nonhyp'

    def rauzy_class_cardinality(self, left_degree=None, reduced=True):
        r"""
        Return the cardinality of Rauzy diagram associated to this non
        hyperelliptic component.

        INPUT:

        - ``left_degree`` - integer

        - ``reduced`` - boolean (default: True)

        EXAMPLES::

            sage: from surface_dynamics import Stratum

        Examples in genus 3::

            sage: c = Stratum([3,3], k=1).non_hyperelliptic_component()
            sage: c.rauzy_class_cardinality()
            15568

            sage: c = Stratum([3,3,0], k=1).non_hyperelliptic_component()
            sage: c.rauzy_class_cardinality()
            173723

            sage: c.rauzy_diagram(left_degree=3)  # long time
            Rauzy diagram with 155680 permutations
            sage: c.rauzy_class_cardinality(left_degree=3)
            155680
            sage: c.rauzy_diagram(left_degree=3, reduced=False)  # not tested
            Rauzy diagram with 311360 permutations
            sage: c.rauzy_class_cardinality(left_degree=3, reduced=False)
            311360

            sage: c.rauzy_diagram(left_degree=0)  # long time
            Rauzy diagram with 18043 permutations
            sage: c.rauzy_class_cardinality(left_degree=0)
            18043
            sage: cc.rauzy_diagram(left_degree=0, reduced=False) # not tested
            Rauzy diagram with 288688 permutations
            sage: c.rauzy_class_cardinality(left_degree=0,reduced=False)
            288688

        When genus growths, the size of the Rauzy diagram becomes very big::

            sage: c = Stratum([5,5], k=1).non_hyperelliptic_component()
            sage: c.rauzy_class_cardinality()
            136116680

            sage: c = Stratum([7,7,0], k=1).non_hyperelliptic_component()
            sage: c.rauzy_class_cardinality()
            88484743236111
            sage: c.rauzy_class_cardinality(left_degree=7, reduced=False)
            334071852804864
        """
        import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rdc

        profile = list(map(lambda x: x+1,self.stratum().signature()))
        hyp = self.stratum().hyperelliptic_component()

        if left_degree is not None:
            assert isinstance(left_degree, (int,Integer)), "if not None, left_degree should be an integer"
            left_degree = int(left_degree) + 1
            assert left_degree in profile, "if not None, the degree should be one of the degree of the stratum"

            if reduced:
                n = Integer(rdc.gamma_irr(profile,left_degree))
                n_hyp = hyp.rauzy_class_cardinality(left_degree-1)

            else:
                return Rational((1,2)) * (Partition(profile).centralizer_size() /
                    ((left_degree) * profile.count(left_degree)) *
                    self.rauzy_class_cardinality(left_degree-1,reduced=True))

        elif reduced is True:
            n = Integer(rdc.gamma_irr(profile))
            n_hyp = hyp.rauzy_class_cardinality()

        else:
            raise NotImplementedError("no formula known for cardinality of  extended labeled Rauzy classes")

        return n - n_hyp

    def standard_permutations_number(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: C = Stratum([3,3], k=1).non_hyperelliptic_component()
            sage: len(C.standard_permutations())  # long time
            275
            sage: C.standard_permutations_number()
            275

            sage: C = Stratum([5,5], k=1).non_hyperelliptic_component()
            sage: C.standard_permutations_number()
            1022399

            sage: C = Stratum([7,7], k=1).non_hyperelliptic_component()
            sage: C.standard_permutations_number()
            19229011199
        """
        import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rdc

        profile = list(map(lambda x: x+1, self.stratum().signature()))
        return rdc.number_of_standard_permutations(profile) - self.stratum().hyperelliptic_component().standard_permutations_number()

    def _cylinder_diagram_iterator(self, ncyls=None):
        r"""
        Return the list of cylinder diagrams (or completely periodic
        configurations) associated to this non-hyperelliptic component.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: cc = Stratum([3,3], k=1).non_hyperelliptic_component()
            sage: it = cc.cylinder_diagram_iterator()
            sage: c0 = next(it); c0
            (0,1,4,6,2,5,3,7)-(0,1,4,5,3,6,2,7)
            sage: c0.stratum_component()
            H_4(3^2)^nonhyp

            sage: it = cc.cylinder_diagram_iterator(4, force_computation=True)
            sage: c0 = next(it); c0  # random
            (0,4,2)-(0,6) (1,5,3)-(1,7) (6)-(4,5) (7)-(2,3)
            sage: c0.stratum_component()
            H_4(3^2)^nonhyp
            sage: c0.ncyls()
            4
            sage: all(cd.stratum_component() == cc and cd.ncyls() == 4 for cd in it)
            True
            sage: sum(1 for _ in cc.cylinder_diagram_iterator(4, force_computation=True))
            184
        """
        return filter(lambda c: not c.is_hyperelliptic(),
                self.stratum().cylinder_diagram_iterator(ncyls, True, True))

NonHypASC = NonHypAbelianStratumComponent

class EvenAbelianStratumComponent(ASC):
    """
    Connected component of Abelian stratum with even spin structure.

    .. warning::

        Internal class! Do not use directly!
    """
    _name = 'even'

    def spin(self):
        r"""
        Return ``0``.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: c = Stratum([4,2], k=1).even_component(); c
            H_4(4, 2)^even
            sage: c.spin()
            0
        """
        return Integer(0)

    def permutation_representative(self, left_degree=None, reduced=True, alphabet=None, relabel=True):
        r"""
        Returns the Zorich representative of this connected component.

        Zorich constructs explicitly interval exchange
        transformations for each stratum in [Zor08]_.

        INPUT:

        - ``reduced`` - boolean (default: True): whether you obtain a reduced or
          labelled permutation

        - ``left_degree`` - integer (optional) - a specified degree of zero at
          the left of the interval.

        - ``alphabet`` - alphabet or None (default: None): whether you want to
          specify an alphabet for your representative

        - ``relabel`` - boolean (default: True) - if False uses Zorich's natural
          numbering otherwise uses 0,1,...

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: c = Stratum([6], k=1).even_component()
            sage: c
            H_4(6)^even
            sage: p = c.permutation_representative(alphabet=range(8))
            sage: p
            0 1 2 3 4 5 6 7
            5 4 3 2 7 6 1 0
            sage: p.stratum_component()
            H_4(6)^even

        ::

            sage: c = Stratum([4,4], k=1).even_component()
            sage: c
            H_5(4^2)^even
            sage: p = c.permutation_representative(alphabet=range(11))
            sage: p
            0 1 2 3 4 5 6 7 8 9 10
            5 4 3 2 6 8 7 10 9 1 0
            sage: p.stratum_component()
            H_5(4^2)^even

        Different markings lead to different Rauzy diagrams::

            sage: c = Stratum([4,2,0], k=1).even_component()
            sage: p = c.permutation_representative(left_degree=4); p
            0 1 2 3 4 5 6 7 8 9
            6 5 4 3 7 9 8 2 0 1
            sage: p.stratum_component()
            H_4(4, 2, 0)^even
            sage: p.marking().left()
            5
            sage: p.rauzy_diagram()   # long time
            Rauzy diagram with 66140 permutations

            sage: p = c.permutation_representative(left_degree=2); p
            0 1 2 3 4 5 6 7 8 9
            7 6 5 4 3 9 8 2 0 1
            sage: p.stratum_component()
            H_4(4, 2, 0)^even
            sage: p.marking().left()
            3
            sage: p.rauzy_diagram()   # long time
            Rauzy diagram with 39540 permutations

            sage: p = c.permutation_representative(left_degree=0); p
            0 1 2 3 4 5 6 7 8 9
            6 4 3 2 7 9 8 1 5 0
            sage: p.stratum_component()
            H_4(4, 2, 0)^even
            sage: p.marking().left()
            1
            sage: p.rauzy_diagram()   # long time
            Rauzy diagram with 11792 permutations
        """
        z = list(m for m in self._stratum.signature() if m)
        n = self._stratum.signature().count(0)
        g = self._stratum.surface_genus()

        if left_degree is not None:
            if not isinstance(left_degree, (int,Integer)):
                raise ValueError("left_degree (=%d) should be one of the degree"%left_degree)
            if left_degree == 0:
                if n == 0:
                    raise ValueError("left_degree (=%d) should be one of the degree" % left_degree)
            elif left_degree not in z:
                raise ValueError("left_degree (=%d) should be one of the degree" % left_degree)
            else:
                z.remove(left_degree)
                z.insert(0,left_degree)

        l0 = list(range(3*g-2))
        l1 = [6, 5, 4, 3, 2, 7, 9, 8]
        for k in range(10, 3*g-4, 3):
            l1 += [k, k+2, k+1]
        l1 += [1, 0]

        k = 4
        for d in z:
            for i in range(d//2 - 1):
                del l0[l0.index(k)]
                del l1[l1.index(k)]
                k += 3
            k += 3

        # if there are marked points we transform 0 in [3g-2, 3g-3, ...]
        if n != 0:
            interval = list(range(3*g-2, 3*g - 2 + n))

            if left_degree == 0:
                k = l0.index(6)
                l0[k:k] = interval
                l1[-1:-1] = interval
            else:
                l0[1:1] = interval
                l1.extend(interval)

        if reduced:
            from surface_dynamics.interval_exchanges.reduced import ReducedPermutationIET
            p = ReducedPermutationIET([l0, l1])

        else:
            from surface_dynamics.interval_exchanges.labelled import LabelledPermutationIET
            p = LabelledPermutationIET([l0, l1])

        if alphabet is not None:
            p.alphabet(alphabet)
        elif relabel:
            p.alphabet(range(len(p)))
        return p

    def rauzy_class_cardinality(self, left_degree=None, reduced=True):
        r"""
        Cardinality of rauzy diagram for even component of a stratum

        INPUT:

        - ``left_degree`` - integer

        - ``reduced`` - boolean


        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: c = Stratum([6], k=1).even_component()
            sage: c.rauzy_diagram()
            Rauzy diagram with 2327 permutations
            sage: c.rauzy_class_cardinality()
            2327

            sage: c = Stratum([4,2,0], k=1).even_component()
            sage: c.rauzy_class_cardinality()
            117472

            sage: c.rauzy_diagram(left_degree=4)  # long time
            Rauzy diagram with 66140 permutations
            sage: c.rauzy_class_cardinality(left_degree=4)
            66140
            sage: c.rauzy_diagram(left_degree=4, reduced=False)  # long time
            Rauzy diagram with 198420 permutations
            sage: c.rauzy_class_cardinality(left_degree=4,reduced=False)
            198420

            sage: c.rauzy_class_cardinality(2)
            39540
            sage: c.rauzy_diagram(left_degree=2)  # long time
            Rauzy diagram with 39540 permutations
            sage: c.rauzy_diagram(left_degree=2, reduced=False)  # long time
            Rauzy diagram with 197700 permutations
            sage: c.rauzy_class_cardinality(left_degree=2, reduced=False)
            197700

            sage: c.rauzy_class_cardinality(0)
            11792
            sage: c.rauzy_diagram(left_degree=0)
            Rauzy diagram with 11792 permutations
            sage: c.rauzy_diagram(left_degree=0, reduced=False)  # long time
            Rauzy diagram with 176880 permutations
            sage: c.rauzy_class_cardinality(left_degree=0, reduced=False)
            176880
        """
        import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rdc

        profile = list(map(lambda x: x+1, self.stratum().signature()))
        if left_degree is not None:
            assert isinstance(left_degree, (int,Integer)), "if not None, left_degree should be an integer"
            left_degree = int(left_degree) + 1
            assert left_degree in profile, "if not None, the degree should be one of the degree of the stratum"

            if reduced is False:
                return (Partition(profile).centralizer_size() /
                        (left_degree * profile.count(left_degree)) *
                        self.rauzy_class_cardinality(left_degree-1, reduced=True))

        elif reduced is False:
            raise NotImplementedError("no formula known for extended labeled Rauzy classes")

        N = Integer(rdc.gamma_irr(profile,left_degree) - rdc.delta_irr(profile,left_degree))/2

        if (self.stratum().number_of_components() == 3 and
            self.stratum().hyperelliptic_component().spin() == 0):
            if left_degree is None:
                hyp_card = self.stratum().hyperelliptic_component().rauzy_class_cardinality()
            else:
                hyp_card = self.stratum().hyperelliptic_component().rauzy_class_cardinality(left_degree-1)

            return N - hyp_card

        return N

    def standard_permutations_number(self):
        r"""
        Return the number of standard permutation of this even component.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

        For strata in genus 3, the number of standard permutations is reasonably
        small and the whole set can be computed::

            sage: C = Stratum([6], k=1).even_component()
            sage: len(C.standard_permutations())  # long time
            44
            sage: C.standard_permutations_number()
            44

            sage: C = Stratum([4,2], k=1).even_component()
            sage: len(C.standard_permutations())   # long time
            136
            sage: C.standard_permutations_number()
            136

            sage: C = Stratum([2,2,2], k=1).even_component()
            sage: len(C.standard_permutations())   # long time
            92
            sage: C.standard_permutations_number()
            92

        For higher genera, this number can be very big::

            sage: C = Stratum([20], k=1).even_component()
            sage: C.standard_permutations_number()
            109398514483439999
        """
        import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rdc

        profile = [x+1 for x in self.stratum().signature()]
        N = Integer(rdc.gamma_std(profile) - rdc.delta_std(profile)) / 2

        if (self.stratum().number_of_components() == 3 and
            self.stratum().hyperelliptic_component().spin() == 0):
            return N - 1

        return N

    def _cylinder_diagram_iterator(self, ncyls=None):
        r"""
        Iterator over cylinder diagram.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: cc = Stratum([4,2], k=1).even_component()
            sage: it = cc.cylinder_diagram_iterator(4)
            sage: next(it).stratum_component()
            H_4(4, 2)^even

            sage: it = cc.cylinder_diagram_iterator(3, force_computation=True)
            sage: next(it).stratum_component()
            H_4(4, 2)^even
        """
        if self.stratum().has_hyperelliptic_component():
            return filter(
                    lambda c: not c.is_hyperelliptic() and c.spin_parity() == 0,
                    self.stratum().cylinder_diagram_iterator(ncyls,True,True))

        return filter(lambda c: c.spin_parity() == 0,
                self.stratum().cylinder_diagram_iterator(ncyls,True,True))

    def single_cylinder_representative(self, alphabet=None, reduced=False):
        r"""
        Returns a single cylinder permutation representative.

        Returns a permutation representative of a square-tiled surface in this
        component having a single vertical cylinder and a single horizontal cylinder.

        Such representatives were constructed for every stratum of Abelian
        differentials by Jeffreys [Jef19]_.

        INPUT:

        - ``alphabet`` -- an optional alphabet for the permutation representative

        - ``reduced`` (boolean, default ``True``) -- whether to return a reduced
          permutation (ie without labels)

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: cc = Stratum([6], k=1).even_component()
            sage: p = cc.single_cylinder_representative(alphabet=Alphabet(name='lower'))
            sage: p
            a b c d e f g h
            c h g f d b e a
            sage: p.stratum_component() == cc
            True

            sage: cc = Stratum([4,4], k=1).even_component()
            sage: p = cc.single_cylinder_representative()
            sage: p
            0 1 2 3 4 5 6 7 8 9 10
            2 10 7 5 8 1 9 6 4 3 0
            sage: p.stratum_component() == cc
            True
        """
        from surface_dynamics.flat_surfaces.single_cylinder import cylinder_concatenation
        from surface_dynamics.flat_surfaces.single_cylinder import (no_two_even,
                one_two_even, two_twos_even, even_twos_even, odd_twos_even)
        from surface_dynamics.interval_exchanges.constructors import GeneralizedPermutation

        zeros = self.stratum().signature()
        real_zeros = [z for z in zeros if z != 0]

        fk_zeros_perm = GeneralizedPermutation([0],[0])
        mk_pt_perm = GeneralizedPermutation([0,1],[1,0])
        for i in range(self.stratum().signature().count(0)):
            fk_zeros_perm = cylinder_concatenation(fk_zeros_perm,mk_pt_perm)

        two_count = real_zeros.count(2)
        if two_count == 0:
            perm = cylinder_concatenation(fk_zeros_perm,no_two_even(real_zeros))
        elif two_count == 1:
            perm = cylinder_concatenation(fk_zeros_perm,one_two_even(real_zeros))
        elif two_count == 2:
            perm = cylinder_concatenation(fk_zeros_perm,two_twos_even(real_zeros))
        elif two_count > 2 and two_count%2 == 0:
            perm = cylinder_concatenation(fk_zeros_perm,even_twos_even(real_zeros,two_count))
        else:
            perm = cylinder_concatenation(fk_zeros_perm,odd_twos_even(real_zeros,two_count))

        if alphabet is not None:
            perm.alphabet(alphabet)
        return perm.reduced() if reduced else perm

EvenASC = EvenAbelianStratumComponent


class OddAbelianStratumComponent(ASC):
    r"""
    Connected component of an Abelian stratum with odd spin parity.
    """
    _name = 'odd'

    def spin(self):
        r"""
        Returns 1 which is, by definition, the spin parity of this stratum component.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: c = Stratum([4], k=1).odd_component(); c
            H_3(4)^odd
            sage: c.spin()
            1
        """
        return 1

    def permutation_representative(self, left_degree=None, reduced=True, alphabet=None, relabel=True):
        """
        Returns the Zorich representative of this connected component.

        A. Zorich constructs explicitly interval exchange
        transformations for each stratum in [Zor08]_.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: a = Stratum([6], k=1).odd_component()
            sage: p = a.permutation_representative()
            sage: p
            0 1 2 3 4 5 6 7
            3 2 5 4 7 6 1 0
            sage: p.stratum_component()
            H_4(6)^odd

        ::

            sage: a = Stratum([4,4], k=1).odd_component()
            sage: p = a.permutation_representative()
            sage: p
            0 1 2 3 4 5 6 7 8 9 10
            3 2 5 4 6 8 7 10 9 1 0
            sage: p.stratum_component()
            H_5(4^2)^odd

        Different markings lead to different Rauzy diagrams::

            sage: c = Stratum([4,2,0], k=1).odd_component()
            sage: p = c.permutation_representative(left_degree=4); p
            0 1 2 3 4 5 6 7 8 9
            4 3 6 5 7 9 8 2 0 1
            sage: p.stratum_component()
            H_4(4, 2, 0)^odd
            sage: p.marking().left()
            5
            sage: p.rauzy_diagram()   # not tested
            Rauzy diagram with 147090 permutations

            sage: p = c.permutation_representative(left_degree=2); p
            0 1 2 3 4 5 6 7 8 9
            4 3 5 7 6 9 8 2 0 1
            sage: p.stratum_component()
            H_4(4, 2, 0)^odd
            sage: p.marking().left()
            3
            sage: p.rauzy_diagram()   # long time
            Rauzy diagram with 87970 permutations

            sage: p = c.permutation_representative(left_degree=0); p
            0 1 2 3 4 5 6 7 8 9
            4 2 6 5 7 9 8 1 3 0
            sage: p.stratum_component()
            H_4(4, 2, 0)^odd
            sage: p.marking().left()
            1
            sage: p.rauzy_diagram()   # long time
            Rauzy diagram with 27754 permutations
        """
        zeros = list(m for m in self.stratum().signature() if m)
        n = self.stratum().signature().count(0)
        g = self.stratum().surface_genus()

        if left_degree is not None:
            if not isinstance(left_degree, (int,Integer)):
                raise ValueError("left_degree (=%d) should be one of the degree" % left_degree)
            if left_degree == 0:
                if n == 0:
                    raise ValueError("left_degree (=%d) should be one of the degree" % left_degree)
            elif left_degree not in zeros:
                raise ValueError("left_degree (=%d) should be one of the degree" % left_degree)
            else:
                zeros.remove(left_degree)
                zeros.insert(0,left_degree)

        z = [x//2 for x in zeros]

        l0 = list(range(3*g-2))
        l1 = [3, 2]
        for k in range(4, 3*g-4, 3):
            l1 += [k, k+2, k+1]
        l1 += [1, 0]

        k = 4
        for d in z:
            for i in range(d-1):
                del l0[l0.index(k)]
                del l1[l1.index(k)]
                k += 3
            k += 3

        # marked points
        if n != 0:
            interval = list(range(3*g-2, 3*g-2+n))

            if left_degree == 0:
                k = l0.index(3)
                l0[k:k] = interval
                l1[-1:-1] = interval
            else:
                l0[1:1] = interval
                l1.extend(interval)

        if reduced:
            from surface_dynamics.interval_exchanges.reduced import ReducedPermutationIET
            p = ReducedPermutationIET([l0, l1])

        else:
            from surface_dynamics.interval_exchanges.labelled import LabelledPermutationIET
            p = LabelledPermutationIET([l0, l1])

        if alphabet is not None:
            p.alphabet(alphabet)
        elif relabel:
            p.alphabet(range(len(p)))
        return p

    def rauzy_class_cardinality(self, left_degree=None, reduced=True):
        r"""
        Cardinality of rauzy diagram for odd component

        INPUT:

        - ``left_degree`` - integer (optional)

        - ``reduced`` - boolean (default: True)

        EXAMPLES::

            sage: from surface_dynamics import Stratum

        The genus must be at least 3 to have an odd component::

            sage: c = Stratum([4], k=1).odd_component()
            sage: c.rauzy_diagram()
            Rauzy diagram with 134 permutations
            sage: c.rauzy_class_cardinality()
            134
            sage: c = Stratum([4,0], k=1).odd_component()
            sage: c.rauzy_diagram()
            Rauzy diagram with 1114 permutations
            sage: c.rauzy_class_cardinality()
            1114

            sage: c = Stratum([2,2], k=1).odd_component()
            sage: c.rauzy_diagram()
            Rauzy diagram with 294 permutations
            sage: c.rauzy_class_cardinality()
            294

            sage: c = Stratum([2,2,0], k=1).odd_component()
            sage: c.rauzy_class_cardinality()
            2723

            sage: c.rauzy_diagram(left_degree=2)
            Rauzy diagram with 2352 permutations
            sage: c.rauzy_class_cardinality(left_degree=2)
            2352
            sage: c.rauzy_diagram(left_degree=2, reduced=False)
            Rauzy diagram with 7056 permutations
            sage: c.rauzy_class_cardinality(left_degree=2, reduced=False)
            7056

            sage: c.rauzy_diagram(left_degree=0)
            Rauzy diagram with 371 permutations
            sage: c.rauzy_class_cardinality(left_degree=0)
            371
            sage: c.rauzy_diagram(left_degree=0, reduced=False)
            Rauzy diagram with 6678 permutations
            sage: c.rauzy_class_cardinality(left_degree=0, reduced=False)
            6678


        Example in higher genus for which an explicit computation of the Rauzy
        diagram would be very long::

            sage: c = Stratum([4,2,0], k=1).odd_component()
            sage: c.rauzy_class_cardinality()
            262814
            sage: c = Stratum([4,4,4], k=1).odd_component()
            sage: c.rauzy_class_cardinality()
            24691288838
            sage: c.rauzy_class_cardinality(left_degree=4, reduced=False)
            1234564441900
        """
        import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rdc

        profile = list(map(lambda x: x+1, self.stratum().signature()))
        if left_degree is not None:
            assert isinstance(left_degree, (int,Integer)), "if not None, left_degree should be an integer"
            left_degree = int(left_degree) + 1
            assert left_degree in profile, "if not None, the degree should be one of the degree of the stratum"

            if reduced is False:
                return (Partition(profile).centralizer_size() /
                        (left_degree * profile.count(left_degree)) *
                        self.rauzy_class_cardinality(left_degree-1, reduced=True))

        elif reduced is False:
            raise NotImplementedError("no formula known for labeled extended Rauzy classes")

        N = sum(rdc.gamma_irr(profile,left_degree) + rdc.delta_irr(profile,left_degree))//2

        if (self.stratum().number_of_components() == 3 and
            self.stratum().hyperelliptic_component().spin() == 1):
            return N - self.stratum().hyperelliptic_component().rauzy_class_cardinality()

        return N

    def standard_permutations_number(self):
        r"""
        Return the number of standard permutation of this even component.


        EXAMPLES::

            sage: from surface_dynamics import Stratum

        In genus 2, there are two strata which contains an odd component::

            sage: C = Stratum([4], k=1).odd_component()
            sage: len(C.standard_permutations())
            7
            sage: C.standard_permutations_number()
            7

            sage: C = Stratum([2,2], k=1).odd_component()
            sage: len(C.standard_permutations())
            11
            sage: C.standard_permutations_number()
            11

        In genus 3, the number of standard permutations is reasonably small and
        the whole set can be computed::

            sage: C = Stratum([6], k=1).odd_component()
            sage: len(C.standard_permutations())   # long time
            135
            sage: C.standard_permutations_number()
            135

            sage: C = Stratum([4,2], k=1).odd_component()
            sage: len(C.standard_permutations())   # long time
            472
            sage: C.standard_permutations_number()
            472

            sage: C = Stratum([2,2,2], k=1).odd_component()
            sage: len(C.standard_permutations())   # long time
            372
            sage: C.standard_permutations_number()
            372

        For higher genera, this number can be very big::

            sage: C = Stratum([8,6,4,2], k=1).odd_component()
            sage: C.standard_permutations_number()
            26596699869748377600
        """
        import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rdc

        profile = [x+1 for x in self.stratum().signature()]
        N = Integer(rdc.gamma_std(profile) + rdc.delta_std(profile)) / 2

        if (self.stratum().number_of_components() == 3 and
            self.stratum().hyperelliptic_component().spin() == 1):
            return N - 1

        return N

    def _cylinder_diagram_iterator(self, ncyls=None):
        r"""
        Return the list of cylinder diagrams associated to this odd connected
        component.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: C = Stratum([4], k=1).odd_component()
            sage: for c in C.cylinder_diagrams(1): print(c)
            (0,2,1,4,3)-(0,4,2,1,3)
            (0,4,1,2,3)-(0,1,3,4,2)
            sage: for c in C.cylinder_diagrams(2): print(c)
            (0,1,2)-(0,3,1,4) (3,4)-(2)
            (0,1,3)-(4) (2,4)-(0,1,2,3)
            (0,2,3)-(2,4) (1,4)-(0,1,3)
            (0,2,3,1)-(0,2,1,4) (4)-(3)
            sage: for c in C.cylinder_diagrams(3): print(c)
            (0,1)-(0,3,4) (2,3)-(1) (4)-(2)
            (0,2)-(0,3) (1,3)-(1,4) (4)-(2)
            (0,2)-(4) (1,4)-(2,3) (3)-(0,1)
            (0,2,1)-(0,3,4) (3)-(1) (4)-(2)
        """
        if self.stratum().has_hyperelliptic_component():
            return filter(
                    lambda c: not c.is_hyperelliptic() and c.spin_parity() == 1,
                    self.stratum().cylinder_diagram_iterator(ncyls,True,True))

        return filter(lambda c: c.spin_parity == 1,
                self.stratum().cylinder_diagram_iterator(ncyls,True,True))

    def single_cylinder_representative(self, alphabet=None, reduced=True):
        r"""
        Returns a single cylinder permutation representative.

        Returns a permutation representative of a square-tiled surface in this
        component having a single vertical cylinder and a single horizontal cylinder.

        Such representatives were constructed for every stratum of Abelian
        differentials by Jeffreys [Jef19]_.

        INPUT:

        - ``alphabet`` -- an optional alphabet for the permutation representative

        - ``reduced`` (boolean, default ``True``) -- whether to return a reduced
          permutation (ie without labels)

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: cc = Stratum([4], k=1).odd_component()
            sage: p = cc.single_cylinder_representative(alphabet=Alphabet(name='upper'))
            sage: p
            A B C D E F
            C F E B D A
            sage: p.stratum_component() == cc
            True

            sage: cc = Stratum([6,2], k=1).odd_component()
            sage: p = cc.single_cylinder_representative()
            sage: p
            0 1 2 3 4 5 6 7 8 9 10
            2 5 4 6 3 8 10 7 1 9 0
            sage: p.stratum_component() == cc
            True
        """
        from surface_dynamics.flat_surfaces.single_cylinder import (cylinder_concatenation,
                no_two_odd, one_two_odd, even_twos_odd, odd_twos_odd)
        from surface_dynamics.interval_exchanges.constructors import GeneralizedPermutation

        zeros = self.stratum().signature()
        real_zeros = [z for z in zeros if z != 0]

        fk_zeros_perm = GeneralizedPermutation([0],[0])
        mk_pt_perm = GeneralizedPermutation([0,1],[1,0])
        for i in range(self.stratum().signature().count(0)):
            fk_zeros_perm = cylinder_concatenation(fk_zeros_perm,mk_pt_perm)

        two_count = real_zeros.count(2)
        if two_count == 0:
            perm = cylinder_concatenation(fk_zeros_perm,no_two_odd(real_zeros))
        elif two_count == 1:
            perm = cylinder_concatenation(fk_zeros_perm,one_two_odd(real_zeros))
        elif two_count >= 2 and two_count % 2 == 0:
            perm = cylinder_concatenation(fk_zeros_perm,even_twos_odd(real_zeros,two_count))
        else:
            perm = cylinder_concatenation(fk_zeros_perm,odd_twos_odd(real_zeros,two_count))

        if alphabet is not None:
            perm.alphabet(alphabet)
        return perm.reduced() if reduced else perm

OddASC = OddAbelianStratumComponent


#
# iterators for Abelian strata with constraints on genus and dimension
#

class AbelianStrata(Strata):
    r"""
    Abelian strata.

    INPUT:

    - ``genus`` - a non negative integer or None

    - ``dimension`` - a non negative integer or None

    - ``fake_zeros`` - boolean

    EXAMPLES::

        sage: from surface_dynamics import *

    Abelian strata with a given genus::

        sage: for s in AbelianStrata(genus=1): print(s)
        H_1(0)

    ::

        sage: for s in AbelianStrata(genus=2): print(s)
        H_2(2)
        H_2(1^2)

    ::

        sage: for s in AbelianStrata(genus=3): print(s)
        H_3(4)
        H_3(3, 1)
        H_3(2^2)
        H_3(2, 1^2)
        H_3(1^4)

    ::

        sage: for s in AbelianStrata(genus=4): print(s)
        H_4(6)
        H_4(5, 1)
        H_4(4, 2)
        H_4(4, 1^2)
        H_4(3^2)
        H_4(3, 2, 1)
        H_4(3, 1^3)
        H_4(2^3)
        H_4(2^2, 1^2)
        H_4(2, 1^4)
        H_4(1^6)

    Get outside of the tests.
    Abelian strata with a given number of intervals

    sage for s in AbelianStrata(dimension=2): print(s)
    H^out([0])

    sage for s in AbelianStrata(dimension=3): print(s)
    H^out([0], 0)

    sage for s in AbelianStrata(dimension=4): print(s)
    H^out([2])
    H^out([0], 0, 0)

    Get outside of tests
    sage  for s in AbelianStrata(dimension=5): print(s)
    H^out(2, [0])
    H^out([2], 0)
    H^out([1], 1)
    H^out([0], 0, 0, 0)
    """
    def __new__(self, genus=None, dimension=None, fake_zeros=None):
        fake_zeros = bool(fake_zeros)
        if dimension is not None:
            dimension = Integer(dimension)
            if dimension < 0:
                raise ValueError("dimension must be a non-negative integer")

        if genus is not None:
            genus = Integer(genus)
            if genus < 0:
                raise ValueError("genus must be a non-negative integer")

        if genus is None:
            if dimension is None:
                cls = AbelianStrata_all
            else:
                cls = AbelianStrata_d
        elif dimension is None:
            cls = AbelianStrata_g
        else:
            cls = AbelianStrata_gd

        S = Strata.__new__(cls, genus, dimension, fake_zeros)
        AbelianStrata.__init__(S, genus, dimension, fake_zeros)
        return S

    def __init__(self, genus=None, dimension=None, fake_zeros=None):
        r"""
        TESTS::

            sage: from surface_dynamics import *

            sage: s = AbelianStrata(genus=3)
            sage: loads(dumps(s)) == s
            True
        """
        if dimension is None and (fake_zeros is True or genus is None):
            category = InfiniteEnumeratedSets()
        else:
            category = FiniteEnumeratedSets()
        Parent.__init__(self, category=category, facade=True)

        self._genus = genus
        self._dimension = dimension
        self._fake_zeros = fake_zeros

    def __eq__(self, other):
        r"""
        Equality test.
        """
        return (isinstance(other, AbelianStrata) and
                (self._dimension == other._dimension) and
                (self._genus == other._genus) and
                (self._fake_zeros == other._fake_zeros))

    def __ne__(self, other):
        r"""
        Difference test.
        """
        return not self.__eq__(other)

    def _repr_(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *

            sage: AbelianStrata()                        # indirect doctest
            Abelian strata

            sage: AbelianStrata(dimension=2)             # indirect doctest
            Abelian strata of dimension 2

            sage: AbelianStrata(genus=3)                 # indirect doctest
            Abelian strata of genus 3 surfaces

            sage: AbelianStrata(genus=2, dimension=4)    # indirect doctest
            Abelian strata of genus 2 surfaces and dimension 4
        """
        s = "Abelian strata"

        l = []
        if self._genus is not None:
            l.append("genus {} surfaces".format(self._genus))
        if self._dimension is not None:
            l.append("dimension {}".format(self._dimension))

        if l:
            return "Abelian strata of " + " and ".join(l)
        else:
            return "Abelian strata"

    def __reduce__(self):
        r"""
        Pickling support.
        """
        return (AbelianStrata, (self._genus, self._dimension, self._fake_zeros))

    def __contains__(self, c):
        r"""
        Containance test

        TESTS::

            sage: from surface_dynamics import Stratum, AbelianStrata

            sage: a = AbelianStrata(genus=3)
            sage: all(s in a for s in a)
            True

            sage: a = AbelianStrata(genus=3,fake_zeros=False)
            sage: all(s in a for s in a)
            True

            sage: a = AbelianStrata(dimension=7,fake_zeros=True)
            sage: all(s in a for s in a)
            True
            sage: Stratum([2,0,0], k=1) in a
            False

            sage: a = AbelianStrata(dimension=7,fake_zeros=False)
            sage: all(s in a for s in a)
            True
            sage: Stratum([4,0], k=1) in a
            False
        """
        if not isinstance(c, AbelianStratum):
            return False

        return ((self._genus is None or c.surface_genus() == self._genus) and
                (self._dimension is None or c.dimension() == self._dimension) and
                (self._fake_zeros is None or self._fake_zeros or not c.signature().count(0)))


class AbelianStrata_g(AbelianStrata):
    r"""
    Stratas of genus g surfaces without fake zeros.

    INPUT:

    - ``genus`` - a non negative integer

    EXAMPLES::

        sage: from surface_dynamics import *

        sage: AbelianStrata(genus=2).list()
        [H_2(2), H_2(1^2)]
        sage: AbelianStrata(genus=3).list()
        [H_3(4), H_3(3, 1), H_3(2^2), H_3(2, 1^2), H_3(1^4)]
        sage: AbelianStrata(genus=4).random_element() #random
        H_4(4, 2)
    """
    def cardinality(self):
        r"""
        Return the number of abelian strata with a given genus.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: AbelianStrata(genus=1).cardinality()
            1
            sage: AbelianStrata(genus=2).cardinality()
            2
            sage: AbelianStrata(genus=3).cardinality()
            5
            sage: AbelianStrata(genus=4).cardinality()
            11
        """
        if self._genus == 0:
            return Integer(0)
        if self._genus == 1:
            return Integer(1)
        return Partitions(2*self._genus-2).cardinality()

    def __iter__(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *

            sage: list(AbelianStrata(genus=1))
            [H_1(0)]
        """
        if self._genus == 0:
            pass
        elif self._genus == 1:
            yield Stratum([0], k=1)
        else:
            for p in Partitions(2*self._genus-2):
                yield Stratum(p._list, k=1)

    def random_element(self):
        r"""
        Return a random stratum.
        """
        if self._genus == 0:
            raise ValueError("No stratum with that genus")
        if self._genus == 1:
            return Stratum([0], k=1)
        return Stratum(Partitions(2*self._genus - 2).random_element(), k=1)

    def first(self):
        r"""
        Return the first element of this list of strata.

        EXAMPLES::

            sage: from surface_dynamics import AbelianStrata

            sage: AbelianStrata(genus=3).first()
            H_3(4)
            sage: AbelianStrata(genus=4).first()
            H_4(6)
        """
        return Stratum([2*self._genus-2], k=1)

    an_element_ = first

    def last(self):
        r"""
        Return the last element of this list of strata.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: AbelianStrata(genus=4).last()
            H_4(1^6)
            sage: AbelianStrata(genus=5).last()
            H_5(1^8)
        """
        return Stratum({1:2*self._genus-2}, k=1)


class AbelianStrata_d(AbelianStrata):
    r"""
    Strata with prescribed dimension.

    INPUT:

    - ``dimension`` - an integer greater than 1

    - ``fake_zeros`` - boolean (default: False) - allows or not fake zeros

    EXAMPLES::

        sage: from surface_dynamics import *

        sage: for a in AbelianStrata(dimension=5,fake_zeros=True):
        ....:     print(a)
        ....:     print(a.permutation_representative())
        H_2(2, 0)
        0 1 2 3 4
        4 1 3 2 0
        H_2(1^2)
        0 1 2 3 4
        4 3 2 1 0
        H_1(0^4)
        0 1 2 3 4
        4 0 1 2 3
    """
    def first(self):
        r"""
        Returns the first stratum.

        EXAMPLES::

            sage: from surface_dynamics import AbelianStrata

            sage: AbelianStrata(dimension=2).first()
            H_1(0)
            sage: AbelianStrata(dimension=3).first()
            H_1(0^2)
            sage: AbelianStrata(dimension=4).first()
            H_2(2)
        """
        n = self._dimension
        if n%2:
            return Stratum([(n-3)//2,(n-3)//2], k=1)
        return Stratum([n-2], k=1)

    an_element = first

    def last(self):
        r"""
        Return the last stratum.

        EXAMPLES::

            sage: from surface_dynamics import AbelianStrata

            sage: AbelianStrata(dimension=9,fake_zeros=True).last()
            H_1(0^8)
            sage: AbelianStrata(dimension=9,fake_zeros=False).last()
            H_3(1^4)

            sage: AbelianStrata(dimension=10,fake_zeros=True).last()
            H_1(0^9)
            sage: AbelianStrata(dimension=10,fake_zeros=False).last()
            H_4(2^3)
        """
        n = self._dimension
        if self._fake_zeros:
            return Stratum({0:n-1}, k=1)
        else:
            if n == 4:
                return Stratum([2], k=1)
            if n == 5:
                return Stratum([1,1], k=1)
            elif n == 6:
                return Stratum([4], k=1)
            else:
                nn = (n-2)%4
                return Stratum({2:3-nn,1:2*((n-10)//4)+2*nn}, k=1)

    def __iter__(self):
        r"""
        Iterator.

        TESTS::

            sage: from surface_dynamics import *

            sage: for a in AbelianStrata(dimension=4,fake_zeros=True): print(a)
            H_2(2)
            H_1(0^3)
        """
        n = self._dimension
        if n < 2:
            pass
        elif self._fake_zeros:
            for s in range(1+n%2, n, 2):
                for p in Partitions(n-1, length=s):
                    yield Stratum([k-1 for k in p], k=1)
        else:
            if n == 2:
                yield Stratum([0], k=1)
            else:
                for s in range(1+n%2, n, 2):
                    for p in Partitions(n-1,length=s,min_part=2):
                        yield Stratum([k-1 for k in p], k=1)

    def cardinality(self):
        r"""
        Return the number of Abelian strata with given dimension.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: AbelianStrata(dimension=5,fake_zeros=True).cardinality()
            3
            sage: AbelianStrata(dimension=5,fake_zeros=False).cardinality()
            1

            sage: AbelianStrata(dimension=6,fake_zeros=True).cardinality()
            4
            sage: AbelianStrata(dimension=6,fake_zeros=False).cardinality()
            1

            sage: AbelianStrata(dimension=7,fake_zeros=True).cardinality()
            6
            sage: AbelianStrata(dimension=7,fake_zeros=False).cardinality()
            2

            sage: AbelianStrata(dimension=12,fake_zeros=True).cardinality()
            29
            sage: AbelianStrata(dimension=12,fake_zeros=False).cardinality()
            7

        TESTS::

            sage: for d in range(1,15):
            ....:   A = AbelianStrata(dimension=d,fake_zeros=True)
            ....:   assert len(A.list()) == A.cardinality()
            ....:   A = AbelianStrata(dimension=d,fake_zeros=False)
            ....:   assert len(A.list()) == A.cardinality()
        """
        n = self._dimension
        if n < 2:
            return Integer(0)

        if self._fake_zeros:
            return sum(Partitions(n-1,length=s).cardinality() for s in range(1+n%2,n,2))
        if n == 2:
            return Integer(1)
        return sum(Partitions(n-1,length=s,min_part=2).cardinality() for s in range(1+n%2,n,2))

class AbelianStrata_gd(AbelianStrata):
    r"""
    Abelian strata of prescribed genus and number of intervals.

    INPUT:

    - ``genus`` - integer: the genus of the surfaces

    - ``dimension`` - integer: the number of intervals

    - ``fake_zeros`` - boolean: whether or not consider fake zeros
    """
    def __iter__(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *

            sage: AbelianStrata(genus=2, dimension=4).list()
            [H_2(2)]

            sage: AbelianStrata(genus=4, dimension=10, fake_zeros=True).list()
            [H_4(6, 0^2), H_4(5, 1, 0), H_4(4, 2, 0), H_4(4, 1^2), H_4(3^2, 0), H_4(3, 2, 1), H_4(2^3)]
            sage: AbelianStrata(genus=4, dimension=10, fake_zeros=False).list()
            [H_4(4, 1^2), H_4(3, 2, 1), H_4(2^3)]

            sage: AbelianStrata(genus=3, dimension=10, fake_zeros=True).list()
            [H_3(4, 0^4), H_3(3, 1, 0^3), H_3(2^2, 0^3), H_3(2, 1^2, 0^2), H_3(1^4, 0)]
        """
        if self._genus == 0:
            pass
        elif self._genus == 1:
            if self._dimension >= 2 and self._fake_zeros:
                yield Stratum([0]*(self._dimension-1), k=1)
        else:
            s = self._dimension - 2*self._genus + 1
            if self._fake_zeros:
                for p in Partitions(2*self._genus - 2 + s, length=s):
                    yield Stratum([k-1 for k in p], k=1)
            else:
                for p in Partitions(2*self._genus - 2, length=s):
                    yield Stratum(p, k=1)

class AbelianStrata_all(AbelianStrata):
    r"""
    Abelian strata.

    INPUT:

    - ``fake_zeros`` - boolean (default: ``False``)

    EXAMPLES::

        sage: from surface_dynamics import *

        sage: A = AbelianStrata()
        sage: it = iter(A)
        sage: for _ in range(10):
        ....:     print(next(it))
        H_1(0)
        H_2(2)
        H_2(1^2)
        H_3(4)
        H_3(3, 1)
        H_3(2^2)
        H_4(6)
        H_3(2, 1^2)
        H_4(5, 1)
        H_4(4, 2)

        sage: A = AbelianStrata(fake_zeros=True)
        sage: it = iter(A)
        sage: for _ in range(10):
        ....:     print(next(it))
        H_1(0)
        H_1(0^2)
        H_2(2)
        H_1(0^3)
        H_2(2, 0)
        H_2(1^2)
        H_1(0^4)
        H_3(4)
        H_2(2, 0^2)
        H_2(1^2, 0)
    """
    def __iter__(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *

            sage: next(iter(AbelianStrata()))  # indirect doctest
            H_1(0)
        """
        from itertools import count
        for d in count(2):
            yield from AbelianStrata(dimension=d, fake_zeros=self._fake_zeros)
