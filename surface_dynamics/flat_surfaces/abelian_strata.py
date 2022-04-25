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

    sage: a = AbelianStratum(1,1)
    sage: a
    H_2(1^2)
    sage: a.genus()
    2
    sage: a.dimension()
    5

::

    sage: a = AbelianStratum(4,3,2,1)
    sage: a
    H_6(4, 3, 2, 1)
    sage: a.genus()
    6
    sage: a.dimension()
    15

By convention, the degrees are always written in decreasing order::

    sage: a1 = AbelianStratum(4,3,2,1)
    sage: a1
    H_6(4, 3, 2, 1)
    sage: a2 = AbelianStratum(2,3,1,4)
    sage: a2
    H_6(4, 3, 2, 1)
    sage: a1 == a2
    True

It is possible to lis strata and their connected components::

    sage: AbelianStratum(10).components()
    [H_6(10)^hyp, H_6(10)^odd, H_6(10)^even]

Get a list of strata with constraints on genus or on the number of intervals
of a representative::

    sage: AbelianStrata(genus=3).list()
    [H_3(4), H_3(3, 1), H_3(2^2), H_3(2, 1^2), H_3(1^4)]

Obtains the connected components of a stratum::

    sage: a = AbelianStratum(0)
    sage: a.components()
    [H_1(0)^hyp]

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
    ....:       for z in set(a.zeros()):
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

        sage: a = AbelianStratum(2, 2)
        sage: a
        H_3(2^2)
        sage: a.components()
        [H_3(2^2)^hyp, H_3(2^2)^odd]

    Get a permutation representative of a connected component::

        sage: a = AbelianStratum(2,2)
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

    @staticmethod
    def __classcall_private__(self, *l):
        if len(l) == 1:
            try:
                Integer(l[0])
            except TypeError:
                l = l[0]
        elif len(l) == 2 and isinstance(l[0], (tuple,list)) and isinstance(1, numbers.Integral):
            l = tuple(l[0]) + (0,) * l[1]

        if isinstance(l, dict):
            l = sum(([v]*e for v,e in iteritems(l)), [])

        zeros = list(map(Integer, filter(lambda x: x, l)))
        if any(z < 0 for z in zeros):
            raise ValueError("the degrees must be non negative")
        nb_fake_zeros = sum(1 for x in l if not x)
        zeros.sort(reverse=True)
        zeros = tuple(zeros)

        s = sum(zeros)
        if s%2:
            raise ValueError("the sum of the degrees must be even")

        if not zeros and not nb_fake_zeros:
            raise ValueError("there must be at least one zero")

        return UniqueRepresentation.__classcall__(AbelianStratum, zeros, nb_fake_zeros)

    def __init__(self, zeros, nb_fake_zeros):
        """
        TESTS::

            sage: from surface_dynamics import *

            sage: s = AbelianStratum(0)
            sage: s == loads(dumps(s))
            True
            sage: s = AbelianStratum(1,1,1,1)
            sage: s == loads(dumps(s))
            True
        """
        self._zeros = zeros
        self._nb_fake_zeros = nb_fake_zeros

        s = sum(self._zeros)
        if s%2:
            raise ValueError("the sum of the degrees must be even")
        genus = s//2 + 1

        if genus == 1:
            self._cc = (HypASC,)

        elif genus == 2:
            self._cc = (HypASC,)

        elif genus == 3:
            if zeros == (2, 2) or zeros == (4,):
                self._cc = (HypASC, OddASC)
            else:
                self._cc = (ASC,)

        elif len(zeros) == 1:
            # just one zeros [2g-2]
            self._cc = (HypASC, OddASC, EvenASC)

        elif zeros == (genus-1, genus-1):
            # two similar zeros [g-1, g-1]
            if genus % 2 == 0:
                self._cc = (HypASC, NonHypASC)

            else:
                self._cc = (HypASC, OddASC, EvenASC)

        elif all(x%2 == 0 for x in zeros):
            # even zeroes [2 l_1, 2 l_2, ..., 2 l_n]
            self._cc = (OddASC, EvenASC)

        else:
            # connected
            self._cc = (ASC, )

    def zeros(self, fake_zeros=True):
        r"""
        Return the multiplicities of the zeros.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: AbelianStratum([1,2,3]).zeros()
            (3, 2, 1)
            sage: AbelianStratum({2:4}).zeros()
            (2, 2, 2, 2)
        """
        if fake_zeros:
            return self._zeros + (0,)*self._nb_fake_zeros
        return self._zeros

    def nb_zeros(self, fake_zeros=True):
        r"""
        Returns the number of zeros of self.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: AbelianStratum(0).nb_zeros()
            1
            sage: AbelianStratum({2:4,3:2}).nb_zeros()
            6
        """
        if fake_zeros:
            return len(self._zeros) + self._nb_fake_zeros
        return len(self._zeros)

    def nb_fake_zeros(self):
        r"""
        Return the number of fake zeros.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: AbelianStratum(0).nb_fake_zeros()
            1
            sage: AbelianStratum(1,1,0,0).nb_fake_zeros()
            2

            sage: QuadraticStratum(0,4,2,2).nb_fake_zeros()
            1
         """
        return self._nb_fake_zeros

    def genus(self):
        r"""
        Return the genus of the stratum.

        OUTPUT:

        integer -- the genus

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: AbelianStratum(0).genus()
            1
            sage: AbelianStratum(1,1).genus()
            2
            sage: AbelianStratum(3,2,1).genus()
            4
        """
        return Integer(sum(self._zeros)//2+1)

    def dimension(self):
        r"""
        Return the complex dimension of this stratum.

        The dimension is `2g-2+s+1` where `g` is the genus of surfaces in the
        stratum, `s` the number of singularities. The complex dimension of a
        stratum is also the number of intervals of any interval exchange
        transformations associated to the strata.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: AbelianStratum(0).dimension()
            2
            sage: AbelianStratum(0,0).dimension()
            3
            sage: AbelianStratum(2).dimension()
            4
            sage: AbelianStratum(1,1).dimension()
            5

        ::

            sage: a = AbelianStratum(4,3,2,1,0)
            sage: p = a.permutation_representative()
            sage: len(p) == a.dimension()
            True
        """
        return 2 * self.genus() + self.nb_zeros() - 1

    def rank(self):
        r"""
        Return the rank of this manifold (half dimension of the absolute part of the tangent space).

        EXAMPLES::

            sage: from surface_dynamics import AbelianStratum

            sage: AbelianStratum(0,0).rank()
            1
            sage: AbelianStratum(2).rank()
            2
            sage: AbelianStratum(2,0,0).rank()
            2
        """
        return self.genus()

    #
    # Connected component
    #

    def has_odd_component(self):
        r"""
        Test whether this stratum has an odd spin component.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: AbelianStratum(2).has_odd_component()
            False
            sage: AbelianStratum(4).has_odd_component()
            True
            sage: AbelianStratum(4).odd_component()
            H_3(4)^odd
        """
        return all(z%2 == 0 for z in self.zeros()) and self.genus() != 2

    def has_even_component(self):
        r"""
        Test whether this stratum has an even spin component.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: AbelianStratum(2,2).has_even_component()
            False
            sage: AbelianStratum(6).has_even_component()
            True
            sage: AbelianStratum(6).even_component()
            H_4(6)^even
        """
        return all(z%2 == 0 for z in self.zeros()) and self.genus() >= 4

    def has_hyperelliptic_component(self):
        r"""
        Test whether this stratum has an hyperelliptic component.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: AbelianStratum(2,1,1).has_hyperelliptic_component()
            False
            sage: AbelianStratum(2,2).has_hyperelliptic_component()
            True
            sage: AbelianStratum(2,2).hyperelliptic_component()
            H_3(2^2)^hyp
        """
        z = self.zeros()
        return len(z) == 1 or (len(z) == 2 and z[0] == z[1])

    def has_non_hyperelliptic_component(self):
        r"""
        Test whether this stratum has a non-hyperelliptic component.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: AbelianStratum(1,1).has_non_hyperelliptic_component()
            False
            sage: AbelianStratum(3,3).has_non_hyperelliptic_component()
            True
            sage: AbelianStratum(3,3).non_hyperelliptic_component()
            H_4(3^2)^nonhyp
        """
        z = self.zeros()
        return len(z) == 2 and z[0] == z[1] and z[0]%2 == 1 and z[0] > 1

    def odd_component(self):
        r"""
        Return the odd component of self (if any).

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: a = AbelianStratum([2,2]); a
            H_3(2^2)
            sage: a.odd_component()
            H_3(2^2)^odd
        """
        if OddASC in self._cc: return OddASC(self)
        raise ValueError("No odd spin component in this stratum")

    def even_component(self):
        r"""
        Return the even component of self (if any)

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: a = AbelianStratum({2:4}); a
            H_5(2^4)
            sage: a.even_component()
            H_5(2^4)^even
        """
        if EvenASC in self._cc: return EvenASC(self)
        raise ValueError("No even spin component in this stratum")

    def hyperelliptic_component(self):
        r"""
        Return the hyperelliptic component of self (if any)

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: a = AbelianStratum(10); a
            H_6(10)
            sage: a.hyperelliptic_component()
            H_6(10)^hyp
        """
        if HypASC in self._cc: return HypASC(self)
        raise ValueError("No hyperelliptic component in this stratum")

    def non_hyperelliptic_component(self):
        r"""
        Return the non hyperelliptic component of self (if any)

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: a = AbelianStratum(3,3); a
            H_4(3^2)
            sage: a.non_hyperelliptic_component()
            H_4(3^2)^nonhyp
        """
        if NonHypASC in self._cc: return NonHypASC(self)
        raise ValueError("No non hyperelliptic component in this stratum")


    #
    # Quadratic cover
    #

    def orientation_quotients(self,fake_zeros=False):
        r"""
        Return the list of quadratic strata such that their orientation cover
        are contained in this stratum.

        If ``fake_zeros`` (default: False) is True we do care about poles which
        becomes a marked zero.

        EXAMPLES::

            sage: from surface_dynamics import *

        The stratum H(2g-2) has one conic singularities of angle `2(2g-1)pi`. The
        only way a surface in H(2g-2) covers a quadratic differential is that
        the quadratic differential has as unique zeros a conical singularity of
        angle `(2g-1) \pi`. The number of poles may vary and give a collection
        of possibilities::

            sage: AbelianStratum(2).orientation_quotients()
            [Q_0(1, -1^5)]
            sage: AbelianStratum(4).orientation_quotients()
            [Q_1(3, -1^3), Q_0(3, -1^7)]
            sage: AbelianStratum(6).orientation_quotients()
            [Q_2(5, -1), Q_1(5, -1^5), Q_0(5, -1^9)]

        A stratum with two zeros may or may not have orientation quotients::

            sage: AbelianStratum(1,1).orientation_quotients()
            [Q_1(2, -1^2), Q_0(2, -1^6)]
            sage: AbelianStratum(2,2).orientation_quotients()
            [Q_1(1^2, -1^2), Q_0(1^2, -1^6), Q_1(4, -1^4), Q_0(4, -1^8)]
            sage: AbelianStratum(3,1).orientation_quotients()
            []

        To impose that covering of poles are fake zeros, switch option
        ``fake_zeros`` to ``True``::

            sage: AbelianStratum(2,2,0,0).orientation_quotients(fake_zeros=True)
            [Q_1(1^2, -1^2)]
        """
        e = {}
        for i in self.zeros(fake_zeros=False):
            if i not in e: e[i] = 0
            e[i] += 1

        # the odd degrees (corresponding to angles 2((2m+1)+1) times pi should
        # be non ramified and hence come by pair.
        if any(e[i]%2 for i in e if i%2):
            return []

        pairings = []
        for d,m in iteritems(e):
            if d%2: # if the degree is odd it is necessarily non ramified
                pairings.append([(d,m//2)])
            else: # if the degree is even ramified and non ramified are possible
                pairings.append([(d,k) for k in range(m//2+1)])

        import itertools
        from .quadratic_strata import QuadraticStratum
        res = []

        for p in itertools.product(*pairings):
            ee = dict((d-1,0) for d in e)
            ee.update((2*d,0) for d in e)
            for d,m in p:
                ee[d-1] += e[d]-2*m
                ee[2*d] += m

            degrees = []
            for d in ee: degrees.extend([d]*ee[d])

            s = sum(degrees)
            for nb_poles in range(s%4,s+5,4):
                q = QuadraticStratum(degrees + [-1]*nb_poles)
                if not q.is_empty() and (not fake_zeros or q.nb_poles() <= self.nb_fake_zeros()):
                    res.append(q)

        return res

    #
    # Separatrix and cylinder diagrams
    #

    def separatrix_diagram_iterator(self, ncyls=None):
        r"""
        Return an iterator over the separatrix diagrams of this stratum.

        For strata of small dimension, it could be faster to use the method
        separatrix_diagrams.

        INPUT:

        - ``ncyls`` -- an optional number of cylinders
        """
        from .separatrix_diagram import separatrix_diagram_iterator
        return separatrix_diagram_iterator([m+1 for m in self.zeros()],ncyls)

    def separatrix_diagrams(self, ncyls=None):
        r"""
        Returns the list of separatrix diagrams that appears in this stratum.

        INPUT:

        - ``database`` - boolean (default: True) - if True, use the
          FlatSurfacesDatabase

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: a = AbelianStratum(2); a
            H_2(2)
            sage: for s in a.separatrix_diagrams(): print(s)
            (0,1,2)-(0,1,2)
            (0)(1,2)-(0,1)(2)

        TESTS::

            sage: from surface_dynamics import *

            sage: for (zeros, ncyl) in [((4,), 3), ((2,2), 4)]:
            ....:     S = AbelianStratum(4).separatrix_diagrams(3)
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

            sage: from surface_dynamics import *

            sage: A = AbelianStratum(4)
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
            if ncyls < 0 or ncyls > self.genus() + self.nb_zeros() - 1:
                raise ValueError("ncyls is not valid")

        if not force_computation:
            for cc in self.components():
                for cd in cc.cylinder_diagram_iterator(ncyls, up_to_symmetry, False):
                    yield cd
        else:
            for sd in self.separatrix_diagram_iterator(ncyls):
                iterator = sd.cylinder_diagram_iterator(up_to_symmetry=True)
                if not up_to_symmetry:
                    iterator = _cylinder_diagrams_with_symmetric(iterator)
                for cd in iterator:
                    yield cd

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

            sage: from surface_dynamics import *

            sage: A = AbelianStratum(2,2)
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
            sage: for c in AbelianStratum(2,1,1).cylinder_diagrams(2):
            ....:     total += 4 // (1 + sum(c.symmetries()))
            sage: total
            61
            sage: len(AbelianStratum(2, 1, 1).cylinder_diagrams(2, up_to_symmetry=False))
            61

        You obtain the same number directly::

            sage: AbelianStratum(2, 1, 1).cylinder_diagrams_number(2, up_to_symmetry=False)
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

            sage: from surface_dynamics import *

            sage: A = AbelianStratum(4)
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
            if ncyls < 0 or ncyls > self.genus() + self.nb_zeros()-1:
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

            sage: from surface_dynamics import *

            sage: a = AbelianStratum(3,2,1); a
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

            sage: from surface_dynamics import *

            sage: H22 = AbelianStratum(2,2)
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

            sage: H31 = AbelianStratum(3,1)
            sage: for d in range(1,5):
            ....:     print("%d %d" %(H31.cylinder_diagrams_number(d, True, False),
            ....:                     H31.cylinder_diagrams_number(d, True, True)))
            2 2
            12 12
            16 16
            4 4

            sage: H211 = AbelianStratum(2,1,1)
            sage: for d in range(1,6):
            ....:     print("%d %d" % (H211.cylinder_diagrams_number(d, True, False),
            ....:               H211.cylinder_diagrams_number(d, True, True)))
            5 5
            29 29
            53 53
            27 27
            8 8
        """
        if ncyls is not None and ncyls > self.genus() + self.nb_zeros() - 1:
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

            sage: from surface_dynamics import *

            sage: C = AbelianStratum(2,0)
            sage: p = C.single_cylinder_representative()
            sage: p
            0 1 2 3 4
            4 3 1 2 0
            sage: p.stratum() == C
            True

            sage: C = AbelianStratum(3,1)
            sage: p = C.single_cylinder_representative(alphabet=Alphabet(name='lower'))
            sage: p
            a b c d e f g
            c f b g e d a
            sage: p.stratum() == C
            True

            sage: C = AbelianStratum(2)
            sage: C.single_cylinder_representative()
            Traceback (most recent call last):
            ...
            ValueError: no 1,1-square-tiled surfaces in this stratum try again with H_2(2, 0)
            sage: C = AbelianStratum(1,1)
            sage: C.single_cylinder_representative()
            Traceback (most recent call last):
            ...
            ValueError: no 1,1-square-tiled surfaces in this stratum try again with H_2(1^2, 0^2)
        """
        genus = self.genus()
        nb_real_zeros = self.nb_zeros()-self.nb_fake_zeros()

        if genus == 2 and nb_real_zeros == 1 and self.nb_fake_zeros() < 1:
            raise ValueError("no 1,1-square-tiled surfaces in this stratum try again with H_2(2, 0)")
        elif genus == 2 and nb_real_zeros == 2 and self.nb_fake_zeros() < 2:
            raise ValueError("no 1,1-square-tiled surfaces in this stratum try again with H_2(1^2, 0^2)")

        return self.one_component().single_cylinder_representative(alphabet, reduced)

    def single_cylinder_origami(self):
        r"""
        Returns an origami associated to a single cylinder permutation representative.

        Returns an origami in this connected component having a single vertical
        cylinder and a single horizontal cylinder.

        Examples::

            sage: from surface_dynamics import *

            sage: C = AbelianStratum(4)
            sage: O = C.single_cylinder_origami()
            sage: O
            (1,2,3,4,5)
            (1,4,3,5,2)
            sage: O.stratum() == AbelianStratum(4)
            True
            sage: C = AbelianStratum(2,0)
            sage: O = C.single_cylinder_origami()
            sage: O
            (1,2,3,4)
            (1,3,2,4)
            sage: O.stratum() == AbelianStratum(2)
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

            sage: from surface_dynamics import *

            sage: c = AbelianStratum([1,1,1,1]).unique_component(); c
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

            sage: from surface_dynamics import *

            sage: c = AbelianStratum(1,1,1,1).unique_component()
            sage: p = c.permutation_representative(alphabet="abcdefghi")
            sage: p
            a b c d e f g h i
            e d c f i h g b a
            sage: p.stratum_component()
            H_3(1^4)^c

            sage: cc = AbelianStratum(3,2,1,0).unique_component()
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

        g = stratum.genus()
        zeros = list(stratum.zeros(fake_zeros=False))
        n = stratum.nb_fake_zeros()

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
            p = ReducedPermutationIET([l0, l1], reduced=True)
        else:
            from surface_dynamics.interval_exchanges.labelled import LabelledPermutationIET
            p = LabelledPermutationIET([l0, l1], reduced=False)

        p.alphabet(alphabet)
        return p

    def lyapunov_exponents_approx(self, **kargs):
        r"""
        Return the approximate Lyapunov exponents of the KZ-cocycle.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: AbelianStratum(2).unique_component().lyapunov_exponents_approx(nb_iterations=2**21)  # abs tol .05
            [1.000, 0.333]

            sage: H4hyp, H4odd = AbelianStratum(4).components()
            sage: H4hyp.lyapunov_exponents_approx(nb_iterations=2**21) # abs tol .05
            [1.000, 0.616, 0.184]
            sage: H4odd.lyapunov_exponents_approx(nb_iterations=2**21) # abs tol .05
            [1.000, 0.418, 0.182]
        """
        perm = self.permutation_representative(reduced=False)
        return perm.lyapunov_exponents_approx(**kargs)

    # TODO
    # def sum_of_lyapunov_exponents(self)
    # TODO
    # def volume(self)
    # TODO
    # def carea(self)

    def random_permutation(self, reduced=True, nsteps=None):
        r"""
        Return a random permutation obtained via a random walk in the extended rauzy diagram.

        INPUT:

        - ``nsteps`` - integer or None - perform nsteps and then stops as soon
          as a Strebel differential is found.

        At each step, with probability 1/3 we perform one of the following
        moves:

        - exchange top,bottom and left,right (proba 1/10)

        - top rauzy move (proba 9/20)

        - bot rauzy move (proba 9/20)

        EXAMPLES:

            sage: from surface_dynamics import *

            sage: C = AbelianStratum(10).hyperelliptic_component()
            sage: p = C.random_permutation(); p   # random
            0 1 2 3 4 5 6 7 8 9 10 11
            11 10 9 8 7 6 5 4 3 2 1 0
            sage: p.stratum_component()
            H_6(10)^hyp


        TESTS::

            sage: from surface_dynamics import *
            sage: from surface_dynamics.interval_exchanges.reduced import ReducedPermutationIET
            sage: from surface_dynamics.interval_exchanges.labelled import LabelledPermutationIET
            sage: for A in [AbelianStratum(6), AbelianStratum(1,1,1,1), AbelianStratum(3,3)]:
            ....:     p = C.random_permutation(reduced=True)
            ....:     assert isinstance(p, ReducedPermutationIET), C
            ....:     assert p._labels is None, C
            ....:     p = C.random_permutation(reduced=False)
            ....:     assert isinstance(p, LabelledPermutationIET), C
            ....:     assert p._labels is not None, C
        """
        import sage.misc.prandom as prandom

        p = self.permutation_representative(reduced=True)
        if nsteps is None:
            nsteps = 64 * self.stratum().dimension()

        for _ in range(nsteps):
            rd = prandom.random()
            if rd < 0.1:   # (inplace) symmetric with proba 1/10
                p._inversed_twin()
                p._reversed_twin()
            elif rd < .55: # (inplace) rauzy move top with proba 9/20
                p._move(1, 0, 1, p._twin[0][0])
            else:          # (inplace) rauzy move bot with proba 9/20
                p._move(0, 0, 0, p._twin[1][0])

        p._check()
        if reduced:
            return p
        if not reduced:
            from surface_dynamics.interval_exchanges.labelled import LabelledPermutationIET
            return LabelledPermutationIET([p[0], p[1]], reduced=False)

    def random_standard_permutation(self, reduced=True, nsteps=None):
        r"""
        Return a random standard permutation representative obtained by performing a random
        walk in the Rauzy diagram.

        INPUT:

        - ``nsteps`` - integer or None - perform nsteps and then stops as soon
          as a Strebel differential is found.

        EXAMPLES:

            sage: from surface_dynamics import *

            sage: C = AbelianStratum(10).hyperelliptic_component()
            sage: p = C.random_standard_permutation(); p   # random
            0 1 2 3 4 5 6 7 8 9 10 11
            11 10 9 8 7 6 5 4 3 2 1 0
            sage: p.stratum_component()
            H_6(10)^hyp

            sage: C = AbelianStratum(6,4,2).odd_component(); C
            H_7(6, 4, 2)^odd
            sage: p = C.random_standard_permutation(); p  # random
            0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
            15 2 14 12 3 11 6 10 8 5 9 13 7 4 1 0
            sage: p.stratum_component()
            H_7(6, 4, 2)^odd

            sage: C = AbelianStratum(2,2,2,2).even_component(); C
            H_5(2^4)^even
            sage: p = C.random_standard_permutation(); p  # random
            0 1 2 3 4 5 6 7 8 9 10 11 12
            12 4 9 11 8 3 7 6 1 10 2 5 0
            sage: p.stratum_component()
            H_5(2^4)^even

            sage: C = AbelianStratum(32).odd_component(); C
            H_17(32)^odd
            sage: p = C.random_standard_permutation(); p  # random
            0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33
            33 30 10 3 32 19 11 28 4 14 24 15 21 20 9 12 25 6 2 29 26 23 27 13 8 1 18 17 16 31 7 22 5 0
            sage: p.stratum_component()
            H_17(32)^odd
        """
        return self.random_permutation(reduced, nsteps).to_standard()

    def rauzy_diagram(self, *args, **kwds):
        r"""
        Returns the extended Rauzy diagram associated to this connected component.

        OUTPUT:

        rauzy diagram -- the Rauzy diagram associated to this stratum

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: c = AbelianStratum(0).components()[0]
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

            sage: a = AbelianStratum({1:4}).unique_component(); a
            H_3(1^4)^c
            sage: a.rauzy_diagram()
            Rauzy diagram with 1255 permutations
            sage: a.rauzy_class_cardinality()
            1255

            sage: cc = AbelianStratum(3,2,1).unique_component()
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

            sage: a = AbelianStratum({1:8}).unique_component(); a
            H_5(1^8)^c
            sage: a.rauzy_class_cardinality()
            55184875

        Cardinalities for labeled Rauzy classes instead of reduced::

            sage: cc=AbelianStratum(2,1,1).unique_component()
            sage: cc.rauzy_diagram(left_degree=2,reduced=False)
            Rauzy diagram with 3676 permutations
            sage: cc.rauzy_class_cardinality(left_degree=2,reduced=False)
            3676

            sage: cc.rauzy_diagram(left_degree=1,reduced=False)
            Rauzy diagram with 3774 permutations
            sage: cc.rauzy_class_cardinality(left_degree=1,reduced=False)
            3774

            sage: cc=AbelianStratum(2,1,1,0).unique_component()
            sage: cc.rauzy_diagram(left_degree=2,reduced=False) # long time
            Rauzy diagram with 33084 permutations
            sage: cc.rauzy_diagram(left_degree=1,reduced=False) # long time
            Rauzy diagram with 33966 permutations
            sage: cc.rauzy_diagram(left_degree=0,reduced=False) # long time
            Rauzy diagram with 30828 permutations

            sage: cc.rauzy_class_cardinality(left_degree=2,reduced=False)
            33084
            sage: cc.rauzy_class_cardinality(left_degree=1,reduced=False)
            33966
            sage: cc.rauzy_class_cardinality(left_degree=0,reduced=False)
            30828
        """
        import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rdc

        profile = list(map(lambda x: x+1, self.stratum().zeros()))
        s = self.stratum().nb_zeros()

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

            sage: from surface_dynamics import *

            sage: cc = AbelianStratum(3,1).unique_component()
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

            sage: cc = AbelianStratum({1:10}).unique_component(); cc
            H_6(1^10)^c
            sage: cc.standard_permutations_number()
            59520825
        """
        import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rdc

        profile = [x+1 for x in self.stratum().zeros()]

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

            sage: from surface_dynamics import *

            sage: C = AbelianStratum(4).odd_component()
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

            sage: from surface_dynamics import *

            sage: A = AbelianStratum(2,2).odd_component()
            sage: c = A.one_cylinder_diagram(); c
            (0,5,1,3,2,4)-(0,5,4,3,2,1)
            sage: c.stratum_component()
            H_3(2^2)^odd

            sage: A = AbelianStratum(3,3).non_hyperelliptic_component()
            sage: c = A.one_cylinder_diagram(); c
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

            sage: from surface_dynamics import *

            sage: A = AbelianStratum(1,1,1,1)
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

            sage: C1 = list(cc.cylinder_diagram_iterator(3, force_computation=False))  # long time
            sage: C2 = list(cc.cylinder_diagram_iterator(3, force_computation=True))   # long time
            sage: assert len(C1) == len(C2)                                            # long time
            sage: isoms = []                                                           # long time
            sage: for c in C1:                                                         # long time
            ....:     isom = []
            ....:     for i,cc in enumerate(C2):
            ....:         if c.is_isomorphic(cc) or \
            ....:            c.is_isomorphic(cc.horizontal_symmetry()) or \
            ....:            c.is_isomorphic(cc.vertical_symmetry()) or \
            ....:            c.is_isomorphic(cc.inverse()):
            ....:              isom.append(i)
            ....:     assert len(isom) == 1, isom
            ....:     isoms.extend(isom)
            sage: assert sorted(isoms) == list(range(len(C1)))                         # long time
        """
        if ncyls is not None:
            if not isinstance(ncyls, (int,Integer)):
                raise TypeError("ncyls should be None or an integer")
            if ncyls < 0 or ncyls > self.stratum().genus() + self.stratum().nb_zeros()-1:
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

            sage: from surface_dynamics import *

            sage: C = AbelianStratum(1,1,1,1).unique_component(); C
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

            sage: from surface_dynamics import *

            sage: C = AbelianStratum(3,1).unique_component()
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

            sage: C = AbelianStratum(6)
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
        if ncyls is not None and ncyls > self.stratum().genus() + self.stratum().nb_zeros() - 1:
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

            sage: from surface_dynamics import *

            sage: a = AbelianStratum(2,2).one_component()
            sage: a.one_origami().stratum()
            H_3(2^2)

            sage: AbelianStratum(3,2,1).unique_component().one_origami().stratum()
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

            sage: from surface_dynamics import *

            sage: cc = AbelianStratum(6).even_component()
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

            sage: from surface_dynamics import *

            sage: H11_hyp = AbelianStratum(1,1).hyperelliptic_component()
            sage: len(H11_hyp.origamis(6))
            88

            sage: T6 = H11_hyp.arithmetic_teichmueller_curves(6)
            sage: len(T6)
            5
            sage: sum(t.veech_group().index() for t in T6)
            88

            sage: H4_odd = AbelianStratum(4).odd_component()
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

            sage: from surface_dynamics import *

            sage: A = AbelianStratum(2).hyperelliptic_component(); A
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

            sage: A = AbelianStratum(1,1).hyperelliptic_component(); A
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

            sage: A = AbelianStratum(4).hyperelliptic_component(); A
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

    def lyapunov_exponents(self, **kargs):
        return(self.permutation_representative(reduced=False).lyapunov_exponents_H_plus(**kargs))

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

            sage: from surface_dynamics import *

            sage: cc = AbelianStratum(1,1,1,1).unique_component()
            sage: p = cc.single_cylinder_representative()
            sage: p
            0 1 2 3 4 5 6 7 8
            2 6 5 3 1 8 4 7 0
            sage: p.stratum_component() == cc
            True

            sage: cc = AbelianStratum(2,1,1).unique_component()
            sage: p = cc.single_cylinder_representative()
            sage: p
            0 1 2 3 4 5 6 7
            2 6 4 1 7 5 3 0
            sage: p.stratum_component() == cc
            True

            sage: cc = AbelianStratum(3,3).non_hyperelliptic_component()
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

        zeros = self.stratum().zeros()
        real_zeros = [z for z in zeros if z != 0]
        odd_zeros = [z for z in real_zeros if z%2 == 1]
        even_zeros = [z for z in real_zeros if z%2 == 0]

        fk_zeros_perm = GeneralizedPermutation([0],[0])
        mk_pt_perm = GeneralizedPermutation([0,1],[1,0])
        for i in range(self.stratum().nb_fake_zeros()):
            fk_zeros_perm = cylinder_concatenation(fk_zeros_perm,mk_pt_perm)

        if even_zeros == [2]:
            perm = only_even_2(odd_zeros)
        elif odd_zeros == [1,1]:
            perm = only_odds_11(even_zeros)
        else:
            if even_zeros:
                even_perm = AbelianStratum(even_zeros).odd_component().single_cylinder_representative()
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

            sage: from surface_dynamics import *

            sage: cc = AbelianStratum(4).odd_component()
            sage: O = cc.single_cylinder_origami()
            sage: O
            (1,2,3,4,5)
            (1,4,3,5,2)
            sage: O.stratum_component() == cc
            True
            sage: cc = AbelianStratum(5,3).unique_component()
            sage: O = cc.single_cylinder_origami()
            sage: O
            (1,2,3,4,5,6,7,8,9,10)
            (1,9,8,10,6,7,4,3,5,2)
            sage: O.stratum_component() == cc
            True
            sage: cc = AbelianStratum(4,2).even_component()
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

            sage: from surface_dynamics import *

        For the strata `H(2g-2)`::

            sage: c = AbelianStratum(0).hyperelliptic_component()
            sage: c.spin()
            1
            sage: p = c.permutation_representative()
            sage: p.arf_invariant()
            1

            sage: c = AbelianStratum(2).hyperelliptic_component()
            sage: c.spin()
            1
            sage: p = c.permutation_representative()
            sage: p.arf_invariant()
            1

            sage: c = AbelianStratum(4).hyperelliptic_component()
            sage: c.spin()
            0
            sage: p = c.permutation_representative()
            sage: p.arf_invariant()
            0

        For the strata `H(g-1,g-1)`::

            sage: c = AbelianStratum(2,2).hyperelliptic_component()
            sage: c.spin()
            0
            sage: p = c.permutation_representative()
            sage: p.arf_invariant()
            0

            sage: c = AbelianStratum(4,4).hyperelliptic_component()
            sage: c.spin()
            1
            sage: p = c.permutation_representative()
            sage: p.arf_invariant()
            1
        """
        z = self.stratum().zeros(fake_zeros=False)
        if not z:
            return Integer(1)
        elif len(z) == 1:
            return Integer(((self.stratum().genus()+1)//2) % 2)
        elif len(z) == 2:
            if z[0] % 2:
                return None
            return Integer(((self.stratum().genus()+1)//2) %2)

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

            sage: from surface_dynamics import *

            sage: c = AbelianStratum(0).hyperelliptic_component()
            sage: p = c.permutation_representative()
            sage: p
            0 1
            1 0
            sage: p.stratum_component()
            H_1(0)^hyp

            sage: c = AbelianStratum(0,0).hyperelliptic_component()
            sage: p = c.permutation_representative(alphabet="abc")
            sage: p
            a b c
            c b a
            sage: p.stratum_component()
            H_1(0^2)^hyp

            sage: c = AbelianStratum(2,2).hyperelliptic_component()
            sage: p = c.permutation_representative(alphabet="ABCDEFGHIJKL")
            sage: p
            A B C D E F G
            G F E D C B A
            sage: c = AbelianStratum(1,1,0).hyperelliptic_component()
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

        TESTS::

            sage: from surface_dynamics import AbelianStratum
            sage: from surface_dynamics.interval_exchanges.reduced import ReducedPermutationIET
            sage: from surface_dynamics.interval_exchanges.labelled import LabelledPermutationIET
            sage: for A in [AbelianStratum(6), AbelianStratum(1,1,1,1), AbelianStratum(3,3)]:
            ....:     for C in A.components():
            ....:         p = C.permutation_representative(reduced=True)
            ....:         assert isinstance(p, ReducedPermutationIET)
            ....:         assert p._labels is None, C
            ....:         p = C.permutation_representative(reduced=False)
            ....:         assert isinstance(p, LabelledPermutationIET)
            ....:         assert p._labels is not None, C
        """
        g = self._stratum.genus()
        n = self._stratum.nb_fake_zeros()
        m = len(self._stratum.zeros(fake_zeros=False))

        if left_degree is not None:
            if not isinstance(left_degree, (int,Integer)) or left_degree not in self.stratum().zeros():
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
            p = ReducedPermutationIET([l0, l1], reduced=True)
        else:
            from surface_dynamics.interval_exchanges.labelled import LabelledPermutationIET
            p = LabelledPermutationIET([l0, l1], reduced=False)
        if alphabet is not None:
            p.alphabet(alphabet)
        elif relabel:
            p.alphabet(range(len(p)))

        p._check()
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

            sage: from surface_dynamics import *

        The case of the torus is a little bit different::

            sage: c = AbelianStratum(0).hyperelliptic_component()
            sage: c.rauzy_diagram()
            Rauzy diagram with 1 permutation
            sage: c.rauzy_class_cardinality()
            1
            sage: c = AbelianStratum(0,0).hyperelliptic_component()
            sage: c.rauzy_diagram()
            Rauzy diagram with 3 permutations
            sage: c.rauzy_class_cardinality()
            3

        Examples in genus 2::

            sage: c = AbelianStratum(2,0).hyperelliptic_component()
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

            sage: c = AbelianStratum(1,1,0,0).hyperelliptic_component()
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

            sage: c = AbelianStratum(12,0,0).hyperelliptic_component()
            sage: c.rauzy_class_cardinality()
            1114200
            sage: c.rauzy_class_cardinality(left_degree=12, reduced=False)
            1965840

            sage: c = AbelianStratum(14).hyperelliptic_component()
            sage: c.rauzy_class_cardinality()
            32767
        """
        from sage.arith.all import binomial

        if left_degree is not None:
            assert isinstance(left_degree, (int,Integer)), "if not None, left_degree should be an integer"
            assert left_degree in self.stratum().zeros(), "if not None, the degree should be one of the degree of the stratum"

        if reduced is False:
            if left_degree is None:
                raise NotImplementedError("no formula known for cardinality of labeled extended Rauzy classes")
            zeros = self.stratum().zeros()
            profile = Partition([x+1 for x in zeros])
            if self.stratum().nb_zeros(fake_zeros=False) == 1:
                epsilon = 1
            else:
                epsilon = Rational((1,self.stratum().genus()))
            return epsilon * (profile.centralizer_size() /
                    ((left_degree+1) * zeros.count(left_degree)) *
                    self.rauzy_class_cardinality(left_degree=left_degree,reduced=True))

        k = self.stratum().nb_fake_zeros()
        dd = self.stratum().dimension()  # it is d+k
        d = dd-k

        if self.stratum().genus() == 1:
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
        if not self.stratum().nb_fake_zeros():
            return self.permutation_representative()

        raise NotImplementedError("not implemented when there are fake zeros")

    def standard_permutations(self, reduced=True):
        r"""
        Return the standard permutations in this hyperelliptic component.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: AbelianStratum(6).hyperelliptic_component().standard_permutations()
            [0 1 2 3 4 5 6 7
             7 6 5 4 3 2 1 0]
        """
        if not self.stratum().nb_fake_zeros():
            d = self.stratum().dimension()
            l0 = list(range(d))
            l1 = list(range(d-1,-1,-1))

            if reduced:
                from surface_dynamics.interval_exchanges.reduced import ReducedPermutationIET
                p = ReducedPermutationIET([l0, l1], reduced=True)
            else:
                from surface_dynamics.interval_exchanges.labelled import LabelledPermutationIET
                p = LabelledPermutationIET([l0, l1], reduced=False)

            return [p]

        raise NotImplementedError("not implemented when there are fake zeros")

    def standard_permutations_number(self):
        r"""
        Return the number of standard permutations in this hyperelliptic
        component.
        """
        if not self.stratum().nb_fake_zeros():
            return Integer(1)

        raise NotImplementedError("not implemented when there are fake zeros")

    def _cylinder_diagram_iterator(self, ncyls=None):
        r"""
        Returns the list of cylinder diagrams associated to this hyperelliptic
        component.

        INPUT:

        - ``ncyls`` -- an optional number of cylinders

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: C = AbelianStratum(2,2).hyperelliptic_component()
            sage: [sum(1 for c in C.cylinder_diagram_iterator(n)) for n in range(1,5)]
            [1, 3, 5, 2]

        When ``ncyls`` is set to ``None``, the iterator can reasonably be used
        with very large data::

            sage: C = AbelianStratum(10,10).hyperelliptic_component()
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

        if stratum.nb_fake_zeros():
            raise ValueError("the stratum has fake zeros")

        z = stratum.zeros()

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

            sage: from surface_dynamics import *

            sage: cc = AbelianStratum(2,0).hyperelliptic_component()
            sage: p = cc.single_cylinder_representative(alphabet=Alphabet(name='upper'))
            sage: p
            A B C D E
            E D B C A
            sage: p.stratum_component() == cc
            True

            sage: cc = AbelianStratum({3:2,0:6}).hyperelliptic_component()
            sage: p = cc.single_cylinder_representative()
            sage: p
            0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
            14 12 13 10 11 8 9 7 5 6 3 4 1 2 0
            sage: p.stratum_component() == cc
            True

            sage: cc = AbelianStratum(2).hyperelliptic_component()
            sage: cc.single_cylinder_representative()
            Traceback (most recent call last):
            ...
            ValueError: no 1,1-square-tiled surfaces in this connected component try again with H_2(2, 0)^hyp
            sage: cc = AbelianStratum({3:2,0:5}).hyperelliptic_component()
            sage: cc.single_cylinder_representative()
            Traceback (most recent call last):
            ...
            ValueError: no 1,1-square-tiled surfaces in this connected component try again with H_4(3^2, 0^6)^hyp
        """
        stratum = self.stratum()
        genus = stratum.genus()
        nb_fk_zeros = stratum.nb_fake_zeros()
        nb_real_zeros = stratum.nb_zeros()-nb_fk_zeros
        add_fk_zeros = nb_fk_zeros - 2*genus+4-nb_real_zeros

        from surface_dynamics.interval_exchanges.constructors import GeneralizedPermutation

        if nb_real_zeros == 1 and add_fk_zeros < 0:
            raise ValueError("no 1,1-square-tiled surfaces in this connected component try again with %s^hyp" %(str(AbelianStratum({2*genus-2:1,0:2*genus-3}))))
        elif nb_real_zeros == 2 and add_fk_zeros < 0:
            raise ValueError("no 1,1-square-tiled surfaces in this connected component try again with %s^hyp" %(str(AbelianStratum({genus-1:2,0:2*genus-2}))))
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

            sage: from surface_dynamics import *

        Examples in genus 3::

            sage: c = AbelianStratum(3,3).non_hyperelliptic_component()
            sage: c.rauzy_class_cardinality()
            15568

            sage: c = AbelianStratum(3,3,0).non_hyperelliptic_component()
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

            sage: c = AbelianStratum(5,5).non_hyperelliptic_component()
            sage: c.rauzy_class_cardinality()
            136116680

            sage: c = AbelianStratum(7,7,0).non_hyperelliptic_component()
            sage: c.rauzy_class_cardinality()
            88484743236111
            sage: c.rauzy_class_cardinality(left_degree=7, reduced=False)
            334071852804864
        """
        import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rdc

        profile = list(map(lambda x: x+1,self.stratum().zeros()))
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

            sage: from surface_dynamics import *

            sage: C = AbelianStratum(3,3).non_hyperelliptic_component()
            sage: len(C.standard_permutations())  # long time
            275
            sage: C.standard_permutations_number()
            275

            sage: C = AbelianStratum(5,5).non_hyperelliptic_component()
            sage: C.standard_permutations_number()
            1022399

            sage: C = AbelianStratum(7,7).non_hyperelliptic_component()
            sage: C.standard_permutations_number()
            19229011199
        """
        import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rdc

        profile = list(map(lambda x: x+1, self.stratum().zeros()))
        return rdc.number_of_standard_permutations(profile) - self.stratum().hyperelliptic_component().standard_permutations_number()

    def _cylinder_diagram_iterator(self, ncyls=None):
        r"""
        Return the list of cylinder diagrams (or completely periodic
        configurations) associated to this non-hyperelliptic component.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: cc = AbelianStratum(3,3).non_hyperelliptic_component()
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

            sage: from surface_dynamics import *

            sage: c = AbelianStratum(4,2).even_component(); c
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

            sage: from surface_dynamics import *

            sage: c = AbelianStratum(6).even_component()
            sage: c
            H_4(6)^even
            sage: p = c.permutation_representative(alphabet=range(8))
            sage: p
            0 1 2 3 4 5 6 7
            5 4 3 2 7 6 1 0
            sage: p.stratum_component()
            H_4(6)^even

        ::

            sage: c = AbelianStratum(4,4).even_component()
            sage: c
            H_5(4^2)^even
            sage: p = c.permutation_representative(alphabet=range(11))
            sage: p
            0 1 2 3 4 5 6 7 8 9 10
            5 4 3 2 6 8 7 10 9 1 0
            sage: p.stratum_component()
            H_5(4^2)^even

        Different markings lead to different Rauzy diagrams::

            sage: c = AbelianStratum(4,2,0).even_component()
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
        z = list(self._stratum.zeros(fake_zeros=False))
        n = self._stratum.nb_fake_zeros()
        g = self._stratum.genus()

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
            for i in range(d/2-1):
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
            p = ReducedPermutationIET([l0, l1], reduced=True)
        else:
            from surface_dynamics.interval_exchanges.labelled import LabelledPermutationIET
            p = LabelledPermutationIET([l0, l1], reduced=False)

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

            sage: from surface_dynamics import *

            sage: c = AbelianStratum(6).even_component()
            sage: c.rauzy_diagram()
            Rauzy diagram with 2327 permutations
            sage: c.rauzy_class_cardinality()
            2327

            sage: c = AbelianStratum(4,2,0).even_component()
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

        profile = list(map(lambda x: x+1, self.stratum().zeros()))
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

            sage: from surface_dynamics import *

        For strata in genus 3, the number of standard permutations is reasonably
        small and the whole set can be computed::

            sage: C = AbelianStratum(6).even_component()
            sage: len(C.standard_permutations())  # long time
            44
            sage: C.standard_permutations_number()
            44

            sage: C = AbelianStratum(4,2).even_component()
            sage: len(C.standard_permutations())   # long time
            136
            sage: C.standard_permutations_number()
            136

            sage: C = AbelianStratum(2,2,2).even_component()
            sage: len(C.standard_permutations())   # long time
            92
            sage: C.standard_permutations_number()
            92

        For higher genera, this number can be very big::

            sage: C = AbelianStratum(20).even_component()
            sage: C.standard_permutations_number()
            109398514483439999
        """
        import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rdc

        profile = [x+1 for x in self.stratum().zeros()]
        N = Integer(rdc.gamma_std(profile) - rdc.delta_std(profile)) / 2

        if (self.stratum().number_of_components() == 3 and
            self.stratum().hyperelliptic_component().spin() == 0):
            return N - 1

        return N

    def _cylinder_diagram_iterator(self, ncyls=None):
        r"""
        Iterator over cylinder diagram.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: cc = AbelianStratum(4,2).even_component()
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

            sage: from surface_dynamics import *

            sage: cc = AbelianStratum(6).even_component()
            sage: p = cc.single_cylinder_representative(alphabet=Alphabet(name='lower'))
            sage: p
            a b c d e f g h
            c h g f d b e a
            sage: p.stratum_component() == cc
            True

            sage: cc = AbelianStratum(4,4).even_component()
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

        zeros = self.stratum().zeros()
        real_zeros = [z for z in zeros if z != 0]

        fk_zeros_perm = GeneralizedPermutation([0],[0])
        mk_pt_perm = GeneralizedPermutation([0,1],[1,0])
        for i in range(self.stratum().nb_fake_zeros()):
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

            sage: from surface_dynamics import *

            sage: c = AbelianStratum(4).odd_component(); c
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

            sage: from surface_dynamics import *

            sage: a = AbelianStratum(6).odd_component()
            sage: p = a.permutation_representative()
            sage: p
            0 1 2 3 4 5 6 7
            3 2 5 4 7 6 1 0
            sage: p.stratum_component()
            H_4(6)^odd

        ::

            sage: a = AbelianStratum(4,4).odd_component()
            sage: p = a.permutation_representative()
            sage: p
            0 1 2 3 4 5 6 7 8 9 10
            3 2 5 4 6 8 7 10 9 1 0
            sage: p.stratum_component()
            H_5(4^2)^odd

        Different markings lead to different Rauzy diagrams::

            sage: c = AbelianStratum(4,2,0).odd_component()
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
        zeros = list(self.stratum().zeros(fake_zeros=False))
        n = self._stratum.nb_fake_zeros()
        g = self._stratum.genus()

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
            p = ReducedPermutationIET([l0, l1], reduced=True)
        else:
            from surface_dynamics.interval_exchanges.labelled import LabelledPermutationIET
            p = LabelledPermutationIET([l0, l1], reduced=False)

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

            sage: from surface_dynamics import *

        The genus must be at least 3 to have an odd component::

            sage: c = AbelianStratum(4).odd_component()
            sage: c.rauzy_diagram()
            Rauzy diagram with 134 permutations
            sage: c.rauzy_class_cardinality()
            134
            sage: c = AbelianStratum(4,0).odd_component()
            sage: c.rauzy_diagram()
            Rauzy diagram with 1114 permutations
            sage: c.rauzy_class_cardinality()
            1114

            sage: c = AbelianStratum(2,2).odd_component()
            sage: c.rauzy_diagram()
            Rauzy diagram with 294 permutations
            sage: c.rauzy_class_cardinality()
            294

            sage: c = AbelianStratum(2,2,0).odd_component()
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

            sage: c = AbelianStratum(4,2,0).odd_component()
            sage: c.rauzy_class_cardinality()
            262814
            sage: c = AbelianStratum(4,4,4).odd_component()
            sage: c.rauzy_class_cardinality()
            24691288838
            sage: c.rauzy_class_cardinality(left_degree=4, reduced=False)
            1234564441900
        """
        import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rdc

        profile = list(map(lambda x: x+1, self.stratum().zeros()))
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

            sage: from surface_dynamics import *

        In genus 2, there are two strata which contains an odd component::

            sage: C = AbelianStratum(4).odd_component()
            sage: len(C.standard_permutations())
            7
            sage: C.standard_permutations_number()
            7

            sage: C = AbelianStratum(2,2).odd_component()
            sage: len(C.standard_permutations())
            11
            sage: C.standard_permutations_number()
            11

        In genus 3, the number of standard permutations is reasonably small and
        the whole set can be computed::

            sage: C = AbelianStratum(6).odd_component()
            sage: len(C.standard_permutations())   # long time
            135
            sage: C.standard_permutations_number()
            135

            sage: C = AbelianStratum(4,2).odd_component()
            sage: len(C.standard_permutations())   # long time
            472
            sage: C.standard_permutations_number()
            472

            sage: C = AbelianStratum(2,2,2).odd_component()
            sage: len(C.standard_permutations())   # long time
            372
            sage: C.standard_permutations_number()
            372

        For higher genera, this number can be very big::

            sage: C = AbelianStratum(8,6,4,2).odd_component()
            sage: C.standard_permutations_number()
            26596699869748377600
        """
        import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rdc

        profile = [x+1 for x in self.stratum().zeros()]
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

            sage: from surface_dynamics import *

            sage: C = AbelianStratum(4).odd_component()
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

            sage: from surface_dynamics import *

            sage: cc = AbelianStratum(4).odd_component()
            sage: p = cc.single_cylinder_representative(alphabet=Alphabet(name='upper'))
            sage: p
            A B C D E F
            C F E B D A
            sage: p.stratum_component() == cc
            True

            sage: cc = AbelianStratum(6,2).odd_component()
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

        zeros = self.stratum().zeros()
        real_zeros = [z for z in zeros if z != 0]

        fk_zeros_perm = GeneralizedPermutation([0],[0])
        mk_pt_perm = GeneralizedPermutation([0,1],[1,0])
        for i in range(self.stratum().nb_fake_zeros()):
            fk_zeros_perm = cylinder_concatenation(fk_zeros_perm,mk_pt_perm)

        two_count = real_zeros.count(2)
        if two_count == 0:
            perm = cylinder_concatenation(fk_zeros_perm,no_two_odd(real_zeros))
        elif two_count == 1 :
            perm = cylinder_concatenation(fk_zeros_perm,one_two_odd(real_zeros))
        elif two_count >= 2 and two_count%2 == 0:
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

            sage: from surface_dynamics import *

            sage: a = AbelianStrata(genus=3)
            sage: all(s in a for s in a)
            True

            sage: a = AbelianStrata(genus=3,fake_zeros=False)
            sage: all(s in a for s in a)
            True

            sage: a = AbelianStrata(dimension=7,fake_zeros=True)
            sage: all(s in a for s in a)
            True
            sage: AbelianStratum(2,0,0) in a
            False

            sage: a = AbelianStrata(dimension=7,fake_zeros=False)
            sage: all(s in a for s in a)
            True
            sage: AbelianStratum(4,0) in a
            False
        """
        if not isinstance(c, AbelianStratum):
            return False

        return ((self._genus is None or c.genus() == self._genus) and
                (self._dimension is None or c.dimension() == self._dimension) and
                (self._fake_zeros is None or self._fake_zeros or not c.nb_fake_zeros()))


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
            yield AbelianStratum([0])
        else:
            for p in Partitions(2*self._genus-2):
                yield AbelianStratum(p._list)

    def random_element(self):
        r"""
        Return a random stratum.
        """
        if self._genus == 0:
            raise ValueError("No stratum with that genus")
        if self._genus == 1:
            return AbelianStratum([0])
        return AbelianStratum(Partitions(2*self._genus - 2).random_element())

    def first(self):
        r"""
        Return the first element of this list of strata.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: AbelianStrata(genus=3).first()
            H_3(4)
            sage: AbelianStrata(genus=4).first()
            H_4(6)
        """
        return AbelianStratum([2*self._genus-2])

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
        return AbelianStratum({1:2*self._genus-2})


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

            sage: from surface_dynamics import *

            sage: AbelianStrata(dimension=2).first()
            H_1(0)
            sage: AbelianStrata(dimension=3).first()
            H_1(0^2)
            sage: AbelianStrata(dimension=4).first()
            H_2(2)
        """
        n = self._dimension
        if n%2:
            return AbelianStratum([(n-3)//2,(n-3)//2])
        return AbelianStratum([n-2])

    an_element = first

    def last(self):
        r"""
        Return the last stratum.

        EXAMPLES::

            sage: from surface_dynamics import *

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
            return AbelianStratum({0:n-1})
        else:
            if n == 4:
                return AbelianStratum([2])
            if n == 5:
                return AbelianStratum([1,1])
            elif n == 6:
                return AbelianStratum([4])
            else:
                nn = (n-2)%4
                return AbelianStratum({2:3-nn,1:2*((n-10)//4)+2*nn})

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
                    yield AbelianStratum([k-1 for k in p])
        else:
            if n == 2:
                yield AbelianStratum([0])
            else:
                for s in range(1+n%2, n, 2):
                    for p in Partitions(n-1,length=s,min_part=2):
                        yield AbelianStratum([k-1 for k in p])

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
                yield AbelianStratum([0]*(self._dimension-1))
        else:
            s = self._dimension - 2*self._genus + 1
            if self._fake_zeros:
                for p in Partitions(2*self._genus - 2 + s, length=s):
                    yield AbelianStratum([k-1 for k in p])
            else:
                for p in Partitions(2*self._genus - 2, length=s):
                    yield AbelianStratum(p)

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
            for stratum in AbelianStrata(dimension=d, fake_zeros=self._fake_zeros):
                yield stratum

