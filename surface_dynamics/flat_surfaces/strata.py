r"""
Strata of differential on Riemann surfaces

This file gather common code used in
:mod:`~surface_dynamics.flat_surfaces.abelian_strata` and
:mod:`~surface_dynamics.flat_surfaces.quadratic_strata`.
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

from functools import total_ordering

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.sage_object import SageObject
from sage.structure.parent import Parent

def list_to_exp_list(l):
    r"""
    Convert list into exponential notation.

    EXAMPLES::

        sage: from surface_dynamics import *

        sage: from surface_dynamics.flat_surfaces.strata import list_to_exp_list
        sage: l = [0,0,2,2,3,2,0,0,0]
        sage: list_to_exp_list(l)
        [(0, 2), (2, 2), (3, 1), (2, 1), (0, 3)]
    """
    d = []
    i = 0
    while i < len(l):
        j = i
        while j < len(l) and l[j] == l[i]:
            j += 1
        d.append((l[i],j-i))
        i = j
    return d

#
# Stratum, stratum component
#

class Stratum(UniqueRepresentation, SageObject):
    r"""
    Generic class for stratum of flat surfaces.

    Assumes there are

    - a method ``.zeros()`` which returns the list of all zeros

    - a method ``.nb_zeros()`` which returns the number of zeros with an option
      fake_zeros which could be true or false

    - a method ``.nb_fake_zeros()`` which returns the number of fake zeros (or
      marked points)

    - a method ``.dimension()`` which returns the dimension of the stratum

    - an attribute ``._cc`` which is a list of classes associated to the
      connected components of self

    There may be

    - an attribute ``._name`` which corresponds to the begining of the string
      representation (default is the empty string)

    - an attribute ``._latex_name`` which corresponds to the begining of the latex
      string representation (uses ``_name`` by default)

    TESTS::

            sage: from surface_dynamics import *

            sage: A = AbelianStratum(2,2)
            sage: B = AbelianStratum(1,1)
            sage: hash(A) == hash(B)
            False
    """
    _name = ''
    _latex_name = ''

    #
    # String representation
    #

    def _flat_zero_str(self):
        r"""
        String representation of the zeros.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: a = AbelianStratum({2:3})
            sage: a._flat_zero_str()
            '2, 2, 2'
        """
        return ', '.join(map(str,self.zeros()))

    def _exp_zero_str(self):
        r"""
        String representation with exponential notation

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: a = AbelianStratum(2,2,2)
            sage: a._exp_zero_str()
            '2^3'
        """
        return ', '.join('%d^%d' %(i,e) if e != 1 else '%d' %i for (i,e) in list_to_exp_list(self.zeros()))

    # this attribute can be switched between _flat_zero_str and _exp_zero_str
    _zero_str = _exp_zero_str

    def _repr_(self):
        """
        TESTS::

            sage: from surface_dynamics import *

            sage: repr(AbelianStratum(1,1))       # indirect doctest
            'H_2(1^2)'
            sage: repr(QuadraticStratum(1,1,1,1)) # indirect doctest
            'Q_2(1^4)'
        """
        return self._name + "_" + str(self.genus()) + "(" + self._zero_str() + ")"

    def _latex_(self):
        r"""
        Latex string representation

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: AbelianStratum(0)._latex_()
            '\\mathcal{H}_1(0)'
            sage: QuadraticStratum({-1:4})._latex_()
            '\\mathcal{Q}_0(-1^4)'
        """
        return self._latex_name + '_' + str(self.genus()) + "(" + self._zero_str() + ")"

    #
    # Equality and comparisons
    #

    def __lt__(self, other):
        r"""
        Comparison

        ALGORITHM:

        First compare the class, then the dimension and then the list of zeros.

        TESTS::

            sage: from surface_dynamics import *

            sage: AbelianStratum(0) == AbelianStratum(0)
            True
            sage: QuadraticStratum(5,-1) == QuadraticStratum(5,-1)
            True
            sage: QuadraticStratum(5,-1) != QuadraticStratum(5,-1)
            False

            sage: AbelianStratum(12) == QuadraticStratum(12)
            False
            sage: QuadraticStratum(12,0,0) == QuadraticStratum(12,0)
            False

            sage: AbelianStratum(2,0) == AbelianStratum(2)
            False
            sage: AbelianStratum(2,0) != AbelianStratum(2)
            True

            sage: AbelianStratum(1,1) < AbelianStratum(1,1,0)
            True
            sage: AbelianStratum(1,1,0) < AbelianStratum(1,1)
            False
            sage: AbelianStratum(1,1,0) < AbelianStratum(1,1,0,0)
            True
            sage: AbelianStratum(2) < AbelianStratum(1,1)
            True
            sage: AbelianStratum(4,0) > AbelianStratum(1,1,1,1)
            False
            sage: AbelianStratum(4,0,0,0) > AbelianStratum(1,1,1,1)
            True

        ::

            sage: QuadraticStratum(2,2) < QuadraticStratum(2,2,0)
            True
            sage: QuadraticStratum(2,2,0) < QuadraticStratum(2,2)
            False
            sage: QuadraticStratum(2,2,0) < QuadraticStratum(2,2,0,0)
            True
            sage: QuadraticStratum(4) < QuadraticStratum(2,2)
            True
            sage: QuadraticStratum(4,0) > QuadraticStratum(1,1,1,1)
            False

            sage: Q1 = QuadraticStratum(4,0,0,0)
            sage: Q2 = QuadraticStratum(1,1,1,1)
            sage: Q1 > Q2
            True
            sage: Q1 >= Q2
            True
            sage: Q1 < Q2
            False
            sage: Q1 <= Q2
            False
        """
        if not isinstance(self, Stratum) or not isinstance(other, Stratum):
            raise TypeError

        if type(self) is type(other):
            # compare the dimension
            if self.dimension() < other.dimension():
                return True
            elif self.dimension() > other.dimension():
                return False

            # compare the list of zeros
            if self.zeros() < other.zeros():
                return True
            elif self.zeros() > other.zeros():
                return False

            # equality
            return False

        else:
            sname = type(self).__name__
            oname = type(other).__name__
            if sname < oname:
                return True
            elif sname > oname:
                return False

            return False
    
    def __ne__(self, other):
        return not self == other

    def __le__(self, other):
        return self == other or self < other

    def __ge__(self, other):
        return not self < other

    def __gt__(self, other):
        return not self <= other

    #
    # Connected components
    #

    def is_connected(self):
        r"""
        Test if the strata is connected.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: AbelianStratum([2]).is_connected()
            True
            sage: AbelianStratum([2,2]).is_connected()
            False
            sage: QuadraticStratum([-1,-1,-1,-1]).is_connected()
            True
            sage: QuadraticStratum([12]).is_connected()
            False
        """
        return len(self._cc) <= 1

    def permutation_representative(self, *args, **kwds):
        r"""
        Return a permutation of interval exchanges associated to this stratum.

        EXAMPLES::

            sage: from surface_dynamics import *

        Examples from Abelian differentials::

            sage: a = AbelianStratum([3,2,1,0,0])
            sage: p = a.permutation_representative()
            sage: p.stratum()
            H_4(3, 2, 1, 0^2)
            sage: a = AbelianStratum([2, 2, 2])
            sage: p = a.permutation_representative()
            sage: p.stratum()
            H_4(2^3)

        Examples from quadratic differentials::

            sage: a = QuadraticStratum([6,-1,-1])
            sage: p = a.permutation_representative()
            sage: p.stratum()
            Q_2(6, -1^2)
            sage: a = QuadraticStratum([-1,-1,-1,-1,0,0])
            sage: p = a.permutation_representative()
            sage: p.stratum()
            Q_0(0^2, -1^4)
        """
        return self.one_component().permutation_representative(*args,**kwds)

    def is_empty(self):
        r"""
        Return True if the stratum is empty

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: AbelianStratum(2).is_empty()
            False
            sage: QuadraticStratum(1,-1).is_empty()
            True
        """
        return len(self._cc) == 0

    def number_of_components(self):
        r"""
        Returns the number of connected components of self

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: AbelianStratum(2).number_of_components()
            1
            sage: AbelianStratum(4).number_of_components()
            2
            sage: AbelianStratum(3,3).number_of_components()
            2
        """
        return len(self._cc)

    def one_component(self):
        r"""
        Returns a connected component of this stratum.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: AbelianStratum(2).one_component()
            H_2(2)^hyp
        """
        if self.components():
            return self.components()[-1]
        from sage.categories.sets_cat import EmptySetError
        raise EmptySetError("The stratum is empty")

    def unique_component(self):
        r"""
        Returns the unique component of self or raise a ValueError.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: a = AbelianStratum(1,1); a
            H_2(1^2)
            sage: a.unique_component()
            H_2(1^2)^hyp

            sage: a = AbelianStratum(3,2,1); a
            H_4(3, 2, 1)
            sage: a.unique_component()
            H_4(3, 2, 1)^c

            sage: QuadraticStratum({1:1, -1:5}).unique_component()
            Q_0(1, -1^5)^c
            sage: QuadraticStratum(3,2,-1).unique_component()
            Q_2(3, 2, -1)^nonhyp

            sage: QuadraticStratum(12).unique_component()
            Traceback (most recent call last):
            ...
            ValueError: several components for this stratum
        """
        if len(self._cc) != 1:
            raise ValueError("several components for this stratum")
        return self._cc[0](self)

    def random_component(self):
        r"""
        Returns a random connected component of this stratum.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: Q = QuadraticStratum(6,6)
            sage: Q.random_component()
            Q_4(6^2)^hyp
            sage: Q.random_component()
            Q_4(6^2)^reg
        """
        if self.components():
            from sage.misc.prandom import choice
            return choice(self.components())
        from sage.categories.sets_cat import EmptySetError
        raise EmptySetError("The stratum is empty")

    def components(self):
        """
        Lists the connected components of the Stratum.

        OUTPUT:

        list -- a list of connected components of stratum

        EXAMPLES::

            sage: from surface_dynamics import *

        Some abelian strata::

            sage: AbelianStratum(0).components()
            [H_1(0)^hyp]
            sage: AbelianStratum(2).components()
            [H_2(2)^hyp]
            sage: AbelianStratum(4).components()
            [H_3(4)^hyp, H_3(4)^odd]
            sage: AbelianStratum(2,2).components()
            [H_3(2^2)^hyp, H_3(2^2)^odd]
            sage: AbelianStratum(1,1,1,1).components()
            [H_3(1^4)^c]

        Some quadratic strata::

            sage: QuadraticStratum(12).components()
            [Q_4(12)^reg, Q_4(12)^irr]
            sage: QuadraticStratum(6,-1,-1).components()
            [Q_2(6, -1^2)^hyp, Q_2(6, -1^2)^nonhyp]
        """
        return list(map(lambda x: x(self), self._cc))

class StratumComponent(SageObject):
    r"""
    Generic class for connected component of a stratum of flat surfaces.

    Assumes there are implemented

    - a method .permutation_representative()

    There may be

    - an attribute ._name

    - an attribute ._latex_name

    """
    _name = ''
    _latex_name = ''

    def __init__(self, stratum):
        r"""
        TEST::

            sage: from surface_dynamics import *

            sage: a = AbelianStratum(4,4).one_component()
            sage: a == loads(dumps(a))
            True
            sage: q = QuadraticStratum(5,5,-1,-1).one_component()
            sage: q == loads(dumps(q))
            True
        """
        self._stratum = stratum

    def __reduce__(self):
        r"""
        Reduce method for pickling

        TESTS::

            sage: from surface_dynamics import *

        Tests for Abelian strata::

            sage: a = AbelianStratum(2,2)
            sage: all(loads(dumps(cc)) == cc for cc in a.components())
            True
            sage: a = AbelianStratum(3,3)
            sage: all(loads(dumps(cc)) == cc for cc in a.components())
            True
            sage: a = AbelianStratum(6)
            sage: all(loads(dumps(cc)) == cc for cc in a.components())
            True
            sage: a = AbelianStratum(1,1,1,1)
            sage: all(loads(dumps(cc)) == cc for cc in a.components())
            True

        Tests for quadratic strata::

            sage: q = QuadraticStratum(-1,-1,-1,-1)
            sage: all(loads(dumps(cc)) == cc for cc in q.components())
            True
            sage: q = QuadraticStratum(12)
            sage: all(loads(dumps(cc)) == cc for cc in q.components())
            True
        """
        return (self.__class__, (self._stratum,))

    def _repr_(self):
        r"""
        String representation

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: a_hyp = AbelianStratum(4).hyperelliptic_component()
            sage: a_hyp._repr_()
            'H_3(4)^hyp'
            sage: a_odd = AbelianStratum(4).odd_component()
            sage: a_odd._repr_()
            'H_3(4)^odd'
        """
        return str(self._stratum) + "^" + self._name

    def __hash__(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *

            sage: A4hyp = AbelianStratum(4).hyperelliptic_component()
            sage: A4odd = AbelianStratum(4).odd_component()
            sage: hash(A4hyp) != hash(A4odd)
            True
            sage: A22hyp = AbelianStratum(2,2).hyperelliptic_component()
            sage: hash(A22hyp) != hash(A4hyp)
            True
        """
        return hash(self._stratum) ^ hash(self._name)

    def stratum(self):
        r"""
        Return the stratum associated to self

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: a = AbelianStratum(4,4)
            sage: all([c.stratum() == a for c in a.components()])
            True
        """
        return self._stratum

    def genus(self):
        r"""
        Return genus of the corresponding stratum

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: a = AbelianStratum(4,4)
            sage: a.one_component().genus()
            5
        """
        return self._stratum.genus()

    def __eq__(self,other):
        r"""
        Equality test

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: c_hyp = AbelianStratum(6).hyperelliptic_component()
            sage: c_odd = AbelianStratum(6).odd_component()
            sage: c_hyp == c_hyp
            True
            sage: c_hyp == c_odd
            False
        """
        if not isinstance(self, StratumComponent) or not isinstance(other, StratumComponent):
            return NotImplemented

        return type(self) == type(other) and self._stratum == other._stratum

    def __lt__(self, other):
        r"""
        Comparison

        TESTS::

            sage: from surface_dynamics import *

            sage: a1 = AbelianStratum(1,1,1,1)
            sage: c1 = a1.components()[0]
            sage: a2 = AbelianStratum(3,1)
            sage: c2 = a2.components()[0]
            sage: c1 == c1
            True
            sage: c1 == c2
            False
            sage: a1 = AbelianStratum(1,1,1,1)
            sage: c1 = a1.components()[0]
            sage: a2 = AbelianStratum(2, 2)
            sage: c2_hyp, c2_odd = a2.components()
            sage: c1 != c1
            False
            sage: c1 != c2_hyp
            True
            sage: c2_hyp != c2_odd
            True
        """
        if not isinstance(self, StratumComponent) or not isinstance(other, StratumComponent):
            return NotImplemented

        if self._stratum < other._stratum:
            return True
        elif self._stratum > other._stratum:
            return False

        if type(self) < type(other):
            return True
        elif type(self) > type(other):
            return False

#
# Strata (family of strata)
#

class Strata(Parent):
    r"""
    Strata of Abelian or Quadratic differentials.
    """
    pass
