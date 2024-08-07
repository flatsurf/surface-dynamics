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

import numbers

from functools import total_ordering

from sage.misc.cachefunc import cached_method
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
    Stratum of holomorphic or meromorphic k-differentials on smooth connected Riemann surfaces.

    INPUT:

    - ``signature`` -- a list of ``n`` integers (determine the angles of
      conical singularities)

    - ``k`` -- optional integer (default ``1``) the order of the differential

    EXAMPLES::

        sage: from surface_dynamics import Stratum
        sage: Stratum((2,), k=1)
        H_2(2)

    TESTS::

        sage: S = Stratum((2,), k=1)
        sage: loads(dumps(S)) == S
        True
        sage: S = Stratum((2, 2, 2), k=3)
        sage: loads(dumps(S)) == S
        True
    """
    @staticmethod
    def __classcall_private__(self, signature, k=1):
        if not isinstance(k, numbers.Integral) or k <= 0:
            raise ValueError('k must be a positive integer')
        k = int(k)
        if isinstance(signature, dict):
            signature = sum(([i] * mult for i, mult in signature.items()), [])
        signature = tuple(signature)
        if not signature:
            raise ValueError('the signature must be non-empty')
        for m in signature:
            if not isinstance(m, numbers.Integral):
                raise ValueError('mu must be a list of integers')

        signature = tuple(sorted(map(int, signature), reverse=True))
        return super().__classcall__(Stratum, signature, k)

    def __new__(cls, signature, k):
        if k == 1:
            from .abelian_strata import AbelianStratum as cls
        elif k == 2:
            from .quadratic_strata import QuadraticStratum as cls
        else:
            cls = Stratum
        return object.__new__(cls)

    def __init__(self, signature, k=1):
        s = sum(signature)
        if s % (2 * k):
            raise ValueError('the sum of orders must be congruent to 2 k = {} got {} (mu={})'.format(2 * k, s, signature))
        if s < - 2 * k:
            raise ValueError('the sum of orders must be at least {} (got {})'.format(- 2 * k, s))
        self._k = k
        self._signature = signature

    def __reduce__(self):
        return Stratum, (self._signature, self._k)

    def surface_differential_order(self):
        r"""
        Return the order of differentials in this stratum.

        EXAMPLES::

            sage: from surface_dynamics import Stratum
            sage: Stratum([1]*6, 3).surface_differential_order()
            3
        """
        return self._k

    # NOTE: called g in diffstrata GeneralisedStratum
    # the name here might be misleading as this is not the genus of the space but of the object in the moduli
    def surface_genus(self):
        r"""
        Return the genus of the surfaces in this stratum.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: Stratum((0,), k=1).surface_genus()
            1
            sage: Stratum((1,1), k=1).surface_genus()
            2
            sage: Stratum((3,2,1), k=1).surface_genus()
            4

            sage: Stratum((-1,-1,-1,-1), k=2).surface_genus()
            0
        """
        return sum(self._signature) // (2 * self._k) + 1

    def genus(self):
        import warnings

        warnings.warn('genus() has been deprecated and will be removed in a future version of surface-dynamics; use surface_genus()')

        return self.surface_genus()

    def surface_has_finite_area(self):
        r"""
        Return whether the k-differentials in this moduli space have finite or
        infinite area.

        EXAMPLES::

            sage: from surface_dynamics import Stratum
            sage: Stratum((3, 2, 1), k=1).surface_has_finite_area()
            True
            sage: Stratum((1, 0, -1), k=1).surface_has_finite_area()
            False

            sage: Stratum([-1]*6, k=3).surface_has_finite_area()
            True
            sage: Stratum([-2]*3, k=3).surface_has_finite_area()
            True
            sage: Stratum([-3]*2, k=3).surface_has_finite_area()
            False
        """
        return self._signature[-1] > -self._k

    def dimension(self):
        r"""
        Return the complex dimension of this stratum.

        The dimension is `2g-2+s+1` where `g` is the genus of surfaces in the
        stratum, `s` the number of singularities. The complex dimension of a
        stratum is also the number of intervals of any interval exchange
        transformations associated to the strata.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: Stratum((0,), k=1).dimension()
            2
            sage: Stratum((0,0), k=1).dimension()
            3
            sage: Stratum((2,), k=1).dimension()
            4
            sage: Stratum((1,1), k=1).dimension()
            5

            sage: Stratum({-1:4}, k=2).dimension()
            2

        ::

            sage: a = Stratum((4,3,2,1,0), k=1)
            sage: p = a.permutation_representative()
            sage: len(p) == a.dimension()
            True
        """
        if self._k == 1 and all(x >= 0 for x in self._signature):
            return 2 * self.surface_genus() + len(self._signature) - 1
        elif self._k == 2 and all(x >= -1 for x in self._signature):
            return 2 * self.surface_genus() + len(self._signature) - 2

        raise NotImplementedError('higher order differential')

    def rank(self):
        r"""
        Return the rank of this GL(2,R)-invariant manifold (half dimension of the
        absolute part of the tangent space).

        EXAMPLES::

            sage: from surface_dynamics import Stratum, QuadraticStrata

            sage: Stratum((0,0), k=1).rank()
            1
            sage: Stratum((2,), k=1).rank()
            2
            sage: Stratum((2,0,0), k=1).rank()
            2

            sage: Stratum({-1: 4}, k=2).rank()
            1
            sage: Stratum({-1:4, 0:5}, k=2).rank()
            1

        Complete list of rank 2 quadratic strata listed by dimension::

            sage: for dim in range(4, 9):
            ....:     quad = [Q for Q in QuadraticStrata(dimension=dim) if Q.rank() == 2]
            ....:     print("%d: %s" % (dim, ", ".join(map(str, quad))))
            4: Q_2(5, -1), Q_1(1^2, -1^2), Q_1(3, -1^3), Q_0(1, -1^5)
            5: Q_3(8), Q_2(2, 1^2), Q_2(4, 1, -1), Q_2(3, 2, -1), Q_2(6, -1^2), Q_1(2, 1, -1^3), Q_1(4, -1^4), Q_0(2, -1^6)
            6: Q_3(6, 2), Q_3(4^2), Q_2(2^2, 1, -1), Q_2(4, 2, -1^2), Q_1(2^2, -1^4)
            7: Q_3(4, 2^2), Q_2(2^3, -1^2)
            8: Q_3(2^4)
        """
        if self._k == 1:
            if all(x >= 0 for x in self._signature):
                return self.surface_genus()
            else:
                # is there a meaning for rank?
                raise NotImplementedError
        elif self._k == 2:
            if all(x >= -1 for x in self._signature):
                return self.surface_genus() + sum(z % 2 for z in self._signature) // 2 - 1
            else:
                # is there a meaning for rank?
                raise NotImplementedError
        elif self._k >= 3:
            raise ValueError('not a GL(2,R)-invariant subvariety')

    def signature(self):
        r"""
        Return the order of zeros with multiplicities.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: Stratum([1, 2, 3], k=1).signature()
            (3, 2, 1)
            sage: Stratum({2: 4}, k=1).signature()
            (2, 2, 2, 2)
            sage: Stratum([-1, 1], k=1).signature()
            (1, -1)

            sage: Stratum({-1: 4}, k=2).signature()
            (-1, -1, -1, -1)
            sage: Stratum({1: 8}, k=2).signature()
            (1, 1, 1, 1, 1, 1, 1, 1)
            sage: Stratum({-2: 2, 0: 1}, k=2).signature()
            (0, -2, -2)
        """
        return self._signature

    def zeros(self, fake_zeros=True, poles=True):
        import warnings

        warnings.warn('zeros() has been deprecated and will be removed in a future version of surface-dynamics; use signature()')

        s = self.signature()
        if not fake_zeros:
            s = tuple(m for m in s if m)
        if not poles:
            s = tuple(m for m in s if m >= 0)
        return s

    def nb_zeros(self, fake_zeros=True, poles=True):
        r"""
        Returns the number of zeros of self.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: Stratum([0], k=1).nb_zeros()
            doctest:warning
            ...
            UserWarning: nb_zero() has been deprecated and will be removed in a future version of surface-dynamics; use signature()
            1
            sage: Stratum({2:4,3:2}, k=1).nb_zeros()
            6

            sage: Stratum({-1:4}, k=2).nb_zeros()
            4
            sage: Stratum({-1:4,1:4}, k=2).nb_zeros()
            8
        """
        import warnings

        warnings.warn('nb_zero() has been deprecated and will be removed in a future version of surface-dynamics; use signature()')

        return sum(int((fake_zeros or m) and (poles or m >= 0)) for m in self.signature())

    def nb_fake_zeros(self):
        r"""
        Return the number of fake zeros.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: Stratum([0], k=1).nb_fake_zeros()
            doctest:warning
            ...
            UserWarning: nb_fake_zeros() has been deprecated and will be removed in a future version of surface-dynamics; use signature()
            1
            sage: Stratum([1,1,0,0], k=1).nb_fake_zeros()
            2

            sage: Stratum([0,4,2,2], k=2).nb_fake_zeros()
            1
        """
        import warnings

        warnings.warn('nb_fake_zeros() has been deprecated and will be removed in a future version of surface-dynamics; use signature()')

        return self._signature.count(0)

    def nb_poles(self):
        r"""
        Return the number of poles of this quadratic stratum.
        """
        import warnings

        warnings.warn('nb_poles() has been deprecated and will be removed in a future version of surface-dynamics; use signature()')

        return sum(int(m < 0) for m in self.signature())

    @cached_method
    def connected_components(self):
        r"""
        Return the connected components of this stratum of differentials.

        - Abelian holomorphic differentials [KonZor03]_
        - Quadratic differentials with at most simple poles [Lan08]_
        """
        zeros = tuple(m for m in self._signature if m)

        if not self.surface_has_finite_area():
            raise NotImplementedError('meromorphic differentials with higher order poles')

        if self._k == 1:
            # Abelian differentials: Kontsevich-Zorich classification
            from .abelian_strata import ASC, HypASC, NonHypASC, OddASC, EvenASC
            s = sum(zeros)
            genus = s // 2 + 1

            if genus == 1:
                return (HypASC(self),)
            elif genus == 2:
                return (HypASC(self),)
            elif genus == 3:
                if zeros == (2, 2) or zeros == (4,):
                    return (HypASC(self), OddASC(self))
                else:
                    return (ASC(self),)
            elif len(zeros) == 1:
                    # just one zeros [2g-2]
                    return (HypASC(self), OddASC(self), EvenASC(self))
            elif zeros == (genus-1, genus-1):
                # two identical zeros [g-1, g-1]
                if genus % 2 == 0:
                    return (HypASC(self), NonHypASC(self))
                else:
                    return (HypASC(self), OddASC(self), EvenASC(self))
            elif all(x%2 == 0 for x in zeros):
                # even zeroes [2 l_1, 2 l_2, ..., 2 l_n]
                return (OddASC(self), EvenASC(self))
            else:
                return (ASC(self),)

        elif self._k == 2:
            # quadratic differentials: Lanneau classification
            from .quadratic_strata import CQSC, HQSC, GTNQSC, GTHQSC, GOQSC, REQSC, IEQSC, GZQSC, NQSC
            nb_poles = zeros.count(-1)
            s = sum(zeros)
            genus = s // 4 + 1

            #TODO: check genus 2 components
            #TODO: in genus 2, decide between GTHQSC/GTNQSC and HQSC/NQSC
            if genus == 0:
                return (GZQSC(self),)

            #TODO: all genus 1 strata are connected, but two are hyperelliptic; give the component a different name then?
            elif genus == 1:
                if zeros == () or zeros == (1, -1):
                    # empty!
                    return ()
                else:
                    return  (GOQSC(self),)

            elif genus == 2:
                if zeros == (4,) or zeros == (3, 1):
                    # empty!
                    return ()
                elif zeros == (6, -1, -1):
                    return (HQSC(self), GTNQSC(self))
                elif zeros == (3, 3, -1, -1):
                    return (HQSC(self), GTNQSC(self))
                elif zeros == (2, 2) or zeros == (2, 1, 1) or zeros == (1, 1, 1, 1):
                    return (GTHQSC(self),)
                else:
                    return (GTNQSC(self),)

            elif genus == 3 and nb_poles == 1 and all((z == -1 or z % 3 == 0) for z in zeros):
                return (REQSC(self), IEQSC(self))

            elif genus == 4 and nb_poles == 0 and all(z % 3 == 0 for z in zeros):
                if zeros == (12,) or zeros == (9, 3):
                    return (REQSC(self), IEQSC(self))
                elif zeros == (6, 6) or zeros == (6, 3, 3) or zeros == (3, 3, 3, 3):
                    return (HQSC(self), REQSC(self), IEQSC(self))
            else:
                if len(zeros) == 2 and zeros[0] % 4 == 2 and zeros[1] % 4 == 2:
                    return (HQSC(self), NQSC(self))
                elif len(zeros) == 4 and zeros[0] == zeros[1] and zeros[2] == zeros[3] and zeros[0] % 2 and zeros[2] % 2:
                    return (HQSC(self), NQSC(self))
                elif len(zeros) == 3 and zeros[0] == zeros[1] and zeros[0] % 2 and zeros[2] % 4 == 2:
                    return (HQSC(self), NQSC(self))
                elif len(zeros) == 3 and zeros[1] == zeros[2] and zeros[1] % 2 and zeros[0] % 4 == 2:
                    return (HQSC(self), NQSC(self))
                else:
                    return (CQSC(self),)

        raise NotImplementedError('higher order differential')

    #
    # String representation
    #

    def _flat_zero_str(self):
        r"""
        String representation of the zeros.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: a = Stratum({2:3}, k=1)
            sage: a._flat_zero_str()
            '2, 2, 2'
        """
        return ', '.join(map(str, self.signature()))

    def _exp_zero_str(self):
        r"""
        String representation with exponential notation

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: a = Stratum({2: 3}, k=1)
            sage: a._exp_zero_str()
            '2^3'
        """
        return ', '.join('%d^%d' %(i,e) if e != 1 else '%d' %i for (i,e) in list_to_exp_list(self.signature()))

    # this attribute can be switched between _flat_zero_str and _exp_zero_str
    _zero_str = _exp_zero_str

    def _repr_(self):
        """
        TESTS::

            sage: from surface_dynamics import Stratum

            sage: repr(Stratum([1,1], k=1))       # indirect doctest
            'H_2(1^2)'
            sage: repr(Stratum([1,1,1,1], k=2)) # indirect doctest
            'Q_2(1^4)'
        """
        try:
            name = self._name
        except AttributeError:
            if self._k == 1:
                name = 'H'
            elif self._k == 2:
                name = 'Q'
            else:
                name = 'H^{{({})}}'.format(self._k)

        return name + "_" + str(self.surface_genus()) + "(" + self._zero_str() + ")"

    def _latex_(self):
        r"""
        Latex string representation

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: Stratum([0], k=1)._latex_()
            '\\mathcal{H}_1(0)'
            sage: Stratum({-1:4}, k=2)._latex_()
            '\\mathcal{Q}_0(-1^4)'
        """
        return self._latex_name + '_' + str(self.surface_genus()) + "(" + self._zero_str() + ")"

    #
    # Equality and comparisons
    #

    def __lt__(self, other):
        r"""
        Comparison

        ALGORITHM:

        First compare the class, then the dimension and then the list of zeros.

        TESTS::

            sage: from surface_dynamics import Stratum

            sage: Stratum([0], k=1) == Stratum([0], k=1)
            True
            sage: Stratum([5,-1], k=2) == Stratum([5,-1], k=2)
            True
            sage: Stratum([5,-1], k=2) != Stratum([5,-1], k=2)
            False

            sage: Stratum([12], k=1) == Stratum([12], k=2)
            False
            sage: Stratum([12,0,0], k=2) == Stratum([12,0], k=2)
            False

            sage: Stratum([2,0], k=1) == Stratum([2], k=1)
            False
            sage: Stratum([2,0], k=1) != Stratum([2], k=1)
            True

            sage: Stratum([1,1], k=1) < Stratum([1,1,0], k=1)
            True
            sage: Stratum([1,1,0], k=1) < Stratum([1,1], k=1)
            False
            sage: Stratum([1,1,0], k=1) < Stratum([1,1,0,0], k=1)
            True
            sage: Stratum([2], k=1) < Stratum([1,1], k=1)
            True
            sage: Stratum([4,0], k=1) > Stratum([1,1,1,1], k=1)
            False
            sage: Stratum([4,0,0,0], k=1) > Stratum([1,1,1,1], k=1)
            True

        ::

            sage: Stratum([2,2], k=2) < Stratum([2,2,0], k=2)
            True
            sage: Stratum([2,2,0], k=2) < Stratum([2,2], k=2)
            False
            sage: Stratum([2,2,0], k=2) < Stratum([2,2,0,0], k=2)
            True
            sage: Stratum([4], k=2) < Stratum([2,2], k=2)
            True
            sage: Stratum([4,0], k=2) > Stratum([1,1,1,1], k=2)
            False

            sage: Q1 = Stratum([4,0,0,0], k=2)
            sage: Q2 = Stratum([1,1,1,1], k=2)
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

        if self._k == other._k:
            # compare the dimension
            if self.dimension() < other.dimension():
                return True
            elif self.dimension() > other.dimension():
                return False

            # compare the list of zeros
            if self.signature() < other.signature():
                return True
            elif self.signature() > other.signature():
                return False

            # equality
            return False

        else:
            return self._k < other._k

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

            sage: from surface_dynamics import Stratum

            sage: Stratum([2], k=1).is_connected()
            True
            sage: Stratum([2,2], k=1).is_connected()
            False

            sage: Stratum([-1,-1,-1,-1], k=2).is_connected()
            True
            sage: Stratum([12], k=2).is_connected()
            False
        """
        return len(self.connected_components()) <= 1

    def permutation_representative(self, *args, **kwds):
        r"""
        Return a permutation of interval exchanges associated to this stratum.

        This method only makes sense for Abelian and quadratic differentials.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

        Examples from Abelian differentials::

            sage: a = Stratum([3,2,1,0,0], k=1)
            sage: p = a.permutation_representative()
            sage: p.stratum()
            H_4(3, 2, 1, 0^2)
            sage: a = Stratum([2, 2, 2], k=1)
            sage: p = a.permutation_representative()
            sage: p.stratum()
            H_4(2^3)

        Examples from quadratic differentials::

            sage: a = Stratum([6,-1,-1], k=2)
            sage: p = a.permutation_representative()
            sage: p.stratum()
            Q_2(6, -1^2)
            sage: a = Stratum([-1,-1,-1,-1,0,0], k=2)
            sage: p = a.permutation_representative()
            sage: p.stratum()
            Q_0(0^2, -1^4)
        """
        if self._k > 2:
            raise ValueError('stratum of higher order differentials')

        return self.one_component().permutation_representative(*args,**kwds)

    def is_empty(self):
        r"""
        Return True if the stratum is empty

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: Stratum([2], k=1).is_empty()
            False
            sage: Stratum([1,-1], k=2).is_empty()
            True
        """
        return len(self.connected_components()) == 0

    def number_of_components(self):
        r"""
        Returns the number of connected components of self

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: Stratum([2], k=1).number_of_components()
            1
            sage: Stratum([4], k=1).number_of_components()
            2
            sage: Stratum([3,3], k=1).number_of_components()
            2
        """
        return len(self.connected_components())

    def one_component(self):
        r"""
        Returns a connected component of this stratum.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: Stratum([2], k=1).one_component()
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

            sage: from surface_dynamics import Stratum

            sage: a = Stratum([1,1], k=1)
            sage: a
            H_2(1^2)
            sage: a.unique_component()
            H_2(1^2)^hyp

            sage: a = Stratum([3,2,1], k=1)
            sage: a
            H_4(3, 2, 1)
            sage: a.unique_component()
            H_4(3, 2, 1)^c

            sage: Stratum({1:1, -1:5}, k=2).unique_component()
            Q_0(1, -1^5)^c
            sage: Stratum([3,2,-1], k=2).unique_component()
            Q_2(3, 2, -1)^nonhyp

            sage: Stratum([12], k=2).unique_component()
            Traceback (most recent call last):
            ...
            ValueError: several components for this stratum
        """
        ccs = self.connected_components()
        if len(ccs) != 1:
            raise ValueError("several components for this stratum")
        return ccs[0]

    def random_component(self):
        r"""
        Returns a random connected component of this stratum.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: Q = Stratum([6,6], k=2)
            sage: Q.random_component() # random
            Q_4(6^2)^hyp
            sage: Q.random_component() # random
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

            sage: from surface_dynamics import Stratum

        Some abelian strata::

            sage: Stratum([0], k=1).components()
            (H_1(0)^hyp,)
            sage: Stratum([2], k=1).components()
            (H_2(2)^hyp,)
            sage: Stratum([4], k=1).components()
            (H_3(4)^hyp, H_3(4)^odd)
            sage: Stratum([2,2], k=1).components()
            (H_3(2^2)^hyp, H_3(2^2)^odd)
            sage: Stratum([1,1,1,1], k=1).components()
            (H_3(1^4)^c,)

        Some quadratic strata::

            sage: Stratum([12], k=2).components()
            (Q_4(12)^reg, Q_4(12)^irr)
            sage: Stratum([6,-1,-1], k=2).components()
            (Q_2(6, -1^2)^hyp, Q_2(6, -1^2)^nonhyp)
        """
        # TODO: deprecate
        return self.connected_components()

    def masur_veech_volume(self, rational=False, method=None):
        r"""
        Return the Masur-Veech volume of this stratum.

        INPUT:

        - ``rational`` (optional, boolean) - if ``False`` (default) return the Masur-Veech volume
          and if ``True`` return the Masur-Veech volume divided by `\zeta(2g)`.

        - ``method`` (optional string) - the method to use to compute the volume either, see
          :func:`~surface_dynamics.flat_surfaces.masur_veech_volumes.masur_veech_volume`

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: Stratum([2], k=1).masur_veech_volume()
            1/120*pi^4
            sage: Stratum([1,1,1,1], k=1).masur_veech_volume()
            1/4860*pi^6
            sage: Stratum([20], k=1).masur_veech_volume()
            1604064377302075061983/792184445986404135075840000000000*pi^22
        """
        from .masur_veech_volumes import masur_veech_volume
        return masur_veech_volume(self, rational, method)


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
        TESTS::

            sage: from surface_dynamics import Stratum

            sage: a = Stratum([4,4], k=1).one_component()
            sage: a == loads(dumps(a))
            True
            sage: q = Stratum([5,5,-1,-1], k=2).one_component()
            sage: q == loads(dumps(q))
            True
        """
        self._stratum = stratum

    def __reduce__(self):
        r"""
        Reduce method for pickling

        TESTS::

            sage: from surface_dynamics import Stratum

        Tests for Abelian strata::

            sage: a = Stratum([2,2], k=1)
            sage: all(loads(dumps(cc)) == cc for cc in a.components())
            True
            sage: a = Stratum([3,3], k=1)
            sage: all(loads(dumps(cc)) == cc for cc in a.components())
            True
            sage: a = Stratum([6], k=1)
            sage: all(loads(dumps(cc)) == cc for cc in a.components())
            True
            sage: a = Stratum([1,1,1,1], k=1)
            sage: all(loads(dumps(cc)) == cc for cc in a.components())
            True

        Tests for quadratic strata::

            sage: q = Stratum([-1,-1,-1,-1], k=2)
            sage: all(loads(dumps(cc)) == cc for cc in q.components())
            True
            sage: q = Stratum([12], k=2)
            sage: all(loads(dumps(cc)) == cc for cc in q.components())
            True
        """
        return self.__class__, (self._stratum,)

    def _repr_(self):
        r"""
        String representation

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: a_hyp = Stratum([4], k=1).hyperelliptic_component()
            sage: a_hyp._repr_()
            'H_3(4)^hyp'
            sage: a_odd = Stratum([4], k=1).odd_component()
            sage: a_odd._repr_()
            'H_3(4)^odd'
        """
        return str(self._stratum) + "^" + self._name

    def __hash__(self):
        r"""
        TESTS::

            sage: from surface_dynamics import Stratum

            sage: A4hyp = Stratum([4], k=1).hyperelliptic_component()
            sage: A4odd = Stratum([4], k=1).odd_component()
            sage: hash(A4hyp) != hash(A4odd)
            True
            sage: A22hyp = Stratum([2,2], k=1).hyperelliptic_component()
            sage: hash(A22hyp) != hash(A4hyp)
            True
        """
        return hash(self._stratum) ^ hash(self._name)

    def stratum(self):
        r"""
        Return the stratum associated to self

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: a = Stratum([4,4], k=1)
            sage: all([c.stratum() == a for c in a.components()])
            True
        """
        return self._stratum

    def surface_genus(self):
        r"""
        Return genus of the corresponding stratum

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: a = Stratum([4,4], k=1)
            sage: a.one_component().surface_genus()
            5
        """
        return self._stratum.surface_genus()

    def surface_has_finite_area(self):
        return self._stratum.surface_has_finite_area()

    def surface_differential_order(self):
        return self._stratum.surface_differential_order()

    def genus(self):
        import warnings

        warnings.warn('genus() has been deprecated and will be removed in a future version of surface-dynamics; use surface_genus()')

        return self.surface_genus()

    def dimension(self):
        r"""
        Return the (complex) dimension of this GL(2,R)-invariant orbifold.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: Stratum([4], k=1).odd_component().dimension()
            6
            sage: Stratum([12], k=2).regular_component().dimension()
            7
        """
        return self._stratum.dimension()

    def rank(self):
        r"""
        Return the rank of this GL(2,R)-invariant orbifold.

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: Stratum([4], k=1).odd_component().rank()
            3
            sage: Stratum([12], k=2).regular_component().rank()
            3
        """
        return self._stratum.rank()

    def __eq__(self,other):
        r"""
        Equality test

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: c_hyp = Stratum([6], k=1).hyperelliptic_component()
            sage: c_odd = Stratum([6], k=1).odd_component()
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

            sage: from surface_dynamics import Stratum

            sage: a1 = Stratum([1,1,1,1], k=1)
            sage: c1 = a1.components()[0]
            sage: a2 = Stratum([3,1], k=1)
            sage: c2 = a2.components()[0]
            sage: c1 == c1
            True
            sage: c1 == c2
            False
            sage: a1 = Stratum([1,1,1,1], k=1)
            sage: c1 = a1.components()[0]
            sage: a2 = Stratum([2, 2], k=1)
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

    def masur_veech_volume(self, rational=False, method=None):
        r"""
        Return the Masur-Veech volume of this stratum component.

        INPUT:

        - ``rational`` (optional, boolean) - if ``False`` (default) return the Masur-Veech volume
          and if ``True`` return the Masur-Veech volume divided by `\zeta(2g)`.

        - ``method`` (optional string) - the method to use to compute the volume either, see
          :func:`~surface_dynamics.flat_surfaces.masur_veech_volumes.masur_veech_volume`

        EXAMPLES::

            sage: from surface_dynamics import Stratum

            sage: Stratum([4], k=1).hyperelliptic_component().masur_veech_volume()
            1/6720*pi^6
            sage: Stratum([6], k=1).even_component().masur_veech_volume()
            32/1913625*pi^8
        """
        from .masur_veech_volumes import masur_veech_volume
        return masur_veech_volume(self, rational, method)

#
# Strata (family of strata)
#

class Strata(Parent):
    r"""
    Strata of Abelian or Quadratic differentials.
    """
    pass
