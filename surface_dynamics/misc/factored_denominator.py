r"""
Factored denominators and associated free modules
"""
# ****************************************************************************
#       Copyright (C) 2019 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import absolute_import, print_function
from six.moves import range, zip

from six import iteritems

from sage.misc.cachefunc import cached_method
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.richcmp import op_LT, op_EQ, op_NE, op_LE, op_GE, op_GT
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.all import Rings
from sage.rings.all import ZZ, QQ
from sage.modules.vector_integer_dense import Vector_integer_dense


def laurent_monomial(R, arg):
    r"""
    EXAMPLES::

        sage: from surface_dynamics.misc.factored_denominator import laurent_monomial
        sage: R = LaurentPolynomialRing(QQ, 'x', 2)
        sage: laurent_monomial(R, (1,2))
        x0*x1^2
        sage: laurent_monomial(R, vector((-1r,3r)))
        x0^-1*x1^3
        sage: laurent_monomial(R, (-1,3)) * laurent_monomial(R, (1,1)) == laurent_monomial(R, (0,4))
        True
    """
    # only in beta 8.1 versions of Sage
    # return R.monomial(*arg)
    arg = tuple(arg)
    if len(arg) != R.ngens():
        raise TypeError("tuple key must have same length as ngens")

    from sage.misc.misc_c import prod
    return prod(x**int(i) for (x, i) in zip(R.gens(), arg))


def vector_to_monomial_string(u, var_names):
    r"""
    EXAMPLES::

        sage: from surface_dynamics.misc.factored_denominator import vector_to_monomial_string
        sage: vector_to_monomial_string((0,4,3), 'x')
        'x1^4*x2^3'
        sage: vector_to_monomial_string((0,4,3), 'hello')
        'hello1^4*hello2^3'
        sage: vector_to_monomial_string((0,4,3), ('x','y','z'))
        'y^4*z^3'
        sage: vector_to_monomial_string((1,0,-1), 'x')
        'x0*x2^-1'
    """
    s = []
    for i, j in enumerate(u):
        if j:
            if isinstance(var_names, str):
                var = '%s%d' % (var_names, i)
            else:
                var = var_names[i]

            if j == 1:
                s.append(var)
            else:
                s.append('%s^%d' % (var, j))

    return '*'.join(s) if s else '1'


def vector_to_linear_form_string(u, var_names):
    r"""
    EXAMPLES::

        sage: from surface_dynamics.misc.factored_denominator import vector_to_linear_form_string
        sage: vector_to_linear_form_string((0,4,3), 'x')
        '4*x1 + 3*x2'
        sage: vector_to_linear_form_string((0,4,3), 'hello')
        '4*hello1 + 3*hello2'
        sage: vector_to_linear_form_string((0,4,3), ('x','y','z'))
        '4*y + 3*z'
        sage: vector_to_linear_form_string((1,0,-1), 'x')
        'x0 - x2'
        sage: vector_to_linear_form_string((2,0,-3), 'x')
        '2*x0 - 3*x2'
        sage: vector_to_linear_form_string((-1,1,0), 'x')
        '-x0 + x1'
        sage: vector_to_linear_form_string((-2,3,0), 'x')
         '-2*x0 + 3*x1'
        sage: vector_to_linear_form_string((0,-1,1,0), 'x')
        '-x1 + x2'
    """
    s = ''
    first = True
    for i, j in enumerate(u):
        if j:
            if isinstance(var_names, str):
                var = '%s%d' % (var_names, i)
            else:
                var = var_names[i]

            if j == 1:
                if first:
                    s += var
                else:
                    s += ' + %s' % var
            elif j == -1:
                if first:
                    s += '-%s' % var
                else:
                    s += ' - %s' % var
            elif first:
                s += '%d*%s' % (j, var)
            elif j < 0:
                s += ' - %d*%s' % (-j, var)
            else:
                s += ' + %d*%s' % (j, var)
            first = False

    return '0' if not s else s


# NOTE: should this be an instance of Factorization?
class FactoredDenominator(object):
    r"""
    Factored denominator

    This class is a simple datastructure to handle factored denominator, that is
    a list of pairs ``(v, m)`` where the ``v`` are integer vectors (of fixed
    dimension) and ``m`` are multiplicities (i.e. positive integers).

    It is used for at least two purposes:

    - (Factored) product of polynomials of the form `(1 - m)^d` where `m` is a
      monomial

    - generalized multiple zeta values, where denominators are products of
      linear forms

    EXAMPLES::

        sage: from surface_dynamics.misc.factored_denominator import FactoredDenominator

        sage: V = ZZ**3
        sage: f1 = FactoredDenominator([((1,0,0), 2)], V)
        sage: f2 = FactoredDenominator([((0,1,2), 3), ((1,1,1), 1)], V)
        sage: f3 = FactoredDenominator([((0,-1,2), 1), ((1,0,0), 1), ((0,0,2), 1)], V)
        sage: f1
        {(1, 0, 0): 2}
        sage: f1 * f2 * f3
        {(1, 0, 0): 3, (0, 1, 2): 3, (1, 1, 1): 1, (0, -1, 2): 1, (0, 0, 2): 1}
        sage: hash(f1)  # random
        9164823176457928306

        sage: FactoredDenominator(f1)
        {(1, 0, 0): 2}
    """
    __slots__ = ['_dict', '_tuple']

    def __init__(self, data, V=None):
        if isinstance(data, (tuple, list)):
            if V is None:
                raise ValueError("a ZZ-free module V must be provided")
            d = V.dimension()
            self._dict = {}
            for mon, mult in data:
                mon = V(mon)
                if mon.is_zero():
                    raise ValueError('zero in denominator')
                for i in range(d):
                    if mon[i]:
                        break
                mult = ZZ(mult)
                mon.set_immutable()
                self._dict[mon] = mult

        elif isinstance(data, dict):
            if V is not None:
                self._dict = {}
                for k, v in data.items():
                    k = V(k)
                    k.set_immutable()
                    self._dict[k] = v
            else:
                self._dict = data

        elif isinstance(data, FactoredDenominator):
            self._dict = data._dict.copy()
            self._tuple = data._tuple[:]
            return

        elif isinstance(data, Vector_integer_dense):
            data = V(data)
            data.set_immutable()
            self._dict = {data: ZZ.one()}
            self._tuple = ((data, ZZ.one()),)
            return

        else:
            raise TypeError('invalid data of type {} to initialized a FactoredDenominator'.format(type(data)))

        self._tuple = tuple(sorted(self._dict.items()))

    def __len__(self):
        r"""
        Return the number of factors (without multiplicities).
        """
        return len(self._tuple)

    def __iter__(self):
        r"""
        Iterates through the pairs ``(monomial, exponent)``.
        """
        return iter(self._tuple)

    def subs(self, m):
        r"""
        Matrix substitution.

        EXAMPLES::

            sage: from surface_dynamics.misc.factored_denominator import FactoredDenominator
            sage: V = ZZ**3
            sage: f = FactoredDenominator([((0,1,2), 3), ((1,1,1), 1)], V)

            sage: m = matrix(3, [1,1,0,0,1,1,1,0,1])
            sage: f.subs(m) == FactoredDenominator([((1,3,2), 3), ((2,2,2), 1)], V)
            True

            sage: m = matrix(3, [1,1,1,1,1,1,1,1,1])
            sage: f.subs(m)
            {(3, 3, 3): 4}

            sage: m = matrix(3, [1,0,0,1,0,0,1,0,0])
            sage: f.subs(m)
            Traceback (most recent call last):
            ...
            ValueError: zero denominator
        """
        if not self._tuple:
            return self
        new_dict = {}
        for mon, mult in self._dict.items():
            new_mon = m * mon
            new_mon.set_immutable()
            if not new_mon:
                raise ValueError('zero denominator')
            if new_mon in new_dict:
                new_dict[new_mon] += mult
            else:
                new_dict[new_mon] = mult

        return FactoredDenominator(new_dict)

    def to_multiplicative_polynomial(self, S, extra_var=False):
        r"""
        Return the product of the term in a given polynomial ring ``S``.

        EXAMPLES::

            sage: from surface_dynamics.misc.factored_denominator import FactoredDenominator

            sage: V = ZZ ** 3
            sage: f = FactoredDenominator([((1,0,0), 2)], V)
            sage: g = FactoredDenominator([((1,0,0), 2), ((1,1,1),1)], V)

            sage: R1 = QQ['x,y,z']
            sage: f.to_multiplicative_polynomial(R1)
            x^2 - 2*x + 1
            sage: g.to_multiplicative_polynomial(R1)
            -x^3*y*z + 2*x^2*y*z - x*y*z + x^2 - 2*x + 1

            sage: f.to_multiplicative_polynomial(R1['E'], extra_var=True)
            x^2*E^2 - 2*x*E + 1
            sage: g.to_multiplicative_polynomial(R1['E'], extra_var=True)
            -x^3*y*z*E^5 + 2*x^2*y*z*E^4 - x*y*z*E^3 + x^2*E^2 - 2*x*E + 1

            sage: g.to_multiplicative_polynomial(QQ['t1,t2,t3'])
            -t1^3*t2*t3 + 2*t1^2*t2*t3 - t1*t2*t3 + t1^2 - 2*t1 + 1
        """
        if extra_var:
            R = S.base_ring()
            T = S.gen()
        else:
            R = S
            T = S.one()

        ans = S.one()
        for a, i in self._tuple:
            ans *= (S.one() - R.monomial(*a) * T ** sum(a)) ** i
        return ans

    def to_additive_polynomial(self, S, extra_var=False):
        r"""
        Return the product of the term in a given polynomial ring ``S``.

        EXAMPLES::

            sage: from surface_dynamics.misc.factored_denominator import FactoredDenominator

            sage: V = ZZ ** 3
            sage: f = FactoredDenominator([((1,0,0), 2)], V)
            sage: g = FactoredDenominator([((1,0,0), 2), ((1,1,1),1)], V)

            sage: R1 = QQ['x,y,z']
            sage: f.to_additive_polynomial(R1)
            x^2
            sage: g.to_additive_polynomial(R1)
            x^3 + x^2*y + x^2*z

            sage: f.to_additive_polynomial(R1['E'], extra_var=True)
            x^2*E^2
            sage: g.to_additive_polynomial(R1['E'], extra_var=True)
            (x^3 + x^2*y + x^2*z)*E^3

            sage: g.to_additive_polynomial(QQ['t1,t2,t3'])
            t1^3 + t1^2*t2 + t1^2*t3
        """
        if extra_var:
            R = S.base_ring()
            T = S.gen()
        else:
            R = S
            T = S.one()

        ans = S.one()
        for a, mult in self._tuple:
            monomial = sum(coeff * R.gen(i) for i, coeff in enumerate(a))
            ans *= (monomial * T) ** mult
        return ans

    def degree(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.misc.factored_denominator import FactoredDenominator

            sage: V = ZZ**3
            sage: FactoredDenominator([], V).degree()
            0
            sage: FactoredDenominator([((1,0,0), 2)], V).degree()
            2
            sage: FactoredDenominator([((0,1,2), 3), ((1,1,1), 1)], V).degree()
            4
            sage: FactoredDenominator([((0,-1,2), 1), ((1,0,0), 1), ((0,0,2), 1)], V).degree()
            3

        TESTS::

            sage: parent(FactoredDenominator([], V).degree())
            Integer Ring
            sage: parent(FactoredDenominator([((1,0,0), 2)], V).degree())
            Integer Ring
        """
        if not self._tuple:
            return ZZ.zero()
        return sum(m for _, m in self._tuple)

    # TODO
    # this method does not really make sense at this level of generality!
    # this can be defined on denominators and extended by linearity
    # we should just provide a method to do this in an abstract parent algebra
    #
    # _image_from_image_of_denominators(function, elt)
    #   return sum(num * function(den) for num,den in elt)
    def logarithmic_minus_derivative(self, j):
        r"""
        Minus the logarithmic derivative -v' / v taken with respect to the ``j``-th variable

        Since this denominator has the form 1 / (1 - m1) (1 - m2) ... (1 - mk) the
        logarithmic derivative is simply

        m1'/(1-m1) + m2'/(1-m2) + ... + mk'/(1-mk)

        INPUT:

        - ``j`` -- index of a variable

        OUTPUT: iterator of triples ``(mult, v, monomial)``

        EXAMPLES::

            sage: from surface_dynamics.misc.factored_denominator import FactoredDenominator

            sage: V = ZZ**3

            sage: f = FactoredDenominator([((1,0,0), 2)], V)
            sage: f
            {(1, 0, 0): 2}
            sage: list(f.logarithmic_minus_derivative(0))
            [(2, (0, 0, 0), (1, 0, 0))]
            sage: list(f.logarithmic_minus_derivative(1))
            []
            sage: list(f.logarithmic_minus_derivative(2))
            []

            sage: f = FactoredDenominator([((1,1,1), 1)], V)
            sage: list(f.logarithmic_minus_derivative(0))
            [(1, (0, 1, 1), (1, 1, 1))]

            sage: f = FactoredDenominator([((1,0,0), 1), ((0,1,0), 1), ((0,0,1), 1)], V)
            sage: list(f.logarithmic_minus_derivative(0))
            [(1, (0, 0, 0), (1, 0, 0))]

            sage: f = FactoredDenominator([((1,0,0), 2), ((1,1,0), 3), ((1,1,1), 1)], V)
            sage: f
            {(1, 0, 0): 2, (1, 1, 0): 3, (1, 1, 1): 1}
            sage: list(f.logarithmic_minus_derivative(0))
            [(2, (0, 0, 0), (1, 0, 0)),
             (3, (0, 1, 0), (1, 1, 0)),
             (1, (0, 1, 1), (1, 1, 1))]
        """
        for mon, mult in self._tuple:
            if mon[j]:
                v = mon.__copy__()
                mult *= v[j]
                v[j] -= 1
                yield (mult, v, mon)

    def __eq__(self, other):
        r"""
        Equality test

        EXAMPLES::

            sage: from surface_dynamics.misc.factored_denominator import FactoredDenominator

            sage: V = ZZ**3

            sage: f1 = FactoredDenominator([((1,0,0), 2)], V)
            sage: f2 = FactoredDenominator([((1,0,0), 2)], V)
            sage: f3 = FactoredDenominator([((1,0,0), 1)], V)
            sage: f4 = FactoredDenominator([((1,0,1), 2)], V)

            sage: f1 == f1 and f1 == f2
            True
            sage: f1 == f3 or f1 == f4
            False
        """
        if type(self) is not type(other):
            raise TypeError
        return self._tuple == other._tuple

    def __ne__(self, other):
        if type(self) is not type(other):
            raise TypeError
        return self._tuple != other._tuple

    def __lt__(self, other):
        if type(self) is not type(other):
            raise TypeError
        return self._tuple < other._tuple

    def __le__(self, other):
        if type(self) is not type(other):
            raise TypeError
        return self._tuple <= other._tuple

    def __gt__(self, other):
        if type(self) is not type(other):
            raise TypeError
        return self._tuple > other._tuple

    def __ge__(self, other):
        if type(self) is not type(other):
            raise TypeError
        return self._tuple >= other._tuple

    def copy(self):
        res = FactoredDenominator.__new__(FactoredDenominator)
        res._dict = self._dict.copy()
        res._tuple = self._tuple
        return res

    def __mul__(self, other):
        r"""
        Multiplication

        TESTS::

            sage: from surface_dynamics.misc.factored_denominator import FactoredDenominator
            sage: V = ZZ**3

            sage: f = FactoredDenominator([((1,0,0), 2)], V)
        """
        if type(self) is not type(other):
            raise TypeError

        new_data = self._dict.copy()
        for i, j in other._dict.items():
            if i in new_data:
                new_data[i] += j
                if new_data[i].is_zero():
                    del new_data[i]
            else:
                new_data[i] = j

        return FactoredDenominator(new_data, None)

    def __truediv__(self, other):
        r"""
        (Partial) division.

        EXAMPLES::

            sage: from surface_dynamics.misc.factored_denominator import FactoredDenominator
            sage: V = ZZ**3

            sage: f1 = FactoredDenominator([((1,0,0), 3)], V)
            sage: f2 = FactoredDenominator([((1,0,0), 1)], V)
            sage: f1/f2
            {(1, 0, 0): 2}
        """
        if type(self) != type(other):
            raise TypeError

        sd = self._dict
        od = other._dict
        nd = sd.copy()
        for i, j in od.items():
            jj = sd.get(i, -1)
            if j > jj:
                raise ArithmeticError
            elif j == jj:
                del nd[i]
            else:
                nd[i] -= j

        return FactoredDenominator(nd, None)

    __div__ = __truediv__

    def __nonzero__(self):
        return bool(self._dict)

    def is_one(self):
        return not self._dict

    def str(self):
        return str(self._dict)

    def __repr__(self):
        return repr(self._dict)

    def __hash__(self):
        return hash(self._tuple)

    def lcm_update(self, other):
        if type(self) is not type(other):
            raise TypeError

        sd = self._dict
        od = other._dict
        for i, j in od.items():
            if j > sd.get(i, -1):
                sd[i] = j

        self._tuple = tuple(sorted(self._dict.items()))

    def lcm(self, other):
        r"""
        Return the lcm of two factored denominators.

        EXAMPLES::

            sage: from surface_dynamics.misc.factored_denominator import FactoredDenominator

            sage: V = ZZ**3

            sage: f1 = FactoredDenominator([((1,0,0), 2),((0,1,0),1)], V)
            sage: f2 = FactoredDenominator([((1,0,0), 1),((0,1,0),2)], V)
            sage: f3 = FactoredDenominator([((1,0,0), 1),((0,1,0),1),((0,0,1),1)], V)
            sage: f4 = FactoredDenominator([], V)

            sage: f1.lcm(f2)
            {(1, 0, 0): 2, (0, 1, 0): 2}
            sage: f1.lcm(f3)
            {(1, 0, 0): 2, (0, 1, 0): 1, (0, 0, 1): 1}
            sage: f2.lcm(f3)
            {(1, 0, 0): 1, (0, 1, 0): 2, (0, 0, 1): 1}

            sage: f1.lcm(f4) == f1 and f4.lcm(f1) == f1
            True
        """
        res = self.copy()
        res.lcm_update(other)
        return res

    def gcd_update(self, other):
        if type(self) is not type(other):
            raise TypeError

        sd = self._dict
        od = other._dict
        for i, j in od.items():
            if j < sd.get(i, j):
                sd[i] = j

        self._tuple = tuple(sorted(self._dict.items()))

    def gcd(self, other):
        r"""
        Return the gcd of two factored denominator.

        EXAMPLES::

            sage: from surface_dynamics.misc.factored_denominator import FactoredDenominator

            sage: V = ZZ**3

            sage: f1 = FactoredDenominator([((1,0,0), 2),((0,1,0),1)], V)
            sage: f2 = FactoredDenominator([((1,0,0), 1),((0,1,0),2)], V)
            sage: f3 = FactoredDenominator([((1,0,0), 1),((0,1,0),1),((0,0,1),1)], V)
            sage: f4 = FactoredDenominator([], V)

            sage: f1.gcd(f2)
            {(1, 0, 0): 1, (0, 1, 0): 1}
            sage: f1.gcd(f2) == f1.gcd(f3) == f2.gcd(f3)
            True

            sage: f4.gcd(f1) == f4 and f1.gcd(f4) == f1
            True
        """
        res = self.copy()
        res.gcd_update(other)
        return res

    # TODO
    # this does not make much sense here
    def str_monomials(self, var_names='x'):
        terms = []
        for mon, mul in self._tuple:
            term = '(1 - %s)' % vector_to_monomial_string(mon, var_names)
            if not mul.is_one():
                term += '^%d' % mul
            terms.append(term)
        return '*'.join(terms)

    # TODO
    # this does not make much sense here
    def str_linear(self, var_names='h'):
        terms = []
        for mon, mul in self._tuple:
            term = '(' + vector_to_linear_form_string(mon, var_names) + ')'
            if not mul.is_one():
                term += '^%d' % mul
            terms.append(term)

        return '*'.join(terms)


# TODO: make it possible to use Laurent polynomials in denominator
# TODO: monomial substitution
# TODO: coefficient expansion Verdoolaege-Woods
# 1 - monomial versus linear form
class AbstractMSum(Element):
    r"""
    An abstract multiple sum
    """
    def __init__(self, parent, data, allow_multiple=False, check=True):
        r"""
        _data is a dictionary: {denominator: polynomial}
        """
        Element.__init__(self, parent)
        V = self.parent().free_module()
        R = self.parent().polynomial_ring()
        self._data = {}

        if not data:
            return

        if isinstance(data, (tuple, list)):
            data = iter(data)

        elif isinstance(data, dict):
            if not check:
                self._data = data
                return
            data = list(data.items())

        elif isinstance(data, AbstractMSum):
            self._data = data.copy()
            return

        for term in data:
            if not isinstance(term, (tuple, list)) or len(term) != 2:
                raise ValueError
            den, num = term

            num = R(num)
            if num.is_zero():
                continue

            if not isinstance(den, FactoredDenominator):
                den = FactoredDenominator(den, V)

            if den in self._data:
                if not allow_multiple:
                    raise ValueError('multiple times the same denominator in the input')
                else:
                    self._data[den] += num
            else:
                self._data[den] = num

    def _den_str(self, den):
        return str(den)

    def _repr_(self):
        r"""
        TESTS::

            sage: from surface_dynamics.misc.multiplicative_multivariate_generating_series import MultiplicativeMultivariateGeneratingSeriesRing

            sage: M = MultiplicativeMultivariateGeneratingSeriesRing('x', 2)
            sage: R = M.polynomial_ring()
            sage: x0, x1 = R.gens()
            sage: M.term(-x0 + x1^3, [((1,1), 1), ((1,2),2)])
            (x1^3 - x0)/((1 - x0*x1)*(1 - x0*x1^2)^2)

            sage: M = MultiplicativeMultivariateGeneratingSeriesRing('t,u')
            sage: t, u = M.polynomial_ring().gens()
            sage: M.term(-t + u^3, [((1,1), 1), ((1,2),2)])
            (u^3 - t)/((1 - t*u)*(1 - t*u^2)^2)
        """
        if not self._data:
            return '0'

        terms = []
        for den in sorted(self._data):
            num = self._data[den]
            fraction = '(' + str(num) + ')'
            if not den.is_one():
                fraction += '/'
                fraction += '(' + self._den_str(den) + ')'

            terms.append(fraction)

        return ' + '.join(terms)

    def is_trivial_zero(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.misc.multiplicative_multivariate_generating_series import MultiplicativeMultivariateGeneratingSeriesRing

            sage: M = MultiplicativeMultivariateGeneratingSeriesRing('x', 2)
            sage: M.term(0, []).is_trivial_zero()
            True
            sage: M.term(0, [((1,1),1)]).is_trivial_zero()
            True
            sage: M.term(1, []).is_trivial_zero()
            False
        """
        return not self._data

    def is_trivial_one(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.misc.multiplicative_multivariate_generating_series import MultiplicativeMultivariateGeneratingSeriesRing

            sage: M = MultiplicativeMultivariateGeneratingSeriesRing('x', 2)
            sage: M.term(0, []).is_trivial_one()
            False
            sage: M.term(0, [((1,1),1)]).is_trivial_one()
            False
            sage: M.term(1, []).is_trivial_one()
            True
            sage: M.term(1, [((1,1),1)]).is_trivial_one()
            False

            sage: f = M.term(1, [((1,0),1)]) + M.term(1, [((0,1),1)])
            sage: f.is_trivial_one()
            False
        """
        if len(self._data) != 1:
            return False
        (den, num) = next(iteritems(self._data))
        return num.is_one() and den.is_one()

    def _add_(self, other):
        r"""
        TESTS::

            sage: from surface_dynamics.misc.multiplicative_multivariate_generating_series import MultiplicativeMultivariateGeneratingSeriesRing

            sage: M = MultiplicativeMultivariateGeneratingSeriesRing('x', 3)
            sage: m1 = M.term(1, [((1,1,0),1)])
            sage: m2 = M.term(-2, [((1,0,0),1),((0,1,0),2)])
            sage: m1 + m2
            (-2)/((1 - x1)^2*(1 - x0)) + (1)/((1 - x0*x1))
        """
        if self.is_trivial_zero():
            return other
        if other.is_trivial_zero():
            return self

        ans_data = self._data.copy()
        for den, num in other._data.items():
            if den in ans_data:
                ans_data[den] += num
                if ans_data[den].is_zero():
                    del ans_data[den]
            else:
                ans_data[den] = num

        R = self.parent()
        return R.element_class(R, ans_data)

    def _mul_(self, other):
        r"""
        TESTS::

            sage: from surface_dynamics.misc.multiplicative_multivariate_generating_series import MultiplicativeMultivariateGeneratingSeriesRing

            sage: M = MultiplicativeMultivariateGeneratingSeriesRing('x', 3)
            sage: m1 = M.term(1, [((1,1,0),1)])
            sage: m2 = M.term(-2, [((1,0,0),1),((0,1,0),2)])
            sage: m1 * m2
            (-2)/((1 - x1)^2*(1 - x0)*(1 - x0*x1))

            sage: R = M.polynomial_ring()
            sage: M.one() * R.one()
            (1)

            sage: R = M.polynomial_ring()
            sage: x0, x1, x2 = R.gens()
            sage: m1 = M.term(x0 - 2*x2 + 1, [((1,1,1),2)])
            sage: m2 = M.term(x0*x1 + 1, [((0,1,2),1)])
            sage: m3 = M.term(1 + x0 + x1, [((1,0,1),1)])
            sage: m4 = M.term(x0 - 1, [((1,0,0),2),((1,0,1),1)])
            sage: (m1 + m2) * (m3 + m4) - m1 * m3 - m2 * m3 - m1 * m4 - m2 * m4
            0
        """
        ans_data = {}
        for den1, num1 in self._data.items():
            for den2, num2 in other._data.items():
                den = den1 * den2
                num = num1 * num2
                if den in ans_data:
                    ans_data[den] += num
                    if not ans_data[den]:
                        del ans_data[den]
                else:
                    ans_data[den] = num

        M = self.parent()
        return M.element_class(M, ans_data)

    def __neg__(self):
        r"""
        TESTS::

            sage: from surface_dynamics.misc.multiplicative_multivariate_generating_series import MultiplicativeMultivariateGeneratingSeriesRing

            sage: M = MultiplicativeMultivariateGeneratingSeriesRing('x', 3)
            sage: m1 = M.term(1, [((1,1,0),1)])
            sage: m1
            (1)/((1 - x0*x1))
            sage: -m1
            (-1)/((1 - x0*x1))
            sage: m1 - m1
            0
        """
        ans_data = {}
        for den, num in self._data.items():
            ans_data[den] = -num
        M = self.parent()
        return M.element_class(M, ans_data)

    def _richcmp_(self, other, op):
        r"""
        TESTS::

            sage: from surface_dynamics.misc.multiplicative_multivariate_generating_series import MultiplicativeMultivariateGeneratingSeriesRing

            sage: M = MultiplicativeMultivariateGeneratingSeriesRing('x', 2)
            sage: x, y = M.polynomial_ring().gens()
            sage: M.term(x*y+1, [((1,1),1)]) == M.term(x*y+1, [((1,1),1)])
            True
        """
        if op != op_EQ and op != op_NE:
            raise TypeError('no comparison available for multivariate series')

        if self._data == other._data:
            return op == op_EQ

        difference = self - other
        return bool(difference.factor()._data) == (op == op_NE)

    def degrees(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.misc.multiplicative_multivariate_generating_series import MultiplicativeMultivariateGeneratingSeriesRing
            sage: from surface_dynamics.misc.additive_multivariate_generating_series import AdditiveMultivariateGeneratingSeriesRing

            sage: M = MultiplicativeMultivariateGeneratingSeriesRing('x', 2)
            sage: f = M.term(1, [((1,0), 1), ((1,1),2)]) + M.term(2, [((1,0), 2), ((1,2),3)])
            sage: f.degrees()
            {3, 5}

            sage: A = AdditiveMultivariateGeneratingSeriesRing('x', 2)
            sage: f = A.term(1, [((1,0), 1), ((1,1),2)]) + A.term(2, [((1,0), 2), ((1,2),3)])
            sage: f.degrees()
            {3, 5}
        """
        res = set()
        for den, num in self._data.items():
            if num(*self.parent()._critical_point()).is_zero():
                raise ValueError('singular numerator')
            res.add(den.degree())
        return res


# TODO: the UniqueRepresentation should normalize the Laurent polynomial
class AbstractMSumRing(UniqueRepresentation, Parent):
    r"""
    TESTS::

        sage: from surface_dynamics.misc.multiplicative_multivariate_generating_series import MultiplicativeMultivariateGeneratingSeriesRing
        sage: MultiplicativeMultivariateGeneratingSeriesRing(5, 'x') is MultiplicativeMultivariateGeneratingSeriesRing(['x0','x1','x2','x3','x4'])
        True
    """
    def __init__(self, poly_ring):
        from sage.modules.free_module import FreeModule
        self._polynomial_ring = poly_ring
        dim = ZZ(poly_ring.ngens())
        self._free_module = FreeModule(ZZ, dim)

        # univariate extension of the polynomial ring
        # (needed in several algorithms)
        self._polynomial_ring_extra_var = self._polynomial_ring['EXTRA_VAR']

        Parent.__init__(self, category=Rings(), base=poly_ring.base_ring())

    def _critical_point(self):
        raise NotImplementedError

    def ngens(self):
        r"""
        Return the number of generators.

        EXAMPLES::

            sage: from surface_dynamics.misc.multiplicative_multivariate_generating_series import MultiplicativeMultivariateGeneratingSeriesRing
            sage: from surface_dynamics.misc.additive_multivariate_generating_series import AdditiveMultivariateGeneratingSeriesRing

            sage: MultiplicativeMultivariateGeneratingSeriesRing('x', 3).ngens()
            3
            sage: AdditiveMultivariateGeneratingSeriesRing('x', 3).ngens()
            3
        """
        return self.polynomial_ring().ngens()

    def free_module(self):
        return self._free_module

    def polynomial_ring(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.misc.multiplicative_multivariate_generating_series import MultiplicativeMultivariateGeneratingSeriesRing
            sage: from surface_dynamics.misc.additive_multivariate_generating_series import AdditiveMultivariateGeneratingSeriesRing


            sage: MultiplicativeMultivariateGeneratingSeriesRing('x', 3).polynomial_ring()
            Multivariate Laurent Polynomial Ring in x0, x1, x2 over Rational Field
            sage: AdditiveMultivariateGeneratingSeriesRing('x', 3).polynomial_ring()
            Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
        """
        return self._polynomial_ring

    def polynomial_ring_extra_var(self):
        return self._polynomial_ring_extra_var

    def with_extra_var(self):
        return self['EXTRA_VAR']

    @cached_method
    def zero(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.misc.multiplicative_multivariate_generating_series import MultiplicativeMultivariateGeneratingSeriesRing

            sage: M = MultiplicativeMultivariateGeneratingSeriesRing('x', 2)
            sage: M.zero()
            0
            sage: M.zero().parent() is M
            True
            sage: M.zero().is_zero()
            True
        """
        return self._element_constructor_(QQ.zero())

    @cached_method
    def one(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.misc.multiplicative_multivariate_generating_series import MultiplicativeMultivariateGeneratingSeriesRing

            sage: M = MultiplicativeMultivariateGeneratingSeriesRing('x', 2)
            sage: M.zero()
            0
            sage: M.zero().parent() is M
            True
            sage: M.one().is_one()
            True
        """
        return self._element_constructor_(QQ.one())

    def term(self, num, den):
        r"""
        Return the term ``num / den``.

        INPUT:

        - ``num`` - a Laurent polynomial

        - ``den`` - a list of pairs ``(vector, power)`` or a dictionary
           whose keys are the vectors and the values the powers. The
           vector ``v = (v_0, v_1, \ldots)`` with power ``n`` corresponds
           to the factor `(1 - x_0^{v_0} x_1^{v_1} \ldots x_k^{v_k})^n`.

        EXAMPLES::

            sage: from surface_dynamics.misc.multiplicative_multivariate_generating_series import MultiplicativeMultivariateGeneratingSeriesRing
            sage: from surface_dynamics.misc.additive_multivariate_generating_series import AdditiveMultivariateGeneratingSeriesRing


            sage: A = AdditiveMultivariateGeneratingSeriesRing('x', 3)
            sage: M = MultiplicativeMultivariateGeneratingSeriesRing('x', 3)
            sage: M.term(1, [([1,1,0],1),([1,0,-1],2)])
            (1)/((1 - x0*x2^-1)^2*(1 - x0*x1))
            sage: M.term(1, {(1,1,0): 1, (1,0,2): 2})
            (1)/((1 - x0*x2^2)^2*(1 - x0*x1))

            sage: A.term(1, [([1,1,0],1),([1,0,-1],2)])
            (1)/((x0 - x2)^2*(x0 + x1))
            sage: A.term(1, {(1,1,0): 1, (1,0,2): 2})
            (1)/((x0 + 2*x2)^2*(x0 + x1))

        Also works if the denominator is already a factored denominator::

            sage: from surface_dynamics.misc.factored_denominator import FactoredDenominator
            sage: M = MultiplicativeMultivariateGeneratingSeriesRing('x', 3)
            sage: den = FactoredDenominator({(1,0,0): 2, (1,1,1):1})
            sage: M.term(1, den)
            (1)/((1 - x0)^2*(1 - x0*x1*x2))

            sage: A.term(1, den)
             (1)/((x0)^2*(x0 + x1 + x2))

        TESTS::

            sage: from surface_dynamics.misc.multiplicative_multivariate_generating_series import MultiplicativeMultivariateGeneratingSeriesRing
            sage: M = MultiplicativeMultivariateGeneratingSeriesRing('x', 2)
            sage: R = M.polynomial_ring()
            sage: x0, x1 = R.gens()
            sage: d = {}
            sage: v0 = vector(QQ, (1,0)); v0.set_immutable()
            sage: v1 = vector(QQ, (1,1)); v1.set_immutable()
            sage: d[v1] = 1
            sage: d[v0] = 1
            sage: M.term(1*x0^2*x1, d)
            (x0^2*x1)/((1 - x0)*(1 - x0*x1))
        """
        return self.element_class(self, [(den, num)], self.free_module())

    def _element_constructor_(self, arg):
        r"""
        TESTS::

            sage: from surface_dynamics.misc.multiplicative_multivariate_generating_series import MultiplicativeMultivariateGeneratingSeriesRing
            sage: from surface_dynamics.misc.additive_multivariate_generating_series import AdditiveMultivariateGeneratingSeriesRing

            sage: M = MultiplicativeMultivariateGeneratingSeriesRing('x', 2)
            sage: A = AdditiveMultivariateGeneratingSeriesRing('x', 2)
            sage: M(1)
            (1)
            sage: A(1)
            (1)

            sage: R = M.polynomial_ring()
            sage: M(R.0)
            (x0)
            sage: A(R.0)
            (x0)
        """
        num = self._polynomial_ring(arg)
        return self.element_class(self, [([], num)], self.free_module())
