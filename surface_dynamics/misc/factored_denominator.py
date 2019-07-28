r"""
Factored denominators
"""
#*****************************************************************************
#       Copyright (C) 2019 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, print_function

from sage.rings.integer_ring import ZZ
from sage.modules.vector_integer_dense import Vector_integer_dense
from sage.structure.element import parent

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
    for i,j in enumerate(u):
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
        'x0 + -x2'
    """
    s = []
    for i,j in enumerate(u):
        if j:
            if isinstance(var_names, str):
                var = '%s%d' % (var_names, i)
            else:
                var = var_names[i]

            if j == 1:
                s.append(var)
            elif j == -1:
                s.append('-' + var)
            else:
                s.append('%d*%s' % (j, var))

    return ' + '.join(s) if s else '0'



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
                if mult <= 0:
                    raise ValueError('non-positive multiplicity in denominator')
                mon.set_immutable()
                self._dict[mon] = mult

        elif isinstance(data, dict):
            if V is not None:
                self._dict = {}
                for k,v in data.items():
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

    def degree(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.misc.factored_denominator import FactoredDenominator

            sage: V = ZZ**3
            sage: FactoredDenominator([((1,0,0), 2)], V).degree()
            2
            sage: FactoredDenominator([((0,1,2), 3), ((1,1,1), 1)], V).degree()
            4
            sage: FactoredDenominator([((0,-1,2), 1), ((1,0,0), 1), ((0,0,2), 1)], V).degree()
            3
        """
        return sum(x[1] for x in self._tuple)

    def to_polynomial(self, S, extra_var=False):
        r"""
        Return the product of the term in a given polynomial ring ``S``.

        EXAMPLES::

            sage: from surface_dynamics.misc.factored_denominator import FactoredDenominator

            sage: V = ZZ ** 3
            sage: f = FactoredDenominator([((1,0,0), 2)], V)
            sage: g = FactoredDenominator([((1,0,0), 2), ((1,1,1),1)], V)

            sage: R1 = QQ['x,y,z']
            sage: f.to_polynomial(R1)
            x^2 - 2*x + 1
            sage: g.to_polynomial(R1)
            -x^3*y*z + 2*x^2*y*z - x*y*z + x^2 - 2*x + 1

            sage: f.to_polynomial(R1['E'], extra_var=True)
            x^2*E^2 - 2*x*E + 1
            sage: g.to_polynomial(R1['E'], extra_var=True)
            -x^3*y*z*E^5 + 2*x^2*y*z*E^4 - x*y*z*E^3 + x^2*E^2 - 2*x*E + 1

            sage: g.to_polynomial(QQ['t1,t2,t3'])
            -t1^3*t2*t3 + 2*t1^2*t2*t3 - t1*t2*t3 + t1^2 - 2*t1 + 1
        """
        if extra_var:
            R = S.base_ring()
            T = S.gen()
        else:
            R = S
            T = R.one()

        ans = S.one()
        for a,i in self._tuple:
            ans *= (S.one() - R.monomial(*a) * T ** sum(a)) ** i
        return ans

    def inverse_series_trunc(self, S, prec):
        r"""
        Return a truncation of the inverse series of this element.

        INPUT:

        - ``S`` - polynomial ring in one variable whose base ring is a
          multivariate polynomial ring with the same number of variables
          as this factor.
        
        - ``prec`` - precision of the Taylor expansion (in terms of total degree)

        EXAMPLES::

            sage: from surface_dynamics.misc.factored_denominator import FactoredDenominator

            sage: V = ZZ**2
            sage: R = QQ['t,u']
            sage: f = FactoredDenominator([((2,3), 2), ((4,1), 1), ((1,1), 1)], V)
            sage: f.inverse_series_trunc(R['X'], 10)
            (t^6*u^3 + 2*t^4*u^5)*X^9 +
            t^4*u^4*X^8 +
            (t^5*u^2 + 2*t^3*u^4)*X^7 +
            t^3*u^3*X^6 +
            (t^4*u + 2*t^2*u^3)*X^5 +
            t^2*u^2*X^4 +
            t*u*X^2 +
            1
        """
        # this is a very stupid way of doing so
        # we should Taylor expand each monomial and then use _mul_trunc_
        return self.to_polynomial(S, extra_var=True).inverse_series_trunc(prec)

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
        return sum(m for _,m in self._tuple)

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

        OUTPUT: itertor of triples ``(mult, v, monomial)``

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
        for i,j in other._dict.items():
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
        for i,j in od.items():
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
        for i,j in od.items():
            if j > sd.get(i, -1):
                sd[i] = j

        self._tuple = tuple(sorted(self._dict.items()))

    def lcm(self, other):
        r"""
        Return the lcm of two factored denomiator.

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
        for i,j in od.items():
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

