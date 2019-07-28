r"""
Rational multivariate quasi-polynomial functions

Given a rational function `F` in the variables `x1, \ldots, xm` with no
pole at zero can be expanded as a (Laurent) power series `\sum f(i0, i1, ...,
im) x0^{i0} ... xm^{im}`. This modules deal with the case when `f` are quasi-polynomial
are equivalently when `F` is a finite sum of rational functions of the form

.. MATH::

    \frac{P}{\prod (1 - m_i)^{d_i}}

where the `m_i` are monomials and `d_i` are positive integers.

REFERENCES:

DH:03

[LHTY] De Loera, R. Hemmecke, J. Tauzer, R. Yoshida
 "Effective lattice point counting in rational convex polytopes"
 (LattE)

[Ba] Barvinok

[St] Stanley
"Enumerative Combinatorics"
Volume I
Cambridge Studies in Advanced Mathematics 49 (1997)

.. TODO::

 - make it barvinok compatible, the data-structure is
   a term is a ``struct short_rat`` while the generating
   series is a ``struct gen_fun``
"""
#*****************************************************************************
#       Copyright (C) 2019 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function, absolute_import
from six.moves import range, map, filter, zip
from six import iteritems

from .factored_denominator import FactoredDenominator

import numbers

from sage.misc.cachefunc import cached_method

from sage.structure.sage_object import SageObject
from sage.structure.parent import Parent
from sage.structure.element import Element, parent

try:
    from sage.structure.richcmp import op_LT, op_EQ, op_NE, op_LE, op_GE, op_GT, op_LT
except ImportError:
    # sage version <= 7.6
    from sage.structure.sage_object import op_LT, op_EQ, op_NE, op_LE, op_GE, op_GT, op_LT


from sage.categories.all import CommutativeAlgebras, Rings

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.modules.vector_integer_dense import Vector_integer_dense



def laurent_monomial(R, arg):
    r"""
    EXAMPLES::

        sage: from surface_dynamics.misc.multivariate_generating_series import laurent_monomial
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
    return prod(x**int(i) for (x,i) in zip(R.gens(), arg))

# custom latte count
def latte_generating_series(L, M=None):
    r"""
    EXAMPLES::

        sage: from surface_dynamics.misc.multivariate_generating_series import latte_generating_series
        sage: ieqs = [[0, 1, 0, 0, 0, 0, 0],
        ....:         [0, 0, 1, -1, 1, 0, 0],
        ....:         [0, 0, 0, 0, 1, 0, 0],
        ....:         [0, 0, 0, 1, 0, 0, 0],
        ....:         [0, 0, 1, 0, 0, 0, 0]]
        sage: eqns = [[0, 0, 1, -1, 1, 0, -1], [0, 0, 1, -1, 1, -1, 0]]
        sage: L = Polyhedron(ieqs=ieqs, eqns=eqns)
        sage: latte_generating_series(L)   # optional - latte_int
        (1)/((1 - x2^-1*x4*x5)*(1 - x2*x3)*(1 - x1*x2)*(1 - x0)) + (1)/((1 - x3*x4*x5)*(1 - x2*x4^-1*x5^-1)*(1 - x1*x4*x5)*(1 - x0))
    """
    if M is None:
        M = MultivariateGeneratingSeriesRing('x', L.ambient_dim())
    try:
        from sage.interfaces.latte import count
    except ImportError:
        from sage.version import version
        raise ValueError('your Sage version is too old ({}) to use this function'.format(version))
    ans = count(L.cdd_Hrepresentation(), cdd=True, multivariate_generating_function=True, raw_output=True)
    return parse_latte_generating_series(M, ans)

def parse_latte_generating_series(M, s):
    r"""
    INPUT:

    - ``M`` - a multivariate generating series ring

    - ``s`` - a string as given by LattE

    OUTPUT: multivariate short rational function
    """
    n = M.ngens()
    R = M.laurent_polynomial_ring()
    V = M.free_module()

    terms = s.strip().split('\n + ')
    m = M.zero()
    for term in terms:
        num, den = term.split('/')
        m_den = {}
        for term in den[2:-2].split(')*('):
            mon = V(0)
            assert term.startswith('1-')
            term = term[2:].split('*')
            for v in term:
                if '^' in v:
                    v, mult = v.split('^')
                    if mult.startswith('(') and mult.endswith(')'):
                        mult = mult[1:-1]
                    mult = int(mult)
                else:
                    mult = 1

                assert v.startswith('x[') and v.endswith(']'), (s, term, v)
                i = int(v[2:-1])
                mon[i] += mult
            mon.set_immutable()
            if mon in m_den:
                m_den[mon] += 1
            else:
                m_den[mon] = 1
        m += M.term(R(num), list(m_den.items()))

    return m


def monomial_projection(num, den, i, M):
    r"""
    Return the monomial substitution corresponding to setting the ``i``-th variable
    to ``1`` in the term ``num / den`` of the multivariate series ring ``M``.

    INPUT:

    - ``num`` - multivariate polynomial

    - ``den`` - a factored denominator

    - ``i`` - (integer) the index of a variable

    - ``M`` - multivariate generating series ring


    EXAMPLES::

        sage:
    """
    R = parent(num)
    n = R.ngens()
    x = R.gen(i)

    # formally we replace by (1+t) and let t tend to 0
    B = []   # singular terms ("bad")
    G = []   # regular terms  ("good")
    for v, mult in den:
        if v[i] != 0 and sum(v[j] != 0 for j in range(n)) == 1:
            B.append((v[i], mult))
        else:
            G.append((v,mult))


    # 1. regular case
    if not B:
        new_dict = {}
        for v,mult in G:
            v = v.__copy__()
            v[i] = 0
            v.set_immutable()
            if v in new_dict:
                new_dict[v] += mult
            else:
                new_dict[v] = mult

        num = num.subs({x: 1})
        return M.term(num, den)

    # 2. singular case
    r = sum(x[1] for x in B) + 1   # r-1 is the degree of the pole
    S = M.laurent_polynomial_ring_extra_var()
    SS = M.with_extra_var()
    t = S.gen()
    tt = SS.gen()

    P = num.subs({x: 1+t}).truncate(r)

    # 2.a) each singular term is 1 / (1 - (1+t)^j)^mult
    for j,mult in B:
        q = ((1 - (1 + t)**j) // t).inverse_series_trunc(r)._power_trunc(mult, r)
        P = q._mul_trunc_(P, r)

    # 2.b) each regular term is 1 / (1 - (1+t)^j s)^mult
    for v,mult in G:
        if v[i] == 0:  # no t involved
            P *= M.term(1, [(v,mult)])
        else:          # develop in powers of t
            j = v[i]
            v = v.__copy__()
            v[i] = 0
            v.set_immutable()
            s = laurent_monomial(S, v)
            z = M.term(1, [(1-s), 1])
            q = (z * (1 - ((1+t)**j - 1) * s * z).inverse_series_trunc(r))._power_trunc(mult, r)
            P = q._mul_trunc_(P, r)

    return P[r-1]



#def monomial_substitution(mat, num, den, M):
#    r"""
#    Apply the monomial substitution on ``num`` / ``den`` given by ``mat``.
#
#    We assume that the matrix applies on rows (in particular it should
#    have the same number of rows as the dimension of the monomial in
#    this factored denominator).
#
#    If some monomial is in the kernel, a vector ``c`` needs to be provided
#    and in that case this function return the constant term of the expansion
#    in the direction given by ``c`` (the result might depend on the choice of
#    the vector ``c`` provided).
#
#    INPUT:
#
#    - ``num`` - a polynomial
#
#    - ``den`` - a factored denominator
#
#    - ``M`` - a multivariate series ring
#
#    - ``c`` - an optional vector
#
#    EXAMPLES::
#
#        sage: from surface_dynamics.misc.multivariate_generating_series import FactoredDenominator
#        sage: V = ZZ**3
#        sage: f = FactoredDenominator([((1,0,0), 1),((0,1,0),1),((0,0,1),1)], V)
#        sage: f.monomial_substitution(matrix(3, 2, [1,1,1,0,1,0]))
#        (1 - x0)^2*(1 - x0*x1)
#        sage: f.monomial_substitution(matrix(3, 4, [1,1,1,1,1,1,0,0,0,0,1,1]))
#        (1 - x2*x3)*(1 - x0*x1)*(1 - x0*x1*x2*x3)
#    """
#    B = []
#    G = []
#    for a,i in self._tuple:
#        b = a * mat
#        if b.is_zero():
#            B.append((a,i))
#        else:
#            b.set_immutable()
#            G.append((a,b,i))
#
#    # no singular part
#    if not B:
#        new_dict = {}
#        for a,b,i in G:
#            if b in new_dict:
#                new_dict[b] += i
#            else:
#                new_dict[b] = i
#        return M.term(1, new_dict)
#
#    # case with singular part
#    l = sum(i for a,i in B)  # degree of the singular part
#    # XXX
#    raise NotImplementedError('singular case not implemented')




# This should be a real number!
# ... and we need an algebra to handle all these values properly
class _OLD_TOTO_(object):
    r"""
    A generalized multiple zeta value.

    A *generalized multiple zeta value* is a sum of the form

    .. MATH::

        \sum_{h in ZZ^s_+} chi(h) /  \ell_1(h)^{-1} \ell_2(h)^{-1} ... \ell_d(h)^{-1}

    where `chi(h)` is a character, `\ell_1(h)`, `\ell_2(h)`, ..., `\ell_d(h)`
    are linear form with non-negative coefficients. The integer `d` is called
    the *weight*.

    The version with character is not implemented!
    """
    def gp_multizeta_sum_cmd(self, method='sumnummonien', tabname='T'):
        r"""
        Return a numerical approximation of the corresponding generalized multizeta
        """
        var_names = ['n%d'%n for n in range(1, self.length()+1)]

        # TODO: reorder the term to move out as much as possible the linear forms
        t = self._den._tuple

        if not t:
            return '1'

        cmd = "%s(%%s=1, %%s, %s)" %(method, tabname)
        d = len(t[0][0])
        for i in range(d-1):
            cmd = cmd % (var_names[i], '%s(%%s=1, %%s, %s)'% (method, tabname))

        tot = []
        for mon, mul in self._den:
            l = den.str_linear()
            if mul:
                l = '%s^%d' %(l, mul)
            tot.append(l)

        f = '1 / (%s)' % '*'.join(tot)
        cmd = cmd % (var_names[d-1], f)
        return cmd


# TODO: make it possible to use Laurent polynomials in denominator
# TODO: monomial substitution
# TODO: coefficient expansion Verdoolaege-Woods
class AbstractMSum(Element):
    r"""
    An abstract multiple sums for generating functions
    """
    def __init__(self, parent, data, allow_multiple=False, check=True):
        r"""
        _data is a dictionary: {denominator: polynomial}
        """
        Element.__init__(self, parent)
        V = self.parent().free_module()
        R = self.parent().laurent_polynomial_ring()
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
            if not isinstance(term, (tuple,list)) or len(term) != 2:
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

    def is_trivial_zero(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.misc.multivariate_generating_series import MultivariateGeneratingSeriesRing

            sage: M = MultivariateGeneratingSeriesRing('x', 2)
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

            sage: from surface_dynamics.misc.multivariate_generating_series import MultivariateGeneratingSeriesRing

            sage: M = MultivariateGeneratingSeriesRing('x', 2)
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

            sage: from surface_dynamics.misc.multivariate_generating_series import MultivariateGeneratingSeriesRing

            sage: M = MultivariateGeneratingSeriesRing('x', 3)
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
        for den,num in other._data.items():
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

            sage: from surface_dynamics.misc.multivariate_generating_series import MultivariateGeneratingSeriesRing

            sage: M = MultivariateGeneratingSeriesRing('x', 3)
            sage: m1 = M.term(1, [((1,1,0),1)])
            sage: m2 = M.term(-2, [((1,0,0),1),((0,1,0),2)])
            sage: m1 * m2
            (-2)/((1 - x1)^2*(1 - x0)*(1 - x0*x1))

            sage: R = M.laurent_polynomial_ring()
            sage: M.one() * R.one()
            (1)

            sage: R = M.laurent_polynomial_ring()
            sage: x0, x1, x2 = R.gens()
            sage: m1 = M.term(x0 - 2*x2 + 1, [((1,1,1),2)])
            sage: m2 = M.term(x0*x1 + 1, [((0,1,2),1)])
            sage: m3 = M.term(1 + x0 + x1, [((1,0,1),1)])
            sage: m4 = M.term(x0 - 1, [((1,0,0),2),((1,0,1),1)])
            sage: (m1 + m2) * (m3 + m4) - m1 * m3 - m2 * m3 - m1 * m4 - m2 * m4
            0
        """
        ans_data = {}
        for den1,num1 in self._data.items():
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

            sage: from surface_dynamics.misc.multivariate_generating_series import MultivariateGeneratingSeriesRing

            sage: M = MultivariateGeneratingSeriesRing('x', 3)
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

            sage: from surface_dynamics.misc.multivariate_generating_series import MultivariateGeneratingSeriesRing

            sage: M = MultivariateGeneratingSeriesRing('x', 2)
            sage: x, y = M.laurent_polynomial_ring().gens()
            sage: M.term(x*y+1, [((1,1),1)]) == M.term(x*y+1, [((1,1),1)])
            True
        """
        if op != op_EQ and op != op_NE:
            raise TypeError('no comparison available for multivariate series')

        if self._data == other._data:
            return op == op_EQ

        (sden,snum), = self.factor()._data.items()
        (oden,onum), = other.factor()._data.items()

        if sden == oden and snum == onum:
            return op == op_EQ

        raise NotImplementedError

# TODO: here the elements are not necessarily with the same number of variables
# we want the following to work
# 1/3 * Sum(1/((n0)^4), n0=1..+oo) + 2/3 * Sum(1/((n0)*(n0 + n1)^3), n0=1..+oo, n1=1..+oo) + 1/3 * Sum(1/((n0)^2*(n0 + n1)^2), n0=1..+oo, n1=1..+oo)
class GeneralizedMultiZetaElement(AbstractMSum):
    def simplify(self):
        r"""
        Canonicalize each term by applying a permutation
        """
        from itertools import permutations

        P = self.parent()
        new_data = []
        n = P.ngens()
        for den,num in iteritems(self._data):
            # TODO: this is stupidly innefficient
            L = [list(k) + [v] for k,v in den._tuple]
            L.sort()
            Lmin = L
            for p in permutations(range(n)):
                LL = [[l[p[i]] for i in range(n)] + [l[n]] for l in L]
                LL.sort()
                if LL < Lmin:
                    Lmin = LL

            new_den = [(k[:-1], k[-1]) for k in Lmin]
            new_data.append((new_den, num))
        return P.element_class(P, new_data, allow_multiple=True)

    def expand(self):
        r"""
        Apply 1 / (L1 L2) = 1/(L1+L2) * (1/L1 + 1/L2) to try getting
        multiple zeta values
        """
        raise NotImplementedError

    def weights(self):
        r"""
        Return the set of weights of the generalized multiple zeta values.
        """
        return sorted(set(den.degree() for den in self._data))

    def _repr_(self):
        r"""
        String representation
        """
        terms = []
        var_names = ['n%d'%n for n in range(self.parent().laurent_polynomial_ring().ngens())]
        for den in sorted(self._data.keys()):
            num = self._data[den]
            fraction = str(num)
            if not den.is_one():
                fraction += ' * Sum(1/'
                linterm = den.str_linear(var_names)
                fraction += '(' + den.str_linear(var_names) + '), '
                fraction += ', '.join("%s=1..+oo"%v for v in var_names) + ')'
            terms.append(fraction)

        return ' + '.join(terms)

class MultivariateGeneratingSeries(AbstractMSum):
    def __nonzero__(self):
        return bool(self._data)

    def summands(self):
        r"""
        Return the list of elementary summands

        EXAMPLES::

           sage: from surface_dynamics.misc.multivariate_generating_series import *
           sage: M = MultivariateGeneratingSeriesRing(2, 'x')
           sage: m1 = M.term(1, [((1,-1),1), ((0,1),1)])
           sage: m2 = M.term(1, [((-1,1),1), ((1,0),1)])
           sage: f = (m1 + m2)**2
           sage: for s in f.summands(): print(s)
           (2)/((1 - x0^-1*x1)*(1 - x1)*(1 - x0*x1^-1)*(1 - x0))
           (1)/((1 - x0^-1*x1)^2*(1 - x0)^2)
           (1)/((1 - x1)^2*(1 - x0*x1^-1)^2)
        """
        M = self.parent()
        res = []
        for den in sorted(self._data):
            num = self._data[den]
            res.append(M.term(num,den))
        return res

    # TODO: this is only working term by term. It completely ignores pole cancellation between
    # the terms
    #
    # 1 / (1 - x0 x1^-1) (1 - x1) + 1 / (1 - x0^-1 x1) (1 - x0)
    #
    #     sage: M = MultivariateGeneratingSeriesRing(2)
    #     sage: m1 = M.term(1, [((1,-1),1), ((0,1),1)])
    #     sage: m2 = M.term(1, [((-1,1),1), ((1,0),1)])
    #     sage: f = m1 + m2
    #
    # This is solved by applying the monomial substitution as in
    # Barvinok Woods 2003 "Short rational generating functions for lattice point problems"
    # (possible alternative in Verdoolaege)
    #
#    def monomial_substitution(self, arg, M=None):
#        r"""
#        Perform a monomial substitution
#
#        - ``arg`` - a list of monomials or vectors or a matrix
#
#        - ``M`` - an optional codomain
#
#        EXAMPLES::
#
#            sage: M = MultivariateGeneratingSeriesRing(2)
#            sage: R = M.laurent_polynomial_ring()
#            sage: x0,x1 = R.gens()
#
#            sage: f = M.term(x0, [((1,0),2),((0,1),1)])
#            sage: f
#            (x0)/((1 - x1)*(1 - x0)^2)
#            sage: f.monomial_substitution(matrix(2, [1,1,1,1]))
#            (x0*x1)/((1 - x0*x1)^3)
#
#            sage: C = MultivariateGeneratingSeriesRing(3, 'z')
#            sage: f.monomial_substitution(matrix(2, 3, [1,0,1,0,1,1]), C)
#            (z0*z2)/((1 - z1*z2)*(1 - z0*z2)^2)
#        """
#        M_domain = self.parent()
#        R_domain = M_domain.laurent_polynomial_ring()
#
#        if M is None:
#            M_codomain = M_domain
#            R_codomain = R_domain
#        else:
#            M_codomain = M
#            R_codomain = M.laurent_polynomial_ring()
#
#        if arg.nrows() != M_domain.ngens() or arg.ncols() != M_codomain.ngens():
#            raise ValueError('wrong number of rows or columns')
#
#        s_dict = {R_domain.gen(i): R_codomain.monomial(*row) for i,row in enumerate(arg)}
#
#        # the monomial substitution also depends on the numerator !!
#        new_data = {}
#        for den, num in self._data.items():
#            monomial_substitution(arg, num, den, M_codomain)
#            assert num.parent() is R_domain
#            num = R_codomain(num.subs(s_dict))
#            if den in new_data:
#                new_data[den] += num
#            else:
#                new_data[den] = num
#
#        return M_codomain.element_class(M_codomain, new_data)

    # TODO: this does not make any sense
    # the residue is 1 / L1(h) ... Ld(h) where the Li are linear forms. At
    # no point we are asking for a sum of them!!
    # In other words, this method should just return a function h -> 1 / prod(Li)
    def residue(self):
        r"""
        denominator: each (1 - mon)^k in denom is replaced with -> mon^k
        numerator: evaluate at (1,1,...,1)

        OUTPUT: a pair '(degree, value)`

        EXAMPLES::

            sage: from surface_dynamics.misc.multivariate_generating_series import MultivariateGeneratingSeriesRing

            sage: M = MultivariateGeneratingSeriesRing(2, 'x')
            sage: R = M.laurent_polynomial_ring()
            sage: x0,x1 = R.gens()
            sage: f = M.term(x0, [((1,1),2)])
            sage: f.residue()
            (2, [(1, {(1, 1): 2})])
            sage: f = M.term(x0, [((1,1),2)]) + M.term(1, [((1,0),1),((0,1),1),((1,1),1)])
            sage: r = f.residue()
            sage: r # random
            (3, [(1, {(1, 0): 1, (0, 1): 1, (1, 1): 1})])
            sage: f = M.term(x0, [((1,1),2)]) + M.term(1, [((1,0),1),((1,1),1)])
            sage: r = f.residue()
            sage: r # random
            (2, [(1, {(1, 1): 2}), (1, {(1, 0): 1, (1, 1): 1})])
        """
        R = self.parent().laurent_polynomial_ring()
        one = QQ.one()
        values = {g: one for g in R.gens()}
        ans = []
        d = -1
        for den, num in self._data.items():
            if den.degree() >= d:
                if den.degree() > d:
                    ans = []
                    d = den.degree()
                num = QQ(num.subs(values))
                if num.is_zero():
                    raise NotImplementedError('zero numerator')
                ans.append((num, den))
        return d, ans

    def derivative(self, var):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.misc.multivariate_generating_series import MultivariateGeneratingSeriesRing

            sage: M = MultivariateGeneratingSeriesRing('x', 2)
            sage: R = M.laurent_polynomial_ring()
            sage: x0, x1 = R.gens()

            sage: f = M.term(1, [((1,0),1)])
            sage: f.derivative(0)
            (1)/((1 - x0)^2)
            sage: f.derivative(0).derivative(0)
            (2)/((1 - x0)^3)

            sage: xx0, xx1 = R.polynomial_ring().gens()
            sage: f.taylor(10).derivative(xx0) == f.derivative(0).taylor(9)
            True
            sage: f.derivative(1)
            0

            sage: f = M.term(x0 * x1, [((1,1),2)])
            sage: f.derivative(0)
            (x1)/((1 - x0*x1)^2) + (2*x0*x1^2)/((1 - x0*x1)^3)
            sage: f.derivative(0).taylor(10) - f.taylor(10).derivative(xx0)
            30*x0^5*x1^6

            sage: f = M.term(1, [((1,2),2), ((0,1),2), ((1,0),2), ((1,1),2)])
            sage: min((f.taylor(10).derivative(xx0) - f.derivative(0).taylor(10)).degrees())
            9
            sage: min((f.taylor(10).derivative(xx1) - f.derivative(1).taylor(10)).degrees())
            10

            sage: f = M.term(1, [((1,0),2), ((1,1),1)])
            sage: f.taylor(10).derivative(xx0).derivative(xx1) - f.derivative(0).derivative(1).taylor(9)
            -18*x0^9*x1 - 42*x0^8*x1^2 - ... - 72*x0^5*x1^3 - 25*x0^4*x1^4

            sage: f = M.zero() + M.term(1/1*x0^2 * x1, [((1,0),1), ((1,1),1)])
            sage: f.derivative(0)
            (2*x0*x1)/((1 - x0)*(1 - x0*x1)) + (x0^2*x1^2)/((1 - x0)*(1 - x0*x1)^2) + (x0^2*x1)/((1 - x0)^2*(1 - x0*x1))

        You can indistinctly use integers, strings or polynomial variables for ``var``::

            sage: f = M.term(1, [((1,1),1)])
            sage: f.derivative('x0')
            (x1)/((1 - x0*x1)^2)
            sage: f.derivative(0)
            (x1)/((1 - x0*x1)^2)
            sage: f.derivative(R.gen(0))
            (x1)/((1 - x0*x1)^2)

        Checking errors in the input::

            sage: f.derivative(-1)
            Traceback (most recent call last):
            ...
            ValueError: generator not defined
            sage: f.derivative('q')
            Traceback (most recent call last):
            ...
            ValueError: 'q' not valid as a variable
        """
        M = self.parent()
        R = M.laurent_polynomial_ring()

        if isinstance(var, numbers.Integral):
            j = var
        else:
            try:
                j = R.gens().index(R(var))
            except (ValueError, TypeError):
                raise ValueError('%r not valid as a variable' % var)

        var = R.gen(j)

        ans = M.zero()
        for den, num in self._data.items():
            ans += M.term(num.derivative(var), den)
            s = sum(mul * laurent_monomial(R, v) * M.term(1, mon) for mul, v, mon in den.logarithmic_minus_derivative(j))
            ans += M.term(num, den) * s
        return ans

    def derivative_up_to_lower_order_terms(self, var):
        r"""
        Each term ``u/v`` is replaced by ``-u v'/v^2``.

        This corresponds to one half of the derivative, the other half
        being ``u' / v``. This second part can be ignored when asymptotics
        question are considered.

        EXAMPLES::

            sage: from surface_dynamics.misc.multivariate_generating_series import MultivariateGeneratingSeriesRing

            sage: M = MultivariateGeneratingSeriesRing('x', 2)
            sage: f = M.term(1, [((1,0),2)])
            sage: f
            (1)/((1 - x0)^2)
            sage: f.derivative_up_to_lower_order_terms(0)
            (2)/((1 - x0)^3)

            sage: f = M.term(1, [((1,0),1),((0,1),1),((1,1),1)])
        """
        M = self.parent()
        R = M.laurent_polynomial_ring()

        ans = M.zero()
        for den, num in self._data.items():
            s = M.zero()
            for m,v,d in den.logarithmic_minus_derivative(var):
                s += M.term(m * laurent_monomial(R,v), [(d,1)])
            ans += M.term(num, den) * s

        return ans

    def delta(self):
        r"""
        Take a derivative (up to lower order terms) with respect to each of the variables.
        """
        a = self
        for i in range(self.parent().ngens()):
            a = a.derivative_up_to_lower_order_terms(i)
        return a

    # what does prec means for Laurent polynomials!?
    def taylor(self, prec, R=None):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.misc.multivariate_generating_series import MultivariateGeneratingSeriesRing

            sage: M = MultivariateGeneratingSeriesRing('t,u')
            sage: t,u = M.laurent_polynomial_ring().gens()
            sage: f = M.term(t + u^3, [((1,1), 1), ((1,2),2)])
            sage: f
            (u^3 + t)/((1 - t*u)*(1 - t*u^2)^2)
            sage: f.taylor(10)
            2*t^4*u^8 + 4*t^3*u^9 + t^4*u^7 + 3*t^3*u^8 + ... + t*u^4 + 2*t^2*u^2 + t^2*u + u^3 + t

        TODO: this is only working term by term but for example containing negative
        powers as::

            sage: M = MultivariateGeneratingSeriesRing('x', 2)   # not tested
            sage: m1 = M.term(1, [((1,-1),1), ((0,1),1)])   # not tested
            sage: m2 = M.term(1, [((-1,1),1), ((1,0),1)])   # not tested
            sage: f = m1 + m2                               # not tested
            sage: f.taylor(10)                              # not tested

        TESTS::

            sage: M = MultivariateGeneratingSeriesRing('x', 2)
            sage: R = M.laurent_polynomial_ring()
            sage: x0, x1 = R.gens()
            sage: type(M.term(1, [((1,0),1)]).taylor(10))
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
        """
        M = self.parent()
        R = M.laurent_polynomial_ring().polynomial_ring()
        S = M.laurent_polynomial_ring().polynomial_ring()['EXTRA_VAR']
        ans = S.zero()
        for den, num in self._data.items():
            ans += (S(num) * den.inverse_series_trunc(S, prec)).truncate(prec)
        return R(ans.subs({S.gen(): R.one()}))

    def is_zero(self):
        if self.is_trivial_zero():
            return True
        raise NotImplementedError

    def is_one(self):
        if self.is_trivial_one():
            return True
        raise NotImplementedError

    def degrees(self):
        return (den.degree() for den in self._data.keys())

    def _repr_(self):
        r"""
        TESTS::

            sage: from surface_dynamics.misc.multivariate_generating_series import MultivariateGeneratingSeriesRing

            sage: M = MultivariateGeneratingSeriesRing('x', 2)
            sage: R = M.laurent_polynomial_ring()
            sage: x0, x1 = R.gens()
            sage: M.term(-x0 + x1^3, [((1,1), 1), ((1,2),2)])
            (x1^3 - x0)/((1 - x0*x1)*(1 - x0*x1^2)^2)

            sage: M = MultivariateGeneratingSeriesRing('t,u')
            sage: t, u = M.laurent_polynomial_ring().gens()
            sage: M.term(-t + u^3, [((1,1), 1), ((1,2),2)])
            (u^3 - t)/((1 - t*u)*(1 - t*u^2)^2)
        """
        if not self._data:
            return '0'

        var_names = self.parent().laurent_polynomial_ring().variable_names()
        terms = []
        for den in sorted(self._data):
            num = self._data[den]
            fraction = '(' + str(num) + ')'
            if not den.is_one():
                fraction += '/'
                fraction += '(' + den.str_monomials(var_names) + ')'

            terms.append(fraction)

        return ' + '.join(terms)

    def subs_numerator(self, *args, **kwds):
        r"""
        TESTS::

            sage: from surface_dynamics.misc.multivariate_generating_series import MultivariateGeneratingSeriesRing
            sage: M = MultivariateGeneratingSeriesRing('x', 2)
            sage: R = M.laurent_polynomial_ring()
            sage: x0, x1 = R.gens()
            sage: f = M.term(-x0 + 2*x0*x1^3, [((1,1), 1), ((1,2),2)])
            sage: f.subs_numerator(x0=1, x1=1)
            (1)/((1 - x0*x1)*(1 - x0*x1^2)^2)
        """
        M = self.parent()
        res = M.zero()
        for den, num in self._data.items():
            res += M.term(num.subs(*args, **kwds), den)
        return res

    def factor(self, simplify=False):
        """
        Group all the partial fractions into a unique fraction.

        EXAMPLES::

            sage: from surface_dynamics.misc.multivariate_generating_series import MultivariateGeneratingSeriesRing
            sage: M = MultivariateGeneratingSeriesRing('x', 3)
            sage: f = M.term(1, [([1,3,0],1),([1,0,-1],1)])
            sage: g = M.term(1, [([1,1,0],1),([1,0,-1],2)])
            sage: (f + g).factor()
            (-x0*x1^3 + x0^2*x1*x2^-1 - x0*x1 - x0*x2^-1 + 2)/((1 - x0*x2^-1)^2*(1 - x0*x1)*(1 - x0*x1^3))
        """
        M = self.parent()

        # compute the product of denominators
        V = M.free_module()
        D = FactoredDenominator([], V)
        for den, num in self._data.items():
            D.lcm_update(den)

        R = M.laurent_polynomial_ring()
        N = R.zero()
        for den,num in self._data.items():
            N += (D / den).to_polynomial(R)

        # TODO: look more carefully for simplification
        # note: (1 - x0x1) and (1 - x0^-1 x1^-1) are the same things up
        # to sign and division by x0x1. We can canonicalize as soon as the numerator
        # is gentle enough... we should really introduce the polynomial version!

        gt = False
        for mon,mult in D._dict.items():
            pmon = laurent_monomial(R, mon)
            q,r = N.quo_rem(pmon)
            while mult and not r:
                gt = True
                mult -= 1
                N = q
                q,r = N.quo_rem(pmon)

            D._dict[mon] = mult

        if gt:
            D._tuple = tuple(sorted(self._dict.items()))

        return M.term(N, D)

    def as_symbolic(self):
        from sage.symbolic.ring import SR
        return SR(str(self))


class AbstractMSumRing(Parent):
    def __init__(self, *args):
        if len(args) == 1 and isinstance(args[0], AbstractMSumRing):
            self._laurent_polynomial_ring = args[0]._laurent_polynomial_ring
            self._free_module = args[0]._free_module
            self._laurent_polynomial_ring_extra_var = args[0]._laurent_polynomial_ring_extra_var
        else:
            from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
            from sage.modules.free_module import FreeModule
            self._laurent_polynomial_ring = LaurentPolynomialRing(QQ, *args)
            dim = ZZ(self._laurent_polynomial_ring.ngens())
            self._free_module = FreeModule(ZZ, dim)

            # univariate extension of the polynomial ring
            # (needed in several algorithms)
            self._laurent_polynomial_ring_extra_var = self._laurent_polynomial_ring['EXTRA_VAR']

        Parent.__init__(self, category=Rings())

    def ngens(self):
        return self.laurent_polynomial_ring().ngens()

    def free_module(self):
        return self._free_module

    def polynomial_ring(self):
        raise ValueError

    def polynomial_ring_extra_var(self):
        raise ValueError

    def laurent_polynomial_ring(self):
        return self._laurent_polynomial_ring

    def laurent_polynomial_ring_extra_var(self):
        return self._laurent_polynomial_ring_extra_var

    def with_extra_var(self):
        return self['EXTRA_VAR']

    @cached_method
    def zero(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.misc.multivariate_generating_series import MultivariateGeneratingSeriesRing

            sage: M = MultivariateGeneratingSeriesRing('x', 2)
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

            sage: from surface_dynamics.misc.multivariate_generating_series import MultivariateGeneratingSeriesRing

            sage: M = MultivariateGeneratingSeriesRing('x', 2)
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

            sage: from surface_dynamics.misc.multivariate_generating_series import MultivariateGeneratingSeriesRing

            sage: M = MultivariateGeneratingSeriesRing('x', 3)
            sage: M.term(1, [([1,1,0],1),([1,0,-1],2)])
            (1)/((1 - x0*x2^-1)^2*(1 - x0*x1))
            sage: M.term(1, {(1,1,0): 1, (1,0,2): 2})
            (1)/((1 - x0*x2^2)^2*(1 - x0*x1))

        Also works if the denominator is already a factored denominator::

            sage: from surface_dynamics.misc.factored_denominator import FactoredDenominator
            sage: M = MultivariateGeneratingSeriesRing('x', 3)
            sage: den = FactoredDenominator({(1,0,0): 2, (1,1,1):1})
            sage: M.term(1, den)
            (1)/((1 - x0)^2*(1 - x0*x1*x2))

        TESTS::


            sage: from surface_dynamics.misc.multivariate_generating_series import MultivariateGeneratingSeriesRing
            sage: M = MultivariateGeneratingSeriesRing('x', 2)
            sage: R = M.laurent_polynomial_ring()
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

            sage: from surface_dynamics.misc.multivariate_generating_series import MultivariateGeneratingSeriesRing

            sage: M = MultivariateGeneratingSeriesRing('x', 2)
            sage: M(1)
            (1)

            sage: R = M.laurent_polynomial_ring()
            sage: M(R.0)
            (x0)
        """
        num = self._laurent_polynomial_ring(arg)
        return self.element_class(self, [([], num)], self.free_module())

class GeneralizedMultiZetaRing(AbstractMSumRing):
    Element = GeneralizedMultiZetaElement

    def _coerce_map_from_(self, other):
        if self.laurent_polynomial_ring().base_ring().has_coerce_map_from(other):
            return True


# TODO: this should actually be an algebra over QQ[x0, x1, ..., xn]
# TODO: we should have two versions, as an algebra over the polynomial
#       ring QQ[x0, x1, ..., xn] and as an algebra over the Laurent
#       polynomial ring QQ[x0, x0^-1, x1, x1^-1, ..., xn, xn^-1]
class MultivariateGeneratingSeriesRing(AbstractMSumRing):
    r"""
    EXAMPLES::

        sage: from surface_dynamics.misc.multivariate_generating_series import MultivariateGeneratingSeriesRing

        sage: MultivariateGeneratingSeriesRing('x', 3)
        Multivariate quasi-polynomial generating series on x0, x1, x2

        sage: MultivariateGeneratingSeriesRing(['x', 'y', 'z'])
        Multivariate quasi-polynomial generating series on x, y, z

        sage: MultivariateGeneratingSeriesRing('a,b,c')
        Multivariate quasi-polynomial generating series on a, b, c

        sage: MultivariateGeneratingSeriesRing('y', 2)
        Multivariate quasi-polynomial generating series on y0, y1

        sage: M = MultivariateGeneratingSeriesRing('x', 3)
        sage: M.zero()
        0
        sage: M.one()
        (1)
        sage: m1 = M.term(1, [((1,1,0),1)])
        sage: m1
        (1)/((1 - x0*x1))
        sage: m2 = M.term(-2, [((1,0,0),1),((0,1,0),2)])
        sage: m2
        (-2)/((1 - x1)^2*(1 - x0))
    """
    Element = MultivariateGeneratingSeries

    @cached_method
    def residue_ring(self):
        return GeneralizedMultiZetaRing(self)

    def _repr_(self):
        vars_string = ', '.join(self.laurent_polynomial_ring().variable_names())
        return 'Multivariate quasi-polynomial generating series on ' + vars_string

    def _coerce_map_from_(self, other):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.misc.multivariate_generating_series import MultivariateGeneratingSeriesRing

            sage: M = MultivariateGeneratingSeriesRing('x', 2)

            sage: M.has_coerce_map_from(ZZ)
            True
            sage: M.coerce_map_from(ZZ)
            Co...ion map:
              From: Integer Ring
              To:   Multivariate quasi-polynomial generating series on x0, x1

            sage: M.has_coerce_map_from(QQ)
            True

            sage: M.has_coerce_map_from(M.laurent_polynomial_ring())
            True
        """
        if self.laurent_polynomial_ring().has_coerce_map_from(other):
            return True
