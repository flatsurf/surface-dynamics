r"""
Generalized multiple zeta values

Most functions in this module have a pair ``(n, den_tuple)`` as
arguments. ``n`` stands for the dimension (number of variables)
and ``den_tuple`` is a list of pairs ``(v, p)`` where ``v``
is a vector of length ``n`` with entries `\{0, 1\}` corresponding
to a linear form and ``p`` is a positive integer corresponding
to a power. For example, the standard mzv

.. MATH::

    \sum_{x,y,z \ge 1} \frac{1}{x (x + y) (x + y + z)^2}

would be encoded as ``[(1, (1, 0, 0)), (1, (1, 1, 0)), (2, (1, 1, 1))]``.
The ordering of the pairs ``(v, p)`` is irrelevant.
"""
#*****************************************************************************
#       Copyright (C) 2019-2023 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

import itertools
from collections import defaultdict

from sage.misc.cachefunc import cached_method, cached_function
from sage.all import ZZ, QQ, matrix, bernoulli_polynomial, prod, FreeModule
from sage.arith.all import binomial, factorial
from sage.sets.disjoint_set import DisjointSet
from sage.combinat.permutation import Permutations


from surface_dynamics.misc.linalg import linearly_independent_vectors
from .options import VERBOSE, DIVERGENT_MZV

import cypari2.handle_error

try:
    from sage.modular.multiple_zeta import Multizetas
except (ImportError, cypari2.handle_error.PariError):
    def Multizetas(*args, **kwds):
        raise ValueError('your sage version does not support multiple zeta values')

from surface_dynamics.misc.linalg import disjoint_vectors
from .gmzv_two_variables import Z2
from .gmzv_three_variables import Z3, Z3_sort_abc

# NOTE: ordering of variables (= inverse base 2)
#  1      2      3      4      5      6      7  ...
#  x      y      x+y    z      x+z    y+z    x+y+z  t      x+t   (y+t) (x+y+t) (z+t) ...
# (1000) (0100) (1100) (0010) (1010) (0110) (1110) (0001) (1001)
def linear_forms(d):
    r"""
    Return the `2^d - 1` linear forms with coefficients `\{0, 1\}` in ``d`` variables.

    EXAMPLES::

        sage: from surface_dynamics.generalized_multiple_zeta_values.generalized_multiple_zeta_values import linear_forms

        sage: list(linear_forms(1))
        [(1)]

        sage: list(linear_forms(2))
         [(1, 0), (0, 1), (1, 1)]

        sage: L1 = list(linear_forms(3))
        sage: L2 = L1[:]
        sage: L2.sort(key = lambda x: x[::-1])
        sage: assert L1 == L2
    """
    F = FreeModule(ZZ,d)
    for n in range(1, 2**d):
        v = F(ZZ(n).digits(2, padto=d))
        v.set_immutable()
        yield v


class DivergentZetaError(Exception):
    pass


def handle_term(n, den_tuple):
    r"""
    Return a linear combination of standard mzv equivalent to the gmzv ``(n, den_tuple)``.

    EXAMPLES::

        sage: from surface_dynamics.generalized_multiple_zeta_values.generalized_multiple_zeta_values import handle_term, is_convergent

        sage: M = Multizetas(QQ) # optional: mzv

        sage: V1 = FreeModule(ZZ, 1)
        sage: v = V1((1,)); v.set_immutable()
        sage: dt = ((v,3),)
        sage: assert is_convergent(1, dt) and handle_term(1, dt) == M((3,)) # optional: mzv


        sage: V2 = FreeModule(ZZ, 2)
        sage: va = V2((1,0)); va.set_immutable()
        sage: vb = V2((0,1)); vb.set_immutable()
        sage: vc = V2((1,1)); vc.set_immutable()
        sage: dt = ((va,2), (vc,3))
        sage: assert is_convergent(2, dt) and handle_term(2, dt) == M((2,3)) # optional: mzv
        sage: dt1 = ((va,2),(vb,3))
        sage: dt2 = ((va,3),(vb,2))
        sage: assert is_convergent(2,dt1) and is_convergent(2,dt2) # optional: mzv
        sage: assert handle_term(2, ((va,2), (vb,3))) == handle_term(2, ((va,3), (vb,2))) == M((2,)) * M((3,)) # optional: mzv

        sage: V3 = FreeModule(ZZ, 3)
        sage: va = V3((1,0,0)); va.set_immutable()
        sage: vb = V3((0,1,0)); vb.set_immutable()
        sage: vc = V3((0,0,1)); vc.set_immutable()
        sage: vd = V3((1,1,0)); vd.set_immutable()
        sage: ve = V3((1,0,1)); ve.set_immutable()
        sage: vf = V3((0,1,1)); vf.set_immutable()
        sage: vg = V3((1,1,1)); vg.set_immutable()
        sage: assert handle_term(3, ((va,2), (vd,3), (vg,4))) == M((2,3,4)) # optional: mzv
        sage: assert handle_term(3, ((va,2), (vb,3), (vc,4))) == handle_term(3, ((va,3), (vb,2), (vc,4))) # optional mzv
        sage: assert handle_term(3, ((va,2), (vb,3), (vc,4))) == M((2,)) * M((3,)) * M((4,)) # optional: mzv
        sage: assert handle_term(3, ((va,1), (vc,2), (vd,3))) == handle_term(3, ((va,1), (vb,2), (ve,3))) # optional: mzv
        sage: assert handle_term(3, ((va,1), (vc,2), (vd,3))) == handle_term(3, ((va,2), (vb,1), (vf,3))) # optional: mzv
        sage: assert handle_term(3, ((va,1), (vc,2), (vd,3))) ==  M((2,)) * M((1,3)) # optional: mzv
    """
    if n == 1:
        M = Multizetas(QQ)
        ans = M.zero()
        for v, p in den_tuple:
            # sum (v[0]x)^p
            ans += QQ.one() / v[0] **p * M((p,))
        return ans

    elif n == 2:
        dat = to_Z2(den_tuple)
        if dat is None:
            raise NotImplementedError("generalized mzv {}".format(den_tuple))
        return Z2(*dat)

    elif n == 3:
        dat = to_Z3(den_tuple)
        if dat is None:
            raise NotImplementedError("generalized mzv {}".format(den_tuple))
        return Z3(*dat)

    return _handle_term(n, den_tuple)


@cached_function
def _handle_term(n, den_tuple):
    if VERBOSE:
        print("handle term({}, {})".format(n, den_tuple))
    assert all(len(v) == n for v,p in den_tuple), (n, den_tuple)

    if any(x < 0 or x > 1 for v,p in den_tuple for x in v):
        raise ValueError("unhandled zeta values {}".format(den_tuple))

    # 0. check convergence
    if not DIVERGENT_MZV and not is_convergent(n, den_tuple):
        raise DivergentZetaError("{} Z({})".format(n, den_tuple))

    # 1. multizeta
    Z = is_multizeta(n, den_tuple)
    if Z is not None:
        if VERBOSE:
            print("standard multizeta")
        return Z

    # 2. apply stuffle
    P = is_stufflisable(n, den_tuple)
    if P is not None:
        if VERBOSE:
            print("stuffle")
        return sum(handle_term(nn, dd) for nn, dd in P)

    # 3. equal rows
    # HOW?

    # 4. diminish the number of linear forms without creating convergence problem
    data = has_term_sum_of_smaller_terms(n, den_tuple)
    if data is not None:
        if VERBOSE:
            print("relation between linear forms: L_i = sum L_j with L_i > L_j")
        return sum(coeff * handle_term(n, new_den_tuple) for coeff, new_den_tuple in kill_relation(n, den_tuple, data[0], data[1]))

    # 5. make "big terms"
    data = is_reducible(n, den_tuple)
    if data is not None:
        if VERBOSE:
            print("reduction")
        return sum(coeff * handle_term(n, new_den_tuple) for coeff, new_den_tuple in data)

    # 5. 3-variables
    if n == 3:
        dat = to_Z3(den_tuple)
        if dat is None:
            raise NotImplementedError("generalized mzv {}".format(den_tuple))
        return Z3(*dat)

    raise NotImplementedError("unhnandled generalized multiple zeta value {}".format(den_tuple))


def clean_term(n, den_tuple):
    D = {}
    for den, p in den_tuple:
        if den in D:
            D[den] += p
            if not D[den]:
                del D[den]
        else:
            D[den] = p
    return tuple(sorted(D.items()))


def to_Z2(den_tuple):
    r"""
    Converts ``den_tuple`` to arguments ``(a, b, c)`` for the function `Z2`.
    """
    a = b = c = 0
    for v,p in den_tuple:
        v = tuple(v)
        if len(v) != 2:
            raise ValueError
        if v == (1,0):
            a += p
        elif v == (0,1):
            b += p
        elif v == (1,1):
            c += p
        else:
            return
    return a,b,c


def to_Z3(den_tuple, sort=True):
    r"""
    Converts ``den_tuple`` to arguments ``(a, b, c, d, e, f, g)`` for the function `Z3`.
    """
    a = b = c = d = e = f = g = 0
    for v,p in den_tuple:
        v = tuple(v)
        if len(v) != 3:
            raise ValueError
        if v == (1,0,0):
            a += p
        elif v == (0,1,0):
            b += p
        elif v == (0,0,1):
            c += p
        elif v == (1,1,0):
            d += p
        elif v == (1,0,1):
            e += p
        elif v == (0,1,1):
            f += p
        elif v == (1,1,1):
            g += p
        else:
            return
    # now normalize
    if sort:
        a,b,c,d,e,f,g = Z3_sort_abc(a,b,c,d,e,f,g)
        if b == c:
            if a == b:
                d,e,f = sorted([d,e,f], reverse=True)
            else:
                d,e = sorted([d,e], reverse=True)
    return a,b,c,d,e,f,g


def linear_form(R, v):
    return sum(R.gen(i) for i,j in enumerate(v) if j)


def negative_rays(n):
    l = []
    v = [0]*n
    for i in range(n):
        v[i] = -1
        l.append(v[:])
        v[i] = 0
    return l


def is_convergent(n, den_tuple):
    r"""
    Test whether the GMZV is convergent.

    TESTS::

        sage: import itertools
        sage: from surface_dynamics.generalized_multiple_zeta_values.generalized_multiple_zeta_values import is_convergent
        sage: V = FreeModule(ZZ, 3)
        sage: va = V((1,0,0)); va.set_immutable()
        sage: vb = V((0,1,0)); vb.set_immutable()
        sage: vc = V((0,0,1)); vc.set_immutable()
        sage: vd = V((1,1,0)); vd.set_immutable()
        sage: ve = V((1,0,1)); ve.set_immutable()
        sage: vf = V((0,1,1)); vf.set_immutable()
        sage: vg = V((1,1,1)); vg.set_immutable()
        sage: gens = [va,vb,vc,vd,ve,vf,vg]
        sage: N = 0
        sage: for p in itertools.product([0,1,2], repeat=7): # optional: mzv
        ....:     if sum(map(bool,p)) == 3 and is_convergent(3, list(zip(gens,p))):
        ....:         print(p)
        ....:         N += 1
        (0, 0, 0, 0, 1, 1, 2)
        (0, 0, 0, 0, 1, 2, 1)
        (0, 0, 0, 0, 1, 2, 2)
        (0, 0, 0, 0, 2, 1, 1)
        (0, 0, 0, 0, 2, 1, 2)
        (0, 0, 0, 0, 2, 2, 1)
        (0, 0, 0, 0, 2, 2, 2)
        (0, 0, 0, 1, 0, 1, 2)
        ...
        (2, 0, 2, 0, 0, 2, 0)
        (2, 0, 2, 2, 0, 0, 0)
        (2, 1, 0, 0, 0, 0, 2)
        (2, 1, 0, 0, 0, 2, 0)
        (2, 2, 0, 0, 0, 0, 2)
        (2, 2, 0, 0, 0, 2, 0)
        (2, 2, 0, 0, 2, 0, 0)
        (2, 2, 2, 0, 0, 0, 0)
        sage: print(N) # optional: mzv
        125
    """
    from sage.geometry.polyhedron.constructor import Polyhedron
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

    assert all(len(v) == n for v, p in den_tuple), (n, den_tuple)

    # TODO: fast code path

    R = PolynomialRing(QQ, 'x', n)
    x = R.gens()
    den_poly = prod(linear_form(R, v)**p for v, p in den_tuple)
    newton_polytope = Polyhedron(vertices=den_poly.exponents(), rays=negative_rays(n))
    V = newton_polytope.intersection(Polyhedron(rays=[[1]*n])).vertices()
    r = max(max(v.vector()) for v in V)
    return r > 1


def convergent_multizeta(t):
    r"""
    Multizeta value at a convergent index ``t``.

    TESTS::

        sage: from surface_dynamics.generalized_multiple_zeta_values.generalized_multiple_zeta_values import convergent_multizeta
        sage: assert all(convergent_multizeta(t) == Multizeta(*t) for t in [(2,),(3,),(1,2),(3,2),(1,1,2)]) # optional: mzv

        sage: convergent_multizeta((0,3)) # optional: mzv
        ζ(2) - ζ(3)
        sage: convergent_multizeta((0,2,2)) # optional: mzv
        ζ(1,2) - ζ(2,2)
        sage: convergent_multizeta((1,0,3)) # optional: mzv
        ζ(1,2) - ζ(1,3) - ζ(2) + ζ(3)
        sage: convergent_multizeta((0,1,3)) # optional: mzv
        -ζ(1,3) + ζ(2) - ζ(3)

        sage: convergent_multizeta((0, 4)) # optional: mzv
        ζ(3) - ζ(4)
        sage: convergent_multizeta((-1, 5)) # optional: mzv
        1/2*ζ(3) - 1/2*ζ(4)
        sage: convergent_multizeta((-2, 5)) # optional: mzv
        1/3*ζ(2) - 1/2*ζ(3) + 1/6*ζ(4)

        sage: convergent_multizeta((-1, 3, 4)) # optional: mzv
        1/2*ζ(1,4) - 1/2*ζ(2,4)
        sage: convergent_multizeta((-1, -1, 8)) # optional: mzv
        1/8*ζ(4) - 5/12*ζ(5) + 3/8*ζ(6) - 1/12*ζ(7)
        sage: convergent_multizeta((-2,-2,10)) # optional: mzv
        1/18*ζ(4) - 4/15*ζ(5) + 31/72*ζ(6) - 1/4*ζ(7) + 1/72*ζ(8) + 1/60*ζ(9)
        sage: convergent_multizeta((-1, -2, 10)) # optional: mzv
        1/10*ζ(5) - 3/8*ζ(6) + 5/12*ζ(7) - 1/8*ζ(8) - 1/60*ζ(9)

        sage: convergent_multizeta((4,-4,10)) # optional: mzv
        -1/3*ζ(1,10) + 1/30*ζ(3,10) + 1/5*ζ(4,5) - 1/2*ζ(4,6) + 1/3*ζ(4,7) - 1/30*ζ(4,9) - 1/10*ζ(8) - 2/5*ζ(9) + 1/2*ζ(10)
        sage: convergent_multizeta((2,-1,8)) # optional: mzv
        -1/2*ζ(1,8) + 1/2*ζ(2,6) - 1/2*ζ(2,7) - 1/2*ζ(7) + 1/2*ζ(8)

        sage: convergent_multizeta((0,3,2)) # optional: mzv
        ζ(2,2) - ζ(3,2)

    Divergent cases::

        sage: convergent_multizeta((0,2)) # optional: mzv
        Traceback (most recent call last):
        ...
        DivergentZetaError: divergent multizeta value (0, 2)
        sage: convergent_multizeta((0,0,3)) # optional: mzv
        Traceback (most recent call last):
        ...
        DivergentZetaError: divergent multizeta value (0, 0, 3)
        sage: convergent_multizeta((1,0,2)) # optional: mzv
        Traceback (most recent call last):
        ...
        DivergentZetaError: divergent multizeta value (1, 0, 2)
        sage: convergent_multizeta((0,1,2)) # optional: mzv
        Traceback (most recent call last):
        ...
        DivergentZetaError: divergent multizeta value (0, 1, 2)
    """
    if VERBOSE:
        print("convergent_multizeta({})".format(t))
    n = len(t)
    for i in range(1,n+1):
        if sum(t[-i:]) <= i:
            raise DivergentZetaError("divergent multizeta value {}".format(t))
    if all(x > 0 for x in t):
        M = Multizetas(QQ)
        W = M.basis().keys()
        return M.term(W(t))
    else:
        # find the first non-positive index and drop the summation on the
        # corresponding variable using Faulhaber's formula
        i = 0
        while t[i] > 0:
            i += 1
        t = list(t)
        tt = t[:i] + t[i+1:]
        x = QQ['x'].gen()
        if t[i] == 0:
            # bug: bernoulli_polynomial(x, 0) is an integer
            P = x - 1
        else:
            P = bernoulli_polynomial(QQ['x'].gen(), -t[i]).integral()
        s = Multizetas(QQ).zero()
        for c,e in zip(P.coefficients(), P.exponents()):
            tt[i] -= e
            s += c * convergent_multizeta(tt[:])
            tt[i] += e
        if i > 0:
            Q = P(x+1)
            for c,e in zip(Q.coefficients(), Q.exponents()):
                tt[i-1] -= e
                s -= c * convergent_multizeta(tt[:])
                tt[i-1] += e
        return s


def is_multizeta(n, den_tuple):
    r"""
    Return the corresponding zeta values if it is one.

    EXAMPLES::

        sage: from surface_dynamics.generalized_multiple_zeta_values.generalized_multiple_zeta_values import linear_forms, is_multizeta
        sage: va, vb, vd, vc, ve, vf, vg = linear_forms(3)

        sage: is_multizeta(3, ((vb, 2), (vf, 5), (vg, 2))) # optional: mzv
        ζ(2,5,2)

        sage: is_multizeta(3, [(vg, 5)]) # optional: mzv
        1/2*ζ(3) - 3/2*ζ(4) + ζ(5)

        sage: assert is_multizeta(3, ((va,3), (vb,3), (vc,3))) is None # optional: mzv
        sage: assert is_multizeta(3, ((vb,2), (ve,5), (vg,1))) is None # optional: mzv
    """
    assert all(len(v) == n for v,p in den_tuple), (n, den_tuple)

    by_support = [None for _ in range(n+1)]
    for vp in den_tuple:
        v,p = vp
        if any(x > 1 for x in v):
            return
        m = sum(vp[0])
        if by_support[m] is None:
            by_support[m] = vp
        else:
            return

    k = 0
    linear_order = [None] * n
    for step, vp in enumerate(by_support):
        if vp is None:
            continue
        v,p = vp
        # discover new variables
        for i in range(n):
            if v[i] and linear_order[i] is None:
                linear_order[i] = k
                k += 1
        if k > step:
            return

    if all(v is not None for v in linear_order[1:]):
        return convergent_multizeta([0 if vp is None else vp[1] for vp in by_support[1:]])


def is_stufflisable(n, den_tuple):
    m = len(den_tuple)
    possibilities = []
    for i in range(n):
        for j in range(i):
            if all(not (v[i] and v[j]) for v,p in den_tuple):
                # see if we can perform a stuffle that will lower rank in next step
                for l in range(n):
                    if all(v[i] + v[j] == v[l] for v,p in den_tuple):
                        return stuffle(n, den_tuple, i, j)
                possibilities.append((i,j))

    # pick the first candidate ?
    # multiple columns shuffling ?
    if possibilities:
        i, j = possibilities[0]
        return stuffle(n, den_tuple, i, j)


def stuffle(n, den_tuple, i, j):
    r"""
    Apply a stuffle relation for the variables i and j.
    """
    m = matrix([v for v,p in den_tuple]).transpose()

    # big subsimplex 1
    mm = m.__copy__()
    mm[i] += mm[j]
    D1 = clean_term(n, tuple((c,p) for c,(v,p) in zip(mm.columns(), den_tuple)))

    # big subsimplex 2
    mm = m.__copy__()
    mm[j] += mm[i]
    D2 = clean_term(n, tuple((c,p) for c,(v,p) in zip(mm.columns(), den_tuple)))

    # small subsimplex
    mm = m.__copy__()
    mm[i] += mm[j]
    mm = mm.delete_rows([j])
    D3 = clean_term(n-1, tuple((c,p) for c,(v,p) in zip(mm.columns(), den_tuple)))

    return [(n,D1),(n,D2),(n-1,D3)]


def is_product(n, den_tuple):
    r"""
    Test whether the GMZV is a product of GMZV of smaller dimensions. If so, return
    a pair of linear combinations whose product is equivalent.

    INPUT:

    - ``n`` - number of variables

    - ``den_tuple`` - tuple of pairs ``(vector, power)``

    TESTS::

        sage: from surface_dynamics.generalized_multiple_zeta_values.generalized_multiple_zeta_values import is_product
        sage: term1, term2 = is_product(3, [((1,0,0),2), ((0,1,0),5), ((1,1,0),1), ((0,0,1),3)])
        sage: term1
        (2, (((0, 1), 5), ((1, 0), 2), ((1, 1), 1)))
        sage: term2
        (1, (((1), 3),))

        sage: is_product(3, [((1,0,0),2), ((0,1,0),3), ((1,0,1),1), ((0,0,1),5)])
        [(2, (((0, 1), 5), ((1, 0), 2), ((1, 1), 1))), (1, (((1), 3),))]

        sage: is_product(3, [((1,1,1),3)]) is None
        True
    """
    D = DisjointSet(n)
    assert all(len(v) == n for v,p in den_tuple), (n, den_tuple)

    # 1. product structure
    for v,_ in den_tuple:
        i0 = 0
        while not v[i0]:
            i0 += 1
        i = i0 + 1
        while i < n:
            if v[i]:
                D.union(i0, i)
            i += 1
        if D.number_of_subsets() == 1:
            # no way to split variables
            return

    # split variables
    Rdict = D.root_to_elements_dict()
    keys = sorted(Rdict.keys())
    key_indices = {k: i for i,k in enumerate(keys)}
    values = [Rdict[k] for k in keys]
    values_indices = [{v:i for i,v in enumerate(v)} for v in values]
    n_list = [len(J) for J in values]
    F = [FreeModule(ZZ, nn) for nn in n_list]
    new_terms = [[] for _ in range(len(Rdict))]
    for v,p in den_tuple:
        i0 = 0
        while not v[i0]:
            i0 += 1
        i0 = D.find(i0)
        assert all(D.find(i) == i0 for i in range(n) if v[i]), (i0, [D.find(i) for i in range(n) if v[i]])
        k = key_indices[i0]
        vv = F[k]()
        for i in range(n):
            if v[i]:
                vv[values_indices[k][i]] = v[i]
        vv.set_immutable()
        new_terms[k].append((vv,p))

    return list(zip(n_list, [tuple(sorted(terms)) for terms in new_terms]))

# TODO: apply a power of a relation at once with multinomial coefficients
def apply_relation(n, den_tuple, i, relation):
    # iterator through the pairs (coeff, den_tuples) obtained by applying relation on the i-th term
    ans = []
    for j in range(n):
        if relation[j] and j != i:
            term = list(den_tuple)
            v, p = term[i]
            term[i] = (v, p+1)
            v, p = term[j]
            term[j] = (v, p-1)
            yield (relation[j], tuple(term))


def kill_relation(n, den_tuple, i, relation):
    r"""
    Make calls to :func:`apply_relation` until we get rid of term.

    EXAMPLES::

        sage: from surface_dynamics.generalized_multiple_zeta_values.generalized_multiple_zeta_values import kill_relation

        sage: V = FreeModule(ZZ, 2)
        sage: va = V((1,0)); va.set_immutable()
        sage: vb = V((0,1)); vb.set_immutable()
        sage: vc = V((1,1)); vc.set_immutable()

        sage: den_tuple = ((va,2), (vb,3), (vc,4))
        sage: kill_relation(2, den_tuple, 2, [1,1,0])
        [(1, (((0, 1), 3), ((1, 1), 6))),
         (2, (((0, 1), 2), ((1, 1), 7))),
         (1, (((1, 0), 2), ((1, 1), 7))),
         (3, (((0, 1), 1), ((1, 1), 8))),
         (3, (((1, 0), 1), ((1, 1), 8)))]
    """
    assert len(relation) == len(den_tuple)
    assert 0 <= i < len(den_tuple)
    s = sum(relation[j] * den_tuple[j][0] for j in range(len(den_tuple)))
    assert s == den_tuple[i][0], (s, den_tuple[i][0])
    D = {den_tuple: QQ.one()}
    todo = [den_tuple]
    while todo:
        den_tuple = todo.pop(0)
        coeff1 = D.pop(den_tuple)
        for coeff, new_den_tuple in apply_relation(n, den_tuple, i, relation):
            coeff *= coeff1
            if new_den_tuple in D:
                D[new_den_tuple] += coeff
            elif any(not new_den_tuple[j][1] for j in range(len(den_tuple))):
                new_den_tuple = tuple(x for x in new_den_tuple if x[1])
                if new_den_tuple not in D:
                    D[new_den_tuple] = coeff
                else:
                    D[new_den_tuple] += coeff
            else:
                todo.append(new_den_tuple)
                D[new_den_tuple] = coeff

    return [(coeff, new_den_tuple) for new_den_tuple, coeff in D.items()]


def try_relation(n, den_tuple):
    r"""
    Assuming that ``den_tuple`` has full rank, make it so that it has only ``n`` columns.

    EXAMPLES::

        sage: from surface_dynamics.generalized_multiple_zeta_values.generalized_multiple_zeta_values import try_relation
        sage: V = FreeModule(ZZ, 2)
        sage: va = V((1,0)); va.set_immutable()
        sage: vb = V((0,1)); vb.set_immutable()
        sage: vc = V((1,1)); vc.set_immutable()

        sage: den_tuple = ((va,2), (vb,2), (vc,2))
        sage: for r in try_relation(2, den_tuple):
        ....:     print(r)
        [(1, (((0, 1), 2), ((1, 1), 4))), (1, (((1, 0), 2), ((1, 1), 4))), (2, (((0, 1), 1), ((1, 1), 5))), (2, (((1, 0), 1), ((1, 1), 5)))]
    """
    # assume it is full rank
    if len(den_tuple) <= n:
        return
    M = matrix([v for v,p in den_tuple])
    for relation in sorted(M.left_kernel().basis(), key=lambda x: sum(bool(cc) for cc in x)):
        if sum(x < 0 for x in relation) > sum(x > 0 for x in relation):
            relation = -relation
        for i in range(len(den_tuple)):
            if relation[i] < 0:
                relation /= -relation[i]
                relation[i] = 0
                yield kill_relation(n, den_tuple, i, relation)


def has_term_sum_of_smaller_terms(n, den_tuple):
    r"""
    Look for a vector v_i and {v_j} with sum(v_j) = v_i

    This is useful only useful when the linear forms have relations between them.

    EXAMPLES::

        sage: from surface_dynamics.generalized_multiple_zeta_values.generalized_multiple_zeta_values import linear_forms, has_term_sum_of_smaller_terms

        sage: va,vb,vc = linear_forms(2)
        sage: has_term_sum_of_smaller_terms(2, ((va,1),(vb,1),(vc,1)))
        [2, [1, 1, 0]]
        sage: assert has_term_sum_of_smaller_terms(2, ((va,1),(vc,1))) is None
        sage: assert has_term_sum_of_smaller_terms(2, ((vb,1),(vc,1))) is None
        sage: assert has_term_sum_of_smaller_terms(2, ((va,1),(vb,1))) is None

        sage: va, vb, vd, vc, ve, vf, vg = linear_forms(3)
        sage: has_term_sum_of_smaller_terms(3, ((va,1),(vb,1),(vc,1),(vg,1)))
        [3, [1, 1, 1, 0]]
        sage: has_term_sum_of_smaller_terms(3, ((va,1),(vb,1),(ve,1),(vf,1),(vg,1)))
        [4, [0, 1, 1, 0, 0]]
        sage: has_term_sum_of_smaller_terms(3, ((vd,1),(ve,1),(vf,1),(vg,1)))
        [3, [1/2, 1/2, 1/2, 0]]
        sage: has_term_sum_of_smaller_terms(3, ((va,1),(vd,1),(ve,1),(vg,1)))
        [3, [-1, 1, 1, 0]]
    """
    l = len(den_tuple)
    for iv in range(l):
        v = den_tuple[iv][0]
        if sum(v) == 1:
            continue

        candidates = [k for k,(u,_) in enumerate(den_tuple) if sum(u) < sum(v) and all(u[i] <= v[i] for i in range(n))]
        M = matrix(QQ, [den_tuple[k][0] for k in candidates])
        try:
            coeffs = M.solve_left(v)
        except ValueError:
            continue
        full = [0] * l
        for k,c in zip(candidates, coeffs):
            full[k] = c
        return [iv, full]


def write_array(level, explanations):
    for i,key in enumerate(sorted(level, key=lambda x: 1000 if level[x] is None else level[x])):
        print("%s & %s \\\\" % (lin_prod(key),'\\infty' if level[key] is None else level[key]))


def is_reducible(n, den_tuple):
    r"""
    Test whether the given GMZV is reducible and if so return an equivalent linear expression of simpler GMZV.

    We call a GMZV is reducible If (x1+x2+...+xd) is not present, use a linear relation to create it. Then
    try to kill using other forms. Then try to write as P(x1,x2,...,x_{d-1}) (x1+x2+...+xd) and recurse.

    For example, this function expresses all Tornheim sum as linear combinations of MZV.

    EXAMPLES::

        sage: from surface_dynamics.generalized_multiple_zeta_values.generalized_multiple_zeta_values import linear_forms, is_reducible
        sage: va, vb, vd, vc, ve, vf, vg = linear_forms(3)
        sage: is_reducible(3, ((va, 3), (vb, 3), (vc, 3)))
        [(1, (((0, 0, 1), 3), ((0, 1, 1), 3), ((1, 1, 1), 3))),
         (1, (((0, 1, 0), 3), ((0, 1, 1), 3), ((1, 1, 1), 3))),
         (3, (((0, 0, 1), 2), ((0, 1, 1), 4), ((1, 1, 1), 3))),
         (3, (((0, 1, 0), 2), ((0, 1, 1), 4), ((1, 1, 1), 3))),
        ...
         (90, (((1, 0, 0), 1), ((1, 0, 1), 1), ((1, 1, 1), 7))),
         (90, (((0, 1, 0), 1), ((1, 1, 0), 1), ((1, 1, 1), 7))),
         (90, (((1, 0, 0), 1), ((1, 1, 0), 1), ((1, 1, 1), 7)))]
    """
    if len(den_tuple) == 1:
        return

    # force the max vector (1,1,...,1) to appear possibly with an exponent 0
    F = FreeModule(ZZ, n)
    vmax = F([1] * n)
    vmax.set_immutable()
    imax = None
    for i, (v, p) in enumerate(den_tuple):
        if v == vmax:
            imax = i
            break
    if imax is None:
        imax = len(den_tuple)
        den_tuple = den_tuple + ((vmax, 0),)
    if imax != len(den_tuple) - 1:
        den_tuple = list(den_tuple)
        den_tuple.append(den_tuple.pop(imax))
        den_tuple = tuple(den_tuple)
        imax = len(den_tuple) - 1

    assert den_tuple[imax][0] == vmax

    M = matrix(QQ, [v for v,p in den_tuple if v != vmax])
    try:
        relation = M.solve_left(vmax)
    except ValueError:
        if den_tuple[-1][1] == 0:
            return
        den_tuple2 = den_tuple[:-1]
        variables = set().union(*[[i for i in range(n) if v[i]] for v,p in den_tuple2])
        if len(variables) == n:
            return
        killed = [(1, den_tuple)]
    else:
        killed = kill_relation(n, den_tuple, imax, list(relation) + [0])

    ans = defaultdict(QQ)
    for coeff, den_tuple2 in killed:
        assert den_tuple2[-1][0] == vmax
        pmax = den_tuple2[-1][1]
        assert pmax
        den_tuple2 = den_tuple2[:-1]

        variables = set().union(*[[i for i in range(n) if v[i]] for v,p in den_tuple2])
        if len(variables) == n:
            # removing the maximal vector (1,1,...,1) is not enough to make a variable disappear
            data = is_reducible(n, den_tuple2)
            if data is None:
                data = [(1, den_tuple2)]
        else:
            # less variables!
            nn = len(variables)
            variables = sorted(variables)
            new_indices = {j:i for i,j in enumerate(variables)}
            G = FreeModule(ZZ, nn)
            new_den_tuple2 = []
            for v,p in den_tuple2:
                vv = G([v[i] for i in variables])
                vv.set_immutable()
                new_den_tuple2.append((vv,p))
            new_den_tuple2 = tuple(new_den_tuple2)
            data = is_reducible(nn, new_den_tuple2)
            if data is None:
                data = [(1, new_den_tuple2)]

            # lift to the n variables version
            new_data = []
            for coeff3, den_tuple3 in data:
                den_tuple3 = [(F([v[new_indices[j]] if j in new_indices else 0 for j in range(n)]), p) for v,p in den_tuple3]
                for v,p in den_tuple3:
                    v.set_immutable()
                new_data.append((coeff3, tuple(sorted(den_tuple3))))
            data = new_data

        # update the answer
        for coeff3, den_tuple3 in data:
            imax = None
            for i,(v,p) in enumerate(den_tuple3):
                if v == vmax:
                    imax = i
                    break
            if imax is None:
                den_tuple3 = den_tuple3 + ((vmax,pmax),)
            else:
                den_tuple3 = den_tuple3[:imax] + ((vmax, den_tuple3[imax][1] + pmax),) + den_tuple3[imax+1:]
            ans[den_tuple3] += coeff * coeff3

    if len(ans) > 1:
        return [(coeff, den_tuple) for den_tuple, coeff in ans.items()]


class GeneralizedMultipleZetaFunction:
    r"""
    Generalized multiple zeta function.

    A generalized multiple zeta function is a series of the form

    .. MATH::

        Z(s_1, \ldots, s_r) = \sum_{x_1, x_2, \ldots, x_d \geq 1} L_1(x)^{-s_1} L_2(x)^{-s_2} \cdots L_r(x)^{-s_r}

    where the sum is over integers $x_i$ and the $L_i$ are linear forms in `x =
    (x_1, x_2, \ldots, x_d)` with `\{0,1\}` coefficients. It could be thought of
    as as generalized polylogarithm in the variable `s = (s_1, \ldots, s_d)`
    evaluated at `z = (1, 1, \ldots, 1)`.

    EXAMPLES::

        sage: from surface_dynamics.generalized_multiple_zeta_values.generalized_multiple_zeta_values import linear_forms, GeneralizedMultipleZetaFunction

        sage: va, vb, vd, vc, ve, vf, vg = linear_forms(3)
        sage: f = GeneralizedMultipleZetaFunction([va, vb, vc, vg])
        sage: f
        GeneralizedMultipleZetaFunction('(0)(1)(2)(0+1+2)')
        sage: f((1, 1, 1, 1)) # optional - mzv
        6*ζ(1,1,2)
    """
    def __init__(self, *args, as_rows=False, reduced=None):
        if len(args) == 1 and isinstance(args[0], str):
            s = args[0].replace(' ', '').replace(')*(', ')(')
            if not s.startswith('(') or not s.endswith(')'):
                raise ValueError
            s = [lin.split('+') for lin in s[1:-1].split(')(')]
            letters = set().union(*s)
            indices = {letter:i for i,letter in enumerate(sorted(letters))}
            F = FreeModule(ZZ, len(indices))
            args = ([sum(F.gen(indices[letter]) for letter in lin) for lin in s],)

        self._entries = matrix(ZZ, *args)
        if not as_rows:
            self._entries = self._entries.transpose()
        self._reduced = reduced

    def children(self):
        r"""
        Compatibility with function with recursive calls.
        """
        return []

    @classmethod
    def from_matrix(cls, m, canonicalize, sort_rows, set_immutable, reduced):
        L = GeneralizedMultipleZetaFunction.__new__(GeneralizedMultipleZetaFunction)
        L._reduced = reduced
        if canonicalize:
            L._entries = m.permutation_normal_form()
            nr = m.nrows()
            nc = m.ncols()
        elif sort_rows:
            L._entries = matrix(ZZ, sorted(m.rows(), reverse=True))
        else:
            L._entries = matrix(ZZ, m)

        if set_immutable:
            L._entries.set_immutable()
        return L

    def __getitem__(self, key):
        return self._entries[key]

    def __eq__(self, other):
        if type(self) is not type(other):
            raise TypeError
        return self._entries == other._entries

    def __ne__(self, other):
        if type(self) is not type(other):
            raise TypeError
        return self._entries != other._entries

    def __call__(self, e):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.generalized_multiple_zeta_values import GeneralizedMultipleZetaFunction, convergent_multizeta
            sage: f = GeneralizedMultipleZetaFunction([[1,0,0],[1,1,1],[1,0,1]])
            sage: f((2,4,3)) # optional - mzv
            ζ(2,3,4)
            sage: f = GeneralizedMultipleZetaFunction([[1,0,0],[1,1,1]])
            sage: f((2,4)) == convergent_multizeta((2,0,4)) # optional - mzv
            True

            sage: G = GeneralizedMultipleZetaFunction([[1,1,1,1],[1,1,1,0],[1,1,0,0]])
            sage: G((2,2,2)) # optional - mzv
            ζ(1,2,2) - ζ(2,2,2)
        """
        M = Multizetas(QQ)
        d = self.nrows()
        n = self.ncols()
        if len(e) != n:
            raise ValueError('invalid argument: generalized mzv with {} columns but got e={}'.format(n, e))
        if d == 1:
            # pure zeta
            return M((sum(e),))
        elif d == 2:
            return Z2(*to_Z2(list(zip(self.columns(), e))))
        elif d == 3:
            return Z3(*to_Z3(list(zip(self.columns(), e))))
        
        dat = self.is_multizeta()
        if dat:
            s = [0] * self.nrows()
            for x, c in zip(e, self.columns()):
                s[sum(c)-1] += x
            return convergent_multizeta(s)

        raise NotImplementedError('generalized mzv {}'.format(self.lin_prod_string()))
    
    def copy(self, immutable=False):
        if immutable and self._entries.is_immutable():
            return self
        L = GeneralizedMultipleZetaFunction.__new__(GeneralizedMultipleZetaFunction)
        L._reduced = self._reduced
        L._entries = self._entries.__copy__()
        if immutable:
            L._entries.set_immutable()
        return L

    def __repr__(self):
        return "GeneralizedMultipleZetaFunction({!r})".format(self.lin_prod_string())

    def lin_prod_string(self): 
        nr = self.nrows()
        return ''.join('(' + '+'.join(str(i) for i in range(nr) if c[i]) + ')' for c in self._entries.columns())    

    def subs(self, substitution_dict, canonicalize=False, set_immutable=False, reduced=False):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.generalized_multiple_zeta_values import GeneralizedMultipleZetaFunction
            sage: L = GeneralizedMultipleZetaFunction([[1,0,0],[0,1,0],[0,0,1]], as_rows=True)
            sage: L.subs({0: [(1,0),(1,1)], 1: [(1,0),(1,1),(1,2)]})
            GeneralizedMultipleZetaFunction('(0+1)(0+1)(1+2)')
        """
        M = self._entries.__copy__()
        for i, lin_comb in substitution_dict.items():
            M[i] = sum(c * self._entries[j] for c,j in lin_comb)
        return GeneralizedMultipleZetaFunction.from_matrix(M, canonicalize, False, set_immutable, reduced)

    def set_immutable(self):
        self._entries.set_immutable()

    @cached_method
    def is_reduced(self):
        if self._reduced is None:
            self._reduced = self._entries.column_module().rank() == self._entries.ncols()
        return self._reduced

    def __hash__(self):
        return hash(self._entries)

    def __eq__(self, other):
        if type(self) is not type(other):
            raise TypeError
        return self._entries == other._entries
    def __ne__(self, other):
        if type(self) is not type(other):
            raise TypeError
        return self._entries == other._entries
    def nrows(self):
        return self._entries.nrows()
    def ncols(self):
        return self._entries.ncols()
    def rows(self):
        return self._entries.rows()
    def columns(self):
        return self._entries.columns()
    def delete_row(self, i):
        self._entries = self._entries.delete_rows([i])
    def delete_column(self, i):
        self._entries = self._entries.delete_columns([i])
    def symmetric(self):
        r"""
        Return all possible symmetric version obtained by permuting variables
        """
        d = self.ncols()
        for p in Permutations(self.columns()):
            yield GeneralizedMultipleZetaFunction.from_matrix(matrix(ZZ, p).transpose(), False, True, True, self._reduced)

    def has_no_zero_row_or_column(self):
        M = self._entries
        return all(any(self._entries[i][j] for j in range(M.ncols())) for i in range(M.nrows())) and all(any(self._entries[i][j] for i in range(M.nrows())) for j in range(M.ncols()))

    def normalized_columns(self):
        f = ZZ(self.nrows()).factorial()
        nc = [f//sum(v) * v for v in self.columns()]
        for v in nc: v.set_immutable()
        return nc

    def is_product(self):
        r"""
        Return ``None`` or a list of matrices.
        """
        d = self.nrows()
        n = self.ncols()
        M = self._entries

        # make group of rows
        D = DisjointSet(d)
        for c in self.columns():
            i0 = 0
            while i0 < d and not c[i0]:
                i0 += 1
            if i0 == d:
                raise ValueError('zero column')
            i = i0 + 1
            while i < d:
                if c[i]:
                    D.union(i0, i)
                i += 1

        if D.number_of_subsets() == 1:
            # no way to split variables
            return

        ans = []
        for rows in D:
            cols = set().union(*[M.nonzero_positions_in_row(i) for i in rows])
            ans.append(GeneralizedMultipleZetaFunction(M.matrix_from_rows_and_columns(rows, sorted(cols)), as_rows=True))
        return ans
           
    def new_columns(self, forward_only=False):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.generalized_multiple_zeta_values.generalized_multiple_zeta_values import GeneralizedMultipleZetaFunction, linear_forms
            sage: va, vb, vd, vc, ve, vf, vg = linear_forms(3)
            sage: L = GeneralizedMultipleZetaFunction([va, vd, vg])
            sage: [x[0] for x in L.new_columns()] == [vb, vc, ve, vf]
            True

            sage: L = GeneralizedMultipleZetaFunction([va, vg])
            sage: [x[0] for x in L.new_columns()] == [vf]
            True

            sage: L = GeneralizedMultipleZetaFunction([vd,ve,vf])
            sage: for c in L.new_columns(False): print(c)
            ((1, 0, 0), (1/2, 1/2, -1/2))
            ((0, 1, 0), (1/2, -1/2, 1/2))
            ((0, 0, 1), (-1/2, 1/2, 1/2))
            ((1, 1, 1), (1/2, 1/2, 1/2))
            sage: for c in L.new_columns(True): print(c)
            ((1, 1, 1), (1/2, 1/2, 1/2))
        """
        d = self._entries.nrows()
        cols = set(self._entries.columns())
        V = FreeModule(ZZ, d)
        for n in range(1, 2**d):
            c = V(ZZ(n).digits(2,padto=d))
            c.set_immutable()
            if c not in cols:
                try:
                    s = self._entries.solve_right(c)
                except ValueError:
                    continue
                if forward_only and any(not c[i] and x[i] for j,x in enumerate(self.columns()) for i in range(d) if s[j]):
                    continue
                yield c, s

    def is_multizeta(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.generalized_multiple_zeta_values.generalized_multiple_zeta_values import GeneralizedMultipleZetaFunction, linear_forms
            sage: va, vb, vd, vc, ve, vf, vg = linear_forms(3)
            sage: GeneralizedMultipleZetaFunction([va,vd,vg]).is_multizeta()
            [0, 1, 2]
            sage: GeneralizedMultipleZetaFunction([vd,vb,vg]).is_multizeta()
            [1, 0, 2]

            sage: GeneralizedMultipleZetaFunction([vb, vg]).is_multizeta()
            [1, 0, 2]
            sage: GeneralizedMultipleZetaFunction([vf, vg, vc]).is_multizeta()
            [2, 1, 0]
            sage: GeneralizedMultipleZetaFunction([va, vf, vg]).is_multizeta()
            False
            sage: GeneralizedMultipleZetaFunction([vf, vg]).is_multizeta()
            [1, 2, 0]
        """
        M = self._entries
        d = M.nrows()
        columns = M.columns()
        columns.sort(key = lambda x: sum(x))
        variables = [False] * M.nrows()
        order = []
        for c in columns:
            for j in range(d):
                if variables[j]:
                    if not c[j]:
                        return False
                elif c[j]:
                    variables[j] = True
                    order.append(j)
        return order

    def canonicalize(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.generalized_multiple_zeta_values import GeneralizedMultipleZetaFunction
            sage: for cols in [[[1,0,1,0],[0,0,1,0],[1,0,0,0],[0,1,0,1]], [[0,1,0],[1,0,1],[1,1,1]], [[1,0,1,0],[0,0,1,0],[0,1,0,1]], [[1,0,1],[1,1,1],[1,0,0]], [[1,0,1],[0,0,1],[0,1,1]], [[0,1],[1,0],[0,1],[1,1],[0,1]]]:
            ....:     L = GeneralizedMultipleZetaFunction(cols)
            ....:     LL = L.copy()
            ....:     rp, cp = LL.canonicalize()
            ....:     assert all(L[rp(i+1)-1,cp(j+1)-1] == LL[i,j] for i in range(L.nrows()) for j in range(L.ncols()))
        """
        assert not self._entries.is_immutable()
        self._entries, (pr, pc) = self._entries.permutation_normal_form(True)
        nr = self._entries.nrows()
        nc = self._entries.ncols()
#        self._entries = self._entries.parent()([[self._entries[nr-i-1, nc-j-1] for j in range(nc)] for i in range(nr)])
#        Sr = pr.parent(); dr = Sr.degree()
#        Sc = pc.parent(); dc = Sc.degree()
#        pr *= Sr([dr-i for i in range(dr)])
#        pc *= Sc([dc-i for i in range(dc)])
        return (pr, pc)

    def is_convergent(self, e):
        r"""
        Return whether the given exponent is convergent

        EXAMPLES::

            sage: from surface_dynamics.generalized_multiple_zeta_values.generalized_multiple_zeta_values import linear_forms, GeneralizedMultipleZetaFunction
            sage: va, vb, vd, vc, ve, vf, vg = linear_forms(3)
            sage: L = GeneralizedMultipleZetaFunction([va,vb,vc,vd,ve,vf,vg])
            sage: L
            GeneralizedMultipleZetaFunction('(0)(1)(2)(0+1)(0+2)(1+2)(0+1+2)')
            sage: L.is_convergent((0,0,0,0,2,1,1))
            True
        """
        if self.ncols() != len(e):
            raise ValueError("generalized mzv with {} variables but got e={}".format(self.ncols(), e))

        from sage.geometry.polyhedron.constructor import Polyhedron

        m = self.ncols()  # number of linear forms
        d = self.nrows()  # number of variables
        ZZd = FreeModule(ZZ, d)
        zero = ZZd.zero()

        U = set([zero])
        V = set()
        W = set()
        for i,c in enumerate(self.columns()):
            V.clear()
            V.add(zero)
            for j,x in enumerate(c):
                if x:
                    v = e[i] * ZZd.gen(j)
                    v.set_immutable()
                    V.add(v)
            W.clear()
            for u in U:
                for v in V:
                    w = u + v
                    w.set_immutable()
                    W.add(w)
            U, W = W, U
        P = Polyhedron(vertices=U)
        vertices = P.intersection(Polyhedron(rays=[[1]*d])).vertices()
        r = max(max(v.vector()) for v in vertices)
        return r > 1

    def stuffle(self, subset, only_highest_weight=False, canonicalize=False, sort_rows=False):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.generalized_multiple_zeta_values import GeneralizedMultipleZetaFunction

            sage: L = GeneralizedMultipleZetaFunction([[1,0,0],[0,1,0],[0,0,1]])
            sage: L((2,3,4)) == sum(LL((2,3,4)) for LL in L.stuffle([0,1]))
            True
            sage: L((2,3,4)) == sum(LL((2,3,4)) for LL in L.stuffle([1,2]))
            True
            sage: L((2,3,4)) == sum(LL((2,3,4)) for LL in L.stuffle([0,1,2]))
            True
        """
        d = self.ncols()
        n = self.nrows()
        rows = self.rows()
        ans = []
        s = sum(rows[i] for i in subset)
        if any(x > 1 for x in s):
            raise ValueError('invalid subset for stuffle')
        for i in subset:
            new_rows = rows[:]
            new_rows[i] = s
            ans.append(GeneralizedMultipleZetaFunction.from_matrix(matrix(ZZ, new_rows), canonicalize, sort_rows, True, self._reduced))
        if only_highest_weight:
            return ans
        for k in range(2, len(subset) + 1):
            for ss in itertools.combinations(subset, k):
                new_rows = [rows[i] for i in range(n) if i not in ss]
                new_rows.append(s)
                ans.append(GeneralizedMultipleZetaFunction.from_matrix(matrix(ZZ, new_rows), canonicalize, sort_rows, True, self._reduced))
        return ans

    def stuffles(self, only_highest_weight=False, canonicalize=False, sort_rows=False):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.generalized_multiple_zeta_values import GeneralizedMultipleZetaFunction, generalized_multiple_zeta_functions
            sage: L = GeneralizedMultipleZetaFunction([[1,0,0],[0,1,0],[0,0,1]])
            sage: for x in L.stuffles():
            ....:     print(x[0], [vv.lin_prod_string() for vv in x[1]])
            [0, 1] ['(0)(0+1)(2)', '(0+1)(1)(2)', '(1)(1)(0)']
            [0, 1, 2] ['(0)(0+1)(0+2)', '(0+1)(1)(1+2)', '(0+2)(1+2)(2)', '(1)(1)(0+1)', '(1)(0+1)(1)', '(0+1)(1)(1)', '(0)(0)(0)']
            [0, 2] ['(0)(1)(0+2)', '(0+2)(1)(2)', '(1)(0)(1)']
            [1, 2] ['(0)(1)(1+2)', '(0)(1+2)(2)', '(0)(1)(1)']

            sage: for L in generalized_multiple_zeta_functions(3, 3):
            ....:     for subset, stuffle in L.stuffles():
            ....:         assert L((2,3,4)) == sum(LL((2,3,4)) for LL in stuffle), (L, stuffle)
            ....:         assert L((2,4,3)) == sum(LL((2,4,3)) for LL in stuffle), (L, stuffle)
            ....:         assert L((3,2,4)) == sum(LL((3,2,4)) for LL in stuffle), (L, stuffle)
        """
        d = self.ncols()
        rows = self.rows()
        ans = []
        for subset, s in disjoint_vectors(rows, min_size=2, max_size=None):
            yield subset, self.stuffle(subset, only_highest_weight, canonicalize, sort_rows)

    def bad_minimals(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.generalized_multiple_zeta_values.generalized_multiple_zeta_values import GeneralizedMultipleZetaFunction
            sage: L = GeneralizedMultipleZetaFunction(matrix([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]]).transpose())
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        m = self.ncols()  # number of linear forms
        d = self.nrows()  # number of variables
        ZZd = FreeModule(ZZ, d)
        ZZm = FreeModule(ZZ, m)
        ZZmd = FreeModule(ZZ, d + m)
        zero = ZZmd.zero()

        rays = []
        for i,c in enumerate(self.columns()):
            for j,x in enumerate(c):
                if x:
                    rays.append(ZZmd(list(ZZd.gen(j)) + list(ZZm.gen(i))))
        C = Polyhedron(rays=rays)
        Q = Polyhedron(vertices=[[1]*d + [0]*m], rays=[[0]*d + list(b) for b in ZZm.basis()])
        bad_minimals = [v[d:] for v in C.intersection(Q).vertices_list()]
        return bad_minimals


def generalized_multiple_zeta_functions(d, r):
    r"""
    Return the set of generalized multiple zeta functions up to the permutation
    action on variables and linear forms.

    INPUT:

    - ``d`` -- dimension

    - ``r`` -- rank

    EXAMPLES::

        sage: from surface_dynamics.generalized_multiple_zeta_values import generalized_multiple_zeta_functions
        sage: for f in generalized_multiple_zeta_functions(3, 3):
        ....:     print(f)
        GeneralizedMultipleZetaFunction('(0)(1)(2)')
        GeneralizedMultipleZetaFunction('(0+1+2)(0)(1)')
        GeneralizedMultipleZetaFunction('(0+1)(0+2)(1)')
        GeneralizedMultipleZetaFunction('(0+1)(0+2)(0)')
        GeneralizedMultipleZetaFunction('(0+1+2)(0+1)(0)')
        GeneralizedMultipleZetaFunction('(0+1)(0)(2)')
        GeneralizedMultipleZetaFunction('(0+1)(0+2)(1+2)')
        GeneralizedMultipleZetaFunction('(0+1+2)(0+1)(0+2)')
    """
    assert 1 <= r <= d
    V = set()
    vectors = list(linear_forms(d))

    if d == r:
        for _, lins in linearly_independent_vectors(vectors, min_size=d, max_size=d):
            V.add(GeneralizedMultipleZetaFunction.from_matrix(lins, True, False, True, True))
    else:
        for _, lins in linearly_independent_vectors(vectors, min_size=r, max_size=r):
            v = GeneralizedMultipleZetaFunction.from_matrix(lins.transpose(), True, False, True, True)
            if v.has_no_zero_row_or_column():
                V.add(v)
    return V
