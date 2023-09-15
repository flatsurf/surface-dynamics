r"""
Generalized multiple zeta values

TODO:

- turns warning for convergence check (ie off/warning/error)

- introduce linear combinations of *divergent* multiple zeta values
  that we possibly simplify at the end.
"""
from __future__ import absolute_import

import itertools
from collections import defaultdict
import cypari2

from sage.all import ZZ, QQ, matrix, bernoulli_polynomial, prod, FreeModule
# Do not use binomial and factorial from sage.all that are slow and broken
from sage.arith.all import binomial, factorial
from sage.sets.disjoint_set import DisjointSet
from sage.misc.cachefunc import cached_method, cached_function

try:
    from sage.modular.multiple_zeta import Multizetas
except (ImportError, cypari2.handle_error.PariError):
    def Multizetas(*args, **kwds):
        raise ValueError('your sage version does not support multiple zeta values')

VERBOSE = False

# pick a new variable and add it in order to the previous ones
#  1      2      3      4      5      6      7  ...
#  x      y      x+y    z      x+z    y+z    x+y+z  t      x+t   (y+t) (x+y+t) (z+t) ...
# (1000) (0100) (1100) (0010) (1010) (0110) (1110) (0001) (1001)
def linear_forms(d):
    r"""
    EXAMPLES::

        sage: from surface_dynamics.misc.generalized_multiple_zeta_values import linear_forms

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
        v = F(ZZ(n).digits(2,padto=d))
        v.set_immutable()
        yield v

class DivergentZetaError(Exception):
    pass

def handle_term(n, den_tuple):
    r"""
    Entry point for reduction of generalized multiple zeta values into standard
    multiple zeta values.

    EXAMPLES::

        sage: from surface_dynamics.misc.generalized_multiple_zeta_values import handle_term, is_convergent

        sage: M = Multizetas(QQ)

        sage: V1 = FreeModule(ZZ, 1)
        sage: v = V1((1,)); v.set_immutable()
        sage: dt = ((v,3),)
        sage: assert is_convergent(1, dt) and handle_term(1, dt) == M((3,))


        sage: V2 = FreeModule(ZZ, 2)
        sage: va = V2((1,0)); va.set_immutable()
        sage: vb = V2((0,1)); vb.set_immutable()
        sage: vc = V2((1,1)); vc.set_immutable()
        sage: dt = ((va,2), (vc,3))
        sage: assert is_convergent(2, dt) and handle_term(2, dt) == M((2,3))
        sage: dt1 = ((va,2),(vb,3))
        sage: dt2 = ((va,3),(vb,2))
        sage: assert is_convergent(2,dt1) and is_convergent(2,dt2)
        sage: assert handle_term(2, ((va,2), (vb,3))) == handle_term(2, ((va,3), (vb,2))) == M((2,)) * M((3,))

        sage: V3 = FreeModule(ZZ, 3)
        sage: va = V3((1,0,0)); va.set_immutable()
        sage: vb = V3((0,1,0)); vb.set_immutable()
        sage: vc = V3((0,0,1)); vc.set_immutable()
        sage: vd = V3((1,1,0)); vd.set_immutable()
        sage: ve = V3((1,0,1)); ve.set_immutable()
        sage: vf = V3((0,1,1)); vf.set_immutable()
        sage: vg = V3((1,1,1)); vg.set_immutable()
        sage: assert handle_term(3, ((va,2), (vd,3), (vg,4))) == M((2,3,4))
        sage: assert handle_term(3, ((va,2), (vb,3), (vc,4))) == handle_term(3, ((va,3), (vb,2), (vc,4)))  # optional mzv
        sage: assert handle_term(3, ((va,2), (vb,3), (vc,4))) == M((2,)) * M((3,)) * M((4,))
        sage: assert handle_term(3, ((va,1), (vc,2), (vd,3))) == handle_term(3, ((va,1), (vb,2), (ve,3)))
        sage: assert handle_term(3, ((va,1), (vc,2), (vd,3))) == handle_term(3, ((va,2), (vb,1), (vf,3)))
        sage: assert handle_term(3, ((va,1), (vc,2), (vd,3))) ==  M((2,)) * M((1,3))
    """
    if n == 1:
        M = Multizetas(QQ)
        ans = M.zero()
        for v,p in den_tuple:
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

    return _handle_term(den_tuple, n)

@cached_function
def _handle_term(den_tuple, n):
    if VERBOSE:
        print("handle term({}, {})".format(n, den_tuple))
    assert all(len(v) == n for v,p in den_tuple), (n, den_tuple)

    if any(x < 0 or x > 1 for v,p in den_tuple for x in v):
        raise ValueError("unhandled zeta values {}".format(den_tuple))

    # 0. check convergence
    if not is_convergent(n, den_tuple):
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
    a = b = c = 0
    for v,p in den_tuple:
        v = tuple(v)
        if len(v) != 2:
            raise ValueError
        if v == (1,0):
            a = p
        elif v == (0,1):
            b = p
        elif v == (1,1):
            c = p
        else:
            return
    return a,b,c

def to_Z3(den_tuple, sort=True):
    a = b = c = d = e = f = g = 0
    for v,p in den_tuple:
        v = tuple(v)
        if len(v) != 3:
            raise ValueError
        if v == (1,0,0):
            a = p
        elif v == (0,1,0):
            b = p
        elif v == (0,0,1):
            c = p
        elif v == (1,1,0):
            d = p
        elif v == (1,0,1):
            e = p
        elif v == (0,1,1):
            f = p
        elif v == (1,1,1):
            g = p
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

def is_convergent(n, den):
    r"""
    TESTS::

        sage: import itertools
        sage: from surface_dynamics.misc.generalized_multiple_zeta_values import is_convergent
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
        sage: for p in itertools.product([0,1,2], repeat=7):
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
        sage: print(N)
        125
    """
    from sage.geometry.polyhedron.constructor import Polyhedron
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

    assert all(len(v) == n for v,p in den), (n, den)

    # TODO: fast code path

    R = PolynomialRing(QQ, 'x', n)
    x = R.gens()
    den_poly = prod(linear_form(R, v)**p for v,p in den)
    newton_polytope = Polyhedron(vertices=den_poly.exponents(), rays=negative_rays(n))
    V = newton_polytope.intersection(Polyhedron(rays=[[1]*n])).vertices()
    r = max(max(v.vector()) for v in V)
    return r > 1

def convergent_multizeta(t):
    r"""
    Multizeta value at a convergent index ``t``.

    TESTS::

        sage: from surface_dynamics.misc.generalized_multiple_zeta_values import convergent_multizeta
        sage: assert all(convergent_multizeta(t) == Multizeta(*t) for t in [(2,),(3,),(1,2),(3,2),(1,1,2)])

        sage: convergent_multizeta((0,3))
        ζ(2) - ζ(3)
        sage: convergent_multizeta((0,2,2))
        ζ(1,2) - ζ(2,2)
        sage: convergent_multizeta((1,0,3))
        ζ(1,2) - ζ(1,3) - ζ(2) + ζ(3)
        sage: convergent_multizeta((0,1,3))
        -ζ(1,3) + ζ(2) - ζ(3)

        sage: convergent_multizeta((0, 4))
        ζ(3) - ζ(4)
        sage: convergent_multizeta((-1, 5))
        1/2*ζ(3) - 1/2*ζ(4)
        sage: convergent_multizeta((-2, 5))
        1/3*ζ(2) - 1/2*ζ(3) + 1/6*ζ(4)

        sage: convergent_multizeta((-1, 3, 4))
        1/2*ζ(1,4) - 1/2*ζ(2,4)
        sage: convergent_multizeta((-1, -1, 8))
        1/8*ζ(4) - 5/12*ζ(5) + 3/8*ζ(6) - 1/12*ζ(7)
        sage: convergent_multizeta((-2,-2,10))
        1/18*ζ(4) - 4/15*ζ(5) + 31/72*ζ(6) - 1/4*ζ(7) + 1/72*ζ(8) + 1/60*ζ(9)
        sage: convergent_multizeta((-1, -2, 10))
        1/10*ζ(5) - 3/8*ζ(6) + 5/12*ζ(7) - 1/8*ζ(8) - 1/60*ζ(9)

        sage: convergent_multizeta((4,-4,10))
        -1/3*ζ(1,10) + 1/30*ζ(3,10) + 1/5*ζ(4,5) - 1/2*ζ(4,6) + 1/3*ζ(4,7) - 1/30*ζ(4,9) - 1/10*ζ(8) - 2/5*ζ(9) + 1/2*ζ(10)
        sage: convergent_multizeta((2,-1,8))
        -1/2*ζ(1,8) + 1/2*ζ(2,6) - 1/2*ζ(2,7) - 1/2*ζ(7) + 1/2*ζ(8)

        sage: convergent_multizeta((0,3,2))
        ζ(2,2) - ζ(3,2)

    Divergent cases::

        sage: convergent_multizeta((0,2))
        Traceback (most recent call last):
        ...
        DivergentZetaError: divergent multizeta value (0, 2)
        sage: convergent_multizeta((0,0,3))
        Traceback (most recent call last):
        ...
        DivergentZetaError: divergent multizeta value (0, 0, 3)
        sage: convergent_multizeta((1,0,2))
        Traceback (most recent call last):
        ...
        DivergentZetaError: divergent multizeta value (1, 0, 2)
        sage: convergent_multizeta((0,1,2))
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
        return Multizetas(QQ)(t)
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

        sage: from surface_dynamics.misc.generalized_multiple_zeta_values import linear_forms, is_multizeta
        sage: va, vb, vd, vc, ve, vf, vg = linear_forms(3)

        sage: is_multizeta(3, ((vb, 2), (vf, 5), (vg, 2)))
        ζ(2,5,2)

        sage: is_multizeta(3, [(vg, 5)])
        1/2*ζ(3) - 3/2*ζ(4) + ζ(5)

        sage: assert is_multizeta(3, ((va,3), (vb,3), (vc,3))) is None
        sage: assert is_multizeta(3, ((vb,2), (ve,5), (vg,1))) is None
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
    INPUT:

    - ``n`` - number of variables

    - ``den_tuple`` - tuple of pairs ``(vector, power)``

    TESTS::

        sage: from surface_dynamics.misc.generalized_multiple_zeta_values import is_product
        sage: is_product(3, [((1,0,0),2), ((0,1,0),5), ((1,1,0),1), ((0,0,1),3)])
        [(2, (((0, 1), 5), ((1, 0), 2), ((1, 1), 1))), (1, (((1), 3),))]
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
    Make calls to `apply_relation` until we get rid of term

    EXAMPLES::

        sage: from surface_dynamics.misc.generalized_multiple_zeta_values import kill_relation

        sage: V = FreeModule(ZZ, 2)
        sage: va = V((1,0)); va.set_immutable()
        sage: vb = V((0,1)); vb.set_immutable()
        sage: vc = V((1,1)); vc.set_immutable()

        sage: den_tuple = ((va,2),(vb,3),(vc,4))
        sage: kill_relation(2, den_tuple, 2, [1,1,0])
        [(1, (((0, 1), 3), ((1, 1), 6))),
         (2, (((0, 1), 2), ((1, 1), 7))),
         (1, (((1, 0), 2), ((1, 1), 7))),
         (3, (((0, 1), 1), ((1, 1), 8))),
         (3, (((1, 0), 1), ((1, 1), 8)))]
    """
    assert len(relation) == len(den_tuple)
    assert 0 <= i < len(den_tuple)
    assert sum(relation[j] * den_tuple[j][0] for j in range(len(den_tuple))) == den_tuple[i][0]
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

def has_term_sum_of_smaller_terms(n, den_tuple):
    r"""
    Look for a vector v_i and {v_j} with sum(v_j) = v_i

    EXAMPLES::

        sage: from surface_dynamics.misc.generalized_multiple_zeta_values import linear_forms, has_term_sum_of_smaller_terms

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
    If (x1+x2+...+xd) is not present, use a linear relation to create it. Then
    try to kill using other forms. Then try to write as P(x1,x2,...,x_{d-1}) (x1+x2+...+xd) and recurse.

    Should solve all Tonrheim

    EXAMPLES::

        sage: from surface_dynamics.misc.generalized_multiple_zeta_values import linear_forms, is_reducible
        sage: va, vb, vd, vc, ve, vf, vg = linear_forms(3)
        sage: is_reducible(3, ((va,3),(vb,3),(vc,3)))
        [(1, (((0, 0, 1), 3), ((0, 1, 1), 3), ((1, 1, 1), 3))),
         (1, (((0, 1, 0), 3), ((0, 1, 1), 3), ((1, 1, 1), 3))),
         (3, (((0, 0, 1), 2), ((0, 1, 1), 4), ((1, 1, 1), 3))),
         (3, (((0, 1, 0), 2), ((0, 1, 1), 4), ((1, 1, 1), 3))),
        ...
         (90, (((1, 0, 0), 1), ((1, 0, 1), 1), ((1, 1, 1), 7))),
         (90, (((0, 1, 0), 1), ((1, 1, 0), 1), ((1, 1, 1), 7))),
         (90, (((1, 0, 0), 1), ((1, 1, 0), 1), ((1, 1, 1), 7)))]
    """
    F = FreeModule(ZZ, n)
    vmax = F([1] * n)
    vmax.set_immutable()

    if len(den_tuple) == 1:
        return

    # force the max vector (1,1,...,1) to appear
    imax = None
    for i,(v,p) in enumerate(den_tuple):
        if v == vmax:
            imax = i
            break
    if imax is None:
        imax = len(den_tuple)
        den_tuple = den_tuple + ((vmax,0),)
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


def Z3_sort_abc(a,b,c,d,e,f,g):
    if b > a:
        if c > b:
            # c > b > a
            a,b,c,d,e,f,g = c,b,a,f,e,d,g
        elif c > a:
            # b >= c > a
            a,b,c,d,e,f,g = b,c,a,f,d,e,g
        else:
            # b > a >= c
            a,b,c,d,e,f,g = b,a,c,d,f,e,g
    elif c > a:
        # c > a >= b
        a,b,c,d,e,f,g = c,a,b,e,f,d,g
    elif c > b:
        # a >= c >= b
        a,b,c,d,e,f,g = a,c,b,e,d,f,g
    else:
        pass
    assert a >= b >= c
    return a,b,c,d,e,f,g

def is_Z2_convergent(a,b,c):
    r"""
    TESTS::

        sage: from surface_dynamics.misc.generalized_multiple_zeta_values import is_Z2_convergent

    Convergent examples::

        sage: assert is_Z2_convergent(1,1,1)
        sage: assert is_Z2_convergent(2,2,0)
        sage: assert is_Z2_convergent(0,0,3)

    Divergent examples::

        sage: assert not is_Z2_convergent(0,0,2)
        sage: assert not is_Z2_convergent(1,2,0)
    """
    from sage.geometry.polyhedron.constructor import Polyhedron
    x, y = ZZ['x,y'].gens()
    poly = x**a * y**b * (x+y)**c
    newton_polytope = Polyhedron(vertices=poly.exponents(), rays=[(-1,0),(0,-1)])
    V = newton_polytope.intersection(Polyhedron(rays=[(1,1)])).vertices()
    r = max(max(v.vector()) for v in V)
    return r > 1

def is_Z3_convergent(a,b,c,d,e,f,g):
    r"""
    TESTS::

        sage: from surface_dynamics.misc.generalized_multiple_zeta_values import is_Z3_convergent

    Convergent examples::

        sage: assert is_Z3_convergent(2,0,2,2,0,0,0)
        sage: assert is_Z3_convergent(0,0,0,0,0,0,4)
        sage: assert is_Z3_convergent(1,0,0,1,0,0,2)
        sage: assert is_Z3_convergent(0,1,0,1,0,0,2)
        sage: assert is_Z3_convergent(0,1,0,0,0,1,2)

    Divergent examples::

        sage: assert not is_Z3_convergent(0,0,0,1,1,1,0)
        sage: assert not is_Z3_convergent(0,0,0,0,0,0,3)

    """
    from sage.geometry.polyhedron.constructor import Polyhedron
    x, y, z = ZZ['x,y,z'].gens()
    poly = x**a * y**b * z**c * (x+y)**d * (x+z)**e * (y+z)**f * (x+y+z)**g
    newton_polytope = Polyhedron(vertices=poly.exponents(), rays=[(-1,0,0),(0,-1,0),(0,0,-1)])
    V = newton_polytope.intersection(Polyhedron(rays=[(1,1,1)])).vertices()
    r = max(max(v.vector()) for v in V)
    return r > 1

def Z2(a, b, c, check_convergence=True):
    M = Multizetas(QQ)
    if a == 0 and b == 0:
        return M((c-1,)) - M((c,))
    else:
        return sum(binomial(a+i-1,i) * M((b-i,c+a+i)) for i in range(b)) + sum(binomial(b+i-1,i) * M((a-i,c+b+i)) for i in range(a))

def Z3(a, b, c, d, e, f, g, check_convergence=True):
    r"""
    The function ``Z3(a,b,c,d,e,f,g)``.

    .. MATH::

        \sum_{x,y,z \geq 1} x^{-a} y^{-b} z^{-c} (x+y)^{-d} (x+z)^{-e} (y+z)^{-f} (x+y+z)^{-g}

    The reduction algorithm was designed by Bill Allombert.

    EXAMPLES::

        sage: from surface_dynamics.misc.generalized_multiple_zeta_values import Z3

        sage: M = Multizetas(QQ)

        sage: Z3(1,1,1,1,1,1,1)
        21/2*ζ(1,1,5) + 9/2*ζ(1,2,4) - 3/2*ζ(1,3,3) - 3/2*ζ(1,4,2) + 9/2*ζ(1,6)
        sage: Z3(3,0,0,0,0,3,0)
        6*ζ(1,4) - 12*ζ(1,5) + 3*ζ(2,3) - 6*ζ(2,4) + ζ(3,2) - 2*ζ(3,3)

        sage: assert Z3(2,3,4,0,0,0,0) == M((2,)) * M((3,)) * M((4,))
        sage: assert Z3(1,0,0,2,0,0,3) == M((1,2,3))

        sage: assert Z3(0,0,0,2,0,1,1) == M((4,)) / 2
        sage: assert Z3(1,0,1,1,0,0,1) == 3 * M((1,1,2))

        sage: assert Z3(0,0,0,0,0,0,4) == 1/2 * M((2,)) - 3/2 * M((3,)) + M((4,))

        sage: assert Z3(0,0,0,2,0,1,1) == 2 * M((1,1,2)) - M((2,2)) - 3 *M((1,3))
    """
    M = Multizetas(QQ)
    CHECK_CONVERGENCE = False

    # x^a y^b z^c (x+y)^d (x+z)^e (y+z)^f (x+y+z)^g
    if VERBOSE:
        print("Z3({},{},{},{},{},{},{})".format(a,b,c,d,e,f,g))
    if a < 0 or b < 0 or c < 0 or d < 0 or e < 0 or f < 0 or g < 0:
        raise ValueError("invalid exponents for Z3: a={} b={} c={} d={} e={} f={} g={}".format(a,b,c,d,e,f,g))

    if check_convergence and not is_Z3_convergent(a,b,c,d,e,f,g):
        raise DivergentZetaError("divergent Z3({},{},{},{},{},{},{})".format(a,b,c,d,e,f,g))

    # step 1: try to get rid of the terms (x+y), (x+z), (y+z)
    if d and e and f:
        if VERBOSE:
            print("reduction (x+y+z) = ((x+y) + (x+z) + (y+z)) / 2")
        return (Z3(a,b,c,d-1,e,f,g+1,CHECK_CONVERGENCE) + Z3(a,b,c,d,e-1,f,g+1,CHECK_CONVERGENCE) + Z3(a,b,c,d,e,f-1,g+1,CHECK_CONVERGENCE)) / 2
    if a and f:
        # x^a (y+z)^f -> (x+y+z)^g
        # additive version Z3(a-1,b,c,d,e,f,g+1) + Z3(a,b,c,d,e,f-1,g+1)
        if VERBOSE:
            print("reduction (x+y+z) = (x) + (y+z)")
        return sum(binomial(a+k-1, k) * Z3(0,b,c,d,e,f-k,g+a+k,CHECK_CONVERGENCE) for k in range(f)) + \
               sum(binomial(f+k-1, k) * Z3(a-k,b,c,d,e,0,g+f+k,CHECK_CONVERGENCE) for k in range(a))
    if b and e:
        # y^b (x+z)^e -> (x+y+z)^g
        # additive version Z3(a,b-1,c,d,e,f,g+1) + Z3(a,b,c,d,e-1,f,g+1)
        if VERBOSE:
            print("reduction (x+y+z) = (y) + (x+z)")
        return sum(binomial(b+k-1, k) * Z3(a,0,c,d,e-k,f,g+b+k,CHECK_CONVERGENCE) for k in range(e)) + \
               sum(binomial(e+k-1, k) * Z3(a,b-k,c,d,0,f,g+e+k,CHECK_CONVERGENCE) for k in range(b))
    if c and d:
        # z^c (x+y)^d -> (x+y+z)^g
        # additive version Z3(a,b,c-1,d,e,f,g+1) + Z3(a,b,c,d-1,e,f,g+1)
        if VERBOSE:
            print("reduction (x+y+z) = (z) + (x+y)")
        return sum(binomial(c+k-1, k) * Z3(a,b,0,d-k,e,f,g+c+k,CHECK_CONVERGENCE) for k in range(d)) + \
               sum(binomial(d+k-1, k) * Z3(a,b,c-k,0,e,f,g+d+k,CHECK_CONVERGENCE) for k in range(c))

    assert d*e*f == a*f == b*e == c*d == 0

    a,b,c,d,e,f,g = Z3_sort_abc(a,b,c,d,e,f,g)

    # step 2: kill c
    if c:
        if VERBOSE:
            print("reduction (x+y+z) = (x) + (y) + (z)")
        return Z3(a-1,b,c,d,e,f,g+1,CHECK_CONVERGENCE) + Z3(a,b-1,c,d,e,f,g+1,CHECK_CONVERGENCE) + Z3(a,b,c-1,d,e,f,g+1,CHECK_CONVERGENCE)

    assert c == 0

    if b:
        # x^a y^b -> (x+y)^d
        # additive version Z3(a-1,b,0,d+1,e,f,g) + Z3(a,b-1,0,d+1,e,f,g)
        if VERBOSE:
            print("reduction (x+y) = (x) + (y)")
        return sum(binomial(a+k-1, k) * Z3(0,b-k,c,d+a+k,e,f,g,CHECK_CONVERGENCE) for k in range(b)) + \
               sum(binomial(b+k-1, k) * Z3(a-k,0,c,d+b+k,e,f,g,CHECK_CONVERGENCE) for k in range(a))

    assert b == c == 0

    assert b == c == d*e*f == 0

    if a and f:
        raise RuntimeError

    assert b == c == d*e*f == a*f == 0

    if a == 0:
        d,e,f = sorted([d,e,f], reverse=True)

    if a and d and e:
        if VERBOSE:
            print("reduction (x+y+z) = (x+y) + (x+z) - (x)")
        return Z3(a,b,c,d-1,e,f,g+1,CHECK_CONVERGENCE) + Z3(a,b,c,d,e-1,f,g+1,CHECK_CONVERGENCE) - Z3(a-1,b,c,d,e,f,g+1,CHECK_CONVERGENCE)

    assert b == c == d*e*f == a*f == a*d*e == 0

    # x^a (x+y)^d (x+y+z)^g
    if e == f == 0:
        return convergent_multizeta((a,d,g))
    # x^a (x+z)^e (x+y+z)^g
    if d == f == 0:
        return convergent_multizeta((a,e,g))
    if a == f == 0:
        if g == 0:
            # (x+y)^d (x+z)^e
            assert d and e
            return M((d-1,e)) - M((d,e)) + M((e-1,d)) - M((e,d)) + M((d+e-1,)) - M((d+e,))
        elif d == e == 1:
            if VERBOSE:
                print("Bill formulas")
            return sum(M((1,g+1-i,i)) for i in range(2,g+1)) - 3*M((1,g+1))
        elif g == 1:
            if VERBOSE:
                print("[0,0,0,d,0,f,1] = [f,0,d-1,1,0,0,1] - [f,d,1] - zeta([d,1,f]) - zeta([1+d,f])")
            f,d = sorted([d,e])
            return Z3(f,0,d-1,1,0,0,1,CHECK_CONVERGENCE) - Z2(f,d,1,CHECK_CONVERGENCE) - M((f,1,d)) - M((f,d+1))
        else:
            # THIS LOOKS WRONG
            # [0,0,0,d,0,f,g] = [f,0,d,g,0,0,0] - [f,d,0,0,0,0,g] - [f,d,g] - zeta([d,g,f]) - zeta([g+d,f])
            if VERBOSE:
                print("[0,0,0,0,e,f,g] = [e,0,f,g,0,0,0] - [e,0,f,0,0,0,g] - Z(e,g,f) - Z(e,f+g) - Z(e+f,g)")
            f,d = sorted([d,e])
            return Z3(f,0,d,g,0,0,0,CHECK_CONVERGENCE) - Z3(f,d,0,0,0,0,g,CHECK_CONVERGENCE) - Z2(f,d,g,CHECK_CONVERGENCE) - M((f,g,d)) - M((f,g+d))

    raise RuntimeError
