r"""
Automatic reduction of generalized multiple zeta values


EXAMPLES::

    sage: from surface_dynamics.generalized_multiple_zeta_values.gmzv_generic_reduction import GMZVSolver
    sage: S = GMZVSolver(3, 3)
    sage: S
    Solve strategy d=3 r=3 (0 known and 8 unknown)
    sage: S.update_multizetas()
    sage: S
    Solve strategy d=3 r=3 (1 known and 7 unknown)
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
import collections

from sage.misc.misc_c import prod
from sage.misc.prandom import random, choice
from sage.misc.cachefunc import cached_method

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.modules.free_module_element import vector
from sage.modules.free_module import VectorSpace, FreeModule
from sage.matrix.constructor import matrix 

from surface_dynamics.misc.linalg import linearly_independent_vectors

from .gmzv_two_variables import Z2
from .gmzv_three_variables import Z3
from .generalized_multiple_zeta_values import linear_forms, convergent_multizeta, DivergentZetaError, GeneralizedMultipleZetaFunction, to_Z2, to_Z3, generalized_multiple_zeta_functions
from .options import DIVERGENT_MZV


class EvalAlgebraic:
    r"""
    """
    def __init__(self, S, v, lin_comb):
        assert v.ncols() == len(lin_comb)
        self.v = v
        self.lin_comb = lin_comb
        self.support = [i for i in range(v.ncols()) if lin_comb[i]]
        columns = v.columns()
        c = sum(x * col for x, col in zip(lin_comb, columns))
        assert all(x == 0 or x == 1 for x in c), c
        self.c = c = vector(ZZ, c)
        self.S = S
        self.new_forms = {}
        for i in self.support:
            new_columns = columns[:]
            new_columns.pop(i)
            new_columns.append(c)
            vv = GeneralizedMultipleZetaFunction(new_columns)
            _, cp = vv.canonicalize()
            vv.set_immutable()
            assert S.strategy[vv] is not None
            self.new_forms[i] = (vv, cp)

    def children(self):
        r"""
        Return the list of zeta used
        """
        return [nf[0] for nf in self.new_forms.values()]

    def __repr__(self):
        return 'algebraic {} -> {}'.format(self.lin_comb, self.c)

    def __call__(self, exponents):
        # TODO: implement multiplicative version
        e = tuple(exponents) + (0,)
        D = {e: ZZ.one()}
        ans = 0
        todo = [e]
        while todo:
            e = todo.pop(0)
            if all(e[j] for j in self.support):
                coeff = D.pop(e)
                if not coeff:
                    continue
                for j in self.support:
                    ee = list(e)
                    ee[j] -= 1
                    ee[-1] += 1
                    ee = tuple(ee)
                    if ee not in D:
                        D[ee] = self.lin_comb[j] * coeff
                        todo.append(ee)
                    else:
                        assert ee in todo
                        D[ee] += self.lin_comb[j] * coeff

        for e, coeff in D.items():
            for j in self.support:
                if e[j] == 0:
                    break
            assert e[j] == 0
            e = e[:j] + e[j+1:]
            vv, cp = self.new_forms[j]
            e = [e[cp(i+1)-1] for i in range(self.v.ncols())]
            ans += coeff * self.S.eval(vv, e)

        return ans


class FunctionWrap:
    def __init__(self, f, groups, exp_position=0, extra_args=None):
        self.f = f
        self.groups = tuple(map(tuple, groups))
        if extra_args is None:
            extra_args=()
        self.extra_args = extra_args
        self.pos = exp_position
    def __hash__(self):
        hash((self.f, self.groups, self.extra_args))
    def __repr__(self):
        gstring = '(' + ','.join('+'.join(map(str, g)) for g in self.groups) + ')'
        frepr = None
        if isinstance(self.f, GeneralizedMultipleZetaFunction):
            dat = self.f.is_multizeta()
            frepr = self.f.lin_prod_string()
        else:
            try:
                is_eval = self.f.__func__ == GMZVSolver.eval
            except AttributeError:
                is_eval = False
            else:
                if is_eval:
                    frepr = self.extra_args[0].lin_prod_string()
        if frepr is None:
            frepr = repr(self.f)

        return frepr + gstring

    def __call__(self, e):
        ee = [sum(e[i] for i in g) for g in self.groups]
        args = self.extra_args[:self.pos] + (ee,) + self.extra_args[self.pos:]
        return self.f(*args)


def stuffle_to_lin_comb(v, rp, cp, stuffle):
    r"""
    Return a linear combination of functions taking list of functions

    INPUT:

    - ``v``

    - ``rp`` -- row permutation

    - ``cp`` -- column permutation

    - ``stuffle`` --

    EXAMPLES::

        sage: from surface_dynamics.generalized_multiple_zeta_values.generalized_multiple_zeta_values import GeneralizedMultipleZetaFunction, linear_forms
        sage: from surface_dynamics.generalized_multiple_zeta_values.gmzv_generic_reduction import stuffle_to_lin_comb
        sage: va, vb, vc, vd, ve, vf, vg = linear_forms(3)
        sage: G = GeneralizedMultipleZetaFunction([vd, va, vb])
        sage: S = SymmetricGroup(3)
        sage: p = S.one()
    """
    vtmp = v.copy()
    _, _ = vtmp.canonicalize()
    assert all(v._entries[rp(i+1)-1, cp(j+1)-1] == vtmp._entries[i][j] for i in range(v.nrows()) for j in range(v.ncols()))

    if isinstance(stuffle, (tuple, list)):
        # the subset is given with respect to the canonical form
        subset, subtasks = stuffle
        subset = [rp(i + 1) - 1 for i in subset]
        stuffle = v.stuffle(subset, only_highest_weight=False)

        main = stuffle[:len(subset)]
        assert len(main) == len(subtasks)
        for vv, (rrp, ccp, task) in zip(main, subtasks):
            assert vv.nrows() == v.nrows()
            assert vv.ncols() == v.ncols()
            for (vvv, rrrp, cccp) in stuffle_to_lin_comb(vv, rrp*rp, ccp*cp, task):
                yield (vvv, rrrp , cccp)
        remainder = stuffle[len(subset):]
        for vv in remainder:
            assert vv.nrows() < v.nrows()
            assert vv.ncols() == v.ncols()
            yield (vv, None, rp.parent().one())
    else:
        # check that we got the correct element
        assert stuffle == vtmp
        yield (stuffle, rp, cp)


class LinearCombination:
    def __init__(self, comb):
        self.comb = [(coeff, f) for coeff, f in comb if coeff]

    def children(self):
        ans = []
        for coeff, f in self.comb:
            if isinstance(f, FunctionWrap):
                try:
                    is_eval = isinstance(f.f.__self__, GMZVSolver)
                except AttributeError:
                    is_eval = False
                if is_eval:
                    ans.append(f.extra_args[0])
                else:
                    ans.append(f.f)
            else:
                ans.append(f)
        return ans

    def __repr__(self):
        return ' + '.join(str(f) if coeff == 1 else '%s %s' % (coeff, f) for coeff, f in self.comb)

    def __call__(self, e):
        return sum(coeff * f(e) for coeff, f in self.comb)


class GMZVSolver:
    r"""
    A solver for finding relations between generalized multiple zeta functions and values.

    At any stage, the solver has some "known" gmzv that can be reduced to linear combinations
    of standard mzv and some "unknown" forms. One can then try different options to look for
    new relations and getting more "known" gmzv.

    Attributes:

    - ``d`` -- dimension

    - ``r`` -- rank

    - ``V`` -- list of multiple zeta functions to be handled

    - ``strategy`` -- dictionary

    - ``what`` -- dictionary

    - ``level`` -- dictionary

    EXAMPLES::

        sage: from surface_dynamics.generalized_multiple_zeta_values.gmzv_generic_reduction import GMZVSolver
        sage: S = GMZVSolver(3, 3)
        sage: S
        Solve strategy d=3 r=3 (0 known and 8 unknown)
    """
    def __init__(self, d, r):
        self.d = d
        self.r = r
        self.V = list(generalized_multiple_zeta_functions(d, r))

        # a list v -> function to evaluate at v
        self.strategy = {v: None for v in self.V}
        self.what = {v: None for v in self.V}
        self.level = {v: None for v in self.V}

    def graphviz_string(self, filename=None):
        if filename is None:
            import sys
            f = sys.stdout
        else:
            f = open(filename, 'w')

        to_index = {v:i for i,v in enumerate(self.V)}
        f.write("digraph G {\n")
        f.write("  rankdir=TB\n")
        for u, s in self.strategy.items():
            i = to_index[u]
            if isinstance(s[0], EvalAlgebraic):
                f.write("  %d [style=filled fillcolor=blue]\n" % i)
            elif isinstance(s[0], LinearCombination):
                f.write("  %d [style=filled fillcolor=red]\n" % i)
            elif isinstance(s[0], GeneralizedMultipleZetaFunction):
                f.write("  %d [style=filled fillcolor=yellow]\n" % i)
            for v in s[0].children():
                if v not in to_index:
                    continue
                f.write("%d -> %d\n" % (i, to_index[v]))

        f.write("}\n")

        if filename is not None:
            f.close()

    def try_small_exponents(self, e_bound=3, s_bound=7, vlist=None):
        if vlist is None:
            vlist = sorted(self.V, key=lambda v: self.level[v])
        for v in vlist:
            for i in IntegerListsLex(length=self.r, min_part=1, max_part=e_bound, max_sum=s_bound):
                if v.is_convergent(i):
                    try:
                        self.eval(v, tuple(i))
                    except DivergentZetaError:
                        print("failure at level {} with\n{}\nand exponent {}".format(self.level[v], v, i))
                    except ValueError:
                        print("called for divergent zeta at level {} with\n{}\nand exponent {}".format(self.level[v], v, i))

    def eval_abstract_sum(self, s):
        return sum(QQ(coeff) * self.eval_den_tuple(den_tuple) for den_tuple, coeff in s._data.items())

    def eval_den_tuple(self, den_tuple):
        columns = []
        exponents = []
        for col, exp in den_tuple._dict.items():
            columns.append(col)
            exponents.append(exp)
        v = GeneralizedMultipleZetaFunction(columns)
        rp, cp = v.canonicalize()
        v.set_immutable()
        exponents = [exponents[cp(i+1)-1] for i in range(len(columns))]
        return self.eval(v, exponents)

    def eval(self, v, exponents):
        if len(exponents) != v.ncols():
            raise ValueError("generalized multiple zeta values with {} linear forms but got e={}".format(v.ncols(), exponents))
        if v not in self.strategy:
            raise NotImplementedError('v={} not in strategy'.format(v.lin_prod_string()))
        if self.strategy[v] is None:
            raise ValueError('no available strategy')
        f, extra_args = self.strategy[v]
        return f(exponents, *extra_args)

    def __repr__(self):
        known = sum(x is not None for x in self.strategy.values())
        unknown = sum(x is None for x in self.strategy.values())
        return 'Solve strategy d={} r={} ({} known and {} unknown)'.format(self.d, self.r, known, unknown)

    def known(self):
        return [v for v,status in self.strategy.items() if status is not None]

    def unknown(self):
        return [v for v,status in self.strategy.items() if status is None]

    def update_multizetas(self):
        r"""
        Set trivial evaluation for multizetas and set their levels to 0
        """
        full_mzv = [[1] * k + [0] * (self.d - k) for k in range(self.d-1, 0, -1)]
        for cols in itertools.combinations(full_mzv, self.r-1):
            v = GeneralizedMultipleZetaFunction(([1]*self.d,) + cols)
            v.set_immutable()
            assert v in self.strategy, v
            self.strategy[v] = (v, ())
            self.what[v] = ('multizeta', None)
            self.level[v] = 0

    def update_algebraic(self, v, lin_comb):
        # take lin_comb of the current columns to generate a new one
        f = EvalAlgebraic(self, v, lin_comb)
        self.strategy[v] = (f, ())
        self.what[v] = ('algebraic', lin_comb)
        self.level[v] = max(self.level[w] for w in f.children()) + 1

    def algebraic(self, forward_only=False):
        r"""
        Iterate through available algebraic relations.

        EXAMPLES::

            sage: from surface_dynamics.generalized_multiple_zeta_values.gmzv_generic_reduction import GMZVSolver
            sage: S = GMZVSolver(3, 3)
            sage: S.update_multizetas()
            sage: S
            Solve strategy d=3 r=3 (1 known and 7 unknown)
            sage: algebraic = list(S.algebraic())
            sage: print([ans[0].lin_prod_string() for ans in algebraic])
            ['(0+1+2)(0)(1)']

        One can make the solver aware of this relation via::

            sage: S.update_algebraic(*algebraic[0])
            sage: S
            Solve strategy d=3 r=3 (2 known and 6 unknown)
        """
        for v in self.strategy:
            assert v.ncols() == self.r
            if self.strategy[v] is not None:
                continue
            for c, lin_comb in v.new_columns(forward_only):
                unknown = False
                for i in range(self.r):
                    if not lin_comb[i]:
                        continue
                    m = v._entries.__copy__()
                    m[:,i] = c
                    vv = GeneralizedMultipleZetaFunction.from_matrix(m, True, False, True, True)
                    if self.strategy[vv] is None:
                        unknown = True
                        break
                if unknown:
                    continue
                yield v, lin_comb

    def update_forward_stuffle(self, v, stuffle):
        d = v.nrows()  # number of variables
        n = v.ncols()

        # stuffle is a tree: each element is either a LFF or a list
        # subset to stuffle, list of triples (rperm, cperm, new node)
        comb = []
        Sn = SymmetricGroup(n)
        level = 0
        for vv, rp, cp in stuffle_to_lin_comb(v, Sn.one(), Sn.one(), stuffle):
            groups = [(cp(i+1) - 1,) for i in range(n)]
            if vv.nrows() < d:
                # assume .eval is working for less rows
                comb.append((1, FunctionWrap(vv, groups)))
            else:
                level = max(level, self.level[vv])
                comb.append((1, FunctionWrap(self.eval, groups, exp_position=1, extra_args=(vv,))))

        self.strategy[v] = (LinearCombination(comb), ())
        self.what[v] = ('forward stuffle', stuffle)
        self.level[v] = level + 1

    def _find_iterated_stuffle(self, v):
        if self.strategy[v] is not None:
            return v
        for subset, stuffle in v.stuffles(only_highest_weight=True):
            has = True
            for i, x in enumerate(stuffle):
                y = x.copy()
                rperm, cperm = y.canonicalize()
                y.set_immutable()
                ans = self._find_iterated_stuffle(y)
                if not ans:
                    has = False
                    break
                stuffle[i] = (rperm, cperm, ans)
            if has:
                return [subset, stuffle]
        return False

    def iterated_forward_stuffles(self):
        r"""
        Return the unknown gmzv that could be reduced via iteration of forward
        stuffles to known gmzv.

        EXAMPLES::

            sage: from surface_dynamics.generalized_multiple_zeta_values.gmzv_generic_reduction import GMZVSolver
            sage: S = GMZVSolver(3, 3)
            sage: S.update_multizetas()
            sage: stuffles = list(S.iterated_forward_stuffles())
            sage: list(ans[0].lin_prod_string() for ans in stuffles)
            ['(0)(1)(2)', '(0+1)(0+2)(0)', '(0+1)(0)(2)']
            sage: S.update_forward_stuffle(*stuffles[1])
            sage: S
            Solve strategy d=3 r=3 (2 known and 6 unknown)
        """
        for v in self.strategy:
            if self.strategy[v] is not None:
                continue
            stuffle = self._find_iterated_stuffle(v)
            if stuffle:
                yield v, stuffle

    def update_multistuffle(self, v, vnsym, lin_comb):
        r"""
        INPUT:

        - ``v``

        - ``vnsym``

        - ``lin_comb`` - list of ``(typ, coeff, vector)``
        """
        D = {}  # our linear combination of full data
        rem = {} # remainders
        level = 0
        for t, coeff, data in lin_comb:
            if t == 0:
                # by design, the highest length values compensate and we ignore them
                vv, subset = data
                for vvv in vv.stuffle(subset, only_highest_weight=False, sort_rows=True):
                    if vvv.nrows() == v.nrows():
                        continue
                    if vvv in D:
                        D[vvv] += coeff
                    else:
                        D[vvv] = coeff
            elif t == 1:
                # known value
                vv = data
                if vv in D:
                    D[vv] += coeff
                else:
                    D[vv] = coeff

        for vv in list(D):
            if D[vv] == 0:
                del D[vv]

        comb = []
        for vv, coeff in D.items():
            vv = vv.copy()
            rp, cp = vv.canonicalize()
            vv.set_immutable()
            groups = [(cp(i+1)-1,) for i in range(self.r)]
            if vv.nrows() < self.d:
                # assume .eval is working for less rows
                comb.append((coeff, FunctionWrap(vv, groups)))
            else:
                comb.append((coeff, FunctionWrap(self.eval, groups, exp_position=1, extra_args=(vv,))))
                level = max(level, self.level[vv])

        self.strategy[v] = (LinearCombination(comb), ())
        self.what[v] = ('multistuffle', lin_comb, D)
        self.level[v] = level + 1

    def multistuffles(self, verbose=True):
        # TODO: we should make all matrices sparse
        known, M_known = self.known_matrix()
        relations, M_rel = self.stuffle_relation_matrix()
        M = matrix(ZZ, M_rel + M_known)
        M = M.sparse_matrix()
        M.set_immutable()
        Vnsym, nsym_to_can = self._all_forms()
        indicesnsym = {v:i for i,v in enumerate(Vnsym)}
        Fnsym = FreeModule(ZZ, len(indicesnsym))
        nrel = len(M_rel)
        for i, v in enumerate(self.V):
            if self.strategy[v] is not None:
                continue
            vnsym = next(v.symmetric())  # ?
            try:
                sol = M.solve_left(Fnsym.gen(indicesnsym[vnsym]))
            except ValueError:
                sol = None
            if sol is not None:
                yield v, vnsym, [(0, x, relations[i]) if i < nrel else (1, x, known[i-nrel]) for i,x in enumerate(sol) if x]

    @cached_method
    def _all_forms(self):
        d = self.d
        r = self.r
        assert 1 <= r <= d
        vectors = list(linear_forms(d))

        if d == r:
            Vnsym = []
            for _, lins in linearly_independent_vectors(vectors, min_size=d, max_size=d):
                Vnsym.append(GeneralizedMultipleZetaFunction.from_matrix(lins, False, True, True, True))
        else:
            # columns are independent
            # rows are sorted
            Vnsym = set()
            for _, lins in linearly_independent_vectors(vectors, min_size=r, max_size=r):
                for p in itertools.permutations(range(r)):
                    v = GeneralizedMultipleZetaFunction.from_matrix(matrix(ZZ, [lins[p[i]] for i in range(r)]).transpose(), False, True, True, True)
                    if v.has_no_zero_row_or_column():
                        Vnsym.add(v)
            Vnsym = list(Vnsym)

        nsym_to_can = {}
        for v in Vnsym:
            vv = v.copy()
            rp, cp = vv.canonicalize()
            vv.set_immutable()
            nsym_to_can[v] = (vv, rp, cp)

        return Vnsym, nsym_to_can

    @cached_method
    def stuffle_relation_matrix(self):
        r"""
        return relations, mat
        """
        Vnsym, nsym_to_can = self._all_forms()
        indices = {v: i for i,v in enumerate(Vnsym)}
        F = VectorSpace(QQ, len(Vnsym))
        assert len(Vnsym) == len(indices)

        relations = []
        mat = []
        for v in Vnsym:
            for subset, stuffle in v.stuffles(only_highest_weight=True, sort_rows=True):
                relations.append((v, subset[:]))
                mat.append(F.gen(indices[v]) - sum(F.gen(indices[vv]) for vv in stuffle))

        return relations, mat

    def known_matrix(self):
        Vnsym, _ = self._all_forms()
        indices = {v: i for i,v in enumerate(Vnsym)}
        F = VectorSpace(QQ, len(Vnsym))
        elts = []
        mat = []
        for v, status in self.strategy.items():
            if status is not None:
                for vv in v.symmetric():
                    mat.append(F.gen(indices[vv]))
                    elts.append(vv)
        return elts, mat

def print_array(level, explanations):
    width = max(len(L.lin_prod_string()) for L in level)

    line = '{L:<{width}} {level:<2} {handling:<10} {neighbors}'
    for L in sorted(level, key=lambda x: 100000 if level[x] is None else level[x]):
        Ls = '{L:<{width}}'.format(L=L.lin_prod_string(), width=width)
        if level[L] is None:
            print(line.format(L=Ls, width=width, level='oo', handling='x', neighbors=''))
            continue

        handling = explanations[L][0]
        N = explanations[L][2]
        Ns = ', '.join(vv.lin_prod_string()+':'+str(level[vv]) for vv in N)
        print(line.format(L=Ls, width=width, level=level[L], handling=handling, neighbors=Ns))


def print_stuffle3():
    gens = V3gens()
    simplices = [c for c in itertools.combinations(gens,3) if matrix(c).rank() == 3]
    todo = []
    for s1,s2,s3 in simplices:
        if s1[0] + s1[1] <= 1 and s2[0] + s2[1] <= 1 and s3[0] + s3[1] <= 1:
            todo.append((s1,s2,s3,0,1,2))
        if s1[0] + s1[2] <= 1 and s2[0] + s2[2] <= 1 and s3[0] + s3[2] <= 1:
            todo.append((s1,s2,s3,0,2,1))
        if s1[1] + s1[2] <= 1 and s2[1] + s2[2] <= 1 and s3[1] + s3[2] <= 1:
            todo.append((s1,s2,s3,1,2,0))

    a,b,c = ZZ['a,b,c'].gens()
    for s1,s2,s3,i,j,k in todo:
        D = ((s1,a),(s2,b),(s3,c))

        # big subsimplex 1
        m = matrix([s1,s2,s3]).transpose()
        m[i] += m[j]
        t1,t2,t3 = m.columns()
        D1 = clean_term(3, ((t1,a),(t2,b),(t3,c)))

        # big subsimplex 2
        m = matrix([s1,s2,s3]).transpose()
        m[j] += m[i]
        t1,t2,t3 = m.columns()
        D2 = clean_term(3, ((t1,a),(t2,b),(t3,c)))

        # join
        m = matrix([s1,s2,s3]).transpose()
        m[i] += m[j]
        m = m.delete_rows([j])
        t1,t2,t3 = m.columns()
        D3 = clean_term(2, ((t1,a),(t2,b),(t3,c)))

        dat = to_Z3(D)
        if dat[0] == 0:
            continue
        mzv = Z3_to_mzv(*dat)
        if mzv is not None:
            LHS =  "ζ{}".fomat(mzv)
        else:
            LHS = "Z3{}".format(dat)

        dat1 = to_Z3(D1)
        mzv1 = Z3_to_mzv(*dat1)
        if mzv1 is not None:
            RHS1 =  "ζ{}".format(mzv1)
        else:
            RHS1 = "Z3{}".format(dat1)

        dat2 = to_Z3(D2)
        mzv2 = Z3_to_mzv(*dat2)
        if mzv2 is not None:
            RHS2 =  "ζ{}".format(mzv2)
        else:
            RHS2 = "Z3{}".format(dat2)

        dat3 = to_Z2(D3)
        mzv3 = Z2_to_mzv(*dat3)
        if mzv3 is not None:
            RHS3 =  "ζ{}".format(mzv3)
        else:
            RHS3 = "Z2{}".format(dat3)

        print("{} = {} + {} + {}".format(LHS, RHS1, RHS2, RHS3))



def try_algebraic_forward(S):
    L = list(S.algebraic(True))
    if not L:
        return 0
    S.update_algebraic(*choice(L))
    return 1


def try_algebraic(S):
    L = list(S.algebraic(False))
    if not L:
        return 0
    S.update_algebraic(*choice(L))
    return 1


def try_forward_stuffle(S):
    L = list(S.iterated_forward_stuffles())
    if not L:
        return 0
    S.update_forward_stuffle(*choice(L))
    return 1

def try_multistuffle(S):
    L = list(S.multistuffles())
    if not L:
        return 0
    S.update_multistuffle(*choice(L))
    return 1


class Solvers:
    @staticmethod
    def random(d, r):
        import random
        S = GMZVSolver(d, r)
        S.update_multizetas()
        update_functions = [try_algebraic, try_algebraic_forward, try_forward_stuffle, try_multistuffle]
        while any(status is None for status in S.strategy.values()):
            random.shuffle(update_functions)
            for F in update_functions:
                status = F(S)
                if status:
                    break
            if status == 0:
                raise ValueError('not completable strategy')
        return S

    @staticmethod
    def multistuffle_first(d, r, verbose=False):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.generalized_multiple_zeta_values import gmzv_solvers, GeneralizedMultipleZetaFunction
            sage: S = gmzv_solvers.multistuffle_first(3, 3)
            sage: l = GeneralizedMultipleZetaFunction([[1,0,0], [0,1,0], [0,0,1]])
            sage: l.set_immutable()
            sage: assert S.eval(l, (2,2,2)) == Multizeta(2)**3 # optional - mzv
        """
        S = GMZVSolver(d, r)
        S.update_multizetas()
        update_functions = [try_multistuffle, try_algebraic_forward, try_algebraic]
        while any(status is None for status in S.strategy.values()):
            if verbose:
                print("  {}/{}".format(len(S.known()), len(S.strategy)))
            for F in update_functions:
                res = F(S)
                if res:
                    break
            if not res:
                print("NOT COMPLETABLE")
                return S
        return S

    @staticmethod
    def reasonable(d, r, verbose=False):
        S = GMZVSolver(d, r)
        S.update_multizetas()
        update_functions = [try_forward_stuffle, try_algebraic_forward, try_multistuffle, try_algebraic]
        while any(status is None for status in S.strategy.values()):
            if verbose:
                print("  {}/{}".format(len(S.known()), len(S.strategy)))
            for F in update_functions:
                res = F(S)
                if res:
                    break
            if not res:
                print("NOT COMPLETABLE")
                return S
        return S

gmzv_solvers = Solvers()
