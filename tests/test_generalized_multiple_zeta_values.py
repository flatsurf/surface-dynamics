#!/usr/bin/env python
#*****************************************************************************
#       Copyright (C) 2020 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************
import pytest
import sage.all
import cypari2

try:
    import sage.modular.multiple_zeta
    from sage.all import Multizetas
    from surface_dynamics.misc.generalized_multiple_zeta_values import linear_forms, handle_term, is_convergent, Z2, to_Z2, Z3, to_Z3, is_Z3_convergent, clean_term, convergent_multizeta
    # put more multiple zeta values in the cache
    from sage.modular.multiple_zeta import Values
    Values.reset(max_weight=12)
except ImportError:
    pytestmark = pytest.mark.skip(reason="sage.modular.multiple_zeta not present")
except cypari2.handle_error.PariError:
    pytestmark = pytest.mark.skip(reason="sage.modular.multiple_zeta not functional")

import itertools
import random

from sage.all import FreeModule, ZZ, QQ, matrix

def V1gens():
    return linear_forms(1)

def V2gens():
    return linear_forms(2)

def V3gens():
    va, vb, vd, vc, ve, vf, vg = linear_forms(3)
    return (va,vb,vc,vd,ve,vf,vg)

def test_convergent_multizeta():
    M = Multizetas(QQ)
    assert convergent_multizeta((0,3)) == M((2,)) - M((3,))
    assert convergent_multizeta((0,4)) == M((3,)) - M((4,))
    # binomal(x-1, 1) = x - 1
    assert convergent_multizeta((0,2,6)) == M((1,6)) - M((2,6))
    # binomial(x-1, 2) = 1/2 x^2 - 3/2 x + 1
    assert convergent_multizeta((0,0,6)) == (M((4,)) - 3*M((5,)) + 2*M((6,))) / 2
    # binomial(x-y-1,1) = x - y - 1
    assert convergent_multizeta((2,0,6)) == M((2,5)) - M((1,6)) - M((2,6))

    assert convergent_multizeta((0,0,0,5)) == 4 * convergent_multizeta((-1,0,0,6))
    assert convergent_multizeta((0,0,0,6)) == 4 * convergent_multizeta((-1,0,0,7))
    assert convergent_multizeta((0,0,0,0,8)) == 5 * convergent_multizeta((-1,0,0,0,9))

    # Stuffle relation for zeta(a,b) * zeta(c,d)
    for [a,b,c,d] in [[1,2,1,2],[-1,4,-1,4],[-1,4,-1,5],[-1,4,-2,5],[-2,5,-2,5],[-2,6,-2,6]]:
        P = convergent_multizeta((a,b)) * convergent_multizeta((c,d))
        S = sum(convergent_multizeta(t) for t in [
            (a,b,c,d),(a,c,b,d),(a,c,d,b),(c,a,b,d),(c,a,d,b),
             (c,d,a,b),(a+c,b,d),(a,b+c,d),(a+c,d,b),(a,c,b+d),
             (c,a+d,b),(c,a,b+d),(a+c,b+d)])
        assert P == S, (a,b,c,d)

def stuffle(A):
    A.set_immutable()
    ans = {A: 1}
    todo = {}
    for i in range(A.nrows()):
        for j in range(i):
            if all(A[i,k] + A[i,j] <= 1 for k in range(A.ncols())):
                do_stuffle(A, i, j)

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

def test_bill_data():
    M = Multizetas(QQ)
    bill_data = [
    ([0,0,2,0,3,2,0],[[[4,3],-3],[[5,2],-21],[[6,1],-83],[[2,3,2],1],[[2,4,1],12],[[3,2,2],11],[[3,3,1],36],[[4,1,2],30],[[4,2,1],66],[[5,1,1],132]]),
    ([0,0,2,0,2,3,0],[[[4,3],-3],[[5,2],-21],[[6,1],-83],[[2,3,2],1],[[2,4,1],12],[[3,2,2],11],[[3,3,1],36],[[4,1,2],30],[[4,2,1],66],[[5,1,1],132]]),
    ([0,0,3,0,2,2,0],[[[5,2],9],[[6,1],58],[[2,2,3],2],[[2,4,1],-12],[[3,1,3],4],[[3,2,2],-6],[[3,3,1],-24],[[4,1,2],-18],[[4,2,1],-30],[[5,1,1],-60]]),
    ([0,1,1,3,2,0,0],[[[4,3],-2],[[5,2],-10],[[6,1],-22],[[2,4,1],1],[[3,2,2],6],[[3,3,1],14],[[4,1,2],15],[[4,2,1],30],[[5,1,1],47]]),
    ([0,1,1,2,3,0,0],[[[4,3],-2],[[5,2],-10],[[6,1],-22],[[2,4,1],1],[[3,2,2],6],[[3,3,1],14],[[4,1,2],15],[[4,2,1],30],[[5,1,1],47]]),
    ([0,2,1,2,2,0,0],[[[5,2],-6],[[6,1],-19],[[2,4,1],2],[[3,2,2],3],[[3,3,1],8],[[4,1,2],8],[[4,2,1],22],[[5,1,1],46]]),
    ([0,2,1,1,3,0,0],[[[4,3],-1],[[5,2],-3],[[6,1],-18],[[2,4,1],3],[[3,3,1],8],[[4,1,2],5],[[4,2,1],22],[[5,1,1],57]]),
    ([0,1,2,3,1,0,0],[[[4,3],-1],[[5,2],-3],[[6,1],-18],[[2,4,1],3],[[3,3,1],8],[[4,1,2],5],[[4,2,1],22],[[5,1,1],57]]),
    ([0,1,2,2,2,0,0],[[[5,2],-6],[[6,1],-19],[[2,4,1],2],[[3,2,2],3],[[3,3,1],8],[[4,1,2],8],[[4,2,1],22],[[5,1,1],46]]),
    ([0,2,2,2,1,0,0],[[[5,2],-1],[[6,1],-9],[[2,4,1],2],[[3,2,2],3],[[3,3,1],8],[[4,1,2],6],[[4,2,1],12],[[5,1,1],22]]),
    ([0,2,2,1,2,0,0],[[[5,2],-1],[[6,1],-9],[[2,4,1],2],[[3,2,2],3],[[3,3,1],8],[[4,1,2],6],[[4,2,1],12],[[5,1,1],22]])]

    for e, S in bill_data:
        z1 = Z3(*e)
        z2 = sum(coeff * M(e[::-1]) for e,coeff in S)
        assert z1 == z2

def test_convergence():
    v, = V1gens()
    assert not is_convergent(1, ((v,1),))
    assert is_convergent(1, ((v,2),))

    va,vb,vc = V2gens()
    assert is_convergent(2, ((va,2),(vb,2)))
    assert is_convergent(2, ((va,1),(vc,2)))
    assert is_convergent(2, ((vc,3),))

    assert not is_convergent(2, ((va,3),))
    assert not is_convergent(2, ((vb,3),))
    assert not is_convergent(2, ((va,1),(vb,1)))
    assert not is_convergent(2, ((va,2),(vb,1)))
    assert not is_convergent(2, ((va,1),(vb,2)))
    assert not is_convergent(2, ((vc,2),))
    assert not is_convergent(2, ((va,1),(vc,1)))

    va,vb,vc,vd,ve,vf,vg = V3gens()

    assert is_convergent(3, ((va,2),(vc,2),(vd,2)))
    assert is_convergent(3, ((vg,4),))
    assert is_convergent(3, ((va,1),(vd,1),(vg,2)))
    assert is_convergent(3, ((vb,1),(vd,1),(vg,2)))
    assert is_convergent(3, ((vb,1),(ve,1),(vg,2)))

    assert not is_convergent(3, ((vd,1),(ve,1),(vf,1)))
    assert not is_convergent(3, ((vg,1),))

def test_zeta_positive3():
    exponents = [
            # weight 4
            [1,0,0,0,0,1,2],[0,0,0,0,0,0,4],
            # weight 5
            [0,0,0,0,1,2,2],[1,0,0,1,0,2,1],[0,0,0,1,2,2,0],
            # weight 6
            [0,0,0,0,1,2,3],[0,0,0,0,2,2,2],[1,1,0,0,1,3,0],
            # weight 7
            [0,0,0,3,2,0,2], [0,0,2,0,3,2,0], [1,1,1,1,1,1,1]]
    for e in exponents:
        v = Z3(*e)
        assert v.n() > 0, (e, v, v.n())

def Z2_to_mzv(a,b,c):
    if a == 0:
        return b,c
    if b == 0:
        return a,c

def test_stuffle3():
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
    a1 = 2
    a2 = 3
    a3 = 4
    for s1,s2,s3,i,j,k in todo:
        D = ((s1,a1),(s2,a2),(s3,a3))
        z = handle_term(3, D)

        # big subsimplex 1
        m = matrix([s1,s2,s3]).transpose()
        m[i] += m[j]
        t1,t2,t3 = m.columns()
        D1 = clean_term(3, ((t1,a1),(t2,a2),(t3,a3)))
        z1 = handle_term(3, D1)

        # big subsimplex 2
        m = matrix([s1,s2,s3]).transpose()
        m[j] += m[i]
        t1,t2,t3 = m.columns()
        D2 = clean_term(3, ((t1,a1),(t2,a2),(t3,a3)))
        z2 = handle_term(3, D2)

        # join
        m = matrix([s1,s2,s3]).transpose()
        m[i] += m[j]
        m = m.delete_rows([j])
        t1,t2,t3 = m.columns()
        D3 = clean_term(2, ((t1,a1),(t2,a2),(t3,a3)))
        z3 = handle_term(2, D3)

        assert z == z1 + z2 + z3

def test_zeta_symmetries3():
    for _ in range(50):
        t = [random.randint(0,3) for _ in range(7)]
        if sum(t) > 10:
            continue
        if not is_Z3_convergent(*t):
            continue
        t0 = t
        z = Z3(*t)
        for _ in range(10):
            a,b,c,d,e,f,g = t
            # x^a y^b z^c (x+y)^d (x+z)^e (y+z)^f (x+y+z)^g
            r = random.randint(0,2)
            if r == 0:
                # x <-> y
                t = b,a,c,d,f,e,g
            elif r == 1:
                # x <-> z
                t = c,b,a,f,e,d,g
            else:
                # y <-> z
                t = a,c,b,e,d,f,g
            assert z == Z3(*t), (t0, t)

def test_zeta_relations3():
    gens = V3gens()
    #      a  b  c  d  e  f  g
    relations = [
         # (x) + (y+z) = (x+y+z)
         [ 1, 0, 0, 0, 0, 1,-1],
         [ 0, 1, 0, 0, 1, 0,-1],
         [ 0, 0, 1, 1, 0, 0,-1],
         # (x) + (y) = (x+y)
         [ 1, 1, 0,-1, 0, 0, 0],
         [ 1, 0, 1, 0,-1, 0, 0],
         [ 0, 1, 1, 0, 0,-1, 0],
         # (x) + (y) + (z) = (x+y+z)
         [ 1, 1, 1, 0, 0, 0, -1],
         # (x+y) + (x+z) + (y+z) = 2(x+y+z)
         [ 0, 0, 0, 1, 1, 1, -2],
         # (x) + (y) + (x+z) + (y+z) = 2(x+y+z)
         [1, 1, 0, 0, 1, 1, -2],
         # (x) + (y+z) = (y) + (x+z)
         [1, -1, 0, 0, -1, 1, 0],
         # (x) + (x+y+z) = (x+y) + (x+z)
         [1, 0, 0, -1, -1, 0, 1]
         ]
    for r in relations:
        assert sum(r[i]*gens[i] for i in range(7)).is_zero()

    def expo_to_den_tuple(e):
        return tuple((gens[i], e[i]) for i in range(7) if e[i])

    exponents = [[1,0,0,0,1,1,2], [1,0,0,1,2,2,0],[1,1,1,1,1,1,1],[1,0,0,0,0,1,4], [0,0,0,1,1,1,0], [0,0,0,0,1,2,3],[1,0,0,1,0,2,2],[0,0,0,1,2,2,1],[1,0,1,0,1,0,2],[0,0,0,0,0,0,4]]
    for e in exponents:
        for r in relations:
            nc = False
            ee = e[:]
            for i in range(7):
                if r[i] and not ee[i]:
                    ee[i] += 1
            z1 = 0
            z2 = 0
            for i in range(7):
                if r[i]:
                    eee = ee[:]
                    if eee[i] == 0:
                        nc = True
                        break
                    eee[i] -= 1
                    expo_to_den_tuple(eee)
                    if not is_Z3_convergent(*eee):
                        nc = True
                        break
                    z1 += r[i] * handle_term(3, expo_to_den_tuple(eee))
                    z2 += r[i] * Z3(*eee)
            if nc:
                print("non-convergent relation: {} {}".format(ee,r))
            else:
                assert z1.is_zero() and z2.is_zero(), (z1, z2, ee, r)
