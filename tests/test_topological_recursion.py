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
import itertools
from sage.all import ZZ, QQ, factorial, multinomial, parent, Partitions, prod

def test_kontsevich():
    from surface_dynamics.topological_recursion import KontsevichTR
    K = KontsevichTR()
    # dim 0
    assert K.F(0, 3, [0,0,0]) == 1

    # dim 1
    assert K.F(1, 1, [0]) == 0
    assert K.F(1, 1, [1]) == QQ((1,8))
    assert K.F(0, 4, [1,0,0,0]) == 3

    # dim 2
    assert K.F(0, 5, [2,0,0,0,0]) == 15
    assert K.F(0, 5, [1,1,0,0,0]) == 18
    assert K.F(1, 2, [2,0]) == QQ((5,8))
    assert K.F(1, 2, [1,1]) == QQ((3,8))

    # dim 3
    assert K.F(2, 1, [4]) == QQ((105,128))

    # dim 4
    assert K.F(2, 2, [5,0]) == QQ((1155,128))
    assert K.F(2, 2, [4,1]) == QQ((945,128))
    assert K.F(2, 2, [3,2]) == QQ((1015,128))

    # genus zero
    assert K.F(0, 3, (0, 0, 0)) == QQ((1, 1))
    assert K.F(0, 4, (1, 0, 0, 0)) == QQ((3, 1))
    assert K.F(0, 5, (2, 0, 0, 0, 0)) == QQ((15, 1))
    assert K.F(0, 5, (1, 1, 0, 0, 0)) == QQ((18, 1))
    assert K.F(0, 6, (3, 0, 0, 0, 0, 0)) == QQ((105, 1))
    assert K.F(0, 6, (2, 1, 0, 0, 0, 0)) == QQ((135, 1))
    assert K.F(0, 6, (1, 1, 1, 0, 0, 0)) == QQ((162, 1))
    assert K.F(0, 7, (4, 0, 0, 0, 0, 0, 0)) == QQ((945, 1))
    assert K.F(0, 7, (3, 1, 0, 0, 0, 0, 0)) == QQ((1260, 1))
    assert K.F(0, 7, (2, 2, 0, 0, 0, 0, 0)) == QQ((1350, 1))
    assert K.F(0, 7, (2, 1, 1, 0, 0, 0, 0)) == QQ((1620, 1))
    assert K.F(0, 7, (1, 1, 1, 1, 0, 0, 0)) == QQ((1944, 1))

    # genus one
    assert K.F(1, 1, (1,)) == QQ((1, 8))
    assert K.F(1, 2, (2, 0)) == QQ((5, 8))
    assert K.F(1, 2, (1, 1)) == QQ((3, 8))
    assert K.F(1, 3, (3, 0, 0)) == QQ((35, 8))
    assert K.F(1, 3, (2, 1, 0)) == QQ((15, 4))
    assert K.F(1, 3, (1, 1, 1)) == QQ((9, 4))
    assert K.F(1, 4, (4, 0, 0, 0)) == QQ((315, 8))
    assert K.F(1, 4, (3, 1, 0, 0)) == QQ((315, 8))
    assert K.F(1, 4, (2, 2, 0, 0)) == QQ((75, 2))
    assert K.F(1, 4, (2, 1, 1, 0)) == QQ((135, 4))
    assert K.F(1, 4, (1, 1, 1, 1)) == QQ((81, 4))

    # genus two
    assert K.F(2, 1, [4]) == QQ((105, 128))
    assert K.F(2, 2, [5, 0]) == QQ((1155, 128))
    assert K.F(2, 2, [4, 1]) == QQ((945, 128))
    assert K.F(2, 2, [3, 2]) == QQ((1015, 128))
    assert K.F(2, 3, [6, 0, 0]) == QQ((15015, 128))
    assert K.F(2, 3, [5, 1, 0]) == QQ((3465, 32))
    assert K.F(2, 3, [4, 2, 0]) == QQ((3465, 32))
    assert K.F(2, 3, [3, 3, 0]) == QQ((7105, 64))
    assert K.F(2, 3, [4, 1, 1]) == QQ((2835, 32))
    assert K.F(2, 3, [3, 2, 1]) == QQ((3045, 32))
    assert K.F(2, 4, [7, 0, 0, 0]) == QQ((225225, 128))
    assert K.F(2, 4, [6, 1, 0, 0]) == QQ((225225, 128))
    assert K.F(2, 4, [5, 2, 0, 0]) == QQ((3465, 2))
    assert K.F(2, 4, [5, 1, 1, 0]) == QQ((51975, 32))
    assert K.F(2, 4, [4, 3, 0, 0]) == QQ((112455, 64))
    assert K.F(2, 4, [4, 2, 1, 0]) == QQ((51975, 32))
    assert K.F(2, 4, [4, 1, 1, 1]) == QQ((42525, 32))
    assert K.F(2, 4, [3, 3, 1, 0]) == QQ((106575, 64))
    assert K.F(2, 4, [3, 2, 2, 0]) == QQ((13125, 8))
    assert K.F(2, 4, [3, 2, 1, 1]) == QQ((45675, 32))
    assert K.F(2, 4, [2, 2, 2, 1]) == QQ((23625, 16))

    # genus three
    assert K.F(3, 1, [7]) == QQ((25025, 1024))
    assert K.F(3, 2, [8, 0]) == QQ((425425, 1024))
    assert K.F(3, 2, [7, 1]) == QQ((375375, 1024))
    assert K.F(3, 2, [6, 2]) == QQ((385385, 1024))
    assert K.F(3, 2, [5, 3]) == QQ((193655, 512))
    assert K.F(3, 2, [4, 4]) == QQ((191205, 512))
    assert K.F(3, 3, [9, 0, 0]) == QQ((8083075, 1024))
    assert K.F(3, 3, [8, 1, 0]) == QQ((3828825, 512))
    assert K.F(3, 3, [7, 2, 0]) == QQ((3828825, 512))
    assert K.F(3, 3, [7, 1, 1]) == QQ((3378375, 512))
    assert K.F(3, 3, [6, 3, 0]) == QQ((7732725, 1024))
    assert K.F(3, 3, [6, 2, 1]) == QQ((3468465, 512))
    assert K.F(3, 3, [5, 4, 0]) == QQ((1923075, 256))
    assert K.F(3, 3, [5, 3, 1]) == QQ((1742895, 256))
    assert K.F(3, 3, [5, 2, 2]) == QQ((883575, 128))
    assert K.F(3, 3, [4, 4, 1]) == QQ((1720845, 256))
    assert K.F(3, 3, [4, 3, 2]) == QQ((1765575, 256))
    assert K.F(3, 3, [3, 3, 3]) == QQ((3570875, 512))

def test_kontsevich_genus0_formula():
    from surface_dynamics.topological_recursion.kontsevich import psi_correlator

    # d1 + ... + dn = n-3
    # equivalently (d1+1) + ... + (dn+1) = 2n-3
    # <tau_{d1} ... tau_{dn}> = (n-3)! / prod d_i!
    for n in range(3, 10):
        for p in Partitions(2*n-3, length=n):
            p = tuple(i-1 for i in p)
            val = QQ((factorial(n-3), prod(factorial(d) for d in p)))
            assert psi_correlator(*p) == val, p

def test_kontsevich_genus1_formula():
    # from Andersen-Borot-Charbonnier-Delecroix-Giacchetto-Lewanski-Wheeler
    # appendix A

    from surface_dynamics.topological_recursion.kontsevich import psi_correlator

    # d1 + ... + dn = n
    # equivalently (d1+1) + ... + (dn+1) = 2n
    for n in range(1, 8):
        for p in Partitions(2*n, length=n):
            p = tuple(i-1 for i in p)
            s = sum(multinomial([i-j for i,j in zip(p,diff)]) * factorial(sum(diff) - 2) for diff in itertools.product([0,1], repeat=n) if sum(diff) >= 2)
            assert 24 * psi_correlator(*p) == multinomial(p) - s, (p, s)

def test_kontsevich_witten_formula():
    from surface_dynamics.topological_recursion.kontsevich import psi_correlator

    # <tau_{3g-2}>_{g,1} = 1 / (24^g g!)
    from surface_dynamics.topological_recursion import KontsevichTR
    K = KontsevichTR()
    for g in range(1, 10):
        assert 24**g * factorial(g) * psi_correlator(3*g-2) == 1, g

def test_masur_veech():
    from surface_dynamics.topological_recursion import MasurVeechTR

    MV = MasurVeechTR()

    for g,n,value in [(0, 4, QQ((2,1))),
            (0, 5, QQ((1,1))), (0, 6, QQ((1, 2))),
            (0, 7, QQ((1, 4))), (1, 1, QQ((2,3))),
            (1, 2, QQ((1, 3))), (1, 3, QQ((11, 60))),
            (1, 4, QQ((1, 10))), (1, 5, QQ((163, 3024))),
            (2, 1, QQ((29, 840))), (2, 2, QQ((337, 18144))),
            (3, 1, QQ((4111, 2223936)))]:
        coeff = 2**(4*g-2+n) * ZZ(4*g-4+n).factorial() / ZZ(6*g-7+2*n).factorial()
        mv = MV.F(g, n, (0,)*n)
        assert coeff * mv == value, (g, n, mv, value)

def test_masur_veech_edge_weight():
    from surface_dynamics.topological_recursion import MasurVeechTR
    from sage.all import PolynomialRing

    R = PolynomialRing(QQ, 't')
    t = R.gen()
    MV = MasurVeechTR(edge_weight=t)
    for g,n,value in [
                (0, 3, R.one()),
                (0, 4, t),
                (0, 5, QQ((4,9))*t + QQ((5,9))*t**2),
                (0, 6, QQ((8,27))*t + QQ((4,9))*t**2 + QQ((7,27))*t**3),
                (1, 1, t),
                (1, 2, QQ((5,9))*t + QQ((4,9))*t**2),
                (1, 3, QQ((4,9))*t + QQ((13,33))*t**2 + QQ((16,99))*t**3),
                (2, 1, QQ((76,261))*t + QQ((125,261))*t**2 + QQ((440,2349))*t**3 + QQ((100,2349))*t**4),
                (2, 2, QQ((296,1011))*t + QQ((19748,45495))*t**2 + QQ((9127,45495))*t**3 + QQ((560,9099))*t**4 + QQ((100,9099))*t**5)]:
        p = MV.F(g, n, (0,)*n)
        assert parent(p) is R
        p /= p(1)
        assert p == value, (g, n, p, value)
