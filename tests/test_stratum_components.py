#*****************************************************************************
#       Copyright (C) 2020 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

import pytest
import random

from surface_dynamics import AbelianStratum
from surface_dynamics import QuadraticStratum

def check_abelian_component(C, repeat):
    for reduced in [True, False]:
        p = C.permutation_representative(reduced=reduced)
        assert p.stratum_component() == C, (p, C)
        path = []
        choices = ['tl','bl','tr','br']
        for _ in range(repeat):
            x = random.choice(choices)
            path.append(x)
            p = p.rauzy_move(*x)
            assert p.stratum_component() == C, (reduced, p, C, ','.join(path))

def check_quadratic_component(C, repeat):
    for reduced in [True, False]:
        p = C.permutation_representative(reduced=reduced)
        assert p.stratum_component() == C, (p, C)
        path = []
        choices = ['t', 'b']
        for _ in range(repeat):
            if not p.has_rauzy_move('t'):
                x = 'b'
            elif not p.has_rauzy_move('b'):
                x = 't'
            else:
                x = random.choice(choices)
            path.append(x)
            p = p.rauzy_move(x)
            assert p.stratum_component() == C, (reduced, p, C, ''.join(path))

def test_H_6():
    H = AbelianStratum(6)
    check_abelian_component(H.hyperelliptic_component(), 100)
    check_abelian_component(H.even_component(), 100)
    check_abelian_component(H.odd_component(), 100)

def test_H_3_3():
    H = AbelianStratum(3,3)
    check_abelian_component(H.hyperelliptic_component(), 100)
    check_abelian_component(H.non_hyperelliptic_component(), 100)

def test_Q_6_2():
    Q = QuadraticStratum(6, 2)
    check_quadratic_component(Q.hyperelliptic_component(), 100)
    check_quadratic_component(Q.non_hyperelliptic_component(), 100)

def test_Q_9_p():
    Q = QuadraticStratum(9, -1)
    check_quadratic_component(Q.regular_component(), 100)
    check_quadratic_component(Q.irregular_component(), 100)

def test_Q_6_3_p():
    Q = QuadraticStratum(6,3,-1)
    check_quadratic_component(Q.regular_component(), 100)
    check_quadratic_component(Q.irregular_component(), 100)

def test_Q_12():
    Q = QuadraticStratum(12)
    check_quadratic_component(Q.irregular_component(), 100)
    check_quadratic_component(Q.regular_component(), 100)

def testQ_6_6():
    Q = QuadraticStratum(6,6)
    check_quadratic_component(Q.hyperelliptic_component(), 100)
    check_quadratic_component(Q.regular_component(), 100)
    check_quadratic_component(Q.irregular_component(), 100)

@pytest.mark.slow
def test_Q_6_3_3():
    Q = QuadraticStratum(6,3,3)
    check_quadratic_component(Q.hyperelliptic_component(), 100)
    check_quadratic_component(Q.regular_component(), 100)
    check_quadratic_component(Q.irregular_component(), 100)

@pytest.mark.slow
def test_Q_3_3_3_3():
    Q = QuadraticStratum(3,3,3,3)
    check_quadratic_component(Q.hyperelliptic_component(), 100)
    check_quadratic_component(Q.regular_component(), 100)
    check_quadratic_component(Q.irregular_component(), 100)
