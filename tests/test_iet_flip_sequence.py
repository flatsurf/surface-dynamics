#!/usr/bin/env python
#*****************************************************************************
#       Copyright (C) 2021 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

import pytest
import itertools
import random

from surface_dynamics import iet
from surface_dynamics.misc.permutation import perm_random

def random_flip_sequence(p, n):
    f = iet.FlipSequence(p)
    for _ in range(n):
        # left or right
        side = random.choice([0, -1])
        if not f._end.has_rauzy_move(0, side):
            winner = 1
        elif not f._end.has_rauzy_move(1, side):
            winner = 0
        else:
            winner = random.choice([0, 1])
        f.rauzy_move(winner, side)
    if random.randint(0, 2):
        f.top_bottom_inverse()
    if random.randint(0, 2):
        f.left_right_inverse()

    f.relabel(perm_random(len(p)))

    return f

@pytest.mark.parametrize("top, bot",
    [
     ('a b', 'b a'),
     ('a b c', 'c b a'),
     ('a b c d', 'd c b a'),
     ('a b c d e f', 'f c b e d a')
     ])
def test_mul(top, bot):
    p = iet.Permutation(top, bot)
    for _ in range(20):
        f1 = random_flip_sequence(p, random.randint(2, 5))
        f2 = random_flip_sequence(f1._end, random.randint(2, 5))
        f3 = random_flip_sequence(f2._end, random.randint(2, 5))

        # path multiplication
        f = (f1 * f2) * f3
        g = f1 * (f2 * f3)
        f._check()
        g._check()
        assert f._start == g._start == f1._start
        assert f._end == g._end == f3._end
        assert f == g

        # matrix multiplication
        m1 = f1.matrix()
        m2 = f2.matrix()
        m3 = f3.matrix()
        m = f.matrix()
        assert m == m1 * m2 * m3

        # orbit substitution multiplication
        s1 = f1.substitution()
        s2 = f2.substitution()
        s3 = f3.substitution()
        s = f.substitution()
        assert s1.incidence_matrix() == m1
        assert s2.incidence_matrix() == m2
        assert s3.incidence_matrix() == m3
        assert s.incidence_matrix() == m

        if not f1._top_bottom_inverse and \
           not f2._top_bottom_inverse and \
           not f3._top_bottom_inverse:
            assert s == s1 * s2 * s3

        # dual substitution multiplication
        s1 = f1.dual_substitution()
        s2 = f2.dual_substitution()
        s3 = f3.dual_substitution()
        s = f.dual_substitution()

        if not f1._left_right_inverse and \
           not f2._left_right_inverse and \
           not f3._left_right_inverse:
               assert s == s3 * s2 * s1
