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
from copy import copy

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


def quadratic_permutation(name):
    from surface_dynamics import Stratum
    if name == 'Q_2(1^4)^hyp':
        C = Stratum([1, 1, 1, 1], k=2).hyperelliptic_component()
    elif name == 'Q_2(2^2)^hyp':
        C = Stratum([2, 2], k=2).hyperelliptic_component()
    elif name == 'Q_2(5,-1)^nonhyp':
        C = Stratum([5, -1], k=2).non_hyperelliptic_component()
    else:
        raise ValueError(name)
    return C.permutation_representative(reduced=False)


def random_valid_flip_sequence(p, n, rng):
    r"""
    Return a flip sequence from ``p`` of at most ``n`` Rauzy moves, each
    chosen at random among the moves (left or right) defined at its start.
    """
    f = iet.FlipSequence(p)
    for _ in range(n):
        moves = [(winner, side) for winner in (0, 1) for side in (-1, 0)
                 if f._end.has_rauzy_move(winner, side)]
        if not moves:
            break
        f.rauzy_move(*rng.choice(moves))
    if rng.randint(0, 1):
        f.top_bottom_inverse()
    if rng.randint(0, 1):
        f.left_right_inverse()
    relabelling = list(range(len(p)))
    rng.shuffle(relabelling)
    f.relabel(relabelling)
    return f


QUADRATIC = ['Q_2(1^4)^hyp', 'Q_2(2^2)^hyp', 'Q_2(5,-1)^nonhyp']


@pytest.mark.parametrize("name", QUADRATIC)
def test_quadratic_check_and_reversed(name):
    rng = random.Random(1)
    p = quadratic_permutation(name)
    for _ in range(20):
        f = random_valid_flip_sequence(p, rng.randint(0, 12), rng)
        f._check()
        forward = [(copy(q), winner, side) for q, winner, side in f]
        backward = [(copy(q), winner, side) for q, winner, side in reversed(f)]
        assert len(forward) == len(f)
        assert backward[::-1] == forward
        assert hash(copy(f)) == hash(f) and copy(f) == f


@pytest.mark.parametrize("name", QUADRATIC)
def test_quadratic_mul(name):
    rng = random.Random(2)
    p = quadratic_permutation(name)
    for _ in range(20):
        f1 = random_valid_flip_sequence(p, rng.randint(2, 5), rng)
        f2 = random_valid_flip_sequence(f1._end, rng.randint(2, 5), rng)
        f3 = random_valid_flip_sequence(f2._end, rng.randint(2, 5), rng)

        f = (f1 * f2) * f3
        g = f1 * (f2 * f3)
        f._check()
        g._check()
        assert f == g
        assert f._start == f1._start and f._end == f3._end
        assert f.matrix() == f1.matrix() * f2.matrix() * f3.matrix()


@pytest.mark.parametrize("name", QUADRATIC)
def test_quadratic_undefined_move(name):
    rng = random.Random(3)
    p = quadratic_permutation(name)
    refused = 0
    for _ in range(50):
        f = random_valid_flip_sequence(p, rng.randint(0, 12), rng)
        for winner in 'tb':
            for side in 'rl':
                if f._end.has_rauzy_move(winner, side):
                    continue
                g = copy(f)
                with pytest.raises(ValueError, match="Rauzy induction is not well defined"):
                    g.rauzy_move(winner, side)
                with pytest.raises(ValueError, match="Rauzy induction is not well defined"):
                    g.rauzy_move(winner, side, iterations=3)
                assert g == f and g._end == f._end
                refused += 1
    assert refused > 0


def test_quadratic_undefined_move_example():
    # Example 2.9 of Boissy-Lanneau: the bottom right move is not defined,
    # and after the top right move no right move is defined
    p = iet.GeneralizedPermutation('A B A', 'B D C C D')
    with pytest.raises(ValueError):
        iet.FlipSequence(p, 'b')
    with pytest.raises(ValueError):
        iet.FlipSequence(p, 'tt')
    f = iet.FlipSequence(p, 't')
    with pytest.raises(ValueError):
        f.rauzy_move('t', iterations=2)
    assert f == iet.FlipSequence(p, 't')


@pytest.mark.parametrize("name", QUADRATIC)
def test_quadratic_abelian_only(name):
    p = quadratic_permutation(name)
    f = random_valid_flip_sequence(p, 5, random.Random(4))
    for method in ['orbit_substitution', 'interval_substitution', 'substitution',
                   'dual_substitution', 'self_similar_interval_exchange_transformation',
                   'self_similar_iet']:
        with pytest.raises(NotImplementedError):
            getattr(f, method)()
