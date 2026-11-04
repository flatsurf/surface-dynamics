#!/usr/bin/env python
#*****************************************************************************
#       Copyright (C) 2026 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

import pytest

from surface_dynamics import iet, Stratum


# name -> (stratum, component); the search takes no left degree on a
# quadratic component
DIAGRAMS = {
    'Q1111hyp': ([1, 1, 1, 1], 'hyp'),
    'Q22hyp': ([2, 2], 'hyp'),
    'Q5m1': ([5, -1], 'unique'),
}


def diagram(name):
    stratum, kind = DIAGRAMS[name]
    S = Stratum(stratum, k=2)
    if kind == 'hyp':
        return S.hyperelliptic_component()
    return S.unique_component()


def counts(walks, N):
    r"""Number of walks of each length 1..N."""
    res = [0] * N
    for f in walks:
        res[len(f) - 1] += 1
    return res


def test_boissy_lanneau_example_2_8():
    # Boissy-Lanneau 2009, Example 2.8 (p. 784)
    p = iet.GeneralizedPermutation('A B C B D', 'D E A C E')
    assert p.rauzy_move('t') == iet.GeneralizedPermutation('A B C B D', 'D E E A C')
    assert p.rauzy_move('b') == iet.GeneralizedPermutation('A B C B', 'D D E A C E')


def test_boissy_lanneau_example_2_9():
    # Boissy-Lanneau 2009, Example 2.9 (p. 784): the bottom move is undefined
    p = iet.GeneralizedPermutation('A B A', 'B D C C D')
    assert p.rauzy_move('t') == iet.GeneralizedPermutation('D A B A', 'B D C C')
    assert p.has_rauzy_move('t')
    assert not p.has_rauzy_move('b')
    with pytest.raises(ValueError):
        p.rauzy_move('b')


def test_boissy_lanneau_example_2_10():
    # Boissy-Lanneau 2009, Example 2.10 (p. 784): neither move is defined
    p = iet.GeneralizedPermutation('A B A', 'B C C')
    assert not p.has_rauzy_move('t')
    assert not p.has_rauzy_move('b')


def test_boissy_lanneau_figure_a2():
    # Boissy-Lanneau 2009, Figure A2 (p. 814): a Rauzy class of
    # Q(-1, -1, -1, 7) has 28884 generalized permutations up to relabelling
    C = Stratum([7, -1, -1, -1], k=2).unique_component()
    assert C.permutation_representative().rauzy_diagram().cardinality() == 28884


# name -> (vertices of the reduced Rauzy diagram, vertices with a single
# defined move, marking of the representative, degree of the singularity at
# the left end); from rauzy_diagram(), has_rauzy_move() and marked_profile()
REDUCED_DIAGRAMS = {
    'Q1111hyp': (380, 112, '3o3', 1),
    'Q22hyp': (73, 26, '4o4', 2),
    'Q5m1': (440, 150, '7|5', 5),
}


@pytest.mark.parametrize('name', list(REDUCED_DIAGRAMS))
def test_reduced_diagrams(name):
    V, one_move, marking, left_degree = REDUCED_DIAGRAMS[name]
    p = diagram(name).permutation_representative()
    D = p.rauzy_diagram()
    assert D.cardinality() == V
    assert sum(1 for q in D.vertices() if q.has_rauzy_move(0) != q.has_rauzy_move(1)) == one_move
    assert all(q.has_rauzy_move(0) or q.has_rauzy_move(1) for q in D.vertices())
    m = p.marked_profile().marking()
    assert str(m) == marking
    # angles in units of pi: a singularity of angle k pi has degree k - 2
    assert m.left() - 2 == left_degree


ALL_COUNTS = {
    # brute force over all words of defined moves from every reduced vertex,
    # agreeing with the trace formula (1/n) sum_{e|n} mu(n/e) tr(A^e) on the
    # adjacency matrix A of the reduced diagram
    'Q1111hyp': [76, 34, 34, 32, 54, 58, 104, 134, 212, 302, 480, 706],
    'Q22hyp': [12, 6, 2, 6, 4, 10, 8, 21, 22, 54, 52, 123, 130, 293, 342, 709],
    'Q5m1': [50, 24, 30, 34, 26, 40, 38, 72, 88, 144, 186, 322, 404, 713],
}


@pytest.mark.xfail(raises=ValueError, strict=True)
@pytest.mark.parametrize('name', list(ALL_COUNTS))
def test_counts_all_walks(name):
    C = diagram(name)
    expected = ALL_COUNTS[name]
    N = len(expected)
    assert counts(iet.rauzy_closed_walks(C, N, complete=False), N) == expected


# brute force over all words of defined moves from every reduced vertex,
# keeping the primitive closed walks whose closed flip sequence is complete;
# up to the shortest length with a complete walk
COMPLETE_COUNTS = {
    'Q1111hyp': [0] * 12 + [6],
    'Q22hyp': [0] * 7 + [2, 4, 6, 8, 29, 42],
    'Q5m1': [0] * 13 + [9],
}


@pytest.mark.xfail(raises=ValueError, strict=True)
@pytest.mark.parametrize('name', list(COMPLETE_COUNTS))
def test_counts_complete(name):
    C = diagram(name)
    expected = COMPLETE_COUNTS[name]
    N = len(expected)
    assert counts(iet.rauzy_closed_walks(C, N), N) == expected


COMPLETE_COUNTS_SLOW = {
    'Q1111hyp': [0] * 12 + [6, 14, 48, 88],
    'Q22hyp': [0] * 7 + [2, 4, 6, 8, 29, 42, 88, 150, 299, 488, 918],
    'Q5m1': [0] * 13 + [9, 44, 129, 216],
}


@pytest.mark.slow
@pytest.mark.xfail(raises=ValueError, strict=True)
@pytest.mark.parametrize('name', list(COMPLETE_COUNTS_SLOW))
def test_counts_complete_slow(name):
    C = diagram(name)
    expected = COMPLETE_COUNTS_SLOW[name]
    N = len(expected)
    assert counts(iet.rauzy_closed_walks(C, N), N) == expected
