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

from sage.all import QQ, polygen

from surface_dynamics import iet, Stratum


def _component(stratum, kind):
    S = Stratum(stratum)
    if kind == 'unique':
        return S.unique_component()
    if kind == 'hyp':
        return S.hyperelliptic_component()
    return S.odd_component()


# name -> (stratum, component, left_degree)
DIAGRAMS = {
    'H2': ([2], 'unique', None),
    'H4hyp': ([4], 'hyp', None),
    'H4odd': ([4], 'odd', None),
    'H22hyp': ([2, 2], 'hyp', None),
    'H22odd': ([2, 2], 'odd', None),
    'H31left3': ([3, 1], 'unique', 3),
    'H31left1': ([3, 1], 'unique', 1),
}


def diagram(name):
    stratum, kind, left_degree = DIAGRAMS[name]
    return _component(stratum, kind), left_degree


def counts(walks, N):
    r"""Number of walks of each length 1..N."""
    res = [0] * N
    for f in walks:
        res[len(f) - 1] += 1
    return res


def test_h2_length_8():
    # the only complete closed walk of length 8 of H(2), the loop ttbtbbtb of
    # the reduced Rauzy diagram, read from its least vertex as btbttbtb
    C = Stratum([2]).unique_component()
    walks = list(iet.rauzy_closed_walks(C, 8, min_length=8))
    assert len(walks) == 1
    f = walks[0]
    assert f.is_closed()
    assert f.is_complete()
    assert f.start() == iet.Permutation([0, 1, 2, 3], [1, 3, 2, 0])
    assert f._simplified_flip_sequence() == 'btbttbtb'
    x = polygen(QQ)
    assert f.matrix().charpoly() == x**4 - 7*x**3 + 13*x**2 - 7*x + 1


ALL_COUNTS = {
    # brute force over all words of moves from every reduced vertex, agreeing
    # with the trace formula (1/n) sum_{e|n} mu(n/e) tr(A^e) on the adjacency
    # matrix A of the reduced diagram
    'H2': [4, 2, 4, 4, 8, 11, 22, 35, 64, 110, 204, 360],
    'H4hyp': [16, 8, 12, 14, 28, 34, 68, 102, 180, 295],
    'H4odd': [44, 22, 24, 35, 56, 75, 120, 190, 300, 514],
    'H22hyp': [32, 16, 24, 28, 54, 68, 134, 200, 352, 574],
    'H22odd': [76, 38, 36, 47, 74, 106, 160, 248, 380, 629],
    'H31left3': [168, 84, 88, 125, 190, 279, 438, 675, 1056, 1763],
    'H31left1': [84, 42, 46, 64, 92, 148, 214, 343, 528, 900],
}


@pytest.mark.parametrize('name', list(ALL_COUNTS))
def test_counts_all_walks(name):
    C, k = diagram(name)
    expected = ALL_COUNTS[name]
    N = len(expected)
    assert counts(iet.rauzy_closed_walks(C, N, complete=False, left_degree=k), N) == expected


# brute force over all words of moves from every reduced vertex, keeping the
# primitive closed walks whose closed flip sequence is complete
COMPLETE_COUNTS = {
    'H2': [0, 0, 0, 0, 0, 0, 0, 1, 6, 14, 36, 89],
    'H22odd': [0] * 12 + [8],
    'H31left1': [0] * 12 + [8],
}


@pytest.mark.parametrize('name', list(COMPLETE_COUNTS))
def test_counts_complete(name):
    C, k = diagram(name)
    expected = COMPLETE_COUNTS[name]
    N = len(expected)
    assert counts(iet.rauzy_closed_walks(C, N, left_degree=k), N) == expected


COMPLETE_COUNTS_SLOW = {
    'H22odd': [0] * 12 + [8, 28],
    'H31left1': [0] * 12 + [8, 16],
    'H4hyp': [0] * 15 + [6],
    'H4odd': [0] * 14,
    'H22hyp': [0] * 14,
    'H31left3': [0] * 14,
}


@pytest.mark.slow
@pytest.mark.parametrize('name', list(COMPLETE_COUNTS_SLOW))
def test_counts_complete_slow(name):
    C, k = diagram(name)
    expected = COMPLETE_COUNTS_SLOW[name]
    N = len(expected)
    assert counts(iet.rauzy_closed_walks(C, N, left_degree=k), N) == expected


# complete walks per length up to 18, from the module before the orbits of a
# base were read off the tree labellings; the relabelling group has several
# orbits (order 3 for H22odd, order 2 for H31left1)
COMPLETE_COUNTS_18 = {
    'H22odd': [0] * 12 + [8, 28, 82, 206, 512, 1200],
    'H31left1': [0] * 12 + [8, 16, 58, 132, 326, 716],
}


@pytest.mark.parametrize('name', list(COMPLETE_COUNTS_18))
def test_counts_complete_18(name):
    C, k = diagram(name)
    expected = COMPLETE_COUNTS_18[name]
    assert sum(expected) == {'H22odd': 2036, 'H31left1': 1256}[name]
    assert counts(iet.rauzy_closed_walks(C, 18, left_degree=k), 18) == expected


def _walk_list(walks):
    return [(str(f.start()), f._simplified_flip_sequence()) for f in walks]


# (diagram, max_length): the walks of complete=False number 828, 20622, 11927
# and 17032, of which 146, 6, 36 and 24 are complete; the relabelling group is
# trivial for H2 and H4hyp, of order 3 for H22odd and 2 for H31left1
@pytest.mark.parametrize('name, N', [('H2', 12), ('H4hyp', 16), ('H22odd', 14), ('H31left1', 14)])
def test_complete_against_filter(name, N):
    C, k = diagram(name)
    filtered = [f for f in iet.rauzy_closed_walks(C, N, complete=False, left_degree=k) if f.is_complete()]
    assert _walk_list(iet.rauzy_closed_walks(C, N, left_degree=k)) == _walk_list(filtered)


def test_h2_visits():
    # counts of the labelled Rauzy diagram iterator over the closed loops
    # through the representative, summed over the complete loops of each length
    C = Stratum([2]).unique_component()
    p = C.permutation_representative(reduced=False)
    pr = p.reduced()
    visits = []
    for n in range(8, 13):
        total = 0
        for f in iet.rauzy_closed_walks(C, n, min_length=n):
            total += sum(1 for q, winner, side in f if q.reduced() == pr)
        visits.append(total)
    assert visits == [2, 12, 28, 76, 204]


def _positive_power(M):
    d = M.nrows()
    return all(x > 0 for x in (M ** ((d - 1)**2 + 1)).list())


def test_walk_properties():
    C = Stratum([2, 2]).odd_component()
    for f in iet.rauzy_closed_walks(C, 13):
        d = len(f.start())
        assert f.is_closed()
        assert f.is_complete()
        assert f.start()._labels[0] == list(range(d))
        assert _positive_power(f.matrix())
    for f in iet.rauzy_closed_walks(C, 8, complete=False):
        assert f.is_complete() == _positive_power(f.matrix())


# name -> number of vertices of the reduced Rauzy diagram
KEY_MOVE_VERTICES = {
    'H2': 7,
    'H4hyp': 31,
    'H4odd': 134,
    'H22hyp': 63,
    'H22odd': 294,
    'H31left3': 509,
    'H31left1': 261,
}


@pytest.mark.parametrize('name', list(KEY_MOVE_VERTICES))
def test_key_moves(name):
    # the moves on keys against the moves of the canonical labellings, on
    # every edge of the diagram
    from surface_dynamics.interval_exchanges.rauzy_walks import ReducedRauzyDiagram
    C, k = diagram(name)
    G = ReducedRauzyDiagram(C.permutation_representative(left_degree=k, reduced=False))
    V = G.vertices()
    assert len(V) == KEY_MOVE_VERTICES[name]
    for key in V:
        q = G._canonical(key)
        for move in (0, 1):
            moved = G._key_move(key, move)
            assert q.has_rauzy_move(move)
            assert moved is not None
            successor, winner, cycle = moved
            r = q.rauzy_move(move)
            assert successor == G._key(r)
            assert winner == q._labels[move][-1]
            assert _expand(cycle, G._d) == G._move_map(r)


def _expand(cycle, d):
    # the relabelling list of a cycle (s, t): m[t] = s, and the ranks from s
    # to t shift by one towards t; None is the identity
    m = list(range(d))
    if cycle is not None:
        s, t = cycle
        assert s != t
        m[t] = s
        if t < s:
            m[t + 1:s + 1] = range(t, s)
        else:
            m[s:t] = range(s + 1, t + 1)
    return m


def _rotations(d, i, j):
    # the two rotations of the ranks i < j: j goes to i, and i goes to j
    head, tail = list(range(i)), list(range(j + 1, d))
    return (head + [j] + list(range(i, j)) + tail,
            head + list(range(i + 1, j + 1)) + [i] + tail)


def _edge_relabelling(G, u, move, lab):
    # the relabelling stored by G for the edge applied to the tuple lab
    getter = (G._top_getters, G._bot_getters)[move][u]
    return lab if getter is None else getter(lab)


@pytest.mark.parametrize('name', list(KEY_MOVE_VERTICES))
def test_edge_relabellings(name):
    # on every edge, the relabelling is the identity or one rotation of an
    # interval of ranks, the identity for every top move, and the data
    # stored for the edge applies it; the getters are shared, one per
    # distinct cycle
    from surface_dynamics.interval_exchanges.rauzy_walks import ReducedRauzyDiagram
    C, k = diagram(name)
    G = ReducedRauzyDiagram(C.permutation_representative(left_degree=k, reduced=False))
    d = G._d
    lab = tuple(range(10, 10 + d))
    cycles = set()
    for u, key in enumerate(G.vertices()):
        for move in (0, 1):
            cycle = G._key_move(key, move)[2]
            m = _expand(cycle, d)
            moved = [b for b in range(d) if m[b] != b]
            if moved:
                assert move == 1
                assert m in _rotations(d, moved[0], moved[-1])
                cycles.add(cycle)
            assert _edge_relabelling(G, u, move, lab) == tuple(lab[b] for b in m)
    getters = {id(g) for g in G._top_getters + G._bot_getters if g is not None}
    assert len(getters) == len(cycles)


def test_arguments():
    C = Stratum([2]).unique_component()
    p = C.permutation_representative()
    with pytest.raises(ValueError):
        list(iet.rauzy_closed_walks(p, 8, left_degree=1))
    assert (sorted(str(f) for f in iet.rauzy_closed_walks(p, 8)) ==
            sorted(str(f) for f in iet.rauzy_closed_walks(C, 8)))

    # a reduced permutation gives the same walks as the labelled one
    q = iet.Permutation('a b c d', 'd c b a')
    r = iet.Permutation('a b c d', 'd c b a', reduced=True)
    assert (sorted(str(f) for f in iet.rauzy_closed_walks(r, 8)) ==
            sorted(str(f) for f in iet.rauzy_closed_walks(q, 8)))

    # the walks keep the alphabet of the permutation, and their words of moves
    # are those of the stratum component
    walks = list(iet.rauzy_closed_walks(q, 9))
    assert walks
    for f in walks:
        assert list(f.start().alphabet()) == ['a', 'b', 'c', 'd']
    assert ([f._simplified_flip_sequence() for f in walks] ==
            [f._simplified_flip_sequence() for f in iet.rauzy_closed_walks(C, 9)])
