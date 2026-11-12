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


# name -> (vertices of the reduced Rauzy diagram of the labelled
# representative, undefined moves)
KEY_MOVE_DIAGRAMS = {
    'Q1111hyp': (380, 112),
    'Q22hyp': (73, 26),
    'Q5m1': (440, 150),
}


@pytest.mark.parametrize('name', list(KEY_MOVE_DIAGRAMS))
def test_key_moves(name):
    # the moves on keys against the moves of the canonical labellings, on
    # every edge of the diagram and every undefined move
    from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
    n, undefined = KEY_MOVE_DIAGRAMS[name]
    G = LazyReducedRauzyDiagram(diagram(name).permutation_representative(reduced=False))
    V = G.vertices()
    assert len(V) == n
    seen = 0
    for key in V:
        q = G._canonical(key)
        for move in (0, 1):
            moved = G._key_move(key, move)
            assert (moved is None) == (not q.has_rauzy_move(move))
            if moved is None:
                seen += 1
                continue
            successor, winner, m = moved
            r = q.rauzy_move(move)
            assert successor == G._key(r)
            assert winner == q._labels[move][-1]
            assert m == G._move_map(r)
    assert seen == undefined


def _rotations(d, i, j):
    # the two rotations of the ranks i < j: j goes to i, and i goes to j
    head, tail = list(range(i)), list(range(j + 1, d))
    return (head + [j] + list(range(i, j)) + tail,
            head + list(range(i + 1, j + 1)) + [i] + tail)


def _edge_relabelling(G, u, move, lab):
    # the relabelling stored by G for the edge applied to the tuple lab
    getter = (G._top_getters, G._bot_getters)[move][u]
    return lab if getter is None else getter(lab)


@pytest.mark.parametrize('name', list(KEY_MOVE_DIAGRAMS))
def test_edge_relabellings(name):
    # on every edge and every undefined move, the relabelling is the
    # identity or one rotation of an interval of ranks, and the data stored
    # for the edge applies it
    from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
    G = LazyReducedRauzyDiagram(diagram(name).permutation_representative(reduced=False))
    d = G._d
    lab = tuple(range(10, 10 + d))
    rotations = 0
    for u, key in enumerate(G.vertices()):
        for move in (0, 1):
            moved = G._key_move(key, move)
            m = list(range(d)) if moved is None else moved[2]
            changed = [b for b in range(d) if m[b] != b]
            if changed:
                rotations += 1
                assert m in _rotations(d, changed[0], changed[-1])
            assert _edge_relabelling(G, u, move, lab) == tuple(lab[b] for b in m)
    assert rotations > 0


ALL_COUNTS = {
    # brute force over all words of defined moves from every reduced vertex,
    # agreeing with the trace formula (1/n) sum_{e|n} mu(n/e) tr(A^e) on the
    # adjacency matrix A of the reduced diagram
    'Q1111hyp': [76, 34, 34, 32, 54, 58, 104, 134, 212, 302, 480, 706],
    'Q22hyp': [12, 6, 2, 6, 4, 10, 8, 21, 22, 54, 52, 123, 130, 293, 342, 709],
    'Q5m1': [50, 24, 30, 34, 26, 40, 38, 72, 88, 144, 186, 322, 404, 713],
}


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
@pytest.mark.parametrize('name', list(COMPLETE_COUNTS_SLOW))
def test_counts_complete_slow(name):
    C = diagram(name)
    expected = COMPLETE_COUNTS_SLOW[name]
    N = len(expected)
    assert counts(iet.rauzy_closed_walks(C, N), N) == expected


def _walk_list(walks):
    return [(str(f.start()), f._simplified_flip_sequence()) for f in walks]


def test_complete_against_filter():
    # Q22hyp up to length 13: 450 walks of complete=False (the sum of
    # ALL_COUNTS up to 13), of which 91 are complete (the sum of
    # COMPLETE_COUNTS); the relabelling group is transitive on the letters
    C = diagram('Q22hyp')
    walks = list(iet.rauzy_closed_walks(C, 13, complete=False))
    assert len(walks) == sum(ALL_COUNTS['Q22hyp'][:13]) == 450
    filtered = [f for f in walks if f.is_complete()]
    assert len(filtered) == sum(COMPLETE_COUNTS['Q22hyp']) == 91
    assert _walk_list(iet.rauzy_closed_walks(C, 13)) == _walk_list(filtered)


def _first_appearance(q):
    r"""The letters of q (indices in its alphabet) by first appearance, top row then bottom row."""
    res = []
    for row in q._labels:
        for a in row:
            if a not in res:
                res.append(a)
    return res


def _is_proper_power(f):
    r"""
    Whether the closed walk f is a proper power of a shorter closed walk: its
    word of moves is a power of a shorter word that returns to the reduced
    vertex of its start.
    """
    from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
    from surface_dynamics.interval_exchanges.flip_sequence import RauzyFlipSequence
    word = f._simplified_flip_sequence()
    n = len(word)
    start = LazyReducedRauzyDiagram._key(f.start())
    for p in range(1, n):
        if n % p == 0 and word == word[:p] * (n // p):
            g = RauzyFlipSequence(f.start(), word[:p])
            if LazyReducedRauzyDiagram._key(g.end()) == start:
                return True
    return False


def _check_walks(walks, complete):
    for f in walks:
        start = f.start()
        assert f.is_closed()
        assert f._check() is None
        if complete:
            assert f.is_complete()
        assert not _is_proper_power(f)
        assert _first_appearance(start) == list(range(len(start)))


def test_walk_properties():
    _check_walks(iet.rauzy_closed_walks(diagram('Q22hyp'), 13), True)

    C = diagram('Q5m1')
    N = 8
    walks = list(iet.rauzy_closed_walks(C, N, complete=False))
    assert counts(walks, N) == ALL_COUNTS['Q5m1'][:N]
    _check_walks(walks, False)
    # with the proper powers, the closed walks of length n number the sum of
    # the primitive ones of length e over the divisors e of n
    c = ALL_COUNTS['Q5m1']
    expected = [sum(c[e - 1] for e in range(1, n + 1) if n % e == 0) for n in range(1, N + 1)]
    powers = list(iet.rauzy_closed_walks(C, N, complete=False, primitive=False))
    assert counts(powers, N) == expected
    assert sum(_is_proper_power(f) for f in powers) == sum(expected) - sum(c[:N])


def test_relabelling_orbits():
    # the orbits read off the tree labellings against those of the group
    # generated by the conjugated generators, at every vertex
    from sage.all import PermutationGroup
    from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
    G = LazyReducedRauzyDiagram(diagram('Q5m1').permutation_representative(reduced=False))
    V = G.vertices()
    assert len(V) == 440
    d = G._d
    for key in V:
        gens = [[i + 1 for i in g] for g in G.relabelling_generators(key)]
        orbits = PermutationGroup(gens, domain=range(1, d + 1)).orbits()
        assert G.relabelling_orbits(key) == sorted(tuple(sorted(i - 1 for i in o)) for o in orbits)


def test_arguments():
    C = diagram('Q22hyp')
    with pytest.raises(NotImplementedError):
        list(iet.rauzy_closed_walks(C, 8, left_degree=2))
    p = C.permutation_representative(reduced=False)
    with pytest.raises(ValueError):
        list(iet.rauzy_closed_walks(p, 8, left_degree=2))

    # the labelled and reduced representatives give the walks of the component
    r = C.permutation_representative()
    walks = [str(f) for f in iet.rauzy_closed_walks(C, 10)]
    assert walks
    assert [str(f) for f in iet.rauzy_closed_walks(p, 10)] == walks
    assert [str(f) for f in iet.rauzy_closed_walks(r, 10)] == walks

    # a linear involution over letters keeps its alphabet, and its words of
    # moves are those of the component
    L = 'abcde'
    top, bot = p.list()
    q = iet.GeneralizedPermutation([L[i] for i in top], [L[i] for i in bot], alphabet=L)
    walks_q = list(iet.rauzy_closed_walks(q, 10))
    assert walks_q
    for f in walks_q:
        assert list(f.start().alphabet()) == list(L)
    assert ([f._simplified_flip_sequence() for f in walks_q] ==
            [f._simplified_flip_sequence() for f in iet.rauzy_closed_walks(C, 10)])
