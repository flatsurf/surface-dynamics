r"""
Closed walks of reduced Rauzy diagrams

The reduced Rauzy diagram of an Abelian permutation has one vertex per reduced
permutation of its Rauzy class and two edges out of each vertex, the right
Rauzy moves ``t`` (top winner) and ``b`` (bottom winner). A *closed walk* is a
path of this diagram whose end is its start, taken up to rotation. It is a
*primitive closed walk* when it is not a proper power of a shorter closed walk
(read on the sequence of pairs (vertex, move), not on the word of moves
alone). It is *complete* when every letter is the winner of some move, up to
the orbits of the relabelling that closes it.

A reduced vertex is represented by its *canonical labelling*: the labelled
permutation, relabelled so that its top row is `0 1 \dots d-1`. Its *key* is
the tuple of the bottom row of its canonical labelling, and vertices are
compared by key.

EXAMPLES::

    sage: from surface_dynamics import Stratum
    sage: from surface_dynamics.interval_exchanges.rauzy_walks import rauzy_closed_walks
    sage: C = Stratum([2]).unique_component()
    sage: for f in rauzy_closed_walks(C, 9):
    ....:     print(f.start().str(sep=' / '), f._simplified_flip_sequence())
    0 1 2 3 / 1 3 0 2 tbtbtttbt
    0 1 2 3 / 1 3 0 2 tbbtttbtb
    0 1 2 3 / 1 3 2 0 bttbttbtb
    0 1 2 3 / 1 3 2 0 btbttbtb
    0 1 2 3 / 1 3 2 0 btbttbbtb
    0 1 2 3 / 1 3 2 0 bbtbtbttb
    0 1 2 3 / 1 3 2 0 bbtbbtbtb
"""
# ***************************************************************************
#       Copyright (C) 2026 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

import operator

from .template import OrientablePermutationIET, FlippedPermutationIET
from .labelled import LabelledPermutationIET
from .flip_sequence import IETFlipSequence


class LazyReducedRauzyDiagram:
    r"""
    The reduced Rauzy diagram of a labelled Abelian permutation, built on
    demand.

    The diagram holds, for each reduced vertex met so far, its key, its
    canonical labelling and its successors by the right Rauzy moves already
    requested. It grows by one entry per reduced vertex visited, up to the
    size of the reduced Rauzy class (7 vertices for `H(2)`, 509 for the
    representative of left degree 3 of `H_3(3, 1)`). It is the graph a search
    walks through, not a store of results.

    INPUT:

    - ``p`` -- an irreducible labelled permutation of interval exchange
      transformation (not flipped)

    EXAMPLES::

        sage: from surface_dynamics import iet
        sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
        sage: p = iet.Permutation('a b c d', 'd c b a')
        sage: G = LazyReducedRauzyDiagram(p)
        sage: G
        Lazy reduced Rauzy diagram on 4 letters (1 vertex held)
        sage: G.start()
        (3, 2, 1, 0)
        sage: len(G.vertices())
        7
        sage: len(G)
        7
    """
    def __init__(self, p):
        r"""
        TESTS::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: LazyReducedRauzyDiagram(iet.Permutation('a b c', 'c b a', reduced=True))
            Traceback (most recent call last):
            ...
            ValueError: p must be a labelled permutation of interval exchange transformation
            sage: LazyReducedRauzyDiagram(iet.Permutation('a b c', 'a c b'))
            Traceback (most recent call last):
            ...
            ValueError: p must be irreducible
        """
        if not isinstance(p, LabelledPermutationIET) or isinstance(p, FlippedPermutationIET):
            raise ValueError('p must be a labelled permutation of interval exchange transformation')
        if not p.is_irreducible():
            raise ValueError('p must be irreducible')
        self._d = len(p)
        self._start = self._key(p)
        self.clear()

    @staticmethod
    def _key(q):
        r"""
        Return the key of the labelled permutation ``q``: the bottom row of
        ``q`` relabelled so that its top row is ``0, 1, ..., d-1``.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: LazyReducedRauzyDiagram._key(iet.Permutation('a c b d', 'd a b c'))
            (3, 0, 2, 1)
        """
        top, bot = q._labels
        r = [None] * len(top)
        for i, a in enumerate(top):
            r[a] = i
        return tuple(r[b] for b in bot)

    def _canonical(self, key):
        r"""
        Return the canonical labelling of the vertex ``key``, as a new
        permutation.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
            sage: G._canonical((2, 0, 1))
            0 1 2
            2 0 1
        """
        from .constructors import Permutation
        return Permutation(list(range(self._d)), list(key), reduced=False)

    def clear(self):
        r"""
        Forget every vertex but the start.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c d', 'd c b a'))
            sage: len(G.vertices())
            7
            sage: G.clear()
            sage: len(G)
            1
        """
        self._labellings = {self._start: self._canonical(self._start)}
        self._successors = {self._start: [None, None]}
        self._explored = False

    def __len__(self):
        r"""
        Return the number of vertices held.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c d', 'd c b a'))
            sage: len(G)
            1
            sage: G.successor(G.start(), 0)
            (3, 0, 2, 1)
            sage: len(G)
            2
        """
        return len(self._labellings)

    def __repr__(self):
        r"""
        TESTS::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c d', 'd c b a'))
            sage: _ = G.vertices()
            sage: G
            Lazy reduced Rauzy diagram on 4 letters (7 vertices held)
        """
        n = len(self)
        return 'Lazy reduced Rauzy diagram on %d letters (%d %s held)' % (self._d, n, 'vertex' if n == 1 else 'vertices')

    def start(self):
        r"""
        Return the key of the vertex of the permutation the diagram was built
        from.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: G = LazyReducedRauzyDiagram(iet.Permutation('b a c', 'c a b'))
            sage: G.start()
            (2, 1, 0)
        """
        return self._start

    def canonical_labelling(self, key):
        r"""
        Return a copy of the canonical labelling of the vertex ``key``.

        The vertex must be held, that is, be the start or have been reached by
        :meth:`successor` or :meth:`vertices`.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c d', 'd c b a'))
            sage: G.canonical_labelling(G.successor(G.start(), 1))
            0 1 2 3
            1 3 2 0
            sage: G.canonical_labelling((1, 0, 2, 3))
            Traceback (most recent call last):
            ...
            KeyError: (1, 0, 2, 3)
        """
        return self._labellings[key].__copy__()

    def successor(self, key, move):
        r"""
        Return the key of the vertex reached from the vertex ``key`` by the
        right Rauzy move ``move`` (``0`` for top, ``1`` for bottom).

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c d', 'd c b a'))
            sage: k = G.start()
            sage: G.successor(k, 0)
            (3, 0, 2, 1)
            sage: G.successor(k, 1)
            (1, 3, 2, 0)
            sage: G.successor(G.successor(k, 1), 1)
            (2, 1, 3, 0)
        """
        succ = self._successors[key]
        s = succ[move]
        if s is None:
            s = self._key(self._labellings[key].rauzy_move(move))
            if s not in self._labellings:
                self._labellings[s] = self._canonical(s)
                self._successors[s] = [None, None]
            succ[move] = s
        return s

    def vertices(self):
        r"""
        Return the sorted list of the keys of the vertices of the reduced
        Rauzy class.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c d', 'd c b a'))
            sage: G.vertices()
            [(1, 3, 0, 2),
             (1, 3, 2, 0),
             (2, 0, 3, 1),
             (2, 1, 3, 0),
             (3, 0, 2, 1),
             (3, 1, 0, 2),
             (3, 2, 1, 0)]
        """
        if not self._explored:
            seen = {self._start}
            todo = [self._start]
            while todo:
                k = todo.pop()
                for move in (0, 1):
                    s = self.successor(k, move)
                    if s not in seen:
                        seen.add(s)
                        todo.append(s)
            self._explored = True
        return sorted(self._labellings)


def _labelled_permutation(x, left_degree):
    r"""
    Return the labelled permutation whose reduced Rauzy diagram is the one of
    ``x``, an Abelian permutation or a component of an Abelian stratum.

    EXAMPLES::

        sage: from surface_dynamics import iet, Stratum
        sage: from surface_dynamics.interval_exchanges.rauzy_walks import _labelled_permutation
        sage: C = Stratum([2]).unique_component()
        sage: _labelled_permutation(C, None)
        0 1 2 3
        3 2 1 0
        sage: _labelled_permutation(C.permutation_representative(), None)
        0 1 2 3
        3 2 1 0
        sage: _labelled_permutation(iet.Permutation('a b', 'b a', flips='a'), None)
        Traceback (most recent call last):
        ...
        ValueError: x must be an Abelian permutation or a component of an Abelian stratum
    """
    from surface_dynamics.flat_surfaces.abelian_strata import AbelianStratumComponent
    from .constructors import Permutation

    if isinstance(x, AbelianStratumComponent):
        return x.permutation_representative(left_degree=left_degree, reduced=False)
    if isinstance(x, OrientablePermutationIET) and not isinstance(x, FlippedPermutationIET):
        if left_degree is not None:
            raise ValueError('left_degree must be None when x is a permutation')
        if not x.is_irreducible():
            raise ValueError('x must be irreducible')
        if isinstance(x, LabelledPermutationIET):
            return x
        return Permutation(x.list(), reduced=False)
    raise ValueError('x must be an Abelian permutation or a component of an Abelian stratum')


def rauzy_closed_walks(x, max_length, min_length=1, complete=True, primitive=True, left_degree=None):
    r"""
    Iterate over the closed walks of the reduced Rauzy diagram of ``x`` with
    length between ``min_length`` and ``max_length``.

    INPUT:

    - ``x`` -- an irreducible Abelian permutation (labelled or reduced), or a
      component of an Abelian stratum

    - ``max_length``, ``min_length`` -- (``min_length`` defaults to ``1``)
      non-negative integers, the bounds on the length of the walks

    - ``complete`` -- (default ``True``) whether to keep only the complete
      walks, those in which every letter is a winner up to the orbits of the
      relabelling that closes the walk

    - ``primitive`` -- (default ``True``) whether to keep only the primitive
      closed walks; if ``False``, the proper powers of closed walks are
      yielded as well

    - ``left_degree`` -- (optional) when ``x`` is a stratum component, the
      left degree of the permutation representative whose diagram is walked
      (see ``permutation_representative``); must be ``None`` when ``x`` is a
      permutation

    OUTPUT:

    An iterator over closed :class:`IETFlipSequence`, each walk once up to
    rotation. A walk starts at the canonical labelling (top row ``0 1 ...
    d-1``) of its least vertex, and its word of moves is the least among the
    rotations starting there (``t < b``); it is closed by
    :meth:`IETFlipSequence.close`, so it carries the relabelling that maps
    its end back to its start. The walks come by increasing least vertex,
    then by increasing sequence of first-return loops at this vertex.

    ALGORITHM:

    For each vertex `v`, the first-return loops at `v` (the words from `v`
    back to `v` whose intermediate vertices are all greater than `v`) form a
    prefix code. The walks whose least vertex is `v`, in their canonical
    rotation, are the sequences of these loops that are necklaces over the
    loops ordered as words, and the primitive closed walks are the Lyndon
    words. They are generated by the prenecklace backtracking of
    Fredricksen, Kessler and Maiorana, with the total length bounded by
    ``max_length``; the loops themselves are found by a depth-first search
    pruned by the distance back to `v`.

    EXAMPLES:

    The only complete closed walk of length `8` of `H(2)`::

        sage: from surface_dynamics import Stratum, iet
        sage: from surface_dynamics.interval_exchanges.rauzy_walks import rauzy_closed_walks
        sage: C = Stratum([2]).unique_component()
        sage: walks = list(rauzy_closed_walks(C, 8, min_length=8))
        sage: len(walks)
        1
        sage: f = walks[0]
        sage: f.start()
        0 1 2 3
        1 3 2 0
        sage: f._simplified_flip_sequence()
        'btbttbtb'
        sage: f.is_complete()
        True
        sage: f.matrix().charpoly()
        x^4 - 7*x^3 + 13*x^2 - 7*x + 1

    All the primitive closed walks of `H(2)` of length at most `3`::

        sage: walks = list(rauzy_closed_walks(C, 3, complete=False))
        sage: [sum(1 for f in walks if len(f) == n) for n in (1, 2, 3)]
        [4, 2, 4]

    The same walks from a permutation of the Rauzy class::

        sage: p = iet.Permutation('a b c d', 'd c b a')
        sage: sorted(str(f) for f in rauzy_closed_walks(p, 3, complete=False)) == sorted(str(f) for f in walks)
        True

    The left degree only applies to a stratum component::

        sage: p = C.permutation_representative()
        sage: rauzy_closed_walks(p, 8, left_degree=1)
        Traceback (most recent call last):
        ...
        ValueError: left_degree must be None when x is a permutation

    TESTS:

    With ``primitive=False``, the closed walks of length `n` number
    `\sum_{e | n} c_e` where `c_e` is the number of primitive closed walks
    of length `e`::

        sage: walks = list(rauzy_closed_walks(C, 4, complete=False, primitive=False))
        sage: [sum(1 for f in walks if len(f) == n) for n in (1, 2, 3, 4)]
        [4, 6, 8, 10]
        sage: all(f.is_closed() for f in walks)
        True

    Bounds::

        sage: list(rauzy_closed_walks(C, 0, complete=False))
        []
        sage: list(rauzy_closed_walks(C, 3, min_length=4, complete=False))
        []
        sage: len(list(rauzy_closed_walks(C, 3, min_length=0, complete=False)))
        10
        sage: rauzy_closed_walks(C, -1)
        Traceback (most recent call last):
        ...
        ValueError: max_length and min_length must be non-negative
        sage: rauzy_closed_walks(iet.Permutation('a b c', 'a c b'), 3)
        Traceback (most recent call last):
        ...
        ValueError: x must be irreducible
        sage: rauzy_closed_walks(iet.GeneralizedPermutation('a a b', 'b c c'), 3)
        Traceback (most recent call last):
        ...
        ValueError: x must be an Abelian permutation or a component of an Abelian stratum
    """
    max_length = operator.index(max_length)
    min_length = operator.index(min_length)
    if max_length < 0 or min_length < 0:
        raise ValueError('max_length and min_length must be non-negative')
    G = LazyReducedRauzyDiagram(_labelled_permutation(x, left_degree))
    return _closed_walks(G, max_length, max(min_length, 1), bool(complete), bool(primitive))


def _distances_to(v, predecessors, max_length):
    r"""
    Return the dictionary ``u -> d`` where ``d`` is the length of a shortest
    path from ``u`` to ``v`` whose vertices but ``v`` are all greater than
    ``v``, for the vertices ``u > v`` at distance less than ``max_length``.

    EXAMPLES::

        sage: from surface_dynamics.interval_exchanges.rauzy_walks import _distances_to
        sage: pred = {0: [1, 2], 1: [0, 2], 2: [1]}
        sage: _distances_to(0, pred, 5)
        {1: 1, 2: 1}
        sage: _distances_to(1, pred, 5)
        {2: 1}
        sage: _distances_to(1, {0: [1], 1: [2], 2: [3], 3: [0]}, 5)
        {2: 1, 3: 2}
    """
    dist = {}
    frontier = [v]
    level = 0
    while frontier and level < max_length:
        level += 1
        new = []
        for w in frontier:
            for u in predecessors[w]:
                if u > v and u not in dist:
                    dist[u] = level
                    new.append(u)
        frontier = new
    return dist


def _first_return_loops(G, v, dist, max_length):
    r"""
    Return the list, sorted as words, of the first-return loops at ``v`` of
    length at most ``max_length``: the words of moves from ``v`` back to
    ``v`` whose intermediate vertices are all greater than ``v``.

    ``dist`` is the output of :func:`_distances_to`.

    EXAMPLES::

        sage: from surface_dynamics import iet
        sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram, _distances_to, _first_return_loops
        sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
        sage: V = G.vertices()
        sage: pred = {k: [] for k in V}
        sage: for k in V:
        ....:     for move in (0, 1):
        ....:         pred[G.successor(k, move)].append(k)
        sage: v = V[0]
        sage: dist = _distances_to(v, pred, 4)
        sage: _first_return_loops(G, v, dist, 4)
        [(0,), (1, 0, 0, 1), (1, 1)]
    """
    loops = []
    word = []

    def search(u, length):
        # u is reached from v after length moves; u is v itself or greater
        # than v, and the pruning ensures that some path of length at most
        # max_length - length leads from u back to v
        for move in (0, 1):
            w = G.successor(u, move)
            word.append(move)
            if w == v:
                loops.append(tuple(word))
            elif w > v and w in dist and length + 1 + dist[w] <= max_length:
                search(w, length + 1)
            word.pop()

    if max_length >= 1:
        search(v, 0)
    return loops


def _loop_necklaces(lengths, max_length, min_length, primitive):
    r"""
    Iterate over the necklaces (Lyndon words if ``primitive``) over the
    alphabet ``0, 1, ..., len(lengths) - 1`` whose weight is between
    ``min_length`` and ``max_length``, the letter ``i`` having weight
    ``lengths[i]``.

    The words are generated by the prenecklace backtracking of Fredricksen,
    Kessler and Maiorana, the subtrees of weight greater than ``max_length``
    being cut. The list yielded is modified afterwards: copy it to keep it.

    EXAMPLES::

        sage: from surface_dynamics.interval_exchanges.rauzy_walks import _loop_necklaces
        sage: [tuple(w) for w in _loop_necklaces([1, 1], 4, 1, True)]
        [(0,), (0, 0, 0, 1), (0, 0, 1), (0, 0, 1, 1), (0, 1), (0, 1, 1), (0, 1, 1, 1), (1,)]
        sage: [tuple(w) for w in _loop_necklaces([1, 1], 4, 4, False)]
        [(0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 1, 1), (0, 1, 0, 1), (0, 1, 1, 1), (1, 1, 1, 1)]
        sage: [tuple(w) for w in _loop_necklaces([2, 1], 4, 1, True)]
        [(0,), (0, 1), (0, 1, 1), (1,)]
    """
    n = len(lengths)
    seq = []

    def extend(p, total):
        # seq is a prenecklace of length m with period p (the length of its
        # longest Lyndon prefix)
        m = len(seq)
        if m and total >= min_length and (p == m if primitive else m % p == 0):
            yield seq
        first = seq[m - p] if m else 0
        for i in range(first, n):
            t = total + lengths[i]
            if t > max_length:
                # the letters are not sorted by weight: a later one may fit
                continue
            seq.append(i)
            yield from extend(p if m and i == first else m + 1, t)
            seq.pop()

    yield from extend(0, 0)


def _closed_walks(G, max_length, min_length, complete, primitive):
    r"""
    Iterate over the closed walks of the reduced Rauzy diagram ``G``, as
    described in :func:`rauzy_closed_walks`.

    EXAMPLES::

        sage: from surface_dynamics import iet
        sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram, _closed_walks
        sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
        sage: for f in _closed_walks(G, 3, 1, False, True):
        ....:     print(f._simplified_flip_sequence())
        t
        tbb
        bb
        tt
        ttb
        b
    """
    vertices = G.vertices()
    predecessors = {k: [] for k in vertices}
    for k in vertices:
        for move in (0, 1):
            predecessors[G.successor(k, move)].append(k)

    for v in vertices:
        dist = _distances_to(v, predecessors, max_length)
        loops = _first_return_loops(G, v, dist, max_length)
        lengths = [len(loop) for loop in loops]
        for seq in _loop_necklaces(lengths, max_length, min_length, primitive):
            word = [move for i in seq for move in loops[i]]
            f = IETFlipSequence(G.canonical_labelling(v), word)
            f.close()
            # completeness is checked on the walk built, replaying its moves
            if complete and not f.is_complete():
                continue
            yield f
