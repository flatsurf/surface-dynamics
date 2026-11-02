r"""
Closed walks of reduced Rauzy diagrams

The reduced Rauzy diagram of an Abelian permutation has one vertex per
permutation of its Rauzy class up to relabelling. Two edges reach out
of each vertex, the right Rauzy moves ``t`` (top winner) and ``b`` (bottom
winner). A *closed walk* is a path in this diagram whose end is its start,
taken up to relabelling. It is a *primitive closed walk* when it is not a proper
power of a shorter closed walk (read on the sequence of pairs (vertex, move),
not on the word of moves alone). It is *complete* when every letter is the
winner of some move, up to the orbits of the relabelling that closes it.

A reduced vertex is represented by its *canonical labelling*: the labelled
permutation, relabelled so that its top row is `0 1 \dots d-1`. Its *key* is
the tuple of the bottom row of its canonical labelling, and vertices are
compared by key.

EXAMPLES::

    sage: from surface_dynamics import Stratum, iet
    sage: C = Stratum([2]).unique_component()
    sage: for f in iet.rauzy_closed_walks(C, 9):
    ....:     print(f)
    FlipSequence('0 1 2 3\n1 3 0 2', 'tbtbtttbt')
    FlipSequence('0 1 2 3\n1 3 0 2', 'tbbtttbtb')
    FlipSequence('0 1 2 3\n1 3 2 0', 'bttbttbtb')
    FlipSequence('0 1 2 3\n1 3 2 0', 'btbttbtb')
    FlipSequence('0 1 2 3\n1 3 2 0', 'btbttbbtb')
    FlipSequence('0 1 2 3\n1 3 2 0', 'bbtbtbttb')
    FlipSequence('0 1 2 3\n1 3 2 0', 'bbtbbtbtb')
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

    Once :meth:`relabelling_generators` or :meth:`relabelling_orbits` has
    been called, the diagram also holds the relabelling group of the class:
    one labelling per reduced vertex (along a spanning tree from the start)
    and one list of generators at the start, at most one per edge of the
    diagram. Both are computed once for the whole class.

    :meth:`clear` drops all of it but the start.

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
        Forget every vertex but the start, and the relabelling group.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c d', 'd c b a'))
            sage: len(G.vertices())
            7
            sage: G.relabelling_generators(G.start())
            []
            sage: G.clear()
            sage: len(G)
            1
            sage: G._tree_labellings is None
            True
        """
        self._labellings = {self._start: self._canonical(self._start)}
        self._successors = {self._start: [None, None]}
        self._explored = False
        self._tree_labellings = None
        self._root_generators = None

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

    def _relabelling_group(self):
        r"""
        Compute, once for the class, the labelling of each vertex along a
        spanning tree from the start and the generators of the relabelling
        group at the start.

        The tree labelling of a vertex ``u`` is stored as the top row ``t_u``
        of the labelled permutation reached from the canonical labelling of
        the start along the tree: that permutation is the canonical labelling
        of ``u`` with each letter ``a`` replaced by ``t_u[a]``. For each edge
        ``u -> w``, the move applied to the tree labelling of ``u`` gives the
        tree labelling of ``w`` relabelled by some ``g``; these Schreier
        generators ``g`` generate the group at the start.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
            sage: G._relabelling_group()
            sage: sorted(G._tree_labellings.items())
            [((1, 2, 0), (0, 2, 1)), ((2, 0, 1), (0, 1, 2)), ((2, 1, 0), (0, 1, 2))]
            sage: G._root_generators
            []
        """
        if self._root_generators is not None:
            return
        self.vertices()
        d = self._d
        identity = tuple(range(d))
        tree = {self._start: identity}
        generators = []
        seen = {identity}
        frontier = [self._start]
        while frontier:
            new = []
            for u in frontier:
                tu = tree[u]
                for move in (0, 1):
                    w = self.successor(u, move)
                    # the move from the canonical labelling of u ends at the
                    # canonical labelling of w relabelled by its top row
                    mu = self._labellings[u].rauzy_move(move)._labels[0]
                    x = tuple(tu[a] for a in mu)
                    if w not in tree:
                        tree[w] = x
                        new.append(w)
                        continue
                    y = tree[w]
                    g = [None] * d
                    for a, b in zip(y, x):
                        g[a] = b
                    g = tuple(g)
                    if g not in seen:
                        seen.add(g)
                        generators.append(g)
            frontier = new
        self._tree_labellings = tree
        self._root_generators = generators

    def relabelling_generators(self, key):
        r"""
        Return a list of generators of the relabelling group at the vertex
        ``key``.

        Each closed walk from the canonical labelling of ``key`` ends at this
        labelling up to a relabelling, a permutation of the letters; these
        relabellings form a group, the *relabelling group* at ``key``. The
        relabellings are tuples ``g`` of the integers ``0, ..., d-1``, the
        letter ``i`` being sent to ``g[i]``, as the relabelling of an
        :class:`IETFlipSequence`. The groups at two vertices are conjugate.

        The generators are the Schreier generators along a spanning tree,
        without the identity and without repetition.

        EXAMPLES:

        The relabelling group is trivial for `H(2)`, of order `3` for
        `H_3(2^2)^{odd}` and of order `2` for the representative of left
        degree `1` of `H_3(3, 1)`::

            sage: from surface_dynamics import Stratum
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: def order(G, k):
            ....:     return PermutationGroup([[i + 1 for i in g] for g in G.relabelling_generators(k)]).order()
            sage: G = LazyReducedRauzyDiagram(Stratum([2]).unique_component().permutation_representative(reduced=False))
            sage: G.relabelling_generators(G.start())
            []
            sage: G = LazyReducedRauzyDiagram(Stratum([2, 2]).odd_component().permutation_representative(reduced=False))
            sage: k = G.start()
            sage: order(G, k), order(G, G.vertices()[0])
            (3, 3)
            sage: C = Stratum([3, 1]).unique_component()
            sage: G = LazyReducedRauzyDiagram(C.permutation_representative(left_degree=1, reduced=False))
            sage: order(G, G.start()), order(G, G.vertices()[0])
            (2, 2)
            sage: all(sorted(g) == list(range(7)) for g in G.relabelling_generators(G.start()))
            True
        """
        self._relabelling_group()
        t = self._tree_labellings[key]
        tinv = [None] * self._d
        for a, b in enumerate(t):
            tinv[b] = a
        # conjugation by the tree labelling of key, a bijection: no identity
        # nor repetition appears
        return [tuple(tinv[g[t[a]]] for a in range(self._d)) for g in self._root_generators]

    def relabelling_orbits(self, key):
        r"""
        Return the orbits of the relabelling group at the vertex ``key`` on
        the letters, as a sorted list of sorted tuples.

        See :meth:`relabelling_generators`.

        EXAMPLES::

            sage: from surface_dynamics import Stratum
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: G = LazyReducedRauzyDiagram(Stratum([2]).unique_component().permutation_representative(reduced=False))
            sage: G.relabelling_orbits(G.start())
            [(0,), (1,), (2,), (3,)]
            sage: G = LazyReducedRauzyDiagram(Stratum([2, 2]).odd_component().permutation_representative(reduced=False))
            sage: G.relabelling_orbits(G.start())
            [(0,), (1, 5, 6), (2,), (3,), (4,)]
        """
        generators = self.relabelling_generators(key)
        orbit_of = [None] * self._d
        orbits = []
        for a in range(self._d):
            if orbit_of[a] is not None:
                continue
            orbit = [a]
            orbit_of[a] = a
            for b in orbit:
                for g in generators:
                    c = g[b]
                    if orbit_of[c] is None:
                        orbit_of[c] = a
                        orbit.append(c)
            orbits.append(tuple(sorted(orbit)))
        orbits.sort()
        return orbits


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

    With ``complete``, completeness is tracked along the backtracking. The
    relabelling that closes a walk at `v` lies in the relabelling group at
    `v` (see :meth:`LazyReducedRauzyDiagram.relabelling_generators`),
    computed once from Schreier generators along a spanning tree. The
    orbits of this group hit by the winners of a loop do not depend on the
    labelling the loop is run from, so a subtree is cut when more orbits
    contain no winner than moves are left. When the group is trivial, this
    is exact pruning by the letters still missing; otherwise the winners
    and the relabelling are combined loop by loop and completeness is
    decided exactly on each walk. A flip sequence is built for the walks
    yielded only.

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
        sage: f
        FlipSequence('0 1 2 3\n1 3 2 0', 'btbttbtb')
        sage: f.is_complete()
        True
        sage: f.matrix().charpoly()
        x^4 - 7*x^3 + 13*x^2 - 7*x + 1

    In genus `3`, complete closed walks are long: the components `H_3(4)^{hyp}`,
    `H_3(4)^{odd}`, `H_3(2^2)^{hyp}` and `H_3(3, 1)^c` have none of length at
    most `14`. The odd component of `H(2, 2)` has `8` of length at most `13`::

        sage: C3 = Stratum([2, 2]).odd_component()
        sage: walks3 = list(rauzy_closed_walks(C3, 13))
        sage: len(walks3)
        8
        sage: sorted(f.matrix().charpoly().number_of_real_roots() for f in walks3)
        [3, 3, 5, 5, 5, 5, 5, 5]

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


def _loop_data(G, v, loops, orbits, label_maps):
    r"""
    Return the winners of the ``loops`` at ``v``, and the orbits and label
    maps they give.

    Each loop is run once from the canonical labelling of ``v``. ``orbits``
    is a list of disjoint tuples of letters covering the letters, those of
    the relabelling group at ``v``.

    OUTPUT: a triple ``(winner_masks, orbit_masks, maps)``, lists indexed by
    the loops:

    - ``winner_masks[j]`` -- the bitmask of the letters of the canonical
      labelling of ``v`` that are winners of the loop ``j``

    - ``orbit_masks[j]`` -- the bitmask of the indices in ``orbits`` of the
      orbits containing a winner of the loop ``j``

    - ``maps[j]`` -- if ``label_maps`` is set, the top row of the labelled
      permutation at which the loop ``j`` ends, a tuple ``m`` such that this
      permutation is the canonical labelling of ``v`` with each letter ``a``
      replaced by ``m[a]``; otherwise ``maps`` is ``None``

    EXAMPLES::

        sage: from surface_dynamics import iet, Stratum
        sage: from surface_dynamics.interval_exchanges.flip_sequence import IETFlipSequence
        sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram, _distances_to, _first_return_loops, _loop_data
        sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c d', 'd c b a'))
        sage: V = G.vertices()
        sage: pred = {k: [] for k in V}
        sage: for k in V:
        ....:     for move in (0, 1):
        ....:         pred[G.successor(k, move)].append(k)
        sage: v = V[0]
        sage: loops = _first_return_loops(G, v, _distances_to(v, pred, 8), 8)
        sage: orbits = G.relabelling_orbits(v)
        sage: winner_masks, orbit_masks, maps = _loop_data(G, v, loops, orbits, True)
        sage: loops[0], bin(winner_masks[0]), maps[0]
        ((0, 0), '0b1000', (0, 1, 2, 3))
        sage: winner_masks == orbit_masks
        True

    The winners agree with those of the flip sequence of the loop::

        sage: all(winner_masks[j] == sum(1 << a for a in set(IETFlipSequence(G.canonical_labelling(v), loop).winners()))
        ....:     for j, loop in enumerate(loops))
        True
    """
    orbit_of = [None] * G._d
    for i, orbit in enumerate(orbits):
        for a in orbit:
            orbit_of[a] = i
    winner_masks = []
    orbit_masks = []
    maps = [] if label_maps else None
    for loop in loops:
        q = G.canonical_labelling(v)
        winners = 0
        for move in loop:
            winners |= 1 << q._labels[move][-1]
            q.rauzy_move(move, inplace=True)
        hit = 0
        a = 0
        w = winners
        while w:
            if w & 1:
                hit |= 1 << orbit_of[a]
            w >>= 1
            a += 1
        winner_masks.append(winners)
        orbit_masks.append(hit)
        if label_maps:
            maps.append(tuple(q._labels[0]))
    return winner_masks, orbit_masks, maps


def _loop_necklaces(lengths, max_length, min_length, primitive, orbit_masks=None, num_orbits=0, winner_masks=None, label_maps=None):
    r"""
    Iterate over the necklaces (Lyndon words if ``primitive``) over the
    alphabet ``0, 1, ..., len(lengths) - 1`` whose weight is between
    ``min_length`` and ``max_length``, the letter ``i`` having weight
    ``lengths[i]``.

    The words are generated by the prenecklace backtracking of Fredricksen,
    Kessler and Maiorana, the subtrees of weight greater than ``max_length``
    being cut. The list yielded is modified afterwards: copy it to keep it.

    If ``orbit_masks`` is given, only the complete words are yielded: the
    letter ``i`` is a loop whose winners hit the orbits of the bitmask
    ``orbit_masks[i]`` among the ``num_orbits`` orbits of the relabelling
    group, and a word is kept when the union of these bitmasks is full. A
    subtree is cut when more orbits are missing than moves are left, each
    move having one winner. If moreover ``label_maps`` is given (the
    relabelling group is nontrivial), the completeness of a word is decided
    exactly from ``winner_masks`` and ``label_maps`` (see
    :func:`_loop_data`): the winners of the letter ``j`` of the word are
    transported by the product of the label maps of the letters before it,
    and they must hit every orbit of the product of all the label maps.

    EXAMPLES::

        sage: from surface_dynamics.interval_exchanges.rauzy_walks import _loop_necklaces
        sage: [tuple(w) for w in _loop_necklaces([1, 1], 4, 1, True)]
        [(0,), (0, 0, 0, 1), (0, 0, 1), (0, 0, 1, 1), (0, 1), (0, 1, 1), (0, 1, 1, 1), (1,)]
        sage: [tuple(w) for w in _loop_necklaces([1, 1], 4, 4, False)]
        [(0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 1, 1), (0, 1, 0, 1), (0, 1, 1, 1), (1, 1, 1, 1)]
        sage: [tuple(w) for w in _loop_necklaces([2, 1], 4, 1, True)]
        [(0,), (0, 1), (0, 1, 1), (1,)]

    With two orbits, each loop hitting one of them::

        sage: [tuple(w) for w in _loop_necklaces([1, 1], 4, 1, True, [0b01, 0b10], 2)]
        [(0, 0, 0, 1), (0, 0, 1), (0, 0, 1, 1), (0, 1), (0, 1, 1), (0, 1, 1, 1)]

    A relabelling group of order `2` exchanging the letters ``0`` and
    ``1``. The first loop has winner ``1`` and ends with the letters
    exchanged, the second one has winner ``0`` and ends with the identity.
    The word ``(1,)`` hits the only orbit of the group but closes with the
    identity and misses the letter ``1``; the word ``(0, 0)``, not a
    Lyndon word, has winners ``1`` then ``0``::

        sage: [tuple(w) for w in _loop_necklaces([1, 1], 3, 1, True, [0b1, 0b1], 1, [0b10, 0b01], [(1, 0), (0, 1)])]
        [(0,), (0, 0, 1), (0, 1), (0, 1, 1)]
        sage: [tuple(w) for w in _loop_necklaces([1, 1], 2, 2, False, [0b1, 0b1], 1, [0b10, 0b01], [(1, 0), (0, 1)])]
        [(0, 0), (0, 1)]
    """
    n = len(lengths)
    if not n:
        return
    seq = []

    if orbit_masks is None:
        def extend(p, total):
            # seq is a prenecklace of length m with period p (the length of
            # its longest Lyndon prefix)
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
        return

    # Python integers for int.bit_count, whatever the input (Sage Integer)
    orbit_masks = [int(h) for h in orbit_masks]
    full = (1 << int(num_orbits)) - 1

    if label_maps is None:
        # trivial relabelling group: the orbits are the letters and every
        # walk closes with the identity, so a full mask is completeness
        def extend(p, total, hit):
            m = len(seq)
            if m and hit == full and total >= min_length and (p == m if primitive else m % p == 0):
                yield seq
            first = seq[m - p] if m else 0
            for i in range(first, n):
                t = total + lengths[i]
                if t > max_length:
                    continue
                h = hit | orbit_masks[i]
                if (full ^ h).bit_count() > max_length - t:
                    continue
                seq.append(i)
                yield from extend(p if m and i == first else m + 1, t, h)
                seq.pop()

        yield from extend(0, 0, 0)
        return

    d = len(label_maps[0])
    winner_letters = [[a for a in range(d) if (w >> a) & 1] for w in winner_masks]

    def is_complete(winners, relabelling):
        # every cycle of the relabelling contains a winner
        seen = 0
        for a in range(d):
            if (seen >> a) & 1:
                continue
            cycle = 0
            b = a
            while not (cycle >> b) & 1:
                cycle |= 1 << b
                b = relabelling[b]
            if not cycle & winners:
                return False
            seen |= cycle
        return True

    def extend(p, total, hit, winners, relabelling):
        # winners: the letters of the canonical labelling of v that are
        # winners of the walk seq run from this labelling; relabelling: the
        # top row of the labelled permutation where it ends
        m = len(seq)
        if (m and hit == full and total >= min_length and (p == m if primitive else m % p == 0)
                and is_complete(winners, relabelling)):
            yield seq
        first = seq[m - p] if m else 0
        for i in range(first, n):
            t = total + lengths[i]
            if t > max_length:
                continue
            h = hit | orbit_masks[i]
            if (full ^ h).bit_count() > max_length - t:
                continue
            w = winners
            for a in winner_letters[i]:
                w |= 1 << relabelling[a]
            seq.append(i)
            yield from extend(p if m and i == first else m + 1, t, h, w,
                              tuple(relabelling[a] for a in label_maps[i]))
            seq.pop()

    yield from extend(0, 0, 0, 0, tuple(range(d)))


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
        if complete:
            orbits = G.relabelling_orbits(v)
            # the group is trivial exactly when its orbits are singletons
            nontrivial = len(orbits) < G._d
            winner_masks, orbit_masks, label_maps = _loop_data(G, v, loops, orbits, nontrivial)
            necklaces = _loop_necklaces(lengths, max_length, min_length, primitive,
                                        orbit_masks, len(orbits), winner_masks, label_maps)
        else:
            necklaces = _loop_necklaces(lengths, max_length, min_length, primitive)
        for seq in necklaces:
            word = [move for i in seq for move in loops[i]]
            f = IETFlipSequence(G.canonical_labelling(v), word)
            f.close()
            yield f
