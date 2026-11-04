r"""
Closed walks in Rauzy diagrams and associated flip sequences

The reduced Rauzy diagram of a permutation (Abelian differential)
or generalized permutation (quadratic differential) has one vertex per
permutation of its Rauzy class up to relabelling. Two edges reach out
of each vertex, the right Rauzy moves ``t`` (top winner) and ``b`` (bottom
winner). A *closed walk* is a path in this diagram whose end is its start.
It is a *primitive closed walk* when it is not a proper power of a shorter
closed walk (read on the sequence of pairs (vertex, move), not on the word of
moves alone). It is *complete* when every letter is the winner of some move, up
to the orbits of the relabelling that closes it.

A reduced vertex is represented by its *canonical labelling*: the labelled
permutation, relabelled so that its top row is the alphabet of the input
permutation in its order (`0 1 \dots d-1` for a component of a stratum). Its
*key* is the tuple of the positions in the top row of the letters of its
bottom row, and vertices are compared by key.

This module holds the class :class:`LazyReducedRauzyDiagram`, whose method
:meth:`~LazyReducedRauzyDiagram.closed_walks` does the work. The entry point
for the user is ``iet.rauzy_closed_walks`` (see example below).

EXAMPLES::

    sage: from surface_dynamics import Stratum, iet
    sage: C = Stratum([2]).unique_component()
    sage: for f in iet.rauzy_closed_walks(C, 9):
    ....:     print(f)
    RauzyFlipSequence('0 1 2 3\n1 3 0 2', 'tbtbtttbt')
    RauzyFlipSequence('0 1 2 3\n1 3 0 2', 'tbbtttbtb')
    RauzyFlipSequence('0 1 2 3\n1 3 2 0', 'bttbttbtb')
    RauzyFlipSequence('0 1 2 3\n1 3 2 0', 'btbttbtb')
    RauzyFlipSequence('0 1 2 3\n1 3 2 0', 'btbttbbtb')
    RauzyFlipSequence('0 1 2 3\n1 3 2 0', 'bbtbtbttb')
    RauzyFlipSequence('0 1 2 3\n1 3 2 0', 'bbtbbtbtb')
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

from surface_dynamics.misc.necklaces import _weighted_necklaces

from .template import OrientablePermutationIET, FlippedPermutationIET
from .flip_sequence import RauzyFlipSequence


class LazyReducedRauzyDiagram:
    r"""
    The reduced Rauzy diagram of an Abelian permutation, built on demand.

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

    The canonical labellings use the alphabet of ``p``, in its order.

    INPUT:

    - ``p`` -- an irreducible permutation of Abelian differential, labelled
      or reduced

    EXAMPLES::

        sage: from surface_dynamics import iet
        sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
        sage: p = iet.Permutation('a b c d', 'd c b a')
        sage: G = LazyReducedRauzyDiagram(p)
        sage: G
        Lazy reduced Rauzy diagram on 4 letters (1 vertex held)
        sage: G.start()
        (3, 2, 1, 0)
        sage: len(G)
        1

    The alphabet of ``p`` is kept by the canonical labellings::

        sage: G.canonical_labelling(G.start())
        a b c d
        d c b a

    A reduced permutation gives the same diagram::

        sage: q = iet.Permutation('a b c d', 'd c b a', reduced=True)
        sage: LazyReducedRauzyDiagram(q).start()
        (3, 2, 1, 0)

    Calling some methods such as ``vertices()`` will compute the whole Rauzy diagram::

        sage: len(G.vertices())
        7
        sage: len(G)
        7
        sage: G
        Lazy reduced Rauzy diagram on 4 letters (7 vertices held)
    """
    def __init__(self, p):
        r"""
        TESTS::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: p = iet.Permutation('a b c', 'c b a')
            sage: G = LazyReducedRauzyDiagram(p)
            sage: H = LazyReducedRauzyDiagram(iet.Permutation('a b c', 'c b a', reduced=True))
            sage: G.start() == H.start()
            True
            sage: LazyReducedRauzyDiagram(iet.Permutation('a b', 'b a', flips='a'))
            Traceback (most recent call last):
            ...
            ValueError: p must be a permutation of interval exchange transformation
            sage: LazyReducedRauzyDiagram(iet.Permutation('a b c', 'a c b'))
            Traceback (most recent call last):
            ...
            ValueError: p must be irreducible
            sage: LazyReducedRauzyDiagram(iet.Permutation('a', 'a'))
            Traceback (most recent call last):
            ...
            ValueError: p must have at least two letters
        """
        if not isinstance(p, OrientablePermutationIET) or isinstance(p, FlippedPermutationIET):
            raise ValueError('p must be a permutation of interval exchange transformation')
        if not p.is_irreducible():
            raise ValueError('p must be irreducible')
        if len(p) < 2:
            raise ValueError('p must have at least two letters')
        self._d = len(p)
        self._alphabet = p.alphabet()
        self._start = self._key(p)
        self.clear()

    @staticmethod
    def _key(q):
        r"""
        Return the key of the permutation ``q``: the top position of each
        letter of its bottom row. It is the bottom row of the canonical
        labelling of ``q``.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: LazyReducedRauzyDiagram._key(iet.Permutation('a c b d', 'd a b c'))
            (3, 0, 2, 1)
            sage: LazyReducedRauzyDiagram._key(iet.Permutation('a c b d', 'd a b c', reduced=True))
            (3, 0, 2, 1)
        """
        return tuple(q._twin[1])

    def _canonical(self, key):
        r"""
        Return the canonical labelling of the vertex ``key``, as a new
        permutation, over the alphabet of the input.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
            sage: G._canonical((2, 0, 1))
            a b c
            c a b
        """
        from .constructors import Permutation
        L = list(self._alphabet)
        return Permutation(L, [L[i] for i in key], alphabet=self._alphabet, reduced=False)

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
            a b c d
            b d c a

        A stratum component keeps the letters ``0, 1, ..., d-1``::

            sage: from surface_dynamics import Stratum
            sage: C = Stratum([2]).unique_component()
            sage: G = LazyReducedRauzyDiagram(C.permutation_representative(reduced=False))
            sage: G.canonical_labelling(G.successor(G.start(), 1))
            0 1 2 3
            1 3 2 0

        An unknown vertex::

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
            s = self._key(self._labellings[key].rauzy_move(move, -1, check=False))
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
                    mu = self._labellings[u].rauzy_move(move, -1, check=False)._labels[0]
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
        :class:`RauzyFlipSequence`. The groups at two vertices are conjugate.

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

    def closed_walks(self, max_length, min_length=1, complete=True, primitive=True):
        r"""
        Iterate over the closed walks of the diagram with length between
        ``min_length`` and ``max_length``.

        INPUT:

        - ``max_length``, ``min_length`` -- (``min_length`` defaults to ``1``)
          non-negative integers, the bounds on the length of the walks

        - ``complete`` -- (default ``True``) whether to keep only the complete
          walks, those in which every letter is a winner up to the orbits of
          the relabelling that closes the walk

        - ``primitive`` -- (default ``True``) whether to keep only the
          primitive closed walks; if ``False``, the proper powers of closed
          walks are yielded as well

        The arguments are checked when the method is called, not at the first
        step of the iteration.

        OUTPUT:

        An iterator over closed :class:`RauzyFlipSequence`, each walk once up
        to rotation. A walk starts at the canonical labelling of its least
        vertex (over the alphabet of the input permutation), and its word of
        moves is the least among the rotations starting there (``t < b``); it
        is closed, as by :meth:`RauzyFlipSequence.close`: it carries the
        relabelling that maps its end back to its start. The walks come by
        increasing least vertex, then by increasing sequence of first-return
        loops at this vertex.

        ALGORITHM:

        For each vertex `v`, the first-return loops at `v` (the words from `v`
        back to `v` whose intermediate vertices are all greater than `v`) form
        a prefix code. The walks whose least vertex is `v`, in their canonical
        rotation, are the sequences of these loops that are necklaces over the
        loops ordered as words, and the primitive closed walks are the Lyndon
        words, the weight of a loop being its length. They are generated by
        :func:`~surface_dynamics.misc.necklaces.weighted_necklaces`, the
        prenecklace backtracking of [FreMai78]_, [FreKes86]_ and
        [CatRusSawSerMie00]_, with the total length bounded by
        ``max_length``, each node iterating only over the loops that fit in
        the length left. The loops themselves are found by a depth-first
        search pruned by the distance back to `v`, which also accumulates the
        winners of each loop and the relabelling at which it ends, so that no
        loop is replayed.

        The relabelling that closes a walk is the inverse of the product of
        the relabellings of its loops, computed along the backtracking, so no
        flip sequence is replayed either. When the relabelling group at `v`
        (see :meth:`relabelling_generators`) is trivial, every walk at `v`
        closes with the identity.

        With ``complete``, completeness is tracked along the backtracking. The
        relabelling that closes a walk at `v` lies in the relabelling group at
        `v`, computed once from Schreier generators along a spanning tree. The
        orbits of this group hit by the winners of a loop do not depend on the
        labelling the loop is run from: they are the bits of a mask per loop,
        and a subtree is cut when more orbits contain no winner than moves
        are left. When the group is trivial, this is exact pruning by the
        letters still missing; otherwise the winners and the relabelling are
        combined loop by loop, as a product in a semidirect product, and
        completeness is decided exactly on each walk: every cycle of its
        relabelling must contain a winner. A flip sequence is built for the
        walks yielded only.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
            sage: for f in G.closed_walks(3, complete=False):
            ....:     print(f._simplified_flip_sequence())
            t
            tbb
            bb
            tt
            ttb
            b

        The walks use the alphabet of the permutation::

            sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c d', 'd c b a'))
            sage: next(G.closed_walks(8))
            RauzyFlipSequence('a b c d\nb d c a', 'btbttbtb')

        The walks are closed, with the relabelling that :meth:`RauzyFlipSequence.close`
        gives::

            sage: from surface_dynamics import Stratum
            sage: from surface_dynamics.interval_exchanges.flip_sequence import RauzyFlipSequence
            sage: G = LazyReducedRauzyDiagram(Stratum([2, 2]).odd_component().permutation_representative(reduced=False))
            sage: walks = list(G.closed_walks(9, complete=False))
            sage: len(walks), sum(1 for f in walks if f._relabelling != list(range(7)))
            (1165, 38)
            sage: def closed(f):
            ....:     g = RauzyFlipSequence(f.start(), f._simplified_flip_sequence())
            ....:     g.close()
            ....:     return g
            sage: all(f == closed(f) and f.is_closed() for f in walks)
            True
            sage: all(f._check() is None for f in walks)
            True

        TESTS::

            sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
            sage: G.closed_walks(-1)
            Traceback (most recent call last):
            ...
            ValueError: max_length and min_length must be non-negative
            sage: G.closed_walks(3, min_length=-1)
            Traceback (most recent call last):
            ...
            ValueError: max_length and min_length must be non-negative
        """
        max_length = operator.index(max_length)
        min_length = operator.index(min_length)
        if max_length < 0 or min_length < 0:
            raise ValueError('max_length and min_length must be non-negative')
        return self._closed_walks(max_length, max(min_length, 1), bool(complete), bool(primitive))

    def _diagram_data(self):
        r"""
        Return the whole reduced Rauzy diagram on integer vertices, with the
        data of its edges, as flat lists.

        The vertex ``u`` is the ``u``-th key of :meth:`vertices`, so that the
        integers compare as the keys. The canonical labelling of each vertex is
        moved once per edge here, instead of once per edge of each walk or loop
        found later. The lists hold Python integers only, which index faster
        than nested containers in the search.

        OUTPUT: a tuple ``(vertices, top_successors, bot_successors,
        bot_winners, top_predecessors, bot_predecessors, top_maps, bot_maps)``,
        where ``d`` is the number of letters and ``n`` the number of vertices:

        - ``vertices`` -- the sorted list of the keys of the vertices

        - ``top_successors[u]``, ``bot_successors[u]`` -- the vertex reached
          from ``u`` by the top move (``0``) and by the bottom move (``1``)

        - ``bot_winners[u]`` -- the winner of the bottom move in the canonical
          labelling of ``u``, as an index in ``0, ..., d-1`` (the last entry of
          the key). There is no list for the top move: the top winner of a
          canonical labelling is always ``d - 1``

        - ``top_predecessors[w]``, ``bot_predecessors[w]`` -- the vertex with
          a top (bottom) edge to ``w``. It is unique: a Rauzy move is invertible
          on labelled permutations and the relabelling group acts freely

        - ``top_maps[d*u + b]``, ``bot_maps[d*u + b]`` -- the entry ``b`` of the
          top row ``m`` of the labelled permutation reached by the top (bottom)
          move from the canonical labelling of ``u``. This permutation is the
          canonical labelling of the successor with each letter ``b`` replaced
          by ``m[b]``; so from the canonical labelling of ``u`` relabelled by
          ``lab`` (each letter ``c`` replaced by ``lab[c]``), the move leads to
          the canonical labelling of the successor relabelled by the tuple
          ``operator.itemgetter(*m)(lab)``

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
            sage: vertices, top_succ, bot_succ, bot_winners, top_pred, bot_pred, top_maps, bot_maps = G._diagram_data()
            sage: vertices
            [(1, 2, 0), (2, 0, 1), (2, 1, 0)]
            sage: top_succ, bot_succ
            ([0, 2, 1], [2, 1, 0])
            sage: bot_winners
            [0, 1, 0]
            sage: top_pred, bot_pred
            ([0, 2, 1], [2, 1, 0])
            sage: bot_maps[0:3]
            [0, 2, 1]

        The winners and the top rows agree with the moves of the canonical
        labellings::

            sage: d = 3
            sage: all(bot_winners[u] == G.canonical_labelling(k)._labels[1][-1]
            ....:     and top_maps[d*u:d*u + d] == list(G.canonical_labelling(k).rauzy_move(0)._labels[0])
            ....:     and bot_maps[d*u:d*u + d] == list(G.canonical_labelling(k).rauzy_move(1)._labels[0])
            ....:     for u, k in enumerate(vertices))
            True

        Each vertex is the successor of its predecessors::

            sage: from surface_dynamics import Stratum
            sage: G = LazyReducedRauzyDiagram(Stratum([2, 2]).odd_component().permutation_representative(reduced=False))
            sage: vertices, top_succ, bot_succ, _, top_pred, bot_pred, _, _ = G._diagram_data()
            sage: all(top_succ[top_pred[w]] == w and bot_succ[bot_pred[w]] == w for w in range(len(vertices)))
            True
        """
        vertices = self.vertices()
        index = {k: j for j, k in enumerate(vertices)}
        n = len(vertices)
        top_successors = []
        bot_successors = []
        bot_winners = []
        top_predecessors = [None] * n
        bot_predecessors = [None] * n
        top_maps = []
        bot_maps = []
        for u, k in enumerate(vertices):
            t = index[self.successor(k, 0)]
            b = index[self.successor(k, 1)]
            top_successors.append(t)
            bot_successors.append(b)
            top_predecessors[t] = u
            bot_predecessors[b] = u
            q = self._labellings[k]
            # the top row of the canonical labelling is 0, 1, ..., d-1 and its
            # bottom row is the key: the bottom winner is its last letter
            bot_winners.append(k[-1])
            top_maps.extend(q.rauzy_move(0, -1, check=False)._labels[0])
            bot_maps.extend(q.rauzy_move(1, -1, check=False)._labels[0])
        return (vertices, top_successors, bot_successors, bot_winners,
                top_predecessors, bot_predecessors, top_maps, bot_maps)

    def _closed_walks(self, max_length, min_length, complete, primitive):
        r"""
        Iterate over the closed walks, the arguments being checked as in
        :meth:`closed_walks`.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
            sage: len(list(G._closed_walks(3, 1, False, True)))
            6
        """
        (vertices, top_successors, bot_successors, bot_winners,
         top_predecessors, bot_predecessors, top_maps, bot_maps) = self._diagram_data()
        d = self._d
        n = len(vertices)
        identity = list(range(d))
        # built once for the whole search: the loops of each vertex share them
        top_getters = [operator.itemgetter(*top_maps[d*u:d*u + d]) for u in range(n)]
        bot_getters = [operator.itemgetter(*bot_maps[d*u:d*u + d]) for u in range(n)]

        for v, key in enumerate(vertices):
            dist = _distances_to(v, top_predecessors, bot_predecessors, max_length)
            orbits = self.relabelling_orbits(key)
            # the group is trivial exactly when its orbits are singletons; then
            # every walk at v closes with the identity, and no label map is needed
            nontrivial = len(orbits) < d
            loops, winner_masks, label_maps = _first_return_loops(
                d, v, top_successors, bot_successors, bot_winners, top_getters, bot_getters,
                dist, max_length, nontrivial)
            lengths = [len(loop) for loop in loops]
            masks = None
            full = None
            if complete:
                full = (1 << len(orbits)) - 1
                orbit_of = [None] * d
                for i, orbit in enumerate(orbits):
                    for a in orbit:
                        orbit_of[a] = i
                masks = []
                for winners in winner_masks:
                    hit = 0
                    a = 0
                    while winners:
                        if winners & 1:
                            hit |= 1 << orbit_of[a]
                        winners >>= 1
                        a += 1
                    masks.append(hit)
            # each walk owns its start and end: _from_data stores them as given
            canonical_labelling = self.canonical_labelling

            if not nontrivial:
                # the orbits are the letters: with complete, the masks of the
                # letters hit decide completeness, and every walk closes with
                # the identity
                # seq is the working list of the search: read it before the
                # next step
                for seq in _weighted_necklaces(lengths, max_length, min_length, primitive,
                                               masks, full, None, None):
                    moves = [(move, -1) for i in seq for move in loops[i]]
                    yield RauzyFlipSequence._from_data(canonical_labelling(key), canonical_labelling(key),
                                                     moves, identity[:])
                continue

            # the value of a walk is (W, r): the walk run from the canonical
            # labelling of v has the winners of the bitmask W and ends at the
            # canonical labelling of v with each letter a replaced by r[a]
            for seq, (winners, product) in _weighted_necklaces(
                    lengths, max_length, min_length, primitive, masks, full,
                    list(zip(winner_masks, label_maps)), _semidirect_product):
                if complete and not _hits_every_cycle(winners, product):
                    continue
                # the inverse of r closes the walk
                relabelling = [None] * d
                for a, b in enumerate(product):
                    relabelling[b] = a
                moves = [(move, -1) for i in seq for move in loops[i]]
                yield RauzyFlipSequence._from_data(canonical_labelling(key), canonical_labelling(key),
                                                 moves, relabelling)


def _distances_to(v, top_predecessors, bot_predecessors, max_length):
    r"""
    Return the distances to the vertex ``v``, in the graph below, along the
    paths whose vertices but ``v`` are all greater than ``v``.

    The graph has the vertices ``0, ..., n-1``, where ``n`` is the length of
    ``top_predecessors`` and of ``bot_predecessors``, and the ``2n`` edges
    ``u -> top_successors[u]``, labelled by the move ``0``, and
    ``u -> bot_successors[u]``, labelled by the move ``1``. The two
    arguments are the same graph read backwards: ``top_predecessors[w]`` is
    the vertex ``u`` with ``top_successors[u] == w`` and
    ``bot_predecessors[w]`` the vertex ``u`` with ``bot_successors[u] == w``
    (each vertex has one edge of each label coming in).

    The distance from ``u`` to ``v`` is the least number of edges of a path
    from ``u`` to ``v`` all of whose vertices but the last (``v``) are
    greater than ``v``, found by a breadth-first search from ``v`` along the
    edges backwards. ``max_length`` bounds the distances sought.

    OUTPUT: the list ``dist`` of length ``n``: ``dist[u]`` is the distance
    from ``u`` to ``v`` if ``u > v`` and it is at most ``max_length``, and
    ``None`` otherwise (in particular for ``u <= v``).

    EXAMPLES:

    The graph of the diagram of ``a b c / c b a``, with the edges
    ``0 -> 0, 1 -> 2, 2 -> 1`` of label ``0`` and ``0 -> 2, 1 -> 1, 2 -> 0``
    of label ``1``::

        sage: from surface_dynamics.interval_exchanges.rauzy_walks import _distances_to
        sage: top_pred = [0, 2, 1]
        sage: bot_pred = [2, 1, 0]
        sage: _distances_to(0, top_pred, bot_pred, 5)
        [None, 2, 1]
        sage: _distances_to(1, top_pred, bot_pred, 5)
        [None, None, 1]
        sage: _distances_to(0, top_pred, bot_pred, 1)
        [None, None, 1]

    A graph with a cycle ``0 -> 1 -> 2 -> 3 -> 0`` of label ``0``, the edges
    of label ``1`` being the loops ``u -> u``::

        sage: _distances_to(0, [3, 0, 1, 2], [0, 1, 2, 3], 5)
        [None, 3, 2, 1]
        sage: _distances_to(0, [3, 0, 1, 2], [0, 1, 2, 3], 2)
        [None, None, 2, 1]
        sage: _distances_to(1, [3, 0, 1, 2], [0, 1, 2, 3], 5)
        [None, None, None, None]
    """
    dist = [None] * len(top_predecessors)
    frontier = [v]
    level = 0
    while frontier and level < max_length:
        level += 1
        new = []
        for w in frontier:
            u = top_predecessors[w]
            if u > v and dist[u] is None:
                dist[u] = level
                new.append(u)
            u = bot_predecessors[w]
            if u > v and dist[u] is None:
                dist[u] = level
                new.append(u)
        frontier = new
    return dist


def _first_return_loops(d, v, top_successors, bot_successors, bot_winners, top_getters, bot_getters, dist, max_length, label_maps):
    r"""
    Return the first-return loops at the vertex ``v`` of length at most
    ``max_length``, with their winners and the relabellings at which they
    end.

    The arguments encode a graph with the vertices ``0, ..., n-1`` and the
    ``2n`` edges ``u -> top_successors[u]``, labelled by the move ``0``, and
    ``u -> bot_successors[u]``, labelled by the move ``1``. A path is read
    as the word of the labels of its edges. The first-return loops at ``v``
    are the words of the paths from ``v`` back to ``v`` whose intermediate
    vertices are all greater than ``v`` (``v`` is not met in between).
    ``max_length`` bounds the number of edges of a loop. ``dist`` is the
    output of :func:`_distances_to` for ``v`` in this graph and prunes the
    search: a path is extended only if it can still return to ``v`` within
    ``max_length`` edges.

    On an edge, a letter is a winner as follows. ``d`` is the number of
    letters. The canonical labelling of the vertex ``u`` has the letters
    ``0, ..., d-1``; the winner of the move ``0`` from it is the letter
    ``d - 1`` and the winner of the move ``1`` the letter ``bot_winners[u]``.
    The labelling the edge leads to is the canonical labelling of its end
    vertex with each letter ``b`` replaced by ``m[b]``, for a row ``m`` given
    by the getter of the edge: if the start is the canonical labelling with
    each letter ``c`` replaced by ``lab[c]``, ``lab`` a tuple, the end is the
    canonical labelling of the end vertex relabelled by
    ``top_getters[u](lab)`` (move ``0``) or ``bot_getters[u](lab)`` (move
    ``1``), the getters being the :func:`operator.itemgetter` of the rows
    ``top_maps[d*u:d*u + d]`` and ``bot_maps[d*u:d*u + d]`` of
    :meth:`LazyReducedRauzyDiagram._diagram_data`. The lists
    ``top_successors``, ``bot_successors`` and ``bot_winners`` are those of
    that method too.

    The loops are found by a depth-first search, which also accumulates
    along the current prefix its winners and the relabelling at which it
    ends, from the data of the edges: no loop is replayed.

    OUTPUT: a triple ``(loops, winner_masks, maps)``, the lists indexed by
    the loops sorted as words:

    - ``loops[j]`` -- the tuple of the moves of the loop ``j``

    - ``winner_masks[j]`` -- the bitmask of the letters of the canonical
      labelling of ``v`` that are winners of the loop ``j`` run from it

    - ``maps[j]`` -- if ``label_maps`` is set, the top row of the labelled
      permutation at which the loop ``j`` ends, a tuple ``m`` such that this
      permutation is the canonical labelling of ``v`` with each letter ``a``
      replaced by ``m[a]``; otherwise ``maps`` is ``None``

    EXAMPLES::

        sage: from surface_dynamics import iet
        sage: import operator
        sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram, _distances_to, _first_return_loops
        sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
        sage: vertices, top_succ, bot_succ, bot_winners, top_pred, bot_pred, top_maps, bot_maps = G._diagram_data()
        sage: top_getters = [operator.itemgetter(*top_maps[3*u:3*u + 3]) for u in range(3)]
        sage: bot_getters = [operator.itemgetter(*bot_maps[3*u:3*u + 3]) for u in range(3)]
        sage: dist = _distances_to(0, top_pred, bot_pred, 4)
        sage: loops, winner_masks, maps = _first_return_loops(3, 0, top_succ, bot_succ, bot_winners, top_getters, bot_getters, dist, 4, True)
        sage: loops
        [(0,), (1, 0, 0, 1), (1, 1)]
        sage: winner_masks
        [4, 3, 1]
        sage: _first_return_loops(3, 0, top_succ, bot_succ, bot_winners, top_getters, bot_getters, dist, 4, False)[2] is None
        True

    The winners and the relabellings agree with those of the flip sequence
    of each loop, including when the relabelling group is nontrivial::

        sage: from surface_dynamics import Stratum
        sage: from surface_dynamics.interval_exchanges.flip_sequence import RauzyFlipSequence
        sage: def replay(G, v, L):
        ....:     d = G._d
        ....:     vertices, top_succ, bot_succ, bot_winners, top_pred, bot_pred, top_maps, bot_maps = G._diagram_data()
        ....:     top_getters = [operator.itemgetter(*top_maps[d*u:d*u + d]) for u in range(len(vertices))]
        ....:     bot_getters = [operator.itemgetter(*bot_maps[d*u:d*u + d]) for u in range(len(vertices))]
        ....:     dist = _distances_to(v, top_pred, bot_pred, L)
        ....:     loops, winner_masks, maps = _first_return_loops(d, v, top_succ, bot_succ, bot_winners, top_getters, bot_getters, dist, L, True)
        ....:     flips = [RauzyFlipSequence(G.canonical_labelling(vertices[v]), loop) for loop in loops]
        ....:     return (len(loops), sum(m != tuple(range(G._d)) for m in maps),
        ....:             winner_masks == [sum(1 << list(G._alphabet).index(a) for a in set(f.winners())) for f in flips],
        ....:             maps == [tuple(f.end()._labels[0]) for f in flips])
        sage: replay(LazyReducedRauzyDiagram(iet.Permutation('a b c', 'c b a')), 0, 4)
        (3, 0, True, True)
        sage: replay(LazyReducedRauzyDiagram(iet.Permutation('a b c d', 'd c b a')), 0, 8)
        (8, 0, True, True)
        sage: G = LazyReducedRauzyDiagram(Stratum([2, 2]).odd_component().permutation_representative(reduced=False))
        sage: v = G.vertices().index((2, 0, 3, 5, 4, 6, 1))
        sage: G.relabelling_orbits(G.vertices()[v])
        [(0,), (1,), (2,), (3,), (4, 5, 6)]
        sage: replay(G, v, 9)
        (11, 4, True, True)
    """
    loops = []
    winner_masks = []
    maps = [] if label_maps else None
    word = []

    def search(u, length, lab, winners):
        # word leads from v to u in length moves, through vertices greater
        # than v; run from the canonical labelling of v, it has the winners
        # of the bitmask winners and ends at the canonical labelling of u
        # with each letter a replaced by lab[a]; the pruning ensures that
        # some path of length at most max_length - length leads from u back
        # to v
        w = top_successors[u]
        word.append(0)
        win = winners | (1 << lab[d - 1])
        if w == v:
            loops.append(tuple(word))
            winner_masks.append(win)
            if label_maps:
                maps.append(top_getters[u](lab))
        elif w > v:
            dw = dist[w]
            if dw is not None and length + 1 + dw <= max_length:
                search(w, length + 1, top_getters[u](lab), win)
        word[-1] = 1
        w = bot_successors[u]
        win = winners | (1 << lab[bot_winners[u]])
        if w == v:
            loops.append(tuple(word))
            winner_masks.append(win)
            if label_maps:
                maps.append(bot_getters[u](lab))
        elif w > v:
            dw = dist[w]
            if dw is not None and length + 1 + dw <= max_length:
                search(w, length + 1, bot_getters[u](lab), win)
        word.pop()

    try:
        if max_length >= 1:
            search(v, 0, tuple(range(d)), 0)
    finally:
        # search refers to itself through its closure: without this del, the
        # cycle keeps the loops alive until the cyclic garbage collector runs
        del search
    return loops, winner_masks, maps


def _semidirect_product(x, y):
    r"""
    Return the product of the pairs ``x = (W1, r1)`` and ``y = (W2, r2)``
    of a bitmask over the letters ``0, ..., d-1`` and a tuple of these
    letters.

    The product is ``(W1 | r1(W2), r1∘r2)``, where ``r1(W2)`` is the bitmask
    of the letters ``r1[a]`` for the bits ``a`` of ``W2`` and ``r1∘r2`` the
    tuple ``(r1[r2[0]], ..., r1[r2[d-1]])``. It is associative.

    For a walk from the canonical labelling of a vertex `v` to `v`, let
    ``W`` be the bitmask of the letters of this labelling that are winners
    of the walk, and let ``r`` be the tuple such that the walk ends at this
    labelling with each letter ``a`` replaced by ``r[a]`` (the pairs of
    :func:`_first_return_loops`). Then the pair of the concatenation of two
    walks at `v` is the product of their pairs: the second walk is run from
    the end of the first one, relabelled by ``r1``.

    EXAMPLES::

        sage: from surface_dynamics.interval_exchanges.rauzy_walks import _semidirect_product
        sage: _semidirect_product((0b001, (1, 0, 2)), (0b101, (0, 2, 1)))
        (7, (1, 2, 0))
        sage: _semidirect_product((0b001, (0, 2, 1)), (0b010, (1, 0, 2)))
        (5, (2, 0, 1))

    With two letters exchanged by the first loop, which has the winner
    ``1``, and a second loop ending with the identity, of winner ``0``: the
    walk ``(0, 0)`` has the winners ``1`` then ``0``, while ``(1, 1)``
    misses the letter ``1``::

        sage: a = (0b10, (1, 0))
        sage: b = (0b01, (0, 1))
        sage: _semidirect_product(a, a)
        (3, (0, 1))
        sage: _semidirect_product(b, b)
        (1, (0, 1))
    """
    W1, r1 = x
    W2, r2 = y
    while W2:
        low = W2 & -W2
        W1 |= 1 << r1[low.bit_length() - 1]
        W2 ^= low
    # d >= 2 for an irreducible permutation, so itemgetter returns a tuple
    return W1, operator.itemgetter(*r2)(r1)


def _hits_every_cycle(winners, relabelling):
    r"""
    Return whether every cycle of the permutation ``relabelling`` (a tuple,
    the image of ``a`` being ``relabelling[a]``) contains a bit of the
    bitmask ``winners``.

    A closed walk whose value under :func:`_semidirect_product` is
    ``(winners, relabelling)`` is complete exactly when this holds: the
    letters that are winners up to the relabelling closing the walk are the
    union of the cycles that contain a winner.

    EXAMPLES::

        sage: from surface_dynamics.interval_exchanges.rauzy_walks import _hits_every_cycle
        sage: _hits_every_cycle(0b001, (1, 0, 2))
        False
        sage: _hits_every_cycle(0b101, (1, 0, 2))
        True
        sage: _hits_every_cycle(0b010, (1, 2, 0))
        True
    """
    seen = 0
    for a in range(len(relabelling)):
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
