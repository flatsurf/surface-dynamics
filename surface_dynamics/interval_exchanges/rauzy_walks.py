r"""
Closed walks in Rauzy diagrams and associated flip sequences

The reduced Rauzy diagram of a permutation (Abelian differential)
or of a linear involution (generalized permutation, quadratic differential)
has one vertex per permutation of its Rauzy class up to relabelling. The edges
reaching out of a vertex are the right Rauzy moves ``t`` (top winner) and
``b`` (bottom winner). Both are defined for a permutation; for a linear
involution, one of them may be undefined (see ``has_rauzy_move``), so that a
vertex has one or two outgoing edges. A *closed walk* is a path in this
diagram whose end is its start. It is a *primitive closed walk* when it is not
a proper power of a shorter closed walk (read on the sequence of pairs
(vertex, move), not on the word of moves alone). It is *complete* when every
letter is the winner of some move, up to the orbits of the relabelling that
closes it.

A reduced vertex is represented by its *canonical labelling*: the labelled
permutation whose letters, read on the top row then on the bottom row, appear
for the first time in the order of the alphabet of the input (`0 1 \dots d-1`
for a component of a stratum). For a permutation, its top row is the alphabet
in its order. The *key* of a vertex encodes its canonical labelling with the
letters replaced by their indices in the alphabet: for a permutation, the
tuple of the bottom row (the positions in the top row of the letters of the
bottom row); for a linear involution, the pair of the tuples of the two rows.
Vertices are compared by key.

This module holds the class :class:`ReducedRauzyDiagram`, which builds the
whole reduced Rauzy diagram and its relabelling group at construction, and
whose method :meth:`~ReducedRauzyDiagram.closed_walks` does the work. The
entry point for the user is ``iet.rauzy_closed_walks`` (see example below).

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

The same for a component of a quadratic stratum, whose reduced Rauzy diagram
walked is the one of its permutation representative::

    sage: C = Stratum([2, 2], k=2).hyperelliptic_component()
    sage: for f in iet.rauzy_closed_walks(C, 8):
    ....:     print(f)
    RauzyFlipSequence('0 1 0 2 3\n3 1 4 2 4', 'tbttbbtb', relabelling=(0,1))
    RauzyFlipSequence('0 1 2 1\n3 4 3 2 0 4', 'tbtbbttb', relabelling=(3,4))
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

from .template import (OrientablePermutationIET, FlippedPermutationIET,
                       OrientablePermutationLI, FlippedPermutationLI, PermutationLI)
from .flip_sequence import RauzyFlipSequence


class ReducedRauzyDiagram:
    r"""
    The reduced Rauzy diagram of a permutation or of a generalized
    permutation.

    The whole diagram is built at construction, since a search of closed
    walks needs all of it: a walk is based at its least vertex. It is built
    once and never grows; it is the graph a search walks through, not a
    store of results. It holds:

    - the sorted list of the keys of the vertices, the vertex ``u`` being
      the ``u``-th key, and a dict from each key to its vertex: one entry
      per reduced vertex;

    - for each vertex and each right Rauzy move, top then bottom, the vertex
      reached, the winner of the move and the vertex from which the same
      move leads to it, as lists of integers indexed by the vertices, and
      the :func:`operator.itemgetter` of the relabelling of the move, or
      ``None`` when it is the identity: one entry per vertex and move, that
      is per edge or undefined move;

    - the relabelling group of the class: one labelling per reduced vertex,
      along a spanning tree from the start, the generators of the group at
      the start, at most one per edge, and its orbits on the letters, at
      most one per letter. The generators and orbits at another vertex are
      their images by the tree labelling of that vertex.

    The vertices are represented by their keys, as described in
    :mod:`~surface_dynamics.interval_exchanges.rauzy_walks`, and the moves
    are run on the keys; a canonical labelling is built only when asked for
    (:meth:`canonical_labelling`), as a new object. For a linear
    involution, a Rauzy move may be undefined: :meth:`successor` returns
    ``None`` for it, and it is no edge of the diagram.

    The canonical labellings use the alphabet of ``p``, in its order.

    INPUT:

    - ``p`` -- an irreducible permutation of Abelian differential, or an
      irreducible linear involution without flips, labelled or reduced

    EXAMPLES::

        sage: from surface_dynamics import iet
        sage: from surface_dynamics.interval_exchanges.rauzy_walks import ReducedRauzyDiagram
        sage: p = iet.Permutation('a b c d', 'd c b a')
        sage: G = ReducedRauzyDiagram(p)
        sage: G
        Reduced Rauzy diagram on 4 letters with 7 vertices
        sage: G.start()
        (3, 2, 1, 0)
        sage: len(G)
        7

    The alphabet of ``p`` is kept by the canonical labellings::

        sage: G.canonical_labelling(G.start())
        a b c d
        d c b a

    A reduced permutation gives the same diagram::

        sage: q = iet.Permutation('a b c d', 'd c b a', reduced=True)
        sage: H = ReducedRauzyDiagram(q)
        sage: H.start()
        (3, 2, 1, 0)
        sage: H.vertices() == G.vertices()
        True

    A linear involution, whose key is the pair of the rows of its canonical
    labelling, the letters being numbered by first appearance::

        sage: q = iet.GeneralizedPermutation('b b a', 'a c c', alphabet='abc')
        sage: G = ReducedRauzyDiagram(q)
        sage: G.start()
        ((0, 0, 1), (1, 2, 2))
        sage: G.canonical_labelling(G.start())
        a a b
        b c c
        sage: G.vertices()
        [((0, 0), (1, 1, 2, 2)),
         ((0, 0, 1), (1, 2, 2)),
         ((0, 0, 1, 1), (2, 2)),
         ((0, 1, 1), (2, 2, 0))]
    """
    def __init__(self, p):
        r"""
        Build the whole reduced Rauzy diagram of ``p`` and its relabelling
        group.

        The key of each vertex is moved once per edge (see
        :meth:`_key_move`), from the start, the vertices being numbered in
        the order of discovery, then renumbered in the order of their keys,
        so that the integers compare as the keys. The data of the edges are
        lists of Python integers indexed by the vertices, which index faster
        than nested containers in the search. With ``d`` the number of
        letters:

        - ``_vertices`` -- the sorted list of the keys of the vertices, and
          ``_index`` the dict from each key to its vertex, kept for the
          methods that take a key

        - ``_top_successors[u]``, ``_bot_successors[u]`` -- the vertex
          reached from ``u`` by the top move (``0``) and by the bottom move
          (``1``), ``-1`` if the move is undefined (for a linear involution
          only)

        - ``_top_winners[u]``, ``_bot_winners[u]`` -- the winner of the top
          (bottom) move in the canonical labelling of ``u``, its last letter
          on the top (bottom) row, as an index in ``0, ..., d-1``. The top
          winner of the canonical labelling of a permutation is always
          ``d - 1``, not that of a linear involution

        - ``_top_predecessors[w]``, ``_bot_predecessors[w]`` -- the vertex
          with a top (bottom) edge to ``w``, ``-1`` if there is none (for a
          linear involution only). It is unique: a Rauzy move is invertible
          on labelled permutations and commutes with the relabellings, so two
          vertices with an edge of the same label to ``w`` have canonical
          labellings that differ by a relabelling, hence are equal

        - ``_top_getters[u]``, ``_bot_getters[u]`` -- ``None`` if the
          relabelling ``m`` (see :meth:`_move_map`) of the labelled
          permutation reached by the top (bottom) move from the canonical
          labelling of ``u`` is the identity, else the
          :func:`operator.itemgetter` of ``m``: this permutation is the
          canonical labelling of the successor with each letter ``b``
          replaced by ``m[b]``. For a permutation, ``m`` is its top row. So
          from the canonical labelling of ``u`` relabelled by the tuple
          ``lab`` (each letter ``c`` replaced by ``lab[c]``), the move leads
          to the canonical labelling of the successor relabelled by the tuple
          ``_top_getters[u](lab)``, or by ``lab`` itself for ``None``, which
          the search then reuses rather than copies. The relabelling of an
          undefined move is the identity. The getters are shared: there is
          one per distinct relabelling, the edges with the same relabelling
          holding the same object.

          The canonical labelling numbers the letters in their order of
          first appearance, reading the top row then the bottom row, and a
          Rauzy move changes the rank of one letter only in that order: it
          moves one occurrence of the loser, so the other letters keep their
          relative order. So ``m`` is the identity or one rotation of an
          interval of ranks ``i < j``, either ``m[i] = j`` and
          ``m[k] = k - 1`` for ``i < k <= j``, or ``m[j] = i`` and
          ``m[k] = k + 1`` for ``i <= k < j``. For a permutation, the top
          moves are the identity, since the top row, which gives the order,
          is unchanged

        - ``_tree_labellings``, ``_root_generators``, ``_root_orbits`` --
          the relabelling group (see :meth:`_relabelling_group`)

        The entries ``-1`` compare below every vertex and are never used as
        an index by the search: an edge ``u -> w`` is followed only when
        ``w`` is ``v`` or greater than ``v``, at a base ``v >= 0``, and a
        predecessor only when it is greater than ``v``.

        TESTS::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import ReducedRauzyDiagram
            sage: p = iet.Permutation('a b c', 'c b a')
            sage: G = ReducedRauzyDiagram(p)
            sage: H = ReducedRauzyDiagram(iet.Permutation('a b c', 'c b a', reduced=True))
            sage: G.start() == H.start()
            True

        The data of the vertices and of the edges::

            sage: G._vertices
            [(1, 2, 0), (2, 0, 1), (2, 1, 0)]
            sage: all(G._index[k] == u for u, k in enumerate(G._vertices))
            True
            sage: G._top_successors, G._bot_successors
            ([0, 2, 1], [2, 1, 0])
            sage: G._top_winners, G._bot_winners
            ([2, 2, 2], [0, 1, 0])
            sage: G._top_predecessors, G._bot_predecessors
            ([0, 2, 1], [2, 1, 0])
            sage: G._bot_getters[0]((0, 1, 2))
            (0, 2, 1)

        The winners and the top rows agree with the moves of the canonical
        labellings::

            sage: def relabel(getter, lab):
            ....:     return lab if getter is None else getter(lab)
            sage: identity = (0, 1, 2)
            sage: all(G._bot_winners[u] == G.canonical_labelling(k)._labels[1][-1]
            ....:     and list(relabel(G._top_getters[u], identity)) == G.canonical_labelling(k).rauzy_move(0)._labels[0]
            ....:     and list(relabel(G._bot_getters[u], identity)) == G.canonical_labelling(k).rauzy_move(1)._labels[0]
            ....:     for u, k in enumerate(G._vertices))
            True
            sage: G._top_getters, G._bot_getters[1]
            ([None, None, None], None)

        Each vertex is the successor of its predecessors::

            sage: from surface_dynamics import Stratum
            sage: G = ReducedRauzyDiagram(Stratum([2, 2]).odd_component().permutation_representative(reduced=False))
            sage: all(G._top_successors[G._top_predecessors[w]] == w and G._bot_successors[G._bot_predecessors[w]] == w
            ....:     for w in range(len(G)))
            True

        For a linear involution, the undefined moves have the successor
        ``-1``; the winners and the relabellings agree with the moves of the
        canonical labellings, and each vertex is the successor of its
        predecessors::

            sage: G = ReducedRauzyDiagram(Stratum([2, 2], k=2).hyperelliptic_component().permutation_representative(reduced=False))
            sage: n, d = len(G), 5
            sage: n, G._top_successors.count(-1), G._bot_successors.count(-1)
            (73, 13, 13)
            sage: def agree(u, move, succ, winners, getters):
            ....:     q = G.canonical_labelling(G._vertices[u])
            ....:     m = relabel(getters[u], tuple(range(d)))
            ....:     if not q.has_rauzy_move(move):
            ....:         return succ[u] == -1 and m == tuple(range(d))
            ....:     r = q.rauzy_move(move)
            ....:     s = G.canonical_labelling(G._vertices[succ[u]])
            ....:     return (winners[u] == q._labels[move][-1]
            ....:             and [[m[b] for b in row] for row in s._labels] == r._labels)
            sage: all(agree(u, 0, G._top_successors, G._top_winners, G._top_getters)
            ....:     and agree(u, 1, G._bot_successors, G._bot_winners, G._bot_getters)
            ....:     for u in range(n))
            True
            sage: all(G._top_successors[G._top_predecessors[w]] == w for w in range(n) if G._top_predecessors[w] != -1)
            True
            sage: all(G._bot_successors[G._bot_predecessors[w]] == w for w in range(n) if G._bot_predecessors[w] != -1)
            True
            sage: G._top_predecessors.count(-1), G._bot_predecessors.count(-1)
            (13, 13)

        Invalid inputs::

            sage: ReducedRauzyDiagram(iet.Permutation('a b', 'b a', flips='a'))
            Traceback (most recent call last):
            ...
            ValueError: p must be a permutation of interval exchange transformation or a linear involution, without flips
            sage: ReducedRauzyDiagram(iet.GeneralizedPermutation('a a b', 'b c c', flips='a'))
            Traceback (most recent call last):
            ...
            ValueError: p must be a permutation of interval exchange transformation or a linear involution, without flips
            sage: ReducedRauzyDiagram(iet.Permutation('a b c', 'a c b'))
            Traceback (most recent call last):
            ...
            ValueError: p must be irreducible
            sage: ReducedRauzyDiagram(iet.GeneralizedPermutation('a b c a', 'b d d c'))
            Traceback (most recent call last):
            ...
            ValueError: p must be irreducible
            sage: ReducedRauzyDiagram(iet.Permutation('a', 'a'))
            Traceback (most recent call last):
            ...
            ValueError: p must have at least two letters
        """
        if isinstance(p, OrientablePermutationIET) and not isinstance(p, FlippedPermutationIET):
            self._quadratic = False
        elif isinstance(p, OrientablePermutationLI) and not isinstance(p, FlippedPermutationLI):
            self._quadratic = True
        else:
            raise ValueError('p must be a permutation of interval exchange transformation or a linear involution, without flips')
        if not p.is_irreducible():
            raise ValueError('p must be irreducible')
        if len(p) < 2:
            raise ValueError('p must have at least two letters')
        self._d = len(p)
        self._alphabet = p.alphabet()
        self._start = self._key(p)

        d = self._d
        key_move = self._key_move
        # the getter of each distinct relabelling cycle, shared by the edges
        # with that cycle (at most d - 2 cycles for a permutation, d (d - 1)
        # for a linear involution); the entries of a getter are letters, which
        # the renumbering of the vertices leaves unchanged
        getters = {}
        # discovery order: the vertex u is keys[u] until the renumbering
        index = {self._start: 0}
        keys = [self._start]
        top_successors = []
        bot_successors = []
        top_winners = []
        bot_winners = []
        top_getters = []
        bot_getters = []
        u = 0
        while u < len(keys):
            k = keys[u]
            for move, successors, winners, edge_getters in ((0, top_successors, top_winners, top_getters),
                                                            (1, bot_successors, bot_winners, bot_getters)):
                moved = key_move(k, move)
                if moved is None:
                    # the winner is read by the search before it checks the
                    # successor: it must be a letter
                    successors.append(-1)
                    winners.append(k[move][-1])
                    edge_getters.append(None)
                    continue
                s, w, cycle = moved
                j = index.get(s)
                if j is None:
                    j = len(keys)
                    index[s] = j
                    keys.append(s)
                successors.append(j)
                winners.append(w)
                if cycle is not None:
                    getter = getters.get(cycle)
                    if getter is None:
                        # the rotation m of the cycle (a, t): m[t] = a, and
                        # the ranks from a to t shift by one towards t
                        a, t = cycle
                        m = list(range(d))
                        m[t] = a
                        if t < a:
                            m[t + 1:a + 1] = range(t, a)
                        else:
                            m[a:t] = range(a + 1, t + 1)
                        getter = getters[cycle] = operator.itemgetter(*m)
                    edge_getters.append(getter)
                else:
                    edge_getters.append(None)
            u += 1
        del getters
        n = len(keys)

        # the sorted keys, the dict rewritten in place to the sorted indices,
        # and rank[u] the sorted index of the vertex u of the discovery order
        vertices = sorted(keys)
        for j, k in enumerate(vertices):
            index[k] = j
        rank = [index[k] for k in keys]
        del keys
        self._vertices = vertices
        self._index = index

        def renumbered(values, renumber):
            res = [-1] * n
            for u, x in enumerate(values):
                res[rank[u]] = rank[x] if renumber and x >= 0 else x
            return res

        self._top_winners = renumbered(top_winners, False)
        self._bot_winners = renumbered(bot_winners, False)
        del top_winners, bot_winners
        self._top_successors = top_successors = renumbered(top_successors, True)
        self._bot_successors = bot_successors = renumbered(bot_successors, True)
        top_predecessors = [-1] * n
        bot_predecessors = [-1] * n
        for u in range(n):
            w = top_successors[u]
            if w >= 0:
                top_predecessors[w] = u
            w = bot_successors[u]
            if w >= 0:
                bot_predecessors[w] = u
        self._top_predecessors = top_predecessors
        self._bot_predecessors = bot_predecessors

        # an identity relabelling has no getter, so that the search reuses the
        # relabelling it holds
        self._top_getters = renumbered(top_getters, False)
        self._bot_getters = renumbered(bot_getters, False)
        del top_getters, bot_getters, rank

        self._relabelling_group()

    @staticmethod
    def _key(q):
        r"""
        Return the key of the permutation or linear involution ``q``.

        For a permutation, it is the top position of each letter of its bottom
        row, which is the bottom row of the canonical labelling of ``q``. For a
        linear involution, it is the pair of the rows of ``q`` with each letter
        replaced by its rank of first appearance, reading the top row then the
        bottom row, which are the rows of the canonical labelling of ``q``.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import ReducedRauzyDiagram
            sage: ReducedRauzyDiagram._key(iet.Permutation('a c b d', 'd a b c'))
            (3, 0, 2, 1)
            sage: ReducedRauzyDiagram._key(iet.Permutation('a c b d', 'd a b c', reduced=True))
            (3, 0, 2, 1)
            sage: ReducedRauzyDiagram._key(iet.GeneralizedPermutation('c a c', 'a d b b d'))
            ((0, 1, 0), (1, 2, 3, 3, 2))
            sage: ReducedRauzyDiagram._key(iet.GeneralizedPermutation('c a c', 'a d b b d', reduced=True))
            ((0, 1, 0), (1, 2, 3, 3, 2))
        """
        if not isinstance(q, PermutationLI):
            return tuple(q._twin[1])
        # the twin data, unlike the labels, exists for a reduced linear
        # involution too
        twin = q._twin
        rows = ([None] * len(twin[0]), [None] * len(twin[1]))
        n = 0
        for i in (0, 1):
            row = rows[i]
            for j in range(len(row)):
                if row[j] is None:
                    row[j] = n
                    i2, j2 = twin[i][j]
                    rows[i2][j2] = n
                    n += 1
        return (tuple(rows[0]), tuple(rows[1]))

    @staticmethod
    def _move_map(r):
        r"""
        Return the relabelling of the labelled permutation or linear
        involution ``r``: the list ``m`` such that ``r`` is its canonical
        labelling with each letter ``b`` replaced by ``m[b]``.

        The letters of ``r`` are read as their indices in its alphabet. The
        list ``m`` holds them in the order of their first appearance, reading
        the top row then the bottom row; for a permutation, it is the top row.
        Applied to the result of a move from a canonical labelling, it gives
        the relabelling at which the move ends.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import ReducedRauzyDiagram
            sage: ReducedRauzyDiagram._move_map(iet.Permutation('a c b d', 'd a b c', alphabet='abcd'))
            [0, 2, 1, 3]
            sage: q = iet.GeneralizedPermutation('a b a', 'b c d d c').rauzy_move(0)
            sage: q
            c a b a
            b c d d
            sage: ReducedRauzyDiagram._move_map(q)
            [2, 0, 1, 3]
        """
        top, bot = r._labels
        if not isinstance(r, PermutationLI):
            return list(top)
        m = []
        seen = 0
        for row in (top, bot):
            for a in row:
                if not (seen >> a) & 1:
                    seen |= 1 << a
                    m.append(a)
        return m

    def _key_move(self, key, move):
        r"""
        Return the key reached from the vertex ``key`` by the right Rauzy move
        ``move``, with the winner and the relabelling of the move, or ``None``
        if the move is undefined.

        It moves the rows of the key directly, so that the search builds no
        permutation object per vertex; the Rauzy moves of ``template.py``
        stay the reference, as the doctests below show.

        INPUT:

        - ``key`` -- the key of a vertex

        - ``move`` -- ``0`` (top) or ``1`` (bottom)

        OUTPUT:

        ``None`` if the move is undefined (for a linear involution only, as
        tested by ``has_rauzy_move``), else a tuple
        ``(successor, winner, cycle)``:

        - ``successor`` -- the key of the vertex reached

        - ``winner`` -- the winner of the move in the canonical labelling of
          ``key``, the last letter of its top (bottom) row, as an index in
          ``0, ..., d-1``

        - ``cycle`` -- ``None`` if the relabelling ``m`` such that the
          labelled permutation reached by the move from the canonical
          labelling of ``key`` is the canonical labelling of ``successor``
          with each letter ``b`` replaced by ``m[b]`` (see
          :meth:`_move_map`) is the identity, else the pair ``(s, t)`` of
          ranks, ``s != t``, of the rotation that ``m`` is: ``m[t] = s``, and
          the ranks strictly between ``t`` and ``s``, and ``s`` itself, are
          shifted by one towards ``t`` (``m[k] = k - 1`` for ``t < k <= s``
          when ``t < s``, ``m[k] = k + 1`` for ``s <= k < t`` when ``s < t``)

        EXAMPLES:

        A permutation::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import ReducedRauzyDiagram
            sage: def expand(cycle, d):
            ....:     m = list(range(d))
            ....:     if cycle is not None:
            ....:         s, t = cycle
            ....:         m[t] = s
            ....:         if t < s:
            ....:             m[t+1:s+1] = range(t, s)
            ....:         else:
            ....:             m[s:t] = range(s+1, t+1)
            ....:     return m
            sage: G = ReducedRauzyDiagram(iet.Permutation('a b c d', 'd c b a', alphabet='abcd'))
            sage: k = G.successor(G.start(), 1)
            sage: G._key_move(k, 1)
            ((2, 1, 3, 0), 0, (3, 1))
            sage: expand((3, 1), 4)
            [0, 3, 1, 2]
            sage: r = G.canonical_labelling(k).rauzy_move(1)
            sage: G._key(r), G._move_map(r)
            ((2, 1, 3, 0), [0, 3, 1, 2])

        A linear involution::

            sage: from surface_dynamics import Stratum
            sage: C = Stratum([2, 2], k=2).hyperelliptic_component()
            sage: G = ReducedRauzyDiagram(C.permutation_representative(reduced=False))
            sage: V = G.vertices()
            sage: V[0]
            ((0, 1, 0), (1, 2, 3, 2, 4, 3, 4))
            sage: G._key_move(V[0], 0)
            (((0, 1, 2, 1), (2, 3, 4, 3, 0, 4)), 0, (4, 0))
            sage: expand((4, 0), 5)
            [4, 0, 1, 2, 3]
            sage: r = G.canonical_labelling(V[0]).rauzy_move(0)
            sage: G._key(r), G._move_map(r)
            (((0, 1, 2, 1), (2, 3, 4, 3, 0, 4)), [4, 0, 1, 2, 3])

        An undefined move, the bottom move of the same vertex::

            sage: G.canonical_labelling(V[0]).has_rauzy_move(1)
            False
            sage: G._key_move(V[0], 1) is None
            True
        """
        if not self._quadratic:
            d = self._d
            if move == 0:
                # the top winner d - 1 stays; the last letter of the bottom
                # row goes right after d - 1 there
                j = key.index(d - 1) + 1
                return key[:j] + key[-1:] + key[j:-1], d - 1, None
            # the last letter d - 1 of the top row goes right after the
            # winner w there; the letters after w shift right by one
            w = key[-1]
            s = tuple(a if a <= w else (w + 1 if a == d - 1 else a + 1) for a in key)
            return s, w, (None if w == d - 2 else (d - 1, w + 1))
        rows = (list(key[0]), list(key[1]))
        win = rows[move]
        lose = rows[1 - move]
        a = win[-1]
        b = lose[-1]
        if a == b:
            return None
        if a in lose:
            lose.pop()
            lose.insert(lose.index(a) + 1, b)
        else:
            # the move is undefined when neither end letter appears in the
            # other row and exactly one letter is repeated in the loser row
            if b not in win and len(lose) - len(set(lose)) == 1:
                return None
            lose.pop()
            win.insert(win.index(a), b)
        rank = [-1] * self._d
        m = []
        for row in rows:
            for c in row:
                if rank[c] < 0:
                    rank[c] = len(m)
                    m.append(c)
        # m is the identity or a rotation of an interval of ranks, from its
        # first to its last entry that moves
        moved = [k for k, c in enumerate(m) if c != k]
        if not moved:
            cycle = None
        else:
            i = moved[0]
            j = moved[-1]
            cycle = (j, i) if m[i] == j else (i, j)
        return ((tuple([rank[c] for c in rows[0]]), tuple([rank[c] for c in rows[1]])), a, cycle)

    def _canonical(self, key):
        r"""
        Return the canonical labelling of the vertex ``key``, as a new
        permutation or linear involution, over the alphabet of the input.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import ReducedRauzyDiagram
            sage: G = ReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
            sage: G._canonical((2, 0, 1))
            a b c
            c a b
            sage: G = ReducedRauzyDiagram(iet.GeneralizedPermutation('a a b', 'b c c'))
            sage: G._canonical(((0, 1, 0), (2, 2, 1)))
            a b a
            c c b
        """
        L = list(self._alphabet)
        if self._quadratic:
            from .constructors import GeneralizedPermutation
            top, bot = key
            return GeneralizedPermutation([L[i] for i in top], [L[i] for i in bot],
                                          alphabet=self._alphabet, reduced=False)
        from .constructors import Permutation
        return Permutation(L, [L[i] for i in key], alphabet=self._alphabet, reduced=False)

    def __len__(self):
        r"""
        Return the number of vertices of the diagram.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import ReducedRauzyDiagram
            sage: G = ReducedRauzyDiagram(iet.Permutation('a b c d', 'd c b a'))
            sage: len(G)
            7
            sage: len(ReducedRauzyDiagram(iet.Permutation('a b', 'b a')))
            1
        """
        return len(self._vertices)

    def __repr__(self):
        r"""
        TESTS::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import ReducedRauzyDiagram
            sage: ReducedRauzyDiagram(iet.Permutation('a b c d', 'd c b a'))
            Reduced Rauzy diagram on 4 letters with 7 vertices
            sage: ReducedRauzyDiagram(iet.Permutation('a b', 'b a'))
            Reduced Rauzy diagram on 2 letters with 1 vertex
        """
        n = len(self)
        return 'Reduced Rauzy diagram on %d letters with %d %s' % (self._d, n, 'vertex' if n == 1 else 'vertices')

    def start(self):
        r"""
        Return the key of the vertex of the permutation the diagram was built
        from.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import ReducedRauzyDiagram
            sage: G = ReducedRauzyDiagram(iet.Permutation('b a c', 'c a b'))
            sage: G.start()
            (2, 1, 0)
        """
        return self._start

    def canonical_labelling(self, key):
        r"""
        Return the canonical labelling of the vertex ``key``, as a new object.

        The vertex must be in the reduced Rauzy class of the diagram. The
        labelling is built at each call: the diagram holds the keys only.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import ReducedRauzyDiagram
            sage: G = ReducedRauzyDiagram(iet.Permutation('a b c d', 'd c b a'))
            sage: G.canonical_labelling(G.successor(G.start(), 1))
            a b c d
            b d c a

        A stratum component keeps the letters ``0, 1, ..., d-1``::

            sage: from surface_dynamics import Stratum
            sage: C = Stratum([2]).unique_component()
            sage: G = ReducedRauzyDiagram(C.permutation_representative(reduced=False))
            sage: G.canonical_labelling(G.successor(G.start(), 1))
            0 1 2 3
            1 3 2 0

        An unknown vertex::

            sage: G.canonical_labelling((1, 0, 2, 3))
            Traceback (most recent call last):
            ...
            KeyError: (1, 0, 2, 3)
        """
        if key not in self._index:
            raise KeyError(key)
        return self._canonical(key)

    def successor(self, key, move):
        r"""
        Return the key of the vertex reached from the vertex ``key`` by the
        right Rauzy move ``move`` (``0`` for top, ``1`` for bottom), or
        ``None`` if this move is undefined (for a linear involution only).

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import ReducedRauzyDiagram
            sage: G = ReducedRauzyDiagram(iet.Permutation('a b c d', 'd c b a'))
            sage: k = G.start()
            sage: G.successor(k, 0)
            (3, 0, 2, 1)
            sage: G.successor(k, 1)
            (1, 3, 2, 0)
            sage: G.successor(G.successor(k, 1), 1)
            (2, 1, 3, 0)

        A move of a linear involution may be undefined::

            sage: G = ReducedRauzyDiagram(iet.GeneralizedPermutation('a a b', 'b c c'))
            sage: k = ((0, 0), (1, 1, 2, 2))
            sage: k in G.vertices()
            True
            sage: q = G.canonical_labelling(k)
            sage: q
            a a
            b b c c
            sage: q.has_rauzy_move(1)
            False
            sage: G.successor(k, 1) is None
            True
            sage: G.successor(k, 0)
            ((0, 1, 1), (2, 2, 0))
        """
        successors = (self._top_successors, self._bot_successors)[move]
        s = successors[self._index[key]]
        return None if s < 0 else self._vertices[s]

    def vertices(self):
        r"""
        Return the sorted list of the keys of the vertices of the reduced
        Rauzy class.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import ReducedRauzyDiagram
            sage: G = ReducedRauzyDiagram(iet.Permutation('a b c d', 'd c b a'))
            sage: G.vertices()
            [(1, 3, 0, 2),
             (1, 3, 2, 0),
             (2, 0, 3, 1),
             (2, 1, 3, 0),
             (3, 0, 2, 1),
             (3, 1, 0, 2),
             (3, 2, 1, 0)]

        For a linear involution, the undefined moves are no edges::

            sage: from surface_dynamics import Stratum
            sage: C = Stratum([2, 2], k=2).hyperelliptic_component()
            sage: G = ReducedRauzyDiagram(C.permutation_representative(reduced=False))
            sage: V = G.vertices()
            sage: len(V)
            73
            sage: V[0]
            ((0, 1, 0), (1, 2, 3, 2, 4, 3, 4))
            sage: sum(G.successor(k, move) is None for k in V for move in (0, 1))
            26
        """
        return list(self._vertices)

    def _relabelling_group(self):
        r"""
        Compute, once for the class, the labelling of each vertex along a
        spanning tree from the start, the generators of the relabelling
        group at the start and its orbits on the letters.

        The tree is walked breadth first on the lists of the successors, the
        relabelling of each edge being given by its getter, or the identity
        for an edge without one, with no Rauzy move run (see
        :meth:`__init__`).

        The tree labelling of a vertex ``u`` is stored as the relabelling
        ``t_u`` of the labelled permutation reached from the canonical
        labelling of the start along the tree: that permutation is the
        canonical labelling of ``u`` with each letter ``a`` replaced by
        ``t_u[a]`` (for a permutation, ``t_u`` is its top row). For each edge
        ``u -> w``, the move applied to the tree labelling of ``u`` gives the
        tree labelling of ``w`` relabelled by some ``g``; these Schreier
        generators ``g`` generate the group at the start. The undefined moves
        of a linear involution are no edges.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import ReducedRauzyDiagram
            sage: G = ReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
            sage: G._vertices
            [(1, 2, 0), (2, 0, 1), (2, 1, 0)]
            sage: G._tree_labellings
            [(0, 2, 1), (0, 1, 2), (0, 1, 2)]
            sage: G._root_generators
            []
            sage: G._root_orbits
            [(0,), (1,), (2,)]

        For a linear involution, the undefined moves are skipped; the group
        of `Q_2(2^2)^{hyp}` acts freely, with order the number of labelled
        vertices over the number of reduced ones, `8760 / 73 = 120`::

            sage: from surface_dynamics import Stratum
            sage: C = Stratum([2, 2], k=2).hyperelliptic_component()
            sage: G = ReducedRauzyDiagram(C.permutation_representative(reduced=False))
            sage: len(G._tree_labellings)
            73
            sage: gens = [[i + 1 for i in g] for g in G._root_generators]
            sage: PermutationGroup(gens).order()
            120
        """
        d = self._d
        identity = tuple(range(d))
        generators = []
        seen = {identity}

        def schreier(y, x):
            # the generator of an edge whose move leads to the tree labelling
            # y of its end relabelled
            g = [None] * d
            for a, b in zip(y, x):
                g[a] = b
            g = tuple(g)
            if g not in seen:
                seen.add(g)
                generators.append(g)

        top_successors = self._top_successors
        bot_successors = self._bot_successors
        top_getters = self._top_getters
        bot_getters = self._bot_getters
        labels = [None] * len(self._vertices)
        s = self._index[self._start]
        labels[s] = identity
        frontier = [s]
        while frontier:
            new = []
            for u in frontier:
                tu = labels[u]
                for w, getter in ((top_successors[u], top_getters[u]),
                                  (bot_successors[u], bot_getters[u])):
                    if w < 0:
                        continue
                    # the move from the canonical labelling of u relabelled
                    # by tu ends at the canonical labelling of w relabelled
                    # by x
                    x = tu if getter is None else getter(tu)
                    if labels[w] is None:
                        labels[w] = x
                        new.append(w)
                    else:
                        schreier(labels[w], x)
            frontier = new
        orbit_of = [None] * d
        orbits = []
        for a in range(d):
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
        self._tree_labellings = labels
        self._root_generators = generators
        self._root_orbits = orbits

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
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import ReducedRauzyDiagram
            sage: def order(G, k):
            ....:     return PermutationGroup([[i + 1 for i in g] for g in G.relabelling_generators(k)]).order()
            sage: G = ReducedRauzyDiagram(Stratum([2]).unique_component().permutation_representative(reduced=False))
            sage: G.relabelling_generators(G.start())
            []
            sage: G = ReducedRauzyDiagram(Stratum([2, 2]).odd_component().permutation_representative(reduced=False))
            sage: k = G.start()
            sage: order(G, k), order(G, G.vertices()[0])
            (3, 3)
            sage: C = Stratum([3, 1]).unique_component()
            sage: G = ReducedRauzyDiagram(C.permutation_representative(left_degree=1, reduced=False))
            sage: order(G, G.start()), order(G, G.vertices()[0])
            (2, 2)
            sage: all(sorted(g) == list(range(7)) for g in G.relabelling_generators(G.start()))
            True
        """
        t = self._tree_labellings[self._index[key]]
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
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import ReducedRauzyDiagram
            sage: G = ReducedRauzyDiagram(Stratum([2]).unique_component().permutation_representative(reduced=False))
            sage: G.relabelling_orbits(G.start())
            [(0,), (1,), (2,), (3,)]
            sage: G = ReducedRauzyDiagram(Stratum([2, 2]).odd_component().permutation_representative(reduced=False))
            sage: G.relabelling_orbits(G.start())
            [(0,), (1, 5, 6), (2,), (3,), (4,)]

        The relabelling group of a linear involution may be transitive on the
        letters, as for `Q_2(2^2)^{hyp}`::

            sage: G = ReducedRauzyDiagram(Stratum([2, 2], k=2).hyperelliptic_component().permutation_representative(reduced=False))
            sage: G.relabelling_orbits(G.start())
            [(0, 1, 2, 3, 4)]

        TESTS:

        At every vertex, the orbits are those of the group generated by
        :meth:`relabelling_generators`, for `H_3(2^2)^{odd}` and
        `Q_2(5, -1)`::

            sage: def orbits(G, k):
            ....:     gens = [[i + 1 for i in g] for g in G.relabelling_generators(k)]
            ....:     O = PermutationGroup(gens, domain=range(1, G._d + 1)).orbits()
            ....:     return sorted(tuple(sorted(i - 1 for i in o)) for o in O)
            sage: G = ReducedRauzyDiagram(Stratum([2, 2]).odd_component().permutation_representative(reduced=False))
            sage: len(G.vertices())
            294
            sage: all(G.relabelling_orbits(k) == orbits(G, k) for k in G.vertices())
            True
            sage: G = ReducedRauzyDiagram(Stratum([5, -1], k=2).unique_component().permutation_representative(reduced=False))
            sage: len(G.vertices())
            440
            sage: all(G.relabelling_orbits(k) == orbits(G, k) for k in G.vertices())
            True
        """
        t = self._tree_labellings[self._index[key]]
        tinv = [None] * self._d
        for a, b in enumerate(t):
            tinv[b] = a
        # the generators at key are a -> tinv[g[t[a]]], so the orbit of a at
        # key is the image by tinv of the orbit of t[a] at the start
        return sorted(tuple(sorted(tinv[a] for a in orbit)) for orbit in self._root_orbits)

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
        loops at this vertex. The walks at a vertex share one start object;
        each owns its end.

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
        search pruned by the distance back to `v`, computed by a backward
        breadth-first search stopped at about half of ``max_length``, the
        vertices beyond counting as farther; the search also accumulates the
        winners of each loop and the relabelling at which it ends, so that no
        loop is replayed. The undefined moves of a linear involution are no
        edges, so no loop goes through them.

        The relabelling that closes a walk is the inverse of the product of
        the relabellings of its loops, computed along the backtracking, so no
        flip sequence is replayed either. When the relabelling group at `v`
        (see :meth:`relabelling_generators`) is trivial, every walk at `v`
        closes with the identity. Without ``complete``, the relabellings alone
        are multiplied, by composition.

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
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import ReducedRauzyDiagram
            sage: G = ReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
            sage: for f in G.closed_walks(3, complete=False):
            ....:     print(f._simplified_flip_sequence())
            t
            tbb
            bb
            tt
            ttb
            b

        The walks use the alphabet of the permutation::

            sage: G = ReducedRauzyDiagram(iet.Permutation('a b c d', 'd c b a'))
            sage: next(G.closed_walks(8))
            RauzyFlipSequence('a b c d\nb d c a', 'btbttbtb')

        The walks are closed, with the relabelling that :meth:`RauzyFlipSequence.close`
        gives::

            sage: from surface_dynamics import Stratum
            sage: from surface_dynamics.interval_exchanges.flip_sequence import RauzyFlipSequence
            sage: G = ReducedRauzyDiagram(Stratum([2, 2]).odd_component().permutation_representative(reduced=False))
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

        The same on a linear involution, whose relabelling group is never
        trivial; the walks of `Q_2(2^2)^{hyp}`::

            sage: G = ReducedRauzyDiagram(Stratum([2, 2], k=2).hyperelliptic_component().permutation_representative(reduced=False))
            sage: walks = list(G.closed_walks(10, complete=False))
            sage: len(walks)
            145
            sage: all(f == closed(f) and f.is_closed() and f._check() is None for f in walks)
            True
            sage: [f for f in G.closed_walks(8)]
            [RauzyFlipSequence('0 1 0 2 3\n3 1 4 2 4', 'tbttbbtb', relabelling=(0,1)),
             RauzyFlipSequence('0 1 2 1\n3 4 3 2 0 4', 'tbtbbttb', relabelling=(3,4))]

        TESTS::

            sage: G = ReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
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

    def _closed_walks(self, max_length, min_length, complete, primitive):
        r"""
        Iterate over the closed walks, the arguments being checked as in
        :meth:`closed_walks`.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import ReducedRauzyDiagram
            sage: G = ReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
            sage: len(list(G._closed_walks(3, 1, False, True)))
            6
        """
        vertices = self._vertices
        top_successors = self._top_successors
        bot_successors = self._bot_successors
        top_winners = self._top_winners
        bot_winners = self._bot_winners
        top_predecessors = self._top_predecessors
        bot_predecessors = self._bot_predecessors
        top_getters = self._top_getters
        bot_getters = self._bot_getters
        d = self._d
        n = len(vertices)
        identity = list(range(d))

        # the lower bounds of the distances back to the base, shared by the
        # bases: each base sets the entries it reaches and resets them after.
        # The bounds (see _distances_to and _first_return_loops) are exact
        # up to the depth h and the vertices beyond count as farther, so the
        # pruning is exact on the last max_length // 2 + 2 moves of each
        # loop, while the backward search visits about as few vertices as
        # the forward one. Deeper, it visits more vertices than the pruning
        # saves; a depth above max_length bounds nothing more.
        h = min(max_length, max_length // 2 + 2)
        dist = [h + 1] * n

        # the orbits at v are the images by the inverse of the tree
        # labelling t of v of those at the start: the letter a lies in the
        # orbit of the start that holds t[a], whose index serves as the bit
        # of the orbit (completeness reads the number of bits only)
        tree = self._tree_labellings
        num_orbits = len(self._root_orbits)
        root_index = [None] * d
        for i, orbit in enumerate(self._root_orbits):
            for a in orbit:
                root_index[a] = i
        # the groups at the vertices are conjugate: the group is trivial at
        # every vertex exactly when its orbits are singletons; then every
        # walk closes with the identity, and no label map is needed
        nontrivial = num_orbits < d
        full = (1 << num_orbits) - 1 if complete else None

        for v, key in enumerate(vertices):
            reached = _distances_to(v, top_predecessors, bot_predecessors, h, dist)
            loops, winner_masks, label_maps = _first_return_loops(
                d, v, top_successors, bot_successors, top_winners, bot_winners,
                top_getters, bot_getters, dist, max_length, nontrivial)
            # the next base needs every entry back at h + 1; resetting the
            # reached entries only costs no O(n) per base
            for u in reached:
                dist[u] = h + 1
            lengths = [len(loop) for loop in loops]
            masks = None
            if complete:
                if not nontrivial:
                    # the orbits are the letters
                    masks = winner_masks
                else:
                    t = tree[v]
                    orbit_masks = [0] * num_orbits
                    for a in range(d):
                        orbit_masks[root_index[t[a]]] |= 1 << a
                    masks = []
                    for winners in winner_masks:
                        hit = 0
                        for i, orbit_mask in enumerate(orbit_masks):
                            if winners & orbit_mask:
                                hit |= 1 << i
                        masks.append(hit)
            # the walks at v share one start: no method of RauzyFlipSequence
            # mutates its start in place, while the end is moved, relabelled
            # and inverted in place, so each walk owns its end, a copy of the
            # start; the start is built at the first walk only, since a vertex
            # may have none
            start = None
            canonical_labelling = self._canonical

            if not nontrivial:
                # the orbits are the letters: with complete, the masks of the
                # letters hit decide completeness, and every walk closes with
                # the identity
                # seq is the working list of the search: read it before the
                # next step
                for seq in _weighted_necklaces(lengths, max_length, min_length, primitive,
                                               masks, full, None, None):
                    moves = [(move, -1) for i in seq for move in loops[i]]
                    if start is None:
                        start = canonical_labelling(key)
                    yield RauzyFlipSequence._from_data(start, start.__copy__(),
                                                     moves, identity[:])
                continue

            if not complete:
                # the winners are not needed: the value of a walk is the tuple
                # r such that the walk run from the canonical labelling of v
                # ends at the canonical labelling of v with each letter a
                # replaced by r[a]
                for seq, product in _weighted_necklaces(
                        lengths, max_length, min_length, primitive, None, None,
                        label_maps, _compose):
                    # the inverse of r closes the walk
                    relabelling = [None] * d
                    for a, b in enumerate(product):
                        relabelling[b] = a
                    moves = [(move, -1) for i in seq for move in loops[i]]
                    if start is None:
                        start = canonical_labelling(key)
                    yield RauzyFlipSequence._from_data(start, start.__copy__(),
                                                     moves, relabelling)
                continue

            # the value of a walk is (W, r): the walk run from the canonical
            # labelling of v has the winners of the bitmask W and ends at the
            # canonical labelling of v with each letter a replaced by r[a]
            for seq, (winners, product) in _weighted_necklaces(
                    lengths, max_length, min_length, primitive, masks, full,
                    list(zip(winner_masks, label_maps)), _semidirect_product):
                if not _hits_every_cycle(winners, product):
                    continue
                # the inverse of r closes the walk
                relabelling = [None] * d
                for a, b in enumerate(product):
                    relabelling[b] = a
                moves = [(move, -1) for i in seq for move in loops[i]]
                if start is None:
                    start = canonical_labelling(key)
                yield RauzyFlipSequence._from_data(start, start.__copy__(),
                                                 moves, relabelling)


def _distances_to(v, top_predecessors, bot_predecessors, depth, dist):
    r"""
    Set the distances to the vertex ``v`` up to ``depth``, in the graph
    below, along the paths whose vertices but ``v`` are all greater than
    ``v``, and return the vertices set.

    The graph has the vertices ``0, ..., n-1``, where ``n`` is the length of
    ``top_predecessors`` and of ``bot_predecessors``, and at most ``2n``
    edges: ``u -> top_successors[u]``, labelled by the move ``0``, and
    ``u -> bot_successors[u]``, labelled by the move ``1``, when these
    successors are defined (an undefined move of a linear involution is no
    edge). The two
    arguments are the same graph read backwards: ``top_predecessors[w]`` is
    the vertex ``u`` with ``top_successors[u] == w`` and
    ``bot_predecessors[w]`` the vertex ``u`` with ``bot_successors[u] == w``
    (each vertex has at most one edge of each label coming in, and ``-1``
    stands for none; see :class:`ReducedRauzyDiagram` for why this holds
    and why ``-1`` is never followed).

    The distance from ``u`` to ``v`` is the least number of edges of a path
    from ``u`` to ``v`` all of whose vertices but the last (``v``) are
    greater than ``v``, found by a breadth-first search from ``v`` along the
    edges backwards, stopped after ``depth`` levels.

    ``dist`` is a list of length ``n`` whose entries are all ``depth + 1``
    on input. Afterwards, ``dist[u]`` is the distance from ``u`` to ``v`` at
    the vertices ``u > v`` at distance at most ``depth``, and stays
    ``depth + 1`` elsewhere: for ``u > v`` this is a lower bound of the
    distance (``u`` is farther, or cannot reach ``v``). The list is filled
    in place rather than built per call: one list serves all the bases
    ``v`` of a search, and the caller resets the returned entries to
    ``depth + 1``, so that a base costs the vertices it reaches and no
    allocation of length ``n``.

    OUTPUT: the list of the vertices ``u`` whose entry was set, in the order
    of the search (by increasing distance).

    EXAMPLES:

    The graph of the diagram of ``a b c / c b a``, with the edges
    ``0 -> 0, 1 -> 2, 2 -> 1`` of label ``0`` and ``0 -> 2, 1 -> 1, 2 -> 0``
    of label ``1``::

        sage: from surface_dynamics.interval_exchanges.rauzy_walks import _distances_to
        sage: top_pred = [0, 2, 1]
        sage: bot_pred = [2, 1, 0]
        sage: dist = [6] * 3
        sage: _distances_to(0, top_pred, bot_pred, 5, dist), dist
        ([2, 1], [6, 2, 1])
        sage: dist = [6] * 3
        sage: _distances_to(1, top_pred, bot_pred, 5, dist), dist
        ([2], [6, 6, 1])
        sage: dist = [2] * 3
        sage: _distances_to(0, top_pred, bot_pred, 1, dist), dist
        ([2], [2, 2, 1])
        sage: dist = [1] * 3
        sage: _distances_to(0, top_pred, bot_pred, 0, dist), dist
        ([], [1, 1, 1])

    A graph with a cycle ``0 -> 1 -> 2 -> 3 -> 0`` of label ``0``, the edges
    of label ``1`` being the loops ``u -> u``; at depth ``2``, the vertex
    ``1`` at distance ``3`` keeps the bound ``3``::

        sage: dist = [6] * 4
        sage: _distances_to(0, [3, 0, 1, 2], [0, 1, 2, 3], 5, dist), dist
        ([3, 2, 1], [6, 3, 2, 1])
        sage: dist = [3] * 4
        sage: _distances_to(0, [3, 0, 1, 2], [0, 1, 2, 3], 2, dist), dist
        ([3, 2], [3, 3, 2, 1])
        sage: dist = [6] * 4
        sage: _distances_to(1, [3, 0, 1, 2], [0, 1, 2, 3], 5, dist), dist
        ([], [6, 6, 6, 6])

    The same cycle with no edge of label ``1``::

        sage: dist = [6] * 4
        sage: _distances_to(0, [3, 0, 1, 2], [-1, -1, -1, -1], 5, dist), dist
        ([3, 2, 1], [6, 3, 2, 1])
    """
    reached = []
    frontier = [v]
    level = 0
    while frontier and level < depth:
        level += 1
        new = []
        for w in frontier:
            # an entry above depth is unset: the set ones are at most depth
            u = top_predecessors[w]
            if u > v and dist[u] > depth:
                dist[u] = level
                new.append(u)
            u = bot_predecessors[w]
            if u > v and dist[u] > depth:
                dist[u] = level
                new.append(u)
        reached.extend(new)
        frontier = new
    return reached


def _first_return_loops(d, v, top_successors, bot_successors, top_winners, bot_winners, top_getters, bot_getters, dist, max_length, label_maps):
    r"""
    Return the first-return loops at the vertex ``v`` of length at most
    ``max_length``, with their winners and the relabellings at which they
    end.

    The arguments encode a graph with the vertices ``0, ..., n-1`` and the
    edges ``u -> top_successors[u]``, labelled by the move ``0``, and
    ``u -> bot_successors[u]``, labelled by the move ``1``, a successor
    ``-1`` standing for no edge (it is below ``v``, so it is never
    followed). A path is read as the word of the labels of its edges. The
    first-return loops at ``v``
    are the words of the paths from ``v`` back to ``v`` whose intermediate
    vertices are all greater than ``v`` (``v`` is not met in between).
    ``max_length`` bounds the number of edges of a loop. ``dist`` is the
    list filled by :func:`_distances_to` for ``v`` in this graph, at any
    depth: at each vertex ``w > v``, a lower bound of the distance from
    ``w`` back to ``v`` (exact up to the depth). It prunes the search: a
    path of ``length`` edges is extended to ``w`` only if
    ``length + 1 + dist[w] <= max_length``. The pruning is sound since the
    bound never exceeds the distance, so no path that can still return to
    ``v`` within ``max_length`` edges is cut; a smaller depth cuts fewer
    paths far from ``v``, at no cost to the result.

    On an edge, a letter is a winner as follows. ``d`` is the number of
    letters. The canonical labelling of the vertex ``u`` has the letters
    ``0, ..., d-1``; the winner of the move ``0`` from it is the letter
    ``top_winners[u]`` and the winner of the move ``1`` the letter
    ``bot_winners[u]``. The labelling the edge leads to is the canonical
    labelling of its end vertex with each letter ``b`` replaced by ``m[b]``,
    for a relabelling ``m`` given by the getter of the edge: if the start is
    the canonical labelling with each letter ``c`` replaced by ``lab[c]``,
    ``lab`` a tuple, the end is the canonical labelling of the end vertex
    relabelled by ``top_getters[u](lab)`` (move ``0``) or
    ``bot_getters[u](lab)`` (move ``1``). The getter is ``None`` when ``m``
    is the identity, and the end is then relabelled by ``lab`` itself, with
    no new tuple. A Rauzy move changes the rank of one letter only in the
    order of first appearance that numbers the letters of a canonical
    labelling, so ``m`` is the identity or one rotation of an interval of
    ranks; for a permutation, every top move has the identity. These lists
    are the attributes of the same names, with a leading underscore, of
    :class:`ReducedRauzyDiagram` (see its ``__init__``).

    The loops are found by a depth-first search, which also accumulates
    along the current prefix its winners and the relabelling at which it
    ends, from the data of the edges: no loop is replayed.

    OUTPUT: a triple ``(loops, winner_masks, maps)``, the lists indexed by
    the loops sorted as words:

    - ``loops[j]`` -- the tuple of the moves of the loop ``j``

    - ``winner_masks[j]`` -- the bitmask of the letters of the canonical
      labelling of ``v`` that are winners of the loop ``j`` run from it.
      They are always computed, but only a search of complete walks reads
      them: :meth:`ReducedRauzyDiagram._closed_walks` with
      ``complete=False`` drops them, the winner bit of each edge being
      computed for nothing there

    - ``maps[j]`` -- if ``label_maps`` is set, the relabelling of the
      labelled permutation at which the loop ``j`` ends, a tuple ``m`` such
      that this
      permutation is the canonical labelling of ``v`` with each letter ``a``
      replaced by ``m[a]``; otherwise ``maps`` is ``None``

    EXAMPLES::

        sage: from surface_dynamics import iet
        sage: from surface_dynamics.interval_exchanges.rauzy_walks import ReducedRauzyDiagram, _distances_to, _first_return_loops
        sage: G = ReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
        sage: dist = [5] * 3
        sage: _ = _distances_to(0, G._top_predecessors, G._bot_predecessors, 4, dist)
        sage: edges = (G._top_successors, G._bot_successors, G._top_winners, G._bot_winners,
        ....:          G._top_getters, G._bot_getters)
        sage: loops, winner_masks, maps = _first_return_loops(3, 0, *edges, dist, 4, True)
        sage: loops
        [(0,), (1, 0, 0, 1), (1, 1)]
        sage: winner_masks
        [4, 3, 1]
        sage: _first_return_loops(3, 0, *edges, dist, 4, False)[2] is None
        True

    The winners and the relabellings agree with those of the flip sequence
    of each loop, including when the relabelling group is nontrivial::

        sage: from surface_dynamics import Stratum
        sage: from surface_dynamics.interval_exchanges.flip_sequence import RauzyFlipSequence
        sage: def replay(G, v, L):
        ....:     dist = [L + 1] * len(G)
        ....:     _ = _distances_to(v, G._top_predecessors, G._bot_predecessors, L, dist)
        ....:     loops, winner_masks, maps = _first_return_loops(
        ....:         G._d, v, G._top_successors, G._bot_successors, G._top_winners, G._bot_winners,
        ....:         G._top_getters, G._bot_getters, dist, L, True)
        ....:     flips = [RauzyFlipSequence(G.canonical_labelling(G._vertices[v]), loop) for loop in loops]
        ....:     return (len(loops), sum(m != tuple(range(G._d)) for m in maps),
        ....:             winner_masks == [sum(1 << list(G._alphabet).index(a) for a in set(f.winners())) for f in flips],
        ....:             maps == [tuple(G._move_map(f.end())) for f in flips])
        sage: replay(ReducedRauzyDiagram(iet.Permutation('a b c', 'c b a')), 0, 4)
        (3, 0, True, True)
        sage: replay(ReducedRauzyDiagram(iet.Permutation('a b c d', 'd c b a')), 0, 8)
        (8, 0, True, True)
        sage: G = ReducedRauzyDiagram(Stratum([2, 2]).odd_component().permutation_representative(reduced=False))
        sage: v = G.vertices().index((2, 0, 3, 5, 4, 6, 1))
        sage: G.relabelling_orbits(G.vertices()[v])
        [(0,), (1,), (2,), (3,), (4, 5, 6)]
        sage: replay(G, v, 9)
        (11, 4, True, True)

    On a linear involution, whose top winner is not always the last letter
    and whose moves may be undefined::

        sage: G = ReducedRauzyDiagram(Stratum([2, 2], k=2).hyperelliptic_component().permutation_representative(reduced=False))
        sage: replay(G, 0, 12)
        (17, 17, True, True)

    The loops, their winners and their relabellings do not depend on the
    depth of the backward search::

        sage: def by_depth(G, v, L, h):
        ....:     dist = [h + 1] * len(G)
        ....:     _ = _distances_to(v, G._top_predecessors, G._bot_predecessors, h, dist)
        ....:     return _first_return_loops(
        ....:         G._d, v, G._top_successors, G._bot_successors, G._top_winners, G._bot_winners,
        ....:         G._top_getters, G._bot_getters, dist, L, True)
        sage: G1 = ReducedRauzyDiagram(Stratum([2, 2]).odd_component().permutation_representative(reduced=False))
        sage: len(G1.vertices())
        294
        sage: all(by_depth(G1, v, 9, h) == by_depth(G1, v, 9, 9) for v in range(20) for h in range(10))
        True
        sage: all(by_depth(G, v, 12, h) == by_depth(G, v, 12, 12) for v in range(20) for h in range(13))
        True
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
        # the lower bound dist[u] of the distance from u back to v is at most
        # max_length - length (u != v)
        w = top_successors[u]
        word.append(0)
        win = winners | (1 << lab[top_winners[u]])
        if w == v:
            loops.append(tuple(word))
            winner_masks.append(win)
            if label_maps:
                getter = top_getters[u]
                maps.append(lab if getter is None else getter(lab))
        elif w > v:
            if length + 1 + dist[w] <= max_length:
                getter = top_getters[u]
                search(w, length + 1, lab if getter is None else getter(lab), win)
        word[-1] = 1
        w = bot_successors[u]
        win = winners | (1 << lab[bot_winners[u]])
        if w == v:
            loops.append(tuple(word))
            winner_masks.append(win)
            if label_maps:
                getter = bot_getters[u]
                maps.append(lab if getter is None else getter(lab))
        elif w > v:
            if length + 1 + dist[w] <= max_length:
                getter = bot_getters[u]
                search(w, length + 1, lab if getter is None else getter(lab), win)
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


def _compose(r1, r2):
    r"""
    Return the composition ``r1∘r2`` of the tuples ``r1`` and ``r2``, the
    tuple ``(r1[r2[0]], ..., r1[r2[d-1]])``, for ``d >= 2``.

    It is the second component of :func:`_semidirect_product`, used when
    the winners are not needed.

    EXAMPLES::

        sage: from surface_dynamics.interval_exchanges.rauzy_walks import _compose
        sage: _compose((1, 0, 2), (0, 2, 1))
        (1, 2, 0)
    """
    # d >= 2 for an irreducible permutation, so itemgetter returns a tuple
    return operator.itemgetter(*r2)(r1)


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
