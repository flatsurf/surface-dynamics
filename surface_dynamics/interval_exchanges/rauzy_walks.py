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


class LazyReducedRauzyDiagram:
    r"""
    The reduced Rauzy diagram of an Abelian permutation or of a linear
    involution, built on demand.

    The diagram holds, for each reduced vertex met so far, its key, its
    canonical labelling and its successors by the right Rauzy moves already
    requested. It grows by one entry per reduced vertex visited, up to the
    size of the reduced Rauzy class (7 vertices for `H(2)`, 509 for the
    representative of left degree 3 of `H_3(3, 1)`; 73 for the
    representative of `Q_2(2^2)^{hyp}`, 380 for `Q_2(1^4)^{hyp}` and 440 for
    `Q_2(5, -1)`). It is the graph a search walks through, not a store of
    results.

    The vertices are represented by their keys and canonical labellings, as
    described in :mod:`~surface_dynamics.interval_exchanges.rauzy_walks`. For
    a linear involution, a Rauzy move may be undefined: :meth:`successor`
    returns ``None`` for it, and it is no edge of the diagram.

    Once :meth:`relabelling_generators` or :meth:`relabelling_orbits` has
    been called, the diagram also holds the relabelling group of the class:
    one labelling per reduced vertex (along a spanning tree from the start),
    one list of generators at the start, at most one per edge of the
    diagram, and the orbits of the group at the start on the letters, at
    most one per letter. They are computed once for the whole class: the
    generators and orbits at another vertex are their images by the tree
    labelling of that vertex.

    :meth:`clear` drops all of it but the start.

    The canonical labellings use the alphabet of ``p``, in its order.

    INPUT:

    - ``p`` -- an irreducible permutation of Abelian differential, or an
      irreducible linear involution without flips, labelled or reduced

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

    A linear involution, whose key is the pair of the rows of its canonical
    labelling, the letters being numbered by first appearance::

        sage: q = iet.GeneralizedPermutation('b b a', 'a c c', alphabet='abc')
        sage: G = LazyReducedRauzyDiagram(q)
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
            ValueError: p must be a permutation of interval exchange transformation or a linear involution, without flips
            sage: LazyReducedRauzyDiagram(iet.GeneralizedPermutation('a a b', 'b c c', flips='a'))
            Traceback (most recent call last):
            ...
            ValueError: p must be a permutation of interval exchange transformation or a linear involution, without flips
            sage: LazyReducedRauzyDiagram(iet.Permutation('a b c', 'a c b'))
            Traceback (most recent call last):
            ...
            ValueError: p must be irreducible
            sage: LazyReducedRauzyDiagram(iet.GeneralizedPermutation('a b c a', 'b d d c'))
            Traceback (most recent call last):
            ...
            ValueError: p must be irreducible
            sage: LazyReducedRauzyDiagram(iet.Permutation('a', 'a'))
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
        self.clear()

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
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: LazyReducedRauzyDiagram._key(iet.Permutation('a c b d', 'd a b c'))
            (3, 0, 2, 1)
            sage: LazyReducedRauzyDiagram._key(iet.Permutation('a c b d', 'd a b c', reduced=True))
            (3, 0, 2, 1)
            sage: LazyReducedRauzyDiagram._key(iet.GeneralizedPermutation('c a c', 'a d b b d'))
            ((0, 1, 0), (1, 2, 3, 3, 2))
            sage: LazyReducedRauzyDiagram._key(iet.GeneralizedPermutation('c a c', 'a d b b d', reduced=True))
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
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: LazyReducedRauzyDiagram._move_map(iet.Permutation('a c b d', 'd a b c', alphabet='abcd'))
            [0, 2, 1, 3]
            sage: q = iet.GeneralizedPermutation('a b a', 'b c d d c').rauzy_move(0)
            sage: q
            c a b a
            b c d d
            sage: LazyReducedRauzyDiagram._move_map(q)
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

    def _canonical(self, key):
        r"""
        Return the canonical labelling of the vertex ``key``, as a new
        permutation or linear involution, over the alphabet of the input.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
            sage: G._canonical((2, 0, 1))
            a b c
            c a b
            sage: G = LazyReducedRauzyDiagram(iet.GeneralizedPermutation('a a b', 'b c c'))
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
        self._root_orbits = None

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
        right Rauzy move ``move`` (``0`` for top, ``1`` for bottom), or
        ``None`` if this move is undefined (for a linear involution only).

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

        A move of a linear involution may be undefined::

            sage: G = LazyReducedRauzyDiagram(iet.GeneralizedPermutation('a a b', 'b c c'))
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
        succ = self._successors[key]
        s = succ[move]
        if s is None:
            q = self._labellings[key]
            if self._quadratic and not q.has_rauzy_move(move):
                # False marks a move known to be undefined, None one not
                # computed yet
                succ[move] = False
                return None
            s = self._key(q.rauzy_move(move, -1, check=False))
            if s not in self._labellings:
                self._labellings[s] = self._canonical(s)
                self._successors[s] = [None, None]
            succ[move] = s
        elif s is False:
            return None
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

        For a linear involution, the undefined moves are no edges::

            sage: from surface_dynamics import Stratum
            sage: C = Stratum([2, 2], k=2).hyperelliptic_component()
            sage: G = LazyReducedRauzyDiagram(C.permutation_representative(reduced=False))
            sage: V = G.vertices()
            sage: len(V)
            73
            sage: V[0]
            ((0, 1, 0), (1, 2, 3, 2, 4, 3, 4))
            sage: sum(G.successor(k, move) is None for k in V for move in (0, 1))
            26
        """
        if not self._explored:
            seen = {self._start}
            todo = [self._start]
            while todo:
                k = todo.pop()
                for move in (0, 1):
                    s = self.successor(k, move)
                    if s is not None and s not in seen:
                        seen.add(s)
                        todo.append(s)
            self._explored = True
        return sorted(self._labellings)

    def _relabelling_group(self):
        r"""
        Compute, once for the class, the labelling of each vertex along a
        spanning tree from the start, the generators of the relabelling
        group at the start and its orbits on the letters.

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
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
            sage: G._relabelling_group()
            sage: sorted(G._tree_labellings.items())
            [((1, 2, 0), (0, 2, 1)), ((2, 0, 1), (0, 1, 2)), ((2, 1, 0), (0, 1, 2))]
            sage: G._root_generators
            []
            sage: G._root_orbits
            [(0,), (1,), (2,)]

        For a linear involution, the undefined moves are skipped; the group
        of `Q_2(2^2)^{hyp}` acts freely, with order the number of labelled
        vertices over the number of reduced ones, `8760 / 73 = 120`::

            sage: from surface_dynamics import Stratum
            sage: C = Stratum([2, 2], k=2).hyperelliptic_component()
            sage: G = LazyReducedRauzyDiagram(C.permutation_representative(reduced=False))
            sage: G._relabelling_group()
            sage: len(G._tree_labellings)
            73
            sage: gens = [[i + 1 for i in g] for g in G._root_generators]
            sage: PermutationGroup(gens).order()
            120
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
                    if w is None:
                        continue
                    # the move from the canonical labelling of u ends at the
                    # canonical labelling of w relabelled by mu
                    mu = self._move_map(self._labellings[u].rauzy_move(move, -1, check=False))
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
        self._tree_labellings = tree
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

        The relabelling group of a linear involution may be transitive on the
        letters, as for `Q_2(2^2)^{hyp}`::

            sage: G = LazyReducedRauzyDiagram(Stratum([2, 2], k=2).hyperelliptic_component().permutation_representative(reduced=False))
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
            sage: G = LazyReducedRauzyDiagram(Stratum([2, 2]).odd_component().permutation_representative(reduced=False))
            sage: len(G.vertices())
            294
            sage: all(G.relabelling_orbits(k) == orbits(G, k) for k in G.vertices())
            True
            sage: G = LazyReducedRauzyDiagram(Stratum([5, -1], k=2).unique_component().permutation_representative(reduced=False))
            sage: len(G.vertices())
            440
            sage: all(G.relabelling_orbits(k) == orbits(G, k) for k in G.vertices())
            True
        """
        self._relabelling_group()
        t = self._tree_labellings[key]
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

        The same on a linear involution, whose relabelling group is never
        trivial; the walks of `Q_2(2^2)^{hyp}`::

            sage: G = LazyReducedRauzyDiagram(Stratum([2, 2], k=2).hyperelliptic_component().permutation_representative(reduced=False))
            sage: walks = list(G.closed_walks(10, complete=False))
            sage: len(walks)
            145
            sage: all(f == closed(f) and f.is_closed() and f._check() is None for f in walks)
            True
            sage: [f for f in G.closed_walks(8)]
            [RauzyFlipSequence('0 1 0 2 3\n3 1 4 2 4', 'tbttbbtb', relabelling=(0,1)),
             RauzyFlipSequence('0 1 2 1\n3 4 3 2 0 4', 'tbtbbttb', relabelling=(3,4))]

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
        top_winners, bot_winners, top_predecessors, bot_predecessors,
        top_maps, bot_maps)``, where ``d`` is the number of letters and ``n``
        the number of vertices:

        - ``vertices`` -- the sorted list of the keys of the vertices

        - ``top_successors[u]``, ``bot_successors[u]`` -- the vertex reached
          from ``u`` by the top move (``0``) and by the bottom move (``1``),
          ``-1`` if the move is undefined (for a linear involution only)

        - ``top_winners[u]``, ``bot_winners[u]`` -- the winner of the top
          (bottom) move in the canonical labelling of ``u``, its last letter
          on the top (bottom) row, as an index in ``0, ..., d-1``. The top
          winner of the canonical labelling of a permutation is always
          ``d - 1``, not that of a linear involution

        - ``top_predecessors[w]``, ``bot_predecessors[w]`` -- the vertex with
          a top (bottom) edge to ``w``, ``-1`` if there is none (for a linear
          involution only). It is unique: a Rauzy move is invertible on
          labelled permutations and commutes with the relabellings, so two
          vertices with an edge of the same label to ``w`` have canonical
          labellings that differ by a relabelling, hence are equal

        - ``top_maps[d*u + b]``, ``bot_maps[d*u + b]`` -- the entry ``b`` of
          the relabelling ``m`` (see :meth:`_move_map`) of the labelled
          permutation reached by the top (bottom) move from the canonical
          labelling of ``u``: this permutation is the canonical labelling of
          the successor with each letter ``b`` replaced by ``m[b]``. For a
          permutation, ``m`` is its top row. So from the canonical labelling of
          ``u`` relabelled by ``lab`` (each letter ``c`` replaced by
          ``lab[c]``), the move leads to the canonical labelling of the
          successor relabelled by the tuple ``operator.itemgetter(*m)(lab)``.
          The entries of an undefined move are those of the identity

        The entries ``-1`` compare below every vertex and are never used as an
        index by the search: an edge ``u -> w`` is followed only when ``w``
        is ``v`` or greater than ``v``, at a base ``v >= 0``, and a
        predecessor only when it is greater than ``v``.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram
            sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
            sage: (vertices, top_succ, bot_succ, top_winners, bot_winners,
            ....:  top_pred, bot_pred, top_maps, bot_maps) = G._diagram_data()
            sage: vertices
            [(1, 2, 0), (2, 0, 1), (2, 1, 0)]
            sage: top_succ, bot_succ
            ([0, 2, 1], [2, 1, 0])
            sage: top_winners, bot_winners
            ([2, 2, 2], [0, 1, 0])
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
            sage: vertices, top_succ, bot_succ, _, _, top_pred, bot_pred, _, _ = G._diagram_data()
            sage: all(top_succ[top_pred[w]] == w and bot_succ[bot_pred[w]] == w for w in range(len(vertices)))
            True

        For a linear involution, the undefined moves have the successor
        ``-1``; the winners and the relabellings agree with the moves of the
        canonical labellings, and each vertex is the successor of its
        predecessors::

            sage: G = LazyReducedRauzyDiagram(Stratum([2, 2], k=2).hyperelliptic_component().permutation_representative(reduced=False))
            sage: (vertices, top_succ, bot_succ, top_winners, bot_winners,
            ....:  top_pred, bot_pred, top_maps, bot_maps) = G._diagram_data()
            sage: n, d = len(vertices), 5
            sage: n, top_succ.count(-1), bot_succ.count(-1)
            (73, 13, 13)
            sage: def agree(u, move, succ, winners, maps):
            ....:     q = G.canonical_labelling(vertices[u])
            ....:     if not q.has_rauzy_move(move):
            ....:         return succ[u] == -1 and maps[d*u:d*u + d] == list(range(d))
            ....:     r = q.rauzy_move(move)
            ....:     m = maps[d*u:d*u + d]
            ....:     s = G.canonical_labelling(vertices[succ[u]])
            ....:     return (winners[u] == q._labels[move][-1]
            ....:             and [[m[b] for b in row] for row in s._labels] == r._labels)
            sage: all(agree(u, 0, top_succ, top_winners, top_maps) and agree(u, 1, bot_succ, bot_winners, bot_maps)
            ....:     for u in range(n))
            True
            sage: all(top_succ[top_pred[w]] == w for w in range(n) if top_pred[w] != -1)
            True
            sage: all(bot_succ[bot_pred[w]] == w for w in range(n) if bot_pred[w] != -1)
            True
            sage: top_pred.count(-1), bot_pred.count(-1)
            (13, 13)
        """
        vertices = self.vertices()
        index = {k: j for j, k in enumerate(vertices)}
        # an undefined move has the successor -1
        index[None] = -1
        n = len(vertices)
        d = self._d
        identity = list(range(d))
        quadratic = self._quadratic
        top_successors = []
        bot_successors = []
        top_winners = []
        bot_winners = []
        top_predecessors = [-1] * n
        bot_predecessors = [-1] * n
        top_maps = []
        bot_maps = []
        for u, k in enumerate(vertices):
            t = index[self.successor(k, 0)]
            b = index[self.successor(k, 1)]
            top_successors.append(t)
            bot_successors.append(b)
            q = self._labellings[k]
            if quadratic:
                # the key holds the two rows of the canonical labelling
                top_winners.append(k[0][-1])
                bot_winners.append(k[1][-1])
            else:
                # the top row of the canonical labelling is 0, 1, ..., d-1 and
                # its bottom row is the key
                top_winners.append(d - 1)
                bot_winners.append(k[-1])
            if t == -1:
                top_maps.extend(identity)
            else:
                top_predecessors[t] = u
                top_maps.extend(self._move_map(q.rauzy_move(0, -1, check=False)))
            if b == -1:
                bot_maps.extend(identity)
            else:
                bot_predecessors[b] = u
                bot_maps.extend(self._move_map(q.rauzy_move(1, -1, check=False)))
        return (vertices, top_successors, bot_successors, top_winners, bot_winners,
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
        (vertices, top_successors, bot_successors, top_winners, bot_winners,
         top_predecessors, bot_predecessors, top_maps, bot_maps) = self._diagram_data()
        d = self._d
        n = len(vertices)
        identity = list(range(d))
        # built once for the whole search: the loops of each vertex share them
        top_getters = [operator.itemgetter(*top_maps[d*u:d*u + d]) for u in range(n)]
        bot_getters = [operator.itemgetter(*bot_maps[d*u:d*u + d]) for u in range(n)]

        # the lower bounds of the distances back to the base, shared by the
        # bases: each base sets the entries it reaches and resets them after
        h = _backward_depth(max_length)
        dist = [h + 1] * n

        for v, key in enumerate(vertices):
            reached = _distances_to(v, top_predecessors, bot_predecessors, h, dist)
            orbits = self.relabelling_orbits(key)
            # the group is trivial exactly when its orbits are singletons; then
            # every walk at v closes with the identity, and no label map is needed
            nontrivial = len(orbits) < d
            loops, winner_masks, label_maps = _first_return_loops(
                d, v, top_successors, bot_successors, top_winners, bot_winners,
                top_getters, bot_getters, dist, max_length, nontrivial)
            # the next base needs every entry back at h + 1; resetting the
            # reached entries only costs no O(n) per base
            for u in reached:
                dist[u] = h + 1
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
            # the walks at v share one start: no method of RauzyFlipSequence
            # mutates its start in place, while the end is moved, relabelled
            # and inverted in place, so each walk owns its end; the start is
            # copied at the first walk only, since a vertex may have none (no
            # vertex of Q_3(1^8) has a complete walk of length at most 14)
            start = None
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
                    if start is None:
                        start = canonical_labelling(key)
                    yield RauzyFlipSequence._from_data(start, canonical_labelling(key),
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
                    yield RauzyFlipSequence._from_data(start, canonical_labelling(key),
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
                yield RauzyFlipSequence._from_data(start, canonical_labelling(key),
                                                 moves, relabelling)


def _backward_depth(max_length):
    r"""
    Return the depth at which the search of closed walks of length at most
    ``max_length`` stops the backward search of the distances to each base.

    A first-return loop of length at most ``max_length`` at the base ``v``
    is found by a forward search from ``v``, pruned by lower bounds of the
    distances back to ``v`` (see :func:`_distances_to` and
    :func:`_first_return_loops`). These bounds are exact up to the depth
    returned, and the vertices beyond it count as farther: the pruning is
    exact on the last ``max_length // 2 + 2`` moves of each loop, while the
    backward search visits about as few vertices as the forward one. A depth
    above ``max_length`` would bound nothing more.

    EXAMPLES::

        sage: from surface_dynamics.interval_exchanges.rauzy_walks import _backward_depth
        sage: [_backward_depth(L) for L in range(8)]
        [0, 1, 2, 3, 4, 4, 5, 5]
        sage: _backward_depth(14), _backward_depth(24)
        (9, 14)
    """
    # complete walks, cold process, time (s) at the depths L//2 - 1, L//2,
    # L//2 + 1, L//2 + 2, L:
    #   H_3(4)^odd 22:      0.81  0.80  0.79  0.81  0.79
    #   H_3(3,1) 22:        1.96  1.95  1.96  1.90  1.89
    #   H_4(6)^odd 24:     23.13 22.36 22.07 22.05 22.04
    #   H_4(5,1) 20:       18.10 16.46 15.71 15.48 18.58
    #   H_4(3,2,1) 16:     13.86 11.90 11.01 10.70 15.40
    #   Q_3(7,1) 24:       14.06 13.70 13.72 13.41 16.43
    #   Q_3(1^8) 14:       68.51 65.64 66.58 63.95 97.91
    # the largest ratio over the cases to the best time of the case is 1.30,
    # 1.11, 1.04, 1.03 and 1.53 at these depths: L//2 + 2 is the best; at L
    # the backward search visits far more vertices than the forward one cuts
    # (Q_3(1^8) 14 is 53 % slower)
    return min(max_length, max_length // 2 + 2)


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
    stands for none: it is below every vertex, so it is never followed).

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
    for a relabelling ``m`` given
    by the getter of the edge: if the start is the canonical labelling with
    each letter ``c`` replaced by ``lab[c]``, ``lab`` a tuple, the end is the
    canonical labelling of the end vertex relabelled by
    ``top_getters[u](lab)`` (move ``0``) or ``bot_getters[u](lab)`` (move
    ``1``), the getters being the :func:`operator.itemgetter` of the
    relabellings ``top_maps[d*u:d*u + d]`` and ``bot_maps[d*u:d*u + d]`` of
    :meth:`LazyReducedRauzyDiagram._diagram_data`. The lists
    ``top_successors``, ``bot_successors``, ``top_winners`` and
    ``bot_winners`` are those of that method too.

    The loops are found by a depth-first search, which also accumulates
    along the current prefix its winners and the relabelling at which it
    ends, from the data of the edges: no loop is replayed.

    OUTPUT: a triple ``(loops, winner_masks, maps)``, the lists indexed by
    the loops sorted as words:

    - ``loops[j]`` -- the tuple of the moves of the loop ``j``

    - ``winner_masks[j]`` -- the bitmask of the letters of the canonical
      labelling of ``v`` that are winners of the loop ``j`` run from it.
      They are always computed, but only a search of complete walks reads
      them: :meth:`LazyReducedRauzyDiagram._closed_walks` with
      ``complete=False`` drops them, the winner bit of each edge being
      computed for nothing there

    - ``maps[j]`` -- if ``label_maps`` is set, the relabelling of the
      labelled permutation at which the loop ``j`` ends, a tuple ``m`` such
      that this
      permutation is the canonical labelling of ``v`` with each letter ``a``
      replaced by ``m[a]``; otherwise ``maps`` is ``None``

    EXAMPLES::

        sage: from surface_dynamics import iet
        sage: import operator
        sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram, _distances_to, _first_return_loops
        sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
        sage: (vertices, top_succ, bot_succ, top_winners, bot_winners,
        ....:  top_pred, bot_pred, top_maps, bot_maps) = G._diagram_data()
        sage: top_getters = [operator.itemgetter(*top_maps[3*u:3*u + 3]) for u in range(3)]
        sage: bot_getters = [operator.itemgetter(*bot_maps[3*u:3*u + 3]) for u in range(3)]
        sage: dist = [5] * 3
        sage: _ = _distances_to(0, top_pred, bot_pred, 4, dist)
        sage: loops, winner_masks, maps = _first_return_loops(3, 0, top_succ, bot_succ, top_winners, bot_winners,
        ....:                                                 top_getters, bot_getters, dist, 4, True)
        sage: loops
        [(0,), (1, 0, 0, 1), (1, 1)]
        sage: winner_masks
        [4, 3, 1]
        sage: _first_return_loops(3, 0, top_succ, bot_succ, top_winners, bot_winners,
        ....:                     top_getters, bot_getters, dist, 4, False)[2] is None
        True

    The winners and the relabellings agree with those of the flip sequence
    of each loop, including when the relabelling group is nontrivial::

        sage: from surface_dynamics import Stratum
        sage: from surface_dynamics.interval_exchanges.flip_sequence import RauzyFlipSequence
        sage: def replay(G, v, L):
        ....:     d = G._d
        ....:     (vertices, top_succ, bot_succ, top_winners, bot_winners,
        ....:      top_pred, bot_pred, top_maps, bot_maps) = G._diagram_data()
        ....:     top_getters = [operator.itemgetter(*top_maps[d*u:d*u + d]) for u in range(len(vertices))]
        ....:     bot_getters = [operator.itemgetter(*bot_maps[d*u:d*u + d]) for u in range(len(vertices))]
        ....:     dist = [L + 1] * len(vertices)
        ....:     _ = _distances_to(v, top_pred, bot_pred, L, dist)
        ....:     loops, winner_masks, maps = _first_return_loops(d, v, top_succ, bot_succ, top_winners, bot_winners,
        ....:                                                     top_getters, bot_getters, dist, L, True)
        ....:     flips = [RauzyFlipSequence(G.canonical_labelling(vertices[v]), loop) for loop in loops]
        ....:     return (len(loops), sum(m != tuple(range(G._d)) for m in maps),
        ....:             winner_masks == [sum(1 << list(G._alphabet).index(a) for a in set(f.winners())) for f in flips],
        ....:             maps == [tuple(G._move_map(f.end())) for f in flips])
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

    On a linear involution, whose top winner is not always the last letter
    and whose moves may be undefined::

        sage: G = LazyReducedRauzyDiagram(Stratum([2, 2], k=2).hyperelliptic_component().permutation_representative(reduced=False))
        sage: replay(G, 0, 12)
        (17, 17, True, True)

    The loops, their winners and their relabellings do not depend on the
    depth of the backward search::

        sage: def by_depth(G, v, L, h):
        ....:     d = G._d
        ....:     (vertices, top_succ, bot_succ, top_winners, bot_winners,
        ....:      top_pred, bot_pred, top_maps, bot_maps) = G._diagram_data()
        ....:     top_getters = [operator.itemgetter(*top_maps[d*u:d*u + d]) for u in range(len(vertices))]
        ....:     bot_getters = [operator.itemgetter(*bot_maps[d*u:d*u + d]) for u in range(len(vertices))]
        ....:     dist = [h + 1] * len(vertices)
        ....:     _ = _distances_to(v, top_pred, bot_pred, h, dist)
        ....:     return _first_return_loops(d, v, top_succ, bot_succ, top_winners, bot_winners,
        ....:                                top_getters, bot_getters, dist, L, True)
        sage: G1 = LazyReducedRauzyDiagram(Stratum([2, 2]).odd_component().permutation_representative(reduced=False))
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
                maps.append(top_getters[u](lab))
        elif w > v:
            if length + 1 + dist[w] <= max_length:
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
            if length + 1 + dist[w] <= max_length:
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
