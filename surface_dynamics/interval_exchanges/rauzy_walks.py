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
    IETFlipSequence('0 1 2 3\n1 3 0 2', 'tbtbtttbt')
    IETFlipSequence('0 1 2 3\n1 3 0 2', 'tbbtttbtb')
    IETFlipSequence('0 1 2 3\n1 3 2 0', 'bttbttbtb')
    IETFlipSequence('0 1 2 3\n1 3 2 0', 'btbttbtb')
    IETFlipSequence('0 1 2 3\n1 3 2 0', 'btbttbbtb')
    IETFlipSequence('0 1 2 3\n1 3 2 0', 'bbtbtbttb')
    IETFlipSequence('0 1 2 3\n1 3 2 0', 'bbtbbtbtb')
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
from bisect import bisect_left

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
      transformation (not flipped) on at least two letters

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
            sage: LazyReducedRauzyDiagram(iet.Permutation('a', 'a'))
            Traceback (most recent call last):
            ...
            ValueError: p must have at least two letters
        """
        if not isinstance(p, LabelledPermutationIET) or isinstance(p, FlippedPermutationIET):
            raise ValueError('p must be a labelled permutation of interval exchange transformation')
        if not p.is_irreducible():
            raise ValueError('p must be irreducible')
        # a permutation of one letter has no Rauzy move, hence no Rauzy
        # diagram; the moves below are applied without check, and the edge
        # data (itemgetter of a row) assumes at least two letters
        if len(p) < 2:
            raise ValueError('p must have at least two letters')
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
    rotations starting there (``t < b``); it is closed, as by
    :meth:`IETFlipSequence.close`: it carries the relabelling that maps
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
    ``max_length``, each node iterating only over the loops that fit in the
    length left. The loops themselves are found by a depth-first search
    pruned by the distance back to `v`, which also accumulates the winners
    of each loop and the relabelling at which it ends, so that no loop is
    replayed.

    The relabelling that closes a walk is the inverse of the product of the
    relabellings of its loops, so no flip sequence is replayed either. When
    the relabelling group at `v` (below) is trivial, every walk at `v`
    closes with the identity.

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
        IETFlipSequence('0 1 2 3\n1 3 2 0', 'btbttbtb')
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


def _diagram_data(G):
    r"""
    Return the reduced Rauzy diagram ``G`` on integer vertices, with the data
    of its edges.

    The vertex ``j`` is the ``j``-th key of ``G.vertices()``, so that the
    integers compare as the keys. The canonical labelling of each vertex is
    moved once per edge here, instead of once per edge of each walk or loop
    found later.

    OUTPUT: a tuple ``(vertices, successors, predecessors, edges)``:

    - ``vertices`` -- the sorted list of the keys of the vertices

    - ``successors[u]`` -- the pair of the vertices reached from ``u`` by the
      top and by the bottom move

    - ``predecessors[u]`` -- the list of the vertices with an edge to ``u``,
      once per edge

    - ``edges[u][move]`` -- a pair ``(a, g)``: ``a`` is the winner of the
      move in the canonical labelling of ``u``, and ``g`` the
      :func:`operator.itemgetter` of the top row ``m`` of the labelled
      permutation reached by the move from this labelling. This permutation
      is the canonical labelling of the successor with each letter ``b``
      replaced by ``m[b]``; so from the canonical labelling of ``u`` relabelled
      by ``lab`` (each letter ``c`` replaced by ``lab[c]``), the move leads to
      the canonical labelling of the successor relabelled by ``g(lab)``.

    EXAMPLES::

        sage: from surface_dynamics import iet
        sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram, _diagram_data
        sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
        sage: vertices, successors, predecessors, edges = _diagram_data(G)
        sage: vertices
        [(1, 2, 0), (2, 0, 1), (2, 1, 0)]
        sage: successors
        [(0, 2), (2, 1), (1, 0)]
        sage: predecessors
        [[0, 2], [1, 2], [0, 1]]
        sage: a, g = edges[0][1]
        sage: a, g((0, 1, 2))
        (0, (0, 2, 1))

    The winners and the top rows agree with the moves of the canonical
    labellings::

        sage: all(edges[u][move][0] == G.canonical_labelling(k)._labels[move][-1]
        ....:     and edges[u][move][1](tuple(range(3))) == tuple(G.canonical_labelling(k).rauzy_move(move)._labels[0])
        ....:     for u, k in enumerate(vertices) for move in (0, 1))
        True
    """
    vertices = G.vertices()
    index = {k: j for j, k in enumerate(vertices)}
    d = G._d
    successors = []
    predecessors = [[] for _ in vertices]
    edges = []
    for u, k in enumerate(vertices):
        succ = (index[G.successor(k, 0)], index[G.successor(k, 1)])
        successors.append(succ)
        for w in succ:
            predecessors[w].append(u)
        q = G._labellings[k]
        # the top row of the canonical labelling is 0, 1, ..., d-1 and its
        # bottom row is the key: the winners are their last letters
        edges.append(((d - 1, operator.itemgetter(*q.rauzy_move(0, -1, check=False)._labels[0])),
                      (k[-1], operator.itemgetter(*q.rauzy_move(1, -1, check=False)._labels[0]))))
    return vertices, successors, predecessors, edges


def _distances_to(v, predecessors, max_length):
    r"""
    Return the list of the lengths of the shortest paths to the vertex ``v``
    whose vertices but ``v`` are all greater than ``v``.

    The vertices are the integers ``0, ..., n-1`` where ``n`` is the length
    of ``predecessors``. The entry ``u`` of the list is the length of such a
    path from ``u``, for the vertices ``u > v`` at distance at most
    ``max_length``, and ``None`` for the other vertices.

    EXAMPLES::

        sage: from surface_dynamics.interval_exchanges.rauzy_walks import _distances_to
        sage: pred = [[1, 2], [0, 2], [1]]
        sage: _distances_to(0, pred, 5)
        [None, 1, 1]
        sage: _distances_to(1, pred, 5)
        [None, None, 1]
        sage: _distances_to(1, [[1], [2], [3], [0]], 5)
        [None, None, 1, 2]
        sage: _distances_to(1, [[1], [2], [3], [0]], 1)
        [None, None, 1, None]
    """
    dist = [None] * len(predecessors)
    frontier = [v]
    level = 0
    while frontier and level < max_length:
        level += 1
        new = []
        for w in frontier:
            for u in predecessors[w]:
                if u > v and dist[u] is None:
                    dist[u] = level
                    new.append(u)
        frontier = new
    return dist


def _first_return_loops(d, v, successors, edges, dist, max_length, label_maps):
    r"""
    Return the first-return loops at the vertex ``v`` of length at most
    ``max_length``, with their winners and the relabellings at which they
    end.

    The first-return loops at ``v`` are the words of moves from ``v`` back
    to ``v`` whose intermediate vertices are all greater than ``v``. ``d``
    is the number of letters, ``successors`` and ``edges`` are as in
    :func:`_diagram_data`, and ``dist`` is the output of :func:`_distances_to`.

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
        sage: from surface_dynamics.interval_exchanges.rauzy_walks import LazyReducedRauzyDiagram, _diagram_data, _distances_to, _first_return_loops
        sage: G = LazyReducedRauzyDiagram(iet.Permutation('a b c', 'c b a'))
        sage: vertices, successors, predecessors, edges = _diagram_data(G)
        sage: dist = _distances_to(0, predecessors, 4)
        sage: loops, winner_masks, maps = _first_return_loops(3, 0, successors, edges, dist, 4, True)
        sage: loops
        [(0,), (1, 0, 0, 1), (1, 1)]
        sage: winner_masks
        [4, 3, 1]
        sage: _first_return_loops(3, 0, successors, edges, dist, 4, False)[2] is None
        True

    The winners and the relabellings agree with those of the flip sequence
    of each loop, including when the relabelling group is nontrivial::

        sage: from surface_dynamics import Stratum
        sage: from surface_dynamics.interval_exchanges.flip_sequence import IETFlipSequence
        sage: def replay(G, v, L):
        ....:     vertices, successors, predecessors, edges = _diagram_data(G)
        ....:     loops, winner_masks, maps = _first_return_loops(G._d, v, successors, edges, _distances_to(v, predecessors, L), L, True)
        ....:     flips = [IETFlipSequence(G.canonical_labelling(vertices[v]), loop) for loop in loops]
        ....:     return (len(loops), sum(m != tuple(range(G._d)) for m in maps),
        ....:             winner_masks == [sum(1 << a for a in set(f.winners())) for f in flips],
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
        su = successors[u]
        eu = edges[u]
        for move in (0, 1):
            w = su[move]
            a, g = eu[move]
            word.append(move)
            win = winners | (1 << lab[a])
            if w == v:
                loops.append(tuple(word))
                winner_masks.append(win)
                if label_maps:
                    maps.append(g(lab))
            elif w > v:
                dw = dist[w]
                if dw is not None and length + 1 + dw <= max_length:
                    search(w, length + 1, g(lab), win)
            word.pop()

    try:
        if max_length >= 1:
            search(v, 0, tuple(range(d)), 0)
    finally:
        # search refers to itself through its closure: without this del, the
        # cycle keeps the loops alive until the cyclic garbage collector runs
        del search
    return loops, winner_masks, maps


def _loop_necklaces(lengths, max_length, min_length, primitive, orbit_masks=None, num_orbits=0, winner_masks=None, label_maps=None):
    r"""
    Iterate over the necklaces (Lyndon words if ``primitive``) over the
    alphabet ``0, 1, ..., len(lengths) - 1`` whose weight is between
    ``min_length`` and ``max_length``, the letter ``i`` having weight
    ``lengths[i]``.

    The words are generated by the prenecklace backtracking of Fredricksen,
    Kessler and Maiorana, the subtrees of weight greater than ``max_length``
    being cut. A node of weight ``w`` iterates only over the letters of
    weight at most ``max_length - w``, from an index built once: most
    letters do not fit at the deep nodes, where most of the nodes are.

    Each word comes in a pair ``(word, relabelling)``. The list ``word`` is
    modified afterwards: copy it to keep it. If ``label_maps`` is given,
    ``relabelling`` is the product of the label maps of the letters of the
    word (see below); otherwise it is ``None``.

    If ``orbit_masks`` is given, only the complete words are yielded: the
    letter ``i`` is a loop whose winners hit the orbits of the bitmask
    ``orbit_masks[i]`` among the ``num_orbits`` orbits of the relabelling
    group, and a word is kept when the union of these bitmasks is full. A
    subtree is cut when more orbits are missing than moves are left, each
    move having one winner. If moreover ``label_maps`` is given (the
    relabelling group is nontrivial), the completeness of a word is decided
    exactly from ``winner_masks`` and ``label_maps`` (see
    :func:`_first_return_loops`): the winners of the letter ``j`` of the
    word are transported by the product of the label maps of the letters
    before it, and they must hit every orbit of the product of all the label
    maps.

    EXAMPLES::

        sage: from surface_dynamics.interval_exchanges.rauzy_walks import _loop_necklaces
        sage: [tuple(w) for w, _ in _loop_necklaces([1, 1], 4, 1, True)]
        [(0,), (0, 0, 0, 1), (0, 0, 1), (0, 0, 1, 1), (0, 1), (0, 1, 1), (0, 1, 1, 1), (1,)]
        sage: [tuple(w) for w, _ in _loop_necklaces([1, 1], 4, 4, False)]
        [(0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 1, 1), (0, 1, 0, 1), (0, 1, 1, 1), (1, 1, 1, 1)]
        sage: [tuple(w) for w, _ in _loop_necklaces([2, 1], 4, 1, True)]
        [(0,), (0, 1), (0, 1, 1), (1,)]

    The letters need not be sorted by weight::

        sage: [tuple(w) for w, _ in _loop_necklaces([1, 3, 1], 3, 1, True)]
        [(0,), (0, 0, 2), (0, 2), (0, 2, 2), (1,), (2,)]

    With label maps, the product of those of the word comes with it::

        sage: [(tuple(w), r) for w, r in _loop_necklaces([1, 2], 3, 1, True, label_maps=[(1, 0, 2), (0, 2, 1)])]
        [((0,), (1, 0, 2)), ((0, 1), (1, 2, 0)), ((1,), (0, 2, 1))]

    With two orbits, each loop hitting one of them::

        sage: [tuple(w) for w, _ in _loop_necklaces([1, 1], 4, 1, True, [0b01, 0b10], 2)]
        [(0, 0, 0, 1), (0, 0, 1), (0, 0, 1, 1), (0, 1), (0, 1, 1), (0, 1, 1, 1)]

    A relabelling group of order `2` exchanging the letters ``0`` and
    ``1``. The first loop has winner ``1`` and ends with the letters
    exchanged, the second one has winner ``0`` and ends with the identity.
    The word ``(1,)`` hits the only orbit of the group but closes with the
    identity and misses the letter ``1``; the word ``(0, 0)``, not a
    Lyndon word, has winners ``1`` then ``0``::

        sage: [tuple(w) for w, _ in _loop_necklaces([1, 1], 3, 1, True, [0b1, 0b1], 1, [0b10, 0b01], [(1, 0), (0, 1)])]
        [(0,), (0, 0, 1), (0, 1), (0, 1, 1)]
        sage: [tuple(w) for w, _ in _loop_necklaces([1, 1], 2, 2, False, [0b1, 0b1], 1, [0b10, 0b01], [(1, 0), (0, 1)])]
        [(0, 0), (0, 1)]
    """
    n = len(lengths)
    if not n:
        return
    seq = []

    # by_budget[r]: the increasing list of the letters of weight at most r;
    # the candidates of a node are those of the full scan that fit, in the
    # same order
    by_budget = [[] for _ in range(max_length + 1)]
    for i, length in enumerate(lengths):
        for r in range(length, max_length + 1):
            by_budget[r].append(i)

    if label_maps is not None:
        d = len(label_maps[0])
        identity = tuple(range(d))
        # getters[i](r) is tuple(r[a] for a in label_maps[i]); d >= 2 for
        # an irreducible permutation, so itemgetter returns a tuple
        getters = [operator.itemgetter(*m) for m in label_maps]

    if orbit_masks is None:
        if label_maps is None:
            def extend(p, total):
                # seq is a prenecklace of length m with period p (the length
                # of its longest Lyndon prefix)
                m = len(seq)
                if m and total >= min_length and (p == m if primitive else m % p == 0):
                    yield seq, None
                first = seq[m - p] if m else 0
                cand = by_budget[max_length - total]
                for i in cand[bisect_left(cand, first):]:
                    seq.append(i)
                    yield from extend(p if m and i == first else m + 1, total + lengths[i])
                    seq.pop()

            try:
                yield from extend(0, 0)
            finally:
                # extend refers to itself through its closure: break the cycle
                # so that the index is freed with this generator
                del extend
            return

        def extend(p, total, relabelling):
            m = len(seq)
            if m and total >= min_length and (p == m if primitive else m % p == 0):
                yield seq, relabelling
            first = seq[m - p] if m else 0
            cand = by_budget[max_length - total]
            for i in cand[bisect_left(cand, first):]:
                seq.append(i)
                yield from extend(p if m and i == first else m + 1, total + lengths[i],
                                  getters[i](relabelling))
                seq.pop()

        try:
            yield from extend(0, 0, identity)
        finally:
            del extend
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
                yield seq, None
            first = seq[m - p] if m else 0
            cand = by_budget[max_length - total]
            for i in cand[bisect_left(cand, first):]:
                t = total + lengths[i]
                h = hit | orbit_masks[i]
                if (full ^ h).bit_count() > max_length - t:
                    continue
                seq.append(i)
                yield from extend(p if m and i == first else m + 1, t, h)
                seq.pop()

        try:
            yield from extend(0, 0, 0)
        finally:
            del extend
        return

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
            yield seq, relabelling
        first = seq[m - p] if m else 0
        cand = by_budget[max_length - total]
        for i in cand[bisect_left(cand, first):]:
            t = total + lengths[i]
            h = hit | orbit_masks[i]
            if (full ^ h).bit_count() > max_length - t:
                continue
            w = winners
            for a in winner_letters[i]:
                w |= 1 << relabelling[a]
            seq.append(i)
            yield from extend(p if m and i == first else m + 1, t, h, w, getters[i](relabelling))
            seq.pop()

    try:
        yield from extend(0, 0, 0, 0, identity)
    finally:
        del extend


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

    The walks are closed, with the relabelling that :meth:`IETFlipSequence.close`
    gives::

        sage: from surface_dynamics import Stratum
        sage: from surface_dynamics.interval_exchanges.flip_sequence import IETFlipSequence
        sage: G = LazyReducedRauzyDiagram(Stratum([2, 2]).odd_component().permutation_representative(reduced=False))
        sage: walks = list(_closed_walks(G, 9, 1, False, True))
        sage: len(walks), sum(1 for f in walks if f._relabelling != list(range(7)))
        (1165, 38)
        sage: def closed(f):
        ....:     g = IETFlipSequence(f.start(), f._simplified_flip_sequence())
        ....:     g.close()
        ....:     return g
        sage: all(f == closed(f) and f.is_closed() for f in walks)
        True
        sage: all(f._check() is None for f in walks)
        True
    """
    vertices, successors, predecessors, edges = _diagram_data(G)
    d = G._d
    identity = list(range(d))

    for v, key in enumerate(vertices):
        dist = _distances_to(v, predecessors, max_length)
        orbits = G.relabelling_orbits(key)
        # the group is trivial exactly when its orbits are singletons; then
        # every walk at v closes with the identity, and no label map is needed
        nontrivial = len(orbits) < d
        loops, winner_masks, label_maps = _first_return_loops(d, v, successors, edges, dist, max_length, nontrivial)
        lengths = [len(loop) for loop in loops]
        if complete:
            orbit_of = [None] * d
            for i, orbit in enumerate(orbits):
                for a in orbit:
                    orbit_of[a] = i
            orbit_masks = []
            for winners in winner_masks:
                hit = 0
                a = 0
                while winners:
                    if winners & 1:
                        hit |= 1 << orbit_of[a]
                    winners >>= 1
                    a += 1
                orbit_masks.append(hit)
            necklaces = _loop_necklaces(lengths, max_length, min_length, primitive,
                                        orbit_masks, len(orbits), winner_masks, label_maps)
        else:
            necklaces = _loop_necklaces(lengths, max_length, min_length, primitive, label_maps=label_maps)
        for seq, product in necklaces:
            # the walk ends at the canonical labelling of v with each letter
            # a replaced by product[a]: the inverse of product closes it
            if product is None:
                relabelling = identity[:]
            else:
                relabelling = [None] * d
                for a, b in enumerate(product):
                    relabelling[b] = a
            moves = [(move, -1) for i in seq for move in loops[i]]
            yield IETFlipSequence._from_data(G.canonical_labelling(key), G.canonical_labelling(key),
                                             moves, relabelling)
