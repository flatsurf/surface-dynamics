r"""
Flip sequence of permutation of iet

A concatenation of:
- left rauzy moves (top or bottom)
- right rauzy moves (top or bottom)
- left_right_inverse, top_bottom_inverse, symmetric
- relabeled
"""
#*****************************************************************************
#       Copyright (C) 2021 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function, absolute_import

import numbers
from copy import copy

from sage.structure.sage_object import SageObject
from sage.rings.all import ZZ
from sage.modules.free_module import FreeModule
from sage.matrix.matrix_space import MatrixSpace
import sage.matrix.matrix0
from sage.combinat.words.words import FiniteWords
from sage.combinat.words.morphism import WordMorphism

from .template import interval_conversion, side_conversion
from .labelled import LabelledPermutation
from .constructors import GeneralizedPermutation
from surface_dynamics.misc.permutation import (perm_init, perm_check,
        perm_compose, perm_one, perm_cycle_string, perm_orbit,
        perm_on_list_inplace, perm_invert)


class IETFlipSequence(SageObject):
    r"""
    A flip sequence of interval exchange transformations.

    Each step is a Rauzy move either on the left or on the right.

    EXAMPLES::

        sage: from surface_dynamics import iet

    Making substitutions from a flip sequence::

        sage: p = iet.Permutation([['a','b','c','d'],['d','c','b','a']])
        sage: g = iet.FlipSequence(p, ['t', 't', 'b', 't', 'b', 'b', 't', 'b'])
        sage: assert g.is_complete()
        sage: s1 = g.orbit_substitution()
        sage: print(s1)
        a->adbd, b->adbdbd, c->adccd, d->adcd
        sage: s2 = g.interval_substitution()
        sage: print(s2)
        a->abcd, b->bab, c->cdc, d->dcbababcd
        sage: s1.incidence_matrix() == s2.incidence_matrix().transpose()
        True

    The minimal dilatations in H(2g-2)^hyp from Boissy-Lanneau::

        sage: def gamma(n, k):
        ....:     p = iet.Permutation(list(range(n)), list(range(n-1,-1,-1)))
        ....:     for _ in range(k): p.rauzy_move('t', inplace=True)
        ....:     f = iet.FlipSequence(p, ['b'] * (n-1-k) + ['t'] * (n-1-2*k), True, True)
        ....:     f.close()
        ....:     assert f.is_complete()
        ....:     return f

        sage: x = polygen(QQ)
        sage: gamma(4, 1).matrix().charpoly()
        x^4 - x^3 - x^2 - x + 1
        sage: gamma(6, 1).matrix().charpoly()
        x^6 - x^5 - x^4 - x^3 - x^2 - x + 1
        sage: gamma(6, 2).matrix().charpoly()
        x^6 - x^5 - x^4 + x^3 - x^2 - x + 1
        sage: gamma(8, 1).matrix().charpoly()
        x^8 - x^7 - x^6 - x^5 - x^4 - x^3 - x^2 - x + 1
        sage: gamma(8, 2).matrix().charpoly()
        x^8 - x^7 - x^6 - x^5 + x^4 - x^3 - x^2 - x + 1
        sage: gamma(8, 3).matrix().charpoly()
        x^8 - x^7 - x^6 + x^5 - x^4 + x^3 - x^2 - x + 1
    """
    def __init__(self, p, rauzy_moves=None, top_bottom_inverse=False, left_right_inverse=False, relabelling=None):
        r"""
        INPUT:

        - ``p`` - permutation or generalized permutation

        - ``rauzy_moves`` - (optional; default empty) a sequence of Rauzy moves

        - ``top_bottom_inverse`` - (optional; default ``False``) whether a
          top-bottom inverse is performed at the end of the flip sequence

        - ``left_right_inverse`` - (optional; default ``False``) whether a
          left-right inverse is performed at the end of the flip sequence

        - ``relabelling`` - (optional; default identity permutation) the
          relabelling to be performed at the end of the flip sequence
        """
        if not isinstance(p, LabelledPermutation):
            p = GeneralizedPermutation(p)
        self._start = p.__copy__()
        self._end = p.__copy__()
        self._rauzy_moves = []
        self._top_bottom_inverse = False
        self._left_right_inverse = False
        self._relabelling = perm_one(len(self._start))

        if rauzy_moves is not None:
            for r in rauzy_moves:
                if isinstance(r, numbers.Integral):
                    self.rauzy_move(r)
                else:
                    self.rauzy_move(*r)

        if top_bottom_inverse:
            self.top_bottom_inverse()
        if left_right_inverse:
            self.left_right_inverse()
        if relabelling is not None:
            self.relabel(relabelling)

    def __eq__(self, other):
        if type(self) != type(other):
            raise TypeError('incomparable')

        return self._start == other._start and \
               self._rauzy_moves == other._rauzy_moves and \
               self._top_bottom_inverse == other._top_bottom_inverse and \
               self._left_right_inverse == other._left_right_inverse and \
               self._relabelling == other._relabelling

    def __ne__(self, other):
        return not (self == other)

    def _check(self):
        r"""
        TESTS::

            sage: from surface_dynamics import iet

            sage: p = iet.Permutation('a b c d e f', 'f e d c b a')
            sage: q = p.rauzy_move('t').rauzy_move('b').rauzy_move('t')
            sage: f = iet.FlipSequence(q, ['t', 't', 'b', 'b', 'b', 't', 'b', 't'],
            ....:                      top_bottom_inverse=False,
            ....:                      left_right_inverse=False)
            sage: f._check()
            sage: f.top_bottom_inverse()
            sage: f.left_right_inverse()
            sage: f.close()
            sage: f._check()
        """
        p = copy(self._start)
        for winner, side in self._rauzy_moves:
            p.rauzy_move(winner, side, inplace=True)
        if self._top_bottom_inverse:
            p.top_bottom_inverse(inplace=True)
        if self._left_right_inverse:
            p.left_right_inverse(inplace=True)
        p.relabel(self._relabelling)
        if self._end != p:
            raise ValueError('inconsistent flip sequence: end=%s while p=%s' % (self._end.str('/'), p.str('/')))

    def __copy__(self):
        res = IETFlipSequence.__new__(IETFlipSequence)
        res._start = copy(self._start)
        res._end = copy(self._end)
        res._rauzy_moves = self._rauzy_moves[:]
        res._top_bottom_inverse = self._top_bottom_inverse
        res._left_right_inverse = self._left_right_inverse
        res._relabelling = self._relabelling
        return res

    def __mul__(self, other):
        if type(self) != type(other):
            raise TypeError('incompatible types for multiplication')
        if self._end != other._start:
            raise ValueError('the left path should end with the start of the right path')

        res = IETFlipSequence.__new__(IETFlipSequence)
        res._start = copy(self._start)
        res._end = copy(other._end)
        res._rauzy_moves = self._rauzy_moves[:]

        if self._left_right_inverse and self._top_bottom_inverse:
            res._rauzy_moves.extend((1-winner, -side-1) for winner, side in other._rauzy_moves)
        elif self._left_right_inverse and not self._top_bottom_inverse:
            res._rauzy_moves.extend((winner, -side-1) for winner, side in other._rauzy_moves)
        elif not self._left_right_inverse and self._top_bottom_inverse:
            res._rauzy_moves.extend((1-winner, side) for winner, side in other._rauzy_moves)
        else:
            assert not self._left_right_inverse
            assert not self._top_bottom_inverse
            res._rauzy_moves.extend((winner, side) for winner, side in other._rauzy_moves)

        res._top_bottom_inverse = self._top_bottom_inverse ^ other._top_bottom_inverse
        res._left_right_inverse = self._left_right_inverse ^ other._left_right_inverse
        res._relabelling = perm_compose(self._relabelling, other._relabelling)
        return res

    def __iter__(self):
        p = copy(self._start)
        for winner, side in self._rauzy_moves:
            yield (p, winner, side)
            p.rauzy_move(winner, side, inplace=True)

    def __reversed__(self):
        r"""
        TESTS::

            sage: import itertools
            sage: from surface_dynamics import iet
            sage: p = iet.Permutation('a b c d e','d a e b c')
            sage: for tb, lr in itertools.product([False, True], repeat=2):
            ....:     f = iet.FlipSequence(p, ['t','t','b','t','b'],  tb, lr, "(0,1,3)")
            ....:     l1 = [(copy(p), winner, side) for p, winner, side in f]
            ....:     l2 = [(copy(p), winner, side) for p, winner, side in reversed(f)]
            ....:     assert l2 == l1[::-1]
        """
        p = copy(self._end)
        if self._top_bottom_inverse:
            p.top_bottom_inverse(inplace=True)
        if self._left_right_inverse:
            p.left_right_inverse(inplace=True)
        p.relabel(perm_invert(self._relabelling))
        for winner, side in reversed(self._rauzy_moves):
            p.backward_rauzy_move(winner, side, inplace=True)
            yield (p, winner, side)

    def _simplified_flip_sequence(self):
        l = []
        for winner, side in self._rauzy_moves:
            winner = 'tb'[winner]
            if side == -1:
                l.append(winner)
            else:
                l.append(winner + 'l')
        return l

    def _repr_(self):
        return "FlipSequence({!r}, {!r}, {!r}, {!r}, {!r})".format(
                self._start.str(),
                self._simplified_flip_sequence(),
                self._top_bottom_inverse,
                self._left_right_inverse,
                perm_cycle_string(self._relabelling, singletons=False))

    def rauzy_move(self, winner, side='right', iterations=1):
        winner = interval_conversion(winner)
        side = side_conversion(side)

        for _ in range(iterations):
            self._end.rauzy_move(winner, side, inplace=True)

        if self._top_bottom_inverse:
            winner = 1 - winner
        if self._left_right_inverse:
            side = 0 if side == -1 else -1
        self._rauzy_moves.extend([(winner, side)]*iterations)

    def left_right_inverse(self):
        self._left_right_inverse = True
        self._end.left_right_inverse(inplace=True)

    def top_bottom_inverse(self):
        self._top_bottom_inverse = True
        self._end.top_bottom_inverse(inplace=True)

    def find_closure(self):
        r"""
        Return the unique permutation of the labels that allows to make a
        closed loop out of this path. If no such permutation exists, return
        ``None``.

        EXAMPLES::

            sage: from surface_dynamics import iet

            sage: p = iet.Permutation('a b c d e f', 'f e d c b a')
            sage: q = p.rauzy_move('t').rauzy_move('b').rauzy_move('t')
            sage: f = iet.FlipSequence(q, ['t', 't', 'b', 'b', 'b', 't', 'b', 't'],
            ....:                      top_bottom_inverse=True,
            ....:                      left_right_inverse=True)
            sage: f.find_closure()
            [2, 3, 1, 0, 5, 4]
        """
        if self._start._twin != self._end._twin:
            return None

        l0 = self._start._labels
        l1 = self._end._labels
        p = [None] * len(self._start)
        assert len(l0[0]) == len(l1[0]) and len(l0[1]) == len(l1[1])

        for i,j in zip(l0[0] + l0[1], l1[0] + l1[1]):
            if p[j] is None:
                p[j] = i
            else:
                assert p[j] == i
        return p

    def relabel(self, p):
        if not perm_check(p, len(self._start)):
            p = perm_init(p, n=len(self._start))
        self._end.relabel(p)
        self._relabelling = perm_compose(self._relabelling, p)

    def close(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics import iet

            sage: p = iet.Permutation('a b c d e f', 'f e d c b a')
            sage: q = p.rauzy_move('t').rauzy_move('b').rauzy_move('t')
            sage: f = iet.FlipSequence(q, ['t', 't', 'b', 'b', 'b', 't', 'b', 't'],
            ....:                      top_bottom_inverse=True,
            ....:                      left_right_inverse=True)
            sage: f.close()
            sage: f
            FlipSequence('a b f c d e\nf a e b d c', ['t', 't', 'b', 'b', 'b', 't', 'b', 't'], True, True, '(0,2,1,3)(4,5)')
        """
        if self._start == self._end:
            return
        p = self.find_closure()
        if p is None:
            raise ValueError('can not close this loop with a relabelling')
        self._end.relabel(p)
        self._relabelling = perm_compose(self._relabelling, p)

    def is_loop(self):
        return self._start == self._end

    def winners_losers(self):
        r"""
        Return the list of winners.
        """
        winners = []
        losers = []
        for p, winner, side in self:
            winners.append(p[winner][side])
            losers.append(p[1-winner][side])
        return (winners, losers)

    def winners(self):
        r"""
        Return the list of winners along this flip sequence.

        EXAMPLES::

            sage: from surface_dynamics import iet

            sage: p = iet.Permutation('a b','b a')
            sage: iet.FlipSequence(p).winners()
            []
            sage: iet.FlipSequence(p, [0]).winners()
            ['b']
            sage: iet.FlipSequence(p, [1]).winners()
            ['a']
        """
        return self.winners_losers()[0]

    def losers(self):
        r"""
        Return the list of the loosers along this flip sequence.

        EXAMPLES::

            sage: from surface_dynamics import iet

            sage: p = iet.Permutation('a b c','c b a')
            sage: iet.FlipSequence(p, ['t','b','t']).losers()
            ['a', 'c', 'b']
            sage: iet.FlipSequence(p, ['b','t','b']).losers()
            ['c', 'a', 'b']
        """
        return self.winners_losers()[1]

    def is_complete(self):
        r"""
        Return whether that all winners appear.

        EXAMPLES::

            sage: from surface_dynamics import iet

            sage: p = iet.Permutation('a b c d e f', 'f e d c b a')
            sage: q = p.rauzy_move('t').rauzy_move('b').rauzy_move('t')
            sage: f = iet.FlipSequence(q, ['t', 't', 'b', 'b', 'b', 't', 'b', 't'],
            ....:                      top_bottom_inverse=True,
            ....:                      left_right_inverse=True)
            sage: f.close()
            sage: f.is_complete()
            True
        """
        if not self._start == self._end:
            raise ValueError('not a loop')
        n = len(self._start)
        moved = set()
        for p, winner, side in self:
            moved.update(perm_orbit(self._relabelling, p._labels[winner][side]))
            if len(moved) == n:
                return True
        return False

    is_full = is_complete

    def matrix(self):
        V = FreeModule(ZZ, len(self._start))
        columns = [copy(v) for v in V.basis()]
        for p, winner, side in self:
            winner_letter = p._labels[winner][side]
            loser_letter = p._labels[1-winner][side]
            columns[loser_letter] += columns[winner_letter]
        m = MatrixSpace(ZZ, len(p))(columns).transpose()
        perm_on_list_inplace(self._relabelling, m, swap=sage.matrix.matrix0.Matrix.swap_columns)
        return m

    def _int_list_to_substitutions(self, words):
        A = self._start._alphabet
        W = FiniteWords(A)
        s = {}
        for i, w in enumerate(words):
            s[A.unrank(i)] = [A.unrank(j) for j in w]
        return WordMorphism(s, domain=W, codomain=W)

    def orbit_substitution(self, as_plain_list=False):
        r"""
        Return the substitution on the orbit.

        .. TODO::

            When there is a top/bottom version one should use inverses in the
            free group.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b','b a')
            sage: g0 = iet.FlipSequence(p, ['t'])
            sage: s0 = g0.orbit_substitution()
            sage: print(s0)
            a->ab, b->b
            sage: g1 = iet.FlipSequence(p, ['b'])
            sage: s1 = g1.orbit_substitution()
            sage: print(s1)
            a->a, b->ab
            sage: (g0 * g1).orbit_substitution() == s0 * s1
            True
            sage: (g1 * g0).orbit_substitution() == s1 * s0
            True
        """
        words = [[i] for i in range(len(self._start))]
        for p, winner, side in self:
            top_letter = p._labels[0][side]
            bottom_letter = p._labels[1][side]
            loser_letter = p._labels[1-winner][side]
            words[loser_letter] = words[bottom_letter] + words[top_letter]
        perm_on_list_inplace(self._relabelling, words)

        return words if as_plain_list else self._int_list_to_substitutions(words)

    def interval_substitution(self, as_plain_list=False):
        r"""
        Return the substitution of intervals.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b','b a')
            sage: p0 = iet.FlipSequence(p, ['t'])
            sage: s0 = p0.interval_substitution()
            sage: print(s0)
            a->a, b->ba
            sage: p1 = iet.FlipSequence(p, ['b'])
            sage: s1 = p1.interval_substitution()
            sage: print(s1)
            a->ab, b->b
            sage: (p0 * p1).interval_substitution() == s1 * s0
            True
            sage: (p1 * p0).interval_substitution() == s0 * s1
            True

            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: f0 = iet.FlipSequence(p, ['t','b','t'], relabelling="(0,1,2)")
            sage: f1 = iet.FlipSequence(f0._end, ['b','t','b'], relabelling="(1,3,2)")
            sage: f2 = iet.FlipSequence(f1._end, ['t','b','t'], relabelling="(0,3,1,2)")
            sage: f = f0 * f1 * f2
            sage: f.interval_substitution()
            WordMorphism: a->ba, b->cdbadbadbc, c->dbadbcdbadb, d->adbcba
            sage: f.interval_substitution(as_plain_list=True)
            [[1, 0],
             [2, 3, 1, 0, 3, 1, 0, 3, 1, 2],
             [3, 1, 0, 3, 1, 2, 3, 1, 0, 3, 1],
             [0, 3, 1, 2, 1, 0]]
            sage: s0 = f0.interval_substitution()
            sage: s1 = f1.interval_substitution()
            sage: s2 = f2.interval_substitution()
            sage: f.interval_substitution() == s2 * s1 * s0
            True
        """
        r = self._relabelling
        words = [[r[i]] for i in range(len(self._start))]

        for p, winner, side in reversed(self):
            winner_letter = p._labels[winner][side]
            loser_letter = p._labels[1-winner][side]
            top_letter = p._labels[0][side]
            bottom_letter = p._labels[1][side]
            if side == 0:
                words[winner_letter] = words[loser_letter] + words[winner_letter]
            else:
                words[winner_letter] = words[winner_letter] + words[loser_letter]
        return words if as_plain_list else self._int_list_to_substitutions(words)

    substitution = orbit_substitution  # standard name
    dual_substitution = interval_substitution  # standard name

    def self_similar_interval_exchange_transformation(self):
        r"""
        Return the dilatation and the self-similar interval exchange
        transformation associated to this path.

        EXAMPLES::

            sage: from surface_dynamics import iet

        The golden rotation::

            sage: p = iet.Permutation('a b', 'b a')
            sage: g = iet.FlipSequence(p, ['t', 'b'])
            sage: a, T = g.self_similar_iet()
            sage: T.lengths().parent()
            Vector space of dimension 2 over Number Field ...
            sage: T.lengths().n()
            (1.00000000000000, 1.61803398874989)

        An example from Do-Schmidt::

            sage: code = [1,0,1,0,1,0,0,0,1,0,0,1,1,1,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,0]
            sage: p = iet.Permutation([0,1,2,3,4,5,6],[6,5,4,3,2,1,0])
            sage: g = iet.FlipSequence(p, code)
            sage: a, T = g.self_similar_iet()
            sage: T.sah_arnoux_fathi_invariant()
            (0, 0, 0)
        """
        if not self.is_complete():
            raise ValueError('not a complete loop')
        from sage.rings.qqbar import AA
        m = self.matrix()
        poly = m.charpoly()
        rho = max(poly.roots(AA, False))
        try:
            K, a, _ = rho.as_number_field_element(minimal=True, embedded=True)
        except TypeError:
            # NOTE: the embedded option appeared in sage 8.8
            from sage.rings.number_field.number_field import NumberField
            L, b, _ = rho.as_number_field_element(minimal=True)
            K = NumberField(L.polynomial(), str(L.gen()), embedding=rho)
            a = K(b.polynomial())

        lengths = (m - a).right_kernel().basis()[0]
        if lengths[0] < 0:
            lengths = -lengths
        if any(x <= 0 for x in lengths):
            raise RuntimeError("wrong Perron-Frobenius eigenvector: {}".format(lengths))

        # NOTE: the above code makes "lengths" with parent being the right kernel (that is
        # a submodule of R^d)
        lengths = lengths.parent().ambient_vector_space()(lengths)

        from .iet import IntervalExchangeTransformation
        return a, IntervalExchangeTransformation(self._start, lengths)

    self_similar_iet = self_similar_interval_exchange_transformation

