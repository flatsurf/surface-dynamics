r"""
Flip sequences of permutations of interval exchanges and linear involutions

A concatenation of:
- left rauzy moves (top or bottom)
- right rauzy moves (top or bottom)
- left_right_inverse, top_bottom_inverse, symmetric
- relabeled
"""
# ***************************************************************************
#       Copyright (C) 2021-2022 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************
import numbers
from copy import copy

from sage.structure.sage_object import SageObject
from sage.rings.all import ZZ
from sage.modules.free_module import FreeModule
from sage.matrix.matrix_space import MatrixSpace
import sage.matrix.matrix0
from sage.combinat.words.words import FiniteWords
from sage.combinat.words.morphism import WordMorphism

from .template import (interval_conversion, side_conversion,
        OrientablePermutationIET, PermutationLI)
from .labelled import LabelledPermutation
from .constructors import GeneralizedPermutation
from surface_dynamics.misc.permutation import (perm_init, perm_check,
        perm_compose, perm_one, perm_cycle_string, perm_orbit,
        perm_on_list_inplace, perm_invert, perm_is_one)


class RauzyFlipSequence(SageObject):
    r"""
    A flip sequence of interval exchange transformations or of linear
    involutions.

    Each step is a Rauzy move either on the left or on the right. The
    permutations are labelled, either permutations of interval exchanges
    (``PermutationIET``) or generalized permutations of linear involutions
    (``PermutationLI``), flipped or not.

    EXAMPLES::

        sage: from surface_dynamics import iet

    Making substitutions from a flip sequence::

        sage: p = iet.Permutation([['a','b','c','d'],['d','c','b','a']])
        sage: g = iet.FlipSequence(p, 'ttbtbbtb')
        sage: assert g.is_complete()
        sage: s1 = g.orbit_substitution()
        sage: print(s1)
        a->adbd, b->adbdbd, c->adccd, d->adcd
        sage: s2 = g.interval_substitution()
        sage: print(s2)
        a->abcd, b->bab, c->cdc, d->dcbababcd
        sage: s1.incidence_matrix() == s2.incidence_matrix().transpose()
        True

    If large powers appear in the path it might be more convenient to use
    power notations as in the following example::

        sage: p = iet.Permutation([['a','b','c','d'],['d','c','b','a']])
        sage: fs = iet.FlipSequence(p, 't^5b^3tbt^3b')
        sage: fs == iet.FlipSequence(p, 'tttttbbbtbtttb')
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

    The path of non-orientable flipped iet from [BasLop]_ that gives non-uniquely
    ergodic examples. Since :class:`RauzyFlipSequence` uses labelled permutations, the
    path needs a relabelling to be closed::

        sage: p0 = iet.Permutation([1,2,3,4,5,6,7,8,9,10], [9,10,1,2,3,4,5,6,7,8], flips=[1,2,3,4,5,6,7,10])
        sage: fs1 = iet.FlipSequence(p0, 't^7btbbtbtbtb')
        sage: p16 = fs1.end()
        sage: fs2 = iet.FlipSequence(p16, 'b')
        sage: assert fs2.is_closed()
        sage: fs3 = iet.FlipSequence(p16, 'tb^2t^3b^4t^5b^5')
        sage: p37 = fs3.end()
        sage: fs4 = iet.FlipSequence(p37, 't')
        sage: assert fs4.is_closed()
        sage: fs5 = iet.FlipSequence(p37, 'bt^7btb')
        sage: p49 = fs5.end()
        sage: assert p49.reduced() == p0.reduced()
        sage: fs6 = iet.FlipSequence(p49, [], relabelling=[7,1,2,3,4,5,6,9,8,0])
        sage: assert fs6.end() == p0
        sage: (fs1 * fs2 * fs3 * fs4 * fs5 * fs6).matrix()
        [ 2  2  2  2  2  2  2  2  3  2]
        [ 2  4  2  2  2  2  1  0  0  0]
        [ 0  0  2  2  2  1  0  0  0  0]
        [ 0  0  0  3  2  0  0  0  0  0]
        [ 0  0  1  2  2  0  0  0  0  0]
        [ 1  2  2  2  2  2  0  0  0  0]
        [ 2  3  2  2  2  2  2  0  0  0]
        [ 0  0  0  0  0  0  0  1  1  1]
        [ 1  1  1  1  1  1  1  1  2  1]
        [ 6 10 10 14 13  8  4  2  3  3]

    A symbolic matrix (with polynomial coefficients) where powers of the paths
    ``fs2`` and ``fs4`` are variables can be obtained with the function
    ``symbolic_matrix_power``::

        sage: from surface_dynamics.misc.linalg import symbolic_matrix_power
        sage: QQrs = QQ['r,s']
        sage: m1 = fs1.matrix()
        sage: m2 = symbolic_matrix_power(fs2.matrix(), QQrs.gen(0))
        sage: m3 = fs3.matrix()
        sage: m4 = symbolic_matrix_power(fs4.matrix(), QQrs.gen(1))
        sage: m5 = fs5.matrix()
        sage: m6 = fs6.matrix()
        sage: print(m1 * m2 * m3 * m4 * m5 * m6)
        [      2       2       2       2       2       2       2       2       3       2]
        [    2*s 2*s + 2       2       2       2       2       1       0       0       0]
        [      0       0       2       2       2       1       0       0       0       0]
        [      0       0       0   r + 2   r + 1       0       0       0       0       0]
        [      0       0       1       2       2       0       0       0       0       0]
        [      s   s + 1       2       2       2       2       0       0       0       0]
        [  s + 1   s + 2       2       2       2       2       2       0       0       0]
        [      0       0       0       0       0       0       0       1       1       1]
        [      1       1       1       1       1       1       1       1       2       1]
        [4*s + 2 4*s + 6      10  r + 13  r + 12       8       4       2       3       3]

    Flip sequences of generalized permutations (linear involutions) work the
    same way. A complete closed walk in the Rauzy diagram of the hyperelliptic
    component of `Q_2(2^2)`::

        sage: from surface_dynamics import Stratum
        sage: p = Stratum([2, 2], k=2).hyperelliptic_component().permutation_representative(reduced=False)
        sage: p
        0 1 2 1 3
        4 3 4 2 0
        sage: f = iet.FlipSequence(p, 'ttbtbtbb')
        sage: f.close()
        sage: f
        RauzyFlipSequence('0 1 2 1 3\n4 3 4 2 0', 'ttbtbtbb', relabelling=(3,4))
        sage: f.is_complete()
        True
        sage: f.matrix()
        [1 2 1 0 0]
        [0 1 0 1 0]
        [1 1 2 0 0]
        [2 3 3 0 1]
        [0 0 0 1 1]
        sage: f.matrix().charpoly()
        x^5 - 5*x^4 + 4*x^3 + 4*x^2 - 5*x + 1
    """
    def __init__(self, p, rauzy_moves=None, top_bottom_inverse=False, left_right_inverse=False, relabelling=None, check=True):
        r"""
        INPUT:

        - ``p`` - permutation or generalized permutation

        - ``rauzy_moves`` - (optional; default empty) a sequence of Rauzy moves. It can either
          be a list of valid Rauzy moves or a string with letters ``'t'``, ``'b'``, ``'T'`` or
          ``'B'`` (meaning respectively top-right, bottom-right, top-left and bottom-left Rauzy
          inductions)

        - ``top_bottom_inverse`` - (optional; default ``False``) whether a
          top-bottom inverse is performed at the end of the flip sequence

        - ``left_right_inverse`` - (optional; default ``False``) whether a
          left-right inverse is performed at the end of the flip sequence

        - ``relabelling`` - (optional; default identity permutation) the
          relabelling to be performed at the end of the flip sequence

        - ``check`` - (optional; default ``True``) when ``rauzy_moves`` is a
          list, whether its moves are converted and checked to be defined;
          with ``False``, each move is either a winner ``0`` (top) or ``1``
          (bottom) for a right move, or a pair ``(winner, side)`` with
          ``side`` either ``-1`` (right) or ``0`` (left), and the moves are
          trusted (see :meth:`rauzy_move`); a string ``rauzy_moves`` is
          always checked

        TESTS::

            sage: from surface_dynamics import iet
            sage: p = iet.Permutation('a b', 'b a')
            sage: iet.FlipSequence(p, 'tbcf')
            Traceback (most recent call last):
            ...
            ValueError: invalid flip sequence
            sage: iet.FlipSequence(p, 't^t')
            Traceback (most recent call last):
            ...
            ValueError: invalid flip sequence

        A Rauzy move of a generalized permutation may be undefined
        (Example 2.9 of [BoiLan09]_), and is then refused::

            sage: q = iet.GeneralizedPermutation('A B A', 'B D C C D')
            sage: q.has_rauzy_move('t'), q.has_rauzy_move('b')
            (True, False)
            sage: iet.FlipSequence(q, 't')
            RauzyFlipSequence('A B A\nB D C C D', 't')
            sage: iet.FlipSequence(q, 'b')
            Traceback (most recent call last):
            ...
            ValueError: Rauzy induction is not well defined
            sage: iet.FlipSequence(q, ['b'])
            Traceback (most recent call last):
            ...
            ValueError: Rauzy induction is not well defined

        With ``check=False`` the moves are trusted::

            sage: iet.FlipSequence(q, [1], check=False)
            RauzyFlipSequence('A B A\nB D C C D', 'b')

        Moves on both sides and the final inverses on a generalized
        permutation::

            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: f = iet.FlipSequence(p, [('t', 0), ('b', 0)])
            sage: f
            RauzyFlipSequence('0 1 2 1 3\n4 3 4 2 0', 'TB')
            sage: f.end()
            1 2 1 0 3
            3 4 2 4 0
            sage: f._check()
            sage: f = iet.FlipSequence(p, 'ttttb', top_bottom_inverse=True)
            sage: f.end()
            4 3 0 4 2
            0 1 2 3 1
            sage: f._check()
            sage: f = iet.FlipSequence(p, 'ttttb', left_right_inverse=True)
            sage: f.end()
            1 3 2 1 0
            2 4 0 3 4
            sage: f._check()
            sage: iet.FlipSequence(p, 'tbtt')
            Traceback (most recent call last):
            ...
            ValueError: Rauzy induction is not well defined
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
            if isinstance(rauzy_moves, str):
                    i = 0
                    while i < len(rauzy_moves):
                        r = rauzy_moves[i]
                        if r in ['T', 'B']:
                            s = 0
                            r = r.lower()
                        else:
                            if r not in ['t', 'b']:
                                raise ValueError('invalid flip sequence')
                            s = -1
                        if i + 1 < len(rauzy_moves) and rauzy_moves[i+1] == '^':
                            k = i + 2
                            if len(rauzy_moves) <= k or rauzy_moves[k] == '0' or not rauzy_moves[k].isdigit():
                                raise ValueError('invalid flip sequence')
                            k += 1
                            while k < len(rauzy_moves) and rauzy_moves[k].isdigit():
                                k += 1
                            power = int(rauzy_moves[i+2:k])
                            self.rauzy_move(r, s, iterations=power)
                            i = k
                        else:
                            self.rauzy_move(r, s)
                            i += 1
            else:
                for r in rauzy_moves:
                    if isinstance(r, numbers.Integral):
                        self.rauzy_move(r, check=check)
                    else:
                        self.rauzy_move(*r, check=check)

        if top_bottom_inverse:
            self.top_bottom_inverse()
        if left_right_inverse:
            self.left_right_inverse()
        if relabelling is not None:
            self.relabel(relabelling)

    def start(self):
        r"""
        Return the start of the path.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: iet.FlipSequence(p, 'ttb').start()
            0 1 2 1 3
            4 3 4 2 0
        """
        return self._start.__copy__()

    def end(self):
        r"""
        Return the end of the path.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: iet.FlipSequence(p, 'tttb').end()
            0 3 1 2 1
            4 3 4 2 0
        """
        return self._end.__copy__()

    def __hash__(self):
        r"""
        TESTS::

            sage: from surface_dynamics import iet
            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: f = iet.FlipSequence(p, 'ttt')
            sage: hash(f) == hash(iet.FlipSequence(p, 't^3'))
            True
            sage: len({f, iet.FlipSequence(p, 'ttt'), iet.FlipSequence(p, 'ttb')})
            2
        """
        return hash((self._start, tuple(self._rauzy_moves), self._top_bottom_inverse, self._left_right_inverse, tuple(self._relabelling)))

    def __eq__(self, other) -> bool:
        r"""
        TESTS::

            sage: from surface_dynamics import iet
            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: iet.FlipSequence(p, 'ttt') == iet.FlipSequence(p, 't^3')
            True
            sage: iet.FlipSequence(p, 'ttt') == iet.FlipSequence(p, 'ttb')
            False
            sage: iet.FlipSequence(p, 'ttt') != iet.FlipSequence(p, 'ttb')
            True
        """
        if type(self) != type(other):
            raise TypeError('incomparable')

        return self._start == other._start and \
               self._rauzy_moves == other._rauzy_moves and \
               self._top_bottom_inverse == other._top_bottom_inverse and \
               self._left_right_inverse == other._left_right_inverse and \
               self._relabelling == other._relabelling

    def __ne__(self, other) -> bool:
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

        On a generalized permutation::

            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: f = iet.FlipSequence(p, 'ttttb')
            sage: f.relabel([4, 3, 2, 1, 0])
            sage: f.top_bottom_inverse()
            sage: f.left_right_inverse()
            sage: f._check()
            sage: f._end = p
            sage: f._check()
            Traceback (most recent call last):
            ...
            ValueError: inconsistent flip sequence: end=0 1 2 1 3/4 3 4 2 0 while p=2 0 4 1 0/3 1 2 3 4
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
        r"""
        TESTS::

            sage: from surface_dynamics import iet
            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: f = iet.FlipSequence(p, 'ttt')
            sage: g = copy(f)
            sage: g == f
            True
            sage: g.rauzy_move('b')
            sage: g == f
            False
            sage: f
            RauzyFlipSequence('0 1 2 1 3\n4 3 4 2 0', 'ttt')
        """
        res = type(self).__new__(type(self))
        res._start = copy(self._start)
        res._end = copy(self._end)
        res._rauzy_moves = self._rauzy_moves[:]
        res._top_bottom_inverse = self._top_bottom_inverse
        res._left_right_inverse = self._left_right_inverse
        res._relabelling = self._relabelling
        return res

    @classmethod
    def _from_data(cls, start, end, rauzy_moves, relabelling):
        r"""
        Return the flip sequence from ``start`` to ``end`` with the given Rauzy
        moves and final relabelling, without replaying the moves.

        Nothing is checked and nothing is copied: the caller passes objects it
        owns, consistent with each other (see :meth:`_check`).

        INPUT:

        - ``start``, ``end`` -- labelled permutations, stored as given

        - ``rauzy_moves`` -- a list of pairs ``(winner, side)`` with ``winner``
          either ``0`` or ``1`` and ``side`` either ``-1`` or ``0``, stored as
          given

        - ``relabelling`` -- a list, the relabelling performed at the end of
          the flip sequence, stored as given

        OUTPUT: an ``RauzyFlipSequence`` with neither top-bottom nor left-right
        inverse

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: from surface_dynamics.interval_exchanges.flip_sequence import RauzyFlipSequence
            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: g = iet.FlipSequence(p, 'ttbtbbtb')
            sage: g.close()
            sage: f = RauzyFlipSequence._from_data(p.__copy__(), p.__copy__(),
            ....:                                [(0, -1), (0, -1), (1, -1), (0, -1), (1, -1), (1, -1), (0, -1), (1, -1)],
            ....:                                g._relabelling[:])
            sage: f == g
            True
            sage: f._check()

        A closed walk with a nontrivial relabelling::

            sage: q = iet.Permutation([0, 1, 2, 3, 4, 5, 6], [2, 0, 3, 5, 4, 6, 1])
            sage: g = iet.FlipSequence(q, 'bbtbtbb')
            sage: g.close()
            sage: g
            RauzyFlipSequence('0 1 2 3 4 5 6\n2 0 3 5 4 6 1', 'bbtbtbb', relabelling=(4,6,5))
            sage: moves = [(1, -1), (1, -1), (0, -1), (1, -1), (0, -1), (1, -1), (1, -1)]
            sage: f = RauzyFlipSequence._from_data(q.__copy__(), q.__copy__(), moves, [0, 1, 2, 3, 6, 4, 5])
            sage: f == g and f.is_closed()
            True
            sage: f._check()

        On a generalized permutation::

            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: moves = [(0, -1), (0, -1), (1, -1), (0, -1), (1, -1), (0, -1), (1, -1), (1, -1)]
            sage: f = RauzyFlipSequence._from_data(p.__copy__(), p.__copy__(), moves, [0, 1, 2, 4, 3])
            sage: f
            RauzyFlipSequence('0 1 2 1 3\n4 3 4 2 0', 'ttbtbtbb', relabelling=(3,4))
            sage: f._check()
        """
        res = cls.__new__(cls)
        res._start = start
        res._end = end
        res._rauzy_moves = rauzy_moves
        res._top_bottom_inverse = False
        res._left_right_inverse = False
        res._relabelling = relabelling
        return res

    def __mul__(self, other):
        r"""
        Return the concatenation of this flip sequence and ``other``.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: f = iet.FlipSequence(p, 'tt')
            sage: g = iet.FlipSequence(f.end(), 'b')
            sage: f * g
            RauzyFlipSequence('0 1 2 1 3\n4 3 4 2 0', 'ttb')
            sage: (f * g).matrix() == f.matrix() * g.matrix()
            True
            sage: g * f
            Traceback (most recent call last):
            ...
            ValueError: the left path should end with the start of the right path
        """
        if type(self) != type(other):
            raise TypeError('incompatible types for multiplication')
        if self._end != other._start:
            raise ValueError('the left path should end with the start of the right path')

        res = type(self).__new__(type(self))
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

    def __len__(self):
        r"""
        Return the number of Rauzy moves of this flip sequence.

        TESTS::

            sage: from surface_dynamics import iet
            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: len(iet.FlipSequence(p, 'ttbtbtbb'))
            8
        """
        return len(self._rauzy_moves)

    def __iter__(self):
        r"""
        Iterate over the triples ``(p, winner, side)`` of the Rauzy moves of
        this flip sequence, with ``p`` the permutation the move starts from.

        The permutation ``p`` is modified in place along the iteration.

        TESTS::

            sage: from surface_dynamics import iet
            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: for q, winner, side in iet.FlipSequence(p, 'tTb'):
            ....:     print(q.str('/'), winner, side)
            0 1 2 1 3/4 3 4 2 0 0 -1
            0 1 2 1 3/4 3 0 4 2 0 0
            0 1 2 1 3/3 4 0 4 2 1 -1
        """
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

        The same on a generalized permutation, and on flipped permutations::

            sage: q = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: for tb, lr in itertools.product([False, True], repeat=2):
            ....:     f = iet.FlipSequence(q, 'ttTtbB', tb, lr, "(0,4)(1,3)")
            ....:     l1 = [(copy(p), winner, side) for p, winner, side in f]
            ....:     l2 = [(copy(p), winner, side) for p, winner, side in reversed(f)]
            ....:     assert l2 == l1[::-1]
            sage: p = iet.Permutation('a b c', 'c a b', flips=['b'])
            sage: f = iet.FlipSequence(p, 'tbB')
            sage: l1 = [(copy(p), winner, side) for p, winner, side in f]
            sage: l1 == [(copy(p), winner, side) for p, winner, side in reversed(f)][::-1]
            True
            sage: p = iet.GeneralizedPermutation('a b c b', 'd c d a', flips='abcd')
            sage: f = iet.FlipSequence(p, 'tb')
            sage: l1 = [(copy(p), winner, side) for p, winner, side in f]
            sage: l1 == [(copy(p), winner, side) for p, winner, side in reversed(f)][::-1]
            True
        """
        p = copy(self._end)
        if self._top_bottom_inverse:
            p.top_bottom_inverse(inplace=True)
        if self._left_right_inverse:
            p.left_right_inverse(inplace=True)
        p.relabel(perm_invert(self._relabelling))
        # Only the backward move of OrientablePermutationIET works in place;
        # those of OrientablePermutationLI, FlippedPermutationIET and
        # FlippedPermutationLI take no ``inplace`` and return a new permutation.
        inplace = type(p).backward_rauzy_move is OrientablePermutationIET.backward_rauzy_move
        for winner, side in reversed(self._rauzy_moves):
            if inplace:
                p.backward_rauzy_move(winner, side, inplace=True, check=False)
            else:
                p = p.backward_rauzy_move(winner, side, check=False)
            yield (p, winner, side)

    def _simplified_flip_sequence(self):
        l = []
        for winner, side in self._rauzy_moves:
            winner = 'tb'[winner]
            if side == -1:
                l.append(winner)
            else:
                l.append(winner.upper())
        return ''.join(l)

    def _repr_(self) -> str:
        args = [repr(self._start.str()), repr(self._simplified_flip_sequence())]
        if self._top_bottom_inverse:
            args.append('top_bottom_inverse=True')
        if self._left_right_inverse:
            args.append('left_right_inverse=True')
        if not perm_is_one(self._relabelling):
            args.append('relabelling={}'.format(perm_cycle_string(self._relabelling, singletons=False)))
        return "RauzyFlipSequence({})".format(', '.join(arg for arg in args))

    def rauzy_move(self, winner, side=-1, iterations=1, check=True):
        r"""
        Add a Rauzy move to this flip sequence.

        INPUT:

        - ``winner`` -- either ``'t'`` (or ``0``) for top induction or ``'b'`` (or
          ``1``) for bottom induction
        - ``side`` -- either ``'l'`` (or ``0``) for left or ``'r'`` (or ``-1``) for right
        - ``iterations`` -- the number of times the induction is applied
        - ``check`` -- boolean (default ``True``) whether ``winner`` and ``side``
          are converted and the moves are checked to be defined (see
          ``has_rauzy_move`` of the permutations); with ``False``, ``winner``
          must already be ``0`` or ``1``, ``side`` must be ``-1`` or ``0``,
          the moves must be defined, and nothing is checked

        A move that is not defined raises a ``ValueError`` and leaves this
        flip sequence unchanged.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: f = iet.FlipSequence(p)
            sage: f.rauzy_move('t')
            sage: f.rauzy_move(1, -1, check=False)
            sage: f
            RauzyFlipSequence('a b c d\nd c b a', 'tb')
            sage: f == iet.FlipSequence(p, [0, 1], check=False)
            True

        On a generalized permutation, a Rauzy move may be undefined
        (Example 2.9 of [BoiLan09]_)::

            sage: q = iet.GeneralizedPermutation('A B A', 'B D C C D')
            sage: f = iet.FlipSequence(q)
            sage: f.rauzy_move('b')
            Traceback (most recent call last):
            ...
            ValueError: Rauzy induction is not well defined
            sage: f.rauzy_move('t', iterations=2)
            Traceback (most recent call last):
            ...
            ValueError: Rauzy induction is not well defined
            sage: f
            RauzyFlipSequence('A B A\nB D C C D', '')
            sage: f.end()
            A B A
            B D C C D
            sage: f.rauzy_move('t')
            sage: f
            RauzyFlipSequence('A B A\nB D C C D', 't')
            sage: f.end()
            D A B A
            B D C C

        With ``check=False``, the move is trusted::

            sage: f = iet.FlipSequence(q)
            sage: f.rauzy_move(1, -1, check=False)
            sage: f
            RauzyFlipSequence('A B A\nB D C C D', 'b')
        """
        if check:
            winner = interval_conversion(winner)
            side = side_conversion(side)
            # A single move is tested before it changes the end in place. With
            # several, a move refused after the first one must leave the end
            # unchanged, so they run on a copy.
            end = self._end if iterations == 1 else copy(self._end)
            for _ in range(iterations):
                if not end.has_rauzy_move(winner, side):
                    raise ValueError("Rauzy induction is not well defined")
                end.rauzy_move(winner, side, inplace=True, check=False)
            self._end = end
        else:
            for _ in range(iterations):
                self._end.rauzy_move(winner, side, inplace=True, check=False)

        if self._top_bottom_inverse:
            winner = 1 - winner
        if self._left_right_inverse:
            side = 0 if side == -1 else -1
        self._rauzy_moves.extend([(winner, side)]*iterations)

    def left_right_inverse(self):
        r"""
        Add a left-right inverse at the end of this flip sequence.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: f = iet.FlipSequence(p, 'ttttb')
            sage: f.left_right_inverse()
            sage: f
            RauzyFlipSequence('0 1 2 1 3\n4 3 4 2 0', 'ttttb', left_right_inverse=True)
            sage: f.end()
            1 3 2 1 0
            2 4 0 3 4
            sage: f._check()
        """
        self._left_right_inverse = True
        self._end.left_right_inverse(inplace=True)

    def top_bottom_inverse(self):
        r"""
        Add a top-bottom inverse at the end of this flip sequence.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: f = iet.FlipSequence(p, 'ttttb')
            sage: f.top_bottom_inverse()
            sage: f
            RauzyFlipSequence('0 1 2 1 3\n4 3 4 2 0', 'ttttb', top_bottom_inverse=True)
            sage: f.end()
            4 3 0 4 2
            0 1 2 3 1
            sage: f._check()
        """
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

        On a generalized permutation::

            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: iet.FlipSequence(p, 'ttt').find_closure()
            [0, 1, 2, 3, 4]
            sage: iet.FlipSequence(p, 'ttbtbtbb').find_closure()
            [0, 1, 2, 4, 3]
            sage: iet.FlipSequence(p, 'tt').find_closure() is None
            True
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
        r"""
        Add a relabelling to this flip sequence.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: f = iet.FlipSequence(p, 'ttttb')
            sage: f.relabel([4, 3, 2, 1, 0])
            sage: f
            RauzyFlipSequence('0 1 2 1 3\n4 3 4 2 0', 'ttttb', relabelling=(0,4)(1,3))
            sage: f.end()
            4 3 2 1 3
            0 1 4 0 2
            sage: f._check()
        """
        if not perm_check(p, len(self._start)):
            p = perm_init(p, n=len(self._start))
        self._end.relabel(p)
        self._relabelling = perm_compose(self._relabelling, p)

    def close(self):
        r"""
        Close this path with the unique relabelling that identifies its end to its start.

        EXAMPLES::

            sage: from surface_dynamics import iet

            sage: p = iet.Permutation('a b c d e f', 'f e d c b a')
            sage: q = p.rauzy_move('t').rauzy_move('b').rauzy_move('t')
            sage: f = iet.FlipSequence(q, ['t', 't', 'b', 'b', 'b', 't', 'b', 't'],
            ....:                      top_bottom_inverse=True,
            ....:                      left_right_inverse=True)
            sage: f.close()
            sage: f
            RauzyFlipSequence('a b f c d e\nf a e b d c', 'ttbbbtbt', top_bottom_inverse=True, left_right_inverse=True, relabelling=(0,2,1,3)(4,5))

        On a generalized permutation::

            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: f = iet.FlipSequence(p, 'ttbtbtbb')
            sage: f.is_closed()
            False
            sage: f.close()
            sage: f
            RauzyFlipSequence('0 1 2 1 3\n4 3 4 2 0', 'ttbtbtbb', relabelling=(3,4))
            sage: f.is_closed()
            True
            sage: iet.FlipSequence(p, 'tt').close()
            Traceback (most recent call last):
            ...
            ValueError: can not close this loop with a relabelling
        """
        if self._start == self._end:
            return
        p = self.find_closure()
        if p is None:
            raise ValueError('can not close this loop with a relabelling')
        self._end.relabel(p)
        self._relabelling = perm_compose(self._relabelling, p)

    def is_closed(self) -> bool:
        r"""
        Return whether the path is closed, that is whether its start and end coincide.

        EXAMPLES::

            sage: from surface_dynamics import iet

            sage: p = iet.Permutation('a b c d e f', 'f e d c b a')
            sage: iet.FlipSequence(p, 't^5').is_closed()
            True
            sage: iet.FlipSequence(p, 't^4').is_closed()
            False

        On a generalized permutation::

            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: iet.FlipSequence(p, 'ttt').is_closed()
            True
            sage: iet.FlipSequence(p, 'tt').is_closed()
            False

        TESTS::

            sage: iet.FlipSequence(p, 'ttt').is_loop()
            True
        """
        return self._start == self._end

    is_loop = is_closed

    def winners_losers(self):
        r"""
        Return the pair of list of winner letters and list of loser letters along the path.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: iet.FlipSequence(p, 'ttbtbtbb').winners_losers()
            ([3, 3, 4, 1, 0, 2, 0, 0], [0, 2, 3, 4, 1, 0, 2, 1])
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

        On a generalized permutation::

            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: iet.FlipSequence(p, 'ttt').winners()
            [3, 3, 3]
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

        On a generalized permutation::

            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: iet.FlipSequence(p, 'ttt').losers()
            [0, 2, 4]
        """
        return self.winners_losers()[1]

    def is_complete(self) -> bool:
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

        On a generalized permutation::

            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: iet.FlipSequence(p, 'ttt').is_complete()
            False
            sage: f = iet.FlipSequence(p, 'ttbtbtbb')
            sage: f.close()
            sage: f.is_complete()
            True

        TESTS::

            sage: f.is_full()
            True
            sage: iet.FlipSequence(p, 'tt').is_complete()
            Traceback (most recent call last):
            ...
            ValueError: not a loop
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
        r"""
        Return the Rauzy matrix of this path.

        TESTS::

            sage: from surface_dynamics import iet
            sage: p = iet.Permutation('a b', 'b a')
            sage: iet.FlipSequence(p, '').matrix()
            [1 0]
            [0 1]

        On a generalized permutation::

            sage: p = iet.GeneralizedPermutation([0, 1, 2, 1, 3], [4, 3, 4, 2, 0])
            sage: iet.FlipSequence(p, 'ttb').matrix()
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [1 0 1 1 0]
            [0 0 0 1 1]
            sage: iet.FlipSequence(p, 'ttt').matrix().charpoly()
            x^5 - 5*x^4 + 10*x^3 - 10*x^2 + 5*x - 1
        """
        V = FreeModule(ZZ, len(self._start))
        columns = [copy(v) for v in V.basis()]
        for p, winner, side in self:
            winner_letter = p._labels[winner][side]
            loser_letter = p._labels[1-winner][side]
            columns[loser_letter] += columns[winner_letter]
        m = MatrixSpace(ZZ, len(self._start))(columns).transpose()
        perm_on_list_inplace(self._relabelling, m, swap=sage.matrix.matrix0.Matrix.swap_columns)
        return m

    def _assert_abelian(self, method, reason):
        r"""
        Raise a ``NotImplementedError`` naming ``method`` and ``reason`` if
        this flip sequence is made of generalized permutations (linear
        involutions).

        TESTS::

            sage: from surface_dynamics import iet
            sage: p = iet.GeneralizedPermutation('0 1 2 1 3', '4 3 4 2 0')
            sage: f = iet.FlipSequence(p, 't')
            sage: f._assert_abelian('foo', 'no reason')
            Traceback (most recent call last):
            ...
            NotImplementedError: foo is only implemented for interval exchanges, not for linear involutions: no reason
            sage: iet.FlipSequence(iet.Permutation('a b', 'b a'), 't')._assert_abelian('foo', 'no reason')
        """
        if isinstance(self._start, PermutationLI):
            raise NotImplementedError('{} is only implemented for interval exchanges, not for linear involutions: {}'.format(method, reason))

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

        TESTS:

        Linear involutions are not supported::

            sage: p = iet.GeneralizedPermutation('0 1 2 1 3', '4 3 4 2 0')
            sage: iet.FlipSequence(p, 'ttt').orbit_substitution()
            Traceback (most recent call last):
            ...
            NotImplementedError: orbit_substitution is only implemented for interval exchanges, not for linear involutions: the letters of the substitution would need signs
            sage: iet.FlipSequence(p, 'ttt').substitution()
            Traceback (most recent call last):
            ...
            NotImplementedError: orbit_substitution is only implemented for interval exchanges, not for linear involutions: the letters of the substitution would need signs
        """
        self._assert_abelian('orbit_substitution', 'the letters of the substitution would need signs')
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

        TESTS:

        Linear involutions are not supported::

            sage: p = iet.GeneralizedPermutation('0 1 2 1 3', '4 3 4 2 0')
            sage: iet.FlipSequence(p, 'ttt').interval_substitution()
            Traceback (most recent call last):
            ...
            NotImplementedError: interval_substitution is only implemented for interval exchanges, not for linear involutions: the letters of the substitution would need signs
            sage: iet.FlipSequence(p, 'ttt').dual_substitution()
            Traceback (most recent call last):
            ...
            NotImplementedError: interval_substitution is only implemented for interval exchanges, not for linear involutions: the letters of the substitution would need signs
        """
        self._assert_abelian('interval_substitution', 'the letters of the substitution would need signs')
        r = self._relabelling
        words = [[r[i]] for i in range(len(self._start))]

        for p, winner, side in reversed(self):
            winner_letter = p._labels[winner][side]
            loser_letter = p._labels[1-winner][side]
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

        TESTS:

        Linear involutions are not supported, even along a complete closed walk::

            sage: p = iet.GeneralizedPermutation('0 1 2 1 3', '4 3 4 2 0')
            sage: f = iet.FlipSequence(p, 'ttbtbtbb')
            sage: f.close()
            sage: f.is_complete()
            True
            sage: f.self_similar_interval_exchange_transformation()
            Traceback (most recent call last):
            ...
            NotImplementedError: self_similar_interval_exchange_transformation is only implemented for interval exchanges, not for linear involutions: there is no linear involution analogue of IntervalExchangeTransformation
            sage: f.self_similar_iet()
            Traceback (most recent call last):
            ...
            NotImplementedError: self_similar_interval_exchange_transformation is only implemented for interval exchanges, not for linear involutions: there is no linear involution analogue of IntervalExchangeTransformation
        """
        self._assert_abelian('self_similar_interval_exchange_transformation',
                             'there is no linear involution analogue of IntervalExchangeTransformation')
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


class IETFlipSequence(RauzyFlipSequence):
    r"""
    Deprecated name of :class:`RauzyFlipSequence`.

    The class stays so that the objects pickled under this name still load
    (unpickling does not call the constructor and does not warn).

    EXAMPLES::

        sage: from surface_dynamics import iet
        sage: from surface_dynamics.interval_exchanges.flip_sequence import IETFlipSequence
        sage: p = iet.Permutation('a b c d', 'd c b a')
        sage: f = IETFlipSequence(p, 'ttbtbbtb')
        doctest:warning
        ...
        UserWarning: IETFlipSequence has been deprecated and will be removed in a future version of surface-dynamics; use RauzyFlipSequence
        sage: f
        RauzyFlipSequence('a b c d\nd c b a', 'ttbtbbtb')
        sage: f.is_complete()
        True
        sage: loads(dumps(f)) == f
        True
    """
    def __init__(self, *args, **kwds):
        import warnings
        warnings.warn('IETFlipSequence has been deprecated and will be removed in a future version of surface-dynamics; use RauzyFlipSequence')
        super().__init__(*args, **kwds)
