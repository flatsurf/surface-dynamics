r"""
Template for permutations of interval exchange transformations

This file define high level operations on permutations (alphabet,
the different rauzy moves, ...) shared by reduced and labeled
permutations.

AUTHORS:

- Vincent Delecroix (2008-12-20): initial version

- Vincent Delecroix (2010-02-11): datatype simplification

TODO:

- disallow access to stratum, stratum component for permutations with flip

- construct dynamic Rauzy graphs and paths

- construct coherent _repr_

"""
#*****************************************************************************
#       Copyright (C) 2008 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function, absolute_import
from six.moves import range, map, filter, zip
from six import iterkeys, iteritems

from functools import total_ordering

from sage.structure.sage_object import SageObject

from copy import copy

from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import identity_matrix, matrix
from sage.misc.nested_class import NestedClassMetaclass

from surface_dynamics.misc.permutation import perms_canonical_labels

def to_fat_graphs(twin):
    lt = len(twin[0])
    lb = len(twin[1])
    assert twin[1][-1] == (0,0)

    if any(i == 1 for i,j in twin[0][1:]):
        # non-separating
        n = lt + lb - 2

        labels = [[None]*lt, [None]*lb]
        labels[0][0] = labels[1][-1] = -1
        k = 0
        ep = [None] * n
        for i in range(2):
            for j in range(len(labels[i])):
                if labels[i][j] is None:
                    ii,jj = twin[i][j]
                    labels[i][j] = k
                    labels[ii][jj] = k+1
                    ep[k] = k+1
                    ep[k+1] = k
                    k += 2
        assert k == n
        labels[0].pop(0)
        labels[1].pop(-1)

        fp = [None] * n
        for i in range(lt-1):
            j = (i-1) % (lt-1)
            fp[labels[0][i]] = labels[0][j]
        for i in range(lb-1):
            j = (i+1) % (lb-1)
            fp[labels[1][i]] = labels[1][j]

        return [(ep, fp)]
    else:
        # separating
        labelstop = [None] * lt
        labelstop[0] = -1
        eptop = [None] * (lt-1)
        k = 0
        for j in range(len(labelstop)):
            if labelstop[j] is None:
                jj = twin[0][j][1]
                labelstop[j] = k
                labelstop[jj] = k+1
                eptop[k] = k+1
                eptop[k+1] = k
                k += 2
        labelstop.pop(0)
        fptop = [None] * (lt-1)
        for i in range(lt-1):
            j = (i-1) % (lt-1)
            fptop[labelstop[i]] = labelstop[j]

        labelsbot = [None] * lb
        labelsbot[-1] = -1
        epbot = [None] * (lb-1)
        k = 0
        for j in range(len(labelsbot)):
            if labelsbot[j] is None:
                jj = twin[1][j][1]
                labelsbot[j] = k
                labelsbot[jj] = k+1
                epbot[k] = k+1
                epbot[k+1] = k
                k += 2
        labelsbot.pop(-1)
        fpbot = [None] * (lb-1)
        for i in range(lb-1):
            j = (i+1) % (lb-1)
            fpbot[labelsbot[i]] = labelsbot[j]

        return [(eptop, fptop), (epbot, fpbot)]

def cylindric_canonical(p):
    r"""
    TESTS::

        sage: from surface_dynamics import Stratum
        sage: from surface_dynamics.interval_exchanges.template import cylindric_canonical

        sage: C = Stratum({1:1, -1:5}, k=2).unique_component()
        sage: R = C.permutation_representative().rauzy_diagram()
        sage: K = [p for p in R if p.is_cylindric()]
        sage: Kcan = set(cylindric_canonical(p) for p in K)
        sage: Kcan
        {((1, 0), (1, 2, 3, 4, 5, 0))}
        sage: C = Stratum({1: 2, -1: 6}, k=2).unique_component()
        sage: R = C.permutation_representative().rauzy_diagram()
        sage: K = [p for p in R if p.is_cylindric()]
        sage: Kcan = set(cylindric_canonical(p) for p in K)
        sage: Kcan
        {((1, 0), (1, 2, 3, 4, 6, 0, 7, 8, 9, 5)),
         ((1, 2, 3, 4, 5, 0), (1, 2, 3, 4, 5, 0))}
    """
    twin = p._twin
    if twin[1][-1] != (0,0):
        lt = len(twin[0])
        lb = len(twin[1])
        p._move(0, lt-1, 0, 0)
        p._move(1, 0, 1, lb)
    fg = [perms_canonical_labels([ep,fp])[0][1] for ep,fp in to_fat_graphs(twin)]
    fg.sort()
    return tuple(map(tuple, fg))

def interval_conversion(interval=None):
    r"""
    Converts the argument in 0 or 1.

    INPUT:

    - ``winner`` - 'top' (or 't' or 0) or bottom (or 'b' or 1)

    OUTPUT:

    integer -- 0 or 1

    TESTS::

        sage: from surface_dynamics import *

        sage: from surface_dynamics.interval_exchanges.template import interval_conversion
        sage: interval_conversion('top')
        0
        sage: interval_conversion('t')
        0
        sage: interval_conversion(0)
        0
        sage: interval_conversion('bottom')
        1
        sage: interval_conversion('b')
        1
        sage: interval_conversion(1)
        1

    .. Non admissible strings raise a ValueError::

        sage: interval_conversion('')
        Traceback (most recent call last):
        ...
        ValueError: the interval can not be the empty string
        sage: interval_conversion('right')
        Traceback (most recent call last):
        ...
        ValueError: 'right' can not be converted to interval
        sage: interval_conversion('top_right')
        Traceback (most recent call last):
        ...
        ValueError: 'top_right' can not be converted to interval
    """
    if isinstance(interval, (int,Integer)):
        if interval != 0 and interval != 1:
            raise ValueError("interval must be 0 or 1")
        else:
            return interval

    elif isinstance(interval,str):
        if not interval: raise ValueError("the interval can not be the empty string")
        elif 'top'.startswith(interval): return 0
        elif 'bottom'.startswith(interval): return 1
        else: raise ValueError("'%s' can not be converted to interval" %(interval))

    else:
        raise TypeError("'%s' is not an admissible type" %(str(interval)))

def side_conversion(side=None):
    r"""
    Converts the argument in 0 or -1.

    INPUT:

    - ``side`` - either 'left' (or 'l' or 0) or 'right' (or 'r' or -1)

    OUTPUT:

    integer -- 0 or -1

    TESTS::

        sage: from surface_dynamics.interval_exchanges.template import side_conversion
        sage: side_conversion('left')
        0
        sage: side_conversion('l')
        0
        sage: side_conversion(0)
        0
        sage: side_conversion('right')
        -1
        sage: side_conversion('r')
        -1
        sage: side_conversion(1)
        -1
        sage: side_conversion(-1)
        -1

    .. Non admissible strings raise a ValueError::

        sage: side_conversion('')
        Traceback (most recent call last):
        ...
        ValueError: no empty string for side
        sage: side_conversion('top')
        Traceback (most recent call last):
        ...
        ValueError: 'top' can not be converted to a side
    """
    if side is None:
        return -1

    elif isinstance(side,str):
        if not side: raise ValueError("no empty string for side")
        if 'left'.startswith(side): return 0
        elif 'right'.startswith(side): return -1
        raise ValueError("'%s' can not be converted to a side" %(side))

    elif isinstance(side, (int,Integer)):
        if side == 0: return 0
        elif side == 1 or side == -1: return -1
        else: raise ValueError("side must be 0 or 1")

    else:
        raise TypeError("'%s' is not an admissible type" %(str(side)))

#
# NICE PRINTING OF FLIPS
#

def labelize_flip(couple):
    r"""
    Returns a string from a 2-uple couple of the form (name, flip).

    TESTS::

        sage: from surface_dynamics.interval_exchanges.template import labelize_flip
        sage: labelize_flip((0,1))
        ' 0'
        sage: labelize_flip((0,-1))
        '-0'
    """
    if couple[1] == -1: return '-' + str(couple[0])
    return ' ' + str(couple[0])

#
# CLASSES FOR PERMUTATIONS
#

@total_ordering
class Permutation(SageObject):
    r"""
    Template for all permutations.

    .. warning::

        Internal class! Do not use directly!

    This class implement generic algorithm (stratum, connected component, ...)
    and unfies all its children.

    It has four attributes

    - ``_alphabet`` -- the alphabet on which the permutation is defined. Be
      careful, it might have a different cardinality as the size of the
      permutation!

     - ``_twin`` -- the permutation

     - ``_labels`` -- None or the list of labels

     - ``_flips`` -- None or the list of flips (each flip is either ``1`` or
       ``-1``)

    The datatype for ``_twin`` differs for IET and LI (TODO: unify).
    """
    _alphabet = None
    _twin = None
    _labels = None
    _flips = None

    def __init__(self, intervals=None, alphabet=None, reduced=False, flips=None):
        r"""
        INPUT:

        - ``intervals`` - the intervals as a list of two lists

        - ``alphabet`` - something that should be converted to an alphabet

        - ``reduced`` - (boolean) whether the permutation is reduced or labeled

        - ``flips`` - (optional) a list of letters of the alphabet to be flipped (in which
          case the permutation corresponds to non-orientable surface)

        TESTS::

            sage: from surface_dynamics.interval_exchanges.labelled import LabelledPermutationIET

            sage: p1 = LabelledPermutationIET([[1,2,3],[3,2,1]])
            sage: p1 == loads(dumps(p1))
            True
            sage: p2 = LabelledPermutationIET([['a', 'b', 'c'], ['c', 'b', 'a']])
            sage: p2 == loads(dumps(p2))
            True
            sage: p3 = LabelledPermutationIET([['1','2','3'],['3','2','1']])
            sage: p3 == loads(dumps(p3))
            True
            sage: from surface_dynamics.interval_exchanges.labelled import LabelledPermutationLI
            sage: p1 = LabelledPermutationLI([[1,2,2],[3,3,1]])
            sage: p1 == loads(dumps(p1))
            True
            sage: p2 = LabelledPermutationLI([['a','b','b'],['c','c','a']])
            sage: p2 == loads(dumps(p2))
            True
            sage: p3 = LabelledPermutationLI([['1','2','2'],['3','3','1']])
            sage: p3 == loads(dumps(p3))
            True

            sage: from surface_dynamics.interval_exchanges.reduced import ReducedPermutationIET
            sage: p = ReducedPermutationIET()
            sage: loads(dumps(p)) == p
            True
            sage: p = ReducedPermutationIET([['a','b'],['b','a']])
            sage: loads(dumps(p)) == p
            True
            sage: from surface_dynamics.interval_exchanges.reduced import ReducedPermutationLI
            sage: p = ReducedPermutationLI()
            sage: loads(dumps(p)) == p
            True
            sage: p = ReducedPermutationLI([['a','a'],['b','b']])
            sage: loads(dumps(p)) == p
            True

            sage: from surface_dynamics.interval_exchanges.labelled import FlippedLabelledPermutationIET
            sage: p = FlippedLabelledPermutationIET([['a','b'],['a','b']],flips='a')
            sage: p == loads(dumps(p))
            True
            sage: p = FlippedLabelledPermutationIET([['a','b'],['b','a']],flips='ab')
            sage: p == loads(dumps(p))
            True

            sage: from surface_dynamics.interval_exchanges.labelled import FlippedLabelledPermutationLI
            sage: p = FlippedLabelledPermutationLI([['a','a','b'],['b','c','c']],flips='a')
            sage: p == loads(dumps(p))
            True
            sage: p = FlippedLabelledPermutationLI([['a','a'],['b','b','c','c']],flips='ac')
            sage: p == loads(dumps(p))
            True

            sage: from surface_dynamics import iet
            sage: p = iet.Permutation('a b','b a',reduced=True,flips='a')
            sage: p == loads(dumps(p))
            True
            sage: p = iet.Permutation('a b','b a',reduced=True,flips='b')
            sage: p == loads(dumps(p))
            True
            sage: p = iet.Permutation('a b','b a',reduced=True,flips='ab')
            sage: p == loads(dumps(p))
            True
            sage: p = iet.GeneralizedPermutation('a a','b b',reduced=True,flips='a')
            sage: p == loads(dumps(p))
            True
            sage: p = iet.GeneralizedPermutation('a a','b b',reduced=True,flips='b')
            sage: p == loads(dumps(p))
            True
            sage: p = iet.GeneralizedPermutation('a a','b b',reduced=True,flips='ab')
            sage: p == loads(dumps(p))
            True
        """
        # this constructor assumes that several methods are present
        #  _init_twin(intervals)
        #  _set_alphabet(alphabet)
        #  _init_alphabet(intervals)
        #  _init_flips(intervals, flips)

        # setting twins
        if intervals is None:
            self._twin = [[], []]
        else:
            self._init_twin(intervals)

        # setting alphabet
        if alphabet is not None:
            self._set_alphabet(alphabet)
        elif intervals is not None:
            self._init_alphabet(intervals)

        # optionally setting labels
        if intervals is not None and not reduced:
            self._labels = [
                list(map(self._alphabet.rank, intervals[0])),
                list(map(self._alphabet.rank, intervals[1]))]

        # optionally setting flips
        if flips is not None:
            self._init_flips(intervals, flips)

    def _init_flips(self,intervals,flips):
        r"""
        Initialize the flip list

        TESTS::

            sage: from surface_dynamics import *

            sage: iet.Permutation('a b','b a',flips='a').flips() #indirect doctest
            ['a']
            sage: iet.Permutation('a b','b a',flips='b').flips() #indirect doctest
            ['b']
            sage: iet.Permutation('a b','b a',flips='ab').flips() #indirect doctest
            ['a', 'b']

        ::

            sage: iet.GeneralizedPermutation('a a','b b',flips='a').flips()
            ['a']
            sage: iet.GeneralizedPermutation('a a','b b',flips='b').flips()
            ['b']
            sage: iet.GeneralizedPermutation('a a','b b',flips='ab').flips()
            ['a', 'b']
        """
        self._flips = [[1]*self.length_top(), [1]*self.length_bottom()]
        for interval in (0,1):
            for i,letter in enumerate(intervals[interval]):
                if letter in flips:
                    self._flips[interval][i] = -1

    def __eq__(self,other):
        r"""
        Tests equality

        TESTS::

            sage: from surface_dynamics import *

            sage: p1 = iet.Permutation('a b', 'a b', reduced=True, alphabet='ab')
            sage: p2 = iet.Permutation('a b', 'a b', reduced=True, alphabet='ba')
            sage: q1 = iet.Permutation('a b', 'b a', reduced=True, alphabet='ab')
            sage: q2 = iet.Permutation('a b', 'b a', reduced=True, alphabet='ba')
            sage: p1 == p2 and p2 == p1 and q1 == q2 and q2 == q1
            True
            sage: p1 == q1 or p2 == q1 or q1 == p1 or q1 == p2
            False

            sage: p1 = iet.Permutation('a b', 'b a', alphabet='ab')
            sage: p2 = iet.Permutation('a b', 'b a', alphabet='ba')
            sage: q1 = iet.Permutation('b a', 'a b', alphabet='ab')
            sage: q2 = iet.Permutation('b a', 'a b', alphabet='ba')
            sage: p1 == p2 and p2 == p1
            True
            sage: p1 == q1 or q1 == p1
            False
            sage: p1 == q2 or q2 == p1
            False
            sage: p2 == q1 or q1 == p2
            False
            sage: p2 == q2 or q2 == p2
            False
            sage: q1 == q2 or q2 == q1
            True

        ::

            sage: p1 = iet.GeneralizedPermutation('a a', 'b b', alphabet='ab')
            sage: p2 = iet.GeneralizedPermutation('a a', 'b b', alphabet='ba')
            sage: q1 = iet.GeneralizedPermutation('b b', 'a a', alphabet='ab')
            sage: q2 = iet.GeneralizedPermutation('b b', 'a a', alphabet='ba')
            sage: p1 == p2 and p2 == p1
            True
            sage: p1 == q1 or q1 == p1
            False
            sage: p1 == q2 or q2 == p1
            False
            sage: p2 == q1 or q1 == p2
            False
            sage: p2 == q2 or q2 == p2
            False
            sage: q1 == q2 and q2 == q1
            True

        ::

            sage: p = iet.GeneralizedPermutation('a b b', 'c c a', reduced=True)
            sage: q = iet.GeneralizedPermutation('b a a', 'c c b', reduced=True)
            sage: r = iet.GeneralizedPermutation('t s s', 'w w t', reduced=True)
            sage: p == q
            True
            sage: p == r
            True

        ::

            sage: p = iet.Permutation('a b', 'a b', reduced=True, flips='a')
            sage: q = copy(p)
            sage: q.alphabet([0,1])
            sage: p == q
            True
            sage: l0 = ['a b','a b']
            sage: l1 = ['a b','b a']
            sage: l2 = ['b a', 'a b']
            sage: p0 = iet.Permutation(l0, reduced=True, flips='ab')
            sage: p1 = iet.Permutation(l1, reduced=True, flips='a')
            sage: p2 = iet.Permutation(l2, reduced=True, flips='b')
            sage: p3 = iet.Permutation(l1, reduced=True, flips='ab')
            sage: p4 = iet.Permutation(l2 ,reduced=True,flips='ab')
            sage: p0 == p1 or p0 == p2 or p0 == p3 or p0 == p4
            False
            sage: p1 == p2 and p3 == p4
            True
            sage: p1 == p3 or p1 == p4 or p2 == p3 or p2 == p4
            False

        ::

            sage: a0 = [0,0,1]
            sage: a1 = [1,2,2]
            sage: p = iet.GeneralizedPermutation(a0, a1, reduced=True, flips=[0])
            sage: q = copy(p)
            sage: q.alphabet("abc")
            sage: p == q
            True
            sage: b0 = [1,0,0]
            sage: b1 = [2,2,1]
            sage: r = iet.GeneralizedPermutation(b0, b1, reduced=True, flips=[0])
            sage: p == r or q == r
            False

        ::

            sage: p1 = iet.Permutation('a b c', 'c b a', flips='a')
            sage: p2 = iet.Permutation('a b c', 'c b a', flips='b')
            sage: p3 = iet.Permutation('d e f', 'f e d', flips='d')
            sage: p1 == p1 and p2 == p2 and p3 == p3
            True
            sage: p1 == p2
            False
            sage: p1 == p3
            False
            sage: p1.reduced() == p3.reduced()
            True

            sage: p1 = iet.Permutation('a b', 'a b', reduced=True, alphabet='ab')
            sage: p2 = iet.Permutation('a b', 'a b', reduced=True, alphabet='ba')
            sage: q1 = iet.Permutation('a b', 'b a', reduced=True, alphabet='ab')
            sage: q2 = iet.Permutation('a b', 'b a', reduced=True, alphabet='ba')
            sage: p1 != p2 or p2 != p1 or q1 != q2 or q2 != q1
            False
            sage: p1 != q1 and p2 != q1 and q1 != p1 and q1 != p2
            True

            sage: p1 = iet.Permutation('a b', 'b a', alphabet='ab')
            sage: p2 = iet.Permutation('a b', 'b a', alphabet='ba')
            sage: q1 = iet.Permutation('b a', 'a b', alphabet='ab')
            sage: q2 = iet.Permutation('b a', 'a b', alphabet='ba')
            sage: p1 != p2 or p2 != p1
            False
            sage: p1 != q1 and q1 != p1
            True
            sage: p1 != q2 and q2 != p1
            True
            sage: p2 != q1 and q1 != p2
            True
            sage: p2 != q2 and q2 != p2
            True
            sage: q1 != q2 or q2 != q1
            False

        ::

            sage: p1 = iet.GeneralizedPermutation('a a', 'b b', alphabet='ab')
            sage: p2 = iet.GeneralizedPermutation('a a', 'b b', alphabet='ba')
            sage: q1 = iet.GeneralizedPermutation('b b', 'a a', alphabet='ab')
            sage: q2 = iet.GeneralizedPermutation('b b', 'a a', alphabet='ba')
            sage: p1 != p2 or p2 != p1
            False
            sage: p1 != q1 and q1 != p1
            True
            sage: p1 != q2 and q2 != p1
            True
            sage: p2 != q1 and q1 != p2
            True
            sage: p2 != q2 and q2 != p2
            True
            sage: q1 != q2 or q2 != q1
            False

        ::

            sage: p = iet.GeneralizedPermutation('a b b', 'c c a', reduced=True)
            sage: q = iet.GeneralizedPermutation('b b a', 'c c a', reduced=True)
            sage: r = iet.GeneralizedPermutation('i j j', 'k k i', reduced=True)
            sage: p != q
            True
            sage: p != r
            False

        ::

            sage: p = iet.Permutation('a b','a b',reduced=True,flips='a')
            sage: q = copy(p)
            sage: q.alphabet([0,1])
            sage: p != q
            False
            sage: l0 = ['a b','a b']
            sage: l1 = ['a b','b a']
            sage: l2 = ['b a', 'a b']
            sage: p0 = iet.Permutation(l0, reduced=True, flips='ab')
            sage: p1 = iet.Permutation(l1, reduced=True, flips='a')
            sage: p2 = iet.Permutation(l2, reduced=True, flips='b')
            sage: p3 = iet.Permutation(l1, reduced=True, flips='ab')
            sage: p4 = iet.Permutation(l2 ,reduced=True,flips='ab')
            sage: p0 != p1 and p0 != p2 and p0 != p3 and p0 != p4
            True
            sage: p1 != p2 or p3 != p4
            False
            sage: p1 != p3 and p1 != p4 and p2 != p3 and p2 != p4
            True

        ::

            sage: a0 = [0,0,1]
            sage: a1 = [1,2,2]
            sage: p = iet.GeneralizedPermutation(a0,a1,reduced=True,flips=[0])
            sage: q = copy(p)
            sage: q.alphabet("abc")
            sage: p != q
            False
            sage: b0 = [1,0,0]
            sage: b1 = [2,2,1]
            sage: r = iet.GeneralizedPermutation(b0,b1,reduced=True,flips=[0])
            sage: p != r and q != r
            True

        ::

            sage: p1 = iet.Permutation('a b c','c b a',flips='a')
            sage: p2 = iet.Permutation('a b c','c b a',flips='b')
            sage: p3 = iet.Permutation('d e f','f e d',flips='d')
            sage: p1 != p1 or p2 != p2 or p3 != p3
            False
            sage: p1 != p2
            True
            sage: p1 != p3
            True
            sage: p1.reduced() != p3.reduced()
            False
        """
        if (type(self) != type(other) or
            self._twin != other._twin or
            self._flips != other._flips):
            return False

        if self._labels is not None and other._alphabet is not None:
            if self._alphabet is other._alphabet:
                return self._labels == other._labels
            else:
                # (slower) comparison over different alphabets using letters
                return self.list() == other.list()

        return True

    def __lt__(self, other):
        r"""
        Defines a natural lexicographic order.

        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('a b','a b',reduced=True)
            sage: q = copy(p)
            sage: q.alphabet([0,1])
            sage: p == q
            True
            sage: p0 = iet.GeneralizedPermutation('a b', 'a b', reduced=True)
            sage: p1 = iet.GeneralizedPermutation('a b', 'b a', reduced=True)
            sage: p0 < p1 and p1 > p0
            True
            sage: q0 = iet.GeneralizedPermutation('a b c', 'a b c', reduced=True)
            sage: q1 = iet.GeneralizedPermutation('a b c', 'a c b', reduced=True)
            sage: q2 = iet.GeneralizedPermutation('a b c', 'b a c', reduced=True)
            sage: q3 = iet.GeneralizedPermutation('a b c', 'c a b', reduced=True)
            sage: q4 = iet.GeneralizedPermutation('a b c', 'b c a', reduced=True)
            sage: q5 = iet.GeneralizedPermutation('a b c', 'c b a', reduced=True)
            sage: p0 < q0 and q0 > p0 and p1 < q0 and q0 > p1
            True
            sage: q0 < q1 and q1 > q0
            True
            sage: q1 < q2 and q2 > q1
            True
            sage: q2 < q3 and q3 > q2
            True
            sage: q3 < q4 and q4 > q3
            True
            sage: q4 < q5 and q5 > q4
            True

            sage: p0 = iet.Permutation('1 2', '1 2')
            sage: p1 = iet.Permutation('1 2', '2 1')
            sage: p0 != p0
            False
            sage: (p0 == p0) and (p0 < p1)
            True
            sage: (p1 > p0) and (p1 == p1)
            True

            sage: p0 = iet.GeneralizedPermutation('0 0', '1 1 2 2')
            sage: p1 = iet.GeneralizedPermutation('0 0', '1 2 1 2')
            sage: p2 = iet.GeneralizedPermutation('0 0', '1 2 2 1')
            sage: p3 = iet.GeneralizedPermutation('0 0 1 1', '2 2')
            sage: p4 = iet.GeneralizedPermutation('0 0 1', '1 2 2')
            sage: p5 = iet.GeneralizedPermutation('0 1 0 1', '2 2')
            sage: p6 = iet.GeneralizedPermutation('0 1 1 0', '2 2')
            sage: p0 == p0 and p0 < p1 and p0 < p2 and p0 < p3 and p0 < p4
            True
            sage: p0 < p5 and p0 < p6 and p1 < p2 and p1 < p3 and p1 < p4
            True
            sage: p1 < p5 and p1 < p6 and p2 < p3 and p2 < p4 and p2 < p5
            True
            sage: p2 < p6 and p3 < p4 and p3 < p5 and p3 < p6 and p4 < p5
            True
            sage: p4 < p6 and p5 < p6 and p0 == p0 and p1 == p1 and p2 == p2
            True
            sage: p3 == p3 and p4 == p4 and p5 == p5 and p6 == p6
            True

            sage: p = iet.Permutation('a b', 'a b', reduced=True, flips='a')
            sage: q = copy(p)
            sage: q.alphabet([0,1])
            sage: p == q
            True
            sage: l0 = ['a b','a b']
            sage: l1 = ['a b','b a']
            sage: p1 = iet.Permutation(l1, reduced=True, flips='b')
            sage: p2 = iet.Permutation(l1, reduced=True, flips='a')
            sage: p3 = iet.Permutation(l1, reduced=True, flips='ab')
            sage: p2 > p3 and p3 < p2
            True
            sage: p1 > p2 and p2 < p1
            True
            sage: p1 > p3 and p3 < p1
            True
            sage: q1 = iet.Permutation(l0, reduced=True, flips='a')
            sage: q2 = iet.Permutation(l0, reduced=True, flips='b')
            sage: q3 = iet.Permutation(l0, reduced=True, flips='ab')
            sage: q2 > q1 and q2 > q3 and q1 < q2 and q3 < q2
            True
            sage: q1 > q3
            True
            sage: q3 < q1
            True
            sage: r = iet.Permutation('a b c', 'a b c', reduced=True, flips='a')
            sage: r > p1 and r > p2 and r > p3
            True
            sage: p1 < r and p2 < r and p3 < r
            True

            sage: p1 = iet.Permutation('a b', 'b a', alphabet='ab')
            sage: p2 = iet.Permutation('a b', 'b a', alphabet='ba')
            sage: p1 < p2
            Traceback (most recent call last):
            ...
            ValueError: comparison of permutations over different alphabets
        """
        if type(self) != type(other):
            raise TypeError("Permutations must be of the same type")

        if self._labels is not None and other._labels is not None:
            if self._alphabet is not other._alphabet:
                raise ValueError('comparison of permutations over different alphabets')

        if len(self) < len(other):
            return True
        elif len(self) > len(other):
            return False

        if self._twin < other._twin:
            return True
        elif self._twin > other._twin:
            return False

        if self._labels is not None and other._labels is not None:
            if self._labels < other._labels:
                return True
            elif self._labels > other._labels:
                return False

        if self._flips is not None:
            if self._flips < other._flips:
                return True
            elif self._flips > other._flips:
                return False

        # equality
        return False

    def _check(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: p = iet.Permutation('a b c d', 'd a c b', flips=['a','b'])
            sage: p._check()
        """
        if set(self.letters()) != set().union(*self.list()):
            raise RuntimeError("letters are bad")

        for i in range(2):
            for p in range(self.length(i)):
                j,q = self.twin(i,p)
                if self.twin(j,q) != (i,p):
                    raise RuntimeError("self.twin is not an involution: i={} p={} j={} q={}".format(i,p,j,q))
                if self[i][p] != self[j][q]:
                    raise RuntimeError("uncoherent getitem i={} p={} j={} q={}".format(i,p,j,q))
                if self._labels is not None:
                    if self._labels[i][p] != self._labels[j][q]:
                        raise RuntimeError("self._labels wrong i={} p={} j={} q={}".format(i,p,j,q))
                if self._flips is not None:
                    if self._flips[i][p] != self._flips[j][q]:
                        raise RuntimeError("self._flips wrong i={} p={} j={} q={}".format(i,p,j,q))

    def letters(self):
        r"""
        Returns the list of letters of the alphabet used for representation.

        The letters used are not necessarily the whole alphabet (for example if
        the alphabet is infinite).

        OUTPUT: a list of labels

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation([1,2],[2,1])
            sage: p.alphabet(Alphabet(name="NN"))
            sage: p
            0 1
            1 0
            sage: p.letters()
            [0, 1]

            sage: p = iet.GeneralizedPermutation('a a b','b c c')
            sage: p.letters()
            ['a', 'b', 'c']
            sage: p.alphabet(range(10))
            sage: p.letters()
            [0, 1, 2]
            sage: p._remove_interval(0, 2)
            sage: p.letters()
            [0, 2]

            sage: p = iet.GeneralizedPermutation('a a b', 'b c c', reduced=True)
            sage: p.letters()
            ['a', 'b', 'c']

        For permutations with flips, the letters appear as pairs made of an element
        of the alphabet and the flip::

            sage: p = iet.Permutation('A B C D', 'D C A B', flips='AC')
            sage: p.letters()
            ['A', 'B', 'C',  'D']

            sage: p = iet.GeneralizedPermutation('A A B', 'B C C', flips='B')
            sage: p.letters()
            ['A', 'B', 'C']
        """
        unrank = self._alphabet.unrank
        if self._labels is not None:
            labels = set(self._labels[0]).union(self._labels[1])
            return [unrank(l) for l in sorted(labels)]
        else:
            return [unrank(i) for i in range(len(self))]

    def _repr_(self):
        r"""
        Representation method of self.

        Apply the function str to _repr_type(_repr_options) if _repr_type is
        callable and _repr_type else.

        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c','c b a')
            sage: p._repr_type = 'str'
            sage: p._repr_options = ('\n',)
            sage: p   #indirect doctest
            a b c
            c b a
            sage: p._repr_options = (' / ',)
            sage: p   #indirect doctest
            a b c / c b a

        ::

            sage: p._repr_type = '_twin'
            sage: p   #indirect doctest
            [[2, 1, 0], [2, 1, 0]]
        """
        if self._repr_type is None:
            return ''

        elif self._repr_type == 'reduced':
            return ''.join(map(str,self[1]))

        else:
            f = getattr(self, self._repr_type)
            if callable(f):
                return str(f(*self._repr_options))
            else:
                return str(f)

    def __hash__(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: q = iet.Permutation('a b c', 'c b a', flips=['a'])
            sage: r = iet.Permutation('a b c', 'c b a', reduced=True)
            sage: s = iet.Permutation('a b c', 'c b a', flips=['a'], reduced=True)
            sage: len(set([hash(p), hash(q), hash(r), hash(s)]))
            4

            sage: P = [iet.Permutation(t, b, flips=f, reduced=False) \
            ....:        for t in ['a b', 'b a'] \
            ....:        for b in ['a b', 'b a'] \
            ....:        for f in [None, 'a', 'b', 'ab']]
            sage: P.extend(iet.Permutation('a b', b, flips=f, reduced=True) \
            ....:        for b in ['a b', 'b a'] \
            ....:        for f in [None, 'a', 'b', 'ab'])
            sage: len(set(hash(x) for x in P))
            24
        """
        h = hash(tuple(self._twin[0] + self._twin[1]))
        if self._labels is not None:
            h *= 7461864723258187525
            h ^= hash(tuple(self._labels[0] + self._labels[1]))
            h += hash(tuple(self.letters()))
        if self._flips is not None:
            h *= 5566797465546889505
            h ^= hash(tuple(self._flips[0] + self._flips[1]))
        return h

    def str(self, sep= "\n", align=None):
        r"""
        A string representation of the generalized permutation.

        INPUT:

        - ``sep`` - (default: '\n') a separator for the two intervals

        OUTPUT:

        string -- the string that represents the permutation


        EXAMPLES::

            sage: from surface_dynamics import *

        For permutations of iet::

            sage: p = iet.Permutation('a b c','c b a')
            sage: p.str()
            'a b c\nc b a'
            sage: p.str(sep=' | ')
            'a b c | c b a'

        The permutation can be rebuilt from the standard string::

            sage: p == iet.Permutation(p.str())
            True

        For permutations of li::

            sage: p = iet.GeneralizedPermutation('a b b','c c a')
            sage: p.str()
            'a b b\nc c a'
            sage: p.str(sep=' | ')
            'a b b | c c a'

            sage: iet.GeneralizedPermutation([0,0], [1,2,3,2,1,3])
            0 0
            1 2 3 2 1 3
            sage: iet.GeneralizedPermutation([0,1,2,1,0,2], [3,3])
            0 1 2 1 0 2
            3 3

        Again, the generalized permutation can be rebuilt from the standard string::

            sage: p == iet.GeneralizedPermutation(p.str())
            True

        With flips::

            sage: p = iet.GeneralizedPermutation('a a','b b',flips='a')
            sage: print(p.str())
            -a -a
             b  b
             sage: print(p.str('/'))
             -a -a/ b  b

        Alignment with alphabet of different sizes::

            sage: p = iet.Permutation('aa b ccc d', 'd b ccc aa')
            sage: print(p.str())
            aa b ccc d
            d b ccc aa
            sage: print(p.str(align='left'))
            aa b ccc d
            d  b ccc aa
            sage: print(p.str(align='right'))
            aa b ccc  d
             d b ccc aa

            sage: p = iet.GeneralizedPermutation('aa fff b ccc b fff d', 'eee d eee ccc aa')
            sage: print(p.str())
            aa fff b ccc b fff d
            eee d eee ccc aa
            sage: print(p.str(align='left'))
            aa  fff b   ccc b  fff d
            eee d   eee ccc aa
            sage: print(p.str(align='right'))
             aa fff   b ccc  b fff d
            eee   d eee ccc aa
        """
        s = []
        if self._flips is None:
            l0, l1 = self.list()
            formatter = str
        else:
            l0, l1 = self.list(flips=True)
            formatter = labelize_flip

        n = min(len(l0), len(l1))
        for i in range(n):
            l0[i] = s0 = formatter(l0[i])
            l1[i] = s1 = formatter(l1[i])
            if align is None:
                continue
            elif len(s0) < len(s1):
                if align == 'right':
                    l0[i] = ' ' * (len(s1) - len(s0)) + l0[i]
                elif align == 'left':
                    l0[i] = l0[i] + ' ' * (len(s1) - len(s0))
            elif len(s0) > len(s1):
                if align == 'right':
                    l1[i] = ' ' * (len(s0) - len(s1)) + l1[i]
                elif align == 'left':
                    l1[i] = l1[i] + ' ' * (len(s0) - len(s1))
        for i in range(n, len(l0)):
            l0[i] = formatter(l0[i])
        for i in range(n, len(l1)):
            l1[i] = formatter(l1[i])

        return sep.join([' '.join(l0), ' '.join(l1)])

    def flips(self):
        r"""
        Returns the list of flips.

        If the permutation is not a flipped permutations then ``None`` is returned.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: iet.Permutation('a b c', 'c b a').flips()
            []
            sage: iet.Permutation('a b c', 'c b a', flips='ac').flips()
            ['a', 'c']
            sage: iet.GeneralizedPermutation('a a', 'b b', flips='a').flips()
            ['a']
            sage: iet.GeneralizedPermutation('a a','b b', flips='b', reduced=True).flips()
            ['b']
        """
        if self._flips is None:
            return []

        res = []
        l = self.list(flips=False)
        letters = []
        for i,f in enumerate(self._flips[0]):
            if f == -1 and l[0][i] not in letters:
                res.append(l[0][i])
                letters.append(l[0][i])
        for i,f in enumerate(self._flips[1]):
            if f == -1 and l[1][i] not in letters:
                res.append(l[1][i])
                letters.append(l[1][i])
        return letters

    def __copy__(self) :
        r"""
        Returns a copy of self.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c', 'c b a', reduced=True)
            sage: q = copy(p)
            sage: p == q
            True
            sage: p is q
            False
            sage: p = iet.GeneralizedPermutation('a b b','c c a', reduced=True)
            sage: q = copy(p)
            sage: p == q
            True
            sage: p is q
            False

        TESTS::

            sage: p1 = iet.Permutation('a b c', 'c b a', reduced=True)
            sage: p2 = iet.Permutation('a b c', 'c b a', reduced=False)
            sage: p3 = iet.Permutation('a b c', 'c b a', flips='a', reduced=True)
            sage: p4 = iet.Permutation('a b c', 'c b a', flips='a', reduced=False)

            sage: p5 = iet.GeneralizedPermutation('a a b', 'b c c', reduced=True)
            sage: p6 = iet.GeneralizedPermutation('a a b', 'b c c', reduced=False)
            sage: p7 = iet.GeneralizedPermutation('a a b', 'b c c', flips='a', reduced=True)
            sage: p8 = iet.GeneralizedPermutation('a a b', 'b c c', flips='a', reduced=False)

            sage: for p in (p1,p2,p3,p4,p5,p6,p7,p8):
            ....:     q = p.__copy__()
            ....:     assert type(p) is type(q)
            ....:     assert p == q
        """
        q = self.__class__.__new__(self.__class__)

        q._twin = [self._twin[0][:], self._twin[1][:]]
        q._alphabet = self._alphabet
        q._repr_type = self._repr_type
        q._repr_options = self._repr_options

        if self._labels is not None:
            q._labels = [self._labels[0][:], self._labels[1][:]]

        if self._flips is not None:
            q._flips = [self._flips[0][:], self._flips[1][:]]

        return q

    _repr_type = 'str'
    _repr_options = ("\n",)

    def __len__(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b','b a',reduced=True)
            sage: len(p)
            2
            sage: p = iet.Permutation('a b','b a',reduced=False)
            sage: len(p)
            2
            sage: p = iet.GeneralizedPermutation('a a b','b c c',reduced=True)
            sage: len(p)
            3
            sage: p = iet.GeneralizedPermutation('a a b','b c c',reduced=False)
            sage: len(p)
            3
            sage: p = iet.GeneralizedPermutation('a a','b b c c',reduced=True)
            sage: len(p)
            3
            sage: p = iet.GeneralizedPermutation('a a','b b c c',reduced=False)
            sage: len(p)
            3
        """
        return (len(self._twin[0]) + len(self._twin[1])) // 2

    def length_top(self):
        r"""
        Returns the number of intervals in the top segment.

        OUTPUT:

        integer -- the length of the top segment

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c','c b a', reduced=True)
            sage: p.length_top()
            3
            sage: p = iet.Permutation('a b c','c b a', reduced=False)
            sage: p.length_top()
            3
            sage: p = iet.GeneralizedPermutation('a a b','c d c b d',reduced=True)
            sage: p.length_top()
            3
            sage: p = iet.GeneralizedPermutation('a a b','c d c d b',reduced=False)
            sage: p.length_top()
            3
            sage: p = iet.GeneralizedPermutation('a b c b d c d', 'e a e',reduced=True)
            sage: p.length_top()
            7
            sage: p = iet.GeneralizedPermutation('a b c d b c d', 'e a e', reduced=False)
            sage: p.length_top()
            7
        """
        return len(self._twin[0])

    def length_bottom(self):
        r"""
        Returns the number of intervals in the bottom segment.

        OUTPUT:

        integer -- the length of the bottom segment

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c','c b a',reduced=True)
            sage: p.length_bottom()
            3
            sage: p = iet.Permutation('a b c','c b a',reduced=False)
            sage: p.length_bottom()
            3
            sage: p = iet.GeneralizedPermutation('a a b','c d c b d',reduced=True)
            sage: p.length_bottom()
            5
            sage: p = iet.GeneralizedPermutation('a a b','c d c b d',reduced=False)
            sage: p.length_bottom()
            5
        """
        return len(self._twin[1])

    def length(self, interval=None):
        r"""
        Returns the 2-uple of lengths.

        p.length() is identical to (p.length_top(), p.length_bottom())
        If an interval is specified, it returns the length of the specified
        interval.

        INPUT:

        - ``interval`` - None, 'top' (or 't' or 0) or 'bottom' (or 'b' or 1)

        OUTPUT:

        integer or 2-uple of integers -- the corresponding lengths

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c','c b a',reduced=False)
            sage: p.length()
            (3, 3)
            sage: p = iet.Permutation('a b c','c b a',reduced=True)
            sage: p.length()
            (3, 3)

            sage: p = iet.GeneralizedPermutation('a a b','c d c b d',reduced=False)
            sage: p.length()
            (3, 5)
            sage: p = iet.GeneralizedPermutation('a a b','c d c b d',reduced=True)
            sage: p.length()
            (3, 5)
        """
        if interval is None :
            return len(self._twin[0]),len(self._twin[1])
        else :
            interval = interval_conversion(interval)
            return len(self._twin[interval])

    def _set_alphabet(self, alphabet):
        r"""
        Sets the alphabet of self.

        TESTS:

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('a a','b b')
            sage: p.alphabet([0,1])   #indirect doctest
            sage: p.alphabet() == Alphabet([0,1])
            True
            sage: p
            0 0
            1 1
            sage: p.alphabet("cd")   #indirect doctest
            sage: p.alphabet() == Alphabet(['c','d'])
            True
            sage: p
            c c
            d d

        Tests with reduced permutations::

            sage: p = iet.Permutation('a b','b a',reduced=True)
            sage: p.alphabet([0,1])   #indirect doctest
            sage: p.alphabet() == Alphabet([0,1])
            True
            sage: p
            0 1
            1 0
            sage: p.alphabet("cd")   #indirect doctest
            sage: p.alphabet() == Alphabet(['c','d'])
            True
            sage: p
            c d
            d c

        ::

            sage: p = iet.GeneralizedPermutation('a a','b b',reduced=True)
            sage: p.alphabet([0,1])   #indirect doctest
            sage: p.alphabet() == Alphabet([0,1])
            True
            sage: p
            0 0
            1 1
            sage: p.alphabet("cd")   #indirect doctest
            sage: p.alphabet() == Alphabet(['c','d'])
            True
            sage: p
            c c
            d d
        """
        from sage.combinat.words.alphabet import build_alphabet
        alphabet = build_alphabet(alphabet)
        if alphabet.cardinality() < len(self):
            raise ValueError("not enough letters in alphabet")
        self._alphabet = alphabet

    def alphabet(self, data=None):
        r"""
        Manages the alphabet of self.

        If there is no argument, the method returns the alphabet used. If the
        argument could be converted to an alphabet, this alphabet will be used.

        INPUT:

        - ``data`` - None or something that could be converted to an alphabet


        OUTPUT:

        - either None or the current alphabet


        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b','a b')
            sage: p.alphabet([0,1])
            sage: p.alphabet() == Alphabet([0,1])
            True
            sage: p
            0 1
            0 1
            sage: p.alphabet("cd")
            sage: p.alphabet() == Alphabet(['c','d'])
            True
            sage: p
            c d
            c d
        """
        if data is None:
            return self._alphabet
        else:
            self._set_alphabet(data)

    def left_right_inverse(self, inplace=False):
        r"""
        Return the left-right inverse.

        The left-right inverse of a permutation, is the permutation obtained by
        reversing the order of the underlying ordering.

        You can also use the shorter .lr_inverse()

        There are two other symmetries of permutation which are accessible via
        the methods
        :meth:`Permutation.top_bottom_inverse` and
        :meth:`Permutation.symmetric`.

        OUTPUT: a permutation

        EXAMPLES::

            sage: from surface_dynamics import *

        For labelled permutations::

            sage: p = iet.Permutation('a b c','c a b')
            sage: p.left_right_inverse()
            c b a
            b a c
            sage: p = iet.Permutation('a b c d','c d a b')
            sage: p.left_right_inverse()
            d c b a
            b a d c

        for reduced permutations::

            sage: p = iet.Permutation('a b c','c a b',reduced=True)
            sage: p.left_right_inverse()
            a b c
            b c a
            sage: p = iet.Permutation('a b c d','c d a b',reduced=True)
            sage: p.left_right_inverse()
            a b c d
            c d a b

        for labelled quadratic permutations::

             sage: p = iet.GeneralizedPermutation('a a','b b c c')
             sage: p.left_right_inverse()
             a a
             c c b b

        for reduced quadratic permutations::

             sage: p = iet.GeneralizedPermutation('a a','b b c c',reduced=True)
             sage: p.left_right_inverse() == p
             True
        """
        res = self if inplace else copy(self)
        res._reversed_twin()

        if res._flips is not None:
            res._flips[0].reverse()
            res._flips[1].reverse()

        if res._labels is not None:
            res._labels[0].reverse()
            res._labels[1].reverse()

        return res

    lr_inverse = left_right_inverse
    vertical_inverse = left_right_inverse

    def top_bottom_inverse(self, inplace=False):
        r"""
        Return the top-bottom inverse.

        You can use also use the shorter .tb_inverse().

        There are two other symmetries of permutation which are accessible via
        the methods
        :meth:`Permutation.left_right_inverse` and
        :meth:`Permutation.symmetric`.

        OUTPUT: a permutation

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b','b a')
            sage: p.top_bottom_inverse()
            b a
            a b
            sage: p = iet.Permutation('a b','b a',reduced=True)
            sage: p.top_bottom_inverse() == p
            True

        ::

            sage: p = iet.Permutation('a b c d','c d a b')
            sage: p.top_bottom_inverse()
            c d a b
            a b c d

        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b','a b')
            sage: p == p.top_bottom_inverse()
            True
            sage: p is p.top_bottom_inverse()
            False
            sage: p = iet.GeneralizedPermutation('a a','b b',reduced=True)
            sage: p == p.top_bottom_inverse()
            True
            sage: p is p.top_bottom_inverse()
            False
        """
        res = self if inplace else copy(self)
        res._inversed_twin()

        if res._flips is not None:
            res._flips = [res._flips[1], res._flips[0]]

        if res._labels is not None:
            res._labels = [res._labels[1], res._labels[0]]

        return res

    tb_inverse = top_bottom_inverse
    horizontal_inverse = top_bottom_inverse

    def symmetric(self):
        r"""
        Returns the symmetric permutation.

        The symmetric permutation is the composition of the top-bottom
        inversion and the left-right inversion (which are geometrically
        orientation reversing).

        There are two other symmetries of permutation which are accessible via
        the methods
        :meth:`Permutation.left_right_inverse` and
        :meth:`Permutation.top_bottom_inverse`.

        OUTPUT: a permutation

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c','c b a',reduced=True)
            sage: p.symmetric() == p
            True
            sage: p = iet.Permutation('a b c d','c a d b',reduced=True)
            sage: q = p.symmetric()
            sage: q
            a b c d
            b d a c
            sage: q1 = p.tb_inverse().lr_inverse()
            sage: q2 = p.lr_inverse().tb_inverse()
            sage: q == q1 and q == q2
            True

        It works for any type of permutations::

            sage: p = iet.GeneralizedPermutation('a b b','c c a',flips='ab')
            sage: p
            -a -b -b
             c  c -a
            sage: p.symmetric()
            -a  c  c
            -b -b -a

        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('a a b','b c c',reduced=True)
            sage: q = p.symmetric()
            sage: q1 = p.tb_inverse().lr_inverse()
            sage: q2 = p.lr_inverse().tb_inverse()
            sage: q == q1 and q == q2
            True
        """
        res = copy(self)
        res._inversed_twin()
        res._reversed_twin()

        if res._flips is not None:
            res._flips[0].reverse()
            res._flips[1].reverse()
            res._flips.reverse()

        if res._labels is not None:
            res._labels[0].reverse()
            res._labels[1].reverse()
            res._labels.reverse()

        return res

    def _canonical_signs(self):
        r"""
        Return label positions and orientation in some compatible way.

        OUTPUT:

        - a dictionary whose keys are the labels of this permutation and the
          value associated to a label `a` is a pair `((ip,jp), (im,jm))` made
          of the two positions of `a`.

        - a list of two lists that give the signs of each subinterval

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('d a e a b b','e c c d')
            sage: twins, orientation = p._canonical_signs()
            sage: twins
            {'a': [(0, 1), (0, 3)],
             'b': [(0, 4), (0, 5)],
             'c': [(1, 2), (1, 1)],
             'd': [(0, 0), (1, 3)],
             'e': [(0, 2), (1, 0)]}
            sage: orientation
            [[1, 1, 1, -1, 1, -1], [-1, -1, 1, -1]]

            sage: for a, ((ip,jp),(im,jm)) in twins.items():
            ....:     assert p[ip][jp] == p[im][jm] == a
            ....:     assert orientation[ip][jp] == +1
            ....:     assert orientation[im][jm] == -1

            sage: p = iet.Permutation('a b c', 'c b a', flips='ac')
            sage: twins, orientation = p._canonical_signs()
            sage: twins
            {'a': [(0, 0), (1, 2)],
             'b': [(0, 1), (1, 1)],
             'c': [(0, 2), (1, 0)]}
            sage: orientation
            [[1, 1, 1], [1, -1, 1]]

            sage: p = iet.GeneralizedPermutation('a a', 'b b c c', flips='b')
            sage: twins, orientation = p._canonical_signs()
            sage: twins
            {'a': [(0, 0), (0, 1)],
             'b': [(1, 1), (1, 0)],
             'c': [(1, 3), (1, 2)]}
            sage: orientation
            [[1, -1], [1, 1, -1, 1]]
        """
        flips = self._flips
        letters = set(self.letters())

        label_to_twins = {}
        sign_top = []
        sign_bot = []
        for i,a in enumerate(self[0]):
            if a in letters:  # first time seen
                sign_top.append(1)
                label_to_twins[a] = [(0,i)]
                letters.remove(a)
            else:             # already seen
                if flips is not None and flips[0][i] == -1:
                    sign_top.append(1)
                else:
                    sign_top.append(-1)
                label_to_twins[a].append((0,i))

        n = len(self[1])
        for i,a in enumerate(reversed(self[1])):
            i = n - 1 - i
            if a in letters:  # first time seen
                sign_bot.append(1)
                letters.remove(a)
                label_to_twins[a] = [(1,i)]
            else:             # already seen
                if flips is not None and flips[1][i] == -1:
                    sign_bot.append(1)
                else:
                    sign_bot.append(-1)
                label_to_twins[a].append((1,i))
        sign_bot.reverse()

        return (label_to_twins, [sign_top, sign_bot])

    def interval_diagram(self, glue_ends=True, sign=False):
        r"""
        Return the interval diagram of self.

        The result is in random order.

        INPUT:

        - ``glue_ends`` - bool (default: True)

        - ``sign`` - bool (default: False)

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c','c b a')
            sage: p.interval_diagram()
            [['b', ('c', 'a')], ['b', ('c', 'a')]]

            sage: p = iet.Permutation('a b c','c a b')
            sage: p.interval_diagram()
            [['a', 'b'], [('c', 'b'), ('c', 'a')]]

            sage: p = iet.GeneralizedPermutation('a a','b b c c')
            sage: p.interval_diagram()
            [['a'], [('b', 'a', 'c')], ['b'], ['c']]

            sage: p = iet.GeneralizedPermutation('a a b b','c c')
            sage: p.interval_diagram()
            [['a'], [('b', 'c', 'a')], ['b'], ['c']]
            sage: p.interval_diagram(sign=True)
            [[('a', -1)], [(('b', 1), ('c', 1), ('a', 1))], [('b', -1)], [('c', -1)]]

            sage: p = iet.GeneralizedPermutation((0,1,0,2),(3,2,4,1,4,3))
            sage: p.interval_diagram()
            [[0, 1, 4, 1], [2, 3, 4, (2, 3, 0)]]
            sage: p.interval_diagram(sign=True)
            [[(0, -1), (1, 1), (4, -1), (1, -1)],
             [(2, 1), (3, -1), (4, 1), ((2, -1), (3, 1), (0, 1))]]

            sage: p = iet.GeneralizedPermutation('a b c d b', 'e d f e a f c', flips='bdf')
            sage: p.interval_diagram()
            [['a', 'b', 'd', 'e', 'f', 'c', ('b', 'c'), 'd', 'f'], [('e', 'a')]]

        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('0 1 2 3 2','4 3 4 1 0')
            sage: p.interval_diagram(sign=True)
            [[('1', 1), (('4', 1), ('0', 1)), ('1', -1), (('2', 1), ('0', -1))],
             [('2', -1), ('3', 1), ('4', -1), ('3', -1)]]

            sage: p = iet.Permutation('a b c', 'c b a', flips='a')
            sage: p.interval_diagram(glue_ends=False, sign=True)
            [[('a', -1), ('b', -1), ('c', 1), ('a', 1), ('c', -1), ('b', 1)]]
            sage: p.interval_diagram(glue_ends=True, sign=False)
            [['a', 'b', ('c', 'a', 'c'), 'b']]

            sage: iet.Permutation('a b c', 'c b a', flips='b').interval_diagram(glue_ends=False, sign=True)
            [[('a', -1), ('b', 1), ('a', 1), ('c', 1), ('b', -1), ('c', -1)]]
            sage: iet.Permutation('a b c', 'c b a', flips='c').interval_diagram(glue_ends=False, sign=True)
            [[('a', -1), ('b', 1), ('c', 1), ('b', -1), ('a', 1), ('c', -1)]]

            sage: iet.Permutation('a b c', 'c b a', flips='bc').interval_diagram(glue_ends=False, sign=True)
            [[('a', -1), ('b', 1), ('a', 1), ('c', -1)], [('b', -1), ('c', 1)]]

            sage: iet.Permutation('a b c', 'c b a', flips='ac').interval_diagram(glue_ends=False, sign=True)
            [[('a', -1), ('b', -1), ('c', 1), ('b', 1)], [('a', 1), ('c', -1)]]

            sage: iet.Permutation('a b c', 'c b a', flips='ab').interval_diagram(glue_ends=False, sign=True)
            [[('a', -1), ('b', 1)], [('a', 1), ('c', -1), ('b', -1), ('c', 1)]]

            sage: iet.Permutation('a b c', 'c b a', flips='abc').interval_diagram(glue_ends=False, sign=True)
            [[('a', -1), ('b', 1)], [('a', 1), ('c', -1)], [('b', -1), ('c', 1)]]
        """
        # NOTE: each interval comes with an orientation (obtained from the
        # method ._canonical_signs(). We label each side of each interval
        # by either +1 or -1 as follows
        #
        #   o---------------->--------------o
        #  +1 (for source)       -1 (for target)
        #
        flips = self._flips
        twin = self.twin_list()
        labels = self.list()
        letters = set((label,j) for label in self.letters() for j in (+1,-1))
        label_to_twins, orientation = self._canonical_signs()
        m0 = len(labels[0])
        m1 = len(labels[1])

        singularities = []

        just_glued = False     # True iff the last elt in singularity is paired
        glued_at_begin = False # True iff the 1st elt in singularity is paired
        flip = 1
        while letters:
            # pick a remaining letter
            try:
                # try picking the minimum
                t = min(letters)
                letters.remove(t)
                label, j = t
            except TypeError:
                # pick a random one
                print('failure')
                label,j = letters.pop()
            (i0,p0),(i1,p1) = label_to_twins[label]  # its two positions in the interval

            # try to see if one of (i0, p0) fits for anti-clockwise order
            # otherwise start with clockwise direction
            if j == -1:
                if orientation[i0][p0] == -1:
                    i,p = i0,p0
                    flip = 1
                elif orientation[i1][p1] == -1:
                    i,p = i1,p1
                    flip = 1
                else:
                    i,p = i0,p0
                    flip = -1
            else:
                if orientation[i0][p0] == +1:
                    i,p = i0,p0
                    flip = 1
                elif orientation[i1][p1] == +1:
                    i,p = i1,p1
                    flip = 1
                else:
                    i,p = i0,p0
                    flip = -1


            if sign:
                singularity = [(label,j)]
            else:
                singularity = [label]
            while True:
                i,p = twin[i][p]
                if flips is not None and flips[i][p] == -1:
                    flip *= -1
                if i == 0:      # twin on top
                    p += flip   # next interval
                    if flip == 1 and p == m0:     # at the right end?
                        i = 1
                        p = m1-1
                    elif flip == -1 and p == -1:  # at the left end?
                        i = 1
                        p = 0
                else:           # twin on bot
                    p -= flip   # next interval
                    if flip == 1 and p == -1:     # at the left end?
                        i = 0
                        p = 0
                    elif flip == -1 and p == m1:  # at the right end?
                        i = 0
                        p = m0 - 1

                label = labels[i][p]
                j = flip * orientation[i][p]


                if (label,j) not in letters:
                    if (glue_ends and
                        ((i == 1 and p == m1-1 and flip == +1) or (i == 0 and p == 0 and flip == +1) or
                         (i == 0 and p == m1-1 and flip == -1) or (i == 1 and p == 0 and flip == -1))):
                        sg2 = singularity.pop(0)
                        sg1 = singularity.pop(-1)
                        if glued_at_begin:
                            singularity.append((sg1,) + sg2)
                        elif just_glued:
                            singularity.append(sg1 + (sg2,))
                        else:
                            singularity.append((sg1, sg2))
                    break
                letters.remove((label,j))
                if (glue_ends and (
                    (i == 1 and p == m1-1 and flip == +1) or (i == 0 and p == 0 and flip == 1) or
                    (i == 0 and p == m1-1 and flip == -1) or (i == 1 and p == 0 and flip == -1))):
                    sg1 = singularity.pop(-1)
                    if sign:
                        sg2 = (label,j)
                    else:
                        sg2 = label
                    if len(singularity) == 0:
                        glued_at_begin = True
                    if just_glued:
                        singularity.append(sg1 + (sg2,))
                    else:
                        singularity.append((sg1,sg2))
                    just_glued = True
                else:
                    if sign:
                        singularity.append((label,j))
                    else:
                        singularity.append(label)
                    just_glued = False

            singularities.append(singularity)
            just_glued = False
            glued_at_begin = False

        return singularities

    def _remove_interval(self, i, pos):
        r"""
        Remove the letter in the interval ``i`` and position ``pos`` (and its
        twin).

        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c d', 'd a b c')
            sage: p._remove_interval(0, 1); p
            a c d
            d a c
            sage: p._check()
            sage: p.twin_list()
            [[(1, 1), (1, 2), (1, 0)], [(0, 2), (0, 0), (0, 1)]]
            sage: p._remove_interval(0, 2); p
            a c
            a c
            sage: p._check()
            sage: p.twin_list()
            [[(1, 0), (1, 1)], [(0, 0), (0, 1)]]
            sage: p._remove_interval(0, 0); p
            c
            c
            sage: p.twin_list()
            [[(1, 0)], [(0, 0)]]

            sage: p = iet.Permutation('a b c d', 'd c b a', reduced=True)
            sage: p._remove_interval(0, 1); p
            a b c
            c b a
            sage: p._check()
            sage: p._remove_interval(1, 0); p
            a b
            b a
            sage: p._check()

            sage: p = iet.Permutation('a b c d', 'd a c b', flips=['a','b'])
            sage: p._remove_interval(0, 0); p
            -b  c  d
             d  c -b
            sage: p._check()
            sage: p._remove_interval(1, 0); p
            -b  c
             c -b
            sage: p._check()

            sage: p = iet.Permutation('a b c e d', 'e b d a c')
            sage: p._remove_interval(0, 3); p
            a b c d
            b d a c
            sage: p._check()
            sage: p._remove_interval(1, 3); p
            a b d
            b d a
            sage: p._check()
            sage: p._remove_interval(1, 0); p
            a d
            d a
            sage: p._check()

            sage: p = iet.GeneralizedPermutation('a a b c b e d e', 'f c d f g g')
            sage: p._remove_interval(0, 0)
            sage: p._check()
            sage: p._remove_interval(0, 4)
            sage: p._check()
            sage: p._remove_interval(1, 4)
            sage: p._check()
            sage: p
            b c b e e
            f c f
        """
        assert i == 0 or i == 1
        assert 0 <= pos < self.length(i)

        ii, ppos = self.twin(i, pos)
        if i == ii and pos < ppos:
            pos, ppos = ppos, pos

        for j,t in enumerate(self.twin_list()):
            for q,(jj,qq) in enumerate(t):
                dec = (jj == i and qq > pos) + (jj == ii and qq > ppos)
                if dec:
                    self._set_twin(j, q, jj, qq - dec)

        del self._twin[i][pos]
        del self._twin[ii][ppos]

        if self._flips is not None:
            del self._flips[i][pos]
            del self._flips[ii][ppos]

        if self._labels is not None:
            del self._labels[i][pos]
            del self._labels[ii][ppos]

    def _identify_intervals(self, side):
        r"""
        Identify the two intervals on the right (if ``side=-1``) or on the left
        (if ``side=0``).

        This is needed for the generalized Rauzy induction when the length on
        top and bottom are equal.

        The label of the top interval disappears.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c', 'c a b')
            sage: p._identify_intervals(-1); p
            a b
            b a
            sage: p._check()

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: p._identify_intervals(0); p
            b c
            b c
            sage: p._check()

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: p._identify_intervals(0); p
            b c
            b c
            sage: p._check()

            sage: p = iet.Permutation('a b c d', 'c a d b')
            sage: p._identify_intervals(-1); p
            a b c
            c a b
            sage: p._check()

            sage: p = iet.Permutation('a b c d', 'd a c b')
            sage: p._identify_intervals(-1); p
            a b c
            b a c
            sage: p._check()

            sage: p = iet.GeneralizedPermutation('a d e e a b', 'b c c d')
            sage: p._identify_intervals(0); p
            d e e b b
            c c d
            sage: p._check()
            sage: p._identify_intervals(-1); p
            d e e d
            c c
            sage: p._check()

            sage: p = iet.GeneralizedPermutation('a b b', 'a c c')
            sage: p._identify_intervals(0); p
            b b
            c c
            sage: p._check()
        """
        if self._flips:
            raise ValueError("not implemented if there is a flip")

        if side == 0:
            pos0 = pos1 = 0
        elif side == -1:
            pos0 = self.length_top()-1
            pos1 = self.length_bottom()-1
        twin = self.twin(0, pos0)
        if self.twin(0, pos0) != (1, pos1):
            self._move(1, pos1, twin[0], twin[1])
        self._remove_interval(0, pos0)

    def cover(self, perms, as_tuple=False):
        r"""
        Return a covering of this permutation.

        INPUT:

        - ``perms`` - a list of permutations that describe the gluings

        - ``as_tuple`` - whether permutations need to be considered as `1`-based
          (default) or `0`-based.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: p = iet.Permutation('a b', 'b a')
            sage: p.cover(['(1,2)', '(1,3)'])
            Covering of degree 3 of the permutation:
            a b
            b a

            sage: p.cover([[1,0,2], [2,1,0]], as_tuple=True)
            Covering of degree 3 of the permutation:
            a b
            b a

            sage: p = iet.GeneralizedPermutation('a a b b','c c')
            sage: q = p.cover(['(0,1)', [], [2,1,0]], as_tuple=True)
            sage: q
            Covering of degree 3 of the permutation:
            a a b b
            c c
            sage: q.covering_data('a')
            (1,2)
            sage: q.covering_data('b')
            ()
            sage: q.covering_data('c')
            (1,3)
        """
        if len(perms) != len(self):
            raise ValueError("wrong number of permutations")

        from surface_dynamics.misc.permutation import perm_init, equalize_perms
        if as_tuple:
            perms = [perm_init(p) for p in perms]
        else:
            try:
                # Trac #28652: Rework the constructor of PermutationGroupElement
                from sage.groups.perm_gps.constructor import PermutationGroupElement
            except ImportError:
                from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
            perms = [PermutationGroupElement(p,check=True) for p in perms]
            perms = [[i-1 for i in p.domain()] for p in perms]

        d = equalize_perms(perms)
        from .cover import PermutationCover
        return PermutationCover(self, d, perms)

    def regular_cover(self, G, elts):
        r"""
        Return a regular (or normal) cover of this permutation.

        EXAMPLES::

            sage: from surface_dynamics import iet

            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: G = SymmetricGroup(4)
            sage: ga = G('(1,2)')
            sage: gb = G('(1,3,4)')
            sage: gc = G('(1,3)')
            sage: gd = G('()')
            sage: pp = p.regular_cover(G, [ga, gb, gc, gd])
            sage: pp
            Regular cover of degree 24 with group Symmetric group of order 4! as a permutation group of the permutation:
            a b c d
            d c b a
            sage: pp.genus()
            33

        Additive groups are converted to multiplicative groups::

            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: G = Zmod(5)
            sage: p.regular_cover(G, [0, 1, 2])
            Regular cover of degree 5 with group Multiplicative Abelian group isomorphic to C5 of the permutation:
            a a b
            b c c
        """
        from .cover import RegularCover
        return RegularCover(self, G, elts)

class PermutationIET(Permutation):
    def _init_twin(self, a):
        r"""
        Initializes the twin list.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: iet.Permutation('a b','b a',reduced=True) #indirect doctest
            a b
            b a
            sage: iet.Permutation('a b c','c a b',reduced=True) #indirect doctest
            a b c
            c a b
        """
        if a is None:
            self._twin = [[],[]]

        else:
            self._twin = [[0]*len(a[0]), [0]*len(a[1])]
            for i in range(len(a[0])) :
                j = a[1].index(a[0][i])
                self._twin[0][i] = j
                self._twin[1][j] = i

    def twin(self, i, pos):
        r"""
        Return the twin of the interval in the interval ``i`` at position
        ``pos``.

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: p = iet.Permutation('a b c e d', 'e b d a c')
            sage: p.twin(0,0)
            (1, 3)
            sage: p.twin(0,1)
            (1, 1)

            sage: twin_top = [p.twin(0,i) for i in range(p.length_top())]
            sage: twin_bot = [p.twin(1,i) for i in range(p.length_bottom())]
            sage: p.twin_list() == [twin_top, twin_bot]
            True
        """
        return (1-i, self._twin[i][pos])

    def _set_twin(self, i, p, j, q):
        r"""
        """
        assert i == 0 or i == 1
        assert j == 0 or j == 1
        assert i == 1 - j
        assert 0 <= p < self.length(i)
        assert 0 <= q < self.length(j)
        self._twin[i][p] = q

    def twin_list(self):
        r"""
        Returns the twin list of self.

        The twin list is the involution without fixed point associated to that
        permutation seen as two lines of symbols. As the domain is two lines,
        the position are 2-tuples `(i,j)` where `i` specifies the line and `j`
        the position in the line.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c','c b a')
            sage: p.twin_list()[0]
            [(1, 2), (1, 1), (1, 0)]
            sage: p.twin_list()[1]
            [(0, 2), (0, 1), (0, 0)]

        We may check that it is actually an involution without fixed point::

            sage: t = p.twin_list()
            sage: all(t[i][j] != (i,j) for i in range(2) for j in range(len(t[i])))
            True
            sage: all(t[t[i][j][0]][t[i][j][1]] == (i,j) for i in range(2) for j in range(len(t[i])))
            True
        """
        twin0 = [(1,i) for i in self._twin[0]]
        twin1 = [(0,i) for i in self._twin[1]]
        return [twin0,twin1]

    def _init_alphabet(self,a) :
        r"""
        Initializes the alphabet from intervals.

        INPUT:

        - ``a`` - the two intervals as lists

        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c d','d c a b')   #indirect doctest
            sage: p.alphabet() == Alphabet(['a', 'b', 'c', 'd'])
            True
            sage: p = iet.Permutation([0,1,2],[1,0,2],reduced=True)   #indirect doctest
            sage: p.alphabet() == Alphabet([0,1,2])
            True
        """
        from sage.combinat.words.alphabet import build_alphabet
        self._alphabet = build_alphabet(a[0])

    def _inversed_twin(self):
        r"""
        Inverses the twin of the permutation.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c','c a b',reduced=True)
            sage: p.left_right_inverse() # indirect doc test
            a b c
            b c a

        ::

            sage: p = iet.Permutation('a b c d','d a b c',reduced=True)
            sage: p.left_right_inverse() # indirect doctest
            a b c d
            b c d a
        """
        self._twin = [self._twin[1], self._twin[0]]

    def _reversed_twin(self):
        r"""
        Reverses the twin of the permutation.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c','c a b',reduced=True)
            sage: p.top_bottom_inverse() # indirect doctest
            a b c
            b c a
            sage: p = iet.Permutation('a b c d','d a b c',reduced=True)
            sage: p.top_bottom_inverse() # indircet doctest
            a b c d
            b c d a
        """
        tmp = [self._twin[0][:], self._twin[1][:]]

        n = self.length_top()
        for i in (0,1):
            for j in range(n):
                tmp[i][n- 1 - j] = n - 1 - self._twin[i][j]

        self._twin = tmp

    def _move(self, interval, position, interval_to, position_to):
        r"""
        Moves the element at (interval,position) to (interval, position_to)

        INPUT:

        - ``interval`` - 0 or 1

        - ``position`` - a position in interval

        - ``interval_to`` - 0 or 1

        - ``position_to`` - a position in interval

        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c d','d c b a')
            sage: p._move(0,0,0,2)
            sage: p
            b a c d
            d c b a
            sage: p._move(0,1,0,3)
            sage: p
            b c a d
            d c b a
            sage: p._move(0,2,0,4)
            sage: p
            b c d a
            d c b a
            sage: p._move(1,3,1,1)
            sage: p
            b c d a
            d a c b
        """
        assert interval == interval_to

        if position < position_to:
            if self._flips is not None:
                self._flips[interval].insert(
                        position_to-1,
                        self._flips[interval].pop(position))

            if self._labels is not None:
                self._labels[interval].insert(
                        position_to-1,
                        self._labels[interval].pop(position))

            elti = 1-interval
            eltp = self._twin[interval][position]
            k = self._twin[interval][position+1:position_to]

            # decrement the twin between position and position to
            for j,pos in enumerate(k):
                self._twin[elti][pos] = position+j

            # modify twin of the moved element
            self._twin[elti][eltp] = position_to-1

            # move
            self._twin[interval].insert(
                    position_to-1,
                    self._twin[interval].pop(position))

        elif position_to < position:
            if self._flips is not None:
                self._flips[interval].insert(
                        position_to,
                        self._flips[interval].pop(position))

            if self._labels is not None:
                self._labels[interval].insert(
                        position_to,
                        self._labels[interval].pop(position))

            elti = 1-interval
            eltp = self._twin[interval][position]
            k = self._twin[interval][position_to:position]

            # increment the twin between position and position to
            for j,pos in enumerate(k):
                self._twin[1-interval][pos] = position_to+j+1

            # modify twin of the moved element
            self._twin[elti][eltp] = position_to

            # move
            self._twin[interval].insert(
                    position_to,
                    self._twin[interval].pop(position))

    def has_rauzy_move(self, winner, side='right'):
        r"""
        Test if a Rauzy move can be performed on this permutation.

        EXAMPLES::

            sage: from surface_dynamics import *

        for labelled permutations::

            sage: p = iet.Permutation('a b c','a c b',reduced=False)
            sage: p.has_rauzy_move(0,'right')
            True
            sage: p.has_rauzy_move(0,'left')
            False
            sage: p.has_rauzy_move(1,'right')
            True
            sage: p.has_rauzy_move(1,'left')
            False

        for reduced permutations::

            sage: p = iet.Permutation('a b c','a c b',reduced=True)
            sage: p.has_rauzy_move(0,'right')
            True
            sage: p.has_rauzy_move(0,'left')
            False
            sage: p.has_rauzy_move(1,'right')
            True
            sage: p.has_rauzy_move(1,'left')
            False
        """
        side = side_conversion(side)
        winner = interval_conversion(winner)

        return self._twin[winner][side] % len(self) != side % len(self)

    def is_irreducible(self, return_decomposition=False) :
        r"""
        Test irreducibility.

        A permutation p = (p0,p1) is reducible if:
        set(p0[:i]) = set(p1[:i]) for an i < len(p0)

        OUTPUT:

        - a boolean

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: p.is_irreducible()
            True

            sage: p = iet.Permutation('a b c', 'b a c')
            sage: p.is_irreducible()
            False

            sage: p = iet.Permutation('a b c', 'c b a', flips=['a'])
            sage: p.is_irreducible()
            True

        """
        s0, s1 = 0, 0
        for i in range(len(self)-1) :
            s0 += i
            s1 += self._twin[0][i]
            if s0 == s1 :
                if return_decomposition :
                    return False, (self[0][:i+1], self[0][i+1:], self[1][:i+1], self[1][i+1:])
                return False
        if return_decomposition:
            return True, (self[0],[],self[1],[])
        return True


class PermutationLI(Permutation):
    def _init_twin(self, a):
        r"""
        Initializes the _twin attribute

        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('a a','b b',reduced=True)   #indirect doctest
            sage: p._twin
            [[(0, 1), (0, 0)], [(1, 1), (1, 0)]]
        """
        if a is None:
            self._twin = [[],[]]

        else:
            twin = [
                [(0,j) for j in range(len(a[0]))],
                [(1,j) for j in range(len(a[1]))]]

            for i in (0, 1):
                for j in range(len(twin[i])) :
                    if twin[i][j] == (i,j) :
                        if a[i][j] in a[i][j+1:] :
                            # two up or two down
                            j2 = (a[i][j+1:]).index(a[i][j]) + j + 1
                            twin[i][j] = (i,j2)
                            twin[i][j2] = (i,j)
                        else:
                            # one up, one down (here i=0)
                            j2 = a[1].index(a[i][j])
                            twin[0][j] = (1,j2)
                            twin[1][j2] = (0,j)

            self._twin = twin

    def _init_alphabet(self, intervals) :
        r"""
        Initialization procedure of the alphabet of self from intervals list

        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('a a','b b')   #indirect doctest
            sage: p.alphabet()
            {'a', 'b'}

            sage: p = iet.GeneralizedPermutation('b b','a a')  #indirect doctest
            sage: p.alphabet()
            {'b', 'a'}
        """
        tmp_alphabet = []
        for letter in intervals[0] + intervals[1] :
            if letter not in tmp_alphabet :
                tmp_alphabet.append(letter)

        from sage.combinat.words.alphabet import build_alphabet
        self._alphabet = build_alphabet(tmp_alphabet)

    def _inversed_twin(self):
        r"""
        Inverses the twin of the permutation.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('a a','b b',reduced=True)
            sage: p.left_right_inverse() #indirect doctest
            a a
            b b
            sage: p = iet.GeneralizedPermutation('a a b','b c c',reduced=True)
            sage: p.left_right_inverse() #indirect doctest
            a b b
            c c a
            sage: p = iet.GeneralizedPermutation('a a','b b c c',reduced=True)
            sage: p.left_right_inverse() #indirect doctest
            a a b b
            c c
        """
        self._twin = [self._twin[1], self._twin[0]]

        for interval in (0,1):
            for j in range(self.length(interval)):
                self._twin[interval][j] = (
                    1-self._twin[interval][j][0],
                    self._twin[interval][j][1])

    def _reversed_twin(self):
        r"""
        Reverses the twin of the permutation.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('a b b','c c a',reduced=True)
            sage: p.top_bottom_inverse() #indirect doctest
            a a b
            b c c
            sage: p = iet.GeneralizedPermutation('a a','b b c c',reduced=True)
            sage: p.top_bottom_inverse() #indirect doctest
            a a
            b b c c
        """
        tmp = [self._twin[0][:], self._twin[1][:]]

        n = self.length()

        for i in (0,1):
            for j in range(n[i]):
                interval, position = self._twin[i][j]
                tmp[i][n[i] - 1 - j] = (
                    interval,
                    n[interval] - 1 - position)

        self._twin = tmp

    def _move(self, interval, position, interval_to, position_to):
        r"""
        _move the element at (interval,position) to (interval_to, position_to)

        INPUT:

        - ``interval`` - 0 or 1

        - ``position`` - a position in the corresponding interval

        - ``interval_to`` - 0 or 1

        - ``position_to`` - a position in the corresponding interval

        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('a b c c','d b d a',flips='a')
            sage: p._move(0,0,0,2)
            sage: p
             b -a  c  c
             d  b  d -a
            sage: p._move(0,1,0,0)
            sage: p
            -a  b  c  c
             d  b  d -a
            sage: p._move(1,1,1,4)
            sage: p
            -a  b  c  c
             d  d -a  b
            sage: p._move(1,3,1,1)
            sage: p
            -a  b  c  c
             d  b  d -a
        """
        if interval != interval_to:
            if self._flips is not None:
                self._flips[interval_to].insert(
                        position_to,
                        self._flips[interval].pop(position))

            if self._labels is not None:
                self._labels[interval_to].insert(
                    position_to,
                    self._labels[interval].pop(position))

            elti,eltp = self._twin[interval][position] # the element to move
            k1 = self._twin[interval_to][position_to:] # interval to shift
            k2 = self._twin[interval][position+1:] # interval to unshift

            # increment the twin after the position_to
            for j,(tw,pos) in enumerate(k1):
                self._twin[tw][pos] = (interval_to, position_to+j+1)

            # decrement the twin after the position
            for j,(tw,pos) in enumerate(k2):
                self._twin[tw][pos] = (interval, position+j)

            # modify twin of the moved interval
            self._twin[elti][eltp] = (interval_to, position_to)

            # move
            self._twin[interval_to].insert(
                    position_to,
                    self._twin[interval].pop(position))

        else: # interval == interval_to (just one operation !)
            if position < position_to:
                if self._flips is not None:
                    self._flips[interval].insert(
                            position_to-1,
                            self._flips[interval].pop(position))

                if self._labels is not None:
                    self._labels[interval].insert(
                            position_to-1,
                            self._labels[interval].pop(position))

                elti, eltp = self._twin[interval][position]
                k = self._twin[interval][position+1:position_to]

                # decrement the twin between position and position to
                for j,(tw,pos) in enumerate(k):
                    self._twin[tw][pos] = (interval,position+j)

                # modify twin of the moved element
                self._twin[elti][eltp] = (interval_to, position_to-1)

                # move
                self._twin[interval].insert(
                        position_to-1,
                        self._twin[interval].pop(position))

            elif position_to < position:
                if self._flips is not None:
                    self._flips[interval].insert(
                            position_to,
                            self._flips[interval].pop(position))

                if self._labels is not None:
                    self._labels[interval].insert(
                            position_to,
                            self._labels[interval].pop(position))

                elti, eltp = self._twin[interval][position]
                k = self._twin[interval][position_to:position]

                # increment the twin between position and position to
                for j,(tw,pos) in enumerate(k):
                    self._twin[tw][pos] = (interval,position_to+j+1)

                # modify twin of the moved element
                self._twin[elti][eltp] = (interval_to,position_to)

                # move
                self._twin[interval].insert(
                        position_to,
                        self._twin[interval].pop(position))


    def _set_twin(self, i, p, j, q):
        assert i == 0 or i == 1
        assert j == 0 or j == 1
        assert 0 <= p < self.length(i)
        assert 0 <= q < self.length(j)
        self._twin[i][p] = (j,q)

    def twin(self, i, pos):
        r"""
        Return the twin of the letter in interval ``i`` at position ``pos``

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: p = iet.GeneralizedPermutation('a a b c', 'c e b e')
            sage: p.twin(0,0)
            (0, 1)
            sage: p.twin(0,1)
            (0, 0)

            sage: twin_top = [p.twin(0,i) for i in range(p.length_top())]
            sage: twin_bot = [p.twin(1,i) for i in range(p.length_bottom())]
            sage: p.twin_list() == [twin_top, twin_bot]
            True
        """
        return self._twin[i][pos]

    def twin_list(self):
        r"""
        Returns the twin list of self.

        The twin list is the involution without fixed point which defines it. As the
        domain is naturally split into two lines we use a 2-tuple (i,j) to
        specify the element at position j in line i.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('a a b','b c c')
            sage: p.twin_list()[0]
            [(0, 1), (0, 0), (1, 0)]
            sage: p.twin_list()[1]
            [(0, 2), (1, 2), (1, 1)]

        And we may check that it is actually an involution without fixed point::

            sage: t = p.twin_list()
            sage: all(t[i][j] != (i,j) for i in range(2) for j in range(len(t[i])))
            True
            sage: all(t[t[i][j][0]][t[i][j][1]] == (i,j) for i in range(2) for j in range(len(t[i])))
            True

        A slightly more complicated example::

            sage: q = iet.GeneralizedPermutation('a b c a','d e f e g c b g d f')
            sage: q.twin_list()[0]
            [(0, 3), (1, 6), (1, 5), (0, 0)]
            sage: q.twin_list()[1]
            [(1, 8), (1, 3), (1, 9), (1, 1), (1, 7), (0, 2), (0, 1), (1, 4), (1, 0), (1, 2)]

        ::

            sage: t = q.twin_list()
            sage: all(t[t[i][j][0]][t[i][j][1]] == (i,j) for i in range(2) for j in range(len(t[i])))
            True
        """
        return [self._twin[0][:],self._twin[1][:]]

    def is_irreducible(self, return_decomposition=False):
        r"""
        Test of reducibility

        A quadratic (or generalized) permutation is *reducible* if there exists
        a decomposition

        .. math::

           A1 u B1 | ... | B1 u A2

           A1 u B2 | ... | B2 u A2

        where no corners is empty, or exactly one corner is empty
        and it is on the left, or two and they are both on the
        right or on the left. The definition is due to [BoiLan09]_ where they prove
        that the property of being irreducible is stable under Rauzy induction.

        INPUT:

        -  ``return_decomposition`` - boolean (default: False) - if True, and
           the permutation is reducible, returns also the blocks A1 u B1, B1 u
           A2, A1 u B2 and B2 u A2 of a decomposition as above.

        OUTPUT:

        If return_decomposition is True, returns a 2-uple
        (test,decomposition) where test is the preceding test and
        decomposition is a 4-uple (A11,A12,A21,A22) where:

        A11 = A1 u B1
        A12 = B1 u A2
        A21 = A1 u B2
        A22 = B2 u A2

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: GP = iet.GeneralizedPermutation

            sage: GP('a a','b b').is_irreducible()
            False
            sage: GP('a a b','b c c').is_irreducible()
            True
            sage: GP('1 2 3 4 5 1','5 6 6 4 3 2').is_irreducible()
            True

        TESTS::

            sage: from surface_dynamics import *

        Test reducible permutations with no empty corner::

            sage: GP('1 4 1 3','4 2 3 2').is_irreducible(True)
            (False, (['1', '4'], ['1', '3'], ['4', '2'], ['3', '2']))

        Test reducible permutations with one left corner empty::

            sage: GP('1 2 2 3 1','4 4 3').is_irreducible(True)
            (False, (['1'], ['3', '1'], [], ['3']))
            sage: GP('4 4 3','1 2 2 3 1').is_irreducible(True)
            (False, ([], ['3'], ['1'], ['3', '1']))

        Test reducible permutations with two left corner empty::

            sage: GP('1 1 2 3','4 2 4 3').is_irreducible(True)
            (False, ([], ['3'], [], ['3']))

        Test reducible permutations with two right corner empty::

            sage: GP('1 2 2 3 3','1 4 4').is_irreducible(True)
            (False, (['1'], [], ['1'], []))
            sage: GP('1 2 2','1 3 3').is_irreducible(True)
            (False, (['1'], [], ['1'], []))
            sage: GP('1 2 3 3','2 1 4 4 5 5').is_irreducible(True)
            (False, (['1', '2'], [], ['2', '1'], []))

        A ``NotImplementedError`` is raised when there are flips::

            sage: p = iet.GeneralizedPermutation('a b c e b','d c d a e', flips='abcd', reduced=True)
            sage: p.is_irreducible()
            Traceback (most recent call last):
            ...
            NotImplementedError: irreducibility test not implemented for generalized permutations with flips
            sage: p = iet.GeneralizedPermutation('a b c e b','d c d a e', flips='abcd', reduced=False)
            sage: p.is_irreducible()
            Traceback (most recent call last):
            ...
            NotImplementedError: irreducibility test not implemented for generalized permutations with flips
        """
        if self._flips is not None:
            raise NotImplementedError('irreducibility test not implemented for generalized permutations with flips')

        l0 = self.length_top()
        l1 = self.length_bottom()
        s0,s1 = self.list()

        # testing two corners empty on the right (i12 = i22 = 0)
        A11, A21, A12, A22 = [],[],[],[]

        for i11 in range(1, l0):
            if s0[i11-1] in A11:
                break
            A11 = s0[:i11]

            for i21 in range(1, l1):
                if s1[i21-1] in A21:
                    break
                A21 = s1[:i21]

                if sorted(A11)  == sorted(A21):
                    if return_decomposition:
                        return False,(A11,A12,A21,A22)
                    return False
            A21 = []

        # testing no corner empty but one or two on the left
        t11 = t21 = False
        A11, A12, A21, A22 = [], [], [], []
        for i11 in range(0, l0):
            if i11 > 0 and s0[i11-1] in A11:
                break
            A11 = s0[:i11]

            for i21 in range(0, l1) :
                if i21 > 0 and s1[i21-1] in A21:
                    break
                A21 = s1[:i21]

                for i12 in range(l0 - 1, i11 - 1, -1) :
                    if s0[i12] in A12 or s0[i12] in A21:
                        break
                    A12 = s0[i12:]

                    for i22 in range(l1 - 1, i21 - 1, -1) :
                        if s1[i22] in A22 or s1[i22] in A11:
                            break
                        A22 = s1[i22:]

                        if sorted(A11 + A22) == sorted(A12 + A21) :
                            if return_decomposition :
                                return False, (A11,A12,A21,A22)
                            return False
                    A22 = []
                A12 = []
            A21 = []


        if return_decomposition:
            return True, ()
        return True


    def to_cylindric(self):
        r"""
        Return a cylindric permutation in the same extended Rauzy class

        A generalized permutation is *cylindric* if the first letter in the top
        interval is the same as the last letter in the bottom interval or if the
        laster letter of the top interval is the same as the fist letter of the
        bottom interval.

        EXAMPLES::

            sage: from surface_dynamics import iet

            sage: p = iet.GeneralizedPermutation('a b d a c','c e b e d')
            sage: p.is_irreducible()
            True
            sage: p.to_cylindric().is_cylindric()
            True

            sage: p = iet.GeneralizedPermutation([0,1,2,1,2,3,4,5,6], [6,0,5,4,7,3,7], reduced=True)
            sage: p.to_cylindric()
            0 1 2 1 2 3 4 5 6
            6 0 5 4 7 3 7

        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation([[0,1,1],[2,2,0]], reduced=True)
            sage: p.to_cylindric()
            0 1 1
            2 2 0

        ALGORITHM:

        The algorithm is naive. It computes the extended Rauzy class until it
        finds a cylindric permutation.
        """
        if self.is_cylindric():
            return self

        wait = [self]
        rauzy_class = set([self])
        while wait:
            q = wait.pop()
            if q.has_rauzy_move('t'): # top rauzy move
                qq = q.rauzy_move('t')
                if qq not in rauzy_class:
                    if qq._twin[1][-1] == (0,0) or qq._twin[0][-1] == (1,0):
                        return qq
                    wait.append(qq)
                    rauzy_class.add(qq)
            if q.has_rauzy_move('b'): # bot rauzy move
                qq = q.rauzy_move('b')
                if qq not in rauzy_class:
                    if qq._twin[1][-1] == (0,0) or qq._twin[0][-1] == (1,0):
                        return qq
                    wait.append(qq)
                    rauzy_class.add(qq)
            qq = q.symmetric() # symmetric
            if qq not in rauzy_class:
                if qq._twin[1][-1] == (0,0) or qq._twin[0][-1] == (1,0):
                    return qq
                wait.append(qq)
                rauzy_class.add(qq)

        raise RuntimeError("no cylindric permutation in the extended Rauzy class")

    def is_cylindric(self):
        r"""
        Test if the permutation is cylindric

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: q = iet.GeneralizedPermutation('a b b','c c a')
            sage: q.is_cylindric()
            True
            sage: q = iet.GeneralizedPermutation('a a b b','c c')
            sage: q.is_cylindric()
            False
        """
        return self._twin[0][-1] == (1,0) or self._twin[1][-1] == (0,0)

    def profile(self):
        r"""
        Returns the ``profile`` of self.

        The *profile* of a generalized permutation is the list `(d_1, \ldots,
        d_k)` where `(d_1 \pi, \ldots, d_k \pi)` is the list of angles of any
        suspension of that generalized permutation.

        See also :meth:`marked_profile`.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p1 = iet.GeneralizedPermutation('a a b','b c c')
            sage: p1.profile()
            [1, 1, 1, 1]
            sage: all(p.profile() == [1, 1, 1, 1] for p in p1.rauzy_diagram())
            True

            sage: p2 = iet.GeneralizedPermutation('0 1 2 1 3','4 3 4 2 0')
            sage: p2.profile()
            [4, 4]
            sage: all(p.profile() == [4,4] for p in p2.rauzy_diagram())
            True

            sage: p3 = iet.GeneralizedPermutation('0 1 2 3 3','2 1 4 4 0')
            sage: p3.profile()
            [3, 3, 1, 1]
            sage: all(p.profile() == [3, 3, 1, 1] for p in p3.rauzy_diagram())
            True
        """
        from sage.combinat.partition import Partition
        s = self.interval_diagram(sign=False,glue_ends=True)
        return Partition(sorted((len(x) for x in s),reverse=True))

    def genus(self):
        r"""
        Returns the genus of any suspension of self.

        The genus `g` can be deduced from the profile (see :meth:`profile`)
        `p=(p_1,\ldots,p_k)` of self by the formula:
        `4g-4 = \sum_{i=1}^k (p_i - 2)`.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: iet.GeneralizedPermutation('a a b','b c c').genus()
            0
            sage: iet.GeneralizedPermutation((0,1,2,1,3),(4,3,4,2,0)).genus()
            2
        """
        p = self.profile()
        return Integer((sum(p)-2*len(p)) // 4 + 1)

    def marking(self):
        r"""
        Return the marking induced by the two sides of the interval

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('0 1 2 3 4 3 5 6 7','1 6 8 4 2 7 5 8 0')
            sage: p.marking()
            8|7

            sage: p = iet.GeneralizedPermutation('0 1 2 3 4 3 5 6 7','1 6 8 4 2 7 8 0 5')
            sage: p.marking()
            8o8
        """
        return self.marked_profile().marking()

    def marked_profile(self):
        r"""
        Returns the marked profile of self.

        The *marked profile* of a generalized permutation is an integer
        partition and some additional data associated to the angles of conical
        singularities in the suspension. The partition, called the
        *profile*, is the list of angles divided by `2\pi` (see
        :meth:`profile`). The additional is called the *marking* and may be of
        two different types.

        If the left endpoint and the right endpoint of the interval associated
        to the permutation coincides, then the marking is of *type 1* and the
        additional data consists of a couple `(m,a)` such that `m` is the
        angle of the conical singularity and `a` is the angle between the
        outgoing separatrix associated to the left endpoint and the incoming
        separatrix associated to the right endpoint. A marking of type one is
        denoted `m | a`.

        If the left endpoint and the right endpoint are two different conical
        singularities in the suspension, then the marking is of *type 2* and the
        data consists in a couple `(m_l,m_r)` where `m_l` (resp. `m_r`) is
        the conical angle of the singularity at the left endpoint (resp. right
        endpoint). A marking of type two is denoted `m_l \circ m_r`

        EXAMPLES::

            sage: from surface_dynamics import *

        All possible markings for the profile [1, 1, 1, 1]::

            sage: p = iet.GeneralizedPermutation('a a b','b c c')
            sage: p.marked_profile()
            1o1 [1, 1, 1, 1]
            sage: p = iet.GeneralizedPermutation('a a','b b c c')
            sage: p.marked_profile()
            1|0 [1, 1, 1, 1]

        All possible markings for the profile [4, 4]::

            sage: p = iet.GeneralizedPermutation('0 1 2 1 3','3 4 0 4 2')
            sage: p.marked_profile()
            4o4 [4, 4]

            sage: p = iet.GeneralizedPermutation('0 1 2 1 3','4 3 2 0 4')
            sage: p.marked_profile()
            4|0 [4, 4]

            sage: p = iet.GeneralizedPermutation('0 1 0 2 3 2','4 3 4 1')
            sage: p.marked_profile()
            4|1 [4, 4]

            sage: p = iet.GeneralizedPermutation('0 1 2 3 2','4 3 4 1 0')
            sage: p.marked_profile()
            4|2 [4, 4]

            sage: p = iet.GeneralizedPermutation('0 1 0 1','2 3 2 4 3 4')
            sage: p.marked_profile()
            4|3 [4, 4]
        """
        if self._flips:
            raise ValueError('not available on permutations with flips')

        from .marked_partition import MarkedPartition

        if len(self) == 1:
            return MarkedPartition([],2,(0,0))

        g = self.interval_diagram(glue_ends=True,sign=True)
        signs = self._canonical_signs()[1]
        p = sorted(map(lambda x: len(x), g),reverse=True)

        left1 = ((self[1][0], -signs[1][0]), (self[0][0], signs[0][0]))
        left2 = (left1[1], left1[0])
        right1 = ((self[0][-1], -signs[0][-1]), (self[1][-1], signs[1][-1]))
        right2 = (right1[1], right1[0])
        if len(set(left1+right1)) == 3:
            if left1[0] == right1[0]:
                lr1 = (left1[1], left1[0],right1[1])
                lr2 = (right1[1], left1[0], left1[1])
            elif left1[0] == right1[1]:
                lr1 = (left1[1], left1[0], right1[0])
                lr2 = (right1[0], left1[0], left1[1])
            elif left1[1] == right1[0]:
                lr1 = (left1[0], left1[1], right1[1])
                lr2 = (right1[1], left1[1], left1[0])
            elif left1[1] == right1[1]:
                lr1 = (left1[0], left1[1], right1[0])
                lr2 = (right1[0], left1[1], left1[0])
            for c in g:
                if lr1 in c or lr2 in c:
                    break
            return MarkedPartition(p, 1, (len(c), 0))
        else:
            c_left = c_right = None
            for c in g:
                if left1 in c or left2 in c: c_left = c
                if right1 in c or right2 in c: c_right = c

        if c_left == c_right:
            mm = len(c_left)
            if right1 in c_right:
                r = c_right.index(right1)
            else:
                r = c_right.index(right2)
            if left1 in c_left:
                l = c_left.index(left1)
            else:
                l = c_left.index(left2)
            a = ((r-l)%mm)
            return MarkedPartition(p, 1, (mm, a))

        else:
            m_l = len(c_left)
            m_r = len(c_right)
            return MarkedPartition(p, 2, (m_l,m_r))

    def erase_marked_points(self):
        r"""
        Return a permutation without marked points.

        This method is not implemented for generalized permutations.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('a a b','b c c')
            sage: p.stratum()
            Q_0(-1^4)
            sage: p.erase_marked_points()
            a a b
            b c c
            sage: p = iet.GeneralizedPermutation('a d d a b','b c c')
            sage: p.stratum()
            Q_0(0, -1^4)
            sage: p.erase_marked_points()
            Traceback (most recent call last):
            ...
            NotImplementedError: not yet implemented! Do it!
        """
        if 0 in self.stratum().signature():
            raise NotImplementedError("not yet implemented! Do it!")
        else:
            return self

    def is_hyperelliptic(self, verbose=False):
        r"""
        Test if this permutation is in an hyperelliptic connected component.

        EXAMPLES::

            sage: from surface_dynamics import *

        An example of hyperelliptic permutation::

            sage: p = iet.GeneralizedPermutation([0,1,2,0,6,5,3,1,2,3],[4,5,6,4])
            sage: p.is_hyperelliptic()
            True

        Check for the correspondence::

            sage: q = Stratum([6,6], k=2)
            sage: c_hyp, c_reg, c_irr = q.components()

            sage: p_hyp = c_hyp.permutation_representative()
            sage: p_hyp
            0 1 2 3 4 1 5 6 7
            7 6 5 8 4 3 2 8 0
            sage: p_hyp.is_hyperelliptic()
            True

            sage: p_reg = c_reg.permutation_representative()
            sage: p_reg
            0 1 2 3 4 5 2 6 7 5
            1 4 6 8 7 8 3 0
            sage: p_reg.is_hyperelliptic()
            False

            sage: p_irr = c_irr.permutation_representative()
            sage: p_irr
            0 1 2 3 4 3 5 6 7
            1 6 8 4 2 7 5 8 0
            sage: p_irr.is_hyperelliptic()
            False

            sage: q = Stratum([3,3,2], k=2)
            sage: c_hyp, c_non_hyp = q.components()
            sage: p_hyp = c_hyp.permutation_representative()
            sage: p_hyp.is_hyperelliptic()
            True
            sage: p_non_hyp = c_non_hyp.permutation_representative()
            sage: p_non_hyp.is_hyperelliptic()
            False
            sage: q = Stratum([5,5,2], k=2)
            sage: c_hyp, c_non_hyp = q.components()
            sage: p_hyp = c_hyp.permutation_representative()
            sage: p_hyp.is_hyperelliptic()
            True
            sage: p_non_hyp = c_non_hyp.permutation_representative()
            sage: p_non_hyp.is_hyperelliptic()
            False
            sage: q = Stratum([3,3,1,1], k=2)
            sage: c_hyp, c_non_hyp = q.components()
            sage: p_hyp = c_hyp.permutation_representative()
            sage: p_hyp.is_hyperelliptic()
            True
            sage: p_non_hyp = c_non_hyp.permutation_representative()
            sage: p_non_hyp.is_hyperelliptic()
            False
        """
        p = self.erase_marked_points()
        s = p.stratum()
        zeros = s.signature()

        if not s.has_hyperelliptic_component():
            return False

        q = p.to_cylindric()

        if q[0][0] == q[1][-1]:
            l0 = []
            q0 = q[0][1:]
            q1 = q[1][:-1]
            for i,j in q._twin[0][1:]:
                if i == 0: l0.append((0,j-1))
                else: l0.append((1,j))
            l1 = []
            for i,j in q._twin[1][:-1]:
                if i == 0: l1.append((0,j-1))
                else: l1.append((1,j))
        else:
            l0 = []
            q0 = q[0][:-1]
            q1 = q[1][1:]
            for i,j in q._twin[0][:-1]:
                if i == 1: l0.append((1,j-1))
                else: l0.append((0,j))
            l1 = []
            for i,j in q._twin[1][1:]:
                if i ==1: l1.append((1,j-1))
                else: l1.append((0,j))

        if verbose:
            print("found Jenkins-Strebel")
            print(q)
            print(l0)
            print(l1)

        if any(x[0] == 1 for x in l0):
            if verbose: print("potential form 1")
            i0 = []
            i1 = []
            for i in range(len(l0)):
                if l0[i][0] == 0:
                    i0.append(i)
            for i in range(len(l1)):
                if l1[i][0] == 1:
                    i1.append(i)
            if len(i0) != 2 or len(i1) != 2:
                if verbose: print("no repetition twice in intervals")
                return False

            q0_0 = q0[i0[0]+1:i0[1]]
            q0_1 = q0[i0[1]+1:] + q0[:i0[0]]
            q0_0.reverse()
            q0_1.reverse()

            q1_0 = q1[i1[0]+1:i1[1]]
            q1_1 = q1[i1[1]+1:] + q1[:i1[0]]

            if verbose:
                print(q0_0, q0_1)
                print(q1_0, q1_1)

            return (q0_0 == q1_0 and q0_1 == q1_1) or (q0_0 == q1_1 and q0_1 == q1_0)

        else:
            if verbose: print("potential form 2")
            if any(i==1 for i,_ in l0) or any(i==0 for i,_ in l1):
                return False
            j = len(l0) // 2
            for i in range(j):
                if l0[i][1] != j+i:
                    return False

            j = len(l1) // 2
            for i in range(j):
                if l1[i][1] != j+i:
                    return False

            return True

    def stratum_component(self):
        r"""
        Return the connected component of stratum in which self belongs to.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('a b b','c c a')
            sage: p.stratum_component()
            Q_0(-1^4)^c

        Test the exceptional strata in genus 3::

            sage: Q = Stratum([9,-1], k=2)
            sage: p = Q.regular_component().permutation_representative()
            sage: p.stratum_component()
            Q_3(9, -1)^reg
            sage: p = Q.irregular_component().permutation_representative()
            sage: p.stratum_component()
            Q_3(9, -1)^irr

            sage: Q = Stratum([6,3,-1], k=2)
            sage: p = Q.regular_component().permutation_representative()
            sage: p.stratum_component()
            Q_3(6, 3, -1)^reg
            sage: p = Q.irregular_component().permutation_representative()
            sage: p.stratum_component()
            Q_3(6, 3, -1)^irr

            sage: Q = Stratum([3,3,3,-1], k=2)
            sage: p = Q.regular_component().permutation_representative()
            sage: p.stratum_component()
            Q_3(3^3, -1)^reg
            sage: p = Q.irregular_component().permutation_representative()
            sage: p.stratum_component()
            Q_3(3^3, -1)^irr

        Test the exceptional strata in genus 4::

            sage: Q = Stratum([12], k=2)
            sage: p = Q.regular_component().permutation_representative()
            sage: p.stratum_component()
            Q_4(12)^reg
            sage: p = Q.irregular_component().permutation_representative()
            sage: p.stratum_component()
            Q_4(12)^irr

            sage: Q = Stratum([9,3], k=2)
            sage: p = Q.regular_component().permutation_representative()
            sage: p.stratum_component()
            Q_4(9, 3)^reg
            sage: p = Q.irregular_component().permutation_representative()
            sage: p.stratum_component()
            Q_4(9, 3)^irr

            sage: Q = Stratum([6,6], k=2)
            sage: p = Q.hyperelliptic_component().permutation_representative()
            sage: p.stratum_component()
            Q_4(6^2)^hyp
            sage: p = Q.regular_component().permutation_representative()
            sage: p.stratum_component()
            Q_4(6^2)^reg
            sage: p = Q.irregular_component().permutation_representative()
            sage: p.stratum_component()
            Q_4(6^2)^irr

            sage: Q = Stratum([6,3,3], k=2)
            sage: p = Q.regular_component().permutation_representative()
            sage: p.stratum_component()
            Q_4(6, 3^2)^reg
            sage: p = Q.irregular_component().permutation_representative()
            sage: p.stratum_component()
            Q_4(6, 3^2)^irr

            sage: Q = Stratum([3,3,3,3], k=2)
            sage: p = Q.hyperelliptic_component().permutation_representative()
            sage: p.stratum_component()
            Q_4(3^4)^hyp
            sage: p = Q.regular_component().permutation_representative()
            sage: p.stratum_component()
            Q_4(3^4)^reg
            sage: p = Q.irregular_component().permutation_representative()
            sage: p.stratum_component()
            Q_4(3^4)^irr
        """
        stratum = self.stratum()
        cc = stratum.components()

        if len(cc) == 1: # connected
            return cc[0]

        elif stratum.has_hyperelliptic_component(): # hyp / nonhyp
            if self.is_hyperelliptic():
                return stratum.hyperelliptic_component()
            elif len(cc) == 2:
                return stratum.non_hyperelliptic_component()

        # reg / irr
        from surface_dynamics.databases.flat_surfaces import IrregularComponentTwins
        D = IrregularComponentTwins()
        if len(D.list_strata()) != 8:
            raise NotImplementedError("database of irregular twins not available")

        p = self.erase_marked_points().to_cylindric()
        if cylindric_canonical(p) in D.get(stratum):
            return stratum.irregular_component()
        return stratum.regular_component()

    def has_rauzy_move(self, winner, side='right'):
        r"""
        Test of Rauzy movability (with an eventual specified choice of winner)

        A quadratic (or generalized) permutation is rauzy_movable type
        depending on the possible length of the last interval. It's
        dependent of the length equation.

        INPUT:

        - ``winner`` - the integer 'top' or 'bottom'

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('a a','b b')
            sage: p.has_rauzy_move('top','right')
            False
            sage: p.has_rauzy_move('top','left')
            False
            sage: p.has_rauzy_move('bottom','right')
            False
            sage: p.has_rauzy_move('bottom','left')
            False

        ::

            sage: p = iet.GeneralizedPermutation('a a b','b c c')
            sage: p.has_rauzy_move('top','right')
            True
            sage: p.has_rauzy_move('bottom','right')
            True
            sage: p.has_rauzy_move('top','left')
            True
            sage: p.has_rauzy_move('bottom','left')
            True

        ::

            sage: p = iet.GeneralizedPermutation('a a','b b c c')
            sage: p.has_rauzy_move('top','right')
            True
            sage: p.has_rauzy_move('bottom','right')
            False
            sage: p.has_rauzy_move('top','left')
            True
            sage: p.has_rauzy_move('bottom','left')
            False

        ::

            sage: p = iet.GeneralizedPermutation('a a b b','c c')
            sage: p.has_rauzy_move('top','right')
            False
            sage: p.has_rauzy_move('bottom','right')
            True
            sage: p.has_rauzy_move('top','left')
            False
            sage: p.has_rauzy_move('bottom','left')
            True
        """
        winner = interval_conversion(winner)
        side = side_conversion(side)

        loser = 1 - winner

        if side == -1:
            # the same letter at the right-end (False)
            if self._twin[0][-1] == (1, self.length_bottom()-1):
                return False

            # winner or loser letter is repeated on the other interval (True)
            if self._twin[0][-1][0] == 1: return True
            if self._twin[1][-1][0] == 0: return True

            # the loser letter is the only letter repeated in
            # the loser interval (False)
            if [i for i,_ in self._twin[loser]].count(loser) == 2:
                return False

            return True

        elif side == 0:
            # the same letter at the left-end (False)
            if (self._twin[0][0] == (1,0)):
                return False

            # winner or loser repeated on the other interval (True)
            if self._twin[0][0][0] == 1: return True
            if self._twin[1][0][0] == 0: return True

            # the loser letter is the only letter repeated in
            # the loser interval (False)
            if [i for i,_ in self._twin[loser]].count(loser) == 2:
                return False

            return True

    def orientation_cover(self):
        r"""
        Return the orientation cover of this permutation.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: c = p.orientation_cover()
            sage: c
            Covering of degree 2 of the permutation:
            a a b
            b c c
            sage: c.stratum()
            H_1(0^4)

            sage: C = Stratum([3,2,2,1], k=2).unique_component()
            sage: p = C.permutation_representative()
            sage: c = p.orientation_cover()
            sage: c.stratum()
            H_6(4, 2, 1^4)
        """
        rank = self.alphabet().rank
        p0 = set(map(rank, self[0]))
        p1 = set(map(rank, self[1]))
        inv_letters = p0.symmetric_difference(p1)
        permut_cover = [[1,0] if i in inv_letters else [0,1] for i in range(len(self))]
        return self.cover(permut_cover, as_tuple=True)

class OrientablePermutationIET(PermutationIET):
    """
    Template for permutation of Interval Exchange Transformation.

    .. warning::

        Internal class! Do not use directly!

    AUTHOR:

    - Vincent Delecroix (2008-12-20): initial version

    """
    def is_identity(self):
        r"""
        Returns True if self is the identity.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: iet.Permutation("a b","a b",reduced=False).is_identity()
            True
            sage: iet.Permutation("a b","a b",reduced=True).is_identity()
            True
            sage: iet.Permutation("a b","b a",reduced=False).is_identity()
            False
            sage: iet.Permutation("a b","b a",reduced=True).is_identity()
            False
        """
        return all(self._twin[0][i] == i for i in range(len(self)))

    #TODO: change the name
    def decompose(self):
        r"""
        Returns the decomposition as a concatenation of irreducible permutations.

        OUTPUT:

        a list of permutations

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c','c b a').decompose()[0]
            sage: p
            a b c
            c b a

        ::

            sage: p1,p2,p3 = iet.Permutation('a b c d e','b a c e d').decompose()
            sage: p1
            a b
            b a
            sage: p2
            c
            c
            sage: p3
            d e
            e d
        """
        l = []
        test, t = self.is_irreducible(return_decomposition=True)
        l.append(self.__class__((t[0],t[2])))

        while not test:
            q = self.__class__((t[1],t[3]))
            test, t = q.is_irreducible(return_decomposition=True)
            l.append(self.__class__((t[0],t[2])))

        return l

    def intersection_matrix(self, ring=None):
        r"""
        Returns the intersection matrix.

        This `d*d` antisymmetric matrix is given by the rule :

        .. math::

            m_{ij} = \begin{cases}
                1 & \text{$i < j$ and $\pi(i) > \pi(j)$} \\
                -1 & \text{$i > j$ and $\pi(i) < \pi(j)$} \\
                0 & \text{else}
                \end{cases}

        OUTPUT:

        - a matrix

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c d','d c b a')
            sage: p.intersection_matrix()
            [ 0  1  1  1]
            [-1  0  1  1]
            [-1 -1  0  1]
            [-1 -1 -1  0]

        ::

            sage: p = iet.Permutation('1 2 3 4 5','5 3 2 4 1')
            sage: p.intersection_matrix()
            [ 0  1  1  1  1]
            [-1  0  1  0  1]
            [-1 -1  0  0  1]
            [-1  0  0  0  1]
            [-1 -1 -1 -1  0]

        ::

            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: R = p.rauzy_diagram()
            sage: g = R.path(p, *'tbt')
            sage: m = g.matrix()
            sage: q = g.end()
            sage: q.intersection_matrix() == m.transpose() * p.intersection_matrix() * m
            True
        """
        if ring is None:
            ring = ZZ
        n = self.length_top()
        m = matrix(ring,n)
        # NOTE: because of the extended Rauzy inductions, the labels are just a subset of
        # {0, 1, ...} not necessarily the first n integers.
        use_labels = self._labels is not None and \
                     min(self._labels[0]) == 0 and \
                     max(self._labels[0]) == n - 1
        for i in range(n):
            ii = self._labels[0][i] if use_labels else i
            for j in range(i,n):
                jj = self._labels[0][j] if use_labels else j
                if self._twin[0][i] > self._twin[0][j]:
                    m[ii,jj] = 1
                    m[jj,ii] = -1
        return m

    def attached_out_degree(self):
        r"""
        Returns the degree of the singularity at the left of the interval.

        OUTPUT:

        - a positive integer


        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p1 = iet.Permutation('a b c d e f g','d c g f e b a')
            sage: p2 = iet.Permutation('a b c d e f g','e d c g f b a')
            sage: p1.attached_out_degree()
            3
            sage: p2.attached_out_degree()
            1
        """
        left_corner = (self[1][0], +1)
        for s in self.interval_diagram(glue_ends=False,sign=True):
            if left_corner in s:
                return len(s)//2 - 1

    def attached_in_degree(self):
        r"""
        Returns the degree of the singularity at the right of the interval.

        OUTPUT:

        - a positive integer


        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p1 = iet.Permutation('a b c d e f g','d c g f e b a')
            sage: p2 = iet.Permutation('a b c d e f g','e d c g f b a')
            sage: p1.attached_in_degree()
            1
            sage: p2.attached_in_degree()
            3
        """
        right_corner = (self[1][-1], -1)
        for s in self.interval_diagram(glue_ends=False,sign=True):
            if right_corner in s:
                return len(s)//2 - 1

    def profile(self):
        r"""
        Returns the profile of the permutation

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: iet.Permutation('a b c d','d c b a').profile()
            [3]
            sage: iet.Permutation('a b c d e','e d c b a').profile()
            [2, 2]
        """
        from sage.combinat.partition import Partition
        s = self.interval_diagram(glue_ends=True,sign=False)
        return Partition(sorted((len(x)/2 for x in s),reverse=True))

    def marking(self):
        r"""
        Return the marking induced by the two sides of the interval

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c d e f','f a e b d c')
            sage: p.marking()
            5|0
            sage: p = iet.Permutation('0 1 2 3 4 5 6','3 2 4 6 5 1 0')
            sage: p.marking()
            3o3
        """
        return self.marked_profile().marking()

    def marked_profile(self):
        r"""
        Returns the marked profile of the permutation

        The marked profile of a permutation corresponds to the integer partition
        associated to the angles of conical singularities in the suspension
        together with a data associated to the endpoint called marking.

        If the left endpoint and the right endpoint of the interval associated
        to the permutation, then the marking is of type one and consists in a
        couple ``(m,a)`` such that ``m`` is the angle of the conical singularity
        and ``a`` is the angle between the outgoing separatrix associated to the
        left endpoint and the incoming separatrix associated to the right
        endpoint. A marking of type one is denoted ``(m|a)``.

        If the left endpoint and the right endpoint are two different conical
        singularities in the suspension the marking is of type two and
        consists in a couple ``(m_l,m_r)`` where ``m_l`` (resp. ``m_r``) is the
        conical angle of the singularity at the left endpoint (resp. right
        endpoint). A marking of type two is denoted ``m_l o m_r``

        EXAMPLES::

            sage: from surface_dynamics import *

        The irreducible permutation on 1 interval has marked profile of type 2
        with data `(0,0)`::

            sage: p = iet.Permutation('a','a')
            sage: p.marked_profile()
            0o0 []

        Permutations in H(3,1) with all possible profiles::

            sage: p = iet.Permutation('a b c d e f g','b g a c f e d')
            sage: p.interval_diagram()
            [['a', ('b', 'a'), ('g', 'd'), 'e', 'f', 'g', 'b', 'c'], ['c', 'd', 'e', 'f']]
            sage: p.marked_profile()
            4|0 [4, 2]

            sage: p = iet.Permutation('a b c d e f g','c a g d f b e')
            sage: p.interval_diagram()
            [['a', 'b', 'f', 'g'], ['c', 'd', ('g', 'e'), 'f', 'd', 'e', 'b', ('c', 'a')]]
            sage: p.marked_profile()
            4|1 [4, 2]

            sage: p = iet.Permutation('a b c d e f g','e b d g c a f')
            sage: p.interval_diagram()
            [['a', 'b', 'e', 'f'], ['c', 'd', 'b', 'c', ('g', 'f'), 'g', 'd', ('e', 'a')]]
            sage: p.marked_profile()
            4|2 [4, 2]

            sage: p = iet.Permutation('a b c d e f g', 'e c g b a f d')
            sage: p.interval_diagram()
            [['a', 'b', ('g', 'd'), ('e', 'a'), 'b', 'c', 'e', 'f'], ['c', 'd', 'f', 'g']]
            sage: p.marked_profile()
            4|3 [4, 2]

            sage: p = iet.Permutation('a b c d e f g', 'f d c a g e b')
            sage: p.interval_diagram()
            [['a', 'b', 'e', ('f', 'a'), 'c', 'd', 'f', 'g'], ['c', 'd', 'e', ('g', 'b')]]
            sage: p.marked_profile()
            4o2 [4, 2]
        """
        from .marked_partition import MarkedPartition

        if len(self) == 1:
            return MarkedPartition([],2,(0,0))

        g = self.interval_diagram(glue_ends=True,sign=True)
        p = sorted(map(lambda x: len(x)//2, g),reverse=True)
        left = ((self[1][0], +1), (self[0][0], +1))
        right = ((self[0][-1], -1), (self[1][-1], -1))
        for c in g:
            if left in c: c_left = c
            if right in c: c_right = c

        if c_left == c_right:
            mm = len(c_left)
            a = ((c_right.index(right)-c_left.index(left)-1) %mm) // 2
            return MarkedPartition(p, 1, (mm//2, a))

        else:
            m_l = len(c_left) // 2
            m_r = len(c_right) // 2
            return MarkedPartition(p, 2, (m_l,m_r))

    def stratum(self):
        r"""
        Returns the strata in which any suspension of this permutation lives.

        OUTPUT:

        - a stratum of Abelian differentials

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: p.stratum()
            H_1(0^2)

            sage: p = iet.Permutation('a b c d', 'd a b c')
            sage: p.stratum()
            H_1(0^3)

            sage: p = iet.Permutation(list(range(9)), [8,5,2,7,4,1,6,3,0])
            sage: p.stratum()
            H_3(1^4)

            sage: a = 'a b c d e f g h i j'
            sage: b3 = 'd c g f e j i h b a'
            sage: b2 = 'd c e g f j i h b a'
            sage: b1 = 'e d c g f h j i b a'
            sage: p3 = iet.Permutation(a, b3)
            sage: p3.stratum()
            H_4(3, 2, 1)
            sage: p2 = iet.Permutation(a, b2)
            sage: p2.stratum()
            H_4(3, 2, 1)
            sage: p1 = iet.Permutation(a, b1)
            sage: p1.stratum()
            H_4(3, 2, 1)

        AUTHORS:

        - Vincent Delecroix (2008-12-20)
        """
        from surface_dynamics.flat_surfaces.abelian_strata import Stratum

        if not self.is_irreducible():
            return list(map(lambda x: x.stratum(), self.decompose()))

        if len(self) == 1:
            return Stratum([], k=1)

        singularities = [x - 1 for x in self.profile()]

        return Stratum(singularities, k=1)

    def genus(self) :
        r"""
        Returns the genus corresponding to any suspension of self.

        The genus can be deduced from the profile (see :meth:`profile`)
        `p = (p_1,\ldots,p_k)` of self by the formula:
        `2g-2 = \sum_{i=1}^k (p_i-1)`.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: p.genus()
            1

            sage: p = iet.Permutation('a b c d','d c b a')
            sage: p.genus()
            2
        """
        p = self.profile()
        return Integer((sum(p)-len(p))//2+1)

    def arf_invariant(self):
        r"""
        Returns the Arf invariant of the permutation.

        To a permutation `\pi` is associated a quadratic form on the field with
        2 elements. The *Arf invariant* is the total invariant of linear
        equivalence class of quadratic form of given rank.

        Let `V` be a vector space on the field with two elements `\FF_2`.  `V`
        there are two equivalence classes of non degenerate quadratic forms.  A
        complete invariant for quadratic forms is the *Arf invariant*.

        For non zero degenerate quadratic forms there are three equivalence
        classes. If `B` denotes the bilinear form associated to `q` then the
        three classes are as follows

        - the restriction of `q` to `ker(B)` is non zero

        - the restriction of `q` to `ker(B)` is zero and the spin parity of `q`
          on the quotient `V/ker(B)` is 0

        - the restriction of `q` to `ker(B)` is zero and the spin parity of `q`
          on the quotient `V/ker(B)` is 1

        The function returns respectively `None`, `0` or `1` depending on the
        three alternatives above.

        EXAMPLES::

            sage: from surface_dynamics import *

        Permutations from the odd and even component of H(2,2,2)::

            sage: a = list(range(10))
            sage: b1 = [3,2,4,6,5,7,9,8,1,0]
            sage: b0 = [6,5,4,3,2,7,9,8,1,0]
            sage: p1 = iet.Permutation(a,b1)
            sage: p1.arf_invariant()
            1
            sage: p0 = iet.Permutation(a,b0)
            sage: p0.arf_invariant()
            0

        Permutations from the odd and even component of H(4,4)::

            sage: a = list(range(11))
            sage: b1 = [3,2,5,4,6,8,7,10,9,1,0]
            sage: b0 = [5,4,3,2,6,8,7,10,9,1,0]
            sage: p1 = iet.Permutation(a,b1)
            sage: p1.arf_invariant()
            1
            sage: p0 = iet.Permutation(a,b0)
            sage: p0.arf_invariant()
            0
        """
        if any((z+1)%2 for z in self.profile()):
            return None

        from sage.rings.finite_rings.finite_field_constructor import GF
        GF2 = GF(2)

        M = self.intersection_matrix(GF2)
        F, C = M.symplectic_form()

        g = F.rank() // 2
        n = F.ncols()

        s = GF2(0)
        for i in range(g):
            a = C.row(i)

            a_indices = [k for k in range(n) if a[k]]
            t_a = GF2(len(a_indices))
            for j1 in range(len(a_indices)):
                for j2 in range(j1+1,len(a_indices)):
                    t_a += M[a_indices[j1], a_indices[j2]]

            b = C.row(g+i)
            b_indices = [k for k in range(n) if b[k]]
            t_b = GF2(len(b_indices))
            for j1 in range(len(b_indices)):
                for j2 in range(j1+1,len(b_indices)):
                    t_b += M[b_indices[j1],b_indices[j2]]

            s += t_a * t_b

        return s

    def stratum_component(self):
        r"""
        Returns a connected components of a stratum.

        EXAMPLES::

            sage: from surface_dynamics import *

        Permutations from the stratum H(6)::

            sage: a = list(range(8))
            sage: b_hyp = [7,6,5,4,3,2,1,0]
            sage: b_odd = [3,2,5,4,7,6,1,0]
            sage: b_even = [5,4,3,2,7,6,1,0]
            sage: p_hyp = iet.Permutation(a, b_hyp)
            sage: p_odd = iet.Permutation(a, b_odd)
            sage: p_even = iet.Permutation(a, b_even)
            sage: p_hyp.stratum_component()
            H_4(6)^hyp
            sage: p_odd.stratum_component()
            H_4(6)^odd
            sage: p_even.stratum_component()
            H_4(6)^even

        Permutations from the stratum H(4,4)::

            sage: a = list(range(11))
            sage: b_hyp = [10,9,8,7,6,5,4,3,2,1,0]
            sage: b_odd = [3,2,5,4,6,8,7,10,9,1,0]
            sage: b_even = [5,4,3,2,6,8,7,10,9,1,0]
            sage: p_hyp = iet.Permutation(a,b_hyp)
            sage: p_odd = iet.Permutation(a,b_odd)
            sage: p_even = iet.Permutation(a,b_even)
            sage: p_hyp.stratum() == Stratum([4,4], k=1)
            True
            sage: p_hyp.stratum_component()
            H_5(4^2)^hyp
            sage: p_odd.stratum() == Stratum([4,4], k=1)
            True
            sage: p_odd.stratum_component()
            H_5(4^2)^odd
            sage: p_even.stratum() == Stratum([4,4], k=1)
            True
            sage: p_even.stratum_component()
            H_5(4^2)^even

        As for stratum you can specify that you want to attach the singularity
        on the left of the interval using the option marked_separatrix::

            sage: a = list(range(1,10))
            sage: b_odd = [4,3,6,5,7,9,8,2,1]
            sage: b_even = [6,5,4,3,7,9,8,2,1]
            sage: p_odd = iet.Permutation(a,b_odd)
            sage: p_even = iet.Permutation(a,b_even)
            sage: p_odd.stratum_component()
            H_4(4, 2)^odd
            sage: p_even.stratum_component()
            H_4(4, 2)^even
        """
        # TODO: we reimplement the very same logic in flat_surfaces.separatrix_diagram.stratum_component()
        if not self.is_irreducible():
            return list(map(lambda x: x.stratum_component(), self.decompose()))

        stratum = self.stratum()

        if stratum.is_connected():
            return stratum.unique_component()

        if stratum.has_hyperelliptic_component():
            if self.is_hyperelliptic():
                return stratum.hyperelliptic_component()
            # TODO: we assume that the first entry is the hyperelliptic one
            cc = stratum.connected_components()
            if len(cc) == 2:
                return cc[1]

        # if we still have several components, they must be spin
        spin = self.arf_invariant()
        if spin == 0:
            return stratum.even_component()
        else:
            return stratum.odd_component()

    def order_of_rauzy_action(self, winner, side=None):
        r"""
        Returns the order of the action of a Rauzy move.

        INPUT:

        - ``winner`` - string ``'top'`` or ``'bottom'``

        - ``side`` - string ``'left'`` or ``'right'``

        OUTPUT:

        An integer corresponding to the order of the Rauzy action.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c d','d a c b')
            sage: p.order_of_rauzy_action('top', 'right')
            3
            sage: p.order_of_rauzy_action('bottom', 'right')
            2
            sage: p.order_of_rauzy_action('top', 'left')
            1
            sage: p.order_of_rauzy_action('bottom', 'left')
            3
        """
        winner = interval_conversion(winner)
        side = side_conversion(side)

        if side == -1:
            return self.length(winner) - self._twin[winner][-1] - 1
        elif side == 0:
            return self._twin[winner][0]

    def rauzy_move(self, winner, side='right', inplace=False):
        r"""
        Returns the permutation after a Rauzy move.

        INPUT:

        - ``winner`` - 'top' or 'bottom' interval

        - ``side`` - 'right' or 'left' (default: 'right') corresponding
          to the side on which the Rauzy move must be performed.

        - ``inplace`` - (default ``False``) whether the Rauzy move is
          performed inplace (to be used with care since permutations
          are hashable, set to ``True`` if you are sure to know what
          you are doing)

        OUTPUT:

        - a permutation

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b','b a')
            sage: p.rauzy_move(winner='top', side='right') == p
            True
            sage: p.rauzy_move(winner='bottom', side='right') == p
            True
            sage: p.rauzy_move(winner='top', side='left') == p
            True
            sage: p.rauzy_move(winner='bottom', side='left') == p
            True

        The options winner can be shortened to 't', 'b' and  'r', 'l'. As you
        can see in the following example::

            sage: p = iet.Permutation('a b c','c b a')
            sage: p.rauzy_move(winner='t', side='r')
            a b c
            c a b
            sage: p.rauzy_move(winner='b', side='r')
            a c b
            c b a
            sage: p.rauzy_move(winner='t', side='l')
            a b c
            b c a
            sage: p.rauzy_move(winner='b', side='l')
            b a c
            c b a

        This works as well for reduced permutations::

            sage: p = iet.Permutation('a b c d','d b c a',reduced=True)
            sage: p.rauzy_move('t')
            a b c d
            d a b c

        If Rauzy induction is not well defined, an error is raised::

            sage: p = iet.Permutation('a b', 'a b')
            sage: p.rauzy_move('t')
            Traceback (most recent call last):
            ...
            ValueError: Rauzy induction is not well defined


        Test the inplace option::

            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: q = p.rauzy_move('t', inplace=True)
            sage: assert q is p
            sage: p
            a b c d
            d a c b
            sage: q = p.rauzy_move('b', inplace=True)
            sage: assert q is p
        """
        winner = interval_conversion(winner)
        side = side_conversion(side)
        loser = 1 - winner

        if self._twin[0][side] == side or self._twin[0][side] == len(self)+side:
            raise ValueError("Rauzy induction is not well defined")

        if inplace:
            res = self
        else:
            res = self.__copy__()

        wtp = res._twin[winner][side]

        if side == -1:
            res._move(loser, len(self._twin[loser])-1, loser, wtp+1)

        if side == 0:
            res._move(loser, 0, loser, wtp)

        return res

    def backward_rauzy_move(self, winner, side='right', inplace=False):
        r"""
        Returns the permutation before a Rauzy move.

        INPUT:

        - ``winner`` - 'top' or 'bottom' interval

        - ``side`` - 'right' or 'left' (default: 'right') corresponding
          to the side on which the Rauzy move must be performed.

        - ``inplace`` - (default ``False``) whether the Rauzy move is
          performed inplace (to be used with care since permutations
          are hashable, set to ``True`` if you are sure to know what
          you are doing)

        OUTPUT:

        - a permutation

        TESTS::

            sage: from surface_dynamics import *

        Testing the inversion on labelled permutations::

            sage: p = iet.Permutation('a b c d','d c b a')
            sage: for pos,side in [('t','r'),('b','r'),('t','l'),('b','l')]:
            ....:     q = p.rauzy_move(pos,side)
            ....:     print(q.backward_rauzy_move(pos,side) == p)
            ....:     q = p.backward_rauzy_move(pos,side)
            ....:     print(q.rauzy_move(pos,side) == p)
            True
            True
            True
            True
            True
            True
            True
            True

        Testing the inversion on reduced permutations::

            sage: p = iet.Permutation('a b c d','d c b a',reduced=True)
            sage: for pos,side in [('t','r'),('b','r'),('t','l'),('b','l')]:
            ....:     q = p.rauzy_move(pos,side)
            ....:     print(q.backward_rauzy_move(pos,side) == p)
            ....:     q = p.backward_rauzy_move(pos,side)
            ....:     print(q.rauzy_move(pos,side) == p)
            True
            True
            True
            True
            True
            True
            True
            True

        Test the inplace option::

            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: q = p.backward_rauzy_move('t', inplace=True)
            sage: assert q is p
            sage: p
            a b c d
            d b a c
            sage: q = p.backward_rauzy_move('t', inplace=True)
            sage: q = p.backward_rauzy_move('b', inplace=True)
            sage: assert q is p
            sage: q = p.rauzy_move('b', inplace=True)
            sage: q = p.rauzy_move('t', inplace=True)
            sage: q = p.rauzy_move('t', inplace=True)
            sage: p
            a b c d
            d c b a
        """
        winner = interval_conversion(winner)
        side = side_conversion(side)

        loser = 1 - winner
        winner_twin = self._twin[winner][side]
        d = len(self)

        if inplace:
            res = self
        else:
            res = copy(self)

        if side == -1:
            if self._labels is not None:
                res._labels[loser].append(res._labels[loser].pop(winner_twin+1))

            # move the element
            res._twin[loser].append(res._twin[loser].pop(winner_twin+1))

            # correction for the moved element
            res._twin[winner][res._twin[loser][-1]] = d

            # shift twins that are after the moved element
            for j in range(winner_twin + 1, d):
                res._twin[winner][res._twin[loser][j]] -= 1

        elif side == 0:
            if res._labels is not None:
                res._labels[loser].insert(
                        0,
                        res._labels[loser].pop(winner_twin-1))

            # move the element
            res._twin[loser].insert(
                    0,
                    res._twin[loser].pop(winner_twin-1))

            # correction for the moved element
            res._twin[winner][res._twin[loser][0]] = 0

            # unshift elements before the moved element
            for j in range(1, winner_twin):
                res._twin[winner][res._twin[loser][j]] += 1

        return res

    def erase_marked_points(self):
        r"""
        Returns a permutation equivalent to self but without marked points.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b','b a')
            sage: p.erase_marked_points()
            a b
            b a
            sage: p = iet.Permutation('a b1 b2 c d', 'd c b1 b2 a')
            sage: p.erase_marked_points()
            a b1 c d
            d c b1 a
            sage: p = iet.Permutation('a0 a1 b0 b1 c0 c1 d0 d1','d0 d1 c0 c1 b0 b1 a0 a1')
            sage: p.erase_marked_points()
            a0 b0 c0 d0
            d0 c0 b0 a0
            sage: p = iet.Permutation('a b y0 y1 x0 x1 c d','c x0 x1 a d y0 y1 b')
            sage: p.erase_marked_points()
            a b c d
            c a d b
            sage: p = iet.Permutation('a x y z b','b x y z a')
            sage: p.erase_marked_points()
            a b
            b a
            sage: p = iet.Permutation("0 1 2 3 4 5 6","6 0 3 2 4 1 5")
            sage: p.stratum()
            H_3(4, 0)
            sage: p.erase_marked_points().stratum()
            H_3(4)
        """
        if len(self) == 1:
            return self

        if not self.is_irreducible():
            raise ValueError("the permutation must be irreducible")

        tops = [True]*len(self)  # true if we keep and false if not
        bots = [True]*len(self)

        # remove the zeros which are not at the endpoints
        i = 0
        while i < len(self):
            i += 1
            while i < len(self) and self._twin[0][i] == self._twin[0][i-1]+1:
                tops[i] = False
                bots[self._twin[0][i]] = False
                i += 1

        # remove the fake zero on the left
        i0 = self._twin[1][0]-1
        i1 = self._twin[0][0]-1
        while i0>0 and i1>0 and self._twin[0][i0] == i1:
            tops[i0] = False
            bots[i1] = False
            i0 -= 1
            i1 -= 1

        # remove the fake zero on the right
        i0 = self._twin[1][-1]+1
        i1 = self._twin[0][-1]+1
        n = len(self)
        while i0<n and i1<n and self._twin[0][i0] == i1:
            tops[i0] = False
            bots[i1] = False
            i0 += 1
            i1 += 1


        top_labs = self[0]
        bot_labs = self[1]
        top = []
        bot = []
        for i in range(len(self)):
            if tops[i]:
                top.append(top_labs[i])
            if bots[i]:
                bot.append(bot_labs[i])

        # remove the fake zero on the left-right
        if len(top)>2:
            if top[-1] == bot[0] and bot[-1] != top[0]:
                if bot[1] == top[0] and bot[-1] == top[-2]:
                    del bot[-1]
                    del bot[1]
                    del top[-2]
                    bot.append(top[0])

            elif top[-1] != bot[0] and bot[-1] == top[0]:
                    if top[1] == bot[0] and top[-1] == bot[-2]:
                        del top[-1]
                        del top[1]
                        del bot[-2]
                        top.append(bot[0])

            else:
                i0 = top.index(bot[-1])
                i1 = bot.index(top[-1])
                if bot[i1+1] == top[0] and top[i0+1] == bot[0]:
                    del top[i0+1]
                    del bot[i1]
                    del bot[0]
                    bot.insert(0, top[-1])

        return self.__class__((top,bot))

    def is_hyperelliptic(self):
        r"""
        Returns True if the permutation is in the class of the symmetric
        permutations (with eventual marked points).

        This is equivalent to say that the suspension lives in an hyperelliptic
        stratum of Abelian differentials H_hyp(2g-2) or H_hyp(g-1, g-1) with
        some marked points.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: iet.Permutation('a b c d','d c b a').is_hyperelliptic()
            True
            sage: iet.Permutation('0 1 2 3 4 5','5 2 1 4 3 0').is_hyperelliptic()
            False
        """
        test = self.erase_marked_points()

        n = test.length_top()
        cylindric = test.to_standard()
        return cylindric._twin[0] == list(range(n-1,-1,-1))

    def to_cylindric(self):
        r"""
        Returns a cylindric permutation in the same Rauzy class.

        A permutation is *cylindric* if the first letter in the top interval is
        also the last letter of the bottom interval or if the last letter of the
        top interval is the first letter of the bottom interval.

        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c','c b a')
            sage: p.to_cylindric() == p
            True
            sage: p = iet.Permutation('a b c d','b d a c')
            sage: q = p.to_cylindric()
            sage: q[0][0] == q[1][-1] or q[1][0] == q[1][0]
            True
        """
        tmp = copy(self)
        n = self.length(0)

        a0 = tmp._twin[0][-1]
        a1 = tmp._twin[1][-1]
        p_min = min(a0,a1)

        while p_min > 0:
            if p_min == a0:
                k_min = min(tmp._twin[1][a0+1:])
                k = n - tmp._twin[1].index(k_min) - 1

                for j in range(k):
                    tmp = tmp.rauzy_move(0)

            else:
                k_min = min(tmp._twin[0][a1+1:])
                k = n - tmp._twin[0].index(k_min) - 1

                for j in range(k):
                    tmp = tmp.rauzy_move(1)

            a0 = tmp._twin[0][-1]
            a1 = tmp._twin[1][-1]
            p_min = min(a0,a1)

        return tmp

    def is_cylindric(self):
        r"""
        Returns True if the permutation is cylindric

        A permutation `\pi` is cylindric if `\pi(1) = n` or `\pi(n) = 1`. The
        name cylindric comes from geometry. A cylindric permutation has a
        suspension which is a flat surface with a completely periodic horizontal
        direction which is made of only one cylinder.


        EXAMPLES::

            sage: from surface_dynamics import *

            sage: iet.Permutation('1 2 3','3 2 1').is_cylindric()
            True
            sage: iet.Permutation('1 2 3','3 1 2').is_cylindric()
            True
            sage: iet.Permutation('1 2 3 4','3 1 2 4').is_cylindric()
            False
        """
        return self._twin[0][-1] == 0 or self._twin[1][-1] == 0

    def to_standard(self):
        r"""
        Returns a standard permutation in the same Rauzy class.

        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c','c b a')
            sage: p.to_standard() == p
            True
            sage: p = iet.Permutation('a b c d','b d a c')
            sage: q = p.to_standard()
            sage: q[0][0] == q[1][-1]
            True
            sage: q[1][0] == q[1][0]
            True
        """
        tmp = self.to_cylindric()
        n = len(self)

        a0 = tmp._twin[0][-1]
        a1 = tmp._twin[1][-1]
        p_min = min(a0,a1)

        if a0 == 0:
            for j in range(n - tmp._twin[1].index(0) - 1):
                tmp = tmp.rauzy_move(0)

        else:
            for j in range(n - tmp._twin[0].index(0) - 1):
                tmp = tmp.rauzy_move(1)

        return tmp

    def is_standard(self):
        r"""
        Test if the permutation is standard

        A permutation `\pi` is standard if '\pi(n) = 1` and `\pi(1) = n`.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c d','d c b a')
            sage: p.is_standard()
            True
            sage: p = p.rauzy_move('top')
            sage: p.is_standard()
            False
        """
        return self._twin[0][-1] == 0 and self._twin[1][-1] == 0

    def to_permutation(self):
        r"""
        Returns the permutation as an element of the symmetric group.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c','c b a')
            sage: p.to_permutation()
            [3, 2, 1]

        ::

            sage: p = Permutation([2,4,1,3])
            sage: q = iet.Permutation(p)
            sage: q.to_permutation() == p
            True
        """
        from sage.combinat.permutation import Permutation
        return Permutation(list(map(lambda x: x+1,self._twin[1])))

    def suspension_cone(self, winner=None):
        r"""
        Return the cone of suspension data.

        A suspension data `\tau` for a permutation `(\pi_{top}, \pi_{bot})`
        on the alphabet `\mathcal{A}` is a real vector in `RR^\mathcal{A}`
        so that

        .. MATH::

            \forall 1 \leq k < d,\,
            \sum_{\beta: \pi_{top}(\beta) \leq k} \tau_\beta > 0
            \quad \text{and} \quad
            \sum_{\beta: \pi_{bot}(\beta) \leq k} \tau_\beta < 0.

        A suspension data determines half of a zippered rectangle construction.
        The other half is the length data that is a positive vector in
        `\RR^\mathcal{A}`.

        INPUT:

        - ``winner`` - (optional) either ``None``, ``"top"`` or ``"bottom"``. If
          not ``None`` , then return only half of the suspension cone corresponding
          to data that either comes from a top or bottom Rauzy induction.

        .. SEEALSO::

            :meth:`heights_cone`

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: p = iet.Permutation('a b c d e f', 'e c b f d a')
            sage: H = p.suspension_cone()
            sage: H.dimension()
            6
            sage: rays = [r.vector() for r in H.rays()]
            sage: r = sum(randint(1,5)*ray for ray in rays)
            sage: r[0]>0 and r[0]+r[1] > 0 and r[0]+r[1]+r[2] > 0
            True
            sage: r[0]+r[1]+r[2]+r[3]>0
            True
            sage: r[0]+r[1]+r[2]+r[3]+r[4]>0
            True
            sage: r[4]<0 and r[4]+r[2]<0 and r[4]+r[2]+r[1] < 0
            True
            sage: r[4]+r[2]+r[1]+r[5]<0
            True
            sage: r[4]+r[2]+r[1]+r[5]+r[3]<0
            True

        The construction also works with reduced permutations (ie not carrying
        labels)::

            sage: p = iet.Permutation('a b c d', 'd c b a', reduced=True)
            sage: H = p.suspension_cone()
            sage: r = sum(r.vector() for r in H.rays())
            sage: r[0] > 0 and r[0]+r[1] > 0 and r[0]+r[1]+r[2] > 0
            True
            sage: r[3] < 0 and r[3]+r[2] < 0 and r[3]+r[2]+r[1] < 0
            True
        """
        n = len(self)
        ieqs = []

        if self._labels is not None:
            labels = self._labels
        else:
            labels = [list(range(n)), self._twin[1]]

        for i in range(1,len(self)):
            ieq = [0]*(n+1)
            for j in range(i):
                ieq[labels[0][j]+1] = 1
            ieqs.append(ieq)

            ieq = [0]*(n+1)
            for j in range(i):
                ieq[labels[1][j]+1] = -1
            ieqs.append(ieq)

        if winner is not None:
            winner = interval_conversion(winner)
            if winner == 0:
                # sum of heights is <= 0
                ieqs.append([0] + [-1] * len(self))
            elif winner == 1:
                # sum of heights is >= 0
                ieqs.append([0] + [1] * len(self))

        from sage.geometry.polyhedron.constructor import Polyhedron
        return Polyhedron(ieqs=ieqs)

    def heights_cone(self, side=None):
        r"""
        Return the cone of heights data.

        .. SEEALSO::

            :meth:`suspension_cone`

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: C = p.heights_cone()
            sage: C
            A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 1 vertex and 5 rays
            sage: C.rays_list()
            [[0, 0, 1, 1], [0, 1, 1, 0], [0, 1, 1, 1], [1, 1, 0, 0], [1, 1, 1, 0]]

            sage: p.heights_cone('top').rays_list()
            [[0, 0, 1, 1], [0, 1, 1, 0], [1, 1, 0, 0], [1, 1, 1, 0]]
            sage: p.heights_cone('bot').rays_list()
            [[0, 0, 1, 1], [0, 1, 1, 0], [0, 1, 1, 1], [1, 1, 0, 0]]
        """
        I = self.intersection_matrix()
        C = self.suspension_cone(side)

        from sage.geometry.polyhedron.constructor import Polyhedron
        return Polyhedron(rays=[-I*c.vector() for c in C.rays()])

    def invariant_density_rauzy(self, winner=None, var='x'):
        r"""
        Return the invariant density for the Rauzy induction.

        Goes via the zippered rectangle construction of [Vee1982]_.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: f = iet.Permutation('a b c d', 'd c b a').invariant_density_rauzy()
            sage: f
            (1)/((x2 + x3)*(x1 + x2)*(x1 + x2 + x3)*(x0 + x1)) + (1)/((x2 + x3)*(x1 + x2)*(x0 + x1)*(x0 + x1 + x2))

            sage: f_top = iet.Permutation('a b c d', 'd c b a').invariant_density_rauzy('top')
            sage: f_top
            (1)/((x2 + x3)*(x1 + x2)*(x0 + x1)*(x0 + x1 + x2))
            sage: f_bot = iet.Permutation('a b c d', 'd c b a').invariant_density_rauzy('bot')
            sage: f_bot
            (1)/((x2 + x3)*(x1 + x2)*(x1 + x2 + x3)*(x0 + x1))

            sage: f == f_bot + f_top
            True
        """
        from surface_dynamics.misc.additive_multivariate_generating_series import AdditiveMultivariateGeneratingSeriesRing
        from surface_dynamics.misc.linalg import cone_triangulate

        d = len(self)
        S = self.suspension_cone(winner=winner)
        Omega = self.intersection_matrix()
        M = AdditiveMultivariateGeneratingSeriesRing(var, d)

        ans = M.zero()
        hyperplane = sum(Omega.columns())
        fac = 1 / ZZ(d).factorial()
        for t in cone_triangulate(S, hyperplane):
            heights = [r * Omega for r in t]
            for h in heights: h.set_immutable()
            d = {}
            for h in heights:
                if h not in d: d[h] = ZZ.one()
                else: d[h] += ZZ.one()
            ans += M.term(ZZ.one(), d)

        return ans

    def to_origami(self):
        r"""
        Return the origami associated to a cylindric permutation.

        EXAMPLES::

            sage: from surface_dynamics import iet
            sage: p = iet.Permutation('a b', 'b a')
            sage: p.to_origami()
            (1)
            (1)

            sage: p = iet.Permutation('a b c e d f g', 'f e b g d c a')
            sage: p.to_origami()
            (1,2,3,4,5,6)
            (1,3,2,6,4,5)
            sage: assert p.stratum_component() == p.to_origami().stratum_component()
        """
        n = len(self._twin[0])
        if self._twin[1][-1] != 0:
            raise ValueError("to_origami is only valid for cylindric permutation")

        from surface_dynamics.misc.permutation import perm_invert
        from surface_dynamics.flat_surfaces.origamis.origami import Origami
        r = list(range(1, n-1))
        r.append(0)
        u = perm_invert([self._twin[1][i]-1 for i in range(n-1)])
        return Origami(r, u, as_tuple=True)

    def _masur_polygon_helper(self, lengths, heights):
        n = len(self)

        from sage.structure.sequence import Sequence

        s = Sequence(list(lengths) + list(heights))
        lengths = s[:len(lengths)]
        heights = s[len(lengths):]
        base_ring = s.universe()

        if self._labels is None:
            lengths = [lengths[:n]] + [lengths[j] for j in self._twins[1]]
            heights = [heights[:n]] + [heights[j] for j in self._twins[1]]
        else:
            lengths = [[lengths[j] for j in self._labels[i]] for i in [0, 1]]
            heights = [[heights[j] for j in self._labels[i]] for i in [0, 1]]

        try:
            zero = base_ring.zero()
        except AttributeError:
            zero = base_ring(0)

        # build the polygon in counter-clockwise order
        Ltop = [(zero,zero)]
        for i,dx,dy in zip(range(n), lengths[0], heights[0]):
            x, y = Ltop[-1]
            if dx <= 0 or (y <= 0 and i != 0):
                raise ValueError('invalid suspension data dx={} y={} at i={} on top'.format(dx, y, i))
            Ltop.append((x+dx, y+dy))
        Lbot = [(zero,zero)]
        for i,dx,dy in zip(range(n), lengths[1], heights[1]):
            x, y = Lbot[-1]
            if dx <= 0 or (y >= 0 and i != 0):
                raise ValueError('invalid suspension data dx={} y={} at i={} on bot'.format(dx, y, i))
            Lbot.append((x+dx, y+dy))

        assert Ltop[0] == Lbot[0] and Ltop[-1] == Lbot[-1], (Ltop, Lbot)
        Ltop.pop(0)
        Lbot.pop(0)
        endpoint = Ltop.pop(-1)
        Lbot.pop(-1)

        from flatsurf.geometry.polygon import Polygon

        ptop = Ltop[0]
        pbot = Lbot[0]
        triangles = [Polygon(vertices=[(zero,zero), pbot, ptop])]
        tops = [(0,2)]
        bots = [(0,0)]
        mids = [(0,1)]
        itop = 1
        ibot = 1
        k = 1
        while itop < len(Ltop) or ibot < len(Lbot):
            xtop = Ltop[itop][0] if itop < len(Ltop) else None
            xbot = Lbot[ibot][0] if ibot < len(Lbot) else None
            if xbot is None or (xtop is not None and xtop <= xbot):
                # add a triangle with a new vertex on top
                pptop = Ltop[itop]
                itop += 1
                triangles.append(Polygon(vertices=[ptop,pbot,pptop]))
                tops.append((k,2))
                mids.append((k,0))
                mids.append((k,1))
                ptop = pptop
            else:
                ppbot = Lbot[ibot]
                ibot += 1
                triangles.append(Polygon(vertices=[ptop,pbot,ppbot]))
                bots.append((k,1))
                mids.append((k,0))
                mids.append((k,2))
                pbot = ppbot
            k += 1

        triangles.append(Polygon(vertices=[ptop, pbot, endpoint]))
        tops.append((k, 2))
        bots.append((k, 1))
        mids.append((k, 0))

        assert len(tops) == len(bots) == n, (n, tops, bots)
        assert len(triangles) == 2*n-2, (n, len(triangles))

        return base_ring, triangles, tops, bots, mids

    def masur_polygon(self, lengths, heights):
        r"""
        Return the Masur polygon for the given ``lengths`` and ``heights``

        EXAMPLES::

            sage: from surface_dynamics import iet

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: S = p.masur_polygon([1,4,2], [2,0,-1])  # optional: sage_flatsurf
            sage: stratum = S.stratum()                   # optional: sage_flatsurf # random
            sage: stratum                                 # optional: sage_flatsurf
            H_1(0^2)
            sage: S                                       # optional: sage_flatsurf
            Translation Surface in H_1(0^2) built from 2 isosceles triangles and 2 triangles

        Generic construction using suspension cone::

            sage: p = iet.Permutation('a b c d e f g h i', 'b g f c a d i e h')
            sage: x = polygen(QQ)
            sage: poly = x^3 - x - 1
            sage: emb = AA.polynomial_root(poly, RIF(1.3,1.4))
            sage: K = NumberField(poly, 'a', embedding=emb)
            sage: a = K.gen()
            sage: R = [r.vector() for r in p.suspension_cone().rays()]
            sage: C = [1, a, a+1, 2-a, 2, 1, a, a, 1, a-1, 1]
            sage: H = sum(c*r for c,r in zip(C,R))
            sage: H
            (a + 2, -2, 2, -2*a + 2, 3*a, -a - 4, 0, -a + 1, -2*a - 1)
            sage: L = [1+a**2, 2*a**2-1, 1, 1, 1+a, a**2, a-1, a-1, 2]
            sage: S = p.masur_polygon(L, H)   # optional: sage_flatsurf
            sage: S                           # optional: sage_flatsurf
            Translation Surface in H_3(1^4) built from 15 triangles and a right triangle

        TESTS::

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: for L in [[1,4,2],[2,4,1],[5,1,1],[1,5,1],[1,1,5]]:  # optional: sage_flatsurf
            ....:     S = p.masur_polygon(L, [2,0,-1])
            ....:     assert S.stratum() == p.stratum()
        """
        base_ring, triangles, tops, bots, mids = self._masur_polygon_helper(lengths, heights)

        from flatsurf import MutableOrientedSimilaritySurface
        S = MutableOrientedSimilaritySurface(base_ring)
        for t in triangles:
            S.add_polygon(t)
        for i in range(len(self)):
            p1, e1 = tops[i]
            p2, e2 = bots[self._twin[0][i]]
            S.glue((p1, e1), (p2, e2))
        for i in range(0,len(mids),2):
            p1, e1 = mids[i]
            p2, e2 = mids[i+1]
            S.glue((p1, e1), (p2, e2))
        S.set_immutable()
        return S


class OrientablePermutationLI(PermutationLI):
    r"""
    Template for quadratic permutation.

    .. warning::

        Internal class! Do not use directly!

    AUTHOR:

    - Vincent Delecroix (2008-12-20): initial version

    """
    def rauzy_move(self, winner, side='right', inplace=False):
        r"""
        Returns the permutation after a Rauzy move.

        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('a a b','b c c',reduced=True)
            sage: p.rauzy_move(0)
            a a b
            b c c
            sage: p.rauzy_move(1)
            a a
            b b c c

        ::

            sage: p = iet.GeneralizedPermutation('a a b','b c c',reduced=True)
            sage: p.rauzy_move(0)
            a a b
            b c c
            sage: p.rauzy_move(1)
            a a
            b b c c

        ::

            sage: p = iet.GeneralizedPermutation('a a b','b c c',reduced=True)
            sage: pp = p.rauzy_move(0, inplace=True)
            sage: p
            a a b
            b c c
            sage: pp is p
            True
        """
        winner = interval_conversion(winner)
        side = side_conversion(side)
        loser = 1 - winner

        if inplace:
            res = self
        else:
            res = self.__copy__()

        wti, wtp = res._twin[winner][side]

        if side == -1:
            d = len(res._twin[loser])
            if wti == loser:
                res._move(loser, d-1, loser, wtp+1)
            else:
                res._move(loser, d-1, winner, wtp)

        if side == 0:
            if wti == loser:
                res._move(loser, 0, loser, wtp)
            else:
                res._move(loser, 0, winner, wtp+1)

        return res

    def backward_rauzy_move(self, winner, side='top'):
        r"""
        Return the permutation before the Rauzy move.

        TESTS::

            sage: from surface_dynamics import *

        Tests the inversion on labelled generalized permutations::

            sage: p = iet.GeneralizedPermutation('a a b b','c c d d')
            sage: for pos,side in [('t','r'),('b','r'),('t','l'),('b','l')]:
            ....:     q = p.rauzy_move(pos,side)
            ....:     print(q.backward_rauzy_move(pos,side) == p)
            ....:     q = p.backward_rauzy_move(pos,side)
            ....:     print(q.rauzy_move(pos,side) == p)
            True
            True
            True
            True
            True
            True
            True
            True

        Tests the inversion on reduced generalized permutations::

            sage: p = iet.GeneralizedPermutation('a a b b','c c d d',reduced=True)
            sage: for pos,side in [('t','r'),('b','r'),('t','l'),('b','l')]:
            ....:     q = p.rauzy_move(pos,side)
            ....:     print(q.backward_rauzy_move(pos,side) == p)
            ....:     q = p.backward_rauzy_move(pos,side)
            ....:     print(q.rauzy_move(pos,side) == p)
            True
            True
            True
            True
            True
            True
            True
            True
        """
        winner = interval_conversion(winner)
        side = side_conversion(side)
        loser = 1 - winner

        res = copy(self)

        wti, wtp = res._twin[winner][side]


        if side == -1:
            d = len(self._twin[loser])
            if wti == loser:
                res._move(loser, wtp+1, loser, d)
            else:
                res._move(winner, wtp-1, loser, d)

        if side == 0:
            if wti == loser:
                res._move(loser, wtp-1, loser, 0)
            else:
                res._move(winner, wtp+1, loser, 0)

        return res

    def stratum(self):
        r"""
        Returns the stratum associated to self

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('a b b','c c a')
            sage: p.stratum()
            Q_0(-1^4)
        """
        if self.is_irreducible():
            from surface_dynamics.flat_surfaces.quadratic_strata import Stratum
            return Stratum([x-2 for x in self.profile()], k=2)
        raise ValueError("stratum is well defined only for irreducible permutations")


FlippedPermutation = Permutation


class FlippedPermutationIET(PermutationIET):
    r"""
    Template for flipped Abelian permutations.

    .. warning::

        Internal class! Do not use directly!
    """
    def rauzy_move(self, winner, side='right', inplace=False):
        r"""
        Returns the permutation after a Rauzy move.

        TESTS::

            sage: from surface_dynamics import iet

            sage: p = iet.Permutation('a b c', 'c a b', flips=['b'], reduced=True)
            sage: p.rauzy_move('t','r')
             a -b  c
             c -b  a
            sage: p.rauzy_move('b','r')
             a -b -c
            -b  a -c
            sage: p.rauzy_move('t','l')
             a -b  c
             c  a -b
            sage: p.rauzy_move('b','l')
            -a  b  c
             c  b -a

            sage: p = iet.GeneralizedPermutation('a b c d','d a b c',flips='abcd')
            sage: p
            -a -b -c -d
            -d -a -b -c
            sage: p.rauzy_move('top','right')
            -a -b  c -d
             c -d -a -b
            sage: p.rauzy_move('bottom','right')
            -a -b  d -c
             d -a -b -c
            sage: p.rauzy_move('top','left')
            -a -b -c  d
            -a  d -b -c
            sage: p.rauzy_move('bottom','left')
            -b -c -d  a
            -d  a -b -c

            sage: p = iet.GeneralizedPermutation('a b c d','d a b c',flips='abcd')
            sage: pp = p.rauzy_move('top', 'right', inplace=True)
            sage: pp is p
            True
            sage: p
            -a -b  c -d
             c -d -a -b
        """
        winner = interval_conversion(winner)
        side = side_conversion(side)
        loser = 1 - winner

        if inplace:
            res = self
        else:
            res = self.__copy__()

        wtp = res._twin[winner][side]
        flip = self._flips[winner][side]
        if flip == -1:
            res._flips[loser][side] *= -1
            res._flips[winner][res._twin[loser][side]] *= -1
            flip = 1
        else:
            flip = 0

        if side == -1:
            d = len(res._twin[loser])
            res._move(loser, d-1, loser, wtp+1-flip)

        if side == 0:
            res._move(loser, 0, loser, wtp+flip)

        return res

    def backward_rauzy_move(self, winner, side=-1):
        r"""
        Returns the permutation before a Rauzy move.

        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('a b c d e','d a b e c', flips='abcd')
            sage: for pos,side in [('t','r'),('b','r'),('t','l'),('b','l')]:
            ....:     q = p.rauzy_move(pos,side)
            ....:     print(q.backward_rauzy_move(pos,side) == p)
            ....:     q = p.backward_rauzy_move(pos,side)
            ....:     print(q.rauzy_move(pos,side) == p)
            True
            True
            True
            True
            True
            True
            True
            True

        Testing the inversion on reduced permutations::

            sage: p = iet.Permutation('f a b c d e','d f c b e a', flips='abcd', reduced=True)
            sage: for pos,side in [('t','r'),('b','r'),('t','l'),('b','l')]:
            ....:     q = p.rauzy_move(pos,side)
            ....:     print(q.backward_rauzy_move(pos,side) == p)
            ....:     q = p.backward_rauzy_move(pos,side)
            ....:     print(q.rauzy_move(pos,side) == p)
            True
            True
            True
            True
            True
            True
            True
            True
        """
        winner = interval_conversion(winner)
        side = side_conversion(side)
        loser = 1 - winner

        res = copy(self)

        wtp = res._twin[winner][side]
        flip = self._flips[winner][side]

        if side == -1:
            d = len(self._twin[loser])

            if flip == -1:
                res._flips[loser][wtp-1] *= -1
                res._flips[winner][res._twin[loser][wtp-1]] *= -1
                res._move(loser, wtp-1, loser, d)
            else:
                res._move(loser, wtp+1, loser, d)

        if side == 0:
            if flip == -1:
                res._flips[loser][wtp+1] *= -1
                res._flips[winner][res._twin[loser][wtp+1]] *= -1
                res._move(loser, wtp+1, loser, 0)
            else:
                res._move(loser, wtp-1, loser, 0)

        return res


class FlippedPermutationLI(PermutationLI):
    r"""
    Template for flipped quadratic permutations.

    .. warning::

        Internal class! Do not use directly!

    AUTHORS:

    - Vincent Delecroix (2008-12-20): initial version

    """
    def rauzy_move(self, winner, side='right', inplace=False):
        r"""
        Rauzy move

        TESTS::

            sage: from surface_dynamics import iet

            sage: p = iet.GeneralizedPermutation('a b c b','d c d a',flips='abcd')
            sage: p
            -a -b -c -b
            -d -c -d -a
            sage: p.rauzy_move('top','right')
             a -b  a -c -b
            -d -c -d
            sage: p.rauzy_move('bottom','right')
             b -a  b -c
            -d -c -d -a
            sage: p.rauzy_move('top','left')
            -a -b -c -b
            -c  d -a  d
            sage: p.rauzy_move('bottom','left')
            -b -c -b
            -d -c  a -d  a

            sage: p = iet.GeneralizedPermutation('a b c b','d c d a',flips='abcd')
            sage: pp = p.rauzy_move('top', 'right', inplace=True)
            sage: p
             a -b  a -c -b
            -d -c -d
            sage: pp is p
            True
        """
        winner = interval_conversion(winner)
        side = side_conversion(side)
        loser = 1 - winner

        if inplace:
            res = self
        else:
            res = self.__copy__()

        wti, wtp = res._twin[winner][side]
        flip = res._flips[winner][side]
        if flip == -1:
            res._flips[loser][side] *= -1
            lti,ltp = res._twin[loser][side]
            res._flips[lti][ltp] *= -1
            flip = 1
        else:
            flip = 0

        if side == -1:
            d = len(res._twin[loser])
            if wti == loser:
                res._move(loser, d-1, loser, wtp+1-flip)
            else:
                res._move(loser, d-1, winner, wtp+flip)

        if side == 0:
            if wti == loser:
                res._move(loser, 0, loser, wtp+flip)
            else:
                res._move(loser, 0, winner, wtp+1-flip)

        return res

    def backward_rauzy_move(self, winner, side=-1):
        r"""
        Rauzy move

        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.GeneralizedPermutation('a b c e b','d c d a e',flips='abcd')
            sage: for pos,side in [('t','r'),('b','r'),('t','l'),('b','l')]:
            ....:     q = p.rauzy_move(pos,side)
            ....:     print(q.backward_rauzy_move(pos,side) == p)
            ....:     q = p.backward_rauzy_move(pos,side)
            ....:     print(q.rauzy_move(pos,side) == p)
            True
            True
            True
            True
            True
            True
            True
            True

        Testing the inversion on reduced permutations::

            sage: p = iet.GeneralizedPermutation('a b c e b','d c d a e',flips='abcd',reduced=True)
            sage: for pos,side in [('t','r'),('b','r'),('t','l'),('b','l')]:
            ....:     q = p.rauzy_move(pos,side)
            ....:     print(q.backward_rauzy_move(pos,side) == p)
            ....:     q = p.backward_rauzy_move(pos,side)
            ....:     print(q.rauzy_move(pos,side) == p)
            True
            True
            True
            True
            True
            True
            True
            True
        """
        winner = interval_conversion(winner)
        side = side_conversion(side)
        loser = 1 - winner

        res = copy(self)

        wti, wtp = res._twin[winner][side]
        flip = self._flips[winner][side]

        if side == -1:
            d = len(self._twin[loser])
            if wti == loser:
                if flip == -1:
                    res._flips[loser][wtp-1] *= -1
                    lti,ltp = res._twin[loser][wtp-1]
                    res._flips[lti][ltp] *= -1
                    res._move(loser, wtp-1, loser, d)
                else:
                    res._move(loser, wtp+1, loser, d)

            else:
                if flip == -1:
                    res._flips[winner][wtp+1] *= -1
                    lti,ltp = res._twin[winner][wtp+1]
                    res._flips[lti][ltp] *= -1
                    res._move(winner, wtp+1, loser, d)
                else:
                    res._move(winner, wtp-1, loser, d)


        if side == 0:
            if wti == loser:
                if flip == -1:
                    res._flips[loser][wtp+1] *= -1
                    lti,ltp = res._twin[loser][wtp+1]
                    res._flips[lti][ltp] *= -1
                    res._move(loser, wtp+1, loser, 0)
                else:
                    res._move(loser, wtp-1, loser, 0)
            else:
                if flip == -1:
                    res._flips[winner][wtp-1] *= -1
                    lti,ltp = res._twin[winner][wtp-1]
                    res._flips[lti][ltp] *= -1
                    res._move(winner, wtp-1, loser, 0)
                else:
                    res._move(winner, wtp+1, loser, 0)

        return res

class RauzyDiagram(SageObject):
    r"""
    Template for Rauzy diagrams.

    .. warning:

        Internal class! Do not use directly!

    AUTHORS:

    - Vincent Delecroix (2008-12-20): initial version

    """
    # TODO: pickle problem of Path (it does not understand what is its parent)
    __metaclass__ = NestedClassMetaclass

    class Path(SageObject):
        r"""
        Path in Rauzy diagram.

            A path in a Rauzy diagram corresponds to a subsimplex of the simplex of
            lengths. This correspondence is obtained via the Rauzy induction. To a
            idoc IET we can associate a unique path in a Rauzy diagram. This
            establishes a correspondence between infinite full path in Rauzy diagram
            and equivalence topologic class of IET.
        """
        def __init__(self, parent, *data):
            r"""
            Constructor of the path.

            TESTS::

                sage: from surface_dynamics import *

                sage: p = iet.Permutation('a b c', 'c b a')
                sage: r = p.rauzy_diagram()
                sage: g = r.path(p, 0, 1, 0); g
                Path of length 3 in a Rauzy diagram

            Check for trac ticket 8388::

                sage: loads(dumps(g)) == g
                True
            """
            self._parent = parent

            if len(data) == 0:
                raise ValueError("No empty data")

            start = data[0]
            if start not in self._parent:
                raise ValueError("Starting point not in this Rauzy diagram")

            self._start = self._parent._permutation_to_vertex(start)

            cur_vertex = self._start
            self._edge_types = []

            n = len(self._parent._edge_types)
            for i in data[1:]:
                if not isinstance(i, (int,Integer)): # try parent method
                    i = self._parent.edge_types_index(i)

                if i < 0 or i > n:
                    raise ValueError("indices must be integer between 0 and %d"%(n))
                neighbours = self._parent._succ[cur_vertex]
                if neighbours[i] is None:
                    raise ValueError("Invalid path")

                cur_vertex = neighbours[i]
                self._edge_types.append(i)

            self._end = cur_vertex

        def _repr_(self):
            r"""
            Returns a representation of the path.

            TESTS::

                sage: from surface_dynamics import *

                sage: p = iet.Permutation('a b','b a')
                sage: r = p.rauzy_diagram()
                sage: r.path(p)   #indirect doctest
                Path of length 0 in a Rauzy diagram
                sage: r.path(p,'top')   #indirect doctest
                Path of length 1 in a Rauzy diagram
                sage: r.path(p,'bottom')   #indirect doctest
                Path of length 1 in a Rauzy diagram
            """
            return "Path of length %d in a Rauzy diagram" %(len(self))

        def start(self):
            r"""
            Returns the first vertex of the path.

            EXAMPLES::

                sage: from surface_dynamics import *

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: g = r.path(p, 't', 'b')
                sage: g.start() == p
                True
            """
            return self._parent._vertex_to_permutation(self._start)

        def end(self):
            r"""
            Returns the last vertex of the path.

            EXAMPLES::

                sage: from surface_dynamics import *

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: g1 = r.path(p, 't', 'b', 't')
                sage: g1.end() == p
                True
                sage: g2 = r.path(p, 'b', 't', 'b')
                sage: g2.end() == p
                True
            """
            return self._parent._vertex_to_permutation(self._end)

        def edge_types(self):
            r"""
            Returns the edge types of the path.

            EXAMPLES::

                sage: from surface_dynamics import *

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: g = r.path(p, 0, 1)
                sage: g.edge_types()
                [0, 1]
            """
            return copy(self._edge_types)

        def __eq__(self, other):
            r"""
            Tests equality

            TESTS::

                sage: from surface_dynamics import *

                sage: p1 = iet.Permutation('a b','b a')
                sage: r1 = p1.rauzy_diagram()
                sage: p2 = p1.reduced()
                sage: r2 = p2.rauzy_diagram()
                sage: r1.path(p1,0,1) == r2.path(p2,0,1)
                False
                sage: r1.path(p1,0) == r1.path(p1,0)
                True
                sage: r1.path(p1,1) == r1.path(p1,0)
                False
            """
            return (
                type(self) == type(other) and
                self._parent == other._parent and
                self._start == other._start and
                self._edge_types == other._edge_types)

        def __ne__(self,other):
            r"""
            Tests inequality

            TESTS::

                sage: from surface_dynamics import *

                sage: p1 = iet.Permutation('a b','b a')
                sage: r1 = p1.rauzy_diagram()
                sage: p2 = p1.reduced()
                sage: r2 = p2.rauzy_diagram()
                sage: r1.path(p1,0,1) != r2.path(p2,0,1)
                True
                sage: r1.path(p1,0) != r1.path(p1,0)
                False
                sage: r1.path(p1,1) != r1.path(p1,0)
                True
            """
            return not self == other

        def __copy__(self):
            r"""
            Returns a copy of the path.

            TESTS::

                sage: from surface_dynamics import *

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: g1 = r.path(p,0,1,0,0)
                sage: g2 = copy(g1)
                sage: g1 is g2
                False
            """
            res = self.__class__(self._parent, self.start())
            res._edge_types = copy(self._edge_types)
            res._end = copy(self._end)
            return res

        def pop(self):
            r"""
            Pops the queue of the path

            OUTPUT:

            a path corresponding to the last edge

            EXAMPLES::

                sage: from surface_dynamics import *

                sage: p = iet.Permutation('a b','b a')
                sage: r = p.rauzy_diagram()
                sage: g = r.path(p,0,1,0)
                sage: g0,g1,g2,g3 = g[0], g[1], g[2], g[3]
                sage: g.pop() == r.path(g2,0)
                True
                sage: g == r.path(g0,0,1)
                True
                sage: g.pop() == r.path(g1,1)
                True
                sage: g == r.path(g0,0)
                True
                sage: g.pop() == r.path(g0,0)
                True
                sage: g == r.path(g0)
                True
                sage: g.pop() == r.path(g0)
                True
            """
            if len(self) == 0:
                return self._parent.path(self.start())

            else:
                e = self._edge_types.pop()
                self._end = self._parent._pred[self._end][e]
                return self._parent.path(self.end(),e)

        def append(self, edge_type):
            r"""
            Append an edge to the path.

            EXAMPLES::

                sage: from surface_dynamics import *

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: g = r.path(p)
                sage: g.append('top')
                sage: g
                Path of length 1 in a Rauzy diagram
                sage: g.append('bottom')
                sage: g
                Path of length 2 in a Rauzy diagram
            """
            if not isinstance(edge_type, (int,Integer)):
                edge_type = self._parent.edge_types_index(edge_type)

            elif edge_type < 0 or edge_type >= len(self._parent._edge_types):
                raise ValueError("invalid edge type")

            if self._parent._succ[self._end][edge_type] is None:
                raise ValueError("%d is not a valid edge" %(edge_type))

            self._edge_types.append(edge_type)
            self._end = self._parent._succ[self._end][edge_type]

        def _fast_append(self, edge_type):
            r"""
            Append an edge to the path without verification.

            EXAMPLES::

                sage: from surface_dynamics import *

                sage: p = iet.GeneralizedPermutation('a a','b b c c')
                sage: r = p.rauzy_diagram()

            .. try to add 1 with append::

                sage: g = r.path(p)
                sage: r[p][1] is None
                True
                sage: g.append(1)
                Traceback (most recent call last):
                ...
                ValueError: 1 is not a valid edge

            .. the same with fast append::

                sage: g = r.path(p)
                sage: r[p][1] is None
                True
                sage: g._fast_append(1)
            """
            self._edge_types.append(edge_type)
            self._end = self._parent._succ[self._end][edge_type]

        def extend(self, path):
            r"""
            Extends self with another path.

            EXAMPLES::

                sage: from surface_dynamics import *

                sage: p = iet.Permutation('a b c d','d c b a')
                sage: r = p.rauzy_diagram()
                sage: g1 = r.path(p,'t','t')
                sage: g2 = r.path(p.rauzy_move('t').rauzy_move('t'),'b','b')
                sage: g = r.path(p,'t','t','b','b')
                sage: g == g1 + g2
                True
                sage: g = copy(g1)
                sage: g.extend(g2)
                sage: g == g1 + g2
                True
            """
            if self._parent != path._parent:
                raise ValueError("not on the same Rauzy diagram")

            if self._end != path._start:
                raise ValueError("the end of the first path must the start of the second")

            self._edge_types.extend(path._edge_types)
            self._end = path._end

        def _fast_extend(self, path):
            r"""
            Extension with no verification.

            EXAMPLES::

                sage: from surface_dynamics import *

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: p0, p1 = r[p]
                sage: g = r.path(p)
                sage: g._fast_extend(r.path(p0))
                sage: g
                Path of length 0 in a Rauzy diagram
                sage: g._fast_extend(r.path(p1))
                sage: g
                Path of length 0 in a Rauzy diagram
            """
            self._edge_types.extend(path._edge_types)
            self._end = path._end

        def __len__(self):
            r"""
            Returns the length of the path.

            TESTS::

                sage: from surface_dynamics import *

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: len(r.path(p))
                0
                sage: len(r.path(p,0))
                1
                sage: len(r.path(p,1))
                1
            """
            return len(self._edge_types)

        def __getitem__(self, i):
            r"""
            TESTS::

                sage: from surface_dynamics import *

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: g = r.path(p,'t','b')
                sage: g[0] == p
                True
                sage: g[1] == p.rauzy_move('t')
                True
                sage: g[2] == p.rauzy_move('t').rauzy_move('b')
                True
                sage: g[-1] == g[2]
                True
                sage: g[-2] == g[1]
                True
                sage: g[-3] == g[0]
                True
            """
            if i > len(self) or i < -len(self)-1:
                raise IndexError("path index out of range")

            if i == 0: return self.start()
            if i < 0: i = i + len(self) + 1

            v = self._start
            for k in range(i):
                v = self._parent._succ[v][self._edge_types[k]]
            return self._parent._vertex_to_permutation(v)

        def __add__(self, other):
            r"""
            Concatenation of paths.

            EXAMPLES::

                sage: from surface_dynamics import *

                sage: p = iet.Permutation('a b','b a')
                sage: r = p.rauzy_diagram()
                sage: r.path(p) + r.path(p,'b') == r.path(p,'b')
                True
                sage: r.path(p,'b') + r.path(p) == r.path(p,'b')
                True
                sage: r.path(p,'t') + r.path(p,'b') == r.path(p,'t','b')
                True
            """
            if self._end != other._start:
                raise ValueError("the end of the first path is not the start of the second")

            res = copy(self)
            res._fast_extend(other)
            return res

        def __mul__(self, n):
            r"""
            Multiple of a loop.

            EXAMPLES::

                sage: from surface_dynamics import *

                sage: p = iet.Permutation('a b','b a')
                sage: r = p.rauzy_diagram()
                sage: l = r.path(p,'b')
                sage: l * 2 == r.path(p,'b','b')
                True
                sage: l * 3 == r.path(p,'b','b','b')
                True
            """
            if not self.is_loop():
                raise ValueError("Must be a loop to have multiple")

            if not isinstance(n, (Integer,int)):
                raise TypeError("The multiplier must be an integer")

            if n < 0:
                raise ValueError("The multiplier must be non negative")

            res = copy(self)
            for i in range(n-1):
                res += self
            return res

        def is_loop(self):
            r"""
            Tests whether the path is a loop (start point = end point).

            EXAMPLES::

                sage: from surface_dynamics import *

                sage: p = iet.Permutation('a b','b a')
                sage: r = p.rauzy_diagram()
                sage: r.path(p).is_loop()
                True
                sage: r.path(p,0,1,0,0).is_loop()
                True
            """
            return self._start == self._end

        def winners(self):
            r"""
            Returns the winner list associated to the edge of the path.

            EXAMPLES::

                sage: from surface_dynamics import *

                sage: p = iet.Permutation('a b','b a')
                sage: r = p.rauzy_diagram()
                sage: r.path(p).winners()
                []
                sage: r.path(p,0).winners()
                ['b']
                sage: r.path(p,1).winners()
                ['a']
            """
            return self.composition(
                self._parent.edge_to_winner,
                list.__add__)

        def losers(self):
            r"""
            Returns a list of the loosers on the path.

            EXAMPLES::

                sage: from surface_dynamics import *

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: g0 = r.path(p,'t','b','t')
                sage: g0.losers()
                ['a', 'c', 'b']
                sage: g1 = r.path(p,'b','t','b')
                sage: g1.losers()
                ['c', 'a', 'b']
            """
            return self.composition(
                self._parent.edge_to_loser,
                list.__add__)

        def __iter__(self):
            r"""
            Iterator over the permutations of the path.

            EXAMPLES::

                sage: from surface_dynamics import *

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: g = r.path(p)
                sage: for q in g:
                ....:     print(p)
                a b c
                c b a
                sage: g = r.path(p, 't', 't')
                sage: for q in g:
                ....:     print("%s\n*****" % q)
                a b c
                c b a
                *****
                a b c
                c a b
                *****
                a b c
                c b a
                *****
                sage: g = r.path(p,'b','t')
                sage: for q in g:
                ....:     print("%s\n*****" % q)
                a b c
                c b a
                *****
                a c b
                c b a
                *****
                a c b
                c b a
                *****
            """
            i = self._start

            for edge_type in self._edge_types:
                yield self._parent._vertex_to_permutation(i)
                i = self._parent._succ[i][edge_type]

            yield self.end()

        def composition(self, function, composition = None):
            r"""
            Compose an edges function on a path

            INPUT:

            - ``path`` - either a Path or a tuple describing a path

            - ``function`` - function must be of the form

            - ``composition`` - the composition function

            AUTHOR:

            - Vincent Delecroix (2009-09-29)

            EXAMPLES::

                sage: from surface_dynamics import *

                sage: p = iet.Permutation('a b','b a')
                sage: r = p.rauzy_diagram()
                sage: def f(i,t):
                ....:     if t is None: return []
                ....:     return [t]
                sage: g = r.path(p)
                sage: g.composition(f,list.__add__)
                []
                sage: g = r.path(p,0,1)
                sage: g.composition(f, list.__add__)
                [0, 1]
            """
            result = function(None,None)
            cur_vertex = self._start
            p = self._parent._element

            if composition is None: composition = result.__class__.__mul__

            for i in self._edge_types:
                self._parent._set_element(cur_vertex)
                result = composition(result, function(p,i))
                cur_vertex = self._parent._succ[cur_vertex][i]

            return result

        def right_composition(self, function, composition = None) :
            r"""
            Compose an edges function on a path

            INPUT:

            - ``function`` - function must be of the form (indice,type) ->
              element. Moreover function(None,None) must be an identity element
              for initialization.

            - ``composition`` - the composition function for the function. * if None (default: None)

            TESTS::

                sage: from surface_dynamics import *

                sage: p = iet.Permutation('a b','b a')
                sage: r = p.rauzy_diagram()
                sage: def f(i,t):
                ....:     if t is None: return []
                ....:     return [t]
                sage: g = r.path(p)
                sage: g.right_composition(f,list.__add__)
                []
                sage: g = r.path(p,0,1)
                sage: g.right_composition(f, list.__add__)
                [1, 0]
            """
            result = function(None,None)
            p = self._parent._element
            cur_vertex = self._start

            if composition is None: composition = result.__class__.__mul__

            for i in self._edge_types:
                self._parent._set_element(cur_vertex)
                result = composition(function(p,i),result)
                cur_vertex = self._parent._succ[cur_vertex][i]

            return result

    def __init__(self, p,
        right_induction=True,
        left_induction=False,
        left_right_inversion=False,
        top_bottom_inversion=False,
        symmetric=False):
        r"""
        self._succ contains successors
        self._pred contains predecessors

        self._element_class is the class of elements of self
        self._element is an instance of this class (hence contains the alphabet,
        the representation mode, ...). It is used to store data about property
        of permutations and also as a fast iterator.

         INPUT:

         - ``right_induction`` - boolean or 'top' or 'bottom': consider the
           right induction

         - ``left_induction`` - boolean or 'top' or 'bottom': consider the
           left induction

         - ``left_right_inversion`` - consider the left right inversion

         - ``top_bottom_inversion`` - consider the top bottom inversion

         - ``symmetric`` - consider the symmetric

        TESTS::

            sage: from surface_dynamics import *

            sage: r1 = iet.RauzyDiagram('a b','b a')
            sage: r2 = loads(dumps(r1))
        """
        self._edge_types = []
        self._index = {}

        if right_induction is True:
            self._index['rt_rauzy'] = len(self._edge_types)
            self._edge_types.append(('rauzy_move',(0,-1)))
            self._index['rb_rauzy'] = len(self._edge_types)
            self._edge_types.append(('rauzy_move',(1,-1)))

        elif isinstance(right_induction,str):
            if right_induction == '':
                raise ValueError("right_induction can not be empty string")

            elif 'top'.startswith(right_induction):
                self._index['rt_rauzy'] = len(self._edge_types)
                self._edge_types.append(('rauzy_move',(0,-1)))

            elif 'bottom'.startswith(right_induction):
                self._index['rb_rauzy'] = len(self._edge_types)
                self._edge_types.append(('rauzy_move',(1,-1)))

            else:
                raise ValueError("%s is not valid for right_induction" %(right_induction))

        if left_induction is True:
            self._index['lt_rauzy'] = len(self._edge_types)
            self._edge_types.append(('rauzy_move',(0,0)))
            self._index['lb_rauzy'] = len(self._edge_types)
            self._edge_types.append(('rauzy_move',(1,0)))

        elif isinstance(left_induction,str):
            if left_induction == '':
                raise ValueError("left_induction can not be empty string")

            elif 'top'.startswith(left_induction):
                self._index['lt_rauzy'] = len(self._edge_types)
                self._edge_types.append(('rauzy_move',(0,0)))

            elif 'bottom'.startswith(left_induction):
                self._index['lb_rauzy'] = len(self._edge_types)
                self._edge_types.append(('rauzy_move',(1,0)))

            else:
                raise ValueError("%s is not valid for left_induction" %(right_induction))

        if left_right_inversion is True:
            self._index['lr_inverse'] = len(self._edge_types)
            self._edge_types.append(('left_right_inverse',()))

        if top_bottom_inversion is True:
            self._index['tb_inverse'] =  len(self._edge_types)
            self._edge_types.append(('top_bottom_inverse',()))

        if symmetric is True:
            self._index['symmetric'] = len(self._edge_types)
            self._edge_types.append(('symmetric',()))

        self._n = len(p)
        self._element_class = p.__class__
        self._element = copy(p)
        self._alphabet = self._element._alphabet

        self._pred = {}
        self._succ = {}

        self.complete(p)

    def __eq__(self, other):
        r"""
        Tests equality.

        TESTS::

            sage: from surface_dynamics import *

            sage: iet.RauzyDiagram('a b','b a') == iet.RauzyDiagram('a b c','c b a')
            False
            sage: r = iet.RauzyDiagram('a b c','c b a')
            sage: r1 = iet.RauzyDiagram('a c b','c b a', alphabet='abc')
            sage: r2 = iet.RauzyDiagram('a b c','c a b', alphabet='abc')
            sage: r == r1
            True
            sage: r == r2
            True
            sage: r1 == r2
            True

        ::

            sage: r = iet.RauzyDiagram('a b c d','d c b a')
            sage: for p in r:
            ....:     p.rauzy_diagram() == r
            True
            True
            True
            True
            True
            True
            True
        """
        return type(self) == type(other) and \
               self._edge_types == other._edge_types and \
               next(iterkeys(self._succ)) in other._succ

    def __ne__(self, other):
        r"""
        Tests difference.

        TESTS::

            sage: from surface_dynamics import *

            sage: iet.RauzyDiagram('a b','b a') != iet.RauzyDiagram('a b c','c b a')
            True
            sage: r = iet.RauzyDiagram('a b c','c b a')
            sage: r1 = iet.RauzyDiagram('a c b','c b a', alphabet='abc')
            sage: r2 = iet.RauzyDiagram('a b c','c a b', alphabet='abc')
            sage: r != r1
            False
            sage: r != r2
            False
            sage: r1 != r2
            False
        """
        return not self == other

    def vertices(self):
        r"""
        Returns a list of the vertices.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: for p in r.vertices(): print(p)
            a b
            b a
        """
        return list(map(
            lambda x: self._vertex_to_permutation(x),
            self._succ.keys()))

    def vertex_iterator(self):
        r"""
        Returns an iterator over the vertices

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: for p in r.vertex_iterator(): print(p)
            a b
            b a

        ::

            sage: r = iet.RauzyDiagram('a b c d','d c b a')
            sage: r_1n = filter(lambda x: x.is_standard(), r)
            sage: for p in r_1n: print(p)
            a b c d
            d c b a
        """
        return map(
            lambda x: self._vertex_to_permutation(x),
            self._succ.keys())

    def edges(self,labels=True):
        r"""
        Returns a list of the edges.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: len(r.edges())
            2
        """
        return list(self.edge_iterator())

    def edge_iterator(self):
        r"""
        Returns an iterator over the edges of the graph.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b','b a')
            sage: r = p.rauzy_diagram()
            sage: for e in r.edge_iterator():
            ....:     print('%s --> %s' %(e[0].str(sep='/'), e[1].str(sep='/')))
            a b/b a --> a b/b a
            a b/b a --> a b/b a
        """
        for x in self._succ.keys():
            for i,y in enumerate(self._succ[x]):
                if y is not None:
                    yield(
                        (self._vertex_to_permutation(x),
                         self._vertex_to_permutation(y),
                         i))

    def edge_types_index(self, data):
        r"""
        Try to convert the data as an edge type.

        INPUT:

        - ``data`` - a string

        OUTPUT:

        integer

        EXAMPLES:

            sage: from surface_dynamics import *

        For a standard Rauzy diagram (only right induction) the 0 index
        corresponds to the 'top' induction and the index 1 corresponds to the
        'bottom' one::

            sage: p = iet.Permutation('a b c','c b a')
            sage: r = p.rauzy_diagram()
            sage: r.edge_types_index('top')
            0
            sage: r[p][0] == p.rauzy_move('top')
            True
            sage: r.edge_types_index('bottom')
            1
            sage: r[p][1] == p.rauzy_move('bottom')
            True

        The special operations (inversion and symmetry) always appears after the
        different Rauzy inductions::

            sage: p = iet.Permutation('a b c','c b a')
            sage: r = p.rauzy_diagram(symmetric=True)
            sage: r.edge_types_index('symmetric')
            2
            sage: r[p][2] == p.symmetric()
            True

        This function always try to resolve conflictuous name. If it's
        impossible a ValueError is raised::

            sage: p = iet.Permutation('a b c','c b a')
            sage: r = p.rauzy_diagram(left_induction=True)
            sage: r.edge_types_index('top')
            Traceback (most recent call last):
            ...
            ValueError: left and right inductions must be differentiated
            sage: r.edge_types_index('top_right')
            0
            sage: r[p][0] == p.rauzy_move(0)
            True
            sage: r.edge_types_index('bottom_left')
            3
            sage: r[p][3] == p.rauzy_move('bottom', 'left')
            True

        ::

            sage: p = iet.Permutation('a b c','c b a')
            sage: r = p.rauzy_diagram(left_right_inversion=True,top_bottom_inversion=True)
            sage: r.edge_types_index('inversion')
            Traceback (most recent call last):
            ...
            ValueError: left-right and top-bottom inversions must be differentiated
            sage: r.edge_types_index('lr_inverse')
            2
            sage: p.lr_inverse() == r[p][2]
            True
            sage: r.edge_types_index('tb_inverse')
            3
            sage: p.tb_inverse() == r[p][3]
            True

        Short names are accepted::

            sage: p = iet.Permutation('a b c','c b a')
            sage: r = p.rauzy_diagram(right_induction='top',top_bottom_inversion=True)
            sage: r.edge_types_index('top_rauzy_move')
            0
            sage: r.edge_types_index('t')
            0
            sage: r.edge_types_index('tb')
            1
            sage: r.edge_types_index('inversion')
            1
            sage: r.edge_types_index('inverse')
            1
            sage: r.edge_types_index('i')
            1
        """
        if not isinstance(data,str):
            raise ValueError("the edge type must be a string")

        if ('top_rauzy_move'.startswith(data) or
            't_rauzy_move'.startswith(data)):
            if 'lt_rauzy' in self._index:
                if 'rt_rauzy' in self._index:
                    raise ValueError("left and right inductions must be differentiated")
                return self._index['lt_rauzy']

            if 'rt_rauzy' in self._index:
                return self._index['rt_rauzy']

            raise ValueError("no top induction in this Rauzy diagram")

        if ('bottom_rauzy_move'.startswith(data) or
            'b_rauzy_move'.startswith(data)):
            if 'lb_rauzy' in self._index:
                if 'rb_rauzy' in self._index:
                    raise ValueError("left and right inductions must be differentiated")
                return self._index['lb_rauzy']

            if 'rb_rauzy' in self._index:
                return self._index['rb_rauzy']

            raise ValueError("no bottom Rauzy induction in this diagram")

        if ('left_rauzy_move'.startswith(data) or
            'l_rauzy_move'.startswith(data)):
            if 'lt_rauzy' in self._index:
                if 'lb_rauzy' in self._index:
                    raise ValueError("top and bottom inductions must be differentiated")
                return self._index['lt_rauzy']

            if 'lb_rauzy' in self._index:
                return self._index('lb_rauzy')

            raise ValueError("no left Rauzy induction in this diagram")

        if ('lt_rauzy_move'.startswith(data) or
            'tl_rauzy_move'.startswith(data) or
            'left_top_rauzy_move'.startswith(data) or
            'top_left_rauzy_move'.startswith(data)):
            if 'lt_rauzy' not in self._index:
                raise ValueError("no top-left Rauzy induction in this diagram")
            else:
                return self._index['lt_rauzy']

        if ('lb_rauzy_move'.startswith(data) or
            'bl_rauzy_move'.startswith(data) or
            'left_bottom_rauzy_move'.startswith(data) or
            'bottom_left_rauzy_move'.startswith(data)):
            if 'lb_rauzy' not in self._index:
                raise ValueError("no bottom-left Rauzy induction in this diagram")
            else:
                return self._index['lb_rauzy']

        if 'right'.startswith(data):
            raise ValueError("ambiguity with your edge name: %s" %(data))

        if ('rt_rauzy_move'.startswith(data) or
            'tr_rauzy_move'.startswith(data) or
            'right_top_rauzy_move'.startswith(data) or
            'top_right_rauzy_move'.startswith(data)):
            if 'rt_rauzy' not in self._index:
                raise ValueError("no top-right Rauzy induction in this diagram")
            else:
                return self._index['rt_rauzy']

        if ('rb_rauzy_move'.startswith(data) or
            'br_rauzy_move'.startswith(data) or
            'right_bottom_rauzy_move'.startswith(data) or
            'bottom_right_rauzy_move'.startswith(data)):
            if 'rb_rauzy' not in self._index:
                raise ValueError("no bottom-right Rauzy induction in this diagram")
            else:
                return self._index['rb_rauzy']

        if 'symmetric'.startswith(data):
            if 'symmetric' not in self._index:
                raise ValueError("no symmetric in this diagram")
            else:
                return self._index['symmetric']

        if 'inversion'.startswith(data) or data == 'inverse':
            if 'lr_inverse' in self._index:
                if 'tb_inverse' in self._index:
                    raise ValueError("left-right and top-bottom inversions must be differentiated")
                return self._index['lr_inverse']

            if 'tb_inverse' in self._index:
                return self._index['tb_inverse']

            raise ValueError("no inversion in this diagram")

        if ('lr_inversion'.startswith(data) or
            data == 'lr_inverse' or
            'left_right_inversion'.startswith(data) or
            data == 'left_right_inverse'):
            if 'lr_inverse' not in self._index:
                raise ValueError("no left-right inversion in this diagram")
            else:
                return self._index['lr_inverse']

        if ('tb_inversion'.startswith(data) or
            data == 'tb_inverse' or
            'top_bottom_inversion'.startswith(data)
            or data == 'top_bottom_inverse'):
            if 'tb_inverse' not in self._index:
                raise ValueError("no top-bottom inversion in this diagram")
            else:
                return self._index['tb_inverse']

        raise ValueError("this edge type does not exist: %s" %(data))

    def edge_types(self):
        r"""
        Print information about edges.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: r = iet.RauzyDiagram('a b', 'b a')
            sage: r.edge_types()
            0: rauzy_move(0, -1)
            1: rauzy_move(1, -1)

        ::

            sage: r = iet.RauzyDiagram('a b', 'b a', left_induction=True)
            sage: r.edge_types()
            0: rauzy_move(0, -1)
            1: rauzy_move(1, -1)
            2: rauzy_move(0, 0)
            3: rauzy_move(1, 0)

        ::

            sage: r = iet.RauzyDiagram('a b',' b a',symmetric=True)
            sage: r.edge_types()
            0: rauzy_move(0, -1)
            1: rauzy_move(1, -1)
            2: symmetric()
        """
        for i,(edge_type,t) in enumerate(self._edge_types):
            print(str(i) + ": " + edge_type + str(t))

    def alphabet(self, data=None):
        r"""
        TESTS::

            sage: from surface_dynamics import *

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: r.alphabet() == Alphabet(['a','b'])
            True
            sage: r = iet.RauzyDiagram([0,1],[1,0])
            sage: r.alphabet() == Alphabet([0,1])
            True
        """
        if data is None:
            return self._element._alphabet
        else:
            self._element._set_alphabet(data)

    def letters(self):
        r"""
        Returns the letters used by the RauzyDiagram.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: r.alphabet()
            {'a', 'b'}
            sage: r.letters()
            ['a', 'b']
            sage: r.alphabet('ABCDEF')
            sage: r.alphabet()
            {'A', 'B', 'C', 'D', 'E', 'F'}
            sage: r.letters()
            ['A', 'B']
        """
        return self._element.letters()

    def _vertex_to_permutation(self,data=None):
        r"""
        Converts the (vertex) data to a permutation.

        TESTS:

            sage: from surface_dynamics import *

            sage: r = iet.RauzyDiagram('a b','b a')   #indirect doctest
        """
        if data is not None:
            self._set_element(data)
            return copy(self._element)

    def edge_to_matrix(self, p=None, edge_type=None):
        r"""
        Return the corresponding matrix

        INPUT:

        - ``p`` - a permutation

        - ``edge_type`` - 0 or 1 corresponding to the type of the edge

        OUTPUT:

        A matrix

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c','c b a')
            sage: d = p.rauzy_diagram()
            sage: print(d.edge_to_matrix(p,1))
            [1 0 1]
            [0 1 0]
            [0 0 1]
        """
        if p is None and edge_type is None:
            return identity_matrix(self._n)

        function_name = self._edge_types[edge_type][0] + '_matrix'

        if not hasattr(self._element_class,function_name):
            return identity_matrix(self._n)

        arguments = self._edge_types[edge_type][1]

        return getattr(p,function_name)(*arguments)

    def edge_to_winner(self, p=None, edge_type=None):
        r"""
        Return the corresponding winner

        TESTS::

            sage: from surface_dynamics import *

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: r.edge_to_winner(None,None)
            []
        """
        if p is None and edge_type is None:
            return []

        function_name = self._edge_types[edge_type][0] + '_winner'

        if not hasattr(self._element_class, function_name):
            return [None]

        arguments = self._edge_types[edge_type][1]

        return [getattr(p,function_name)(*arguments)]

    def edge_to_loser(self, p=None, edge_type=None):
        r"""
        Return the corresponding loser

        TESTS::

            sage: from surface_dynamics import *

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: r.edge_to_loser(None,None)
            []
        """
        if p is None and edge_type is None:
            return []

        function_name = self._edge_types[edge_type][0] + '_loser'

        if not hasattr(self._element_class, function_name):
            return [None]

        arguments = self._edge_types[edge_type][1]

        return [getattr(p,function_name)(*arguments)]

    def _all_npath_extension(self, g, length=0):
        r"""
        Returns an iterator over all extension of fixed length of p.

        INPUT:

        - ``p`` - a path

        - ``length`` - a non-negative integer

        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b','b a')
            sage: r = p.rauzy_diagram()
            sage: g0 = r.path(p)
            sage: for g in r._all_npath_extension(g0,0):
            ....:     print(g)
            Path of length 0 in a Rauzy diagram
            sage: for g in r._all_npath_extension(g0,1):
            ....:     print(g)
            Path of length 1 in a Rauzy diagram
            Path of length 1 in a Rauzy diagram
            sage: for g in r._all_npath_extension(g0,2):
            ....:     print(g)
            Path of length 2 in a Rauzy diagram
            Path of length 2 in a Rauzy diagram
            Path of length 2 in a Rauzy diagram
            Path of length 2 in a Rauzy diagram

        ::

            sage: p = iet.GeneralizedPermutation('a a','b b c c',reduced=True)
            sage: r = p.rauzy_diagram()
            sage: g0 = r.path(p)
            sage: len(list(r._all_npath_extension(g0,0)))
            1
            sage: len(list(r._all_npath_extension(g0,1)))
            1
            sage: len(list(r._all_npath_extension(g0,2)))
            2
            sage: len(list(r._all_npath_extension(g0,3)))
            3
            sage: len(list(r._all_npath_extension(g0,4)))
            5
        """
        if length == 0:
            yield g

        else:
            for i in range(len(self._edge_types)):
                if self._succ[g._end][i] is not None:
                    g._fast_append(i)
                    for h in self._all_npath_extension(g,length-1): yield h
                    g.pop()

    def _all_path_extension(self, g, length=0):
        r"""
        Returns an iterator over all path extension of p.

        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b','b a')
            sage: r = p.rauzy_diagram()
            sage: g0 = r.path(p)
            sage: for g in r._all_path_extension(g0,0):
            ....:     print(g)
            Path of length 0 in a Rauzy diagram
            sage: for g in r._all_path_extension(g0, 1):
            ....:     print(g)
            Path of length 0 in a Rauzy diagram
            Path of length 1 in a Rauzy diagram
            Path of length 1 in a Rauzy diagram

        ::

            sage: p = iet.GeneralizedPermutation('a a','b b c c')
            sage: r = p.rauzy_diagram()
            sage: g0 = r.path(p)
            sage: len(list(r._all_path_extension(g0,0)))
            1
            sage: len(list(r._all_path_extension(g0,1)))
            2
            sage: len(list(r._all_path_extension(g0,2)))
            4
            sage: len(list(r._all_path_extension(g0,3)))
            7
        """
        yield g

        if length > 0:
            for i in range(len(self._edge_types)):
                if self._succ[g._end][i] is not None:
                    g._fast_append(i)
                    for h in self._all_path_extension(g,length-1): yield h
                    g.pop()

    def __iter__(self):
        r"""
        Iterator over the permutations of the Rauzy diagram.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: for p in r: print(p)
            a b
            b a
            sage: r = iet.RauzyDiagram('a b c','c b a')
            sage: for p in r: print(p.stratum())
            H_1(0^2)
            H_1(0^2)
            H_1(0^2)
        """
        for data in iterkeys(self._succ):
            yield self._vertex_to_permutation(data)

    def __contains__(self, element):
        r"""
        Containance test.

        TESTS::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c d','d c b a',reduced=True)
            sage: q = iet.Permutation('a b c d','d b c a',reduced=True)
            sage: r = p.rauzy_diagram()
            sage: s = q.rauzy_diagram()
            sage: p in r
            True
            sage: p in s
            False
            sage: q in r
            False
            sage: q in s
            True
        """
        for p in iterkeys(self._succ):
            if self._vertex_to_permutation(p) == element:
                return True

        return False

    def _repr_(self):
        r"""
        Returns a representation of self

        TESTS::

            sage: from surface_dynamics import *

            sage: iet.RauzyDiagram('a b','b a')   #indirect doctest
            Rauzy diagram with 1 permutation
            sage: iet.RauzyDiagram('a b c','c b a')   #indirect doctest
            Rauzy diagram with 3 permutations
        """
        if len(self._succ) == 0:
            return "Empty Rauzy diagram"
        elif len(self._succ) == 1:
            return "Rauzy diagram with 1 permutation"
        else:
            return "Rauzy diagram with %d permutations" %(len(self._succ))

    def __getitem__(self,p):
        r"""
        Returns the neighbors of p.

        Just use the function vertex_to_permutation that must be defined
        in each child.

        INPUT:

        - ``p`` - a permutation in the Rauzy diagram

        TESTS::

            sage: from surface_dynamics import *


            sage: p = iet.Permutation('a b c','c b a')
            sage: p0 = iet.Permutation('a b c','c a b',alphabet="abc")
            sage: p1 = iet.Permutation('a c b','c b a',alphabet="abc")
            sage: r = p.rauzy_diagram()
            sage: r[p] == [p0, p1]
            True
            sage: r[p1] == [p1, p]
            True
            sage: r[p0] == [p, p0]
            True
        """
        if not isinstance(p, self._element_class):
            raise ValueError("Your element does not have the good type")

        perm = self._permutation_to_vertex(p)
        return list(map(lambda x: self._vertex_to_permutation(x),
            self._succ[perm]))

    def __len__(self):
        r"""
        Deprecated use cardinality.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: r.cardinality()
            1
        """
        return self.cardinality()

    def cardinality(self):
        r"""
        Returns the number of permutations in this Rauzy diagram.

        OUTPUT:

        - `integer` - the number of vertices in the diagram

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: r.cardinality()
            1
            sage: r = iet.RauzyDiagram('a b c','c b a')
            sage: r.cardinality()
            3
            sage: r = iet.RauzyDiagram('a b c d','d c b a')
            sage: r.cardinality()
            7
        """
        return len(self._succ)

    def complete(self, p):
        r"""
        Completion of the Rauzy diagram.

        Add to the Rauzy diagram all permutations that are obtained by
        successive operations defined by edge_types(). The permutation must be
        of the same type and the same length as the one used for the creation.

        INPUT:

        - ``p`` - a permutation of Interval exchange transformation

        Rauzy diagram is the reunion of all permutations that could be
        obtained with successive rauzy moves. This function just use the
        functions __getitem__ and has_rauzy_move and rauzy_move which must
        be defined for child and their corresponding permutation types.

        TESTS::

            sage: from surface_dynamics import *

            sage: r = iet.RauzyDiagram('a b c','c b a')   #indirect doctest
            sage: r = iet.RauzyDiagram('a b c','c b a',left_induction=True) #indirect doctest
            sage: r = iet.RauzyDiagram('a b c','c b a',symmetric=True)   #indirect doctest
            sage: r = iet.RauzyDiagram('a b c','c b a',lr_inversion=True)   #indirect doctest
            sage: r = iet.RauzyDiagram('a b c','c b a',tb_inversion=True)   #indirect doctest
        """
        if p.__class__ is not self._element_class:
            raise TypeError("wrong permutation type")

        if len(p) != self._n:
            raise ValueError("wrong length")

        pred = self._pred
        succ = self._succ
        p = self._permutation_to_vertex(p)
        perm = self._element
        l = []

        if p not in succ:
            succ[p] = [None] * len(self._edge_types)
            pred[p] = [None] * len(self._edge_types)
            l.append(p)

        while l:
            p = l.pop()
            self._set_element(p)

            for t,edge in enumerate(self._edge_types):
                if (not hasattr(perm, 'has_'+edge[0]) or
                  getattr(perm, 'has_'+edge[0])(*(edge[1]))):
                    q = getattr(perm,edge[0])(*(edge[1]))
                    q = self._permutation_to_vertex(q)
                    if q not in succ:
                        succ[q] = [None] * len(self._edge_types)
                        pred[q] = [None] * len(self._edge_types)
                        l.append(q)
                    succ[p][t] = q
                    pred[q][t] = p

    def path(self, *data):
        r"""
        Returns a path over this Rauzy diagram.

        INPUT:

        - ``initial_vertex`` - the initial vertex (starting point of the path)

        - ``data`` - a sequence of edges

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c','c b a')
            sage: r = p.rauzy_diagram()
            sage: g = r.path(p, 'top', 'bottom')
        """
        if len(data) == 0:
            raise TypeError("Must be non empty")
        elif len(data) == 1 and isinstance(data[0], self.Path):
                return copy(data[0])
        return self.Path(self,*data)

    def graph(self):
        r"""
        Returns the Rauzy diagram as a Graph object

        The graph returned is more precisely a DiGraph (directed graph) with
        loops and multiedges allowed.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: r = iet.RauzyDiagram('a b c','c b a')
            sage: r
            Rauzy diagram with 3 permutations
            sage: r.graph()
            Looped digraph on 3 vertices
        """
        from sage.graphs.digraph import DiGraph

        if len(self._element) != 2:
            G = DiGraph(loops=True,multiedges=False)
        else:
            G = DiGraph(loops=True,multiedges=True)

        for p,neighbours in iteritems(self._succ):
            p = self._vertex_to_permutation(p)
            for i,n in enumerate(neighbours):
                if n is not None:
                    q = self._vertex_to_permutation(n)
                    G.add_edge(p,q,i)

        return G

class FlippedRauzyDiagram(RauzyDiagram):
    r"""
    Template for flipped Rauzy diagrams.

    .. warning:

        Internal class! Do not use directly!

    AUTHORS:

    - Vincent Delecroix (2009-09-29): initial version

    """
    def complete(self, p, reducible=False):
        r"""
        Completion of the Rauzy diagram

        Add all successors of p for defined operations in edge_types. Could be
        used for generating non (strongly) connected Rauzy diagrams. Sometimes,
        for flipped permutations, the maximal connected graph in all
        permutations is not strongly connected. Finding such components needs to
        call most than once the .complete() method.

        INPUT:

        - ``p`` - a permutation

        - ``reducible`` - put or not reducible permutations

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c','c b a',flips='a')
            sage: d = p.rauzy_diagram()
            sage: d
            Rauzy diagram with 3 permutations
            sage: p = iet.Permutation('a b c','c b a',flips='b')
            sage: d.complete(p)
            sage: d
            Rauzy diagram with 8 permutations
            sage: p = iet.Permutation('a b c','c b a',flips='a')
            sage: d.complete(p)
            sage: d
            Rauzy diagram with 8 permutations
        """
        if p.__class__ is not self._element_class:
            raise ValueError("your permutation is not of good type")

        if len(p) != self._n:
            raise ValueError("your permutation has not the good length")

        pred = self._pred
        succ = self._succ
        p = self._permutation_to_vertex(p)
        l = []

        if p not in succ:
            succ[p] = [None] * len(self._edge_types)
            pred[p] = [None] * len(self._edge_types)
            l.append(p)

        while l:
            p = l.pop()

            for t,edge_type in enumerate(self._edge_types):
                perm = self._vertex_to_permutation(p)

                if (not hasattr(perm,'has_' + edge_type[0]) or
                  getattr(perm, 'has_' + edge_type[0])(*(edge_type[1]))):
                    q = perm.rauzy_move(t)
                    q = self._permutation_to_vertex(q)
                    if reducible is True or perm.is_irreducible():
                        if q not in succ:
                            succ[q] = [None] * len(self._edge_types)
                            pred[q] = [None] * len(self._edge_types)
                            l.append(q)

                        succ[p][t] = q
                        pred[q][t] = p
