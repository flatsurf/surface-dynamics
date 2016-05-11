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

- construct dynamic graphs

- construct coherent _repr_

"""
#*****************************************************************************
#       Copyright (C) 2008 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject

from copy import copy

from sage.rings.integer import Integer
from sage.matrix.constructor import identity_matrix, matrix
from sage.misc.nested_class import NestedClassMetaclass


def interval_conversion(interval=None):
    r"""
    Converts the argument in 0 or 1.

    INPUT:

    - ``winner`` - 'top' (or 't' or 0) or bottom (or 'b' or 1)

    OUTPUT:

    integer -- 0 or 1

    TESTS::

        sage: from surface_dynamics.all import *

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

class Permutation(SageObject):
    r"""
    Template for all permutations.

    .. warning::

        Internal class! Do not use directly!

    This class implement generic algorithm (stratum, connected component, ...)
    and unfies all its children.

    It has three attributes

     - ``_twin`` -- the permutation

     - ``_labels`` -- None or the list of labels

     - ``_flips`` -- None or the list of flips (each flip is either ``1`` or
       ``-1``)

    The datatype for ``_twin`` differs for IET and LI.
    """
    _twin = None
    _labels = None
    _flips = None

    def _repr_(self):
        r"""
        Representation method of self.

        Apply the function str to _repr_type(_repr_options) if _repr_type is
        callable and _repr_type else.

        TESTS::

            sage: from surface_dynamics.all import *

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

            sage: p._repr_type = 'interval_diagram'
            sage: p._repr_options = (False,)
            sage: p   #indirect doctest
            [['a', 'b', 'c'], ['a', 'b', 'c']]
            sage: p._repr_options = (True,)
            sage: p
            [['b', ('c', 'a')], ['b', ('c', 'a')]]

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

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: q = iet.Permutation('a b c', 'c b a', flips=['a'])
            sage: r = iet.Permutation('a b c', 'c b a', reduced=True)
            sage: s = iet.Permutation('a b c', 'c b a', flips=['a'], reduced=True)
            sage: len(set([hash(p), hash(q), hash(r), hash(s)]))
            4
        """
        return hash(tuple(self._twin[0])) ^ \
               hash(tuple(self._twin[1])) ^ \
               (hash(None if self._labels is None else tuple(self._labels[0])) | \
               hash(None if self._labels is None else tuple(self._labels[1]))) ^ \
               (hash(None if self._flips is None else tuple(self._flips[0])) | \
               hash(None if self._flips is None else tuple(self._flips[1])))

    def str(self, sep= "\n"):
        r"""
        A string representation of the generalized permutation.

        INPUT:

        - ``sep`` - (default: '\n') a separator for the two intervals

        OUTPUT:

        string -- the string that represents the permutation


        EXAMPLES::

            sage: from surface_dynamics.all import *

        For permutations of iet::

            sage: p = iet.Permutation('a b c','c b a')
            sage: p.str()
            'a b c\nc b a'
            sage: p.str(sep=' | ')
            'a b c | c b a'

        ..the permutation can be rebuilt from the standard string::

            sage: p == iet.Permutation(p.str())
            True

        For permutations of li::

            sage: p = iet.GeneralizedPermutation('a b b','c c a')
            sage: p.str()
            'a b b\nc c a'
            sage: p.str(sep=' | ')
            'a b b | c c a'

        ..the generalized permutation can be rebuilt from the standard string::

            sage: p == iet.GeneralizedPermutation(p.str())
            True
        """
        l = self.list()
        s = []
        s.append(' '.join(map(str,l[0])))
        s.append(' '.join(map(str,l[1])))
        return sep.join(s)
 
    def __copy__(self) :
        r"""
        Returns a copy of self.

        EXAMPLES::

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *

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
        return (len(self._twin[0]) + len(self._twin[1])) / 2

    def length_top(self):
        r"""
        Returns the number of intervals in the top segment.

        OUTPUT:

        integer -- the length of the top segment

        EXAMPLES::

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *

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
            raise ValueError("Your alphabet has not enough letters")
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

            sage: from surface_dynamics.all import *

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


    def left_right_inverse(self):
        r"""
        Returns the left-right inverse.

        The left-right inverse of a permutation, is the permutation obtained by
        reversing the order of the underlying ordering.

        You can also use the shorter .lr_inverse()

        There are two other symmetries of permutation which are accessible via
        the methods
        :meth:`Permutation.top_bottom_inverse` and
        :meth:`Permutation.symmetric`.

        OUTPUT: a permutation

        EXAMPLES::

            sage: from surface_dynamics.all import *

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
        res = copy(self)
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

    def top_bottom_inverse(self):
        r"""
        Returns the top-bottom inverse.

        You can use also use the shorter .tb_inverse().

        There are two other symmetries of permutation which are accessible via
        the methods
        :meth:`Permutation.left_right_inverse` and
        :meth:`Permutation.symmetric`.

        OUTPUT: a permutation

        EXAMPLES::

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *

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
        res = copy(self)
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

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *

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

    def interval_diagram(self, glue_ends=True, sign=False):
        r"""
        Return the interval diagram of self

        Convention, the first letter is always the left hand side.

        INPUT:

        - ``glue_ends`` - bool (default: True)

        - ``sign`` - bool (default: False)

        EXAMPLES::

            sage: from surface_dynamics.all import *

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

            sage: p = iet.GeneralizedPermutation((0,1,0,2),(3,2,4,1,4,3))
            sage: p.interval_diagram()
            [[0, 1, 4, 1], [2, 3, 4, (2, 3, 0)]]
        """
        if self._flips is not None:
            raise ValueError("not implemented for flipped permutation")
        twin = self.twin_list()
        labels = self.list()
        letters = set((label,j) for label in self.letters() for j in xrange(2))
        label_to_twins = dict((label,[]) for label in self.letters())
               # position of each label
        twin_to_index = [] # list of lists of {0,1} of the same size as self.
                           # Each pairs of letters is exactly mapped to a 0 and
                           # a 1 in a canonic way which is compatible with
                           # label_to_twins.
        for i in xrange(2):
            twin_to_index.append([])
            line = labels[i]
            for p in xrange(len(line)):
                twin_to_index[-1].append(len(label_to_twins[line[p]]))
                label_to_twins[line[p]].append((i,p))
        m0 = len(labels[0])
        m1 = len(labels[1])

        singularities = []

        just_glued = False     # True iff the last elt in singularity is paired
        glued_at_begin = False # True iff the 1st elt in singularity is paired
        while letters:
            label,j = letters.pop()
            i,p = label_to_twins[label][j]
            if sign:
                singularity = [(label,j)]
            else:
                singularity = [label]
            while True:
                i,p = twin[i][p]
                if i == 0:
                    p += 1
                    if p == m0:
                        i = 1
                        p = m1-1
                        j = 1
                else:
                    p-=1
                    if p == -1:
                        i = 0
                        p = 0
                        j = 0
                label = labels[i][p]
                j = twin_to_index[i][p]
                if (label,j) not in letters:
                    if (glue_ends and
                        ((i == 1 and p == m1-1) or (i == 0 and p == 0))):
                        sg2 = singularity.pop(0)
                        sg1 = singularity.pop(-1)
                        if glued_at_begin:
                            singularity.append((sg1,) + sg2)
                        elif just_glued:
                            singularity.append(sg1 +(sg2,))
                        else:
                            singularity.append((sg1,sg2))
                    break
                letters.remove((label,j))
                if (glue_ends and (
                    (i == 1 and p==m1-1) or (i == 0 and p == 0))):
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

    def cover(self, permutations, as_tuple=False):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p = iet.Permutation('a b', 'b a')
            sage: p.cover(['(1,2)', '(1,3)'])
            Covering of degree 3 of the permutation:
            a b
            b a

            sage: p.cover([[1,0,2], [2,1,0]], as_tuple=True)
            Covering of degree 3 of the permutation:
            a b
            b a
        """
        from cover import PermutationCover

        if len(permutations) != len(self):
            raise ValueError

        if not as_tuple:
            from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
            permutations = [PermutationGroupElement(p,check=True) for p in permutations]
            permutations = [[i-1 for i in p.domain()] for p in permutations]

            d = max(len(p) for p in permutations)

            for p in permutations: p.extend(xrange(len(p),d))

        else:
            d = len(permutations[0])

        return PermutationCover(self, d, permutations)

    def orientation_cover(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: c = p.orientation_cover()
            sage: c
            Covering of degree 2 of the permutation:
            a a b
            b c c
            sage: c.stratum()
            H_1(0^4)
        """
        rank = self.alphabet().rank
        p0 = set(map(rank, self[0]))
        p1 = set(map(rank, self[1]))
        inv_letters = p0.symmetric_difference(p1)
        permut_cover = [[1,0] if i in inv_letters else [0,1] for i in range(len(self.alphabet()))]
        return self.cover(permut_cover, as_tuple=True)

class PermutationIET(Permutation):
    def _init_twin(self, a):
        r"""
        Initializes the twin list.

        EXAMPLES::

            sage: from surface_dynamics.all import *

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

    def twin_list(self):
        r"""
        Returns the twin list of self.

        The twin list is the involution without fixed point associated to that
        permutation seen as two lines of symbols. As the domain is two lines,
        the position are 2-tuples `(i,j)` where `i` specifies the line and `j`
        the position in the line.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c','c b a')
            sage: p.twin_list()[0]
            [(1, 2), (1, 1), (1, 0)]
            sage: p.twin_list()[1]
            [(0, 2), (0, 1), (0, 0)]

        We may check that it is actually an involution without fixed point::

            sage: t = p.twin_list()
            sage: all(t[i][j] != (i,j) for i in xrange(2) for j in xrange(len(t[i])))
            True
            sage: all(t[t[i][j][0]][t[i][j][1]] == (i,j) for i in xrange(2) for j in xrange(len(t[i])))
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

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *

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

    def _move(self, interval, position, position_to):
        r"""
        Moves the element at (interval,position) to (interval, position_to)

        INPUT:

        - ``interval`` - 0 or 1

        - ``position`` - a position in interval

        - ``position_to`` - a position in interval

        TESTS::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c d','d c b a')
            sage: p._move(0,0,2)
            sage: p
            b a c d
            d c b a
            sage: p._move(0,1,3)
            sage: p
            b c a d
            d c b a
            sage: p._move(0,2,4)
            sage: p
            b c d a
            d c b a
            sage: p._move(1,3,1)
            sage: p
            b c d a
            d a c b
        """
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

    def _remove_interval(self, pos):
        r"""
        Remove the interval at ``pos`` in the top interval

        TESTS::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c d', 'd a b c')
            sage: p._remove_interval(1)
            sage: p
            a c d
            d a c
            sage: p._twin
            [[1, 2, 0], [2, 0, 1]]
            sage: p._remove_interval(2)
            sage: p
            a c
            a c
            sage: p._twin
            [[0, 1], [0, 1]]
            sage: p._remove_interval(0)
            sage: p
            c
            c
            sage: p._twin
            [[0], [0]]

            sage: p = iet.Permutation('a b c d', 'd c b a', reduced=True)
            sage: p._remove_interval(1)
            sage: p
            a b c
            c b a
            sage: p._remove_interval(2)
            sage: p
            a b
            b a

            sage: p = iet.Permutation('a b c d', 'd a c b', flips=['a','b'])
            sage: p._remove_interval(0)
            sage: p
            -b  c  d
             d  c -b
            sage: p._remove_interval(2)
            sage: p
            -b  c
             c -b
        """
        assert 0 <= pos < len(self)

        twin = self._twin[0][pos]

        if self._flips is not None:
            del self._flips[0][pos]
            del self._flips[1][twin]

        if self._labels is not None:
            del self._labels[0][pos]
            del self._labels[1][twin]

        del self._twin[0][pos]
        del self._twin[1][twin]

        for i,j in enumerate(self._twin[0]):
            if j > twin:
                self._twin[0][i] = j-1
        for i,j in enumerate(self._twin[1]):
            if j > pos:
                self._twin[1][i] = j-1

    def _identify_intervals(self, pos1, pos2):
        r"""
        Identify the interval ``pos1`` in top with ``pos2`` in bottom.

        The label of the top interval at position ``pos1`` disappears.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c', 'c a b')
            sage: p._identify_intervals(2,2)
            sage: p
            a b
            b a

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: p._identify_intervals(0,0)
            sage: p
            b c
            b c
            sage: p.letters()
            ['b', 'c']

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: p._identify_intervals(0,1)
            sage: p
            b c
            c b

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: p._identify_intervals(1,0)
            sage: p
            a c
            c a

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: p._identify_intervals(1,2)
            sage: p
            a c
            c a 

            sage: p = iet.Permutation('a b c d', 'c a d b')
            sage: p._identify_intervals(3,1)
            sage: p
            a b c 
            c a b

            sage: p = iet.Permutation('a b c d', 'd a c b')
            sage: p._identify_intervals(1,0)
            sage: p
            a c d
            a c d
            sage: p._twin
            [[0, 1, 2], [0, 1, 2]]
            sage: p = iet.Permutation('a b c d', 'd a c b')
            sage: p._identify_intervals(1,1)
            sage: p
            a c d
            d c a
            sage: p._twin
            [[2, 1, 0], [2, 1, 0]]

            sage: p = iet.Permutation('a b c d', 'd a c b')
            sage: p._identify_intervals(1,2)
            sage: p
            a c d
            d a c

            sage: p = iet.Permutation('a b c d', 'd a c b')
            sage: p._identify_intervals(1,3)
            Traceback (most recent call last):
            ...
            AssertionError
        """
        assert 0 <= pos1 < len(self)
        assert 0 <= pos2 < len(self)
        twin1 = self._twin[0][pos1]
        assert pos2 != twin1
        if self._flips:
            raise ValueError("not implemented if there is a flip")
        self._remove_interval(pos1)
        self._move(1, pos2 if pos2 < twin1 else pos2-1, twin1)

    def has_rauzy_move(self, winner, side='right'):
        r"""
        Test if a Rauzy move can be performed on this permutation.

        EXAMPLES::

            sage: from surface_dynamics.all import *

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

    def letters(self):
        r"""
        Returns the list of letters of the alphabet used for representation.

        The letters used are not necessarily the whole alphabet (for example if
        the alphabet is infinite).

        OUTPUT: a list of labels

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation([1,2],[2,1])
            sage: p.alphabet(Alphabet(name="NN"))
            sage: p
            0 1
            1 0
            sage: p.letters()
            [0, 1]
        """
        if self._labels is not None:
            return map(self._alphabet.unrank, sorted(self._labels[0]))
        else:
            return map(self._alphabet.unrank, range(len(self)))


class PermutationLI(Permutation):
    def _init_twin(self, a):
        r"""
        Initializes the _twin attribute

        TESTS::

            sage: from surface_dynamics.all import *

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

            for i in (0,1):
              for j in range(len(twin[i])) :
                  if twin[i][j] == (i,j) :
                    if a[i][j] in a[i][j+1:] :
                    # two up or two down
                        j2 = (a[i][j+1:]).index(a[i][j]) + j + 1
                        twin[i][j] = (i,j2)
                        twin[i][j2] = (i,j)
                    else :
                        # one up, one down (here i=0)
                        j2 = a[1].index(a[i][j])
                        twin[0][j] = (1,j2)
                        twin[1][j2] = (0,j)

            self._twin = twin

    def _init_alphabet(self, intervals) :
        r"""
        Intialization procedure of the alphabet of self from intervals list

        TEST::

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

            sage: from surface_dynamics.all import *

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
            for j in xrange(self.length(interval)):
                self._twin[interval][j] = (
                    1-self._twin[interval][j][0],
                    self._twin[interval][j][1])

    def _reversed_twin(self):
        r"""
        Reverses the twin of the permutation.

        EXAMPLES::

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *

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

    def letters(self):
        r"""
        Returns the list of letters of the alphabet used for representation.

        The letters used are not necessarily the whole alphabet (for example if
        the alphabet is infinite).

        OUTPUT: a list of labels

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation([1,2],[2,1])
            sage: p.alphabet(Alphabet(name="NN"))
            sage: p
            0 1
            1 0
            sage: p.letters()
            [0, 1]
        """
        unrank = self._alphabet.unrank
        if self._labels is not None:
            return map(unrank, sorted(set(self._labels[0]+self._labels[1])))
        else:
            return map(unrank, range(len(self)))

    def twin_list(self):
        r"""
        Returns the twin list of self.

        The twin list is the involution without fixed point which defines it. As the
        domain is naturally split into two lines we use a 2-tuple (i,j) to
        specify the element at position j in line i.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('a a b','b c c')
            sage: p.twin_list()[0]
            [(0, 1), (0, 0), (1, 0)]
            sage: p.twin_list()[1]
            [(0, 2), (1, 2), (1, 1)]

        And we may check that it is actually an involution without fixed point::

            sage: t = p.twin_list()
            sage: all(t[i][j] != (i,j) for i in xrange(2) for j in xrange(len(t[i])))
            True
            sage: all(t[t[i][j][0]][t[i][j][1]] == (i,j) for i in xrange(2) for j in xrange(len(t[i])))
            True

        A slightly more complicated example::

            sage: q = iet.GeneralizedPermutation('a b c a','d e f e g c b g d f')
            sage: q.twin_list()[0]
            [(0, 3), (1, 6), (1, 5), (0, 0)]
            sage: q.twin_list()[1]
            [(1, 8), (1, 3), (1, 9), (1, 1), (1, 7), (0, 2), (0, 1), (1, 4), (1, 0), (1, 2)]

        ::

            sage: t = q.twin_list()
            sage: all(t[t[i][j][0]][t[i][j][1]] == (i,j) for i in xrange(2) for j in xrange(len(t[i])))
            True
        """
        return [self._twin[0][:],self._twin[1][:]]

    def to_cylindric(self):
        r"""
        Return a cylindric permutation in the same extended Rauzy class

        A generalized permutation is *cylindric* if the first letter in the top
        interval is the same as the last letter in the bottom interval.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('a b d a c','c e b e d')
            sage: p.is_irreducible()
            True
            sage: p.to_cylindric().is_cylindric()
            True

        ALGORITHM:

        The algorithm is naive. It computes the extended Rauzy class until it
        finds a cylindric permutation.
        """
        wait = []
        rauzy_class = set([self])
        q = self
        while True:
            if q.has_rauzy_move('t'): # top rauzy move
                qq = q.rauzy_move('t')
                if qq not in rauzy_class:
                    if qq._twin[1][-1] == (0,0):
                        return qq
                    wait.append(qq)
                    rauzy_class.add(qq)
            if q.has_rauzy_move('b'): # bot rauzy move
                qq = q.rauzy_move('b')
                if qq not in rauzy_class:
                    if qq._twin[1][-1] == (0,0):
                        return qq
                    wait.append(qq)
                    rauzy_class.add(qq)
            qq = q.symmetric() # symmetric
            if qq not in rauzy_class:
                if qq._twin[1][-1] == (0,0):
                    return qq
                wait.append(qq)
                rauzy_class.add(qq)
            q = wait.pop()

        raise ValueError, "no cylindric permutation in the extended Rauzy class"

    def is_cylindric(self):
        r"""
        Test if the permutation is cylindric

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: q = iet.GeneralizedPermutation('a b b','c c a')
            sage: q.is_cylindric()
            True
            sage: q = iet.GeneralizedPermutation('a a b b','c c')
            sage: q.is_cylindric()
            False
        """
        return self._twin[0][-1] == (1,0) or self._twin[1][-1] == (0,0)

    def stratum(self):
        r"""
        Returns the stratum associated to self

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('a b b','c c a')
            sage: p.stratum()
            Q_0(-1^4)
        """
        from surface_dynamics.flat_surfaces.quadratic_strata import QuadraticStratum
        if self.is_irreducible():
            return QuadraticStratum([x-2 for x in self.profile()])
        raise ValueError("stratum is well defined only for irreducible permutations")

    def profile(self):
        r"""
        Returns the ``profile`` of self.

        The *profile* of a generalized permutation is the list `(d_1, \ldots,
        d_k)` where `(d_1 \pi, \ldots, d_k \pi)` is the list of angles of any
        suspension of that generalized permutation.

        See also :meth:`marked_profile`.

        EXAMPLES::

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *

            sage: iet.GeneralizedPermutation('a a b','b c c').genus()
            0
            sage: iet.GeneralizedPermutation((0,1,2,1,3),(4,3,4,2,0)).genus()
            2
        """
        p = self.profile()
        return Integer((sum(p)-2*len(p))/4+1)

    def marking(self):
        r"""
        Return the marking induced by the two sides of the interval

        EXAMPLES::

            sage: from surface_dynamics.all import *

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
        singularities in the the suspension. The partition, called the
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

            sage: from surface_dynamics.all import *

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
        from marked_partition import MarkedPartition

        if len(self) == 1:
            return MarkedPartition([],2,(0,0))

        g = self.interval_diagram(glue_ends=True,sign=True)
        p = sorted(map(lambda x: len(x), g),reverse=True)

        if self._twin[1][0][0] == 0:
            left1 = ((self[1][0],0),(self[0][0],0))
        else:
            left1 = ((self[1][0],1),(self[0][0],0))
        left2 = (left1[1],left1[0])
        if self._twin[0][-1][0] == 1:
            right1 = ((self[0][-1],1),(self[1][-1],1))
        else:
            right1 = ((self[0][-1],0),(self[1][-1],1))
        right2 = (right1[1],right1[0])
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

            sage: from surface_dynamics.all import *

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
            NotImplementedError: Not yet implemented! Do it!
        """
        if self.stratum().nb_fake_zeros():
            raise NotImplementedError, "Not yet implemented! Do it!"
        else:
            return self

    def is_hyperelliptic(self,verbose=False):
        r"""
        Test if this permutation is in an hyperelliptic connected component.

        EXAMPLES::

            sage: from surface_dynamics.all import *

        An example of hyperelliptic permutation::

            sage: p = iet.GeneralizedPermutation([0,1,2,0,6,5,3,1,2,3],[4,5,6,4])
            sage: p.is_hyperelliptic()
            True

        Check for the corresondance::

            sage: q = QuadraticStratum(6,6)
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

            sage: q = QuadraticStratum(3,3,2)
            sage: c_hyp, c_non_hyp = q.components()
            sage: p_hyp = c_hyp.permutation_representative()
            sage: p_hyp.is_hyperelliptic()
            True
            sage: p_non_hyp = c_non_hyp.permutation_representative()
            sage: p_non_hyp.is_hyperelliptic()
            False
            sage: q = QuadraticStratum(5,5,2)
            sage: c_hyp, c_non_hyp = q.components()
            sage: p_hyp = c_hyp.permutation_representative()
            sage: p_hyp.is_hyperelliptic()
            True
            sage: p_non_hyp = c_non_hyp.permutation_representative()
            sage: p_non_hyp.is_hyperelliptic()
            False
            sage: q = QuadraticStratum(3,3,1,1)
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
        zeros = s.zeros()

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
            print "found Jenkins-Strebel"
            print q
            print l0
            print l1

        if any(x[0] == 1 for x in l0):
            if verbose: print "potential form 1"
            i0 = []; i1 = []
            for i in xrange(len(l0)):
                if l0[i][0] == 0:
                    i0.append(i)
            for i in xrange(len(l1)):
                if l1[i][0] == 1:
                    i1.append(i)
            if len(i0) != 2 or len(i1) != 2:
                if verbose: print "no repetition twice in intervals"
                return False

            q0_0 = q0[i0[0]+1:i0[1]]
            q0_1 = q0[i0[1]+1:] + q0[:i0[0]]
            q0_0.reverse()
            q0_1.reverse()

            q1_0 = q1[i1[0]+1:i1[1]]
            q1_1 = q1[i1[1]+1:] + q1[:i1[0]]

            if verbose:
                print q0_0, q0_1
                print q1_0, q1_1

            return (q0_0 == q1_0 and q0_1 == q1_1) or (q0_0 == q1_1 and q0_1 == q1_0)

        else:
            if verbose: print "potential form 2"
            if any(i==1 for i,_ in l0) or any(i==0 for i,_ in l1):
                return False
            j = len(l0) // 2
            for i in xrange(j):
                if l0[i][1] != j+i:
                    return False

            j = len(l1) // 2
            for i in xrange(j):
                if l1[i][1] != j+i:
                    return False

            return True

    def stratum_component(self):
        r"""
        Return the connected component of stratum in which self belongs to.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('a b b','c c a')
            sage: p.stratum_component()
            Q_0(-1^4)^c

        Test the exceptionnal strata in genus 3::

            sage: Q = QuadraticStratum(9,-1)
            sage: p = Q.regular_component().permutation_representative()
            sage: p.stratum_component()
            Q_3(9, -1)^reg
            sage: p = Q.irregular_component().permutation_representative()
            sage: p.stratum_component()
            Q_3(9, -1)^irr

            sage: Q = QuadraticStratum(6,3,-1)
            sage: p = Q.regular_component().permutation_representative()
            sage: p.stratum_component()
            Q_3(6, 3, -1)^reg
            sage: p = Q.irregular_component().permutation_representative()
            sage: p.stratum_component()
            Q_3(6, 3, -1)^irr

            sage: Q = QuadraticStratum(3,3,3,-1)
            sage: p = Q.regular_component().permutation_representative()
            sage: p.stratum_component()
            Q_3(3^3, -1)^reg
            sage: p = Q.irregular_component().permutation_representative()
            sage: p.stratum_component()
            Q_3(3^3, -1)^irr

        Test the exceptionnal strata in genus 4::

            sage: Q = QuadraticStratum(12)
            sage: p = Q.regular_component().permutation_representative()
            sage: p.stratum_component()
            Q_4(12)^reg
            sage: p = Q.irregular_component().permutation_representative()
            sage: p.stratum_component()
            Q_4(12)^irr

            sage: Q = QuadraticStratum(9,3)
            sage: p = Q.regular_component().permutation_representative()
            sage: p.stratum_component()  # long time - 1.5sec
            Q_4(9, 3)^reg
            sage: p = Q.irregular_component().permutation_representative()
            sage: p.stratum_component()  # long time - 2sec
            Q_4(9, 3)^irr

            sage: Q = QuadraticStratum(6,6)
            sage: p = Q.hyperelliptic_component().permutation_representative()
            sage: p.stratum_component()
            Q_4(6^2)^hyp
            sage: p = Q.regular_component().permutation_representative()
            sage: p.stratum_component()  # long time - 1sec
            Q_4(6^2)^reg
            sage: p = Q.irregular_component().permutation_representative()
            sage: p.stratum_component()  # long time - 1sec
            Q_4(6^2)^irr

            sage: Q = QuadraticStratum(6,3,3)
            sage: p = Q.regular_component().permutation_representative()
            sage: p.stratum_component()  # long time - 3sec
            Q_4(6, 3^2)^reg
            sage: p = Q.irregular_component().permutation_representative()
            sage: p.stratum_component()  # long time - 3sec
            Q_4(6, 3^2)^irr

            sage: Q = QuadraticStratum(3,3,3,3)
            sage: p = Q.hyperelliptic_component().permutation_representative()
            sage: p.stratum_component()
            Q_4(3^4)^hyp
            sage: p = Q.regular_component().permutation_representative()
            sage: p.stratum_component()  # long time - 5sec
            Q_4(3^4)^reg
            sage: p = Q.irregular_component().permutation_representative()
            sage: p.stratum_component()  # long time - 5sec
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
            raise NotImplementedError, "database of irregular twins not available"
        p = self.erase_marked_points().to_cylindric()

        if p._twin in D.get(stratum):
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

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *

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

    def is_irreducible(self, return_decomposition=False) :
        r"""
        Tests the irreducibility.

        An abelian permutation p = (p0,p1) is reducible if:
        set(p0[:i]) = set(p1[:i]) for an i < len(p0)

        OUTPUT:

        - a boolean

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: p.is_irreducible()
            True

            sage: p = iet.Permutation('a b c', 'b a c')
            sage: p.is_irreducible()
            False
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

    #TODO: change the name
    def decompose(self):
        r"""
        Returns the decomposition as a concatenation of irreducible permutations.

        OUTPUT:

        a list of permutations

        EXAMPLES::

            sage: from surface_dynamics.all import *

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

    def intersection_matrix(self,ring=None):
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

            sage: from surface_dynamics.all import *

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
        """
        if ring is None:
            from sage.rings.integer_ring import ZZ
            ring = ZZ
        n = self.length_top()
        m = matrix(ring,n)
        for i in range(n):
            for j in range(i,n):
                if self._twin[0][i] > self._twin[0][j]:
                    m[i,j] = 1
                    m[j,i] = -1
        return m

    def attached_out_degree(self):
        r"""
        Returns the degree of the singularity at the left of the interval.

        OUTPUT:

        - a positive integer


        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p1 = iet.Permutation('a b c d e f g','d c g f e b a')
            sage: p2 = iet.Permutation('a b c d e f g','e d c g f b a')
            sage: p1.attached_out_degree()
            3
            sage: p2.attached_out_degree()
            1
        """
        left_corner = (self[1][0],0)
        for s in self.interval_diagram(glue_ends=False,sign=True):
            if left_corner in s:
                return len(s)//2 - 1

    def attached_in_degree(self):
        r"""
        Returns the degree of the singularity at the right of the interval.

        OUTPUT:

        - a positive integer


        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p1 = iet.Permutation('a b c d e f g','d c g f e b a')
            sage: p2 = iet.Permutation('a b c d e f g','e d c g f b a')
            sage: p1.attached_in_degree()
            1
            sage: p2.attached_in_degree()
            3
        """
        right_corner = (self[1][-1],1)

        for s in self.interval_diagram(glue_ends=False,sign=True):
            if right_corner in s:
                return len(s)//2 - 1

    def profile(self):
        r"""
        Returns the profile of the permutation

        EXAMPLES::

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *

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
        associated to the angles of conical singularities in the the suspension
        together with a data associated to the endpoint called marking.

        If the left endpoint and the right endpoint of the interval associated
        to the permutation, then the marking is of type one and consists in a
        couple ``(m,a)`` such that ``m`` is the angle of the conical singularity
        and ``a`` is the angle between the outgoing separatrix associated to the
        left endpoint and the incoming separatrix associated to the right
        endpoint. A marking of type one is denoted ``(m|a)``.

        If the left endpoint and the right endpoint are two different conical
        singularities in the suspension the the marking is of type two and
        consists in a couple ``(m_l,m_r)`` where ``m_l`` (resp. ``m_r``) is the
        conical angle of the singularity at the left endpoint (resp. right
        endpoint). A marking of type two is denoted ``m_l o m_r``

        EXAMPLES::

            sage: from surface_dynamics.all import *

        The irreducible permutation on 1 interval has marked profile of type 2
        with data `(0,0)`::

            sage: p = iet.Permutation('a','a')
            sage: p.marked_profile()
            0o0 []

        Permutations in H(3,1) with all possible profiles::

            sage: p = iet.Permutation('a b c d e f g','b g a c f e d')
            sage: p.interval_diagram()
            [['a', ('b', 'a'), ('g', 'd'), 'e', 'f', 'g', 'b', 'c'], ['f', 'c', 'd', 'e']]
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
            [['a', 'b', ('g', 'd'), ('e', 'a'), 'b', 'c', 'e', 'f'], ['f', 'g', 'c', 'd']]
            sage: p.marked_profile()
            4|3 [4, 2]

            sage: p = iet.Permutation('a b c d e f g', 'f d c a g e b')
            sage: p.interval_diagram()
            [['a', 'b', 'e', ('f', 'a'), 'c', 'd', 'f', 'g'], [('g', 'b'), 'c', 'd', 'e']]
            sage: p.marked_profile()
            4o2 [4, 2]
        """
        from marked_partition import MarkedPartition

        if len(self) == 1:
            return MarkedPartition([],2,(0,0))

        g = self.interval_diagram(glue_ends=True,sign=True)
        p = sorted(map(lambda x: len(x)//2, g),reverse=True)
        left = ((self[1][0],0),(self[0][0],0))
        right = ((self[0][-1],1),(self[1][-1],1))
        for c in g:
            if left in c: c_left = c
            if right in c: c_right = c

        if c_left == c_right:
            mm = len(c_left)
            a = ((c_right.index(right)-c_left.index(left)-1) %mm) // 2
            return MarkedPartition(p, 1, (mm//2, a))

        else:
            m_l = len(c_left) // 2
            m_r = len(c_right) //2
            return MarkedPartition(p, 2, (m_l,m_r))

    def stratum(self):
        r"""
        Returns the strata in which any suspension of this permutation lives.

        OUTPUT:

        - a stratum of Abelian differentials

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: print p.stratum()
            H_1(0^2)

            sage: p = iet.Permutation('a b c d', 'd a b c')
            sage: print p.stratum()
            H_1(0^3)

            sage: p = iet.Permutation(range(9), [8,5,2,7,4,1,6,3,0])
            sage: print p.stratum()
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
        from surface_dynamics.flat_surfaces.abelian_strata import AbelianStratum

        if not self.is_irreducible():
            return map(lambda x: x.stratum(), self.decompose())

        if len(self) == 1:
            return AbelianStratum([])

        singularities = [x - 1 for x in self.profile()]

        return AbelianStratum(singularities)

    def genus(self) :
        r"""
        Returns the genus corresponding to any suspension of self.

        The genus can be deduced from the profile (see :meth:`profile`)
        `p = (p_1,\ldots,p_k)` of self by the formula:
        `2g-2 = \sum_{i=1}^k (p_i-1)`.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: p.genus()
            1

            sage: p = iet.Permutation('a b c d','d c b a')
            sage: p.genus()
            2

        REFERENCES:

            Veech, 1982
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

            sage: from surface_dynamics.all import *

        Permutations from the odd and even component of H(2,2,2)::

            sage: a = range(10)
            sage: b1 = [3,2,4,6,5,7,9,8,1,0]
            sage: b0 = [6,5,4,3,2,7,9,8,1,0]
            sage: p1 = iet.Permutation(a,b1)
            sage: print p1.arf_invariant()
            1
            sage: p0 = iet.Permutation(a,b0)
            sage: print p0.arf_invariant()
            0

        Permutations from the odd and even component of H(4,4)::

            sage: a = range(11)
            sage: b1 = [3,2,5,4,6,8,7,10,9,1,0]
            sage: b0 = [5,4,3,2,6,8,7,10,9,1,0]
            sage: p1 = iet.Permutation(a,b1)
            sage: print p1.arf_invariant()
            1
            sage: p0 = iet.Permutation(a,b0)
            sage: print p0.arf_invariant()
            0

        REFERENCES:

        [Jo80] D. Johnson, "Spin structures and quadratic forms on surfaces", J.
        London Math. Soc (2), 22, 1980, 365-373

        [KoZo03] M. Kontsevich, A. Zorich "Connected components of the moduli
        spaces of Abelian differentials with prescribed singularities",
        Inventiones Mathematicae, 153, 2003, 631-678
        """
        if any((z+1)%2 for z in self.profile()):
            return None

        from sage.rings.finite_rings.constructor import GF
        GF2 = GF(2)

        M = self.intersection_matrix(GF2)
        F, C = M.symplectic_form()

        g = F.rank()/2
        n = F.ncols()

        s = GF2(0)
        for i in range(g):
            a = C.row(i)

            a_indices = [k for k in xrange(n) if a[k]]
            t_a = GF2(len(a_indices))
            for j1 in xrange(len(a_indices)):
                for j2 in xrange(j1+1,len(a_indices)):
                    t_a += M[a_indices[j1], a_indices[j2]]

            b = C.row(g+i)
            b_indices = [k for k in xrange(n) if b[k]]
            t_b = GF2(len(b_indices))
            for j1 in xrange(len(b_indices)):
                for j2 in xrange(j1+1,len(b_indices)):
                    t_b += M[b_indices[j1],b_indices[j2]]

            s += t_a * t_b

        return s

    def stratum_component(self):
        r"""
        Returns a connected components of a stratum.

        EXAMPLES::

            sage: from surface_dynamics.all import *

        Permutations from the stratum H(6)::

            sage: a = range(8)
            sage: b_hyp = [7,6,5,4,3,2,1,0]
            sage: b_odd = [3,2,5,4,7,6,1,0]
            sage: b_even = [5,4,3,2,7,6,1,0]
            sage: p_hyp = iet.Permutation(a, b_hyp)
            sage: p_odd = iet.Permutation(a, b_odd)
            sage: p_even = iet.Permutation(a, b_even)
            sage: print p_hyp.stratum_component()
            H_4(6)^hyp
            sage: print p_odd.stratum_component()
            H_4(6)^odd
            sage: print p_even.stratum_component()
            H_4(6)^even

        Permutations from the stratum H(4,4)::

            sage: a = range(11)
            sage: b_hyp = [10,9,8,7,6,5,4,3,2,1,0]
            sage: b_odd = [3,2,5,4,6,8,7,10,9,1,0]
            sage: b_even = [5,4,3,2,6,8,7,10,9,1,0]
            sage: p_hyp = iet.Permutation(a,b_hyp)
            sage: p_odd = iet.Permutation(a,b_odd)
            sage: p_even = iet.Permutation(a,b_even)
            sage: p_hyp.stratum() == AbelianStratum(4,4)
            True
            sage: print p_hyp.stratum_component()
            H_5(4^2)^hyp
            sage: p_odd.stratum() == AbelianStratum(4,4)
            True
            sage: print p_odd.stratum_component()
            H_5(4^2)^odd
            sage: p_even.stratum() == AbelianStratum(4,4)
            True
            sage: print p_even.stratum_component()
            H_5(4^2)^even

        As for stratum you can specify that you want to attach the singularity
        on the left of the interval using the option marked_separatrix::

            sage: a = range(1,10)
            sage: b_odd = [4,3,6,5,7,9,8,2,1]
            sage: b_even = [6,5,4,3,7,9,8,2,1]
            sage: p_odd = iet.Permutation(a,b_odd)
            sage: p_even = iet.Permutation(a,b_even)
            sage: p_odd.stratum_component()
            H_4(4, 2)^odd
            sage: p_even.stratum_component()
            H_4(4, 2)^even
        """
        from surface_dynamics.flat_surfaces.abelian_strata import (ASC, HypASC, NonHypASC, OddASC, EvenASC)

        if not self.is_irreducible():
            return map(lambda x: x.stratum_component(), self.decompose())

        stratum = self.stratum()
        cc = stratum._cc

        if len(cc) == 1:
            return stratum.components()[0]

        if HypASC in cc:
            if self.is_hyperelliptic():
                return HypASC(stratum)
            else:
                cc = cc[1:]

        if len(cc) == 1:
            return cc[0](stratum)

        else:
            spin = self.arf_invariant()
            if spin == 0:
                return EvenASC(stratum)
            else:
                return OddASC(stratum)

    def order_of_rauzy_action(self, winner, side=None):
        r"""
        Returns the order of the action of a Rauzy move.

        INPUT:

        - ``winner`` - string ``'top'`` or ``'bottom'``

        - ``side`` - string ``'left'`` or ``'right'``

        OUTPUT:

        An integer corresponding to the order of the Rauzy action.

        EXAMPLES::

            sage: from surface_dynamics.all import *

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

    def rauzy_move(self, winner, side='right'):
        r"""
        Returns the permutation after a Rauzy move.

        INPUT:

        - ``winner`` - 'top' or 'bottom' interval

        - ``side`` - 'right' or 'left' (defaut: 'right') corresponding
          to the side on which the Rauzy move must be performed.

        OUTPUT:

        - a permutation

        EXAMPLES::

            sage: from surface_dynamics.all import *

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
        """
        winner = interval_conversion(winner)
        side = side_conversion(side)
        loser = 1 - winner

        if self._twin[0][side] == side or self._twin[0][side] == len(self)+side:
            raise ValueError("Rauzy induction is not well defined")

        res = copy(self)

        wtp = res._twin[winner][side]

        if side == -1:
            res._move(loser, len(self._twin[loser])-1, wtp+1)

        if side == 0:
            res._move(loser, 0, wtp)

        return res

    def backward_rauzy_move(self, winner, side='right'):
        r"""
        Returns the permutation before a Rauzy move.

        INPUT:

        - ``winner`` - 'top' or 'bottom' interval

        - ``side`` - 'right' or 'left' (defaut: 'right') corresponding
          to the side on which the Rauzy move must be performed.

        OUTPUT:

        - a permutation

        TESTS::

            sage: from surface_dynamics.all import *

        Testing the inversion on labelled permutations::

            sage: p = iet.Permutation('a b c d','d c b a')
            sage: for pos,side in [('t','r'),('b','r'),('t','l'),('b','l')]:
            ...    q = p.rauzy_move(pos,side)
            ...    print q.backward_rauzy_move(pos,side) == p,
            ...    q = p.backward_rauzy_move(pos,side)
            ...    print q.rauzy_move(pos,side) == p,
            True True True True True True True True

        Testing the inversion on reduced permutations::

            sage: p = iet.Permutation('a b c d','d c b a',reduced=True)
            sage: for pos,side in [('t','r'),('b','r'),('t','l'),('b','l')]:
            ...    q = p.rauzy_move(pos,side)
            ...    print q.backward_rauzy_move(pos,side) == p,
            ...    q = p.backward_rauzy_move(pos,side)
            ...    print q.rauzy_move(pos,side) == p,
            True True True True True True True True
        """
        winner = interval_conversion(winner)
        side = side_conversion(side)

        loser = 1 - winner
        winner_twin = self._twin[winner][side]
        d = len(self)

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
            if self._labels is not None:
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

            sage: from surface_dynamics.all import *

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
            raise ValueError, "the permutation must be irreducible"

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
        for i in xrange(len(self)):
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

            sage: from surface_dynamics.all import *

            sage: iet.Permutation('a b c d','d c b a').is_hyperelliptic()
            True
            sage: iet.Permutation('0 1 2 3 4 5','5 2 1 4 3 0').is_hyperelliptic()
            False

        REFERENCES:

        Gerard Rauzy, "Echanges d'intervalles et transformations induites",
        Acta Arith. 34, no. 3, 203-212, 1980

        M. Kontsevich, A. Zorich "Connected components of the moduli space
        of Abelian differentials with prescripebd singularities" Invent. math.
        153, 631-678 (2003)
        """
        test = self.erase_marked_points()

        n = test.length_top()
        cylindric = test.to_standard()
        return cylindric._twin[0] == range(n-1,-1,-1)

    def to_cylindric(self):
        r"""
        Returns a cylindric permutation in the same Rauzy class.

        A permutation is *cylindric* if the first letter in the top interval is
        also the last letter of the bottom interval or if the last letter of the
        top interval is the first letter of the bottom interval.

        TESTS::

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *

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
        Returns the permutation as an element of the symetric group.

        EXAMPLES::

            sage: from surface_dynamics.all import *

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
        return Permutation(map(lambda x: x+1,self._twin[1]))

class OrientablePermutationLI(PermutationLI):
    r"""
    Template for quadratic permutation.

    .. warning::

        Internal class! Do not use directly!

    AUTHOR:

    - Vincent Delecroix (2008-12-20): initial version

    """
    def rauzy_move(self, winner, side=-1):
        r"""
        Returns the permutation after a Rauzy move.

        TESTS::

            sage: from surface_dynamics.all import *

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
        """
        winner = interval_conversion(winner)
        side = side_conversion(side)
        loser = 1 - winner

        res = copy(self)

        wti, wtp = res._twin[winner][side]

        if side == -1:
            d = len(self._twin[loser])
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

            sage: from surface_dynamics.all import *

        Tests the inversion on labelled generalized permutations::

            sage: p = iet.GeneralizedPermutation('a a b b','c c d d')
            sage: for pos,side in [('t','r'),('b','r'),('t','l'),('b','l')]:
            ...    q = p.rauzy_move(pos,side)
            ...    print q.backward_rauzy_move(pos,side) == p,
            ...    q = p.backward_rauzy_move(pos,side)
            ...    print q.rauzy_move(pos,side) == p,
            True True True True True True True True


        Tests the inversion on reduced generalized permutations::

            sage: p = iet.GeneralizedPermutation('a a b b','c c d d',reduced=True)
            sage: for pos,side in [('t','r'),('b','r'),('t','l'),('b','l')]:
            ...    q = p.rauzy_move(pos,side)
            ...    print q.backward_rauzy_move(pos,side) == p,
            ...    q = p.backward_rauzy_move(pos,side)
            ...    print q.rauzy_move(pos,side) == p,
            True True True True True True True True
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
        right or on the left. The definition is due to [BL08]_ where they prove
        that the property of being irreducible is stable under Rauzy induction.

        INPUT:

        -  ``return_decomposition`` - boolean (default: False) - if True, and
           the permutation is reducible, returns also the blocs A1 u B1, B1 u
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

            sage: from surface_dynamics.all import *

            sage: GP = iet.GeneralizedPermutation

            sage: GP('a a','b b').is_irreducible()
            False
            sage: GP('a a b','b c c').is_irreducible()
            True
            sage: GP('1 2 3 4 5 1','5 6 6 4 3 2').is_irreducible()
            True

        TESTS::

            sage: from surface_dynamics.all import *

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

        AUTHORS:

        - Vincent Delecroix (2008-12-20)
        """
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

            for i21 in xrange(0, l1) :
                if i21 > 0 and s1[i21-1] in A21:
                    break
                A21 = s1[:i21]

                for i12 in xrange(l0 - 1, i11 - 1, -1) :
                    if s0[i12] in A12 or s0[i12] in A21:
                        break
                    A12 = s0[i12:]

                    for i22 in xrange(l1 - 1, i21 - 1, -1) :
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

class FlippedPermutation(Permutation):
    r"""
    Template for flipped generalized permutations.

    .. warning::

        Internal class! Do not use directly!

    AUTHORS:

    - Vincent Delecroix (2008-12-20): initial version

    """
    def _init_flips(self,intervals,flips):
        r"""
        Initialize the flip list

        TESTS::

            sage: from surface_dynamics.all import *

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

    def str(self, sep="\n"):
        r"""
        String representation.

        TESTS::

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('a a','b b',flips='a')
            sage: print p.str()
            -a -a
             b  b
             sage: print p.str('/')
             -a -a/ b  b
        """
        l = self.list(flips=True)
        return (' '.join(map(labelize_flip, l[0]))
                + sep
                + ' '.join(map(labelize_flip, l[1])))

        return s


class FlippedPermutationIET(FlippedPermutation, PermutationIET):
    r"""
    Template for flipped Abelian permutations.

    .. warning::

        Internal class! Do not use directly!
    """
    def rauzy_move(self, winner, side=-1):
        r"""
        Returns the permutation after a Rauzy move.

        TESTS::

            sage: from surface_dynamics.all import *

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
        """
        winner = interval_conversion(winner)
        side = side_conversion(side)
        loser = 1 - winner

        res = copy(self)

        wtp = res._twin[winner][side]
        flip = self._flips[winner][side]
        if flip == -1:
            res._flips[loser][side] *= -1
            res._flips[winner][res._twin[loser][side]] *= -1
            flip = 1
        else:
            flip = 0

        if side == -1:
            d = len(self._twin[loser])
            res._move(loser, d-1, wtp+1-flip)

        if side == 0:
            res._move(loser, 0, wtp+flip)

        return res

    def backward_rauzy_move(self, winner, side=-1):
        r"""
        Returns the permutation before a Rauzy move.

        TESTS::

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('a b c d e','d a b e c',flips='abcd')
            sage: for pos,side in [('t','r'),('b','r'),('t','l'),('b','l')]:
            ...    q = p.rauzy_move(pos,side)
            ...    print q.backward_rauzy_move(pos,side) == p,
            ...    q = p.backward_rauzy_move(pos,side)
            ...    print q.rauzy_move(pos,side) == p,
            True True True True True True True True

        Testing the inversion on reduced permutations::

            sage: p = iet.Permutation('f a b c d e','d f c b e a',flips='abcd', reduced=True)
            sage: for pos,side in [('t','r'),('b','r'),('t','l'),('b','l')]:
            ...    q = p.rauzy_move(pos,side)
            ...    print q.backward_rauzy_move(pos,side) == p,
            ...    q = p.backward_rauzy_move(pos,side)
            ...    print q.rauzy_move(pos,side) == p,
            True True True True True True True True
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
                res._move(loser, wtp-1, d)
            else:
                res._move(loser, wtp+1, d)

        if side == 0:
            if flip == -1:
                res._flips[loser][wtp+1] *= -1
                res._flips[winner][res._twin[loser][wtp+1]] *= -1
                res._move(loser, wtp+1, 0)
            else:
                res._move(loser, wtp-1, 0)

        return res

    def flips(self):
        r"""
        Returns the list of flips.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c','c b a',flips='ac')
            sage: p.flips()
            ['a', 'c']
        """
        result = []
        l = self.list(flips=False)
        for i,f in enumerate(self._flips[0]):
            if f == -1:
                result.append(l[0][i])
        return result



class FlippedPermutationLI(FlippedPermutation, PermutationLI):
    r"""
    Template for flipped quadratic permutations.

    .. warning::

        Internal class! Do not use directly!

    AUTHORS:

    - Vincent Delecroix (2008-12-20): initial version

    """
    def flips(self):
        r"""
        Returns the list of flipped intervals.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('a a','b b',flips='a')
            sage: p.flips()
            ['a']
            sage: p = iet.GeneralizedPermutation('a a','b b',flips='b',reduced=True)
            sage: p.flips()
            ['b']
        """
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

    def rauzy_move(self, winner, side=-1):
        r"""
        Rauzy move

        TESTS::

            sage: from surface_dynamics.all import *

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
        """
        winner = interval_conversion(winner)
        side = side_conversion(side)
        loser = 1 - winner

        res = copy(self)

        wti, wtp = res._twin[winner][side]
        flip = self._flips[winner][side]
        if flip == -1:
            res._flips[loser][side] *= -1
            lti,ltp = res._twin[loser][side]
            res._flips[lti][ltp] *= -1
            flip = 1
        else:
            flip = 0

        if side == -1:
            d = len(self._twin[loser])
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

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('a b c e b','d c d a e',flips='abcd')
            sage: for pos,side in [('t','r'),('b','r'),('t','l'),('b','l')]:
            ...    q = p.rauzy_move(pos,side)
            ...    print q.backward_rauzy_move(pos,side) == p,
            ...    q = p.backward_rauzy_move(pos,side)
            ...    print q.rauzy_move(pos,side) == p,
            True True True True True True True True

        Testing the inversion on reduced permutations::

            sage: p = iet.GeneralizedPermutation('a b c e b','d c d a e',flips='abcd',reduced=True)
            sage: for pos,side in [('t','r'),('b','r'),('t','l'),('b','l')]:
            ...    q = p.rauzy_move(pos,side)
            ...    print q.backward_rauzy_move(pos,side) == p,
            ...    q = p.backward_rauzy_move(pos,side)
            ...    print q.rauzy_move(pos,side) == p,
            True True True True True True True True
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
            lengths. This correspondance is obtained via the Rauzy induction. To a
            idoc IET we can associate a unique path in a Rauzy diagram. This
            establishes a correspondance between infinite full path in Rauzy diagram
            and equivalence topologic class of IET.
        """
        def __init__(self, parent, *data):
            r"""
            Constructor of the path.

            TEST::

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

            TEST::

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

                sage: from surface_dynamics.all import *

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

                sage: from surface_dynamics.all import *

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

                sage: from surface_dynamics.all import *

                
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

            TEST::

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
                isinstance(self, type(other)) and
                self._parent == other._parent and
                self._start == other._start and
                self._edge_types == other._edge_types)

        def __ne__(self,other):
            r"""
            Tests inequality

            TEST::

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
            return (
                not isinstance(self, type(other)) or
                self._parent != other._parent or
                self._start != other._start or
                self._edge_types != other._edge_types)

        def __copy__(self):
            r"""
            Returns a copy of the path.

            TESTS::

                sage: from surface_dynamics.all import *

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

                sage: from surface_dynamics.all import *

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

                sage: from surface_dynamics.all import *

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
                raise ValueError, "Edge type not valid"

            if self._parent._succ[self._end][edge_type] is None:
                raise ValueError, "%d is not a valid edge" %(edge_type)

            self._edge_types.append(edge_type)
            self._end = self._parent._succ[self._end][edge_type]

        def _fast_append(self, edge_type):
            r"""
            Append an edge to the path without verification.

            EXAMPLES::

                sage: from surface_dynamics.all import *

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

                sage: from surface_dynamics.all import *

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
                raise ValueError, "Not on the same Rauzy diagram"

            if self._end != path._start:
                raise ValueError, "The end of the first path must the start of the second"

            self._edge_types.extend(path._edge_types)
            self._end = path._end

        def _fast_extend(self, path):
            r"""
            Extension with no verification.

            EXAMPLES::

                sage: from surface_dynamics.all import *

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

            TEST::


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

                sage: from surface_dynamics.all import *

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
                raise IndexError, "path index out of range"

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

                sage: from surface_dynamics.all import *

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
                raise ValueError, "The end of the first path is not the start of the second"

            res = copy(self)
            res._fast_extend(other)
            return res

        def __mul__(self, n):
            r"""
            Multiple of a loop.

            EXAMPLES::

                sage: from surface_dynamics.all import *

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

                sage: from surface_dynamics.all import *

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

                sage: from surface_dynamics.all import *

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

                sage: from surface_dynamics.all import *

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

                sage: from surface_dynamics.all import *

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: g = r.path(p)
                sage: for q in g:
                ...       print p
                a b c
                c b a
                sage: g = r.path(p, 't', 't')
                sage: for q in g:
                ...       print q, "\n*****"
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
                ...       print q, "\n*****"
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

                sage: from surface_dynamics.all import *

                sage: p = iet.Permutation('a b','b a')
                sage: r = p.rauzy_diagram()
                sage: def f(i,t):
                ...       if t is None: return []
                ...       return [t]
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

            - ``composition`` - the composition function for the function. * if None (defaut None)

            TEST::

                sage: p = iet.Permutation('a b','b a')
                sage: r = p.rauzy_diagram()
                sage: def f(i,t):
                ...       if t is None: return []
                ...       return [t]
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

            sage: from surface_dynamics.all import *

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
                raise ValueError, "right_induction can not be empty string"

            elif 'top'.startswith(right_induction):
                self._index['rt_rauzy'] = len(self._edge_types)
                self._edge_types.append(('rauzy_move',(0,-1)))

            elif 'bottom'.startswith(right_induction):
                self._index['rb_rauzy'] = len(self._edge_types)
                self._edge_types.append(('rauzy_move',(1,-1)))

            else:
                raise ValueError, "%s is not valid for right_induction" %(right_induction)

        if left_induction is True:
            self._index['lt_rauzy'] = len(self._edge_types)
            self._edge_types.append(('rauzy_move',(0,0)))
            self._index['lb_rauzy'] = len(self._edge_types)
            self._edge_types.append(('rauzy_move',(1,0)))

        elif isinstance(left_induction,str):
            if left_induction == '':
                raise ValueError, "left_induction can not be empty string"

            elif 'top'.startswith(left_induction):
                self._index['lt_rauzy'] = len(self._edge_types)
                self._edge_types.append(('rauzy_move',(0,0)))

            elif 'bottom'.startswith(left_induction):
                self._index['lb_rauzy'] = len(self._edge_types)
                self._edge_types.append(('rauzy_move',(1,0)))

            else:
                raise ValueError, "%s is not valid for left_induction" %(right_induction)

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

            sage: from surface_dynamics.all import *

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
            ...       p.rauzy_diagram() == r
            True
            True
            True
            True
            True
            True
            True
        """
        return (
            isinstance(self, type(other)) and
            self._edge_types == other._edge_types and
            self._succ.keys()[0] in other._succ)

    def __ne__(self, other):
        r"""
        Tests difference.


        TEST::

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
        return (
            not isinstance(self, type(other)) or
            self._edge_types != other._edge_types or
            self._succ.keys()[0] not in other._succ)

    def vertices(self):
        r"""
        Returns a list of the vertices.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: for p in r.vertices(): print p
            a b
            b a
        """
        return map(
            lambda x: self._vertex_to_permutation(x),
            self._succ.keys())

    def vertex_iterator(self):
        r"""
        Returns an iterator over the vertices

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: for p in r.vertex_iterator(): print p
            a b
            b a

        ::

            sage: r = iet.RauzyDiagram('a b c d','d c b a')
            sage: from itertools import ifilter
            sage: r_1n = ifilter(lambda x: x.is_standard(), r)
            sage: for p in r_1n: print p
            a b c d
            d c b a
        """
        from itertools import imap
        return imap(
            lambda x: self._vertex_to_permutation(x),
            self._succ.keys())

    def edges(self,labels=True):
        r"""
        Returns a list of the edges.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: len(r.edges())
            2
        """
        return list(self.edge_iterator())

    def edge_iterator(self):
        r"""
        Returns an iterator over the edges of the graph.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b','b a')
            sage: r = p.rauzy_diagram()
            sage: for e in r.edge_iterator():
            ...    print e[0].str(sep='/'), '-->', e[1].str(sep='/')
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

            sage: from surface_dynamics.all import *

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
            raise ValueError, "the edge type must be a string"

        if ('top_rauzy_move'.startswith(data) or
            't_rauzy_move'.startswith(data)):
            if self._index.has_key('lt_rauzy'):
                if self._index.has_key('rt_rauzy'):
                    raise ValueError, "left and right inductions must be differentiated"
                return self._index['lt_rauzy']

            if self._index.has_key('rt_rauzy'):
                return self._index['rt_rauzy']

            raise ValueError, "no top induction in this Rauzy diagram"

        if ('bottom_rauzy_move'.startswith(data) or
            'b_rauzy_move'.startswith(data)):
            if self._index.has_key('lb_rauzy'):
                if self._index.has_key('rb_rauzy'):
                    raise ValueError, "left and right inductions must be differentiated"
                return self._index['lb_rauzy']

            if self._index.has_key('rb_rauzy'):
                return self._index['rb_rauzy']

            raise ValueError, "no bottom Rauzy induction in this diagram"

        if ('left_rauzy_move'.startswith(data) or
            'l_rauzy_move'.startswith(data)):
            if self._index.has_key('lt_rauzy'):
                if self._index.has_key('lb_rauzy'):
                    raise ValueError, "top and bottom inductions must be differentiated"
                return self._index['lt_rauzy']

            if self._index.has_key('lb_rauzy'):
                return self._index('lb_rauzy')

            raise ValueError, "no left Rauzy induction in this diagram"

        if ('lt_rauzy_move'.startswith(data) or
            'tl_rauzy_move'.startswith(data) or
            'left_top_rauzy_move'.startswith(data) or
            'top_left_rauzy_move'.startswith(data)):
            if not self._index.has_key('lt_rauzy'):
                raise ValueError, "no top-left Rauzy induction in this diagram"
            else:
                return self._index['lt_rauzy']

        if ('lb_rauzy_move'.startswith(data) or
            'bl_rauzy_move'.startswith(data) or
            'left_bottom_rauzy_move'.startswith(data) or
            'bottom_left_rauzy_move'.startswith(data)):
            if not self._index.has_key('lb_rauzy'):
                raise ValueError, "no bottom-left Rauzy induction in this diagram"
            else:
                return self._index['lb_rauzy']

        if 'right'.startswith(data):
            raise ValueError, "ambiguity with your edge name: %s" %(data)

        if ('rt_rauzy_move'.startswith(data) or
            'tr_rauzy_move'.startswith(data) or
            'right_top_rauzy_move'.startswith(data) or
            'top_right_rauzy_move'.startswith(data)):
            if not self._index.has_key('rt_rauzy'):
                raise ValueError, "no top-right Rauzy induction in this diagram"
            else:
                return self._index['rt_rauzy']

        if ('rb_rauzy_move'.startswith(data) or
            'br_rauzy_move'.startswith(data) or
            'right_bottom_rauzy_move'.startswith(data) or
            'bottom_right_rauzy_move'.startswith(data)):
            if not self._index.has_key('rb_rauzy'):
                raise ValueError, "no bottom-right Rauzy induction in this diagram"
            else:
                return self._index['rb_rauzy']

        if 'symmetric'.startswith(data):
            if not self._index.has_key('symmetric'):
                raise ValueError, "no symmetric in this diagram"
            else:
                return self._index['symmetric']

        if 'inversion'.startswith(data) or data == 'inverse':
            if self._index.has_key('lr_inverse'):
                if self._index.has_key('tb_inverse'):
                    raise ValueError, "left-right and top-bottom inversions must be differentiated"
                return self._index['lr_inverse']

            if self._index.has_key('tb_inverse'):
                return self._index['tb_inverse']

            raise ValueError, "no inversion in this diagram"

        if ('lr_inversion'.startswith(data) or
            data == 'lr_inverse' or
            'left_right_inversion'.startswith(data) or
            data == 'left_right_inverse'):
            if not self._index.has_key('lr_inverse'):
                raise ValueError, "no left-right inversion in this diagram"
            else:
                return self._index['lr_inverse']

        if ('tb_inversion'.startswith(data) or
            data == 'tb_inverse' or
            'top_bottom_inversion'.startswith(data)
            or data == 'top_bottom_inverse'):
            if not self._index.has_key('tb_inverse'):
                raise ValueError, "no top-bottom inversion in this diagram"
            else:
                return self._index['tb_inverse']

        raise ValueError, "this edge type does not exist: %s" %(data)

    def edge_types(self):
        r"""
        Print information about edges.

        EXAMPLES::

            sage: from surface_dynamics.all import *

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
            print str(i) + ": " + edge_type + str(t)

    def alphabet(self, data=None):
        r"""
        TESTS::

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c','c b a')
            sage: d = p.rauzy_diagram()
            sage: print d.edge_to_matrix(p,1)
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

        TEST::

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

        TEST::

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

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b','b a')
            sage: r = p.rauzy_diagram()
            sage: g0 = r.path(p)
            sage: for g in r._all_npath_extension(g0,0):
            ...       print g
            Path of length 0 in a Rauzy diagram
            sage: for g in r._all_npath_extension(g0,1):
            ...       print g
            Path of length 1 in a Rauzy diagram
            Path of length 1 in a Rauzy diagram
            sage: for g in r._all_npath_extension(g0,2):
            ...       print g
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

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b','b a')
            sage: r = p.rauzy_diagram()
            sage: g0 = r.path(p)
            sage: for g in r._all_path_extension(g0,0):
            ...       print g
            Path of length 0 in a Rauzy diagram
            sage: for g in r._all_path_extension(g0, 1):
            ...       print g
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

            sage: from surface_dynamics.all import *

            sage: r = iet.RauzyDiagram('a b','b a')
            sage: for p in r: print p
            a b
            b a
            sage: r = iet.RauzyDiagram('a b c','c b a')
            sage: for p in r: print p.stratum()
            H_1(0^2)
            H_1(0^2)
            H_1(0^2)
        """
        for data in self._succ.iterkeys():
            yield self._vertex_to_permutation(data)

    def __contains__(self, element):
        r"""
        Containance test.

        TESTS::

            sage: from surface_dynamics.all import *

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
        for p in self._succ.iterkeys():
            if self._vertex_to_permutation(p) == element:
                return True

        return False

    def _repr_(self):
        r"""
        Returns a representation of self

        TEST::

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

            sage: from surface_dynamics.all import *


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
            raise ValueError, "Your element does not have the good type"

        perm = self._permutation_to_vertex(p)
        return map(lambda x: self._vertex_to_permutation(x),
            self._succ[perm])

    def __len__(self):
        r"""
        Deprecated use cardinality.

        EXAMPLES::

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *

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

        TEST::

            sage: r = iet.RauzyDiagram('a b c','c b a')   #indirect doctest
            sage: r = iet.RauzyDiagram('a b c','c b a',left_induction=True) #indirect doctest
            sage: r = iet.RauzyDiagram('a b c','c b a',symmetric=True)   #indirect doctest
            sage: r = iet.RauzyDiagram('a b c','c b a',lr_inversion=True)   #indirect doctest
            sage: r = iet.RauzyDiagram('a b c','c b a',tb_inversion=True)   #indirect doctest
        """
        if p.__class__ is not self._element_class:
            raise ValueError, "your permutation is not of good type"

        if len(p) != self._n:
            raise ValueError, "your permutation has not the good length"

        pred = self._pred
        succ = self._succ
        p = self._permutation_to_vertex(p)
        perm = self._element
        l = []

        if not succ.has_key(p):
            succ[p] = [None] * len(self._edge_types)
            pred[p] = [None] * len(self._edge_types)
            l.append(p)

        while(l != []):
            p = l.pop()
            self._set_element(p)

            for t,edge in enumerate(self._edge_types):
                if (not hasattr(perm, 'has_'+edge[0]) or
                  getattr(perm, 'has_'+edge[0])(*(edge[1]))):
                    q = getattr(perm,edge[0])(*(edge[1]))
                    q = self._permutation_to_vertex(q)
                    if not succ.has_key(q):
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

            sage: from surface_dynamics.all import *

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

            sage: from surface_dynamics.all import *

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

        for p,neighbours in self._succ.iteritems():
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

            sage: from surface_dynamics.all import *

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
            raise ValueError, "your permutation is not of good type"

        if len(p) != self._n:
            raise ValueError, "your permutation has not the good length"

        pred = self._pred
        succ = self._succ
        p = self._permutation_to_vertex(p)
        l = []

        if not succ.has_key(p):
            succ[p] = [None] * len(self._edge_types)
            pred[p] = [None] * len(self._edge_types)
            l.append(p)

        while(l != []):
            p = l.pop()

            for t,edge_type in enumerate(self._edge_types):
                perm = self._vertex_to_permutation(p)

                if (not hasattr(perm,'has_' + edge_type[0]) or
                  getattr(perm, 'has_' + edge_type[0])(*(edge_type[1]))):
                    q = perm.rauzy_move(t)
                    q = self._permutation_to_vertex(q)
                    if reducible == True or perm.is_irreducible():
                        if not succ.has_key(q):
                            succ[q] = [None] * len(self._edge_types)
                            pred[q] = [None] * len(self._edge_types)
                            l.append(q)

                        succ[p][t] = q
                        pred[q][t] = p

