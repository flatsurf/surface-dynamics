r"""
Interval Exchange Transformations and Linear Involution

An interval exchage transformation is a map defined on an interval (see
help(iet.IntervalExchangeTransformation) for a more complete help.

EXAMPLES::

    sage: from surface_dynamics.all import *

Initialization of a simple iet with integer lengths::

    sage: T = iet.IntervalExchangeTransformation(Permutation([3,2,1]), [3,1,2])
    sage: print T
    Interval exchange transformation of [0, 6[ with permutation
    1 2 3
    3 2 1

Rotation corresponds to iet with two intervals::

    sage: p = iet.Permutation('a b', 'b a')
    sage: T = iet.IntervalExchangeTransformation(p, [1, (sqrt(5)-1)/2])
    sage: print T.in_which_interval(0)
    a
    sage: print T.in_which_interval(T(0))
    a
    sage: print T.in_which_interval(T(T(0)))
    b
    sage: print T.in_which_interval(T(T(T(0))))
    a

There are two plotting methods for iet::

    sage: p = iet.Permutation('a b c','c b a')
    sage: T = iet.IntervalExchangeTransformation(p, [1, 2, 3])

.. plot the domain and the range of T::

    sage: T.plot_two_intervals()   # not tested (problem with matplotlib font cache)
    Graphics object consisting of 12 graphics primitives

.. plot T as a function::

    sage: T.plot_function()  # not tested (problem with matplotlib font cache)
    Graphics object consisting of 3 graphics primitives
"""
from copy import copy
from sage.modules.free_module_element import free_module_element
from sage.rings.all import ZZ, QQ

from template import side_conversion, interval_conversion


def wedge(u, v):
    r"""
    Return the wedge of the vectors ``u`` and ``v``
    """
    d = len(u)
    R = u.base_ring()
    assert len(u) == len(v) and v.base_ring() == R
    return free_module_element(R, d*(d-1)//2, [(u[i]*v[j] - u[j]*v[i]) for i in range(d-1) for j in range(i+1,d)])

class IntervalExchangeTransformation(object):
    r"""
    Interval exchange transformation

    INPUT:

    - ``permutation`` - a permutation (LabelledPermutationIET)

    - ``lengths`` - the list of lengths

    EXAMPLES::

        sage: from surface_dynamics.all import *

    Direct initialization::

        sage: p = iet.IET(('a b c','c b a'),{'a':1,'b':1,'c':1})
        sage: p.permutation()
        a b c
        c b a
        sage: p.lengths()
        (1, 1, 1)

    Initialization from a iet.Permutation::

        sage: perm = iet.Permutation('a b c','c b a')
        sage: l = vector([0.5,1,1.2])
        sage: t = iet.IET(perm,l)
        sage: t.permutation() == perm
        True
        sage: t.lengths() == l
        True

    Initialization from a Permutation::

        sage: p = Permutation([3,2,1])
        sage: iet.IET(p, [1,1,1])
        Interval exchange transformation of [0, 3[ with permutation
        1 2 3
        3 2 1

    If it is not possible to convert lengths to real values an error is raised::

        sage: iet.IntervalExchangeTransformation(('a b','b a'),['e','f'])
        Traceback (most recent call last):
        ...
        TypeError: unable to convert x (='e') into a real number

    The value for the lengths must be positive::

        sage: iet.IET(('a b','b a'),[-1,-1])
        Traceback (most recent call last):
        ...
        ValueError: lengths must be positive
    """
    def __init__(self, permutation=None, lengths=None, base_ring=None):
        r"""
        INPUT:

        - ``permutation`` - a permutation (LabelledPermutationIET)

        - ``lengths`` - the list of lengths

        TEST::

            sage: from surface_dynamics.all import *

            sage: p = iet.IntervalExchangeTransformation(('a','a'), [1])
            sage: p == loads(dumps(p))
            True
        """
        from labelled import LabelledPermutationIET
        from sage.modules.free_module_element import free_module_element as vector

        if permutation is None or lengths is None:
            self._permutation = None
            self._lengths = []
            self._base_ring = None
            self._vector_space = None
        else:
            self._permutation = copy(permutation)
            self._lengths = vector(lengths)
            V = self._lengths.parent()
            self._base_ring = self._lengths.base_ring()
            self._vector_space = V.ambient_module()

    def vector_space(self):
        return self._vector_space

    def base_ring(self):
        r"""
        Return the base ring over which the lengths are defined

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p = iet.Permutation('a b', 'b a')

            sage: T = iet.IntervalExchangeTransformation(p, [3, 12])
            sage: T.base_ring()
            Integer Ring

            sage: T = iet.IntervalExchangeTransformation(p, [3, 12/5])
            sage: T.base_ring()
            Rational Field

            sage: T = iet.IntervalExchangeTransformation(p, [3, AA(12).sqrt()])
            sage: T.base_ring()
            Algebraic Real Field

            sage: sqrt2 = QuadraticField(2).gen()
            sage: T = iet.IntervalExchangeTransformation(p, [1, sqrt2])
            sage: T.base_ring()
            Number Field in a with defining polynomial x^2 - 2
        """
        return self._base_ring

    def permutation(self):
        r"""
        Returns the permutation associated to this iet.

        OUTPUT:

        permutation -- the permutation associated to this iet

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: perm = iet.Permutation('a b c','c b a')
            sage: p = iet.IntervalExchangeTransformation(perm,(1,2,1))
            sage: p.permutation() == perm
            True
        """
        return copy(self._permutation)

    def length(self):
        r"""
        Returns the total length of the interval.

        OUTPUT:

        real -- the length of the interval

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: t = iet.IntervalExchangeTransformation(('a b','b a'),[1,1])
            sage: t.length()
            2
        """
        return sum(self._lengths)

    def lengths(self):
        r"""
        Returns the list of lengths associated to this iet.

        OUTPUT:

        vector -- the list of lengths of subinterval (the order of the entries
                  correspond to the alphabet)

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.IntervalExchangeTransformation(('a b','b a'),[1,3])
            sage: p.lengths()
            (1, 3)
        """
        return self._lengths.__copy__()

    def translations(self):
        r"""
        Return the vector of translations operated on each intervals.

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p = iet.Permutation('a b c', 'c b a')
            sage: T = iet.IntervalExchangeTransformation(p, [5,1,3])
            sage: T.translations()
            (4, -2, -6)

        The order of the entries correspond to the alphabet::

            sage: p = iet.Permutation('a c d b', 'b d c a', alphabet='abcd')
            sage: T = iet.IntervalExchangeTransformation(p, [1, 1, 1, 1])
            sage: T.translations()
            (3, -3, 1, -1)

        This vector is covariant with respect to the Rauzy matrices::

            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: R = p.rauzy_diagram()
            sage: g = R.path(p, *'ttbtbtbtbb')
            sage: T = g.self_similar_iet()
            sage: for i in range(12):
            ....:     S, code = T.zorich_move(iterations=i, data=True)
            ....:     gg = R.path(p, *code)
            ....:     m = gg.matrix()
            ....:     assert m * S.lengths() == T.lengths()
            ....:     assert m.transpose() * T.translations() == S.translations()
        """
        dom_sg = self.domain_singularities()
        im_sg = self.range_singularities()

        p = self._permutation._labels
        top_twin = self._permutation._twin[0]
        top = p[0]
        bot = p[1]

        translations = self.vector_space()()
        for i0,j in enumerate(top):
            i1 = top_twin[i0]
            translations[j] = im_sg[i1] - dom_sg[i0]

        return translations

    def sah_arnoux_fathi_invariant(self):
        r"""
        Return the Sah-Arnoux-Fathi invariant

        The interval exchange needs to be defined over a number field. The output
        is then a vector with rational entries of dimension `d (d-1) / 2` where
        `d` is the degree of the field.

        EXAMPLES:

        The golden rotation::

            sage: from surface_dynamics.all import *
            sage: p = iet.Permutation('a b','b a')
            sage: R = p.rauzy_diagram()
            sage: g = R.path(p, 't', 'b')
            sage: T = g.self_similar_iet()
            sage: T.sah_arnoux_fathi_invariant()
            (2)

        The Sah-Arnoux-Fathi invariant is not changed under Rauzy
        (or Zorich) induction::

            sage: S = T.zorich_move(iterations=100)
            sage: S.sah_arnoux_fathi_invariant()
            (2)
            sage: (T.length().n(), S.length().n())
            (2.61803398874989, 0.000000000000000)

        An other rotation::

            sage: g = R.path(p, 't', 'b', 'b')
            sage: T = g.self_similar_iet()
            sage: T.sah_arnoux_fathi_invariant()
            (1)
            sage: T.rauzy_move().sah_arnoux_fathi_invariant()
            (1)

        Arnoux-Yoccoz examples in genus 4::

            sage: x = polygen(ZZ)
            sage: poly = x^4 - x^3 - x^2 - x - 1
            sage: l = max(poly.roots(AA, False))
            sage: K.<a> = NumberField(poly, embedding=l)
            sage: top = 'A1l A1r A2 B1 B2 C1 C2 D1 D2'
            sage: bot = 'A1r B2 B1 C2 C1 D2 D1 A2 A1l'
            sage: p = iet.Permutation(top, bot)
            sage: lengths = vector((a**4-a**3, 2*a**3-a**4, a**3, a**2, a**2, a, a, 1, 1))
            sage: T = iet.IntervalExchangeTransformation(p, lengths)
            sage: T.sah_arnoux_fathi_invariant()
            (0, 0, 0, 0, 0, 0)
            sage: T.zorich_move(iterations=10).sah_arnoux_fathi_invariant()
            (0, 0, 0, 0, 0, 0)
        """
        if self.base_ring() is ZZ:
            return free_module_element(ZZ, [])
        elif self.base_ring() is QQ:
            return free_module_element(QQ, [])

        try:
            K, from_V, to_V = self.base_ring().vector_space()
        except (AttributeError, ValueError):
            raise ValueError("the interval exchange needs to be define over a number field")

        return sum(wedge(to_V(u), to_V(v)) for u,v in zip(self.lengths(), self.translations()))

    def normalize(self, total=1, inplace=False):
        r"""
        Returns a interval exchange transformation of normalized lengths.

        The normalization consist in consider a constant homothetic value for
        each lengths in such way that the sum is given (default is 1).

        INPUT:

        - ``total`` - (default: 1) The total length of the interval

        OUTPUT:

        iet -- the normalized iet

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: t = iet.IntervalExchangeTransformation(('a b','b a'), [1,3])
            sage: t.length()
            4
            sage: s = t.normalize(2)
            sage: s.length()
            2
            sage: s.lengths()
            (1/2, 3/2)
        """
        try:
            y = float(total)
        except ValueError:
            raise TypeError("unable to convert x (='%s') into a real number"  %(str(x)))

        if total <= 0:
           raise ValueError("the total length must be positive")

        if inplace:
            res = self
        else:
            res = copy(self)

        res._lengths *= total / res.length()
        if not inplace:
            return res

    def __repr__(self):
        r"""
        A representation string.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: a = iet.IntervalExchangeTransformation(('a','a'),[1])
            sage: a   #indirect doctest
            Interval exchange transformation of [0, 1[ with permutation
            a
            a
        """
        interval = "[0, %s["%self.length()
        s = "Interval exchange transformation of %s "%interval
        s += "with permutation\n%s"%self._permutation
        return s

    def is_identity(self):
        r"""
        Returns True if self is the identity.

        OUTPUT:

        boolean -- the answer

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation("a b","b a")
            sage: q = iet.Permutation("c d","d c")
            sage: s = iet.IET(p, [1,5])
            sage: t = iet.IET(q, [5,1])
            sage: (s*t).is_identity()
            True
            sage: (t*s).is_identity()
            True
        """
        return self._permutation.is_identity()

    def inverse(self):
        r"""
        Returns the inverse iet.

        OUTPUT:

        iet -- the inverse interval exchange transformation

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation("a b","b a")
            sage: s = iet.IET(p, [1,sqrt(2)-1])
            sage: t = s.inverse()
            sage: t.permutation()
            b a
            a b
            sage: t.lengths()
            (1, sqrt(2) - 1)
            sage: t*s
            Interval exchange transformation of [0, sqrt(2)[ with permutation
            aa bb
            aa bb

        We can verify with the method .is_identity()::

            sage: p = iet.Permutation("a b c d","d a c b")
            sage: s = iet.IET(p, [1, sqrt(2), sqrt(3), sqrt(5)])
            sage: (s * s.inverse()).is_identity()
            True
            sage: (s.inverse() * s).is_identity()
            True
        """
        res = copy(self)
        res._permutation = self._permutation.top_bottom_inverse()
        return res

    def __mul__(self, other):
        r"""
        Composition of iet.

        The domain (i.e. the length) of the two iets must be the same). The
        alphabet choosen depends on the permutation.

        TESTS:

        ::

            sage: p = iet.Permutation("a b", "a b")
            sage: t = iet.IET(p, [1,1])
            sage: r = t*t
            sage: r.permutation()
            aa bb
            aa bb
            sage: r.lengths()
            [1, 1]

        ::

            sage: p = iet.Permutation("a b","b a")
            sage: t = iet.IET(p, [1,1])
            sage: r = t*t
            sage: r.permutation()
            ab ba
            ab ba
            sage: r.lengths()
            [1, 1]

        ::

            sage: p = iet.Permutation("1 2 3 4 5","5 4 3 2 1")
            sage: q = iet.Permutation("a b","b a")
            sage: s = iet.IET(p, [1]*5)
            sage: t = iet.IET(q, [1/2, 9/2])
            sage: r = s*t
            sage: r.permutation()
            a5 b1 b2 b3 b4 b5
            b5 a5 b4 b3 b2 b1
            sage: r.lengths()
            [1/2, 1, 1, 1, 1, 1/2]
            sage: r = t*s
            sage: r.permutation()
            1b 2b 3b 4b 5a 5b
            5b 4b 3b 2b 1b 5a
            sage: r.lengths()
            [1, 1, 1, 1, 1/2, 1/2]
            sage: t = iet.IET(q, [3/2, 7/2])
            sage: r = s*t
            sage: r.permutation()
            a4 a5 b1 b2 b3 b4
            a5 b4 a4 b3 b2 b1
            sage: r.lengths()
            [1/2, 1, 1, 1, 1, 1/2]
            sage: t = iet.IET(q, [5/2,5/2])
            sage: r = s*t
            sage: r.permutation()
            a3 a4 a5 b1 b2 b3
            a5 a4 b3 a3 b2 b1
            sage: r = t*s
            sage: r.permutation()
            1b 2b 3a 3b 4a 5a
            3b 2b 1b 5a 4a 3a

        ::

            sage: p = iet.Permutation("a b","b a")
            sage: s = iet.IET(p, [4,2])
            sage: q = iet.Permutation("c d","d c")
            sage: t = iet.IET(q, [3, 3])
            sage: r1 = t * s
            sage: r1.permutation()
            ac ad bc
            ad bc ac
            sage: r1.lengths()
            [1, 3, 2]
            sage: r2 = s * t
            sage: r2.permutation()
            ca cb da
            cb da ca
            sage: r2.lengths()
            [1, 2, 3]
        """
        assert(
            isinstance(other, IntervalExchangeTransformation) and
            self.length() == other.length())

        from labelled import LabelledPermutationIET
        from sage.combinat.words.words import Words

        other_sg = other.range_singularities()[1:]
        self_sg = self.domain_singularities()[1:]

        n_other = len(other._permutation)
        n_self = len(self._permutation)

        interval_other = other._permutation._labels[1]
        interval_self = self._permutation._labels[0]

        d_other = dict([(i,[]) for i in interval_other])
        d_self = dict([(i,[]) for i in interval_self])

        i_other = 0
        i_self = 0

        x = self.base_ring().zero()
        l_lengths = []
        while i_other < n_other and i_self < n_self:
            j_other = interval_other[i_other]
            j_self = interval_self[i_self]

            d_other[j_other].append(j_self)
            d_self[j_self].append(j_other)

            if other_sg[i_other] < self_sg[i_self]:
                l = other_sg[i_other] - x
                x = other_sg[i_other]
                i_other += 1
            elif other_sg[i_other] > self_sg[i_self]:
                l = self_sg[i_self] - x
                x = self_sg[i_self]
                i_self += 1
            else:
                l = self_sg[i_self] - x
                x = self_sg[i_self]
                i_other += 1
                i_self += 1

            l_lengths.append(((j_other,j_self),l))

        alphabet_other = other._permutation.alphabet()
        alphabet_self = self._permutation.alphabet()

        d_lengths = dict(l_lengths)

        l_lengths = []
        top_interval = []
        for i in other._permutation._labels[0]:
            for j in d_other[i]:
                a = alphabet_other.unrank(i)
                b = alphabet_self.unrank(j)
                top_interval.append(str(a)+str(b))
                l_lengths.append(d_lengths[(i,j)])

        bottom_interval = []
        for i in self._permutation._labels[1]:
            for j in d_self[i]:
                a = alphabet_other.unrank(j)
                b = alphabet_self.unrank(i)
                bottom_interval.append(str(a)+str(b))

        p = LabelledPermutationIET((top_interval,bottom_interval))
        return IntervalExchangeTransformation(p,l_lengths)

    def __eq__(self, other):
        r"""
        Tests equality

        TESTS::

            sage: from surface_dynamics.all import *

            sage: t = iet.IntervalExchangeTransformation(('a b','b a'),[1,1])
            sage: t == t
            True
        """
        return (
            isinstance(self, type(other)) and
            self._permutation == other._permutation and
            self._lengths == other._lengths)

    def __ne__(self, other):
        r"""
        Tests difference

        TESTS::

            sage: from surface_dynamics.all import *

            sage: t = iet.IntervalExchangeTransformation(('a b','b a'),[1,1])
            sage: t != t
            False
        """
        return (
            not isinstance(self, type(other)) or
            self._permutation != other._permutation or
            self._lengths != other._lengths)

    def recoding(self, n):
        r"""
        Recode this interval exchange transformation on the words of length
        ``n``.

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p = iet.Permutation('a d c b', 'b c a d', alphabet='abcd')
            sage: T = iet.IntervalExchangeTransformation(p, [119,213,82,33])
            sage: T.recoding(2)
            Interval exchange transformation of [0, 447[ with permutation
            ab db cc cb ba bd bc
            ba bd bc cc cb ab db
            sage: T.recoding(3)
            Interval exchange transformation of [0, 447[ with permutation
            aba abd abc dbc ccb cba bab bdb bcc bcb
            cba aba abd abc dbc bcc bcb ccb bab bdb
        """
        length = len(self._permutation)
        lengths = self._lengths
        A = self.permutation().alphabet()
        unrank = A.unrank

        top = self._permutation._labels[0]
        bot = self._permutation._labels[1]
        bot_twin = self._permutation._twin[1]

        sg_top = self.domain_singularities()[1:-1]  # iterates of the top singularities
        cuts = [[[i], lengths[i]] for i in bot]     # the refined bottom interval

        translations = self.translations()


        for step in range(n-1):
            i = 0
            y = self.base_ring().zero()
            new_sg_top = []
            limits = [0]
            for j,x in enumerate(sg_top):
                while y < x:
                    cuts[i][0].append(top[j])
                    y += cuts[i][1]
                    i += 1

                limits.append(i)
                if y != x:
                    cuts.insert(i, [cuts[i-1][0][:-1], cuts[i-1][1]])
                    cuts[i-1][1] -= y-x
                    cuts[i][0].append(top[j+1])
                    cuts[i][1] = y-x
                    i += 1
            while i < len(cuts):
                cuts[i][0].append(top[j+1])
                i += 1
            limits.append(len(cuts))

            # now we reorder cuts with respect to T according to the cut at
            # limits
            if step != n-2:
                new_cuts = []
                for j in bot_twin:
                    new_cuts.extend(cuts[limits[j]:limits[j+1]])
                cuts = new_cuts

        # build the new interval exchange transformations
        itop = [None]*(max(top)+1)
        for i,j in enumerate(top):
            itop[j] = i
        key = lambda x: [itop[i] for i in x[0]]
        top = sorted(cuts, key=key)
        lengths = [y for x,y in top]
        top = [''.join(str(unrank(i)) for i in x) for x,y in top]
        bot = cuts
        bot = [''.join(str(unrank(i)) for i in x) for x,y in bot]

        from labelled import LabelledPermutationIET
        p = LabelledPermutationIET((top,bot))
        return IntervalExchangeTransformation(p,lengths)

    def in_which_interval(self, x, interval=0):
        r"""
        Returns the letter for which x is in this interval.

        INPUT:

        - ``x`` - a positive number

        - ``interval`` - (default: 'top') 'top' or 'bottom'


        OUTPUT:

        label -- a label corresponding to an interval

        TESTS::

            sage: from surface_dynamics.all import *

            sage: t = iet.IntervalExchangeTransformation(('a b c','c b a'),[1,1,1])
            sage: t.in_which_interval(0)
            'a'
            sage: t.in_which_interval(0.3)
            'a'
            sage: t.in_which_interval(1)
            'b'
            sage: t.in_which_interval(1.9)
            'b'
            sage: t.in_which_interval(2)
            'c'
            sage: t.in_which_interval(2.1)
            'c'
            sage: t.in_which_interval(3)
            Traceback (most recent call last):
            ...
            ValueError: your value does not lie in [0; 3[

        .. and for the bottom interval::

            sage: t.in_which_interval(0,'bottom')
            'c'
            sage: t.in_which_interval(1.2,'bottom')
            'b'
            sage: t.in_which_interval(2.9,'bottom')
            'a'
        """
        interval = interval_conversion(interval)

        if x < 0 or x >= self.length():
            raise ValueError("your value does not lie in [0; {}[".format(self.length()))

        i = 0

        while x >= 0:
            x -= self._lengths[self._permutation._labels[interval][i]]
            i += 1

        i -= 1
        x += self._lengths[self._permutation._labels[interval][i]]

        j = self._permutation._labels[interval][i]
        return self._permutation._alphabet.unrank(j)

    def singularities(self):
        r"""
        The list of singularities of 'T' and 'T^{-1}'.

        OUTPUT:

        list -- two lists of positive numbers which corresponds to extremities
            of subintervals

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: t = iet.IntervalExchangeTransformation(('a b','b a'),[1/2,3/2])
            sage: t.singularities()
            [[0, 1/2, 2], [0, 3/2, 2]]
        """
        return [self.domain_singularities(), self.range_singularities()]

    def domain_singularities(self):
        r"""
        Returns the list of singularities of T

        OUTPUT:

        list -- positive reals that corresponds to singularities in the top
            interval

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: t = iet.IET(("a b","b a"), [1, sqrt(2)])
            sage: t.domain_singularities()
            [0, 1, sqrt(2) + 1]
        """
        l = [self.base_ring().zero()]
        for j in self._permutation._labels[0]:
            l.append(l[-1] + self._lengths[j])
        return l

    def range_singularities(self):
        r"""
        Returns the list of singularities of `T^{-1}`

        OUTPUT:

        list -- real numbers that are singular for `T^{-1}`


        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: t = iet.IET(("a b","b a"), [1, sqrt(2)])
            sage: t.range_singularities()
            [0, sqrt(2), sqrt(2) + 1]
        """
        l = [self.base_ring().zero()]
        for j in self._permutation._labels[1]:
            l.append(l[-1] + self._lengths[j])
        return l

    def __call__(self, value):
        r"""
        Return the image of value by this transformation

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: t = iet.IntervalExchangeTransformation(('a b','b a'),[1/2,3/2])
            sage: t(0)
            3/2
            sage: t(1/2)
            0
            sage: t(1)
            1/2
            sage: t(3/2)
            1
        """
        assert(value >= 0 and value < self.length())

        dom_sg = self.domain_singularities()
        im_sg = self.range_singularities()

        a = self.in_which_interval(value)

        i0 = self._permutation[0].index(a)
        i1 = self._permutation[1].index(a)

        return value - dom_sg[i0] + im_sg[i1]

    def rauzy_move(self, side='right', iterations=1, data=False, error_on_saddles=True):
        r"""
        Performs a Rauzy move.

        INPUT:

        - ``side`` - 'left' (or 'l' or 0) or 'right' (or 'r' or 1)

        - ``iterations`` - integer (default :1) the number of iteration of Rauzy
           moves to perform

        - ``data`` - whether to return also the paths and composition of towers

        - ``error_on_saddles`` - (default: ``True``) whether to stop when a saddle
          is encountered

        OUTPUT:

        - ``iet`` -- the Rauzy move of self

        - ``path`` -- (if ``data=True``) a list of 't' and 'b'

        - ``towers`` -- (if ``data=True``) the towers of the Rauzy induction as a word morphism

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: phi = QQbar((sqrt(5)-1)/2)
            sage: t1 = iet.IntervalExchangeTransformation(('a b','b a'),[1,phi])
            sage: t2 = t1.rauzy_move().normalize(t1.length())
            sage: l2 = t2.lengths()
            sage: l1 = t1.lengths()
            sage: l2[0] == l1[1] and l2[1] == l1[0]
            True

            sage: tt,path,sub = t1.rauzy_move(iterations=3, data=True)
            sage: tt
            Interval exchange transformation of [0, 0.3819660112501051?[ with
            permutation
            a b
            b a
            sage: path
            ['b', 't', 'b']
            sage: sub
            WordMorphism: a->aab, b->aabab

        The substitution can also be recovered from the Rauzy diagram::

            sage: p = t1.permutation()
            sage: p.rauzy_diagram().path(p, *path).substitution() == sub
            True

        An other examples involving 3 intervals::

            sage: t = iet.IntervalExchangeTransformation(('a b c','c b a'),[1,1,3])
            sage: t
            Interval exchange transformation of [0, 5[ with permutation
            a b c
            c b a
            sage: t1 = t.rauzy_move()
            sage: t1
            Interval exchange transformation of [0, 4[ with permutation
            a b c
            c a b
            sage: t2 = t1.rauzy_move()
            sage: t2
            Interval exchange transformation of [0, 3[ with permutation
            a b c
            c b a
            sage: t2.rauzy_move()
            Traceback (most recent call last):
            ...
            ValueError: saddle connection found
            sage: t2.rauzy_move(error_on_saddles=False)
            Interval exchange transformation of [0, 2[ with permutation
            a b
            a b

        Degenerate cases::

            sage: p = iet.Permutation('a b', 'b a')
            sage: T = iet.IntervalExchangeTransformation(p, [1,1])
            sage: T.rauzy_move(error_on_saddles=False)
            Interval exchange transformation of [0, 1[ with permutation
            a
            a
        """
        if data:
            towers = {a:[a] for a in self._permutation.letters()}
            path = []

        side = side_conversion(side)

        res = copy(self)
        for i in range(iterations):
            winner,(a,b,c) = res._rauzy_move(side, error_on_saddles=error_on_saddles)
            if data:
                if winner is None:
                    raise ValueError("does not handle degenerate situations")
                towers[a] = towers[b] + towers[c]
                path.append(winner)

        if data:
            from sage.combinat.words.morphism import WordMorphism
            return res, path, WordMorphism(towers)
        else:
            return res

    def zorich_move(self, side='right', iterations=1, data=False):
        r"""
        Performs a Rauzy move.

        INPUT:

        - ``side`` - 'left' (or 'l' or 0) or 'right' (or 'r' or 1)

        - ``iterations`` - integer (default :1) the number of iteration of Rauzy
           moves to perform

        - ``data`` - whether to return also the path

        OUTPUT:

        - ``iet`` -- the Rauzy move of self

        - ``path`` -- (if ``data=True``) a list of 't' and 'b'

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: T = iet.IntervalExchangeTransformation(p, [12, 35, 67])
            sage: T.zorich_move()
            Interval exchange transformation of [0, 55[ with permutation
            a b c
            c a b
            sage: assert T.permutation() == p and T.lengths() == vector((12,35,67))

        A self similar example in genus 2::

            sage: p = iet.Permutation('a b c d', 'd a c b')
            sage: R = p.rauzy_diagram()
            sage: code = 'b'*4 + 't'*1 + 'b'*3 + 't'*1 + 'b'*3 + 't'*1 + 'b'*1 + 't'*1 + 'b'*4 + 't'*1 + 'b'*2 + 't'*7
            sage: g = R.path(p, *code)
            sage: m = g.matrix()
            sage: poly = m.charpoly()
            sage: l = max(poly.roots(AA, False))
            sage: K.<a> = NumberField(poly, embedding=l)
            sage: lengths = (m - a).right_kernel().basis()[0]
            sage: T = iet.IntervalExchangeTransformation(p, lengths)
            sage: T.normalize(a, inplace=True)
            sage: T
            Interval exchange transformation of [0, a[ with permutation
            a b c d
            d a c b
            sage: T2, path = T.zorich_move(iterations=12, data=True)
            sage: a*T2.lengths() == T.lengths()
            True
            sage: path == code
            True

        Saddle connection detection::

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: T = iet.IntervalExchangeTransformation(p, [41, 22, 135])
            sage: T.zorich_move(iterations=100)
            Traceback (most recent call last):
            ...
            ValueError: saddle connection found
            sage: p = iet.Permutation('a b c d e f', 'f c b e d a')
            sage: T = iet.IntervalExchangeTransformation(p, [41, 132, 22, 135, 55, 333])
            sage: T.zorich_move(iterations=100)
            Traceback (most recent call last):
            ...
            ValueError: saddle connection found
        """
        if data:
            path = ''

        side = side_conversion(side)

        res = copy(self)
        for i in range(iterations):
            winner, m = res._zorich_move(side)
            if data:
                path += winner * m

        if data:
            return res, path
        else:
            return res


    def _rauzy_move(self, side=-1, error_on_saddles=True):
        r"""
        Perform a Rauzy move inplace

        INPUT:

        - ``side`` - must be 0 or -1 (no verification)

        - ``error_on_saddles`` - (default ``True``) whether an error is raised
          when a saddle connections is encountered

        OUTPUT:

        - ``T`` - the iet after Rauzy induction

        - ``winner`` - either ``'t'`` or ``'b'`` (and possibly ``None`` if
          ``error_on_saddles`` is ``False``)

        - ``(a,b,c)`` - a triple of letter such that the towers are obtained by
          applying the substitution `a \mapsto bc`.

        TEST::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: K.<sqrt2> = QuadraticField(2)
            sage: T = iet.IntervalExchangeTransformation(p, [1,1,1,sqrt2])
            sage: T._rauzy_move()
            ('t', ('a', 'a', 'd'))
            sage: T._rauzy_move()
            ('b', ('d', 'b', 'd'))
            sage: T._rauzy_move()
            ('t', ('b', 'b', 'c'))
            sage: T._rauzy_move()
            ('b', ('c', 'b', 'c'))
            sage: T._rauzy_move()
            ('t', ('b', 'b', 'd'))
            sage: T._rauzy_move()
            ('b', ('d', 'c', 'd'))
            sage: T._rauzy_move()
            ('t', ('c', 'c', 'd'))
            sage: T._rauzy_move()
            ('b', ('d', 'a', 'd'))
            sage: T
            Interval exchange transformation of [0, -4*sqrt2 + 7[ with permutation
            a d b c
            d c b a

        Saddle connections::

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: T = iet.IntervalExchangeTransformation(p, [1,2,1])
            sage: T._rauzy_move()
            Traceback (most recent call last):
            ...
            ValueError: saddle connection found
            sage: T._rauzy_move(error_on_saddles=False)
            (None, ('a', 'a', 'c'))
            sage: T.lengths()
            (1, 2, 0)
            sage: T.permutation()
            a b
            a b
        """
        top = self._permutation._labels[0][side]
        bot = self._permutation._labels[1][side]

        unrank = self._permutation.alphabet().unrank
        top_letter = unrank(top)
        bot_letter = unrank(bot)

        length_top = self._lengths[top]
        length_bot = self._lengths[bot]

        if length_top > length_bot:
            winner = 0 # TODO: this value is ignored
            winner_interval = top
            loser_interval = bot
            abc = (bot_letter, bot_letter, top_letter)
            winner = 't'
        elif length_top < length_bot:
            winner = 1 # TODO: this value is ignored
            winner_interval = bot
            loser_interval = top
            abc = (top_letter, bot_letter, top_letter)
            winner = 'b'
        elif error_on_saddles:
            raise ValueError("saddle connection found")
        else:
            # we do a pseudo-top Rauzy induction and remove the bot interval
            # (the iet get one interval less)
            p = self._permutation
            top = p._labels[0][-1]
            bot = p._labels[1][-1]
            p._identify_intervals(-1)
            self._lengths[top] = 0
            return None, (bot_letter,bot_letter,top_letter)

        self._permutation = self._permutation.rauzy_move(winner=winner, side=side, inplace=True)
        self._lengths[winner_interval] -= self._lengths[loser_interval]

        return winner, abc

    def _zorich_move(self, side=-1):
        r"""
        Performs a Zorich move (acceleration of Rauzy) inplace

        INPUT:

        - ``side`` - must be 0 or -1 (no verification)

        OUTPUT:

        - ``winner_letter`` - whether top or bot was done

        - ``m`` -- number of Rauzy made on the same side

        TEST::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b', 'b a')
            sage: K.<sqrt3> = QuadraticField(3)
            sage: t = iet.IntervalExchangeTransformation(p, [1, sqrt3 - 1])
            sage: t._zorich_move()
            ('b', 1)
            sage: t._zorich_move()
            ('t', 2)
            sage: t._zorich_move()
            ('b', 1)
            sage: t._zorich_move()
            ('t', 2)
            sage: t._zorich_move()
            ('b', 1)
            sage: continued_fraction(sqrt3 - 1)
            [0; (1, 2)*]

            sage: x = polygen(ZZ)
            sage: K.<a> = NumberField(x^3 - 2, embedding=AA(2)**(1/3))
            sage: t = iet.IntervalExchangeTransformation(p, [1, a])
            sage: [t._zorich_move()[1] for _ in range(10)]
            [1, 3, 1, 5, 1, 1, 4, 1, 1, 8]
            sage: continued_fraction(a)
            [1; 3, 1, 5, 1, 1, 4, 1, 1, 8, 1, 14, 1, 10, 2, 1, 4, 12, 2, 3, ...]

        A self-similar example in genus 2::

            sage: x = polygen(ZZ)
            sage: poly = x^4 - 11*x^3 + 21*x^2 - 11*x + 1
            sage: l = max(poly.roots(AA, False))
            sage: K.<a> = NumberField(poly, embedding=l)
            sage: lengths = vector([-3*a^3 + 33*a^2 - 63*a + 33, -4*a^3 + 42*a^2 - 63*a + 14,
            ....:     7*a^3 - 74*a^2 + 116*a - 34, -a^2 + 10*a - 10])
            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: t = iet.IntervalExchangeTransformation(p, lengths)
            sage: [t._zorich_move() for _ in range(6)]
            [('t', 1), ('b', 2),  ('t', 2), ('b', 1), ('t', 2), ('b', 2)]
            sage: lengths == a * t.lengths()
            True
            sage: [t._zorich_move(side=0) for _ in range(6)]
            [('b', 1), ('t', 2), ('b', 2), ('t', 1), ('b', 2), ('t', 2)]
            sage: lengths == a^2 * t.lengths()
            True

        Saddle connection detection::

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: T = iet.IntervalExchangeTransformation(p, [41, 1, 5])
            sage: T._zorich_move()
            Traceback (most recent call last):
            ...
            ValueError: saddle connection found
        """
        top = self._permutation._labels[0][side]
        bot = self._permutation._labels[1][side]

        length_top = self._lengths[top]
        length_bot = self._lengths[bot]

        if length_top > length_bot:
            winner = 0
            winner_letter = 't'
            winner_interval = top
            loser_interval = bot
        elif length_top < length_bot:
            winner = 1
            winner_letter = 'b'
            winner_interval = bot
            loser_interval = top
        else:
            raise ValueError("saddle connection found")

        # number of full loops
        lwin = self._lengths[winner_interval]
        loser = 1 - winner
        loser_row = self._permutation._labels[loser]
        llos = 0
        k = 0
        if side == -1:
            # right induction
            j = -1
            while loser_row[j] != winner_interval:
                llos += self._lengths[loser_row[j]]
                j -= 1
                k += 1
        else:
            # left induction
            j = 0
            while loser_row[j] != winner_interval:
                llos += self._lengths[loser_row[j]]
                j += 1
                k += 1

        # remove the full loops
        m = (lwin / llos).floor()
        self._lengths[winner_interval] -= m*llos

        # remaining steps
        r = 0
        while self._lengths[winner_interval] >= self._lengths[loser_interval]:
            self._lengths[winner_interval] -= self._lengths[loser_interval]
            self._permutation.rauzy_move(winner=winner, side=side, inplace=True)
            loser_interval = self._permutation._labels[loser][side]
            r += 1

        if self._lengths[winner_interval].is_zero():
            raise ValueError("saddle connection found")

        return winner_letter, k*m + r

    def __copy__(self):
        r"""
        Returns a copy of this interval exchange transformation.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: t = iet.IntervalExchangeTransformation(('a b','b a'),[1,1])
            sage: s = copy(t)
            sage: s == t
            True
            sage: s is t
            False
        """
        res = self.__class__()
        res._base_ring = self._base_ring
        res._vector_space = self._vector_space
        res._permutation = copy(self._permutation)
        res._lengths = copy(self._lengths)
        return res

    def plot_function(self,**d):
        r"""
        Return a plot of the interval exchange transformation as a
        function.

        INPUT:

        - Any option that is accepted by line2d

        OUTPUT:

        2d plot -- a plot of the iet as a function

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: t = iet.IntervalExchangeTransformation(('a b c d','d a c b'),[1,1,1,1])
            sage: t.plot_function(rgbcolor=(0,1,0))    # not tested (problem with matplotlib font cache)
            Graphics object consisting of 4 graphics primitives
        """
        from sage.plot.plot import Graphics
        from sage.plot.plot import line2d

        G = Graphics()
        l = self.singularities()
        t = self._permutation._twin

        for i in range(len(self._permutation)):
            j = t[0][i]
            G += line2d([(l[0][i],l[1][j]),(l[0][i+1],l[1][j+1])],**d)

        return G

    def plot_two_intervals(
        self,
        position=(0,0),
        vertical_alignment = 'center',
        horizontal_alignment = 'left',
        interval_height=0.1,
        labels_height=0.05,
        fontsize=14,
        labels=True,
        colors=None):
        r"""
        Returns a picture of the interval exchange transformation.

        INPUT:

        - ``position`` - a 2-uple of the position

        - ``horizontal_alignment`` - left (defaut), center or right

        - ``labels`` - boolean (defaut: True)

        - ``fontsize`` - the size of the label


        OUTPUT:

        2d plot -- a plot of the two intervals (domain and range)

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: t = iet.IntervalExchangeTransformation(('a b','b a'),[1,1])
            sage: t.plot_two_intervals()  # not tested (problem with matplotlib font cache)
            Graphics object consisting of 8 graphics primitives
        """
        from sage.plot.plot import Graphics
        from sage.plot.plot import line2d
        from sage.plot.plot import text
        from sage.plot.colors import rainbow

        G = Graphics()

        lengths = map(float,self._lengths)
        total_length = sum(lengths)

        if colors is None:
            colors = rainbow(len(self._permutation), 'rgbtuple')

        if horizontal_alignment == 'left':
            s = position[0]
        elif horizontal_alignment == 'center':
            s = position[0] - total_length / 2
        elif horizontal_alignment == 'right':
            s = position[0] - total_length
        else:
            raise ValueError, "horizontal_alignement must be left, center or right"

        top_height = position[1] + interval_height
        for i in self._permutation._labels[0]:
            G += line2d([(s,top_height),(s+lengths[i],top_height)],
                rgbcolor=colors[i])
            if labels == True:
                G += text(str(self._permutation._alphabet.unrank(i)),
                    (s+float(lengths[i])/2,top_height+labels_height),
                    horizontal_alignment='center',
                    rgbcolor=colors[i],
                    fontsize=fontsize)

            s += lengths[i]

        if horizontal_alignment == 'left':
            s = position[0]
        elif horizontal_alignment == 'center':
            s = position[0] - total_length / 2
        elif horizontal_alignment == 'right':
            s = position[0] - total_length
        else:
            raise ValueError, "horizontal_alignement must be left, center or right"

        bottom_height = position[1] - interval_height
        for i in self._permutation._labels[1]:
            G += line2d([(s,bottom_height), (s+lengths[i],bottom_height)],
                rgbcolor=colors[i])
            if labels == True:
                G += text(str(self._permutation._alphabet.unrank(i)),
                    (s+float(lengths[i])/2,bottom_height-labels_height),
                    horizontal_alignment='center',
                    rgbcolor=colors[i],
                    fontsize=fontsize)
            s += lengths[i]

        return G

    def plot_towers(self, iterations, position=(0,0), colors=None):
        """
        Plot the towers of this interval exchange obtained from Rauzy induction.

        INPUT:

        - ``nb_iterations`` -- the number of steps of Rauzy induction

        - ``colors`` -- (optional) colors for the towers

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('A B', 'B A')
            sage: T = iet.IntervalExchangeTransformation(p, [0.41510826, 0.58489174])
            sage: T.plot_towers(iterations=5)   # not tested (problem with matplotlib font cache)
            Graphics object consisting of 65 graphics primitives
        """
        px,py = map(float, position)

        T,_,towers = self.rauzy_move(iterations=iterations,data=True)
        pi = T.permutation()
        A = pi.alphabet()
        lengths = map(float, T.lengths())

        if colors is None:
            from sage.plot.colors import rainbow
            colors = {a:z for a,z in zip(A, rainbow(len(A)))}

        from sage.plot.graphics import Graphics
        from sage.plot.line import line2d
        from sage.plot.polygon import polygon2d
        from sage.plot.text import text
        G = Graphics()
        x = px
        for letter in pi[0]:
            y = x + lengths[A.rank(letter)]
            tower = towers.image(letter)
            h = tower.length()
            G += line2d([(x,py),(x,py+h)], color='black')
            G += line2d([(y,py),(y,py+h)], color='black')
            for i,a in enumerate(tower):
                G += line2d([(x,py+i),(y,py+i)], color='black')
                G += polygon2d([(x,py+i),(y,py+i),(y,py+i+1),(x,py+i+1)], color=colors[a], alpha=0.4)
                G += text(a, ((x+y)/2, py+i+.5), color='darkgray')
            G += line2d([(x,py+h),(y,py+h)], color='black', linestyle='dashed')
            G += text(letter, ((x+y)/2, py+h+.5), color='black', fontsize='large')
            x = y
        x = px
        G += line2d([(px,py-.5),(px+sum(lengths),py-.5)], color='black')
        for letter in pi[1]:
            y = x + lengths[A.rank(letter)]
            G += line2d([(x,py-.7),(x,py-.3)], color='black')
            G += text(letter, ((x+y)/2, py-.5), color='black', fontsize='large')
            x = y
        G += line2d([(x,py-.7),(x,py-.3)], color='black')

        return G

    plot = plot_two_intervals

    def show(self):
        r"""
        Shows a picture of the interval exchange transformation

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: phi = QQbar((sqrt(5)-1)/2)
            sage: t = iet.IntervalExchangeTransformation(('a b','b a'),[1,phi])
            sage: t.show() # not tested (problem with matplotlib font cache)
        """
        self.plot_two_intervals().show(axes=False)




