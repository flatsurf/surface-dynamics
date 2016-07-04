r"""
Reduced permutations

A reduced (generalized) permutation is better suited to study strata of Abelian
(or quadratic) holomorphic forms on Riemann surfaces. The Rauzy diagram is an
invariant of such a component. Corentin Boissy proved the identification of
Rauzy diagrams with connected components of stratas. But the geometry of the
diagram and the relation with the strata is not yet totally understood.

AUTHORS:

- Vincent Delecroix (2000-09-29): initial version

TESTS::

    sage: from surface_dynamics.interval_exchanges.reduced import ReducedPermutationIET
    sage: ReducedPermutationIET([['a','b'],['b','a']])
    a b
    b a
    sage: ReducedPermutationIET([[1,2,3],[3,1,2]])
    1 2 3
    3 1 2
    sage: from surface_dynamics.interval_exchanges.reduced import ReducedPermutationLI
    sage: ReducedPermutationLI([[1,1],[2,2,3,3,4,4]])
    1 1
    2 2 3 3 4 4
    sage: ReducedPermutationLI([['a','a','b','b','c','c'],['d','d']])
    a a b b c c
    d d
    sage: from surface_dynamics.interval_exchanges.reduced import FlippedReducedPermutationIET
    sage: FlippedReducedPermutationIET([[1,2,3],[3,2,1]],flips=[1,2])
    -1 -2  3
     3 -2 -1
    sage: FlippedReducedPermutationIET([['a','b','c'],['b','c','a']],flips='b')
     a -b  c
    -b  c  a
    sage: from surface_dynamics.interval_exchanges.reduced import FlippedReducedPermutationLI
    sage: FlippedReducedPermutationLI([[1,1],[2,2,3,3,4,4]], flips=[1,4])
    -1 -1
     2  2  3  3 -4 -4
    sage: FlippedReducedPermutationLI([['a','a','b','b'],['c','c']],flips='ac')
    -a -a  b  b
    -c -c
    sage: from surface_dynamics.interval_exchanges.reduced import ReducedRauzyDiagram
    sage: p = ReducedPermutationIET([[1,2,3],[3,2,1]])
    sage: d = ReducedRauzyDiagram(p)
"""
#*****************************************************************************
#       Copyright (C) 2008 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject

from copy import copy

from sage.combinat.words.alphabet import Alphabet
from sage.rings.integer import Integer

from template import OrientablePermutationIET, OrientablePermutationLI   # permutations
from template import FlippedPermutationIET, FlippedPermutationLI         # flipped permutations
from template import RauzyDiagram, FlippedRauzyDiagram

from template import interval_conversion, side_conversion

class ReducedPermutation(SageObject) :
    r"""
    Template for reduced objects.

    .. warning::

        Internal class! Do not use directly!
   """
    def __init__(self, intervals=None, alphabet=None):
        r"""
        INPUT:

        - ``intervals`` - a list of two lists of labels

        - ``alphabet`` - (default: None) alphabet

        TESTS::

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
        """
        self._hash = None

        if intervals is None:
            self._twin = [[],[]]
            self._alphabet = alphabet

        else:
            self._init_twin(intervals)

            if alphabet is None:
                self._init_alphabet(intervals)
            else:
                alphabet = Alphabet(alphabet)
                if alphabet.cardinality() < len(self):
                    raise TypeError("the alphabet is too short")
                self._alphabet = alphabet

    def __getitem__(self, i):
        r"""
        TESTS::

            sage: p = iet.Permutation('a b', 'b a', reduced=True)
            sage: print p[0]
            ['a', 'b']
            sage: print p[1]
            ['b', 'a']
            sage: p.alphabet([0,1])
            sage: print p[0]
            [0, 1]
            sage: print p[1]
            [1, 0]
        """
        return self.list()[i]

def ReducedPermutationsIET_iterator(
    nintervals=None,
    irreducible=True,
    alphabet=None):
    r"""
    Returns an iterator over reduced permutations

    INPUT:

    - ``nintervals`` - integer or None

    - ``irreducible`` - boolean

    - ``alphabet`` - something that should be converted to an alphabet
      of at least nintervals letters

    TESTS::

        sage: from surface_dynamics.all import *

        sage: for p in iet.Permutations_iterator(3,reduced=True,alphabet="abc"):
        ...    print p  #indirect doctest
        a b c
        b c a
        a b c
        c a b
        a b c
        c b a
    """
    from itertools import imap,ifilter
    from sage.combinat.permutation import Permutations

    if irreducible is False:
        if nintervals is None:
            raise NotImplementedError, "choose a number of intervals"
        else:
            assert(isinstance(nintervals,(int,Integer)))
            assert(nintervals > 0)

            a0 = range(1,nintervals+1)
            f = lambda x: ReducedPermutationIET([a0,list(x)],
                alphabet=alphabet)
            return imap(f, Permutations(nintervals))
    else:
        return ifilter(lambda x: x.is_irreducible(),
        ReducedPermutationsIET_iterator(nintervals,False,alphabet))

class ReducedPermutationIET(ReducedPermutation, OrientablePermutationIET):
    """
    Reduced permutation from iet

    Permutation from iet without numerotation of intervals. For initialization,
    you should use GeneralizedPermutation which is the class factory for all
    permutation types.

    EXAMPLES::

        sage: from surface_dynamics.all import *

    Equality testing (no equality of letters but just of ordering)::

        sage: p = iet.Permutation('a b c', 'c b a', reduced = True)
        sage: q = iet.Permutation('p q r', 'r q p', reduced = True)
        sage: p == q
        True

    Reducibility testing::

        sage: p = iet.Permutation('a b c', 'c b a', reduced = True)
        sage: p.is_irreducible()
        True

    ::

        sage: q = iet.Permutation('a b c d', 'b a d c', reduced = True)
        sage: q.is_irreducible()
        False


    Rauzy movability and Rauzy move::

        sage: p = iet.Permutation('a b c', 'c b a', reduced = True)
        sage: p.has_rauzy_move(1)
        True
        sage: print p.rauzy_move(1)
        a b c
        b c a

    Rauzy diagrams::

        sage: p = iet.Permutation('a b c d', 'd a b c')
        sage: p_red = iet.Permutation('a b c d', 'd a b c', reduced = True)
        sage: d = p.rauzy_diagram()
        sage: d_red = p_red.rauzy_diagram()
        sage: p.rauzy_move(0) in d
        True
        sage: print d.cardinality(), d_red.cardinality()
        12 6
    """
    def list(self):
        r"""
        Returns a list of two list that represents the permutation.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('a b','b a',reduced=True)
            sage: p.list() == [['a', 'b'], ['b', 'a']]
            True
            sage: p = iet.GeneralizedPermutation('a a','b b',reduced=True)
            sage: p.list() == [['a', 'a'], ['b', 'b']]
            True
        """
        return [
            map(lambda x: self._alphabet.unrank(x), range(len(self._twin[0]))),
            map(lambda x: self._alphabet.unrank(x), self._twin[1])]

    def rauzy_move_relabel(self, winner, side='right'):
        r"""
        Returns the relabelization obtained from this move.

        EXAMPLE::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c d','d c b a')
            sage: q = p.reduced()
            sage: p_t = p.rauzy_move('t')
            sage: q_t = q.rauzy_move('t')
            sage: s_t = q.rauzy_move_relabel('t')
            sage: print s_t
            a->a, b->b, c->c, d->d
            sage: map(s_t, p_t[0]) == map(Word, q_t[0])
            True
            sage: map(s_t, p_t[1]) == map(Word, q_t[1])
            True
            sage: p_b = p.rauzy_move('b')
            sage: q_b = q.rauzy_move('b')
            sage: s_b = q.rauzy_move_relabel('b')
            sage: print s_b
            a->a, b->d, c->b, d->c
            sage: map(s_b, q_b[0]) == map(Word, p_b[0])
            True
            sage: map(s_b, q_b[1]) == map(Word, p_b[1])
            True
        """
        from surface_dynamics.interval_exchanges.labelled import LabelledPermutationIET
        from sage.combinat.words.morphism import WordMorphism

        winner = interval_conversion(winner)
        side = side_conversion(side)

        p = LabelledPermutationIET(self.list())

        l0_q = p.rauzy_move(winner, side).list()[0]

        d = dict([(self._alphabet[i],l0_q[i]) for i in range(len(self))])

        return WordMorphism(d)

    def rauzy_diagram(self, extended=False, **kwds):
        r"""
        Returns a Rauzy diagram associated to this permutation

        INPUT:

        - ``extended`` - boolean (default: False) - if True return extended
          Rauzy diagram

        OUTPUT:

        A Rauzy diagram

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c d', 'd a b c',reduced=True)
            sage: d = p.rauzy_diagram()
            sage: p.rauzy_move(0) in d
            True
            sage: p.rauzy_move(1) in d
            True

        For more information, try help RauzyDiagram
        """
        options = kwds
        if extended:
            options['symmetric'] = True
        return ReducedRauzyDiagram(self,**options)

    def rauzy_class_cardinality(self, extended=False):
        r"""
        Cardinality of Rauzy diagram

        As proved in [Del10], there exists a closed formula for the cardinality
        of Rauzy diagrams. This function uses the formula without any
        explicit computation of the rauzy class.

        INPUT:

        - ``extended`` - boolean (default: False) - return cardinality for
          extended Rauzy diagrams

        EXAMPLES::

            sage: from surface_dynamics.all import *

        Examples for permutations such that the suspensions are tori::

            sage: p=iet.Permutation('a b','b a',reduced=True)
            sage: p.stratum()
            H_1(0)
            sage: p.rauzy_diagram()
            Rauzy diagram with 1 permutation
            sage: p.rauzy_class_cardinality()
            1

            sage: p = iet.Permutation('a 1 b','b 1 a',reduced=True)
            sage: p.stratum()
            H_1(0^2)
            sage: p.rauzy_diagram()
            Rauzy diagram with 3 permutations
            sage: p.rauzy_class_cardinality()
            3

            sage: p = iet.Permutation('a 1 2 b','b 1 2 a',reduced=True)
            sage: p.stratum()
            H_1(0^3)
            sage: p.rauzy_diagram()
            Rauzy diagram with 6 permutations
            sage: p.rauzy_class_cardinality()
            6

            sage: p = iet.Permutation('a 1 2 3 b','b 1 2 3 a',reduced=True)
            sage: p.rauzy_class_cardinality()
            10
            sage: p = iet.Permutation('a 1 2 3 4 b','b 1 2 3 4 a',reduced=True)
            sage: p.rauzy_class_cardinality()
            15

        You should have recognize the sequence 1, 3, 6, 10, 15... which is the
        sequence with general term the binomial ``n(n+1)/2``.

        An example of extended Rauzy diagram which is different from Rauzy
        diagram::

            sage: p = iet.Permutation('a b c d e f g','g c b f e d a',reduced=True)
            sage: p.marked_profile()
            2o4 [4, 2]
            sage: pp = p.left_right_inverse()
            sage: pp.marked_profile()
            4o2 [4, 2]
            sage: p.rauzy_class_cardinality()
            261
            sage: pp.rauzy_class_cardinality()
            509
            sage: p.rauzy_class_cardinality(extended=True)
            770
            sage: 261 + 509 == 770
            True

        And one can check that this the algorithm for cardinality is True::

            sage: p.rauzy_diagram()
            Rauzy diagram with 261 permutations
            sage: pp.rauzy_diagram()
            Rauzy diagram with 509 permutations
            sage: p.rauzy_diagram(extended=True)
            Rauzy diagram with 770 permutations

        And we end by an example of Rauzy diagram associated to an hyperelliptic
        component::

            sage: p = iet.Permutation('a b c d e 0 f','f e d c b 0 a',reduced=True)
            sage: p.rauzy_diagram()
            Rauzy diagram with 37 permutations
            sage: p.rauzy_class_cardinality()
            37
            sage: p.rauzy_diagram(extended=True)
            Rauzy diagram with 254 permutations
            sage: p.rauzy_class_cardinality(extended=True)
            254
        """
        from rauzy_class_cardinality import gamma_irr,delta_irr
        from sage.arith.all import binomial

        s = self.stratum()
        mp = self.marked_profile()
        p = mp.partition()

        if extended:
            return self.stratum_component().rauzy_class_cardinality()

        if s.is_connected():  # easy stuff
            return gamma_irr(p,mp.left())

        cc = self.stratum_component()
        zeros = s.zeros(fake_zeros=False)
        g = s.genus()
        if zeros == [2*g-2] or zeros == [g-1,g-1]:  # hyp component + others
            if cc.spin() == s.hyperelliptic_component().spin():
                d = len(self) # number of intervals of self
                if g == 1: # genus 1 is particular
                    nb_hyp = binomial(d,2)
                else:
                    k = s.nb_fake_zeros()
                    d -= k
                    if extended:
                        nb_hyp = binomial(d+k,k) * (2*d-1) + binomial(d+k-1,k-1) * d
                    elif mp.left() == 1:
                        nb_hyp = binomial(d+k,k-1) * (2**(d-1)-1+d)
                    elif mp.left() == d-1:
                        nb_hyp = binomial(d+k,k) * (2**(d-1)-1)
                    else:
                        raise ValueError, "this should not happen"

                if cc == s.hyperelliptic_component(): # hyp
                    return nb_hyp

                if s.number_of_components() == 2: # other is alone
                    N = gamma_irr(p,mp.left())

                else: # others are odd and even
                    t = (-1)**(1-cc.spin())
                    N = (gamma_irr(p,mp.left()) + t * delta_irr(p,mp.left()))/2

                return N - nb_hyp

            else: # cc.spin() != s.hyperelliptic_component().spin()
                t = (-1)**(1-cc.spin())
                return (gamma_irr(p,mp.left()) + t * delta_irr(p,mp.left()))/2

        else:   # odd and even components
            t = (-1)**(1-cc.spin())  # 1 if odd and -1 if even
            return (gamma_irr(p,mp.left()) + t * delta_irr(p,mp.left()))/2

class ReducedPermutationLI(ReducedPermutation, OrientablePermutationLI):
    r"""
    Reduced quadratic (or generalized) permutation.

    EXAMPLES::

        sage: from surface_dynamics.all import *

    Reducibility testing::

        sage: p = iet.GeneralizedPermutation('a b b', 'c c a', reduced = True)
        sage: p.is_irreducible()
        True

    ::

        sage: p = iet.GeneralizedPermutation('a b c a', 'b d d c', reduced = True)
        sage: p.is_irreducible()
        False
        sage: test, decomposition = p.is_irreducible(return_decomposition = True)
        sage: test
        False
        sage: decomposition
        (['a'], ['c', 'a'], [], ['c'])

    Rauzy movavability and Rauzy move::

        sage: p = iet.GeneralizedPermutation('a b b', 'c c a', reduced = True)
        sage: p.has_rauzy_move(0)
        True
        sage: p.rauzy_move(0)
        a a b b
        c c
        sage: p.rauzy_move(0).has_rauzy_move(0)
        False
        sage: p.rauzy_move(1)
        a b b
        c c a

    Rauzy diagrams::

        sage: p_red = iet.GeneralizedPermutation('a b b', 'c c a', reduced = True)
        sage: d_red = p_red.rauzy_diagram()
        sage: d_red.cardinality()
        4
    """
    def list(self) :
        r"""
        The permutations as a list of two lists.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('a b b', 'c c a', reduced = True)
            sage: list(p)
            [['a', 'b', 'b'], ['c', 'c', 'a']]
        """
        i_a = 0
        l = [[False]*len(self._twin[0]),[False]*len(self._twin[1])]
        # False means empty here
        for i in range(2) :
            for j in range(len(l[i])) :
                if  l[i][j] is False :
                    l[i][j] = self._alphabet[i_a]
                    l[self._twin[i][j][0]][self._twin[i][j][1]] = self._alphabet[i_a]
                    i_a += 1
        return l

    def rauzy_diagram(self, **kargs):
        r"""
        Returns the associated Rauzy diagram.

        The Rauzy diagram of a permutation corresponds to all permutations
        that we could obtain from this one by Rauzy move. The set obtained
        is a labelled Graph. The label of vertices being 0 or 1 depending
        on the type.

        OUTPUT:

        Rauzy diagram -- the graph of permutations obtained by rauzy induction

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c d', 'd a b c')
            sage: d = p.rauzy_diagram()
        """
        return ReducedRauzyDiagram(self, **kargs)

def labelize_flip(couple):
    r"""
    Returns a string from a 2-uple couple of the form (name, flip).

    TESTS::

        sage: from surface_dynamics.interval_exchanges.reduced import labelize_flip
        sage: labelize_flip((4,1))
        ' 4'
        sage: labelize_flip(('a',-1))
        '-a'
    """
    if couple[1] == -1: return '-' + str(couple[0])
    return ' ' + str(couple[0])

class FlippedReducedPermutation(ReducedPermutation):
    r"""
    Flipped Reduced Permutation.

    .. warning::

        Internal class! Do not use directly!

    INPUT:

    - ``intervals`` - a list of two lists

    - ``flips`` - the flipped letters

    - ``alphabet`` - an alphabet
    """
    def __init__(self, intervals=None, flips=None, alphabet=None):
        r"""
        TESTS::

            sage: from surface_dynamics.all import *

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
        self._hash = None

        if intervals is None:
            self._twin = [[],[]]
            self._flips = [[],[]]
            self._alphabet = None

        else:
            if flips is None: flips = []

            if alphabet is None : self._init_alphabet(intervals)
            else : self._alphabet = Alphabet(alphabet)

            self._init_twin(intervals)
            self._init_flips(intervals, flips)

            self._hash = None

class FlippedReducedPermutationIET(
    FlippedReducedPermutation,
    FlippedPermutationIET,
    ReducedPermutationIET):
    r"""
    Flipped Reduced Permutation from iet

    EXAMPLES::

        sage: from surface_dynamics.all import *

        sage: p = iet.Permutation('a b c', 'c b a', flips=['a'], reduced=True)
        sage: p.rauzy_move(1)
        -a -b  c
        -a  c -b

    TESTS::

        sage: p = iet.Permutation('a b','b a',flips=['a'])
        sage: p == loads(dumps(p))
        True
    """
    def list(self, flips=False):
        r"""
        Returns a list representation of self.

        INPUT:

        - ``flips`` - boolean (default: False) if True the output contains
           2-uple of (label, flip)

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b','b a',reduced=True,flips='b')
            sage: p.list(flips=True)
            [[('a', 1), ('b', -1)], [('b', -1), ('a', 1)]]
            sage: p.list(flips=False)
            [['a', 'b'], ['b', 'a']]
            sage: p.alphabet([0,1])
            sage: p.list(flips=True)
            [[(0, 1), (1, -1)], [(1, -1), (0, 1)]]
            sage: p.list(flips=False)
            [[0, 1], [1, 0]]

        One can recover the initial permutation from this list::

            sage: p = iet.Permutation('a b','b a',reduced=True,flips='a')
            sage: iet.Permutation(p.list(), flips=p.flips(), reduced=True) == p
            True
        """
        if flips:
            a0 = zip(map(self.alphabet().unrank, range(0,len(self))), self._flips[0])
            a1 = zip(map(self.alphabet().unrank, self._twin[1]), self._flips[1])

        else:
            a0 = map(self.alphabet().unrank, range(0,len(self)))
            a1 = map(self.alphabet().unrank, self._twin[1])

        return [a0,a1]

    def rauzy_diagram(self, **kargs):
        r"""
        Returns the associated Rauzy diagram.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b','b a',reduced=True,flips='a')
            sage: r = p.rauzy_diagram()
            sage: p in r
            True
        """
        return FlippedReducedRauzyDiagram(self, **kargs)

class FlippedReducedPermutationLI(
    FlippedReducedPermutation,
    FlippedPermutationLI,
    ReducedPermutationLI):
    r"""
    Flipped Reduced Permutation from li

    EXAMPLES:

            sage: from surface_dynamics.all import *

    Creation using the GeneralizedPermutation function::

        sage: p = iet.GeneralizedPermutation('a a b', 'b c c', reduced=True, flips='a')
    """
    def list(self, flips=False):
        r"""
        Returns a list representation of self.

        INPUT:

        - ``flips`` - boolean (default: False) return the list with flips

        EXAMPLES:

            sage: from surface_dynamics.all import *

        ::

            sage: p = iet.GeneralizedPermutation('a a','b b',reduced=True,flips='a')
            sage: p.list(flips=True)
            [[('a', -1), ('a', -1)], [('b', 1), ('b', 1)]]
            sage: p.list(flips=False)
            [['a', 'a'], ['b', 'b']]

            sage: p = iet.GeneralizedPermutation('a a b','b c c',reduced=True,flips='abc')
            sage: p.list(flips=True)
            [[('a', -1), ('a', -1), ('b', -1)], [('b', -1), ('c', -1), ('c', -1)]]
            sage: p.list(flips=False)
            [['a', 'a', 'b'], ['b', 'c', 'c']]

        one can rebuild the permutation from the list::

            sage: p = iet.GeneralizedPermutation('a a b','b c c',flips='a',reduced=True)
            sage: iet.GeneralizedPermutation(p.list(),flips=p.flips(),reduced=True) == p
            True
        """
        i_a = 0
        l = [[False]*len(self._twin[0]),[False]*len(self._twin[1])]

        if flips:
            for i in range(2):  # False means empty here
                for j in range(len(l[i])):
                   if  l[i][j] is False:
                        l[i][j] = (self._alphabet.unrank(i_a), self._flips[i][j])
                        l[self._twin[i][j][0]][self._twin[i][j][1]] = l[i][j]
                        i_a += 1

        else:
            for i in range(2):  # False means empty here
                for j in range(len(l[i])):
                   if  l[i][j] is False:
                        l[i][j] = self._alphabet.unrank(i_a)
                        l[self._twin[i][j][0]][self._twin[i][j][1]] = l[i][j]
                        i_a += 1
        return l

    def rauzy_diagram(self, **kargs):
        r"""
        Returns the associated Rauzy diagram.

        For more explanation and a list of arguments try help(iet.RauzyDiagram)

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('a a b','c c b',reduced=True)
            sage: r = p.rauzy_diagram()
            sage: p in r
            True
        """
        return FlippedReducedRauzyDiagram(self, **kargs)

class ReducedRauzyDiagram(RauzyDiagram):
    r"""
    Rauzy diagram of reduced permutations
    """
    def _permutation_to_vertex(self, p):
        r"""
        The necessary data to store the permutation.

        TESTS::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c','c b a',reduced=True)   #indirect doctest
            sage: r = p.rauzy_diagram()
            sage: p in r
            True
        """
        return (tuple(p._twin[0]), tuple(p._twin[1]))

    def _set_element(self, data=None):
        r"""
        Sets self._element with data.

        TESTS::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c','c b a',reduced=True)
            sage: r = p.rauzy_diagram()
            sage: p in r   #indirect doctest
            True
        """
        self._element._twin = [list(data[0]), list(data[1])]

class FlippedReducedRauzyDiagram(FlippedRauzyDiagram, ReducedRauzyDiagram):
    r"""
    Rauzy diagram of flipped reduced permutations.
    """
    def _permutation_to_vertex(self, p):
        r"""
        TESTS::

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('a b b','c c a',flips='a',reduced=True)
            sage: r = p.rauzy_diagram()
            sage: p in r   #indirect doctest
            True
        """
        return ((tuple(p._twin[0]), tuple(p._twin[1])),
                (tuple(p._flips[0]), tuple(p._flips[1])))

    def _set_element(self, data=None):
        r"""
        Sets self._element with data.

        TESTS::

            sage: from surface_dynamics.all import *

            sage: r = iet.RauzyDiagram('a b c','c b a',flips='b',reduced=True)   #indirect doctest
        """
        self._element._twin = [list(data[0][0]), list(data[0][1])]
        self._element._flips = [list(data[1][0]), list(data[1][1])]

