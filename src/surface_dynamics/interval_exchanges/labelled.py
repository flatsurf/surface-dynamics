r"""
Labelled permutations

A labelled (generalized) permutation is better suited to study the dynamic of
a translation surface than a reduced one (see the module
:mod:`surface_dynamics.interval_exchanges.reduced`). The latter is more adapted to the
study of strata. This kind of permutation was introduced by Yoccoz [Yoc05]_
(see also [MMY03]_).

In fact, there is a geometric counterpart of labelled permutations. They
correspond to translation surface with marked outgoing separatrices (i.e. we fi
a label for each of them).

Remarks that Rauzy diagram of reduced objects are significantly smaller than
the one for labelled object (for the permutation a b d b e / e d c a c the
labelled Rauzy diagram contains 8760 permutations, and the reduced only 73).
But, as it is in geometrical way, the labelled Rauzy diagram is a covering of
the reduced Rauzy diagram.

AUTHORS:

- Vincent Delecroix (2009-09-29) : initial version

- Vincent Delecroix (2010-02-11) : correction and simplification of datatypes

TESTS::

    sage: from surface_dynamics.interval_exchanges.labelled import LabelledPermutationIET
    sage: LabelledPermutationIET([['a','b','c'],['c','b','a']])
    a b c
    c b a
    sage: LabelledPermutationIET([[1,2,3,4],[4,1,2,3]])
    1 2 3 4
    4 1 2 3
    sage: from surface_dynamics.interval_exchanges.labelled import LabelledPermutationLI
    sage: LabelledPermutationLI([[1,1],[2,2,3,3,4,4]])
    1 1
    2 2 3 3 4 4
    sage: LabelledPermutationLI([['a','a','b','b','c','c'],['d','d']])
    a a b b c c
    d d
    sage: from surface_dynamics.interval_exchanges.labelled import FlippedLabelledPermutationIET
    sage: FlippedLabelledPermutationIET([[1,2,3],[3,2,1]],flips=[1,2])
    -1 -2  3
     3 -2 -1
    sage: FlippedLabelledPermutationIET([['a','b','c'],['b','c','a']],flips='b')
     a -b  c
    -b  c  a
    sage: from surface_dynamics.interval_exchanges.labelled import FlippedLabelledPermutationLI
    sage: FlippedLabelledPermutationLI([[1,1],[2,2,3,3,4,4]], flips=[1,4])
    -1 -1
     2  2  3  3 -4 -4
    sage: FlippedLabelledPermutationLI([['a','a','b','b'],['c','c']],flips='ac')
    -a -a  b  b
    -c -c
    sage: from surface_dynamics.interval_exchanges.labelled import LabelledRauzyDiagram
    sage: p = LabelledPermutationIET([[1,2,3],[3,2,1]])
    sage: d1 = LabelledRauzyDiagram(p)
    sage: p = LabelledPermutationIET([['a','b'],['b','a']])
    sage: d = p.rauzy_diagram()
    sage: g1 = d.path(p, 'top', 'bottom')
    sage: g1.matrix()
    [1 1]
    [1 2]
    sage: g2 = d.path(p, 'bottom', 'top')
    sage: g2.matrix()
    [2 1]
    [1 1]
    sage: p = LabelledPermutationIET([['a','b','c','d'],['d','c','b','a']])
    sage: d = p.rauzy_diagram()
    sage: g = d.path(p, 't', 't', 'b', 't', 'b', 'b', 't', 'b')
    sage: g
    Path of length 8 in a Rauzy diagram
    sage: g.is_loop()
    True
    sage: g.is_full()
    True
    sage: s1 = g.orbit_substitution()
    sage: print s1
    a->adbd, b->adbdbd, c->adccd, d->adcd
    sage: s2 = g.interval_substitution()
    sage: print s2
    a->abcd, b->bab, c->cdc, d->dcbababcd
    sage: s1.incidence_matrix() == s2.incidence_matrix().transpose()
    True

REFERENCES:

.. [Yoc05] Jean-Cristophe Yoccoz "Echange d'Intervalles", Cours au college de
   France

.. [MMY03] Jean-Cristophe Yoccoz, Stefano Marmi and Pierre Moussa "On the
   cohomological equation for interval exchange maps", arXiv:math/0304469v1
"""
#*****************************************************************************
#       Copyright (C) 2008 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.misc.lazy_attribute import lazy_attribute

from copy import copy

import time
import surface_dynamics.interval_exchanges.lyapunov_exponents as lyapunov_exponents  # the cython bindings

from sage.combinat.words.alphabet import Alphabet, OrderedAlphabet
from sage.combinat.words.morphism import WordMorphism

from sage.matrix.constructor import Matrix, identity_matrix
from sage.rings.integer import Integer
from sage.combinat.words.alphabet import Alphabet
from sage.rings.infinity import Infinity

from template import OrientablePermutationIET, OrientablePermutationLI
from template import FlippedPermutationIET, FlippedPermutationLI
from template import RauzyDiagram, FlippedRauzyDiagram
from template import interval_conversion, side_conversion

def mean_and_std_dev(l):
    r"""
    Return the mean and standard deviation of the floatting point numbers in
    the list l.

    The implementation is very naive and should not be used for large list
    (>1000) of numbers.

    .. NOTE::
    mean and std are implemented in Sage but are quite buggy!
    """
    from math import sqrt
    m = sum(l) / len(l)
    if len(l) == 1:
        d = 0
    else:
        d = sum((x-m)**2 for x in l) / (len(l)-1)
    return m,sqrt(d)


class LabelledPermutation(SageObject):
    r"""
    General template for labelled objects.

    .. warning::
       Internal class! Do not use directly!
    """
    def __init__(self, intervals=None, alphabet=None):
        r"""
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
        """
        self._hash = None

        if intervals is None:
            self._twin = [[], []]
            self._labels = [[],[]]
            self._alphabet = None

        else:
            self._init_twin(intervals)

            if alphabet is not None:
                self._set_alphabet(alphabet)
            else:
                self._init_alphabet(intervals)

            self._labels = [
                map(self._alphabet.rank, intervals[0]),
                map(self._alphabet.rank, intervals[1])]

    def __eq__(self, other):
        r"""
        Equality test.

        TESTS::

            sage: from surface_dynamics.all import *

            sage: p1 = iet.Permutation('a b', 'b a', alphabet='ab')
            sage: p2 = iet.Permutation('a b', 'b a', alphabet='ba')
            sage: q1 = iet.Permutation('b a', 'a b', alphabet='ab')
            sage: q2 = iet.Permutation('b a', 'a b', alphabet='ba')
            sage: p1 == p2 or p2 == p1
            False
            sage: p1 == q1 or q1 == p1
            False
            sage: p1 == q2 or q2 == p1
            False
            sage: p2 == q1 or q1 == p2
            False
            sage: p2 == q2 or q2 == p2
            False
            sage: q1 == q2 or q2 == q1
            False

        ::

            sage: p1 = iet.GeneralizedPermutation('a a','b b',alphabet='ab')
            sage: p2 = iet.GeneralizedPermutation('a a','b b',alphabet='ba')
            sage: q1 = iet.GeneralizedPermutation('b b','a a',alphabet='ab')
            sage: q2 = iet.GeneralizedPermutation('b b','a a',alphabet='ba')
            sage: p1 == p2 or p2 == p1
            False
            sage: p1 == q1 or q1 == p1
            False
            sage: p1 == q2 or q2 == p1
            False
            sage: p2 == q1 or q1 == p2
            False
            sage: p2 == q2 or q2 == p2
            False
            sage: q1 == q2 or q2 == q1
            False
        """
        return (
            type(self) == type(other) and
            self._alphabet == other._alphabet and
            self._labels == other._labels)

    def __ne__(self, other):
        r"""
        Non equality test.

        TESTS::

            sage: from surface_dynamics.all import * 
            sage: p1 = iet.Permutation('a b', 'b a', alphabet='ab')
            sage: p2 = iet.Permutation('a b', 'b a', alphabet='ba')
            sage: q1 = iet.Permutation('b a', 'a b', alphabet='ab')
            sage: q2 = iet.Permutation('b a', 'a b', alphabet='ba')
            sage: p1 != p2 and p2 != p1
            True
            sage: p1 != q1 and q1 != p1
            True
            sage: p1 != q2 and q2 != p1
            True
            sage: p2 != q1 and q1 != p2
            True
            sage: p2 != q2 and q2 != p2
            True
            sage: q1 != q2 and q2 != q1
            True

        ::

            sage: p1 = iet.GeneralizedPermutation('a a','b b',alphabet='ab')
            sage: p2 = iet.GeneralizedPermutation('a a','b b',alphabet='ba')
            sage: q1 = iet.GeneralizedPermutation('b b','a a',alphabet='ab')
            sage: q2 = iet.GeneralizedPermutation('b b','a a',alphabet='ba')
            sage: p1 != p2 or p2 != p1
            True
            sage: p1 != q1 or q1 != p1
            True
            sage: p1 != q2 or q2 != p1
            True
            sage: p2 != q1 or q1 != p2
            True
            sage: p2 != q2 or q2 != p2
            True
            sage: q1 != q2 or q2 != q1
            True
        """
        return (
            type(self) != type(other) or
            self._alphabet != other._alphabet or
            self._labels != other._labels)

    def __getitem__(self,i):
        r"""
        TESTS::

            sage: p = iet.Permutation([0,1,2,3],[3,2,1,0])
            sage: p[0][0]
            0
            sage: p[1][2]
            1
            sage: p = iet.Permutation('a b c','c b a')
            sage: p[0][1]
            'b'
            sage: p[1][2]
            'a'
        """
        return map(self._alphabet.unrank, self._labels[i])

        def change(x,index):
            return "Interval(%s,%i)"%(str(x), 1 if index else -1)

        for i in range(2):
            for j in range(len(labels[i])):
                labels[i][j] = change(labels[i][j], twin_to_index[i][j])

        s = "IntExchange([[" + ', '.join(labels[0]) + "], [" + ', '.join(labels[1]) +"]])"        
        print s

    def list(self):
        r"""
        Returns a list of two lists corresponding to the intervals.

        OUTPUT:

        list -- two lists of labels

        EXAMPLES::

            sage: from surface_dynamics.all import *

        The list of an permutation from iet::

            sage: p1 = iet.Permutation('1 2 3', '3 1 2')
            sage: p1.list()
            [['1', '2', '3'], ['3', '1', '2']]
            sage: p1.alphabet("abc")
            sage: p1.list()
            [['a', 'b', 'c'], ['c', 'a', 'b']]

        Recovering the permutation from this list (and the alphabet)::

            sage: q1 = iet.Permutation(p1.list(),alphabet=p1.alphabet())
            sage: p1 == q1
            True

        The list of a quadratic permutation::

            sage: p2 = iet.GeneralizedPermutation('g o o', 'd d g')
            sage: p2.list()
            [['g', 'o', 'o'], ['d', 'd', 'g']]

        Recovering the permutation::

            sage: q2 = iet.GeneralizedPermutation(p2.list(),alphabet=p2.alphabet())
            sage: p2 == q2
            True
        """
        a0 = map(self._alphabet.unrank, self._labels[0])
        a1 = map(self._alphabet.unrank, self._labels[1])
        return [a0, a1]

    def rauzy_move_matrix(self, winner=None, side='right'):
        r"""
        Returns the Rauzy move matrix.

        This matrix corresponds to the action of a Rauzy move on the vector of
        lengths. By convention (to get a positive matrix), the matrix is define
        as the inverse transformation on the length vector.

        OUTPUT:

        matrix -- a square matrix of positive integers

        EXAMPLES:

            sage: from surface_dynamics.all import *

        ::

            sage: p = iet.Permutation('a b','b a')
            sage: p.rauzy_move_matrix('t')
            [1 0]
            [1 1]
            sage: p.rauzy_move_matrix('b')
            [1 1]
            [0 1]

        ::

            sage: p = iet.Permutation('a b c d','b d a c')
            sage: q = p.left_right_inverse()
            sage: m0 = p.rauzy_move_matrix(winner='top',side='right')
            sage: n0 = q.rauzy_move_matrix(winner='top',side='left')
            sage: m0 == n0
            True
            sage: m1 = p.rauzy_move_matrix(winner='bottom',side='right')
            sage: n1 = q.rauzy_move_matrix(winner='bottom',side='left')
            sage: m1 == n1
            True
        """
        if winner is None and side is None:
            return identity_matrix(len(self))

        winner = interval_conversion(winner)
        side = side_conversion(side)

        winner_letter = self._labels[winner][side]
        loser_letter = self._labels[1-winner][side]

        m = copy(identity_matrix(len(self)))
        m[winner_letter, loser_letter] = 1

        return m

    def rauzy_move_winner(self,winner=None,side=None):
        r"""
        Returns the winner of a Rauzy move.

        INPUT:

        - ``winner`` - either 'top' or 'bottom' ('t' or 'b' for short)

        - ``side`` - either 'left' or 'right' ('l' or 'r' for short)

        OUTPUT:

        -- a label

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c d','b d a c')
            sage: p.rauzy_move_winner('top','right')
            'd'
            sage: p.rauzy_move_winner('bottom','right')
            'c'
            sage: p.rauzy_move_winner('top','left')
            'a'
            sage: p.rauzy_move_winner('bottom','left')
            'b'

        ::

            sage: p = iet.GeneralizedPermutation('a b b c','d c a e d e')
            sage: p.rauzy_move_winner('top','right')
            'c'
            sage: p.rauzy_move_winner('bottom','right')
            'e'
            sage: p.rauzy_move_winner('top','left')
            'a'
            sage: p.rauzy_move_winner('bottom','left')
            'd'
        """
        if winner is None and side is None:
            return None

        winner = interval_conversion(winner)
        side = side_conversion(side)

        return self[winner][side]

    def rauzy_move_loser(self,winner=None,side=None):
        r"""
        Returns the loser of a Rauzy move

        INPUT:

        - ``winner`` - either 'top' or 'bottom' ('t' or 'b' for short)

        - ``side`` - either 'left' or 'right' ('l' or 'r' for short)

        OUTPUT:

        -- a label

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c d','b d a c')
            sage: p.rauzy_move_loser('top','right')
            'c'
            sage: p.rauzy_move_loser('bottom','right')
            'd'
            sage: p.rauzy_move_loser('top','left')
            'b'
            sage: p.rauzy_move_loser('bottom','left')
            'a'
        """
        if winner is None and side is None:
            return None

        winner = interval_conversion(winner)
        side = side_conversion(side)

        return self[1-winner][side]

    def lyapunov_exponents_H_plus(self, nb_vectors=None, nb_experiments=100,
                                  nb_iterations=32768, return_speed=False, 
                                  verbose=False, output_file=None, lengths=None, 
                                  sigma=None, isotypic_decomposition=False):
        r"""
        Compute the H^+ Lyapunov exponents in  the covering locus.

        It calls the C-library lyap_exp interfaced with Cython. The computation
        might be significantly faster if ``nb_vectors=1`` (or if it is not
        provided but genus is 1).

        INPUT:
 
        - ``nb_vectors`` -- the number of exponents to compute. The number of
          vectors must not exceed the dimension of the space!

         - ``nb_experiments`` -- the number of experiments to perform. It might
           be around 100 (default value) in order that the estimation of
           confidence interval is accurate enough.

         - ``nb_iterations`` -- the number of iteration of the Rauzy-Zorich
           algorithm to perform for each experiments. The default is 2^15=32768
           which is rather small but provide a good compromise between speed and
           quality of approximation.

        - ``verbose`` -- if ``True`` provide additional informations rather than
          returning only the Lyapunov exponents (i.e. ellapsed time, confidence
          intervals, ...)

        - ``output_file`` -- if provided (as a file object or a string) output
          the additional information in the given file rather than on the
          standard output.

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: QuadraticStratum([1,1,-1,-1]).one_component().lyapunov_exponents_H_plus() #abs tol .1
            [2./3]

            sage: from surface_dynamics.all import *
            sage: q = QuadraticStratum([1,1,-1,-1]).one_component()
            sage: p = q.permutation_representative().orientation_cover()
            sage: p.lyapunov_exponents_H_plus() #abs tol .1
            [1., 2./3, 1./3]
        """
        n = len(self)

        if nb_vectors is None:
            nb_vectors = self.stratum().genus()

        if output_file is None:
            from sys import stdout
            output_file = stdout
        elif isinstance(output_file, str):
            output_file = open(output_file, "w")

        if sigma == None or sigma == []: sigma = [0]*len(self)        #look at the trivial cover
        if isinstance(sigma, dict):
            sigma = map(lambda x : sigma[self._alphabet.unrank(x)], range(n))
        if isinstance(sigma, list) and isinstance(sigma[0], list):
            sigma = reduce(lambda x,y: x+y, sigma)

        nb_vectors = int(nb_vectors)
        nb_experiments = int(nb_experiments)
        nb_iterations = int(nb_iterations)

        if verbose:
            output_file.write("Stratum : " + str(self.stratum()))
            output_file.write("\n")

        if nb_vectors < 0 :     raise ValueError("the number of vectors must be positive")
        if nb_vectors == 0:     return []
        if nb_experiments <= 0: raise ValueError("the number of experiments must be positive")
        if nb_iterations <= 0 : raise ValueError("the number of iterations must be positive")
        if len(sigma) %n != 0 : raise ValueError("you must give a permutation for each interval")
        
        #Translate our structure to the C structure"
        k = len(self[0])
        def convert((i,j)):
            return(j + i*k)

        gp, twin = range(2*n), range(2*n)

        for i in range(2):
            for j in range(len(self[i])):
                gp[convert((i,j))] = int(self._alphabet.rank(self[i][j]))
                if isinstance(self._twin[i][j], int):
                       twin[convert((i,j))] = int(convert((1-i, self._twin[i][j])))
                else:
                       twin[convert((i,j))] = int(convert(self._twin[i][j]))


        if lengths is not None:
            lengths = map(int, lengths)
        sigma = map(int, sigma)

        projections = None
        if sigma and isotypic_decomposition:
            projections = [self.isotypic_projection_matrix(i) for i in range(self.n_characters())]

        t0 = time.time()
        res = lyapunov_exponents.lyapunov_exponents_H_plus_cover(
            gp, int(k), twin, sigma, 
            nb_experiments, nb_iterations,
            [nb_vectors], projections, None, verbose)
        t1 = time.time()

        res_final = []

        m,d = mean_and_std_dev(res[0])
        lexp = m

        if verbose:
            from math import log, floor, sqrt
            output_file.write("sample of %d experiments\n"%nb_experiments)
            output_file.write("%d iterations (~2^%d)\n"%(
                    nb_iterations,
                    floor(log(nb_iterations) / log(2))))
            output_file.write("ellapsed time %s\n"%time.strftime("%H:%M:%S",time.gmtime(t1-t0)))
            output_file.write("Lexp Rauzy-Zorich: %f (std. dev. = %f, conf. rad. 0.01 = %f)\n"%(
                    m,d, 2.576*d/sqrt(nb_experiments)))
        for i in xrange(1,nb_vectors+1):
            m,d = mean_and_std_dev(res[i])
            if verbose:
                output_file.write("theta%d           : %f (std. dev. = %f, conf. rad. 0.01 = %f)\n"%(
                    i,m,d, 2.576*d/sqrt(nb_experiments)))
            res_final.append(m)

        if return_speed: return (lexp, res_final)
        else: return res_final


def LabelledPermutationsIET_iterator(
    nintervals=None,
    irreducible=True,
    alphabet=None):
    r"""
    Returns an iterator over labelled permutations.

    INPUT:

    - ``nintervals`` - integer or None

    - ``irreducible`` - boolean (default: True)

    - ``alphabet`` - something that should be converted to an alphabet of at least nintervals letters

    OUTPUT:

    iterator -- an iterator over permutations

    TESTS::

        sage: for p in iet.Permutations_iterator(2, alphabet="ab"):
        ...       print p, "\n****"   #indirect doctest
        a b
        b a
        ****
        b a
        a b
        ****
        sage: for p in iet.Permutations_iterator(3, alphabet="abc"):
        ...       print p, "\n*****"   #indirect doctest
        a b c
        b c a
        *****
        a b c
        c a b
        *****
        a b c
        c b a
        *****
        a c b
        b a c
        *****
        a c b
        b c a
        *****
        a c b
        c b a
        *****
        b a c
        a c b
        *****
        b a c
        c a b
        *****
        b a c
        c b a
        *****
        b c a
        a b c
        *****
        b c a
        a c b
        *****
        b c a
        c a b
        *****
        c a b
        a b c
        *****
        c a b
        b a c
        *****
        c a b
        b c a
        *****
        c b a
        a b c
        *****
        c b a
        a c b
        *****
        c b a
        b a c
        *****
    """
    from itertools import imap, ifilter, product
    from sage.combinat.permutation import Permutations

    if irreducible is False:
        if nintervals is None:
            raise ValueError, "choose a number of intervals"
        else:
            assert(isinstance(nintervals,(int,Integer)))
            assert(nintervals > 0)

            f = lambda x: LabelledPermutationIET([list(x[0]),list(x[1])],alphabet=alphabet)

            alphabet = Alphabet(alphabet)
            g = lambda x: [alphabet.unrank(k-1) for k in x]
            P = map(g, Permutations(nintervals))
            return imap(f,product(P,P))
    else:
        return ifilter(
            lambda x: x.is_irreducible(),
            LabelledPermutationsIET_iterator(nintervals,False,alphabet))

class LabelledPermutationIET(LabelledPermutation, OrientablePermutationIET):
    """
    Labelled permutation for iet

    EXAMPLES::

        sage: from surface_dynamics.all import *

    Reducibility testing::

        sage: p = iet.Permutation('a b c', 'c b a')
        sage: p.is_irreducible()
        True

        sage: q = iet.Permutation('a b c d', 'b a d c')
        sage: q.is_irreducible()
        False

    Rauzy movability and Rauzy move::

        sage: p = iet.Permutation('a b c', 'c b a')
        sage: p.has_rauzy_move('top')
        True
        sage: print p.rauzy_move('bottom')
        a c b
        c b a
        sage: p.has_rauzy_move('top')
        True
        sage: print p.rauzy_move('top')
        a b c
        c a b

    Rauzy diagram::

        sage: p = iet.Permutation('a b c', 'c b a')
        sage: d = p.rauzy_diagram()
        sage: p in d
        True
    """
    def __cmp__(self, other):
        r"""
        ALGORITHM:

        The order is lexicographic on intervals[0] + intervals[1]

        TESTS::
            sage: list_of_p2 = []
            sage: p0 = iet.Permutation('1 2', '1 2')
            sage: p1 = iet.Permutation('1 2', '2 1')
            sage: p0 != p0
            False
            sage: (p0 == p0) and (p0 < p1)
            True
            sage: (p1 > p0) and (p1 == p1)
            True
        """
        if not isinstance(self, type(other)):
            return -1

        n = len(self)
        if n != len(other):
            return n - len(other)

        i, j = 0,0
        while (self._labels[i][j] == other._labels[i][j]):
            j += 1
            if j == n:
                if i == 1: return 0
                i = 1
                j = 0
        return self._labels[i][j] - other._labels[i][j]

    def reduced(self):
        r"""
        Returns the associated reduced abelian permutation.

        OUTPUT:

        a reduced permutation -- the underlying reduced permutation


        EXAMPLES:

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation("a b c d","d c a b")
            sage: q = iet.Permutation("a b c d","d c a b",reduced=True)
            sage: p.reduced() == q
            True
        """
        from reduced import ReducedPermutationIET

        return ReducedPermutationIET(self.list(),alphabet=self._alphabet)

    def rauzy_move_interval_substitution(self,winner=None,side=None):
        r"""
        Returns the interval substitution associated.

        INPUT:

        - ``winner`` - the winner interval ('top' or 'bottom')

        - ``side`` - (default: 'right') the side ('left' or 'right')

        OUTPUT:

        WordMorphism -- a substitution on the alphabet of the permutation

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b','b a')
            sage: print p.rauzy_move_interval_substitution('top','right')
            a->a, b->ba
            sage: print p.rauzy_move_interval_substitution('bottom','right')
            a->ab, b->b
            sage: print p.rauzy_move_interval_substitution('top','left')
            a->ba, b->b
            sage: print p.rauzy_move_interval_substitution('bottom','left')
            a->a, b->ab
        """
        d = dict([(letter,[letter]) for letter in self.letters()])

        if winner is None and side is None:
            return WordMorphism(d)

        winner = interval_conversion(winner)
        side = side_conversion(side)

        winner_letter = self.rauzy_move_winner(winner,side)
        loser_letter = self.rauzy_move_loser(winner,side)

        if side == 0:
            d[winner_letter] = [loser_letter,winner_letter]
        else:
            d[winner_letter] = [winner_letter,loser_letter]

        return WordMorphism(d)

    def rauzy_move_orbit_substitution(self,winner=None,side=None):
        r"""
        Return the action fo the rauzy_move on the orbit.

        INPUT:

        - ``i`` - integer

        - ``winner`` - the winner interval ('top' or 'bottom')

        - ``side`` - (default: 'right') the side ('right' or 'left')

        OUTPUT:

        WordMorphism -- a substitution on the alphabet of self

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b','b a')
            sage: print p.rauzy_move_orbit_substitution('top','right')
            a->ab, b->b
            sage: print p.rauzy_move_orbit_substitution('bottom','right')
            a->a, b->ab
            sage: print p.rauzy_move_orbit_substitution('top','left')
            a->a, b->ba
            sage: print p.rauzy_move_orbit_substitution('bottom','left')
            a->ba, b->b

        TESTS::

            sage: p = iet.Permutation('a1 a2', 'a2 a1')
            sage: p.rauzy_move_orbit_substitution('top','right').codomain().alphabet()
            {'a1', 'a2'}
        """
        d = dict([(letter,[letter]) for letter in self.letters()])

        if winner is None and side is None:
            return WordMorphism(d)

        winner = interval_conversion(winner)
        side = side_conversion(side)

        loser_letter = self.rauzy_move_loser(winner,side)

        top_letter = self.alphabet().unrank(self._labels[0][side])
        bottom_letter = self.alphabet().unrank(self._labels[1][side])

        d[loser_letter] = [bottom_letter,top_letter]

        return WordMorphism(d)

    def rauzy_diagram(self, **args):
        """
        Returns the associated Rauzy diagram.

        For more information try help(iet.RauzyDiagram).

        OUTPUT:

        Rauzy diagram -- the Rauzy diagram of the permutation

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c', 'c b a')
            sage: d = p.rauzy_diagram()
        """
        return LabelledRauzyDiagram(self, **args)

class LabelledPermutationLI(LabelledPermutation, OrientablePermutationLI):
    r"""
    Labelled quadratic (or generalized) permutation

    EXAMPLES::

            sage: from surface_dynamics.all import *

    Reducibility testing::

        sage: p = iet.GeneralizedPermutation('a b b', 'c c a')
        sage: p.is_irreducible()
        True

    Reducibility testing with associated decomposition::

        sage: p = iet.GeneralizedPermutation('a b c a', 'b d d c')
        sage: p.is_irreducible()
        False
        sage: test, decomposition = p.is_irreducible(return_decomposition = True)
        sage: print test
        False
        sage: print decomposition
        (['a'], ['c', 'a'], [], ['c'])

    Rauzy movability and Rauzy move::

        sage: p = iet.GeneralizedPermutation('a a b b c c', 'd d')
        sage: p.has_rauzy_move(0)
        False
        sage: p.has_rauzy_move(1)
        True
        sage: q = p.rauzy_move(1)
        sage: print q
        a a b b c
        c d d
        sage: q.has_rauzy_move(0)
        True
        sage: q.has_rauzy_move(1)
        True

    Rauzy diagrams::

        sage: p = iet.GeneralizedPermutation('0 0 1 1','2 2')
        sage: r = p.rauzy_diagram()
        sage: p in r
        True
    """
    def __cmp__(self, other):
        r"""
        ALGORITHM:

        Order is lexicographic on length of intervals and on intervals.

        TESTS::
            sage: p0 = iet.GeneralizedPermutation('0 0','1 1 2 2')
            sage: p1 = iet.GeneralizedPermutation('0 0','1 2 1 2')
            sage: p2 = iet.GeneralizedPermutation('0 0','1 2 2 1')
            sage: p3 = iet.GeneralizedPermutation('0 0 1','1 2 2')
            sage: p4 = iet.GeneralizedPermutation('0 0 1 1','2 2')
            sage: p5 = iet.GeneralizedPermutation('0 1 0 1','2 2')
            sage: p6 = iet.GeneralizedPermutation('0 1 1 0','2 2')
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
        """
        if not isinstance(self, type(other)):
            return -1

        n = len(self)

        if n != len(other): return n - len(other)

        l0 = self._labels[0]
        l1 = other._labels[0]

        n = len(self._labels[0])

        if n != len(other._labels[0]): return n - len(other._labels[0])

        i = 0
        while (i < n) and (l0[i] == l1[i]):
            i += 1

        if i != n:
            return l0[i] - l1[i]

        l0 = self._labels[1]
        l1 = other._labels[1]
        n = len(self._labels[1])

        i = 0
        while (i < n) and (l0[i] == l1[i]):
            i += 1

        if i != n:
            return l0[i] - l1[i]

        return 0

    def has_right_rauzy_move(self, winner):
        r"""
        Test of Rauzy movability with a specified winner)

        A quadratic (or generalized) permutation is rauzy_movable type
        depending on the possible length of the last interval. It's
        dependent of the length equation.

        INPUT:

        - ``winner`` - 'top' (or 't' or 0) or 'bottom' (or 'b' or 1)

        OUTPUT:

        bool -- True if self has a Rauzy move

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('a a','b b')
            sage: p.has_right_rauzy_move('top')
            False
            sage: p.has_right_rauzy_move('bottom')
            False

        ::

            sage: p = iet.GeneralizedPermutation('a a b','b c c')
            sage: p.has_right_rauzy_move('top')
            True
            sage: p.has_right_rauzy_move('bottom')
            True

        ::

            sage: p = iet.GeneralizedPermutation('a a','b b c c')
            sage: p.has_right_rauzy_move('top')
            True
            sage: p.has_right_rauzy_move('bottom')
            False

        ::

            sage: p = iet.GeneralizedPermutation('a a b b','c c')
            sage: p.has_right_rauzy_move('top')
            False
            sage: p.has_right_rauzy_move('bottom')
            True
        """
        winner = interval_conversion(winner)
        loser = self._labels[1-winner][-1]

        # the same letter at the right-end (False)
        if self._labels[0][-1] == self._labels[1][-1] :
            return False

        # the winner (or loser) letter is repeated on the other interval (True)
        if self._labels[0][-1] in self._labels[1]: return True
        if self._labels[1][-1] in self._labels[0]: return True

        # the loser letters is the only letter repeated in the loser
        # interval (False)
        for i,c in enumerate((self._labels[1-winner])):
            if c != loser and c in self._labels[1-winner][i+1:]:
                return True

        return False

    def right_rauzy_move(self, winner):
        r"""
        Perform a Rauzy move on the right (the standard one).

        INPUT:

        - ``winner`` - 'top' (or 't' or 0) or 'bottom' (or 'b' or 1)

        OUTPUT:

        boolean -- True if self has a Rauzy move

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('a a b','b c c')
            sage: p.right_rauzy_move(0)
            a a b
            b c c
            sage: p.right_rauzy_move(1)
            a a
            b b c c

        ::

            sage: p = iet.GeneralizedPermutation('a b b','c c a')
            sage: p.right_rauzy_move(0)
            a a b b
            c c
            sage: p.right_rauzy_move(1)
            a b b
            c c a

        TESTS::

            sage: p = iet.GeneralizedPermutation('a a b','b c c')
            sage: q = p.top_bottom_inverse()
            sage: q = q.right_rauzy_move(0)
            sage: q = q.top_bottom_inverse()
            sage: q == p.right_rauzy_move(1)
            True
            sage: q = p.top_bottom_inverse()
            sage: q = q.right_rauzy_move(1)
            sage: q = q.top_bottom_inverse()
            sage: q == p.right_rauzy_move(0)
            True
            sage: p = p.left_right_inverse()
            sage: q = q.left_rauzy_move(0)
            sage: q = q.left_right_inverse()
            sage: q == p.right_rauzy_move(0)
            True
            sage: q = p.left_right_inverse()
            sage: q = q.left_rauzy_move(1)
            sage: q = q.left_right_inverse()
            sage: q == p.right_rauzy_move(1)
            True
        """
        result = copy(self)

        winner_letter = result._labels[winner][-1]
        loser_letter = result._labels[1-winner].pop(-1)

        if winner_letter in result._labels[winner][:-1]:
            loser_to = result._labels[winner].index(winner_letter)
            result._labels[winner].insert(loser_to, loser_letter)
        else:
            loser_to = result._labels[1-winner].index(winner_letter) + 1
            result._labels[1-winner].insert(loser_to, loser_letter)

        return result

    def left_rauzy_move(self, winner):
        r"""
        Perform a Rauzy move on the left.

        INPUT:

        - ``winner`` - 'top' or 'bottom'

        OUTPUT:

        permutation -- the Rauzy move of self

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('a a b','b c c')
            sage: p.left_rauzy_move(0)
            a a b b
            c c
            sage: p.left_rauzy_move(1)
            a a b
            b c c

        ::

            sage: p = iet.GeneralizedPermutation('a b b','c c a')
            sage: p.left_rauzy_move(0)
            a b b
            c c a
            sage: p.left_rauzy_move(1)
            b b
            c c a a


        TESTS::

            sage: p = iet.GeneralizedPermutation('a a b','b c c')
            sage: q = p.top_bottom_inverse()
            sage: q = q.left_rauzy_move(0)
            sage: q = q.top_bottom_inverse()
            sage: q == p.left_rauzy_move(1)
            True
            sage: q = p.top_bottom_inverse()
            sage: q = q.left_rauzy_move(1)
            sage: q = q.top_bottom_inverse()
            sage: q == p.left_rauzy_move(0)
            True
            sage: q = p.left_right_inverse()
            sage: q = q.right_rauzy_move(0)
            sage: q = q.left_right_inverse()
            sage: q == p.left_rauzy_move(0)
            True
            sage: q = p.left_right_inverse()
            sage: q = q.right_rauzy_move(1)
            sage: q = q.left_right_inverse()
            sage: q == p.left_rauzy_move(1)
            True
        """
        result = copy(self)

        winner_letter = result._labels[winner][0]
        loser_letter = result._labels[1-winner].pop(0)

        if winner_letter in result._labels[winner][1:]:
            loser_to = result._labels[winner][1:].index(winner_letter)+2
            result._labels[winner].insert(loser_to, loser_letter)

        else:
            loser_to = result._labels[1-winner].index(winner_letter)
            result._labels[1-winner].insert(loser_to, loser_letter)

        return result

    def reduced(self):
        r"""
        Returns the associated reduced quadratic permutations.

        OUTPUT:

        permutation -- the underlying reduced permutation

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('a a','b b c c')
            sage: q = p.reduced()
            sage: q
            a a
            b b c c
            sage: p.rauzy_move(0).reduced() == q.rauzy_move(0)
            True
        """
        from reduced import ReducedPermutationLI

        return ReducedPermutationLI(self.list(),alphabet=self._alphabet)

    def rauzy_diagram(self, **kargs):
        r"""
        Returns the associated RauzyDiagram.

        OUTPUT:

        Rauzy diagram -- the Rauzy diagram of the permutation

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('a b c b', 'c d d a')
            sage: d = p.rauzy_diagram()
            sage: p in d
            True

        For more information, try help(iet.RauzyDiagram)
        """
        return LabelledRauzyDiagram(self, **kargs)

class FlippedLabelledPermutation(LabelledPermutation):
    r"""
    General template for labelled objects

    .. warning::
       Internal class! Do not use directly!
    """
    def __init__(self, intervals=None, alphabet=None, flips=None):
        r"""
        INPUT:

        - `intervals` - the intervals as a list of two lists

        - `alphabet` - something that should be converted to an alphabe

        - `flips` - a list of letters of the alphabet

        TESTS:

        ::

            sage: from surface_dynamics.interval_exchanges.labelled import FlippedLabelledPermutationIET
            sage: p = FlippedLabelledPermutationIET([['a','b'],['a','b']],flips='a')
            sage: p == loads(dumps(p))
            True
            sage: p = FlippedLabelledPermutationIET([['a','b'],['b','a']],flips='ab')
            sage: p == loads(dumps(p))
            True

        ::

            sage: from surface_dynamics.interval_exchanges.labelled import FlippedLabelledPermutationLI
            sage: p = FlippedLabelledPermutationLI([['a','a','b'],['b','c','c']],flips='a')
            sage: p == loads(dumps(p))
            True
            sage: p = FlippedLabelledPermutationLI([['a','a'],['b','b','c','c']],flips='ac')
            sage: p == loads(dumps(p))
            True
        """
        if intervals is None: intervals=[[],[]]
        if flips is None: flips = []

        super(FlippedLabelledPermutation, self).__init__(intervals, alphabet)
        self._init_flips(intervals, flips)

    def __eq__(self,other):
        r"""
        Test of equality

        ALGORITHM:

        not considering the alphabet used for the representation but just the
        order

        TESTS::

            sage: p1 = iet.Permutation('a b c','c b a',flips='a')
            sage: p2 = iet.Permutation('a b c','c b a',flips='b')
            sage: p3 = iet.Permutation('d e f','f e d',flips='d')
            sage: p1 == p1 and p2 == p2 and p3 == p3
            True
            sage: p1 == p2
            False
            sage: p1 == p3
            True
        """
        return (
            type(self) == type(other) and
            self._labels == other._labels and
            self._flips == other._flips)

    def __ne__(self,other):
        r"""
        Test of difference

        ALGORITHM:

        not considering the alphabet used for the representation

        TESTS::

            sage: p1 = iet.Permutation('a b c','c b a',flips='a')
            sage: p2 = iet.Permutation('a b c','c b a',flips='b')
            sage: p3 = iet.Permutation('d e f','f e d',flips='d')
            sage: p1 != p1 or p2 != p2 or p3 != p3
            False
            sage: p1 != p2
            True
            sage: p1 != p3
            False
        """
        return (
            type(self) != type(other) or
            self._labels != other._labels or
            self._flips != other._flips)

    def list(self, flips=False):
        r"""
        Returns a list associated to the permutation.

        INPUT:

        - ``flips`` - boolean (default: False)

        OUTPUT:

        list -- two lists of labels

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('0 0 1 2 2 1', '3 3', flips='1')
            sage: p.list(flips=True)
            [[('0', 1), ('0', 1), ('1', -1), ('2', 1), ('2', 1), ('1', -1)], [('3', 1), ('3', 1)]]
            sage: p.list(flips=False)
            [['0', '0', '1', '2', '2', '1'], ['3', '3']]

        The list can be used to reconstruct the permutation

        ::

            sage: p = iet.Permutation('a b c','c b a',flips='ab')
            sage: p == iet.Permutation(p.list(), flips=p.flips())
            True

        ::

            sage: p = iet.GeneralizedPermutation('a b b c','c d d a',flips='ad')
            sage: p == iet.GeneralizedPermutation(p.list(),flips=p.flips())
            True
        """
        if flips:
            a0 = zip(map(self._alphabet.unrank, self._labels[0]), self._flips[0])
            a1 = zip(map(self._alphabet.unrank, self._labels[1]), self._flips[1])
        else:
            a0 = map(self._alphabet.unrank, self._labels[0])
            a1 = map(self._alphabet.unrank, self._labels[1])

        return [a0,a1]

    def __getitem__(self,i):
        r"""
        Get labels and flips of specified interval.

        The result is a 2-uple (letter, flip) where letter is the name of the
        sub-interval and flip is a number corresponding to the presence of flip
        as following: 1 (no flip) and -1 (a flip).

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b', 'b a', flips='a')
            sage: print p[0]
            [('a', -1), ('b', 1)]
            sage: p = iet.GeneralizedPermutation('c p p', 't t c', flips='ct')
            sage: print p[1]
            [('t', -1), ('t', -1), ('c', -1)]
        """
        if not isinstance(i, (Integer, int)):
            raise TypeError("Must be an integer")
        if i != 0 and i != 1:
            raise IndexError("The integer must be 0 or 1")

        letters = map(self._alphabet.unrank, self._labels[i])
        flips = self._flips[i]

        return zip(letters,flips)

class FlippedLabelledPermutationIET(
    FlippedLabelledPermutation,
    FlippedPermutationIET,
    LabelledPermutationIET):
    r"""
    Flipped labelled permutation from iet.

    EXAMPLES:

            sage: from surface_dynamics.all import *

    Reducibility testing (does not depends of flips)::

        sage: p = iet.Permutation('a b c', 'c b a',flips='a')
        sage: p.is_irreducible()
        True
        sage: q = iet.Permutation('a b c d', 'b a d c', flips='bc')
        sage: q.is_irreducible()
        False

    Rauzy movability and Rauzy move::

        sage: p = iet.Permutation('a b c', 'c b a',flips='a')
        sage: print p
        -a  b  c
         c  b -a
        sage: print p.rauzy_move(1)
        -c -a  b
        -c  b -a
        sage: print p.rauzy_move(0)
        -a  b  c
         c -a  b

    Rauzy diagrams::

        sage: d = iet.RauzyDiagram('a b c d','d a b c',flips='a')

    AUTHORS:

    - Vincent Delecroix (2009-09-29): initial version
    """
    def reduced(self):
        r"""
        The associated reduced permutation.

        OUTPUT:

        permutation -- the associated reduced permutation

        EXAMPLE::

            sage: p = iet.Permutation('a b c','c b a',flips='a')
            sage: q = iet.Permutation('a b c','c b a',flips='a',reduced=True)
            sage: p.reduced() == q
            True
        """
        from surface_dynamics.interval_exchanges.reduced import FlippedReducedPermutationIET

        return FlippedReducedPermutationIET(
            intervals=self.list(flips=False),
            flips=self.flips(),
            alphabet=self.alphabet())

    def __hash__(self):
        r"""
        ALGORITHM:

        Uses hash of string

        TESTS::

            sage: p =[]
            sage: p.append(iet.Permutation('a b','a b',flips='a'))
            sage: p.append(iet.Permutation('a b','a b',flips='b'))
            sage: p.append(iet.Permutation('a b','a b',flips='ab'))
            sage: p.append(iet.Permutation('a b','b a',flips='a'))
            sage: p.append(iet.Permutation('a b','b a',flips='b'))
            sage: p.append(iet.Permutation('a b','b a',flips='ab'))
            sage: h = map(hash, p)
            sage: for i in range(len(h)-1):
            ...      if h[i] == h[i+1]:
            ...          print "You choose a bad hash!"
        """
        if self._hash is None:
            f = self._flips
            i = self._labels
            l = []
            l.extend([str(j*(1+k)) for j,k in zip(f[0],i[0])])
            l.extend([str(-j*(1+k)) for j,k in zip(f[1],i[1])])
            self._hash = hash(''.join(l))

        return self._hash

    def rauzy_diagram(self, **kargs):
        r"""
        Returns the Rauzy diagram associated to this permutation.

        For more information, try help(iet.RauzyDiagram)

        OUTPUT:

        RauzyDiagram -- the Rauzy diagram of self

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c', 'c b a',flips='a')
            sage: p.rauzy_diagram()
            Rauzy diagram with 3 permutations
        """
        return FlippedLabelledRauzyDiagram(self, **kargs)

class FlippedLabelledPermutationLI(
    FlippedLabelledPermutation,
    FlippedPermutationLI,
    LabelledPermutationLI):
    r"""
    Flipped labelled quadratic (or generalized) permutation.

    EXAMPLES:

            sage: from surface_dynamics.all import *

    Reducibility testing::

        sage: p = iet.GeneralizedPermutation('a b b', 'c c a', flips='a')
        sage: p.is_irreducible()
        True

    Reducibility testing with associated decomposition::

        sage: p = iet.GeneralizedPermutation('a b c a', 'b d d c', flips='ab')
        sage: p.is_irreducible()
        False
        sage: test, decomp = p.is_irreducible(return_decomposition = True)
        sage: print test
        False
        sage: print decomp
        (['a'], ['c', 'a'], [], ['c'])

    Rauzy movability and Rauzy move::

        sage: p = iet.GeneralizedPermutation('a a b b c c', 'd d', flips='d')
        sage: p.has_rauzy_move(0)
        False
        sage: p.has_rauzy_move(1)
        True
        sage: p = iet.GeneralizedPermutation('a a b','b c c',flips='c')
        sage: p.has_rauzy_move(0)
        True
        sage: p.has_rauzy_move(1)
        True
    """
    def reduced(self):
        r"""
        The associated reduced permutation.

        OUTPUT:

        permutation -- the associated reduced permutation

        EXAMPLE::

            sage: p = iet.GeneralizedPermutation('a a','b b c c',flips='a')
            sage: q = iet.GeneralizedPermutation('a a','b b c c',flips='a',reduced=True)
            sage: p.reduced() == q
            True
        """
        from surface_dynamics.interval_exchanges.reduced import FlippedReducedPermutationLI

        return FlippedReducedPermutationLI(
            intervals=self.list(flips=False),
            flips=self.flips(),
            alphabet=self.alphabet())

    def right_rauzy_move(self, winner):
        r"""
        Perform a Rauzy move on the right (the standard one).

        INPUT:

        - ``winner`` - either 'top' or 'bottom' ('t' or 'b' for short)

        OUTPUT:

        permutation -- the Rauzy move of self

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('a a b','b c c',flips='c')
            sage: p.right_rauzy_move(0)
             a  a  b
            -c  b -c
            sage: p.right_rauzy_move(1)
             a  a
            -b -c -b -c

        ::

            sage: p = iet.GeneralizedPermutation('a b b','c c a',flips='ab')
            sage: p.right_rauzy_move(0)
             a -b  a -b
             c  c
            sage: p.right_rauzy_move(1)
             b -a  b
             c  c -a
        """
        result = copy(self)

        winner_letter = result._labels[winner][-1]
        winner_flip = result._flips[winner][-1]

        loser_letter = result._labels[1-winner].pop(-1)
        loser_flip = result._flips[1-winner].pop(-1)

        if loser_letter in result._labels[winner]:
            loser_twin = result._labels[winner].index(loser_letter)
            result._flips[winner][loser_twin] = loser_flip*winner_flip
        else:
            loser_twin = result._labels[1-winner].index(loser_letter)
            result._flips[1-winner][loser_twin] = loser_flip*winner_flip

        if winner_letter in result._labels[winner][:-1]:
            loser_to = result._labels[winner].index(winner_letter)
            if winner_flip == -1: loser_to += 1
            result._labels[winner].insert(loser_to, loser_letter)
            result._flips[winner].insert(loser_to, loser_flip*winner_flip)
        else:
            loser_to = result._labels[1-winner].index(winner_letter)
            if loser_flip == 1: loser_to += 1
            result._labels[1-winner].insert(loser_to, loser_letter)
            result._flips[1-winner].insert(loser_to, loser_flip*winner_flip)

        return result

    def left_rauzy_move(self, winner):
        r"""
        Perform a Rauzy move on the left.

        INPUT:

        - ``winner`` - either 'top' or 'bottom' ('t' or 'b' for short)

        OUTPUT:

        -- a permutation

        EXAMPLES:

            sage: from surface_dynamics.all import *

        ::

            sage: p = iet.GeneralizedPermutation('a a b','b c c')
            sage: p.left_rauzy_move(0)
            a a b b
            c c
            sage: p.left_rauzy_move(1)
            a a b
            b c c

        ::

            sage: p = iet.GeneralizedPermutation('a b b','c c a')
            sage: p.left_rauzy_move(0)
            a b b
            c c a
            sage: p.left_rauzy_move(1)
            b b
            c c a a
        """
        result = copy(self)

        winner_letter = result._labels[winner][0]
        loser_letter = result._labels[1-winner].pop(0)

        if winner_letter in result._labels[winner][1:]:
            loser_to = result._labels[winner][1:].index(winner_letter)+2
            result._labels[winner].insert(loser_to, loser_letter)

        else:
            loser_to = result._labels[1-winner].index(winner_letter)
            result._labels[1-winner].insert(loser_to, loser_letter)

        return result

    def rauzy_diagram(self, **kargs):
        r"""
        Returns the associated Rauzy diagram.

        For more information, try help(RauzyDiagram)

        OUTPUT :

        -- a RauzyDiagram

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.GeneralizedPermutation('a b b a', 'c d c d')
            sage: d = p.rauzy_diagram()
        """
        return FlippedLabelledRauzyDiagram(self, **kargs)

class LabelledRauzyDiagram(RauzyDiagram):
    r"""
    Template for Rauzy diagrams of labelled permutations.

        ...DO NOT USE...
    """
    class Path(RauzyDiagram.Path):
        r"""
        Path in Labelled Rauzy diagram.
        """
        def matrix(self):
            r"""
            Returns the matrix associated to a path.

            The matrix associated to a Rauzy induction, is the linear
            application that allows to recover the lengths of self from the
            lengths of the induced.

            OUTPUT:

            matrix -- a square matrix of integers

            EXAMPLES::

                sage: from surface_dynamics.all import *

                sage: p = iet.Permutation('a1 a2','a2 a1')
                sage: d = p.rauzy_diagram()
                sage: g = d.path(p,'top')
                sage: g.matrix()
                [1 0]
                [1 1]
                sage: g = d.path(p,'bottom')
                sage: g.matrix()
                [1 1]
                [0 1]

            ::

                sage: p = iet.Permutation('a b c','c b a')
                sage: d = p.rauzy_diagram()
                sage: g = d.path(p)
                sage: g.matrix() == identity_matrix(3)
                True
                sage: g = d.path(p,'top')
                sage: g.matrix()
                [1 0 0]
                [0 1 0]
                [1 0 1]
                sage: g = d.path(p,'bottom')
                sage: g.matrix()
                [1 0 1]
                [0 1 0]
                [0 0 1]
            """
            return self.composition(self._parent.edge_to_matrix)

        def interval_substitution(self):
            r"""
            Returns the substitution of intervals obtained.

            OUTPUT:

            WordMorphism -- the word morphism corresponding to the interval

            EXAMPLES::

                sage: from surface_dynamics.all import *

                sage: p = iet.Permutation('a b','b a')
                sage: r = p.rauzy_diagram()
                sage: p0 = r.path(p,0)
                sage: s0 = p0.interval_substitution()
                sage: print s0
                a->a, b->ba
                sage: p1 = r.path(p,1)
                sage: s1 = p1.interval_substitution()
                sage: print s1
                a->ab, b->b
                sage: (p0 + p1).interval_substitution() == s1 * s0
                True
                sage: (p1 + p0).interval_substitution() == s0 * s1
                True
            """
            return self.right_composition(self._parent.edge_to_interval_substitution)

        def orbit_substitution(self):
            r"""
            Returns the substitution on the orbit of the left extremity.

            OUTPUT:

            WordMorhpism -- the word morphism corresponding to the orbit

            EXAMPLES::

                sage: from surface_dynamics.all import *

                sage: p = iet.Permutation('a b','b a')
                sage: d = p.rauzy_diagram()
                sage: g0 = d.path(p,'top')
                sage: s0 = g0.orbit_substitution()
                sage: print s0
                a->ab, b->b
                sage: g1 = d.path(p,'bottom')
                sage: s1 = g1.orbit_substitution()
                sage: print s1
                a->a, b->ab
                sage: (g0 + g1).orbit_substitution() == s0 * s1
                True
                sage: (g1 + g0).orbit_substitution() == s1 * s0
                True
            """
            return self.composition(self._parent.edge_to_orbit_substitution)

        substitution = orbit_substitution  # standard name
        dual_substitution = interval_substitution  # standard name

        def is_full(self):
            r"""
            Tests the fullness.

            A path is full if all intervals win at least one time.

            OUTPUT:

            boolean -- True if the path is full and False else

            EXAMPLE::

                sage: p = iet.Permutation('a b c','c b a')
                sage: r = p.rauzy_diagram()
                sage: g0 = r.path(p,'t','b','t')
                sage: g1 = r.path(p,'b','t','b')
                sage: g0.is_full()
                False
                sage: g1.is_full()
                False
                sage: (g0 + g1).is_full()
                True
                sage: (g1 + g0).is_full()
                True
            """
            return set(self._parent.letters()) == set(self.winners())

    def edge_to_interval_substitution(self, p=None, edge_type=None):
        r"""
        Returns the interval substitution associated to an edge

        OUTPUT:

        WordMorphism -- the WordMorphism corresponding to the edge

        EXAMPLE::

            sage: p = iet.Permutation('a b c','c b a')
            sage: r = p.rauzy_diagram()
            sage: print r.edge_to_interval_substitution(None,None)
            a->a, b->b, c->c
            sage: print r.edge_to_interval_substitution(p,0)
            a->a, b->b, c->ca
            sage: print r.edge_to_interval_substitution(p,1)
            a->ac, b->b, c->c
        """
        if p is None and edge_type is None:
            return WordMorphism(dict((a,[a]) for a in self.letters()))

        function_name = self._edge_types[edge_type][0] + '_interval_substitution'
        if not hasattr(self._element_class,function_name):
            return WordMorphism(dict((a,[a]) for a in self.letters()))

        arguments = self._edge_types[edge_type][1]

        return getattr(p,function_name)(*arguments)

    def edge_to_orbit_substitution(self, p=None, edge_type=None):
        r"""
        Returns the interval substitution associated to an edge

        OUTPUT:

        WordMorphism -- the word morphism corresponding to the edge

        EXAMPLE::

            sage: p = iet.Permutation('a b c','c b a')
            sage: r = p.rauzy_diagram()
            sage: print r.edge_to_orbit_substitution(None,None)
            a->a, b->b, c->c
            sage: print r.edge_to_orbit_substitution(p,0)
            a->ac, b->b, c->c
            sage: print r.edge_to_orbit_substitution(p,1)
            a->a, b->b, c->ac

        TESTS::

            sage: from surface_dynamics.all import *

            sage: pi0 = iet.Permutation('A1 A2 B', 'B A1 A2')
            sage: G = pi0.rauzy_diagram()
            sage: s1 = G.edge_to_orbit_substitution(pi0,0)
            sage: print s1.domain().alphabet()
            {'A1', 'A2', 'B'}
            sage: print s1.codomain().alphabet()
            {'A1', 'A2', 'B'}
        """
        if p is None and edge_type is None:
            return WordMorphism(dict((a,[a]) for a in self.letters()))

        function_name = self._edge_types[edge_type][0] + '_orbit_substitution'

        if not hasattr(self._element_class,function_name):
            return WordMorphism(dict((a,[a]) for a in self.letters()))

        arguments = self._edge_types[edge_type][1]
        return getattr(p,function_name)(*arguments)

    def full_loop_iterator(self, start=None, max_length=1):
        r"""
        Returns an iterator over all full path starting at start.

        INPUT:

        - ``start`` - the start point

        - ``max_length`` - a limit on the length of the paths

        OUTPUT:

        iterator -- iterator over full loops

        EXAMPLE::

            sage: p = iet.Permutation('a b','b a')
            sage: r = p.rauzy_diagram()
            sage: for g in r.full_loop_iterator(p,2):
            ...       print g.matrix(), "\n*****"
            [1 1]
            [1 2]
            *****
            [2 1]
            [1 1]
            *****
        """
        from itertools import ifilter, imap

        g = self.path(start)

        ifull = ifilter(
            lambda x: x.is_loop() and x.is_full(),
            self._all_path_extension(g,max_length))

        return imap(copy,ifull)

    def full_nloop_iterator(self, start=None, length=1):
        r"""
        Returns an iterator over all full loops of given length.

        INPUT:

        - ``start`` - the initial permutation

        - ``length`` - the length to consider

        OUTPUT:

        iterator -- an iterator over the full loops of given length

        EXAMPLE::

            sage: p = iet.Permutation('a b','b a')
            sage: d = p.rauzy_diagram()
            sage: for g in d.full_nloop_iterator(p,2):
            ...       print g.matrix(), "\n*****"
            [1 1]
            [1 2]
            *****
            [2 1]
            [1 1]
            *****
        """
        from itertools import ifilter, imap

        g = self.path(start)

        ifull = ifilter(
            lambda x: x.is_loop() and x.is_full(),
            self._all_npath_extension(g,length))

        return imap(copy, ifull)

    def _permutation_to_vertex(self, p):
        r"""
        Translation of a labelled permutation to a vertex

        INPUT:

        - ``p`` - a labelled Permutation

        TESTS::

            sage: p = iet.Permutation('a b c','c b a')
            sage: r = p.rauzy_diagram()
            sage: p in r   #indirect doctest
            True
        """
        return (
        tuple(p._labels[0]),tuple(p._labels[1]),
        tuple(p._twin[0]),tuple(p._twin[1]))

    def _set_element(self,data):
        r"""
        Sets self._element with data

        TESTS::

            sage: p = iet.Permutation('a b c','c b a')
            sage: r = p.rauzy_diagram()
            sage: r[p][0] == p.rauzy_move(0)   #indirect doctest
            True
            sage: r[p][1] == p.rauzy_move(1)   #indirect doctest
            True
        """
        self._element._labels = [list(data[0]), list(data[1])]
        self._element._twin = [list(data[2]), list(data[3])]

class FlippedLabelledRauzyDiagram(FlippedRauzyDiagram, LabelledRauzyDiagram):
    r"""
    Rauzy diagram of flipped labelled permutations
    """
    def _permutation_to_vertex(self, p):
        r"""
        Returns what must be stored from p.

        INPUT:

        - ``p`` - a Flipped labelled permutation

        TESTS::

            sage: p = iet.Permutation('a b c','c b a',flips='a')
            sage: r = p.rauzy_diagram()
            sage: p in r   #indirect doctest
            True
        """
        return ((tuple(p._labels[0]),tuple(p._labels[1])),
                (tuple(p._twin[0]), tuple(p._twin[1])),
                (tuple(p._flips[0]), tuple(p._flips[1])))

    def _set_element(self, data):
        r"""
        Returns what the vertex i as a permutation.

        TESTS::

            sage: p = iet.Permutation('a b','b a',flips='a')
            sage: r = p.rauzy_diagram()
            sage: p in r   #indirect doctest
            True
        """
        self._element._labels = [list(data[0][0]), list(data[0][1])]
        self._element._twin = [list(data[1][0]), list(data[1][1])]
        self._element._flips = [list(data[2][0]), list(data[2][1])]

