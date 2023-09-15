# encoding=utf-8
r"""
Separatrix diagrams and cylinder diagrams

A separatrix diagram is a couple of permutation ``(bot,top)`` that have the same
number of cycles in their cycle decompositions. A cylinder diagram is a
separatrix diagram together with a bijection between the cycles of ``bot`` and
``top``.

A cylinder diagram encodes the combinatorics of cylinder decomposition of a
completely periodic direction in a translation surface. If we adjoin coordinates
to this combinatorial datum, we have a complete description of the underlying
surface. In the case of arithmetic curves, the coordinates can be taken to be
rational numbers.

This representation of a surface is used in various constructions:

- square tiled surfaces

- Thurston-Veech construction of pseudo-Anosov diffeomorphism

- description of the cusp of Teichmueller curves

.. TODO::

    - We need a more general structure to encode configurations of structure of
      saddle connections (which need not be completely periodic directions (see
      [EskMasZor03]_)

    - Gray code for conjugacy classes of permutation in order to optimize the
      generation of separatrix and cylinder diagrams.

    - Implement separatrix diagram as a list of ribbon graphs with a pairing
      of its faces. In the Abelian case, the maps are endowed with bicoloration
      of the faces and the pairing is a pair of faces with different colors.
      (alternative for quadratic case: orientable + involution)

    - enumerators for fat graphs of genus g with prescribed vertex degrees
      (possibly up to isomorphism)

    - enumerators for fat graphs of genus g with prescribed vertex degrees
      (possibly up to isomorphism)

EXAMPLES::

    sage: from surface_dynamics import *

Separatrix diagrams::

    sage: s = SeparatrixDiagram('(0,1,2)(3,4)(5,6,7)','(0,4,1,2)(3,7)(5,6)')
    sage: s
    (0,1,2)(3,4)(5,6,7)-(0,4,1,2)(3,7)(5,6)
    sage: s.bot_cycle_tuples()
    [(0, 1, 2), (3, 4), (5, 6, 7)]
    sage: s.top_cycle_tuples()
    [(0, 4, 1, 2), (3, 7), (5, 6)]

Cylinder diagrams::

    sage: c = CylinderDiagram([((0,),(4,)),((1,2),(0,1,3)),((3,4),(2,))])
    sage: print(c)
    (0)-(4) (1,2)-(0,1,3) (3,4)-(2)
    sage: print(c.separatrix_diagram())
    (0)(1,2)(3,4)-(0,1,3)(2)(4)

They can also be built from separatrix diagram::

    sage: s = SeparatrixDiagram('(0,1,2)(3,4)(5,6,7)','(0,4,1,2)(3,7)(5,6)')
    sage: s
    (0,1,2)(3,4)(5,6,7)-(0,4,1,2)(3,7)(5,6)
    sage: s.to_cylinder_diagram([(0,1),(1,0),(2,2)])
    (0,1,2)-(3,7) (3,4)-(0,4,1,2) (5,6,7)-(5,6)
"""
#*****************************************************************************
#       Copyright (C) 2010-2019 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function, absolute_import
from six.moves import range, map, filter, zip
from six import iteritems

from functools import total_ordering

import collections
import numbers

from sage.structure.sage_object import SageObject

import itertools
import sage.arith.misc as arith
from sage.rings.integer import Integer

from sage.graphs.digraph import DiGraph

from surface_dynamics.misc.permutation import (perm_check, equalize_perms, perm_init,
        perm_cycles, perm_cycle_type, perm_compose_i,
        perm_invert, perms_canonical_labels,
        perms_transitive_components, canonical_perm, canonical_perm_i,
        argmin)
from surface_dynamics.misc.linalg import cone_triangulate

from sage.misc.decorators import rename_keyword

#
# Abelian and quadratic Separatrix Diagram
#
def two_non_connected_perms_canonical_labels(bot, top):
    r"""
    EXAMPLES::

        sage: from surface_dynamics.flat_surfaces.separatrix_diagram import two_non_connected_perms_canonical_labels
        sage: two_non_connected_perms_canonical_labels([3,2,1,0],[0,1,2,3])
        ([1, 0, 3, 2], [0, 1, 2, 3])
    """
    n = len(bot)

    cs_type_nb = {} # (bot,top) -> nb of them
    c_inv = [None] * n

    for c in perms_transitive_components([bot,top]):
        for i,j in enumerate(c):
            c_inv[j] = i
        cbot = [None]*len(c)
        ctop = [None]*len(c)
        for i in c:
            cbot[c_inv[i]] = c_inv[bot[i]]
            ctop[c_inv[i]] = c_inv[top[i]]

        (cbot,ctop), _ = perms_canonical_labels([cbot,ctop])

        bt = (tuple(cbot),tuple(ctop))

        if bt in cs_type_nb:
            cs_type_nb[bt] += 1
        else:
            cs_type_nb[bt] = 1

    shift = 0
    bot = []
    top = []
    keys = sorted(cs_type_nb.keys())
    for key in keys:
        for _ in range(cs_type_nb[key]):
            bot.extend(shift+i for i in key[0])
            top.extend(shift+i for i in key[1])
            shift += len(key[0])

    return bot,top


# main class
@total_ordering
class SeparatrixDiagram(SageObject):
    r"""
    Separatrix diagram of oriented foliation.

    A separatrix diagram is a 2-tuple of permutations ``(bot,top)`` such that
    ``bot`` and ``top`` share the same number of cycles.

    bot (resp. top) has to be thought a bottom (resp. top) of a potential face
    as in the following::

            -- bot -->
        -------------------
           <-- top --

    The order for bot and top is chosen in such a way that it corresponds to
    the orientation of a face.

    EXAMPLES::

        sage: from surface_dynamics import *

        sage: s = SeparatrixDiagram('(0,2)(1,3,4)','(0,4)(2,1,3)')
        sage: print(s)
        (0,2)(1,3,4)-(0,4)(1,3,2)
        sage: print(s.stratum())
        H_3(4)
    """
    def __init__(self, data, top=None, check=True, copy=True):
        r"""
        TESTS::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0,1)(2,3,4)','(0,2,4)(1,3)')
            sage: s == loads(dumps(s))
            True
            sage: SeparatrixDiagram(s)
            (0,1)(2,3,4)-(0,2,4)(1,3)
            sage: s == SeparatrixDiagram(s)
            True
            sage: SeparatrixDiagram(str(s))
            (0,1)(2,3,4)-(0,2,4)(1,3)
            sage: s == SeparatrixDiagram(str(s))
            True
        """
        if copy:
            if top is None:
                if isinstance(data,SeparatrixDiagram):
                    bot = data.bot()
                    top = data.top()
                elif isinstance(data,(list,tuple)) and len(data) == 2:
                    bot,top = data
                elif isinstance(data,str):
                    bot,top = data.split('-')
                else:
                    raise ValueError("the argument data is not valid")
            else:
                bot = data

            bot = perm_init(bot)
            top = perm_init(top)
            equalize_perms([bot, top])
        else:
            bot = data

        self._bot = bot
        self._top = top
        n = len(bot)

        bot_seen = [True] * n
        top_seen = [True] * n
        bot_to_cycle = [None]*n
        top_to_cycle = [None]*n
        bot_cycles = []
        top_cycles = []
        for i in range(n):
            if bot_seen[i]:
                c=[]
                k = len(bot_cycles)
                while bot_seen[i]:
                    bot_to_cycle[i] = k
                    bot_seen[i] = False
                    c.append(i)
                    i = self._bot[i]
                bot_cycles.append(tuple(c))
                k += 1
            if top_seen[i]:
                c=[]
                k = len(top_cycles)
                while top_seen[i]:
                    top_to_cycle[i] = k
                    top_seen[i] = False
                    c.append(i)
                    i = self._top[i]
                top_cycles.append(tuple(c))
                k += 1

        self._bot_cycles = bot_cycles
        self._top_cycles = top_cycles
        self._bot_to_cycle = bot_to_cycle
        self._top_to_cycle = top_to_cycle

        if check:
            self._check()

    def _check(self):
        r"""
        Check that the data of self is valid, i.e.

          * self._bot, self._top are valid permutations
          * the number of cycles of bot and top are the same

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: c = SeparatrixDiagram('(0,1)(2,3)','(0,2,3,1)') #indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: bot has 2 cylinders whereas top has 1
        """
        perm_check(self._bot)
        perm_check(self._top)

        p_bot = self.bot_cycle_tuples()
        p_top = self.top_cycle_tuples()
        if len(p_top) != len(p_bot):
            raise ValueError("bot has %d cylinders whereas top has %d"%(len(p_bot),len(p_top)))

    def _sage_input_(self, sib, coerced):
        r"""
        Sage input support.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0,3,2)(1,4)','(0,1)(2,3,4)')
            sage: sage_input(s)
            SeparatrixDiagram('(0,3,2)(1,4)-(0,1)(2,3,4)')

        We can check that evaluating the code actually gives back the same
        object::

            sage: t = eval(str(sage_input(s)))
            sage: t == s
            True
        """
        return sib.name('SeparatrixDiagram')(str(self))

    def to_directed_graph(self):
        r"""
        Return a graph that encodes this separatrix diagram.

        The vertices correspond to separatrix and the edges are of two types

        - 'b' neighbor corresponds to the right neighbors on the bottom
          permutation

        - 't' edges correspond to the neighbor of the top permutation

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: S = SeparatrixDiagram('(0,1)(2,3,4)','(0,3,2)(1,4)')
            sage: G = S.to_directed_graph(); G
            Looped multi-digraph on 5 vertices
            sage: G.vertices(sort=True)
            [0, 1, 2, 3, 4]
            sage: G.edges(sort=True)
            [(0, 1, 'b'), (0, 3, 't'), (1, 0, 'b'), (1, 4, 't'), (2, 0, 't'), (2, 3, 'b'), (3, 2, 't'), (3, 4, 'b'), (4, 1, 't'), (4, 2, 'b')]
        """
        G = DiGraph(multiedges=True, loops=True)
        G.add_edges((i, self._top[i], 't') for i in range(self.nseps()))
        G.add_edges((i, self._bot[i], 'b') for i in range(self.nseps()))
        return G

    def _repr_(self):
        r"""
        String representation of self

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: d = SeparatrixDiagram('(0,1)(2)','(0)(1,2)')
            sage: repr(d) #indirect doctest
            '(0,1)(2)-(0)(1,2)'
        """
        return self.bot_cycle_string() + "-" + self.top_cycle_string()

    #TODO
    def _latex_(self):
        r"""
        Latex representation

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: print("to be done")
            to be done
        """
        n = self._n
        if len(self.vertices_out().cycle_type()) == 1:
            v = self.vertices()[0]
            m = 360. / (2*n)
            d = dict([i,(v.index(i),v.index(-i))] for i in range(1,self._n+1))
            s = "\\begin{tikzpicture}\n"
            for i,(vout,vin) in iteritems(d):
                s += "    \\draw [-triangle 45] (0,0) -- (%f:0.8cm);\n" %(vout*m)
                s += "    \\draw (%f:0.8cm) -- (%f:1cm);\n" %(vout*m,vout*m)
                s += "    \\draw (%f:1cm) \n" %(vout*m)
                v1 = Integer(vout+vin)/2
                if vout+vin > 2*n:
                    v2 = v1 - n
                else:
                    v2 = v1 + n
                d1 = min([abs(vout-v1), abs(vout-v1-2*n), abs(vout-v1+2*n)])
                d2 = min([abs(vout-v2), abs(vout-v2-2*n), abs(vout-v2+2*n)])
                if d1 < d2:
                    vint = v1
                    d = d1
                else:
                    vint = v2
                    d = d2

                dint = '%fcm' %(1.5+d/2.)
                ct1 = '%fcm' %(d/2.)

                if cyclic_direction(vout,vint,vin) == 1:
                    ct2 = '%fcm' %(-d/2.)
                    ct3 = '%fcm' %(d/2.)
                else:
                    ct2 = '%fcm' %(d/2.)
                    ct3 = '%fcm' %(-d/2.)

                s += "    ..controls +(%f:%s) and +(%f:%s) ..\n" %(vout*m,ct1,(vint+n/2.)*m,ct2)
                s += "    (%f:%s)\n" %(vint*m,dint)
                s += "    ..controls +(%f:%s) and +(%f:%s) ..\n" %((vint+n/2.)*m,ct3,vin*m,ct1)
                s += "    (%f:1cm);\n" %(vin*m)
                s += "    \\draw [-open triangle 45] (%f:1cm) -- (%f:0.6cm);\n" %(vin*m,vin*m)
                s += "    \\draw (%f:0.6cm) -- (0,0);\n" %(vin*m)
            s += "\\end{tikzpicture}"
            return s
        else:
            return ""

    #
    # Comparisons and canonic labels
    #

    def __eq__(self, other):
        r"""
        Equality test

        TESTS::

            sage: from surface_dynamics import *

            sage: d1 = SeparatrixDiagram('(0)','(0)')
            sage: d2 = SeparatrixDiagram('(0,1)(2)','(0,1)(2)')
            sage: d3 = SeparatrixDiagram('(0,1)(2)','(0,2)(1)')
            sage: d1 == d1 and d2 == d2 and d3 == d3
            True
            sage: d1 == d2 or d1 == d3 or d2 == d3 or d3 == d2
            False

            sage: d1 = SeparatrixDiagram('(0)','(0)')
            sage: d2 = SeparatrixDiagram('(0,1)(2)','(0,1)(2)')
            sage: d3 = SeparatrixDiagram('(0,1)(2)','(0,2)(1)')
            sage: d1 != d1 or d2 != d2 or d3 != d3
            False
            sage: d1 != d2 and d1 != d3 and d2 != d3 and d3 != d2
            True
        """
        if type(self) != type(other):
            raise TypeError

        return self._bot == other._bot and self._top == other._top

    def __ne__(self, other):
        if type(self) is not type(other):
            raise TypeError
        return not (self == other)

    def __lt__(self, other):
        r"""
        Test whether self is lesser or equal than other

        TESTS::

            sage: from surface_dynamics import *

            sage: s = ['(0,1,2)-(0,1,2)',
            ....:      '(0,2,1)-(0,1,2)',
            ....:      '(0,2,1)-(0,2,1)',
            ....:      '(0)(1,2)-(0)(1,2)',
            ....:      '(0)(1,2)-(0,1)(2)',
            ....:      '(0)(1,2)-(0,2)(1)',
            ....:      '(0,1)(2)-(0)(1,2)',
            ....:      '(0,1)(2)-(0,1)(2)',
            ....:      '(0,1)(2)-(0,2)(1)',
            ....:      '(0,2)(1)-(0)(1,2)',
            ....:      '(0,2)(1)-(0,1)(2)',
            ....:      '(0,2)(1)-(0,2)(1)',
            ....:      '(0)(1)(2)-(0)(1)(2)']
            sage: s0 = list(map(SeparatrixDiagram, s))
            sage: s1 = s0[:]
            sage: for _ in range(10):
            ....:     shuffle(s1)
            ....:     assert sorted(s1) == s0
        """
        if type(self) != type(other):
            raise TypeError("only separatrix diagram can be compared to separatrix diagrams")

        if self.nseps() < other.nseps():
            return True
        elif self.nseps() > other.nseps():
            return False

        if self.ncyls() < other.ncyls():
            return True
        elif self.ncyls() > other.ncyls():
            return False

        if self._bot_cycles < other._bot_cycles:
            return True
        elif self._bot_cycles > other._bot_cycles:
            return False

        if self._top_cycles < other._top_cycles:
            return True
        elif self._top_cycles > other._top_cycles:
            return False

        # equal
        return False

    def is_isomorphic(self, other, return_map=False):
        r"""
        Test whether this separatrix diagram is isomorphic to ``other``.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: bot = [1,2,0,3]
            sage: top = [1,0,3,2]
            sage: s = SeparatrixDiagram(bot,top); s
            (0,1,2)(3)-(0,1)(2,3)
            sage: m = [3,0,1,2]
            sage: bot2 = [0]*4
            sage: top2 = [0]*4
            sage: for i in range(4):
            ....:     bot2[m[i]] = m[bot[i]]
            ....:     top2[m[i]] = m[top[i]]
            sage: ss = SeparatrixDiagram(bot2,top2)
            sage: s.is_isomorphic(ss)
            True
            sage: m = [1,2,0,3]
            sage: for i in range(4):
            ....:   bot2[m[i]] = m[bot[i]]
            ....:   top2[m[i]] = m[top[i]]
            sage: ss = SeparatrixDiagram(bot2,top2)
            sage: s.is_isomorphic(ss)
            True
        """
        if type(self) is not type(other):
            raise TypeError("other must be a separatrix diagram")
        if return_map:
            raise NotImplementedError
        return self._get_normal_perms() == other._get_normal_perms()

    def relabel(self, perm, inplace=False):
        r"""
        Relabel self according to the permutation ``perm``.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0)(2,3,4)','(0,3,2)(1)')
            sage: s
            (0)(1)(2,3,4)-(0,3,2)(1)(4)
            sage: s.relabel(perm=[1,0,2,3,4])
            (0)(1)(2,3,4)-(0)(1,3,2)(4)
            sage: s.relabel(perm=[1,2,0,3,4])
            (0,3,4)(1)(2)-(0,1,3)(2)(4)
        """
        n = self.degree()
        perm.extend(range(len(perm),n))
        bot = [None] * self.degree()
        top = [None] * self.degree()

        for i in range(n):
            bot[perm[i]] = perm[self._bot[i]]
            top[perm[i]] = perm[self._top[i]]

        S = SeparatrixDiagram(bot,top)
        if inplace:
            self._bot = S._bot
            self._top = S._top
            self._bot_cycles = S._bot_cycles
            self._top_cycles = S._top_cycles
            self._bot_to_cycle = S._bot_to_cycle
            self._top_to_cycle = S._top_to_cycle

        return S

    def canonical_label(self, inplace=False):
        r"""
        Relabel self according to some canonical labels.

        The result is cached.

        INPUT:

        - ``inplace`` - boolean (default: ``True``) - if True modify self if not
          return a new separatrix diagram.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: bot = '(0,1,3,6,7,5)(2,4)(8)(9)'
            sage: top = '(0)(1,2)(3,4,5)(6,7,8,9)'
            sage: s = SeparatrixDiagram(bot,top)
            sage: s.canonical_label()
            (0)(1)(2,3,4,5,6,7)(8,9)-(0,1,2,3)(4,7,9)(5)(6,8)

        TESTS::

            sage: from surface_dynamics import *

            sage: bot = [3,2,4,0,1]
            sage: top = [1,0,3,4,2]
            sage: b = [None]*5; t = [None]*5
            sage: for p in Permutations([0,1,2,3,4]):
            ....:     for i in range(5):
            ....:         b[p[i]] = p[bot[i]]
            ....:         t[p[i]] = p[top[i]]
            ....:     s = SeparatrixDiagram(b,t)
            ....:     print(s.canonical_label())
            (0,1)(2,3,4)-(0,2,4)(1,3)
            (0,1)(2,3,4)-(0,2,4)(1,3)
            (0,1)(2,3,4)-(0,2,4)(1,3)
            (0,1)(2,3,4)-(0,2,4)(1,3)
            (0,1)(2,3,4)-(0,2,4)(1,3)
            (0,1)(2,3,4)-(0,2,4)(1,3)
            ...
            (0,1)(2,3,4)-(0,2,4)(1,3)
            (0,1)(2,3,4)-(0,2,4)(1,3)
            (0,1)(2,3,4)-(0,2,4)(1,3)
            (0,1)(2,3,4)-(0,2,4)(1,3)
        """
        try:
            sep = self._normal_form
        except AttributeError:
            bot, top = self._get_normal_perms()
            self._normal_form = SeparatrixDiagram(bot, top, check=False, copy=False)
            self._normal_form._normal_form = self._normal_form
            sep = self._normal_form

        if inplace:
            other = self._normal_form
            self._bot = other._bot
            self._top = other._top
            self._bot_cycles = other._bot_cycles
            self._top_cycles = other._top_cycles
            self._bot_to_cycle = other._bot_to_cycle
            self._top_to_cycle = other._top_to_cycle
            return

        return sep

    def horizontal_symmetry(self):
        r"""
        Return the horizontal symmetric of this separatrix diagram.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0,1,2,3)(4,5)','(1,2,3)(4,5,0)')
            sage: sh = s.horizontal_symmetry()
            sage: print(sh)
            (0,5,4)(1,3,2)-(0,3,2,1)(4,5)

            sage: c = s.cylinder_diagrams(up_to_symmetry=False)[0].horizontal_symmetry()
            sage: ch = sh.cylinder_diagrams(up_to_symmetry=False)[0]
            sage: c.is_isomorphic(ch)
            True
        """
        return SeparatrixDiagram(tuple(t[::-1] for t in self._top_cycles),
                                 tuple(b[::-1] for b in self._bot_cycles))

    def vertical_symmetry(self):
        r"""
        Return the vertical symmetric of this separatrix diagram.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0,1,2,3)(4,5)','(1,2,3)(4,5,0)')
            sage: sv = s.vertical_symmetry()
            sage: print(sv)
            (0,3,2,1)(4,5)-(0,5,4)(1,3,2)

            sage: c = sv.cylinder_diagrams(up_to_symmetry=False)[0].vertical_symmetry()
            sage: cv = sv.cylinder_diagrams(up_to_symmetry=False)[0]
            sage: c.is_isomorphic(cv)
            True
        """
        return SeparatrixDiagram(tuple(b[::-1] for b in self._bot_cycles),
                                 tuple(t[::-1] for t in self._top_cycles))

    def inverse(self):
        r"""
        Return the inverse of this separatrix diagram, that is the one we obtain
        after application of `-Id`.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0,1,2)(3,4,5,6,7,8)-(0,1,3,5,7)(2,4,6,8)')
            sage: s.inverse()
            (0,1,3,5,7)(2,4,6,8)-(0,1,2)(3,4,5,6,7,8)
            sage: s.horizontal_symmetry().vertical_symmetry() == s.inverse()
            True
            sage: s.vertical_symmetry().horizontal_symmetry() == s.inverse()
            True
        """
        return SeparatrixDiagram(tuple(self._top_cycles), tuple(self._bot_cycles))

    def _get_normal_perms(self):
        r"""
        Returns the normal form of the permutations bot top defining self.

        Note that the result is cached.

        ALGORITHM:

        1) compute the orbit of G = <bot,top>

        2) for each of the connected component compute a normal form

        3) sort the list of (normal_top, normal_bot) and concatenate them

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0,3,2)(1,4)','(0,1)(2,3,4)')
            sage: s._get_normal_perms()
            ([1, 0, 3, 4, 2], [2, 3, 4, 1, 0])

            sage: s = SeparatrixDiagram('(0,5,2)(1,3,4)(6,7,8)','(0,3,7,8)(1,5)(2,4,6)')
            sage: s._get_normal_perms()
            ([1, 2, 0, 4, 5, 3, 7, 8, 6], [1, 3, 5, 6, 8, 7, 0, 2, 4])

            (0,5,2)(1,3,4)(6,7,8)-(0,3,7,8)(1,5)(2,4,6)
            sage: s.canonical_label() #indirect doctest
            (0,1,2)(3,4,5)(6,7,8)-(0,1,3,6)(2,5,7)(4,8)
        """
        try:
            return self._normal_bot, self._normal_top
        except AttributeError:
            pass

        self._normal_bot, self._normal_top = two_non_connected_perms_canonical_labels(self._bot, self._top)

        return self._normal_bot, self._normal_top

    def _get_sym_perms(self):
        r"""
        Return the four symmetric version of self.
        """
        n = len(self._top)
        bot, top = self._get_normal_perms()

        # compute the inverses
        ibot = [None]*n
        itop = [None]*n
        for i in range(n):
            ibot[bot[i]] = i
            itop[top[i]] = i

        hbot, htop = two_non_connected_perms_canonical_labels(itop, ibot)
        vbot, vtop = two_non_connected_perms_canonical_labels(ibot, itop)
        sbot, stop = two_non_connected_perms_canonical_labels(top, bot)

        return (bot,top),(hbot,htop),(vbot,vtop),(sbot,stop)

    def symmetries(self):
        r"""
        Return a triple of boolean ``(horiz_sym, vert_sym, inverse_sym)`` which
        correspond to the symmetry of ``self``.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0,1,2)(3,4,5)-(0,1)(2,3,4,5)')
            sage: s.symmetries()
            (False, True, False)
            sage: s.horizontal_symmetry().is_isomorphic(s)
            False
            sage: s.vertical_symmetry().is_isomorphic(s)
            True
            sage: s.inverse().is_isomorphic(s)
            False

            sage: s = SeparatrixDiagram('(0,1,3,5)(2,4)-(0,4,1,5)(2,3)')
            sage: s.symmetries()
            (True, False, False)
            sage: s.horizontal_symmetry().is_isomorphic(s)
            True
            sage: s.vertical_symmetry().is_isomorphic(s)
            False
            sage: s.inverse().is_isomorphic(s)
            False

            sage: s = SeparatrixDiagram('(0,1,3,5)(2,4)-(0,3,2,1)(5,4)')
            sage: s.symmetries()
            (False, False, True)
            sage: s.horizontal_symmetry().is_isomorphic(s)
            False
            sage: s.vertical_symmetry().is_isomorphic(s)
            False
            sage: s.inverse().is_isomorphic(s)
            True

            sage: s = SeparatrixDiagram('(0)(1,2,3,4,5)-(0,1,2,5,3)(4)')
            sage: s.symmetries()
            (False, False, False)
            sage: s.horizontal_symmetry().is_isomorphic(s)
            False
            sage: s.vertical_symmetry().is_isomorphic(s)
            False
            sage: s.inverse().is_isomorphic(s)
            False

        TESTS::

            sage: sym = lambda s: (s.horizontal_symmetry().is_isomorphic(s),
            ....:                  s.vertical_symmetry().is_isomorphic(s),
            ....:                  s.inverse().is_isomorphic(s))
            sage: from surface_dynamics.flat_surfaces.separatrix_diagram import separatrix_diagram_iterator
            sage: for s in separatrix_diagram_iterator((2,2,2,2)):
            ....:     assert s.symmetries() == sym(s)
            sage: for s in separatrix_diagram_iterator((4,2)):
            ....:     assert s.symmetries() == sym(s)
        """
        n = len(self._top)
        bot, top = self._get_normal_perms()

        # compute the inverses
        ibot = [None]*n
        itop = [None]*n
        for i in range(n):
            ibot[bot[i]] = i
            itop[top[i]] = i

        # horiz
        bot1, top1 = two_non_connected_perms_canonical_labels(itop, ibot)
        horiz_sym = bot == bot1 and top == top1

        # vert
        bot1, top1 = two_non_connected_perms_canonical_labels(ibot, itop)
        vert_sym = bot == bot1 and top == top1

        # inv
        if horiz_sym and vert_sym:  # got the two
            inverse_sym = True
        elif horiz_sym^vert_sym:    # got exactly one
            inverse_sym = False
        else:                       # none of them
            bot1, top1 = two_non_connected_perms_canonical_labels(top, bot)
            inverse_sym = bot == bot1 and top == top1

        return (horiz_sym, vert_sym, inverse_sym)

    def is_in_normal_form(self):
        r"""
        Test normal form

        Return True if self is in normal form and False otherwise.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0,1,2)(3,4,5)(6,7,8)','(0,3,7,8)(1,5)(2,4,6)')
            sage: s.is_in_normal_form()
            False
            sage: s.canonical_label().is_in_normal_form()
            True
        """
        return (self._bot, self._top) == self._get_normal_perms()

    #
    # Attributes access
    #

    def degree(self):
        r"""
        Return the degree (number of separatrices) of this separatrix diagram.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: S = SeparatrixDiagram('(0,1)(2,3)','(1,3,2)(0)')
            sage: S.degree()
            4
        """
        return len(self._top)

    nseps = degree
    num_edges = degree

    def ncyls(self):
        r"""
        Return the number of cylinders of this separatrix diagram.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: S = SeparatrixDiagram('(0,1)(2,3)','(1,3,2)(0)')
            sage: S.ncyls()
            2
        """
        return len(self._top_cycles)

    def profile(self):
        r"""
        Return the angles around each vertex

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: a = AbelianStratum(1,1,0)
            sage: s = a.separatrix_diagrams()[0]
            sage: s.profile()
            [2, 2, 1]
        """
        p = perm_cycle_type(self.outgoing_edges_perm())

        from sage.combinat.partition import Partition
        return Partition(p)

    def euler_characteristic(self):
        r"""
        Return the Euler characteristic

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: SeparatrixDiagram('(0)','(0)').euler_characteristic()
            0

            sage: CylinderDiagram([((0,),(0,))]).euler_characteristic()
            0
            sage: CylinderDiagram([((0,1),(0,2)), ((2,),(1,))]).euler_characteristic()
            -2
        """
        p = self.profile()
        return Integer(len(p)-sum(p))

    def genus(self):
        r"""
        Return the genus

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: CylinderDiagram([((0,),(0,))]).genus()
            1
            sage: CylinderDiagram([((0,1),(0,1))]).genus()
            1
            sage: CylinderDiagram([((0,1,2),(0,1,2))]).genus()
            2
            sage: CylinderDiagram([((0,1,2,3),(0,1,2,3))]).genus()
            2
            sage: CylinderDiagram([((0,1,2,3,4),(0,1,2,3,4))]).genus()
            3
        """
        return Integer(1 - self.euler_characteristic()//2)

    def stratum(self):
        r"""
        Return the Abelian stratum this separatrix diagram belongs to.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: SeparatrixDiagram('(0)(1)(2)','(0)(1)(2)').stratum()
            H_1(0^3)
            sage: SeparatrixDiagram('(0,1)(2)','(0,2)(1)').stratum()
            H_2(2)
        """
        from .abelian_strata import AbelianStratum

        return AbelianStratum([i-1 for i in self.profile()])

    def bot(self):
        r"""
        The bot permutation as a list from 0 to nseps-1

        Warning: the output list should not be modified

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0)(1,2)','(0,1)(2)')
            sage: s.bot()
            [0, 2, 1]
        """
        return list(self._bot)

    def bot_perm(self):
        r"""
        Return the bot as a permutation (element of a group)

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0)(1,2)','(0,1)(2)')
            sage: s.bot_perm()
            (2,3)
        """
        try:
            # Trac #28652: Rework the constructor of PermutationGroupElement
            from sage.groups.perm_gps.constructor import PermutationGroupElement
        except ImportError:
            from sage.groups.perm_gps.permgroup_element import PermutationGroupElement

        return PermutationGroupElement([i+1 for i in self._bot])

    def bot_orbit(self, i):
        r"""
        Return the orbit of i under the bot permutation

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0,1)(2,5)(3,4,6)','(0,1,5)(2,3,6)(4)')
            sage: s.bot_orbit(0)
            (0, 1)
            sage: s.bot_orbit(4)
            (3, 4, 6)
        """
        return self._bot_cycles[self._bot_to_cycle[i]]

    def bot_cycle_tuples(self):
        r"""
        Return the cycles of the bottom permutation as a list of tuples.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: S = SeparatrixDiagram('(0,2)(3,4)','(0)(1,2,3)')
            sage: S.bot_cycle_tuples()
            [(0, 2), (1,), (3, 4)]
        """
        return self._bot_cycles

    def bot_cycle_string(self):
        r"""
        Return the cycles of the top permutation as a string.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: S = SeparatrixDiagram('(0,2)(3,4)','(0)(1,2,3)')
            sage: S.bot_cycle_string()
            '(0,2)(1)(3,4)'
        """
        return ''.join('(' + ','.join(map(str,c)) +')' for c in self.bot_cycle_tuples())

    def top(self):
        r"""
        Return the top permutation of self as a list.

        Warning: the output should not be modified

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0,1,3)(2,4)','(0,4)(1,2,3)')
            sage: s.top()
            [4, 2, 3, 1, 0]
        """
        return self._top

    def top_perm(self):
        r"""
        Return the top as a permutation

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0)(1,2)','(1)(0,2)')
            sage: s.top_perm()
            (1,3)
        """
        try:
            # Trac #28652: Rework the constructor of PermutationGroupElement
            from sage.groups.perm_gps.constructor import PermutationGroupElement
        except ImportError:
            from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
        return PermutationGroupElement([i+1 for i in self._top])

    def top_orbit(self,i):
        r"""
        Return the orbit of ``i`` under the top permutation.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0,1)(2,5)(3,4,6)','(0,1,5)(2,3,6)(4)')
            sage: s.top_orbit(0)
            (0, 1, 5)
            sage: s.top_orbit(6)
            (2, 3, 6)
        """
        return self._top_cycles[self._top_to_cycle[i]]

    def top_cycle_tuples(self):
        r"""
        Return the cycle of the top permutation as a list of tuples.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: S = SeparatrixDiagram('(0,2)(3,4)','(0)(1,2,3)')
            sage: S.top_cycle_tuples()
            [(0,), (1, 2, 3), (4,)]
        """
        return self._top_cycles

    def top_cycle_string(self):
        r"""
        Return the cycle of the top permutation as a string.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: S = SeparatrixDiagram('(0,2)(3,4)','(0)(1,2,3)')
            sage: S.top_cycle_string()
            '(0)(1,2,3)(4)'
        """
        return ''.join('(' + ','.join(map(str,c)) + ')' for c in self.top_cycle_tuples())

    def automorphism_group(self, implementation='graph'):
        r"""
        Return the automorphism group of self.

        That is the centralizer of the permutations top and bottom.

        INPUT:

        - ``implementation`` - either graph or gap

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: S = SeparatrixDiagram('(0,3,1,4,2)','(0,1,2,3,4)')
            sage: G1 = S.automorphism_group(implementation='graph'); G1
            Permutation Group with generators [(0,1,2,3,4)]
            sage: G2 = S.automorphism_group(implementation='gap'); G2
            Subgroup ...
            sage: G1.is_isomorphic(G2)
            True
        """
        if implementation == 'graph':
            return self.to_directed_graph().automorphism_group(edge_labels=True)

        elif implementation == 'gap':
            from sage.groups.perm_gps.permgroup import PermutationGroup
            from sage.groups.perm_gps.permgroup_named import SymmetricGroup

            return SymmetricGroup(self.nseps()).centralizer(PermutationGroup([self.top_perm(),self.bot_perm()]))

        else:
            raise ValueError("implementation should be either 'graph' or 'gap'")

    def homological_dimension_of_cylinders(self):
        r"""
        Returns the dimension in the first homology group of the span of waist
        curves of horizontal cylinders.

        EXAMPLES::

            sage: from surface_dynamics import *

        Homological dimension in the stratum H(2)::

            sage: c = CylinderDiagram('(0,1,2)-(0,1,2)')
            sage: c.stratum()
            H_2(2)
            sage: c.homological_dimension_of_cylinders()
            1
            sage: c = CylinderDiagram('(0,1)-(1,2) (2)-(0)')
            sage: c.stratum()
            H_2(2)
            sage: c.homological_dimension_of_cylinders()
            2

        Homological dimensions for cylinder diagrams in H(1,1)::

            sage: c = CylinderDiagram('(0,1,2,3)-(0,1,2,3)')
            sage: c.stratum()
            H_2(1^2)
            sage: c.homological_dimension_of_cylinders()
            1
            sage: c = CylinderDiagram('(0,1)-(0,2) (2,3)-(1,3)')
            sage: c.stratum()
            H_2(1^2)
            sage: c.homological_dimension_of_cylinders()
            2
            sage: c = CylinderDiagram('(0,1,2)-(1,2,3) (3)-(0)')
            sage: c.stratum()
            H_2(1^2)
            sage: c.homological_dimension_of_cylinders()
            2
            sage: c = CylinderDiagram('(0,1)-(2,3) (2)-(0) (3)-(1)')
            sage: c.stratum()
            H_2(1^2)
            sage: c.homological_dimension_of_cylinders()
            2
        """
        return Integer(self.ncyls() - SeparatrixDiagram.to_directed_graph(self).connected_components_number() + 1)

    def homologous_cylinders(self):
        r"""
        Return the list of homologous cylinders.

        OUTPUT: a list of lists. Each sublist is an equivalence class of > 1
        homologous cylinders.

        EXAMPLES::

            sage: from surface_dynamics import CylinderDiagram
            sage: c = CylinderDiagram('(0,7,1,2)-(3,6,4,5) (3,6,4,5)-(0,7,1,2)')
            sage: c.homologous_cylinders()
            [[0, 1]]

            sage: c = CylinderDiagram('(0,1)-(2,3) (2)-(0) (3)-(1)')
            sage: c.homologous_cylinders()
            []

            sage: c = CylinderDiagram('(0,2,1)-(9) (3,6,4,5)-(7,10,8) (7,9,8)-(3,6,4,5) (10)-(0,2,1)')
            sage: c.homologous_cylinders()
            [[0, 3], [1, 2]]
        """
        from .twist_space import TwistSpace
        return TwistSpace(self).homologous_cylinders()

    #
    # Vertices of the separatrix diagram
    #

    def outgoing_edges_perm(self):
        r"""
        Permutation associated to turning around vertices in trigonometric
        order.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0,1)','(2,3)')
            sage: s.outgoing_edges_perm()
            [1, 0, 3, 2]

            sage: s = SeparatrixDiagram('(0,5,2)(1,3,4)(6,7,8)','(0,3,7,8)(1,5)(2,4,6)')
            sage: s.outgoing_edges_perm()
            [6, 2, 1, 5, 0, 8, 7, 4, 3]

            sage: c = CylinderDiagram('(0,5,2,1,4,3)-(0,4,2,1,5,3)')
            sage: c.outgoing_edges_perm()
            [5, 4, 1, 0, 2, 3]
        """
        return perm_compose_i(self._bot,self._top)

    def incoming_edges_perm(self):
        r"""
        Permutation associated to turning around vertices in trigonometric
        order.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0,1)','(2,3)')
            sage: s.incoming_edges_perm()
            [1, 0, 3, 2]

            sage: s = SeparatrixDiagram('(0,5,2)(1,3,4)(6,7,8)','(0,3,7,8)(1,5)(2,4,6)')
            sage: s.incoming_edges_perm()
            [7, 0, 8, 2, 5, 4, 3, 1, 6]

            sage: c = CylinderDiagram('(0,5,2,1,4,3)-(0,4,2,1,5,3)')
            sage: c.incoming_edges_perm()
            [4, 5, 1, 0, 3, 2]
        """
        return perm_compose_i(self._top,self._bot)

    def vertices(self):
        r"""
        Return a list of pairs ``(vout, vin)`` where each pair represent a vertex and
        ``vout``, ``vin`` are respectively the labels of outgoing/incoming separatrices
        at this vertex.

        EXAMPLES::


            sage: from surface_dynamics import CylinderDiagram
            sage: from surface_dynamics.misc.permutation import perm_cycles

            sage: c = CylinderDiagram('(0,1)-(0,2,6,7,4,5,3,8) (2,7,3,4,8,6,5)-(1)')
            sage: c.vertices()
            [([0, 1, 8, 7], [0, 4, 2, 1]), ([2, 4, 5], [3, 6, 5]), ([3, 6], [7, 8])]
            sage: perm_cycles(c.outgoing_edges_perm())
            [[0, 1, 8, 7], [2, 4, 5], [3, 6]]
            sage: perm_cycles(c.incoming_edges_perm())
            [[0, 4, 2, 1], [3, 6, 5], [7, 8]]

            sage: c = CylinderDiagram('(0,1,3,4,2)-(0,1,5,7,6) (5,6)-(4) (7)-(2,3)')
            sage: perm_cycles(c.outgoing_edges_perm())
            [[0, 3], [1, 6], [2, 4], [5, 7]]
            sage: perm_cycles(c.incoming_edges_perm())
            [[0, 5], [1, 2], [3, 4], [6, 7]]
        """
        top = perm_invert(self._top)
        bot = perm_invert(self._bot)

        seen = [False] * self.nseps()
        res = []
        for i in range(self.nseps()):
            if seen[i]:
                continue
            cout = []
            cin = []
            while seen[i] is False:
                seen[i] = True
                cout.append(i)
                cin.append(bot[i])
                i = top[bot[i]]
            j = argmin(cin)
            res.append((cout, cin[j:] + cin[:j]))
        return res

    #
    # to cylinder diagram
    #

    def to_cylinder_diagram(self, pairing):
        r"""
        Return a cylinder diagram with the given pairing

        The pairing should be a list of 2-tuples of integer.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0,1,3)(2,4)','(0,2)(1,4,3)'); s
            (0,1,3)(2,4)-(0,2)(1,4,3)

            sage: s.to_cylinder_diagram([(0,0),(1,1)])
            (0,1,3)-(0,2) (2,4)-(1,4,3)
            sage: s.to_cylinder_diagram([(1,1),(0,0)])
            (0,1,3)-(0,2) (2,4)-(1,4,3)

            sage: s.to_cylinder_diagram([(0,1),(1,0)])
            (0,1,3)-(1,4,3) (2,4)-(0,2)
            sage: s.to_cylinder_diagram([(1,0),(0,1)])
            (0,1,3)-(1,4,3) (2,4)-(0,2)
        """
        from copy import copy

        other = copy(self)
        other.__class__ = CylinderDiagram

        bots = self.bot_cycle_tuples()
        tops = self.top_cycle_tuples()

        other._bot_to_cyl = [None]*self.nseps()
        other._top_to_cyl = [None]*self.nseps()

        for i in range(len(pairing)):
            b = bots[pairing[i][0]]
            t = tops[pairing[i][1]]
            cyl = (b[0],t[0])

            for j in b:
                other._bot_to_cyl[j] = cyl
            for j in t:
                other._top_to_cyl[j] = cyl

        return other

    @rename_keyword(deprecation=666, up_to_isomorphism='up_to_symmetry')
    def cylinder_diagram_iterator(self, connected=True, up_to_symmetry=True):
        r"""
        Construct all cylinder diagrams from given separatrix diagram (i.e. a pair
        of permutations).

        INPUT:

        - ``connected`` - boolean (default: True) - if true, returns only
          connected cylinder diagrams.

        - ``up_to_symmetry`` - boolean (default: True) take care of the horizontal
          and vertical symmetries.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0,1)(2,3)(4,5)','(1,2)(3,4)(5,0)')
            sage: for c in s.cylinder_diagram_iterator(): print(c)  # random
            (0,5)-(0,4) (1,4)-(1,3) (2,3)-(2,5)
            (0,3)-(0,5) (1,2)-(1,4) (4,5)-(2,3)
            (0,5)-(3,4) (1,4)-(0,2) (2,3)-(1,5)
            sage: sum(1 for _ in s.cylinder_diagram_iterator(up_to_symmetry=False))
            6
            sage: sum(1 for _ in s.cylinder_diagram_iterator(up_to_symmetry=True))
            3

        Here is an example with some symmetry::

            sage: s = SeparatrixDiagram('(0)(1)(2,3)(4,5,6)-(0,1)(2,4)(3,5)(6)')
            sage: s.vertical_symmetry().canonical_label() == s
            True
            sage: C1 = [CylinderDiagram('(0,1)-(4) (2,4,3)-(5,6) (5)-(0,2) (6)-(1,3)'),
            ....:       CylinderDiagram('(0,3,1)-(0,6) (2,6)-(4,5) (4)-(1) (5)-(2,3)'),
            ....:       CylinderDiagram('(0,1)-(0,4) (2,3,4)-(5,6) (5)-(2) (6)-(1,3)')]
            sage: C2 = s.cylinder_diagrams()
            sage: assert len(C1) == len(C2)
            sage: for (c1, c2) in zip(C1, C2): assert c1.is_isomorphic(c2)

        TESTS::

            sage: from surface_dynamics import SeparatrixDiagram
            sage: s = SeparatrixDiagram('(0,3)(1,4,5)(2)','(0)(1,2)(3,4,5)')
            sage: C2 = s.cylinder_diagram_iterator(up_to_isomorphism=True)
            doctest:warning
            ...
            DeprecationWarning: use the option 'up_to_symmetry' instead of 'up_to_isomorphism'
            ...
        """
        cbot = self.bot_cycle_tuples()
        ctop0 = self.top_cycle_tuples()

        connected = not connected

        if up_to_symmetry:
            # NOTE: here we should only consider symmetries of the cylinder diagrams
            # only when this underlying separatrix diagrams has some. But the
            # canonical labels of cylinder diagrams and separatrix diagrams are
            # not compatible!!

            s = set([])
            hsym, vsym, isym = self.symmetries()
            for ctop in itertools.permutations(ctop0):
                c = CylinderDiagram(zip(cbot,ctop), check=False)
                c.canonical_label(inplace=True)
                if c in s:
                    continue

                cc = [c]
                if hsym:
                    c1 = c.horizontal_symmetry()
                    c1.canonical_label(inplace=True)
                    cc.append(c1)
                if vsym:
                    c1 = c.vertical_symmetry()
                    c1.canonical_label(inplace=True)
                    cc.append(c1)
                if isym:
                    c1 = c.inverse()
                    c1.canonical_label(inplace=True)
                    cc.append(c1)
                s.update(cc)

                if (connected or c.is_connected()) and c.smallest_integer_lengths():
                    yield min(cc)

        else:
            for ctop in itertools.permutations(ctop0):
                c = CylinderDiagram(zip(cbot,ctop),check=False)
                if (connected or c.is_connected()) and c.smallest_integer_lengths():
                    yield c

    @rename_keyword(deprecation=666, up_to_isomorphism='up_to_symmetry')
    def cylinder_diagrams(self, connected=True, up_to_symmetry=True):
        r"""
        Return the list of cylinder diagrams associated to this separatrix
        diagram.

        We warn that the cylinder diagram may be renumeroted in the output list
        (in order to prevent repetitions). If you care about numerotation the
        option ``up_to_symmetry`` should be set to False.

        INPUT:

        - ``connected`` - boolean (default: True)

        - ``up_to_symmetry`` - boolean (default: True)

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0)(1)(2)','(0)(1)(2)')
            sage: for c in s.cylinder_diagrams(connected=True): print(c)
            (0)-(2) (1)-(0) (2)-(1)
            sage: for c in s.cylinder_diagrams(connected=False): print(c)
            (0)-(0) (1)-(1) (2)-(2)
            (0)-(1) (1)-(0) (2)-(2)
            (0)-(2) (1)-(0) (2)-(1)

            sage: s = SeparatrixDiagram('(0,1)(2)','(0)(1,2)')
            sage: C1 = [CylinderDiagram('(0,1)-(0,2) (2)-(1)')]
            sage: C2 = s.cylinder_diagrams()
            sage: assert len(C1) == len(C2)
            sage: for (c1, c2) in zip(C1, C2): assert c1.is_isomorphic(c2)

        In the example below, there is no isomorphism problem for the cylinder
        diagram generation as the separatrix diagram admit no automorphism::

            sage: s = SeparatrixDiagram('(0,3)(1,4,5)(2)','(0)(1,2)(3,4,5)')
            sage: s.automorphism_group()
            Permutation Group with generators [()]
            sage: C1 = [CylinderDiagram('(0,3,1)-(5) (2,5)-(3,4) (4)-(0,2,1)'),
            ....:       CylinderDiagram('(0,1,2)-(0,1,5) (3,5)-(2,4) (4)-(3)'),
            ....:       CylinderDiagram('(0,2,3)-(2,5) (1,4)-(0,1,3) (5)-(4)')]
            sage: C2 = s.cylinder_diagrams()
            sage: C3 = s.cylinder_diagrams(up_to_symmetry=False)
            sage: assert len(C1) == len(C2) == len(C3)
            sage: for (c1, c2, c3) in zip(C1, C2, C3): assert c1.is_isomorphic(c2) and c1.is_isomorphic(c3)

        TESTS::

            sage: from surface_dynamics import SeparatrixDiagram
            sage: s = SeparatrixDiagram('(0,3)(1,4,5)(2)','(0)(1,2)(3,4,5)')
            sage: C2 = s.cylinder_diagrams(up_to_isomorphism=True)
            doctest:warning
            ...
            DeprecationWarning: use the option 'up_to_symmetry' instead of 'up_to_isomorphism'
            ...
        """
        return list(self.cylinder_diagram_iterator(
            connected=connected,
            up_to_symmetry=up_to_symmetry))

    def saddle_connections_graph(self, mutable=False):
        r"""
        Return the fat graph (or ribbon graph) made by the saddle connections.

        The return graph is a
        :class:`~surface_dynamics.topology.fat_graph.FatGraph`. The saddle
        connection labelled `i` on this diagram gets labels `2i` and `2i+1` in
        the graph (there one label per half-edge in the fat graph). The even
        labels correspond to half-edges in the bottom of cylinders while the
        odd ones correspond to the top.

        EXAMPLES::

            sage: from surface_dynamics import AbelianStratum
            sage: H11 = AbelianStratum(1,1).unique_component()
            sage: for cd in H11.cylinder_diagrams():
            ....:     fg = cd.saddle_connections_graph()
            ....:     print(cd.ncyls(), [comp.genus() for comp in fg.connected_components()])
            1 [1]
            2 [0]
            2 [0]
            3 [0, 0]

        TESTS::

            sage: from surface_dynamics import AbelianStrata
            sage: for g in (2, 3):
            ....:     for H in AbelianStrata(genus=g):
            ....:         for C in H.components():
            ....:             for cd in C.cylinder_diagrams():
            ....:                 fg = cd.saddle_connections_graph()
            ....:                 assert fg.num_faces() == 2 * cd.ncyls()
            ....:                 assert fg.euler_characteristic() == cd.euler_characteristic() + 2 * cd.ncyls()
            ....:                 assert len(fg.connected_components()) == 1 + cd.ncyls() - cd.homological_dimension_of_cylinders()
        """
        n = len(self._bot)
        fp = [-1] * (2 * n)
        for i in range(n):
            fp[2 * i] = 2 * self._bot[i]
            fp[2 * i + 1] = 2 * self._top[i] + 1

        from surface_dynamics.topology.fat_graph import FatGraph
        return FatGraph(fp=fp, mutable=mutable)


def cyclic_direction(x,y,z):
    r"""
    Returns 1 or -1 depending on the cyclic ordering of (x,y,z)

    TESTS::

        sage: from surface_dynamics.flat_surfaces.separatrix_diagram import cyclic_direction
        sage: cyclic_direction(0,1,2)
        1
        sage: cyclic_direction(1,2,0)
        1
        sage: cyclic_direction(2,0,1)
        1
        sage: cyclic_direction(2,1,0)
        -1
        sage: cyclic_direction(1,0,2)
        -1
        sage: cyclic_direction(0,2,1)
        -1
    """
    if (x < y < z) or (y < z < x) or (z < x < y): return 1
    else: return -1


# iterators

def separatrix_diagram_fast_iterator(profile,ncyls=None):
    r"""
    Iterator over separatrix diagram with given ``profile``

    Return a list of 3-tuples ``[bot, top, s]`` where ``bot`` and ``top`` are
    list on 0, ..., nseps-1 that corresponds to a separatrix diagram with
    profile ``profile`` while ``s`` is the element conjugacy class corresponding
    to the profile which equals ``bot * top``.

    If ncyls is not None, it should be a list of integers from which the number
    of cylinders is considered.

    Warning: each isomorphism class of separatrix diagram is output more than
    once in general. If you want a unique representative in each isomorphism
    class you may consider the method separatrix_diagram_iterator instead.

    Uses bounds from [Nav08]_.

    EXAMPLES::

        sage: from surface_dynamics import *

        sage: from surface_dynamics.flat_surfaces.separatrix_diagram import separatrix_diagram_fast_iterator
        sage: for s in sorted(separatrix_diagram_fast_iterator([3])): print(s)
        ([0, 2, 1], [1, 0, 2], [(0, 1, 2)])
        ([1, 2, 0], [1, 2, 0], [(0, 2, 1)])
        ([2, 1, 0], [1, 0, 2], [(0, 2, 1)])
        sage: for s in sorted(separatrix_diagram_fast_iterator([2,2])): print(s)
        ([0, 1, 3, 2], [1, 0, 2, 3], [(0, 1), (2, 3)])
        ([0, 2, 3, 1], [1, 2, 0, 3], [(0, 1), (2, 3)])
        ([1, 2, 3, 0], [1, 2, 3, 0], [(0, 2), (1, 3)])
        ([1, 3, 2, 0], [1, 2, 0, 3], [(0, 2), (1, 3)])
        ([2, 3, 0, 1], [1, 0, 3, 2], [(0, 3), (1, 2)])
        ([3, 1, 0, 2], [1, 2, 0, 3], [(0, 3), (1, 2)])
        ([3, 2, 1, 0], [1, 0, 3, 2], [(0, 2), (1, 3)])
    """
    from sage.combinat.partition import Partition,Partitions
    from sage.groups.perm_gps.symgp_conjugacy_class import conjugacy_class_iterator

    part = Partition(profile)
    n = sum(part)
    d = (n+len(part))//2  # the maximum number of cylinders is known
                          # to be g+s-1 from a theorem of Y. Naveh

    tops = [[]]
    if ncyls is None:
        ncyls = list(range(1,d+1))
    else:
        if isinstance(ncyls, numbers.Integral):
            ncyls = set([int(ncyls)])
        else:
            ncyls = set(map(int,ncyls))
        for i in ncyls:
            if i < 1 or i > d:
                raise ValueError("%d is not possible as number of cylinders"%i)

    # build the list of admissible tops up to conjugacy class
    for k in range(1,d+1):
        tops.append([])
        if k in ncyls:
            for p in Partitions(n,length=k):
                tops[-1].append((canonical_perm(p),canonical_perm_i(p)))

    for s in conjugacy_class_iterator(part, list(range(n))):
        for k in range(len(tops)):
            for top,top_i in tops[k]:
                bot = list(range(len(top_i)))
                for cycle in s:
                    for i in range(len(cycle)-1):
                        bot[cycle[i]] = top_i[cycle[i+1]]
                    bot[cycle[-1]] = top_i[cycle[0]]

                seen = [True]*len(bot)
                nb_cycles = 0
                for i in range(len(bot)):
                    if seen[i]:
                        seen[i] = False
                        nb_cycles += 1
                        if nb_cycles > k:
                            break
                        j = bot[i]
                        while seen[j]:
                            seen[j] = False
                            j = bot[j]

                if nb_cycles == k:
                    yield (bot,top,s)

def separatrix_diagram_iterator(profile, ncyls=None):
    r"""
    Iterator over separatrix diagram with given ``profile`` and number of
    cylinders.

    Warning: to prevent isomorphism class to be output twice the function
    implement a cache mechanism. If you intend to iterate through a huge
    class of separatrix_diagram and do not care about isomorphism problem use
    separatrix_diagram_fast_iterator instead.

    EXAMPLES::

        sage: from surface_dynamics import *

        sage: from surface_dynamics.flat_surfaces.separatrix_diagram import separatrix_diagram_iterator

        sage: for s in sorted(separatrix_diagram_iterator([1,1])): print(s)
        (0,1)-(0,1)
        (0)(1)-(0)(1)

        sage: for s in sorted(separatrix_diagram_iterator([3])): print(s)
        (0,1,2)-(0,1,2)
        (0)(1,2)-(0,1)(2)

        sage: for s in sorted(separatrix_diagram_iterator([2,2])): print(s)
        (0,1,2,3)-(0,1,2,3)
        (0)(1,2,3)-(0,1,2)(3)
        (0,1)(2,3)-(0,2)(1,3)
        (0)(1)(2,3)-(0,1)(2)(3)

        sage: sum(1 for s in separatrix_diagram_iterator([3,2,2]))
        64
    """
    res = set([])
    for bot,top,_ in separatrix_diagram_fast_iterator(profile,ncyls):
        bot, top = two_non_connected_perms_canonical_labels(bot, top)
        s_perm = tuple(bot+top)
        if s_perm not in res:
            s = SeparatrixDiagram(bot, top, check=False, copy=False)
            syms = s._get_sym_perms()
            s = SeparatrixDiagram(*min(syms), check=False, copy=False)
            res.update(tuple(bot+top) for bot,top in syms)
            yield s

#
# Cylinder diagram
#  (or completely periodic decomposition)
#

def string_to_cycle(s):
    r"""
    TESTS::

        sage: from surface_dynamics.flat_surfaces.separatrix_diagram import string_to_cycle
        sage: string_to_cycle('(3,1,2)')
        (3, 1, 2)
    """
    if len(s) < 2:
        raise ValueError("Wrong syntax")
    if s[0] != '(':
        raise ValueError("A cycle string should start with an opening parenthesis")
    if s[-1] != ')':
        raise ValueError("A cycle string should end with a closing parenthesis")
    return tuple(int(i) for i in s[1:-1].split(','))


def orientation_cover(alpha,phi,a,verbose=0):
    r"""
    Build the cylinder diagram of Abelian differentials that double covers it.

    A quadratic differrential separatrix diagram is given by three permutations

    - sigma: the permutation of 1/2-separatrices around vertices
    - alpha: the permutation of 1/2-separatrices that describe the separatrices
       (it is a fixed point free involution)
    - phi: the permutation of 1/2-separatrices that describe the cycles.

    INPUT:

    - ``alpha`` -- permutation

    - ``phi`` -- permutation

    - ``a`` -- number of half separatrices

    EXAMPLES::

        sage: from surface_dynamics import *

        sage: from surface_dynamics.flat_surfaces.separatrix_diagram import orientation_cover
        sage: alpha = [3, 2, 1, 0, 5, 4, 7, 6]
        sage: phi = [3, 1, 0, 2, 5, 4, 7, 6]
        sage: orientation_cover(alpha,phi,3)
        (0,2)-(0,1) (1)-(2)
    """
    if verbose: print(" orientation cover")
    cyls = []
    todo = [True]*a

    for i in range(a):
        if todo[i]:
            todo[i] = False
            b = [i]
            if alpha[i] >= a:
                t = [i]
            else:
                t = [alpha[i]]
            if verbose: print("  top from %d,  bot from %d" % (i, t[0]))
            j = phi[i]
            if j >= a:
                j = phi[j]
            while j != i:
                todo[j] = False
                b.append(j)

                if alpha[j] >= a:
                    t.append(j)
                else:
                    t.append(alpha[j])
                if verbose: print("  add %d to bot,  add %d to top" % (j, b[-1]))

                j = phi[j]
                if j >= a:
                    j = phi[j]

            cyls.append((b,t))

    return CylinderDiagram(cyls)

#TODO: do something less stupid for symmetries
def hyperelliptic_cylinder_diagram_iterator(a,verbose=False):
    r"""
    Return an iterator over cylinder diagrams of Abelian differentials that
    double covers Q((a-2), -1^(a+2)).

    The generator is up to isomorphism.

    TODO:

    - An optimization could be obtained by considering the generation of
      k-subsets of {1,...,n} up to the cyclic symmetry of the tree.

    INPUT:

    - ``a`` - integer - angle of the conical singularity of the quadratic
      differential.

    - ``verbose`` - integer (default: 0) - output various information during the
      iteration (mainly for debug).

    EXAMPLES::

        sage: from surface_dynamics import *

        sage: from surface_dynamics.flat_surfaces.separatrix_diagram import hyperelliptic_cylinder_diagram_iterator
        sage: it = hyperelliptic_cylinder_diagram_iterator(3)
        sage: c = next(it); c.is_isomorphic(CylinderDiagram('(0,1)-(0,2) (2)-(1)'))
        True
        sage: c.stratum_component()
        H_2(2)^hyp

        sage: hyp = AbelianStratum(2,2).hyperelliptic_component()
        sage: all(c.stratum_component() == hyp for c in hyperelliptic_cylinder_diagram_iterator(6))
        True
    """
    from surface_dynamics.misc.plane_tree import admissible_plane_tree_iterator
    from sage.combinat.gray_codes import combinations

    cyl_diags = set([])
    if verbose is True: verbose=1

    B = [False]*(2*a+2)  # open loops indicator
                         # if B[k] is not False, it is where loop k starts
    sigma = list(range(1,a)) + [0] + list(range(a,2*a+2))
    for t,n,l in admissible_plane_tree_iterator(a):
        # Build the initial tree
        L = []
        p = 2*n-a               # the number of poles
        ll = 0                  # leaf counter
        s = 0                   # 1/2-separatrix counter
        sp = a                  # pole counter
        alpha = [None]*(2*a+2)  # edge permutation
        phi = [None]*(2*a+2)    # face permutation
        if verbose:
            print("n = %d,  l = %d,  p = %d" % (n, l, p))
            print("t =", t)
        for k in range(1,n+2):
            if verbose: print(" k = %d" % k)
            for kk in range(t[k-1],t[k]-1,-1): # close current loops
                if B[kk] is not False:
                    if verbose:
                        print(" close loop from %d to %d" % (B[kk], s))
                    alpha[B[kk]] = s
                    alpha[s] = B[kk]
                    phi[s] = (B[kk]-1)%a
                    phi[B[kk]] = (s-1)%a
                    s += 1
                    if verbose > 2:
                        print(" alpha =", alpha)
                        print(" phi   =", phi)
            if ll < p and t[k] >= t[k+1]:
                L.append(s) # store the leaf
                # t[k] is a pole
                if verbose: print(" pole at %d" % s)
                alpha[s] = sp
                alpha[sp] = s
                phi[s] = sp
                phi[sp] = (s-1)%a
                s += 1
                sp += 1
                ll += 1
                B[t[k]] = False
                if verbose > 2:
                    print(" alpha =", alpha)
                    print(" phi   =", phi)

            elif k != n+1: # not at the end -> open a loop
                if t[k] >= t[k+1]: # store the leaf
                    L.append(s)
                if verbose: print(" open loop at %d" % s)
                B[t[k]] = s
                s += 1
                if verbose > 2:
                    print(" alpha =", alpha)
                    print(" phi   =", phi)

        if verbose:
            print(" tree is over")
            print(" alpha =", alpha)
            print(" phi =", phi)

        for pp in range(a+p,2*a+2,2):
            if verbose: print(" paired poles (%d,%d)" % (pp, pp+1))
            alpha[pp] = phi[pp] = pp+1
            alpha[pp+1] = phi[pp+1] = pp
            if verbose > 1:
                print(" alpha =", alpha)
                print(" phi   =", phi)

        assert len(L) == l, "This may not happen"

        # yield the canonical sepx. diag
        if verbose:
            print(" =" * (3*a+7))
            print(" sigma =", sigma)
            print(" alpha =", alpha)
            print(" phi   =", phi)
            print(" =" * (3*a+7))

        c = orientation_cover(alpha,phi,a,verbose=verbose)
        c.canonical_label(inplace=True)
        if c not in cyl_diags:
            c_sym = [c]
            cc = c.horizontal_symmetry()
            cc.canonical_label(inplace=True)
            c_sym.append(cc)

            cc = c.vertical_symmetry()
            cc.canonical_label(inplace=True)
            c_sym.append(cc)
            cyl_diags.update(c_sym)
            yield c

        # Make the poles vary among the leaves
        #TODO: optimization when tree has nontrivial cyclic symmetry
        if p != 0 and p != l:
            if verbose:
                print(" start revolving door(%d,%d)" % (l, p))
                print(" leaves are at separatrices", L)
            for i,j in combinations(l,p):
                i = L[i]
                j = L[j]
                if verbose > 1:
                    print(" revolve i=%d j=%d" % (i, j))
                a_i = alpha[i]
                a_j = alpha[j]
                s_i = sigma[i]
                s_a_j = sigma[a_j]
                ss_i = phi[alpha[i]]  # sigma^-1(i)
                ss_j = phi[alpha[j]]  # sigma^-1(j)
                a_s_i = alpha[s_i]
                a_s_a_j = alpha[s_a_j]

                assert sigma[j] == a_j, "sigma[%d] = %d != alpha[%d]"%(j,sigma[j],a_j)
                assert phi[i] == a_i, "phi[%d] != alpha[i]"%(i,)
                assert phi[a_s_i] == i, "phi[%d] != %d"%(a_s_i,i)
                assert phi[j] == j, "phi[%d] + %d"%(j,j)
                assert phi[a_s_a_j] == a_j, "phi[%d] != alpha[%d]"%(a_s_a_j,j)

                alpha[i]     = a_j
                alpha[a_j]   = i
                alpha[j]     = a_i
                alpha[a_i]   = j

                sigma[i]     = a_j    # old_sigma[j]
                sigma[a_j]   = s_i    # old_sigma[i]
                sigma[j]     = s_a_j  # old_sigma[a_j]

                phi[i]       = i      # old_phi[a_s_i]
                phi[j]       = a_i    # old_phi[i]

                if s_i != j: # and a_s_i == a_j
                    phi[a_s_i]   = a_j    # old_phi[a_s_a_j]
                    phi[a_i]     = ss_j   # old_phi[a_j]
                else:
                    phi[a_i] = a_j

                if s_a_j != i:
                    phi[a_j]     = ss_i   # old_phi[a_i]
                    phi[a_s_a_j] = j      # old_phi[j]
                else:
                    phi[a_j] = j

                if verbose:
                    print(" =" * (3*a+7))
                    print(" sigma =", sigma)
                    print(" alpha =", alpha)
                    print(" phi   =", phi)
                    print(" =" * (3*a+7))

                for i in range(2*a+2):
                    ii = phi[alpha[sigma[i]]]
                    assert ii == i, "f_a_s(%d) == %d != %d"%(i,ii,i)
                for i in range(a,2*a+2):
                    assert sigma[i] == i, "sigma[%d] = %d != %d"%(i,sigma[i],i)
                c = orientation_cover(alpha,phi,a,verbose=verbose)
                c.canonical_label(inplace=True)
                if c not in cyl_diags:
                    c_sym = [c]
                    cc = c.horizontal_symmetry()
                    cc.canonical_label(inplace=True)
                    c_sym.append(cc)

                    cc = c.vertical_symmetry()
                    cc.canonical_label(inplace=True)
                    c_sym.append(cc)
                    cyl_diags.update(c_sym)
                    yield c

            # reinitialize sigma
            sigma = list(range(1,a)) + [0] + list(range(a,2*a+2))

@total_ordering
class CylinderDiagram(SeparatrixDiagram):
    r"""
    Separatrix diagram with pairing.

    Each cylinder is stored as a couple (bot,top) for which the orientation is
    as follows::

        +--------------------+
        |     <-- top --     |
        |                    |
        |                    |
        |      -- bot -->    |
        +--------------------+

    INPUT:

    - ``data`` - list of 2-tuples - matching of bottom-top pairs

    EXAMPLES::

        sage: from surface_dynamics import *


    We first build the simplest cylinder diagram which corresponds to a torus::

        sage: CylinderDiagram([((0,),(0,))])
        (0)-(0)

    The same initialized from a string::

        sage: CylinderDiagram('(0)-(0)')
        (0)-(0)

    The following initialize a cylinder diagram with two cylinder which gives a
    surface of genus 2 with one singularity of degree 2::

        sage: CylinderDiagram([((0,1),(0,2)),((2,),(1,))])
        (0,1)-(0,2) (2)-(1)

    ALGORITHM:

    A cylinder is represented by a couple (i,j) where i is the min in bot and j
    is the min in top. The data _top_to_cyl and _bot_to_cyl corresponds to the
    association of a separatrix to the corresponding 2-tuple. The twist
    coordinate correspond to the shift between those two indices.
    """
    def __init__(self, data, check=True):
        r"""
        TESTS::

            sage: from surface_dynamics import *

            sage: c = CylinderDiagram([((0,),(0,))])
            sage: CylinderDiagram(str(c)) == c
            True
            sage: loads(dumps(c)) == c
            True

            sage: CylinderDiagram('(0,1)-(2,1) (2)-(3)')
            Traceback (most recent call last):
            ...
            ValueError: each separatrix must appear once on bot and once on top
        """
        bot = []
        top = []

        if isinstance(data,str):
            data = [(string_to_cycle(b),string_to_cycle(t)) for b,t in (w.split('-') for w in data.split(' '))]
        else:
            data = list(data)

        for b,t in data:
            bot.append(tuple(b))
            top.append(tuple(t))

        SeparatrixDiagram.__init__(self, tuple(bot), tuple(top), check=False)

        if sum(len(x) for x,y in data) != self.nseps() or\
           sum(len(y) for x,y in data) != self.nseps():
            raise ValueError('each separatrix must appear once on bot and once on top')

        b2c = [None] * self.nseps() # bot separatrix -> cylinder (bot_min_index, top_min_index)
        t2c = [None] * self.nseps() # top separatrix -> cylinder (bot_min_index, top_min_index)
        for b,t in data:
            cyl = (min(b),min(t))
            for j in b: b2c[j] = cyl
            for j in t: t2c[j] = cyl

        self._bot_to_cyl = b2c
        self._top_to_cyl = t2c

        self._check()

    def __hash__(self):
        r"""
        Hash value: This is bad since we can modify it inplace!!
        """
        return hash(tuple(self._bot_to_cyl + self._top_to_cyl +
                          self._bot_cycles + self._top_cycles))

    def _check(self):
        SeparatrixDiagram(self)._check()
        b2c = self._bot_to_cyl
        t2c = self._top_to_cyl

        for j in range(self.nseps()):
            bo = self.bot_orbit(b2c[j][0])
            if j not in bo or b2c[j][0] != min(bo):
                raise ValueError("invalid data: j={} bo={} b2c[j]={}".format(j, bo, b2c[j]))
            to = self.top_orbit(t2c[j][1])
            if j not in to or t2c[j][1] != min(to):
                raise ValueError("invalid data: j={} to={} t2c[j]={}".format(j, to, t2c[j]))

    def weighted_adjacency_matrix(self, canonicalize=True):
        r"""
        Return the weighted adjacency matrix of this cylinder diagram.

        There is a vertex per cylinder and the weight from cyl1 to cyl2
        is the number of saddle connections which are on top of cyl1 and
        in the bottom of cyl2.

        EXAMPLES::

            sage: from surface_dynamics import AbelianStratum
            sage: for cd in AbelianStratum(1,1).cylinder_diagrams():
            ....:     print(cd.weighted_adjacency_matrix())
            [4]
            [0 1]
            [1 2]
            [1 1]
            [1 1]
            [0 0 1]
            [0 0 1]
            [1 1 0]
        """
        from sage.rings.all import ZZ
        from sage.matrix.constructor import matrix

        A = matrix(ZZ, self.ncyls())
        cylinder_indices = {}
        for key in self._bot_to_cyl:
            if key not in cylinder_indices:
                cylinder_indices[key] = len(cylinder_indices)
        for i in range(self.nseps()):
            bot = cylinder_indices[self._bot_to_cyl[i]]
            top = cylinder_indices[self._top_to_cyl[i]]
            A[bot,top] += 1

        if canonicalize:
            from sage.graphs.digraph import DiGraph
            D = DiGraph(A, format='weighted_adjacency_matrix')
            A = D.canonical_label(edge_labels=True).weighted_adjacency_matrix()

        A.set_immutable()

        return A

    def lengths_minimal_solution(self):
        r"""
        Return an integral solution for the lengths equation with minimal sum
        """
        from sage.numerical.mip import MixedIntegerLinearProgram
        M = MixedIntegerLinearProgram(maximization=False)
        x = M.new_variable(integer=True)
        M.set_objective(M.sum(x[i] for i in range(self.nseps())))
        for i in range(self.nseps()):
            M.add_constraint(x[i] >= 1)
        for bot,top in self.cylinders():
            M.add_constraint(M.sum(x[i] for i in bot) == M.sum(x[i] for i in top))
        M.solve()
        v = M.get_values(x)
        return [v[i] for i in range(self.nseps())]


    def tikz_string(self, cyl_height=0.8, cyl_sep=1.3, cyl_pos=None, sep_labels=None):
        r"""
        INPUT:

        - ``cyl_height`` - (optional) height of cylinders

        - ``cyl_sep`` - (optional) vertical space between cylinders

        - ``cyl_pos`` - (optional) position of left corners of cylinders. If provided, ``cyl_sep``
          is ignored.

        - ``sep_labels`` - labels of the separatrices
        """
        if cyl_pos is None:
            cyl_pos = [(0, i*(cyl_height+cyl_sep)) for i in range(self.nseps())]
        if sep_labels is None:
            sep_labels = range(self.nseps())

        v = self.lengths_minimal_solution()
        s = []
        for i,(bot,top) in enumerate(self.cylinders()):
            w = sum(v[i] for i in bot)
            assert w == sum(v[i] for i in top)

            xl, yl = cyl_pos[i]

            # fill the cylinder
            s.append('\\fill[gray!20] ({},{}) rectangle ({},{});'.format(
                xl,yl,xl+w,yl+cyl_height))
            s.append('\\draw ({},{}) rectangle ({},{});'.format(
                xl,yl,xl+w,yl+cyl_height))

            # singularities
            x = xl
            for i in bot:
                s.append('\\fill ({},{}) circle (2pt);'.format(x, yl))
                s.append('\\node[below=-2] at ({},{}) {{{}}};'.format(x+v[i]/2, yl, sep_labels[i]))
                x += v[i]
            s.append('\\fill ({},{}) circle (2pt);'.format(x, yl))

            x = xl+w
            for i in top:
                s.append('\\fill ({},{}) circle (2pt);'.format(x, yl+cyl_height))
                s.append('\\node[above=-2] at ({},{}) {{{}}};'.format(x-v[i]/2, yl+cyl_height, sep_labels[i]))
                x -= v[i]
            s.append('\\fill ({},{}) circle (2pt);'.format(x, yl+cyl_height))

        return '\n'.join(s)

    def _repr_(self):
        r"""
        String representation

        TESTS::

            sage: from surface_dynamics import *

            sage: c = CylinderDiagram([((0,1),(1,2)),((2,),(0,))])
            sage: repr(c) #indirect doctest
            '(0,1)-(1,2) (2)-(0)'
        """
        l = []
        for b,t in self.cylinders():
            l.append('(' + ','.join(map(str,b)) + ')-(' + ','.join(map(str,t)) + ')')
        return ' '.join(l)

    def __eq__(self, other):
        r"""
        TESTS::

            sage: from surface_dynamics import AbelianStratum, CylinderDiagram

            sage: c1 = CylinderDiagram('(0,5)-(0,4) (1,4)-(1,3) (2,3)-(2,5)')
            sage: c2 = CylinderDiagram('(0,5)-(1,3) (1,4)-(0,4) (2,3)-(2,5)')
            sage: c3 = CylinderDiagram('(0,5)-(2,5) (1,4)-(1,3) (2,3)-(0,4)')
            sage: c1 == c2 or c2 == c3 or c1 == c3
            False
            sage: c1 != c2 and c2 != c3 and c1 != c3
            True

            sage: from operator import not_
            sage: C = AbelianStratum(4).cylinder_diagrams()
            sage: for c1 in C:
            ....:     assert c1 == c1
            ....:     assert not_(c1 != c1)
            ....:     for c2 in C:
            ....:         assert (c1 == c2) == (c2 == c1) == not_(c2 != c1) == not_(c1 != c2)
        """
        if type(self) is not type(other):
            raise TypeError
        return SeparatrixDiagram.__eq__(self, other) and self._bot_to_cyl == other._bot_to_cyl

    def __ne__(self, other):
        if type(self) is not type(other):
            raise TypeError
        return not (self == other)

    def __lt__(self, other):
        r"""
        Comparison

        TESTS::

            sage: from surface_dynamics import AbelianStratum
            sage: from operator import not_
            sage: C = AbelianStratum(4).cylinder_diagrams()
            sage: for c1 in C:
            ....:     for c2 in C:
            ....:         assert (c1 < c2) == (c2 > c1) == not_(c1 >= c2) == not_(c2 <= c1)
            ....:         assert (c1 <= c2) == (c2 >= c1) == not_(c1 > c2) == not_(c2 < c1)
        """
        if type(self) is not type(other):
            raise TypeError

        if SeparatrixDiagram.__lt__(self, other):
            return True
        if SeparatrixDiagram.__lt__(other, self):
            return False

        if self._bot_to_cyl < other._bot_to_cyl:
            return True
        elif self._bot_to_cyl > other._bot_to_cyl:
            return False

        if self._top_to_cyl < other._top_to_cyl:
            return True
        elif self._top_to_cyl > other._top_to_cyl:
            return False

        # equality
        return False

    #
    # access to attribute
    #


    def to_directed_graph(self):
        r"""
        Return a labeled directed graph that encodes the cylinder diagram.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: c = CylinderDiagram('(0,1,5)-(2,5) (2)-(0,1,3) (3,4)-(4)'); c
            (0,1,5)-(2,5) (2)-(0,1,3) (3,4)-(4)
            sage: G = c.to_directed_graph(); G
            Looped multi-digraph on 6 vertices
            sage: G.edges(sort=True)
            [(0, 1, 'b'), (0, 1, 't'), (0, 2, 'c'), (0, 5, 'c'), (1, 2, 'c'), (1, 3, 't'), (1, 5, 'b'), (1, 5, 'c'), (2, 0, 'c'), (2, 1, 'c'), (2, 2, 'b'), (2, 3, 'c'), (2, 5, 't'), (3, 0, 't'), (3, 4, 'b'), (3, 4, 'c'), (4, 3, 'b'), (4, 4, 'c'), (4, 4, 't'), (5, 0, 'b'), (5, 2, 'c'), (5, 2, 't'), (5, 5, 'c')]
        """
        G = SeparatrixDiagram.to_directed_graph(self)

        for cb,ct in self.cylinders():
            for i in cb:
                for j in ct:
                    G.add_edge(i,j,'c')

        return G

    def is_isomorphic(self, other, return_map=False):
        r"""
        Test whether this cylinder diagram is isomorphic to ``other``.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: c1 = CylinderDiagram('(0,2,1)-(3,4,5) (3)-(1) (4)-(2) (5)-(0)')
            sage: c2 = CylinderDiagram('(0,2,1)-(3,5,4) (3)-(1) (4)-(2) (5)-(0)')
            sage: c1.is_isomorphic(c2)
            False

            sage: c3 = CylinderDiagram('(1,3,5)-(2,4,0) (2)-(5) (4)-(3) (0)-(1)')
            sage: c1.is_isomorphic(c3)
            True
            sage: ans, m = c1.is_isomorphic(c3, return_map=True)
            sage: assert ans is True and c1.relabel(m) == c3

            sage: c4 = CylinderDiagram('(4,2,3)-(1,5,0) (1)-(3) (0)-(2) (5)-(4)')
            sage: c2.is_isomorphic(c4)
            True
            sage: ans, m = c2.is_isomorphic(c4, return_map=True)
            sage: assert ans is True and c2.relabel(m) == c4

            sage: c1 = CylinderDiagram('(0,5)-(3,4) (1,4)-(2,5) (2)-(0) (3)-(1)')
            sage: c2 = CylinderDiagram('(0,5)-(3,4) (1,4)-(2,5) (2)-(1) (3)-(0)')
            sage: c1.is_isomorphic(c2)
            False

            sage: c1 = CylinderDiagram('(0,1)-(4,5) (2,4)-(0,3) (3,5)-(1,2)')
            sage: c2 = CylinderDiagram('(0,1)-(2,3) (2,4)-(1,4) (3,5)-(0,5)')
            sage: c3 = CylinderDiagram('(0,5)-(2,5) (1,4)-(0,4) (2,3)-(1,3)')
            sage: c1.is_isomorphic(c2)
            False
            sage: c1.is_isomorphic(c3)
            False
            sage: c2.is_isomorphic(c3)
            False

        TESTS::

            sage: from surface_dynamics import CylinderDiagram
            sage: c = CylinderDiagram('(0,3)-(3,7) (1,6,4)-(0,2,4,6,8) (2,8,7,5)-(1,5)')
            sage: c1 = c.relabel([4,7,0,3,2,1,6,8,5])
            sage: c2 = c.relabel([6,3,1,8,2,0,4,5,7])
            sage: assert c.is_isomorphic(c1) and c.is_isomorphic(c2) and c1.is_isomorphic(c2)
            sage: _, m1 = c1.is_isomorphic(c, return_map=True)
            sage: assert c1.relabel(m1) == c
            sage: _, m2 = c2.is_isomorphic(c, return_map=True)
            sage: assert c2.relabel(m2) == c
        """
        if type(self) is not type(other):
            raise TypeError
        if return_map:
            s, sm = self.canonical_label(return_map=True)
            o, om = other.canonical_label(return_map=True)
            if s != o:
                return (False, None)
            else:
                ominv = {j:i for i,j in om.items()}
                return (True, {i: ominv[sm[i]] for i in sm})
        else:
            return self.canonical_label() == other.canonical_label()

    def relabel(self, perm, inplace=False):
        r"""
        EXAMPLES::

            sage: from surface_dynamics import CylinderDiagram
            sage: c = CylinderDiagram('(0,3)-(3,7) (1,6,4)-(0,2,4,6,8) (2,8,7,5)-(1,5)')
            sage: c.relabel([3,1,2,6,4,5,0,7,8])
            (0,4,1)-(0,8,3,2,4) (2,8,7,5)-(1,5) (3,6)-(6,7)
            sage: c
            (0,3)-(3,7) (1,6,4)-(0,2,4,6,8) (2,8,7,5)-(1,5)
            sage: c.relabel([3,1,2,6,4,5,0,7,8], inplace=True)
            (0,4,1)-(0,8,3,2,4) (2,8,7,5)-(1,5) (3,6)-(6,7)
            sage: c
            (0,4,1)-(0,8,3,2,4) (2,8,7,5)-(1,5) (3,6)-(6,7)

        TESTS::

            sage: c = CylinderDiagram('(0,3)-(3,7) (1,6,4)-(0,2,4,6,8) (2,8,7,5)-(1,5)')
            sage: c1 = c.relabel([3,1,2,6,4,5,0,7,8], inplace=False)
            sage: c2 = c.relabel([3,1,2,6,4,5,0,7,8], inplace=True)
            sage: assert c is c2 and c1 == c2
            sage: c1._check()
            sage: c2._check()

            sage: c1 = c.relabel([1,5,0,6,8,3,2,7,4], inplace=False)
            sage: c2 = c.relabel([1,5,0,6,8,3,2,7,4], inplace=True)
            sage: assert c is c2 and c1 == c2
            sage: c1._check()
            sage: c2._check()
        """
        if inplace:
            n = self.nseps()

            # find new (bot min, top min) in each cylinder
            new_mins = {}
            for b, t in self.cylinders():
                bb = [perm[i] for i in b]
                tt = [perm[i] for i in t]

                new_mins[(min(b), min(t))] = (min(bb), min(tt))

            b2c = [None] * n
            t2c = [None] * n
            for i in range(n):
                b2c[perm[i]] = new_mins[self._bot_to_cyl[i]]
                t2c[perm[i]] = new_mins[self._top_to_cyl[i]]

            self._bot_to_cyl = b2c
            self._top_to_cyl = t2c

            SeparatrixDiagram.relabel(self, perm, True)

            return self

        else:
            return CylinderDiagram([(tuple(perm[i] for i in b), tuple(perm[i] for i in t)) for b,t in self.cylinders()])

    def canonical_label(self, inplace=False, return_map=False):
        r"""
        Return a cylinder diagram with canonical labels.

        .. NOTE::

            The canonical label might change depending on your Sage version and the
            optional packages available.

        EXAMPLES::

            sage: from surface_dynamics import CylinderDiagram
            sage: c1 = CylinderDiagram('(0,5,4)-(0,3,2,1) (1,3,2)-(4,5)')
            sage: c1.canonical_label()   # random
            (0,4,5)-(1,3,2,5) (1,3,2)-(0,4)
            sage: c2 = c1.relabel([2,4,5,1,3,0])
            sage: c3 = c1.relabel([5,3,0,1,4,2])
            sage: c1.canonical_label() == c2.canonical_label() == c3.canonical_label()
            True

        TESTS::

            sage: import itertools
            sage: from surface_dynamics import CylinderDiagram
            sage: c = CylinderDiagram('(0)-(1) (1,2)-(0,3) (3)-(2)')
            sage: can = c.canonical_label()
            sage: for p in itertools.permutations([0,1,2,3]):
            ....:    cc = c.relabel(p)
            ....:    ccan, m = cc.canonical_label(return_map=True)
            ....:    assert cc.relabel(m) == can == ccan

            sage: c = CylinderDiagram('(0,4)-(0,3) (1,3)-(1,5) (2,5)-(2,4)')
            sage: can = c.canonical_label()
            sage: for p in itertools.permutations([0,1,2,3,4,5]):
            ....:    cc = c.relabel(p)
            ....:    ccan, m = cc.canonical_label(return_map=True)
            ....:    assert cc.relabel(m) == can == ccan

            sage: c = CylinderDiagram('(0,1)-(0,2) (3,5,4)-(1,4,6) (2,6)-(3,5)')
            sage: c is c.canonical_label()
            False
            sage: c.canonical_label() is c.canonical_label()
            True
            sage: c.canonical_label().canonical_label() is c.canonical_label()
            True
        """
        if not hasattr(self, '_normal_form'):
            G = self.to_directed_graph()
            _, m = G.canonical_label(certificate=True, edge_labels=True)
            # m = [m[i] for i in range(self.nseps())]
            # GG the new digraph
            # m from the digraph to its canonic labels
            cyls = []
            for b,t in self.cylinders():
                cyls.append((tuple(m[i] for i in b), tuple(m[i] for i in t)))

            self._normal_form = CylinderDiagram(cyls, check=False)
            self._normal_labels = m

            self._normal_form._normal_form = self._normal_form
            self._normal_form._normal_labels = range(self.nseps())

        if inplace:
            self.__dict__ = self._normal_form.__dict__

        if return_map:
            return self._normal_form, self._normal_labels
        elif not inplace:
            return self._normal_form

    def separatrix_diagram(self):
        r"""
        Return the underlying separatrix diagram

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: s = SeparatrixDiagram('(0,1)(2,3,4)','(0,3)(1,4,2)'); s
            (0,1)(2,3,4)-(0,3)(1,4,2)
            sage: c = s.to_cylinder_diagram([(0,1),(1,0)]); c
            (0,1)-(1,4,2) (2,3,4)-(0,3)
            sage: c.separatrix_diagram() == s
            True
        """
        return SeparatrixDiagram(self._bot,self._top,check=False)

    def lengths_cone(self):
        r"""
        Return the rational polyhedron corresponding to the set of length with
        the given fixed heights.

        -> one can obtain ehrhart series for each of them! It tells us that we
        have a nice asymptotics... and the asymptotics is simply given by the
        volume of this polytope (up to the ignored twists parameters)!
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        from sage.rings.integer_ring import ZZ

        n = self.nseps()
        eqns = []
        ieqs = []

        # lengths are non-negative
        e = [ZZ.zero()] * (n+1)
        for i in range(n):
            e[i+1] = ZZ.one()
            ieqs.append(e[:])
            e[i+1] = ZZ.zero()

        # for each cylinder, length top = length bot
        for i,(bot,top) in enumerate(self.cylinders()):
            e = [ZZ.zero()] * (n+1)
            for j in bot:
                e[j+1] += ZZ.one()
            for j in top:
                e[j+1] -= ZZ.one()
            eqns.append(e)

        return Polyhedron(ieqs=ieqs, eqns=eqns)

    def widths_generating_series(self, var='w'):
        r"""
        Generating series of the number of lengths solutions given the widths.

        WARNING: when a triangulation is involved, the generating series ignore
        some lower dimensional polytopes!!

        INPUT:

        - ``var`` - (optional) an optional name for the variable

        EXAMPLES::

            sage: from surface_dynamics import CylinderDiagram

            sage: c = CylinderDiagram('(0,1)-(0,2) (2)-(1)')
            sage: c.widths_generating_series()  # optional: latte_int
            (1)/((1 - w0)*(1 - w0*w1))

            sage: c = CylinderDiagram('(0)-(2) (1,2,3)-(4,5) (4)-(3) (5)-(0,1)')
            sage: c.widths_generating_series()  # optional: latte_int
            (1)/((1 - w1*w3)*(1 - w1*w2)*(1 - w0*w1*w3))

            sage: c = CylinderDiagram('(0,1,3)-(0,2,5) (2,4)-(1,3) (5)-(4)')
            sage: c.widths_generating_series()  # optional: latte_int
            (1)/((1 - w0)*(1 - w0*w1)*(1 - w0*w1*w2)^2) + (1)/((1 - w0)*(1 - w0*w1)^2*(1 - w0*w1*w2))
        """
        from sage.matrix.constructor import matrix
        from surface_dynamics.misc.multiplicative_multivariate_generating_series import MultiplicativeMultivariateGeneratingSeriesRing
        from sage.geometry.polyhedron.constructor import Polyhedron

        from sage.misc.misc_c import prod

        import re
        re_var = re.compile('w(\\d+)')

        sub1 = [None] * self.nseps()
        sub2 = [None] * self.nseps()
        for i,(top,bot) in enumerate(self.cylinders()):
            for j in top:
                sub1[j] = i
            for j in bot:
                sub2[j] = i

        M = MultiplicativeMultivariateGeneratingSeriesRing(self.ncyls(), 'w')
        V = M.free_module()
        ans = M.zero()
        dim = self.nseps() - self.ncyls() + 1
        for rays in cone_triangulate(self.lengths_cone()):
            assert len(rays) == dim, (self, dim)
            num = simplex_count(rays)

            den1 = {}
            den2 = {}
            for r in rays:
                l1 = [0] * self.ncyls()
                l2 = [0] * self.ncyls()
                for j,c in enumerate(r):
                    l1[sub1[j]] += c
                    l2[sub2[j]] += c
                v1 = V(l1); v1.set_immutable()
                v2 = V(l2); v2.set_immutable()
                if v1 in den1:
                    den1[v1] += 1
                else:
                    den1[v1] = 1
                if v2 in den2:
                    den2[v2] += 1
                else:
                    den2[v2] = 1

            f1 = M.term(num, list(den1.items()))
            f2 = M.term(num, list(den2.items()))

            assert f1 == f2, (f1, f2)
            degrees = set(f1.degrees())
            assert degrees == set([dim]), (self, degrees, dim)


            ans += f1

        return ans

    def volume_contribution(self):
        r"""
        Return the volume contribution of this cylinder diagram as a
        generalized multiple zeta values.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: c0, c1 = AbelianStratum(2).cylinder_diagrams()
            sage: v0 = c0.volume_contribution()  # optional: latte_int
            sage: v0                             # optional: latte_int
            (1/3)/((w)^4)
            sage: v0.integral_sum_as_mzv()       # optional: latte_int
            1/3*ζ(4)
            sage: v1 = c1.volume_contribution()  # optional: latte_int
            sage: v1                             # optional: latte_int
            (2/3)/((w1)*(w0 + w1)^3) + (1/3)/((w1)^2*(w0 + w1)^2)
            sage: v1.integral_sum_as_mzv()       # optional: latte_int
            2/3*ζ(1,3) + 1/3*ζ(2,2)

            sage: for c in AbelianStratum(1,1).cylinder_diagrams():  # optional: latte_int
            ....:     print(c, c.volume_contribution().integral_sum_as_mzv())
            (0,3,1,2)-(0,3,1,2) 1/6*ζ(5)
            (0)-(1) (1,2,3)-(0,2,3) 1/3*ζ(2,3) + 1/3*ζ(3,2)
            (0,3)-(1,3) (1,2)-(0,2) ζ(1,4) + 1/3*ζ(2,3)
            (0,1)-(2,3) (2)-(1) (3)-(0) 1/3*ζ(1,3) + 1/3*ζ(2,2) - 1/3*ζ(2,3) - 1/3*ζ(3,2) + 1/3*ζ(4) - 1/3*ζ(5)

            sage: sum(c.volume_contribution() for c in AbelianStratum(2,1,1).cylinder_diagrams(1)).integral_sum_as_mzv()  # optional: latte_int
            7/180*ζ(8)

        Detailed contribution of 2 cylinder diagrams::

            sage: cyls = AbelianStratum(2,1,1).cylinder_diagrams(2)
            sage: sum(cyls[k].volume_contribution() for k in [2,7,8,21,22]).integral_sum_as_mzv()  # optional: latte_int
            13/630*ζ(5,3) + 13/252*ζ(6,2)
            sage: sum(cyls[k].volume_contribution() for k in [0,11,19,20]).integral_sum_as_mzv()  # optional: latte_int
            1/21*ζ(4,4) + 4/63*ζ(5,3)
            sage: sum(cyls[k].volume_contribution() for k in [3,10]).integral_sum_as_mzv()  # optional: latte_int
            2/35*ζ(3,5) + 3/70*ζ(4,4)
            sage: sum(cyls[k].volume_contribution() for k in [13,23,24]).integral_sum_as_mzv()  # optional: latte_int
            2/21*ζ(2,6) + 4/105*ζ(3,5)
            sage: sum(cyls[k].volume_contribution() for k in [9]).integral_sum_as_mzv()  # optional: latte_int
            1/21*ζ(1,7) + 1/126*ζ(2,6)
            sage: sum(cyls[k].volume_contribution() for k in [5]).integral_sum_as_mzv()  # optional: latte_int
            2/105*ζ(3,5) + 1/70*ζ(4,4)
            sage: sum(cyls[k].volume_contribution() for k in [6,26]).integral_sum_as_mzv()  # optional: latte_int
            1/7*ζ(1,7) + 1/14*ζ(2,6) + 1/21*ζ(3,5) + 1/28*ζ(4,4) + 2/105*ζ(5,3)
            sage: sum(cyls[k].volume_contribution() for k in [1,14]).integral_sum_as_mzv()  # optional: latte_int
            2/7*ζ(1,7) + 1/7*ζ(2,6) + 2/21*ζ(3,5) + 3/70*ζ(4,4)
            sage: sum(cyls[k].volume_contribution() for k in [17,27]).integral_sum_as_mzv()  # optional: latte_int
            2/7*ζ(1,7) + 1/7*ζ(2,6) + 4/105*ζ(3,5)
            sage: sum(cyls[k].volume_contribution() for k in [4,25]).integral_sum_as_mzv()  # optional: latte_int
            1/7*ζ(1,7) + 1/42*ζ(2,6)
            sage: sum(cyls[k].volume_contribution() for k in [12,28]).integral_sum_as_mzv()  # optional: latte_int
            4/21*ζ(1,7) + 2/21*ζ(2,6) + 8/315*ζ(3,5)
            sage: sum(cyls[k].volume_contribution() for k in [15,16]).integral_sum_as_mzv()  # optional: latte_int
            4/21*ζ(1,7) + 2/21*ζ(2,6) + 8/315*ζ(3,5)
            sage: sum(cyls[k].volume_contribution() for k in [18]).integral_sum_as_mzv()  # optional: latte_int
            1/36*ζ(7) - 1/36*ζ(8)
        """
        from surface_dynamics.misc.additive_multivariate_generating_series import AdditiveMultivariateGeneratingSeriesRing

        from sage.misc.misc_c import prod
        from sage.rings.integer_ring import ZZ


        # automorphism group
        aut_size = self.automorphism_group().cardinality()

        # equivalent sep diag in the component
        sym = ZZ(4) // (1 + sum(self.symmetries()))

        # numbering zeros
        mult = {}
        for i in perm_cycles(self.outgoing_edges_perm(), singletons=True):
            i = len(i)
            if i in mult:
                mult[i] += 1
            else:
                mult[i] = 1
        numbering = prod(ZZ(m).factorial() for m in mult.values())

        F = self.widths_generating_series()
        deg, res = F.delta().residue()
        assert deg == self.nseps() + 1, (deg, self.nseps(), F, res)

        A = F.parent().residue_ring()
        res = sum(A.term(num, den) for num,den in res)

        m = ZZ(2) * numbering * sym / ZZ(self.nseps()).factorial() / aut_size
        return m * res

    def cylinders(self):
        r"""
        Return the cylinders as a list of pairs ``(bot, top)``.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: c = CylinderDiagram('(0,2,4)-(1,3,5) (1,5)-(0) (3)-(2,4)')
            sage: c
            (0,2,4)-(1,3,5) (1,5)-(0) (3)-(2,4)
            sage: c.cylinders()
            [((0, 2, 4), (1, 3, 5)), ((1, 5), (0,)), ((3,), (2, 4))]
        """
        return [(b,self.top_orbit(self._bot_to_cyl[b[0]][1])) for b in self.bot_cycle_tuples()]

    def bot_to_cyl(self, j):
        r"""
        Return the cylinder above the ``j``-th separatrix.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: c = CylinderDiagram('(0,2,4)-(1,3,5) (1,5)-(0) (3)-(2,4)')
            sage: c
            (0,2,4)-(1,3,5) (1,5)-(0) (3)-(2,4)
            sage: c.bot_to_cyl(0)
            ((0, 2, 4), (1, 3, 5))
            sage: c.bot_to_cyl(1)
            ((1, 5), (0,))
            sage: c.bot_to_cyl(3)
            ((3,), (2, 4))
        """
        jb,jt = self._bot_to_cyl[j]
        return self.bot_orbit(jb), self.top_orbit(jt)

    def top_to_cyl(self, j):
        r"""
        Return the cylinder below the ``j``-th separatrix.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: c = CylinderDiagram('(0,2,4)-(1,3,5) (1,5)-(0) (3)-(2,4)')
            sage: c.top_to_cyl(0)
            ((1, 5), (0,))
            sage: c.top_to_cyl(2)
            ((3,), (2, 4))
        """
        jb,jt = self._top_to_cyl[j]
        return self.bot_orbit(jb), self.top_orbit(jt)

    #
    # properties
    #

    def is_connected(self):
        r"""
        Test the connectedness of this cylinder diagram.

        TESTS::

            sage: from surface_dynamics import *

            sage: CylinderDiagram('(0)-(1) (1)-(0)').is_connected()
            True
            sage: CylinderDiagram('(0,1)-(0) (2)-(1,2)').is_connected()
            True

            sage: CylinderDiagram('(0)-(0) (1)-(1)').is_connected()
            False
            sage: CylinderDiagram('(0,1)-(3) (2)-(2) (3)-(0,1)').is_connected()
            False
        """
        from sage.graphs.graph import Graph
        G = Graph(multiedges=True, loops=True)
        for b,t in self.cylinders():
            G.add_edges((b[0],b[j]) for j in range(1,len(b)))
            G.add_edges((t[0],t[j]) for j in range(1,len(t)))
            G.add_edge(b[0],t[0])
        return G.num_verts() == self.nseps() and G.is_connected()

    #
    # symmetries
    #

    def inverse(self):
        r"""
        Return the inverse cylinder diagram.

        The inverse of a cylinder diagram is the cylinder diagram in which all
        cylinders have been reversed. It corresponds to the multiplication by
        `-1` on the underlying Abelian differential.

        Combinatorially the operation is b0-t0 ... bk-tk becomes t0-b0 ... tk-bk

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: c = CylinderDiagram('(0,1)-(0,2) (3,5,4)-(1,4,6) (2,6)-(3,5)')
            sage: c
            (0,1)-(0,2) (2,6)-(3,5) (3,5,4)-(1,4,6)
            sage: c.inverse()
            (0,2)-(0,1) (1,4,6)-(3,5,4) (3,5)-(2,6)

        The operation can also be defined at the level of the separatrix
        diagrams and the two operation commutes::

            sage: c.separatrix_diagram().inverse() == c.inverse().separatrix_diagram()
            True

        The inversion can also be seen as the composition of the horizontal and
        vertical symmetries::

            sage: c.horizontal_symmetry().vertical_symmetry() == c.inverse()
            True
            sage: c.vertical_symmetry().horizontal_symmetry() == c.inverse()
            True

        The inversion is an involution on cylinder diagrams::

            sage: all(cc.inverse().inverse() == cc for cc in AbelianStratum(4).cylinder_diagrams()) # long time
            True
        """
        return CylinderDiagram([(t,b) for (b,t) in self.cylinders()])

    def vertical_symmetry(self):
        r"""
        Return the cylinder diagram obtained by reflecting the cylinder
        configuration along the vertical axis.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: c = CylinderDiagram('(0,3,4)-(0,3,5) (1,2,5)-(1,2,4)')
            sage: c.vertical_symmetry()
            (0,4,3)-(0,5,3) (1,5,2)-(1,4,2)

            sage: c.separatrix_diagram().vertical_symmetry() == c.vertical_symmetry().separatrix_diagram()
            True

            sage: A = AbelianStratum(2,2)
            sage: all(c.vertical_symmetry().stratum() == A for c in A.cylinder_diagrams())
            True
        """
        return CylinderDiagram(tuple((b[::-1],t[::-1]) for b,t in self.cylinders()))

    def horizontal_symmetry(self):
        r"""
        Return the cylinder diagram obtained by reflecting the cylinder
        configuration along the horizontal axis.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: c = CylinderDiagram('(0,3,4)-(0,3,5) (1,2,5)-(1,2,4)')
            sage: c.horizontal_symmetry()
            (0,5,3)-(0,4,3) (1,4,2)-(1,5,2)

            sage: c.separatrix_diagram().horizontal_symmetry() == c.horizontal_symmetry().separatrix_diagram()
            True

            sage: A = AbelianStratum(2,2)
            sage: all(c.horizontal_symmetry().stratum() == A for c in A.cylinder_diagrams())
            True
        """
        return CylinderDiagram(tuple((t[::-1],b[::-1]) for b,t in self.cylinders()))

    def symmetries(self):
        r"""
        Return a triple ``(horiz_sym, vert_sym, inv_sym)``

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: c = CylinderDiagram('(0,1)-(2,3,5) (2,3,4)-(1) (5)-(0,4)')
            sage: c.symmetries()
            (False, True, False)

            sage: c.horizontal_symmetry().is_isomorphic(c)
            False
            sage: c.vertical_symmetry().is_isomorphic(c)
            True
            sage: c.inverse().is_isomorphic(c)
            False

            sage: CylinderDiagram('(0,2,1,4,3)-(0,4,2,1,3)').symmetries()
            (True, True, True)
            sage: CylinderDiagram('(0,3)-(0,4,5) (1,4,2)-(1,3) (5)-(2)').symmetries()
            (False, False, True)
            sage: CylinderDiagram('(0,2,3,1)-(0,2,1,4) (4)-(3)').symmetries()
            (True, False, False)
            sage: CylinderDiagram('(0,1,2)-(0,3,1,4) (3,4)-(2)').symmetries()
            (False, True, False)
            sage: CylinderDiagram('(0,1)-(0,3,4) (2,3)-(1) (4)-(2)').symmetries()
            (False, False, False)
        """
        n = len(self._top)

        # we first consider the separatrix diagram as it is much faster
        bot, top = two_non_connected_perms_canonical_labels(self._top, self._bot)

        # compute the inverses
        ibot = [None]*n
        itop = [None]*n
        for i in range(n):
            ibot[bot[i]] = i
            itop[top[i]] = i

        # horiz
        bot1, top1 = two_non_connected_perms_canonical_labels(itop, ibot)
        sep_horiz_sym = bot == bot1 and top == top1

        # vert
        bot1, top1 = two_non_connected_perms_canonical_labels(ibot, itop)
        sep_vert_sym = bot == bot1 and top == top1

        # inv
        if sep_horiz_sym and sep_vert_sym:  # got the two
            sep_inverse_sym = True
        elif sep_horiz_sym^sep_vert_sym:    # got exactly one
            sep_inverse_sym = False
        else:                       # none of them
            bot1, top1 = two_non_connected_perms_canonical_labels(top, bot)
            sep_inverse_sym = bot == bot1 and top == top1

        # next we check the cylinder diagram if needed
        if sep_horiz_sym:
            c1 = self.canonical_label(inplace=False)
            c2 = self.horizontal_symmetry().canonical_label(inplace=False)
            horiz_sym = c1 == c2
        else:
            horiz_sym = False

        if sep_vert_sym:
            c1 = self.canonical_label(inplace=False)
            c2 = self.vertical_symmetry().canonical_label(inplace=False)
            vert_sym = c1 == c2
        else:
            vert_sym = False

        if horiz_sym and vert_sym:  # got the two
            inverse_sym = True
        elif horiz_sym^vert_sym:    # got exactly one
            inverse_sym = False
        else:                       # none of them
            c1 = self.canonical_label(inplace=False)
            c2 = self.inverse().canonical_label(inplace=False)
            inverse_sym = c1 == c2

        return (horiz_sym, vert_sym, inverse_sym)


    def automorphism_group(self, order=False):
        r"""
        Return the automorphism group

        INPUT:

        - ``order`` - boolean (default: False) - whether or not return the order
          of the group

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: cyl = CylinderDiagram('(0,1)-(0,2) (2,3)-(1,3)')
            sage: cyl.automorphism_group()
            Permutation Group with generators [(0,3)(1,2)]
        """
        return self.to_directed_graph().automorphism_group(edge_labels=True, order=order)

    def is_hyperelliptic(self,verbose=False):
        r"""
        Test of hyperellipticity

        Each stratum of Abelian differentials as up to three connected
        components. For the strata H(2g-2) and H(g-1,g-1) there is a special
        component called *hyperelliptic* in which all translation surfaces
        `(X,\omega)` in that component are such that `X` is hyperelliptic.

        This function returns True if and only if the cylinder diagrams
        correspond to a decomposition of a surface associated to the
        hyperelliptic components in H(2g-2) or H(g-1,g-1).

        EXAMPLES::

            sage: from surface_dynamics import *

        In genus 2, strata H(2) and H(1,1), all surfaces are hyperelliptic::

            sage: for c in AbelianStratum(2).cylinder_diagrams():
            ....:     print("%d %s" % (c.ncyls(), c.is_hyperelliptic()))
            1 True
            2 True

            sage: for c in AbelianStratum(1,1).cylinder_diagrams():
            ....:     print("%d %s" % (c.ncyls(), c.is_hyperelliptic()))
            1 True
            2 True
            2 True
            3 True

        In higher genera, some of them are, some of them are not::

            sage: C = AbelianStratum(4).cylinder_diagrams()
            sage: len(C)
            15
            sage: sum(c.is_hyperelliptic() for c in C)
            5

            sage: C = AbelianStratum(2,2).cylinder_diagrams()
            sage: len(C)
            41
            sage: sum(c.is_hyperelliptic() for c in C)
            11
        """
        z = self.stratum().zeros()
        if z == [0] or z == [2] or z == [1,1]: return True
        if 0 in z:
            raise NotImplementedError("is_hyperelliptic method not implemented for cylinder diagrams with fake zeros")
        ns = self.nseps()
        if len(z) == 1: # minimal stratum H(2g-2)
            for cy in self.cylinders():
                if len(cy[0]) != len(cy[1]): return False
            b = self.bot()
            t = self.top()
            # build list of seps in cyclic order around zero, starting by outgoing sep 0
            lout = [0]
            lin = []
            for _ in range(ns):
                lin.append(t[lout[-1]])
                lout.append(b[lin[-1]])
            if verbose:
                print('lin  ', lin)
                print('lout', lout)
            # build involution on separatrices
            p = [None]*ns
            for a in range(ns):
                p[lout[a]] = lin[(a+ns//2)%ns]
            if verbose: print("involution on seps", p)
            # wsep = counter of sepatrices with a wpt
            wsep = 0
            for cy in self.cylinders():
                for k in cy[0]:
                    # check that p(k) is on the top of the cyl that has k on its bottom
                    if p[k] not in cy[1]: return False
                    # check that if k is on bot and top of cyl, then p(k) = k
                    if k in cy[1]:
                        if k != p[k]: return False
                        wsep += 1
            if verbose: print("wsep", wsep)
            # check number of w pts
            if wsep + 2*self.ncyls() != z[0] + 3: return False
            # check that cylinders are stable under involution
            if self != CylinderDiagram(
                [(cy[0],tuple(map(lambda x: p[x],cy[0]))) for cy in self.cylinders()]):
                return False
            return True
        elif len(z) == 2: # should be stratum H(g-1,g-1)
            if z[0] != z[1]: return False
            for cy in self.cylinders():
                if len(cy[0]) != len(cy[1]): return False
            b = self.bot()
            t = self.top()
            # build list of seps in cyclic order around first zero, starting by outgoing sep 0
            lout = [0]
            lin = []
            for _ in range(ns//2):
                lin.append(t[lout[-1]])
                lout.append(b[lin[-1]])
            if verbose:
                print('lin  ', lin)
                print('lout', lout)
            # build list of seps in cyclic order around the other zero
            a = 0
            while a in lout:
                a += 1
            llout = [a]
            llin = []
            for _ in range(ns//2):
                llin.append(t[llout[-1]])
                llout.append(b[llin[-1]])
            if verbose:
                print('llin  ', llin)
                print('llout', llout)
            # now, try each way the involution could send lout to llout
            for j in range(ns//2):
                test = True
                # build involution on separatrices
                p = [None]*ns
                for a in range(ns//2):
                    p[lout[a]] = llin[(j+a)%(ns//2)]
                    p[llout[a]] = lin[(a-j-1)%(ns//2)]
                if verbose: print("involution on seps", p)
                wsep = 0
                for cy in self.cylinders():
                    for k in cy[0]:
                        # check that p(k) is on the top of the cyl that has k on its bottom
                        if p[k] not in cy[1]:
                            test = False
                            break
                        # check that if k is on bot and top of cyl, then p(k) = k
                        if k in cy[1]:
                            if k != p[k]:
                                test = False
                                break
                            wsep += 1
                    if test is False: break
                if test is False: continue # try next j
                if verbose: print("wsep", wsep)
                # check number of w pts
                if wsep + 2*self.ncyls() != 2*z[0] + 4:
                    continue # try next j
                # check that cylinders are stable under involution
                if self != CylinderDiagram(
                    [(cy[0],tuple(map(lambda x: p[x],cy[0]))) for cy in self.cylinders()]):
                    continue # try next j
                return True
            return False

        else:
            return False

    #
    # construction
    #

    def dual_graph(self):
        r"""
        The dual graph of the stable curve at infinity in the horizontal
        direction.

        This graph is defines as follows. Cut each horizontal cylinder along a
        circumference, then the vertices are the equivalence class of half
        cylinder modulo the relation "linked by a saddle connection" and the
        edges are the circumferences.

        EXAMPLES::

            sage: from surface_dynamics import *

        We consider the three diagrams of the stratum H(1,1)::

            sage: c1 = CylinderDiagram('(0,1,2,3)-(0,1,2,3)')
            sage: c1.stratum()
            H_2(1^2)
            sage: c1.dual_graph()
            Looped multi-graph on 1 vertex
            sage: c2 = CylinderDiagram('(0,1)-(1,2) (2,3)-(0,3)')
            sage: c2.stratum()
            H_2(1^2)
            sage: c2.dual_graph()
            Looped multi-graph on 1 vertex
            sage: c3 = CylinderDiagram('(0,1)-(2,3) (2)-(0) (3)-(1)')
            sage: c3.stratum()
            H_2(1^2)
            sage: c3.dual_graph()
            Looped multi-graph on 2 vertices
        """
        from sage.graphs.graph import Graph
        cb = self.bot_cycle_tuples()
        ct = self.top_cycle_tuples()

        # first compute the equivalence class of half cylinders (i.e. gives vertices)
        V = Graph()
        V.add_vertices('%db' %c[0] for c in cb)
        V.add_vertices('%dt' %c[0] for c in ct)

        for i in range(self.nseps()):
            V.add_edge(
                ('%db' %self._bot_to_cyl[i][0]),
                ('%dt' %self._top_to_cyl[i][1]))

        # the dual graph
        G = Graph(loops=True,multiedges=True)
        cc = map(tuple,V.connected_components())
        hc2cc = {} # half-cyl to conn comp
        for c in cc:
            for e in c:
                hc2cc[e] = c
        for c in self.cylinders():
            G.add_edge(hc2cc['%db' %c[0][0]],hc2cc['%dt' %c[1][0]],(c[0][0],c[1][0]))
        return G

    def matrix_relation(self):
        r"""
        Return the matrix of relation on the lengths of the separatrices.

        The output matrix has size `ncyls \times nseps`.

        EXAMPLES::

            sage: from surface_dynamics import *

        For a one cylinder diagram, there is no relations::

            sage: cyl = CylinderDiagram('(0,1,2,3)-(0,1,2,3)')
            sage: cyl.matrix_relation()
            [0 0 0 0]

        Here is an example in the stratum H(2)::

            sage: cyl = CylinderDiagram('(0,1)-(0,2) (2)-(1)')
            sage: cyl.stratum()
            H_2(2)
            sage: cyl.matrix_relation()
            [ 0  1 -1]
            [ 0 -1  1]
        """
        from sage.matrix.constructor import matrix

        m = matrix(self.ncyls(),self.nseps(),sparse=True)
        for i,(top,bot) in enumerate(self.cylinders()):
            for t in top:
                m[i,t] = 1
            for b in bot:
                m[i,b] += -1
        return m

    #
    # Abelian differentials / coordinates
    #

    def stratum_component(self):
        r"""
        Return the connected component of stratum of ``self``.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: CylinderDiagram('(0,1)-(0,2) (2)-(1)').stratum_component()
            H_2(2)^hyp

            sage: c = CylinderDiagram('(0,3,2,1)-(1,4,3,2) (4,7,6,5)-(0,7,6,5)')
            sage: c.stratum_component()
            H_4(3^2)^hyp
            sage: c = CylinderDiagram('(0,1,4)-(1,6,7) (2,5,3)-(0,2,4) (6)-(5) (7)-(3)')
            sage: c.stratum_component()
            H_4(3^2)^nonhyp

            sage: c = CylinderDiagram('(0,6)-(1,7) (1,5,4,3,2)-(2,6,5,4,3) (7,9,8)-(0,9,8)')
            sage: c.stratum_component()
            H_5(4^2)^hyp
            sage: c = CylinderDiagram('(0,2,6,1)-(0,8,1,9,2,5,7,4) (3,7,4,8,9,5)-(3,6)')
            sage: c.stratum_component()
            H_5(4^2)^even
            sage: c = CylinderDiagram('(3,7,4,8,9,5)-(0,8,1,9,2,5,7,4) (0,2,6,1)-(3,6)')
            sage: c.stratum_component()
            H_5(4^2)^odd
        """
        stratum = self.stratum()
        cc = stratum._cc
        if len(cc) == 1:
            return cc[0](stratum)

        from .abelian_strata import HypASC
        if cc[0] is HypASC:
            if self.is_hyperelliptic():
                return HypASC(stratum)
            elif len(cc) == 2:
                return cc[1](stratum)

        if self.spin_parity() == 0:
            from .abelian_strata import EvenASC
            return EvenASC(stratum)
        else:
            from .abelian_strata import OddASC
            return OddASC(stratum)

    def smallest_integer_lengths(self):
        r"""
        Check if there is a integer solution that satisfy the cylinder
        conditions.

        If there is a solution, the function returns a list a pair
        ``(total_length, list_of_lengths)`` that consists of the sum of the
        length of the separatrices together with the list of lengths. Otherwise,
        returns False.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: c = CylinderDiagram('(0,1)-(0,2) (2,3)-(1,3)')
            sage: c.smallest_integer_lengths()
            (4, [1, 1, 1, 1])
            sage: c = CylinderDiagram('(0,1,2)-(3) (3)-(0) (4)-(1,2,4)')
            sage: c.smallest_integer_lengths()
            False

            sage: c = CylinderDiagram('(0,1)-(0,5) (2)-(3) (3,6)-(8) (4,8)-(6,7) (5)-(2,4) (7)-(1)')
            sage: c.smallest_integer_lengths()
            (13, [1, 2, 1, 1, 1, 2, 1, 2, 2])
        """
        if self.ncyls() == 1:
            return (self.nseps(), [1] * self.nseps())

        from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException

        n = self.nseps()
        bot = self.bot_cycle_tuples()
        top = [self.top_orbit(self._bot_to_cyl[b[0]][1]) for b in bot]

        p = MixedIntegerLinearProgram(maximization=False)
        scl = p.new_variable(nonnegative=True)
        p.set_objective(sum(scl[i] for i in range(n)))
        for i in range(n):
            p.add_constraint(scl[i],min=1)
        for b,t in zip(bot,top):
            p.add_constraint(
                    p.sum(scl[i] for i in set(b).difference(t)) ==
                    p.sum(scl[i] for i in set(t).difference(b))
                    )

        try:
            total = Integer(p.solve())
            lengths = [Integer(p.get_values(scl[i])) for i in range(n)]
            return total, lengths
        except MIPSolverException:
            return False

    #
    # homology
    #

    def to_ribbon_graph(self):
        r"""
        Return a ribbon graph

        A *ribbon graph* is a graph embedded in an oriented surface such that
        its complement is a union of topological discs. To a cylinder diagram we
        associate the graph which consists of separatrices together with a
        choice of one vertical edge in each cylinder.

        The edges of the ribbon graph are labeled by ``(i,nseps+i)`` for
        separatrices and by ``(2(nseps+j),2(nseps+j)+1)`` for vertical in
        cylinders.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: C = CylinderDiagram([((0,1),(0,2)),((2,),(1,))])
            sage: C.stratum()
            H_2(2)
            sage: R = C.to_ribbon_graph(); R
            Ribbon graph with 1 vertex, 5 edges and 2 faces
            sage: l,m = R.cycle_basis(intersection=True)
            sage: m.rank() == 2 * C.genus()
            True

        TESTS::

            sage: f = lambda c: c.to_ribbon_graph().cycle_basis(intersection=True)[1]

            sage: a = AbelianStratum(2)
            sage: all(f(c).rank() == 4 for c in a.cylinder_diagrams())
            True
            sage: a = AbelianStratum(1,1)
            sage: all(f(c).rank() == 4 for c in a.cylinder_diagrams())
            True
        """
        from .homology import RibbonGraphWithAngles

        n = self.nseps()
        m = self.ncyls()

        edges = [(i,n+i) for i in range(n)] + [(2*(n+i),2*(n+i)+1) for i in range(m)]
        faces = []
        angles = [1] * (2*(n+m))
        half = Integer(1)/Integer(2)
        for j,(b,t) in enumerate(self.cylinders()):
            face = [i for i in b] + [2*(n+j)] + [n+i for i in t[1:]+t[:1]] + [2*(n+j)+1]
            faces.append(tuple(face))
            t1 = t[1] if len(t) > 1 else t[0]
            angles[b[0]] = angles[n+t1] = angles[2*(n+j)] = angles[2*(n+j)+1] = half

        return RibbonGraphWithAngles(edges=edges,faces=faces,angles=angles)

    def to_ribbon_graph_with_holonomies(self, lengths, heights, twists):
        from .homology import RibbonGraphWithHolonomies

        n = self.nseps()
        m = self.ncyls()

        edges = [(i,n+i) for i in range(n)] + [(2*(n+i),2*(n+i)+1) for i in range(m)]
        faces = []

        for j,(b,t) in enumerate(self.cylinders()):
            face = [i for i in b] + [2*(n+j)] + [n+i for i in t[1:]+t[:1]] + [2*(n+j)+1]
            faces.append(tuple(face))

        holonomies = [None] * (2*(n+m))
        for i in range(n):
            holonomies[i]   = (lengths[i],0)
            holonomies[n+i] = (-lengths[i],0)
        for i in range(m):
            holonomies[2*(n+i)]   = (twists[i], heights[i])
            holonomies[2*(n+i)+1] = (-twists[i], -heights[i])

        return RibbonGraphWithHolonomies(edges=edges,faces=faces,holonomies=holonomies)



    def spin_parity(self):
        r"""
        Return the spin parity of any surface that is built from this cylinder
        diagram.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: c = CylinderDiagram('(0,1,2,3,4)-(0,1,2,3,4)')
            sage: c.spin_parity()
            0
            sage: c = CylinderDiagram('(0,1,2,3,4)-(0,1,4,2,3)')
            sage: c.spin_parity()
            1

            sage: c = CylinderDiagram('(0,2,6,1)-(0,8,1,9,2,5,7,4) (3,7,4,8,9,5)-(3,6)')
            sage: c.spin_parity()
            0
            sage: c = CylinderDiagram('(3,7,4,8,9,5)-(0,8,1,9,2,5,7,4) (0,2,6,1)-(3,6)')
            sage: c.spin_parity()
            1
        """
        if any(z%2 for z in self.stratum().zeros()):
            return None
        return self.to_ribbon_graph().spin_parity()

#    def circumferences_of_cylinders(self,ring=None):
#        r"""
#        Return the set of circumferences of cylinders as cycles in the chain
#        space.
#        """
#        from sage.modules.free_module import FreeModule
#        from copy import copy
#
#        if ring is None:
#            from sage.rings.integer_ring import ZZ
#            ring = ZZ
#
#        g = self.to_ribbon_graph()
#        C = g.chain_complex(ring)
#        C1 = C.chain_space(1)
#        Z1 = C.cycle_space(1)
#        n = g.num_edges()
#
#        V = FreeModule(ring, n)
#
#        l = []
#        for (b,t) in self.cylinders():
#            v = copy(V.zero())
#            for i in b:
#                v[g.dart_to_edge(i,orientation=True)] = 1
#            l.append(Z1(V(v)))
#        return l

    #
    # build one or many origamis
    #

    def an_origami(self):
        r"""
        Return one origami with this diagram cylinder if any.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: cyl = CylinderDiagram('(0,1)-(0,2) (2,3)-(1,3)')
            sage: cyl.an_origami()
            (1,2)(3,4)
            (1,3,4,2)
        """
        res = self.smallest_integer_lengths()
        if res is False:
            return False
        m,lengths = res

        widths = [sum(lengths[i] for i in bot) for bot in self.bot_cycle_tuples()]
        areas = [widths[i] for i in range(self.ncyls())]

        v = [0]
        for a in areas:
            v.append(v[-1] + a)

        # initialization of bottom squares: sep_i -> bottom position
        sep_bottom_pos = [None] * self.nseps()
        for i,(bot,_) in enumerate(self.cylinders()):
            w = 0
            for j in bot:
                sep_bottom_pos[j] = v[i] + w
                w += lengths[j]

        # initialization of sigma_h which remains constant
        lx = list(range(1, v[-1]+1))
        for i in range(self.ncyls()):
            for j in range(v[i], v[i+1], widths[i]):
                lx[j+widths[i]-1] = j

        # initialization of y except the top
        ly = []
        for i in range(self.ncyls()):
            ly.extend([None]*widths[i])

        # build the top interval without twist
        for i,(_,top_seps) in enumerate(self.cylinders()):
            top = []
            for k in reversed(top_seps):
                top.extend(range(sep_bottom_pos[k],sep_bottom_pos[k]+lengths[k]))
            ly[v[i+1]-widths[i]:v[i+1]] = top

        # yield the origami without twist
        from surface_dynamics.flat_surfaces.origamis.origami import Origami_dense_pyx
        return Origami_dense_pyx(tuple(lx), tuple(ly))

    def origami_iterator(self,n):
        r"""
        Iteration over all origamis with n squares.

        INPUT:

        - ``n`` - positive integer - the number of squares

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: cyl = CylinderDiagram('(0,1,2)-(3,1,2) (3)-(0)')
            sage: for o in cyl.origami_iterator(4):
            ....:     print(o)
            ....:     print(o.stratum())
            ....:     print(o.nb_squares())
            (1,2,3)(4)
            (1,4)(2,3)
            H_2(1^2)
            4
            (1,2,3)(4)
            (1,2,4)(3)
            H_2(1^2)
            4
            (1,2,3)(4)
            (1,3,4)(2)
            H_2(1^2)
            4
        """
        for w, h in self.widths_and_heights_iterator(n):
            yield from self.cylcoord_to_origami_iterator(w, h)

    def origamis(self,n=None):
        r"""
        Return the set of origamis having ``n`` squares.

        If ``n`` is None then return the origamis with less number of squares.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: cyl = CylinderDiagram('(0,1,2)-(0,1,3) (3)-(2)')
            sage: o5 = cyl.origamis(5)
            sage: o5[0]
            (1,2,3,4)(5)
            (1,5,4,2,3)
            sage: o5[1].nb_squares()
            5
            sage: o5[2].stratum_component()
            H_2(1^2)^hyp

        TESTS::

            sage: c1 = CylinderDiagram("(0,1)-(0,3,4,5) (2,3,5)-(1) (4)-(2)")
            sage: c2 = CylinderDiagram("(0,1)-(0,5) (2)-(4) (3,4)-(1) (5)-(2,3)")
            sage: c3 = CylinderDiagram("(0,3)-(5) (1)-(0) (2,5)-(3,4) (4)-(1,2)")
            sage: len(c1.origamis(8))
            12
            sage: len(c2.origamis(8))
            12
            sage: len(c3.origamis(8))
            12
        """
        if n is None:
            res = self.smallest_integer_lengths()
            if res is False:
                return False
            n = res[0]

        return list(self.origami_iterator(n))

    def widths_and_heights_iterator(self, n, height_one=False):
        """
        Iterate over the possible integer widths and heights of the cylinders
        for which the corresponding translation surface has area ``n``.

        At each iteration, the output is a pair of ``(lengths, heights)``. You
        can then use :meth:`cylcoord_to_origami` to build the corresponding
        origami.

        INPUT:

        - ``height_one`` -- (boolean default ``False``) whether to return only
          coordinates with cylinder height one.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: cyl = CylinderDiagram([((0,1),(0,2)),((2,),(1,))])
            sage: cyl
            (0,1)-(0,2) (2)-(1)

            sage: it = cyl.widths_and_heights_iterator(10)
            sage: l,h = next(it)
            sage: print(l)
            (2, 1, 1)
            sage: print(h)
            [3, 1]
            sage: cyl.cylcoord_to_origami(l,h)
            (1,2,3)(4,5,6)(7,8,9)(10)
            (1,4,7)(2,5,8)(3,6,9,10)

            sage: it = cyl.widths_and_heights_iterator(10, height_one=True)
            sage: l,h = next(it)
            sage: print(l)
            (8, 1, 1)
            sage: print(h)
            [1, 1]
            sage: cyl.cylcoord_to_origami(l,h)
            (1,2,3,4,5,6,7,8,9)(10)
            (1)(2)(3)(4)(5)(6)(7)(8)(9,10)

        TESTS::

            sage: c1 = CylinderDiagram("(0,1)-(0,3,4,5) (2,3,5)-(1) (4)-(2)")
            sage: c2 = CylinderDiagram("(0,1)-(0,5) (2)-(4) (3,4)-(1) (5)-(2,3)")
            sage: c3 = CylinderDiagram("(0,3)-(5) (1)-(0) (2,5)-(3,4) (4)-(1,2)")
            sage: list(c1.widths_and_heights_iterator(8))
            [((1, 3, 1, 1, 1, 1), [1, 1, 1])]
            sage: list(c2.widths_and_heights_iterator(8))
            [((1, 2, 1, 1, 1, 2), [1, 1, 1, 1])]
            sage: list(c3.widths_and_heights_iterator(8))
            [((1, 1, 1, 1, 2, 2), [1, 1, 1, 1])]

            sage: for n in range(8, 12):
            ....:     for c in c1,c2,c3:
            ....:         L1 = [wh for wh in c.widths_and_heights_iterator(n) if all(h == 1 for h in wh[1])]
            ....:         L2 = list(c.widths_and_heights_iterator(n, height_one=True))
            ....:         assert L1 == L2, (L1, L2)
        """
        from sage.combinat.integer_lists import IntegerListsLex
        from sage.rings.integer_ring import ZZ
        from sage.modules.free_module import FreeModule
        from copy import copy

        V = FreeModule(ZZ,self.nseps())

        m = self.matrix_relation()

        min_lengths = [1] * self.nseps()
        for i in range(self.ncyls()):
            pos = m.nonzero_positions_in_row(i)
            pos_m = list(filter(lambda j: m[i,j] == -1, pos))
            pos_p = list(filter(lambda j: m[i,j] == 1, pos))
            if len(pos_m) == 1:
                min_lengths[pos_m[0]] = max(min_lengths[pos_m[0]], len(pos_p))
            if len(pos_p) == 1:
                min_lengths[pos_p[0]] = max(min_lengths[pos_m[0]], len(pos_m))

        min_widths = []
        for bot,top in self.cylinders():
            min_widths.append(max(
                sum(min_lengths[j] for j in top),
                sum(min_lengths[j] for j in bot)))

        for a in filter(
              lambda x: all(x[i] >= min_widths[i] for i in range(self.ncyls())),
              IntegerListsLex(n=n, length=self.ncyls(), min_part=1)):
            if height_one:
                area_div = [[x] for x in a]
            else:
                area_div = tuple(list(filter(lambda d: d >= min_widths[i],arith.divisors(a[i]))) for i in range(self.ncyls()))
            for w in itertools.product(*area_div):
                h = [Integer(a[i]/w[i]) for i in range(self.ncyls())]

                # from here the resolution becomes linear and convex ...
                #TODO: program a linear and convex solution
                seps_b = [c[0] for c in self.cylinders()]
                nseps_b = list(map(len, seps_b))
                lengths = tuple(IntegerListsLex(n=w[i], length=nseps_b[i], min_part=1) for i in range(self.ncyls()))
                for l_by_cyl in itertools.product(*lengths):
                    l = copy(V.zero())
                    for i in range(self.ncyls()):
                        for j in range(nseps_b[i]):
                            l[seps_b[i][j]] = l_by_cyl[i][j]
                    if not m*l:
                        yield l,h

    def cylcoord_to_origami_iterator(self, lengths, heights):
        r"""
        Convert coordinates of the cylinders into an origami.

        INPUT:

        - ``lengths`` - lengths of the separatrices

        - ``heights`` - heights of the cylinders

        OUTPUT:

        - iterator over all possible origamis with those lengths and heights...

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: cyl = CylinderDiagram('(0,1,2)-(3,1,2) (3)-(0)')
            sage: for o in cyl.cylcoord_to_origami_iterator((1,1,1,1),(1,1)):
            ....:     print(o)
            (1,2,3)(4)
            (1,4)(2,3)
            (1,2,3)(4)
            (1,2,4)(3)
            (1,2,3)(4)
            (1,3,4)(2)

        The number of origamis generated is just the product of the widths::

            sage: sum(1 for _ in cyl.cylcoord_to_origami_iterator((2,1,1,2),(3,2)))
            8
        """
        from surface_dynamics.flat_surfaces.origamis.origami_dense import Origami_dense_pyx
        from sage.combinat.gray_codes import product

        widths = [sum(lengths[i] for i in bot) for bot in self.bot_cycle_tuples()]
        areas = [heights[i]*widths[i] for i in range(self.ncyls())]

        # initialization of partial volumes: the set of squares in cylinder i is range(v[i],v[i+1])
        v = [0]
        for a in areas:
            v.append(v[-1] + a)

        # initialization of bottom squares: sep_i -> bottom position
        sep_bottom_pos = [None] * self.nseps()
        for i,(bot,_) in enumerate(self.cylinders()):
            w = 0
            for j in bot:
                sep_bottom_pos[j] = v[i] + w
                w += lengths[j]

        # initialization of sigma_h which remains constant
        lx = list(range(1, v[-1]+1))
        for i in range(self.ncyls()):
            for j in range(v[i], v[i+1], widths[i]):
                lx[j+widths[i]-1] = j

        # initialization of y except the top
        ly = []
        for i in range(self.ncyls()):
            ly.extend(range(v[i]+widths[i],v[i+1]))
            ly.extend([None]*widths[i])

        # build the top interval without twist
        for i,(_,top_seps) in enumerate(self.cylinders()):
            top = []
            for k in reversed(top_seps):
                top.extend(range(sep_bottom_pos[k],sep_bottom_pos[k]+lengths[k]))
            ly[v[i+1]-widths[i]:v[i+1]] = top

        # yield the one without twist
        yield Origami_dense_pyx(tuple(lx),tuple(ly))

        # yield the others using a Gray code
        for i,o in product(widths):
            if o == 1:
                ly.insert(v[i+1]-widths[i],ly.pop(v[i+1]-1))
            else:
                ly.insert(v[i+1]-1,ly.pop(v[i+1]-widths[i]))
            yield Origami_dense_pyx(tuple(lx),tuple(ly))

    def cylcoord_to_origami(self, lengths, heights, twists=None):
        r"""
        Convert coordinates of the cylinders into an origami.

        INPUT:

        - ``lengths`` - lengths of the separatrices

        - ``heights`` - heights of the cylinders

        - ``twists`` - twists for cylinders

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: c = CylinderDiagram([((0,1),(1,2)),((2,),(0,))])
            sage: c.stratum()
            H_2(2)
            sage: c.cylcoord_to_origami([1,1,1],[1,1]).stratum()
            H_2(2)
            sage: o1 = c.cylcoord_to_origami([2,1,2],[1,1],[1,0])
            sage: o1 = o1.relabel()
            sage: o2 = c.cylcoord_to_origami([2,1,2],[1,1],[0,1])
            sage: o2 = o2.relabel()
            sage: o3 = c.cylcoord_to_origami([2,1,2],[1,1],[1,1])
            sage: o3 = o3.relabel()
            sage: all(o.stratum() == AbelianStratum(2) for o in [o1,o2,o3])
            True
            sage: o1 == o2 or o1 == o3 or o3 == o1
            False

        If the lengths are not compatible with the cylinder diagram a ValueError
        is raised::

            sage: c.cylcoord_to_origami([1,2,3],[1,1])
            Traceback (most recent call last):
            ...
            ValueError: lengths are not compatible with cylinder equations

        TESTS::

            sage: c = CylinderDiagram([((0,),(1,)), ((1,2,3),(0,2,3))])
            sage: c
            (0)-(1) (1,2,3)-(0,2,3)
            sage: lengths = [1,1,1,1]
            sage: heights = [1,1]
            sage: c.cylcoord_to_origami(lengths,heights,[0,0])
            (1)(2,3,4)
            (1,2)(3,4)
            sage: c.cylcoord_to_origami(lengths,heights,[0,1])
            (1)(2,3,4)
            (1,2,3)(4)
            sage: c.cylcoord_to_origami(lengths,heights,[0,2])
            (1)(2,3,4)
            (1,2,4)(3)
        """
        from sage.rings.integer_ring import ZZ
        from surface_dynamics.flat_surfaces.origamis.origami_dense import Origami_dense_pyx

        if len(lengths) != self.nseps():
            raise ValueError("the 'lengths' vector has wrong length")
        else:
            lengths = [ZZ.coerce(l) for l in lengths]

        if len(heights) != self.ncyls():
            raise ValueError("the 'heights' vector has wrong length")
        else:
            heights = [ZZ.coerce(h) for h in heights]

        widths = [sum(lengths[i] for i in bot) for bot,_ in self.cylinders()]

        if widths != [sum(lengths[i] for i in top) for _,top in self.cylinders()]:
            raise ValueError("lengths are not compatible with cylinder equations")

        if twists is None:
            twists = [0] * len(widths)
        elif len(twists) != self.ncyls():
            raise ValueError("the 'twists' vector has wrong length")
        else:
            twists = [(-twists[i])%widths[i] for i in range(len(widths))]

        areas = [heights[i]*widths[i] for i in range(self.ncyls())]

        # initialization of partial volumes: the set of squares in cylinder i is range(v[i],v[i+1])
        v = [0]
        for a in areas:
            v.append(v[-1] + a)

        # initialization of bottom squares: sep_i -> bottom position
        sep_bottom_pos = [None] * self.nseps()
        for i,(bot,_) in enumerate(self.cylinders()):
            w = 0
            for j in bot:
                sep_bottom_pos[j] = v[i] + w
                w += lengths[j]

        # build the permutation r
        lx = list(range(1, v[-1]+1))
        for i in range(self.ncyls()):
            for j in range(v[i], v[i+1], widths[i]):
                lx[j+widths[i]-1] = j

        # build permutation u with the given twists
        ly = []
        for i,(_,top_seps) in enumerate(self.cylinders()):
            # everything excepted the top
            ly.extend(range(v[i]+widths[i],v[i+1]))

            # the top
            k = top_seps[0]
            top = list(range(sep_bottom_pos[k],sep_bottom_pos[k]+lengths[k]))
            for k in reversed(top_seps[1:]):
                top.extend(range(sep_bottom_pos[k],sep_bottom_pos[k]+lengths[k]))
            ly.extend(top[twists[i]:] + top[:twists[i]])

        # yield the one without twist
        return Origami_dense_pyx(tuple(lx), tuple(ly))

    def cylinder_graph(self):
        """
        Return the cylinder graph.

        The cylinder graph is the graph whose vertex set are the cylinders
        and for each saddle connection there is a directed edge from the
        adjacent cylinders. The multiplicities are encoded in the labels.

        EXAMPLES::

            sage: from surface_dynamics import CylinderDiagram

            sage: c = CylinderDiagram('(0,1,5)-(2,5) (2)-(0,1,3) (3,4)-(4)')
            sage: c.cylinder_graph()
            Looped digraph on 3 vertices
            sage: c.cylinder_graph().edges(sort=True)
            [(0, 0, 1), (0, 1, 1), (1, 0, 2), (1, 2, 1), (2, 2, 1)]

            sage: c = CylinderDiagram('(0,1,3,5)-(2,5,3) (2,4)-(0,4,1)')
            sage: c.cylinder_graph().edges(sort=True)
            [(0, 0, 2), (0, 1, 1), (1, 0, 2), (1, 1, 1)]
        """
        bot_to_cyl = [None] * self.degree()
        top_to_cyl = [None] * self.degree()
        for i, (bot, top) in enumerate(self.cylinders()):
            for saddle in bot:
                bot_to_cyl[saddle] = i
            for saddle in top:
                top_to_cyl[saddle] = i

        edges = collections.defaultdict(int)
        for u, v in zip(top_to_cyl, bot_to_cyl):
            edges[u, v] += 1

        G = DiGraph(self.ncyls(), loops=True, multiedges=False, weighted=True)
        G.add_edges([(u, v, m) for ((u, v), m) in edges.items()])
        return G

    #TODO
#    def chain_complex_dual(self, ring=None):
#        r"""
#        Return a chain complex for the cylinder diagram
#
#        The vertices are in bijection with the cylinder of self
#        The edges are in bijection with separatrices and cylinders
#
#        """
#        from homology import TranslationSurfaceChainComplex
#        from sage.rings.integer import Integer
#
#        if ring is None:
#            from sage.rings.integer_ring import IntegerRing
#            ring = IntegerRing()
#
#        vertices = []   # list of list of vertex = integers from 0 to ncyls
#                        # (in or out, edge)
#        edges = {}      # label -> (start vertex,end)
#        angles = {}     # (in or out,label) -> angle to next edge
#
#        cyls = self.cylinders()
#        t2c = [None]*self.nseps()
#        b2c = [None]*self.nseps()
#
#        for k,(cb,ct) in enumerate(cyls):
#            for i in cb: b2c[i] = k
#            for i in ct: t2c[i] = k
#
#        for k,(cb,ct) in enumerate(cyls):
#            vertex = []
#            e = 'c%d' %k
#            edges[e] = (k,k)
#            angles[(1,e)] = Integer(1)/Integer(2)
#            angles[(-1,e)] = Integer(1)/Integer(2)
#            # the incoming edges from bottom
#            for i in cb:
#                e = 's%d' %i
#                edges[e] = (t2c[i],k)
#                vertex.append((-1,e))
#                angles[(-1,e)] = Integer(0)
#            angles[(-1,e)] = Integer(1)/Integer(2)
#
#            # the central edge (outgoing)
#            vertex.append((1, 'c%d' %k))
#
#            # the outgoing edges from top
#            for i in ct:
#                e = 's%d' %i
#                vertex.append((1,e))
#                angles[(1,e)] = Integer(0)
#            angles[(1,e)] = Integer(1)/Integer(2)
#
#            # the central edge (incoming)
#            vertex.append((-1, 'c%d' %k))
#
#            vertices.append(vertex)
#
#        return TranslationSurfaceChainComplex(ring,vertices,edges,angles)


def move_forward(i, v, g01, g23):
    r"""
    Helper function to build pillowcase covers

    EXAMPLES::

        sage: from surface_dynamics.flat_surfaces.separatrix_diagram import move_forward
        sage: g23 = [1,2,0]
        sage: g01 = [2,0,1]
        sage: i,v = 0,0
        sage: for _ in range(6):
        ....:     i,v = move_forward(i, v, g01, g23)
        ....:     print("%d %d" % (i,v))
        0 1
        1 0
        1 1
        2 0
        2 1
        0 0

        sage: i,v = 0,2
        sage: for _ in range(6):
        ....:     i,v = move_forward(i, v, g01, g23)
        ....:     print("%d %d" % (i,v))
        0 3
        2 2
        2 3
        1 2
        1 3
        0 2
    """
    if v == 0:
        v = 1
    elif v == 2:
        v = 3
    elif v == 1:
        v = 0
        i = g23[i]
    elif v == 3:
        v = 2
        i = g01[i]
    else:
        raise ValueError

    return i,v

def move_backward(i, v, g01, g23):
    r"""
    Helper function to build pillowcase covers

    EXAMPLES::

        sage: from surface_dynamics.flat_surfaces.separatrix_diagram import move_backward
        sage: g23 = [1,2,0]
        sage: g01 = [2,0,1]
        sage: i,v = 0,0
        sage: for _ in range(6):
        ....:     i,v = move_backward(i, v, g01, g23)
        ....:     print("%d %d" % (i,v))
        2 1
        2 0
        1 1
        1 0
        0 1
        0 0

        sage: i,v = 0,2
        sage: for _ in range(6):
        ....:     i,v = move_backward(i, v, g01, g23)
        ....:     print("%d %d" % (i,v))
        1 3
        1 2
        2 3
        2 2
        0 3
        0 2
    """
    if v == 0:
        v = 1
        i = g01[i]
    elif v == 1:
        v = 0
    elif v == 2:
        v = 3
        i = g23[i]
    elif v == 3:
        v = 2
    else:
        raise ValueError

    return i,v

def simplex_count(rays):
    r"""
    EXAMPLES::

        sage: from surface_dynamics.flat_surfaces.separatrix_diagram import simplex_count

        sage: rays = [(0,1,1), (1,0,1), (1,1,0)]
        sage: simplex_count(rays)
        1

        sage: rays = [(0,1,1), (1,0,1), (1,1,0)]
        sage: simplex_count(rays)
        1

        sage: rays = [(0,1,1,1),(1,0,1,1),(1,1,0,1),(1,1,1,0)]
        sage: simplex_count(rays)
        1
    """
    from sage.geometry.polyhedron.constructor import Polyhedron
    d = len(rays[0])
    return Polyhedron([[0]*d] + list(rays)).integral_points_count() - len(rays)

class QuadraticCylinderDiagram(SageObject):
    r"""
    Cylinder diagram for quadratic differentials

    Cylinder diagrams are encoded as a Ribbon graph together with a pairing of
    faces (in particular the number of faces must be even).

    EXAMPLES::

        sage: from surface_dynamics import *

    If you start with strings, the cylinders are preserved but the names of
    saddle connections are changed::

        sage: QuadraticCylinderDiagram('(4,4,5)-(6,6,1) (2,3,2,0)-(1,0,5,3)')
        (0,0,1)-(2,2,3) (4,5,4,6)-(3,6,1,5)
    """
    def __init__(self, arg1, arg2=None):
        r"""
        INPUT: there are two input formats

        - with a unique argument

        - with two arguments

        - ``g`` -- a ribbon graph

        - ``pairing`` -- the pairing of faces of g (= permutation)
        """
        from .homology import RibbonGraph

        if arg2 is None:
            if isinstance(arg1, str):
                data = [(string_to_cycle(b),string_to_cycle(t)) for b,t in (w.split('-') for w in arg1.split(' '))]
            else:
                data = arg1
            N = 0
            for i,pair in enumerate(data):
                if not isinstance(pair, (tuple, list)) or len(pair) != 2:
                    raise TypeError('input must be a list of pairs')
                data[i] = [[int(x) for x in pair[0]], [int(x) for x in pair[1]]]
                N += len(pair[0]) + len(pair[1])

            if N%2:
                raise ValueError('each symbol must appear exactly twice')
            n = N//2
            for b,t in data:
                if any(i < 0 or i >= n for i in b) or \
                   any(i < 0 or i >= n for i in t):
                        raise ValueError('symbol out of range')

            k = 0
            faces = []
            edges = [None] * N
            seen = [None] * n
            p = []
            for i, bt in enumerate(data):
                # constructing p (= face matching)
                p.append(2*i+1)
                p.append(2*i)

                # constructing edges and faces
                for f in bt:
                    faces.append(tuple(range(k, k+len(f))))

                    for j in range(len(f)):
                        e = f[j]
                        if seen[e] is not None:
                            if seen[e] == -1:
                                raise ValueError('number %d appears more than twice' % e)
                            kk = seen[e]
                            edges[kk] = k + j
                            edges[k + j] = kk
                            seen[e] = -1
                        else:
                            seen[e] = k + j

                    k += len(f)

            g = RibbonGraph(edges=edges, faces=faces, connected=False)

        else:
            g = arg1
            p = arg2
            if not isinstance(g, RibbonGraph):
                raise ValueError

            p = perm_init(p)

        self._g = g
        self._p = p

        self._check()

        k = 0
        fc = self._face_to_cylinder_index = [None] * self._g.num_faces()
        for i in range(len(self._p)):
            if fc[i] is None:
                fc[i] = k
                fc[self._p[i]] = k
                k += 1

    def _check(self):
        self._g._check()

        f = self._g.num_faces()
        if f % 2:
            raise ValueError('the number of faces of the fatgraph must be even')
        if len(self._p) != f:
            raise ValueError('the pairing has wrong length')
        for i,j in enumerate(self._p):
            if i == j or self._p[j] != i:
                raise ValueError('the pairing is not an involution')

        from sage.graphs.graph import Graph
        G = Graph(loops=True, multiedges=True)
        for i in range(self._g._total_darts):
            if self._g._active_darts[i]:
                G.add_edge(i,self._g._vertices[i])
                G.add_edge(i,self._g._edges[i])
                G.add_edge(i,self._g._faces[i])
        for i,j in enumerate(self._p):
            k1 = self._g._face_cycles[i][0]
            k2 = self._g._face_cycles[j][0]
            G.add_edge(k1, k2)

        if not G.is_connected():
            raise ValueError("the graph is not connected")

    def num_darts(self):
        r"""
        Number of darts
        """
        return self._g.num_darts()

    def num_edges(self):
        r"""
        Number of edges
        """
        return self._g.num_edges()

    nseps = num_edges

    def edges(self):
        r"""
        The set of edges.
        """
        return self._g.edges()

    def num_cylinders(self):
        r"""
        Number of cylinders.
        """
        return len(self._p) // 2

    ncyls = num_cylinders

    def stratum(self):
        r"""
        Return the stratum of quadratic differentials associated to this cylinder diagram

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: QuadraticCylinderDiagram('(0,0)-(1,1)').stratum()
            Q_0(-1^4)

            sage: QuadraticCylinderDiagram('(0,0)-(1,1,2,2,3,3)').stratum()
            Q_0(1, -1^5)

            sage: QuadraticCylinderDiagram('(0,2,3,2)-(1,0,1,3)').stratum()
            Q_2(2^2)

            sage: QuadraticCylinderDiagram('(0,1)-(2,3) (0)-(4,4) (1)-(5,5) (2)-(6,6) (3)-(7,7)').stratum()
            Q_0(2^2, -1^8)
        """
        from surface_dynamics.flat_surfaces.quadratic_strata import QuadraticStratum
        return QuadraticStratum([len(t)-2 for t in self._g._vertex_cycles])

    def cylinders(self, dart=False):
        r"""
        Cylinders of self

        Return a list of pairs ``(bot, top)`` where, by convention the bottom
        corresponds to the face with smaller index in the list of faces.

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: from surface_dynamics.flat_surfaces.separatrix_diagram import QuadraticCylinderDiagram
            sage: rg = RibbonGraph(edges='(0,1)(2,3)(4,5)(6,7)', faces='(0,1)(2,4,5)(3)(6,7)', connected=False)
            sage: q = QuadraticCylinderDiagram(rg, '(0,1)(2,3)')
            sage: q.cylinders()
            [((0, 0), (1, 2, 2)), ((1,), (3, 3))]
            sage: q.cylinders(True)
            [((0, 1), (2, 4, 5)), ((3,), (6, 7))]
        """
        ans = []
        g = self._g
        for i,j in enumerate(self._p):
            if i > j:
                continue
            bot = g._face_cycles[i]
            top = g._face_cycles[j]
            if dart:
                bot = tuple(bot)
                top = tuple(top)
            else:
                bot = tuple(g._dart_to_edge_index[x] for x in bot)
                top = tuple(g._dart_to_edge_index[x] for x in top)
            ans.append((bot,top))

        return ans

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics import *
            sage: from surface_dynamics.flat_surfaces.separatrix_diagram import QuadraticCylinderDiagram
            sage: rg = RibbonGraph(edges='(0,1)(2,3)(4,5)(6,7)', faces='(0,1)(2,4,5)(3)(6,7)', connected=False)
            sage: q = QuadraticCylinderDiagram(rg, '(0,1)(2,3)')
            sage: q
            (0,0)-(1,2,2) (1)-(3,3)
        """
        s = []
        for bot,top in self.cylinders():
            bot = ','.join(str(x) for x in bot)
            top = ','.join(str(x) for x in top)
            s.append('({})-({})'.format(bot,top))
        return ' '.join(s)

    def lengths_cone(self):
        r"""
        Return the polytope of admissible lengths.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: rg = RibbonGraph(edges='(0,1)(2,3)(4,5)(6,7)', faces='(0,1)(2,4,5)(3)(6,7)', connected=False)
            sage: q = QuadraticCylinderDiagram(rg, '(0,1)(2,3)')
            sage: q.stratum()
            Q_0(1, -1^5)
            sage: q.lengths_cone().rays_list()
            [[1, 0, 1, 0], [1, 2, 0, 1]]

            sage: rg = RibbonGraph(edges='(0,1)(2,3)(4,5)(6,7)', faces='(0,2,3)(1)(4,6,7)(5)', connected=False)
            sage: q = QuadraticCylinderDiagram(rg, '(0,2)(1,3)')

            sage: rg = RibbonGraph(edges='(1,2)(3,4)(5,6)(7,8)(9,10)(11,12)', faces='(1,2)(3,4,6)(5)(8)(7,9,10)(11,12)', connected=False)
            sage: q = QuadraticCylinderDiagram(rg, '(0,1)(2,4)(3,5)')
            sage: q.stratum()
            Q_0(1^2, -1^6)
            sage: L = q.lengths_cone()
            sage: L
            A 3-dimensional polyhedron in QQ^6 defined as the convex hull of 1 vertex and 3 rays
        """
        from sage.rings.integer_ring import ZZ

        n = self.num_edges()

        ieqs = []
        eqns = []

        # lengths are non-negative
        e = [ZZ.zero()] * (n+1)
        for i in range(n):
            e[i+1] = ZZ.one()
            ieqs.append(e[:])
            e[i+1] = ZZ.zero()

        # for each cylinder, length top = length bot
        for i,(bot,top) in enumerate(self.cylinders()):
            e = [ZZ.zero()] * (n+1)
            for i in bot:
                e[i+1] += ZZ.one()
            for i in top:
                e[i+1] -= ZZ.one()
            eqns.append(e)

        from sage.geometry.polyhedron.constructor import Polyhedron
        return Polyhedron(eqns=eqns, ieqs=ieqs)

    def widths_generating_series(self, var='w'):
        r"""
        Generating series of the number of saddle connection lengths of this quadratic
        differential cylinder diagram.

        .. WARNING::

            When a triangulation is involved, the generating series ignore
            some lower dimensional polytopes that are counted twice!

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: q = QuadraticCylinderDiagram('(0,1,2,3,3)-(0,4,4,2,1)')
            sage: q.widths_generating_series()  # optional -- latte_int
            (1)/((1 - w)^3*(1 - w^2))

            sage: q = QuadraticCylinderDiagram('(0,0,1,1,2,2)-(3,3,4,4)')
            sage: q.widths_generating_series()  # optional -- latte_int
            (3)/((1 - w^2)^4)

            sage: q = QuadraticCylinderDiagram('(0,0,1)-(2,2,3) (1,4)-(3,4)')
            sage: q.widths_generating_series()  # optional -- latte_int
            (1)/((1 - w1)*(1 - w0*w1)*(1 - w0^2))

            sage: q = QuadraticCylinderDiagram('(0,0,1,2,3)-(1,4,4,5,6) (2,5,7,7,8)-(3,6,8,9,9)')
            sage: F = q.widths_generating_series()  # optional -- latte_int
            sage: q = QuadraticCylinderDiagram('(0,0,1,2,3)-(1,4,4,5,6) (2,5,7,8,8)-(3,6,7,9,9)')
            sage: F = q.widths_generating_series()  # optional -- latte_int
        """
        from sage.matrix.constructor import matrix
        from surface_dynamics.misc.multiplicative_multivariate_generating_series import MultiplicativeMultivariateGeneratingSeriesRing
        from sage.geometry.polyhedron.constructor import Polyhedron

        sub1 = [[] for i in range(self.nseps())]
        sub2 = [[] for i in range(self.nseps())]
        for i,(top,bot) in enumerate(self.cylinders()):
            for j in top:
                sub1[j].append(i)
            for j in bot:
                sub2[j].append(i)

        dim = self.nseps() - self.ncyls()
        M = MultiplicativeMultivariateGeneratingSeriesRing(self.ncyls(), 'w')
        V = M.free_module()
        ans = M.zero()
        for rays in cone_triangulate(self.lengths_cone()):
            assert len(rays) == dim, (self, rays, dim)
            num = simplex_count(rays)

            den1 = {}
            den2 = {}
            for r in rays:
                l1 = [0] * self.ncyls()
                l2 = [0] * self.ncyls()
                for j,c in enumerate(r):
                    for k in sub1[j]:
                        l1[k] += c
                    for k in sub2[j]:
                        l2[k] += c
                v1 = V(l1)
                v1.set_immutable()
                v2 = V(l2)
                v2.set_immutable()
                if v1 in den1:
                    den1[v1] += 1
                else:
                    den1[v1] = 1
                if v2 in den2:
                    den2[v2] += 1
                else:
                    den2[v2] = 1

            f1 = M.term(num, list(den1.items()))
            f2 = M.term(num, list(den2.items()))

            assert f1 == f2, (f1, f2)

            ans += f1

        return ans

    def _lengths_to_dart_x_coords(self, lengths):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.flat_surfaces.separatrix_diagram import QuadraticCylinderDiagram
            sage: from surface_dynamics import *
            sage: r = RibbonGraph(vertices='(0,3,1,2)(4,7,5,6)',
            ....:                 edges='(0,1)(2,3)(4,5)(6,7)',
            ....:                 faces='(0,3,1,2)(4,7,5,6)', connected=False)
            sage: q = QuadraticCylinderDiagram(r, [1,0])
            sage: q._lengths_to_dart_x_coords(lengths=[2,2,2,2])
            [0, 4, 6, 2, 0, 4, 6, 2]
            sage: q._lengths_to_dart_x_coords(lengths=[1,2,3,4])
            [0, 3, 4, 1, 0, 7, 10, 3]
        """
        G = self._g
        dart_x_coords = [None] * G.num_darts()
        for f in G.faces():
            dart_x_coords[f[0]] = 0
            w = lengths[G._dart_to_edge_index[f[0]]]
            for i in range(1, len(f)):
                dart_x_coords[f[i]] = w
                w += lengths[G._dart_to_edge_index[f[i]]]
        return dart_x_coords

    def _cylcoord_dart_to_corner(self, dart_x_coords, heights, twists, verbose=False):
        r"""
        TESTS::

            sage: from surface_dynamics.flat_surfaces.separatrix_diagram import QuadraticCylinderDiagram
            sage: from surface_dynamics import *

        The pillow::

            sage: q = QuadraticCylinderDiagram('(0,0)-(1,1)')
            sage: x = q._lengths_to_dart_x_coords([1, 1])
            sage: q._cylcoord_dart_to_corner(x, [1], [0])
            [0, 1, 3, 2]
            sage: q._cylcoord_dart_to_corner(x, [1], [1])
            [0, 1, 2, 3]
            sage: x = q._lengths_to_dart_x_coords([2, 2])
            sage: q._cylcoord_dart_to_corner(x, [1], [0])
            [0, 0, 3, 3]
            sage: q._cylcoord_dart_to_corner(x, [1], [1])
            [0, 0, 2, 2]

        An example in `Q(2^2)` (everything is authorized)::

            sage: q = QuadraticCylinderDiagram('(0,1,0,1)-(2,3,2,3)')
            sage: x = q._lengths_to_dart_x_coords([2,2,2,2])
            sage: q._cylcoord_dart_to_corner(x, [1], [0])
            [0, 0, 0, 0, 3, 3, 3, 3]
            sage: q._cylcoord_dart_to_corner(x, [1], [1])
            [0, 0, 0, 0, 2, 2, 2, 2]
            sage: q._cylcoord_dart_to_corner(x, [2], [0])
            [0, 0, 0, 0, 0, 0, 0, 0]
            sage: q._cylcoord_dart_to_corner(x, [2], [1])
            [0, 0, 0, 0, 1, 1, 1, 1]

        Another example in `Q(2^2)`::

            sage: q = QuadraticCylinderDiagram('(0,2,0,3)-(1,2,1,3)')
            sage: x = q._lengths_to_dart_x_coords([1,2,2,1])
            sage: q._cylcoord_dart_to_corner(x, [2], [0])
            [0, 1, 1, 0, 0, 1, 1, 0]

        An example in `Q(3, -1^3)`::

            sage: q = QuadraticCylinderDiagram('(0,1,1)-(0,2,2,3,3)')
            sage: x = q._lengths_to_dart_x_coords([2, 1, 1, 1])
            sage: q._cylcoord_dart_to_corner(x, [2], [0])
            [0, 0, 1, 0, 0, 1, 0, 1]
            sage: q._cylcoord_dart_to_corner(x, [1], [0])
            Traceback (most recent call last):
            ...
            ValueError: invalid lengths/heights/twists
            sage: q._cylcoord_dart_to_corner(x, [1], [1])
            Traceback (most recent call last):
            ...
            ValueError: invalid lengths/heights/twists

        An example in `Q(4, 2, -1^2)`::

            sage: q = QuadraticCylinderDiagram('(0,1,0,2)-(1,3,4,3) (2,4)-(5,5)')
            sage: x = q._lengths_to_dart_x_coords([1, 2, 2, 1, 2, 2])
            sage: q._cylcoord_dart_to_corner(x, [2, 1], [1, 0])
            [0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 3, 3]
            sage: q._cylcoord_dart_to_corner(x, [2, 1], [1, 1])
            [0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 2, 2]
            sage: q._cylcoord_dart_to_corner(x, [2, 1], [0, 0])
            Traceback (most recent call last):
            ...
            ValueError: invalid lengths/heights/twists
        """
        G = self._g
        # pillowcase vertices
        dart_to_corner = [None] * G.num_darts()  # which vertex of the pillowcase the attached
        dart_to_corner[0] = 0
        todo = set([0])  # dart from which we should complete vertices and faces
        cross = set([0]) # dart from which we should cross cylinders
        while todo or cross:

            # 1. completing vertices and faces
            while todo:
                i = todo.pop()
                if verbose:
                    print('** new loop **')
                    print('  i = {}'.format(i))
                    print('  dart_to_corner = {}'.format(dart_to_corner))

                v = dart_to_corner[i]
                assert v is not None

                for j in G.vertex_orbit(i):
                    if dart_to_corner[j] is None:
                        dart_to_corner[j] = v
                        todo.add(j)
                    elif dart_to_corner[j] != dart_to_corner[i]:
                        if verbose:
                            print('i = {} j = {} dart_to_corner[{}] = {} dart_to_corner[{}] = {}'.format(
                                i, j, i, dart_to_corner[i], j, dart_to_corner[j]))
                        raise ValueError('invalid lengths/heights/twists')
                if verbose:
                    print('done completing vertex')
                    print('  dart_to_corner = {}'.format(dart_to_corner))

                f = G.face_orbit(i)
                m = dart_x_coords[i] % 2
                for j in f:
                    if dart_x_coords[j] % 2 == m:
                        if dart_to_corner[j] is None:
                            dart_to_corner[j] = v
                            todo.add(j)
                        elif dart_to_corner[j] != dart_to_corner[i]:
                            if verbose:
                                print('i = {} j = {} dart_to_corner[{}] = {} dart_to_corner[{}] = {}'.format(
                                i, j, i, dart_to_corner[i], j, dart_to_corner[j]))
                            raise ValueError('invalid lengths/heights/twists')
                    else:
                        # 0 <-> 1 and 2 <-> 3
                        vv = v + 1 - 2 * (v % 2)
                        if dart_to_corner[j] is None:
                            dart_to_corner[j] = vv
                            todo.add(j)
                        elif dart_to_corner[j] != vv:
                            if verbose:
                                print('j = {} dart_to_corner[{}] = {}'.format(j,j,dart_to_corner[j]))
                            raise ValueError('invalid lengths/heights/twists')

                cross.add(min(f))

                if verbose:
                    print('done completing face')
                    print('  dart_to_corner = {}'.format(dart_to_corner))

            # 2. crossing cylinders
            while cross:
                # i: dart index
                # j: face index
                i = cross.pop()
                v = dart_to_corner[i]
                j = G._dart_to_face_index[i]
                k = self._face_to_cylinder_index[j]
                jj = self._p[j]
                ii = G._face_cycles[jj][0]

                vv = v
                if heights[k] % 2:
                    # 0 <-> 3 and 1 <-> 2
                    vv = 3 - vv
                if twists[k] % 2:
                    # 0 <-> 1 and 2 <-> 3
                    vv = vv + 1 - 2 * (vv % 2)

                if verbose:
                    print('crossing cylinder from i={} to ii={}'.format(i, ii))
                    print('  dart {} is above pillow vertex {}'.format(i, v))
                    print('  dart {} is above pillow vertex {}'.format(ii, vv))

                if dart_to_corner[ii] is None:
                    dart_to_corner[ii] = vv
                    todo.add(ii)
                elif dart_to_corner[ii] != vv:
                    if verbose:
                        print('ii = {} dart_to_corner[{}] = {}'.format(ii, ii, dart_to_corner[ii]))
                    raise ValueError('invalid lengths/heights/twists')

                if verbose:
                    print('done crossing cylinder {} from dart {}'.format(j, i))
                    print('  dart_to_corner = {}'.format(dart_to_corner))

        if verbose:
            print('dart_to_corner: {}'.format(dart_to_corner))

        return dart_to_corner

    # figure out whether this awfull method is really useful...
    # square tiled quadratic surface would be much more simple to handle
    def cylcoord_to_pillowcase_cover(self, lengths, heights, twists=None, verbose=False):
        r"""
        Convert coordinates of the cylinders into a pillowcase cover.

        The pillow is considered as made of two 1 x 1 squares in order to avoid denominators
        in the input ``lengths``, ``heights`` and ``twists``.

        INPUT:

        - ``lengths`` - positive integers - lengths of the separatrices

        - ``heights`` - positive integers - heights of the cylinders

        - ``twists`` - (optional) non-negative integers - twists. The twist
          is measured as the difference in horizontal coordinates between
          the smallest element in bottom and the smallest in element in top.

        OUTPUT: a pillowcase cover

        EXAMPLES::

            sage: from surface_dynamics import *

        Some pillows in `Q(-1^4)`::

            sage: q = QuadraticCylinderDiagram('(0,0)-(1,1)')
            sage: q.cylcoord_to_pillowcase_cover([1,1],[1],[0])
            g0 = (1)
            g1 = (1)
            g2 = (1)
            g3 = (1)
            sage: q.cylcoord_to_pillowcase_cover([1,1],[1],[1])
            g0 = (1)
            g1 = (1)
            g2 = (1)
            g3 = (1)
            sage: q.cylcoord_to_pillowcase_cover([2,2],[1],[0])
            g0 = (1)(2)
            g1 = (1,2)
            g2 = (1,2)
            g3 = (1)(2)
            sage: q.cylcoord_to_pillowcase_cover([2,2],[1],[1])
            g0 = (1)(2)
            g1 = (1,2)
            g2 = (1)(2)
            g3 = (1,2)
            sage: q.cylcoord_to_pillowcase_cover([2,2],[1],[2]) == q.cylcoord_to_pillowcase_cover([2,2],[1],[2])
            True
            sage: q.cylcoord_to_pillowcase_cover([3,3],[1],[0])
            g0 = (1)(2,3)
            g1 = (1,3)(2)
            g2 = (1,3)(2)
            g3 = (1)(2,3)
            sage: q.cylcoord_to_pillowcase_cover([3,3],[1],[1])
            g0 = (1)(2,3)
            g1 = (1,3)(2)
            g2 = (1)(2,3)
            g3 = (1,2)(3)

            sage: q.cylcoord_to_pillowcase_cover([1,1],[2],[0])
            g0 = (1)(2)
            g1 = (1)(2)
            g2 = (1,2)
            g3 = (1,2)
            sage: q.cylcoord_to_pillowcase_cover([1,1],[2],[1]) == q.cylcoord_to_pillowcase_cover([1,1],[2],[0])
            True

            sage: q.cylcoord_to_pillowcase_cover([2,2],[2],[1])
            g0 = (1)(2)(3,4)
            g1 = (1,2)(3)(4)
            g2 = (1,3)(2,4)
            g3 = (1,4)(2,3)

        Two one cylinder examples in `Q(2^2)`::

            sage: q1 = QuadraticCylinderDiagram('(0,1,0,1)-(2,3,2,3)')
            sage: q2 = QuadraticCylinderDiagram('(0,1,0,2)-(3,1,3,2)')

            sage: q1.cylcoord_to_pillowcase_cover([2,2,2,2],[1],[0])
            g0 = (1,2,3,4)
            g1 = (1,3)(2,4)
            g2 = (1,3)(2,4)
            g3 = (1,4,3,2)
            sage: q1.cylcoord_to_pillowcase_cover([2,2,2,2],[1],[1])
            g0 = (1,2,3,4)
            g1 = (1,3)(2,4)
            g2 = (1,4,3,2)
            g3 = (1,3)(2,4)
            sage: p = q1.cylcoord_to_pillowcase_cover([2,6,4,4],[3],[1])
            sage: p.stratum()
            Q_2(2^2)
            sage: p.nb_pillows()
            24

        One two cylinders example in `Q(2^2)`::

            sage: q = QuadraticCylinderDiagram('(0,1)-(2,3) (0,3)-(1,2)')
            sage: p = q.cylcoord_to_pillowcase_cover([1,1,1,1], [2,2], [0,1])
            sage: p
            g0 = (1,4,2,3)
            g1 = (1,3,2,4)
            g2 = (1,2)(3,4)
            g3 = (1,2)(3,4)
            sage: p.stratum()
            Q_2(2^2)
            sage: q.cylcoord_to_pillowcase_cover([1,3,1,3], [2,2], [0,1]).stratum()
            Q_2(2^2)

        TESTS::

            sage: from surface_dynamics import *
            sage: q = QuadraticCylinderDiagram('(0,0)-(1,1)')
            sage: q.cylcoord_to_pillowcase_cover([1,2],[1])
            Traceback (most recent call last):
            ...
            ValueError: sum of lengths on top and bottom differ
        """
        from sage.rings.integer_ring import ZZ

        G = self._g

        if len(lengths) != G.num_edges():
            raise ValueError("'lengths' has wrong length")
        if len(heights) != self.num_cylinders():
            raise ValueError("'heights' has wrong length")

        lengths = [ZZ.coerce(l) for l in lengths]
        heights = [ZZ.coerce(h) for h in heights]
        ZZ_zero = ZZ.zero()
        for l in lengths:
            if l <= ZZ_zero:
                raise ValueError("each length should be positive")
        for h in heights:
            if h <= ZZ_zero:
                raise ValueError("each height should be positive")

        widths = []
        for bot, top in self.cylinders(dart=False):
            w1 = sum(lengths[i] for i in bot)
            w2 = sum(lengths[i] for i in top)
            if w1 != w2:
                raise ValueError('sum of lengths on top and bottom differ')
            widths.append(w1)

        areas = [heights[i] * widths[i] for i in range(self.ncyls())]
        N = sum(areas)  # number of squares
        if N % 2:
            raise ValueError('got an odd area')
        N = N // 2
        if verbose:
            print('area = {}'.format(N))

        if twists is None:
            twists = [ZZ_zero] * self.num_cylinders()
        elif len(twists) != len(widths):
            raise ValueError("the 'twists' vector has wrong length")
        else:
            twists = [ZZ.coerce(t) % w for (t,w) in zip(twists, widths)]

        # now glue the boundaries
        # (distance with respect to the minimum dart in the face)
        x = self._lengths_to_dart_x_coords(lengths)
        dart_to_corner = self._cylcoord_dart_to_corner(x, heights, twists, verbose)

        # building dart_to_pillow
        if verbose:
            print('Constructing dart_to_pillow')
        dart_to_pillow = [None] * G.num_darts()
        n = 0
        for w,h,t,(bot,top) in zip(widths, heights, twists, self.cylinders(True)):
            w = w // 2

            if verbose:
                print('new cylinder ({},{}) w={} h={} t={}'.format(bot,top,w,h,t))
                print('pillows in [n, n + h*w[ = [{}, {}['.format(n, n+h*w))
            # pillows in the bottom
            dart_to_pillow[bot[0]] = n
            nn = n
            for i in range(len(bot) - 1):
                j = bot[i]
                l = lengths[G._dart_to_edge_index[j]]
                if l % 2:
                    # changing vertex nature (0 <-> 1, 2 <-> 3)
                    v = dart_to_corner[j]
                    if verbose:
                        print('  odd length {} for dart {} at ({},{})'.format(l,j,nn,v))
                    if v == 0 or v == 2:
                        nn += (l - 1) // 2
                    else:
                        nn += (l + 1) // 2
                else:
                    # keeping vertex nature
                    if verbose:
                        print('  even length {} for dart {}'.format(l,j))
                    nn += l // 2

                if nn >= n + w:
                    nn -= w
                if verbose:
                    print('  n={} w={} nn = {}'.format(n,w,nn))
                assert n <= nn < n + w
                dart_to_pillow[bot[i+1]] = nn

            if verbose:
                print('  dart to pillow done on bot')
                print('  dart_to_pillow {}'.format(dart_to_pillow))

            # pillow of top[0] (depends on twist)
            if t % 2 == 0:
                tadj = t // 2
            else:
                v = dart_to_corner[bot[0]]
                if v % 2:
                    tadj = (t+1) // 2
                else:
                    tadj = (t-1) // 2

            nn = dart_to_pillow[bot[0]] + (h-1)*w + tadj
            if nn >= n + h * w:
                nn -= w
            assert n + (h-1)*w <= nn < n + h * w
            dart_to_pillow[top[0]] = nn

            if verbose:
                print('  crossing cylinder: dart {} at ({},{})'.format(
                         top[0], dart_to_pillow[top[0]], dart_to_corner[top[0]]))

            # other pillows on top
            for i in range(len(top) - 1):
                j = top[i]
                l = lengths[G._dart_to_edge_index[j]]
                if l % 2:
                    # changing vertex nature (0 <-> 1, 2 <-> 3)
                    v = dart_to_corner[j]
                    if verbose:
                        print('  odd length {} for dart {} at ({},{})'.format(l,j,nn,v))
                    if v == 0 or v == 2:
                        nn -= (l - 1) // 2
                    else:
                        nn -= (l + 1) // 2
                    if verbose:
                        print('  ends at {}'.format(nn))
                else:
                    # keeping vertex nature
                    if verbose:
                        print('  even length {} for dart {}'.format(l,j))
                    nn -= l // 2

                if nn < n + (h-1)*w:
                    nn += w

                assert n <= nn < n + h * w
                dart_to_pillow[top[i+1]] = nn

            n += h * w

        if verbose:
            print('  dart to pillow done on top')
            print('  dart_to_pillow: {}'.format(dart_to_pillow))
        # safety check
        assert all(i is None or 0 <= i < N for i in dart_to_pillow)

        # fill the cylinders and record the number at a dart
        # (g01 and g23 correspond to horizontal displacements)
        g0 = [None] * N
        g1 = [None] * N
        g2 = [None] * N
        g3 = [None] * N
        g = [g0, g1, g2, g3]
        g01 = [None] * N    # redundant information (will be used to glue edges)
        g23 = [None] * N    # redundant information (will be used to glue edges)
        n = 0
        if verbose:
            print('Filling permutations inside cylinders')
        for w,h,t,(bot,top) in zip(widths, heights, twists, self.cylinders(True)):
            w = w // 2
            if verbose:
                print('** new cylinder **')
                print(' w = {} h = {} t = {} n = {}'.format(w,h,t,n))

            u = dart_to_corner[bot[0]]
            if u < 2:
                h0 = g0
                h1 = g1
                h2 = g2
                h3 = g3
                h01 = g01
                h23 = g23
            else:
                h0 = g2
                h1 = g3
                h2 = g0
                h3 = g1
                h01 = g23
                h23 = g01


            # filling all rows but the top one
            for k in range(h - 1):
                if verbose:
                    print('  filling row k={} from n={}'.format(k, n))
                for i in range(n, n + w):
                    h2[i] = i + w
                    h2[i + w] = i
                    h3[i] = i + w - 1
                    h3[i + w - 1] = i
                    h23[i] = i + 1
                    h01[i + 1] = i

                # adjustments at the gluings (left and right of the cylinders)
                h3[n] = n + 2*w - 1
                h3[n + 2*w - 1] = n
                h23[n + w - 1] = n
                h01[n] = n + w - 1

                # 3 2  -> 1 0
                # 0 1     2 3
                h0, h1, h2, h3 = h2, h3, h0, h1
                h01, h23 = h23, h01
                n += w
                if verbose:
                    print('  done')
                    print('    g0 = {}'.format(g0))
                    print('    g1 = {}'.format(g1))
                    print('    g2 = {}'.format(g2))
                    print('    g3 = {}'.format(g3))
                    print('    g01 = {}'.format(g01))
                    print('    g23 = {}'.format(g23))

            # top row
            if verbose:
                print('  filling top row from n={}'.format(n))
            for i in range(n, n + w - 1):
                h23[i] = i + 1
                h01[i + 1] = i
            h23[n + w - 1] = n
            h01[n] = n + w - 1
            n += w
            if verbose:
                print('  done')
                print('    g0 = {}'.format(g0))
                print('    g1 = {}'.format(g1))
                print('    g2 = {}'.format(g2))
                print('    g3 = {}'.format(g3))
                print('    g01 = {}'.format(g01))
                print('    g23 = {}'.format(g23))

        if verbose:
            print('done filling cylinders:')
            print('  g0 : {}'.format(g0))
            print('  g1 : {}'.format(g1))
            print('  g2 : {}'.format(g2))
            print('  g3 : {}'.format(g3))
            print('  g01: {}'.format(g01))
            print('  g23: {}'.format(g23))

        # safety check
        assert all(i is None or 0 <= i < N for i in g0)
        assert all(i is None or 0 <= i < N for i in g1)
        assert all(i is None or 0 <= i < N for i in g2)
        assert all(i is None or 0 <= i < N for i in g3)
        assert all(i is not None and 0 <= i < N for i in g01)
        assert all(i is not None and 0 <= i < N for i in g23)
        assert all(g0.count(i) < 2 for i in range(N))
        assert all(g1.count(i) < 2 for i in range(N))
        assert all(g2.count(i) < 2 for i in range(N))
        assert all(g3.count(i) < 2 for i in range(N))

        # gluing corners
        for c in G._vertex_cycles:
            i = dart_to_pillow[c[0]]
            u = dart_to_corner[c[0]]
            if verbose:
                print('gluing vertex c = {}'.format(c))
                print('  c[0] = {} at ({},{})'.format(c[0], i, u))

            for k in range(len(c)-1):
                j = dart_to_pillow[c[k+1]]
                v = dart_to_corner[c[k+1]]
                if verbose:
                    print('  c[{}] = {} at ({},{})'.format(k+1, c[k+1], j, v))
                assert u == v
                # safety check
                g[u][i] = j
                i = j
            j = dart_to_pillow[c[0]]
            v = dart_to_corner[c[0]]
            g[u][i] = j

        if verbose:
            print('done gluing corners')
            print('  g0: {}'.format(g0))
            print('  g1: {}'.format(g1))
            print('  g2: {}'.format(g2))
            print('  g3: {}'.format(g3))

        # safety check
        assert all(i is None or 0 <= i < N for i in g0)
        assert all(i is None or 0 <= i < N for i in g1)
        assert all(i is None or 0 <= i < N for i in g2)
        assert all(i is None or 0 <= i < N for i in g3)
        assert all(i is not None and 0 <= i < N for i in g01)
        assert all(i is not None and 0 <= i < N for i in g23)
        assert all(g0.count(i) < 2 for i in range(N))
        assert all(g1.count(i) < 2 for i in range(N))
        assert all(g2.count(i) < 2 for i in range(N))
        assert all(g3.count(i) < 2 for i in range(N))

        # gluing edges
        for s1, s2 in G._edge_cycles:
            t1 = G._faces[s1] # endpoint of first edge
            t2 = G._faces[s2] # endpoint of second edge
            i1 = dart_to_pillow[s1]
            u1 = dart_to_corner[s1]
            j1 = dart_to_pillow[t1]
            v1 = dart_to_corner[t1]
            i2 = dart_to_pillow[s2]
            u2 = dart_to_corner[s2]
            j2 = dart_to_pillow[t2]
            v2 = dart_to_corner[t2]
            if verbose:
                print('gluing edge {} ({},{}) - {} ({},{}) to {} ({},{}) - {} ({},{})'.format(
                         s1,i1,u1,
                         t1,j1,v1,
                         s2,i2,u2,
                         t2,j2,v2))

            # safety check
            assert u1 == v2 and u2 == v1

            # gluing the central part
            j2, v2 = move_backward(j2, v2, g01, g23)
            if verbose:
                print('MOVE BACKWARD ({},{})'.format(j2, v2))
            i1, u1 = move_forward(i1, u1, g01, g23)
            if verbose:
                print('MOVE FORWARD ({},{})'.format(i1, u1))
            while (j2,v2) != (i2,u2):
                if verbose:
                    print('  glue (i1,u1) = ({},{}) and (j2,v2) = ({},{})'.format(i1, u1, j2, v2))
                # safety check
                assert u1 == v2
                g[v2][j2] = i1
                g[v2][i1] = j2
                j2, v2 = move_backward(j2, v2, g01, g23)
                i1, u1 = move_forward(i1, u1, g01, g23)

        if verbose:
            print('done gluing edges')
            print('  (i1,u1) = ({},{})'.format(i1,u1))
            print('  (j1,v1) = ({},{})'.format(j1,v1))
            print('  (i2,u2) = ({},{})'.format(i2,u2))
            print('  (j2,v2) = ({},{})'.format(j2,v2))
            print('  g0: {}'.format(g0))
            print('  g1: {}'.format(g1))
            print('  g2: {}'.format(g2))
            print('  g3: {}'.format(g3))

        # safety check
        assert (i1,u1) == (j1,v1)

        from surface_dynamics.flat_surfaces.origamis.pillowcase_cover import PillowcaseCover
        return PillowcaseCover(g0, g1, g2, g3, as_tuple=True)
