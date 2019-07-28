r"""
Fat graph.

This module is experimental.
"""
#*****************************************************************************
#       Copyright (C) 2019 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, print_function
from six.moves import range, map, zip

from sage.misc.cachefunc import cached_function
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from array import array
from collections import deque
from surface_dynamics.misc.permutation import *

###########################
# Miscellaneous functions #
###########################

def num_and_weighted_num(it):
    from sage.rings.integer_ring import ZZ
    from sage.rings.rational_field import QQ
    s = QQ.zero()
    n = ZZ.zero()
    for _,aut in it:
        n += ZZ.one()
        if aut is None:
            s += QQ.one()
        else:
            s += QQ((1,aut.group_cardinality()))
    return n,s

def list_extrems(l, n):
    if not n:
        raise ValueError
    vdmin = vdmax = l[0]
    for i in range(1, n):
        if l[i] > vdmax:
            vdmax = l[i]
        if l[i] < vdmin:
            vdmin = l[i]
    return (vdmin, vdmax)

#####################
# Fat graph #
#####################
#
# For Abelian strata we should use constellations (= bipartite stuff)

class FatGraph(object):
    r"""
    EXAMPLES:

    The once punctured torus::

        sage: from surface_dynamics import FatGraph

        sage: vp = '(0,2,1,3)'
        sage: ep = '(0,1)(2,3)'
        sage: fp = '(0,2,1,3)'
        sage: FatGraph(vp, ep, fp)
        FatGraph('(0,2,1,3)', '(0,1)(2,3)', '(0,2,1,3)')

    Actually it is enough to specify 2 of the 3 permutations::

        sage: vp = '(0,3,1)(4)(2,5,6,7)'
        sage: ep = '(0,1)(2,3)(4,5)(6,7)'
        sage: fp = '(0,3,7,5,4,2)(1)(6)'
        sage: F0 = FatGraph(vp=vp, ep=ep, fp=fp)
        sage: F1 = FatGraph(ep=ep, fp=fp)
        sage: F2 = FatGraph(vp=vp, fp=fp)
        sage: F3 = FatGraph(vp=vp, ep=ep)
        sage: F0 == F1 and F0 == F2 and F0 == F3
        True
    """
    __slots__ = ['_n',  # number of darts (non-negative integer)
                 '_vp', # vertex permutation (array of length _n)
                 '_ep', # edge permutation (array of length _n)
                 '_fp', # face permutation (array of length _n)
                 # labels
                 # TODO: add _el and care about folded edges!!
                 '_vl', # vertex labels (array of length _n)
                 '_fl', # face labels (array of length _n)
                 # numbers
                 # TODO: add _ne
                 '_nv', # number of vertices (non-negative integer)
                 '_nf', # number of faces (non-negative integer)
                 # degrees
                 # TODO: add _ed
                 '_vd', # vertex degrees (array of length _nv)
                 '_fd'] # face degrees (array of length _nf)

    def __init__(self, vp=None, ep=None, fp=None, max_num_dart=None, check=True):
        vp, ep, fp = constellation_init(vp, ep, fp)
        self._vp = vp
        self._ep = ep
        self._fp = fp
        if len(vp) != len(ep) or len(vp) != len(fp):
            raise ValueError("invalid permutations")
        self._n = len(vp)   # number of darts
        self._nf = 0        # number of faces

        self._vl, self._vd = perm_dense_cycles(vp, self._n)
        self._nv = len(self._vd) # number of vertices
        self._fl, self._fd = perm_dense_cycles(fp, self._n)
        self._nf = len(self._fd) # number of faces

        if max_num_dart is not None:
            if max_num_dart < self._n:
                raise ValueError
            self._realloc(max_num_dart)

        if check:
            self._check()

    def __hash__(self):
        raise TypeError("FatGraph not hashable")

    def _realloc(self, max_num_dart):
            if max_num_dart < self._n:
                return
            self._vp.extend([-1] * (max_num_dart - self._n))
            self._ep.extend([-1] * (max_num_dart - self._n))
            self._fp.extend([-1] * (max_num_dart - self._n))
            self._vl.extend([-1] * (max_num_dart - self._n))
            self._fl.extend([-1] * (max_num_dart - self._n))
            self._vd.extend([-1] * (max_num_dart - self._nv))
            self._fd.extend([-1] * (max_num_dart - self._nf))

    def integral_points(self, b):
        r"""
        Return the edge lengths solution to the face lengths constraint.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: vp = '(0,4,2,5,1,3)'
            sage: ep = '(0,1)(2,3)(4,5)'
            sage: fp = '(0,5)(1,3,4,2)'
            sage: cm = FatGraph(vp, ep, fp)
            sage: cm.integral_points((2,4))
            ((1, 1, 1, 1, 1, 1),)
            sage: cm.integral_points((3,5))
            ((1, 1, 1, 1, 2, 2), (2, 2, 1, 1, 1, 1))
            sage: cm.integral_points((5,11))
            ((1, 1, 3, 3, 4, 4),
             (2, 2, 3, 3, 3, 3),
             (3, 3, 3, 3, 2, 2),
             (4, 4, 3, 3, 1, 1))
        """
        # we want to count solutions in integral edges of the face stuff
        if self.num_faces() != len(b):
            raise ValueError

        ep = self._ep
        fp = self._fp

        ieqs = []
        # positivity
        for i in range(self._n):
            l = [0] * self._n
            l[i] = 1
            ieqs.append([-1] + l)

        eqns = []
        # half edge equations
        for i in range(self._n):
            if ep[i] < i:
                continue
            l = [0] * self._n
            l[i] = 1
            l[ep[i]] = -1
            eqns.append([0] + l)

        # face equations
        for c,f in zip(b,self.faces()):
            l = [0] * self._n
            for i in f:
                l[i] = 1
            eqns.append([-c] + l)

        from sage.geometry.polyhedron.constructor import Polyhedron
        return Polyhedron(ieqs=ieqs, eqns=eqns).integral_points()

    def copy(self):
        """
        EXAMPLES::

            sage: from surface_dynamics import FatGraph
            sage: F = FatGraph.from_unicellular_word([0,1,0,2,3,4,1,4,3,2])
            sage: G = F.copy()
            sage: G._check()
        """
        F = FatGraph.__new__(FatGraph)
        F._vp = self._vp[:]
        F._ep = self._ep[:]
        F._fp = self._fp[:]
        F._n = self._n
        F._nf = self._nf
        F._nv = self._nv
        F._vl = self._vl[:]
        F._vd = self._vd[:]
        F._fl = self._fl[:]
        F._fd = self._fd[:]
        return F

    @staticmethod
    def from_unicellular_word(X):
        r"""
        Build a fat graph from a word on the letters {0, ..., n-1} where
        each letter appears exactly twice.

        EXAMPLES::

            sage: from surface_dynamics import FatGraph
            sage: FatGraph.from_unicellular_word([0,1,0,2,3,4,1,4,3,2])
            FatGraph('(0,3)(1,2,6,7)(4,9)(5,8)', '(0,2)(1,6)(3,9)(4,8)(5,7)', '(0,1,2,3,4,5,6,7,8,9)')
            sage: FatGraph.from_unicellular_word([0,1,2,0,3,2,4,1,3,4])
            FatGraph('(0,6,2,7,9,4)(1,3,5,8)', '(0,3)(1,7)(2,5)(4,8)(6,9)', '(0,1,2,3,4,5,6,7,8,9)')
        """
        n = len(X)
        m = n // 2
        ep = [None] * n
        vp = [None] * n
        fp = list(range(1,n)) + [0]
        symb_to_pos = [None] * m
        for i,k in enumerate(X):
            j = symb_to_pos[k]
            if j is not None:
                ep[i] = j
                ep[j] = i
                vp[(j + 1) % n] = i
                vp[(i + 1) % n] = j
            else:
                symb_to_pos[k] = i
        return FatGraph(vp, ep, fp)

    @staticmethod
    def from_string(s):
        r"""
        Build a fat graph from a serialized string.

        See also :meth:`to_string`.

        EXAMPLES::

            sage: from surface_dynamics import FatGraph
            sage: s = '20_i23017546b98jchedfag_2301547698badcfehgji_12346758ab9igdhfejc0'
            sage: F = FatGraph.from_string(s)
            sage: F.to_string() == s
            True
            sage: FatGraph.from_string('0___')
            FatGraph('()', '()', '()')
        """
        if not isinstance(s, str) or s.count('_') != 3:
            raise ValueError("invalid input")
        n, vp, ep, fp = s.split('_')
        n = int(n)
        vp = perm_from_base64_str(vp, n)
        ep = perm_from_base64_str(ep, n)
        fp = perm_from_base64_str(fp, n)
        return FatGraph(vp, ep, fp)

    def to_string(self):
        r"""
        Serialization to string.

        EXAMPLES::

            sage: from surface_dynamics import FatGraph
            sage: FatGraph.from_unicellular_word([0,1,0,2,3,4,1,4,3,2]).to_string()
            '10_3260987154_2609871543_1234567890'
            sage: FatGraph('', '', '').to_string()
            '0___'
        """
        n = self._n
        return str(n) + "_" + \
               perm_base64_str(self._vp, n) + "_" + \
               perm_base64_str(self._ep, n) + "_" + \
               perm_base64_str(self._fp, n)

    def _check(self, error=RuntimeError):
        vp = self._vp
        vl = self._vl
        vd = self._vd

        ep = self._ep

        fp = self._fp
        fl = self._fl
        fd = self._fd

        n = self._n
        nf = self._nf
        nv = self._nv

        m = sum(vp[i] != -1 for i in range(n))

        if not perm_check(vp, n):
            raise ValueError("invalid vertex permutation: %s" % vp)
        if not perm_check(ep, n):
            raise ValueError("invalid edge permutation: %s" % ep)
        if not perm_check(fp, n):
            raise ValueError("invalid face permutation: %s" % fp)

        if perm_num_cycles(vp, n) != self._nv:
            raise error("wrong number of vertices")
        if perm_num_cycles(fp, n) != self._nf:
            raise error("wrong number of faces")

        if len(vl) < n or len(fl) < n or len(vd) < nv or len(fd) < nf:
               raise error("inconsistent lengths")

        if any(x < 0 or x > n for x in vd[:nv]) or sum(vd[:nv]) != m:
            raise error("invalid vertex degrees")
        if any(x < 0 or x > n for x in fd[:nf]) or sum(fd[:nf]) != m:
            raise error("invalid face degrees")

        ffd = [0] * nf
        vvd = [0] * nv

        for i in range(n):
            if vp[i] == -1:
                if ep[i] != -1 or fp[i] != -1:
                    raise ValueError("inconsistent dart activity for i={}".format(i))
                continue
            elif ep[i] == -1 or fp[i] == -1:
                raise ValueError("inconsistent dart activity for i={}".format(i))

            if fp[ep[vp[i]]] != i:
                raise error("fp[ep[vp[%d]]] = %d" % (i, fp[ep[vp[i]]]))
            if fl[i] < 0 or fl[i] >= nf:
                raise error("face label out of range: fl[%d] = %d" % (i, fl[i]))
            if vl[i] < 0 or vl[i] >= nv:
                raise error("vertex label out of range: vl[%d] = %d" % (i, vl[i]))

            if fl[fp[i]] != fl[i]:
                raise error("fl[fp[%d]] = %d while fl[%d] = %d" %(i, fl[fp[i]], i, fl[i]))
            
            if vl[vp[i]] != vl[i]:
                raise error("vl[vp[%d]] = vl[%d] = %d while vl[%d] = %d" %(i, vp[i], vl[vp[i]], i, vl[i]))

            ffd[fl[i]] += 1
            vvd[vl[i]] += 1

        if vvd != vd[:nv]:
            raise error("inconsistent face labels/degrees, got %s instead of vd = %s" % (vvd, vd[:nv]))
        if ffd != fd[:nf]:
            raise error("inconsistent vertex labels/degrees, got %s instead of fd = %s" %(ffd, fd[:nf]))

    def is_face_bipartite(self):
        r"""
        Test whether the faces admit a bi-coloring.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import *
            sage: vp = '(0,2,1,3)'
            sage: ep = '(0,1)(2,3)'
            sage: fp = '(0,2,1,3)'
            sage: F = FatGraph(vp, ep, fp, 6)
            sage: F.is_face_bipartite()
            False
            sage: F.split_face(0,1)
            sage: F.is_face_bipartite()
            True

            sage: vp = '(0,5,1,2,3,4)'
            sage: ep = '(0,2)(1,3)(4,5)'
            sage: fp = '(0,1,2,4)(3,5)'
            sage: FatGraph(vp, ep, fp).is_face_bipartite()
            False

            sage: from surface_dynamics.topology.fat_graph_exhaustive_generation import FatGraphs_g_nf_nv
            sage: F = FatGraphs_g_nf_nv(1, 3, 3, vertex_min_degree=3)
            sage: F.cardinality_and_weighted_cardinality(filter=lambda x,a: x.is_face_bipartite())
            (3, 5/3)
        """
        # trivial cases
        if self._nf == 0:
            return True
        elif self._nf == 1:
            return False

        n = self._n
        ep = self._ep
        fp = self._fp
        fl = self._fl
        nf = self._nf
        colors = [-1] * nf
        edge_seen = [False] * n
        to_test = perm_orbit(fp, 0)
        colors[self._fl[0]] = 1
        while to_test:
            e1 = to_test.pop()
            if edge_seen[e1]:
                continue
            e2 = ep[e1]
            f1 = fl[e1]
            f2 = fl[e2]
            if colors[f1] == -1:
                raise RuntimeError
            elif colors[f2] == -1:
                # discover a new face
                colors[f2] = 1 - colors[f1]
                to_test.extend(perm_orbit(fp, e2))
            elif colors[f1] == colors[f2]:
                # contradiction in colors
                return False

            edge_seen[e1] = edge_seen[e2] = True

        return True

    def __copy__(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph

            sage: vp = '(0,2,1,3)'
            sage: ep = '(0,1)(2,3)'
            sage: fp = '(0,2,1,3)'
            sage: cm = FatGraph(vp, ep, fp)
            sage: cm2 = cm.__copy__()
            sage: cm2._check()
        """
        cm = FatGraph.__new__(FatGraph)
        cm._vp = self._vp[:]
        cm._ep = self._ep[:]
        cm._fp = self._fp[:]
        cm._n = self._n
        cm._nf = self._nf
        cm._nv = self._nv
        cm._vl = self._vl[:]
        cm._fl = self._fl[:]
        cm._vd = self._vd[:]
        cm._fd = self._fd[:]
        return cm

    def __repr__(self):
        n = self._n
        fd = self._fd[:self._nf]
        vd = self._vd[:self._nv]
        fd.sort(reverse=True)
        vd.sort(reverse=True)
        return "FatGraph('%s', '%s', '%s')" % (perm_cycle_string(self._vp, True, n),
                                               perm_cycle_string(self._ep, True, n),
                                               perm_cycle_string(self._fp, True, n))

    def __eq__(self, other):
        r"""
        TESTS::

            sage: from surface_dynamics.topology.fat_graph import FatGraph

            sage: vp = '(0,2,1,3)'
            sage: ep = '(0,1)(2,3)'
            sage: fp = '(0,2,1,3)'
            sage: cm1 = FatGraph(vp, ep, fp)
            sage: cm2 = FatGraph(vp, ep, fp, 100)
            sage: cm1 == cm2
            True
        """
        if type(self) != type(other):
            raise TypeError

        if self._n != other._n or self._nf != other._nf or self._nv != other._nv:
            return False

        for i in range(self._n):
            if self._vp[i] != other._vp[i] or \
               self._ep[i] != other._ep[i] or \
               self._fp[i] != other._fp[i]:
                   return False

        # here we ignore the vertex and face labels...
        return True

    def __ne__(self, other):
        return not self == other

    def vertex_permutation(self, copy=True):
        if copy:
            return self._vp[:self._n]
        else:
            return self._vp

    def edge_permutation(self, copy=True):
        if copy:
            return self._ep[:self._n]
        else:
            return self._ep

    def face_permutation(self, copy=True):
        if copy:
            return self._fp[:self._n]
        else:
            return self._fp

    def vertex_profile(self):
        return perm_cycle_type(self._vp, self._n)

    def edge_profile(self):
        return perm_cycle_type(self._ep, self._n)

    def face_profile(self):
        return perm_cycle_type(self._fp, self._n)

    def profile(self):
        return (self.vertex_profile(), self.edge_profile(), self.face_profile())

    def num_darts(self):
        return self._n

    def num_folded_edges(self):
        return sum(self._ep[i] == 1 for i in range(self._n))

    def num_faces(self):
        return self._nf

    def num_vertices(self):
        return self._nv

    def vertices(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: vp = '(0,4,2,5,1,3)'
            sage: ep = '(0,1)(2,3)(4,5)'
            sage: fp = '(0,5)(1,3,4,2)'
            sage: cm = FatGraph(vp, ep, fp)
            sage: cm.vertices()
            [[0, 4, 2, 5, 1, 3]]
        """
        return perm_cycles(self._vp, True, self._n)

    def vertex_degrees(self):
        return self._vd[:self._nv]

    def vertex_degree_extrems(self):
        return list_extrems(self._vd, self._nv)

    def vertex_degree_min(self):
        return list_extrems(self._vd, self._nv)[0]
    
    def vertex_degree_max(self):
        return list_extrems(self._vd, self._nv)[1]

    def face_degree_extrems(self):
        return list_extrems(self._fd, self._nf)

    def face_degree_min(self):
        return list_extrems(self._fd, self._nf)[0]

    def face_degree_max(self):
        return list_extrems(Self._fd, self._nf)[1]

    def face_degree_min(self):
        return s

    def edges(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: vp = '(0,4,2,5,1,3)'
            sage: ep = '(0,1)(2,3)(4,5)'
            sage: fp = '(0,5)(1,3,4,2)'
            sage: cm = FatGraph(vp, ep, fp)
            sage: cm.edges()
            [[0, 1], [2, 3], [4, 5]]
        """
        return perm_cycles(self._ep, True, self._n)

    def num_edges(self):
        return self._n // 2

    def faces(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: vp = '(0,4,2,5,1,3)'
            sage: ep = '(0,1)(2,3)(4,5)'
            sage: fp = '(0,5)(1,3,4,2)'
            sage: cm = FatGraph(vp, ep, fp)
            sage: cm.faces()
            [[0, 5], [1, 3, 4, 2]]
        """
        return perm_cycles(self._fp, True, self._n)

    def face_degrees(self):
        return self._fd[:self._nf]

    def euler_characteristic(self):
        return self._nf - self._n//2 + self._nv

    def dual(self):
        r"""
        Return the dual fat graph.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: F = FatGraph(vp=None,ep='(0,1)',fp='(0)(1)')
            sage: F.dual()
            sage: F
            FatGraph('(0)(1)', '(0,1)', '(0,1)')
            sage: F._check()
            sage: s = '20_i23017546b98jchedfag_2301547698badcfehgji_12346758ab9igdhfejc0'
            sage: F = FatGraph.from_string(s)
            sage: F.dual()
            sage: F._check()
            sage: F.dual()
            sage: F._check()
            sage: F == FatGraph.from_string(s)
            True
        """
        # TODO: invert in place !!!!
        self._vp, self._fp = perm_invert(self._fp, self._n), perm_invert(self._vp, self._n)
        self._nv, self._nf = self._nf, self._nv
        self._vl, self._fl = self._fl, self._vl
        self._vd, self._fd = self._fd, self._vd

    def kontsevich_volume_rational_function(self, R=None):
        r"""
        This is not under an appropriate form...
        """
        raise NotImplementedError
        print('This is not quite the form under which we would like it... it should remains factorized')
        from sage.rings.rational_field import QQ
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        nf = self._nf
        fl = self._fl
        n = self._n
        ep = self._ep
        if R is None:
            R = PolynomialRing(QQ, 'b', nf)
        gens = R.gens()
        res = R.one()
        for i in range(self._n):
            j = ep[i]
            if j < i:
                continue
            res *= 1 / self.automorphism_group().group_cardinality() / (gens[fl[i]] + gens[fl[j]])
        return res

    ##############################
    # Augmentation and reduction #
    ##############################

    def _check_alloc(self, n, nv, nf):
        if len(self._vp) < n or \
           len(self._ep) < n or \
           len(self._fp) < n or \
           len(self._vl) < n or \
           len(self._fl) < n or \
           len(self._fd) < nf or \
           len(self._vd) < nv:
               raise TypeError("reallocation needed")

    def split_face(self, i, j):
        r"""
        Insert an edge between the darts ``i`` and ``j`` to split the face.

        One of the face will contains i, fp[i], ...,  (the x-face) and the 
        other one will contain j, fp[j], ... In the special case i=j, a
        monogon (= face with only one edge) is created.

        The converse operation is implemented in :meth:`remove_edge`.

        EXAMPLES:

        The once punctured torus::

            sage: from surface_dynamics.topology.fat_graph import FatGraph

            sage: vp = '(0,2,1,3)'
            sage: ep = '(0,1)(2,3)'
            sage: fp = '(0,2,1,3)'
            sage: eps = '(0,1)(2,3)(4,5)'

            sage: vp20 = '(0,4,2,5,1,3)'
            sage: fp20 = '(0,5)(1,3,4,2)'
            sage: cm = FatGraph(vp, ep, fp, 6)
            sage: cm.split_face(2,0)
            sage: cm == FatGraph(vp20, eps, fp20)
            True

            sage: vp10 = '(0,4,2,1,5,3)'
            sage: fp10 = '(0,2,5)(1,3,4)'
            sage: cm = FatGraph(vp, ep, fp, 6)
            sage: cm.split_face(1,0)
            sage: cm == FatGraph(vp10, eps, fp10)
            True

            sage: vp30 = '(0,4,2,1,3,5)'
            sage: fp30 = '(0,2,1,5)(3,4)'
            sage: cm = FatGraph(vp, ep, fp, 6)
            sage: cm.split_face(3,0)
            sage: cm == FatGraph(vp30, eps, fp30)
            True

            sage: vp00 = '(0,5,4,2,1,3)'
            sage: fp00 = '(0,2,1,3,4)(5)'
            sage: cm = FatGraph(vp, ep, fp, 6)
            sage: cm.split_face(0,0)
            sage: cm == FatGraph(vp00, eps, fp00)
            True

            sage: vp22 = '(0,2,5,4,1,3)'
            sage: fp22 = '(0,4,2,1,3)(5)'
            sage: cm = FatGraph(vp, ep, fp, 6)
            sage: cm.split_face(2,2)
            sage: cm == FatGraph(vp22, eps, fp22)
            True

        A genus 2 surface::

            sage: vp = '(0,3,6,8)(1,10,9,12,5)(2,7,4,11)(13)'
            sage: ep = '(0,1)(2,3)(4,5)(6,7)(8,9)(10,11)(12,13)'
            sage: fp = '(0,5,7,3,11,1,8,10,4,12,13,9,6,2)'
            sage: cm = FatGraph(vp, ep, fp, 21)
            sage: cm.split_face(0,1); cm._check()
            sage: cm.split_face(4,13); cm._check()
            sage: cm.split_face(5,14); cm._check()
            sage: cm.remove_edge(18); cm._check()
            sage: cm.remove_edge(16); cm._check()
            sage: cm.remove_edge(14); cm._check()
            sage: cm == FatGraph(vp, ep, fp, 21)
            True
        """
        vp = self._vp
        vl = self._vl
        vd = self._vd

        ep = self._ep

        fp = self._fp
        fl = self._fl
        fd = self._fd

        n = self._n
        nf = self._nf
        nv = self._nv

        i = int(i)
        j = int(j)
        if i < 0 or i >= self._n or j < 0 or j >= n or fl[i] != fl[j]:
            raise ValueError("invalid darts i=%d and j=%d for face splitting" %(i, j))

        self._check_alloc(n + 2, nv, nf + 1)

        x = self._n
        y = self._n + 1
        ii = ep[vp[i]] # = fp^-1(i)
        jj = ep[vp[j]] # = fp^-1(j)

        ep[x] = y
        ep[y] = x

        self._n += 2
        self._nf += 1

        if i == j:
            # add a monogon
            # fp (i A) -> (i A x)(y)
            # vp (... i ...) -> (... i y x ...)
            fp[ii] = x
            fp[x] = i
            fp[y] = y
            vp[x] = vp[i]
            vp[y] = x
            vp[i] = y


            vl[x] = vl[y] = vl[i]
            fl[x] = fl[i]
            fl[y] = nf

            fd[fl[i]] += 1
            fd[nf] = 1

            vd[vl[i]] += 2

        else:
            # general case
            # update permutations:
            # fp  (i A j B)               -> (i A x) (j B y)
            # ep                          -> (x y)
            # vp  (... i ...) (... j ...) -> (... i y ...) (... j x ...)
            fp[jj] = x
            fp[x] = i
            fp[ii] = y
            fp[y] = j

            vp[y] = vp[i]
            vp[i] = y
            vp[x] = vp[j]
            vp[j] = x

            # update labels and degrees
            vl[x] = vl[j]
            vl[y] = vl[i]

            fl[x] = fl[i]
            fl[y] = fl[j]

            dfy = 0  # degree of the y-face
            while fl[y] != nf:
                fl[y] = nf
                y = fp[y]
                dfy += 1
            dfx = fd[fl[x]] + 2 - dfy

            fd[fl[x]] = dfx
            fd[fl[y]] = dfy
            vd[vl[x]] += 1
            vd[vl[y]] += 1

    def remove_edge(self, i):
        r"""
        Remove an edge.

        If the edge has the same face on both sides, then the genus drops by 1.
        Inverse operation of :meth:`split_face` or :meth:`trisect_face`.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph

            sage: vp = '(0,2,1,3)'
            sage: ep = '(0,1)(2,3)'
            sage: fp = '(0,2,1,3)'
            sage: cm = FatGraph(vp, ep, fp)

            sage: eps = '(0,1)(2,3)(4,5)'
            sage: vp20 = '(0,5,4,2,1,3)'
            sage: fp20 = '(0,2,1,3,4)(5)'
            sage: cm2 = FatGraph(vp20, eps, fp20, 6)
            sage: cm2.remove_edge(4)
            sage: cm2 == cm
            True
            sage: cm2 = FatGraph(vp20, eps, fp20, 6)
            sage: cm2.remove_edge(5)
            sage: cm2 == cm
            True

            sage: vp10 = '(0,4,2,5,1,3)'
            sage: fp10 = '(0,5)(1,3,4,2)'
            sage: cm2 = FatGraph(vp10, eps, fp10)
            sage: cm2.remove_edge(4)
            sage: cm2 == cm
            True
            sage: cm2 = FatGraph(vp10, eps, fp10)
            sage: cm2.remove_edge(5)
            sage: cm2 == cm
            True

            sage: vp30 = '(0,4,2,1,5,3)'
            sage: fp30 = '(0,2,5)(1,3,4)'
            sage: cm2 = FatGraph(vp30, eps, fp30)
            sage: cm2.remove_edge(4)
            sage: cm2 == cm
            True
            sage: cm2 = FatGraph(vp30, eps, fp30)
            sage: cm2.remove_edge(5)
            sage: cm2 == cm
            True

            sage: vp00 = '(0,5,4,2,1,3)'
            sage: fp00 = '(0,2,1,3,4)(5)'
            sage: cm2 = FatGraph(vp00, eps, fp00)
            sage: cm2.remove_edge(4)
            sage: cm2 == cm
            True

            sage: vp22 = '(0,2,5,4,1,3)'
            sage: fp22 = '(0,4,2,1,3)(5)'
            sage: cm2 = FatGraph(vp00, eps, fp00)
            sage: cm2.remove_edge(4)
            sage: cm2 == cm
            True
        """
        vp = self._vp
        ep = self._ep
        fp = self._fp

        vl = self._vl
        fl = self._fl

        vd = self._vd
        fd = self._fd

        n = self._n
        nf = self._nf
        nv = self._nv

        i = int(i)
        if i < 0 or i >= self._n:
            raise ValueError("dart index out of range")
        j = ep[i]

        fi = fl[i]
        fj = fl[j]
        if fi == fj:
            raise ValueError("i=%d and j=%d on the same face" %(i,j))
        fmin = min(fi, fj)
        if i < n - 2 or j < n - 2 or max(fi, fj) != nf-1:
            raise NotImplementedError

        ii = ep[vp[i]]
        jj = ep[vp[j]]
        if fd[fl[i]] == 1:
            # monogon
            assert vp[i] == j
            fp[jj] = fp[j]
            vp[fp[j]] = vp[j]
        elif fd[fl[j]] == 1:
            # monogon
            assert vp[j] == i
            fp[ii] = fp[i]
            vp[fp[i]] = vp[i]
        else:
            # none of them are monogons
            fp[ii] = fp[j]
            fp[jj] = fp[i]

            vp[fp[j]] = vp[i]
            vp[fp[i]] = vp[j]

        # update vertex and face degrees
        vd[vl[i]] -= 1
        vd[vl[j]] -= 1

        d = fd[fl[i]] + fd[fl[j]] - 2
        fd[fmin] = d

        # update face labels
        k = fp[i]
        while fl[k] != fmin:
            fl[k] = fmin
            k = fp[k]
        k = fp[j]
        while fl[k] != fmin:
            fl[k] = fmin
            k = fp[k]

        self._n -= 2
        self._nf -= 1

    def split_vertex(self, i, j):
        r"""
        Insert a new edge to split the vertex located at the darts i and j.
        
        This operation keeps the genus constant. The inverse operation is implemented
        in :meth:`contract_edge`.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph

            sage: vp = '(0,2,1,3)'
            sage: ep = '(0,1)(2,3)'
            sage: fp = '(0,2,1,3)'
            sage: eps = '(0,1)(2,3)(4,5)'

            sage: vp02 = '(0,4,1,3)(2,5)'
            sage: fp02 = '(0,4,2,1,3,5)'
            sage: cm = FatGraph(vp, ep, fp, 6)
            sage: cm.split_vertex(0,2)
            sage: cm == FatGraph(vp02, eps, fp02)
            True

            sage: vp01 = '(0,4,3)(1,5,2)'
            sage: fp01 = '(0,2,4,1,3,5)'
            sage: cm = FatGraph(vp, ep, fp, 6)
            sage: cm.split_vertex(0,1)
            sage: cm == FatGraph(vp01, eps, fp01)
            True

            sage: vp03 = '(0,4)(1,3,5,2)'
            sage: fp03 = '(0,2,1,4,3,5)'
            sage: cm = FatGraph(vp, ep, fp, 6)
            sage: cm.split_vertex(0,3)
            sage: cm == FatGraph(vp03, eps, fp03)
            True
        """
        vp = self._vp
        ep = self._ep
        fp = self._fp
        vl = self._vl
        fl = self._fl
        n = self._n
        nf = self._nf
        nv = self._nv
        vd = self._vd
        fd = self._fd

        i = int(i)
        j = int(j)
        if i < 0 or i >= self._n or j < 0 or j >= self._n or vl[i] != vl[j]:
            raise ValueError("invalid darts i=%d and j=%d for vertex splitting" %(i, j))

        self._check_alloc(n + 2, nv + 1, nf)

        x = self._n
        y = self._n + 1
        ii = vp[i]
        jj = vp[j]
        ep[x] = y
        ep[y] = x
        self._n += 2
        self._nv += 1

        if i == j:
            # introduce a vertex of degree 1
            vp[x] = ii
            vp[i] = x
            vp[y] = y

            fp[y] = i
            fp[x] = y
            fp[ep[ii]] = x

            fl[x] = fl[y] = fl[i]
            vl[x] = vl[i]
            vl[y] = nv

            vd[vl[x]] += 1
            vd[vl[y]] = 1

            fd[fl[x]] += 2

        else:
            # general case
            # update permutations
            # fp (... i ...) (... j ...) -> (... y i ...) (... x j ...)
            # ep                         -> (x y)
            # vp (A i B j)               -> (A i x) (B j y)
            vp[x] = jj
            vp[i] = x
            vp[y] = ii
            vp[j] = y

            fp[ep[jj]] = x
            fp[x] = j
            fp[ep[ii]] = y
            fp[y] = i 

            # update labels and degrees

            fl[x] = fl[j]
            fl[y] = fl[i]

            vl[x] = vl[i]
            vl[y] = vl[j]

            dvy = 0
            while vl[y] != nv:
                vl[y] = nv
                y = vp[y]
                dvy += 1
            dvx = vd[vl[x]] + 2 - dvy

            vd[vl[x]] = dvx
            vd[vl[y]] = dvy
            fd[fl[x]] += 1
            fd[fl[y]] += 1

    def contract_edge(self, i):
        r"""
        Contract an edge between two distinct zeros.

        Inverse operation of :meth:`split_vertex` except that here we allow
        vertices of degree one.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph

            sage: vp = '(0,2,1,3)'
            sage: ep = '(0,1)(2,3)'
            sage: fp = '(0,2,1,3)'
            sage: eps = '(0,1)(2,3)(4,5)'

            sage: vp02 = '(0,4,1,3)(2,5)'
            sage: fp02 = '(0,4,2,1,3,5)'
            sage: cm = FatGraph(vp02, eps, fp02)
            sage: cm.contract_edge(4)
            sage: cm == FatGraph(vp, ep, fp)
            True
            sage: cm = FatGraph(vp02, eps, fp02)
            sage: cm.contract_edge(5)
            sage: cm == FatGraph(vp, ep, fp)
            True

            sage: vp01 = '(0,4,3)(1,5,2)'
            sage: fp01 = '(0,2,4,1,3,5)'
            sage: cm = FatGraph(vp01, eps, fp01)
            sage: cm.contract_edge(4)
            sage: cm == FatGraph(vp, ep, fp)
            True
            sage: cm = FatGraph(vp01, eps, fp01)
            sage: cm.contract_edge(5)
            sage: cm == FatGraph(vp, ep, fp)
            True

            sage: vp03 = '(0,4)(1,3,5,2)'
            sage: fp03 = '(0,2,1,4,3,5)'
            sage: cm = FatGraph(vp03, eps, fp03)
            sage: cm.contract_edge(4)
            sage: cm == FatGraph(vp, ep, fp)
            True
            sage: cm = FatGraph(vp03, eps, fp03)
            sage: cm.contract_edge(5)
            sage: cm == FatGraph(vp, ep, fp)
            True

        Degree 1 vertices::

            sage: cm = FatGraph('(0,2)(1)(3)', '(0,1)(2,3)', '(0,1,2,3)')
            sage: cm.contract_edge(2)
            sage: cm
            FatGraph('(0)(1)', '(0,1)', '(0,1)')
            sage: cm2 = FatGraph('(0,2)(1)(3)', '(0,1)(2,3)', '(0,1,2,3)')
            sage: cm2.contract_edge(3)
            sage: cm == cm2
            True
        """
        vp = self._vp
        ep = self._ep
        fp = self._fp

        vl = self._vl
        fl = self._fl

        vd = self._vd
        fd = self._fd

        n = self._n
        nf = self._nf
        nv = self._nv

        i = int(i)
        if i < 0 or i >= self._n:
            raise ValueError("dart index out of range")
        j = ep[i]
        if vl[i] == vl[j]:
            raise ValueError("i=%d and j=%d on the same vertex" %(i,j))

        vi = vl[i]
        vj = vl[j]
        vmin = min(vi, vj)
        if i < n - 2 or j < n - 2 or max(vi, vj) != nv-1:
            raise NotImplementedError

        ii = ep[vp[i]]
        jj = ep[vp[j]]
        if vd[vl[i]] == 1:
            # vertex of degree one
            assert fp[j] == i
            vp[fp[i]] = vp[j]
            fp[jj] = fp[i]

        elif vd[vl[j]] == 1:
            # vertex of degree one
            assert fp[i] == j
            vp[fp[j]] = vp[i]
            fp[ii] = fp[j]

        else:
            vp[fp[i]] = vp[i]
            vp[fp[j]] = vp[j]
            fp[ii] = fp[i]
            fp[jj] = fp[j]

        # update vertex and face degree
        fd[fl[i]] -= 1
        fd[fl[j]] -= 1

        d = vd[vl[i]] + vd[vl[j]] - 2
        vd[vmin] = d

        # update vertex labels
        k = vp[i]
        while vl[k] != vmin:
            vl[k] = vmin
            k = vp[k]
        k = vp[j]
        while vl[k] != vmin:
            vl[k] = vmin
            k = vp[k]

        self._n -= 2
        self._nv -= 1

    def trisect_face(self, i, j, k):
        r"""
        Insert a bridge

        INPUT:

        - ``i``, ``j``, ``k`` - dart in the same face in counter-clockwise
          order


        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph

            sage: vp = '(0,2,1,3)'
            sage: ep = '(0,1)(2,3)'
            sage: fp = '(0,2,1,3)'
            sage: cm = FatGraph(vp, ep, fp, 8)

            sage: vp021 = '(0,7,2,6,5,1,4,3)'
            sage: ep021 = '(0,1)(2,3)(4,5)(6,7)'
            sage: fp021 = '(0,5,1,3,7,2,4,6)'
            sage: cm021 = FatGraph(vp021, ep021, fp021)
            sage: cm.trisect_face(0, 2, 1)
            sage: cm == cm021
            True

            sage: cm = FatGraph(vp, ep, fp, 10)
            sage: cm.trisect_face(0, 0, 3)

            sage: cm = FatGraph(vp, ep, fp, 10)
            sage: cm.trisect_face(0, 3, 3)

            sage: cm = FatGraph(vp, ep, fp, 10)
            sage: cm.trisect_face(0, 3, 0)

            sage: cm = FatGraph(vp, ep, fp, 10)
            sage: cm.trisect_face(0, 0, 0)
        """
        vp = self._vp
        ep = self._ep
        fp = self._fp

        vl = self._vl
        fl = self._fl

        vd = self._vd
        fd = self._fd

        n = self._n
        nf = self._nf
        nv = self._nv

        i = int(i)
        j = int(j)
        k = int(k)

        if i < 0 or i >= n or j < 0 or j >= n or k < 0 or k >= n:
            raise ValueError("dart index out of range")
        if fl[i] != fl[j] or fl[i] != fl[k]:
            raise ValueError("darts in distinct faces")

        self._check_alloc(n + 4, nv, nf)
        self._n += 4

        ii = ep[vp[i]] # = fp^-1(i) at the end of B
        jj = ep[vp[j]] # = fp^-1(j) at the end of A
        kk = ep[vp[k]] # = fp^-1(k) at the end of C

        x = n
        y = n + 1
        xx = n + 2
        yy = n + 3

        ep[x] = y
        ep[y] = x
        ep[xx] = yy
        ep[yy] = xx

        fl[x] = fl[y] = fl[xx] = fl[yy] = fl[i]
        vl[x] = vl[k]
        vl[xx] = vl[y] = vl[j]
        vl[yy] = vl[i]

        fd[fl[i]] += 4
        vd[vl[i]] += 1
        vd[vl[j]] += 2
        vd[vl[k]] += 1

        if i == j == k:
            # face: -> (x xx y yy j C) 
            # (j C kk)
            vp[x] = vp[j]
            vp[yy] = x
            vp[y] = yy
            vp[xx] = y
            vp[j] = xx

            fp[kk] = x
            fp[x] = xx
            fp[xx] = y
            fp[y] = yy
            fp[yy] = j
        elif i == j:
            # face: -> (x xx y k B yy j C)
            # (j C kk) (k B ii)
            vp[yy] = vp[j]
            vp[y] = yy
            vp[xx] = y
            vp[j] = xx
            vp[x] = vp[k]
            vp[k] = x

            fp[ii] = yy
            fp[yy] = j
            fp[kk] = x
            fp[x] = xx
            fp[xx] = y
            fp[y] = k
        elif j == k:
            # face: -> (x xx i A y k B yy)
            # (i A jj) (k B ii)
            vp[yy] = vp[i]
            vp[i] = yy
            vp[y] = vp[k]
            vp[xx] = y
            vp[x] = xx
            vp[k] = x

            fp[ii] = yy
            fp[yy] = x
            fp[x] = xx
            fp[xx] = i
            fp[jj] = y
            fp[y] = k
        elif k == i:
            # face: -> (x xx i A y yy j C)
            # (i A jj) (j C kk)
            vp[y] = vp[j]
            vp[xx] = y
            vp[j] = xx
            vp[x] = vp[i]
            vp[yy] = x
            vp[i] = yy

            fp[kk] = x
            fp[x] = xx
            fp[xx] = i
            fp[jj] = y
            fp[y] = yy
            fp[yy] = j
        else:
            # general case
            # vertex: (...i...)(...j...)(...k...) -> (...i yy...)(...j xx y...)(...k x...)
            # edge  : add (x y) (xx yy)
            # face  : (i A j C k B) -> (x xx i A y k B yy j C)
            #
            # (i A jj) (j C kk) (k B ii)
            vp[yy] = vp[i]
            vp[i] = yy
            vp[y] = vp[j]
            vp[xx] = y
            vp[j] = xx
            vp[x] = vp[k]
            vp[k] = x

            fp[kk] = x
            fp[x] = xx
            fp[xx] = i
            fp[jj] = y
            fp[y] = k
            fp[ii] = yy
            fp[yy] = j

    def remove_face_trisection(self, x):
        ep = self._ep
        fp = self._fp
        fl = self._fl
        fd = self._fd
        vp = self._vp
        vl = self._vl
        vd = self._vd

        x = int(x)
        xx = fp[x]
        y = ep[x]
        yy = ep[xx]

        if fl[x] != fl[y] or fl[x] != fl[xx] or fl[x] != fl[yy]:
            raise ValueError("not a trisection")

        # face: (x xx i A y k B yy j C) -> (i A j C k B)
        #    -> (i A jj) (j C kk) (k B ii)
        i = fp[xx]
        k = fp[y]
        j = fp[yy]
        ii = ep[vp[yy]] # = fp^-1(yy)
        jj = ep[vp[y]]  # = fp^-1(y)
        kk = ep[vp[x]]  # = fp^-1(x)

        if fp[xx] == y and fp[y] == yy:
            # vertex (... j xx y yy x ...) -> (... j ...)
            # face (x xx y yy j C) -> (j C)
            # (j C kk)
            assert vp[j] == xx
            assert vp[xx] == y
            assert vp[yy] == x
            vp[j] = vp[x]
            fp[kk] = j
        elif fp[xx] == y:
            # face: (x xx y k B yy j C) -> (j C k B)
            # (j C kk) (k B ii)
            assert fp[y] != yy and fp[yy] != x
            assert vp[j] == xx
            assert vp[xx] == y
            assert vp[y] == yy
            vp[j] = vp[yy]
            assert vp[k] == x
            vp[k] = vp[x]
            fp[kk] = k
            fp[ii] = j
        elif fp[yy] == x:
            # face: (x xx i A y k B yy) -> (i A k B)
            # (i A jj) (k B ii)
            assert fp[xx] != y and fp[y] != yy
            assert vp[i] == yy
            vp[i] = vp[yy]
            assert vp[k] == x
            assert vp[x] == xx
            assert vp[xx] == y
            vp[k] = vp[y]
            fp[ii] = i
            fp[jj] = k
        elif fp[y] == yy:
            # face: (x xx i A y yy j C) -> (i A j C)
            # (i A jj) (j C kk)
            assert fp[xx] != y and fp[yy] != x
            assert vp[i] == yy
            assert vp[yy] == x
            vp[i] = vp[x]
            assert vp[j] == xx
            assert vp[xx] == y
            vp[j] = vp[y]
            fp[kk] = i
            fp[jj] = j
        else:
            # face: (x xx i A y k B yy j C) -> (i A j C k B)
            # (i A jj) (j C kk) (k B ii)
            assert fp[xx] != y and fp[y] != yy and fp[yy] != x
            assert vp[i] == yy
            vp[i] = vp[yy]
            assert vp[j] == xx
            assert vp[xx] == y
            vp[j] = vp[y]
            assert vp[k] == x
            vp[k] = vp[x]
            fp[jj] = j
            fp[kk] = k
            fp[ii] = i

        self._n -= 4
        fd[fl[i]] -= 4
        vd[vl[i]] -= 1
        vd[vl[j]] -= 2
        vd[vl[k]] -= 1

        assert perm_check(vp, self._n), vp
        assert perm_check(fp, self._n), fp

    ######################################
    # canonical labels and automorphisms #
    ######################################
    # This is chosen so that augmentation works fast (for the exhaustive generation)
    #   - augment1: trisection (single vertex, single face maps)
    #   - augment2: face split (single vertex)
    #   - augment3: vertex split

    def _good_starts(self, i0=-1):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph

            sage: CM = []
            sage: w = [0,1,2,3,4,5,0,6,7,1,2,5,3,4,6,7]
            sage: CM.append(FatGraph.from_unicellular_word(w))
            sage: vp = '(0,15,3,6,8)(1,14,18,10,9,12,5,19)(2,7,4,17,11)(13,16)'
            sage: ep = '(0,1)(2,3)(4,5)(6,7)(8,9)(10,11)(12,13)(14,15)(16,17)(18,19)'
            sage: fp = '(0,19,14)(1,8,10,17,13,9,6,2,15)(3,11,18,5,7)(4,12,16)'
            sage: CM.append(FatGraph(vp, ep, fp))
            sage: for cm in CM:
            ....:     gs = set(cm._good_starts())
            ....:     assert gs
            ....:     for i in range(cm.num_darts()):
            ....:         ggs = cm._good_starts(i)
            ....:         if i in gs:
            ....:             ggs = cm._good_starts(i)
            ....:             assert ggs and ggs[0] == i and sorted(ggs) == sorted(gs), (gs, ggs, i)
            ....:         else:
            ....:             assert not ggs, (gs, ggs, i)
        """
        n = self._n
        ep = self._ep
        vd = self._vd
        fd = self._fd
        fp = self._fp
        vl = self._vl
        fl = self._fl
        if i0 == -1:
            ans = []
        else:
            ans = [i0]

        if self._nv > 1:
            # consider only edges with distinct start and end.
            # Maximize the degrees of vertices, then whether the
            # adjacent faces are distinct, then the degree of adjacent
            # faces.
            if i0 != -1:
                j0 = ep[i0]
                if vl[i0] == vl[j0]:
                    return None
                best = (vd[vl[i0]], vd[vl[j0]], fl[i0] != fl[j0], fd[fl[i0]], fd[fl[j0]])
            else:
                best = None
            for i in range(self._n):
                if i == i0:
                    continue
                j = ep[i]
                if vl[i] == vl[j]:
                    continue
                cur = (vd[vl[i]], vd[vl[j]], fl[i] != fl[j], fd[fl[i]], fd[fl[j]])
                if best is None:
                    best = cur
                if cur > best:
                    if i0 != -1:
                        return None
                    else:
                        del ans[:]
                        best = cur
                if cur == best:
                    ans.append(i)

        elif self._nf > 1:
            # (we have a single vertex but several faces)
            # consider only edges with distinct faces on their sides.
            # Maximize their degrees.
            if i0 != -1:
                j0 = ep[i0]
                if fl[i0] == fl[j0]:
                    return None
                best = (fd[fl[i0]], fd[fl[j0]])
            else:
                best = None
            for i in range(self._n):
                if i == i0:
                    continue
                j = ep[i]
                if fl[i] == fl[j]:
                    continue
                cur = (fd[fl[i]], fd[fl[j]])
                if best is None:
                    best = cur
                if cur > best:
                    if i0 != -1:
                        return None
                    else:
                        del ans[:]
                        best = cur
                if cur == best:
                    ans.append(i)
        else:
            # (we have a single face and a single vertex)
            # Minimize the face angle between i and ep[i]
            
            # 1. compute the "face angle" between the half edges
            fa = [None] * n
            fa[0] = 0
            i = fp[0]
            j = 1
            while i != 0:
                fa[i] = j
                i = fp[i]
                j += 1

            # 2. first guess
            if i0 != -1:
                j0 = ep[i0]
                best = fa[j0] - fa[i0]
                if best < 0: best += n
            else:
                best = None

            # 3. run accross the edges
            for i in range(self._n):
                if i == i0:
                    continue
                j = ep[i]
                cur = fa[j] - fa[i]
                if cur < 0: cur += n
                if best is None:
                    best = cur
                if cur < best:
                    if i0 != -1:
                        return None
                    else:
                        del ans[:]
                        best = cur
                if cur == best:
                    ans.append(i)

        return ans

    def _canonical_labelling_from(self, i0):
        r"""
        Edges gets relabelled (2i, 2i+1).

        OUTPUT: a triple ``(fc, fd, rel)`` where
        
        - ``fc`` is the list of edges seen along the walk (with respect to the new
          numbering)
        
        - ``fd``: face degrees seen along the walk

        - ``rel``: relabelling map {current labels} -> {canonical labels}

        EXAMPLES::
   
            sage: from surface_dynamics.topology.fat_graph import FatGraph

            sage: vp = '(0,15,3,6,8)(1,14,18,10,9,12,5,19)(2,7,4,17,11)(13,16)'
            sage: ep = '(0,1)(2,3)(4,5)(6,7)(8,9)(10,11)(12,13)(14,15)(16,17)(18,19)'
            sage: fp = '(0,19,14)(1,8,10,17,13,9,6,2,15)(3,11,18,5,7)(4,12,16)'
            sage: cm = FatGraph(vp, ep, fp)
            sage: for i in range(20):
            ....:     fc, fd, rel = cm._canonical_labelling_from(i)
            ....:     assert len(fc) == 20
            ....:     assert sorted(fd, reverse=True) == [9, 5, 3, 3]
            ....:     assert sorted(rel) == list(range(20))
        """
        n = self._n
        ep = self._ep
        fp = self._fp

        fc = []        # faces seen along the walk
        fd = []        # face degrees seen along the walk
        rel = [-1] * n # dart relabeling

        rel[i0] = 0   # first edge is relabelled (0,1)
        fc.append(0)
        c = 2         # current dart number (in the new labelling scheme)

        # walk along faces first, starting from i0
        # along the way, we collect unseen edges
        i = fp[i0]
        wait = deque([ep[i0]])
        cyc = [0]
        d = 1
        while i != i0:
            if rel[i] == -1:
                j = ep[i]
                if rel[j] != -1:
                    assert rel[j] % 2 == 0
                    rel[i] = rel[j] + 1
                else:
                    rel[i] = c
                    c += 2
                    wait.append(j)
            fc.append(rel[i])
            d += 1
            i = fp[i]
        fd.append(d)

        while wait:
            i0 = wait.popleft()
            if rel[i0] != -1:
                continue
            assert rel[ep[i0]] != -1 and rel[ep[i0]] % 2 == 0
            rel[i0] = rel[ep[i0]] + 1
            fc.append(rel[i0])
            i = fp[i0]
            d = 1
            while i != i0:
                if rel[i] == -1:
                    j = ep[i]
                    if rel[j] != -1:
                        assert rel[j] % 2 == 0
                        rel[i] = rel[j] + 1
                    else:
                        rel[i] = c
                        c += 2
                        wait.append(j)
                fc.append(rel[i])
                i = fp[i]
                d += 1
            fd.append(d)

        assert len(fc) == self._n, (fc, fd, rel)
        assert len(fd) == self._nf, (fc, fd, rel)

        return fc, fd, rel

    # TODO: this is a waste! We should implement the proper partial relabeling
    def _canonical_labelling_from_if_better(self, best, i0):
        r"""
        INPUT:

        - ``best`` - a triple ``(fc, fd, rel)`` as given from _canonical_labeling_from

        - ``i0`` - start edge
        """
        fc_best, fd_best, rel_best = best

        n = self._n
        ep = self._ep
        fp = self._fp
        fl = self._fl
        sfd = self._fd

        is_fd_better = 0  # whether the current face degree works better
        is_fc_better = 0  # whether the current relabelling works better
                          #  0 = equal
                          #  1 = better
                          # -1 = worse

        fd = []        # face degrees seen along the walk (want to maximize)
        fc = []        # edges seen along the walk (want to minimize)
        rel = [-1] * n # dart relabeling

        # walk along faces first, starting from i0
        # along the way, we collect unseen edges

        rel[i0] = 0   # first edge is relabelled (0,1)
        c = 2         # current dart number (in the new labelling scheme)

        d = sfd[fl[i0]]
        if d < fd_best[0]:
            return -1, None
        elif d > fd_best[0]:
            is_fd_better = 1
        fd.append(d)

        fc.append(0)
        i = fp[i0]
        wait = deque([ep[i0]])
        cyc = [0]
        while i != i0:
            if rel[i] == -1:
                j = ep[i]
                if rel[j] != -1:
                    assert rel[j] % 2 == 0
                    rel[i] = rel[j] + 1
                else:
                    rel[i] = c
                    c += 2
                    wait.append(j)

            # edge comparison
            ii = rel[i]
            if not is_fd_better and not is_fc_better:
                if ii < fc_best[len(fc)]:
                    is_fc_better = 1
                elif ii > fc_best[len(fc)]:
                    if self._nf == 1:
                        return -1, None
                    is_fc_better = -1
            fc.append(ii)
            i = fp[i]

        while wait:
            i0 = wait.popleft()

            # face already seen?
            if rel[i0] != -1:
                continue

            # face degree comparison
            d = sfd[fl[i0]]
            if not is_fd_better:
                if d < fd_best[len(fd)]:
                    return -1, None
                elif d > fd_best[len(fd)]:
                    is_fd_better = 1
            fd.append(d)

            # root edge comparison
            assert rel[ep[i0]] != -1 and rel[ep[i0]] % 2 == 0
            ii0 = rel[ep[i0]] + 1
            if not is_fd_better and not is_fc_better:
                if ii0 < fc_best[len(fc)]:
                    is_fc_better = 1
                elif ii0 > fc_best[len(fc)]:
                    is_fc_better = -1

            rel[i0] = ii0
            fc.append(ii0)
            i = fp[i0]

            while i != i0:
                # label the i-th edge (if not already)
                if rel[i] == -1:
                    j = ep[i]
                    if rel[j] != -1:
                        assert rel[j] % 2 == 0
                        rel[i] = rel[j] + 1
                    else:
                        rel[i] = c
                        c += 2
                        wait.append(j)

                # edge comparison
                ii = rel[i]
                if not is_fd_better and not is_fc_better:
                    if ii < fc_best[len(fc)]:
                        is_fc_better = 1
                    elif ii > fc_best[len(fc)]:
                        is_fc_better = -1

                # update
                fc.append(ii)
                i = fp[i]

        cur = (fc, fd, rel)

        if is_fd_better:
            return 1, cur
        else:
            return is_fc_better, cur

    def _is_canonical(self, i0):
        r"""
        Return a pair ``(answer, automorphisms)`` where answer is a boolean
        that says whether this map is in canonical form and ``automorphisms`` form
        a generating set of the group of automorphisms.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph

            sage: vp = '(0,15,3,6,8)(1,14,18,10,9,12,5,19)(2,7,4,17,11)(13,16)'
            sage: ep = '(0,1)(2,3)(4,5)(6,7)(8,9)(10,11)(12,13)(14,15)(16,17)(18,19)'
            sage: fp = '(0,19,14)(1,8,10,17,13,9,6,2,15)(3,11,18,5,7)(4,12,16)'
            sage: cm = FatGraph(vp, ep, fp)
            sage: any(cm._is_canonical(i)[0] for i in range(20))
            True

        A genus 1 example with 4 symmetries::

            sage: vp = '(0,8,6,4,3,7)(1,9,11,5,2,10)'
            sage: ep = '(0,1)(2,3)(4,5)(6,7)(8,9)(10,11)'
            sage: fp = '(0,10,9)(2,4,11)(1,7,8)(3,5,6)'
            sage: cm = FatGraph(vp, ep, fp)
            sage: for i in range(12):
            ....:     test, aut_grp = cm._is_canonical(i)
            ....:     if test: print(aut_grp.group_cardinality())
            4
            4
            4
            4
        """
        roots = self._good_starts(i0)
        if roots is None:
            return False, None

        if len(roots) == 1:
            return True, None

        # perform complete relabelling
        P = PermutationGroupOrbit(self._n, [], roots)
        i = next(P)
        assert i == i0
        best = self._canonical_labelling_from(i)
        rel0 = perm_invert(best[2])

        for i in P:
            test, cur = self._canonical_labelling_from_if_better(best, i)

            if test == 1:
                return False, None

            elif test == 0:
                fc, fd, rel = cur
                aut = perm_compose(rel, rel0)
                P.add_generator(aut)

        return True, P

    def automorphism_group(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: from surface_dynamics.misc.permutation import perm_conjugate

        The four unicellular map with 4 edges in genus 2::

            sage: cm0 = FatGraph.from_unicellular_word([0,1,0,1,2,3,2,3])
            sage: cm1 = FatGraph.from_unicellular_word([0,1,0,2,1,3,2,3])
            sage: cm2 = FatGraph.from_unicellular_word([0,1,0,2,3,1,2,3])
            sage: cm3 = FatGraph.from_unicellular_word([0,1,2,3,0,1,2,3])
            sage: for cm in [cm0, cm1, cm2, cm3]:
            ....:     P = cm.automorphism_group()
            ....:     print(P.group_cardinality())
            ....:     vp = cm.vertex_permutation()
            ....:     ep = cm.edge_permutation()
            ....:     fp = cm.face_permutation()
            ....:     for a in P.gens():
            ....:         pp = perm_conjugate(vp, a)
            ....:         assert pp == vp, (vp, pp)
            2
            1
            1
            8

            sage: cm = FatGraph.from_unicellular_word([0,1,2,3,0,4,1,2,3,4])
            sage: cm.automorphism_group().group_cardinality()
            2

        An example with two faces::

            sage: vp = '(0,9,5,6,7,4,8,1,2,3)'
            sage: ep = '(0,2)(1,3)(4,6)(5,7)(8,9)'
            sage: fp = '(0,1,2,3,8)(4,5,6,7,9)'
            sage: cm = FatGraph(vp,ep,fp)
            sage: cm.automorphism_group()
            PermutationGroupOrbit(10, [(0,4)(1,5)(2,6)(3,7)(8,9)])
        """
        roots = self._good_starts()
        P = PermutationGroupOrbit(self._n, [], roots)
        if len(roots) == 1:
            return P
        i0 = next(P)
        best = self._canonical_labelling_from(i0)
        rel0 = perm_invert(best[2])

        for i in P:
            test, cur = self._canonical_labelling_from_if_better(best, i)

            if test == 1:
                rel0 = perm_invert(cur[2])
                best = cur
            elif test == 0:
                fc, fd, rel = cur
                aut = perm_compose(rel, rel0)
                P.add_generator(aut)

        return P

    def relabel(self, r):
        r"""
        Relabel according to the permutation ``p``

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph

            sage: cm = FatGraph.from_unicellular_word([0,1,0,1,2,3,2,3])
            sage: cm.relabel([4,7,0,2,3,5,1,6])
            sage: cm._check()
        """
        n = self._n
        self._vp = perm_conjugate(self._vp, r, n)
        self._ep = perm_conjugate(self._ep, r, n)
        self._fp = perm_conjugate(self._fp, r, n)
        
        vl = [None] * n
        fl = [None] * n
        fa = [None] * n

        for i in range(n):
            j = r[i]
            vl[j] = self._vl[i]
            fl[j] = self._fl[i]

        self._vl = vl
        self._fl = fl
