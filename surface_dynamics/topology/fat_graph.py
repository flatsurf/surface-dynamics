r"""
Fat graph.

"""
# ****************************************************************************
#       Copyright (C) 2019-2023 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import numbers

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from collections import deque
from surface_dynamics.misc.permutation import (perm_compose, perm_conjugate, perm_conjugate_inplace,
                                               perm_dense_cycles, constellation_init,
                                               perm_num_cycles, perm_on_list_inplace, perm_orbit,
                                               perm_from_base64_str, perm_base64_str,
                                               perm_check, perm_cycle_string, perm_hash,
                                               perm_cycles, perm_cycle_type, perm_invert,
                                               perm_is_regular, perm_is_on_list_stabilizer, PermutationGroupOrbit, perms_are_transitive, perms_transitive_components)


###########################
# Miscellaneous functions #
###########################


def num_and_weighted_num(it):
    s = QQ.zero()
    n = ZZ.zero()
    for _, aut in it:
        n += ZZ.one()
        if aut is None:
            s += QQ.one()
        else:
            s += QQ((1, aut.group_cardinality()))
    return n, s


def list_extrems(l, n):
    """
    EXAMPLES::

        sage: from surface_dynamics.topology.fat_graph import list_extrems
        sage: list_extrems([3,5,3,17,7,2,1,19],4)
        (3, 17)
    """
    if not n:
        raise ValueError
    vdmin = vdmax = l[0]
    for i in range(1, n):
        if l[i] > vdmax:
            vdmax = l[i]
        elif l[i] < vdmin:
            vdmin = l[i]
    return (vdmin, vdmax)


class FatGraph(object):
    r"""
    Fat graph or graph embedded in an orientable surface.

    A fat graph is a graph embedded in an orientable surface. It is encoded as
    a pair of permutations ``vp`` (for vertex permutation) and ``fp`` (for face
    permutation). The encoding is done by first associating an integer to each
    oriented edge such that the orientation reversing maps is `0
    \leftrightarrow 1`, `2 \leftrightarrow 3`, etc. Then the vertex permutation
    is the permutation of the oriented edge which corresponds to move
    counterclockwise to the next outgoing edge at each vertex.  The face
    permutation is the permutation of the oriented edge which corresponds to
    follow the arrow inside the face on its left side (so that faces also run
    counterclockwise).

    EXAMPLES:

    The once punctured torus::

        sage: from surface_dynamics import FatGraph

        sage: vp = '(0,2,1,3)'
        sage: fp = '(0,2,1,3)'
        sage: FatGraph(vp, fp)
        FatGraph('(0,2,1,3)', '(0,2,1,3)')

    Actually it is enough to specify one of the vertex permutation or face
    permutation as one determines the other::

        sage: vp = '(0,3,1)(4)(2,5,6,7)'
        sage: fp = '(0,3,7,5,4,2)(1)(6)'
        sage: F0 = FatGraph(vp=vp, fp=fp)
        sage: F1 = FatGraph(fp=fp)
        sage: F2 = FatGraph(vp=vp, fp=fp)
        sage: F3 = FatGraph(vp=vp)
        sage: F0 == F1 and F0 == F2 and F0 == F3
        True
    """
    __slots__ = ['_n',       # number of darts (non-negative integer)
                 '_mutable', # mutability flag
                 # permutations
                 '_vp',  # vertex permutation (array of length _n)
                 '_fp',  # face permutation (array of length _n)
                 # TODO: think whether it is useful to keep any of these...
                 # labels (identify uniquely the vertices)
                 '_vl',  # vertex labels (array of length _n)
                 '_fl',  # face labels (array of length _n)
                 # numbers
                 '_nv',  # number of vertices (non-negative integer)
                 '_nf',  # number of faces (non-negative integer)
                 # degrees
                 '_vd',  # vertex degrees (array of length _nv)
                 '_fd']  # face degrees (array of length _nf)

    def __init__(self, vp=None, fp=None, max_num_dart=None, mutable=False, check=True):
        vp, _, fp = constellation_init(vp, 2, fp)
        self._vp = vp
        self._fp = fp
        if len(vp) != len(fp):
            raise ValueError("invalid permutations")
        self._n = len(vp)   # number of darts
        self._nf = 0        # number of faces

        self._vl, self._vd = perm_dense_cycles(vp, self._n)
        self._nv = len(self._vd)  # number of vertices
        self._fl, self._fd = perm_dense_cycles(fp, self._n)
        self._nf = len(self._fd)  # number of faces

        self._mutable = bool(mutable)

        if max_num_dart is not None:
            if max_num_dart < self._n:
                raise ValueError
            self._realloc(max_num_dart)

        if check:
            self._check()

    @staticmethod
    def from_unicellular_word(w):
        r"""
        Build a fat graph from a word on the letters {0, ..., n-1} where
        each letter appears exactly twice.

        EXAMPLES::

            sage: from surface_dynamics import FatGraph
            sage: FatGraph.from_unicellular_word([0,1,0,2,3,4,1,4,3,2])
            FatGraph('(0,4)(1,3,9,2)(5,6)(7,8)', '(0,2,1,4,6,8,3,9,7,5)')
            sage: FatGraph.from_unicellular_word([0,1,2,0,3,2,4,1,3,4])
            FatGraph('(0,8,4,3,9,6)(1,5,7,2)', '(0,2,4,1,6,5,8,3,7,9)')
        """
        # edge 0 is relabelled (0,1), edge 1 is relabelled (2,3), etc
        n = len(w)
        m = n // 2
        fp = [None] * n
        seen = [0] * m
        previous = 2 * w[-1] + 1
        for k in w:
            if not 0 <= k < n:
                raise ValueError('invalid unicellular word')
            if seen[k] == 0:
                # first time we see an edge it is labelled 0, 2, 4, ...
                fp[previous] = previous = 2*k
                seen[k] = True
            elif seen[k] == 1:
                # twins get labelled 1, 3, 5, ...
                fp[previous] = previous = 2*k + 1
            else:
                raise ValueError('invalid unicellular word')

        # consistency check
        assert previous == 2 * w[-1] + 1

        return FatGraph(fp=fp)

    @staticmethod
    def from_string(s):
        r"""
        Build a fat graph from a serialized string.

        See also :meth:`to_string`.

        EXAMPLES::

            sage: from surface_dynamics import FatGraph
            sage: s = '20_i31027546b98jchedfag_23146758ab9igdhfejc0'
            sage: F = FatGraph.from_string(s)
            sage: F.to_string() == s
            True
            sage: FatGraph.from_string('0__')
            FatGraph('()', '()')
        """
        if not isinstance(s, str) or s.count('_') != 2:
            raise ValueError("invalid input")
        n, vp, fp = s.split('_')
        n = int(n)
        vp = perm_from_base64_str(vp, n)
        fp = perm_from_base64_str(fp, n)
        return FatGraph(vp, fp)

    def _check(self, error=RuntimeError):
        vp = self._vp
        vl = self._vl
        vd = self._vd

        fp = self._fp
        fl = self._fl
        fd = self._fd

        n = self._n
        nf = self._nf
        nv = self._nv

        if any(vp[i] == -1 or fp[i] == -1 for i in range(self._n)):
            raise ValueError('inactive dart not allowed')

        if not perm_check(vp, n):
            raise ValueError("invalid vertex permutation: %s" % vp)
        if not perm_check(fp, n):
            raise ValueError("invalid face permutation: %s" % fp)

        if n and (perm_num_cycles(vp, n) != self._nv):
            raise error("wrong number of vertices")
        if n and (perm_num_cycles(fp, n) != self._nf):
            raise error("wrong number of faces")

        if len(vl) < n or len(fl) < n or len(vd) < nv or len(fd) < nf:
            raise error("inconsistent lengths")

        if any(x < 0 or x > n for x in vd[:nv]):
            raise error("invalid vertex degrees")
        if any(x < 0 or x > n for x in fd[:nf]):
            raise error("invalid face degrees")

        ffd = [0] * nf
        vvd = [0] * nv

        for i in range(n):
            if fp[vp[i] ^ 1] != i:
                raise error("fp[ep[vp[%d]]] = %d" % (i, fp[vp[i] ^ 1]))
            if fl[i] < 0 or fl[i] >= nf:
                raise error("face label out of range: fl[%d] = %d" % (i, fl[i]))
            if vl[i] < 0 or vl[i] >= nv:
                raise error("vertex label out of range: vl[%d] = %d" % (i, vl[i]))
            if fl[fp[i]] != fl[i]:
                raise error("fl[fp[%d]] = %d while fl[%d] = %d" % (i, fl[fp[i]], i, fl[i]))

            if vl[vp[i]] != vl[i]:
                raise error("vl[vp[%d]] = vl[%d] = %d while vl[%d] = %d" % (i, vp[i], vl[vp[i]], i, vl[i]))

            ffd[fl[i]] += 1
            vvd[vl[i]] += 1

        if vvd != vd[:nv]:
            raise error("inconsistent vertex labels/degrees, got %s instead of vd = %s" % (vvd, vd[:nv]))
        if ffd != fd[:nf]:
            raise error("inconsistent face labels/degrees, got %s instead of fd = %s" % (ffd, fd[:nf]))

    def _check_dart(self, i):
        if not isinstance(i, numbers.Integral) or i < 0:
            raise TypeError('dart must be a non-negative integer, got i=%d' % i)
        i = int(i)
        if i >= self._n:
            raise ValueError('dart i=%d out of range' % i)
        return i

    def __hash__(self):
        if self._mutable:
            raise TypeError("mutable FatGraph not hashable")
        return perm_hash(self._vp) ^ (13522761 * perm_hash(self._fp))

    def _realloc(self, max_num_dart):
        if max_num_dart < self._n:
            return
        self._vp.extend([-1] * (max_num_dart - self._n))
        self._fp.extend([-1] * (max_num_dart - self._n))
        self._vl.extend([-1] * (max_num_dart - self._n))
        self._fl.extend([-1] * (max_num_dart - self._n))
        self._vd.extend([-1] * (max_num_dart - self._nv))
        self._fd.extend([-1] * (max_num_dart - self._nf))

    def is_connected(self):
        r"""
        Return whether the graph is connected.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: FatGraph(vp='(0,2)(1,3)').is_connected()
            True
            sage: FatGraph(vp='(0,1)(2,3)').is_connected()
            False
        """
        return perms_are_transitive([self._vp, self._fp], self._n)

    def connected_components(self):
        r"""
        Return the list of connected components.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: FatGraph(vp='(0,2)(1,3)').connected_components()
            [FatGraph('(0,2)(1,3)', '(0,3)(1,2)')]
            sage: FatGraph(vp='(0,1)(2,3)').connected_components()
            [FatGraph('(0,1)', '(0)(1)'), FatGraph('(0,1)', '(0)(1)')]
        """
        ccs = perms_transitive_components([self._vp, self._fp], self._n)
        if len(ccs) == 1 and not self._mutable:
            return [self]

        # build a FatGraph for each connected component
        connected_graphs = []
        for cc in ccs:
            relabel = {j: i for i, j in enumerate(cc)}
            vp = [-1] * len(cc)
            fp = [-1] * len(cc)
            for j in cc:
                vp[relabel[j]] = relabel[self._vp[j]]
                fp[relabel[j]] = relabel[self._fp[j]]
            connected_graphs.append(FatGraph(vp, fp))

        return connected_graphs

    def disjoint_union(self, *args, mutable=False):
        r"""
        Return the union of ``self`` with the graphs provided as arguments.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: fg0 = FatGraph('(0,4,2,5)(1,6)(3,7)')
            sage: fg1 = FatGraph('(0,3,5)(1,2,4)')
            sage: fg2 = FatGraph('(0,1)')
            sage: fg = fg0.disjoint_union(fg1, fg2)
            sage: fg
            FatGraph('(0,4,2,5)(1,6)(3,7)(8,11,13)(9,10,12)(14,15)', '(0,6,3,4,2,7,1,5)(8,12,11,9,13,10)(14)(15)')
            sage: fg.connected_components()
            [FatGraph('(0,4,2,5)(1,6)(3,7)', '(0,6,3,4,2,7,1,5)'),
             FatGraph('(0,3,5)(1,2,4)', '(0,4,3,1,5,2)'),
             FatGraph('(0,1)', '(0)(1)')]
        """
        if not args:
            if not self._mutable and not mutable:
                return self
            else:
                return self.copy(mutable=mutable)

        shifts = [self._n]
        for fg in args:
            shifts.append(shifts[-1] + fg._n)
        n = shifts[-1]
        vp = list(self._vp[:self._n]) + [-1] * (n - self._n)
        fp = list(self._fp[:self._n]) + [-1] * (n - self._n)
        for fg, shift in zip(args, shifts):
            for i in range(fg._n):
                vp[shift + i] = shift + fg._vp[i]
                fp[shift + i] = shift + fg._fp[i]
        return FatGraph(vp, fp)

    def edge_flip(self, i):
        r"""
        Return the half-edge that together with ``i`` makes an edge.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: F = FatGraph('(0,2,1,3)')
            sage: F.edge_flip(0)
            1
            sage: F.edge_flip(1)
            0
            sage: F.edge_flip(4)
            Traceback (most recent call last):
            ...
            ValueError: dart i=4 out of range
        """
        i = self._check_dart(i)
        return int(i) ^ 1

    def darts(self):
        r"""
        Iterator through the list of darts.

        EXAMPLES::

            sage: from surface_dynamics import FatGraph
            sage: list(FatGraph('(0,3,5)(1,4,2)').darts())
            [0, 1, 2, 3, 4, 5]
        """
        return range(self._n)

    def copy(self, mutable=None):
        """
        Return a copy of this fat graph.

        INPUT:

        - ``mutable`` -- (optional) whether the copy must be mutable

        EXAMPLES::

            sage: from surface_dynamics import FatGraph
            sage: F = FatGraph.from_unicellular_word([0,1,0,2,3,4,1,4,3,2])
            sage: G = F.copy()
            sage: G._check()
        """
        if mutable is None:
            mutable = self._mutable
        elif not isinstance(mutable, bool):
            raise TypeError('argument mutable must be a boolean')
        if not mutable and not self._mutable:
            # no need to duplicate immutable graphs
            return self
        F = FatGraph.__new__(FatGraph)
        F._n = self._n
        F._vp = self._vp[:]
        F._fp = self._fp[:]
        F._nf = self._nf
        F._nv = self._nv
        F._vl = self._vl[:]
        F._vd = self._vd[:]
        F._fl = self._fl[:]
        F._fd = self._fd[:]
        F._mutable = mutable
        return F

    __copy__ = copy

    def to_string(self):
        r"""
        Serialization to string.

        EXAMPLES::

            sage: from surface_dynamics import FatGraph
            sage: FatGraph.from_unicellular_word([0,1,0,2,3,4,1,4,3,2]).to_string()
            '10_4319065872_2419608537'
            sage: FatGraph('', '').to_string()
            '0__'
        """
        n = self._n
        return str(n) + "_" + \
            perm_base64_str(self._vp, n) + "_" + \
            perm_base64_str(self._fp, n)

    def _symmetric_group(self):
        r"""
        Underlying symmetric group.
        """
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        return SymmetricGroup([i for i in range(self._n) if self._vp[i] != -1])

    def monodromy_group(self):
        r"""
        Return the group generated by the vertex and face permutations.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph

            sage: r1 = FatGraph('(0)(2)(1,3,4,5)')
            sage: G1 = r1.monodromy_group()
            sage: G1
            Subgroup ...
            sage: G1.is_isomorphic(SymmetricGroup(5))
            True

            sage: r2 = FatGraph('(0)(2)(1,3,4)(5,6,7)')
            sage: G2 = r2.monodromy_group()
            sage: G2
            Subgroup ...
            sage: G2.is_isomorphic(PSL(2,7))
            True
        """
        S = self._symmetric_group()
        v = S([i for i in self._vp if i != -1])
        f = S([i for i in self._fp if i != -1])
        return S.subgroup([v, f])

    def is_face_bipartite(self, certificate=False):
        r"""
        Return whether the faces admit a proper 2-coloring.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import *
            sage: vp = '(0,2,1,3)'
            sage: fp = '(0,2,1,3)'
            sage: F = FatGraph(vp, fp, 6, mutable=True)
            sage: F.is_face_bipartite()
            False
            sage: F.split_face(0, 1)
            sage: F.is_face_bipartite()
            True

            sage: vp = '(0,5,2,1,3,4)'
            sage: fp = '(0,2,1,4)(3,5)'
            sage: FatGraph(vp, fp).is_face_bipartite()
            False

            sage: from surface_dynamics import FatGraphs
            sage: F = FatGraphs(g=1, nf=3, nv=3, vertex_min_degree=3)
            sage: F.cardinality_and_weighted_cardinality(filter=lambda x,a: x.is_face_bipartite())
            (3, 5/3)
        """
        # trivial cases
        if self._nf == 0:
            return (True, []) if certificate else True
        elif self._nf == 1:
            return (False, None) if certificate else False

        n = self._n
        fp = self._fp
        fl = self._fl
        nf = self._nf
        colors = [-1] * nf
        seen = [False] * n
        root = 0
        to_test = perm_orbit(fp, root)
        colors[self._fl[root]] = 1
        while to_test:
            e1 = to_test.pop()
            if seen[e1]:
                continue
            e2 = e1 ^ 1
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
                return (False, None) if certificate else False

            seen[e1] = seen[e2] = True

        return (True, colors) if certificate else True

    def __repr__(self):
        n = self._n
        fd = self._fd[:self._nf]
        vd = self._vd[:self._nv]
        fd.sort(reverse=True)
        vd.sort(reverse=True)
        return "FatGraph('%s', '%s')" % (perm_cycle_string(self._vp, True, n),
                                               perm_cycle_string(self._fp, True, n))

    def __eq__(self, other):
        r"""
        TESTS::

            sage: from surface_dynamics.topology.fat_graph import FatGraph

            sage: vp = '(0,2,1,3)'
            sage: fp = '(0,2,1,3)'
            sage: cm1 = FatGraph(vp, fp)
            sage: cm2 = FatGraph(vp, fp, 100)
            sage: cm1 == cm2
            True
        """
        if type(self) != type(other):
            raise TypeError

        if self._n != other._n or self._nf != other._nf or self._nv != other._nv:
            return False

        # here we ignore the vertex and face labels...
        return all(self._vp[i] == other._vp[i] and
                   self._fp[i] == other._fp[i]
                   for i in range(self._n))

    def __ne__(self, other):
        return not self == other

    def vertex_permutation(self, copy=True):
        r"""
        Return the vertex permutation as a list.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: FatGraph('(0,4,2)(1,5,3)').vertex_permutation()
            [4, 5, 0, 1, 2, 3]
        """
        return self._vp[:self._n] if copy else self._vp

    def face_permutation(self, copy=True):
        r"""
        Return the face permutation as a list.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: FatGraph('(0,4,2)(1,5,3)').face_permutation()
            [3, 2, 5, 4, 1, 0]
        """
        return self._fp[:self._n] if copy else self._fp

    def vertex_degrees(self):
        r"""
        Return the degree of vertices.

        EXAMPLES::

            sage: from surface_dynamics import FatGraph
            sage: FatGraph('(0,3)(1,2,5)(4)').vertex_degrees()
            [2, 3, 1]
            sage: FatGraph('()', '()').vertex_degrees()
            [0]
        """
        return [0] if not self._n else self._vd[:self._nv]

    def vertex_min_degree(self):
        return 0 if not self._n else min(self._vd[i] for i in range(self._nv))

    def vertex_max_degree(self):
        return 0 if not self._n else max(self._vd[i] for i in range(self._nv))

    def face_degrees(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics import FatGraph
            sage: FatGraph('()', '()').face_degrees()
            [0]
        """
        return [0] if not self._n else self._fd[:self._nf]

    def face_min_degree(self):
        return 0 if not self._nf else min(self._fd[i] for i in range(self._nf))

    def face_max_degree(self):
        return 0 if not self._nf else max(self._fd[i] for i in range(self._nf))

    def profile(self):
        r"""
        Return the pair of vertex and face profiles.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: FatGraph('(0,4,2,5,1,3)', '(0,5)(1,3,4,2)').profile()
            ([6], [2, 4])
            sage: FatGraph('()', '()').profile()
            ([0], [0])
        """
        return (self.vertex_degrees(), self.face_degrees())

    def num_vertices(self):
        r"""
        Return the number of vertices.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: FatGraph('(0,4,2,5,1,3)', '(0,5)(1,3,4,2)').num_vertices()
            1
            sage: FatGraph('()', '()').num_vertices()
            1
        """
        return 1 if not self._n else self._nv

    def num_edges(self):
        r"""
        Return the number of edges.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: FatGraph('(0,4,2,5,1,3)', '(0,5)(1,3,4,2)').num_edges()
            3
            sage: FatGraph('()', '()').num_edges()
            0
        """
        return self._n // 2

    def num_faces(self):
        r"""
        Return the number of faces.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: FatGraph('(0,4,2,5,1,3)', '(0,5)(1,3,4,2)').num_faces()
            2
            sage: FatGraph('()', '()').num_faces()
            1
        """
        return 1 if not self._n else self._nf

    def vertices(self):
        r"""
        Return the vertices as a list of lists.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: FatGraph('(0,4,2,5,1,3)', '(0,5)(1,3,4,2)').vertices()
            [[0, 4, 2, 5, 1, 3]]
        """
        return perm_cycles(self._vp, True, self._n)

    def edges(self):
        r"""
        Return the edges.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: FatGraph('(0,4,2,5,1,3)', '(0,5)(1,3,4,2)').edges()
            [[0, 1], [2, 3], [4, 5]]
        """
        return [[i, i ^ 1] for i in range(0, self._n, 2) if self._vp[i] != -1]

    def faces(self):
        r"""
        Return the faces.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: FatGraph('(0,4,2,5,1,3)', '(0,5)(1,3,4,2)').faces()
            [[0, 5], [1, 3, 4, 2]]
        """
        return perm_cycles(self._fp, True, self._n)

    def is_vertex_regular(self, d=None):
        r"""
        Return whether all vertices have the same degree.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: FatGraph('(0,5,2,3)(1,6,7,4)').is_vertex_regular()
            True
        """
        if d is not None:
            if not isinstance(d, numbers.Integral):
                raise TypeError('d must be integral')
            if d < 0:
                raise ValueError('d must be non-negative')
        return perm_is_regular(self._vp, d, self._n)

    def is_face_regular(self, d=None):
        r"""
        Return whether all faces have the same degree.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: FatGraph(fp='(0,5,2,3)(1,6,7,4)').is_face_regular()
            True
        """
        if d is not None:
            if not isinstance(d, numbers.Integral):
                raise TypeError('d must be integral')
            if d < 0:
                raise ValueError('d must be non-negative')
        return perm_is_regular(self._fp, d, self._n)

    def is_cubic(self):
        r"""
        Return whether all vertices have degree 3.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: FatGraph('(0,4,1)(2,3,5)').is_cubic()
            True
            sage: FatGraph('(0,4,1)(2)(3)(5)').is_cubic()
            False
        """
        return self.is_vertex_regular(3)

    def is_triangulation(self):
        r"""
        Return whether the underlying graph is a triangulation.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: FatGraph(fp='(0,4,1)(2,3,5)').is_triangulation()
            True
            sage: FatGraph(fp='(0,4,1)(2)(3)(5)').is_triangulation()
            False
        """
        return self.is_face_regular(3)

    def genus(self):
        r"""
        Return the genus of this graph.

        If the graph is not connected, then the answer is the sum of the
        genera of the components.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: vp0 = '(0,4,2,5,1,3)'
            sage: fp0 = '(0,5)(1,3,4,2)'
            sage: cm0 = FatGraph(vp0, fp0)
            sage: cm0.genus()
            1

            sage: vp1 = '(0,6,5,7,2,1,3,4)'
            sage: fp1 = '(0,2,1,4,6,5,3,7)'
            sage: cm1 = FatGraph(vp1, fp1)
            sage: cm1.genus()
            2

        Non-connected examples::

            sage: cm0.disjoint_union(cm0, cm1).genus()
            4
            sage: cm1.disjoint_union(cm0).genus()
            3
        """
        nc = len(self.connected_components())
        chi = self._nf - self._n // 2 + self._nv
        return nc - chi // 2

    def euler_characteristic(self):
        r"""
        Return the Euler characteristic of the associated surface.

        The *Euler characteristic* of a surface complex is `v - e + f`, where
        `v` is the number of vertices, `e` the number of edges and `f` the
        number of faces.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: vp = '(0,4,2,5,1,3)'
            sage: fp = '(0,5)(1,3,4,2)'
            sage: cm = FatGraph(vp, fp)
            sage: cm.euler_characteristic()
            0

        A non-connected example (Euler characteristic is additive)::

            sage: cm.disjoint_union(cm, cm, cm).euler_characteristic()
            0
        """
        return self._nf - self._n // 2 + self._nv

    # TODO: fix your mind about the meaning of dual!!!
    def dual(self):
        r"""
        Return the dual fat graph.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: F = FatGraph(fp='(0)(1)')
            sage: F.dual()
            sage: F
            FatGraph('(0)(1)', '(0,1)')
            sage: F._check()
            sage: s = '20_i31027546b98jchedfag_23146758ab9igdhfejc0'
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

    def edge_lengths_polytope(self, b, min_length=0):
        r"""
        Return the polytope of edge lengths where the input ``b`` specifies the
        length of faces.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: vp = '(0,4,2,5,1,3)'
            sage: fp = '(0,5)(1,3,4,2)'
            sage: cm = FatGraph(vp, fp)

            sage: cm.edge_lengths_polytope([3, 5], min_length=1).vertices_list()
            [[1, 1, 1, 1, 2, 2], [2, 2, 1, 1, 1, 1]]
            sage: cm.edge_lengths_polytope([3, 5], min_length=0).vertices_list()
            [[0, 0, 1, 1, 3, 3], [3, 3, 1, 1, 0, 0]]

            sage: vp = '(0,3,4)(1,2,6)(5)(7)'
            sage: fp = '(0,6,7,2)(1,4,5,3)'
            sage: cm = FatGraph(vp, fp)
            sage: cm.edge_lengths_polytope([5, 7])
            A 2-dimensional polyhedron in QQ^8 defined as the convex hull of 3 vertices

        TESTS::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: vp = '(0,4,2,5,1,3)'
            sage: fp = '(0,5)(1,3,4,2)'
            sage: cm = FatGraph(vp, fp)

            sage: cm.edge_lengths_polytope([3, 5, 2])
            Traceback (most recent call last):
            ...
            ValueError: the length of b must be the number of faces
        """
        if self.num_faces() != len(b):
            raise ValueError("the length of b must be the number of faces")

        ieqs = []
        # positivity
        for i in range(self._n):
            l = [0] * self._n
            l[i] = 1
            ieqs.append([-min_length] + l)

        eqns = []
        # half edge equations
        for i in range(0, self._n, 2):
            j = i ^ 1
            l = [0] * self._n
            l[i] = 1
            l[j] = -1
            eqns.append([0] + l)

        # face equations
        for bb, f in zip(b, self.faces()):
            l = [0] * self._n
            for i in f:
                l[i] = 1
            eqns.append([-bb] + l)

        from sage.geometry.polyhedron.constructor import Polyhedron
        return Polyhedron(ieqs=ieqs, eqns=eqns)

    def integral_points(self, b, min_length=1):
        r"""
        Return the edge lengths solution to the face lengths constraint ``b``.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: vp = '(0,4,2,5,1,3)'
            sage: fp = '(0,5)(1,3,4,2)'
            sage: cm = FatGraph(vp, fp)
            sage: cm.integral_points((2, 4))
            ((1, 1, 1, 1, 1, 1),)
            sage: cm.integral_points((3, 5))
            ((1, 1, 1, 1, 2, 2), (2, 2, 1, 1, 1, 1))
            sage: cm.integral_points((5, 11))
            ((1, 1, 3, 3, 4, 4),
             (2, 2, 3, 3, 3, 3),
             (3, 3, 3, 3, 2, 2),
             (4, 4, 3, 3, 1, 1))
        """
        return self.edge_lengths_polytope(b, min_length).integral_points()

#    def kontsevich_volume_rational_function(self, R=None):
#        r"""
#        This is not under an appropriate form...
#        """
#        from sage.rings.rational_field import QQ
#        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
#
#        nf = self._nf
#        fl = self._fl
#        n = self._n
#        if R is None:
#            R = PolynomialRing(QQ, 'b', nf)
#        gens = R.gens()
#        res = R.one()
#        for i in range(self._n):
#            j = i ^ 1
#            if j < i:
#                continue
#            res *= 1 / self.automorphism_group().group_cardinality() / (gens[fl[i]] + gens[fl[j]])
#        return res

    ##############################
    # Augmentation and reduction #
    ##############################

    def _check_alloc(self, n, nv, nf):
        if len(self._vp) < n or \
           len(self._fp) < n or \
           len(self._vl) < n or \
           len(self._fl) < n or \
           len(self._fd) < nf or \
           len(self._vd) < nv:
            raise TypeError("reallocation needed")

    def _set_genus0_loop(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: f = FatGraph('()', '()', mutable=True)
            sage: f._realloc(2)
            sage: f._set_genus0_loop()
            sage: f._check()
            sage: f
            FatGraph('(0,1)', '(0)(1)')
            sage: f.remove_edge(0)
            sage: f
            FatGraph('()', '()')
        """
        self._check_alloc(2, 1, 2)
        self._n = 2
        self._nf = 2
        self._nv = 1
        # set face
        self._fp[0] = 0
        self._fp[1] = 1
        self._fl[0] = 0
        self._fl[1] = 1
        self._fd[0] = self._fd[1] = 1
        # set vertex
        self._vp[0] = 1
        self._vp[1] = 0
        self._vl[0] = self._vl[1] = 0
        self._vd[0] = 2

    def _set_genus0_edge(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: f = FatGraph('()', '()', mutable=True)
            sage: f._realloc(2)
            sage: f._set_genus0_edge()
            sage: f._check()
            sage: f
            FatGraph('(0)(1)', '(0,1)')
            sage: f.contract_edge(0)
            sage: f
            FatGraph('()', '()')
        """
        if not self._mutable:
            raise ValueError('immutable graph; use a copy instead')
        self._check_alloc(2, 2, 1)
        self._n = 2
        self._nf = 1
        self._nv = 2
        # set vertex
        self._vp[0] = 0
        self._vp[1] = 1
        self._vl[0] = 0
        self._vl[1] = 1
        self._vd[0] = self._vd[1] = 1
        # set face
        self._fp[0] = 1
        self._fp[1] = 0
        self._fl[0] = self._fl[1] = 0
        self._fd[0] = 2

    def _set_genus1_square(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: f = FatGraph('()', '()', mutable=True)
            sage: f._realloc(4)
            sage: f._set_genus1_square()
            sage: f._check()
            sage: f
            FatGraph('(0,2,1,3)', '(0,2,1,3)')
            sage: f.remove_face_trisection(0)
            sage: f
            FatGraph('()', '()')
            sage: f._check()
        """
        if not self._mutable:
            raise ValueError('immutable graph; use a copy instead')
        self._check_alloc(4, 1, 1)
        self._n = 4
        self._nv = self._nf = 1
        vp = self._vp
        fp = self._fp
        vl = self._vl
        fl = self._fl
        vd = self._vd
        fd = self._fd
        fp[0] = 2; fp[1] = 3; fp[2] = 1; fp[3] = 0
        vp[0] = 2; vp[1] = 3; vp[2] = 1; vp[3] = 0
        fl[0] = fl[1] = fl[2] = fl[3] = 0
        fd[0] = 4
        vl[0] = vl[1] = vl[2] = vl[3] = 0
        vd[0] = 4

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
            sage: fp = '(0,2,1,3)'

            sage: vp20 = '(0,4,2,5,1,3)'
            sage: fp20 = '(0,5)(1,3,4,2)'
            sage: cm = FatGraph(vp, fp, 6, mutable=True)
            sage: cm.split_face(2, 0)
            sage: cm == FatGraph(vp20, fp20)
            True

            sage: vp10 = '(0,4,2,1,5,3)'
            sage: fp10 = '(0,2,5)(1,3,4)'
            sage: cm = FatGraph(vp, fp, 6, mutable=True)
            sage: cm.split_face(1, 0)
            sage: cm == FatGraph(vp10, fp10)
            True

            sage: vp30 = '(0,4,2,1,3,5)'
            sage: fp30 = '(0,2,1,5)(3,4)'
            sage: cm = FatGraph(vp, fp, 6, mutable=True)
            sage: cm.split_face(3, 0)
            sage: cm == FatGraph(vp30, fp30)
            True

            sage: vp00 = '(0,5,4,2,1,3)'
            sage: fp00 = '(0,2,1,3,4)(5)'
            sage: cm = FatGraph(vp, fp, 6, mutable=True)
            sage: cm.split_face(0, 0)
            sage: cm == FatGraph(vp00, fp00)
            True

            sage: vp22 = '(0,2,5,4,1,3)'
            sage: fp22 = '(0,4,2,1,3)(5)'
            sage: cm = FatGraph(vp, fp, 6, mutable=True)
            sage: cm.split_face(2, 2)
            sage: cm == FatGraph(vp22, fp22)
            True

        A genus 2 surface::

            sage: vp = '(0,3,6,8)(1,10,9,12,5)(2,7,4,11)(13)'
            sage: fp = '(0,5,7,3,11,1,8,10,4,12,13,9,6,2)'
            sage: cm = FatGraph(vp, fp, 21, mutable=True); cm._check()
            sage: cm.split_face(0, 1); cm._check()
            sage: cm.split_face(4, 13); cm._check()
            sage: cm.split_face(5, 14); cm._check()
            sage: cm
            FatGraph('(0,15,3,6,8)(1,14,18,10,9,12,5,19)(2,7,4,17,11)(13,16)', '(0,19,14)(1,8,10,17,13,9,6,2,15)(3,11,18,5,7)(4,12,16)')
            sage: cm.remove_edge(18); cm._check()
            sage: cm.remove_edge(16); cm._check()
            sage: cm.remove_edge(14); cm._check()
            sage: cm == FatGraph(vp, fp, 21)
            True
        """
        if not self._mutable:
            raise ValueError('immutable graph; use a mutable copy instead')

        i = self._check_dart(i)
        j = self._check_dart(j)

        vp = self._vp
        vl = self._vl
        vd = self._vd

        fp = self._fp
        fl = self._fl
        fd = self._fd

        n = self._n
        nf = self._nf
        nv = self._nv

        self._check_alloc(n + 2, nv, nf + 1)

        x = self._n
        y = self._n + 1
        ii = vp[i] ^ 1  # = fp^-1(i)
        jj = vp[j] ^ 1  # = fp^-1(j)

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
            sage: fp = '(0,2,1,3)'
            sage: cm = FatGraph(vp, fp)

            sage: vp20 = '(0,5,4,2,1,3)'
            sage: fp20 = '(0,2,1,3,4)(5)'
            sage: cm2 = FatGraph(vp20, fp20, 6, mutable=True)
            sage: cm2.remove_edge(4)
            sage: cm2 == cm
            True
            sage: cm2 = FatGraph(vp20, fp20, 6, mutable=True)
            sage: cm2.remove_edge(5)
            sage: cm2 == cm
            True

            sage: vp10 = '(0,4,2,5,1,3)'
            sage: fp10 = '(0,5)(1,3,4,2)'
            sage: cm2 = FatGraph(vp10, fp10, mutable=True)
            sage: cm2.remove_edge(4)
            sage: cm2 == cm
            True
            sage: cm2 = FatGraph(vp10, fp10, mutable=True)
            sage: cm2.remove_edge(5)
            sage: cm2 == cm
            True

            sage: vp30 = '(0,4,2,1,5,3)'
            sage: fp30 = '(0,2,5)(1,3,4)'
            sage: cm2 = FatGraph(vp30, fp30, mutable=True)
            sage: cm2.remove_edge(4)
            sage: cm2 == cm
            True
            sage: cm2 = FatGraph(vp30, fp30, mutable=True)
            sage: cm2.remove_edge(5)
            sage: cm2 == cm
            True

            sage: vp00 = '(0,5,4,2,1,3)'
            sage: fp00 = '(0,2,1,3,4)(5)'
            sage: cm2 = FatGraph(vp00, fp00, mutable=True)
            sage: cm2.remove_edge(4)
            sage: cm2 == cm
            True

            sage: vp22 = '(0,2,5,4,1,3)'
            sage: fp22 = '(0,4,2,1,3)(5)'
            sage: cm2 = FatGraph(vp00, fp00, mutable=True)
            sage: cm2.remove_edge(4)
            sage: cm2 == cm
            True
        """
        if not self._mutable:
            raise ValueError('immutable graph; use a mutable copy instead')

        i = self._check_dart(i)
        j = i ^ 1

        vp = self._vp
        fp = self._fp

        vl = self._vl
        fl = self._fl

        vd = self._vd
        fd = self._fd

        n = self._n
        nf = self._nf

        fi = fl[i]
        fj = fl[j]
        if fi == fj:
            raise ValueError("i=%d and j=%d on the same face" % (i, j))
        fmin = min(fi, fj)
        if i < n - 2 or j < n - 2 or max(fi, fj) != nf - 1:
            raise NotImplementedError

        ii = vp[i] ^ 1
        jj = vp[j] ^ 1
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
        self._vp[self._n] = self._vp[self._n + 1] = -1
        self._fp[self._n] = self._fp[self._n + 1] = -1

    def split_vertex(self, i, j):
        r"""
        Insert a new edge to split the vertex located at the darts i and j.

        This operation keeps the genus constant. The inverse operation is implemented
        in :meth:`contract_edge`.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph

            sage: vp = '(0,2,1,3)'
            sage: fp = '(0,2,1,3)'

            sage: vp02 = '(0,4,1,3)(2,5)'
            sage: fp02 = '(0,4,2,1,3,5)'
            sage: cm = FatGraph(vp, fp, 6, mutable=True)
            sage: cm.split_vertex(0,2)
            sage: cm == FatGraph(vp02, fp02)
            True

            sage: vp01 = '(0,4,3)(1,5,2)'
            sage: fp01 = '(0,2,4,1,3,5)'
            sage: cm = FatGraph(vp, fp, 6, mutable=True)
            sage: cm.split_vertex(0,1)
            sage: cm == FatGraph(vp01, fp01)
            True

            sage: vp03 = '(0,4)(1,3,5,2)'
            sage: fp03 = '(0,2,1,4,3,5)'
            sage: cm = FatGraph(vp, fp, 6, mutable=True)
            sage: cm.split_vertex(0,3)
            sage: cm == FatGraph(vp03, fp03)
            True
        """
        if not self._mutable:
            raise ValueError('immutable graph; use a mutable copy instead')

        i = self._check_dart(i)
        j = self._check_dart(j)

        vp = self._vp
        fp = self._fp
        vl = self._vl
        fl = self._fl
        n = self._n
        nf = self._nf
        nv = self._nv
        vd = self._vd
        fd = self._fd

        self._check_alloc(n + 2, nv + 1, nf)

        x = self._n
        y = self._n + 1
        ii = vp[i]
        jj = vp[j]
        self._n += 2
        self._nv += 1

        if i == j:
            # introduce a vertex of degree 1
            vp[x] = ii
            vp[i] = x
            vp[y] = y

            fp[y] = i
            fp[x] = y
            fp[ii ^ 1] = x

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

            fp[jj ^ 1] = x
            fp[x] = j
            fp[ii ^ 1] = y
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
            sage: fp = '(0,2,1,3)'

            sage: vp02 = '(0,4,1,3)(2,5)'
            sage: fp02 = '(0,4,2,1,3,5)'
            sage: cm = FatGraph(vp02, fp02, mutable=True)
            sage: cm.contract_edge(4)
            sage: cm == FatGraph(vp, fp)
            True
            sage: cm = FatGraph(vp02, fp02, mutable=True)
            sage: cm.contract_edge(5)
            sage: cm == FatGraph(vp, fp)
            True

            sage: vp01 = '(0,4,3)(1,5,2)'
            sage: fp01 = '(0,2,4,1,3,5)'
            sage: cm = FatGraph(vp01, fp01, mutable=True)
            sage: cm.contract_edge(4)
            sage: cm == FatGraph(vp, fp)
            True
            sage: cm = FatGraph(vp01, fp01, mutable=True)
            sage: cm.contract_edge(5)
            sage: cm == FatGraph(vp, fp)
            True

            sage: vp03 = '(0,4)(1,3,5,2)'
            sage: fp03 = '(0,2,1,4,3,5)'
            sage: cm = FatGraph(vp03, fp03, mutable=True)
            sage: cm.contract_edge(4)
            sage: cm == FatGraph(vp, fp)
            True
            sage: cm = FatGraph(vp03, fp03, mutable=True)
            sage: cm.contract_edge(5)
            sage: cm == FatGraph(vp, fp)
            True

        Degree 1 vertices::

            sage: cm = FatGraph('(0,2)(1)(3)', '(0,1,2,3)', mutable=True)
            sage: cm.contract_edge(2)
            sage: cm
            FatGraph('(0)(1)', '(0,1)')
            sage: cm2 = FatGraph('(0,2)(1)(3)', '(0,1,2,3)', mutable=True)
            sage: cm2.contract_edge(3)
            sage: cm == cm2
            True
        """
        if not self._mutable:
            raise ValueError('immutable graph; use a mutable copy instead')

        i = self._check_dart(i)
        j = i ^ 1

        vp = self._vp
        fp = self._fp

        vl = self._vl
        fl = self._fl

        vd = self._vd
        fd = self._fd

        n = self._n
        nv = self._nv

        if vl[i] == vl[j]:
            raise ValueError("i=%d and j=%d on the same vertex" % (i, j))

        vi = vl[i]
        vj = vl[j]
        vmin = min(vi, vj)
        if i < n - 2 or j < n - 2 or max(vi, vj) != nv - 1:
            raise NotImplementedError

        ii = vp[i] ^ 1
        jj = vp[j] ^ 1
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

    def relocalise(self,e,i): # e is the dart that will be relocated, ep[e] will not and i is the edge after the coin where we put e.
        vp = self._vp
        ep = self._ep
        fp = self._fp

        e2 = ep[e]
        j = ep[vp[i]]
        next_e = fp[e2]
        pre_e = ep[vp[e]]

        vp[next_e] = ep[pre_e]
        vp[i] = e
        vp[e] = ep[j]

        fp[pre_e] = next_e
        fp[j] = e
        fp[e2] = i

        # if necessary, we need to update vertex and face degree

    def disconnect(self,e): #disconnect the edge e from the vertex it comes
        vp = self._vp
        ep = self._ep
        fp = self._fp

        e2 = ep[e]
        pre_e = ep[vp[e]]
        next_e = fp[e2]
        next_e2 = fp[e]

        vp[next_e] = ep[pre_e]
        vp[e] = e

        fp[pre_e] = next_e
        fp[e2] = e

        self._nf += -1
        self._nv += 1

    def split_coin(self,e): #split the coin before e by puting and edge and a vertex, return the index of the dart that cut the coin
        vp = self._vp
        ep = self._ep
        fp = self._fp

        x = self._n
        y = self._n + 1
        self._n += 2
        ep.append(y)
        ep.append(x)
        pre_e = ep[vp[e]]

        vp[e] = x
        vp.append(ep[pre_e])
        vp.append(y)

        fp[pre_e] = x
        fp.append(y)
        fp.append(e)

        self._nv += 1
    
        return x

    def trisect_face(self, i, j, k):
        r"""
        Insert a bridge

        INPUT:

        - ``i``, ``j``, ``k`` - dart in the same face in counter-clockwise
          order


        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph

            sage: vp = '(0,2,1,3)'
            sage: fp = '(0,2,1,3)'
            sage: cm = FatGraph(vp, fp, 8, mutable=True)

            sage: vp021 = '(0,7,2,6,5,1,4,3)'
            sage: fp021 = '(0,5,1,3,7,2,4,6)'
            sage: cm021 = FatGraph(vp021, fp021)
            sage: cm.trisect_face(0, 2, 1)
            sage: cm == cm021
            True

            sage: cm = FatGraph(vp, fp, 10, mutable=True)
            sage: cm.trisect_face(0, 0, 3)

            sage: cm = FatGraph(vp, fp, 10, mutable=True)
            sage: cm.trisect_face(0, 3, 3)

            sage: cm = FatGraph(vp, fp, 10, mutable=True)
            sage: cm.trisect_face(0, 3, 0)

            sage: cm = FatGraph(vp, fp, 10, mutable=True)
            sage: cm.trisect_face(0, 0, 0)
        """
        if not self._mutable:
            raise ValueError('immutable graph; use a mutable copy instead')

        vp = self._vp
        fp = self._fp

        vl = self._vl
        fl = self._fl

        vd = self._vd
        fd = self._fd

        n = self._n
        nf = self._nf
        nv = self._nv

        i = self._check_dart(i)
        j = self._check_dart(j)
        k = self._check_dart(k)

        if fl[i] != fl[j] or fl[i] != fl[k]:
            raise ValueError("darts in distinct faces")

        self._check_alloc(n + 4, nv, nf)
        self._n += 4

        ii = vp[i] ^ 1 # = fp^-1(i) at the end of B
        jj = vp[j] ^ 1 # = fp^-1(j) at the end of A
        kk = vp[k] ^ 1 # = fp^-1(k) at the end of C

        x = n
        y = n + 1
        xx = n + 2
        yy = n + 3

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
        r"""
        Remove the face trisection at ``x``.

        TESTS::

            sage: from surface_dynamics.topology.fat_graph import FatGraph

            sage: vp = '(0,2,1,3)'
            sage: fp = '(0,2,1,3)'
            sage: for i, j, k in [(0, 2, 1), (0, 0, 3), (0, 3, 3), (0, 3, 0), (0, 0, 0)]:
            ....:     cm = FatGraph(vp, fp, 8, mutable=True)
            ....:     cm.trisect_face(i, j, k)
            ....:     cm.remove_face_trisection(4)
            ....:     assert cm == FatGraph(vp, fp)
        """
        fp = self._fp
        fl = self._fl
        fd = self._fd
        vp = self._vp
        vl = self._vl
        vd = self._vd

        x = self._check_dart(x)
        xx = fp[x]
        y = x ^ 1
        yy = xx ^ 1

        if fl[x] != fl[y] or fl[x] != fl[xx] or fl[x] != fl[yy]:
            raise ValueError("not a trisection")

        if (x, y, xx, yy) != (self._n - 4, self._n - 3, self._n - 2, self._n - 1):
            raise NotImplementedError('x={} y={} xx={} yy={}'.format(x, y, xx, yy))

        # face: (x xx i A y k B yy j C) -> (i A j C k B)
        #    -> (i A jj) (j C kk) (k B ii)
        i = fp[xx]
        k = fp[y]
        j = fp[yy]
        ii = vp[yy] ^ 1  # = fp^-1(yy)
        jj = vp[y] ^ 1  # = fp^-1(y)
        kk = vp[x] ^ 1  # = fp^-1(x)

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
        vp[self._n] = vp[self._n + 1] = vp[self._n + 2] = vp[self._n + 3] = -1
        vl[self._n] = vl[self._n + 1] = vl[self._n + 2] = vl[self._n + 3] = -1
        fp[self._n] = fp[self._n + 1] = fp[self._n + 2] = fp[self._n + 3] = -1
        fl[self._n] = fl[self._n + 1] = fl[self._n + 2] = fl[self._n + 3] = -1

        assert perm_check(vp, self._n), vp
        assert perm_check(fp, self._n), fp

    def rotate(self, e):
        r"""
        Rotate the edge ``e`` counterclockwise.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: fg = FatGraph("(0,5,3)(1,4,2)", "(0,2,5,1,3,4)", mutable=True)

            sage: fg.rotate(0)
            sage: fg
            FatGraph('(0,5,1,3)(2,4)', '(0,5,2,1,3,4)')
            sage: fg.rotate(5)
            sage: fg
            FatGraph('(0,5,1,4,3)(2)', '(0,5,1,3,2,4)')

            sage: vp = "(0,8,15)(1,16,12)(2,13,3)(4,11,5)(6,9,7)(10,17,14)"
            sage: fp = "(0,12,2,13,16,10,4,11,14,8,6,9)(1,15,17)(3)(5)(7)"
            sage: fg = FatGraph(vp, fp, mutable=True)
            sage: fg.rotate(0)
            sage: fg.rotate(1)
            sage: fg
            FatGraph('(0,14,10,17)(1,13,3,2)(4,11,5)(6,9,7)(8,15)(12,16)', '(0,2,13,16,10,4,11,14,8,6,9,15)(1,17,12)(3)(5)(7)')
        """
        if not self._mutable:
            raise ValueError('immutable graph; use a mutable copy instead')

        # <----- <----- <-----        <----- <----- <----
        #     c      b ^    a            c  ^     b     a
        #             /                     |
        #            /E                     |E
        #           /                       |
        #          /          -->           |
        #        e/                        e|
        #        /                          |
        # -----> ----->                -----> ----->
        vp = self._vp
        fp = self._fp
        vl = self._vl
        fl = self._fl
        vd = self._vd
        fd = self._fd

        e = self._check_dart(e)
        E = e ^ 1
        if vp[E] == E:
            raise ValueError("can not rotate the given edge")

        b = fp[e]
        B = b ^ 1
        c = fp[b]
        a = vp[E] ^ 1
        A = a ^ 1

        fp[e] = c
        fp[a] = b
        fp[b] = E

        vp[c] = E
        vp[E] = B
        vp[b] = A

        vl[E] = vl[c]
        fl[b] = fl[a]
        vd[vl[b]] -= 1
        vd[vl[c]] += 1
        fd[fl[e]] -= 1
        fd[fl[E]] += 1
        self._check()

    def rotate_dual(self, e):
        r"""
        Rotate the edge ``e`` counterclockwise in the dual graph.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: fg = FatGraph("(0,5,3)(1,4,2)", "(0,2,5,1,3,4)")

            sage: vp = "(0,8,15)(1,16,12)(2,13,3)(4,11,5)(6,9,7)(10,17,14)"
            sage: fp = "(0,12,2,13,16,10,4,11,14,8,6,9)(1,15,17)(3)(5)(7)"
            sage: fg = FatGraph(vp, fp, mutable=True)
            sage: fg.rotate_dual(0)
            sage: fg.rotate_dual(1)
            sage: fg
            FatGraph('(0,15,16)(1,12,8)(2,13,3)(4,11,5)(6,9,7)(10,17,14)', '(0,8,6,9,12,2,13,1,16,10,4,11,14)(3)(5)(7)(15,17)')
        """
        if not self._mutable:
            raise ValueError('immutable graph; use a mutable copy instead')

        #        |b       |c                 \b      |c
        #        |        |                    \     |
        #       B|       C|                     B\   |
        #  a     | e      |           a       e    \ |
        # ------ o ------ o    ====>  ----- o ------ o
        #     A        E  |               A        E |
        #                 |                          |
        e = self._check_dart(e)
        B = self._vp[e]
        b = B ^ 1
        return self.rotate(b)

    ######################################
    # canonical labels and automorphisms #
    ######################################

    # NOTE: this is broken for non-connected surfaces!
    def _good_starts(self, i0=-1, fix_vertices=False, fix_faces=False):
        r"""
        Return a set of half edges invariant under the (yet unknown)
        automorphism group.

        This function runs together with the augmentation algorithm to
        enumerate fat graphs up to isomorphism and must not be changed
        alone. It runs in O(n) as it goes once through every half edge.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph

            sage: CM = []
            sage: w = [0,1,2,3,4,5,0,6,7,1,2,5,3,4,6,7]
            sage: CM.append(FatGraph.from_unicellular_word(w))
            sage: vp = '(0,15,3,6,8)(1,14,18,10,9,12,5,19)(2,7,4,17,11)(13,16)'
            sage: fp = '(0,19,14)(1,8,10,17,13,9,6,2,15)(3,11,18,5,7)(4,12,16)'
            sage: CM.append(FatGraph(vp, fp))
            sage: for cm in CM:
            ....:     gs = set(cm._good_starts())
            ....:     assert gs
            ....:     for i in cm.darts():
            ....:         ggs = cm._good_starts(i)
            ....:         if i in gs:
            ....:             ggs = cm._good_starts(i)
            ....:             assert ggs and ggs[0] == i and sorted(ggs) == sorted(gs), (gs, ggs, i)
            ....:         else:
            ....:             assert not ggs, (gs, ggs, i)
        """
        n = self._n
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
            # Maximize the (degrees of) vertices, then whether the
            # adjacent faces are distinct, then the (degree of) adjacent
            # faces.
            if i0 != -1:
                j0 = i0 ^ 1
                if vl[i0] == vl[j0]:
                    return None
                vi0 = vl[i0]
                vj0 = vl[j0]
                fi0 = fl[i0]
                fj0 = fl[j0]
                if not fix_vertices:
                    vi0 = vd[vi0]
                    vj0 = vd[vj0]
                if not fix_faces:
                    fi0 = fd[fi0]
                    fj0 = fd[fj0]
                best = (vi0, vj0, fl[i0] != fl[j0], fi0, fj0)
            else:
                best = None
            for i in range(self._n):
                if i == i0:
                    continue
                j = i ^ 1
                if vl[i] == vl[j]:
                    continue
                vi = vl[i]
                vj = vl[j]
                fi = fl[i]
                fj = fl[j]
                if not fix_vertices:
                    vi = vd[vi]
                    vj = vd[vj]
                if not fix_faces:
                    fi = fd[fi]
                    fj = fd[fj]
                cur = (vi, vj, fl[i] != fl[j], fi, fj)
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
                j0 = i0 ^ 1
                if fl[i0] == fl[j0]:
                    return None
                fi0 = fl[i0]
                fj0 = fl[j0]
                if not fix_faces:
                    fi0 = fd[fi0]
                    fj0 = fd[fj0]
                best = (fi0, fj0)
            else:
                best = None
            for i in range(self._n):
                if i == i0:
                    continue
                j = i ^ 1
                if fl[i] == fl[j]:
                    continue
                fi = fl[i]
                fj = fl[j]
                if not fix_faces:
                    fi = fd[fi]
                    fj = fd[fj]
                cur = (fi, fj)
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
            # Minimize the face angle between i and i ^ 1

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
                j0 = i0 ^ 1
                best = fa[j0] - fa[i0]
                if best < 0:
                    best += n
            else:
                best = None

            # 3. run across the edges
            for i in range(self._n):
                if i == i0:
                    continue
                j = i ^ 1
                cur = fa[j] - fa[i]
                if cur < 0:
                    cur += n
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
            sage: fp = '(0,19,14)(1,8,10,17,13,9,6,2,15)(3,11,18,5,7)(4,12,16)'
            sage: cm = FatGraph(vp, fp)
            sage: for i in range(20):
            ....:     fc, fd, rel = cm._canonical_labelling_from(i)
            ....:     assert len(fc) == 20
            ....:     assert sorted(fd, reverse=True) == [9, 5, 3, 3]
            ....:     assert sorted(rel) == list(range(20))
        """
        n = self._n
        fp = self._fp

        fc = []         # faces seen along the walk
        fd = []         # face degrees seen along the walk
        rel = [-1] * n  # dart relabeling

        rel[i0] = 0   # first edge is relabelled (0,1)
        fc.append(0)
        c = 2         # current dart number (in the new labelling scheme)

        # walk along faces first, starting from i0
        # along the way, we collect unseen edges
        i = fp[i0]
        wait = deque([i0 ^ 1])
        d = 1
        while i != i0:
            if rel[i] == -1:
                j = i ^ 1
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
            assert rel[i0 ^ 1] != -1 and rel[i0 ^ 1] % 2 == 0
            rel[i0] = rel[i0 ^ 1] + 1
            fc.append(rel[i0])
            i = fp[i0]
            d = 1
            while i != i0:
                if rel[i] == -1:
                    j = i ^ 1
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
        fp = self._fp
        fl = self._fl
        sfd = self._fd

        is_fd_better = 0   # whether the current face degree works better
        is_fc_better = 0   # whether the current relabelling works better
        #  0 = equal
        #  1 = better
        # -1 = worse

        fd = []         # face degrees seen along the walk (want to maximize)
        fc = []         # edges seen along the walk (want to minimize)
        rel = [-1] * n  # dart relabeling

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
        wait = deque([i0 ^ 1])
        while i != i0:
            if rel[i] == -1:
                j = i ^ 1
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
            assert rel[i0 ^ 1] != -1 and rel[i0 ^ 1] % 2 == 0
            ii0 = rel[i0 ^ 1] + 1
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
                    j = i ^ 1
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

    # NOTE: this is broken for non-connected surfaces!
    def _is_canonical(self, i0):
        r"""
        Return a pair ``(answer, automorphisms)`` where answer is a boolean
        that says whether this map is in canonical form and ``automorphisms`` form
        a generating set of the group of automorphisms.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph

            sage: vp = '(0,15,3,6,8)(1,14,18,10,9,12,5,19)(2,7,4,17,11)(13,16)'
            sage: fp = '(0,19,14)(1,8,10,17,13,9,6,2,15)(3,11,18,5,7)(4,12,16)'
            sage: cm = FatGraph(vp, fp)
            sage: any(cm._is_canonical(i)[0] for i in range(20))
            True

        A genus 1 example with 4 symmetries::

            sage: vp = '(0,8,6,4,3,7)(1,9,11,5,2,10)'
            sage: fp = '(0,10,9)(2,4,11)(1,7,8)(3,5,6)'
            sage: cm = FatGraph(vp, fp)
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

    def automorphism_vertex_action(self, g):
        r"""
        Return the action of ``g`` on vertices.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: r = FatGraph('(0,5,4)(1,2,3)', '(0,3,1,4)(2)(5)')
            sage: A = r.automorphism_group()
            sage: r.automorphism_vertex_action(A.gens()[0])
            [1, 0]
        """
        if not perm_check(g, self._n):
            raise ValueError('g must be a permutation')
        vp = self._vp
        vl = self._vl
        p = [-1] * self._nv
        for i in range(self._n):
            if vp[i] == -1:
                continue
            if p[vl[i]] == -1:
                p[vl[i]] = vl[g[i]]
            elif p[vl[i]] != vl[g[i]]:
                raise ValueError('not an automorphism')
        return p

    def automorphism_face_action(self, g):
        r"""
        Return the action of ``g`` on faces.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: r = FatGraph('(0,5,4)(1,2,3)', '(0,3,1,4)(2)(5)')
            sage: A = r.automorphism_group()
            sage: r.automorphism_face_action(A.gens()[0])
            [0, 2, 1]
        """
        if not perm_check(g, self._n):
            raise ValueError('g must be a permutation')
        fp = self._fp
        fl = self._fl
        p = [-1] * self._nf
        for i in range(self._n):
            if fp[i] == -1:
                continue
            if p[fl[i]] == -1:
                p[fl[i]] = fl[g[i]]
            elif p[fl[i]] != fl[g[i]]:
                raise ValueError('not an automorphism')
        return p

    def automorphism_group(self, fix_vertices=False, fix_edges=False, fix_faces=False):
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

            sage: vp = '(0,9,6,5,7,4,8,2,1,3)'
            sage: fp = '(0,2,1,3,8)(4,6,5,7,9)'
            sage: cm = FatGraph(vp, fp)
            sage: cm.automorphism_group()
            PermutationGroupOrbit(10, [(0,4)(1,5)(2,6)(3,7)(8,9)])

        One can compute face stabilizer::

            sage: r = FatGraph('(0,5,4)(1,2,3)', '(0,3,1,4)(2)(5)')
            sage: r.automorphism_group().group_cardinality()
            2
            sage: r.automorphism_group(fix_faces=True).group_cardinality()
            1
            sage: r = FatGraph('(0,4,3)(5,2,1)', '(0,2,4,1,3,5)')
            sage: A = r.automorphism_group()
            sage: A.group_cardinality()
            6
            sage: r.automorphism_group(fix_faces=True) == A
            True

            sage: r = FatGraph('(0,2,1,3)', '(0,2,1,3)')
            sage: r.automorphism_group().group_cardinality()
            4
        """
        # NOTE: alternatively, one can compute the centralizer
        # self._symmetric_group().centralizer(self.monodromy_group())

        if fix_edges:
            raise NotImplementedError

        if not self.is_connected():
            raise NotImplementedError('not implemented for non-connected graphs')

        roots = self._good_starts(fix_vertices=fix_vertices, fix_faces=fix_faces)
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
                if ((not fix_vertices or perm_is_on_list_stabilizer(aut, self._vl, self._n)) and
                    (not fix_faces or perm_is_on_list_stabilizer(aut, self._fl, self._n))):
                    P.add_generator(aut)

        return P

    def relabel(self, r):
        r"""
        Relabel according to the permutation ``r``.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: from surface_dynamics.misc.permutation import perm_on_list_inplace, perm_invert

            sage: cm = FatGraph('(0,11,9,14,7,10,12,4,2,1,3)(5,15,8,6)(13)',
            ....:               '(0,2,1,3,4,6,14,5,12,13,10)(7,8,11)(9,15)')
            sage: r = [4,5,8,9,1,0,11,10,15,14,2,3,7,6,12,13]
            sage: cm.relabel(r)
            sage: cm._check()
        """
        n = self._n
        if any(r[i] ^ 1 != r [i ^ 1] for i in range(n)):
            raise ValueError('invalid relabelling permutation')
        perm_conjugate_inplace(self._vp, r, n)
        perm_conjugate_inplace(self._fp, r, n)
        perm_on_list_inplace(r, self._vl, n)
        perm_on_list_inplace(r, self._fl, n)
