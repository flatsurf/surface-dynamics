r"""
Deprecated module.

TESTS::

    sage: from surface_dynamics import *
    sage: R = RibbonGraph(vertices='(0,1,4,3)(5,2)',edges='(0,3)(1,2)(4,5)')
    doctest:warning
    ...
    DeprecationWarning: RibbonGraph is deprecated; use surface_dynamics.fat_graph.FatGraph instead
    sage: R
    Ribbon graph with 2 vertices, 3 edges and 3 faces
    sage: R.vertices()
    [[0, 1, 4, 3], [2, 5]]
    sage: R.edges()
    [[0, 3], [1, 2], [4, 5]]
    sage: R.faces()
    [[0, 4, 2], [3], [1, 5]]
    sage: R.dart_to_vertex(1)
    0
    sage: 1 in R.vertices()[0]
    True
    sage: R.dart_to_vertex(2)
    1
    sage: 2 in R.vertices()[1]
    True

    sage: R.dart_to_edge(3)
    0
    sage: R.dart_to_edge(4)
    2

    sage: R.dart_to_face(4)
    0
    sage: R.dart_to_face(3)
    1
    sage: e = '(0,1)(2,3)'
    sage: f = '(0,2,1,3)'
    sage: a = [1/2,1/2,1/2,1/2]
    sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
    doctest:warning
    ...
    DeprecationWarning: RibbonGraphWithAngles is deprecated; use surface_dynamics.fat_graph.FatGraph instead
    sage: r.spin_parity()
    1
"""
#*****************************************************************************
#       Copyright (C) 2019-2023 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function, absolute_import
from six.moves import range, map, filter, zip

from surface_dynamics.misc.permutation import (perm_init, constellation_init, perm_cycles, perm_invert,
        perm_check, perm_compose, equalize_perms, perm_orbit)

from sage.misc.cachefunc import cached_method

from sage.modules.free_module import FreeModule
from sage.matrix.constructor import matrix, identity_matrix

from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ


class RibbonGraph(SageObject):
    r"""
    Deprecated class.

    TESTS::

        sage: from surface_dynamics import *
        sage: RibbonGraph([],[],[])
        Ribbon graph with 1 vertex, 0 edge and 1 face
        sage: RibbonGraph('()','(0,1)','(0,1)')
        Ribbon graph with 2 vertices, 1 edge and 1 face
        sage: G = RibbonGraph('(0,3)(1,2)','(0,1)(2,3)','(0,2)(1,3)')
        sage: G
        Ribbon graph with 2 vertices, 2 edges and 2 faces
        sage: G.darts()
        [0, 1, 2, 3]
        sage: G.genus()
        0
        sage: G = RibbonGraph(edges='(0,2)(1,3)(4,6)(5,7)',faces='(0,1,2,3,4,5,6,7)')
        sage: G
        Ribbon graph with 1 vertex, 4 edges and 1 face
        sage: G.darts()
        [0, 1, 2, 3, 4, 5, 6, 7]
        sage: G.genus()
        2
        sage: G = RibbonGraph(vertices='(0,2,3,6)(1,4,5,7)')
        sage: G
        Ribbon graph with 2 vertices, 4 edges and 4 faces
        sage: G.edges()
        [[0, 1], [2, 3], [4, 5], [6, 7]]
        sage: G.genus()
        0
    """
    def __init__(self, vertices=None, edges=None, faces=None, connected=True, check=True):
        r"""
        TESTS::

            sage: from surface_dynamics import RibbonGraph
            sage: RibbonGraph([0r,1r],[1r,0r],[1r,0r])
            Ribbon graph with 2 vertices, 1 edge and 1 face

            sage: RibbonGraph([-1,1,-1,3],[-1,3,-1,1],[-1,3,-1,1])
            Ribbon graph with 2 vertices, 1 edge and 1 face

            sage: RibbonGraph(vertices='(1)(5)', edges='(1,5)', faces='(1,5)')
            Ribbon graph with 2 vertices, 1 edge and 1 face
        """
        from warnings import warn
        warn('RibbonGraph is deprecated; use surface_dynamics.fat_graph.FatGraph instead', DeprecationWarning)
        vertices, edges, faces = constellation_init(vertices, edges, faces, check=check)

        n = len(vertices)
        vp = [-1] * n
        fp = [-1] * n
        mapping = [-1] * n
        inv_mapping = [-1] * n
        num_darts = 0
        for i in range(n):
            ii = edges[i]
            if ii != -1 and mapping[i] == -1:
                mapping[i] = num_darts
                mapping[ii] = num_darts + 1
                inv_mapping[num_darts] = i
                inv_mapping[num_darts + 1] = ii
                num_darts += 2
        for i in range(n):
            if mapping[i] != -1:
                vp[mapping[i]] = mapping[vertices[i]]
                fp[mapping[i]] = mapping[faces[i]]
        vp = vp[:num_darts]
        fp = fp[:num_darts]
        self._fat_graph_map = mapping
        self._fat_graph_inv_map = inv_mapping
        from surface_dynamics.topology.fat_graph import FatGraph
        self._fat_graph = FatGraph(vp=vp, fp=fp, max_num_dart=n, mutable=True, check=True)

        self._total_darts = total_darts = len(vertices)
        num_darts = 0
        self._active_darts = [True] * total_darts
        for i,(jv, je, jf) in enumerate(zip(vertices, edges, faces)):
            if jv == -1:
                if je != -1 or jf != -1:
                    raise ValueError("inconsistent permutations")
                self._active_darts[i] = False
            else:
                if je == -1 or jf == -1:
                    raise ValueError("inconsistent permutations")
                num_darts += 1

        self._num_darts = num_darts

        self._vertices = vertices
        self._vertex_cycles = perm_cycles(vertices, singletons=True)
        self._dart_to_vertex_index = [None] * total_darts
        for i,c in enumerate(self._vertex_cycles):
            for j in c:
                self._dart_to_vertex_index[j] = i

        self._edges = edges
        self._edge_cycles = perm_cycles(edges, singletons=True)
        self._dart_to_edge_index = [None] * total_darts
        for i,edge in enumerate(self._edge_cycles):
            for e in edge:
                self._dart_to_edge_index[e] = i

        self._faces = faces
        self._face_cycles = perm_cycles(faces, singletons=True)
        self._dart_to_face_index = [None] * total_darts
        for i,c in enumerate(self._face_cycles):
            for j in c:
                self._dart_to_face_index[j] = i

        self._vertices_inv = [None] * total_darts
        self._edges_inv = [None] * total_darts
        self._faces_inv = [None] * total_darts
        for i in range(total_darts):
            if self._vertices[i] != -1:
                self._vertices_inv[self._vertices[i]] = i
                self._edges_inv[self._edges[i]] = i
                self._faces_inv[self._faces[i]] = i

        self._connected = bool(connected)
        if check:
            self._check()

    def _repr_(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *

            sage: RibbonGraph(edges='(0,1)(2,3)',faces='(0,1,2,3)')._repr_()
            'Ribbon graph with 3 vertices, 2 edges and 1 face'
            sage: RibbonGraph(edges='(0,1)',faces='(0)(1)')._repr_()
            'Ribbon graph with 1 vertex, 1 edge and 2 faces'
            sage: RibbonGraph(edges='(0,1)(2,3)',faces='(0,2)(1,3)')._repr_()
            'Ribbon graph with 2 vertices, 2 edges and 2 faces'
        """
        n = self._fat_graph.num_vertices()
        if n <= 1:
            vert_str = "%d vertex" %n
        else:
            vert_str = "%d vertices" %n
        n = self._fat_graph.num_edges()
        if n <= 1:
            edge_str = "%d edge" %n
        else:
            edge_str = "%d edges" %n
        n = self._fat_graph.num_faces()
        if n <= 1:
            face_str = "%d face" %n
        else:
            face_str = "%d faces" %n
        return "Ribbon graph with %s, %s and %s" %(vert_str,edge_str,face_str)

    @cached_method
    def _symmetric_group(self):
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        return SymmetricGroup([i for i, j in enumerate(self._fat_graph_map) if j != -1])

    def monodromy_group(self):
        r"""
        TESTS::

            sage: from surface_dynamics import RibbonGraph
            sage: r1 = RibbonGraph(vertices='(0)(2)(1,3,4,5)', edges='(0,1)(2,3)(4,5)')
            sage: G1 = r1.monodromy_group()
            sage: G1
            Subgroup ...
            sage: G1.is_isomorphic(SymmetricGroup(5))
            True
            sage: r2 = RibbonGraph(vertices='(0)(2)(1,3,4)(5,6,7)', edges='(0,1)(2,3)(4,5)(6,7)')
            sage: G2 = r2.monodromy_group()
            sage: G2
            Subgroup ...
            sage: G2.is_isomorphic(PSL(2,7))
            True
            sage: r3 = RibbonGraph(vertices='(1)(5)', edges='(1,5)', faces='(1,5)')
            sage: r3.monodromy_group()
            Subgroup ...
        """
        S = self._symmetric_group()
        v = S([self._fat_graph_inv_map[j] for j in self._fat_graph._vp[:self._fat_graph._n]])
        e = S([(self._fat_graph_inv_map[j], self._fat_graph_inv_map[j + 1]) for j in range(0, self._fat_graph._n, 2)])
        f = S([self._fat_graph_inv_map[j] for j in self._fat_graph._fp[:self._fat_graph._n]])
        assert (v * e * f).is_one()
        return S.subgroup([v, e, f])

    def automorphism_group(self, fix_vertices=False, fix_edges=False, fix_faces=False):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: r = RibbonGraph('(1,6,5)(2,3,4)', '(1,2)(3,4)(5,6)', '(1,4,2,5)(3)(6)')
            sage: r.automorphism_group().cardinality()
            2
            sage: r.automorphism_group(fix_faces=True).cardinality()
            1
            sage: r = RibbonGraph('(1,4,5)(2,6,3)', '(1,2)(3,4)(5,6)', '(1,3)(2,5)(4,6)')
            sage: r.automorphism_group().cardinality()
            6
            sage: r = RibbonGraph('(1,5,4)(6,3,2)', '(1,2)(3,4)(5,6)', '(1,3,5,2,4,6)')
            sage: A = r.automorphism_group()
            sage: A
            Subgroup ...
            sage: A.cardinality()
            6
            sage: r.automorphism_group(fix_faces=True) == A
            True
            sage: r = RibbonGraph('(1,3,2,4)', '(1,2)(3,4)', '(1,3,2,4)')
            sage: r.automorphism_group().cardinality()
            4
        """
        C = self._symmetric_group().centralizer(self.monodromy_group())
        if C.cardinality().is_one():
            return C

        if fix_vertices or fix_edges or fix_faces:
            S = self._symmetric_group()
            B = []
            if fix_vertices:
                B.append(self._vertex_cycles)
            if fix_edges:
                B.append(self._edge_cycles)
            if fix_faces:
                B.append(self._face_cycles)
            for BB in B:
                gens = []
                for s in BB:
                    if len(s) == 1:
                        continue
                    gens.append(S([(s[0], s[1])]))
                    if len(s) > 2:
                        gens.append(S([tuple(s)]))
                C = C.intersection(S.subgroup(gens))

        return C

    # TODO: perform the expansion as quasi-polynomial
    # TODO: see how the recursion formula translates on this generating series
    def length_rational_fraction(self, var='b'):
        return self._fat_graph.lengths_generating_series(variable_name=var)

    def add_extra_darts(self, n):
        m = self._total_darts
        if 2*n > m:
            self._total_darts = int(2*n)
            for p in [self._vertices, self._vertices_inv,
                      self._edges, self._edges_inv,
                      self._faces, self._faces_inv]:
                p.extend([None] * (2*n - m))
            self._active_darts.extend([False] * (2*n - m))
            self._check()

    def _check(self):
        if len(self._active_darts) != self._total_darts:
            raise ValueError("the length of active darts is not total_darts")
        if self._active_darts.count(True) != self._num_darts:
            raise ValueError("the number of darts do not coincide with active darts")

        for i in range(self._total_darts):
            if self._active_darts[i]:
                if self._vertices[i] == -1 or self._vertices_inv[i] is None:
                    raise ValueError("dart %d is active but has no vertex" %i)
                if self._edges[i] == -1 or self._edges_inv[i] is None:
                    raise ValueError("dart %d is active but has no edge" %i)
                if self._faces[i] == -1 or self._faces_inv[i] is None:
                    raise ValueError("dart %d is active but has no face" %i)

                if self._vertices[self._vertices_inv[i]] != i:
                    raise ValueError("vertices is not the inverse of vertices_inv")
                if self._edges[self._edges_inv[i]] != i:
                    raise ValueError("edges is not the inverse of edges_inv")
                if self._faces[self._faces_inv[i]] != i:
                    raise ValueError("faces is not the inverse of faces_inv")

                if self._faces[self._edges[self._vertices[i]]] != i:
                    raise ValueError("the Ribbon graph condition vef=() is not satisfied for %d" %i)

            else:
                if self._vertices[i] != -1 or self._vertices_inv[i] is not None:
                    raise ValueError("dart %d is not active but has a vertex" %i)
                if self._edges[i] != -1 or self._edges_inv[i] is not None:
                    raise ValueError("dart %d is not active but has an edge" %i)
                if self._faces[i] != -1 or self._faces_inv[i] is not None:
                    raise ValueError("dart %d is not active but has a face" %i)

        if self._connected and not self.is_connected():
            raise ValueError("the graph is not connected")

    def is_connected(self, force_computation=False):
        return self._fat_graph.is_connected()

    def relabel(self, perm=None):
        if perm is None:
            perm = [None] * self.num_darts()
            k = 0
            for i in range(self.num_darts()):
                if self._active_darts[i]:
                    perm[i] = k
                    k += 1

        vertices = [None] * self.num_darts()
        edges = [None] * self.num_darts()
        faces = [None] * self.num_darts()
        for i in range(self.num_darts()):
            if self._active_darts[i]:
                vertices[perm[i]] = perm[self._vertices[i]]
                edges[perm[i]] = perm[self._edges[i]]
                faces[perm[i]] = perm[self._faces[i]]

        return RibbonGraph(vertices, edges, faces)

    def num_darts(self):
        return self._num_darts

    def darts(self):
        return [i for i in range(self._total_darts) if self._fat_graph_map[i] != -1]

    def num_vertices(self):
        return self._fat_graph.num_vertices()

    def vertex_perm(self, copy=True):
        return self._vertices

    def vertex_orbit(self, i):
        j = self._fat_graph_map[i]
        if j != -1:
            return [self._fat_graph_inv_map[x] for x in perm_orbit(self._fat_graph._vp, j)]
        return None

    def vertices(self):
        return [[self._fat_graph_inv_map[x] for x in cycle] for cycle in self._fat_graph.vertices()]

    def dart_to_vertex(self,i):
        j = self._fat_graph_map[i]
        if j != -1:
            return self._fat_graph._vl[j]
        raise ValueError("dart %d is not active" %i)

    def num_edges(self):
        return self._fat_graph.num_edges()

    def edge_perm(self):
        return self._edges

    def edge_orbit(self, i):
        j = self._fat_gramp_map[i]
        if j != -1:
            return [[self._fat_graph_inv_map[x], self._fat_graph_inv_map[x + 1]] for x in range(self._fat_graph._n)]
        return None

    def edges(self):
        return [[self._fat_graph_inv_map[x] for x in cycle] for cycle in self._fat_graph.edges()]

    def dart_to_edge(self, i, orientation=False):
        j = self._fat_graph_map[i]
        if j == -1:
            raise ValueError("dart %d is not active" %i)
        return (j // 2, j % 2) if orientation else j // 2

    def num_faces(self):
        return self._fat_graph.num_faces()

    def face_perm(self):
        return self._faces

    def face_orbit(self, i):
        j = self._fat_graph_map[i]
        if j != -1:
            return [self._fat_graph_inv_map[x] for x in perm_orbit(self._fat_graph._fp, j)]
        return None

    def faces(self):
        return [[self._fat_graph_inv_map[x] for x in cycle] for cycle in self._fat_graph.faces()]

    def dart_to_face(self, i):
        j = self._fat_graph_map[i]
        if j == -1:
            raise ValueError("dart {} is not active".format(i))
        return self._fat_graph._fl[j]

    def dual(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: r = RibbonGraph(edges='(0,1)',faces='(0)(1)'); r
            Ribbon graph with 1 vertex, 1 edge and 2 faces
            sage: r.dual()
            Ribbon graph with 2 vertices, 1 edge and 1 face
        """
        return RibbonGraph(
            vertices=perm_invert(self._faces),
            edges=self._edges,
            faces=perm_invert(self._vertices))

    def euler_characteristic(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: r = RibbonGraph(edges='(0,1)(2,3)(4,5)',faces='(0,2,4)(1)(3,5)')
            sage: r.euler_characteristic()
            2
            sage: r = RibbonGraph(edges='(0,1)(2,3)',faces='(0,2,1,3)')
            sage: r.euler_characteristic()
            0
        """
        return self._fat_graph.euler_characteristic()

    def is_plane(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: r = RibbonGraph(vertices='(0)(1)',edges='(0,1)')
            sage: r.is_plane()
            True
            sage: r = RibbonGraph(vertices='(0,1)',edges='(0,1)')
            sage: r.is_plane()
            True
            sage: r = RibbonGraph(edges='(0,1)(2,3)',faces='(0,2)(1,3)')
            sage: r.is_plane()
            True
            sage: r = RibbonGraph(edges='(0,1)(2,3)',faces='(0,2,1,3)')
            sage: r.is_plane()
            False
        """
        return self._fat_graph.euler_characteristic() == 2

    def is_plane_tree(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: r = RibbonGraph(vertices='(0)(1)',edges='(0,1)')
            sage: r.is_plane_tree()
            True
            sage: r = RibbonGraph(vertices='(0)(1,2,4)(3)(5)',edges='(0,1)(2,3)(4,5)')
            sage: r.is_plane_tree()
            True
            sage: r = RibbonGraph(vertices='(0,1)',edges='(0,1)')
            sage: r.is_plane_tree()
            False
            sage: r.is_plane()
            True
        """
        return self.num_faces() == 1 and self.genus() == 0

    def is_triangulated(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: r = RibbonGraph(edges='(0,1)(2,3)(4,5)',faces='(0,2,4)(1,5,3)')
            sage: r.is_triangulated()
            True
            sage: r = RibbonGraph(edges='(0,1)(2,3)',faces='(0,2,1,3)')
            sage: r.is_triangulated()
            False
        """
        return self._fat_graph.is_triangulation()

    def genus(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: R = RibbonGraph(vertices='(1)(2)',edges='(1,2)')
            sage: R.genus()
            0
            sage: e='(1,3)(2,4)'
            sage: f='(1,2,3,4)'
            sage: RibbonGraph(edges=e,faces=f).genus()
            1
            sage: e='(1,3)(2,4)(5,7)(6,8)'
            sage: f='(1,2,3,4,5,6,7,8)'
            sage: RibbonGraph(edges=e,faces=f).genus()
            2
            sage: e='(1,3)(2,4)(5,7)(6,8)(9,11)(10,12)'
            sage: f='(1,2,3,4,5,6,7,8,9,10,11,12)'
            sage: RibbonGraph(edges=e,faces=f).genus()
            3
        """
        return self._fat_graph.genus()

    def spanning_tree(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: R = RibbonGraph('(1,2,3)','(1,2)(3,4)')
            sage: R
            Ribbon graph with 2 vertices, 2 edges and 2 faces
            sage: T,o = R.spanning_tree()
            sage: T
            Digraph on 2 vertices
            sage: T.edges(sort=True)
            [(0, 1, (3, 4))]
            sage: o
            [(1, 2)]
            sage: R = RibbonGraph('(1,2,3)(4,5,6)','(1,2)(3,4)(5,6)')
            sage: R
            Ribbon graph with 2 vertices, 3 edges and 3 faces
            sage: T,o = R.spanning_tree()
            sage: T
            Digraph on 2 vertices
            sage: T.edges(sort=True)
            [(0, 1, (3, 4))]
            sage: o
            [(1, 2), (5, 6)]
            sage: e = '(1,3)(5,7)(2,4)(6,8)'
            sage: f = '(1,2,3,4,5,6,7,8)'
            sage: R = RibbonGraph(edges=e, faces=f)
            sage: T,o = R.spanning_tree()
            sage: T
            Digraph on 1 vertex
            sage: o
            [[1, 3], [2, 4], [5, 7], [6, 8]]
        """
        from sage.graphs.digraph import DiGraph

        d = self.darts()
        v = self.vertices()
        e = self.edge_perm()

        if self.num_darts() == 0:
            return DiGraph(),[]
        if self.num_vertices() == 1:
            return DiGraph({0:[]}),self.edges()

        T = DiGraph()

        v0 = 0
        for root in v[0]:
            v1 = self.dart_to_vertex(e[root])
            if v1 != 0:
                break

        T.add_edge(v0,v1,(root,e[root]))
        o = []
        seen = set([v0,v1])       # seen vertices
        waiting = [e[root],root]  # waiting darts
        cc = []

        while waiting:
            ii = waiting.pop() # this is a dart
            v0 = self.dart_to_vertex(ii)
            seen.add(v0)
            for j in self.vertex_orbit(ii)[1:]:
                v1 = self.dart_to_vertex(e[j])
                if v1 in seen:
                    if j < e[j]:
                        o.append((j,e[j]))
                else:
                    T.add_edge(v0,v1,(j,e[j]))
                    waiting.append(e[j])
                    seen.add(v1)

        return T, sorted(o)

    def collapse(self, spanning_tree=None):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: R = RibbonGraph(vertices='(0,1,2,5)(3,7)(4,10,9)(6,11,12)(8,13)')
            sage: R.genus()
            1
            sage: R.num_vertices()
            5
            sage: R.num_edges()
            7
            sage: R.num_faces()
            2
            sage: R2 = R.collapse()
            sage: R2
            Ribbon graph with 1 vertex, 3 edges and 2 faces
            sage: R
            Ribbon graph with 5 vertices, 7 edges and 2 faces
            sage: R3 = R2.dual().collapse().dual()
            sage: R3
            Ribbon graph with 1 vertex, 2 edges and 1 face
        """
        from copy import deepcopy

        if spanning_tree is None:
            spanning_tree,_ = self.spanning_tree()

        darts_to_kill = set([])
        for v0,v1,e in spanning_tree.edges(sort=True):
            darts_to_kill.add(e[0])
            darts_to_kill.add(e[1])

        new_edges = []
        for e in self.edges():
            if e[0] not in darts_to_kill:
                new_edges.append(e)

        new_faces = []
        for f in self.faces():
            ff = tuple(i for i in f if i not in darts_to_kill)
            if ff:
                new_faces.append(ff)

        return RibbonGraph(edges=tuple(new_edges), faces=tuple(new_faces))

    def boundaries(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: r = RibbonGraph('(1,2,3)(4,5,6)','(1,2)(3,4)(5,6)')
            sage: r.boundaries()
            [[(1, 2)],  [(2, 1), (3, 4), (6, 5), (4, 3)], [(5, 6)]]
            sage: r = RibbonGraph('(1,2,3)(4,5)(6,7,8)',edges='(1,2)(3,4)(5,6)(7,8)')
            sage: r.boundaries()
            [[(1, 2)],  [(2, 1), (3, 4), (5, 6), (8, 7), (6, 5), (4, 3)], [(7, 8)]]
        """
        e = self.edge_perm()
        return sorted([[(i,e[i]) for i in f] for f in self.faces()])

    def cycle_basis(self, intersection=False, verbose=False):
        r"""
        TESTS::

            sage: from surface_dynamics import RibbonGraph
            sage: r = RibbonGraph('(1,2,3)(4,5,6)','(1,2)(3,4)(5,6)')
            sage: r.cycle_basis()
            []
            sage: r = RibbonGraph('(1,2,3)(4,5)(6,7,8)',edges='(1,2)(3,4)(5,6)(7,8)')
            sage: r.cycle_basis()
            []
            sage: r = RibbonGraph('(1,4,5)(2,3)(6,7,8)',edges='(1,2)(3,4)(5,6)(7,8)')
            sage: r.cycle_basis()
            []
            sage: e = '(1,3)(2,4)(5,7)(6,8)'
            sage: f = '(1,2,3,4,5,6,7,8)'
            sage: r = RibbonGraph(edges=e,faces=f)
            sage: r.cycle_basis()
            [[(1, 3)], [(2, 4)], [(5, 7)], [(6, 8)]]
            sage: r.dual().cycle_basis()
            [[(1, 3)], [(2, 4)], [(5, 7)], [(6, 8)]]
            sage: r.cycle_basis(intersection=True)
            (
                                                      [ 0  1  0  0]
                                                      [-1  0  0  0]
                                                      [ 0  0  0  1]
            [[(1, 3)], [(2, 4)], [(5, 7)], [(6, 8)]], [ 0  0 -1  0]
            )
            sage: r.dual().cycle_basis(intersection=True)
            (
                                                      [ 0 -1  0  0]
                                                      [ 1  0  0  0]
                                                      [ 0  0  0 -1]
            [[(1, 3)], [(2, 4)], [(5, 7)], [(6, 8)]], [ 0  0  1  0]
            )
            sage: f = '(0,10,13)(6,17,11)(2,14,7)(15,12,3)(16,20,19)(18,1,9)(4,22,21)(23,8,5)'
            sage: e = tuple((i,i+1) for i in range(0,24,2))
            sage: r = RibbonGraph(edges=e,faces=f); r
            Ribbon graph with 2 vertices, 12 edges and 8 faces
            sage: c,m = r.cycle_basis(intersection=True)
            sage: c
            [[(0, 1), (4, 5)], [(6, 7)], [(14, 15), (1, 0)], [(22, 23), (1, 0)]]
            sage: m
            [ 0  0  0  1]
            [ 0  0  1  0]
            [ 0 -1  0  0]
            [-1  0  0  0]
        """
        if intersection:
            cycles, intersection_matrix = self._fat_graph.cycle_basis(True)
        else:
            cycles = self._fat_graph.cycle_basis(False)
        cycles = [[(self._fat_graph_inv_map[x], self._fat_graph_inv_map[x ^ 1]) for x in cycle] for cycle in cycles]
        return (cycles, intersection_matrix) if intersection else cycles

    def is_cycle(self, c):
        for i in range(len(c)-1):
            if self.dart_to_vertex(c[i][1]) != self.dart_to_vertex(c[i+1][0]):
                return False
        if self.dart_to_vertex(c[-1][1]) != self.dart_to_vertex(c[0][0]):
            return False
        return True


class RibbonGraphWithAngles(RibbonGraph):
    r"""
    Deprecated class.
    """
    def __init__(self, vertices=None, edges=None, faces=None, angles=None):
        from warnings import warn
        warn('RibbonGraphWithAngles is deprecated; use surface_dynamics.fat_graph.FatGraph instead', DeprecationWarning)

        r = RibbonGraph(vertices,edges,faces)
        RibbonGraph.__init__(self, r.vertex_perm(), r.edge_perm(), r.face_perm())

        if len(angles) != self.num_darts():
            raise ValueError("there are %d angles and %d darts" %(len(angles),self.num_darts()))
        self._angles = list(map(QQ,angles))
          # angle between a dart and its vertex-neighbour
          # (rational number as multiple of pi)

        self._total_angle = []
          # total angle around vertices
          # (integer which corresponds to a multiple of pi)

        for v in self.vertices():
            self._total_angle.append(sum(angles[i] for i in v))

        for f in self.faces():
            a = sum(self._angles[i] for i in f)
            if a != len(f)-2:
                raise ValueError("the angle of a face should be (nb_edges - 2) x pi")

    def angle_between_darts(self, d1, d2):
        v = self.vertex_orbit(d1)
        if d2 not in v:
            raise ValueError("d1=%s and d2=%s are not at the same vertex" %(str(d1),str(d2)))

        a = 0
        i = 0
        while v[i] != d2:
            a += self._angles[v[i]]
            i += 1
        return a

    def angle_at_vertex(self,v):
        return self._total_angle[v]

    def angle_at_vertices(self):
        return self._total_angle

    def winding(self, c):
        a = 0
        for i in range(len(c)-1):
            d1 = c[i][1]
            d2 = c[i+1][0]
            if self.dart_to_vertex(d1) != self.dart_to_vertex(d2):
                raise ValueError("c is not a cycle")
            a += self.angle_between_darts(d1,d2)-1
        d1 = c[-1][1]
        d2 = c[0][0]
        if self.dart_to_vertex(d1) != self.dart_to_vertex(d2):
            raise ValueError("c is not a cycle")
        a += self.angle_between_darts(d1,d2)-1

        return a

    def holonomy_representation(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: e = '(0,1)(2,3)'
            sage: f = '(0,2,1,3)'
            sage: a = [1/2,1/2,1/2,1/2]
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.holonomy_representation()
            ([0], [0, 0])
            sage: e = tuple((i,i+1) for i in range(0,24,2))
            sage: f = '(0,20,7,10)(16,22,19,21)(2,9,5,23)(14,3,17,1)(12,8,15,11)(18,4,13,6)'
            sage: a = [1/2]*24
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.holonomy_representation()
            ([3/2, 3/2, 3/2, 3/2, 3/2, 3/2, 3/2, 3/2], [])
            sage: e = '(0,1)(2,3)(4,5)'
            sage: f = '(0,2,4)(1,5,3)'
            sage: a = [1/2,1/6,1/3,1/3,1/6,1/2]
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.holonomy_representation()
            ([1, 1/2, 1/2], [])
            sage: a = [1/3,7/15,1/5,1/5,7/15,1/3]
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.holonomy_representation()
            ([2/3, 2/3, 2/3], [])
        """
        from sage.functions.other import floor

        l1 = []
        for c in range(self.num_vertices()):
            w = self.angle_at_vertex(c)
            l1.append(w - 2*floor(w/2))

        l2 = []
        for c in self.cycle_basis():
            w = self.winding(c)
            l2.append(w - 2*floor(w/2))

        return l1,l2

    def has_trivial_holonomy(self):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: e = '(0,1)(2,3)'
            sage: f = '(0,2,1,3)'
            sage: a = [1/2,1/2,1/2,1/2]
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.has_trivial_holonomy()
            True
            sage: e = '(0,1)(2,3)(4,5)'
            sage: f = '(0,2,4)(1,5,3)'
            sage: a = [1/3,7/15,1/5,1/5,7/15,1/3]
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.has_trivial_holonomy()
            False
        """
        l1,l2 = self.holonomy_representation()
        return all(i==0 for i in l1) and all(i==0 for i in l2)

    def spin_parity(self, check=True, verbose=False):
        r"""
        TESTS::

            sage: from surface_dynamics import *
            sage: e = '(0,1)(2,3)'
            sage: f = '(0,2,1,3)'
            sage: a = [1/2,1/2,1/2,1/2]
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.spin_parity()
            1
            sage: e = '(0,1)(2,3)(4,5)(6,7)'
            sage: f = '(0,2,4,3,6,1,7,5)'
            sage: a = [1/2,1/2,1,1/2,1/2,1,3/2,1/2]
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.spin_parity()
            1
            sage: e = '(0,1)(2,3)(4,5)(6,7)'
            sage: f = '(0,2,4,6,1,3,5,7)'
            sage: a = [1/2,1/2,1,1,1,1,1/2,1/2]
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.spin_parity()
            1
            sage: e = '(0,1)(2,3)(4,5)(6,7)'
            sage: f = '(0,2,4,6,1,3,5,7)'
            sage: a = [3/4]*8
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.spin_parity()
            1
            sage: e = '(0,1)(2,3)(4,5)(6,7)(8,9)(10,11)'
            sage: f1 = '(0,4,6,8,10,2,1,9,11,5,7,3)'
            sage: f2 = '(0,4,6,8,10,2,1,5,7,9,11,3)'
            sage: a = [1/2,1/2,1/2,1/2] + [1]*8
            sage: r1 = RibbonGraphWithAngles(edges=e,faces=f1,angles=a)
            sage: r1.spin_parity()
            1
            sage: r2 = RibbonGraphWithAngles(edges=e,faces=f2,angles=a)
            sage: r2.spin_parity()
            0
        """
        from sage.rings.finite_rings.finite_field_constructor import GF
        # mod F2 we have: q(x+y) = B(x,y) + q(x) + q(y)

        if not self.has_trivial_holonomy():
            raise ValueError("the surface does not have trivial holonomy")
        if any((i+2)%4 for i in self.angle_at_vertices()):
            raise ValueError("each angle should be odd multiple of 2pi")

        GF2 = GF(2)

        c, M = self.cycle_basis(intersection=True)

        winding = []
        for cc in c:
            w = self.winding(cc)
            if w % 2 != 0:
                raise ValueError("fatal error ! each winding should be a multiple of 2")
            winding.append(GF2(w//2))

        if verbose:
            print("cycles with winding")
            for i in range(len(c)):
                print(c[i], winding[i])
            print("intersection matrix on Z")
            print(M)

        # compute a base change to get a symplectic basis
        _,P = M.symplectic_form()
        M = M.change_ring(GF2)
        P = P.change_ring(GF2)
        if verbose:
            print("base change for symplectic basis on GF(2)")
            print(P)

        g = self.genus()

        s = GF2(0)
        for i in range(g):
            # 1. computation of q(P.row(i))
            a = P.row(i)
            a_indices = [j for j in range(2*g) if a[j] != 0]
            ## winding + nb_components
            t_a = sum(winding[i]+1 for i in a_indices)
            ## self intersection
            for j1 in range(len(a_indices)):
                for j2 in range(j1+1,len(a_indices)):
                    t_a += M[a_indices[j1],a_indices[j2]]

            # 2. computation of q(P.row(g+i))
            b = P.row(g+i)
            b_indices = [j for j in range(2*g) if b[j] != 0]
            ## winding + nb_components
            t_b = sum(winding[i]+1 for i in b_indices)
            ## self intersection
            for j1 in range(len(b_indices)):
                for j2 in range(j1+1,len(b_indices)):
                    t_b += M[b_indices[j1],b_indices[j2]]

            # 3. add to s the contribution of the couple
            if verbose:
                print("contribution from %d is %d * %d = %d" % (i, t_a, t_b, t_a * t_b))
            s += t_a*t_b

        return s


def angle(v):
    r"""
    Return the argument of the vector ``v``.
    """
    from math import acos,asin,sqrt,pi
    x = float(v[0])
    y = float(v[1])
    r = sqrt(x**2 + y**2)
    if abs(x) >= abs(y):
        if x >= 0:
            return asin(y / r) / pi
        else:
            return -asin(y / r) / pi
    else:
        if y >= 0:
            return acos(x / r) / pi
        else:
            return -acos(x / r) / pi
