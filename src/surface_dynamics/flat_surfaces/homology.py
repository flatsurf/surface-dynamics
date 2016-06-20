r"""
Simplicial complex, homology of surfaces and translation surfaces

In this module are implemented simple homology computation for translation
surfaces. There are three main classes:

- :class:`RibbonGraph`: decomposition of a surface into polygons. The
  combinatorics is stored as a triple of permutations `v` (vertices), `e`
  (edges), `f` (faces) so that the product `vef` is the identity. The domain of
  the permutations correspond to the half edges or *darts*. The permutation `e` is
  an involution without fixed point so that `e(i)` is the other half of the edge
  starting at `i`. The permutation `v` is obtained by turning around a vertex,
  while `f` turning around a face.

- :class:`RibbonGraphWithAngles`: a ribbon graph with an additional angle
  structure.

- :class:`RibbonGraphWithHolonomies`: a ribbon graph with an additional holonomy
  structure on its edges.

EXAMPLES::

    sage: from surface_dynamics.all import *

To create a ribbon graph you just need to fix two of the permutations `v`, `e`,
`f`::

    sage: R = RibbonGraph(vertices='(0,1,4,3)(5,2)',edges='(0,3)(1,2)(4,5)')
    sage: R
    Ribbon graph with 2 vertices, 3 edges and 3 faces

The vertices, edges and faces are by definition the cycles of the permutation.
Calling the method :meth:`~RibbonGraph.vertices`, :meth:`~RibbonGraph.edges` or
:meth:`~RibbonGraph.faces` gives you access to these cycles::

    sage: R.vertices()
    ([0, 1, 4, 3], [2, 5])
    sage: R.edges()
    ([0, 3], [1, 2], [4, 5])
    sage: R.faces()
    ([0, 4, 2], [1, 5], [3])

Given a half edge (i.e. a dart), you can get the index of the vertex, edge or
face it belongs with the methods :meth:`~RibbonGraph.dart_to_vertex`,
:meth:`~RibbonGraph.dart_to_edge` and :meth:`~RibbonGraph.dart_to_edge`::

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
    2

To initialize a ribbon graph with angles, you have to input the standard data to
initialize a ribbon graph plus a list of positive rational numbers which
corresponds to the angles between darts (more precisely, the number at position
i is the angle between i and v(i))::

    sage: e = '(0,1)(2,3)'
    sage: f = '(0,2,1,3)'
    sage: a = [1/2,1/2,1/2,1/2]
    sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
    sage: r.spin_parity()
    1
"""
from surface_dynamics.misc.permutation import (init_perm, perm_cycle_tuples, perm_invert,
        perm_check, perm_compose, equalize_perms, perm_orbit)

from sage.misc.cachefunc import cached_method

from sage.modules.free_module import FreeModule
from sage.matrix.constructor import matrix, identity_matrix

from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ


def clean_perm_data(vertices, edges, faces, check):
    r"""
    TESTS::

        sage: from surface_dynamics.flat_surfaces.homology import clean_perm_data
        sage: clean_perm_data([2,1,0],[1,2,0],None, check=True)
        ([2, 1, 0], [1, 2, 0], [0, 2, 1])
        sage: clean_perm_data([1,0,2],[2,1,0],None, check=True)
        ([1, None, 2], [2, None, 0], [2, None, 1])
    """
    if vertices is not None:
        vertices = init_perm(vertices)
        if check: perm_check(vertices)

    if edges is not None:
        edges = init_perm(edges)
        if check: perm_check(edges)

    if faces is not None:
        faces = init_perm(faces)
        if check: perm_check(faces)

    if edges is None:
        if vertices is None and faces is None:
            raise ValueError("at least vertices or faces should be not None")
        if not (vertices is None or faces is None):
            equalize_perms([vertices,faces])
            edges = perm_compose(perm_invert(vertices),perm_invert(faces))
        else:
            if vertices is None:
                n = len(faces)
            elif faces is None:
                n = len(vertices)
            if n%2:
                raise ValueError("there should be an even number of darts")
            edges = []
            for i in xrange(n//2):
                edges.append(2*i+1)
                edges.append(2*i)
            if vertices is None:
                vertices = perm_compose(perm_invert(faces),perm_invert(edges))
            elif faces is None:
                faces = perm_compose(perm_invert(edges),perm_invert(vertices))
    elif vertices is None:
        if edges is None or faces is None:
            raise ValueError("at least two of the entries should be not None")
        equalize_perms([edges,faces])
        vertices = perm_compose(perm_invert(faces),perm_invert(edges))
    elif faces is None:
        if vertices is None or edges is None:
            raise ValueError("at least two of the entries should be not None")
        equalize_perms([vertices,edges])
        faces = perm_compose(perm_invert(edges),perm_invert(vertices))
    else:
        equalize_perms([vertices,edges,faces])

    for i in range(len(edges)):
        if edges[i] == i:
            vertices[i] = edges[i] = faces[i] = None

    return (vertices, edges, faces)

class RibbonGraph(SageObject):
    r"""
    Generic class for Ribbon graph.

    A Ribbon graph is a graph embedded in a surface. This class uses
    representation as a triple ``(v,e,f)`` of permutations such that `vef = 1`
    and the action of the group generated by `v,e,f` acts transitvely in the
    the domain. The cycles of ``v`` are considered as vertices, the ones of
    ``e`` are considered as edges and the ones of ``f`` as the faces. Each
    element of the domain is a half-edge which is called a *dart*. A dart is
    also associated to an oriented edge.

    The domain of the permutations must be a subset of [0, ..., N-1] for some N.

    A dense ribbon graph has the following attributes

      - total_darts - non negative integer - the total number darts
      - num_darts - non negative integer - the number of active darts
      - active_darts - bitset - list of lengths _total_darts with True or
        False. The position i is True if i is an active dart.

      - vertices, vertices_inv - list - partial permutations of [0,N] which are
        inverse of each other
      - vertex_cycles - the cycles of the partial permutation vertices
      - dart_to_vertex_index

      - edges, edges_inv - list - partial permutations of [0,N] which are
        inverse of each other
      - edge_cycles - the cycles of the partial permutation edge
      - dart_to_edge_index

      - faces, faces_inv - list - partial permutations of [0,N] which are
        inverse of each other
      - face_cycles - the cycles of the partial permutation faces
      - dart_to_face_index

    EXAMPLES::

        sage: from surface_dynamics.all import *


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
        ([0, 1], [2, 3], [4, 5], [6, 7])
        sage: G.genus()
        0
    """
    def __init__(self, vertices=None, edges=None, faces=None, check=True):

        vertices, edges, faces = clean_perm_data(vertices, edges, faces, check)

        self._total_darts = total_darts = len(vertices)
        num_darts = 0
        self._active_darts = [True] * total_darts
        for i,j in enumerate(edges):
            if j is None or i == j:
                vertices[i] = edges[i] = faces[i] = None
                self._active_darts[i] = False
            else:
                num_darts += 1

        self._num_darts = num_darts

        self._vertices = vertices
        self._vertex_cycles = perm_cycle_tuples(vertices, singletons=True)
        self._dart_to_vertex_index = [None] * total_darts
        for i,c in enumerate(self._vertex_cycles):
            for j in c:
                self._dart_to_vertex_index[j] = i

        self._edges = edges
        self._edge_cycles = perm_cycle_tuples(edges, singletons=True)
        self._dart_to_edge_index = [None] * total_darts
        for i,e in enumerate(self._edge_cycles):
            self._dart_to_edge_index[e[0]] = i
            self._dart_to_edge_index[e[1]] = i

        self._faces = faces
        self._face_cycles = perm_cycle_tuples(faces, singletons=True)
        self._dart_to_face_index = [None] * total_darts
        for i,c in enumerate(self._face_cycles):
            for j in c:
                self._dart_to_face_index[j] = i

        self._vertices_inv = [None] * total_darts
        self._edges_inv = [None] * total_darts
        self._faces_inv = [None] * total_darts
        for i in xrange(total_darts):
            if self._vertices[i] is not None:
                self._vertices_inv[self._vertices[i]] = i
                self._edges_inv[self._edges[i]] = i
                self._faces_inv[self._faces[i]] = i

        if check:
            self._check()

    def _repr_(self):
        r"""
        String representation.

        TESTS::

            sage: from surface_dynamics.all import *

            sage: RibbonGraph(edges='(0,1)(2,3)',faces='(0,1,2,3)')._repr_()
            'Ribbon graph with 3 vertices, 2 edges and 1 face'
            sage: RibbonGraph(edges='(0,1)',faces='(0)(1)')._repr_()
            'Ribbon graph with 1 vertex, 1 edge and 2 faces'
            sage: RibbonGraph(edges='(0,1)(2,3)',faces='(0,2)(1,3)')._repr_()
            'Ribbon graph with 2 vertices, 2 edges and 2 faces'
        """
        n = self.num_vertices()
        if n <= 1:
            vert_str = "%d vertex" %n
        else:
            vert_str = "%d vertices" %n
        n = self.num_edges()
        if n <= 1:
            edge_str = "%d edge" %n
        else:
            edge_str = "%d edges" %n
        n = self.num_faces()
        if n <= 1:
            face_str = "%d face" %n
        else:
            face_str = "%d faces" %n
        return "Ribbon graph with %s, %s and %s" %(vert_str,edge_str,face_str)

    def add_extra_darts(self,n):
        r"""
        Add extra darts to the current vertex in order to support a total of
        ``n`` darts.
        """
        m = self._total_darts
        if n > m:
            self._total_darts = int(n)
            for p in [self._vertices, self._vertices_inv,
                      self._edges,self._edges_inv,
                      self._faces, self._faces_inv]:
                p.extend([None] * (n - m))
            self._active_darts.extend([False] * (n-m))
            self._check()

    def _check(self):
        r"""
        Check that the data of the Ribbon graph is coherent
        """
        from sage.graphs.graph import Graph
        G = Graph()

        if len(self._active_darts) != self._total_darts:
            raise ValueError, "the length of active darts is not total_darts"
        if self._active_darts.count(True) != self._num_darts:
            raise ValueError, "the number of darts do not coincide with active darts"

        for i in xrange(self._total_darts):
            if self._active_darts[i]:
                G.add_edge(i,self._vertices[i])
                G.add_edge(i,self._edges[i])
                G.add_edge(i,self._faces[i])
                if self._vertices[i] is None or self._vertices_inv[i] is None:
                    raise ValueError, "dart %d is active but has no vertex" %i
                if self._edges[i] is None or self._edges_inv[i] is None:
                    raise ValueError, "dart %d is active but has no edge" %i
                if self._faces[i] is None or self._faces_inv[i] is None:
                    raise ValueError, "dart %d is active but has no face" %i

                if self._vertices[self._vertices_inv[i]] != i:
                    raise ValueError, "vertices is not the inverse of vertices_inv"
                if self._edges[self._edges_inv[i]] != i:
                    raise ValueError, "edges is not the inverse of edges_inv"
                if self._faces[self._faces_inv[i]] != i:
                    raise ValueError, "faces is not the inverse of faces_inv"

                if self._faces[self._edges[self._vertices[i]]] != i:
                    raise ValueError, "the Ribbon graph condition vef=() is not satisfied for %d" %i
                if self._edges[i] == i or self._edges[self._edges[i]] != i:
                    raise ValueError, "edges is not an involution without fixed point for %d" %i

            else:
                if self._vertices[i] is not None or self._vertices_inv[i] is not None:
                    raise ValueError, "dart %d is not active but has a vertex" %i
                if self._edges[i] is not None or self._edges_inv[i] is not None:
                    raise ValueError, "dart %d is not active but has an edge" %i
                if self._faces[i] is not None or self._faces_inv[i] is not None:
                    raise ValueError, "dart %d is not active but has a face" %i

        if not G.is_connected():
            raise ValueError, "the graph is not connected"

    def relabel(self, perm=None):
        r"""
        perm is a of range(0,N)

        If ``perm`` is None, relabel the darts on 0,2M keeping the relative
        order of the darts.
        """
        if perm is None:
            perm=[None]*self.num_darts()
            k = 0
            for i in xrange(self.num_darts()):
                if self._active_darts[i]:
                    perm[i] = k
                    k += 1

        vertices = [None] * self.num_darts()
        edges = [None] * self.num_darts()
        faces = [None] * self.num_darts()
        for i in xrange(self.num_darts()):
            if self._active_darts[i]:
                vertices[perm[i]] = perm[self._vertices[i]]
                edges[perm[i]] = perm[self._edges[i]]
                faces[perm[i]] = perm[self._faces[i]]

        return RibbonGraph(vertices,edges,faces)

    #
    # Darts
    #

    def num_darts(self):
        r"""
        Returns the number of darts.
        """
        return self._num_darts

    def darts(self):
        r"""
        Return the list of darts
        """
        return [i for i in xrange(self._total_darts) if self._active_darts[i]]

    def num_vertices(self):
        r"""
        Returns the number of vertices.
        """
        return max(1,len(self._vertex_cycles))

    def vertex_perm(self):
        r"""
        Returns the permutation that define the vertices.
        """
        return self._vertices

    def vertex_orbit(self, i):
        r"""
        Return the orbit of ``i`` under the permutation that define the
        vertices.
        """
        if self._active_darts[i]:
            return perm_orbit(self._vertices,i)
        return None

    def vertices(self):
        r"""
        Return the list of vertices as cycles decomposition of the vertex
        permutation.
        """
        return self._vertex_cycles

    def dart_to_vertex(self,i):
        r"""
        Return the vertex on which the dart ``i`` is attached.
        """
        if self._active_darts[i]:
            return self._dart_to_vertex_index[i]
        raise ValueError, "dart %d is not active" %i

    #
    # Edges
    #

    def num_edges(self):
        r"""
        Returns the number of edges.
        """
        return len(self._edge_cycles)

    def edge_perm(self):
        r"""
        Return the permutation that define the edges.
        """
        return self._edges

    def edge_orbit(self, i):
        r"""
        Return the orbit of the dart ``i`` under the permutation that defines
        the edges.
        """
        if self._active_darts[i]:
            return perm_orbit(self._edges,i)
        return None

    def edges(self):
        r"""
        Return the set of edges.
        """
        return self._edge_cycles

    def dart_to_edge(self, i, orientation=False):
        r"""
        Returns the edge the darts ``i`` belongs to.

        If orientation is set to ``True`` then the output is a `2`-tuple
        ``(e,o)`` where ``e`` is the index of the edge and ``o`` is its
        orientation as ``+1`` or ``-1``.
        """
        if self._active_darts[i]:
            if not orientation:
                return self._dart_to_edge_index[i]
            j = self._dart_to_edge_index[i]
            if i == self._edge_cycles[j][0]:
                return (j,1)
            elif i == self._edge_cycles[j][1]:
                return (j,-1)
            else:
                raise ValueError, "this should not happen!"
        raise ValueError, "dart %d is not active" %i

    #
    # Faces
    #

    def num_faces(self):
        r"""
        Return the number of faces.
        """
        return max(1,len(self._face_cycles))

    def face_perm(self):
        r"""
        Return the permutation that defines the face.
        """
        return self._faces

    def face_orbit(self, i):
        r"""
        Return the orbit of ``i`` under the permutation associated to faces.
        """
        if self._active_darts[i]:
            return perm_orbit(self._faces,i)
        return None

    def faces(self):
        r"""
        Return the list of faces.
        """
        return self._face_cycles

    def dart_to_face(self, i):
        if self._active_darts[i]:
            return self._dart_to_face_index[i]
        raise ValueError("dart {} is not active".format(i))

    def dual(self):
        r"""
        Returns the dual Ribbon graph.

        The *dual* ribbon graph of `(v,e,f)` is `(f^{-1}, e, v^{-1})`.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: r = RibbonGraph(edges='(0,1)',faces='(0)(1)'); r
            Ribbon graph with 1 vertex, 1 edge and 2 faces
            sage: r.dual()
            Ribbon graph with 2 vertices, 1 edge and 1 face
        """
        return RibbonGraph(
            vertices=perm_invert(self._faces),
            edges=self._edges,
            faces=perm_invert(self._vertices))

    #
    # euler characteristic
    #

    def euler_characteristic(self):
        r"""
        Returns the Euler characteristic of the embedded surface.

        The *Euler characteristic* of a surface complex is `V - E + F`, where
        `V` is the number of vertices, `E` the number of edges and `F` the
        number of faces.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: r = RibbonGraph(edges='(0,1)(2,3)(4,5)',faces='(0,2,4)(3,5)')
            sage: r.euler_characteristic()
            2

            sage: r = RibbonGraph(edges='(0,1)(2,3)',faces='(0,2,1,3)')
            sage: r.euler_characteristic()
            0
        """
        return self.num_vertices() - self.num_edges() + self.num_faces()

    def is_plane(self):
        r"""
        Returns true if and only if the ribbon graph belongs in a sphere. In
        other words if it has genus 0.

        EXAMPLES::

            sage: from surface_dynamics.all import *

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
        return self.euler_characteristic() == 2

    def is_plane_tree(self):
        r"""
        Returns True if and only if the ribbon graph is a planar tree. In other
        words, it has genus 0 and only one face.

        EXAMPLES::

            sage: from surface_dynamics.all import *

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
        return (self.num_faces() == 1 and self.genus() == 0)

    def is_triangulated(self):
        r"""
        Returns True if the surface is triangulated. In other words, faces
        consist only of the product of 3-cycles.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: r = RibbonGraph(edges='(0,1)(2,3)(4,5)',faces='(0,2,4)(1,5,3)')
            sage: r.is_triangulated()
            True

            sage: r = RibbonGraph(edges='(0,1)(2,3)',faces='(0,2,1,3)')
            sage: r.is_triangulated()
            False
        """
        return all(len(c) == 3 for c in self.faces())

    def genus(self):
        r"""
        Return the genus of the surface associated to this Ribbon graph.

        EXAMPLES::

            sage: from surface_dynamics.all import *

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
        return 1 - self.euler_characteristic()//2

    #
    # cycles and fundamental group
    #

    def spanning_tree(self):
        r"""
        Return a spanning tree

        OUTPUT:

        - spanning tree as a DiGraph

        - remaining edges as 2-tuples ``(i,e[i])``

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: R = RibbonGraph('(1,2,3)','(1,2)(3,4)')
            sage: R
            Ribbon graph with 2 vertices, 2 edges and 2 faces
            sage: T,o = R.spanning_tree()
            sage: T
            Digraph on 2 vertices
            sage: T.edges()
            [(0, 1, (3, 4))]
            sage: o
            [(1, 2)]

            sage: R = RibbonGraph('(1,2,3)(4,5,6)','(1,2)(3,4)(5,6)')
            sage: R
            Ribbon graph with 2 vertices, 3 edges and 3 faces
            sage: T,o = R.spanning_tree()
            sage: T
            Digraph on 2 vertices
            sage: T.edges()
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
            ([1, 3], [2, 4], [5, 7], [6, 8])
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
        Return a ribbon graph callapsed along a spanning tree.

        The resulting graph is on the same surface as the preceding but has only
        one vertex. It could be used twice to provide a polygonal representation
        with one vertex and one face.

        EXAMPLES::

            sage: from surface_dynamics.all import *

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
        for v0,v1,e in spanning_tree.edges():
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
        Return the list of cycles which are boundaries.

        A cycle is a *boundary* if it bounds a face.

        EXAMPLES::

            sage: from surface_dynamics.all import *

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
        Returns a base of oriented cycles of the Ribbon graph modulo boundaries.

        If ``intersection`` is set to True then the method also returns the
        intersection matrix of the cycles.

        EXAMPLES::

            sage: from surface_dynamics.all import *

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
            [[[1, 3]], [[2, 4]], [[5, 7]], [[6, 8]]]

            sage: f = '(0,10,13)(6,17,11)(2,14,7)(15,12,3)(16,20,19)(18,1,9)(4,22,21)(23,8,5)'
            sage: e = tuple((i,i+1) for i in xrange(0,24,2))
            sage: r = RibbonGraph(edges=e,faces=f); r
            Ribbon graph with 2 vertices, 12 edges and 8 faces
            sage: c,m = r.cycle_basis(intersection=True)
            sage: c
            [[(0, 1), [4, 5]], [[8, 9]], [[12, 13]], [[14, 15], (1, 0)]]
            sage: m
            [ 0  1  0  0]
            [-1  0  0  0]
            [ 0  0  0  1]
            [ 0  0 -1  0]
        """
        T,o = self.spanning_tree()

        # build a Ribbon graph with one vertex and one face
        r = self.collapse(T).dual().collapse().dual()
        if T is None:
            return r.edges()

        if intersection:
            c = r.vertices()[0]
            M = len(c)
            I = []

        cycles = [] # the cycles
        for e in r.edges():
            if verbose:
                print "build cycle from edge %s between vertex v0=%d and v1=%d" %(str(e),self.dart_to_vertex(e[0]),self.dart_to_vertex(e[1]))

            # build the branch to the root from v0
            v0 = self.dart_to_vertex(e[0])
            if verbose:
                print " build branch from v0=%d" %v0
            p0 = []
            while v0 != 0:
                v0,_,e0 = T.incoming_edges(v0)[0] # (v_in,v_out,label)
                p0.append(e0)
                if verbose:
                    print " add %d" %v0
            if verbose:
                print " branch is %s" %str(p0)
            # build the branch to the root from v1

            v1 = self.dart_to_vertex(e[1])
            if verbose:
                print " build branch from v1=%d" %v1
            p1 = []
            while v1 != 0:
                v1,_,e1 = T.incoming_edges(v1)[0]
                p1.append(e1)
                if verbose:
                    print " add %d" %v1
            if verbose:
                print " branch is %s" %str(p1)
            # clean the branches by removing common part
            while p0 and p1 and p0[-1] == p1[-1]:
                if verbose:
                    print "find common element",p0[-1]
                p0.pop(-1)
                p1.pop(-1)

            # add the cycle to the list
            cycles.append((p0,e,p1))

            # compute algebraic intersection with preceding cycles
            if intersection:
                i = []
                for _,ee,_ in cycles:
                    if verbose:
                        print "compute intersection"
                    p_in = c.index(e[1])
                    p_out = (c.index(e[0]) - p_in) % M
                    q_in  = (c.index(ee[1]) - p_in) % M
                    q_out = (c.index(ee[0]) - p_in) % M
                    if verbose:
                        print "  after reduction: p_out = %d, q_in = %d, q_out = %d" %(p_out,q_in,q_out)

                    # compute intersection
                    # p_in = 0 and the others 3 are positive
                    if q_in < p_out and p_out < q_out:
                        i.append(1)
                    elif q_out < p_out and p_out < q_in:
                        i.append(-1)
                    else:
                        i.append(0)

                I.append(i)

        # make cycle as list
        cycles = [p0[::-1]+[e]+[c[::-1] for c in p1] for p0,e,p1 in cycles]

        if intersection:
            m = matrix(len(cycles))
            for j in xrange(len(I)):
                for jj in xrange(len(I[j])):
                    m[j,jj] = I[j][jj]
                    m[jj,j] = -I[j][jj]

            return cycles, m
        return cycles

    def is_cycle(self,c):
        r"""
        Test whether ``c`` is a cycle.

        A *path* is a sequence of oriented edges such that each edge starts
        where the preceding one ends. A *cycle* is a path which starts where it
        ends.
        """
        for i in xrange(len(c)-1):
            if self.dart_to_vertex(c[i][1]) != self.dart_to_vertex(c[i+1][0]):
                return False
        if self.dart_to_vertex(c[-1][1]) != self.dart_to_vertex(c[0][0]):
            return False
        return True

class RibbonGraphWithAngles(RibbonGraph):
    r"""
    A Ribbon graph with angles between edges

    Currently angles can only be *rational* multiples of pi.

    TODO:

    - allows any kind of angles by providing a sum for the total and considering
      each angle as a (projective) portion of the total angle.
    """
    def __init__(self, vertices=None, edges=None, faces=None, angles=None):
        r = RibbonGraph(vertices,edges,faces)
        RibbonGraph.__init__(self,r.vertex_perm(),r.edge_perm(),r.face_perm())

        if len(angles) != self.num_darts():
            raise ValueError, "there are %d angles and %d darts" %(len(angles),self.num_darts())
        self._angles = map(QQ,angles)
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
                raise ValueError, "the angle of a face should be (nb_edges - 2) x pi"

    def angle_between_darts(self, d1, d2):
        r"""
        Return the angle between the darts ``d1`` and ``d2``
        """
        v = self.vertex_orbit(d1)
        if d2 not in v:
            raise ValueError, "d1=%s and d2=%s are not at the same vertex" %(str(d1),str(d2))

        a = 0
        i = 0
        while v[i] != d2:
            a += self._angles[v[i]]
            i += 1
        return a

    def angle_at_vertex(self,v):
        r"""
        Angle at a vertex (coefficient of pi)
        """
        return self._total_angle[v]

    def angle_at_vertices(self):
        r"""
        Return the list of angles at a vertex.
        """
        return self._total_angle

    def winding(self, c):
        r"""
        Return winding number along the cycle ``c``.

        This is NOT well defined because it depends on the way we choose to pass
        on the left or on the right at singularity.
        """
        a = 0
        for i in xrange(len(c)-1):
            d1 = c[i][1]
            d2 = c[i+1][0]
            if self.dart_to_vertex(d1) != self.dart_to_vertex(d2):
                raise ValueError, "c is not a cycle"
            a += self.angle_between_darts(d1,d2)-1
        d1 = c[-1][1]
        d2 = c[0][0]
        if self.dart_to_vertex(d1) != self.dart_to_vertex(d2):
            raise ValueError, "c is not a cycle"
        a += self.angle_between_darts(d1,d2)-1

        return a

    def holonomy_representation(self):
        r"""
        Return the holonomy representation in `SO(2)` as two lists.

        The first list correspond to cycles around vertices, while the second
        correspond to a cycle basis that generate homology.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: e = '(0,1)(2,3)'
            sage: f = '(0,2,1,3)'
            sage: a = [1/2,1/2,1/2,1/2]
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.holonomy_representation()
            ([0], [0, 0])

        The standard cube::

            sage: e = tuple((i,i+1) for i in xrange(0,24,2))
            sage: f = '(0,20,7,10)(16,22,19,21)(2,9,5,23)(14,3,17,1)(12,8,15,11)(18,4,13,6)'
            sage: a = [1/2]*24
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.holonomy_representation()
            ([3/2, 3/2, 3/2, 3/2, 3/2, 3/2, 3/2, 3/2], [])

        Two copies of a triangle::

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
        for c in xrange(self.num_vertices()):
            w = self.angle_at_vertex(c)
            l1.append(w - 2*floor(w/2))

        l2 = []
        for c in self.cycle_basis():
            w = self.winding(c)
            l2.append(w - 2*floor(w/2))

        return l1,l2

    def has_trivial_holonomy(self):
        r"""
        Test whether self has trivial holonomy representation

        EXAMPLES::

            sage: from surface_dynamics.all import *

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

    def spin_parity(self,check=True,verbose=False):
        r"""
        Return the spin parity of the Ribbon graph with angles.

        The surface should be holonomy free and with odd multiple of 2 pi
        angles.

        EXAMPLES:

            sage: from surface_dynamics.all import *

        We first consider the case of the torus::

            sage: e = '(0,1)(2,3)'
            sage: f = '(0,2,1,3)'
            sage: a = [1/2,1/2,1/2,1/2]
            sage: r = RibbonGraphWithAngles(edges=e,faces=f,angles=a)
            sage: r.spin_parity()
            1

        Then the case of genus 2 surface (with an angle of 6pi)::

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

        In genus 3 two spin parities occur for one conical angle 10pi::

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
            raise ValueError, "the surface does not have trivial holonomy"
        if any((i+2)%4 for i in self.angle_at_vertices()):
            raise ValueError, "each angle should be odd multiple of 2pi"

        GF2 = GF(2)

        c,M = self.cycle_basis(intersection=True)

        winding = []
        for cc in c:
            w = self.winding(cc)
            if w % 2 != 0:
                raise ValueError, "fatal error ! each winding should be a multiple of 2"
            winding.append(GF2(w//2))

        if verbose:
            print "cycles with winding"
            for i in xrange(len(c)):
                print c[i], winding[i]
            print "intersection matrix on Z"
            print M

        # compute a base change to get a symplectic basis
        _,P = M.symplectic_form()
        M = M.change_ring(GF2)
        P = P.change_ring(GF2)
        if verbose:
            print "base change for symplectic basis on GF(2)"
            print P

        g = self.genus()

        s = GF2(0)
        for i in xrange(g):
            # 1. computation of q(P.row(i))
            a = P.row(i)
            a_indices = [j for j in xrange(2*g) if a[j] != 0]
            ## winding + nb_components
            t_a = sum(winding[i]+1 for i in a_indices)
            ## self intersection
            for j1 in xrange(len(a_indices)):
                for j2 in xrange(j1+1,len(a_indices)):
                    t_a += M[a_indices[j1],a_indices[j2]]

            # 2. computation of q(P.row(g+i))
            b = P.row(g+i)
            b_indices = [j for j in xrange(2*g) if b[j] != 0]
            ## winding + nb_components
            t_b = sum(winding[i]+1 for i in b_indices)
            ## self intersection
            for j1 in xrange(len(b_indices)):
                for j2 in xrange(j1+1,len(b_indices)):
                    t_b += M[b_indices[j1],b_indices[j2]]

            # 3. add to s the contribution of the couple
            if verbose:
                print "contribution from %d is %d * %d = %d"%(i,t_a,t_b,t_a*t_b)
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

class RibbonGraphWithHolonomies(RibbonGraph):
    r"""
    A Ribbon graph with holonomies.

    For now
    """
    def __init__(self, vertices=None, edges=None, faces=None, holonomies=None):
        r = RibbonGraph(vertices,edges,faces)
        RibbonGraph.__init__(self,r.vertex_perm(),r.edge_perm(),r.face_perm())

        if len(holonomies) != self.num_darts():
            raise ValueError, "there are %d angles and %d darts" %(len(angles),self.num_darts())

        from sage.modules.free_module import FreeModule
        V = FreeModule(ZZ,2)
        self._holonomies = map(V, holonomies)

        #self._angles = map(angle, self._holonomies)




