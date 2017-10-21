"""
Fatgraphs

Fatgraph (or ribbon graph, or combinatorial map) is a map embedded in an orientable
surface. We encode a fatgraph with n edges with a permutation on
`X_n = \{-n, -n+1, \ldots, n-1\}`.  Each element of `X_n` is called a dart (or
half edge). An edge is a pair `(i, -i-1)` identified to its positive part.

The paths along the edges are naturally identified to elements
of the Free group with n generators labeled 0 to n-1.


BIG QUESTION:

- how do we generate fatgraphs with fixed num verts, num edges, num faces?
- is it easier if we only fix num edges, num faces?
"""

from array import array
from surface_dynamics.misc.even_permutation import *


class FatGraph(object):
    r"""
    Fatgraph

    EXAMPLES:

    The sphere obtained from two bigons::

        sage: from surface_dynamics import *
        sage: FatGraph([1,2,0,-2,-1,-3])
        Fatgraph of genus 0 with 3 vertices, 3 edges and 2 faces

    The torus obtained from gluing the opposite sides of a square::

        sage: FatGraph([1,-1,0,-2])
        Fatgraph of genus 1 with 1 vertex, 2 edges and 1 face
    """
    __slots__ = ['_n', '_v_perm', '_f_perm', '_vertices', '_dart_to_vertex', '_faces', '_dart_to_face', '_can']

    def __init__(self, f=None, v=None):
        if v is None and f is None:
            self._n = 0
            self._v_perm = []
            self._f_perm = []
            self._vertices = []
            self._dart_to_vertex = []
            self._faces = []
            self._dart_to_face = []
            self._can = None
            return

        if v is not None:
            self._v_perm = array('l', v)
            even_perm_check(self._v_perm)
            self._n = len(v) // 2
        if f is not None:
            self._f_perm = array('l', f)
            even_perm_check(self._f_perm)
            self._n = len(f) // 2

        e = even_perm_minus(self._n)
        if v is None:
            self._v_perm = even_perm_compose_i(self._f_perm, e)
        elif f is None:
            self._f_perm = even_perm_compose_i(e, self._v_perm)

        assert len(self._v_perm) == len(self._f_perm)
        assert even_perm_compose(even_perm_compose(self._v_perm, e), self._f_perm) == even_perm_identity(self._n)
        assert even_perm_is_transitive(self._v_perm)

        self._vertices, self._dart_to_vertex = even_perm_cycles(self._v_perm)
        self._faces, self._dart_to_face = even_perm_cycles(self._f_perm)
        self._can = None

    def __repr__(self):
        r"""
        """
        m = self.num_verts()
        verts = '{} vert{}'.format(m, 'ex' if m == 1 else 'ices')
        m = self.num_edges()
        edges = '{} edge{}'.format(m, '' if m == 1 else 's')
        m = self.num_faces()
        faces = '{} face{}'.format(m, '' if m == 1 else 's')
        return "Fatgraph of genus {} with {}, {} and {}".format(
                self.genus(), verts, edges, faces)

    def __copy__(self):
        F = FatGraph()
        F._n = self._n
        F._v_perm = self._v_perm[:]
        F._f_perm = self._f_perm[:]
        F._vertices = self._vertices[:]
        F._faces = self._faces[:]
        F._dart_to_vertex = self._dart_to_vertex[:]
        F._dart_to_face = self._dart_to_face[:]
        F._can = self._can

    def _check_composition(self):
        for i in range(-self._n, self._n):
            if self._f_perm[~self._v_perm[i]] != i:
                raise ValueError('problem at i = {}'.format(i))

    def euler_characteristic(self):
        return self.num_verts() - self.num_edges() + self.num_faces()

    def num_verts(self):
        return len(self._vertices)

    def num_edges(self):
        return self._n

    def num_faces(self):
        return len(self._faces)

    def faces(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics import *
            sage: T = FatGraph([1r,2r,0r,~4r,~5r,~3r,~2r,~1r,~0r,4,3,5])
            sage: T.faces()
            [(-6, -3, 4), (-5, -2, 3), (-4, -1, 5), (0, 1, 2)]
        """
        return self._faces

    def genus(self):
        return 1 - self.euler_characteristic()//2

    def vertex_profile(self):
        return map(len, self._vertices)

    def face_profile(self):
        return map(len, self._faces)

    def dual(self):
        r"""
        The dual fatgraph (vertices become faces and vice-versa)
        """
        raise NotImplementedError

    def relabel(self, m):
        r"""
        EXAMPLES::

            sage: from surface_dynamics import *
            sage: T = FatGraph([1r,2r,0r,~4r,~5r,~3r,~2r,~1r,~0r,4,3,5])
            sage: m = [1,0,2,3,4,5,~5r,~4r,~3r,~2r,~0r,~1r]
        """
        if not is_signed_perm(m):
            raise ValueError('m should be a signed permutation')
        even_perm_relabel(self._v_perm, m)
        even_perm_relabel(self._f_perm, m)
        self._vertices, self._dart_to_vertex = even_perm_cycles(self._v_perm)
        self._faces, self._dart_to_face = even_perm_cycles(self._f_perm)

    def canonical_label(self, return_map=False, return_graph=True):
        r"""
        EXAMPLES:

        The tetrahedron::

            sage: from surface_dynamics import *
            sage: T = FatGraph([1r,2r,0r,~4r,~5r,~3r,~2r,~1r,~0r,4,3,5])
        """
        if self._can is not None:
            return self._can

        f,m = even_perm_canonical_lael(self._f_perm)
        self._can = FatGraph(f)
        return (self._can, m) if return_map else g

def fatgraphs(vertex_profile, face_profile):
    r"""
    Run through the list of fatgraphs with the given vertex and face profiles.
    """
    raise NotImplementedError

    from sage.groups.perm_gps.symgp_conjugacy_class import conjugacy_class_iterator
    from surface_dynamics.misc.permutation import canonical_perm, perm_compose

    n = sum(vertex_profile)
    if n%2 or sum(face_profile) != n:
        raise ValueError

    e = canonical_perm([2]*(n//2))

