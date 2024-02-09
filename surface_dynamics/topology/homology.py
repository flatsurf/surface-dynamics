r"""
Fat graph homology.
"""
# ****************************************************************************
#       Copyright (C) 2023 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method

from sage.modules.module import Module
from sage.structure.element import ModuleElement, parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.modules import Modules

from sage.modules.free_module import FreeModule
from sage.matrix.constructor import matrix
from sage.arith.misc import gcd
from sage.rings.integer_ring import ZZ

from surface_dynamics.misc.permutation import perm_orbit


def path_to_edge_coefficients(fat_graph, path, check=True):
    if check and not fat_graph.is_path(path):
        raise ValueError('not a valid path')
    coeffs = [0] * fat_graph.num_edges()
    for e in path:
        if e % 2:
            coeffs[e // 2] -= 1
        else:
            coeffs[e // 2] += 1
    return coeffs


class FatGraphHomologyElement(ModuleElement):
    def __init__(self, parent, v):
        ModuleElement.__init__(self, parent)
        self._v = parent._module(v)
        self._v.set_immutable()

    def _richcmp_(self, other, op):
        return self._v._richcmp_(other._v, op)

    def monomial_coefficients(self):
        return self._v.monomial_coefficients()

    def _add_(self, other):
        P = self.parent()
        return P.element_class(P, self._v + other._v)

    def vector(self):
        return self._v

    def _repr_(self):
        if not self._v:
            return '0'

        parent = self.parent()
        fg = parent.fat_graph()
        edge_coeffs= [parent.base_ring().zero()] * fg.num_edges()
        cycles = self.parent()._cycle_basis()[0]
        for coeff, cycle in zip(self._v, cycles):
            if coeff:
                for e in cycle:
                    if e % 2:
                        edge_coeffs[e // 2] -= coeff
                    else:
                        edge_coeffs[e // 2] += coeff

        return '(' + str(edge_coeffs)[1:-1] + ')'

    def intersection(self, other):
        r"""
        Return the algebraic intersection with ``other``.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: fg = FatGraph('(0,13,10)(1,6,8,5,3,12,9,14)(2,15,7,4,11)', '(0,14,2,5,7,1,10,4,8,12)(3,11,13)(6,15,9)')
            sage: H = fg.homology()
            sage: e1 = H.element_from_path([0, 12])
            sage: e2 = H.element_from_path([0, 3, 11])
            sage: e3 = H.element_from_path([13, 14, 11])
            sage: e1.intersection(e2)
            0
            sage: e1.intersection(e3)
            1
            sage: e2.intersection(e3)
            1
        """
        if parent(self) != parent(other):
            raise ValueError
        P = self.parent()
        u = self.vector()
        v = other.vector()
        I = self.parent()._cycle_basis()[1]
        return u * I * v


class FatGraphAbsoluteHomology(UniqueRepresentation, Module):
    r"""
    Absolute homology of a fat graph.

    A canonical basis of homology (in general not symplectic) is provided by a
    tree cotree decomposition of the underlying fat graph. Namely, each
    homology class has a unique representative that excludes edges of the
    cotree.

    TESTS::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: fg = FatGraph('(0,13,10)(1,6,8,5,3,12,9,14)(2,15,7,4,11)', '(0,14,2,5,7,1,10,4,8,12)(3,11,13)(6,15,9)')
            sage: H = fg.homology()
            sage: TestSuite(H).run()
    """
    Element = FatGraphHomologyElement

    def __init__(self, fat_graph, base_ring, tree_cotree_decomposition):
        self._fat_graph = fat_graph.copy(mutable=False)
        self._tree, self._cotree, self._complementary_edges = tree_cotree_decomposition
        nv = fat_graph.num_vertices()
        ne = fat_graph.num_edges()
        nf = fat_graph.num_faces()
        d = self.dimension()
        self._module = FreeModule(base_ring, d)
        self._edge_module = FreeModule(base_ring, ne)
        self._vertex_module = FreeModule(base_ring, nv)
        self._face_module = FreeModule(base_ring, nf)
        Module.__init__(self, base_ring, category=Modules(base_ring).FiniteDimensional().WithBasis())

    def _repr_(self):
        r"""
        TESTS::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: fg = FatGraph('(0,13,10)(1,6,8,5,3,12,9,14)(2,15,7,4,11)', '(0,14,2,5,7,1,10,4,8,12)(3,11,13)(6,15,9)')
            sage: fg.homology()
            Homology(FatGraph('(0,13,10)(1,6,8,5,3,12,9,14)(2,15,7,4,11)', '(0,14,2,5,7,1,10,4,8,12)(3,11,13)(6,15,9)'); Integer Ring)
        """
        return 'Homology({}; {})'.format(self.fat_graph(), self.base_ring())

    def _cotree_dfs(self):
        r"""
        Return a dfs ordering of the edges in the cotree.

        TESTS::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: fg = FatGraph('(0,13,10)(1,6,8,5,3,12,9,14)(2,15,7,4,11)', '(0,14,2,5,7,1,10,4,8,12)(3,11,13)(6,15,9)')
            sage: fg.homology()._cotree_dfs()
            [3, 15]
        """
        cotree = self._cotree
        fl = self._fat_graph._fl
        nf = self._fat_graph._nf
        children = [[] for _ in range(nf)]
        for f, e in enumerate(cotree):
            if e == -1:
                root = f
            else:
                assert fl[e] == f
                children[fl[e ^ 1]].append(e)
        i = 1
        res = children[root]
        while i < len(res):
            res.extend(children[fl[res[i]]])
            i += 1
        return res

    def base_ring(self):
        r"""
        Return the base ring.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: fg = FatGraph('(0,13,10)(1,6,8,5,3,12,9,14)(2,15,7,4,11)', '(0,14,2,5,7,1,10,4,8,12)(3,11,13)(6,15,9)')
            sage: fg.homology().base_ring()
            Integer Ring

            sage: fg.homology(Zmod(2)).base_ring()
            Ring of integers modulo 2
        """
        return self._module.base_ring()

    def change_ring(self, base_ring):
        r"""
        Change the underlying ring.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: fg = FatGraph('(0,13,10)(1,6,8,5,3,12,9,14)(2,15,7,4,11)', '(0,14,2,5,7,1,10,4,8,12)(3,11,13)(6,15,9)')
            sage: fg.homology().change_ring(Zmod(2)).base_ring()
            Ring of integers modulo 2
        """
        if base_ring == self.base_ring():
            return self
        return FatGraphAbsoluteHomology(self._fat_graph, base_ring, (self._tree, self._cotree, self._complementary_edges))

    def _cycle_basis(self):
        return self._fat_graph.cycle_basis(intersection=True, tree_cotree_decomposition=(self._tree, self._cotree, self._complementary_edges))

    def fat_graph(self):
        r"""
        Return the underlying graph.
        """
        return self._fat_graph

    def an_element(self):
        r"""
        Return an element in this homology group.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: fg = FatGraph('(0,13,10)(1,6,8,5,3,12,9,14)(2,15,7,4,11)', '(0,14,2,5,7,1,10,4,8,12)(3,11,13)(6,15,9)')
            sage: H = fg.homology()
            sage: H.an_element()
            (-1, 0, 1, 0, 0, 1, 0, 0)
        """
        return self.element_class(self, self._module.an_element())

    def random_element(self, *args, **kwds):
        r"""
        Return a random element in this homology group.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: fg = FatGraph('(0,13,10)(1,6,8,5,3,12,9,14)(2,15,7,4,11)', '(0,14,2,5,7,1,10,4,8,12)(3,11,13)(6,15,9)')
            sage: H = fg.homology()
            sage: H.random_element()  # random
            (-2, 0, -1, -1, -12, 0, -2, 0)
        """
        return self.element_class(self, self._module.random_element(*args, **kwds))

    def __iter__(self):
        return (self.element_class(self, v) for v in self._module)

    def dimension(self):
        r"""
        Return the dimension which is twice the genus of the underlying graph.
        """
        return 2 * self._fat_graph.genus()

    def cycle_basis(self):
        r"""
        Return the basis as simple paths in the underlying fat graph.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: fg = FatGraph('(0,13,10)(1,6,8,5,3,12,9,14)(2,15,7,4,11)', '(0,14,2,5,7,1,10,4,8,12)(3,11,13)(6,15,9)')
            sage: H = fg.homology()
            sage: H.cycle_basis()
            [[10, 4, 1], [0, 6, 11], [8], [0, 12]]
        """
        return self._cycle_basis()[0]

    def intersection_matrix(self):
        r"""
        Return the intersection matrix on the canonical basis.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: fg = FatGraph('(0,13,10)(1,6,8,5,3,12,9,14)(2,15,7,4,11)', '(0,14,2,5,7,1,10,4,8,12)(3,11,13)(6,15,9)')
            sage: H = fg.homology()
            sage: H.intersection_matrix()
            [ 0  1  1  0]
            [-1  0  0  0]
            [-1  0  0  1]
            [ 0  0 -1  0]
        """
        return self._cycle_basis()[1]

    def basis(self):
        r"""
        Return the canonical basis built from a tree-cotree decomposition of the underlying fat graph.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: fg = FatGraph('(0,13,10)(1,6,8,5,3,12,9,14)(2,15,7,4,11)', '(0,14,2,5,7,1,10,4,8,12)(3,11,13)(6,15,9)')
            sage: H = fg.homology()
            sage: B = H.basis()
            sage: B
            ((-1, 0, 1, 0, 0, 1, 0, 0),
             (1, 0, 0, 1, 0, -1, 0, 0),
             (0, 0, 0, 0, 1, 0, 0, 0),
             (1, 0, 0, 0, 0, 0, 1, 0))
            sage: matrix(ZZ, 4, [B[i].intersection(B[j]) for j in range(4) for i in range(4)])
            [ 0 -1 -1  0]
            [ 1  0  0  0]
            [ 1  0  0 -1]
            [ 0  0  1  0]
        """
        return tuple(self.element_class(self, b) for b in self._module.basis())

    def symplectic_basis(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: fg = FatGraph('(0,13,10)(1,6,8,5,3,12,9,14)(2,15,7,4,11)', '(0,14,2,5,7,1,10,4,8,12)(3,11,13)(6,15,9)')
            sage: H = fg.homology()
            sage: B = H.symplectic_basis()
            sage: B
            ((-1, 0, 1, 0, 0, 1, 0, 0),
             (-1, 0, 0, -1, 1, 1, 0, 0),
             (1, 0, 0, 1, 0, -1, 0, 0),
             (1, 0, 0, 0, 0, 0, 1, 0))
            sage: matrix(ZZ, 4, [B[i].intersection(B[j]) for j in range(4) for i in range(4)])
            [ 0  0 -1  0]
            [ 0  0  0 -1]
            [ 1  0  0  0]
            [ 0  1  0  0]
        """
        # the rows of C form a symplectic basis
        # F = C * self * C.transpose()
        F, C = self._cycle_basis()[1].symplectic_form()
        return tuple(self.element_class(self, b) for b in C.rows())

    def element_from_path(self, path):
        r"""
        Return an homology element from a path.

        The ``path`` should be given as a sequence of consecurive half edges.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: fg1 = FatGraph('(0,13,10)(1,6,8,5,3,12,9,14)(2,15,7,4,11)', '(0,14,2,5,7,1,10,4,8,12)(3,11,13)(6,15,9)')
            sage: H1 = fg1.homology()
            sage: H1.element_from_path([10, 4, 1])
            (-1, 0, 1, 0, 0, 1, 0, 0)

            sage: fg2 = FatGraph('(0,13,11,9,8,14,2,1,3)(4,15,6,5,12,7,10)', '(0,2,1,3,14,4,6,12)(5,10,13)(7,15,8,11)(9)')
            sage: H2 = fg2.homology()
            sage: H2.element_from_path([10, 13])
            (0, 0, 0, 0, 0, 1, -1, 0)
            sage: H2.element_from_path([0, 11, 4, 15])
            (1, 0, 0, 1, 0, 1, -1, 0)

        Faces are mapped to zero::

            sage: assert all(H1.element_from_path(face).is_zero() for face in fg1.faces())
            sage: assert all(H2.element_from_path(face).is_zero() for face in fg2.faces())
        """
        coeffs = path_to_edge_coefficients(self._fat_graph, path)
        return self.element_from_edge_coefficients(coeffs)

    def element_from_half_edge_coefficients(self, coeffs):
        r"""
        Return an element from half edge coefficients.

        The argument ``coeffs`` must be a list or a vector whose length is the
        number of half edges in the underlying fat graph. The entries are the
        coefficients for each half edge.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: fg = FatGraph('(0,13,10)(1,6,8,5,3,12,9,14)(2,15,7,4,11)', '(0,14,2,5,7,1,10,4,8,12)(3,11,13)(6,15,9)')
            sage: H = fg.homology()
            sage: H.element_from_half_edge_coefficients((1, 1, 3, 0, 1, 3, 0, 0, -2, 0, 1, 0, 1, 0, 0, 0))
            (0, 0, -2, 0, -2, -2, -2, 0)
        """
        if len(coeffs) != 2 * self._fat_graph.num_edges():
            raise ValueError('invalid coefficients')
        edge_coeffs = [self.base_ring().zero()] * self._fat_graph.num_edges()
        for e in range(self._fat_graph.num_edges()):
            edge_coeffs[e] += coeffs[2 * e] - coeffs[2 * e + 1]
        return self.element_from_edge_coefficients(edge_coeffs)

    def boundary(self, coeffs):
        coeffs = self._edge_module(coeffs)
        bdry = [0] * self._fat_graph._nv
        vl = self._fat_graph._vl
        for e, coeff in enumerate(coeffs):
            if coeff:
                bdry[vl[2 * e]] += coeff
                bdry[vl[2 * e + 1]] -= coeff
        return self._vertex_module(bdry)

    def element_from_edge_coefficients(self, coeffs):
        r"""
        Return an element from edge coefficients.

        The argument ``coeffs`` must be a list or a vector whose length is the
        number of edges in the underlying fat graph. The entries are the
        coefficients for each edge.

        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: fg = FatGraph('(0,13,10)(1,6,8,5,3,12,9,14)(2,15,7,4,11)', '(0,14,2,5,7,1,10,4,8,12)(3,11,13)(6,15,9)')
            sage: H = fg.homology()
            sage: H.element_from_edge_coefficients((1, 1, 0, 1, 0, 0, 1, 0))
            (1, 0, 0, 1, 0, -1, 0, 0)
            sage: H.element_from_edge_coefficients((1, 0, 0, 1, 0, -1, 0, 0))
            (1, 0, 0, 1, 0, -1, 0, 0)
        """
        if len(coeffs) != self._fat_graph.num_edges():
            raise ValueError('invalid coefficients')
        coeffs = self._edge_module(coeffs)
        if self.boundary(coeffs):
            raise ValueError('not a cycle')

        # step one: rewrite edge in the cotree using faces
        for e in self._cotree_dfs():
            c = coeffs[e // 2]
            if c:
                coeffs[e // 2] = 0
                face = perm_orbit(self._fat_graph._fp, e)
                assert face[0] == e
                for i in range(1, len(face)):
                    ee = face[i]
                    if ee % 2 == e % 2:
                        coeffs[ee // 2] -= c
                    else:
                        coeffs[ee // 2] += c

        assert all(coeffs[e // 2] == 0 for e in self._cotree if e != -1), (coeffs, self._cotree, self._cotree_dfs())
        assert self.boundary(coeffs) == 0

        # step two: the remaining part is a sum of cycles and one
        # obtain the coefficients by reading the complementary edges
        v = [coeffs[e // 2] for e in self._complementary_edges]
        return self.element_class(self, v)

    def element_from_vector(self, v):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph import FatGraph
            sage: fg = FatGraph('(0,13,10)(1,6,8,5,3,12,9,14)(2,15,7,4,11)', '(0,14,2,5,7,1,10,4,8,12)(3,11,13)(6,15,9)')
            sage: H = fg.homology()
            sage: h = H.element_from_vector([0, 0, 0, 1])
            sage: h
            (1, 0, 0, 0, 0, 0, 1, 0)
            sage: h.vector()
            (0, 0, 0, 1)
        """
        return self.element_class(self, v)
