r"""
Origami

An origami is:

- a couple of permutations

- a covering of the torus ramified over one point

- gluing of squares

- Abelian differential on Riemann surface with rational periods

EXAMPLES::

    sage: from surface_dynamics.all import *

    sage: E = origamis.EierlegendeWollmilchsau()
    sage: E.r()
    (1,2,3,4)(5,6,7,8)
    sage: E.u()
    (1,5,3,7)(2,8,4,6)
    sage: E.stratum_component()
    H_3(1^4)^c
    sage: E.lyapunov_exponents_approx()   # abs tol 1e-3
    [0.0000485630931783940, 0.0000452662733371477]

    sage: o = Origami('(1,2,3,4,5,6)','(1,7)')
    sage: V = o.veech_group()
    sage: V
    Arithmetic subgroup of index 54
    sage: G = V.coset_graph()
    sage: G.diameter()
    16
"""
from surface_dynamics.flat_surfaces.origamis.origami_dense import Origami_dense_pyx

from sage.structure.sage_object import SageObject
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.groups.perm_gps.permgroup import PermutationGroupElement
from sage.misc.cachefunc import cached_method
from copy import copy
from sage.matrix.constructor import matrix, identity_matrix

from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets

from sage.rings.integer import Integer
from sage.modules.free_module import VectorSpace
from sage.rings.rational_field import QQ
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.unique_representation import UniqueRepresentation

from sage.interfaces.gap import gap

from sage.plot.plot import options

def flatten_word(w):
    l = []
    for i,j in w:
        l.append(i*j)
    return ''.join(l)

def permutation_simplicial_action(r,u,n,w):
    r"""
    From a word in 'l','r' return the simplicial action on homology as
    well as the obtained origami.

    INPUT:

    - ``r`` and ``u`` - permutations

    - ``n`` - the degree of the permutations

    - ``w`` - a string in 'l' and 'r' or a list of 2-tuples which are made of
      the letter 'l' or 'r' and a positive integer

    """
    if w is None:
        w = []
    elif isinstance(w,list):
        w = flatten_word(w)

    res = identity_matrix(2*n)

    for letter in reversed(w):
        if letter == 'l':
            u = u*~r
            m = identity_matrix(2*n)
            m[:n,n:] = u.matrix()
            res = m * res
        elif letter == 'r':
            r = r*~u
            m = identity_matrix(2*n)
            m[n:,:n] = r.matrix()
            res = m * res
        else:
            raise ValueError, "does not understand the letter %s" %str(letter)

    return r,u,res

#
# Origami and pillow case cover constructors
#

def Origami(r, u,
        sparse=False,
        check=True,
        as_tuple=False,
        positions=None, name=None):
    r"""
    Constructor for origami

    INPUT:

    - ``r``, ``u`` - two permutations

    - ``sparse`` - boolean (default: False)

    - ``check`` - boolean (default: True) - whether or not check the input

    - ``as_tuple`` - boolean (default: False) - if True, assume that ``r`` and
      ``u`` are tuples on [0,...,N-1] (time efficient)

    - ``positions`` - list of 2-tuples (default: None) - position of the squares
      for drawings

    - ``name`` - an optional name to the origami

    """
    if not as_tuple:
        r = PermutationGroupElement(r, check=check)
        u = PermutationGroupElement(u, check=check)

        r = [i-1 for i in r.domain()]
        u = [i-1 for i in u.domain()]

        N = max(len(r),len(u))
        r.extend(xrange(len(r),N))
        u.extend(xrange(len(u),N))

    elif check:
        sr = set(r)
        su = set(u)
        N = len(r)
        if len(u) != N:
            raise ValueError("the two tuples must be of the same length")
        for i in xrange(N):
            if not i in sr:
                raise ValueError("%d is not in r=%s" %(i,str(r)))
            if not i in su:
                raise ValueError("%d is not in u=%s" %(i,str(u)))

    o = Origami_dense_pyx(tuple(r),tuple(u))

    if check and not o.is_connected():
        print "Warning: the origami is not connected"

    if name is not None:
        o.rename(name)
    if positions is not None:
        o.set_positions(positions)
    return o

#
# Other stuff which should be moved in a more geometrical place (like CW
# complex or embeeded graph or whatever).
#

from sage.categories.action import Action

class ActionOnOrigamiObjects(Action):
    r"""
    Generic action of the automorphism group of an origami.
    """
    def __init__(self, objects):
        import operator
        o = objects.origami()
        # Action.__init__(G,S,is_left,op)
        Action.__init__(
                self,
                o.automorphism_group(),
                objects,
                True,
                operator.mul)

    def _call_(self, g, x):
        r"""
        Action of g on x.
        """
        return x._acted_upon_(g, True)

class OrigamiObjects(Parent):
    def chain_space(self):
        return VectorSpace(QQ, self.cardinality())

    @cached_method
    def cycle_space(self):
        r"""
        The space of cycles.
        """
        return self.derivative().right_kernel()

    @cached_method
    def border_space(self):
        r"""
        The border space.
        """
        return self.derivative().column_space()

class OrigamiFaces(OrigamiObjects):
    def __init__(self, origami):
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._origami = origami

    def cardinality(self):
        return self._origami.nb_squares()

    def derivative(self):
        n = self._origami.nb_squares()
        ri = self._origami.r_inv_tuple()
        ui = self._origami.u_inv_tuple()
        der = matrix(2*n,n,ring=QQ)
        for i in xrange(n):
            if ui[i] != i: # horiz sides
                der[i,i] = -1
                der[ui[i],i] = 1
            if ri[i] != i: # vert sides
                der[n+i,i] = 1
                der[n+ri[i],i] = -1
        return der

class OrigamiEdges(OrigamiObjects):
    r"""
    """
    def __init__(self, origami):
        Parent.__init__(self)
        self._origami = origami
        ri = origami.r_inv_tuple()
        ui = origami.u_inv_tuple()
        n = origami.nb_squares()

        self._starts = [None]*(2*n)
        self._ends = [None]*(2*n)
        for i in xrange(n):
            self._ends[i] = origami.vertex(i+1).index()
            self._starts[i] = origami.vertex(ri[i]+1).index()

            self._ends[n+i] = origami.vertex(i+1).index()
            self._starts[n+i] = origami.vertex(ui[i]+1).index()

    def origami(self):
        return self._origami

    # about elements

    def _element_constructor(self, i):
        if i >= 0 and i < 2*o.nb_squares():
            return o.chain_complex().chain_space(1).gen(i)

    def start(self,i):
        return self._starts[i]

    def end(self,i):
        return self._ends[i]

    # about global derivation

    def cardinality(self):
        return 2*self._origami.nb_squares()

    def derivative(self):
        m = self._origami.nb_vertices(True)
        n = 2*self._origami.nb_squares()
        der = matrix(m,n,ring=QQ)
        for i in xrange(n):
            if self._starts[i] != self._ends[i]:
                der[self._ends[i], i] = 1
                der[self._starts[i], i] = -1
        return der

    # geometric stuff (should be moved in cylinder diagrams)

    def is_simple_closed_curve(self, c):
        r"""
        Test if c is a simple closed curve
        """
        assert(c in self.chain_space())
        d = c.dict()
        for i in d:
            if d[i] != 1 and d[i] != -1:
                return False

        vertices = set([])
        for i in d:
            if d[i] == 1:
                v = self.end(i)
            else:
                v = self.start(i)
            if v in vertices:
                return False

            vertices.add(v)

        return True

    def simple_closed_curve_to_vertex_in_out(self, c):
        r"""
        From a curve c (as a vector) returns two dictionnaries associated to
        input/output view from each visited vertices
        """
        c_in = {}
        c_out = {}
        for i,j in c.dict().iteritems():
            if j == 1:
                c_in[self.start(i)] = i
                c_out[self.end(i)] = i
            elif j == -1:
                c_in[self.end(i)] = i
                c_out[self.start(i)] = i
            else:
                raise ValueError, "not a simple closed curve"
        return c_in, c_out

    def basis_of_simple_closed_curves(self):
        from sage.graphs.digraph import DiGraph
        from sage.all import randint

        n = self._origami.nb_squares()
        C = self.chain_space()
        G = DiGraph(multiedges=True,loops=True,implementation='c_graph')

        for i in xrange(2*n):
            G.add_edge(self._starts[i], self._ends[i], i)

        waiting = [0]
        gens = []
        reps = [None] * G.num_verts()
        reps[0] = C.zero()

        while waiting:
            x = waiting.pop(randint(0,len(waiting)-1))
            for v0,v1,e in G.outgoing_edges(x):
                if reps[v1] is not None:
                    gens.append(reps[v0] + C.gen(e) - reps[v1])
                else:
                    reps[v1] = reps[v0] + C.gen(e)
                    waiting.append(v1)

        return gens

    def intersection(self, c1, c2):
        r"""
        Returns the intersection of c1 and c2 assuming that they are simple
        closed curves
        """
        assert(self.is_simple_closed_curve(c1) and self.is_simple_closed_curve(c2))

        c1_in, c1_out = self.simple_closed_curve_to_vertex_in_out(c1)
        c2_in, c2_out = self.simple_closed_curve_to_vertex_in_out(c2)

        intersection = 0
        for vert in c1_in:
            if vert in c2_in:
                if c1_in[vert] == c2_in[vert]:
                    if c1_out[vert] == c2_out[vert]:
                        pass
                    elif are_cyclically_ordered(
                            c1_in[vert], c1_out[vert], c2_out[vert]):
                        intersection += 1

                elif c1_out[vert] == c2_out[vert]:
                    if are_cyclically_ordered(
                            c1_out[vert], c2_in[vert], c1_in[vert]):
                        intersection -= 1

                elif are_cyclically_ordered(
                        c1_in[vert],c2_in[vert],c1_out[vert]):
                    if are_cyclically_ordered(
                            c1_out[vert],c2_out[vert],c1_in[vert]):
                            intersection += 1

                elif are_cyclically_ordered(
                        c1_in[vert],c2_out[vert],c1_out[vert]):
                    intersection -= 1

        return intersection

    def winding(self, c):
        r"""
        Return the winding of the curve c
        """
        assert(self.is_simple_closed_curve(c))
        pass

class OrigamiVertices(OrigamiObjects):
    r"""
    The set of vertices of an origami.

    It is in bijection with the biclasses H  G  C where H is the stabilizer of
    1 and C is the subgroup generated by the commutator ``r u r**-1 u **-1``.
    """
    def __init__(self, origami, register_automorphism_action=False):
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._origami = origami

        self._vertices = []
        self._vertices_inv = {}
        self._vertex_from_ur = {}
        self._vertex_from_dl = {}

        r = origami.r()
        u = origami.u()
        ru = r*u
        vperm = ru*~r*~u

        for ur in vperm.cycle_tuples(singletons=True):
            dl = map(ru,ur)
            v = OrigamiVertex(self, ur, tuple(dl))
            self._vertices_inv[v] = len(self._vertices)
            self._vertices.append(v)
            for i in ur:
                self._vertex_from_ur[i] = v
            for i in dl:
                self._vertex_from_dl[i] = v

        if register_automorphism_action:
            self.register_action(ActionOnOrigamiObjects(self))

    def origami(self):
        r"""
        Return the underlying origami.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: o = Origami('(1,2,3,4)','(1,5)')
            sage: V = o._vertices()
            sage: V.origami()
            (1,2,3,4)(5)
            (1,5)(2)(3)(4)
            sage: V.origami() is o
            True
        """
        return self._origami

    def cardinality(self):
        r"""
        The number of vertices.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: o = Origami('(1,2,3,4)','(1,5)')
            sage: V = o._vertices()
            sage: V.cardinality()
            3
        """
        return Integer(len(self._vertices))

    def chain_space(self):
        r"""
        Return the space of chain on this finite set of vertices.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: o = Origami('(1,2,3,4)', '(1,5)')
            sage: V = o._vertices()
            sage: V.chain_space()
            Vector space of dimension 3 over Rational Field
        """
        return VectorSpace(QQ, self.cardinality())

    def derivative(self):
        r"""
        Return the derivative matrix.

        That is a matrix from the chain space to `\QQ`.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: o = Origami('(1,2,3,4)', '(1,5)')
            sage: V = o._vertices()
            sage: V.derivative()
            [1 1 1]
        """
        return matrix([1]*len(self._vertices), ring=QQ)

    def _repr_(self):
        return "Vertices of origami\n%s" %self.origami()

    def __contains__(self, v):
        return v in self._vertices

    def up_right_vertex(self, i):
        r"""
        Return the vertex that is in the up right corner of the ``i``-th square.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: o = Origami('(1,2,3,4)', '(1,5)')
            sage: V = o._vertices()
            sage: V.up_right_vertex(2)
            vertex (2,)
        """
        return self._vertex_from_ur[i]

    def down_left_vertex(self, i):
        r"""
        Return the vertex that is in the down left corner of the ``i``-th square.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: o = Origami('(1,2,3,4)', '(1,5)')
            sage: V = o._vertices()
            sage: V.down_left_vertex(2)
            vertex (1, 5, 4)
        """
        return self._vertex_from_dl[i]

    vertex = up_right_vertex
    _element_constructor_ = up_right_vertex

    def vertex_index(self, v):
        r"""
        Return the index of the vertex ``v``.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: o = Origami('(1,2,3,4)', '(1,5)')
            sage: V = o._vertices()
            sage: v0 = V.up_right_vertex(2)
            sage: v0
            vertex (2,)
            sage: V.vertex_index(v0)
            1
            sage: V[1] == v0
            True

            sage: v1 = V.up_right_vertex(5)
            sage: v1
            vertex (1, 5, 4)
            sage: V.vertex_index(v1)
            0
            sage: V[0] == v1
            True
        """
        return self._vertices_inv[v]

    def __getitem__(self, i):
        r"""
        Access to a given element.
        """
        return self._vertices[i]

    def __iter__(self):
        r"""
        Iterator over the vertices.
        """
        return iter(self._vertices)

class OrigamiVertex(Element):
    def __init__(self, parent, ur, dl):
        Element.__init__(self, parent)
        self._ur = ur
        self._dl = dl

    def _repr_(self):
        return "vertex %s" %(str(self._ur))

    def __str__(self):
        return str(self._ur)

    def __hash__(self):
        return hash((self._ur,self._dl))

    def __eq__(self, other):
        return (type(self) == type(other) and
                self.parent() == other.parent() and
                self._ur == other._ur)

    def __ne__(self, other):
        return (type(self) != type(other) or
                self.parent() != other.parent() or
                self._ur != other._ur)

    def _acted_upon_(self, g, on_left):
        r"""
        Return the action of the element g of the automorphism group of the
        origami on self.
        """
        if not on_left:
            return self.parent().up_right_vertex(g(self._ur[0]))
        else:
            return self.parent().up_right_vertex(g(self._ur[0]))

    def index(self):
        return self.parent().vertex_index(self)

    def degree(self):
        return len(self._ur)-1

    def up_right_tuple(self):
        return self._ur

    def down_left_tuple(self):
        return self._dl

    # adjacent vectors

    def adjacent_edge_indices(self):
        r"""
        Returns a 2-tuple (incoming,outgoing)
        """
        return self.incoming_edge_indices(), self.outgoing_edge_indices()

    def outgoing_edge_indices(self):
        n = self.parent().origami().nb_squares()
        ri = ~self.parent().origami().r()
        ui = ~self.parent().origami().u()
        res = []
        for i in self._dl:
            res.append(ui(i)-1)
            res.append(n+ri(i)-1)
        return res

    def incoming_edge_indices(self):
        res = []
        n = self.parent().origami().nb_squares()
        for i in self._ur:
            res.append(i-1)
            res.append(n+i-1)
        return res

    @cached_method
    def edge_positions(self):
        r"""
        The position of the edges.
        """
        d_in = {}
        d_out = {}
        e_in = self.incoming_edge_indices()
        e_out = self.outgoing_edge_indices()
        k = 0
        for i in xrange(0,len(e_in),2):
            d_in[e_in[i]] = k
            d_in[e_in[i+1]] = k+1
            d_out[e_out[i]] = k+2
            d_out[e_out[i+1]] = k+3
            k += 4
        return d_in, d_out

    def incoming_edge_position(self, i):
        return self.edge_positions()[0][i]

    def outgoing_edge_positions(self, i):
        return self.edge_positions()[1][i]

class OrigamiChainComplex(SageObject, UniqueRepresentation):
    r"""
    Chain complex for reduced homology of the origami
    """
    def __init__(self, origami):
        self._origami = origami
        self._objects = [OrigamiVertices(origami), OrigamiEdges(origami), OrigamiFaces(origami)]

    def origami(self):
        r"""
        Return the underlying origami.
        """
        return self._origami

    def chain_space(self, degree):
        r"""
        Chain space.

        INPUT:

        - ``degree`` -- either

        """
        assert(degree > -2 and degree < 3)
        if degree == -1:
            return VectorSpace(QQ,1)
        return self._objects[degree].chain_space()

    def cycle_space(self, degree):
        r"""
        Returns the space of cycles of degree i

        Zi = ker(der: C_i -> C_{i-1})
        """
        assert(degree > -2 and degree < 3)
        if degree == -1:
            return VectorSpace(QQ,1)
        return self._objects[degree].cycle_space()

    def border_space(self, degree):
        r"""
        Returns the space of borders of degree i

        Bi = im(der: C_{i+1} -> Ci)
        """
        assert(degree > -2 and degree < 3)
        if degree == 2:
            return self._objects[2].chain_space().subspace([])
        return self._objects[degree+1].border_space()

    def _repr_(self):
        return "Chain complex of origami\n%s" %(self.origami())
