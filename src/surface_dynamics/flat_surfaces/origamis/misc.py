r"""
Old code that might moved or be removed.
"""
#
#from sage.structure.sage_object import SageObject
#
#@cached_function
#def _fact_and_irr(k):
#    r"""
#    intermediate function to accelerate the computation of
#    number_of_irreducible_permutations
#
#    the length of the result is k-1
#
#    [p(j)*factorial(k-j)]
#
#    TESTS::
#
#        sage: _fact_and_irr(1) == ()
#        True
#        sage: _fact_and_irr(2) == (1,)
#        True
#        sage: _fact_and_irr(3) == (2,1)
#        True
#    """
#    if k == 1: return ()
#    rec = _fact_and_irr(k-1)
#    return (number_of_irreducible_permutations(k-1),) + tuple((j+2)*rec[j] for j in xrange(k-2))
#
#@cached_function
#def number_of_irreducible_permutations(k):
#    r"""
#    Return the number of irreducible permutations on k letters.
#
#    Sloane A003319 (see also A158882, A167894)
#
#    EXAMPLES::
#
#        sage: for i in xrange(6): print i, number_of_irreducible_permutations(i)
#        0 1
#        1 1
#        2 1
#        3 3
#        4 13
#        5 71
#    """
#    if k <= 0: return 0
#    return (factorial(k) - sum(_fact_and_irr(k)))
#
#def rauzy_move((ltop,lbot),l):
#    r"""
#    Start from a generalized permutation given as ((ltop,lbot),l)
#    where ltop and lbot are lists with the top and bottom labels
#    (the labels are assumed to pair up -- no check), and l is a
#    dictionary of lengths for the labels used in ltop, lbot.
#
#    rauzy_move will apply one Rauzy move at the rightmost end:
#
#    - if rightmost top and bottom intervals have distinct lengths,
#      do the usual thing
#
#    - if rightmost top and bottom intervals have same label, error
#
#    - if rightmost top and bottom intervals have same length but
#      not same label, just remove one of the labels
#
#    Note: this code was originally used to compute the connected component of an
#    origami. This is no more ncessary but it might be a good idea to keep for
#    double check
#    """
#    if len(ltop) < 2 and len(lbot) < 2:
#        raise ValueError("cannot induce, too few intervals")
#    ktop = ltop[-1]
#    kbot = lbot[-1]
#    lltop = ltop[:-1]
#    llbot = lbot[:-1]
#    ll = l.copy()
#    if l[ktop] - l[kbot] > 0:
#        # print "bottom interval shorter\n"
#
#        # ltop[-1] is some symbol ktop, which appears twice: once as ltop[-1], and the other
#        # time either as ltop[j] for j < len(ltop) - 1 or as lbot[j] for j < len(lbot) - 1;
#        # as ltop[-1], ktop needs to be replaced by [ltop[-1],lbot[-1]],
#        # the other occurrence of ktop needs to be replaced
#        # - by [ltop[-1],lbot[-1]] if it is lbot[j]
#        # - by [lbot[-1],ltop[-1]] if it is ltop[j]
#        # the dictionary of lengths needs one update: ll[ktop] = l[ktop] - l[kbot]
#
#        if ktop in lltop: lltop.insert(lltop.index(ktop),kbot)
#        else: llbot.insert(llbot.index(ktop)+1,kbot)
#        lltop.append(ktop)
#        ll[ktop] = l[ktop] - l[kbot]
#
#    elif l[ktop] - l[kbot] < 0:
#        # print "top interval shorter\n"
#
#        # lbot[-1] is some symbol kbot, which appears twice: once as lbot[-1], and the other
#        # time either as lbot[j] for j < len(lbot) - 1 or as ltop[j] for j < len(ltop) - 1;
#        # as lbot[-1], kbot needs to be replaced by [lbot[-1],ltop[-1]],
#        # the other occurrence of ktop needs to be replaced
#        # - by [lbot[-1],ltop[-1]] if it is ltop[j]
#        # - by [ltop[-1],lbot[-1]] if it is lbot[j]
#        # the dictionary of lenghts needs one update: ll[kbot] = l[kbot] - l[ktop]
#
#        if kbot in llbot: llbot.insert(llbot.index(kbot),ktop)
#        else: lltop.insert(lltop.index(kbot)+1,ktop)
#        llbot.append(kbot)
#        ll[kbot] = l[kbot] - l[ktop]
#
#    elif ktop == kbot:
#        # in this case the surface has a torus component... weird.
#        raise ValueError("it seems the surface \n %s \n %s \n %s either is disconnected or has genus 1" %s(ltop,lbot,l))
#    else:
#        # ktop != kbot but the intervals labeled ktop and kbot have equal lengths
#        # just remove ktop and kbot at the end of ltop and lbot, and match the
#        # remaining ktop and kbot (keep only one of the names; also in the length dict)
#        if ktop in lltop: lltop[lltop.index(ktop)] = kbot
#        else: llbot[llbot.index(ktop)] = kbot
#        ll.pop(ktop)
#
#    return ((lltop,llbot),ll)
#
#def iet_of_origami(o):
#    r"""
#    x, y are permutations which together represent an origami
#    this origami is assumed to be a torus cover
#    (ie corresponds to an abelian differential)
#
#    First get a multi-iet by considering the first return to I, the union of antidiagonals
#    of all squares, of the flow in direction (1,u), where u is the golden mean.
#    The antidiagonals are parametrized by their horizontal coordinate, so that they appear
#    to have length 1, therefore I has length n. Apply Rauzy moves until I has length 1.
#
#    We could define an iet by using this permutation and this length vector
#    p = Permutation([(lambda k: 2*y((k+1)//2) if k%2 else 2*x(k//2)-1)(j) for j in range(1,n+1)])
#    l = [2-u,u-1] * n
#    """
#    from surface_dynamics.interval_exchanges.all import iet
#    from sage.rings.number_field.all import NumberField
#    from sage.rings.polynomial.all import PolynomialRing
#    from sage.rings.all import RR
#
#    sigma = o.r()
#    tau = o.u()
#    n = o.nb_squares()
#
#    R = PolynomialRing(QQ,'x')
#    x = R.gen()
#    K = NumberField(x**2 - x - 1, 'a',embedding=RR(1.618))
#    phi = K.gen()
#
#    ltop = []
#    for k in xrange(1,n+1):
#        ltop.extend([(k,0),(k,1)])
#    lbot = []
#    for k in xrange(1,n+1):
#        lbot.extend([(sigma.inverse()(k),1),(tau.inverse()(k),0)])
#    l = {}
#    for k in xrange(1,n+1):
#        l[k,0] = 2 - RR(phi)
#        l[k,1] = RR(phi) - 1
#
#    while sum(l.values())>1:
#        ((ltop, lbot), l) = rauzy_move((ltop,lbot),l)
#
#    return iet.Permutation(ltop,lbot),l
#
#    # chain complex of homology
#
#    @cached_method
#    def chain_complex(self,ring=None):
#        r"""
#        Return the chain complex associated to self
#        """
#        return SurfaceChainComplex(self,ring)
#
#

##############################
# Chain complex and homology #
##############################

#class SurfaceChainComplex(SageObject):
#    r"""
#    Surface chain complex associated to a Ribbon graph
#
#    -> intersection of cycles
#
#    EXAMPLES::
#
#        sage: g = RibbonGraph(edges='(0,1)(2,3)(4,5)(6,7)',faces='(0,2,4,6,1,3,5,7)')
#        sage: C = g.chain_complex()
#        sage: C.chain_space(0).rank() == g.num_vertices()
#        True
#        sage: C.chain_space(1).rank() == g.num_edges()
#        True
#        sage: C.chain_space(2).rank() == g.num_faces()
#        True
#
#        sage: C.cycle_space(0).rank() - C.boundary_space(0).rank() == 0
#        True
#        sage: C.cycle_space(1).rank() - C.boundary_space(1).rank() == 2 * g.genus()
#        True
#        sage: C.cycle_space(2).rank() - C.boundary_space(2).rank() == 1
#        True
#    """
#    def __init__(self, g, ring=None):
#        self._ribbon_graph = g
#        if ring is None:
#            ring = ZZ
#        self._base_ring = ring
#        self._chain_spaces = [
#            FreeModule(ring,g.num_vertices()),
#            FreeModule(ring,g.num_edges()),
#            FreeModule(ring,g.num_faces())]
#
#        self._derivatives = []
#        m = matrix(ring, [1]*g.num_vertices())
#        self._derivatives.append(m)
#
#        m = matrix(ring, g.num_vertices(), g.num_edges())
#        for i,e in enumerate(g.edges()):
#            m[g.dart_to_vertex(e[0]),i] -= 1
#            m[g.dart_to_vertex(e[1]),i] += 1
#        self._derivatives.append(m)
#
#        m = matrix(ring, g.num_edges(), g.num_faces())
#        for i,f in enumerate(g.faces()):
#            for e in f:
#                j,o = g.dart_to_edge(e,orientation=True)
#                m[j,i] += o
#        self._derivatives.append(m)
#
#    def parent(self):
#        r"""
#        Returns the ribbon graph from which this chain complex is defined.
#        """
#        return self._ribbon_graph
#
#    def _repr_(self):
#        return "Chain complex over %s" %str(self._base_ring)
#
#    def differential(self, degree=None):
#        r"""
#        Returns a differential of given degree.
#        """
#        if degree is None:
#            return self._derivatives
#        return self._derivatives[degree]
#
#    def chain_space(self, degree=None):
#        if degree is None:
#            return self._chain_spaces
#        return self._chain_spaces[degree]
#
#    @cached_method
#    def cycle_space(self,degree=None):
#        r"""
#        The space of cycles.
#        """
#        if degree is None:
#            return [self.cycle_space(i) for i in xrange(3)]
#        return self._derivatives[degree].right_kernel()
#
#    @cached_method
#    def boundary_space(self,degree=None):
#        r"""
#        Return the boundary space of the given ``degree``.
#        """
#        if degree is None:
#            return [self.boundary_space(i) for i in xrange(3)]
#        if degree == 2:
#            return self.chain_space(2).submodule([])
#        else:
#            return self._derivatives[degree+1].column_space()
#

