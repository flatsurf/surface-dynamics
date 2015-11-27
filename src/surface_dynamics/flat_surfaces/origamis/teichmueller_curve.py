r"""
Teichmueller curves of Origamis.
"""

from sage.structure.sage_object import SageObject
from sage.structure.element import Element
from sage.structure.parent import Parent

from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from origami import Origami

from sage.modular.arithgroup.congroup_sl2z import SL2Z
from sage.modular.arithgroup.arithgroup_perm import Lm,Rm,S2m,S3m,EvenArithmeticSubgroup_Permutation, OddArithmeticSubgroup_Permutation
from copy import copy

class TeichmuellerCurve(SageObject):
    pass

class Cusp(SageObject):
    r"""
    A cusp in a Teichmueller curve.

    - width of the cusp

    - cylinder decomposition

    - lengths

    - heights

    - a representative
    """
    def __init__(self, parent, origami, slope):
        self._parent = parent
        self._origami = origami
        self._slope = slope

        self._width = 1
        L,_ = origami.gl2z_edges()
        oo = L[origami]
        while oo != origami:
            oo = L[oo]
            self._width += 1

        c,l,h = origami.cylinder_diagram(data=True)
        self._cyl_diag = c
        self._lengths = l
        self._heights = h

    def _repr_(self):
        r"""
        String representation
        """
        return "Cusp of %s" %str(parent)

    def width(self):
        r"""
        Return the width of the cusp
        """
        return self._width

    def cylinder_diagram(self,data=False):
        r"""
        Return the cylinder diagram associated to this cusp.
        """
        if data:
            return self._cyl_diag, self._lengths, self._heights
        return self._cyl_diags

    def slope(self):
        r"""
        Return one slope
        """
        return self._slope

def TeichmuellerCurveOfOrigami(origami):
    r"""
    Return the teichmueller curve for an origami

    TESTS::

        sage: from surface_dynamics.all import *

        sage: o = Origami('(1,2)','(1,3)')
        sage: o.teichmueller_curve() #indirect test
        Teichmueller curve of the origami
        (1)(2,3)
        (1,2)(3)
    """
    origami = origami.to_standard_form()

    dL,dR,dS = origami.sl2z_edges()
    n = len(dL)

    # find numerotation for the elements of the SL2Z-orbit
    mapping = set(dL)
    mapping.remove(origami)
    mapping = sorted(mapping)
    mapping.insert(0,origami)
    inv_mapping = dict((mapping[i],i) for i in xrange(len(mapping)))

    # consider dL and dR as lists and
    # compute the action of s2 and s3
    # from the one of l and r
    l_edges = [None]*n
    r_edges = [None]*n
    s2_edges = [None]*n
    s3_edges = [None]*n
    for o in mapping:
        i = inv_mapping[o]
        l_edges[i] = inv_mapping[dL[o]]
        r_edges[i] = inv_mapping[dR[o]]
        s2_edges[i] = inv_mapping[dS[o]]
        s3_edges[r_edges[i]] = l_edges[i]

    return TeichmuellerCurveOfOrigami_class(
                mapping, inv_mapping, l_edges, r_edges, s2_edges, s3_edges)

def TeichmuellerCurvesOfOrigamis(origamis, assume_normal_form=False, limit=0, verbose_level=0):
    r"""
    Return a set of Teichmueller curve from a set of origamis

    INPUT:

    - ``origamis`` - iterable of origamis

    - ``assume_normal_form`` - whether or not assume that each origami is in
      normal form

    - ``limit`` - integer (default: 0) - if positive, then stop the computation
      if the size of an orbit is bigger than ``limit``.
    """
    from surface_dynamics.flat_surfaces.origamis.origami_dense import sl2z_orbits

    tcurves = []

    if not assume_normal_form:
        origamis = set(o.relabel() for o in origamis)
    try:
        n = iter(origamis).next().nb_squares()
    except StopIteration:
        return []

    for dL,dR,dS in sl2z_orbits(origamis,n,limit):
        n = len(dL)

        # find numerotation for the elements of the SL2Z-orbit
        mapping = set(dL)
        mapping = sorted(mapping)
        inv_mapping = dict((mapping[i],i) for i in xrange(len(mapping)))

        # consider dL and dR as lists and
        # compute the action of s2 and s3
        # from the one of l and r
        l_edges = [None]*n
        r_edges = [None]*n
        s2_edges = [None]*n
        s3_edges = [None]*n
        for o in mapping:
            i = inv_mapping[o]
            l_edges[i] = inv_mapping[dL[o]]
            r_edges[i] = inv_mapping[dR[o]]
            s2_edges[i] = inv_mapping[dS[o]]
            s3_edges[r_edges[i]] = l_edges[i]

        tcurves.append(TeichmuellerCurveOfOrigami_class(
                    mapping, inv_mapping, l_edges, r_edges, s2_edges, s3_edges))

    return tcurves

class TeichmuellerCurveOfOrigami_class(TeichmuellerCurve):
    def __init__(self, mapping, inv_mapping, l_edges, r_edges, s2_edges, s3_edges):
        N = len(s2_edges)
        ss2 = [None] * len(s2_edges)
        ss3 = [None] * len(s3_edges)
        ll = [None] * len(l_edges)
        rr = [None] * len(r_edges)
        for i in xrange(N):
            ss2[s2_edges[i]] = i
            ss3[s3_edges[i]] = i
            ll[l_edges[i]] = i
            rr[r_edges[i]] = i

        if mapping[0].is_orientation_cover():
            self._veech_group = EvenArithmeticSubgroup_Permutation(ss2,ss3,ll,rr)
        else:
            self._veech_group = OddArithmeticSubgroup_Permutation(ss2,ss3,ll,rr)
        self._l_edges = l_edges
        self._li_edges = ll
        self._r_edges = r_edges
        self._ri_edges = rr
        self._s2_edges = s2_edges
        self._s2i_edges = ss2
        self._s3_edges = s3_edges
        self._s3i_edges = ss3

        self._mapping = mapping
        self._inv_mapping = inv_mapping

    def __reduce__(self):
        return (TeichmuellerCurveOfOrigami_class,
                (self._mapping, self._inv_mapping,
                    self._l_edges, self._r_edges,
                    self._s2_edges, self._s3_edges))

    def _repr_(self):
        r"""
        string representation of self
        """
        return "Teichmueller curve of the origami\n%s" %self.origami()

    def stratum(self):
        r"""
        Returns the stratum of the Teichmueller curve

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: o = Origami('(1,2)','(1,3)')
            sage: t = o.teichmueller_curve()
            sage: t.stratum()
            H_2(2)
        """
        return self[0].stratum()

    def origami(self):
        r"""
        Returns an origmami in this teichmueller curve.

        Beware that the labels of the initial origami may have changed.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: o = Origami('(1,2)','(1,3)')
            sage: t = o.teichmueller_curve()
            sage: t.origami()
            (1)(2,3)
            (1,2)(3)
        """
        return self._mapping[0]

    an_element = origami

    def orbit_graph(self,
            s2_edges=True, s3_edges=True, l_edges=False, r_edges=False,
            vertex_labels=True):
        r"""
        Return the graph of action of PSL on this origami

        INPUT:

        - ``return_map`` - return the list of origamis in the orbit
        """
        G = self._veech_group.coset_graph(
                    s2_edges=s2_edges,
                    s3_edges=s3_edges,
                    l_edges=l_edges,
                    r_edges=r_edges)
        if vertex_labels:
            G.relabel(dict(enumerate(self._mapping)))
        return G

    def __getitem__(self,i):
        assert(isinstance(i,(int,Integer)))
        if i < 0 or i >= len(self._mapping):
            raise IndexError
        return self._mapping[i]

    def veech_group(self):
        r"""
        Return the veech group of the Teichmueller curve
        """
        return self._veech_group

    def plot_graph(self):
        r"""
        Plot the graph of the veech group.

        The graph corresponds to the action of the generators on the cosets
        determined by the Veech group in PSL(2,ZZ).
        """
        return self.orbit_graph.plot()

#    def plot_fundamental_domain(self):
#        from sage.geometry.hyperbolic_geometry import HH
#        from math import cos,sin,pi
#        from sage.matrix.constructor import matrix
#        from sage.all import CC
#        S=matrix([[0,-1],[1,0]])
#        T=matrix([[1,1],[0,1]])
#        g3=S*T
#        tree,gen,lifts=self._graph.spanning_tree(gq=g3)
#        pol = HH.polygon([Infinity,CC(0,1),CC(0,0),CC(-cos(pi/3),sin(pi/3))])
#        res = pol.plot(face_color='blue')
#        for m in lifts[1:]:
#            res += (m*pol).plot(face_color='green',face_alpha=0.3)
#        return res

    def cusp_representative_iterator(self):
        r"""
        Iterator over the cusp of self.

        Each term is a couple ``(o,w)`` where ``o`` is a representative of the
        cusp (an origami) and ``w`` is the width of the cusp (an integer).
        """
        cusps = []
        n = len(self._mapping)
        seen = [True]*n
        l = self._l_edges
        for i in xrange(n):
            if seen[i]:
                k=1
                seen[i] = False
                j = l[i]
                while j != i:
                    seen[j] = False
                    j = l[j]
                    k += 1
                yield (self._mapping[i],k)

    def cusp_representatives(self):
        return list(self.cusp_representative_iterator())

    def sum_of_lyapunov_exponents(self):
        r"""
        Returns the sum of Lyapunov exponents for this origami

        EXAMPLES::

            sage: from surface_dynamics.all import *

        Let us consider few examples in H(2) for which the sum is independant of
        the origami::

            sage: o = Origami('(1,2)','(1,3)')
            sage: o.stratum()
            H_2(2)
            sage: o.sum_of_lyapunov_exponents()
            4/3
            sage: o = Origami('(1,2,3)(4,5,6)','(1,5,7)(2,6)(3,4)')
            sage: o.stratum()
            H_2(2)
            sage: o.sum_of_lyapunov_exponents()
            4/3

        This is true as well for the stratum H(1,1)::

            sage: o = Origami('(1,2)','(1,3)(2,4)')
            sage: o.stratum()
            H_2(1^2)
            sage: o.sum_of_lyapunov_exponents()
            3/2
            sage: o = Origami('(1,2,3,4,5,6,7)','(2,6)(3,7)')
            sage: o.stratum()
            H_2(1^2)
            sage: o.sum_of_lyapunov_exponents()
            3/2

        ALGORITHM:

        Kontsevich-Zorich formula
        """
        K = Integer(1)/Integer(12) * sum(m*(m+2)/(m+1) for m in self.stratum().zeros())
        KK = 0
        for o in self._mapping:
            KK += sum(w/h for (h,w) in o.widths_and_heights())
        return K + Integer(1)/Integer(len(self._mapping)) * KK


    #TODO: mysterious signs and interversion problems...
    def simplicial_action_generators(self):
        r"""
        Return action of generators of the Veech group on homolgy
        """
        from sage.matrix.constructor import matrix, identity_matrix

        tree,reps,word_reps,gens = self._veech_group._spanning_tree_verrill(on_right=False)
        n = self[0].nb_squares()

        l_mat = identity_matrix(2*n)
        s_mat = matrix(2*n)
        s_mat[:n,n:] = -identity_matrix(n)

        reps = [None]*len(tree)
        reps_o = [None]*len(tree)

        reps[0] = identity_matrix(2*n)
        reps_o[0] = self.origami()

        waiting = [0]

        while waiting:
            x = waiting.pop()
            for (e0,e1,label) in tree.outgoing_edges(x):
                waiting.append(e1)
                o = reps_o[e0]
                if label == 's':
                    m = copy(s_mat)
                    m[n:,:n] = o.y().matrix().transpose()
                    reps[e1] = m * reps[e0]
                    reps_o[e1] = o.S_action()
                elif label == 'l':
                    o = o.L_action()
                    m = copy(l_mat)
                    m[:n,n:] = (~o.y()).matrix().transpose()
                    reps[e1] = m * reps[e0]
                    reps_o[e1] = o
                elif label == 'r':
                    o = o.R_action()
                    m = copy(l_mat)
                    m[n:,:n] = (~o.x()).matrix().transpose()
                    reps[e1] = m * reps[e0]
                    reps_o[e1] = o
                else:
                    raise ValueError, "does know how to do only with slr"

        m_gens = []

        for (e0,e1,label) in gens:
            o = reps_o[e0]
            if label == 's':
                oo = o.S_action()
                m = copy(s_mat)
                m[n:,:n] = o.y().matrix().transpose()
            elif label == 'l':
                oo = o.L_action()
                m = copy(l_mat)
                m[:n,n:] = (~oo.y()).matrix().transpose()
            elif label == 'r':
                oo = o.R_action()
                m = copy(l_mat)
                m[:n,n:] = (~oo.x()).matrix().transpose()
            else:
                raise ValueError, "does know how to do only with slr"

            ooo = reps_o[e1]
            ot, oo_renum = oo.standard_form(True)
            os, ooo_renum = ooo.standard_form(True)

            assert(ot == os)

            m_ren = (oo_renum * ~ooo_renum).matrix().transpose()
            m_s = matrix(2*n)
            m_s[:n,:n] = m_ren;
            m_s[n:,n:] = m_ren

            m_gens.append(~reps[e1] * m_s * m * reps[e0])

        return tree,reps,reps_o,gens,m_gens,word_reps

