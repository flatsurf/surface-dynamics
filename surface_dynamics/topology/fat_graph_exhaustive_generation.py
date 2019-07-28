r"""
Exhaustive generation of fat graphs.

This is done following the McKay canonical augmentation. This module
is experimental.
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

import numbers

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from sage.stats.basic_stats import mean

from .fat_graph import FatGraph

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


##########################
# Augmentation functions #
##########################

# augment1: trisection
def augment1(cm, aut_grp, g, callback):
    r"""
    Given a unicellular map ``cm`` with a single vertex and automorphism group
    ``aut_grp``, iterate through all its canonical extensions that are
    uniface-univertex maps of greater genus.

    This operation inserts two edes.

    This augmentation function is sufficient to iterate through unicellular map.
    """
    n = cm._n
    fd = cm._fd
    fl = cm._fl
    fp = cm._fp
    i = 0
    if aut_grp is None:
        R = range(n)
    else:
        aut_grp.reset_iterator()
        R = aut_grp

    for i in R:
        j = i
        for sj in range(fd[fl[i]]):
            k = j
            for sk in range(fd[fl[i]] - sj + (i != j)):
                cm.trisect_face(i, j, k)
                test, aaut_grp = cm._is_canonical(n)
                callback('augment1', test, cm, aaut_grp, g-1)
                if test and g>1:
                    augment1(cm, aaut_grp, g - 1, callback)

                cm.remove_face_trisection(n)
                k = fp[k]
            j = fp[j]
        i = fp[i]

# augment2: face split
# (essentially the same as augment3)
def augment2(cm, aut_grp, depth, callback):
    r"""
    Given a map ``cm`` with a single vertex and automorphism group ``aut_grp``
    iterate through all its canonical extensions that are obtained by splitting
    one of its faces (by adding a single edge).

    Because of the chosen canonical labellings, we only need to consider the
    faces with maximal degree and split in such way that the secondly created
    face is still at least as big as the second biggest.
    """
    n = cm._n
    nf = cm._nf
    fp = cm._fp
    fd = cm._fd
    fl = cm._fl

    if aut_grp is None:
        R = range(n)
    else:
        aut_grp.reset_iterator()
        R = aut_grp

    # TODO find edges with the highest face degrees on their sides
    # (the split must remain higher than them)
    fdmax0 = 0
    fdmax1 = 0
    for i in range(nf):
        if fd[i] > fdmax0:
            fdmax1 = fdmax0
            fdmax0 = fd[i]

    if callback is not None:
        parent = cm.to_string()

    for i in R:
        j = i
        if fdmax1 <= 1:
            if fd[fl[i]] < fdmax0 - 1:
                continue
            niter = fd[fl[i]]
        elif fd[fl[i]] != fdmax0 or vd[fl[i]] < 2 * fdmax1 - 2:
            continue
        else:
            for _ in range(fdmax1 - 1):
                j = fp[j]
            niter = fd[fl[i]] - 2*fdmax1 + 3

        for _ in range(niter):
            cm.split_face(i, j)
            test, aaut_grp = cm._is_canonical(n)
            callback('augment2', test, cm, aaut_grp, depth-1)
            if test and depth > 1:
                augment2(cm, aaut_grp, depth-1, callback)
            cm.remove_edge(n)
            j = fp[j]

# augment3: vertex split
def augment3(cm, aut_grp, depth, min_degree, callback):
    r"""
    Given a map ``cm``, its automorphism group ``aut_grp`` and a minimum
    degree ``min_degree``, iterate through all the canonical extensions of
    ``cm`` that are obtained by splitting ``depth`` times a vertices into
    two vertices.

    This augmentation add ``depth`` edges to the fat graph.

    In principle, because of the chosen canonical labellings, we only need to
    consider the vertices with maximal degree.
    """
    n = cm._n
    nv = cm._nv
    vp = cm._vp
    vd = cm._vd
    vl = cm._vl

    if aut_grp is None:
        R = range(n)
    else:
        aut_grp.reset_iterator()
        R = aut_grp

    # TODO: find edges with the highest vertex degrees on their ends
    # the split must remain higher than them
    vdmax0 = 0
    vdmax1 = 0
    for i in range(nv):
        if vd[i] > vdmax0:
            vdmax1 = vdmax0
            vdmax0 = vd[i]
    min_degree_loc = max(min_degree, vdmax1)

    for i in R:
        # vertex degrees are split as d -> (d1 + 1, d2 + 1)
        # so, if min_degree > 1 we can ignore the first/last half edges
        # moreover, by the choosen canonical labelling, the degrees of
        # the split vertices must remain larger than any other

        j = i
        if min_degree_loc == 1:
            if vd[vl[i]] < vdmax0 - 1:
                continue
            niter = vd[vl[i]]
        elif vd[vl[i]] != vdmax0 or vd[vl[i]] < 2 * min_degree_loc - 2:
            continue
        else:
            for _ in range(min_degree_loc - 1):
                j = vp[j]
            niter = vd[vl[i]] - 2*min_degree_loc + 3

        for _ in range(niter):
            cm.split_vertex(i, j)
            assert vd[vl[i]] >= min_degree_loc
            assert vd[vl[j]] >= min_degree_loc
            test, aaut_grp = cm._is_canonical(n)
            callback('augment3', test, cm, aaut_grp, depth-1)
            if test and depth > 1:
                augment3(cm, aaut_grp, depth-1, min_degree, callback)
            cm.contract_edge(n)
            j = vp[j]

# TODO
# def augment4(cm):
#    r"""
#    Plant vertices of degree 1,2.
#    """

########################################
# Callbacks for the various map reduce #
########################################

# callback to count elements (behaves somehow as a 2-tuple)
class CountAndWeightedCount(object):
    def __init__(self):
        self.count = ZZ(0)
        self.weighted_count = QQ(0)
    def __repr__(self):
        return "(%s, %s)" % (self.count, self.weighted_count)
    def __len__(self):
        return 2
    def __getitem__(self, i):
        if not isinstance(i, numbers.Integral):
            raise TypeError
        i = int(i)
        if i == -1 or i == 1:
            return self.weighted_count
        elif i == 0:
            return self.count
        else:
            raise IndexError("index out of range")
    def __eq__(self, other):
        if type(self) is not type(other):
            raise TypeError
        return self.count == other.count and self.weighted_count == other.weighted_count
    def __ne__(self, other):
        if type(self) is not type(other):
            raise TypeError
        return self.count != other.count or self.weighted_count != other.weighted_count
    def __call__(self, cm, aut):
        self.count += ZZ(1)
        self.weighted_count += QQ((1, (1 if aut is None else aut.group_cardinality())))

# callback to list elements
class ListCallback(object):
    def __init__(self):
        self._list = []

    def __call__(self, cm, aut):
        self._list.append(cm.copy())

    def list(self):
        return self._list

# TODO: make it work again! This is the most precious piece of information
# to enhance the exhaustive generation...
# Callback for getting a full trace of the execution
grapvhiz_header="""/****************************************************************/
/* Trace execution of fat graphs generation                     */
/*                                                              */
/* root: vp={vp:6} ep={ep:6} fp={fp:6}                          */
/* g = {g:2}                                                       */
/* nf = {nf:2}                                                      */
/* nv = {nv:2}                                                      */
                                                                */
 To compile to a graph in pdf format run                        */
                                                                */
      $ sfpdf -Tpdf -o OUTPUT.pdf INPUT.dot                     */
                                                                */
/****************************************************************/
"""

class FatGraphsTrace(object):
    """
    A class to trace the execution of the fat graphs generation.

    It is mostly used for debugging/profiling/illustration purposes.
    """
    def __init__(self, filename=None, verbosity=0):
        self._verbosity = int(verbosity)
        self._properties = {}
        self._k = 0            # current number of vertices
        self._vdepth = {}      # vertex -> depth
        self._vnum   = {}      # vertex -> apparition in the iteration
        self._edges = {}       # parent -> child
        self._bad_explore = {} # vertex -> number of dead end
        self._vaut = {}        # number of automorphisms

    def __repr__(self):
        return 'FatGraphs trace for {%s}' % (', '.join('%s=%s' % (k,v) for k,v in sorted(self._properties.items())))

    def summary(self, filename=None):
        if filename is None:
            from sys import stdout
            f = stdout
        else:
            f = open(fiename, 'w')

        f.write(repr(self))
        f.write('\n')

        if not self._vdepth:
            if filename is not None:
                f.close()
            return

        count_by_depth = defaultdict(int)
        for v in self._vdepth.values():
            count_by_depth[v] += 1
        max_depth = max(count_by_depth)
        mean_depth = float(sum(k*v for k,v in count_by_depth.items())) / sum(v for v in count_by_depth.values())

        childless_by_depth = defaultdict(int)
        for v,child in self._edges.items():
            if self._vdepth[v] != max_depth and not child:
                childless_by_depth[self._vdepth[v]] += 1

        bad_explore_by_depth = defaultdict(int)
        for v,num in self._bad_explore.items():
            bad_explore_by_depth[self._vdepth[v]] += num

        f.write('depth                   : %d\n' % max_depth)
        f.write('mean depth              : %f\n' % mean_depth)
        f.write('total num fat graphs    : %d\n' % len(self._vdepth))
        f.write('by depth num fat graphs : %s\n' % ' '.join('%d' % count_by_depth[i] for i in range(max_depth + 1)))
        f.write('total bad explore       : %d\n' % sum(bad_explore_by_depth.values()))
        f.write('by depth bad explore    : %s\n' % ' '.join('%d' % bad_explore_by_depth[i] for i in range(max_depth + 1)))
        f.write('total childless         : %d\n' % sum(childless_by_depth.values()))
        f.write('by depth childless      : %s\n' % ' '.join('%d' % childless_by_depth[i] for i in range(max_depth + 1)))

        if filename is not None:
            f.close()

    def __call__(self, cm, aut):
        pass
    def add_vertex(self, s, aut_grp, depth):
        if self._verbosity >= 1:
            print('add_vertex(s={}, depth={})'.format(s, depth))
        if s in self._vdepth:
            raise RuntimeError("already explored vertex!")
        self._vdepth[s] = depth
        self._bad_explore[s] = 0
        self._vnum[s] = self._k
        self._edges[s] = []
        self._k += 1
        self._vaut[s] = 1 if aut_grp is None else aut_grp.group_cardinality()

    def add_edge(self, s0, s1):
        if self._verbosity >= 1:
            print('add_edge(s0={}, s1={})'.format(s0, s1))
        if s0 not in self._edges:
            raise RuntimeError("_edges not properly initialized")
        self._edges[s0].append(s1)

    def root(self, cm, aut_grp):
        if self._k:
            raise RuntimeError("trying to set root in a non-empty trace")
        s = cm.to_string()
        if self._verbosity >= 1:
            print('root(cm={})'.format(s))
        self.add_vertex(s, aut_grp, 0)

    def canonical_edge(self, parent, cm, aut_grp, caller):
        if not isinstance(parent, str):
            raise RuntimeError

        s = cm.to_string()
        if self._verbosity >= 1:
            print('canonical_edge(parent={}, cm={}, caller={})'.format(parent, s, caller))

        if parent not in self._edges:
            raise RuntimeError("_edges not properly initialized at %s" % parent)
        if s in self._edges:
            raise RuntimeError("s already in _edges")

        self.add_vertex(s, aut_grp, self._vdepth[parent] + 1)
        self.add_edge(parent, s)

    def non_canonical_edge(self, parent, cm, caller):
        if not isinstance(parent, str):
            raise RuntimeError
        if self._verbosity >= 1:
            print('non_canonical_edge(parent={}, caller={})'.format(parent, caller))
        if parent not in self._edges:
            raise RuntimeError("_edges not properly initialized at %s" % parent)
        self._bad_explore[parent] += 1

    def grapvhiz_tree(self, filename=None):
        if filename is None:
            from sys import stdout
            output = stdout
        else:
            output = open(filename, 'w')

        if filename is not None:
            output.close()

        f = open(filename, 'w')
        f.write(header.format(vp=vp, ep=ep, fp=fp, g=g, nf=nf, nv=nv))

        f.write('digraph Tree {\n')
        f.write('   rankdir = TB;\n')

        col1 = "#FF0000"
        col2 = "#00FF00"
        col3 = "#0000FF"

        a0 = cm0.automorphism_group()
        s0 = cm0.to_string()
        f.write("""   %s [label="%s"];\n""" % (s0, 0))
        for cm1, a1 in augment1(cm0, a0, g, False):
            s1 = cm1.to_string()
            if s0 != s1:
                f.write("""   %s [label="%s"];\n""" % (s1, 0))
                f.write("""   %s -> %s [color="%s"];\n""" %(s0, s1, col1))
            for cm2, a2 in augment2(cm1, a1, nnf, intermediate):
                s2 = cm2.to_string()
                if s2 != s1:
                    f.write("""   %s [label="%s"];\n""" % (s2, 0))
                    f.write("""   %s -> %s [color="%s"];\n""" %(s1, s2, col2))
                for cm3, a3 in augment3(cm2, a2, nnv, vertex_min_degree, intermediate):
                    s3 = cm3.to_string()
                    if s3 != s2:
                        f.write("""   %s [label="%s"];\n""" % (s3, 0))
                        f.write("""   %s -> %s [color="%s"];\n""" %(s2, s3, col3))
                    yield cm3, a3

#################
# Main iterator #
#################

class StackCallback(object):
    def __init__(self, cm, aut, gmin, gdepth, nfmin, nfdepth, nvmin, nvdepth, min_degree, callback, filter):
        self._gmin = gmin
        self._gdepth = gdepth
        self._nfmin = nfmin 
        self._nfdepth = nfdepth
        self._nvmin = nvmin
        self._nvdepth = nvdepth
        self._cm = cm
        self._aut = aut
        self._callback = callback
        self._min_degree = min_degree
        self._filter = filter
        # print("StackCallback(gmin={}, gdepth={}, nfmin={}, nfdepth={}, nvmin={}, nvdepth={}".format(
        #    self._gmin, self._gdepth, self._nfmin, self._nfdepth, self._nvmin, self._nvdepth))

    def __call__(self, caller, test, cm, aut, depth):
        if test:
            if caller == 'augment1':
                # we know the map has one vertex and one edge and hence 4g = n
                if cm._n >= 4 * self._gmin:
                    if cm._nv >= self._nvmin and cm._nf >= self._nfmin and \
                       (self._filter is None or self._filter(cm,aut)):
                        self._callback(cm, aut)

                    # more faces?
                    if self._nfdepth:
                        augment2(cm, aut, self._nfdepth, self)
                    # more vertices?
                    elif self._nvdepth:
                        augment3(cm, aut, self._nvdepth, self._min_degree, self)

            elif caller == 'augment2':
                if cm._nf >= self._nfmin:
                    if cm._nv >= self._nvmin and \
                       (self._filter is None or self._filter(cm, aut)):
                        self._callback(cm, aut)

                    # more vertices?
                    if self._nvdepth:
                        augment3(cm, aut, self._nvdepth, self._min_degree, self)

            elif caller == 'augment3':
                if cm._nv >= self._nvmin and \
                   (self._filter is None or self._filter(cm, aut)):
                    self._callback(cm, aut)

    def run(self):
        if self._gdepth:
            augment1(self._cm, self._aut, self._gdepth, self)
        elif self._nfdepth:
            augment2(self._cm, self._aut, self._nfdepth, self)
        elif self._nvdepth:
            augment3(self._cm, self._aut, self._nvdepth, self._min_degree, self)
        elif self._filter is None or self._filter(cm, aut):
            self._callback(self._cm, self._aut)

##############
# Main class #
##############

class FatGraphs_g_nf_nv(object):
    r"""
    Isomorphism classes of fat graphs with given genus, number of faces and
    number of vertices.

    EXAMPLES::

        sage: from surface_dynamics.topology.fat_graph_exhaustive_generation import FatGraphs_g_nf_nv

    Trees and their dual (maps with single vertex) in genus zero are counted by
    Catalan numbers::

        sage: for n in range(2, 10):
        ....:     ntrees1 = 2 * (n-1) * FatGraphs_g_nf_nv(0, n, 1).weighted_cardinality()
        ....:     ntrees2 = 2 * (n-1) * FatGraphs_g_nf_nv(0, 1, n).weighted_cardinality()
        ....:     assert catalan_number(n-1) == ntrees1 == ntrees2, n

    Genus zero with same number of vertices and faces::
        
        sage: FatGraphs_g_nf_nv(0, 2, 2).cardinality_and_weighted_cardinality()
        (2, 5/4)
        sage: FatGraphs_g_nf_nv(0, 3, 3).cardinality_and_weighted_cardinality()
        (23, 41/2)
        sage: FatGraphs_g_nf_nv(0, 4, 4).cardinality_and_weighted_cardinality()
        (761, 8885/12)

    Duality checks

        sage: for g,nf,nv in [(0,2,3), (0,2,4), (0,3,4),
        ....:                 (1,1,2), (1,1,3), (1,1,4), (1,1,5), (1,2,3), (1,2,4), (1,3,4),
        ....:                 (2,1,2), (2,1,3)]:
        ....:     n1 = FatGraphs_g_nf_nv(g, nf, nv).cardinality_and_weighted_cardinality()
        ....:     n2 = FatGraphs_g_nf_nv(g, nv, nf).cardinality_and_weighted_cardinality()
        ....:     assert n1 == n2
        ....:     print(g, nf, nv, n1)
        0 2 3 (5, 11/3)
        0 2 4 (14, 93/8)
        0 3 4 (108, 103)
        1 1 2 (3, 5/3)
        1 1 3 (11, 35/4)
        1 1 4 (46, 42)
        1 1 5 (204, 385/2)
        1 2 3 (180, 172)
        1 2 4 (1198, 14065/12)
        1 3 4 (18396, 18294)
        2 1 2 (53, 483/10)
        2 1 3 (553, 539)

    Unicellular map with one vertex in genus 3::

        sage: FatGraphs_g_nf_nv(3, 1, 1).cardinality_and_weighted_cardinality()
        (131, 495/4)

    Minimum vertex degree bounds::

        sage: for k in range(2,5):
        ....:    F = FatGraphs_g_nf_nv(1, 2, 2, vertex_min_degree=1)
        ....:    c1 = F.cardinality_and_weighted_cardinality(lambda cm,_: cm.vertex_degree_min() >= k)
        ....:    G = FatGraphs_g_nf_nv(1, 2, 2, vertex_min_degree=k)
        ....:    c2 = G.cardinality_and_weighted_cardinality()
        ....:    assert c1 == c2
        ....:    print(c1)
        (14, 87/8)
        (8, 47/8)
        (4, 15/8)

        sage: for k in range(2,6):
        ....:    F = FatGraphs_g_nf_nv(0, 5, 2, vertex_min_degree=1)
        ....:    c1 = F.cardinality_and_weighted_cardinality(lambda cm,_: cm.vertex_degree_min() >= k)
        ....:    G = FatGraphs_g_nf_nv(0, 5, 2, vertex_min_degree=k)
        ....:    c2 = G.cardinality_and_weighted_cardinality()
        ....:    assert c1 == c2
        ....:    print(c1)
        (28, 123/5)
        (21, 88/5)
        (13, 48/5)
        (7, 18/5)

    Using ranges for vertex and face numbers::

        sage: F = FatGraphs_g_nf_nv(1, nv_min=2, nv_max=4, nf_min=2, nf_max=4)
        sage: def check(cm,aut):
        ....:    if cm.num_vertices() < 2 or \
        ....:       cm.num_vertices() >= 4 or \
        ....:       cm.num_faces() < 2 or \
        ....:       cm.num_faces() >= 4:
        ....:           raise ValueError(str(cm))
        sage: F.map_reduce(check)
        sage: for nf in [2,3]:
        ....:     for nv in [2,3]:
        ....:         c1 = F.cardinality_and_weighted_cardinality(lambda cm,_: cm.num_vertices() == nv and cm.num_faces() == nf)
        ....:         c2 = FatGraphs_g_nf_nv(1, nf, nv).cardinality_and_weighted_cardinality()
        ....:         assert c1 == c2
        ....:         print(nf, nv, c1)
        2 2 (24, 167/8)
        2 3 (180, 172)
        3 2 (180, 172)
        3 3 (2048, 6041/3)
    """
    def __init__(self, g=None, nf=None, nv=None, vertex_min_degree=1, g_min=None, g_max=None, nf_min=None, nf_max=None, nv_min=None, nv_max=None):
        r"""
       INPUT:

        - ``g``, ``g_min``, ``g_max`` - the genus

        - ``nf``, ``nf_min``, ``nf_max`` - number of faces

        - ``nv``, ``nv_min``, ``nv_max`` - number of vertices

        - ``vertex_min_degree`` - minimal number of vertices (default to ``1``)

        - ``intermediate`` - if set to ``True`` then return all graphs with genus = g,
          number of faces <= nf and number of vertices <= nv that satisfy the
          ``vertex_min_degree`` constraint
        """
        self._gmin, self._gmax = self._get_interval(g, g_min, g_max, 0, 'g')
        if self._gmax != self._gmin + 1:
            raise ValueError("not implemented for genus in an interval")
        self._fmin, self._fmax = self._get_interval(nf, nf_min, nf_max, 1, 'nf')
        self._vmin, self._vmax = self._get_interval(nv, nv_min, nv_max, 1, 'nv')
        self._vertex_min_degree = ZZ(vertex_min_degree)

    def _get_interval(self, v, vmin, vmax, low_bnd, name):
        if v is not None:
            if not isinstance(v, numbers.Integral):
                raise TypeError("%s must be an integral" % name)
            v = ZZ(v)
            if v < low_bnd:
                raise ValueError("%s must be >= %d" % (name, low_bnd))
            return (v, v+1)
        if vmax is None:
            raise ValueError("at least %s or %s_max must be set" % (name, name))
        if not isinstance(vmax, numbers.Integral):
            raise ValueError("%s_max must be integral" % name)
        vmax = ZZ(vmax)
        if vmin is None:
            vmin = low_bnd
        elif not isinstance(vmin, numbers.Integral):
            raise TypeError("%s_min must be integral" % name)
        vmin = ZZ(vmin)
        if vmin < low_bnd:
            raise ValueError("%s_min must be >= %s" % (name, low_bnd))
        return vmin, vmax

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph_exhaustive_generation import FatGraphs_g_nf_nv
            sage: FatGraphs_g_nf_nv(0, 5, 1, vertex_min_degree=3)
            Fat graphs of genus 0, with 5 faces and 1 vertices
            sage: FatGraphs_g_nf_nv(0, nf_min=2, nf_max=5, nv_min=1, nv_max=5, vertex_min_degree=3)
            Fat graphs of genus 0, with between 2 and 5 faces and between 1 and 5 vertices
        """
        if self._gmax == self._gmin + 1:
            genus = '%d' % self._gmin
        else:
            genus = 'between %d and %d' % (self._gmin, self._gmax)

        if self._fmax == self._fmin + 1:
            faces = '%d' % self._fmin
        else:
            faces = 'between %d and %d' % (self._fmin, self._fmax)

        if self._vmax == self._vmin + 1:
            vertices = '%d' % self._vmin
        else:
            vertices = 'between %d and %d' % (self._vmin, self._vmax)

        return "Fat graphs of genus %s, with %s faces and %s vertices" % (genus, faces, vertices)

    def map_reduce(self, callback, filter=None):
        r"""
        EXAMPLES::

            sage: from __future__ import print_function
            sage: from surface_dynamics.topology.fat_graph_exhaustive_generation import FatGraphs_g_nf_nv
            sage: FatGraphs_g_nf_nv(1, 2, 2).map_reduce(lambda x,y: print(x))
            FatGraph('(0,6,5,4,1,2,3)(7)', '(0,2)(1,3)(4,5)(6,7)', '(0,1,2,3,4,6,7)(5)')
            FatGraph('(0,6,1,2,3)(4,7,5)', '(0,2)(1,3)(4,5)(6,7)', '(0,1,2,3,6,4,7)(5)')
            FatGraph('(0,5,4,1,6,2,3)(7)', '(0,2)(1,3)(4,5)(6,7)', '(0,6,7,1,2,3,4)(5)')
            FatGraph('(0,7,2,3)(1,6,5,4)', '(0,2)(1,3)(4,5)(6,7)', '(0,7,1,2,3,4,6)(5)')
            FatGraph('(0,5,4,1,2,6,3)(7)', '(0,2)(1,3)(4,5)(6,7)', '(0,1,6,7,2,3,4)(5)')
            FatGraph('(0,7,3)(1,2,6,5,4)', '(0,2)(1,3)(4,5)(6,7)', '(0,1,7,2,3,4,6)(5)')
            FatGraph('(0,5,4,1,2,3,6)(7)', '(0,2)(1,3)(4,5)(6,7)', '(0,1,2,6,7,3,4)(5)')
            FatGraph('(0,7)(1,2,3,6,5,4)', '(0,2)(1,3)(4,5)(6,7)', '(0,1,2,7,3,4,6)(5)')
            FatGraph('(0,5,4,6,1,2,3)(7)', '(0,2)(1,3)(4,5)(6,7)', '(0,1,2,3,6,7,4)(5)')
            FatGraph('(0,5,4,6,2,3)(1,7)', '(0,2)(1,3)(4,5)(6,7)', '(0,6,1,2,3,7,4)(5)')
            FatGraph('(0,5,6,4,1,2,3)(7)', '(0,2)(1,3)(4,5)(6,7)', '(0,1,2,3,4)(5,6,7)')
            FatGraph('(0,5,6,1,2,3)(4,7)', '(0,2)(1,3)(4,5)(6,7)', '(0,1,2,3,6,4)(5,7)')
            FatGraph('(0,5,6,2,3)(1,7,4)', '(0,2)(1,3)(4,5)(6,7)', '(0,6,1,2,3,4)(5,7)')
            FatGraph('(0,5,6,3)(1,2,7,4)', '(0,2)(1,3)(4,5)(6,7)', '(0,1,6,2,3,4)(5,7)')
            FatGraph('(0,6,5,1,2,4,3)(7)', '(0,2)(1,3)(4,5)(6,7)', '(0,1,4,6,7)(2,3,5)')
            FatGraph('(0,6,1,2,4,3)(5,7)', '(0,2)(1,3)(4,5)(6,7)', '(0,1,4,7)(2,3,6,5)')
            FatGraph('(0,6,4,3)(1,2,7,5)', '(0,2)(1,3)(4,5)(6,7)', '(0,1,4,7)(2,3,5,6)')
            FatGraph('(0,6,5,1,2,3,4)(7)', '(0,2)(1,3)(4,5)(6,7)', '(0,1,2,4,6,7)(3,5)')
            FatGraph('(0,5,1,6,2,3,4)(7)', '(0,2)(1,3)(4,5)(6,7)', '(0,6,7,1,2,4)(3,5)')
            FatGraph('(0,5,1,6,3,4)(2,7)', '(0,2)(1,3)(4,5)(6,7)', '(0,7,1,6,2,4)(3,5)')
            FatGraph('(0,5,1,2,3,6,4)(7)', '(0,2)(1,3)(4,5)(6,7)', '(0,1,2,4)(3,5,6,7)')
            FatGraph('(0,5,1,2,3,6)(4,7)', '(0,2)(1,3)(4,5)(6,7)', '(0,1,2,6,4)(3,5,7)')
            FatGraph('(0,7,4)(1,2,3,6,5)', '(0,2)(1,3)(4,5)(6,7)', '(0,1,2,4,6)(3,5,7)')
            FatGraph('(0,5,7,4)(1,2,3,6)', '(0,2)(1,3)(4,5)(6,7)', '(0,1,2,4)(3,6,5,7)')
        """
        g = self._gmin
        fmin = self._fmin
        fmax = self._fmax
        vmin = self._vmin
        vmax = self._vmax
        if g == 0:
            if fmin == 1 and vmin == 1:
                raise NotImplementedError
            elif fmin > 1:
                # start with face splitting of the bicellular map (g = 0, nv = 1, nf = 2)
                cm0 = FatGraph('(0,1)', '(0,1)', '(0)(1)')
                gshift = 0
                vshift = - 1
                fshift = - 2
            elif vmin > 1:
                # start with vertex splitting of the unicellular map (g = 0, nv = 2, nf = 1)
                cm0 = FatGraph('(0)(1)', '(0,1)', '(0,1)')
                gshift = 0
                vshift = - 2
                fshift = - 1
            else:
                raise RuntimeError("this should not happen")
        else:
            # start with trisection of the trivial map (g = 1, nv = 1, nf = 1)
            cm0 = FatGraph.from_unicellular_word([0,1,0,1])
            gshift = - 1
            vshift = - 1
            fshift = - 1
        cm0._realloc(4 * g + 2 * (fmax + vmax - 2))

        a0 = cm0.automorphism_group()
        StackCallback(cm0, a0,
                g,                  # gmin
                g + gshift,         # depth for trisections
                fmin,               # nfmin
                fmax + fshift - 1,  # depth for face splitting
                vmin,               # nvmin
                vmax + vshift - 1,  # depth for vertex splitting
                self._vertex_min_degree,
                callback,
                filter).run()

    def cardinality_and_weighted_cardinality(self, filter=None):
        N = CountAndWeightedCount()
        self.map_reduce(N, filter)
        return tuple(N)

    def weighted_cardinality(self, filter=None):
        N = CountAndWeightedCount()
        self.map_reduce(N, filter)
        return N[1]

    def list(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.topology.fat_graph_exhaustive_generation import FatGraphs_g_nf_nv
            sage: L21 = FatGraphs_g_nf_nv(0, 2, 1).list()
            sage: L21[0].num_faces()
            2
            sage: L21[0].num_vertices()
            1


            sage: L12 = FatGraphs_g_nf_nv(0, 1, 2).list()
            sage: L12[0].num_faces()
            1
            sage: L12[0].num_vertices()
            2
        """
        L = ListCallback()
        self.map_reduce(L)
        return L.list()

