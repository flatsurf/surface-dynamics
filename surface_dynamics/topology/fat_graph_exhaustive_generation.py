r"""
Exhaustive generation of fat graphs.

This is done following the McKay canonical augmentation. This module
is experimental.
"""
# ****************************************************************************
#       Copyright (C) 2019 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import absolute_import, print_function
from six.moves import range

import numbers
from collections import defaultdict

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from .fat_graph import FatGraph, list_extrems

###########################
# Miscellaneous functions #
###########################

def minmax(l):
    m = M = l[0]
    for i in range(1, len(l)):
        if l[i] < m:
            m = l[i]
        elif l[i] > M:
            M = l[i]
    return m, M


def num_and_weighted_num(it):
    from sage.rings.integer_ring import ZZ
    from sage.rings.rational_field import QQ
    s = QQ.zero()
    n = ZZ.zero()
    for _, aut in it:
        n += ZZ.one()
        if aut is None:
            s += QQ.one()
        else:
            s += QQ((1, aut.group_cardinality()))
    return n, s


##########################
# Augmentation functions #
##########################

# augment1: trisection
def augment1(cm, aut_grp, g, callback):
    r"""
    Given a unicellular map ``cm`` with a single vertex and automorphism group
    ``aut_grp``, iterate through all its canonical extensions that are
    uniface-univertex maps of greater genus.

    This operation inserts two edges.

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

    if cm._n == 0:
        cm._set_genus1_square()
        aaut_grp = cm.automorphism_group()
        callback('augment1', True, cm, aaut_grp, g - 1)
        if g > 1:
            augment1(cm, aaut_grp, g - 1, callback)
        cm.remove_face_trisection(n)
        return

    for i in R:
        j = i
        for sj in range(fd[fl[i]]):
            k = j
            for sk in range(fd[fl[i]] - sj + (i != j)):
                cm.trisect_face(i, j, k)
                test, aaut_grp = cm._is_canonical(n)
                callback('augment1', test, cm, aaut_grp, g - 1)
                if test and g > 1:
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

    if cm._n == 0:
        # trivial map -> loop (1 vertex, 2 faces)
        cm._set_genus0_loop()
        aaut_grp = cm.automorphism_group()
        callback('augment2', True, cm, aaut_grp, depth - 1)
        if depth > 1:
            augment2(cm, aaut_grp, depth - 1, callback)
        cm.remove_edge(0)
        return

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

    # if callback is not None:
    #     parent = cm.to_string()

    for i in R:
        j = i
        if fdmax1 <= 1:
            if fd[fl[i]] < fdmax0 - 1:
                continue
            niter = fd[fl[i]]
        elif fd[fl[i]] != fdmax0 or fd[fl[i]] < 2 * fdmax1 - 2:
            continue
        else:
            for _ in range(fdmax1 - 1):
                j = fp[j]
            niter = fd[fl[i]] - 2 * fdmax1 + 3

        for _ in range(niter):
            cm.split_face(i, j)
            test, aaut_grp = cm._is_canonical(n)
            callback('augment2', test, cm, aaut_grp, depth - 1)
            if test and depth > 1:
                augment2(cm, aaut_grp, depth - 1, callback)
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

    if cm._n == 0:
        # trivial map -> edge (2 vertices, 1 face)
        cm._set_genus0_edge()
        aaut_grp = cm.automorphism_group()
        callback('augment3', True, cm, aaut_grp, depth - 1)
        if depth > 1:
            augment3(cm, aaut_grp, depth - 1, min_degree, callback)
        cm.contract_edge(0)
        return

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
        # moreover, by the chosen canonical labelling, the degrees of
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
            niter = vd[vl[i]] - 2 * min_degree_loc + 3

        for _ in range(niter):
            cm.split_vertex(i, j)
            assert vd[vl[i]] >= min_degree_loc
            assert vd[vl[j]] >= min_degree_loc
            test, aaut_grp = cm._is_canonical(n)
            callback('augment3', test, cm, aaut_grp, depth - 1)
            if test and depth > 1:
                augment3(cm, aaut_grp, depth - 1, min_degree, callback)
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
class CountAndWeightedCount:
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
class ListCallback:
    def __init__(self):
        self._list = []

    def __call__(self, cm, aut):
        self._list.append(cm.copy())

    def list(self):
        return self._list


# TODO: make it work again! This is the most precious piece of information
# to enhance the exhaustive generation...
# Callback for getting a full trace of the execution
graphviz_header = """/****************************************************************/
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

class FatGraphsTrace:
    """
    A class to trace the execution of the fat graphs generation.

    It is mostly used for debugging/profiling/illustration purposes.
    """
    def __init__(self, filename=None, verbosity=0):
        self._verbosity = int(verbosity)
        self._properties = {}
        self._k = 0             # current number of vertices
        self._vdepth = {}       # vertex -> depth
        self._vnum = {}         # vertex -> apparition in the iteration
        self._edges = {}        # parent -> child
        self._bad_explore = {}  # vertex -> number of dead end
        self._vaut = {}         # number of automorphisms

    def __repr__(self):
        return 'FatGraphs trace for {%s}' % (', '.join('%s=%s' % (k, v)
                                                       for k, v in sorted(self._properties.items())))

    def summary(self, filename=None):
        if filename is None:
            from sys import stdout
            f = stdout
        else:
            f = open(filename, 'w')

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
        mean_depth = float(sum(k * v for k, v in count_by_depth.items())) / sum(v for v in count_by_depth.values())

        childless_by_depth = defaultdict(int)
        for v, child in self._edges.items():
            if self._vdepth[v] != max_depth and not child:
                childless_by_depth[self._vdepth[v]] += 1

        bad_explore_by_depth = defaultdict(int)
        for v, num in self._bad_explore.items():
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

    def graphviz_tree(self, filename=None):
        if filename is None:
            from sys import stdout
            output = stdout
        else:
            output = open(filename, 'w')

        if filename is not None:
            output.close()

        f = open(filename, 'w')
        f.write(graphviz_header.format(vp=vp, ep=ep, fp=fp, g=g, nf=nf, nv=nv))

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
                f.write("""   %s -> %s [color="%s"];\n""" % (s0, s1, col1))
            for cm2, a2 in augment2(cm1, a1, nnf, intermediate):
                s2 = cm2.to_string()
                if s2 != s1:
                    f.write("""   %s [label="%s"];\n""" % (s2, 0))
                    f.write("""   %s -> %s [color="%s"];\n""" % (s1, s2, col2))
                for cm3, a3 in augment3(cm2, a2, nnv, vertex_min_degree, intermediate):
                    s3 = cm3.to_string()
                    if s3 != s2:
                        f.write("""   %s [label="%s"];\n""" % (s3, 0))
                        f.write("""   %s -> %s [color="%s"];\n""" % (s2, s3, col3))
                    yield cm3, a3

#################
# Main iterator #
#################

class StackCallback:
    def __init__(self, gmin, gmax, fmin, fmax, emin, emax, vmin, vmax, vertex_min_degree, callback, filter):
        self._gmin = gmin
        self._gmax = gmax
        self._fmin = fmin
        self._fmax = fmax
        self._emin = emin
        self._emax = emax
        self._vmin = vmin
        self._vmax = vmax
        self._callback = callback
        self._vertex_min_degree = vertex_min_degree
        self._filter = filter

    def __call__(self, caller, test, cm, aut, depth):
        # argument test: answers whether this is a canonical map
        nv = cm._nv
        ne = cm._n // 2
        nf = cm._nf
        g = (-nv + ne - nf)/2 + 1
        assert nv < self._vmax and ne < self._emax and nf < self._fmax, (cm, depth)

        if test:
            if g >= self._gmin and \
               nv >= self._vmin and \
               ne >= self._emin and \
               nf >= self._fmin and \
               (self._vertex_min_degree <= 1 or all(d >= self._vertex_min_degree for d in cm.vertex_profile())) and \
               (self._filter is None or self._filter(cm, aut)):
                self._callback(cm, aut)

            if caller == 'augment1':
                # augment1 creates fat graphs with a single vertex and a single face.
                # Here: nv=1, nf=1, 2g=e.
                assert nv == nf == 1
                if ne >= 2 * self._gmin:
                    # more faces?
                    nfdepth = min(self._fmax - 2, self._emax - ne - 1)
                    if nfdepth:
                        augment2(cm, aut, nfdepth, self)
                    # more vertices?
                    if nf >= self._fmin:
                        nvdepth = min(self._vmax - 2, self._emax - ne - 1)
                        if nvdepth:
                            augment3(cm, aut, nvdepth, max(1, self._vertex_min_degree), self)

            elif caller == 'augment2':
                # augment2 performs face splitting
                # Here: nv=1
                if nf >= self._fmin:
                    # more vertices?
                    nvdepth = min(self._vmax - 2, self._emax - ne - 1)
                    if nvdepth:
                        augment3(cm, aut, nvdepth, max(1, self._vertex_min_degree), self)

            elif caller == 'augment3':
                pass

            else:
                raise RuntimeError('unknown caller')

    def run(self):
        # trivial map (g = 0, nv = 1, nf = 1)
        cm = FatGraph('()', '()', '()')
        nv = 1
        ne = 0
        nf = 1
        g = 0

        if self._vmin <= nv < self._vmax and \
           self._emin <= ne < self._emax and \
           self._fmin <= nf < self._fmax and \
           self._gmin <= g < self._gmax and \
           self._vertex_min_degree == 0 and \
           (self._filter is None or self._filter(cm, aut)):
               self._callback(cm, None)

        cm._realloc(2*self._emax - 2)
        if self._gmax > 1:
            augment1(cm, None, self._gmax - 1, self)
        if self._gmin == 0 and self._fmax > 2:
            assert cm._n == 0, cm
            augment2(cm, None, self._fmax - 2, self)
        if self._gmin == 0 and self._fmin == 1 and self._vmax > 2:
            assert cm._n == 0, cm
            augment3(cm, None, self._vmax - 2, max(1, self._vertex_min_degree), self)

##############
# Main class #
##############

class FatGraphs:
    r"""
    Isomorphism classes of fat graphs with topological constraints.

    EXAMPLES::

        sage: from surface_dynamics import FatGraphs

    Trees and their dual (maps with single vertex) in genus zero are counted by
    Catalan numbers::

        sage: for n in range(2, 10):
        ....:     ntrees1 = 2 * (n-1) * FatGraphs(g=0, nf=n, nv=1).weighted_cardinality()
        ....:     ntrees2 = 2 * (n-1) * FatGraphs(g=0, nf=1, nv=n).weighted_cardinality()
        ....:     assert catalan_number(n-1) == ntrees1 == ntrees2, (n, ntrees1, ntrees2)

    Tutte formulas in genus 0 (combinatorial maps counted by number of edges)::

        sage: R.<t> = QQ[]
        sage: F = FatGraphs(g=0, ne_max=8)
        sage: poly = R.zero()
        sage: def update(cm, aut):
        ....:     global poly
        ....:     p = [cm.face_profile()]
        ....:     aut_card = 1 if aut is None else aut.group_cardinality()
        ....:     poly += 2*cm.num_edges() // aut_card * t**cm.num_edges()
        sage: F.map_reduce(update)
        sage: poly
        208494*t^7 + 24057*t^6 + 2916*t^5 + 378*t^4 + 54*t^3 + 9*t^2 + 2*t
        sage: sum(2 * 3**n / (n+2) / (n+1) * binomial(2*n,n) *t**n for n in range(1,8))
        208494*t^7 + 24057*t^6 + 2916*t^5 + 378*t^4 + 54*t^3 + 9*t^2 + 2*t

    Genus zero with same number of vertices and faces::

        sage: FatGraphs(g=0, nf=2, nv=2).cardinality_and_weighted_cardinality()
        (2, 5/4)
        sage: FatGraphs(g=0, nf=3, nv=3).cardinality_and_weighted_cardinality()
        (23, 41/2)
        sage: FatGraphs(g=0, nf=4, nv=4).cardinality_and_weighted_cardinality()
        (761, 8885/12)

    Duality checks

        sage: for g,nf,nv in [(0,2,3), (0,2,4), (0,3,4),
        ....:                 (1,1,2), (1,1,3), (1,1,4), (1,1,5), (1,2,3), (1,2,4), (1,3,4),
        ....:                 (2,1,2), (2,1,3)]:
        ....:     n1 = FatGraphs(g=g, nf=nf, nv=nv).cardinality_and_weighted_cardinality()
        ....:     n2 = FatGraphs(g=g, nf=nv, nv=nf).cardinality_and_weighted_cardinality()
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

        sage: FatGraphs(g=3, nf=1, nv=1).cardinality_and_weighted_cardinality()
        (131, 495/4)

    Minimum vertex degree bounds::

        sage: for k in range(2,5):
        ....:    F = FatGraphs(g=1, nf=2, nv=2, vertex_min_degree=1)
        ....:    c1 = F.cardinality_and_weighted_cardinality(lambda cm,_: cm.vertex_degree_min() >= k)
        ....:    G = FatGraphs(g=1, nf=2, nv=2, vertex_min_degree=k)
        ....:    c2 = G.cardinality_and_weighted_cardinality()
        ....:    assert c1 == c2
        ....:    print(c1)
        (14, 87/8)
        (8, 47/8)
        (4, 15/8)

        sage: for k in range(2,6):
        ....:    F = FatGraphs(g=0, nf=5, nv=2, vertex_min_degree=1)
        ....:    c1 = F.cardinality_and_weighted_cardinality(lambda cm,_: cm.vertex_degree_min() >= k)
        ....:    G = FatGraphs(g=0, nf=5, nv=2, vertex_min_degree=k)
        ....:    c2 = G.cardinality_and_weighted_cardinality()
        ....:    assert c1 == c2
        ....:    print(c1)
        (28, 123/5)
        (21, 88/5)
        (13, 48/5)
        (7, 18/5)

    Using ranges for vertex and face numbers::

        sage: F = FatGraphs(g=1, nv_min=2, nv_max=4, nf_min=2, nf_max=4)
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
        ....:         c2 = FatGraphs(g=1, nf=nf, nv=nv).cardinality_and_weighted_cardinality()
        ....:         assert c1 == c2
        ....:         print(nf, nv, c1)
        2 2 (24, 167/8)
        2 3 (180, 172)
        3 2 (180, 172)
        3 3 (2048, 6041/3)

    TESTS::

        sage: from surface_dynamics.topology.fat_graph_exhaustive_generation import FatGraphs
        sage: FatGraphs(g=0, nf=1, nv=1).list()
        [FatGraph('()', '()', '()')]
        sage: FatGraphs(g=1, nf=1, nv=1).list()
        [FatGraph('(0,1,2,3)', '(0,2)(1,3)', '(0,1,2,3)')]
        sage: FatGraphs(g=0, nf=2, nv=1).list()
        [FatGraph('(0,1)', '(0,1)', '(0)(1)')]
        sage: FatGraphs(g=0, nf=1, nv=2).list()
        [FatGraph('(0)(1)', '(0,1)', '(0,1)')]
        sage: FatGraphs(g=0, ne=1).list()
        [FatGraph('(0,1)', '(0,1)', '(0)(1)'), FatGraph('(0)(1)', '(0,1)', '(0,1)')]

        sage: F3 = FatGraphs(g=0, nf_max=4, vertex_min_degree=3)
        sage: for fg in F3.list(): assert all(d >= 3 for d in fg.vertex_profile()), fg
        sage: F3.cardinality_and_weighted_cardinality()
        (3, 7/6)
        sage: F4 = FatGraphs(g=0, nf_max=4, vertex_min_degree=4)
        sage: for fg in F4.list(): assert all(d >= 4 for d in fg.vertex_profile()), fg
        sage: F4.cardinality_and_weighted_cardinality()
        (1, 1/2)

        sage: FatGraphs(g=1, ne=3).list()
        [FatGraph('(0,5,4,1,2,3)', '(0,2)(1,3)(4,5)', '(0,1,2,3,4)(5)'),
         FatGraph('(0,5,1,2,4,3)', '(0,2)(1,3)(4,5)', '(0,1,4)(2,3,5)'),
         FatGraph('(0,5,1,2,3,4)', '(0,2)(1,3)(4,5)', '(0,1,2,4)(3,5)'),
         FatGraph('(0,4,1,2,3)(5)', '(0,2)(1,3)(4,5)', '(0,1,2,3,4,5)'),
         FatGraph('(0,4,2,3)(1,5)', '(0,2)(1,3)(4,5)', '(0,4,1,2,3,5)'),
         FatGraph('(0,4,3)(1,2,5)', '(0,2)(1,3)(4,5)', '(0,1,4,2,3,5)')]
    """
    def __init__(self, g=None, nf=None, ne=None, nv=None, vertex_min_degree=0, g_min=None, g_max=None, nf_min=None, nf_max=None, ne_min=None, ne_max=None, nv_min=None, nv_max=None):
        r"""
       INPUT:

        - ``g``, ``g_min``, ``g_max`` - the genus

        - ``nf``, ``nf_min``, ``nf_max`` - number of faces

        - ``ne``, `ne_min``, ``ne_max`` - number of edges

        - ``nv``, ``nv_min``, ``nv_max`` - number of vertices

        - ``vertex_min_degree`` - minimal degree of vertices (default to ``1``)
        """
        self._gmin, self._gmax = self._get_interval(g, g_min, g_max, 0, 'g')
        self._fmin, self._fmax = self._get_interval(nf, nf_min, nf_max, 1, 'nf')
        self._vmin, self._vmax = self._get_interval(nv, nv_min, nv_max, 1, 'nv')
        self._emin, self._emax = self._get_interval(ne, ne_min, ne_max, 0, 'ne')

        self._vertex_min_degree = ZZ(vertex_min_degree)
        if self._vertex_min_degree < 0:
            raise ValueError('vertex_min_degree must be non-negative')

        self._adjust_bounds()

    def _get_interval(self, v, vmin, vmax, low_bnd, name):
        if v is not None:
            if not isinstance(v, numbers.Integral):
                raise TypeError("%s must be an integral" % name)
            v = ZZ(v)
            if v < low_bnd:
                raise ValueError("%s must be >= %d" % (name, low_bnd))
            return (v, v + 1)
        if vmax is None:
            pass
        elif not isinstance(vmax, numbers.Integral):
            raise ValueError("%s_max must be integral" % name)
        else:
            vmax = ZZ(vmax)

        if vmin is None:
            vmin = low_bnd
        elif not isinstance(vmin, numbers.Integral):
            raise TypeError("%s_min must be integral" % name)
        else:
            vmin = ZZ(vmin)
        if vmin < low_bnd:
            raise ValueError("%s_min must be >= %s" % (name, low_bnd))
        return vmin, vmax

    def _adjust_bounds(self):
        r"""
        TESTS::

            sage: from surface_dynamics import FatGraphs
            sage: FatGraphs(g=0, nv_max=4, vertex_min_degree=3)
            Traceback (most recent call last):
            ...
            ValueError: infinitely many fat graphs
            sage: FatGraphs(g=0, nf_max=4, vertex_min_degree=3)
            FatGraphs(g=0, nf_min=2, nf_max=4, ne_min=1, ne_max=4, nv_min=1, nv_max=3, vertex_min_degree=3)

            sage: F = FatGraphs(g=0, nv=2, ne=1, nf=2)
            sage: F
            EmptySet
            sage: F.cardinality_and_weighted_cardinality()
            (0, 0)

            sage: FatGraphs(ne=0).list()
            [FatGraph('()', '()', '()')]
            sage: FatGraphs(ne=1).list()
            [FatGraph('(0,1)', '(0,1)', '(0)(1)'),
             FatGraph('(0)(1)', '(0,1)', '(0,1)')]
            sage: FatGraphs(ne=2).list()
            [FatGraph('(0,1,2,3)', '(0,2)(1,3)', '(0,1,2,3)'),
             FatGraph('(0,2,1)(3)', '(0,1)(2,3)', '(0,2,3)(1)'),
             FatGraph('(0,2)(1,3)', '(0,1)(2,3)', '(0,3)(1,2)'),
             FatGraph('(0,3,2,1)', '(0,1)(2,3)', '(0,2)(1)(3)'),
             FatGraph('(0,2)(1)(3)', '(0,1)(2,3)', '(0,1,2,3)')]
        """
        # variable order: v, e, f, g
        from sage.geometry.polyhedron.constructor import Polyhedron
        eqns = [(-2, 1, -1, 1, 2)] # -2 + v - e + f + 2g = 0
        ieqs = [(-self._vmin,  1, 0, 0, 0), # -vim + v >= 0
                (-self._emin,  0, 1, 0, 0), # -emin + e >= 0
                (-self._fmin,  0, 0, 1, 0), # -fmin + f >= 0
                (-self._gmin,  0, 0, 0, 1), # -gmin + g >= 0
                (0, -self._vertex_min_degree, 2, 0, 0)]  # -v_min_degree*v + 2e >= 0
        if self._vmax is not None:
            ieqs.append((self._vmax-1, -1, 0, 0, 0)) # v < vmax
        if self._emax is not None:
            ieqs.append((self._emax-1, 0, -1, 0, 0)) # e < emax
        if self._fmax is not None:
            ieqs.append((self._fmax-1, 0, 0, -1, 0)) # f < fmax
        if self._gmax is not None:
            ieqs.append((self._gmax-1, 0, 0, 0, -1)) # g < gmax
        P = Polyhedron(ieqs=ieqs, eqns=eqns, ambient_dim=4, base_ring=QQ)
        if P.is_empty():
            self._vmin = self._vmax = 1
            self._emin = self._emax = 0
            self._fmin = self._fmax = 1
            self._gmin = self._gmax = 0
            self._vertex_min_degree = 0
            return
        if not P.is_compact():
            raise ValueError('infinitely many fat graphs')

        half = QQ((1,2))
        self._vmin, self._vmax = minmax([v[0] for v in P.vertices_list()])
        self._vmin = self._vmin.floor()
        self._vmax = (self._vmax + half).ceil()
        self._emin, self._emax = minmax([v[1] for v in P.vertices_list()])
        self._emin = self._emin.floor()
        self._emax = (self._emax + half).ceil()
        self._fmin, self._fmax = minmax([v[2] for v in P.vertices_list()])
        self._fmin = self._fmin.floor()
        self._fmax = (self._fmax + half).ceil()
        self._gmin, self._gmax = minmax([v[3] for v in P.vertices_list()])
        self._gmin = self._gmin.floor()
        self._gmax = (self._gmax + half).ceil()

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics import FatGraphs
            sage: FatGraphs(g=0, ne=1)
            FatGraphs(g=0, nf_min=1, nf_max=3, ne=1, nv_min=1, nv_max=3)
            sage: FatGraphs(g=0, nf=5, nv=1, vertex_min_degree=3)
            FatGraphs(g=0, nf=5, ne=4, nv=1, vertex_min_degree=3)
            sage: FatGraphs(g=0, nf_min=2, nf_max=5, nv_min=1, nv_max=5, vertex_min_degree=3)
            FatGraphs(g=0, nf_min=2, nf_max=5, ne_min=1, ne_max=7, nv_min=1, nv_max=5, vertex_min_degree=3)
        """
        if self._vmin == self._vmax:
            return 'EmptySet'

        if self._gmax == self._gmin + 1:
            genus = 'g=%d' % self._gmin
        else:
            genus = 'gmin=%d, gmax=%d' % (self._gmin, self._gmax)

        if self._fmax == self._fmin + 1:
            faces = 'nf=%d' % self._fmin
        else:
            faces = 'nf_min=%d, nf_max=%d' % (self._fmin, self._fmax)

        if self._emax == self._emin + 1:
            edges = 'ne=%d' % self._emin
        else:
            edges = 'ne_min=%d, ne_max=%d' % (self._emin, self._emax)

        if self._vmax == self._vmin + 1:
            vertices = 'nv=%d' % self._vmin
        else:
            vertices = 'nv_min=%d, nv_max=%d' % (self._vmin, self._vmax)

        constraints = []
        if self._vertex_min_degree > 1:
            constraints.append("vertex_min_degree=%d" % self._vertex_min_degree)

        return "FatGraphs({})".format(", ".join([genus, faces, edges, vertices] + constraints))

    def map_reduce(self, callback, filter=None):
        r"""
        EXAMPLES::

            sage: from surface_dynamics import FatGraphs
            sage: FatGraphs(g=1, nf=2, nv=2).map_reduce(lambda x,y: print(x))
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
        StackCallback(
                self._gmin, self._gmax,
                self._fmin, self._fmax,
                self._emin, self._emax,
                self._vmin, self._vmax,
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

            sage: from surface_dynamics import FatGraphs
            sage: L21 = FatGraphs(g=0, nf=2, nv=1).list()
            sage: L21[0].num_faces()
            2
            sage: L21[0].num_vertices()
            1

            sage: L12 = FatGraphs(g=0, nf=1, nv=2).list()
            sage: L12[0].num_faces()
            1
            sage: L12[0].num_vertices()
            2
        """
        L = ListCallback()
        self.map_reduce(L)
        return L.list()

###################
# Deprecated code #
###################

def FatGraphs_g_nf_nv(g=None, nf=None, nv=None, vertex_min_degree=1, g_min=None, g_max=None, nf_min=None, nf_max=None, nv_min=None, nv_max=None):
    r"""
    TESTS::

        sage: from surface_dynamics.topology.fat_graph_exhaustive_generation import FatGraphs_g_nf_nv
        sage: F = FatGraphs_g_nf_nv(g=1, nf=1, nv=1)
        doctest:warning
        ...
        DeprecationWarning: FatGraphs_g_nf_nv is deprecated. Use FatGraphs
        sage: F.list()
        [FatGraph('(0,1,2,3)', '(0,2)(1,3)', '(0,1,2,3)')]
    """
    from warnings import warn
    warn('FatGraphs_g_nf_nv is deprecated. Use FatGraphs', DeprecationWarning)
    return FatGraphs(g=g, nf=nf, nv=nv, vertex_min_degree=vertex_min_degree, g_min=g_min, g_max=g_max, nf_min=nf_min, nf_max=nf_max, nv_min=nv_min, nv_max=nv_max)
