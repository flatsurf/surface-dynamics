#coding: utf-8
r"""
Some utility function for representations of finite groups

Most of the functions are just GAP wrappers.
"""
# *************************************************************************
# Copyright (C) 2015-2016 Charles Fougeron <charlesfougeron@gmail.com>
#                         Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at http://www.gnu.org/licenses/
# *************************************************************************

from sage.libs.gap.libgap import libgap

def real_characters(G):
    r"""
    Return a pair ``(table of characters, character degrees)`` for the
    group ``G``.

    OUPUT:

    - table of characters - a list of characters. Each character is represented
      as the list of its values on conjugacy classes. The order of conjugacy
      classes is the same as the one returned by GAP.

    - degrees - the list of degrees of the characters

    EXAMPLES::

        sage: from surface_dynamics.misc.group_representation import real_characters
        sage: T, deg = real_characters(AlternatingGroup(5))
        sage: T
        [(1, 1, 1, 1, 1),
         (3, -1, 0, -E(5) - E(5)^4, -E(5)^2 - E(5)^3),
         (3, -1, 0, -E(5)^2 - E(5)^3, -E(5) - E(5)^4),
         (4, 0, 1, -1, -1),
         (5, 1, -1, 0, 0)]
        sage: set(parent(x) for chi in T for x in chi)
        {Universal Cyclotomic Field}
        sage: deg
        [1, 3, 3, 4, 5]

        sage: T, deg = real_characters(CyclicPermutationGroup(6))
        sage: T
        [(1, 1, 1, 1, 1, 1),
         (1, -1, 1, -1, 1, -1),
         (2, -1, -1, 2, -1, -1),
         (2, 1, -1, -2, -1, 1)]
        sage: deg
        [1, 1, 1, 1]
    """
    from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField

    G = libgap(G)
    UCF = UniversalCyclotomicField()
    n = G.ConjugacyClasses().Length()
    Tgap = G.CharacterTable().Irr()
    degrees = [chi.Degree().sage() for chi in Tgap]
    Tgap = [tuple(UCF(chi[j]) for j in xrange(n)) for chi in Tgap]

    real_T = []
    real_degrees = []
    seen = set()
    for i,chi in enumerate(Tgap):
        if chi in seen:
            continue

        real_degrees.append(degrees[i])
        if all(x.is_real() for x in chi):
            real_T.append(chi)
            seen.add(chi)
        else:
            seen.add(chi)
            chi_bar = tuple(z.conjugate() for z in chi)
            seen.add(chi_bar)
            real_T.append(tuple(chi[j] + chi[j].conjugate() for j in xrange(n)))

    return (real_T, real_degrees)

def conjugacy_class_matrix(cl, d):
    r"""
    Return the matrix associated to a given conjugacy class of a permutation
    group.

    The result is a `d \times d` matrix that is invariant under conjugation by
    the group action on `\ZZ^d`. It is used to produce the projection matrices
    on the isotypic subspaces.

    EXAMPLES::
        
        sage: from surface_dynamics.misc.group_representation import conjugacy_class_matrix

        sage: G = QuaternionGroup()
        sage: Ggap = libgap(G)
        sage: cls = Ggap.ConjugacyClasses()
        sage: m = conjugacy_class_matrix(cls[2], 8)
        sage: m
        [0 0 1 0 0 0 0 0]
        [0 0 0 1 0 0 0 0]
        [1 0 0 0 0 0 0 0]
        [0 1 0 0 0 0 0 0]
        [0 0 0 0 0 0 1 0]
        [0 0 0 0 0 0 0 1]
        [0 0 0 0 1 0 0 0]
        [0 0 0 0 0 1 0 0]

        sage: for cl in cls:
        ....:      m = conjugacy_class_matrix(cl, 8)
        ....:      for g in G:
        ....:          g = g.matrix()
        ....:          assert g*m*~g == m

        sage: A = AlternatingGroup(5)
        sage: Agap = libgap(A)
        sage: Ggap = Agap.Action(Agap)
        sage: gens = libgap.GeneratorsOfGroup(Ggap).sage()
        sage: G = PermutationGroup(gens)
        sage: cls = Ggap.ConjugacyClasses()
        sage: for cl in cls:
        ....:      m = conjugacy_class_matrix(cls[2], 60)
        ....:      for _ in range(20):
        ....:           g = G.random_element().matrix()
        ....:           assert g*m*~g == m
    """
    res = [[0]*d for _ in range(d)]

    for p in cl.AsList():
        p = [i-1 for i in libgap.ListPerm(p, d).sage()]
        for k in range(d):
            res[k][p[k]] += 1

    from sage.matrix.constructor import matrix
    return matrix(res)

def isotypic_projection_matrix(G, d, chi, deg, conj_mats=None):
    r"""
    Return an isotypic projection matrix

    INPUT:

    - ``G`` -- a permutation group

    - ``d`` -- (integer) the domain of the group is `\{1, 2, \ldots, d\}`

    - ``chi`` -- (tuple) real or complex character

    - ``deg`` -- (integer) degree of the character

    Recall the formula for the projection as given in Theorem 8 in [Ser]_. If
    `G` is a permutation group, then
    
    .. MATH::
    
        \pi_\chi = \sum_{g \in G} \overline_{\chi(g)} g

    REFERENCES::

    .. [Ser] J.-P. Serre, "Représentation des groupes finis."

    EXAMPLES::

        sage: from surface_dynamics.misc.group_representation import real_characters, isotypic_projection_matrix
        sage: G = AlternatingGroup(5)
        sage: T,deg = real_characters(G)
        sage: isotypic_projection_matrix(G, 5, T[0], deg[0])
        [1/5 1/5 1/5 1/5 1/5]
        [1/5 1/5 1/5 1/5 1/5]
        [1/5 1/5 1/5 1/5 1/5]
        [1/5 1/5 1/5 1/5 1/5]
        [1/5 1/5 1/5 1/5 1/5]
        sage: isotypic_projection_matrix(G, 5, T[1], deg[1])
        [0 0 0 0 0]
        [0 0 0 0 0]
        [0 0 0 0 0]
        [0 0 0 0 0]
        [0 0 0 0 0]
        sage: isotypic_projection_matrix(G, 5, T[2], deg[2])
        [0 0 0 0 0]
        [0 0 0 0 0]
        [0 0 0 0 0]
        [0 0 0 0 0]
        [0 0 0 0 0]
        sage: isotypic_projection_matrix(G, 5, T[3], deg[3])
        [ 4/5 -1/5 -1/5 -1/5 -1/5]
        [-1/5  4/5 -1/5 -1/5 -1/5]
        [-1/5 -1/5  4/5 -1/5 -1/5]
        [-1/5 -1/5 -1/5  4/5 -1/5]
        [-1/5 -1/5 -1/5 -1/5  4/5]
        sage: isotypic_projection_matrix(G, 5, T[4], deg[4])
        [0 0 0 0 0]
        [0 0 0 0 0]
        [0 0 0 0 0]
        [0 0 0 0 0]
        [0 0 0 0 0]

        sage: sum(isotypic_projection_matrix(G, 5, T[i], deg[i]) for i in range(5)).is_one()
        True
    """
    from sage.matrix.special import zero_matrix
    res = zero_matrix(d)

    Ggap = libgap(G)

    for t,cl in enumerate(Ggap.ConjugacyClasses()):
        if conj_mats is None:
            m = conjugacy_class_matrix(cl, d)
        else:
            m = conj_mats[t]
        res += chi[t] * conjugacy_class_matrix(cl,d)

    return deg / G.cardinality() * res

def real_isotypic_projection_matrices(G, d):
    r"""
    Return the real projections

    EXAMPLES::

        sage: from surface_dynamics.misc.group_representation import real_isotypic_projection_matrices

        sage: G = AlternatingGroup(5)
        sage: mats = real_isotypic_projection_matrices(G, 6)
        sage: sum(mats).is_one()
        True

        sage: G = CyclicPermutationGroup(6)
        sage: mats = real_isotypic_projection_matrices(G, 6)
        sage: sum(mats).is_one()
        True

        sage: G = QuaternionGroup()
        sage: mats = real_isotypic_projection_matrices(G, 8)
        sage: sum(mats).is_one()
        True
    """
    return [isotypic_projection_matrix(G, d, chi, deg) for chi,deg in zip(*real_characters(G))]

