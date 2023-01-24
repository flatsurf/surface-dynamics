r"""
Some linear algebra routines
"""
#*****************************************************************************
#       Copyright (C) 2019-2023 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.arith.all import gcd, lcm
from sage.arith.misc import binomial
from sage.rings.all import ZZ, QQ
from sage.matrix.constructor import matrix
from sage.geometry.polyhedron.constructor import Polyhedron


def relation_space(v):
    r"""
    Relation space of the given vector ``v``

    This is the sub vector space of `\QQ^d` given as the kernel of the map
    `n \mapsto n \cdot \lambda`. The dimension is `d - rank`.

    EXAMPLES::

        sage: from surface_dynamics.misc.linalg import relation_space
        sage: K.<sqrt2> = QuadraticField(2)
        sage: v3 = vector([sqrt2, 1, 1+sqrt2])
        sage: relation_space(v3)
        Vector space of degree 3 and dimension 1 over Rational Field
        Basis matrix:
        [ 1  1 -1]
        sage: v4 = vector([sqrt2, 1, 1+sqrt2, 1-sqrt2])
        sage: relation_space(v4)
        Vector space of degree 4 and dimension 2 over Rational Field
        Basis matrix:
        [   1    0 -1/2  1/2]
        [   0    1 -1/2 -1/2]

        sage: v = vector([1,2,5,3])
        sage: relation_space(v)
        Vector space of degree 4 and dimension 3 over Rational Field
        Basis matrix:
        [   1    0    0 -1/3]
        [   0    1    0 -2/3]
        [   0    0    1 -5/3]

    The relation space has some covariance relation with respect to matrix
    actions::

        sage: m3 = matrix(3, [1,-1,0,2,-3,4,5,-2,2])
        sage: relation_space(v3 * m3) == relation_space(v3) * ~m3.transpose()
        True
        sage: relation_space(m3 * v3) == relation_space(v3) * ~m3
        True

        sage: m4 = matrix(4, [1,-1,0,1,2,-3,0,4,5,3,-2,2,1,1,1,1])
        sage: relation_space(v4 * m4) == relation_space(v4) * ~m4.transpose()
        True
        sage: relation_space(m4 * v4) == relation_space(v4) * ~m4
        True
    """
    from sage.matrix.constructor import matrix
    try:
        m_lengths = matrix([u.vector() for u in v])
    except AttributeError:
        v = [QQ.coerce(i) for i in v]
        m_lengths = matrix([[i] for i in v])
    return m_lengths.left_kernel()

def deformation_space(lengths):
    r"""
    Deformation space of the given ``lengths``

    This is the smallest vector space defined over `\QQ` that contains
    the vector ``lengths``. Its dimension is `rank`.

    EXAMPLES::

        sage: from surface_dynamics.misc.linalg import deformation_space
        sage: K.<sqrt2> = QuadraticField(2)
        sage: v3 = vector([sqrt2, 1, 1+sqrt2])
        sage: deformation_space(v3)
        Vector space of degree 3 and dimension 2 over Rational Field
        Basis matrix:
        [1 0 1]
        [0 1 1]
        sage: v4 = vector([sqrt2, 1, 1+sqrt2, 1-sqrt2])
        sage: deformation_space(v4)
        Vector space of degree 4 and dimension 2 over Rational Field
        Basis matrix:
        [ 1  0  1 -1]
        [ 0  1  1  1]

        sage: v = vector([1, 5, 2, 9])
        sage: deformation_space(v)
        Vector space of degree 4 and dimension 1 over Rational Field
        Basis matrix:
        [1 5 2 9]

    The deformation space has some covariance relation with respect to matrix
    actions::

        sage: m3 = matrix(3, [1,-1,0,2,-3,4,5,-2,2])
        sage: deformation_space(v3 * m3) == deformation_space(v3) * m3
        True
        sage: deformation_space(m3 * v3) == deformation_space(v3) * m3.transpose()
        True

        sage: m4 = matrix(4, [1,-1,0,1,2,-3,0,4,5,3,-2,2,1,1,1,1])
        sage: deformation_space(v4 * m4) == deformation_space(v4) * m4
        True
        sage: deformation_space(m4 * v4) == deformation_space(v4) * m4.transpose()
        True
    """
    from sage.matrix.constructor import matrix
    try:
        m_lengths = matrix([u.vector() for u in lengths])
    except AttributeError:
        lengths = [QQ.coerce(i) for i in lengths]
        m_lengths = matrix([[i] for i in lengths])

    return m_lengths.column_space()

def deformation_cone(v):
    r"""
    Return the deformation cone of the given vector ``v``

    EXAMPLES::

        sage: from surface_dynamics.misc.linalg import deformation_cone
        sage: K.<sqrt2> = QuadraticField(2)
        sage: v3 = vector([sqrt2, 1, 1+sqrt2])
        sage: P = deformation_cone(v3)
        sage: P
        A 2-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex and 2 rays
        sage: P.rays_list()
        [[1, 0, 1], [0, 1, 1]]
    """
    V = deformation_space(v)
    P = Polyhedron(lines=deformation_space(v).basis())
    B = Polyhedron(rays=(QQ**V.degree()).basis())
    return P.intersection(B)

def cone_triangulate(C, hyperplane=None):
    r"""
    Triangulation of rational cone contained in the positive quadrant.

    EXAMPLES::

        sage: from surface_dynamics.misc.linalg import cone_triangulate
        sage: P = Polyhedron(rays=[(1,0,0),(0,1,0),(1,0,1),(0,1,1)])
        sage: list(cone_triangulate(P)) # random
        [[(0, 1, 1), (0, 1, 0), (1, 0, 0)], [(0, 1, 1), (1, 0, 1), (1, 0, 0)]]
        sage: len(_)
        2

        sage: rays = [(0, 1, 0, -1, 0, 0),
        ....: (1, 0, -1, 0, 0, -1),
        ....: (0, 1, -1, 0, 0, -1),
        ....: (0, 0, 1, 0, 0, 0),
        ....: (0, 0, 0, 1, 0, -1),
        ....: (1, -1, 0, 0, 1, -1),
        ....: (0, 0, 0, 0, 1, -1),
        ....: (0, 0, 1, -1, 1, 0),
        ....: (0, 0, 1, -1, 0, 0),
        ....: (0, 0, 1, 0, -1, 0),
        ....: (0, 0, 0, 1, -1, -1),
        ....: (1, -1, 0, 0, 0, -1),
        ....: (0, 0, 0, 0, 0, -1)]
        sage: P = Polyhedron(rays=rays)
        sage: list(cone_triangulate(P, hyperplane=(1, 2, 3, -1, 0, -5))) # random
        [[(0, 0, 0, 0, 0, -1),
          (0, 0, 0, 0, 1, -1),
          (0, 0, 0, 1, -1, -1),
          (0, 0, 1, 0, 0, 0),
          (0, 1, -1, 0, 0, -1),
          (1, -1, 0, 0, 1, -1)],
          ...
          (0, 0, 1, 0, 0, 0),
          (0, 1, -1, 0, 0, -1),
          (0, 1, 0, -1, 0, 0),
          (1, -1, 0, 0, 1, -1),
          (1, 0, -1, 0, 0, -1)]]
        sage: len(_)
        16
    """
    rays = [r.vector() for r in C.rays()]
    dim = len(rays[0])
    if hyperplane is None:
        hyperplane = [1] * dim
    scalings = [sum(x*h for x,h in zip(r, hyperplane)) for r in rays]
    assert all(s > 0 for s in scalings)
    normalized_rays = [r / s for r,s in zip(rays, scalings)]
    P = Polyhedron(vertices=normalized_rays)
    for t in P.triangulate():
        simplex = [P.Vrepresentation(i).vector() for i in t]
        yield [(r / gcd(r)).change_ring(ZZ) for r in simplex]

def symbolic_matrix_power(M, n):
    r"""
    Return the symbolic power ``M^n`` of the unipotent matrix ``M``.

    EXAMPLES::

        sage: from surface_dynamics.misc.linalg import symbolic_matrix_power
        sage: m = matrix(3, [1,1,1,0,1,1,0,0,1])
        sage: n = polygen(QQ, 'n')
        sage: symbolic_matrix_power(m, n)
        [              1               n 1/2*n^2 + 1/2*n]
        [              0               1               n]
        [              0               0               1]

        sage: m = matrix(2, [2,1,1,1])
        sage: symbolic_matrix_power(m, n)
        Traceback (most recent call last):
        ...
        NotImplementedError: power only implemented for unipotent matrices
    """
    d = M.nrows()
    I = M.parent().identity_matrix()
    N = M - M.parent().identity_matrix()
    char = N.charpoly()
    if any(char[i] for i in range(d)):
        raise NotImplementedError('power only implemented for unipotent matrices')

    result = I
    P = N
    p = 1
    while P:
        result += binomial(n, p) * P
        P *= N
        p += 1
    return result


def disjoint_vectors(vectors, min_size=2, max_size=None):
    r"""
    EXAMPLES::

        sage: from surface_dynamics.misc.linalg import disjoint_vectors
        sage: V = FreeModule(ZZ, 4)
        sage: v0 = V((1,0,0,0))
        sage: v1 = V((1,1,0,0))
        sage: v2 = V((0,0,1,0))
        sage: v3 = V((0,1,1,1))
        sage: v4 = V((0,1,0,1))
        sage: for subset, s in disjoint_vectors((v0,v1,v2,v3,v4)):
        ....:     print(subset, s)
        [0, 2] (1, 0, 1, 0)
        [0, 2, 4] (1, 1, 1, 1)
        [0, 3] (1, 1, 1, 1)
        [0, 4] (1, 1, 0, 1)
        [1, 2] (1, 1, 1, 0)
        [2, 4] (0, 1, 1, 1)
        sage: for subset, s in disjoint_vectors((v0,v1,v2,v3,v4), min_size=2, max_size=2):
        ....:     print(subset, s)
        [0, 2] (1, 0, 1, 0)
        [0, 3] (1, 1, 1, 1)
        [0, 4] (1, 1, 0, 1)
        [1, 2] (1, 1, 1, 0)
        [2, 4] (0, 1, 1, 1)
    """
    if max_size is None:
        max_size = len(vectors)
    F = vectors[0].parent()
    n = len(vectors)
    d = F.dimension()
    stack = []
    i = 0
    c = F.zero()  # current sum
    while True:
        while i < len(vectors) and len(stack) < max_size:
            cc = c + vectors[i]
            if all(x < 2 for x in cc):
                stack.append(i)
                c = cc
                if len(stack) >= min_size:
                    c.set_immutable()
                    yield stack, c
            i += 1
        if not stack:
            return
        i = stack.pop()
        c -= vectors[i]
        i += 1


def linearly_independent_vectors(vectors, min_size=0, max_size=None):
    r"""
    Iterate through the subsets of ``vectors`` made of linearly independent vectors.

    EXAMPLES::

        sage: from surface_dynamics.misc.linalg import linearly_independent_vectors
        sage: V = FreeModule(ZZ, 4)
        sage: vecs = [V((1, 0, 0, 0)), V((1, 1, 0, 0)), V((0, 0, 1, 0)), V((0, 1, 1, 1)), V((0, 1, 0, 1))]
        sage: for subset, s in linearly_independent_vectors(vecs):
        ....:     print(subset, s.rows())
        [0] [(1, 0, 0, 0)]
        [0, 1] [(1, 0, 0, 0), (1, 1, 0, 0)]
        [0, 1, 2] [(1, 0, 0, 0), (1, 1, 0, 0), (0, 0, 1, 0)]
        [0, 1, 2, 3] [(1, 0, 0, 0), (1, 1, 0, 0), (0, 0, 1, 0), (0, 1, 1, 1)]
        [0, 1, 2, 4] [(1, 0, 0, 0), (1, 1, 0, 0), (0, 0, 1, 0), (0, 1, 0, 1)]
        [0, 1, 3] [(1, 0, 0, 0), (1, 1, 0, 0), (0, 1, 1, 1)]
        [0, 1, 3, 4] [(1, 0, 0, 0), (1, 1, 0, 0), (0, 1, 1, 1), (0, 1, 0, 1)]
        [0, 1, 4] [(1, 0, 0, 0), (1, 1, 0, 0), (0, 1, 0, 1)]
        [0, 2] [(1, 0, 0, 0), (0, 0, 1, 0)]
        [0, 2, 3] [(1, 0, 0, 0), (0, 0, 1, 0), (0, 1, 1, 1)]
        [0, 2, 4] [(1, 0, 0, 0), (0, 0, 1, 0), (0, 1, 0, 1)]
        [0, 3] [(1, 0, 0, 0), (0, 1, 1, 1)]
        [0, 3, 4] [(1, 0, 0, 0), (0, 1, 1, 1), (0, 1, 0, 1)]
        [0, 4] [(1, 0, 0, 0), (0, 1, 0, 1)]
        [1] [(1, 1, 0, 0)]
        [1, 2] [(1, 1, 0, 0), (0, 0, 1, 0)]
        [1, 2, 3] [(1, 1, 0, 0), (0, 0, 1, 0), (0, 1, 1, 1)]
        [1, 2, 4] [(1, 1, 0, 0), (0, 0, 1, 0), (0, 1, 0, 1)]
        [1, 3] [(1, 1, 0, 0), (0, 1, 1, 1)]
        [1, 3, 4] [(1, 1, 0, 0), (0, 1, 1, 1), (0, 1, 0, 1)]
        [1, 4] [(1, 1, 0, 0), (0, 1, 0, 1)]
        [2] [(0, 0, 1, 0)]
        [2, 3] [(0, 0, 1, 0), (0, 1, 1, 1)]
        [2, 4] [(0, 0, 1, 0), (0, 1, 0, 1)]
        [3] [(0, 1, 1, 1)]
        [3, 4] [(0, 1, 1, 1), (0, 1, 0, 1)]
        [4] [(0, 1, 0, 1)]

        sage: for subset, s in linearly_independent_vectors(vecs, min_size=2, max_size=2):
        ....:     print(subset)
        [0, 1]
        [0, 2]
        [0, 3]
        [0, 4]
        [1, 2]
        [1, 3]
        [1, 4]
        [2, 3]
        [2, 4]
        [3, 4]
        sage: for subset, s in linearly_independent_vectors(vecs, min_size=3, max_size=3):
        ....:     print(subset)
        [0, 1, 2]
        [0, 1, 3]
        [0, 1, 4]
        [0, 2, 3]
        [0, 2, 4]
        [0, 3, 4]
        [1, 2, 3]
        [1, 2, 4]
        [1, 3, 4]
    """
    F = vectors[0].parent()
    n = len(vectors)
    d = F.dimension()
    if max_size is None:
        max_size = d
    stack = []
    i = 0
    M = matrix(ZZ, max_size, d) # current matrix
    while True:
        while len(stack) < max_size and i < len(vectors):
            M[len(stack)] = vectors[i]
            if M.rank() == len(stack) + 1:
                stack.append(i)
                if len(stack) >= min_size:
                    yield stack, M[:len(stack)]
            else:
                M[len(stack)] = F.zero()
            i += 1
        if not stack:
            return
        i = stack.pop()
        M[len(stack)] = F.zero()
        i += 1
