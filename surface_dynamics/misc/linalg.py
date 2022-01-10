r"""
Some linear algebra routines
"""
#*****************************************************************************
#       Copyright (C) 2019-2021 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.all import ZZ, QQ, vector, Compositions, gcd, lcm, Polyhedron
from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException

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

def pivots(n, length, force_first=False):
    r"""
    Iterator through increasing sequences of given ``length` in ``{0, 1, ..., n-1}``.

    INPUT:

    - ``force_first`` - whether the first pivot is on the first column

    EXAMPLES::

        sage: from surface_dynamics.misc.linalg import pivots
        sage: list(pivots(4, 3, force_first=False))
        [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
        sage: list(pivots(4, 3, force_first=True))
        [[0, 1, 2], [0, 1, 3], [0, 2, 3]]
    """
    if force_first:
        for gaps in reversed(Compositions(n, length=length)):
            gaps = [1] + list(gaps)
            pivots = [gaps[0] - 1]
            for i in gaps[1:-1]:
                pivots.append(pivots[-1] + i)
            yield pivots
    else:
        for gaps in reversed(Compositions(n + 1, length=length+1)):
            pivots = [gaps[0] - 1]
            for i in gaps[1:-1]:
                pivots.append(pivots[-1] + i)
            yield pivots

def isotropic_subspaces(B, d, bound=1, contains_positive_vector=False):
    r"""
    Iterator over rational isotropic subspaces of dimension ``d`` of the
    anti-symmetric bilinear form ``B``.

    INPUT:

    - ``B`` -- anti-symmetric matrix with integral coefficients

    - ``d`` -- dimension of the isotropic subspace

    - ``bound`` (optional, default 1) -- bound on the coefficients for generators (in echelon
      form)

    - ``contains_positive_vector`` (optional, default ``False``) - if ``True``,
      then only iterate through subspaces that contain a positive vector

    EXAMPLES::

        sage: from surface_dynamics.misc.linalg import isotropic_subspaces

        sage: B = matrix(4, [0, 1, 1, 1, -1, 0, 1, 1, -1, -1, 0, 1, -1, -1, -1, 0])
        sage: for vectors in isotropic_subspaces(B, 2, bound=2, contains_positive_vector=True):
        ....:     print(matrix(vectors))
        ....:     print('*' * 15)
        [1 0 0 1]
        [0 1 1 0]
        ***************
        [1 0 0 1]
        [0 1 2 0]
        ***************
        [1 0 0 2]
        [0 1 1 2]
        ***************
        [ 1  0  1  0]
        [ 0  1 -2  1]
        ***************
        [1 0 1 2]
        [0 1 0 1]
        ***************
        [1 0 1 2]
        [0 1 2 2]
        ***************
        [ 1  0  1  2]
        [ 0  1 -2  0]
        ***************
        [ 1  0  2  0]
        [ 0  1 -2  1]
        ***************
        [1 0 2 2]
        [0 1 0 1]
        ***************
        [ 1  0 -2  1]
        [ 0  1  1  2]
        ***************
        [ 1  0 -2  1]
        [ 0  1  2  2]
        ***************
        [ 1  0 -2  2]
        [ 0  1  1  0]
        ***************
        [ 1  0 -2  2]
        [ 0  1  2 -1]
        ***************

    Some countings::

        sage: B = matrix(4, [0,1,1,1,-1,0,1,1,-1,-1,0,1,-1,-1,-1,0])
        sage: sum(1 for _ in isotropic_subspaces(B, 2, bound=1, contains_positive_vector=True))
        1
        sage: sum(1 for _ in isotropic_subspaces(B, 2, bound=2, contains_positive_vector=True))
        13
        sage: sum(1 for _ in isotropic_subspaces(B, 2, bound=1, contains_positive_vector=False))
        10
        sage: sum(1 for _ in isotropic_subspaces(B, 2, bound=2, contains_positive_vector=False))
        66

        sage: B = matrix(5, [0,1,1,1,1,-1,0,1,1,1,-1,-1,0,1,1,-1,-1,-1,0,1,-1,-1,-1,-1,0])
        sage: sum(1 for _ in isotropic_subspaces(B, 2, bound=1, contains_positive_vector=True))
        17
        sage: sum(1 for _ in isotropic_subspaces(B, 3, bound=1, contains_positive_vector=True))
        0
        sage: sum(1 for _ in isotropic_subspaces(B, 2, bound=2, contains_positive_vector=True))
        204
        sage: sum(1 for _ in isotropic_subspaces(B, 2, bound=1, contains_positive_vector=False))
        175
        sage: sum(1 for _ in isotropic_subspaces(B, 3, bound=1, contains_positive_vector=False))
        2
        sage: sum(1 for _ in isotropic_subspaces(B, 2, bound=2, contains_positive_vector=False))
        1649
    """
    assert B.is_square()
    n = B.nrows()

    for piv in pivots(n, d, force_first=contains_positive_vector):
        positions_to_be_filled = []
        for i in range(d):
            positions_to_be_filled.append(list(enumerate(sorted(set(range(piv[i]+1, n)).difference(piv[i+1:])))))

        # backtracking
        # the ordering on entries is 0, 1, -1, 2, -2, 3, -3, ....
        vectors = [vector(ZZ, n) for _ in range(d)]
        for i, j in enumerate(piv):
            vectors[i][j] = 1
        products = [vectors[i] * B for i in range(d)]

        if (all(products[i].dot_product(vectors[j]).is_zero() for i in range(d) for j in range(i)) and
            (not contains_positive_vector or has_positive_linear_combination(vectors, n))):
            yield vectors

        i = 1
        while True:
            # find the next vector at position i
            v = vectors[i]
            found_next = False
            for jpos, j in reversed(positions_to_be_filled[i]):
                if v[j] == -bound:
                    v[j] = 0
                elif v[j] > 0:
                    v[j] = -v[j]
                    found_next = True
                    break
                else:
                    v[j] = -v[j] + 1
                    found_next = True
                    break

            # backtrack ?
            if not found_next:
                i -= 1
                if i == -1:
                    break
                continue

            # go for the next vector if span(vectors) has positive elements and v[i] is orthogonal to the i-1 first vectors
            products[i] = v * B
            if all(products[i].dot_product(vectors[j]).is_zero() for j in range(i)):
                if i == d - 1:
                    if not contains_positive_vector or has_positive_linear_combination(vectors, n):
                        yield vectors
                else:
                    i = i + 1

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

def has_positive_linear_combination(vectors, d, solver="PPL"):
    r"""
    Test whether the given ``vectors`` admit a linear combination that is positive.

    EXAMPLES::

        sage: from surface_dynamics.misc.linalg import has_positive_linear_combination
        sage: has_positive_linear_combination([(1,0), (0,1)], 2)
        True
        sage: has_positive_linear_combination([(1,-1)], 2)
        False
        sage: vectors = [(1, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, -1), (0, 0, 1, 0, -1, 0)]
        sage: has_positive_linear_combination(vectors, 3)
        True
        sage: has_positive_linear_combination(vectors, 4)
        False
    """
    M = MixedIntegerLinearProgram(solver=solver)
    x = M.new_variable()
    for i in range(d):
        M.add_constraint(M.sum(x[j] * v[i] for j, v in enumerate(vectors)) >= 1)
        try:
            M.solve(objective_only=True)
        except MIPSolverException:
            return False
    return True
