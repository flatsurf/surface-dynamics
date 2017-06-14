r"""
Some linear algebra routines
"""

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
        from sage.rings.rational_field import QQ
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
        from sage.rings.rational_field import QQ
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
    from sage.rings.all import QQ
    from sage.geometry.polyhedron.constructor import Polyhedron
    V = deformation_space(v)
    P = Polyhedron(lines=deformation_space(v).basis())
    B = Polyhedron(rays=(QQ**V.degree()).basis())
    return P.intersection(B)

