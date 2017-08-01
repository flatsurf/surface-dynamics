r"""
Some extra ppl utilities

ppl is a library to deal with rational polytopes. pplpy is a thin
Python wrapper. This file provides some helper functions to use it.
"""

import ppl

def ppl_check_non_negative_cone(C):
    r"""
    Check whether the input cone ``C`` is a subcone of the positive cone

    EXAMPLES::

        sage: from surface_dynamics.misc.ppl_utils import ppl_check_non_negative_cone # optional - pplpy
        sage: import ppl    # optional - pplpy

        sage: gs = ppl.Generator_System()     # optional - pplpy
        sage: gs.insert(ppl.point())          # optional - pplpy
        sage: gs.insert(ppl.ray(ppl.Variable(0) + ppl.Variable(2))) # optional - pplpy
        sage: ppl_check_non_negative_cone(ppl.C_Polyhedron(gs)) # optional - pplpy

        sage: gs.insert(ppl.ray(-ppl.Variable(0) - ppl.Variable(2))) # optional - pplpy
        sage: ppl_check_non_negative_cone(ppl.C_Polyhedron(gs)) # optional - pplpy
        Traceback (most recent call last):
        ...
        ValueError: C must be a subcone of the non-negative cone

        sage: gs = ppl.Generator_System() # optional - pplpy
        sage: gs.insert(ppl.point(ppl.Variable(2))) # optional - pplpy
        sage: ppl_check_non_negative_cone(ppl.C_Polyhedron(gs)) # optional - pplpy
        Traceback (most recent call last):
        ...
        ValueError: the cone does not contain zero
        sage: gs.insert(ppl.point()) # optional - pplpy
        sage: ppl_check_non_negative_cone(ppl.C_Polyhedron(gs)) # optional - pplpy
        Traceback (most recent call last):
        ...
        ValueError: should have only zero as vertices
    """
    if not isinstance(C, ppl.C_Polyhedron):
        raise ValueError('C must be a polyhedron in the right ambient space')

    n = C.space_dimension()
    if not C.contains(ppl.C_Polyhedron(ppl_zero_point(n))):
        raise ValueError('the cone does not contain zero')
    if not ppl_positive_cone(n).contains(C):
        raise ValueError('C must be a subcone of the non-negative cone')
    for g in C.generators():
        if g.is_point() and not g.is_equivalent_to(ppl_zero_point(n)):
            raise ValueError('should have only zero as vertices'.format(g))
        if g.is_line():
            raise ValueError('the cone contains a line')

def ppl_zero_point(n):
    r"""
    Return the origin in R^n

    EXAMPLES::

        sage: from surface_dynamics.misc.ppl_utils import ppl_zero_point # optional - pplpy
        sage: ppl_zero_point(3)   # optional - pplpy
        point(0/1, 0/1, 0/1)
        sage: ppl_zero_point(5)   # optional - pplpy
        point(0/1, 0/1, 0/1, 0/1, 0/1)
    """
    z = ppl.point()
    z.set_space_dimension(n)
    return z

def ppl_positive_cone(n):
    r"""
    Return the positive cone in R^n

    EXAMPLES::

        sage: from surface_dynamics.misc.ppl_utils import ppl_positive_cone # optional - pplpy
        sage: C = ppl_positive_cone(3)  # optional - pplpy
        sage: C.minimized_generators()  # optional - pplpy
        Generator_System {point(0/1, 0/1, 0/1), ray(0, 0, 1), ray(0, 1, 0), ray(1, 0, 0)}
    """
    gs = ppl.Generator_System(ppl_zero_point(n))
    l = [0]*n
    for i in range(n):
        gs.insert(ppl.ray(ppl.Variable(i)))
    return ppl.C_Polyhedron(gs)

def ppl_cone(rays):
    r"""
    Convert the list ``rays`` into a ppl cone

    EXAMPLES::

        sage: from surface_dynamics.misc.ppl_utils import ppl_cone  # optional - pplpy
        sage: C = ppl_cone([(0,1,2),(1,1,1),(1,0,1)])   # optional - pplpy
        sage: C.minimized_generators()   # optional - pplpy
        Generator_System {point(0/1, 0/1, 0/1), ray(0, 1, 2), ray(1, 0, 1), ray(1, 1, 1)}
    """
    n = len(rays[0])
    gs = ppl.Generator_System(ppl_zero_point(n))
    for r in rays:
        gs.insert(ppl.ray(sum(j * ppl.Variable(i) for i,j in enumerate(r))))
    return ppl.C_Polyhedron(gs)


def ppl_convert(P):
    r"""
    Convert a Sage polyhedron to a ppl polyhedron

    EXAMPLES::

        sage: from surface_dynamics.misc.ppl_utils import ppl_convert  # optional - pplpy
        sage: P = ppl_convert(Polyhedron(vertices=[(0,1,0),(1,0,1)], rays=[(0,0,1),[3,2,1]]))  # optional - pplpy
        sage: P.minimized_generators()  # optional - pplpy
        Generator_System {ray(0, 0, 1), point(0/1, 1/1, 0/1), point(1/1, 0/1, 1/1), ray(3, 2, 1)}
    """
    if isinstance(P, ppl.C_Polyhedron):
        return P
    gs = ppl.Generator_System()
    for v in P.vertices_list():
        gs.insert(ppl.point(sum(j * ppl.Variable(i) for i,j in enumerate(v))))
    for r in P.rays_list():
        gs.insert(ppl.ray(sum(j * ppl.Variable(i) for i,j in enumerate(r))))
    for l in P.lines_list():
        gs.insert(ppl.line(sum(j * ppl.Variable(i) for i,j in enumerate(l))))
    return ppl.C_Polyhedron(gs)

