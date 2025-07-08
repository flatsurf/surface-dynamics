# ****************************************************************************
#       Copyright (C) 2011-2019 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cdef gl2z_orbits(origamis, int n, int limit)
cpdef sl2z_orbits(origamis, int n, int limit)

cdef class Origami_dense_pyx:
    cdef int * _r
    cdef int * _u
    cdef int _n
    cdef object _l_edges
    cdef object _i_edges
    cdef object _name
    cdef list _pos
    cdef set _rl_frontiers
    cdef set _tb_frontiers

    cdef Origami_dense_pyx _new_c(self, int *rr_and_uu)

    cpdef _set_standard_form(self, return_map=*)

    cpdef inverse(self)
    cpdef mirror(self)
    cpdef horizontal_twist(self, width=*, cylinder=*)
    cpdef vertical_twist(self, width=*, cylinder=*)

    cdef _compute_gl2z_edges(self)

    cdef pyx_lyapunov_exponents_approx(self,
                                       nb_iterations, nb_experiments,
                                       nb_vectors, only_mean)

    cdef pyx_lyapunov_exponents_approx_with_involution(
        self, involution, nb_iterations, nb_experiments,
        nb_vectors_p, nb_vectors_m, only_mean)


cdef class PillowcaseCover_dense_pyx:
    cdef int *_g
    cdef int _n
    cdef object _l_edges
    cdef object _i_edges

    cdef PillowcaseCover_dense_pyx _new_c(self, int *g)
