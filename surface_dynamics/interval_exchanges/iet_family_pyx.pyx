# distutils: language=c++
# distutils: extra_compile_args=-std=c++11
# distutils: libraries=ppl gmp
#*****************************************************************************
#       Copyright (C) 2019 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

##########################
# C declarations
##########################

cdef extern from "gmp.h":
    ctypedef unsigned long mp_limb_t
    ctypedef long mp_size_t

    ctypedef struct __mpz_struct:
        pass
    ctypedef __mpz_struct *mpz_srcptr
    ctypedef __mpz_struct mpz_t[1]

    void mpz_init (mpz_t integer) nogil
    void mpz_clear (mpz_t integer) nogil

    void mpz_set (mpz_t rop, mpz_t op)
    void mpz_set_ui (mpz_t rop, unsigned long int op) nogil
    void mpz_set_si (mpz_t rop, signed long int op) nogil
    void mpz_add (mpz_t rop, mpz_t op1, mpz_t op2) nogil
    int mpz_cmp (mpz_t op1, mpz_t op2) nogil
    int mpz_sgn (mpz_t op)

    size_t mpz_size (mpz_t op)
    mp_limb_t mpz_getlimbn (mpz_t op, mp_size_t n)

########################
# Cython imports
########################

from libc.stdlib cimport malloc, calloc, realloc, free, qsort
from libc.string cimport memcpy

from cpython.object cimport Py_EQ, Py_NE
from cython.operator cimport dereference as deref

from cysignals.memory cimport sig_malloc, sig_calloc, sig_realloc, sig_free, check_malloc, check_calloc, check_realloc

from ppl.ppl_decl cimport PPL_Variable, PPL_Generator, PPL_Generator_System, PPL_gs_iterator
from ppl.linear_algebra cimport Variable
from ppl.generator cimport Generator_System
from ppl.polyhedron cimport C_Polyhedron

from sage.ext.stdsage cimport PY_NEW
from sage.rings.integer cimport Integer
from sage.rings.rational_field import QQ

#from gmpy2 cimport import_gmpy2, MPZ_Object, MPZ, GMPy_MPZ_New
#import_gmpy2()


########################
# Beginning of the code
########################

cdef int compare(const void * _left, const void * _right) nogil:
    cdef mpz_t * left = (<mpz_t **> _left)[0]
    cdef mpz_t * right = (<mpz_t **> _right)[0]
    cdef int c = mpz_cmp(left[0], right[0])
    while c == 0:
        left += 1
        right += 1
        c = mpz_cmp(left[0], right[0])
    return c

cdef inline int mpz_vec_any_smaller(mpz_t * x, mpz_t * y, size_t length):
    r"""
    Check whether any x[i] < y[i]
    """
    cdef size_t k
    for k in range(length):
        if mpz_cmp(x[k], y[k]) < 0:
            return 1
    return 0

cdef inline int mpz_vec_equal(mpz_t * x, mpz_t * y, size_t length):
    cdef size_t k
    for k in range(length):
        if mpz_cmp(x[k], y[k]): return 0
    return 1

cdef Py_hash_t mpz_pythonhash(mpz_srcptr z):
    """
    Hash an ``mpz``, where the hash value is the same as the hash value
    of the corresponding Python ``long``, except that we do not replace
    -1 by -2 (the Cython wrapper for ``__hash__`` does that).
    """
    # Add all limbs, adding 1 for every carry
    cdef mp_limb_t h1 = 0
    cdef mp_limb_t h0
    cdef size_t i, n
    n = mpz_size(z)
    for i in range(n):
        h0 = h1
        h1 += mpz_getlimbn(z, i)
        # Add 1 on overflow
        if h1 < h0: h1 += 1

    cdef Py_hash_t h = h1
    if mpz_sgn(z) < 0:
        return -h
    return h

cdef inline void rauzy_top(int * perm, size_t dim):
    cdef int a = perm[dim - 1]
    cdef size_t i, j
    for i in range(dim, 2*dim):
        if perm[i] == a:
            break
    a = perm[2 * dim - 1]
    for j in range(2 * dim - 1, i + 1, -1):
        perm[j] = perm[j - 1]
    perm[i + 1] = a

cdef inline void rauzy_bot(int * perm, size_t dim):
    cdef int a = perm[2*dim - 1]
    cdef size_t i, j
    for i in range(dim):
        if perm[i] == a:
            break
    a = perm[dim - 1]
    for j in range(dim - 1, i + 1, -1):
        perm[j] = perm[j - 1]
    perm[i + 1] = a

# TODO: here we would better store a linear space rather than a positive cone
# (of course, it will be useful to know about the intersection of this space
# with the positive cone)
cdef class IETFamily_pyx:
    r"""
    A linear family of interval exchange transformations

    This class consists of two attributes:

    - a permutation (must be labeled)

    - a cone of lengths (a subcone of the positive cone)

    EXAMPLES::

        sage: from surface_dynamics import *

        sage: p = iet.Permutation('a b c d', 'd c b a')
        sage: F = iet.IETFamily(p, (ZZ**4).basis())
        sage: F
        Linear iet family of dimension 4 in RR^4
        top a b c d
        bot d c b a
        0 0 0 1
        0 0 1 0
        0 1 0 0
        1 0 0 0

        sage: iet.IETFamily(p, [(2,3,0,0), (0,1,1,1)])
        Linear iet family of dimension 2 in RR^4
        top a b c d
        bot d c b a
        0 1 1 1
        2 3 0 0

        sage: iet.IETFamily(p, [(1,0,0,0), (1,-1,1,1)])
        Traceback (most recent call last):
        ...
        ValueError: C must be a subcone of the non-negative cone
        sage: iet.IETFamily(p, Polyhedron(vertices=[(0,0,0,0), (1,0,0,0),(0,1,0,0),(1,1,1,1)]))
        Traceback (most recent call last):
        ...
        ValueError: should have only zero as vertices
    """
    cdef int * perm         # top perm
    cdef size_t dim         # dimension (= number of letters)
    cdef mpz_t * entries    # ray coefficients
    cdef mpz_t ** rows      # ray pointers to entries (beginning of vectors)
    cdef size_t alloc       # ray allocation size
    cdef size_t length      # number of rays
    cdef C_Polyhedron poly  # the polyhedron (that we should get rid of)
    cdef _free_module

    def __cinit__(self):
        self.perm = NULL
        self.entries = NULL
        self.rows = NULL
        self.alloc = 0
        self.dim = 0
        self.length = 0

    def __dealloc__(self):
        cdef size_t i
        for i in range(self.dim * self.alloc):
            mpz_clear(self.entries[i])
        free(self.perm)
        free(self.entries)
        free(self.rows)

    def _new(self):
        return IETFamily_pyx.__new__(IETFamily_pyx)

    def __init__(self, p, C):
        # convert and check input
        from surface_dynamics.interval_exchanges.labelled import LabelledPermutationIET
        if not isinstance(p, LabelledPermutationIET):
            raise ValueError('p must be a labelled permutation')

        if not isinstance(C, C_Polyhedron):
            try:
                C = C._ppl_polyhedron
            except AttributeError:
                from sage.geometry.polyhedron.constructor import Polyhedron
                try:
                    C = Polyhedron(C, base_ring=QQ)
                    C = C._ppl_polyhedron
                except Exception:
                    raise TypeError('invalid input')

        cdef C_Polyhedron CC = <C_Polyhedron> C

        from surface_dynamics.misc.ppl_utils import ppl_check_non_negative_cone
        ppl_check_non_negative_cone(C)
        if C.space_dimension() != len(p):
            raise ValueError('dimension mismatch')

        # fill C data from input
        self.dim = len(p)
        self.perm = <int *> check_malloc(2 * self.dim * sizeof(int))
        cdef size_t i
        cdef int j
        for i,j in enumerate(p._labels[0]):
            self.perm[i] = j
        for i,j in enumerate(p._labels[1]):
            self.perm[self.dim + i] = j

        self.poly = <C_Polyhedron> C
        self.fill_rays(self.poly.thisptr.minimized_generators())

    def free_module(self):
        if self._free_module is None:
            from sage.modules.free_module import FreeModule
            from sage.rings.rational_field import QQ
            self._free_module = FreeModule(QQ, self.ambient_dimension())
        return self._free_module

    def ray_coefficient(self, size_t nray, size_t i):
        r"""
        Return the ``i``-th coefficient of the ``nray``-th ray

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: F = iet.IETFamily(p, [(2,3,0,0), (0,1,1,1)])
            sage: F.ray_coefficient(0, 2)
            1
            sage: F.ray_coefficient(1, 0)
            2
            sage: F.ray_coefficient(3, 0)
            Traceback (most recent call last):
            ...
            IndexError: row index out of range
            sage: F.ray_coefficient(1, 5)
            Traceback (most recent call last):
            ...
            IndexError: column index out of range
        """
        if nray >= self.length:
            raise IndexError('row index out of range')
        if i >= self.dim:
            raise IndexError('column index out of range')
        cdef Integer z= PY_NEW(Integer)
        mpz_set(z.value, self.rows[nray][i])
        return z

    def _debug_info(self):
        import sys
        sys.stdout.write('dim      %lu\n' % self.dim)
        sys.stdout.write('length   %lu\n' % self.length)
        sys.stdout.write('alloc    %lu\n' % self.alloc)
        sys.stdout.flush()
        sys.stdout.write('entries  %lu\n' % <size_t> self.entries)
        sys.stdout.flush()
        sys.stdout.write('rows     %lu\n' % <size_t> self.rows)
        sys.stdout.flush()
        cdef size_t i
        for i in range(self.alloc):
            sys.stdout.write('rows[%d]  %lu\n' %(i, <size_t> (self.rows[i] - self.entries)))
            sys.stdout.flush()


    cdef void fit_length(self, size_t length):
        r"""
        Manage allocation

        After this function is called, it is ensured that there is enough space
        for ``length`` rays.
        """
        cdef size_t i
        cdef mpz_t * old_entries
        cdef ptrdiff_t offset
        if length > self.alloc:
            length = max(length, 2 * self.alloc)

            old_entries = self.entries
            self.entries = <mpz_t *> realloc(self.entries, self.dim * length * sizeof(mpz_t))
            offset = self.entries - old_entries

            for i in range(self.dim * self.alloc, self.dim * length):
                mpz_init(self.entries[i])

            self.rows = <mpz_t **> realloc(self.rows, length * sizeof(mpz_t *))
            if offset:
                for i in range(self.alloc):
                    self.rows[i] += offset
            for i in range(self.alloc, length):
                self.rows[i] = self.entries + (i * self.dim)
            self.alloc = length

    cdef void fill_rays(self, const PPL_Generator_System& gs):
        r"""
        Set the attributes related to the cone
        """
        cdef PPL_gs_iterator * gsi_ptr = new PPL_gs_iterator(gs.begin())
        cdef size_t j
        while gsi_ptr[0] != gs.end():
            if deref(gsi_ptr[0]).is_ray():
                self.fit_length(self.length + 1)
                for j in range(self.dim):
                    mpz_set(self.rows[self.length][j],
                            deref(gsi_ptr[0]).coefficient(PPL_Variable(j)).get_mpz_t())
                self.length += 1
            gsi_ptr[0].inc(1)
        del gsi_ptr
        qsort(self.rows, self.length, sizeof(mpz_t *), & compare)

    def dimension(self):
        return self.poly.affine_dimension()

    def ambient_dimension(self):
        return self.poly.space_dimension()

    def n_rays(self):
        return self.length

    def permutation(self):
        top = [self.perm[i] for i in range(self.dim)]
        bot = [self.perm[self.dim+i] for i in range(self.dim)]
        return [top, bot]

    def polytope(self):
        return self.poly

    def __hash__(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics import *
            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: F = iet.IETFamily(p, Polyhedron(rays=(ZZ**4).basis()))
            sage: hash(F)
            4696861030007372737
        """
        cdef Py_hash_t x, y, mult
        cdef size_t i, j

        mult = 3
        x = 0x345678
        for i in range(self.dim):
            x = (x ^ (self.perm[i])) * mult
            mult += 1345
            x += 35

        for i in range(self.length):
            for j in range(self.dim):
                y = mpz_pythonhash(self.rows[i][j])
                x = (x ^ y) * mult
                mult += 82520
            x += 97531

        return x

    def __richcmp__(IETFamily_pyx self, IETFamily_pyx other, int op):
        r"""
        TESTS::

            sage: from surface_dynamics import *

            sage: p1 = iet.Permutation('a b c d', 'd c b a')
            sage: p2 = iet.Permutation('a b c d', 'd a b c')
            sage: rays1 = [[0,1,0,1],[1,0,1,0],[1,1,0,0],[0,0,1,1]]
            sage: rays2 = [[1,0,0,0],[1,0,1,0],[1,1,0,0],[0,0,1,1]]
            sage: rays3 = [[1,0,1,0],[1,1,0,0],[0,0,1,1]]

            sage: iet.IETFamily(p1, rays1) == iet.IETFamily(p1, rays1)
            True
            sage: iet.IETFamily(p1, rays1) != iet.IETFamily(p1, rays1)
            False

            sage: iet.IETFamily(p1, rays1) == iet.IETFamily(p2, rays1)
            False
            sage: iet.IETFamily(p1, rays1) != iet.IETFamily(p2, rays1)
            True

            sage: iet.IETFamily(p1, rays1) == iet.IETFamily(p2, rays1)
            False
            sage: iet.IETFamily(p1, rays1) != iet.IETFamily(p2, rays1)
            True

            sage: iet.IETFamily(p1, rays1) == iet.IETFamily(p1, rays3)
            False
            sage: iet.IETFamily(p1, rays1) != iet.IETFamily(p1, rays3)
            True

            sage: iet.IETFamily(p1, rays1) == iet.IETFamily(p2, rays3)
            False
            sage: iet.IETFamily(p1, rays1) != iet.IETFamily(p2, rays3)
            True
        """
        if op != Py_EQ and op != Py_NE:
            raise TypeError('only equality and difference tests available')
        if self.dim != other.dim or self.length != other.length:
            return op == Py_NE

        cdef size_t i,j
        for i in range(2 * self.dim):
            if self.perm[i] != other.perm[i]:
                return op == Py_NE
        for i in range(self.length):
            for j in range(self.dim):
                if mpz_cmp(self.rows[i][j], other.rows[i][j]):
                    return op == Py_NE

        return op == Py_EQ

    def has_connection(self):
        from warnings import warn
        warn('has_connection is deprecated; use has_zero_connection instead')
        return self.has_zero_connection()

    def has_zero_connection(self, bint certificate=False):
        r"""
        Check whether this family has a connection of length zero.

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: p = iet.Permutation('A B C D E', 'E D C B A')
            sage: C0 = Polyhedron(rays=(ZZ**5).basis())
            sage: iet.IETFamily(p, C0).has_zero_connection()
            False

            sage: C = Polyhedron(rays=[(1,2,3,2,3),(1,1,0,0,2)])
            sage: iet.IETFamily(p, C).has_zero_connection()
            True

            sage: C = Polyhedron(rays=[(1,1,1,1,4),(2,0,1,2,5)])
            sage: iet.IETFamily(p, C).has_zero_connection()
            True

            sage: C = Polyhedron(rays=[(1,0,0,1,0),(1,0,0,0,1),(0,1,0,1,0),(0,1,0,0,1)])
            sage: iet.IETFamily(p, C).has_zero_connection()
            True

            sage: C = Polyhedron(rays=[(1,2,2,2,2)])
            sage: iet.IETFamily(p, C).has_zero_connection()
            False

            sage: FF = iet.IETFamily(iet.Permutation('a b c d', 'd a c b'), [[2,1,0,3],[4,0,1,4]])
            sage: FF.has_zero_connection()
            True
            sage: FF.has_zero_connection(certificate=True)
            (True, (1, 2))
        """
        cdef bint ans
        ans_certif = None

        # x runs through the top singularities
        # y runs through the bottom singularities
        cdef mpz_t * x = <mpz_t *> check_malloc(2 * self.length * sizeof(mpz_t))
        cdef mpz_t * y = x + self.length
        cdef size_t i, j, k
        for i in range(self.length):
            mpz_init(x[i])
            mpz_set(x[i], self.rows[i][self.perm[0]])
            mpz_init(y[i])
            mpz_set(y[i], self.rows[i][self.perm[self.dim]])

        j = 1
        for i in range(1, self.dim):
            while j < self.dim and mpz_vec_any_smaller(y, x, self.length):
                for k in range(self.length):
                    mpz_add(y[k], y[k], self.rows[k][self.perm[self.dim+j]])
                j += 1
            if mpz_vec_equal(x, y, self.length):
                ans = True
                if certificate:
                    ans_certif = (j, i)
                break
            for k in range(self.length):
                mpz_add(x[k], x[k], self.rows[k][self.perm[i]])
        else:
            ans = False

        for i in range(2 * self.length):
            mpz_clear(x[i])
        sig_free(x)

        return (ans, ans_certif) if certificate else ans

    def children(self):
        r"""
        Compute the children of the pair ``(p, C)`` composed of a permutation and a
        slice of a Rauzy simplex.

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: F = iet.IETFamily(p, Polyhedron(rays=(ZZ**4).basis()))
            sage: top, bot = F.children()
            sage: top[0]
            't'
            sage: top[1]
            Linear iet family of dimension 4 in RR^4
            top a b c d
            bot d a c b
            0 0 0 1
            0 0 1 0
            0 1 0 0
            1 0 0 0
            sage: bot[0]
            'b'
            sage: bot[1]
            Linear iet family of dimension 4 in RR^4
            top a d b c
            bot d c b a
            0 0 0 1
            0 0 1 0
            0 1 0 0
            1 0 0 0

            sage: p = iet.Permutation('a b c d e f', 'e b f c a d')
            sage: F = iet.IETFamily(p, Polyhedron(rays=[(1,1,1,1,1,2),(0,1,0,1,0,1)]))
            sage: c = F.children()
            sage: len(c)
            1
            sage: c[0][0]
            't'
            sage: c[0][1]
            Linear iet family of dimension 2 in RR^6
            top a b c d e f
            bot e b f d c a
            0 1 0 1 0 0
            1 1 1 1 1 1

            sage: F = iet.IETFamily(p, Polyhedron(rays=[(0,0,0,1,0,0),(1,0,0,0,0,0),(1,1,1,1,1,1)]))
            sage: c = F.children()
            sage: len(c)
            1
            sage: c[0][0]
            'b'
            sage: c[0][1]
            Linear iet family of dimension 3 in RR^6
            top a b c d f e
            bot e b f c a d
            0 0 0 1 0 0
            1 0 0 0 0 0
            1 1 1 0 1 1
        """
        cdef int itop = self.perm[self.dim - 1]
        cdef int ibot = self.perm[2*self.dim - 1]
        cdef list ans = []
        cdef size_t d = self.poly.affine_dimension()
        cdef IETFamily_pyx FF

        cdef C_Polyhedron Ctop = C_Polyhedron(self.poly)
        cdef C_Polyhedron ineq = C_Polyhedron(Variable(itop) >= Variable(ibot))
        ineq.add_space_dimensions_and_embed(self.dim - max(ibot, itop) - 1)
        Ctop.intersection_assign(ineq)
        if Ctop.affine_dimension() == d:
            Ctop.affine_image(Variable(itop), Variable(itop) - Variable(ibot))
            FF = self._new()
            FF.dim = self.dim
            FF.perm = <int *> check_malloc(2 * self.dim * sizeof(int))
            memcpy(FF.perm, self.perm, 2 * sizeof(int) * FF.dim)
            rauzy_top(FF.perm, FF.dim)
            FF.fill_rays(Ctop.thisptr.minimized_generators())
            FF.poly = Ctop
            ans.append(('t', FF))

        cdef C_Polyhedron Cbot = C_Polyhedron(self.poly)
        ineq = C_Polyhedron(Variable(itop) <= Variable(ibot))
        ineq.add_space_dimensions_and_embed(self.dim - max(ibot, itop) - 1)
        Cbot.intersection_assign(ineq)
        if Cbot.affine_dimension() == d:
            Cbot.affine_image(Variable(ibot), Variable(ibot) - Variable(itop))
            FF = self._new()
            FF.dim = self.dim
            FF.perm = <int *> check_malloc(2 * self.dim * sizeof(int))
            memcpy(FF.perm, self.perm, 2 * sizeof(int) * FF.dim)
            rauzy_bot(FF.perm, FF.dim)
            FF.fill_rays(Cbot.thisptr.minimized_generators())
            FF.poly = Cbot
            ans.append(('b', FF))

        return ans

    def random_element(self, ring, *args, **kwds):
        r"""
        Return a random point in the cone

        EXAMPLES::

            sage: from surface_dynamics import *
            sage: from surface_dynamics.misc.linalg import deformation_space
            sage: q = iet.Permutation([0,1,2,3,4,5],[5,3,2,1,0,4])
            sage: rays = [[0, 0, 0, 1, 1, 0], [3, 1, 0, 1, 0, 2], [5, 0, 1, 2, 0, 3]]
            sage: F = iet.IETFamily(q, rays)
            sage: F
            Linear iet family of dimension 3 in RR^6
            top 0 1 2 3 4 5
            bot 5 3 2 1 0 4
            0 0 0 1 1 0
            3 1 0 1 0 2
            5 0 1 2 0 3
            sage: x = polygen(ZZ)
            sage: K.<a> = NumberField(x^10 - 3, embedding=AA(3)**(1/10))
            sage: T = F.random_element(K)
            sage: V = deformation_space(T.lengths())
            sage: (QQ**6).subspace(F.rays()) == V
            True
        """
        from .labelled import LabelledPermutationIET as Permutation
        from .iet import IntervalExchangeTransformation as IET

        top = [self.perm[i] for i in range(self.dim)]
        bot = [self.perm[i] for i in range(self.dim,2*self.dim)]
        p = Permutation((top, bot), alphabet=range(self.dim))

        coeffs = []
        for _ in range(self.length):
            c = ring.random_element(*args, **kwds)
            while c.is_zero():
                c = ring.random_element(*args, **kwds)
            coeffs.append(c.abs())
        lengths = [sum(coeffs[i] * self.ray_coefficient(i, j) for i in range(self.length)) for j in range(self.dim)]

        return IET(p, lengths)

    # TODO: find a better name and output!!!
    def random_element_statistics(self, ring, num_exp=100, num_iterations=200, intervalxt=None, *args, **kwds):
        r"""
        Return a triple ``(num_minimals, num_saddles, num_unknowns)`` of number
        of interval exchange transformations found to be respectively minimal, with
        a saddle connection or "unknown" by mean of at most ``num_iterations`` steps
        of Zorich induction over a sample of ``num_exp`` random interval exchange
        transformations in this family.

        INPUT:

        - ``ring`` -- the ring in which is taken the lengths (usually a number field)

        - ``num_exp`` -- number of experiments

        - ``num_iterations`` -- the number of Zorich move to be performed on the
          interval exchange transformations to find a saddle connection

        - ``intervalxt`` -- whether to use pyintervalxt for iteration (much faster for
          large number of induction steps)

        - extra argument are transfered to the method :meth:`random_element`

        EXAMPLES::

            sage: from surface_dynamics import *

            sage: kwds = {'num_bound': 10000, 'den_bound': 1}

            sage: q = iet.Permutation([0,1,2,3,4,5],[5,3,2,1,0,4])
            sage: rays = [[0, 0, 0, 1, 1, 0], [3, 1, 0, 1, 0, 2], [5, 0, 1, 2, 0, 3]]
            sage: F = iet.IETFamily(q, rays)
            sage: x = polygen(QQ)
            sage: K.<cbrt3> = NumberField(x^3 - 3, embedding=AA(3)**(1/3))
            sage: F.random_element_statistics(K, intervalxt=False, **kwds)
            (0, 100, 0)
            sage: F.random_element_statistics(K, intervalxt=True, **kwds) # optional - gmpxxyy pyeantic pyintervalxt
            (0, 100, 0)

            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: F = iet.IETFamily(p, Polyhedron(rays=[(1,2,0,1), (3,1,1,0)]))
            sage: num_minimals, num_saddles, num_unknowns = F.random_element_statistics(K, **kwds)
            sage: num_saddles
            0

        A conjectural counterexample to Dynnikov-Skripshenko conjecture in genus 4 (stratum
        component is H(3,3)^hyp)::

            sage: path = 'bbbbbtbttbttbbttbtttbbbbttbttbtbbtbbtbbtt'
            sage: p = iet.Permutation('a b c d e f g h i', 'i h g f e d c b a')
            sage: R = p.rauzy_diagram()
            sage: T = R.path(p, *path).self_similar_iet()
            sage: T.sah_arnoux_fathi_invariant()
            (0, 0, 0, 0, 0, 0)
            sage: F = iet.IETFamily(T)
            sage: F
            Linear iet family of dimension 4 in RR^9
            top a b c d e f g h i
            bot i h g f e d c b a
            71 67 70 107 0 0 19 2 0
            95 63 21 44 35 0 0 0 43
            105 94 88 147 0 0 0 15 19
            253 310 0 243 34 101 68 0 0
            425 395 404 611 14 0 129 0 0
            524 703 0 628 0 245 51 68 0
            619 462 0 229 238 63 0 0 272
            701 910 0 823 0 308 0 119 51
            sage: F.random_element_statistics(K, num_exp=20, num_iterations=256, intervalxt=False, **kwds)
            (0, 0, 20)
        """
        if intervalxt is None:
            try:
                import gmpxxyy
                import pyintervalxt
                import pyeantic
            except ImportError:
                intervalxt = False
            else:
                intervalxt = True

        num_minimals = 0
        num_saddles = 0
        num_unknowns = 0
        bads = []

        if intervalxt:
            from pyintervalxt import intervalxt
            Result = intervalxt.InductionStep.Result
            LIMIT_REACHED = [Result.LIMIT_REACHED]
            SADDLE = [Result.CYLINDER, Result.SEPARATING_CONNECTION, Result.NON_SEPARATING_CONNECTION]
            MINIMAL = [Result.WITHOUT_PERIODIC_TRAJECTORY_BOSHERNITZAN, Result.WITHOUT_PERIODIC_TRAJECTORY_AUTO_SIMILAR]

            from .conversion import iet_to_pyintervalxt
            for i in range(num_exp):
                T = self.random_element(ring, *args, **kwds)
                U = iet_to_pyintervalxt(T)
                res = U.induce(int(num_iterations))
                if res.result in LIMIT_REACHED:
                    num_unknowns += 1
                elif res.result in SADDLE:
                    num_saddles += 1
                elif res.result in MINIMAL:
                    num_minimals += 1
                else:
                    raise RuntimeError("unknown return code from pyintervalxt {}".format(res.result))
        else:
            for i in range(num_exp):
                T = self.random_element(ring, *args, **kwds)
                try:
                    T.zorich_move(iterations=num_iterations)
                except ValueError:
                    num_saddles += 1
                else:
                    num_unknowns += 1

        return (num_minimals, num_saddles, num_unknowns)

    def random_integer_point_statistics(self, num_exp=100):
        r"""
        OUTPUT: a dictionary whose keys are integers ``i`` and values are the
        number of integer points found so that the number of cylinders is
        ``i``.
        """
        raise NotImplementedError
