# distutils: language=c++
# distutils: libraries=ppl gmp
r"""
Linear families of interval exchange transformations
"""

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

from cpython.object cimport Py_EQ, Py_NE

from cysignals.memory cimport sig_malloc, sig_calloc, sig_realloc, sig_free, check_malloc, check_calloc, check_realloc
from libc.stdlib cimport malloc, calloc, realloc, free, qsort
from libc.string cimport memcpy

from ppl.ppl_decl cimport PPL_Variable, PPL_Generator, PPL_Generator_System, PPL_gs_iterator
from ppl.linear_algebra cimport Variable
from ppl.generator cimport Generator_System
from ppl.polyhedron cimport C_Polyhedron

from gmpy2 cimport import_gmpy2, MPZ_Object, MPZ, GMPy_MPZ_New

from cython.operator cimport dereference as deref

import_gmpy2()

########################
# Python imports
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

cdef class IETFamily(object):
    r"""
    A linear family of interval exchange transformations

    This class consists of two attributes:

    - a permutation (must be labeled)

    - a cone of lengths (a subcone of the positive cone)

    EXAMPLES::

        sage: from surface_dynamics.all import *

        sage: p = iet.Permutation('a b c d', 'd c b a')
        sage: F = iet.IETFamily(p, (ZZ**4).basis())
        sage: F
        top 0 1 2 3
        bot 3 2 1 0
        0 0 0 1
        0 0 1 0
        0 1 0 0
        1 0 0 0

        sage: iet.IETFamily(p, [(2,3,0,0), (0,1,1,1)])
        top 0 1 2 3
        bot 3 2 1 0
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
    cdef mpz_t ** rows      # ray pointers to entries (begining of vectors)
    cdef size_t alloc       # ray allocation size
    cdef size_t length      # number of rays
    cdef C_Polyhedron poly  # the polyhedron (that we should get rid of)

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

    def __init__(self, p, C_Polyhedron C):
        # convert and check input
        from surface_dynamics.interval_exchanges.labelled import LabelledPermutationIET
        if not isinstance(p, LabelledPermutationIET):
            raise ValueError('p must be a labelled permutation')

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

    def ray_coefficient(self, size_t nray, size_t i):
        r"""
        Return the ``i``-th coefficient of the ``nray``-th ray

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: F = iet.IETFamily(p, [(2,3,0,0), (0,1,1,1)])
            sage: F.ray_coefficient(0, 2)
            mpz(1)
            sage: F.ray_coefficient(1, 0)
            mpz(2)
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
        cdef MPZ_Object * z = GMPy_MPZ_New(NULL)
        mpz_set(MPZ(z), self.rows[nray][i])
        return <object> z

    def __repr__(self):
        r"""
        TESTS::

            sage: from surface_dynamics.all import *
            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: F = iet.IETFamily(p, [(2,3,0,0), (0,1,1,1)])
            sage: repr(F)  # indirect doctest
            'top 0 1 2 3\nbot 3 2 1 0\n0 1 1 1\n2 3 0 0'
        """
        cdef size_t i, j
        s = []
        s.append('top ' + ' '.join(str(self.perm[i]) for i in range(self.dim)))
        s.append('bot ' + ' '.join(str(self.perm[self.dim+i]) for i in range(self.dim)))
        for i in range(self.length):
            s.append(' '.join(str(self.ray_coefficient(i, j)) for j in range(self.dim)))
        return '\n'.join(s)

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

    def permutation(self):
        r"""
        Return the permutation as a list of integers
        """
        cdef size_t i
        return [[self.perm[i] for i in range(self.dim)], [self.perm[self.dim+i] for i in range(self.dim)]]

    def __hash__(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.all import *
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

    def __richcmp__(IETFamily self, IETFamily other, int op):
        r"""
        TESTS::

            sage: from surface_dynamics.all import *

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
        r"""
        Check whether this family has a connection.

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p = iet.Permutation('A B C D E', 'E D C B A')
            sage: C0 = Polyhedron(rays=(ZZ**5).basis())
            sage: iet.IETFamily(p, C0).has_connection()
            False

            sage: C = Polyhedron(rays=[(1,2,3,2,3),(1,1,0,0,2)])
            sage: iet.IETFamily(p, C).has_connection()
            True

            sage: C = Polyhedron(rays=[(1,1,1,1,4),(2,0,1,2,5)])
            sage: iet.IETFamily(p, C).has_connection()
            True

            sage: C = Polyhedron(rays=[(1,0,0,1,0),(1,0,0,0,1),(0,1,0,1,0),(0,1,0,0,1)])
            sage: iet.IETFamily(p, C).has_connection()
            True

            sage: C = Polyhedron(rays=[(1,2,2,2,2)])
            sage: iet.IETFamily(p, C).has_connection()
            False
        """
        cdef bint ans

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
                break
            for k in range(self.length):
                mpz_add(x[k], x[k], self.rows[k][self.perm[i]])
        else:
            ans = False

        for i in range(2 * self.length):
            mpz_clear(x[i])
        sig_free(x)

        return ans

    def children(self):
        r"""
        Compute the children of the pair ``(p, C)`` composed of a permutation and a
        slice of a Rauzy simplex.

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: F = iet.IETFamily(p, Polyhedron(rays=(ZZ**4).basis()))
            sage: top, bot = F.children()
            sage: top[0]
            't'
            sage: top[1]
            top 0 1 2 3
            bot 3 0 2 1
            0 0 0 1
            0 0 1 0
            0 1 0 0
            1 0 0 0
            sage: bot[0]
            'b'
            sage: bot[1]
            top 0 3 1 2
            bot 3 2 1 0
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
            top 0 1 2 3 4 5
            bot 4 1 5 3 2 0
            0 1 0 1 0 0
            1 1 1 1 1 1

            sage: F = iet.IETFamily(p, Polyhedron(rays=[(0,0,0,1,0,0),(1,0,0,0,0,0),(1,1,1,1,1,1)]))
            sage: c = F.children()
            sage: len(c)
            1
            sage: c[0][0]
            'b'
            sage: c[0][1]
            top 0 1 2 3 5 4
            bot 4 1 5 2 0 3
            0 0 0 1 0 0
            1 0 0 0 0 0
            1 1 1 0 1 1
        """
        cdef int itop = self.perm[self.dim - 1]
        cdef int ibot = self.perm[2*self.dim - 1]
        cdef list ans = []
        cdef size_t d = self.poly.affine_dimension()
        cdef IETFamily FF

        cdef C_Polyhedron Ctop = C_Polyhedron(self.poly)
        cdef C_Polyhedron ineq = C_Polyhedron(Variable(itop) >= Variable(ibot))
        ineq.add_space_dimensions_and_embed(self.dim - max(ibot, itop) - 1)
        Ctop.intersection_assign(ineq)
        if Ctop.affine_dimension() == d:
            Ctop.affine_image(Variable(itop), Variable(itop) - Variable(ibot))
            FF = IETFamily.__new__(IETFamily)
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
            FF = IETFamily.__new__(IETFamily)
            FF.dim = self.dim
            FF.perm = <int *> check_malloc(2 * self.dim * sizeof(int))
            memcpy(FF.perm, self.perm, 2 * sizeof(int) * FF.dim)
            rauzy_bot(FF.perm, FF.dim)
            FF.fill_rays(Cbot.thisptr.minimized_generators())
            FF.poly = Cbot
            ans.append(('b', FF))

        return ans

    def tree(self, max_depth=5, verbose=False):
        r"""
        The tree should be smarter in several ways:

        - once a saddle connection is found it should check whether it is stable and
          whether it corresponds to a splitting, a cylinder decomposition, etc

        - implement a C version!

        - understand what the parabolic are doing!

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p = iet.Permutation([0,1,2,3,4,5],[5,4,3,2,1,0])
            sage: rays = [[5, 1, 0, 0, 3, 8], [2, 1, 0, 3, 0, 5], [1, 0, 1, 2, 0, 3], [3, 0, 1, 0, 2, 5]]
            sage: F = iet.IETFamily(p, rays)
        """
        s = ''
        f = self
        branch = [[(s, f)]]
        seen = set([f])
        while True:
            if verbose:
                print("branch:")
                for ss,ff in branch[-1]:
                    print(ss)
                    print(ff)
                    print()
            while len(branch) < max_depth:
                branch.append([])
                for ss, ff in f.children():
                    if verbose:
                        print("looking children", ss)
                        print("cone", ff)
                    # check for saddles among the children
                    if ff.has_connection():
                        yield 'saddle', s+ss, ff

                    # check for auto-simlarity
                    elif ff in seen:
                        yield 'autosim', s+ss, ff

                    else:
                        branch[-1].append((s+ss, ff))

                if not branch[-1]:
                    branch.pop(-1)
                    break
                else:
                    s, f = branch[-1][-1]
                    assert f not in seen
                    seen.add(f)

            # backtrack
            while branch and len(branch[-1]) == 1:
                s, f = branch.pop()[0]
                if f not in seen:
                    raise RuntimeError("s = {}\nf =\n{}".format(s, f))
                seen.remove(f)

            if not branch:
                return

            s, f = branch[-1].pop(-1)
            seen.remove(f)
            s, f = branch[-1][-1]
            seen.add(f)

