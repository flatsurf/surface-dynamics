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

#from gmpy2 cimport import_gmpy2, MPZ_Object, MPZ, GMPy_MPZ_New
#import_gmpy2()


########################
# Begining of the code
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
        sage: F = iet.IETFamily(p, (ZZ**4).basis())    # optional - pplpy
        sage: F  # optional - pplpy
        top 0 1 2 3
        bot 3 2 1 0
        0 0 0 1
        0 0 1 0
        0 1 0 0
        1 0 0 0

        sage: iet.IETFamily(p, [(2,3,0,0), (0,1,1,1)]) # optional - pplpy
        top 0 1 2 3
        bot 3 2 1 0
        0 1 1 1
        2 3 0 0

        sage: iet.IETFamily(p, [(1,0,0,0), (1,-1,1,1)]) # optional - pplpy
        Traceback (most recent call last):
        ...
        ValueError: C must be a subcone of the non-negative cone
        sage: iet.IETFamily(p, Polyhedron(vertices=[(0,0,0,0), (1,0,0,0),(0,1,0,0),(1,1,1,1)])) # optional - pplpy
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
            sage: F = iet.IETFamily(p, [(2,3,0,0), (0,1,1,1)]) # optional - pplpy
            sage: F.ray_coefficient(0, 2) # optional - pplpy
            1
            sage: F.ray_coefficient(1, 0) # optional - pplpy
            2
            sage: F.ray_coefficient(3, 0) # optional - pplpy
            Traceback (most recent call last):
            ...
            IndexError: row index out of range
            sage: F.ray_coefficient(1, 5) # optional - pplpy
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

    def rays(self):
        r"""
        Return the rays as vectors

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: F = iet.IETFamily(p, [(2,3,0,0), (0,1,1,1)]) # optional - pplpy
            sage: F.rays() # optional - pplpy
            [(0, 1, 1, 1), (2, 3, 0, 0)]
        """
        from sage.modules.free_module import FreeModule
        from sage.rings.rational_field import QQ
        F = FreeModule(QQ, self.dim)
        return [F([self.ray_coefficient(i, j) for j in range(self.dim)]) for i in range(self.length)]

    def __repr__(self):
        r"""
        TESTS::

            sage: from surface_dynamics.all import *
            sage: p = iet.Permutation('a b c d', 'd c b a')
            sage: F = iet.IETFamily(p, [(2,3,0,0), (0,1,1,1)]) # optional - pplpy
            sage: repr(F)  # indirect doctest # optional - pplpy
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
            sage: F = iet.IETFamily(p, Polyhedron(rays=(ZZ**4).basis())) # optional - pplpy
            sage: hash(F) # optional - pplpy
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

            sage: iet.IETFamily(p1, rays1) == iet.IETFamily(p1, rays1) # optional - pplpy
            True
            sage: iet.IETFamily(p1, rays1) != iet.IETFamily(p1, rays1) # optional - pplpy
            False

            sage: iet.IETFamily(p1, rays1) == iet.IETFamily(p2, rays1) # optional - pplpy
            False
            sage: iet.IETFamily(p1, rays1) != iet.IETFamily(p2, rays1) # optional - pplpy
            True

            sage: iet.IETFamily(p1, rays1) == iet.IETFamily(p2, rays1) # optional - pplpy
            False
            sage: iet.IETFamily(p1, rays1) != iet.IETFamily(p2, rays1) # optional - pplpy
            True

            sage: iet.IETFamily(p1, rays1) == iet.IETFamily(p1, rays3) # optional - pplpy
            False
            sage: iet.IETFamily(p1, rays1) != iet.IETFamily(p1, rays3) # optional - pplpy
            True

            sage: iet.IETFamily(p1, rays1) == iet.IETFamily(p2, rays3) # optional - pplpy
            False
            sage: iet.IETFamily(p1, rays1) != iet.IETFamily(p2, rays3) # optional - pplpy
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

    # this is completely useless and checks many times the same thing
    # we should just look at left/right Rauzy induction
    def has_connection(self):
        r"""
        Check whether this family has a connection.

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p = iet.Permutation('A B C D E', 'E D C B A')
            sage: C0 = Polyhedron(rays=(ZZ**5).basis())
            sage: iet.IETFamily(p, C0).has_connection() # optional - pplpy
            False

            sage: C = Polyhedron(rays=[(1,2,3,2,3),(1,1,0,0,2)])
            sage: iet.IETFamily(p, C).has_connection() # optional - pplpy
            True

            sage: C = Polyhedron(rays=[(1,1,1,1,4),(2,0,1,2,5)])
            sage: iet.IETFamily(p, C).has_connection() # optional - pplpy
            True

            sage: C = Polyhedron(rays=[(1,0,0,1,0),(1,0,0,0,1),(0,1,0,1,0),(0,1,0,0,1)])
            sage: iet.IETFamily(p, C).has_connection() # optional - pplpy
            True

            sage: C = Polyhedron(rays=[(1,2,2,2,2)])
            sage: iet.IETFamily(p, C).has_connection() # optional - pplpy
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
            sage: F = iet.IETFamily(p, Polyhedron(rays=(ZZ**4).basis())) # optional - pplpy
            sage: top, bot = F.children() # optional - pplpy
            sage: top[0] # optional - pplpy
            't'
            sage: top[1] # optional - pplpy
            top 0 1 2 3
            bot 3 0 2 1
            0 0 0 1
            0 0 1 0
            0 1 0 0
            1 0 0 0
            sage: bot[0] # optional - pplpy
            'b'
            sage: bot[1] # optional - pplpy
            top 0 3 1 2
            bot 3 2 1 0
            0 0 0 1
            0 0 1 0
            0 1 0 0
            1 0 0 0

            sage: p = iet.Permutation('a b c d e f', 'e b f c a d')
            sage: F = iet.IETFamily(p, Polyhedron(rays=[(1,1,1,1,1,2),(0,1,0,1,0,1)])) # optional - pplpy
            sage: c = F.children() # optional - pplpy
            sage: len(c) # optional - pplpy
            1
            sage: c[0][0] # optional - pplpy
            't'
            sage: c[0][1] # optional - pplpy
            top 0 1 2 3 4 5
            bot 4 1 5 3 2 0
            0 1 0 1 0 0
            1 1 1 1 1 1

            sage: F = iet.IETFamily(p, Polyhedron(rays=[(0,0,0,1,0,0),(1,0,0,0,0,0),(1,1,1,1,1,1)])) # optional - pplpy
            sage: c = F.children() # optional - pplpy
            sage: len(c) # optional - pplpy
            1
            sage: c[0][0] # optional - pplpy
            'b'
            sage: c[0][1] # optional - pplpy
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
            sage: F = iet.IETFamily(p, rays) # optional - pplpy
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

    def random_element(self, ring, *args, **kwds):
        r"""
        Return a random point in the cone

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: from surface_dynamics.misc.linalg import deformation_space
            sage: q = iet.Permutation([0,1,2,3,4,5],[5,3,2,1,0,4])
            sage: rays = [[0, 0, 0, 1, 1, 0], [3, 1, 0, 1, 0, 2], [5, 0, 1, 2, 0, 3]]
            sage: F = iet.IETFamily(q, rays) # optional - pplpy
            sage: F # optional - pplpy
            top 0 1 2 3 4 5
            bot 5 3 2 1 0 4
            0 0 0 1 1 0
            3 1 0 1 0 2
            5 0 1 2 0 3
            sage: x = polygen(ZZ)
            sage: K.<a> = NumberField(x^10 - 3, embedding=AA(3)**(1/10))
            sage: T = F.random_element(K) # optional - pplpy
            sage: V = deformation_space(T.lengths()) # optional - pplpy
            sage: (QQ**6).subspace(F.rays()) == V # optional - pplpy
            True
        """
        from constructors import Permutation, IET

        top = [self.perm[i] for i in range(self.dim)]
        bot = [self.perm[i] for i in range(self.dim,2*self.dim)]
        p = Permutation(top, bot, alphabet=range(self.dim))

        coeffs = []
        for _ in range(self.length):
            c = ring.random_element(*args, **kwds)
            while c.is_zero():
                c = ring.random_element(*args, **kwds)
            coeffs.append(c.abs())
        lengths = [sum(coeffs[i] * self.ray_coefficient(i, j) for i in range(self.length)) for j in range(self.dim)]

        return IET(p, lengths)

    # TODO: find a better name and output!!!
    def random_element_statistics(self, ring, num_exp=100, num_iterations=200, *args, **kwds):
        r"""
        INPUT:

        - ``ring`` -- the ring in which is taken the lengths (usually a number field)

        - ``num_exp`` -- number of experiments

        - ``num_iterations`` -- the number of Zorich move to be performed on the
          interval exchange transformations to find a saddle connection

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: q = iet.Permutation([0,1,2,3,4,5],[5,3,2,1,0,4])
            sage: rays = [[0, 0, 0, 1, 1, 0], [3, 1, 0, 1, 0, 2], [5, 0, 1, 2, 0, 3]]
            sage: F = iet.IETFamily(q, rays)  # optional - pplpy
            sage: F.random_element_statistics(K)  # optional - pplpy
            100 saddles found on 100 random point

        A conjectural counterexample to Dynnikov-Skripshenko conjecture::

            sage: path = 'bbbbbtbttbttbbttbtttbbbbttbttbtbbtbbtbbtt'
            sage: p = iet.Permutation('a b c d e f g h i', 'i h g f e d c b a')
            sage: R = p.rauzy_diagram()
            sage: T = R.path(p, *path).self_similar_iet()
            sage: T.sah_arnoux_fathi_invariant()
            (0, 0, 0, 0, 0, 0)
            sage: F = iet.IETFamily(T)  # optional - pplpy
            sage: F  # optional - pplpy
            top 0 1 2 3 4 5 6 7 8
            bot 8 7 6 5 4 3 2 1 0
            71 67 70 107 0 0 19 2 0
            95 63 21 44 35 0 0 0 43
            105 94 88 147 0 0 0 15 19
            253 310 0 243 34 101 68 0 0
            425 395 404 611 14 0 129 0 0
            524 703 0 628 0 245 51 68 0
            619 462 0 229 238 63 0 0 272
            701 910 0 823 0 308 0 119 51
            sage: K.<cbrt3> = NumberField(x^3 - 3, embedding=AA(3)**(1/3))
            sage: F.random_element_statistics(K, num_exp=20, num_iterations=200)  # optional - pplpy # random
            ... saddles found on 20 random points
        """
        num_saddles = 0
        bads = []
        for i in range(num_exp):
            T = self.random_element(ring, *args, **kwds)
            try:
                T.zorich_move(iterations=num_iterations)
            except ValueError:
                num_saddles += 1
            else:
                bads.append(T)
        print('%d saddles found on %d random points' % (num_saddles, num_exp))

    def random_integer_point_statistics(self, num_exp=100):
        r"""
        OUTPUT: a dictionary whose keys are integers ``i`` and values are the
        number of integer points found so that the number of cylinders is
        ``i``.
        """
        # TODO: make iet_interface part of flatsurf
        import iet_interface
        from sage.rings.integer_ring import ZZ
        top, bot = self.permutation()
        dim = len(top)
        cdef dict cyls = {}
        for _ in range(num_exp):
            lengths = sum(ZZ.random_element(0,100000000) * r for r in self.rays())
            ncyls = iet_interface.number_of_cylinders(top, bot, lengths)
            if ncyls not in cyls:
                cyls[ncyls] = 0
            cyls[ncyls] += 1
        return cyls

