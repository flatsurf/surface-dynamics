from __future__ import print_function

from cysignals.memory cimport check_malloc, sig_free
from libc.string cimport memcpy

def _test_iterator1(int n, int k):
    r"""
    TESTS::

        sage: from surface_dynamics.interval_exchanges.integer_iet import _test_iterator1

        sage: list(_test_iterator1(0, 0))
        [[]]
        sage: list(_test_iterator1(0, 1))
        []
        sage: list(_test_iterator1(1, 0))
        []

        sage: list(_test_iterator1(1, 1))
        [[1L]]
        sage: list(_test_iterator1(2, 1))
        [[2L]]

        sage: list(_test_iterator1(2, 2))
        [[1L, 1L]]
        sage: list(_test_iterator1(3, 2))
        [[2L, 1L], [1L, 2L]]
        sage: list(_test_iterator1(4, 2))
        [[3L, 1L], [2L, 2L], [1L, 3L]]

        sage: list(_test_iterator1(6, 3))
        [[4L, 1L, 1L],
         [3L, 2L, 1L],
         [2L, 3L, 1L],
         [1L, 4L, 1L],
         [3L, 1L, 2L],
         [2L, 2L, 2L],
         [1L, 3L, 2L],
         [2L, 1L, 3L],
         [1L, 2L, 3L],
         [1L, 1L, 4L]]

        sage: list(_test_iterator1(3, 4))
        []
    """
    cdef uint64_t * x

    x = <uint64_t *> check_malloc(k * sizeof(uint64_t))

    if not int_vector_first(x, n, k):
        sig_free(x)
        return
    yield [x[i] for i in range(k)]
    while int_vector_next(x, n, k):
        yield [x[i] for i in range(k)]

    sig_free(x)

def _test_iterator2(int n, int kfree, int ktop, int kbot):
    r"""
    TESTS::

        sage: from surface_dynamics.interval_exchanges.integer_iet import _test_iterator2

        sage: list(_test_iterator2(3, 1, 1, 1))
        [[1L, 1L, 1L]]
        sage: list(_test_iterator2(4, 1, 1, 1))
        [[2L, 1L, 1L]]
        sage: list(_test_iterator2(5, 1, 1, 1))
        [[3L, 1L, 1L], [1L, 2L, 2L]]

        sage: list(_test_iterator2(3, 0, 1, 1))
        []
        sage: list(_test_iterator2(4, 0, 2, 1))
        [[1L, 1L, 2L]]
        sage: list(_test_iterator2(1, 1, 0, 0))
        [[1L]]

        sage: list(_test_iterator2(5, 2, 1, 1))
        [[2L, 1L, 1L, 1L], [1L, 2L, 1L, 1L]]
        sage: list(_test_iterator2(7, 2, 2, 1))
        [[2L, 1L, 1L, 1L, 2L], [1L, 2L, 1L, 1L, 2L]]
        sage: list(_test_iterator2(8, 2, 2, 1))
        [[3L, 1L, 1L, 1L, 2L],
         [2L, 2L, 1L, 1L, 2L],
         [1L, 3L, 1L, 1L, 2L],
         [1L, 1L, 2L, 1L, 3L],
         [1L, 1L, 1L, 2L, 3L]]
        sage: list(_test_iterator2(8, 2, 1, 2))
        [[3L, 1L, 2L, 1L, 1L],
         [2L, 2L, 2L, 1L, 1L],
         [1L, 3L, 2L, 1L, 1L],
         [1L, 1L, 3L, 2L, 1L],
         [1L, 1L, 3L, 1L, 2L]]

         sage: for a in _test_iterator2(5, 1, 1, 1):
         ....:     assert a[0] + 2*a[1] == a[0] + 2*a[2] == 5
    """
    cdef li_vector_iterator_t v

    int_li_vector_init(v, n, kfree, ktop, kbot)

    int_li_vector_prefirst(v)
    while int_li_vector_first_or_next(v):
        yield [(v.x)[i] for i in range(kfree+ktop+kbot)]

    int_li_vector_clear(v)

def perm_to_twin(p):
    r"""
    Return the twin associated to the permutation ``p``.

    EXAMPLES::

        sage: from surface_dynamics.interval_exchanges.integer_iet import perm_to_twin
        sage: perm_to_twin([0,0,1,1,2,2])
        [1, 0, 3, 2, 5, 4]
        sage: perm_to_twin([0,1,2,0])
        Traceback (most recent call last):
        ...
        ValueError: p must contain integer only between 0 and 1
        sage: perm_to_twin([0,0,1,0])
        Traceback (most recent call last):
        ...
        ValueError: each label must appear exactly twice
        sage: perm_to_twin([1,0,1])
        Traceback (most recent call last):
        ...
        ValueError: p must have even length
    """
    n = len(p)
    if n%2 == 1:
        raise ValueError("p must have even length")
    n = n//2
    seen = [None] * len(p)
    t = [None] * len(p)
    for i,j in enumerate(p):
        j = int(j)
        if j < 0 or j >= n:
            raise ValueError("p must contain integer only between 0 and {}".format(n-1))

        ii = seen[j]
        if ii is not None:
            if t[i] is not None or t[ii] is not None:
                raise ValueError("each label must appear exactly twice")
            t[i] = ii
            t[ii] = i
        else:
            seen[j] = i

    return t

cdef set_int_iet(int_iet_t t, top, bot, lengths):
    "Initialize and set t with the data from top, bot, lengths"
    cdef int * labels
    cdef int * twin
    cdef uint64_t * clengths
    cdef int k = len(top)
    cdef int n = (len(top)+len(bot))/2
    cdef int i,j

    # checks and twin construction
    lengths = map(int, lengths)
    if len(lengths) != n:
        raise ValueError('invalid lengths')
    p = map(int, top + bot)
    if len(p) != 2*n:
        raise ValueError('invalid top, bot')
    python_twin = perm_to_twin(top + bot)

    assert all(l > 0 for l in lengths)
    if sum(lengths[i] for i in top) != sum(lengths[i] for i in bot):
        raise ValueError('different lengths for top and bot')


    # fill C data
    labels = <int *> check_malloc(2*n*sizeof(int));
    twin = <int *> check_malloc(2*n*sizeof(int));
    clengths = <uint64_t *> check_malloc(n*sizeof(uint64_t));

    for i in range(2*n):
        labels[i] = p[i]
        twin[i] = python_twin[i]
    for j in range(n):
        clengths[j] = lengths[j]

    int_iet_init(t, n)
    int_iet_set_labels_and_twin(t, labels, twin, k)
    int_iet_set_lengths(t, clengths)
    sig_free(twin)
    sig_free(labels)
    sig_free(clengths)

    if int_iet_check(t):
        int_iet_clear(t)
        raise RuntimeError("invalid iet")

# we want an iterator over lengths

def _relabel(top, bot):
    r"""
    TESTS::

        sage: from surface_dynamics.interval_exchanges.integer_iet import _relabel
        sage: _relabel([1,0,1],[2,2,0])
        ([1, 0, 1], [2, 2, 0], 1, 1, 1)

        sage: _relabel([1,1,2,2], [0,0])
        ([0, 0, 1, 1], [2, 2], 0, 2, 1)

        sage: _relabel([3,1,3,0,0,4],[4,2,2,1])
        ([2, 0, 2, 3, 3, 1], [1, 4, 4, 0], 2, 2, 1)

        sage: _relabel([3, 2, 0, 1], [1, 2, 0, 3])
        ([0, 1, 2, 3], [3, 1, 2, 0], 4, 0, 0)
    """
    n = (len(top) + len(bot)) / 2
    k = len(top)
    j = 0
    p = top + bot
    twin = perm_to_twin(top + bot)
    for i in range(k):
        if twin[i] >= k:
            p[i] = p[twin[i]] = j
            j += 1
    nfree = j
    for i in range(k):
        if i < twin[i] < k:
            p[i] = p[twin[i]] = j
            j += 1
    ntop = j - nfree
    for i in range(k, 2*n):
        if twin[i] >= k and twin[i] > i:
            p[i] = p[twin[i]] = j
            j += 1
    nbot = j - ntop - nfree
    return p[:k], p[k:], nfree, ntop, nbot

# there are several interesting statistics one might want to compute
# such as
#  - number of cylinders
#  - sum of heights of cylinders
#  - mean ratio in view of Siegel Veech constants (not available for now) as
#    we would need the circumferences of cylinders as well
def cylinder_number_statistics(top, bot, uint64_t L):
    r"""
    Return the statistics of the number of cylinders for a given total length

    INPUT:

    - ``top`` -- (list) composition of the top of the cylinder

    - ``bot`` -- (list) composition of the bottom of the cylinder
    
    - ``L`` -- (positive integer) the length we are considering

    EXAMPLES::

        sage: from surface_dynamics.interval_exchanges.integer_iet import cylinder_number_statistics

    The case of H(0,0)::

        sage: top = [0, 1]
        sage: bot = [1, 0]
        sage: for n in range(2,20):
        ....:     s = cylinder_number_statistics(top, bot, n)
        ....:     print "%2d : %3d %3d %3d" %(n, s[1], s[2], s[1] + s[2])
         2 :   1   1   2
         3 :   4   2   6
         4 :   7   5  12
         5 :  16   4  20
         6 :  15  15  30
         7 :  36   6  42
         8 :  35  21  56
         9 :  52  20  72
        10 :  53  37  90
        11 : 100  10 110
        12 :  65  67 132
        13 : 144  12 156
        14 : 115  67 182
        15 : 132  78 210
        16 : 155  85 240
        17 : 256  16 272
        18 : 165 141 306
        19 : 324  18 342

    Which are the same numbers as Q(0,-1^4)::

        sage: top1 = [0, 1]
        sage: bot1 = [1, 0]
        sage: top2 = [0, 0]
        sage: bot2 = [1, 1, 2, 2]
        sage: for _ in range(20):
        ....:     n = randint(2, 100)
        ....:     s1 = cylinder_number_statistics(top1, bot1, n)
        ....:     s2 = cylinder_number_statistics(top2, bot2, 2*n)
        ....:     assert sorted(s1.keys()) == sorted(s2.keys())
        ....:     assert all(2*s1[i] == s2[i] for i in s1)

    Q(1, -1^5)

        sage: top = [1,1]
        sage: bot = [2,2,3,3,4,4]


    """
    cdef int * labels1
    cdef int * labels2
    cdef int * twin1
    cdef int * twin2
    cdef uint64_t * clengths
    cdef int k1 = len(top)
    cdef int k2 = k1 + 1
    cdef int n = (len(top)+len(bot))/2
    cdef int i,j
    cdef uint64_t twist
    cdef int_iet_t t1, t2
    cdef li_vector_iterator_t v

    from collections import defaultdict
    statistics = defaultdict(int)

    top, bot, kfree, ktop, kbot = _relabel(top, bot)

    p1 = top + bot
    p2 = [n] + p1 + [n]
    python_twin1 = perm_to_twin(p1)
    python_twin2 = perm_to_twin(p2)

    labels1 = <int *> check_malloc(2*n*sizeof(int))
    twin1 = <int *> check_malloc(2*n*sizeof(int))
    labels2 = <int *> check_malloc(2*(n+1)*sizeof(int))
    twin2 = <int *> check_malloc(2*(n+1)*sizeof(int))

    clengths = <uint64_t *> check_malloc((n+1)*sizeof(uint64_t));

    for i in range(2*n):
        labels1[i] = p1[i]
        twin1[i] = python_twin1[i]
    int_iet_init(t1, n)

    for i in range(2*(n+1)):
        labels2[i] = p2[i]
        twin2[i] = python_twin2[i]
    int_iet_init(t2, n+1)

    int_li_vector_init(v, L, kfree, ktop, kbot)
    int_li_vector_prefirst(v)
    while int_li_vector_first_or_next(v):
        memcpy(clengths, v.x, n * sizeof(uint64_t))

        # twist 0 case
        int_iet_set_labels_and_twin(t1, labels1, twin1, k1)
        int_iet_set_lengths(t1, clengths)
        statistics[int_iet_num_cylinders(NULL, t1)] += 1

        # positive twists case
        for twist in range(1, L):
            clengths[n] = twist
            int_iet_set_labels_and_twin(t2, labels2, twin2, k2)
            int_iet_set_lengths(t2, clengths)
            statistics[int_iet_num_cylinders(NULL, t2)] += 1

    sig_free(labels1)
    sig_free(twin1)
    sig_free(labels2)
    sig_free(twin2)
    sig_free(clengths)
    int_iet_clear(t1)
    int_iet_clear(t2)
    int_li_vector_clear(v)

    return dict(statistics)

def cylinder_number(top, bot, lengths):
    r"""
    Return the number of cylinders and of a given interval exchange

    EXAMPLES::

        sage: from surface_dynamics.interval_exchanges.integer_iet import cylinder_number

        sage: top = [0, 1]
        sage: bot = [1, 0]
        sage: cylinder_number(top, bot, [13, 24])
        1

        sage: top = [0, 0]
        sage: bot = [1, 1, 2, 2, 3, 3, 4, 4]
        sage: cylinder_number(top, bot, [10, 4, 3, 2, 1])
        2

        sage: top = [0, 0, 1, 1, 2, 2]
        sage: bot = [3, 3, 4, 4, 5, 5, 6, 6]
        sage: cylinder_number(top, bot, [8, 7, 10, 8, 10, 1, 6])
        3
    """
    cdef int_iet_t t
    cdef uint64_t * widths
    set_int_iet(t, top, bot, lengths)
    cdef int res = int_iet_num_cylinders(NULL, t)
    int_iet_clear(t)
    return res

def cylinder_widths(top, bot, lengths):
    r"""
    Return the widths of cylinders of the interval exchange determined by ``top``, ``bot``, ``lengths``.

    EXAMPLES::

        sage: from surface_dynamics.interval_exchanges.integer_iet import cylinder_widths

        sage: top = [0, 1]
        sage: bot = [1, 0]
        sage: cylinder_widths(top, bot, [13, 24])
        [1L]

        sage: top = [0, 0]
        sage: bot = [1, 1, 2, 2, 3, 3, 4, 4]
        sage: cylinder_widths(top, bot, [10, 4, 3, 2, 1])
        [1L, 2L]

    Genus 0 examples::

        sage: top = [0, 0, 1, 1, 2, 2]
        sage: bot = [3, 3, 4, 4, 5, 5, 6, 6]
        sage: cylinder_widths(top, bot, [8, 7, 10, 8, 10, 1, 6])
        [1L, 2L, 8L]
        sage: cylinder_widths(top, bot, [1298, 814, 23, 140, 1034, 819, 142])
        [2L, 3L]
        sage: cylinder_widths(top, bot, [78745, 39773, 19984, 21665, 31239, 437, 85161])
        [1L, 1L, 1L]
        sage: cylinder_widths(top, bot, [57453, 9483, 110211, 100249, 36644, 9853, 30401])
        [1L, 1L, 3L]

        sage: top = [0, 0, 1, 1]
        sage: bot = [2, 2, 3, 3, 4, 4, 5, 5]
        sage: cylinder_widths(top, bot, [746825, 30952, 206252, 348, 538642, 32535])
        [1L]
        sage: cylinder_widths(top, bot, [219989, 91122, 70542, 187806, 5742, 47021])
        [1L, 2L]
        sage: cylinder_widths(top, bot, [49622, 73834, 4390, 46944, 43358, 28764])
        [2L, 10L]
        sage: cylinder_widths(top, bot, [25489, 109646, 25473, 7492, 11490, 90680])
        [1L, 2L]

    Genus 3 example::

        sage: top = [0, 1, 2, 3, 4, 5, 3, 4, 6]
        sage: bot = [7, 8, 1, 7, 8, 2, 5, 6, 0]
        sage: cylinder_widths(top, bot, [8, 10, 8, 12, 7, 3, 10, 11, 8])
        [1L]
        sage: cylinder_widths(top, bot, [9, 10, 8, 1, 6, 3, 5, 4, 3])
        [1L, 1L]
        sage: cylinder_widths(top, bot, [11, 10, 8, 12, 7, 3, 10, 11, 8])
        [1L, 1L, 2L]
        sage: cylinder_widths(top, bot, [948, 803, 775, 12, 7, 951, 10, 11, 8])
        [1L]
    """
    cdef int_iet_t t
    cdef uint64_t * widths
    set_int_iet(t, top, bot, lengths)
    widths = <uint64_t *> check_malloc((len(bot) + len(top)) / 2 * sizeof(uint64_t))
    cdef int res = int_iet_num_cylinders(widths, t)
    int_iet_clear(t)
    output =  [widths[i] for i in range(res)]
    sig_free(widths)
    output.sort()
    return output
