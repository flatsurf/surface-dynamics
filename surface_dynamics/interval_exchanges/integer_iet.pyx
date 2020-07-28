#*****************************************************************************
#       Copyright (C) 2019 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

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
    lengths = list(map(int, lengths))
    if len(lengths) != n:
        raise ValueError('invalid lengths')
    p = list(map(int, top + bot))
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

cdef void get_stat(statistics, bint flat, int_iet_t t, uint64_t * widths, uint64_t * heights, int kind):
    r"""
    - ``kind=0``: number of cylinders
    - ``kind=1``: number of components
    - ``kind=2``: tuple of cylinder widths (kind=0 is the length and kind=1 is
      the sum)
    - ``kind=3``: tuple of heights
    - ``kind=4``: tuple of pairs (width, height)
    """
    if kind == 0:
        j = int_iet_num_cylinders(NULL, NULL, t)
        key = j
    elif kind == 1:
        if widths == NULL:
            raise ValueError("widths should be allocated")
        j = int_iet_num_cylinders(widths, NULL, t)
        s = 0
        for i in range(j):
            s += widths[i]
        key = s
    elif kind == 2:
        if widths == NULL:
            raise ValueError("widths should be allocated")
        j = int_iet_num_cylinders(widths, NULL, t)
        key = tuple(sorted(widths[i] for i in range(j)))
    elif kind == 3:
        if heights == NULL:
            raise ValueError("heights should be allocated")
        j = int_iet_num_cylinders(NULL, heights, t)
        key = tuple(sorted(heights[i] for i in range(j)))
    elif kind == 4:
        if widths == NULL or heights == NULL:
            raise ValueError("widths and heights should be allocated")
        j = int_iet_num_cylinders(widths, heights, t)
        key = tuple(sorted((widths[i], heights[i]) for i in range(j)))
    else:
        raise ValueError("unknown kind %d; must be one of 0, 1, 2, 3, 4" % kind)

    if flat:
        (<list> statistics).append(key)
    else:
        if key not in (<dict> statistics):
            (<dict> statistics)[key] = 1
        else:
            (<dict> statistics)[key] += 1

def interval_exchange_statistics(top, bot, uint64_t L, int kind=0, bint flat=False):
    r"""
    Return the statistics about cylinder decomposition of all the integral iet
    with given permutation ``top`` and ``bot`` and total length ``L``.

    An integral interval exchange transformation (or iet for short)
    decomposes into periodic components. We want to study these
    components. Equivalently, such iet corresponds to a permutation
    and we want to study its cycle decomposition.

    INPUT:

    - ``top``, ``bot`` -- top and bottom permutation (as list of
      numbers starting from ``0``)

    - ``L`` (integer) -- the total length of the integral iet to be
      tested

    - ``kind`` -- if ``0`` (default) return statistics about the number of
      cylinders, if ``1`` return statistics about the number of components, if
      ``2`` the keys are tuples of cylinder heights

    OUTPUT: a dictionary whose keys are integers and the value associated to a
    key ``k`` is the number of length data with sum ``L`` which corresponds to
    an iet with data ``top`` and ``bot`` and whose cylinder decomposition has
    ``k`` cylinders.

    EXAMPLES:

    Rotations always give a unique cylinder::

        sage: from surface_dynamics.interval_exchanges.integer_iet import interval_exchange_statistics
        sage: interval_exchange_statistics([0,1], [1,0], 10, 0)
        {1: 9}

     Though, the widths of this unique cylinder is the gcd of the two lengths::

        sage: interval_exchange_statistics([0,1], [1,0], 10, 1)
        {1L: 4, 2L: 4, 5L: 1}
        sage: [gcd(k, 10-k) for k in range(1,10)]
        [1, 2, 1, 2, 5, 2, 1, 2, 1]

    Complete information about the heights of cylinders are obtained by setting ``kind=2``::

        sage: s2 = interval_exchange_statistics([0,1,2,3], [3,2,1,0], 10, 2)
        sage: for k in sorted(s2):
        ....:     print("%-6s: %s" % ('(' + ', '.join('%d' % i for i in k) + ')', s2[k]))
        (1)   : 34
        (1, 1): 16
        (1, 2): 16
        (1, 3): 4
        (1, 4): 6
        (1, 6): 2
        (2, 2): 4
        (2, 3): 2

    You can recover ``kind=0`` and ``kind=1`` by looking and lengths and sums::

        sage: from collections import defaultdict
        sage: s0 = defaultdict(int)
        sage: s1 = defaultdict(int)
        sage: for k, v in s2.items():
        ....:     s0[len(k)] += v
        ....:     s1[sum(k)] += v

        sage: dict(s0) == interval_exchange_statistics([0,1,2,3], [3,2,1,0], 10, 0)
        True
        sage: dict(s1) == interval_exchange_statistics([0,1,2,3], [3,2,1,0], 10, 1)
        True

    To get statistics for the heights, consider using ``kind=3``::

        sage: s3 = interval_exchange_statistics([0,1,2,3], [3,2,1,0], 10, 3)
        sage: for k in sorted(s3):
        ....:     print("%-6s: %s" % ('(' + ', '.join('%d' % i for i in k) + ')', s3[k]))
        (1, 4): 4
        (1, 6): 4
        (1, 8): 4
        (2, 2): 4
        (2, 3): 2
        (2, 4): 8
        (2, 6): 6
        (2, 8): 8
        (3, 4): 2
        (4, 6): 8
        (10)  : 34

    And for both lengths and heights, ``kind=4``::

        sage: s4 = interval_exchange_statistics([0,1,2,3], [3,2,1,0], 10, 4)
        sage: for key in s4:
        ....:     assert sum(l*h for l,h in key) == 10
    """
    cdef int * labels
    cdef int * twin
    cdef uint64_t * widths
    cdef uint64_t * heights
    cdef uint64_t s
    cdef int k = len(top)
    cdef int n = (len(top)+len(bot))/2
    cdef int i
    cdef int_iet_t t
    cdef li_vector_iterator_t v

    if flat:
        statistics = []
    else:
        statistics = {}

    top, bot, kfree, ktop, kbot = _relabel(top, bot)

    p = top + bot
    python_twin = perm_to_twin(p)

    labels = <int *> check_malloc(2 * n * sizeof(int))
    twin = <int *> check_malloc(2 * n * sizeof(int))
    widths = <uint64_t *> check_malloc(n * sizeof(uint64_t))
    heights = <uint64_t *> check_malloc(n * sizeof(uint64_t))

    for i in range(2*n):
        labels[i] = p[i]
        twin[i] = python_twin[i]
    int_iet_init(t, n)

    int_li_vector_init(v, L, kfree, ktop, kbot)
    int_li_vector_prefirst(v)
    while int_li_vector_first_or_next(v):
        int_iet_set_labels_and_twin(t, labels, twin, k)
        int_iet_set_lengths(t, v.x)
        get_stat(statistics, flat, t, widths, heights, kind)

    sig_free(labels)
    sig_free(twin)
    sig_free(widths)
    sig_free(heights)
    int_iet_clear(t)
    int_li_vector_clear(v)

    return statistics

def interval_exchange_statistics_sample(top, bot, uint64_t L, uint64_t sample_size, int kind=0, bint flat=False):
    r"""
    Return the statistics of cylinder decomposition of a random sample of ``sample_size``
    integral iet with given permutation ``top`` and ``bot``.

    INPUT:

    - ``top``, ``bot`` -- top and bottom permutation (as list of
      numbers starting from ``0``)

    - ``L`` (integer) -- a magnitude order for the total length of the integral
      iet to be generated

    - ``kind`` -- if ``0`` (default) return statistics about the number of
      cylinders, if ``1`` return statistics about the number of components, if
      ``2`` the keys are tuples of cylinder heights

    OUTPUT: a dictionary whose keys are integers and the value associated to a
    key ``k`` is the number of length data with sum ``L`` which corresponds to
    an iet with data ``top`` and ``bot`` and whose cylinder decomposition has
    ``k`` cylinders.

    EXAMPLES:

    Rotations always give a unique cylinder::

        sage: from surface_dynamics.interval_exchanges.integer_iet import interval_exchange_statistics
        sage: interval_exchange_statistics([0,1],[1,0],10,0)
        {1: 9}

     Though, the widths of this unique cylinder is the gcd of the two lengths::

        sage: interval_exchange_statistics([0,1],[1,0],10,1)
        {1L: 4, 2L: 4, 5L: 1}
        sage: [gcd(k, 10-k) for k in range(1,10)]
        [1, 2, 1, 2, 5, 2, 1, 2, 1]

    Complete information about the heights of cylinders are obtained by setting ``kind=2``::

        sage: s2 = interval_exchange_statistics([0,1,2,3], [3,2,1,0], 10, 2)
        sage: for k in sorted(s2):
        ....:     print("%-6s: %s" % ('(' + ', '.join('%d' % i for i in k) + ')', s2[k]))
        (1)   : 34
        (1, 1): 16
        (1, 2): 16
        (1, 3): 4
        (1, 4): 6
        (1, 6): 2
        (2, 2): 4
        (2, 3): 2

    You can recover ``kind=0`` and ``kind=1`` by looking and lengths and sums::

        sage: from collections import defaultdict
        sage: s0 = defaultdict(int)
        sage: s1 = defaultdict(int)
        sage: for k, v in s2.items():
        ....:     s0[len(k)] += v
        ....:     s1[sum(k)] += v

        sage: dict(s0) == interval_exchange_statistics([0,1,2,3], [3,2,1,0], 10, 0)
        True
        sage: dict(s1) == interval_exchange_statistics([0,1,2,3], [3,2,1,0], 10, 1)
        True
    """
    cdef int * labels
    cdef int * twin
    cdef uint64_t * widths
    cdef uint64_t * heights
    cdef uint64_t s
    cdef int k = len(top)
    cdef int n = (len(top)+len(bot))/2
    cdef int i
    cdef int count
    cdef int_iet_t t

    if flat:
        statistics = []
    else:
        statistics = {}

    top, bot, kfree, ktop, kbot = _relabel(top, bot)

    p = top + bot
    python_twin = perm_to_twin(p)

    labels = <int *> check_malloc(2 * n * sizeof(int))
    twin = <int *> check_malloc(2 * n * sizeof(int))
    widths = <uint64_t *> check_malloc(n * sizeof(uint64_t))
    heights = <uint64_t *> check_malloc(n * sizeof(uint64_t))

    for i in range(2*n):
        labels[i] = p[i]
        twin[i] = python_twin[i]
    int_iet_init(t, n)

    for count in range(sample_size):
        int_iet_set_labels_and_twin(t, labels, twin, k)
        int_iet_set_random_lengths(t, L)
        get_stat(statistics, flat, t, widths, heights, kind)

    sig_free(labels)
    sig_free(twin)
    sig_free(widths)
    sig_free(heights)
    int_iet_clear(t)

    return statistics


# there are several interesting statistics one might want to compute
# such as
#  - number of cylinders
#  - sum of heights of cylinders
#  - mean ratio in view of Siegel Veech constants (not available for now) as
#    we would need the circumferences of cylinders as well
def cylinder_statistics(top, bot, uint64_t L, int kind=0, bint flat=False):
    r"""
    Return the statistics of the number of cylinders for a given total length

    INPUT:

    - ``top`` -- (list) composition of the top of the cylinder

    - ``bot`` -- (list) composition of the bottom of the cylinder
    
    - ``L`` -- (positive integer) the length we are considering

    - ``kind`` -- kind of statistics to gather
      - ``0`` number of cylinders
      - ``1`` number of components
      - ``2`` tuple of widths
      - ``3`` tuple of heights
      - ``4`` tuple of pairs (width, height)

    EXAMPLES::

        sage: from surface_dynamics.interval_exchanges.integer_iet import cylinder_statistics

    The case of H(0,0)::

        sage: top = [0, 1]
        sage: bot = [1, 0]
        sage: for n in range(2,20):
        ....:     s = cylinder_statistics(top, bot, n)
        ....:     print("%2d : %3d %3d %3d" %(n, s[1], s[2], s[1] + s[2]))
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
        ....:     s1 = cylinder_statistics(top1, bot1, n)
        ....:     s2 = cylinder_statistics(top2, bot2, 2*n)
        ....:     assert sorted(s1.keys()) == sorted(s2.keys())
        ....:     assert all(2*s1[i] == s2[i] for i in s1)
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
    cdef uint64_t * widths
    cdef uint64_t * heights
    cdef int_iet_t t1, t2
    cdef li_vector_iterator_t v

    if flat:
        statistics = []
    else:
        statistics = {}

    top, bot, kfree, ktop, kbot = _relabel(top, bot)

    p1 = top + bot
    p2 = [n] + p1 + [n]
    python_twin1 = perm_to_twin(p1)
    python_twin2 = perm_to_twin(p2)

    labels1 = <int *> check_malloc(2*n*sizeof(int))
    twin1 = <int *> check_malloc(2*n*sizeof(int))
    labels2 = <int *> check_malloc(2*(n+1)*sizeof(int))
    twin2 = <int *> check_malloc(2*(n+1)*sizeof(int))
    widths = <uint64_t *> check_malloc((n+1) * sizeof(uint64_t))
    heights = <uint64_t *> check_malloc((n+1) * sizeof(uint64_t))

    clengths = <uint64_t *> check_malloc((n+1)*sizeof(uint64_t))

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
        get_stat(statistics, flat, t1, widths, heights, kind)

        # positive twists case
        for twist in range(1, L):
            clengths[n] = twist
            int_iet_set_labels_and_twin(t2, labels2, twin2, k2)
            int_iet_set_lengths(t2, clengths)
            get_stat(statistics, flat, t2, widths, heights, kind)

    sig_free(labels1)
    sig_free(twin1)
    sig_free(labels2)
    sig_free(twin2)
    sig_free(clengths)
    sig_free(widths)
    sig_free(heights)
    int_iet_clear(t1)
    int_iet_clear(t2)
    int_li_vector_clear(v)

    return statistics

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
    set_int_iet(t, top, bot, lengths)
    cdef int res = int_iet_num_cylinders(NULL, NULL, t)
    int_iet_clear(t)
    return res

def cylinder_widths_and_heights(top, bot, lengths):
    r"""
    Return the list of pairs ``(width, height)`` of each periodic component of the
    iet determined by ``top``, ``bot`` and ``lengths``.

    EXAMPLES::

        sage: from surface_dynamics.interval_exchanges.integer_iet import cylinder_widths_and_heights

        sage: cylinder_widths_and_heights([0,1],[1,0],[5,2])
        [(1L, 7L)]
        sage: cylinder_widths_and_heights([0,1],[1,0],[8,12])
        [(4L, 5L)]

        sage: cylinder_widths_and_heights([0,1,2],[2,1,0],[5,3,2])
        [(1L, 10L)]

        sage: cylinder_widths_and_heights(list(range(8)), list(range(7,-1,-1)), [1386,924,660,495,385,308,252,210])
        [(2L, 930L), (3L, 920L)]
    """
    cdef int_iet_t t
    cdef uint64_t * widths
    cdef uint64_t * heights
    set_int_iet(t, top, bot, lengths)
    widths = <uint64_t *> check_malloc((len(bot) + len(top)) / 2 * sizeof(uint64_t))
    heights = <uint64_t *> check_malloc((len(bot) + len(top)) / 2 * sizeof(uint64_t))
    cdef int res = int_iet_num_cylinders(widths, heights, t)
    int_iet_clear(t)
    output =  [(widths[i], heights[i]) for i in range(res)]
    sig_free(widths)
    sig_free(heights)
    output.sort()
    return output

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
    cdef int res = int_iet_num_cylinders(widths, NULL, t)
    int_iet_clear(t)
    output =  [widths[i] for i in range(res)]
    sig_free(widths)
    output.sort()
    return output
