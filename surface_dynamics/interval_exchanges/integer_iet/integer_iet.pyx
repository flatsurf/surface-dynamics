from cysignals.memory cimport check_malloc, sig_free

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
