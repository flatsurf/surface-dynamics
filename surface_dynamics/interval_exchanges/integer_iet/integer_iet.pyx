from libc.stdlib cimport malloc, free

def perm_to_twin(p):
    r"""
    Return the twin associated to the permutation ``p``.

    EXAMPLES::

        sage: perm_to_twin([0,0,1,1,2,2])
        [1, 0, 3, 2, 5, 4]
    """
    n = len(p)
    if n%2 == 1:
        raise ValueError("p must have even length")
    n = n/2
    d = {}
    for i,j in enumerate(p):
        j = int(j)
        if j < 0 or j >= n:
            raise ValueError("p must contain integer only between 0 and {}".format(n-1))
        if j not in d:
            d[j] = []
        d[j].append(i)

    t = [None] * len(p)
    for l in d.itervalues():
        if len(l) != 2:
            raise ValueError("each label must appear twice")
        k0,k1 = l
        t[k0] = k1
        t[k1] = k0

    return t

cdef set_int_iet(int_iet_t t, top, bot, lengths):
    cdef int * labels
    cdef int * twin
    cdef uint64_t * clengths
    cdef int k = len(top)
    cdef int n = (len(top)+len(bot))/2
    cdef int i,j

    # quick checks
    lengths = map(int, lengths)
    assert len(lengths) == n
    p = map(int, top + bot)
    assert len(p) == 2*n
    for i in range(n):
        assert p.count(i) == 2
    assert all(l > 0 for l in lengths)
    assert sum(lengths[i] for i in top) == sum(lengths[i] for i in bot)

    python_twin = perm_to_twin(top + bot)

    # fill C data
    labels = <int *> malloc(2*n*sizeof(int));
    twin = <int *> malloc(2*n*sizeof(int));
    clengths = <uint64_t *> malloc(n*sizeof(uint64_t));
    if(labels == NULL or twin == NULL or clengths == NULL):
        free(labels)
        free(twin)
        free(clengths)
        raise MemoryError("not able to allocate memory!")

    for i in range(2*n):
        labels[i] = p[i]
        twin[i] = python_twin[i]
    for j in range(n):
        clengths[j] = lengths[j]

    int_iet_init(t, n)
    int_iet_set_labels_and_twin(t, labels, twin, k)
    int_iet_set_lengths(t, clengths)
    free(twin)
    free(labels)
    free(clengths)

    if int_iet_check(t):
        int_iet_clear(t)
        raise ValueError("problem with input data")


def number_of_cylinders(top, bot, lengths):
    r"""
    Return the number of cylinders of a given interval exchange
    """
    cdef int_iet_t t
    set_int_iet(t, top, bot, lengths)
    cdef int res = int_iet_num_cylinders(t)
    int_iet_clear(t)
    return res
