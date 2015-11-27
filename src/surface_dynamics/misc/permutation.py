r"""
Permutation in Sage works by default on `{1, 2, ..., n}` but it might be much
more convenient to work on `{0, 1, ..., n-1}`. This module provide simple
functions for the latter representation.
"""

def permutation_to_perm(p):
    r"""
    Returns a list on `[0, n-1]` from a permutation on `[1, n]`

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import permutation_to_perm
        sage: permutation_to_perm(PermutationGroupElement([3,1,2]))
        [2, 0, 1]
    """
    return map(lambda x: x-1, p.domain())


def perm_to_permutation(l):
    r"""
    Returns a permutation on `[1, n]` from a list on `[0, n-1]`

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_to_permutation
        sage: perm_to_permutation([2,1,0])
        (1,3)
    """
    from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
    return PermutationGroupElement(map(lambda x: x+1, l))


def cycles_to_list(t):
    r"""
    Returns a permutation on `[0, n-1]` from a list of cycles on `[0, n-1]`

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import cycles_to_list
        sage: cycles_to_list([[1,3,5],[0,2,4],[6]])
        [2, 3, 4, 5, 0, 1, 6]

        sage: cycles_to_list([])
        []
        sage: cycles_to_list([[],[]])
        []
    """
    if not any(tt for tt in t):
        return []

    res = range(max(map(max, t))+1)

    for c in t:
        for j in xrange(len(c)-1):
            res[c[j]] = int(c[j+1])
        res[c[-1]] = int(c[0])

    return res


def str_to_cycles(s):
    """
    Returns a list of cycles from a string

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import str_to_cycles
        sage: str_to_cycles('(0,1)')
        [[0, 1]]
        sage: str_to_cycles('(0,1)(3,2)')
        [[0, 1], [3, 2]]

        sage: str_to_cycles('()(0,1)()(2,3)')
        [[0, 1], [2, 3]]
    """
    r = []
    for c_str in s[1:-1].split(')('):
        if not c_str:
            continue
        r.append(map(int, c_str.replace(' ', '').split(',')))
    return r


def init_perm(data):
    """
    Returns a permutation from different kinds of data

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import init_perm
        sage: init_perm([3,2,1,4])
        [3, 2, 1, 4]
        sage: init_perm(([2,1],[3,4,0]))
        [3, 2, 1, 4, 0]
        sage: init_perm('(0,1)(3,2)')
        [1, 0, 3, 2]
        sage: init_perm([3,1,None,0])
        [3, 1, None, 0]
        sage: init_perm([[2,1],[3,4,0]])
        [3, 2, 1, 4, 0]
        sage: init_perm([])
        []
        sage: init_perm([[]])
        []
    """
    if isinstance(data, (tuple,list)):
        if not data:
            return []
        if isinstance(data[0], (tuple,list)):
            return cycles_to_list(data)
        else:
            return [x if (x is None or isinstance(x,int)) else int(x) for x in data]

    if isinstance(data, str):
        c = str_to_cycles(data)
        return cycles_to_list(c)

    raise TypeError("The input must be list, tuple or string")


def equalize_perms(l):
    """
    Ensures that permutations have the same size ?

    INPUT:

    a list of permutations

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import equalize_perms
        sage: li = [[0, 1], [2, 1, 0]]
        sage: equalize_perms(li)
        3
    """
    n = max(map(len, l))
    for p in l:
        p.extend(xrange(len(p), n))
    return n


def perm_check(l):
    r"""
    Checks that `l` is a permutation of `[0, n-1]` for some `n`

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_check

    Good permutations::

        sage: perm_check([1r, 0r, 3r, 2r])
        sage: perm_check([None, None])
        sage: perm_check([None, 3r, None, 1r])

    Bad permutations::

        sage: perm_check([1, 0, 3, 2])
        Traceback (most recent call last):
        ...
        TypeError: entries must be int or None, not Integer
        sage: perm_check([2r, 0r])
        Traceback (most recent call last):
        ...
        ValueError: permutation value 2 out of range
        sage: perm_check([1r,0r,1r])
        Traceback (most recent call last):
        ...
        ValueError: 1 is repeated
    """
    n = len(l)
    seen = [False]*n
    for i in xrange(n):
        if l[i] is None:
            continue
        if type(l[i]) != int:
            raise TypeError("entries must be int or None, not {}".format(type(l[i]).__name__))
        if l[i] < 0 or l[i] >= n:
            raise ValueError("permutation value {} out of range".format(l[i]))
        if seen[l[i]]:
            raise ValueError("{} is repeated".format(l[i]))
        seen[l[i]] = True


def perm_invert(l):
    r"""
    Returns the inverse of the permutation `l`

    TESTS::

        sage: from itertools import permutations
        sage: from surface_dynamics.misc.permutation import perm_invert, perm_compose
        sage: all(perm_compose(perm_invert(p),p) == range(3) for p in permutations(range(3)))
        True
        sage: all(perm_compose(p,perm_invert(p)) == range(3) for p in permutations(range(3)))
        True

        sage: perm_invert([2, None, 5, 0, None, 3])
        [3, None, 0, 5, None, 2]
    """
    res = [0]*len(l)
    for i in xrange(len(l)):
        if l[i] is None:
            res[i] = None
        else:
            res[l[i]] = i
    return res


def perm_compose(l1, l2):
    r"""
    Returns the product `p_1 * p_2` where `p_1` and `p_2` are the permutations
    associated to `l_1` and `l_2`

    INPUT:

    two lists that are permutations on `[0, n-1]`.

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_compose
        sage: perm_compose([0,2,1],[0,2,1])
        [0, 1, 2]
        sage: perm_compose([None,2,3,1],[None,2,1,3])
        [None, 1, 3, 2]
    """
    r = [None] * len(l1)
    for i in xrange(len(l1)):
        if l1[i] is not None and l1[i] < len(l2):
            r[i] = l2[l1[i]]
    return r

def perm_compose_i(l1, l2):
    r"""
    Returns the product `p_1^{-1} * p_2^{-1}` where `p_1` and `p_2` are the
    permutations associated to the list `l_1` and `l_2`

    INPUT:

    two lists that are permutations on `[0, n-1]`.

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_compose_i
        sage: perm_compose_i([0,1,2],[1,2,0])
        [2, 0, 1]
    """
    assert(len(l1) == len(l2))

    res = [None]*len(l1)
    for i in xrange(len(l1)):
        res[l2[l1[i]]] = i

    return res


def perm_orbit(p, i):
    r"""
    Returns the orbit of an integer `i` under the permutation `p`

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_orbit
        sage: perm_orbit([0,3,1,2],2)
        [2, 1, 3]
    """
    res = [i]
    j = p[i]
    while j != i:
        res.append(j)
        j = p[j]
    return res


def perm_cycle_tuples(p, singletons=False):
    r"""
    Return the cycle decomposition of `p`

    INPUT:

    - ``p`` -- the permutation

    - ``singletons`` -- bool (default: ``False``) - return or not the singletons

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_cycle_tuples
        sage: perm_cycle_tuples([0,2,1])
        ([1, 2],)
        sage: perm_cycle_tuples([0,2,1],True)
        ([0], [1, 2])

        sage: perm_cycle_tuples([2,None,0])
        ([0, 2],)
    """
    seen = [1]*len(p)
    res = []

    for i in xrange(len(p)):
        if seen[i] and p[i] is not None:
            cycle = []
            j = i
            while seen[j]:
                seen[j] = 0
                cycle.append(j)
                j = p[j]
            if singletons or len(cycle) > 1:
                res.append(cycle)

    return tuple(res)


def perm_cycle_string(p, singletons=False):
    r"""
    Returns a string representing the cycle decomposition of `p`

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_cycle_string
        sage: perm_cycle_string([0,2,1])
        '(1,2)'
        sage: perm_cycle_string([0,2,1],True)
        '(0)(1,2)'
    """
    return ''.join(map(lambda x: '('+','.join(map(str, x))+')',
                       perm_cycle_tuples(p, singletons)))


def canonical_perm(part,i=0):
    r"""
    Return the canonical permutation with the given part

    EXAMPLES::

        sage: from surface_dynamics.flat_surfaces.separatrix_diagram import canonical_perm
        sage: canonical_perm([3,2])
        [1, 2, 0, 4, 3]
        sage: canonical_perm([2,2,2])
        [1, 0, 3, 2, 5, 4]
        sage: canonical_perm([1,1,3])
        [0, 1, 3, 4, 2]
    """
    res = []
    for p in part:
        res.extend(xrange(i+1,i+p))
        res.append(i)
        i += p
    return res


def canonical_perm_i(part,i=0):
    r"""
    Return the canonical permutation reversed

    EXAMPLES::

        sage: from surface_dynamics.flat_surfaces.separatrix_diagram import canonical_perm_i
        sage: canonical_perm_i([3,2])
        [2, 0, 1, 4, 3]
        sage: canonical_perm_i([2,2,2])
        [1, 0, 3, 2, 5, 4]
        sage: canonical_perm_i([1,1,3])
        [0, 1, 4, 2, 3]
    """
    res = []
    for p in part:
        res.append(i+p-1)
        res.extend(xrange(i,i+p-1))
        i += p
    return res


def perm_switch(p1, p2, i, j):
    """
    Exchanges the values at positions `i` and `j` in two permutations `p_1` and `p_2`

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_switch
        sage: a = [0,1,2]
        sage: b = [1,2,0]
        sage: perm_switch(a,b,0,1)
        sage: a;b
        [1, 0, 2]
        [2, 1, 0]
    """
    i1 = p1[i]
    j1 = p1[j]
    p1[i] = j1
    p1[j] = i1

    i2 = p2[i1]
    j2 = p2[j1]
    p2[i1] = j2
    p2[j1] = i2


def perms_transitive_components(p):
    r"""
    Return the list of transitive components of the subgroup generated by the
    permutations ``p``.

    INPUT:

    - ``p`` -- a list of permutations given as lists

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perms_transitive_components

        sage: perms_transitive_components([[1,0,2,3],[0,1,3,2]])
        [(0, 1), (2, 3)]

        sage: perms_transitive_components([[2,3,0,1]])
        [(0, 2), (1, 3)]

        sage: perms_transitive_components([[3,1,2,0], [0,3,2,1], [0,1,3,2]])
        [(0, 1, 2, 3)]
    """
    n = len(p[0])
    seen = [-1] * n
    cc_num = 0
    for i in range(n):
        if seen[i] != -1:
            continue

        todo = [i]
        seen[i] = cc_num
        while todo:
            j = todo.pop()
            for pp in p:
                k = pp[j]
                if seen[k] == -1:
                    todo.append(k)
                    seen[k] = cc_num
        cc_num += 1

    return [tuple(i for i in range(n) if seen[i] == j) for j in range(cc_num)]


def perms_are_transitive(p):
    """
    Test whether the group generated by the permutations in ``p`` is transitive.

    INPUT:

    - ``p`` - a list of permutations of `[0, n-1]`

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perms_are_transitive
        sage: perms_are_transitive([[0,1,2],[0,2,1]])
        False
        sage: perms_are_transitive([[0,1,2],[1,2,0]])
        True

        sage: p0 = [0,2,3,1,7,5,6,4]
        sage: p1 = [7,1,2,3,4,5,6,0]
        sage: p2 = [6,1,2,3,4,5,0,7]
        sage: p3 = [1,0,2,3,4,5,6,7]
        sage: p4 = [0,1,4,5,2,3,6,7]
        sage: perms_are_transitive([p0,p1,p2,p3,p4])
        True
        sage: perms_are_transitive([p0,p1,p2,p3])
        False
    """
    if not p:
        raise ValueError("empty list")

    n = len(p[0])

    # compute the connected component of 0
    cc0 = [False] * n
    todo = [0]
    cc0[0] = True
    while todo:
        j = todo.pop()
        for pp in p:
            k = pp[j]
            if cc0[k] is False:
                todo.append(k)
                cc0[k] = True

    return all(j is True for j in cc0)


def perms_relabel(p, m):
    """
    `p` is a list of permutations

    `m` is a permutation

    The result is the relabelling of the permutations in `p` according to `m`.

    EXAMPLES::

        sage: from surface_dynamics.misc.constellation import perms_relabel
        sage: perms_relabel([[0,1,2],[0,2,1]],[2,1,0])
        [[0, 1, 2], [1, 0, 2]]
    """
    q = [k[:] for k in p]
    for i in xrange(len(m)):
        for j in xrange(len(p)):
            q[j][m[i]] = m[p[j][i]]
    return q


def perms_canonical_labels_from(x, y, j0):
    r"""
    Return canonical labels for ``x``, ``y`` that starts at ``j0``

    .. WARNING:

        The group generated by ``x`` and the elements of ``y`` should be
        transitive.

    INPUT:

    - ``x`` -- list - a permutation of `[0, ..., n]` as a list

    - ``y`` -- list of permutations of `[0, ..., n]` as a list of lists

    - ``j0`` -- an index in [0, ..., n]

    OUTPUT:

    mapping: a permutation that specify the new labels

    EXAMPLES::

        sage: from surface_dynamics.misc.constellation import perms_canonical_labels_from
        sage: perms_canonical_labels_from([0,1,2],[[1,2,0]],0)
        [0, 1, 2]
        sage: perms_canonical_labels_from([1,0,2],[[2,0,1]],0)
        [0, 1, 2]
        sage: perms_canonical_labels_from([1,0,2],[[2,0,1]],1)
        [1, 0, 2]
        sage: perms_canonical_labels_from([1,0,2],[[2,0,1]],2)
        [2, 1, 0]
    """
    n = len(x)

    k = 0
    mapping = [None] * n
    waiting = [[] for i in xrange(len(y))]

    while k < n:
        # initialize at j0
        mapping[j0] = k
        waiting[0].append(j0)
        k += 1
        # complete x cycle from j0
        j = x[j0]
        while j != j0:
            mapping[j] = k
            waiting[0].append(j)
            k += 1
            j = x[j]

        # find another guy
        l = 0
        while l < len(waiting):
            i = 0
            while i < len(waiting[l]):
                j1 = waiting[l][i]
                if mapping[y[l][j1]] is None:
                    break
                i += 1

            if i == len(waiting[l]):  # not found: go further in waiting
                if l < len(waiting)-1:
                    waiting[l+1].extend(waiting[l])
                waiting[l] = []
                l += 1
                i = 0

            else:  # found: complete cycle from new guy
                j0 = y[l][j1]
                if l < len(waiting)-1:
                    waiting[l+1].extend(waiting[l][:i+1])
                del waiting[l][:i+1]
                break

    return mapping


def perms_canonical_labels(p, e=None):
    assert(len(p) > 1)
    n = len(p[0])

    c_win = None
    m_win = range(n)

    x = p[0]
    y = p[1:]

    if e is None:
        e = range(n)

    # get canonical label from i in to_test and compare
    while e:
        i = e.pop()
        m_test = perms_canonical_labels_from(x, y, i)
        c_test = perms_relabel(p, m_test)
        if c_win is None or c_test < c_win:
            c_win = c_test
            m_win = m_test

    return c_win, m_win


