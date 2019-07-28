r"""
Even permutations

We call even permutation, a permutation on the set of even size
`X_n = \{-n, -n+1, \ldots, n-1\}`. The set `X_n` comes with a canonical
involution without fixed points `i \mapsto -i-1` (note that `-i-1` is
the bit complement to `i`). The commutator of this canonical involution
are the signed permutations.
"""
#*****************************************************************************
#       Copyright (C) 2019 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, print_function
from array import array

def is_signed_perm(p):
    if not isinstance(p, array):
        p = array('l', p)
    if len(p) % 2:
        return False
    n = len(p)//2
    for i in range(-n, n):
        if p[i] != ~p[~i]:
            return False
    return True

def signed_permutations(n):
    r"""
    EXAMPLES::

        sage: from surface_dynamics.misc.even_permutation import signed_permutations
        sage: for p in signed_permutations(2):
        ....:     assert p[0r] == ~p[~0r]
        ....:     assert p[1r] == ~p[~1r]
    """
    from itertools import permutations, combinations

    R = list(range(n))
    for p in permutations(R):
        for k in range(n+1):
            for s in combinations(R, k):
                q = array('l', p + tuple(~i for i in reversed(p)))
                for i in s:
                    q[i] = ~q[i]
                    q[~i] = ~q[~i]
                yield q

def even_permutations(n):
    from itertools import permutations
    for p in permutations(range(-n, n)):
        yield array('l', p)

def even_perm_init(data):
    if isinstance(data, (tuple, list)):
        if not data:
            return []
        if isinstance(data[0], (tuple,list)):
            return even_perm_from_cycles(data)
        else:
            return [int(x) for x in data]
    if isinstance(data,str):
        raise NotImplementedError

    raise TypeError("The input must be list, tuple or string")

def even_perm_check(p):
    r"""
    EXAMPLES::

        sage: from surface_dynamics.misc.even_permutation import even_perm_check
        sage: even_perm_check([0, 1, 2, -3, -2, -1])
        sage: even_perm_check([0])
        Traceback (most recent call last):
        ...
        ValueError: odd length

        sage: even_perm_check([-2,0])
        Traceback (most recent call last):
        ...
        ValueError: element -2 at position 0 out of range

        sage: even_perm_check([0,1,-2,0])
        Traceback (most recent call last):
        ...
        ValueError: element 0 repeated
    """
    if not isinstance(p, (array, tuple, list)):
        raise TypeError('invalid type for permutation')
    if len(p) % 2:
        raise ValueError('odd length')
    seen = [False] * len(p)
    n = len(p)//2
    for i in range(-n, n):
        j = p[i]
        if not (-n <= j < n):
            raise ValueError('element {} at position {} out of range'.format(j, i))
        if seen[j]:
            raise ValueError('element {} repeated'.format(j))
        seen[j] = True

def even_perm_identity(n):
    r"""
    Return the identity

    EXAMPLES::

        sage: from surface_dynamics.misc.even_permutation import even_perm_identity

        sage: even_perm_identity(3)
        array('l', [0, 1, 2, -3, -2, -1])
        sage: even_perm_identity(0)
        array('l')
    """
    return array('l', list(range(n)) + list(range(-n,0)))

def even_perm_minus(n):
    r"""
    EXAMPLES::

        sage: from surface_dynamics.misc.even_permutation import even_perm_minus
        sage: even_perm_minus(2)
        array('l', [-1, -2, 1, 0])
    """
    return array('l', [~i for i in list(range(n))] + [~i for i in list(range(-n,0))])

def even_perm_is_orientable(p):
    r"""
    Test whether the permutation is orientable (that is do not mix non-negatives with negatives)

    EXAMPLES::

        sage: from surface_dynamics.misc.even_permutation import even_perm_is_orientable
        sage: even_perm_is_orientable([1,0,-1,-2])
        True
        sage: even_perm_is_orientable([0,1,-2,-1])
        True
        sage: even_perm_is_orientable([1,-1,0,-2])
        False
    """
    n = len(p) // 2
    for i in range(n):
        if p[i] < 0:
            return False
    return True

def even_perm_cycles(p):
    r"""
    EXAMPLES::

        sage: from surface_dynamics.misc.even_permutation import even_perm_cycles

        sage: even_perm_cycles([-1,0])
        ([(-1, 0)], [0, 0])
        sage: even_perm_cycles([-1,0,-2,1])
        ([(-2,), (-1, 1, 0)], [1, 1, 0, 1])
    """
    ans = []
    n = len(p)//2
    cycles = [-1] * (2*n)
    k = 0
    t = []
    for i in range(-n, n):
        if cycles[i] != -1:
            continue
        del t[:]
        while cycles[i] == -1:
            cycles[i] = k
            t.append(i)
            i = p[i]
        ans.append(tuple(t))
        k += 1

    return ans, cycles

def even_perm_invert(p):
    r"""
    Return the inverse of ``p``.

    EXAMPLES::

        sage: from surface_dynamics.misc.even_permutation import (even_perm_invert,
        ....:    even_perm_compose, even_perm_identity)

        sage: p = [-2,0,-1,1]
        sage: q = even_perm_invert(p)
        sage: q
        array('l', [1, -1, 0, -2])
        sage: even_perm_compose(p,q)
        array('l', [0, 1, -2, -1])
        sage: even_perm_compose(q,p)
        array('l', [0, 1, -2, -1])

        sage: for _ in range(10):
        ....:     n = randint(1, 20)
        ....:     p = list(range(-n, n))
        ....:     shuffle(p)
        ....:     q = even_perm_invert(p)
        ....:     assert even_perm_compose(p, q) == even_perm_compose(q, p) == even_perm_identity(n)
    """
    n = len(p) // 2
    q = array('l', p)
    for i in range(-n, n):
        q[p[i]] = i
    return q

def even_perm_compose(p1, p2):
    r"""
    Return the composition of ``p1`` and ``p2``

    EXAMPLES::

        sage: from surface_dynamics.misc.even_permutation import even_perm_compose

        sage: p1 = [0, -1, 2, -3, 1, -2]
        sage: p2 = [-1, 0, -2, 1, -3, 2]
        sage: q = even_perm_compose(p1, p2)
        sage: q
        array('l', [-1, 2, -2, 1, 0, -3])
        sage: p2[p1[3]] == q[3]
        True
    """
    assert len(p1) == len(p2)
    n = len(p1) // 2
    q = array('l', p1)
    for i in range(-n, n):
        q[i] = p2[p1[i]]
    return q

def even_perm_compose_i(p1, p2):
    r"""
    Return the compositions of the inverses of ``p1`` and ``p2``

    EXAMPLES::

        sage: from surface_dynamics.misc.even_permutation import (even_perm_compose,
        ....:      even_perm_invert, even_perm_compose_i)
        sage: p1 = [0, -2, 2, 1, -1, -3]
        sage: p2 = [1, -1, -2, -3, 2, 0]
        sage: even_perm_compose_i(p1, p2)
        array('l', [-1, -3, -2, 1, 0, 2])
        sage: even_perm_compose(even_perm_invert(p1), even_perm_invert(p2))
        array('l', [-1, -3, -2, 1, 0, 2])
    """
    assert len(p1) == len(p2)
    n = len(p1) // 2
    q = array('l', p1)
    for i in range(-n, n):
        q[p1[p2[i]]] = i
    return q

def even_perm_is_transitive(p):
    r"""
    EXAMPLES::

        sage: from array import array
        sage: from surface_dynamics.misc.even_permutation import even_perm_is_transitive

        sage: p = array('l', [0, 1, -2, -1])
        sage: even_perm_is_transitive(p)
        False
        sage: p = array('l', [0, 1, -1, -2])
        sage: even_perm_is_transitive(p)
        True

        sage: from itertools import permutations
        sage: for p in permutations((-2,-1,0,1)):
        ....:     print("%s %s" % (p, even_perm_is_transitive(array('l',p))))
        (-2, -1, 0, 1) True
        (-2, -1, 1, 0) True
        (-2, 0, -1, 1) True
        (-2, 0, 1, -1) True
        (-2, 1, -1, 0) True
        (-2, 1, 0, -1) True
        (-1, -2, 0, 1) True
        (-1, -2, 1, 0) False
        (-1, 0, -2, 1) True
        (-1, 0, 1, -2) True
        (-1, 1, -2, 0) False
        (-1, 1, 0, -2) True
        (0, -2, -1, 1) True
        (0, -2, 1, -1) False
        (0, -1, -2, 1) True
        (0, -1, 1, -2) True
        (0, 1, -2, -1) False
        (0, 1, -1, -2) True
        (1, -2, -1, 0) True
        (1, -2, 0, -1) True
        (1, -1, -2, 0) True
        (1, -1, 0, -2) True
        (1, 0, -2, -1) True
        (1, 0, -1, -2) True

    """
    assert isinstance(p, array)
    n = len(p) // 2

    # compute the connected component of 0
    cc0 = [False] * n
    todo = [0]
    cc0[0] = True
    num = 1
    while todo and num < n:
        j = todo.pop()
        k = p[j]
        if k < 0: k = ~k
        if not cc0[k]:
            todo.append(k)
            cc0[k] = True
            num += 1
        k = p[~j]
        if k < 0: k = ~k
        if not cc0[k]:
            todo.append(k)
            cc0[k] = True
            num += 1
    return num == n

def even_perm_tilde(p):
    r"""
    EXAMPLES::

        sage: from surface_dynamics.misc.even_permutation import even_perm_tilde
        sage: even_perm_tilde([1r, 0r, -1r, -2r])
        array('l', [1, 0, -1, -2])
        sage: even_perm_tilde([1r, -1r, -2r, 0r])
        array('l', [-1, 1, 0, -2])

        sage: from surface_dynamics.misc.even_permutation import even_permutations, is_signed_perm
        sage: all((p == even_perm_tilde(p)) == is_signed_perm(p) for p in even_permutations(3))
        True
    """
    q = array('l', p)
    for i in range(len(p)):
        q[~i] = ~p[i]
    return q

def even_perm_relabel(p, m):
    """
    Relabel (= conjugate) the permutation ``p`` according to ``m``.

    INPUT:

    - ``p`` - a list of even permutations

    - ``m`` - a permutation

    EXAMPLES::

        sage: from surface_dynamics.misc.even_permutation import (even_perm_relabel,
        ....:      even_perm_identity, even_perm_tilde)
        sage: even_perm_relabel([-1, 0, -2, 1], [1, 0, -2, -1])
        array('l', [1, -1, -2, 0])

        sage: from array import array
        sage: p = array('l', [1, -2, 0, -1])
        sage: even_perm_relabel(p, even_perm_identity(2)) == p
        True

        sage: from surface_dynamics.misc.even_permutation import signed_permutations
        sage: p = array('l', [-2,0,1,2,-1,-3])
        sage: for s in signed_permutations(3):
        ....:     p1 = even_perm_relabel(even_perm_tilde(p), s)
        ....:     p2 = even_perm_tilde(even_perm_relabel(p, s))
        ....:     assert p1 == p2
    """
    assert len(p) == len(m)
    q = array('l', p)
    for i in range(len(m)):
        q[m[i]] = m[p[i]]
    return q

# TODO
# the automorphism group of a fatgraph is completely determined
# by the image of a single edge! Hence at most 2n elements.
# Below, we basically run through # all possible relabelings
# which is somehow a big waste # since we relabel in order
# (edge 0 is determined first, then edge number 1, etc),
# this should be much faster
def even_perm_canonical_label_from(p, i):
    r"""
    Return canonical labels for ``p`` that starts at ``i``.

    INPUT:

    - ``p`` -- even permutation

    - ``i`` -- integer in `X_n = {-n, -n+1, \ldots, n-1\}`

    OUTPUT: a signed permutation that specifies the new labels

    EXAMPLES::

        sage: from array import array
        sage: from surface_dynamics.misc.even_permutation import (even_perm_canonical_label_from,
        ....:     is_signed_perm)

        sage: p = array('l', [-1, 2, 0, 1, -3, -2])
        sage: even_perm_canonical_label_from(p, 0)
        array('l', [0, -2, -3, 2, 1, -1])

        sage: p = array('l', [1, 2, 0, -2, -1, -3])
        sage: is_signed_perm(even_perm_canonical_label_from(p, 0))
        True

        sage: p = array('l', [2,1,0,-1,-3,-2])
        sage: even_perm_canonical_label_from(p, 1)
        array('l', [2, 0, 1, -2, -1, -3])
    """
    assert isinstance(p, array)

    n = len(p) // 2
    mapping = array('l', [n] * len(p))
    i = int(i)
    waiting = [i, ~i]
    mapping[i] = 0
    mapping[~i] = -1
    k = int(1)

    while waiting:
        i = waiting.pop(0)
        assert mapping[i] != n

        j = p[i]
        s = mapping[i] >= 0
        while j != i:
            if mapping[j] == n:
                if s:
                    mapping[j] = k
                    mapping[~j] = ~k
                else:
                    mapping[j] = ~k
                    mapping[~j] = k
                waiting.append(~j)
                k += 1
            j = p[j]

    return mapping

def even_perm_canonical_label(p):
    r"""
    Canonical label for the permutation ``p`` assuming that
    together with the canonical involution it acts transitively
    on `X_n = \{-n, -n+1, \ldots, n-1\}`.

    For fatgraphs, we use canonical labelings for the faces, so that
    canonical labels detect the orientable fatgraphs.

    For the possible starting points, we only need to consider the
    ones that have one side in the largest and the other side in the
    smallest possible. For now, we run through all possibilities.
    
    OUTPUT: a pair ``(canonical_p, map)``

    EXAMPLES::

        sage: from surface_dynamics.misc.even_permutation import (even_perm_check,
        ....:    even_perm_canonical_label, even_perm_is_transitive, even_perm_relabel,
        ....:    signed_permutations, even_perm_tilde)
        sage: from array import array

        sage: def test_perm(p):
        ....:     p = array('l', p)
        ....:     even_perm_check(p)
        ....:     n = len(p) // 2
        ....:     assert even_perm_is_transitive(p)
        ....:     c = even_perm_canonical_label(p)[0]
        ....:     q = even_perm_tilde(p)
        ....:     cc = even_perm_canonical_label(q)[0]
        ....:     assert c == cc, (c, cc)
        ....:     for q in signed_permutations(n):
        ....:         pp = even_perm_relabel(p, q)
        ....:         cc, mm = even_perm_canonical_label(pp)
        ....:         assert cc == even_perm_relabel(pp, mm)
        ....:         assert c == cc, (q, pp, cc)

        sage: test_perm([-2,0,1,-1,2,-3])
        sage: test_perm([-2,0,1,2,-1,-3])
        sage: test_perm([2,1,0,-1,-3,-2])
        sage: test_perm([2,0,1,-1,-3,-2])

        sage: test_perm([2,0,-4,1,-1,-3,-2,3])
        sage: test_perm([0,2,1,3,-2,-1,-4,-3])
    """
    p = array('l', p)

    if not p:
        return p, []

    assert len(p) %2 == 0
    n = len(p) // 2

    m_win = p_win = None

    for i in range(-n,n):
        m_test = even_perm_canonical_label_from(p, i)
        p_test = even_perm_relabel(p, m_test)
        if p_win is None or p_test < p_win:
            p_win = p_test
            m_win = m_test

    return p_win, m_win
