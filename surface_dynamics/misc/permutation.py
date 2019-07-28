r"""
Partial permutation on `\{0, 1, ..., n-1\}` as lists.

Permutations are implemented as lists or arrays. When
the image is undefined, the number -1 is used.
"""
#*****************************************************************************
#       Copyright (C) 2019 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function, absolute_import
from six.moves import range

from math import log

import numbers

try:
    import sage.all
except ImportError:
    from random import shuffle, randint
else:
    from sage.misc.prandom import shuffle, randint
    from sage.arith.functions import lcm

def argmin(l):
    r"""
    Return the position of the minimal element in the list ``l``.

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import argmin
        sage: argmin([3,0,1,2])
        1
        sage: argmin([-1,3,5,-2,50])
        3
    """
    if not(l):
        raise ValueError('empty list')
    imin = 0
    jmin = l[0]
    for i,j in enumerate(l):
        if j < jmin:
            jmin = j
            imin = i
    return imin

#####################################################################
# Initialization
#####################################################################

def permutation_to_perm(p):
    r"""
    Return a list on `[0, n-1]` from a permutation on `[1, n]`

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import permutation_to_perm
        sage: permutation_to_perm(PermutationGroupElement([3,1,2]))
        [2, 0, 1]
    """
    return list(map(lambda x: x-1, p.domain()))

def perm_to_permutation(l):
    r"""
    Returns a permutation on `[1, n]` from a list on `[0, n-1]`

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_to_permutation
        sage: perm_to_permutation([2,1,0])
        (1,3)
    """
    from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
    return PermutationGroupElement(list(map(lambda x: x+1, l)))


def perm_init(data, n=None, partial=False):
    """
    Returns a permutation from ``data``.

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_init
        sage: perm_init([3,2,1,0])
        [3, 2, 1, 0]
        sage: perm_init([3,2,1,0], 5)
        [3, 2, 1, 0, 4]
        sage: perm_init(([2,1],[3,4,0]))
        [3, 2, 1, 4, 0]
        sage: perm_init('(0,1)(3,2)')
        [1, 0, 3, 2]
        sage: perm_init([3,1,-1,0])
        [3, 1, -1, 0]
        sage: perm_init([[2,1],[3,4,0]])
        [3, 2, 1, 4, 0]

        sage: perm_init([0r])
        [0]

    Zerology::

        sage: perm_init([])
        []
        sage: perm_init([], 4)
        [0, 1, 2, 3]
        sage: perm_init([[]])
        []
        sage: perm_init([[]], 4)
        [0, 1, 2, 3]
    """
    if isinstance(data, (tuple,list)):
        if not data:
            return list(range(n)) if n is not None else []
        if isinstance(data[0], (tuple,list)):
            return cycles_to_list(data, n, partial)
        else:
            p = [int(x) for x in data]
            if n is not None:
                p.extend(range(len(p),n))
            return p

    if isinstance(data, str):
        c = str_to_cycles(data)
        return cycles_to_list(c, n, partial)

    raise TypeError("The input must be list, tuple or string")

def init_perm(l):
    from warnings import warn
    warn('[surface_dynamics] init_perm is deprecated use perm_init instead', Warning, 3)
    return perm_init(l)

def equalize_perms(l):
    """
    Extend the permutations in ``l`` to have the same lengths.

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import equalize_perms
        sage: li = [[0, 1], [2, 1, 0]]
        sage: equalize_perms(li)
        3
    """
    n = max(map(len, l))
    for p in l:
        p.extend(range(len(p), n))
    return n


def perm_check(l, n=None):
    r"""
    Checks that ``l`` is a partial permutation of `\{0, 1, ..., n-1\}`.

    INPUT:

    - ``n`` - integer (optional)

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_check

    Good permutations::

        sage: perm_check([1, 0, 3, 2], 4)
        True
        sage: perm_check([-1], 1)
        True
        sage: perm_check([-1, 3, -1, 1], 4)
        True

        sage: perm_check([1,0,-1,-1,-1], 2)
        True

        sage: perm_check([0,-1])
        True
        sage: perm_check([-1,1])
        True

    Bad permutations::

        sage: perm_check([1, 0, 3, 2], 3)
        False
        sage: perm_check([2, 0])
        False
        sage: perm_check([1, 0, 1])
        False

        sage: perm_check([1,-1])
        False
        sage: perm_check([-1,0])
        False
    """
    if not isinstance(l, list):
        return False
    if n is None:
        n = len(l)
    else:
        n = int(n)

    im_seen = [False] * n
    ra_seen = [False] * n
    for i in range(n):
        if l[i] == -1:
            continue
        if not isinstance(l[i], numbers.Integral):
            raise TypeError("entries must be integers >= -1, not {}".format(type(l[i]).__name__))
        if l[i] < 0 or l[i] >= n or im_seen[l[i]]:
            return False
        ra_seen[i] = True
        im_seen[l[i]] = True
    return ra_seen == im_seen

def perm_is_one(l, n=None):
    r"""
    Test whether ``l`` is the identity on its domain.

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_is_one
        sage: perm_is_one([])
        True
        sage: perm_is_one([-1])
        True
        sage: perm_is_one([0])
        True
        sage: perm_is_one([1, 0])
        False

        sage: perm_is_one([0, -1])
        True
        sage: perm_is_one([-1, 1])
        True

        sage: perm_is_one([0, 2, 1], 1)
        True
    """
    if n is None:
        n = len(l)
    for i in range(n):
        if l[i] != -1 and l[i] != i:
            return False
    return True

#####################################################################
# Conversion
#####################################################################

def cycles_to_list(t, n=None, partial=False):
    r"""
    Returns a permutation on `[0, n-1]` from a list of cycles on `[0, n-1]`

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import cycles_to_list

        sage: cycles_to_list([[1,3,5]])
        [0, 3, 2, 5, 4, 1]
        sage: cycles_to_list([[1,3,5]], partial=True)
        [-1, 3, -1, 5, -1, 1]

        sage: cycles_to_list([[1,3,5],[0,2,4],[6]])
        [2, 3, 4, 5, 0, 1, 6]

        sage: cycles_to_list([])
        []
        sage: cycles_to_list([[],[]])
        []

        sage: cycles_to_list([[1,4,2,5],[3],[6]], partial=True)
        [-1, 4, 5, 3, 2, 1, 6]

        sage: cycles_to_list([[0,5]], 3)
        Traceback (most recent call last):
        ...
        ValueError: cycle value out of range
    """
    if not any(tt for tt in t):
        return list(range(n)) if n is not None else []

    if n is None:
        n = max(map(max, t)) + 1
    if partial:
        res = [-1] * n
    else:
        res = list(range(n))

    for c in t:
        for j in range(len(c)-1):
            if not isinstance(c[j], numbers.Integral) or c[j] < 0 or c[j] >= n:
                raise ValueError("cycle values out of range")
            res[c[j]] = int(c[j+1])
        if not isinstance(c[-1], numbers.Integral) or c[-1] < 0 or c[-1] >= n:
            raise ValueError("cycle value out of range")
        res[c[-1]] = int(c[0])

    return res


def str_to_cycles(s):
    """
    Returns a list of cycles from a string.

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
        r.append(list(map(int, c_str.replace(' ', '').split(','))))
    return r


def perm_cycles(p, singletons=False, n=None):
    r"""
    Return the cycle decomposition of ``p``

    INPUT:

    - ``p`` -- the permutation

    - ``singletons`` -- bool (default: ``False``) - return or not the singletons

    - ``n`` - (optional) only use the first ``n`` elements of the permutation ``p``

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_cycles, perm_init
        sage: perm_cycles([0,2,1])
        [[1, 2]]
        sage: perm_cycles([0,2,1],True)
        [[0], [1, 2]]

        sage: perm_cycles([2,-1,0])
        [[0, 2]]

        sage: perm_cycles([2,0,1,-1,-1], n=3)
        [[0, 2, 1]]

        sage: perm_cycles([-1,3,-1,1], singletons=True)
        [[1, 3]]

        sage: perm_cycles(perm_init('(1,3)', partial=True))
        [[1, 3]]
    """
    if n is None:
        n = len(p)
    seen = [1] * n
    res = []

    for i in range(n):
        if seen[i] and p[i] != -1:
            cycle = []
            j = i
            while seen[j]:
                seen[j] = 0
                cycle.append(j)
                j = p[j]
            if singletons or len(cycle) > 1:
                res.append(cycle)

    return res

def perm_cycle_tuples(*args, **kwds):
    print("WARNING: perm_cycle_tuples is deprecated, use perm_cycles")
    return perm_cycles(*args, **kwds)

def perm_num_cycles(p, n=None):
    r"""
    Return the number of cycles of the permutation ``p``.

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_num_cycles
        sage: perm_num_cycles([1,2,3,0])
        1
        sage: perm_num_cycles([0,2,3,1])
        2
        sage: perm_num_cycles([3,2,1,0])
        2
        sage: perm_num_cycles([3,1,2,0])
        3
        sage: perm_num_cycles([0,1,2,3])
        4
    """
    if n is None:
        n = len(p)
    seen = [False] * n
    ans = 0
    for i in range(n):
        if seen[i] or p[i] == -1:
            continue
        ans += 1
        j = i
        while not seen[j]:
            seen[j] = True
            j = p[j]
    return ans

def perm_cycle_type(p, n=None):
    r"""
    Return the lengths of the cycles of the permutation ``p`` in size of
    decreasing order.

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_cycle_type
        sage: perm_cycle_type([1,2,3,0])
        [4]
        sage: perm_cycle_type([0,2,3,1])
        [3, 1]
        sage: perm_cycle_type([3,2,1,0])
        [2, 2]
        sage: perm_cycle_type([3,1,2,0])
        [2, 1, 1]
        sage: perm_cycle_type([0,1,2,3])
        [1, 1, 1, 1]
    """
    if n is None:
        n = len(p)
    seen = [False] * n
    c = []
    for i in range(n):
        if seen[i] or p[i] == -1:
            continue
        k = 0
        j = i
        while not seen[j]:
            seen[j] = True
            k += 1
            j = p[j]
        c.append(k)
    c.sort(reverse=True)
    return c

def perm_dense_cycles(p, n=None):
    r"""
    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_dense_cycles

        sage: perm_dense_cycles([1,2,0])
        ([0, 0, 0], [3])

        sage: perm_dense_cycles([0,2,1])
        ([0, 1, 1], [1, 2])

        sage: perm_dense_cycles([2,1,0])
        ([0, 1, 0], [2, 1])
    """
    if n is None:
        n = len(p)
    deg = []
    res = [-1] * n
    k = 0
    for i in range(n):
        if res[i] != -1:
            continue
        d = 0
        while res[i] == -1:
            res[i] = k
            i = p[i]
            d += 1
        k += 1
        deg.append(d)
    return res, deg


def perm_dense_cycles_and_angles(p, n=None):
    if n is None:
        n = len(p)

    lab = [-1] * n  # labels
    ang = [-1] * n  # angle
    deg = []          # degrees
    k = 0
    for i in range(n):
        if p[i] == -1 or lab[i] != -1:
            continue
        d = 0
        while lab[i] == -1:
            lab[i] = k
            ang[i] = d
            d += 1
            i = p[i]
        k += 1
        deg.append(d)
    return lab, ang, deg

def perm_cycle_string(p, singletons=False, n=None):
    r"""
    Returns a string representing the cycle decomposition of `p`

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_cycle_string
        sage: perm_cycle_string([0,2,1])
        '(1,2)'
        sage: perm_cycle_string([0,2,1], True)
        '(0)(1,2)'
        sage: perm_cycle_string([0,1,2,-1], False, 3)
        '()'
    """
    c = perm_cycles(p, singletons, n)
    if not c:
        return '()'
    else:
        return ''.join(map(lambda x: '('+','.join(map(str, x))+')', c))

def _canonical_reg_perm(n, k):
    r"""
    Return a standard product of disjoint k-cycles on {0, 1, ..., n-1}

    TESTS::

        sage: from surface_dynamics.misc.permutation import _canonical_reg_perm
        sage: _canonical_reg_perm(6, 2)
        [1, 0, 3, 2, 5, 4]
        sage: _canonical_reg_perm(6, 3)
        [1, 2, 0, 4, 5, 3]
    """
    if not isinstance(n, numbers.Integral) or not isinstance(k, numbers.Integral):
        raise ValueError
    if n <= 0 or k <= 0 or n % k:
        raise ValueError

    p = []
    for i in range(0,n,k):
        p.extend(range(i+1, i+k))
        p.append(i)
    return p


def constellation_init(vertices, edges, faces, n=None, domain=None, check=True):
    r"""
    Each of ``vertices``, ``edges`` or ``faces can be ``None``, an
    integer or an object to initialize a (partial) permutation.

    INPUT:

    - ``vertices``, ``edges``, ``faces`` - permutations given as strings or lists

    - ``n`` - (optional) number of darts

    - ``check`` - boolean default ``True``)

    TESTS::

        sage: from surface_dynamics.misc.permutation import constellation_init

        sage: constellation_init([0,1], [1,0], [1,0])
        [[0, 1], [1, 0], [1, 0]]
        sage: constellation_init([0r,1r],[1r,0r], [1r,0r])
        [[0, 1], [1, 0], [1, 0]]

        sage: constellation_init([2,1,0], [1,2,0], None)
        [[2, 1, 0], [1, 2, 0], [0, 2, 1]]

        sage: constellation_init(3, 2, None, 6)
        [[1, 2, 0, 4, 5, 3], [1, 0, 3, 2, 5, 4], [0, 2, 5, 1, 4, 3]]

        sage: constellation_init('(0,1)', '(0,1)', '()')
        [[1, 0], [1, 0], [0, 1]]

        sage: constellation_init('(0,2,3,6)(1,4,5,7)', None, None)
        [[2, 4, 3, 6, 5, 7, 0, 1], [1, 0, 3, 2, 5, 4, 7, 6], [7, 6, 2, 0, 4, 1, 5, 3]]

    Each number (including fixed point) should be specified otherwise they are just
    ignored::

        sage: constellation_init('(1,6,5)(2,3,4)', '(1,2)(3,4)(5,6)', '(1,4,2,5)(3)(6)')
        [[-1, 6, 3, 4, 2, 1, 5], [-1, 2, 1, 4, 3, 6, 5], [-1, 4, 5, 3, 2, 1, 6]]

    TESTS::

        sage: constellation_init(None, '(0,1)(2,3)(4,5)', '(0,2,4)(1)(3,5)')
        [[5, 0, 1, 4, 3, 2], [1, 0, 3, 2, 5, 4], [2, 1, 4, 5, 0, 3]]
    """
    nones = [p is None for p in [vertices, edges, faces]]
    if sum(nones) > 1:
        if edges is None:
            edges = 2
            nones[1] = False
        else:
            raise ValueError("at least two of the three data should be provided")

    nums = [isinstance(p, numbers.Integral) for p in [vertices, edges, faces]]
    if sum(nums) == 3 and n is None:
        raise ValueError("numbers not enough to construct a constellation")

    if not (nones[0] or nums[0]):
        vertices = perm_init(vertices, n, partial=True)
        if check and not perm_check(vertices):
            raise ValueError("invalid vertex permutation")

    if not (nones[1] or nums[1]):
        edges = perm_init(edges, n, partial=True)
        if check and not perm_check(edges):
            raise ValueError("invalid edge permutation")

    if not (nones[2] or nums[2]):
        faces = perm_init(faces, n, partial=True)
        if check and not perm_check(faces):
            raise ValueError("invalid face permutation")

    P = [vertices, edges, faces]
    perms = [p for i,p in enumerate(P) if not (nones[i] or nums[i])]
    if perms:
        equalize_perms(perms)
        j = 0
        for i in range(3):
            if not (nones[i] or nums[i]):
                P[i] = perms[j]
                j += 1
        if n is None:
            n = len(perms[0])

    for i in range(3):
        if nums[i]:
            P[i] = _canonical_reg_perm(n, P[i])

    for i in range(3):
        if nones[i]:
            j = (i - 1) % 3
            k = (i + 1) % 3
            P[i] = perm_compose(perm_invert(P[j]), perm_invert(P[k]))

    if check:
        if set(P[0]) != set(P[1]) or set(P[0]) != set(P[2]):
            raise ValueError("permutations defined on different domains")
        if not perm_is_one(perm_compose(perm_compose(P[0], P[1]), P[2])):
            raise ValueError("product is not identity")

    return P

#####################################################################
# Group operations
#####################################################################

def perm_invert(l, n=None):
    r"""
    Returns the inverse of the permutation ``l``.

    TESTS::

        sage: from itertools import permutations
        sage: from surface_dynamics.misc.permutation import perm_init, perm_invert, perm_compose, perm_is_one
        sage: all(perm_is_one(perm_compose(perm_invert(p),p)) for p in permutations(range(3)))
        True
        sage: all(perm_is_one(perm_compose(p,perm_invert(p))) for p in permutations(range(3)))
        True

        sage: perm_invert([2, -1, 5, 0, -1, 3])
        [3, -1, 0, 5, -1, 2]

        sage: p = perm_init('(1,3,5)(2,9)', partial=True)
        sage: q = perm_invert(p)
        sage: q
        [-1, 5, 9, 1, -1, 3, -1, -1, -1, 2]
        sage: perm_is_one(perm_compose(p, q))
        True

        sage: perm_invert([2,1,0,-1,-1], 3)
        [2, 1, 0]
    """
    if n is None:
        n = len(l)
    res = [-1] * n
    for i in range(n):
        if l[i] != -1:
            res[l[i]] = i
    return res

def perm_invert_inplace(p, n=None):
    r"""
    Inverse in place the permutation p
    """
    if n is None:
        n = len(l)
    seen = [0] * n
    for i in range(n):
        if seen[i]:
            continue
        j = p[i]
        while not seen[i]:
            seen[i] = 1
            i, p[j] = p[j], i

def perm_compose(p1, p2):
    r"""
    Returns the product ``p1 p2``.

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_init, perm_compose
        sage: perm_compose([0,2,1],[0,2,1])
        [0, 1, 2]
        sage: perm_compose([-1,2,3,1],[-1,2,1,3])
        [-1, 1, 3, 2]

        sage: p1 = perm_init('(1,3,5)', partial=True)
        sage: p2 = perm_init('(1,5)(3)', partial=True)
        sage: perm_compose(p1, p2)
        [-1, 3, -1, 1, -1, 5]
    """
    r = [-1] * len(p1)
    for i in range(len(p1)):
        if p1[i] != -1 and p1[i] < len(p2):
            r[i] = p2[p1[i]]
    return r


def perm_compose_i(p1, p2):
    r"""
    Returns the product `p1^{-1} p2^{-1}`.

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_compose_i
        sage: perm_compose_i([0,1,2],[1,2,0])
        [2, 0, 1]

        sage: from surface_dynamics.misc.permutation import perm_invert, perm_compose
        sage: from itertools import permutations
        sage: for p1 in permutations(range(4)):
        ....:     for p2 in permutations(range(4)):
        ....:         assert perm_compose_i(p1, p2) == perm_compose(perm_invert(p1), perm_invert(p2))
    """
    assert(len(p1) == len(p2))

    res = [None]*len(p1)
    for i in range(len(p1)):
        res[p1[p2[i]]] = i

    return res

# can we do that inplace?
def perm_conjugate(p1, p2, n=None):
    r"""
    Conjugate ``p1`` by ``p2``.

    Let call ``res`` the output of this function. If ``p1`` was
    mapping ``a`` to ``b`` then ``res`` will map ``p2[a]``
    to ``p2[b]``.

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_random, perm_conjugate

        sage: p1 = perm_random(23)
        sage: p2 = perm_random(23)
        sage: res = perm_conjugate(p1, p2)
        sage: res[p2[14]] == p2[p1[14]]
        True
        sage: res[p2[19]] == p2[p1[19]]
        True
    """
    if n is None:
        n = len(p1)
    res = [-1] * n
    for i in range(n):
        res[p2[i]] = p2[p1[i]]
    return res

def perm_conjugate_inplace(p1, p2, n=None):
    r"""
    we want

    p1[p2[i]] <- p2[p1[i]]

    save tmp = p1[p2[0]]
    set  p1[p2[0]] = p2[p1[0]]

    """
    if n is None:
        n = len(p1)
    unseen = [True] * n
    for i in range(n):
        if not unseen[i]:
            continue
        unseen[i] = True

def perm_on_list_inplace(p, a, n=None, swap=None):
    r"""
    Inplace action of permutation on list-like objects.

    INPUT:

    - ``p`` - permutation

    - ``a`` - list, array

    - ``n`` - (optional) size of permutation

    - ``swap`` - (optional) a swap function

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import *
        sage: l = [0,1,2,3,4]
        sage: p = [4,2,3,0,1]
        sage: perm_on_list_inplace(p,l)
        sage: l
        [3, 4, 1, 2, 0]

    Permutation action on matrix rows::

        sage: m1 = matrix(ZZ, 5, 5, 1)
        sage: m2 = matrix(ZZ, 5, 5, 1)
        sage: m = matrix(ZZ, 5, 5, 1)
        sage: p1 = perm_init([4,1,3,2,0])
        sage: p2 = perm_init([1,0,3,4,2])
        sage: perm_on_list_inplace(p1, m1, swap=sage.matrix.matrix0.Matrix.swap_rows)
        sage: perm_on_list_inplace(p2, m2, swap=sage.matrix.matrix0.Matrix.swap_rows)
        sage: perm_on_list_inplace(perm_compose(p1, p2), m, swap=sage.matrix.matrix0.Matrix.swap_rows)
        sage: m == m2 * m1
        True
    """
    if n is None:
        n = len(p)
    seen = [False] * n
    for i in range(n):
        if seen[i]:
            continue
        seen[i] = True
        j = p[i]
        while seen[j] is False:
            if swap:
                swap(a, i, j)
            else:
                tmp = a[i]
                a[i] = a[j]
                a[j] = tmp
            seen[j] = True
            j = p[j]

################################################################
# Various permutation constructors (including randomized ones) #
################################################################

def perm_one(n):
    r"""
    The identity permutation

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_one
        sage: perm_one(0)
        []
        sage: perm_one(3)
        [0, 1, 2]
    """
    return list(range(n))

def perm_random(n):
    r"""
    Return a random permutation.

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_random, perm_check
        sage: perm_check(perm_random(13), 13)
        True
    """
    r = list(range(n))
    shuffle(r)
    return r

def perm_random_centralizer(p):
    r"""
    Return a random permutation in the centralizer of ``p``.

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_random, perm_compose, perm_random_centralizer
        sage: p = perm_random(10)
        sage: q = perm_random_centralizer(p)
        sage: perm_compose(p, q) == perm_compose(q, p)
        True
    """
    if not p:
        return p

    cyc = perm_cycles(p)
    cyc.sort(key = lambda c: len(c))
    i = 0
    ans = [-1] * len(p)
    while i < len(cyc):
        j = i + 1
        s = len(cyc[i])
        while j < len(cyc) and len(cyc[j]) == s:
            j += 1

        # permutation of the cycles i, i+1, ..., j-1
        m = perm_random(j - i)

        for ii in range(i, j):
            jj = i + m[ii - i]
            shift = randint(0, s - 1)  # random shift
            for k in range(len(cyc[i])):
                ans[cyc[ii][k]] = cyc[jj][(k + shift) % s]

        # next round
        i = j

    return ans

def perm_random_conjugacy_class(c):
    r"""
    Return a random permutation with given conjugacy class ``c``.

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import (perm_random_conjugacy_class,
        ....:     perm_cycle_type)

        sage: p = perm_random_conjugacy_class([5,2])
        sage: perm_cycle_type(p)
        [5, 2]

        sage: p = perm_random_conjugacy_class([7,3,3,1])
        sage: perm_cycle_type(p)
        [7, 3, 3, 1]
    """
    n = sum(c)
    r = list(range(n))
    shuffle(r)
    p = [-1] * n
    i = 0
    for k in c:
        # add a k-cycle following the list r
        for j in range(i, i+k-1):
            p[r[j]] = r[j+1]
        p[r[i+k-1]] = r[i]
        i += k
    return p

#####################################################################
# Actions
#####################################################################

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


def perm_on_list(p, t):
    r"""
    Action of the permutation ``p`` on the list ``t``.

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_on_list
        sage: perm_on_list([2,1,3,0], [2,1,2,0])
        [3, 1, 3, 2]
    """
    return [p[i] for i in t]


def perm_on_cyclic_list(p, t):
    r"""
    Action of the permutation ``p`` on the list ``t`` up to cyclic order.

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_on_cyclic_list
        sage: perm_on_cyclic_list([0,1,2], [2,1,2])
        [1, 2, 2]
        sage: perm_on_cyclic_list([0,1], [0,1,0,0,1,0,0,0,1,1])
        [0, 0, 0, 1, 1, 0, 1, 0, 0, 1]

        sage: a = [1, 0, 3, 2, 5, 4]
        sage: perm_on_cyclic_list(a, [0, 5, 3])
        [1, 4, 2]
        sage: perm_on_cyclic_list(a, [1, 4, 2])
        [0, 5, 3]

        sage: a1 = [0, 1, 4, 2, 3, 5]
        sage: a2 = [1, 5, 3, 4, 2, 0]
        sage: a3 = [2, 3, 1, 5, 0, 4]
        sage: a4 = [5, 0, 2, 3, 4, 1]
        sage: t1 = [0, 5, 1]
        sage: t2 = [2, 4, 3]
        sage: perm_on_cyclic_list(a1, t1) == perm_on_cyclic_list(a2, t1) == perm_on_cyclic_list(a4, t1) == t1
        True
        sage: perm_on_cyclic_list(a3, t1) == t2
        True
        sage: perm_on_cyclic_list(a3, t2) == t1
        True
    """
    res = r = [p[i] for i in t]
    for i in range(1, len(r)):
        rr = r[i:] + r[:i]
        if rr < res:
            res = rr
    return res


def canonical_perm(part, i=0):
    r"""
    Return the canonical permutation with the given part.

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
        res.extend(range(i+1,i+p))
        res.append(i)
        i += p
    return res


def canonical_perm_i(part, i=0):
    r"""
    Return the canonical permutation reversed.

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
        res.extend(range(i,i+p-1))
        i += p
    return res


def perm_switch(p1, p2, i, j):
    """
    Exchanges the values at positions ``i`` and ``j`` in two permutations
    ``p_1`` and ``p_2``.

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
    Relabel the list of permutations ``p`` according to ``m``.

    INPUT:

    - ``p`` - a list of permutations

    - ``m`` - the relabeling permutation

    EXAMPLES::

        sage: from surface_dynamics.misc.constellation import perms_relabel
        sage: perms_relabel([[0,1,2],[0,2,1]],[2,1,0])
        [[0, 1, 2], [1, 0, 2]]
    """
    q = [k[:] for k in p]
    for i in range(len(m)):
        for j in range(len(p)):
            q[j][m[i]] = m[p[j][i]]
    return q


def perms_canonical_labels_from(x, y, j0):
    r"""
    Return canonical labels for ``x``, ``y`` that starts at ``j0``.

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
    waiting = [[] for i in range(len(y))]

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
    m_win = list(range(n))

    x = p[0]
    y = p[1:]

    if e is None:
        e = list(range(n))

    # get canonical label from i in to_test and compare
    while e:
        i = e.pop()
        m_test = perms_canonical_labels_from(x, y, i)
        c_test = perms_relabel(p, m_test)
        if c_win is None or c_test < c_win:
            c_win = c_test
            m_win = m_test

    return c_win, m_win


#######################################################################
# Dynamical permutation groups
#######################################################################
def perm_order(p, n=None):
    r"""
    Return the multiplicative order of the permutation ``p``.

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_init, perm_order
        sage: p = perm_init('(1,3)(2,4,6)(5)')
        sage: perm_order(p)
        6
    """
    return lcm(perm_cycle_type(p))

# TODO: redesign
# when adding a generator we need to check the orbit under what was already computed
class PermutationGroupOrbit(object):
    r"""
    Dynamical orbit generation.

    This is to be used when computing the automorphism group of an object (e.g.
    a ribbon graph). It allows to iterate through orbit representatives
    dynamically (i.e. group generators can be added on the fly).

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import PermutationGroupOrbit
        sage: P1 = PermutationGroupOrbit(7, ['(0,1)', '(0,1,2,3,4,5,6)'])
        sage: P2 = PermutationGroupOrbit(4, ['(0,1)', '(2,3)'])
    """
    __slots__ = ['_n', '_gens', '_S', '_s', '_seen', '_elts']

    def __init__(self, n, gens, S=None, check=True):
        if check:
            gens = [perm_init(g, n=n) for g in gens]
        if S is None:
            S = list(range(n))

        self._n = n              # underlying set is {0, 1, ..., n-1}
        self._gens = gens        # group generators
        self._S = S              # the set we are interested in
        self._s = 0              # index in the set S
        self._seen = [False] * n # the elements that we already visited (dense version)
        self._elts = []          # the elements that we already visited (sparse version)

    def __iter__(self):
        return self

    def libgap_group(self):
        from sage.libs.gap.libgap import libgap
        if not self._gens:
            return libgap.Group([libgap.PermList([])])
        gens = []
        for a in self._gens:
            gens.append(libgap.PermList([a[i]+1 for i in range(self._n)]))
        return libgap.Group(gens)

    def group_cardinality(self):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.misc.permutation import PermutationGroupOrbit
            sage: PermutationGroupOrbit(7, []).group_cardinality()
            1
            sage: PermutationGroupOrbit(5, ['(0,2)(3,4)']).group_cardinality()
            2
            sage: PermutationGroupOrbit(7, ['(0,1)', '(0,1,2,3,4,5,6)']).group_cardinality()
            5040
        """
        if not self._gens:
            return 1
        elif len(self._gens) == 1:
            return perm_order(self._gens[0])
        return int(self.libgap_group().Size())

    def __repr__(self):
        n = self._n
        return "PermutationGroupOrbit({}, [{}])".format(
                n,
                ', '.join(perm_cycle_string(g, singletons=False, n=n) for g in self._gens))

    def num_gens(self):
        return len(self._gens)

    def gens(self):
        return self._gens

    def reset_iterator(self, S=None):
        r"""
        Reset the self we iterate over with.

        EXAMPLES::

            sage: from surface_dynamics.misc.permutation import PermutationGroupOrbit

            sage: O = PermutationGroupOrbit(4, [])
            sage: list(O)
            [0, 1, 2, 3]
            sage: list(O)
            []
            sage: O.reset_iterator()
            sage: list(O)
            [0, 1, 2, 3]
            sage: list(O)
            []

            sage: O.add_generator('(0,1)(2,3)')
            sage: list(O)
            []
            sage: O.reset_iterator()
            sage: list(O)
            [0, 2]
        """
        n = self._n
        if S is None:
            S = list(range(n))
        self._S = S
        self._seen = [False] * n
        del self._elts[:]
        self._s = 0

    def __next__(self):
        S = self._S        # candidates
        s = self._s        # current index in S
        if s == len(S):
            raise StopIteration

        seen = self._seen  # already visited elements (dense)
        i = S[s]
        s += 1
        while seen[i] and s < len(S):
            i = S[s]
            s += 1

        # update s and remove the orbit of i
        self._s = s
        if seen[i]:
            raise StopIteration

        elts = self._elts  # already visited elements (sparse)
        k = len(elts)
        elts.append(i)
        while k < len(elts):
            u = elts[k]
            for g in self._gens:
                v = g[u]
                if not seen[v]:
                    seen[v] = True
                    elts.append(v)
            k += 1

        return i

    next = __next__  # Python2 support

    def add_generator(self, g, check=True):
        r"""
        EXAMPLES::

            sage: from surface_dynamics.misc.permutation import PermutationGroupOrbit

            sage: O = PermutationGroupOrbit(4, [])
            sage: next(O)
            0
            sage: O.add_generator('(0,1)(2,3)')
            sage: next(O)
            2
            sage: next(O)
            Traceback (most recent call last):
            ...
            StopIteration
        """
        if check:
            g = perm_init(g, n=self._n)
        self._gens.append(g)
        elts = self._elts
        seen = self._seen

        k = 0
        K = len(elts)
        while k < K:
            u = elts[k]
            v = g[u]
            if not seen[v]:
                seen[v] = True
                elts.append(v)
            k += 1

        while k < len(elts):
            u = elts[k]
            for g in self._gens:
                v = g[u]
                if not seen[v]:
                    seen[v] = True
                    elts.append(v)
            k += 1

#################################3
# Serialization

# the first 64 characters are used for integer encoding
# the last one is used for "undefined" (corresponds to -1 in arrays)
CHARS = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ+-.'
CHARS_INV = {j:i for i,j in enumerate(CHARS)}

def uint_base64_str(n, l=None):
    r"""
    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import uint_base64_str

        sage: uint_base64_str(15)
        'f'
        sage: uint_base64_str(15, 3)
        '00f'

        sage: uint_base64_str(-1)
        '.'
        sage: uint_base64_str(-1, 3)
        '...'
    """
    n = int(n)
    if n == -1:
        if l is None:
            return CHARS[64]
        else:
            return CHARS[64] * l
    elif n < -1:
        raise ValueError("invalid negative integer")
    s = ''
    while n:
        s = CHARS[n % 64] + s
        n //= 64
    if l is not None:
        if len(s) > l:
            raise ValueError
        else:
            s = CHARS[0] * (l - len(s)) + s
    return s

def uint_from_base64_str(s):
    r"""
    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import uint_from_base64_str, uint_base64_str

        sage: uint_from_base64_str('mqb')
        91787
        sage: uint_base64_str(91787)
        'mqb'

        sage: uint_from_base64_str('00mqb')
        91787

        sage: uint_from_base64_str('...')
        -1
    """
    if not s or not isinstance(s, str):
        raise ValueError
    if s[0] == CHARS[64]:
        assert all(i == CHARS[64] for i in s)
        return -1
    n = 0
    d = 1
    for c in reversed(s):
        n += CHARS_INV[c] * d
        d *= 64
    return n

def perm_base64_str(p, n=None):
    r"""
    Make a canonical ASCII string out of ``p``.

    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_base64_str, perm_from_base64_str
        sage: from array import array

        sage: perm_base64_str([])
        ''
        sage: perm_base64_str([3,1,0,2])
        '3102'
        sage: s = perm_base64_str(range(2048))
        sage: s
        '00010203...v+v-'
        sage: perm_from_base64_str(s, 2048) == list(range(2048))
        True
    """
    if n is None:
        n = len(p)
    if not n:
        return ''
    l = int(log(n, 64)) + 1 # number of digits used for each entry
    return ''.join(uint_base64_str(p[i], l) for i in range(n))

def perm_from_base64_str(s, n):
    r"""
    EXAMPLES::

        sage: from surface_dynamics.misc.permutation import perm_from_base64_str, perm_base64_str
        sage: from array import array

        sage: p = [3,0,2,1]
        sage: s = perm_base64_str(p)
        sage: perm_from_base64_str(s, 4) == p
        True

        sage: perm_base64_str([-1,-1])
        '..'
        sage: perm_from_base64_str('..', 2)
        [-1, -1]

        sage: p = list(range(3000))
        sage: shuffle(p)
        sage: perm_from_base64_str(perm_base64_str(p), 3000) == p
        True
        sage: p[18] = p[2003] = -1
        sage: perm_from_base64_str(perm_base64_str(p), 3000) == p
        True
    """
    if not n:
        if s:
            raise ValueError("invalid input")
        return []
    l = int(log(n, 64)) + 1 # number of digits used for each entry
    if len(s) != n * l:
        raise ValueError('wrong length')
    return [uint_from_base64_str(s[i:i+l]) for i in range(0,len(s),l)]
