r"""
Cardinality of Rauzy classes

The following functions implement algorithms relative to article [Del2010]_ and
[Boi2010]_ where are given formulas for the cardinality of Rauzy classes of
permutations.

  - ``c``: number of standard labeled permutations
  - ``d``: spin difference of standard labeled permutations
  - ``gamma_std``: number of standard permutations (with given profile and
    marking)
  - ``gamma_irr``: number of irreducible permutations (with given profile and
    marking)
  - ``delta_std``: spin difference for standard permutations (with given profile
    and marking)
  - ``delta_irr``: spin difference for irreducible permutations (with given
    profile and marking)

AUTHOR:

Vincent Delecroix

REFERENCES:

.. [Boi2010] Boissy 2010

.. [Del2010] Delecroix 2010

.. [Vee1982] W. Veech, "Gauss measures for transformations on the space of
   interval exchange maps", Ann. of Math., vol. 115, no. 2 (1982), pp. 201-242.

"""

from sage.misc.cachefunc import cached_function

from sage.rings.integer import Integer

from sage.combinat.partition import Partition
from sage.arith.all import factorial, binomial

#########################
# PROFILES AND MARKINGS #
#########################

def marking_iterator(profile,left=None,standard=False):
    r"""
    Returns the marked profile associated to a partition

    EXAMPLES::

        sage: import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rcc
        sage: p = Partition([3,2,2])
        sage: list(rcc.marking_iterator(p))
        [(1, 2, 0),
         (1, 2, 1),
         (1, 3, 0),
         (1, 3, 1),
         (1, 3, 2),
         (2, 2, 2),
         (2, 2, 3),
         (2, 3, 2)]
    """
    e = Partition(sorted(profile,reverse=True)).to_exp_dict()

    if left is not None:
        assert(left in e)

    if left is not None: keys = [left]
    else: keys = e.keys()

    for m in keys:
        if standard: angles = range(1,m-1)
        else: angles = range(0,m)
        for a in angles:
            yield (1,m,a)

    for m_l in keys:
        for m_r in e:
            if m_l != m_r or e[m_l] > 1:
                yield (2,m_l,m_r)

def split(p,k,i=0):
    r"""
    Splits the i-th term of p into two parts of size k and n-k-1

    There is a symmetry split(p, k, i) = split(p, p[i]-k-1, i)

    INPUT:

    - ``p`` - a partition

    - ``k`` - an integer between 2 and p[i]

    - ``i`` - integer - the index of the element to split

    OUTPUT: a partition

    EXAMPLES::

        sage: import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rcc
        sage: p = Partition([5,1])
        sage: rcc.split(p,1,0)
        [3, 1, 1]
        sage: rcc.split(p,2,0)
        [2, 2, 1]
        sage: rcc.split(p,3,0)
        [3, 1, 1]
    """
    l = list(p)
    n = l.pop(i)
    l.append(n-k-1)
    l.append(k)
    l.sort(reverse=True)
    return Partition(l)

def collapse(p,i,j):
    r"""
    Collapses the i-th term and the j-th term of a permutation

    INPUT:

    - ``p`` - a partition

    - ``i,j`` - two different indices of p

    OUTPUT:

    - a partition

    EXAMPLES::

        sage: import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rcc
        sage: p = Partition([4,2,1])
        sage: rcc.collapse(p, 0, 1)
        [5, 1]
        sage: rcc.collapse(p, 0, 2)
        [4, 2]
        sage: rcc.collapse(p, 1, 2)
        [4, 2]
    """
    assert(i != j)

    l = list(p)

    n_i = l[i]
    n_j = l[j]

    l[i] = n_i+n_j-1
    del l[j]

    l.sort(reverse=True)
    return Partition(l)

def check_std_marking(p, marking):
    r"""
    Tiny internal function that checks the validity of ``marking`` on the
    partition ``p``.

    EXAMPLES::

        sage: import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rcc
        sage: p = Partition([3,2,2])
        sage: rcc.check_std_marking(p, (1,3,1))
        (1, 3, 1)
        sage: rcc.check_std_marking(p, (1,3,0))
        Traceback (most recent call last):
        ...
        ValueError: marking[2] is not good
    """
    if len(marking) != 3:
        raise ValueError("marking must be a 3-tuple (t,i,j)")

    if marking[0] == 1:
        if marking[1] not in p:
            raise ValueError("marking[1] not in p")
        if marking[2] < 1 or marking[2] > marking[1]-2:
            raise ValueError("marking[2] is not good")

    elif marking[0] == 2:
        if marking[1] == marking[2]:
            if not p.to_exp(marking[1]) > 1:
                raise ValueError("wrong marking type 2")
        elif marking[1] not in p or marking[2] not in p:
            raise ValueError("marking not in p")
    else:
        raise ValueError("marking[0] must be 1 or 2")

    return tuple(marking)

def check_marking(p, marking):
    r"""
    Tiny internal function that checks that ``marking`` is compatible with ``p``.

    OUTPUT:

        A 3-tuple.

    EXAMPLES::

        sage: import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rcc
        sage: p = Partition([3,2,2])
        sage: rcc.check_marking(p, (1,3,1))
        (1, 3, 1)
        sage: rcc.check_marking(p, (1,3,2))
        (1, 3, 2)
        sage: rcc.check_marking(p, (1,3,3))
        Traceback (most recent call last):
        ...
        ValueError: marking[2] is not good
        sage: rcc.check_marking(p, (1,3,-1))
        Traceback (most recent call last):
        ...
        ValueError: marking[2] is not good
        sage: rcc.check_marking(p, (2,3,2))
        (2, 3, 2)
        sage: rcc.check_marking(p, (2,3,3))
        Traceback (most recent call last):
        ...
        ValueError: wrong marking type 2
    """
    if len(marking) != 3:
        raise ValueError("marking must be a 3-tuple (t,i,j)")

    if marking[0] == 1:
        if marking[1] not in p:
            raise ValueError("marking[1] not in p")
        if marking[2] < 0 or marking[2] > marking[1]-1:
            raise ValueError("marking[2] is not good")

    elif marking[0] == 2:
        if marking[1] == marking[2]:
            if not p.to_exp_dict()[marking[1]] > 1:
                raise ValueError("wrong marking type 2")
        elif marking[1] not in p or marking[2] not in p:
            raise ValueError("marking not in p")
    else:
        raise ValueError("marking[0] must be 1 or 2")

    return tuple(marking)

def bidecompositions(p):
    r"""
    Iterator through the pair of partitions ``(q1,q2)`` such that the union of
    the parts of ``q1`` and ``q2`` equal ``p``.

    EXAMPLES::

        sage: import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rcc
        sage: list(rcc.bidecompositions(Partition([3,1])))
        [([], [3, 1]), ([3], [1]), ([1], [3]), ([3, 1], [])]
        sage: list(rcc.bidecompositions(Partition([2,1,1])))
        [([], [2, 1, 1]),
         ([2], [1, 1]),
         ([1], [2, 1]),
         ([2, 1], [1]),
         ([1, 1], [2]),
         ([2, 1, 1], [])]
    """
    from itertools import product

    exp = p.to_exp()
    for i in product(*tuple(xrange(i+1) for i in exp)):
        p1 = Partition(exp=i)
        p2 = Partition(exp=[exp[j]-i[j] for j in xrange(len(exp))])
        yield p1,p2

#########################
# STANDARD PERMUTATIONS #
#########################

# number of permutations

@cached_function
def _c_rec(p):
    r"""
    Recurrence function that is called by :func:`c`.
    """
    if p[0] == 1: return factorial(len(p)-1)

    return (
         sum(_c_rec(split(p,k)) for k in xrange(1,p[0]-1)) +
         sum(p[i]*_c_rec(collapse(p,0,i)) for i in xrange(1,len(p))))

def c(p):
    r"""
    Number of labeled standard permutations with given profile

    There is an explicit formula for this number

    .. MATH::

        c(p) = \frac{2 (n-1)!}{n+1} \left( \sum_{q \subset (p_2,p_3,\ldots,p_k)} (-1)^{s(q)-l(q)} \binom{n}{s(q)}^{-1} \right).

    Though, for huge partition `p` this is not very useful. This function
    implements an induction formula to compute `c(p)`.

    EXAMPLES::

        sage: import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rcc

    Partition of length 1::

        sage: n = 7
        sage: rcc.c([n]) == 2 * factorial(n-1) / (n+1)
        True
        sage: all(rcc.c([n]) == 2 * factorial(n-1) / (n+1) for n in xrange(11,18,2))
        True

    Partitions of length 2 with two odd numbers::

        sage: p = [5,3]
        sage: n = sum(p)
        sage: b = binomial(n,p[0])
        sage: rcc.c(p) == 2 * factorial(n-1) * (1 + 1 / binomial(n,p[0])) / (n+1)
        True

        sage: p = [13,5]
        sage: n = sum(p)
        sage: b = binomial(n,p[0])
        sage: rcc.c(p) == 2 * factorial(n-1) * (1 + 1 / binomial(n,p[0])) / (n+1)
        True

    Partitions of length 2 with even numbers::

        sage: p = [4,4]
        sage: n = sum(p)
        sage: b = binomial(n,p[0])
        sage: rcc.c(p) == 2 * factorial(n-1) * (1 - 1 / binomial(n,p[0])) / (n+1)
        True

        sage: p = [10,2]
        sage: n = sum(p)
        sage: b = binomial(n,p[0])
        sage: rcc.c(p) == 2 * factorial(n-1) * (1 - 1 / binomial(n,p[0])) / (n+1)
        True

    Add marked points to an integer partition::

        sage: p = [3,2,2]
        sage: n = sum(p)
        sage: all(rcc.c(p + [1]*k) == factorial(n+k-1) / factorial(n-1) * rcc.c(p) for k in xrange(1,6))
        True
    """
    return _c_rec(Partition(p))


def gamma_std(profile, marking=None):
    r"""
    Return the number of standard permutations of given profile

    INPUT:

    - ``profile`` - an integer partition such that the its sum plus its length
      is congruent to 0 modulo 2

    - ``marking`` - either None, an element of the profile or a 3-tuple ``(1, n1, a)`` or ``(2, n1, n2)``

    EXAMPLES::

        sage: import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rcc

    A ValueError is raised if the partition does not satisfy the requirement::

        sage: rcc.gamma_std([5,2])
        Traceback (most recent call last):
        ...
        ValueError: the sum of the profile (=[5, 2]) plus its length must be congruent to 0 modulo 2

    The Rauzy classes associated to connected strata in genus 3::

        sage: from surface_dynamics.all import AbelianStratum
        sage: cc = AbelianStratum(1,1,1,1).unique_component()
        sage: d = cc.rauzy_diagram()
        sage: d
        Rauzy diagram with 1255 permutations
        sage: len(filter(lambda x: x.is_standard(), d)) == rcc.gamma_std([2,2,2,2])
        True

        sage: cc = AbelianStratum(2,1,1).unique_component()
        sage: d = cc.rauzy_diagram()
        sage: d
        Rauzy diagram with 2177 permutations
        sage: len(filter(lambda x: x.is_standard(), d)) == rcc.gamma_std([3,2,2])
        True

        sage: cc = AbelianStratum(3,1).unique_component()
        sage: d = cc.rauzy_diagram()
        sage: d
        Rauzy diagram with 770 permutations
        sage: len(filter(lambda x: x.is_standard(), d)) == rcc.gamma_std([4,2])
        True

    The non connected strata in genus 3::

        sage: cc_odd = AbelianStratum(2,2).odd_component()
        sage: cc_hyp = AbelianStratum(2,2).hyperelliptic_component()
        sage: d_odd = cc_odd.rauzy_diagram()
        sage: d_hyp = cc_hyp.rauzy_diagram()
        sage: d_odd
        Rauzy diagram with 294 permutations
        sage: d_hyp
        Rauzy diagram with 63 permutations
        sage: n_odd = len(filter(lambda x: x.is_standard(), d_odd))
        sage: n_hyp = len(filter(lambda x: x.is_standard(), d_hyp))
        sage: n_odd + n_hyp == rcc.gamma_std([3,3])
        True

        sage: cc_odd = AbelianStratum(4).odd_component()
        sage: cc_hyp = AbelianStratum(4).hyperelliptic_component()
        sage: d_odd = cc_odd.rauzy_diagram()
        sage: d_hyp = cc_hyp.rauzy_diagram()
        sage: d_odd
        Rauzy diagram with 134 permutations
        sage: d_hyp
        Rauzy diagram with 31 permutations
        sage: n_odd = len(filter(lambda x: x.is_standard(), d_odd))
        sage: n_hyp = len(filter(lambda x: x.is_standard(), d_hyp))
        sage: n_odd + n_hyp == rcc.gamma_std([5])
        True
    """
    p = Partition(sorted(profile,reverse=True))
    if (sum(p) + len(p)) % 2 != 0:
        raise ValueError, "the sum of the profile (=%s) plus its length must be congruent to 0 modulo 2" %p

    if len(p) == 0 and marking==(1,0,0):
        return 1

    if marking is None:
        return sum(gamma_std(p,m) for m in marking_iterator(p,left=None,standard=True))

    elif isinstance(marking, (int,Integer)):
        return sum(gamma_std(p,m) for m in marking_iterator(p,left=marking,standard=True))

    marking = check_std_marking(p, marking)
    l = list(p)

    if marking[0] == 1: # marking of type 1
        i = l.index(marking[1])
        del l[i]
        return _c_rec(split(p,marking[2],i)) / Partition(l).centralizer_size()

    else: # marking of type 2
        i0 = l.index(marking[1])
        if marking[2] <= marking[1]:
            i1 = i0 + 1 + l[i0+1:].index(marking[2])
        else:
            i1 = l.index(marking[2])
            i0,i1 = i1,i0
        del l[i1]
        del l[i0]

        return _c_rec(collapse(p,i0,i1)) / Partition(l).centralizer_size()

number_of_standard_permutations = gamma_std

# spin difference

def d(p):
    r"""
    Difference between the number of odd spin parity and even spin parity
    standard labeled permutations with given profile

    There is an explicit formula

    .. MATH::

        d(p) = \frac{(n-1)!}{2^{(n-k)/2}}

    where `n` is the sum of the partition `p` and `k` is its length.

    EXAMPLES::

        sage: import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rcc
        sage: p = [3,3,1]
        sage: rcc.d([3,3,1]) == factorial(6) / 2**2
        True
        sage: rcc.d([13]) == factorial(12) / 2**6
        True

     Adding marked points::

        sage: p = [5,3,3]
        sage: n = sum(p)
        sage: all(rcc.d(p + [1]*k) == factorial(n+k-1) / factorial(n-1) * rcc.d(p) for k in xrange(1,6))
        True
    """
    n = sum(p)
    k = len(p)
    return factorial(n-1) / 2**((n-k)/2)


def delta_std(profile, marking=None):
    r"""
    Return the difference odd-even in the given stratum

    INPUT:

    - ``p`` - partition with odd terms

    - ``marking`` - a 3-tuple ``(1, n1, a)`` or ``(2, n1, n2)``

    EXAMPLES::

        sage: import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rcc


    A ValueError is raised if the partition does not fullfill the requirement::

        sage: rcc.delta_std([5,2])
        Traceback (most recent call last):
        ...
        ValueError: the profile (=[5, 2]) must contain only odd numbers

    Non connected strata in genus 3 has two connected components distinguished
    by their spin parity::

        sage: from surface_dynamics.all import AbelianStratum
        sage: cc_odd = AbelianStratum(2,2).odd_component()
        sage: cc_hyp = AbelianStratum(2,2).hyperelliptic_component()
        sage: d_odd = cc_odd.rauzy_diagram()
        sage: d_hyp = cc_hyp.rauzy_diagram()
        sage: d_odd
        Rauzy diagram with 294 permutations
        sage: d_hyp
        Rauzy diagram with 63 permutations
        sage: n_odd = len(filter(lambda x: x.is_standard(), d_odd))
        sage: n_hyp = len(filter(lambda x: x.is_standard(), d_hyp))
        sage: n_odd - n_hyp == rcc.delta_std([3,3])
        True

        sage: cc_odd = AbelianStratum(4).odd_component()
        sage: cc_hyp = AbelianStratum(4).hyperelliptic_component()
        sage: d_odd = cc_odd.rauzy_diagram()
        sage: d_hyp = cc_hyp.rauzy_diagram()
        sage: d_odd
        Rauzy diagram with 134 permutations
        sage: d_hyp
        Rauzy diagram with 31 permutations
        sage: n_odd = len(filter(lambda x: x.is_standard(), d_odd))
        sage: n_hyp = len(filter(lambda x: x.is_standard(), d_hyp))
        sage: n_odd - n_hyp == rcc.delta_std([5])
        True
    """
    p = Partition(sorted(profile,reverse=True))

    if any(not x%2 for x in p):
        raise ValueError, "the profile (=%s) must contain only odd numbers"%p

    if marking is None:
        return sum(delta_std(p,m) for m in marking_iterator(p,left=None,standard=True))

    elif isinstance(marking, (int,Integer)):
        return sum(delta_std(p,m) for m in marking_iterator(p,left=marking,standard=True))

    marking = check_std_marking(p, marking)
    l = list(p)

    if marking[0] == 1: # marking of type 1
        if marking[2]%2 == 0: return 0
        i = l.index(marking[1])
        del l[i]
        return d(split(p,marking[2],i)) / Partition(l).centralizer_size()

    else: # marking of type 2
        i0 = l.index(marking[1])
        if marking[2] <= marking[1]:
            i1 = i0 + 1 + l[i0+1:].index(marking[2])
        else:
            i1 = l.index(marking[2])
            i0,i1 = i1,i0
        del l[i1]
        del l[i0]
        return d(collapse(p,i0,i1)) / Partition(l).centralizer_size()

spin_difference_for_standard_permutations = delta_std

#####################################
# FROM STANDARD PERMUTATIONS TO ALL #
#####################################

@cached_function
def _gamma_irr_rec(p, marking):
    r"""
    Internal recursive function called by :func:`gamma_irr`
    """
    if len(p) == 0:
        return 1

    if marking[0] == 1:
        m = marking[1]
        a = marking[2]
        i = p.index(m)
        pp = Partition(p._list[:i]+p._list[i+1:]) # the partition p'

        N = gamma_std(pp._list + [m+2],(1,m+2,m-a))

        for m1 in xrange(1,m-1):
            m2 = m-m1-1
            for a1 in xrange(max(0,a-m2),min(a,m1)):
                a2 = a - a1 - 1
                for p1,p2 in bidecompositions(pp):
                    l1 = sorted([m1]+p1._list,reverse=True)
                    l2 = sorted([m2+2]+p2._list,reverse=True)
                    if (sum(l1)+len(l1)) % 2 == 0 and (sum(l2)+len(l2)) % 2 == 0:
                        N -= (_gamma_irr_rec(Partition(l1), (1,m1,a1)) *
                              gamma_std(Partition(l2),(1,m2+2,m2-a2)))
        return N


    elif marking[0] == 2:
        m1 = marking[1]
        m2 = marking[2]
        i1 = p.index(m1)
        i2 = p.index(m2)
        if m1 == m2: i2 += 1
        if i2 < i1: i1,i2 = i2,i1
        pp = Partition(p._list[:i1] + p._list[i1+1:i2] + p._list[i2+1:])

        N = gamma_std(pp._list + [m1+1,m2+1],(2,m1+1,m2+1))

        for p1,p2 in bidecompositions(pp):
            for k1 in xrange(1,m1): # remove (m'_1|.) (m''_1 o m_2)
                k2 = m1-k1-1
                l1 = sorted(p1._list+[k1],reverse=True)
                l2 = sorted(p2._list+[k2+1,m2+1],reverse=True)
                if (sum(l1)+len(l1)) %2 == 0 and (sum(l2)+len(l2)) %2 == 0:
                    for a in xrange(k1): # a is an angle
                        N -= (_gamma_irr_rec(Partition(l1), (1,k1,a))*
                              gamma_std(Partition(l2),(2,k2+1,m2+1)))

            for k1 in xrange(1,m2): # remove (m_1 o m'_2) (m''_2|.)
                k2 = m2-k1-1
                l1 = sorted(p1._list+[m1,k1],reverse=True)
                l2 = sorted(p2._list+[k2+2],reverse=True)
                if (sum(l1)+len(l1)) %2 == 0 and (sum(l2)+len(l2)) %2 == 0:
                    for a in xrange(1,k2+1): # a is an angle for standard perm
                        N -= (_gamma_irr_rec(Partition(l1), (2,m1,k1)) *
                              gamma_std(Partition(l2),(1,k2+2,a)))

        for m in pp.to_exp_dict(): # remove (m_1, k_1) (k_2, m_2) for k1+k2+1 an other zero
            q = pp._list[:]
            del q[q.index(m)]
            for p1,p2 in bidecompositions(Partition(q)):
                for k1 in xrange(1,m):
                    k2 = m-k1-1
                    l1 = sorted(p1._list+[m1,k1],reverse=True)
                    l2 = sorted(p2._list+[k2+1,m2+1],reverse=True)
                    if (sum(l1)+len(l1))%2 == 0 and (sum(l2)+len(l2))%2 == 0:
                        N -= (_gamma_irr_rec(Partition(l1), (2,m1,k1)) *
                              gamma_std(Partition(l2),(2,k2+1,m2+1)))

        return N

    else:
        raise ValueError, "marking must be a 3-tuple of the form (1,m,a) or (2,m1,m2)"

def gamma_irr(profile=None, marking=None):
    r"""
    Number of permutations for the given profile and marking

    INPUT:

    - ``profile`` - an integer partition such that its sum plus its length is
      congruent to 0 modulo 2

    - ``markings`` - None, an element of the profile or a 3-tuple

    EXAMPLES::

        sage: import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rcc

    The connected strata in genus 3::

        sage: from surface_dynamics.all import AbelianStratum

        sage: c = AbelianStratum(1,1,1,1).unique_component()
        sage: c.rauzy_diagram()
        Rauzy diagram with 1255 permutations
        sage: rcc.gamma_irr([2,2,2,2])
        1255

        sage: c = AbelianStratum(2,1,1).unique_component()
        sage: c.rauzy_diagram()
        Rauzy diagram with 2177 permutations
        sage: rcc.gamma_irr([3,2,2])
        2177

        sage: c = AbelianStratum(3,1).unique_component()
        sage: c.rauzy_diagram()
        Rauzy diagram with 770 permutations
        sage: rcc.gamma_irr([4,2])
        770

    The non connecte strata in genus 3::

        sage: c_odd = AbelianStratum(2,2).odd_component()
        sage: c_hyp = AbelianStratum(2,2).hyperelliptic_component()
        sage: c_odd.rauzy_diagram()
        Rauzy diagram with 294 permutations
        sage: c_hyp.rauzy_diagram()
        Rauzy diagram with 63 permutations
        sage: rcc.gamma_irr([3,3]) == 294 + 63
        True

        sage: c_odd = AbelianStratum(4).odd_component()
        sage: c_hyp = AbelianStratum(4).hyperelliptic_component()
        sage: c_odd.rauzy_diagram()
        Rauzy diagram with 134 permutations
        sage: c_hyp.rauzy_diagram()
        Rauzy diagram with 31 permutations
        sage: rcc.gamma_irr([5]) == 134 + 31
        True
    """
    p = Partition(sorted(profile,reverse=True))

    if marking is None:
        return sum(gamma_irr(profile,m) for m in marking_iterator(profile,left=None))

    elif isinstance(marking, (int,Integer)):
        return sum(gamma_irr(profile,m) for m in marking_iterator(profile,left=marking))

    marking = check_marking(p, marking)
    return _gamma_irr_rec(p, marking)

number_of_irreducible_permutations = gamma_irr

@cached_function
def _delta_irr_rec(p, marking):
    r"""
    Internal recursive function called by :func:`delta_irr`.
    """
    if len(p) == 0:
        return 0

    if marking[0] == 1:
        m = marking[1]
        a = marking[2]
        i = p.index(m)
        pp = Partition(p._list[:i]+p._list[i+1:]) # the partition p'

        N = (-1)**a* delta_std(
                pp._list + [m+2],
                (1,m+2,m-a))

        for m1 in xrange(1,m-1,2):
            m2 = m-m1-1
            for a1 in xrange(max(0,a-m2),min(a,m1)):
                a2 = a - a1 - 1
                for p1,p2 in bidecompositions(pp):
                    l1 = sorted([m1]+p1._list,reverse=True)
                    l2 = sorted([m2+2]+p2._list,reverse=True)
                    N += (-1)**a2*(_delta_irr_rec(Partition(l1),(1,m1,a1)) *
                        spin_difference_for_standard_permutations(Partition(l2),(1,m2+2,m2-a2)))
        return N

    elif marking[0] == 2:
        m1 = marking[1]
        m2 = marking[2]
        i1 = p.index(m1)
        i2 = p.index(m2)
        if m1 == m2: i2 += 1
        if i2 < i1: i1,i2 = i2,i1
        pp = Partition(p._list[:i1] + p._list[i1+1:i2] + p._list[i2+1:])

        N = d(Partition(sorted(pp._list+[m1+m2+1],reverse=True))) /  pp.centralizer_size()
        # nb of standard permutations that corrresponds to extension of good
        # guys

        for p1,p2 in bidecompositions(Partition(pp)):
            for k1 in xrange(1,m1,2): # remove (k1|.) (k2 o m_2)
                k2 = m1-k1-1
                q1 = Partition(sorted(p1._list+[k1],reverse=True))
                q2 = Partition(sorted(p2._list+[k2+m2+1],reverse=True))
                for a in xrange(k1): # a is a angle
                    N += _delta_irr_rec(q1, (1,k1,a)) * d(q2) / p2.centralizer_size()

            for k1 in xrange(1,m2,2): # remove (m_1 o k1) (k2|.)
                k2 = m2-k1-1
                l1 = sorted(p1._list+[m1,k1],reverse=True)
                l2 = sorted(p2._list+[k2+2],reverse=True)
                for a in xrange(1,k2+1): # a is an angle for standard perm
                    N += (_delta_irr_rec(Partition(l1), (2,m1,k1)) *
                        spin_difference_for_standard_permutations(Partition(l2), (1,k2+2,a)))

        for m in pp.to_exp_dict(): # remove (m_1 o k_1) (k_2 o m_2) for k1+k2+1 an other zero
            q = pp._list[:]
            del q[q.index(m)]
            for p1,p2 in bidecompositions(Partition(q)):
                for k1 in xrange(1,m,2):
                    k2 = m-k1-1
                    q1 = Partition(sorted(p1._list+[m1,k1],reverse=True))
                    q2 = Partition(sorted(p2._list+[k2+m2+1],reverse=True))
                    N += _delta_irr_rec(q1, (2,m1,k1)) * d(q2) / p2.centralizer_size()

        return N

def delta_irr(profile, marking=None):
    r"""
    Spin difference for the given profile and marking

    EXAMPLES::

        sage: import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rcc

    The non connecte strata in genus 3::

        sage: from surface_dynamics.all import AbelianStratum

        sage: c_odd = AbelianStratum(2,2).odd_component()
        sage: c_hyp = AbelianStratum(2,2).hyperelliptic_component()
        sage: c_odd.rauzy_diagram()
        Rauzy diagram with 294 permutations
        sage: c_hyp.rauzy_diagram()
        Rauzy diagram with 63 permutations
        sage: rcc.delta_irr([3,3]) == 294 - 63
        True

        sage: c_odd = AbelianStratum(4).odd_component()
        sage: c_hyp = AbelianStratum(4).hyperelliptic_component()
        sage: c_odd.rauzy_diagram()
        Rauzy diagram with 134 permutations
        sage: c_hyp.rauzy_diagram()
        Rauzy diagram with 31 permutations
        sage: rcc.delta_irr([5]) == 134 - 31
        True

    A non connected strata in genus 4::

        sage: import surface_dynamics.interval_exchanges.rauzy_class_cardinality as rdc
        sage: a = AbelianStratum(6)
        sage: c_hyp = a.hyperelliptic_component()
        sage: c_odd = a.odd_component()
        sage: c_even = a.even_component()
        sage: c_hyp.rauzy_diagram()
        Rauzy diagram with 127 permutations
        sage: c_hyp.rauzy_class_cardinality()
        127
        sage: c_odd.rauzy_diagram()
        Rauzy diagram with 5209 permutations
        sage: c_odd.rauzy_class_cardinality()
        5209
        sage: c_even = a.even_component()
        sage: c_even.rauzy_diagram()
        Rauzy diagram with 2327 permutations
        sage: c_even.rauzy_class_cardinality()
        2327
        sage: 5209 - 2327 - 127
        2755
        sage: rdc.delta_irr([7])
        2755

    An example with a very big Rauzy class::

        sage: c = AbelianStratum(6,6).odd_component()
        sage: c.rauzy_class_cardinality()
        11609364656
    """
    p = Partition(sorted(profile,reverse=True))

    if marking is None:
        return sum(delta_irr(profile,m) for m in marking_iterator(profile,left=None))

    elif isinstance(marking, (int,Integer)):
        return sum(delta_irr(profile,m) for m in marking_iterator(profile,left=marking))

    marking = check_marking(p, marking)
    return _delta_irr_rec(p, marking)

spin_difference_for_irreducible_permutations = delta_irr

