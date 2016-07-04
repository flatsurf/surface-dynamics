r"""
Plane trees.

Generation of plane trees
-------------------------

This file contains generators for level sequences of

 * rooted (free) trees
 * (free) trees
 * rooted plane trees
 * plane trees

REFERENCES:

  * Beyer-Hedetniemi "Constant time generation of rooted trees", 1980
  Generation of rooted trees (a root but no cyclic order at vertices). They use
  a successor function on level sequences.
  (good explanation of the algorithm)
  A regular tree is a tree such that for all vertices, its sequence of subtrees
  is in decreasing order.
  How do PREV and SAVE work ?

  * Wright-Richmond-Odlyzko-McKay "Constant time generation of free trees", 1986
  TODO: look at

.. [Nak02]   Nakano "Efficient generation of plane trees", 2002
  For him, rooted means left-right ordering of children of each node.
  Enumerating rooted plane trees with at most `n` vertices is elementary (a
  depth-first search in the tree of trees).
"""

from sage.rings.integer import Integer

###################################
# Rooted trees (Beyer-Hedetniemi) #
###################################

def _TMP_successor(l):
    r"""
    given a level sequence generate the next (Beyer-Hedetniemi)
    """
    # find p the last index in l which is not 1
    p = len(l)-1
    while p > 0 and l[p] == 1:
        p -= 1
    if p == 0:
        return False
    # find q the parent of p
    q = p-1
    while q > 0 and l[q] != l[p]-1:
        q -= 1
    # build the new level sequence
    # l[:p] does not change and then copy as many times as you can the level
    # sequence l[q:p].
    diff = p-q
    ll = l[:p] + l[q:p]*((len(l)-p)/diff)
    ll.extend(l[q:q+len(l)-len(ll)])

    return ll

def _TMP_rooted_trees(n):
    r"""
    Return the list of rooted trees using the preceding function
    """
    l = range(n)
    res = [l]
    ll = _TMP_successor(res[-1])
    while ll:
        res.append(ll)
        ll = _TMP_successor(res[-1])
    return res

def rooted_tree_iterator(n,verbose=False):
    r"""
    Iterator through regular level sequences of rooted trees.
    (only works for n >= 3)

    EXAMPLES::

        sage: from surface_dynamics.misc.plane_tree import rooted_tree_iterator
        sage: for t in rooted_tree_iterator(4): print t
        [0, 1, 2, 3]
        [0, 1, 2, 2]
        [0, 1, 2, 1]
        [0, 1, 1, 1]
        sage: for t in rooted_tree_iterator(5): print t
        [0, 1, 2, 3, 4]
        [0, 1, 2, 3, 3]
        [0, 1, 2, 3, 2]
        [0, 1, 2, 3, 1]
        [0, 1, 2, 2, 2]
        [0, 1, 2, 2, 1]
        [0, 1, 2, 1, 2]
        [0, 1, 2, 1, 1]
        [0, 1, 1, 1, 1]

        sage: for t in rooted_tree_iterator(5,verbose=True): pass
          p =    4
          prev = [0, 1, 2, 3, 4]
          save = [0, 0, 0, 0, 0]
        [0, 1, 2, 3, 4]
          p =    4
          prev = [0, 1, 2, 3, 4]
          save = [0, 0, 0, 0, 0]
        [0, 1, 2, 3, 3]
          p =    4
          prev = [0, 1, 2, 3, 4]
          save = [0, 0, 0, 0, 0]
        [0, 1, 2, 3, 2]
          p =    3
          prev = [0, 1, 2, 0, 4]
          save = [0, 0, 0, 0, 0]
        [0, 1, 2, 3, 1]
          p =    4
          prev = [0, 1, 3, 0, 4]
          save = [0, 0, 0, 2, 0]
        [0, 1, 2, 2, 2]
          p =    3
          prev = [0, 1, 2, 0, 4]
          save = [0, 0, 0, 2, 0]
        [0, 1, 2, 2, 1]
          p =    4
          prev = [0, 3, 2, 0, 4]
          save = [0, 0, 0, 1, 0]
        [0, 1, 2, 1, 2]
          p =    2
          prev = [0, 1, 0, 0, 4]
          save = [0, 0, 0, 1, 0]
        [0, 1, 2, 1, 1]
          p =    0
          prev = [0, 0, 0, 0, 4]
          save = [0, 0, 0, 1, 0]
        [0, 1, 1, 1, 1]
    """
    assert n >= 3

    l = range(n)
    prev = range(n)  # function: level -> ?
    save = [0]*n
    p = n-1

    if verbose:
        print "  p =   ", p
        print "  prev =", prev
        print "  save =", save
        print l
    yield l

    while p > 0:
        l[p] = l[p] - 1
        if p < n and (l[p] != 1 or l[p-1] != 1):
            diff = p - prev[l[p]] # = p-q
            while p < n-1:
                save[p] = prev[l[p]]
                prev[l[p]] = p
                p += 1
                l[p] = l[p-diff]
        while l[p] == 1:
            p -= 1
            prev[l[p]] = save[p]

        if verbose:
            print "  p =   ", p
            print "  prev =", prev
            print "  save =", save
            print l
        yield l

###############################
# Rooted plane trees (Nakano) #
###############################

def rooted_plane_tree_iterator(nmin,nmax=None, verbose=False):
    r"""
    Iterator through rooted plane trees with at least ``nmin`` and at most
    ``nmax`` vertices. If ``nmax`` is not furnished then it is automatically set
    to ``nmin``.

    The output are list which corresponds to level sequences. The algorithm,
    which corresponds to a depth first search in a tree of trees is due to
    Nakano [Nak02]_.

    EXAMPLES::

        sage: from surface_dynamics.misc.plane_tree import rooted_plane_tree_iterator

        sage: for t in rooted_plane_tree_iterator(4): print t
        [0, 1, 2, 3]
        [0, 1, 2, 2]
        [0, 1, 2, 1]
        [0, 1, 1, 2]
        [0, 1, 1, 1]
        sage: for t in rooted_plane_tree_iterator(1,3): print t
        [0]
        [0, 1]
        [0, 1, 2]
        [0, 1, 1]
    """
    nmin = int(nmin)
    if nmax is None:
        nmax = nmin
    else:
        nmax = int(nmax)
        if nmin > nmax:
            raise ValueError, "nmin (=%d) should be lesser than nmax (=%d)"%(nmin,nmax)
    l = [0]
    b = [] # current branch
    while True:
        if len(l) >= nmin:
            yield l
        if len(l) < nmax: # go down
            if verbose: print "go down"
            if verbose: print " append: %d"%(l[-1]+1)
            l.append(l[-1]+1)
            b.append(l[-1])
        else: # go up
            if verbose: print "go up"
            while b and b[-1] == 1:
                if verbose: print " pop"
                b.pop(-1)
                l.pop(-1)
            if b == []:
                break
            if verbose: print " modify last to become %d"%b[-1]
            b[-1] -= 1
            l[-1] = b[-1]

def half_turn(t,i=None):
    r"""
    Half turn of canonical form of a tree t which has an odd diameter.

    INPUT:

    - ``s`` - position of 1 in t

    - ``i`` - the guy of depth m-1 as an index in s

    EXAMPLES::

        sage: from surface_dynamics.misc.plane_tree import half_turn

        sage: half_turn([0,1,2,2,2,1])
        [0, 1, 2, 1, 1, 1]
        sage: half_turn([0,1,2,1,1,1])
        [0, 1, 2, 2, 2, 1]

        sage: half_turn([0,1,2,3,2,1,1,2,1])
        [0, 1, 2, 2, 3, 2, 1, 2, 1]
        sage: half_turn([0,1,2,2,3,2,1,2,1])
        [0, 1, 2, 3, 2, 1, 1, 2, 1]
    """
    if i is None:
        i = 2
        while i < len(t) and t[i] != 1:
            i += 1
    return [0,1] + [j+1 for j in t[i:]] + [j-1 for j in t[2:i]]

def cmp_halves(t,i=None):
    r"""
    Compare the subtrees [2:i] and [i:].

    Used in the case of a central edge in the tree.
    """
    t1 = [j-1 for j in t[2:i]]
    t2 = t[i:]
    return cmp(t1,t2)

def cmp_subtree(t1,d1,t2,d2):
    test = cmp(d1,d2)
    if test: return test

    return cmp(t1,t2)

def is_lyndon(t,s=None,d=None,verbose=False):
    r"""
    Returns
    Lyndon -> True
    periodic -> period
    otherwise -> False

    ALGORITHM

    Lothaire, Combinatorics on words

    INPUT:

    - ``t`` - the tree

    - ``s`` - the positions of 1 in t (limit of subtrees)

    - ``d`` - the depths of subtrees

    EXAMPLES::

        sage: from surface_dynamics.misc.plane_tree import is_lyndon

        sage: is_lyndon([1,2,3,3,1,2])
        True
        sage: is_lyndon([1,2,3,2,3,4,1,2,1,2,3,4])
        False

        sage: is_lyndon([1,1,2,1,2])
        False
        sage: is_lyndon([1,2,1,1,2])
        False
        sage: is_lyndon([1,2,1,2,1])
        True

        sage: is_lyndon([1,1,2,3,2,1,2,3,2,1,2,3,2])
        False
        sage: is_lyndon([1,2,3,2,1,1,2,3,2,1,2,3,2])
        False
        sage: is_lyndon([1,2,3,2,1,2,3,2,1,1,2,3,2])
        False
        sage: is_lyndon([1,2,3,2,1,2,3,2,1,2,3,2,1])
        True

        sage: is_lyndon([1,2,1,2])
        1

        sage: is_lyndon([1,2,3,1,2,1,2,3,1,2]) # lyndon periodic
        2
        sage: is_lyndon([1,2,1,2,3,1,2,1,2,3]) # not lyndon periodic
        False
    """
    if s is None:
        s = [i for i in xrange(len(t)) if t[i] == 1]
    n = len(s)
    if n <= 1:
        return True
    w = [t[s[i]+1:s[i+1]] for i in xrange(len(s)-1)]
    w.append(t[s[-1]+1:])
    if d is None:
        d = [max([1]+ww) for ww in w]
    if verbose:
        print "w=",w
    i,j=0,1
    while j < n:
        if verbose:
            print " i,j =",i,j
            print " w[i] =",w[i]
            print " w[j] =",w[j]
        c = cmp_subtree(w[i],d[i],w[j],d[j])
        if c == 0:
            i += 1
            j += 1
        elif c > 0:
            i = 0
            j += 1
        else:
            return False
    # i = 0: Lyndon, j-i = i: periodic ?
    if verbose:
        print " out", i,j
    if i == 0:
        return True
    elif n%(j-i) == 0:
        return j-i
    return False

def unrooted_plane_tree_iterator(nmin,nmax=None,verbose=False, check=False):
    r"""
    Iterator through plane trees with at most nmax vertices. A representative of
    each equivalence class of rooted tree is returned. The order is decreasing
    among level sequences.

    One can choose a canonical representative by choosing either the middle
    vertex as a root (and the maximal cyclic order among subtrees) or the middle
    edge as the left-most edge (and we have to check between two).

    Sloane's A002995
    1, 1, 1, 1, 2, 3, 6, 14, 34, 95, 280, 854, 2694, 8714, 28640, 95640, 323396,
    1105335, 3813798, 13269146

    TODO:

    one more restriction for go down: if the max depth is attained only at first
    position and there is not enough available vertices to go deeply enough.

    EXAMPLES::

        sage: from surface_dynamics.misc.plane_tree import unrooted_plane_tree_iterator
        sage: for t in unrooted_plane_tree_iterator(4): print t
        [0, 1, 2, 1]
        [0, 1, 1, 1]
        sage: for t in unrooted_plane_tree_iterator(5): print t
        [0, 1, 2, 2, 1]
        [0, 1, 2, 1, 2]
        [0, 1, 1, 1, 1]
        sage: for t in unrooted_plane_tree_iterator(6): print t
        [0, 1, 2, 3, 1, 2]
        [0, 1, 2, 2, 2, 1]
        [0, 1, 2, 2, 1, 2]
        [0, 1, 2, 2, 1, 1]
        [0, 1, 2, 1, 2, 1]
        [0, 1, 1, 1, 1, 1]

        sage: for t in unrooted_plane_tree_iterator(1,5): print t
        [0]
        [0, 1]
        [0, 1, 2, 2, 1]
        [0, 1, 2, 1]
        [0, 1, 2, 1, 2]
        [0, 1, 1]
        [0, 1, 1, 1]
        [0, 1, 1, 1, 1]

        sage: [sum(1 for _ in unrooted_plane_tree_iterator(i)) for i in xrange(1,13)]
        [1, 1, 1, 2, 3, 6, 14, 34, 95, 280, 854, 2694]
    """
    nmin = int(nmin)
    if nmax is None:
        nmax = nmin
    else:
        nmax = int(nmax)
        if nmin > nmax:
            raise ValueError, "nmin (=%d) should be lesser than nmax (=%d)"%(nmin,nmax)
    nnmax = nmax//2

    #TODO: start at the right place and remove trivial cases in recursion
    t = [0]    # current tree
    b = []     # current branch in the tree of trees
    s = []     # positions of 1 in t (or limit of subtrees)
    d = []     # depths of subtrees of t
    d_occ = [] #TODO: nb of occurrences in d
    while True:
        if verbose:
            print "  t =", t
            print "  d =", d
            print "  s =",s
        # verifications
        if check:
            assert len(d) == len(s)
            assert all(t[i] == 1 for i in s)
            assert (len(d) == 0 and t == [0]) or (max(d) == max(t))
        # do we need to output t
        if len(t) >= nmin:
            if t == [0]: #trivial case
                yield t
            else:
                #TODO: remove the loop and use d_occ
                m = d[0]       # depth of first subtree
                m_nb = 1       # nb of subtrees with depth m
                has_mm = False # is there a subtree with depth m-1
                for i in xrange(1,len(d)):
                    if d[i] > m: # the deepest branch is not the first one
                        break
                    elif d[i] == m:
                        m_nb += 1
                    elif d[i] == m-1:
                        has_mm = True
                else: # we did not break the loop
                    if ((m == 1) or
                        (m_nb == 1 and has_mm and cmp_halves(t,s[1]) >= 0) or
                        (m_nb > 1 and is_lyndon(t,s,d))):
                        if verbose: print "YIELD"
                        yield t
                    elif verbose:
                        print "NOT YIELD"
        # do we go down
        #TODO: remove len(s) <= 1 (corresponds to the leftmost branch in the
        #      tree of trees)
        if (len(t) < nmax and
            ((len(s) <= 1 and t[-1] <= nnmax) or
             (len(s) > 1 and t[-1] <= d[0]))):
            if verbose: print "go down"
            if verbose: print " append: %d"%(t[-1]+1)
            #TODO: append the right thing which is not necessarily t[-1]+1
            t.append(t[-1]+1)
            b.append(t[-1])
            if t[-1] == 1: # trivial case (t was [0] before)
                if check:
                    assert len(s) == 0
                    assert len(d) == 0
                s.append(len(t)-1)
                d.append(1)
            elif t[-1] > d[-1]:
                d[-1] = t[-1]
        # or go up
        else:
            if verbose: print "go up"
            while b and b[-1] == 1:
                if verbose: print " pop"
                b.pop(-1)
                t.pop(-1)
                if check:
                    assert d[-1] == 1
                    assert s[-1] == len(t)
                d.pop(-1)
                s.pop(-1)
            if b == []:
                break
            if verbose:
                print " modify last to become %d"%(b[-1]-1)
            b[-1] -= 1
            t[-1] = b[-1]
            if d[-1] == t[-1]+1: # update the depth of the rightmost branch or
                                 # the one before if t[-1] == 1
                i = 1
                d[-1] = 1
                for i in xrange(s[-1]+1,len(t)-1):
                    if t[i] > d[-1]:
                        d[-1] = t[i]
                if verbose: print " updated length of rightmost which is now %d"%d[-1]
            if t[-1] == 1: # new branch
                d.append(1)
                s.append(len(t)-1)

##########################################################################
# Modified version of Nakano's algorithm to generate separatrix diagrams #
##########################################################################

def _TMP_admissible_plane_tree_iterator(a, verbose=False):
    r"""
    Iterator through plane trees with at most nmax vertices.

    INPUT:

    - ``a`` - the angle at the vertices divided by pi.

    OUTPUT:

    - ``(t,n,l)`` - ``t`` is a tree, ``n`` its number of nodes and ``l`` its
      number of leaves.

    ALGORITHM:

    variables in the algorithm

    t : the current tree
    b : the current branch explored in the tree of trees
    n : number of nodes of t except the root
    l : number of leaves of t

    the tree is valid if 2n >= a and 2n-l <= a
    """
    t = [0] # the tree
    b = []  # current branch
    n = 0   # number of nodes (=len(t)-1)
    l = 1   # number of leaves
    while True:
        if verbose:
            print t
            print "n = %2d,   l = %2d" %(n,l)
        if (2*n >= a) and (2*n-l <= a):
            if verbose: print "yield"
            t.append(0)
            yield t, n, l
            t.pop(-1)

        if n <= a and 2*n-l < a: # could add more nodes -> go down
            if verbose:
                print "go down"
                print "  append: %d"%(t[-1]+1)
            t.append(t[-1]+1)
            b.append(t[-1])
            n += 1
            # l is unchanged

        else: # we are too far -> go up
            if verbose: print "go up"
            # update l and n, then pop the rightmost node
            while b and b[-1] == 1:
                if verbose: print " pop"
                b.pop(-1)
                t.pop(-1)
                n -= 1
                l -= 1
            if b == []:
                break
            if verbose:
                print " modify last"
            if len(t) > 1 and t[-2] == t[-1] - 1:
                l += 1
            b[-1] -= 1
            t[-1] = b[-1]

        if verbose: print

def admissible_plane_tree_iterator(a,verbose=False):
    r"""
    Plane tree for separatrix diagram of hyperelliptic strata. Return triple
    (t,n,l) where t is a tree, n its number of nodes and l its number of leaves.

    A plane tree with n nodes and l leaves is *a-admissible* if it satisfies
    `n <= a` and `2n-l < a`

    EXAMPLES::

        sage: from surface_dynamics.misc.plane_tree import admissible_plane_tree_iterator
        sage: for t,n,l in admissible_plane_tree_iterator(5): print t,"\n",n,l
        [0, 1, 2, 2, 1, 0]
        4 3
        [0, 1, 2, 1, 0]
        3 2
        [0, 1, 1, 1, 0]
        3 3
        [0, 1, 1, 1, 1, 0]
        4 4
        [0, 1, 1, 1, 1, 1, 0]
        5 5
    """
    aa = (a+1)//2

    #TODO: start at the right place and remove trivial cases in recursion
    t = [0]    # current tree
    b = []     # current branch in the tree of trees
    s = []     # positions of 1 in t (or limit of subtrees)
    d = []     # depths of subtrees of t
    n = 0      # number of nodes (=len(t)-1)
    l = 1      # number of leaves
    d_occ = [] #TODO: nb of occurrences in d
    while True:
        if verbose:
            print "  t =", t
            print "  d =", d
            print "  s =",s
        # do we need to output t
        if (2*n >= a) and (2*n-l <= a):
            #TODO: remove the loop and use d_occ
            m = d[0]       # depth of first subtree
            m_nb = 1       # nb of subtrees with depth m
            has_mm = False # is there a subtree with depth m-1
            for i in xrange(1,len(d)):
                if d[i] > m: # the deepest branch is not the first one
                    break
                elif d[i] == m:
                    m_nb += 1
                elif d[i] == m-1:
                    has_mm = True
            else: # we did not break the loop
                if ((m == 1) or
                    (m_nb == 1 and has_mm and cmp_halves(t,s[1]) >= 0) or
                    (m_nb > 1 and is_lyndon(t,s,d))):
                    if verbose: print "YIELD"
                    t.append(0)
                    yield t,n,l
                    t.pop(-1)
                elif verbose:
                    print "NOT YIELD"
        # do we go down
        #TODO: remove len(s) <= 1 (corresponds to the leftmost branch in the
        #      tree of trees)
        if ((n <= a and 2*n-l < a) and
           ((len(s) <= 1 and t[-1] <= aa) or (len(s) > 1 and t[-1] <= d[0]))):
            if verbose: print "go down"
            if verbose: print " append: %d"%(t[-1]+1)
            #TODO: append the right thing which is not necessarily t[-1]+1
            t.append(t[-1]+1)
            b.append(t[-1])
            n += 1
            if t[-1] == 1: # trivial case (t was [0] before)
                s.append(len(t)-1)
                d.append(1)
            elif t[-1] > d[-1]:
                d[-1] = t[-1]
        # or go up
        else:
            if verbose: print "go up"
            while b and b[-1] == 1:
                if verbose: print " pop"
                b.pop(-1)
                t.pop(-1)
                d.pop(-1)
                s.pop(-1)
                n -= 1
                l -= 1
            if b == []:
                break
            if verbose:
                print " modify last to become %d"%(b[-1]-1)
            if len(t) > 1 and t[-2] == t[-1] - 1:
                l += 1
            b[-1] -= 1
            t[-1] = b[-1]
            if d[-1] == t[-1]+1: # update the depth of the rightmost branch or
                                 # the one before if t[-1] == 1
                i = 1
                d[-1] = 1
                for i in xrange(s[-1]+1,len(t)-1):
                    if t[i] > d[-1]:
                        d[-1] = t[i]
                if verbose: print " updated length of rightmost which is now %d"%d[-1]
            if t[-1] == 1: # new branch
                d.append(1)
                s.append(len(t)-1)

