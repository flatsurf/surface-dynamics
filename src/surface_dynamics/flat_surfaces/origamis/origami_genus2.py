#encoding=utf-8
r"""
Origamis in genus 2

In genus 2, especially for the stratum H(2), the orbits of origamis under the
action of SL(2,ZZ) are well known

REFERENCES:
.. [HuLe] Hubert Lelièvre
.. [McMul] Mc Mullen
"""
from origami import Origami
from sage.combinat.partition import Partitions,OrderedPartitions
from sage.combinat.permutation import CyclicPermutationsOfPartition
from sage.arith.misc import gcd,divisors
from sage.rings.integer import Integer

################
# stratum H(2) #
################

#TODO:
# add options (different orbits, reducible, ...)
def H2_cardinality_irreducible(n):
    return 3/8 * n^2 * (n-2) * prod([1-1/p^2 for p in prime_divisors(n)])

#TODO:
# optimize generation of lists corresponding to permutation of the origamis from
# the coordinates (i.e. do not call origami_H2_1cyl...)

def origami_H2_1cyl_iterator(n, reducible=False, output='coordinates'):
    r"""
    INPUT:

    - ``n`` - positive integer - the number of squares

    - ``reducible`` - bool (default: False) - if True returns also the reducible
      origamis.

    - ``output`` - either "coordinates" or "origami" or "standard origami"
    """
    if reducible:
        hlist = divisors(n)
    else:
        hlist = [1]

    for h in hlist:  
     for p in Partitions(Integer(n/h),length=3):
      for l in CyclicPermutationsOfPartition([p]):
        l1,l2,l3 = l[0]
        if l1 == l2 and l2 == l3:
          if reducible or l1 == 1:
            for t in xrange(l1):
              if output == "coordinates":
                yield l1,l2,l3,h,t
              elif output == "origami":
                yield origami_H2_1cyl(l1,l2,l3,h,t)
              elif output == "standard_origami":
                yield origami_H2_1cyl(l1,l2,l3,h,t).standard_form()

        elif reducible or gcd([l1,l2,l3]) == 1:
          for t in xrange(n):
            if output == "coordinates":
              yield l1,l2,l3,h,t
            elif output == "origami":
              yield origami_H2_1cyl(l1,l2,l3,h,t)
            elif output == "standard_origami":
              yield origami_H2_1cyl(l1,l2,l3,h,t).standard_form()

#TODO:
# optimize generation of lists corresponding to permutation of the origamis from
# the coordinates (i.e. do not call origami_H2_2cyl...)
def origami_H2_2cyl_iterator(n, reducible=False, output="coordinates"):
    r"""
    INPUT:

    - ``n`` - positive integer - the number of squares

    - ``reducible`` - bool (default: False) - if True returns also the reducible
      origamis.

    - ``output`` - either "coordinates" or "origami" or "standard origami"
    """
    for n1,n2 in OrderedPartitions(n,k=2):
        for w1 in divisors(n1):
          h1 = n1//w1
          for w2 in filter(lambda x: x < w1, divisors(n2)):
            h2 = n2//w2
            if reducible or gcd(h1,h2) == 1:
              d = gcd(w1,w2)
              if d == 1:
                for t1 in xrange(w1):
                  for t2 in xrange(w2):
                    if output == "coordinates":
                      yield w1,h1,t1,w2,h2,t2
                    elif output == "origami":
                      yield origami_H2_2cyl(w1,h1,t1,w2,h2,t2)
                    elif output == "standard_origami":
                      yield origami_H2_2cyl(w1,h1,t1,w2,h2,t2).standard_form()

              else:
                for t1 in xrange(w1):
                  for t2 in xrange(w2):
                      if reducible or gcd(d,h2*t1-h1*t2) == 1:
                        if output == "coordinates":
                          yield w1,h1,t1,w2,h2,t2
                        elif output == "origami":
                          yield origami_H2_2cyl(w1,h1,t1,w2,h2,t2)
                        elif output == "standard_origami":
                          yield origami_H2_2cyl(w1,h1,t1,w2,h2,t2).standard_form()

def origami_H2_iterator(n, reducible=False, output="coordinates"):
    r"""
    Iterator over l shape origami with n squares.

    INPUT:

    - ``n`` - positive integer - the number of squares

    - ``irreducible`` - bool (default: True) - return only irreducible origamis
    """
    from itertools import chain
    return chain(
        origami_H2_1cyl_iterator(n,reducible,output),
        origami_H2_2cyl_iterator(n,reducible,output))

def origamis_H2(n,reducible=False,output="coordinates"):
    return list(origami_H2_iterator(n,reducible,output))

def is_irreducible_1cyl(l1,l2,l3,h,t):
    r"""
    Test of irreducibility.
    """
    return h == 1 and gcd([l1,l2,l3]) == 1

def is_irreducible_2cyl(w1,h1,t1,w2,h2,t2):
    r"""
    Test of irreducibility.
    """
    print "testing", w1, h1, t1, w2, h2, t2
    print "  result heights gcd(%d,%d) = %d" %(h1,h2,gcd(h1,h2)), "and", gcd([w1, w2, h2*t1-h1*t2]) == 1
    return gcd(h1,h2) == 1 and gcd([w1, w2, h2*t1-h1*t2]) == 1

def origami_H2_1cyl(l1,l2,l3,h,t):
    r"""
    Returns the origami in H(2) with one cylinder.

    INPUT:

    - ``l1,l2,l3`` - three positive integers corresponding to the lengths of
      horizontal parameters

    - ``h`` - positive intger - the height of the cylinder

    - ``t`` - the twist parameter
    """
    l = l1 + l2 + l3
    z = (h-1)*l+1
    x = [None] + range(2,h*l+2)
    for i in xrange(0,h*l,l):
        x[i+l] = i+1

    y = [None] + range(l+1,l*h+1) + [None]*l
    for i in xrange(l3):
        y[z + (t+i)%l] = l1+l2+i+1
    for i in xrange(l2):
        y[z + (l3+t+i)%l] = l1+i+1
    for i in xrange(l1):
        y[z + (l3+l2+t+i)%l] = i+1

    return Origami(x[1:],y[1:])

def origami_H2_2cyl(w1,h1,t1,w2,h2,t2):
    r"""
    Returns the origami in H(2) with two cylinder

    INPUT:

    - ``w_1, h_1, t_1`` - the width, height and twist of the first cylinder

    - ``w_2, h_2, t_2``- the width, height and twist of the second cylinder
    """
    assert((w2 < w1) and (t1 < w1) and (t2 < w2))

    # v for volumes and z for z
    v1 = h1*w1
    v2 = h2*w2
    z1 = (h1-1)*w1 + 1
    z2 = v1 + (h2-1)*w2 + 1

    # the horizontal permutation
    x = [None] + range(2,v1+v2+1) + [1]
    for i in range(h1):
        x[(i+1)*w1] = i*w1 + 1
    for i in range(h2):
        x[v1 + (i+1)*w2] = v1 + i*w2 + 1

    # the vertical permutation
    y = ([None] +
        range(w1+1,v1+1) + [None]*w1 +
        range(v1+w2+1,v1+v2+1) + [None]*w2)

    for i in range(w2):
        # up-left of the first cylinder
        # print "U1L) z1 + (t1+i)%w1 -> 1+v1+i: ", z1+(t1+i)%w1, 1+v1+i
        y[z1+(t1+i)%w1] = 1+v1+i
    for i in range(w2):
        # up of the second cylinder
        # print "U2) z2+(t2+i)%w2 -> 1 + (t1+i)%w1: ", z2+(t2+i)%w2, 1+(t1+i)%w1
        y[z2+(t2+i)%w2] = 1+i
    for i in range(w1-w2):
        # up-right of the first cylinder
        # print "U1R) z1+w2+(t1+i) -> 1+i: ", z1+(w2+t1+i)%w1, 1+w2+i
        y[z1+(w2+t1+i)%w1] = 1+w2+i

    return Origami(x[1:],y[1:])

##################
# stratum H(1,1) #
##################

# This times there are four diagrams

