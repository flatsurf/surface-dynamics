#!/usr/bin/env python
r"""
Check the Carrell-Chapuy formula "Simple recurrence formulas to count maps on orientable surfaces" (2015)
"""
#*****************************************************************************
#       Copyright (C) 2020-2021 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.all import ZZ, cached_function

class CarrellChapuyCount:
    r"""
    Weighted sum of combinatorial maps with weight x^#{faces of m}
    """
    def __init__(self, ne):
        self.R = ZZ['x']
        self.x = self.R.gen()
        self.counter = self.R.zero()
        self.ne = ne
    def __call__(self, cm, aut):
        # NOTE: in order to count rooted maps, we multiply by 2 ne / |Aut|
        self.counter += 2 * self.ne // (1 if aut is None else aut.group_cardinality()) * self.x ** cm.num_faces()

@cached_function
def carrell_chapuy_polynomial(g, n):
    if n == 0:
        R = ZZ['x']
        if g == 0:
            return R.gen()
        else:
            return R.zero()

    from surface_dynamics import FatGraphs
    F = FatGraphs(g, ne=n)
    C = CarrellChapuyCount(n)
    F.map_reduce(C)
    return C.counter

def check_carrell_chapuy(g, n):
    R = ZZ['x']
    x = R.gen()

    lhs = (n+1) / ZZ(6) * carrell_chapuy_polynomial(g, n)

    if n > 0:
        rhs1 = (1 + x) * (2*n - 1) / ZZ(3) * carrell_chapuy_polynomial(g, n-1)
    else:
        rhs1 = R.zero()

    if n > 1 and g > 0:
        rhs2 = (2*n - 3) * (2*n - 2) * (2*n - 1) / ZZ(12) * carrell_chapuy_polynomial(g-1, n-2)
    else:
        rhs2 = R.zero()

    # 1 <= k <= n-1
    # 0 <= i <= g
    rhs3 = 1 / ZZ(2) * sum((2*k - 1) * (2*(n-k) - 1) * carrell_chapuy_polynomial(i, k-1) * carrell_chapuy_polynomial(g-i, n-k-1) for k in range(1,n) for i in range(g+1))

    assert lhs == rhs1 + rhs2 + rhs3

def test_carrell_chapuy_initial_condition():
    x = ZZ['x'].gen()
    assert carrell_chapuy_polynomial(0, 0) == x
    assert carrell_chapuy_polynomial(0, 1) == x**2 + x
    assert carrell_chapuy_polynomial(0, 2) == 2*x**3 + 5*x**2 + 2*x

def test_carrell_chapuy():
    for n in range(1,7):
        check_carrell_chapuy(0, n)
        check_carrell_chapuy(1, n)
        check_carrell_chapuy(2, n)
        check_carrell_chapuy(3, n)
