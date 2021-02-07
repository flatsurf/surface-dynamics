#!/usr/bin/env python
#*****************************************************************************
#       Copyright (C) 2021 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

import pytest

pytest.importorskip('pyintervalxt')

from sage.all import ZZ, QQ, AA, polygen, NumberField, vector
from surface_dynamics import iet
from surface_dynamics.interval_exchanges.conversion import iet_to_pyintervalxt, iet_from_pyintervalxt

def test_back_and_forth_int():
    p = iet.Permutation("a b c", "c b a")
    T = iet.IntervalExchangeTransformation(p, [12, 5, 9])
    T2 = iet_to_pyintervalxt(T)
    T3 = iet_from_pyintervalxt(T2)
    assert T == T3

    p = iet.Permutation("a b c", "c a b")
    T = iet.IntervalExchangeTransformation(p, [12, 5, 9])
    T2 = iet_to_pyintervalxt(T)
    T3 = iet_from_pyintervalxt(T2)
    assert T == T3

    p = iet.Permutation("a b c", "a c b")
    T = iet.IntervalExchangeTransformation(p, [12, 5, 9])
    T2 = iet_to_pyintervalxt(T)
    T3 = iet_from_pyintervalxt(T2)
    assert T == T3

def test_back_and_forth_rat():
    p = iet.Permutation("a b c", "c b a")
    T = iet.IntervalExchangeTransformation(p, [QQ((12,5)), QQ((2,7)), QQ((1,21))])
    T2 = iet_to_pyintervalxt(T)
    T3 = iet_from_pyintervalxt(T2)
    assert T == T3

def test_back_and_forth_nf():
    p = iet.Permutation("a b c", "c b a")
    x = polygen(QQ)
    K = NumberField(x**2 - 2, 'sqrt2', embedding=AA(2).sqrt())
    sqrt2 = K.gen()
    T = iet.IntervalExchangeTransformation(p, [1, sqrt2, sqrt2-1])
    T2 = iet_to_pyintervalxt(T)
    T3 = iet_from_pyintervalxt(T2)
    assert T == T3

    r1 = vector(ZZ, (0, 0, 0, 1, 1, 0))
    r2 = vector(ZZ, (3, 1, 0, 1, 0, 2))
    r3 = vector(ZZ, (5, 0, 1, 2, 0, 3))
    K = NumberField(x**3 - 3, 'cbrt3', embedding=AA(3)**QQ((1,3)))
    cbrt3 = K.gen()
    p = iet.Permutation(["a","b","c","d","e","f"], ["f","d","c","b","a","e"])
    l = r1 + (cbrt3 - 1) * r2 + cbrt3**2 * r3
    T = iet.IntervalExchangeTransformation(p, l)
    assert T.sah_arnoux_fathi_invariant().is_zero()
    T2 = iet_to_pyintervalxt(T)
    assert list(T2.safInvariant()) == [0,0,0]
    T3 = iet_from_pyintervalxt(T2)
    assert T == T3
