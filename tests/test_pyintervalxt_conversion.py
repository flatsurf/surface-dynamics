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

from sage.all import QQ, AA, polygen, NumberField
from surface_dynamics import iet
from surface_dynamics.interval_exchanges.conversion import iet_to_pyintervalxt, iet_from_pyintervalxt

def test_back_and_forth_int():
    p = iet.Permutation("a b c", "c b a")
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
