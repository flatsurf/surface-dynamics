#!/usr/bin/env python
#*****************************************************************************
#       Copyright (C) 2026 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

import operator
from functools import reduce

import pytest

from sage.all import LyndonWords, divisors, euler_phi

from surface_dynamics.misc.necklaces import weighted_necklaces


def compose(r1, r2):
    return tuple(r1[a] for a in r2)


@pytest.mark.parametrize('k, n', [(1, 1), (1, 3), (2, 1), (2, 7), (3, 5), (4, 4)])
def test_lyndon_words(k, n):
    expected = sorted(tuple(a - 1 for a in w) for w in LyndonWords(k, n))
    assert list(weighted_necklaces([1] * k, n, n, primitive=True)) == expected


@pytest.mark.parametrize('k, n', [(1, 4), (2, 6), (2, 9), (3, 6), (4, 5)])
def test_necklace_counts(k, n):
    # the number of necklaces of length n on k letters, by Burnside's lemma
    expected = sum(euler_phi(d) * k**(n // d) for d in divisors(n)) // n
    words = list(weighted_necklaces([1] * k, n, n))
    assert len(words) == expected
    assert all(w <= w[i:] + w[:i] for w in words for i in range(n))


def test_unit_weights():
    assert sum(1 for _ in weighted_necklaces([1, 1, 1], 6, primitive=True)) == 196
    assert sum(LyndonWords(3, k).cardinality() for k in range(1, 7)) == 196
    assert sum(1 for _ in weighted_necklaces([1, 1, 1], 6, 6)) == 130


def test_lexicographic_order():
    words = list(weighted_necklaces([2, 1, 3], 9))
    assert words == sorted(words)
    assert len(set(words)) == len(words)


def test_weights():
    words = list(weighted_necklaces([2, 1, 3], 7, primitive=True))
    assert len(words) == 25
    assert all(1 <= sum([2, 1, 3][i] for i in w) <= 7 for w in words)
    assert list(weighted_necklaces([2, 1, 3], 7, 5, primitive=True)) == [
        w for w in words if sum([2, 1, 3][i] for i in w) >= 5]


def test_masks():
    words = list(weighted_necklaces([1, 1, 1], 5, primitive=True, masks=[1, 2, 4], num_bits=3))
    assert len(words) == 41
    assert words == [w for w in weighted_necklaces([1, 1, 1], 5, primitive=True) if len(set(w)) == 3]


def test_masks_weights():
    # each letter of weight w carries at most w bits, so the cut loses no word
    weights = [1, 2, 1, 3]
    masks = [0b0001, 0b0110, 0b0010, 0b1001]
    words = list(weighted_necklaces(weights, 8, primitive=True, masks=masks, num_bits=4))
    assert words == [w for w in weighted_necklaces(weights, 8, primitive=True)
                     if reduce(operator.or_, (masks[i] for i in w)) == 0b1111]


def test_permutation_monoid():
    elements = [(1, 0, 2, 3), (0, 2, 3, 1), (3, 2, 1, 0)]
    pairs = list(weighted_necklaces([1, 2, 1], 7, primitive=True,
                                    monoid_elements=elements, product=compose))
    assert [w for w, _ in pairs] == list(weighted_necklaces([1, 2, 1], 7, primitive=True))
    assert all(r == reduce(compose, (elements[i] for i in w)) for w, r in pairs)
    assert sum(1 for _, r in pairs if r == (0, 1, 2, 3)) == 4
    # sum over n <= 7 of (1/n) sum_{d | n} mu(d) p_{n/d}, where p_m = a^m + b^m
    # for the roots a, b of z^2 - 2z - 1 (two letters of weight 1, one of weight 2)
    assert len(pairs) == 129


def test_monoid_masks():
    pairs = list(weighted_necklaces([1, 1, 1], 4, 4, masks=[1, 2, 4], num_bits=3,
                                    monoid_elements=['x', 'y', 'z'], product=operator.add))
    assert [w for w, _ in pairs] == [(0, 0, 1, 2), (0, 0, 2, 1), (0, 1, 0, 2), (0, 1, 1, 2),
                                     (0, 1, 2, 1), (0, 1, 2, 2), (0, 2, 1, 1), (0, 2, 1, 2),
                                     (0, 2, 2, 1)]
    assert all(s == ''.join('xyz'[i] for i in w) for w, s in pairs)


def test_empty():
    assert list(weighted_necklaces([], 5)) == []
    assert list(weighted_necklaces([1, 1], 0)) == []
    assert list(weighted_necklaces([3], 2)) == []
    assert list(weighted_necklaces([1], 2, -1)) == [(0,), (0, 0)]


@pytest.mark.parametrize('weights', [[1, 0], [-1], [1, 1.5]])
def test_bad_weights(weights):
    with pytest.raises(ValueError):
        weighted_necklaces(weights, 3)


def test_monoid_arguments():
    with pytest.raises(ValueError):
        weighted_necklaces([1, 1], 3, monoid_elements=['x', 'y'])
    with pytest.raises(ValueError):
        weighted_necklaces([1, 1], 3, product=operator.add)
