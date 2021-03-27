#!/usr/bin/env python
#*****************************************************************************
#       Copyright (C) 2020 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

import pytest
import itertools
import random

from surface_dynamics import AbelianStratum, OrigamiDatabase

def representative(cd):
    return min([cd.canonical_label(), cd.horizontal_symmetry().canonical_label(), cd.vertical_symmetry().canonical_label(), cd.inverse().canonical_label()])

def match_list_up_to_symmetry(A, k):
    cds_compute = A.cylinder_diagrams(k, True, True)
    cds_database = A.cylinder_diagrams(k, True, False)

    assert len(cds_compute) == len(cds_database)
    assert len(cds_compute) == A.cylinder_diagrams_number(k, True, True)
    assert len(cds_compute) == A.cylinder_diagrams_number(k, True, False)

    cds_compute = [representative(cd) for cd in cds_compute]
    cds_database = [representative(cd) for cd in cds_database]
    cds_compute.sort()
    cds_database.sort()

    assert cds_compute == cds_database

    assert len(set(cds_compute)) == len(cds_database)

def match_list_no_symmetry(A, k):
    cds_compute = A.cylinder_diagrams(k, False, True)
    cds_database = A.cylinder_diagrams(k, False, False)

    assert len(cds_compute) == len(cds_database)
    assert len(cds_compute) == A.cylinder_diagrams_number(k, False, True)
    assert len(cds_compute) == A.cylinder_diagrams_number(k, False, False)

    cds_compute = [cd.canonical_label(inplace=False) for cd in cds_compute]
    cds_database = [cd.canonical_label(inplace=False) for cd in cds_database]
    cds_compute.sort()
    cds_database.sort()

    assert cds_compute == cds_database

    assert len(set(cds_compute)) == len(cds_database)

def cylinder_diagrams_testing(A):
    ccs = A.components()
    d = A.genus() + A.nb_zeros() - 1

    for k in range(1, d + 1):
        match_list_up_to_symmetry(A, k)
        match_list_no_symmetry(A, k)

        for cc in A.components():
            match_list_up_to_symmetry(cc, k)
            match_list_no_symmetry(cc, k)

def origami_check(C, sample_size, repeat):
    D = OrigamiDatabase()
    q = D.query(stratum=C.stratum(), component=C._name)
    cyl_diags = set(representative(c) for c in C.cylinder_diagrams())
    for o in itertools.islice(q, sample_size):
        for _ in range(repeat):
            o = o.horizontal_twist(random.randint(1,30)).mirror()
            assert representative(o.cylinder_diagram()) in cyl_diags, o

def test_H2():
    A = AbelianStratum(2)
    cylinder_diagrams_testing(A)
    origami_check(A.unique_component(), 10, 10)

def test_H11():
    A = AbelianStratum(1,1)
    cylinder_diagrams_testing(A)
    origami_check(A.unique_component(), 10, 10)

def test_H4():
    A = AbelianStratum(4)
    cylinder_diagrams_testing(A)
    origami_check(A.hyperelliptic_component(), 10, 10)
    origami_check(A.odd_component(), 10, 10)

@pytest.mark.slow
def test_H22():
    A = AbelianStratum(2,2)
    cylinder_diagrams_testing(A)
    origami_check(A.hyperelliptic_component(), 10, 10)
    origami_check(A.odd_component(), 10, 10)

def test_H31():
    A = AbelianStratum(3,1)
    cylinder_diagrams_testing(A)
    origami_check(A.unique_component(), 10, 10)

@pytest.mark.slow
def test_H6():
    A = AbelianStratum(6)
    cylinder_diagrams_testing(A)
    origami_check(A.hyperelliptic_component(), 10, 10)
    origami_check(A.odd_component(), 10, 10)
    origami_check(A.even_component(), 10, 10)

@pytest.mark.slow
def test_H211():
    A = AbelianStratum(2,1,1)
    cylinder_diagrams_testing(A)
    origami_check(A.unique_component(), 5, 10)
