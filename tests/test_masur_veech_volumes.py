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
import sage.all

pytest.importorskip('sage.modular.multiple_zeta')

from sage.all import QQ, Multizeta as M
from surface_dynamics import AbelianStratum
from surface_dynamics.flat_surfaces.masur_veech_volumes import masur_veech_volume

def as_pure_zeta(vol):
    support = set(sum(u) for u in vol.support())
    support = [d for d in support if not vol.homogeneous_component(d).is_zero()]
    if len(support) != 1:
        raise ValueError("not homogeneous")
    d = support.pop()
    vol = vol.homogeneous_component(d)
    u = vol.phi_as_vector()
    zeta = M(d)
    v = zeta.phi_as_vector()
    for i in range(len(u)):
        if u[i]:
            break
    result = u[i]/v[i] * zeta
    if result != vol:
        raise ValueError("not pure zeta phi(vol)={} vs phi(Z({})) = {}".format(u, d, v))
    return result

def test_H2_cylinder_diagram_contributions():
    A = AbelianStratum(2).unique_component()
    c1 = sum(c.volume_contribution() for c in A.cylinder_diagrams(1)).to_mzv()
    c2 = sum(c.volume_contribution() for c in A.cylinder_diagrams(2)).to_mzv()
    assert c1 == M(4)/3
    assert c2 == (2*M(1,3) + M(2,2)) / 3

    vol = c1 + c2
    assert vol == masur_veech_volume(A, rational=True) * M(2*A.stratum().genus())

def test_H11_cylinder_diagram_contributions():
    A = AbelianStratum(1,1).hyperelliptic_component()
    c1 = sum(c.volume_contribution() for c in A.cylinder_diagrams(1)).to_mzv()
    c2 = sum(c.volume_contribution() for c in A.cylinder_diagrams(2)).to_mzv()
    c3 = sum(c.volume_contribution() for c in A.cylinder_diagrams(3)).to_mzv()

    assert c1 == M(5)/6
    assert c2 == M(2)*M(3)/3 - M(5)/6
    assert c3 == (2*M(4) - M(2)*M(3))/3

    vol = c1 + c2 + c3
    assert vol == masur_veech_volume(A, rational=True) * M(2*A.stratum().genus())

def test_H4hyp_cylinder_diagram_contributions():
    A = AbelianStratum(4).hyperelliptic_component()
    c1 = sum(c.volume_contribution() for c in A.cylinder_diagrams(1)).to_mzv()
    c2 = sum(c.volume_contribution() for c in A.cylinder_diagrams(2)).to_mzv()
    c3 = sum(c.volume_contribution() for c in A.cylinder_diagrams(3)).to_mzv()

    vol = c1 + c2 + c3
    assert vol == masur_veech_volume(A, rational=True) * M(2*A.stratum().genus())

def test_H4odd_cylinder_diagram_contributions():
    A = AbelianStratum(4).odd_component()
    c1 = sum(c.volume_contribution() for c in A.cylinder_diagrams(1)).to_mzv()
    c2 = sum(c.volume_contribution() for c in A.cylinder_diagrams(2)).to_mzv()
    c3 = sum(c.volume_contribution() for c in A.cylinder_diagrams(3)).to_mzv()

    vol = c1 + c2 + c3
    assert vol == masur_veech_volume(A, rational=True) * M(2*A.stratum().genus())

def test_H31_cylinder_diagram_contributions():
    # values from Table 3 of A. Zorich "Square tiled surfaces
    # and Teichmueller volumes of the moduli spaces of Abelian
    # differentials" (2002)
    A = AbelianStratum(3,1).unique_component()
    c1 = sum(c.volume_contribution() for c in A.cylinder_diagrams(1)).to_mzv()
    c2 = sum(c.volume_contribution() for c in A.cylinder_diagrams(2)).to_mzv()
    c3 = sum(c.volume_contribution() for c in A.cylinder_diagrams(3)).to_mzv()

    assert c1 == M(7)/15
    assert c2 == (55*M(1,6) + 29*M(2,5) + 15*M(3,4) + 8*M(4,3) + 4*M(5,2)) / 45

    vol = 16 * M(6) / 45
    print (c1.n() / vol.n())
    print (c2.n() / vol.n())
    print (c3.n() / vol.n())

def test_H211_cylinder_diagram_contributions():
    A = AbelianStratum(2,1,1).unique_component()
    c1 = sum(c.volume_contribution() for c in A.cylinder_diagrams(1)).to_mzv()
    c2 = sum(c.volume_contribution() for c in A.cylinder_diagrams(2)).to_mzv()
    c3 = sum(c.volume_contribution() for c in A.cylinder_diagrams(3)).to_mzv()

    assert c1 == 7* M(8)/180
    assert c2 == (1620*M(1,7) + 850*M(2,6) + 436*M(3,5) + 231*M(4,4) + 130*M(5,3) +
                 65*M(6,2) + 35*M(7) - 35*M(8)) / 1260


