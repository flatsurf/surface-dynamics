#!/usr/bin/env python
#*****************************************************************************
#       Copyright (C) 2020 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

import sys
import pytest

from sage.all import ZZ, QQ, factorial, parent

class CarrellChapuyCount():
    def __init__(self, R):
        self.R = R
        self.count = R.zero()
    def __call__(self, cm, aut):
        x,y = R.gens()
        self.count += x**cm.num_vertices() * y**cm.num_edges() * z**cm.num_faces()

def test_carrell_chapuy():
    return
    # Q_g^{n, f}  (ne=n, nf=f)
    # 2 - 2g = nv - n + f  =>  nv = 2 - 2g - ne + nf is fixed >= 1
    # 2 - 2g + n = nv + f
    # bound on n: n = 2g-2 + nv + f
    from surface_dynamics.topology.fat_graph_exhaustive_generation import FatGraphs_g_nf_nv

    R = QQ['x,y,z']

    numbers = {}
    for g in [1,2]:
        for n in range(2*g, 2*g+4):
            # here we would like to bound nv+nf instead of nv and nf
            # (equivalently we would like to fix ne)
            print(g, n)
            nv_min = nf_min = 1
            nv_max = nf_max = 2-2*g+n-1
            print("FatGraphs_g_nf_nv(g=%d, nv_min=%d, nv_max=%d, nf_min=%d, nf_max=%d)" % (g, nv_min, nv_max, nf_min, nf_max))
            F = FatGraphs_g_nf_nv(g=g, nv_min=nv_min, nv_max=nv_max,
                                       nf_min=nf_min, nf_max=nf_max)
            C = CarrellChapuyCount(R)
            F.map_reduce(C)
            print(g, n, F.count)

