#!/usr/bin/env python
#*****************************************************************************
#       Copyright (C) 2021 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

def test_substitution_additive():
    from surface_dynamics.misc.additive_multivariate_generating_series import AdditiveMultivariateGeneratingSeriesRing

    A = AdditiveMultivariateGeneratingSeriesRing('x', 3)
    R = A.polynomial_ring()
    x0, x1, x2 = R.gens()

    f1 = A.term(x0 - x1, [((1,0,0), 1), ((1,1,0),2)]) + A.term(2*x0**2, [((0,1,0), 2), ((1,0,2),3)])
    f2 = A.term(1, [((1,1,0),1), ((1,1,1),1)]) + A.term(x0 - x1 - x2, [((1,0,0),2)])

    sdict1 = {x0: x0+x1, x1:x0-x1, x2: 2*x2}
    sdict2 = {x0: x1+x2}
    for f in [f1, f2]:
        for sdict in [sdict2, sdict2]:
            assert f.to_rational_function().subs(sdict) == \
                   f.subs(sdict).to_rational_function() == \
                   f.subs(**{str(k): v for k,v in sdict.items()}).to_rational_function()
