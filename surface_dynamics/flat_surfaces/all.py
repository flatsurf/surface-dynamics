# ****************************************************************************
#       Copyright (C) 2009-2019 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import absolute_import

from .abelian_strata import AbelianStrata, AbelianStratum
from .quadratic_strata import QuadraticStrata, QuadraticStratum
from .homology import RibbonGraph, RibbonGraphWithAngles
from .separatrix_diagram import SeparatrixDiagram, CylinderDiagram, QuadraticCylinderDiagram
from .origamis.all import *

del absolute_import
