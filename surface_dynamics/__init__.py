#*****************************************************************************
#       Copyright (C) 2019 Vincent Delecroix <20100.delecroix@gmail.com>
#                     2024 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from .version import version as __version__

import warnings
warnings.filterwarnings('default',
    r'\[surface_dynamics].*')

# Make sure that sage's imports are going to resolve in the correct order
import sage.all

from surface_dynamics.flat_surfaces.all import *
from surface_dynamics.interval_exchanges.all import *
from surface_dynamics.misc.constellation import Constellation, Constellations
from surface_dynamics.topological_recursion import *
from surface_dynamics.topology.all import *
