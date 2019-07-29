r"""
Flat surfaces
"""

from __future__ import absolute_import

from .version import version as __version__

import warnings
warnings.filterwarnings('default',
    r'\[surface_dynamics].*')

# Make sure that sage's imports are going to resolve in the correct order
import sage.all

from .topology.all import *
from .flat_surfaces.all import *
from .interval_exchanges.all import *
from .misc.constellation import Constellation, Constellations

del absolute_import
