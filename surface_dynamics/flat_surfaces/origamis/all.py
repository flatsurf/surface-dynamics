# ****************************************************************************
#       Copyright (C) 2011-2019 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import absolute_import

from .origami import Origami
from .pillowcase_cover import PillowcaseCover

from .teichmueller_curve import TeichmuellerCurvesOfOrigamis
from .generators import origamis

from .origami_database import OrigamiDatabase

del absolute_import
