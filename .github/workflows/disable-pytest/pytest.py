r"""
Pretend that pytest is not installed.

`sage -t` tries to invoke pytest on our files which does
not always work. We use this file to pretend that pytest
is not installed, see
https://trac.sagemath.org/ticket/31103#comment:49
"""
#*********************************************************************
#  This file is part of surface-dynamics.
#
#        Copyright (C) 2021-2022 Julian Rüth
#
#  surface-dynamics is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  surface-dynamics is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with surface-dynamics. If not, see <https://www.gnu.org/licenses/>.
#*********************************************************************
raise ModuleNotFoundError
