#!/bin/sh
#*********************************************************************
#  This file is part of surface-dynamics.
#
#        Copyright (C) 2026 Vincent Delecroix
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

# Run the SageMath doctester with the given arguments.
#
# The conda-forge builds of sagelib with the meson layout ship neither
# sage-runtests nor "sage -t", and their sage.doctest imports pytest at
# module level, so they run the doctester as a module without the
# disable-pytest shim. The older builds keep "sage -t" with the shim, since
# SageMath 10.4 and earlier have no module sage.doctest.__main__.
#
# sage-runtests is looked up next to the python of the environment rather
# than on the PATH, which may still hold another environment's bin/.
if [ -x "$(dirname "$(command -v python)")/sage-runtests" ]; then
    PYTHONPATH="$(dirname "$0")/disable-pytest${PYTHONPATH:+:$PYTHONPATH}" exec sage -t "$@"
else
    exec python -m sage.doctest "$@"
fi
