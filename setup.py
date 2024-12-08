"""
Setup script for the surface-dynamics package.

This script defines the compilation of the .pyx files into extension modules
which cannot be configured in pyproject.toml yet.
"""
# ****************************************************************************
#       Copyright (C) 2024 Julian Rüth
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import numpy
from Cython.Build import cythonize
from setuptools import setup, Extension

# NOTE: without this option, tab-completion and documentation are mostly broken within sage
# See https://trac.sagemath.org/ticket/31632
import Cython.Compiler.Options
Cython.Compiler.Options.embed_pos_in_docstring = True

# NOTE: When adding files here, make sure to duplicate the list in the inputs
# of the build task in pyproject.toml.
extensions = [
    Extension(
        "surface_dynamics.flat_surfaces.origamis.origami_dense",
        sources=[
            "surface_dynamics/flat_surfaces/origamis/origami_dense.pyx",
            "surface_dynamics/flat_surfaces/origamis/normal_form.c",
            "surface_dynamics/flat_surfaces/origamis/lyapunov_exponents.c",
        ],
        include_dirs=[numpy.get_include(), "surface_dynamics/flat_surfaces/origamis"],
    ),
    Extension(
        "surface_dynamics.interval_exchanges.lyapunov_exponents",
        sources=[
            "surface_dynamics/interval_exchanges/lyapunov_exponents.pyx",
            "surface_dynamics/interval_exchanges/generalized_permutation.c",
            "surface_dynamics/interval_exchanges/lin_alg.c",
            "surface_dynamics/interval_exchanges/quad_cover.c",
            "surface_dynamics/interval_exchanges/random.c",
            "surface_dynamics/interval_exchanges/permutation.c",
        ],
        include_dirs=[numpy.get_include()],
    ),
    Extension(
        "surface_dynamics.interval_exchanges.integer_iet",
        sources=[
            "surface_dynamics/interval_exchanges/integer_iet.pyx",
            "surface_dynamics/interval_exchanges/int_iet.c",
            "surface_dynamics/interval_exchanges/int_vector.c",
        ],
        include_dirs=[numpy.get_include()],
    )
]

try:
    import ppl as _
except ImportError:
    import sys
    print('Warning: pplpy not installed. Will not compile iet_family.', file=sys.stderr)
else:
    extensions.append(Extension(
        "surface_dynamics.interval_exchanges.iet_family",
        sources=[
            "surface_dynamics/interval_exchanges/iet_family.pyx",
        ],
        include_dirs=[numpy.get_include()],
    ))

# Work around changes in SageMath 9.7, see
# https://trac.sagemath.org/wiki/ReleaseTours/sage-9.7#Packagessagesage.rings...arenownamespaces,
# i.e., work around import errors with Cython <3. This is also the reason we
# cannot currently build with meson.
try:
    from sage.misc.package_dir import cython_namespace_package_support
except (ImportError, ModuleNotFoundError):
    from contextlib import nullcontext as cython_namespace_package_support

with cython_namespace_package_support():
    extensions = cythonize(extensions, compiler_directives={"language_level": "3"})


setup(ext_modules=extensions)
