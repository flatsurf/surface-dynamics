# encoding=utf-8
r"""
Installation script for the flatsurf module
"""

import sys
import os
import numpy as np
import cypari2

from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize

# NOTE: without this option, tab-completion and documentation are mostly broken within sage
# See https://trac.sagemath.org/ticket/31632
import Cython.Compiler.Options
Cython.Compiler.Options.embed_pos_in_docstring = True


try:
    import ppl
except ImportError:
    sys.stderr.write('Warning: pplpy not installed. Will not compile iet_family\n')
    WITH_PPL = False
else:
    WITH_PPL = True


extensions_data = {
    'origamis': {
        'name': 'surface_dynamics.flat_surfaces.origamis.origami_dense',
        'dir': os.path.join('flat_surfaces', 'origamis'),
        'sources': ['origami_dense.pyx', 'normal_form.c', 'lyapunov_exponents.c'],
        'headers': ['origami_dense.pxd', 'lyapunov_exponents.h', 'normal_form.h']
        },

    'lyapunov_exponents': {
        'name': 'surface_dynamics.interval_exchanges.lyapunov_exponents',
        'dir': 'interval_exchanges',
        'sources': ['lyapunov_exponents.pyx', 'generalized_permutation.c' , 'lin_alg.c', 'quad_cover.c', 'random.c', 'permutation.c'],
        'headers': ['lyapunov_exponents.h']
        },

    # this will not compile on sage < 8.0 (troubles with cysignals)
    'integer_iet': {
        'name': 'surface_dynamics.interval_exchanges.integer_iet',
        'dir': 'interval_exchanges',
        'sources': ['integer_iet.pyx', 'int_iet.c', 'int_vector.c'],
        'headers': ['integer_iet.pxd', 'int_iet.h'],
        },

    # build iet_family only if pplpy is available
    'iet_family': {
        'name': 'surface_dynamics.interval_exchanges.iet_family',
        'dir': 'interval_exchanges',
        'sources': ['iet_family.pyx'],
        'headers': [],
        'condition': WITH_PPL
        }
}

extensions = []
source_files = []

for name, data in extensions_data.items():
    if data.get('condition', True):
        print('Adding extension {}:\n  sources = {}\n  headers = {}'.format(data['name'], data['sources'], data['headers']))

        full_dir = os.path.join('surface_dynamics', data['dir'])
        sources = [os.path.join(full_dir, src) for src in data['sources']]
        headers = [os.path.join(full_dir, data['dir'], head) for head in data['headers']]
        #extensions.append(Extension(data['name'], sources=sources, include_dirs=[full_dir, np.get_include()]))

        sources = [os.path.join(data['dir'], src) for src in data['sources']]
        headers = [os.path.join(data['dir'], head) for head in data['headers']]
        source_files.extend(sources)
        source_files.extend(headers)

# Work around changes in SageMath 9.7, see https://trac.sagemath.org/wiki/ReleaseTours/sage-9.7#Packagessagesage.rings...arenownamespaces
try:
    from sage.misc.package_dir import cython_namespace_package_support
except (ImportError, ModuleNotFoundError):
    from contextlib import nullcontext as cython_namespace_package_support

with cython_namespace_package_support():
    extensions = cythonize(extensions)

setup(ext_modules=extensions)
