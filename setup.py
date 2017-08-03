#!/usr/bin/env bash
r"""
Installation script for the flatsurf module

It depends on distutils
"""

try:
    from sage.env import SAGE_SRC
except ImportError:
    raise ValueError("this package currently installs only inside SageMath (http://www.sagemath.org)")

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import sys, os

with open("surface_dynamics/version.py") as f:
    version = f.read().strip()
    prefix = "version='"
    suffix = "'"
    assert version.startswith(prefix) and version.endswith(suffix)
    version = version[len(prefix):len(version)-len(suffix)]
with open("README") as f:
    long_description = f.read()

ORIGAMIS_DIR = os.path.join('surface_dynamics', 'flat_surfaces', 'origamis')
LYAPUNOV_DIR = os.path.join('surface_dynamics', 'interval_exchanges', 'lyapunov_exponents')
INTEGER_IET_DIR = os.path.join('surface_dynamics', 'interval_exchanges', 'integer_iet')

extensions = [
    Extension('surface_dynamics.flat_surfaces.origamis.origami_dense',
            sources = [
            os.path.join(ORIGAMIS_DIR, filename) for filename in ('origami_dense.pyx', 'normal_form.c', 'lyapunov_exponents.c')],
            include_dirs = [SAGE_SRC, ORIGAMIS_DIR] + sys.path,
            libraries = ['m'],
            ),

    Extension('surface_dynamics.interval_exchanges.lyapunov_exponents',
            sources = [
                os.path.join(LYAPUNOV_DIR, filename) for filename in ('lyapunov_exponents.pyx', 'generalized_permutation.c' , 'lin_alg.c', 'quad_cover.c', 'random.c', 'permutation.c')],
            include_dirs = [SAGE_SRC, LYAPUNOV_DIR] + sys.path,
            depends = [os.path.join(LYAPUNOV_DIR, 'lyapunov_exponents.h')]),

    Extension('surface_dynamics.interval_exchanges.integer_iet',
        sources = [os.path.join(INTEGER_IET_DIR, 'int_iet.c'),
                   os.path.join(INTEGER_IET_DIR, 'int_vector.c'),
                   os.path.join(INTEGER_IET_DIR, 'integer_iet.pyx')],
        include_dirs = [SAGE_SRC, INTEGER_IET_DIR] + sys.path,
        depends = [os.path.join(INTEGER_IET_DIR, 'int_iet.h')])]

# build the iet family only if pplpy is available
try:
    import ppl
except ImportError:
    sys.stderr.write('Warning: pplpy not installed. Will not compile iet_family\n')
else:
    extensions.append(
    Extension('surface_dynamics.interval_exchanges.iet_family',
            sources = [os.path.join('surface_dynamics', 'interval_exchanges', 'iet_family.pyx')],
            include_dirs = [SAGE_SRC] + sys.path)
    )

setup(name='surface_dynamics',
      version=version,
      description="Dynamics on surfaces",
      long_description=long_description,
      author='Vincent Delecroix',
      author_email='vincent.delecroix@u-bordeaux.fr',
      url='http://www.labri.fr/perso/vdelecro/',
      license="GPL v3",
      packages=['surface_dynamics',
                'surface_dynamics/misc',
                'surface_dynamics/flat_surfaces',
                'surface_dynamics/databases',
                'surface_dynamics/flat_surfaces/origamis',
                'surface_dynamics/interval_exchanges'],
      package_data={
          'surface_dynamics/interval_exchanges': [
              'iet_family.pyx',
              'lyapunov_exponents/*.h',
              'lyapunov_exponents/*.c',
              'lyapunov_exponents/*.pyx',
              'lyapunov_exponents/*.pxd',
              'integer_iet/*.h',
              'integer_iet/*.c',
              'integer_iet/*.pyx',
              'integer_iet/*.pxd'],
          'surface_dynamics/databases': [
              'cylinder_diagrams/cyl_diags*',
              'generalized_permutation_twins/twins*'],
          'surface_dynamics/flat_surfaces/origamis': [
              '*.pyx',
              '*.pxd',
              'normal_form.h',
              'normal_form.c',
              'lyapunov_exponents.c',
              'lyapunov_exponents.h',
              'origamis.db'],
          },
      ext_modules=cythonize(extensions),
    classifiers=[
      'Development Status :: 4 - Beta',
      'Intended Audience :: Science/Research',
      'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
      'Operating System :: OS Independent',
      'Programming Language :: C',
      'Programming Language :: C++',
      'Programming Language :: Python',
      'Programming Language :: Cython',
      'Topic :: Scientific/Engineering :: Mathematics',
    ],
    keywords='surfaces, dynamics, geometry, flat surfaces, Abelian differentials, quadratic differentials, Riemann surfaces',
)
