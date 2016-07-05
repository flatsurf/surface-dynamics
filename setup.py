#!/usr/bin/env bash
r"""
Installation script for the flatsurf module

It depends on distutils
"""

try:
    from sage.env import SAGE_LOCAL, SAGE_SHARE, SAGE_DOC, SAGE_SRC, SAGE_VERSION
except ImportError:
    raise ValueError("this package currently installs only inside SageMath (http://www.sagemath.org)")

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os

from surface_dynamics.version import version
with open("README") as f:
    long_description = f.read()

ORIGAMIS_DIR = os.path.join('surface_dynamics', 'flat_surfaces', 'origamis')
LYAPUNOV_DIR = os.path.join('surface_dynamics', 'interval_exchanges', 'lyapunov_exponents')

extensions = [
    Extension('surface_dynamics.flat_surfaces.origamis.origami_dense',
            sources = [
            os.path.join(ORIGAMIS_DIR, filename) for filename in ('origami_dense.pyx', 'normal_form.c', 'lyapunov_exponents.c')],
            libraries = ['m'],
            include_dirs = [ORIGAMIS_DIR, SAGE_SRC],
            ),

    Extension('surface_dynamics.interval_exchanges.lyapunov_exponents',
        sources = [
            os.path.join(LYAPUNOV_DIR, filename) for filename in ('lyapunov_exponents.pyx', 'generalized_permutation.c' , 'lin_alg.c', 'quad_cover.c', 'random.c', 'permutation.c')],
        depends = [os.path.join(LYAPUNOV_DIR, 'lyapunov_exponents.h')])
    ]


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
              'lyapunov_exponents/lyapunov_exponents.h'],
          'surface_dynamics/databases': [
              'cylinder_diagrams/*',
              'generalized_permutation_twins/*'],
          'surface_dynamics/flat_surfaces/origamis': [
              'origamis.db',
              'origami_dense.pxd',
              'normal_form.h',
              'lyapunov_exponents.h'],
          },
      ext_modules=extensions,
    cmdclass={'build_ext': build_ext},
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
