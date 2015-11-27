#!/usr/bin/env bash
# setup file for the Sage module of translation surfaces
# this is called from the spkg-install script

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import Cython.Compiler.Options
import os

Cython.Compiler.Options.old_style_globals = True


from sage.env import SAGE_LOCAL, SAGE_SRC

setup(name='surface_dynamics',
      version='0.1',
      description="Dynamics on surfaces",
      author='Vincent Delecroix',
      author_email='vincent.delecroix@labri.fr',
      url='http://www.labri.fr/perso/vdelecro/',
      license = "GPL v3",
      packages=['surface_dynamics',
                'surface_dynamics/misc',
                'surface_dynamics/flat_surfaces',
                'surface_dynamics/databases',
                'surface_dynamics/flat_surfaces/origamis',
                'surface_dynamics/interval_exchanges'],
      ext_modules=[
        Extension('surface_dynamics.flat_surfaces.origamis.origami_dense',
                sources = [
                os.path.join('surface_dynamics','flat_surfaces','origamis', filename) for filename in ('origami_dense.pyx', 'normal_form.c', 'lyapunov_exponents.c')],
                libraries = ['m'],
                include_dirs = [
                      os.path.join('surface_dynamics','flat_surfaces','origamis'),
                      SAGE_SRC,
                      os.path.join(SAGE_SRC, 'sage', 'ext')
                      ],
    )], 
      cmdclass = {'build_ext': build_ext}
)

