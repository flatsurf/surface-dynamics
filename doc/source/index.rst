surface-dynamics: Dynamics with Surfaces
========================================

The surface-dynamics package for `SageMath <https://www.sagemath.org>`_
provides functionality related to interval exchange transformations,
translation surfaces, mapping classes and more. It is based on SageMath and
relies heavily on:

* `GMP <https://gmplib.org/>`_ for arbitrary precision arithmetic
* `PARI/GP <https://pari.math.u-bordeaux.fr/>`_ for number field computations
* `GAP <https://www.gap-system.org/>`_ for finite groups representation and permutation groups
* `PPL <https://www.bugseng.com/ppl>`_ (Parma Polyhedra Library) and
  `LattE <https://www.math.ucdavis.edu/~latte/>`_ (Lattice point Enumeration) for
  polytope computations

The development of surface-dynamics is coordinated on GitHub at https://github.com/flatsurf/surface-dynamics.

Installation
------------

The preferred way to install software should be to use your package manager
(e.g. ``apt`` on Debian or Ubuntu, ``pacman`` on Arch Linux, ``brew`` on macOS,
etc). However, as of this writing, surface-dynamics has not been picked up by `any
of the major distributions yet <https://repology.org/project/python:surface-dynamics/packages>`_.

We therefore recommend to install the entire flatsurf suite from the
sage-flatsurf `Releases Page
<https://github.com/flatsurf/sage-flatsurf/releases>`_.

Detailed installation instructions:

* :external+sage-flatsurf:ref:`Install the flatsurf suite with the pixi tarball <installation-tarball>` for Linux and macOS (recommended)
* :external+sage-flatsurf:ref:`Install the flatsurf suite with the installer <installation-installer>` for Windows (recommended)
* :external+sage-flatsurf:ref:`Install with Conda <installation-conda>`
* :ref:`Install into an existing SageMath source build <installation-sagemath>`
* :ref:`Install with pip <installation-pip>`

If you are planning to work on the surface-dynamics source code, you can also
build surface-dynamics from source. For this, please have a look at our
:ref:`Developer's Guide <developers-guide>`.

.. toctree::
   :hidden:

   install
   developer

A Tour of surface-dynamics
--------------------------

Demos of some of the capabilities of surface-dynamics:

.. toctree::
   :maxdepth: 1
   :class: bullet-toc

   examples/interval_exchanges
   examples/square_tiled_surfaces
   examples/rank2_genus3_classification

Module Reference
----------------

The surface-dynamics source code is split into several Python modules. The
links below lead to the documentation for these modules with more usage
examples and the source code for all the classes and functions implemented in
surface-dynamics. (Note that you can also access this documentation from an
interactive SageMath session with |help|_.)

.. |help| replace:: ``?`` and ``??``
.. _help: https://ipython.readthedocs.io/en/stable/interactive/python-ipython-diff.html#accessing-help

.. toctree::
   :maxdepth: 1
   :class: bullet-toc

   databases/flat_surfaces
   flat_surfaces
   flat_surfaces/origamis
   interval_exchanges
   misc
   topological_recursion
   topology

Citation
--------

To cite surface-dynamics, please follow the instructions on our `zenodo website <https://zenodo.org/badge/latestdoi/347440823>`_.

See our :doc:`references` for everything that is cited in this documentation.

.. toctree::
   :hidden:

   references
