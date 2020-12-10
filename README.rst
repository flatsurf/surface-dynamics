================
Surface dynamics
================

The ``surface_dynamics`` package for SageMath adds functionality related to
interval exchange transformations, translation surfaces, mapping classes
and more. It is based on `SageMath <http://www.sagemath.org>`_ and relies
heavily on:

* gmp or mpir for arbitrary precision arithmetic
* PARI/GP for number field computations
* GAP for finite groups representation and permutation groups
* PPL (Parma Polyhedra Library) and LattE (Lattice point Enumeration)
  for polytope computations

Prerequisites
-------------

Installing ``surface_dynamics`` requires a working Sage installation (with
Cython and gcc). Installing the optional SageMath packages ``gap_packages``,
and ``latte_int`` is recommended and will improve or extend the functionality
in ``surface_dynamics``. The optional package ``database_gap`` is also
recommended if using SageMath < 8.6 (in SageMath 8.6 it was merged partly
into the ``gap`` and partly into the ``gap_packages`` packages).

Installation
------------

The module is distributed on PyPI and is easily installed through the
Python package manager ``pip``. If you downloaded a binary from the SageMath
website (including the Cygwin version running on Windows) or compiled
from source, run the following command::

    $ sage -pip install surface_dynamics [--user]

The ``--user`` option is optional and allows to install the module in your
user space (and does not require administrator rights).

If you use Debian or Ubuntu and you installed Sage through the operating
system's package manager (that is, the package ``sagemath``), run these
two commands::

    $ source /usr/share/sagemath/bin/sage-env
    $ pip install surface_dynamics --user

If you use Arch Linux, you need to install from source (see next section).

Install and use source version
------------------------------

This section provides detailed instructions on how to download, modify
and install the development version of ``surface_dynamics``. In all commands,

* ``PIP`` has to be replaced by either ``pip``, ``pip2``, or ``sage -pip``
* ``PYTHON`` has to be replaced by either ``python``, ``python2`` or ``sage -python``

If you are an Arch Linux user with the ``sagemath`` package installed, use
``PIP=pip2`` and ``PYTHON=python2``. If you downloaded SageMath as a tarball
or installed it from source use ``PIP='sage -pip'`` and ``PYTHON='sage -python'``.

You can install the latest development version in one line with::

    $ PIP install git+https://gitlab.com/videlec/surface_dynamics [--user]

As before, the ``--user`` option is optional and when specified will
install the module in your user space.

You can also perform a two stage installation that will allow you to
modify the source code. The first step is to clone the repository::

    $ git clone https://gitlab.com/videlec/surface_dynamics

The above command creates a repository ``surface_dynamics`` with the source code,
documentation and miscellaneous files. You can then change to the directory
thus created and install the surface dynamics module with::

    $ cd surface_dynamics
    $ PIP install . [--user]

Do not forget the ``.`` that refers to the current directory.

When you don't want to install the package or you are testing some
modifications to the source code, a more convenient way of using
surface dynamics is to do everything locally. To do so, you need
to compile the module in place via::

    $ PYTHON setup.py build_ext --inplace

Once done, you can import the ``surface_dynamics`` module. To check that you
are actually using the right module (i.e. the local one) you can do in a
SageMath session::

    sage: import surface_dynamics
    sage: surface_dynamics.__path__        # random
    ['/home/you/surface_dynamics/surface_dynamics/']

The result of the command must correspond to the path of the repository
created by the command ``git clone`` given above. The compilation step
``PYTHON setup.py build_ext`` has to be redone each time you modify
a C or Cython source file (i.e. with ``.c``, ``.h``, ``.pxd`` or ``.pyx``
extension). In other words, it is not needed if you only
modify or create Python files (i.e. ``.py`` files).

If you wish to install your custom version of ``surface_dynamics``
just use ``PIP`` as indicated before.

Documentation
-------------

* short tutorial: http://www.labri.fr/perso/vdelecro/flatsurf.html
* complete module documentation: http://www.labri.fr/perso/vdelecro/surface-dynamics/

Check
-----

After installing ``surface_dynamics``, check that it works by launching Sage
and typing the following commands. You should get the same
output as below. ::

    sage: from surface_dynamics.all import *
    sage: o = Origami('(1,2)', '(1,3)')
    sage: o
    (1,2)(3)
    (1,3)(2)
    sage: o.sum_of_lyapunov_exponents()
    4/3
    sage: o.lyapunov_exponents_approx()    # abs tol 0.05
    [0.33441823619678734]
    sage: o.veech_group()
    Arithmetic subgroup with permutations of right cosets
     S2=(2,3)
     S3=(1,2,3)
     L=(1,2)
     R=(1,3)
    sage: q = QuadraticStratum(1, 1, 1, 1)
    sage: q.orientation_cover()
    H_5(2^4)
    sage: q.components()
    [Q_2(1^4)^hyp]
    sage: c = q.components()[0]
    sage: c
    Q_2(1^4)^hyp
    sage: c.orientation_cover_component()
    H_5(2^4)^odd

    sage: AbelianStrata(genus=3).list()
    [H_3(4), H_3(3, 1), H_3(2^2), H_3(2, 1^2), H_3(1^4)]

    sage: O = OrigamiDatabase()
    sage: q = O.query(("stratum", "=", AbelianStratum(2)), ("nb_squares", "=", 5))
    sage: q.number_of()
    2
    sage: for o in q:
    ....:     print("%s\n- - - - - - - -" % o)
    (1)(2)(3)(4,5)
    (1,2,3,4)(5)
    - - - - - - - -
    (1)(2)(3,4,5)
    (1,2,3)(4)(5)
    - - - - - - - -
    sage: Q12_reg = QuadraticStratum(12).regular_component()
    sage: Q12_reg.lyapunov_exponents_H_plus(nb_iterations=2**20)   # abs tol 0.05
    [0.6634, 0.4496, 0.2305, 0.0871]
    sage: Q12_reg.lyapunov_exponents_H_minus(nb_iterations=2**20)  # abs tol 0.05
    [1.0000, 0.3087, 0.1192]

Installing development version - source code
--------------------------------------------

The development webpage is

* https://gitlab.com/videlec/surface_dynamics

Assuming you have the program ``git`` on your computer, you can install the
development version with the command::

    $ sage -pip install git+https://gitlab.com/videlec/surface_dynamics [--user]

Contact
-------

Your comments and help are welcome: vincent.delecroix@labri.fr

For problems with macOS: samuel.lelievre@gmail.com

Authors
-------

* Vincent Delecroix: maintainer
* Samuel Lelièvre: origami and permutation representatives for quadratic strata
* Charles Fougeron: Lyapunov exponents for strata coverings
* Luke Jeffreys: single cylinder representatives for strata of Abelian
  differentials

Citation
--------

To cite the software, use the following Bibtex entry::

    @manual{ Sdyn,
      Author = { Delecroix, V. et al. },
      Month  = { March },
      Year   = { 2019 },
      Title  = { surface_dynamics - SageMath package, Version 0.4.1 },
      Doi    = { 10.5281/zenodo.3237923 },
      Url    = { https://doi.org/10.5281/zenodo.3237923 }
    }

Versions
--------

* flatsurf 0.4.4 was released on 2020-01-31 (as a Python package on PyPI)
* flatsurf 0.4.3 was released on 2019-07-28 (as a Python package on PyPI)
* flatsurf 0.4.2 was released on 2019-06-21 (as a Python package on PyPI)
* flatsurf 0.4.1 was released on 2019-03-26 (as a Python package on PyPI)
* flatsurf 0.4.0 was released on 2018-05-14 (as a Python package on PyPI)
* flatsurf 0.3 was released on 2017-08-11 (as a Python package on PyPI)
* flatsurf 0.2 was released on 2015-11-15 (as a Sage spkg)
* flatsurf 0.1 was released on 2015-07-30 (as a Sage spkg)