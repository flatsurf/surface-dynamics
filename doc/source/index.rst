Welcome to surface_dynamic's documentation!
===========================================

Installation
------------

The basic way to install the surface_dynamic package is via the pip
Python installer that comes with your `SageMath <http://www.sagemath.org/>`_
installation.

If you installed Sage from a downloaded binary (included the Cygwin version
on Windows) or compiled Sage from sources then run::

    $ sage -pip install surface_dynamics

You might want (or need) to provide the option ``--user`` to this command in order
to install surface dynamics in your home directory and keep intact your SageMath
installation.

In case you use the sagemath package in a Debian or Ubuntu system you need to perform
a two step installation::

    $ source /usr/share/sagemath/bin/sage-env
    $ sage -pip install surface_dynamics --user

If you use the sagemath package from Archlinux run::

    $ pip2 install git+https://gitlab.com/videlec/surface_dynamics --user

Other sources of information:

- there is a short tutorial at http://www.labri.fr/perso/vdelecro/flatsurf.html

- development webpage (including source code) is available at
  https://gitlab.com/videlec/surface_dynamics

Module documentation
--------------------

.. toctree::
    :maxdepth: 2

    strata
    surface_topology
    origamis
    interval_exchanges
    database

Citation
--------

To cite the library, follow the following Bibtex entry::

    @manual{ Sdyn,
       Author = { Delecroix, V. et al. },
       Month  = { March },
       Year   = { 2019 },
       Title  = { surface_dynamics - SageMath package, Version 0.4.1 },
       Doi    = { 10.5281/zenodo.3237923 },
       Url    = { https://doi.org/10.5281/zenodo.3237923 }
    }

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

