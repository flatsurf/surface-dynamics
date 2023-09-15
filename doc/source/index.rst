Welcome to the surface-dynamics documentation!
==============================================

Installation
------------

The basic way to install the surface-dynamics package is via the pip
Python installer that comes with your `SageMath <http://www.sagemath.org/>`_
installation.

If you installed Sage from a downloaded binary (including the Cygwin version
on Windows) or if you compiled Sage from source, then run::

    $ sage -pip install surface-dynamics

You might want (or need) to provide the option ``--user`` to this command in order
to install surface-dynamics in your home directory and keep intact your SageMath
installation.

In case you use the sagemath package in a Debian or Ubuntu system you need to perform
a two step installation::

    $ source /usr/share/sagemath/bin/sage-env
    $ sage -pip install surface-dynamics --user

If you use the sagemath package from Archlinux, run::

    $ pip3 install git+https://github.com/flatsurf/surface-dynamics --user

Other sources of information:

- short tutorial: http://www.labri.fr/perso/vdelecro/flatsurf.html

- development webpage (including source code):
  https://github.com/flatsurf/surface-dynamics

Walk through
------------

.. toctree::
    :maxdepth: 1

    examples/interval_exchanges
    examples/square_tiled_surfaces
    examples/rank2_genus3_classification

Module documentation
--------------------

.. toctree::
    :maxdepth: 2

    strata
    surface_topology
    origamis
    interval_exchanges/index
    database
    topological_recursion/index
    misc/index

References
----------

.. toctree::
    :maxdepth: 2

    references

Citation
--------

To cite the library, follow the instructions on our `zenodo website <https://zenodo.org/badge/latestdoi/347440823>`_.

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
