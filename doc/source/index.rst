Welcome to the surface_dynamics documentation!
==============================================

Installation
------------

The basic way to install the surface_dynamics package is via the pip
Python installer that comes with your `SageMath <http://www.sagemath.org/>`_
installation.

If you installed Sage from a downloaded binary (including the Cygwin version
on Windows) or if you compiled Sage from source, then run::

    $ sage -pip install surface_dynamics

You might want (or need) to provide the option ``--user`` to this command in order
to install surface_dynamics in your home directory and keep intact your SageMath
installation.

In case you use the sagemath package in a Debian or Ubuntu system you need to perform
a two step installation::

    $ source /usr/share/sagemath/bin/sage-env
    $ sage -pip install surface_dynamics --user

If you use the sagemath package from Archlinux, run::

    $ pip3 install git+https://github.com/flatsurf/surface_dynamics --user

Other sources of information:

- short tutorial: http://www.labri.fr/perso/vdelecro/flatsurf.html

- development webpage (including source code):
  https://github.com/flatsurf/surface_dynamics

Walk through
------------

.. toctree::
    :maxdepth: 1

    demo_interval_exchanges

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
