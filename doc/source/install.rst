Installing surface-dynamics
===========================

There are several different ways to install sage-flatsurf on your machine. You can

* extract :external+sage-flatsurf:ref:`our pixi tarball <installation-tarball>` for Linux and macOS (recommended),
* use :external+sage-flatsurf:ref:`our executable installer <installation-installer>` on Windows (recommended),
* create a :external+sage-flatsurf:ref:`conda environment <installation-conda>` on Linux and macOS,
* install surface-dynamics `into an existing source build of SageMath <#installation-sagemath>`_,
* or `pip install <#installation-pip>`_ surface-dynamics.

If you are having trouble with this or are on another operating system, please
`contact us <https://flatsurf.github.io>`_. We're thrilled about any new user
of surface-dynamics and we're very happy to help and in the process improve
these installation instructions.

.. _installation-sagemath:

Install into SageMath
---------------------

If you are using a `source build of SageMath
<https://doc.sagemath.org/html/en/installation/source.html>`_ or if you
downloaded a SageMath binary, you can install surface-dynamics into SageMath.
Note that this does not install all the optional dependencies of
surface-dynamics so some computations might fail in this setup::

        sage -pip install surface-dynamics

To uninstall surface-dynamics again later::

        sage -pip uninstall surface-dynamics

.. _installation-pip:

Install from PyPI
-----------------

You can install surface-dynamics from `PyPI
<https://pypi.org/project/surface-dynamics/>`_ if you installed sagelib as a
Python package. Again, this does not come with the optional dependencies, so
some computations might fail in this setup::

        pip install --user surface-dynamics

To uninstall surface-dynamics again later::

        pip uninstall --user surface-dynamics
