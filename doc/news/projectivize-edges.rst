**Fixed:**

* Fixed compilation with Cython 3.3, which rejected a list assigned to a
  variable declared as a ``dict`` in ``origami_dense.pyx``.
* Fixed ``Origami.pgl2z_edges()`` and ``Origami.psl2z_edges()`` which raised a
  ``TypeError`` whenever ``-Id`` does not fix the origami.

**Changed:**

* Restricted the build requirement to ``Cython<3.4`` so that a new Cython
  release cannot break the build of a released version.
