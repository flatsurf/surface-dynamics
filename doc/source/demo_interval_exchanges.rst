.. -*- coding: utf-8 *-*
.. linkall

Interval exchange transformations
=================================

.. contents::
   :depth: 1

Introduction
------------

This file is a good entry point to learn how to use interval exchange
transformations in `surface-dynamics
<https://github.com/flatsurf/surface-dynamics>`_. Recall that an interval
exchange transformation is a piecewise translation of an interval. It can be
encoded by a pair :math:`(\pi, \lambda)` where :math:`\pi` is a permutation and
:math:`\lambda` is a vector of positive real numbers.  These are respectively
called the *combinatorial datum* and the *length datum* of the interval
exchange transformation.

Building an interval exchange transformation and simple manipulations
---------------------------------------------------------------------

Here is an example of interval exchange transformation on 4 subintervals
with rational lengths::

    sage: from surface_dynamics import iet
    sage: perm = iet.Permutation('a b c d', 'd c b a')
    sage: length = [1/3, 2/7, 4/13, 5/19]
    sage: t = iet.IntervalExchangeTransformation(perm, length)
    sage: print(t)
    Interval exchange transformation of [0, 6172/5187[ with permutation
    a b c d
    d c b a

It can be visualized with::

    sage: G = t.plot_function()
    sage: G.set_aspect_ratio(1)
    sage: G # random (matplotlib warning in conda)
    Graphics object consisting of 4 graphics primitives
    sage: t.plot_two_intervals()
    Graphics object consisting of 16 graphics primitives

Given a point in the interval :math:`[0, 6172/5187[` it is possible to compute
its image under the map :math:`t`::

    sage: x = 1/12
    sage: t(x)
    19501/20748

To know which translation has been applied to the point :math:`x` you can
use::

    sage: t.in_which_interval(x)
    'a'
    sage: t.translations()
    (1481/1729, 176/741, -142/399, -253/273)
    sage: x + t.translations()[0] == t(x)
    True

Can you compute the first 20 elements of the orbit of :math:`x`, that is the
sequence :math:`(x, t(x), t^2(x), \ldots, t^{19}(x))`? ::

    sage: # edit here

Can you determine whether the orbit of :math:`x` is periodic, that is whether
there exists :math:`n > 0` such that :math:`t^n(x) = x`? ::

    sage: # edit here

To construct an interval exchange on other numbers than rationals you need
to manipulate exact real numbers. The current sage version (9.2) only supports
computation with algebraic numbers such as ::

    sage: x = polygen(QQ)
    sage: K.<sqrt2> = NumberField(x^2 - 2, embedding=AA(2).sqrt())
    sage: length2 = [1, sqrt2, sqrt2**2, sqrt2**3]
    sage: t2 = iet.IntervalExchangeTransformation(perm, length2)
    sage: print(t2)
    Interval exchange transformation of [0, 3*sqrt2 + 3[ with permutation
    a b c d
    d c b a

Rauzy induction and self-similar iet
------------------------------------

Rauzy induction is a map from the space of interval exchange transformations to itself.
The image :math:`\mathcal{R}(t)` of an iet :math:`t` is an induced map.::

    sage: t3 = t2.rauzy_move()
    sage: r = max(set(flatten(t2.singularities())) - set([0,3*sqrt2+3]))
    sage: G = t2.plot_two_intervals() + t3.plot_two_intervals(position=(0,-.8))
    sage: G += line2d([(r, -1.), (r, .2)], color='black', linestyle='dotted')
    sage: G.axes(False)
    sage: G
    Graphics object consisting of 33 graphics primitives

For the next demo, we need a bit of setup to silence an error in automated
testing. You don't need to run these commands unless you are using a very old
version of SageMath::

    sage: import numpy
    sage: numpy.float = float

From a flip sequence given combinatorially you can build the associated self-similar
interval exchange transformation::

    sage: R = perm.rauzy_diagram()
    sage: G = R.graph()
    sage: P = G.plot(color_by_label=True, edge_labels=True, vertex_size=800)
    sage: P
    Graphics object consisting of 42 graphics primitives

::

    sage: seq = ['t', 'b', 't', 'b', 't', 't', 'b', 'b', 't', 'b']
    sage: f = iet.FlipSequence(perm, seq)
    sage: print(f.is_loop())
    True
    sage: print(f.is_full())
    True

::

    sage: dilatation, t4 = f.self_similar_iet()
    sage: print(dilatation, '~', dilatation.n())
    3*a + 2 ~ 6.85410196624968

Above ``dilatation`` is the expansion of the self-similarity and ``t4`` is the self-similar
exchange transformation associated to the flip sequence ``f``::

    sage: t5 = t4.rauzy_move(iterations=len(seq))
    sage: G = t4.plot() + t5.plot(position=(0,-.5))
    sage: G.axes(False)
    sage: G
    Graphics object consisting of 32 graphics primitives

::

    sage: print(t4.lengths())
    (1, 6/5*a + 2/5, 3/5*a + 1/5, 6/5*a + 2/5)
    sage: print(t5.lengths())
    (-3*a + 5, 6/5*a - 8/5, 3/5*a - 4/5, 6/5*a - 8/5)
    sage: dilatation * t5.lengths() == t4.lengths()
    True

The command below checks that ``t4`` is indeed self induced::

    sage: t4.normalize() == t5.normalize()
    True

Suspension
----------

`sage-flatsurf <https://flatsurf.github.io/sage-flatsurf/>`_ is a Python library for translation
surfaces (and more generally similarity surfaces). One can build Masur polygons via::

    sage: height = [1, 0, 0, -1]
    sage: S = perm.masur_polygon(length2, height) # optional - sage_flatsurf
    sage: S # optional - sage_flatsurf
    TranslationSurface built from 6 polygons
    sage: S.stratum() # optional - sage_flatsurf
    H_2(2)

Could you construct a self-similar translation surface from the flip sequence ``f``? (in other words
a translation surface that admits a pseudo-Anosov preserving the horizontal and vertical
foliations)::

    sage: # edit here

Using pyintervalxt
------------------

`intervalxt <https://github.com/flatsurf/intervalxt>`_ is a C++ library with a Python interface
that implements optimized routines to deal with interval exchange
transformations. If ``intervalxt`` is part of your installation you can convert
interval exchange transformations back and forth between ``surface-dynamics``
and ``pyintervalxt``::

    sage: from surface_dynamics.interval_exchanges.conversion import iet_to_pyintervalxt, iet_from_pyintervalxt # optional - gmpxxyy, pyintervalxt
    sage: u2 = iet_to_pyintervalxt(t2) # optional - gmpxxyy, pyintervalxt
    sage: print(u2) # optional - gmpxxyy, pyintervalxt
    [a: 1] [b: (sqrt2 ~ 1.4142136)] [c: 2] [d: (2*sqrt2 ~ 2.8284271)] / [d] [c] [b] [a]
    sage: v2 = iet_from_pyintervalxt(u2) # optional - gmpxxyy, pyintervalxt
    sage: print(v2) # optional - gmpxxyy, pyintervalxt
    Interval exchange transformation of [0, 3*sqrt2 + 3[ with permutation
    a b c d
    d c b a

One feature of ``intervalxt`` is that it can certify that an iet has no periodic trajectory::

    sage: u2.boshernitzanNoPeriodicTrajectory() # optional - gmpxxyy, pyintervalxt
    True

Other features
--------------

This short tour did not exhaust all the possibilities of ``surface-dynamics``, in particular

- iet statistics :mod:`~surface_dynamics.interval_exchanges.integer_iet`

- linear families of iet :mod:`~surface_dynamics.interval_exchanges.iet_family`

- coverings and Lyapunov exponents of the Kontsevich-Zorich cocycle

These topics might be included in later versions of this document.

Missing features
----------------

- generalizations (linear involution associated to generalized permutations,
  interval exchange transformations with flips, affine iet, system of isometries)

- Veech zippered rectangle construction

- constructing the self-similar surface (aka pseudo-Anosov) associated to a
  flip sequence

If you are interested in developing any of these or have any request, get in
touch with us at https://github.com/flatsurf/surface-dynamics
