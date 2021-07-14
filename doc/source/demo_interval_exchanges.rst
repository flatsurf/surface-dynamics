.. linkall:

Interval exchange transformations
=================================

This file is a good entry point to learn how to use interval exchange
transformations in `surface_dynamics
<https://github.com/flatsurf/surface_dynamics>`_. Recall that an interval
exchange transformation is a piecewise translation of an interval. It can be
encoded by a pair `(\pi, \lambda)` where `\pi` is a permutation and `\lambda`
is a vector of positive real numbers.  These are respectively called the
*combinatorial datum* and the *length datum* of the interval echange
transformation.

Building an interval exchange transformation and simple manipulations
---------------------------------------------------------------------

Here is a simple interval exchange transformation on 4 subintervals
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

    sage: t.plot_function().show(aspect_ratio=1)
    sage: t.plot_two_intervals().show(axes=False)

Given a point in the interval `[0, 6172/5187[` it is possible to compute
its image under the map `t`::

    sage: x = 1/12
    sage: t(x)

Can you compute the first 20 elements of the orbit of `x`, that is the
list `[x, t(x), t^2(x), \ldots, t^{19}(x)]`? ::

    sage: # edit here

Can you determine whether the orbit of `x` is periodic, that is whether
there exists `n > 0` such that `t^n(x) = x`? ::

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
The image `\mathcal{R}(t)` of an iet `t` is an induced map.::

    sage: t3 = t2.rauzy_move()
    sage: r = max(set(flatten(t2.singularities())) - set([0,3*sqrt2+3]))
    sage: G = t2.plot_two_intervals() + t3.plot_two_intervals(position=(0,-.8))
    sage: G += line2d([(r, -1.), (r, .2)], color='black', linestyle='dotted')
    sage: G.show(axes=False)

From a flip sequence given combinatorially you can build the associated self-similar
interval exchange transformation::

    sage: R = perm.rauzy_diagram()
    sage: R.graph().plot(color_by_label=True, edge_labels=True, vertex_size=800).show(figsize=12)

::

    sage: seq = ['t', 'b', 't', 'b', 't', 't', 'b', 'b', 't', 'b']
    sage: f = iet.FlipSequence(perm, seq)
    sage: print(f.is_loop())
    True
    sage: print(f.is_full())
    True

::

    sage: a, t4 = f.self_similar_iet()
    sage: t5 = t4.rauzy_move(iterations=len(seq))
    sage: a * t5.length() == t4.length()
    True

::

    sage: (t4.plot() + t5.plot(position=(0,-.5))).show(axes=False)

The command below checks that it is indeed self induced::

    sage: t4.normalize() == t5.normalize()
    True

Suspension
----------

`sage-flatsurf <https://flatsurf.github.io/sage-flatsurf/>`_ is a Python library for translation
surfaces (and more generally similarity surfaces). ::

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
that implements optimized routines to deal with interval exchange transformations. If it is part
of your installation you can convert interval exchange transformations back and forth between
``surface_dynamics`` and ``pyintervalxt``::

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

Missing features
----------------

- generalizations (linear involution associated to generalized permutations,
  interval exchange transformations with flips, affine iet)

- Veech zippered rectangle construction

- constructing the self-similar surface (aka pseudo-Anosov) associated to a
  flip sequence

If you are interested in developing any of these, get in touch with us at
https://github.com/flatsurf/surface_dynamics
