**Changed:**

* Changed the class of the flip sequences built by ``iet.FlipSequence`` to
  ``RauzyFlipSequence``. It serves generalized permutations (linear
  involutions) as well as permutations of interval exchanges; the
  substitutions and ``self_similar_iet`` raise ``NotImplementedError`` on
  linear involutions.

* Changed the Rauzy moves of a flip sequence: with ``check=True`` (the
  default), a move that is not defined raises a ``ValueError`` and leaves
  the flip sequence unchanged.

**Deprecated:**

* Deprecated ``IETFlipSequence``: building one warns and returns a
  ``RauzyFlipSequence``; flip sequences pickled under this name load as
  ``RauzyFlipSequence``.

**Fixed:**

* Fixed ``reversed`` on flip sequences of linear involutions and of
  flipped permutations.
