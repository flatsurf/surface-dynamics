r"""
Conversion to and from pyintervalxt
"""

from sage.all import NumberField, QQ, ZZ
from sage.rings.number_field.number_field_base import is_NumberField

from . import constructors as iet
from .iet import IntervalExchangeTransformation

def iet_to_pyintervalxt(T):
    r"""
    Return the pyintervalxt version of the interval exchange transformation ``T``.

    EXAMPLES::

        sage: from surface_dynamics import iet
        sage: from surface_dynamics.interval_exchanges.conversion import iet_to_pyintervalxt
        sage: p = iet.Permutation("a b c", "c b a")

    Over integers::

        sage: T = iet.IntervalExchangeTransformation(p, [12, 5, 9])
        sage: T.base_ring()
        Integer Ring
        sage: iet_to_pyintervalxt(T) # optional - gmpxxyy, pyintervalxt
        [a: 12] [b: 5] [c: 9] / [c] [b] [a]

    Over rationals::

        sage: T = iet.IntervalExchangeTransformation(p, [12/5, 2/7, 1/21])
        sage: T.base_ring()
        Rational Field
        sage: iet_to_pyintervalxt(T) # optional - gmpxxyy, pyintervalxt
        [a: 12/5] [b: 2/7] [c: 1/21] / [c] [b] [a]

    Over number fields::

        sage: x = polygen(QQ)
        sage: K.<sqrt2> = NumberField(x^2 - 2, embedding=AA(2).sqrt())
        sage: T = iet.IntervalExchangeTransformation(p, [1, sqrt2, sqrt2-1])
        sage: iet_to_pyintervalxt(T) # optional - pyeantic, pyintervalxt
        [a: 1] [b: (sqrt2 ~ 1.4142136)] [c: (sqrt2-1 ~ 0.41421356)] / [c] [b] [a]
    """
    try:
        # instantiate some voudou magic
        import gmpxxyy
    except ImportError:
        pass

    if not isinstance(T, IntervalExchangeTransformation):
        raise TypeError
    K = T._base_ring
    if K == ZZ:
        from gmpxxyy import mpz
        lengths = tuple(mpz(x) for x in T._lengths)
    elif K == QQ:
        from gmpxxyy import mpq
        lengths = tuple(mpq(x.numerator(), x.denominator()) for x in T._lengths)
    elif is_NumberField(K):
        from pyeantic.real_embedded_number_field import RealEmbeddedNumberField
        L = RealEmbeddedNumberField(K)
        lengths = tuple(L(x).renf_elem for x in T._lengths)
    else:
        raise NotImplementedError("unknown base ring K={}".format(K))

    perm = tuple(T._permutation._twin[1])

    import pyintervalxt
    return pyintervalxt.IntervalExchangeTransformation(lengths, perm)


def lengths_to_sage(lengths):
    try:
        # NOTE: the module is not used anywhere but importing it triggers some
        # needed coercions
        import gmpxxyy
    except ImportError:
        pass

    if isinstance(lengths[0], int):
        return [ZZ(x) for x in lengths]

    try:
        from cppyy.gbl import mpz_class
    except ImportError:
        pass
    else:
        if isinstance(lengths[0], mpz_class):
            return [ZZ(x) for x in lengths]

    try:
        from cppyy.gbl import mpq_class
    except ImportError:
        pass
    else:
        if isinstance(lengths[0], mpq_class):
            return [QQ(x) for x in lengths]

    try:
        from cppyy.gbl import eantic
    except ImportError:
        pass
    else:
        if isinstance(lengths[0], eantic.renf_elem_class):
            from pyeantic.real_embedded_number_field import RealEmbeddedNumberField
            K = RealEmbeddedNumberField(lengths[0].parent())
            K_sage = K.number_field
            return [K_sage(K(x)) for x in lengths]

    raise NotImplementedError('unknown length type {}'.format(type(lengths[0])))

def iet_from_pyintervalxt(T, alphabet=None):
    r"""
    Conversion from pyintervalxt

    EXAMPLES::

        sage: import pyintervalxt # optional - gmpxxyy, pyintervalxt
        sage: from surface_dynamics.interval_exchanges.conversion import iet_from_pyintervalxt
        sage: perm = (int(1), int(0))
        sage: T = pyintervalxt.IntervalExchangeTransformation((int(18), int(3)), perm) # optional - gmpxxyy, pyintervalxt
        sage: iet_from_pyintervalxt(T) # optional - gmpxxyy, pyintervalxt
        Interval exchange transformation of [0, 21[ with permutation
        a b
        b a

        sage: from cppyy.gbl import mpz_class # optional - gmpxxyy, pyintervalxt
        sage: l1 = mpz_class(12384758375127328356724597182479485) # optional - gmpxxyy, pyintervalxt
        sage: l2 = mpz_class(571349847513463874558781940004928374) # optional - gmpxxyy, pyintervalxt
        sage: T = pyintervalxt.IntervalExchangeTransformation((l1, l2), perm) # optional - gmpxxyy, pyintervalxt
        sage: iet_from_pyintervalxt(T) # optional - gmpxxyy, pyintervalxt
        Interval exchange transformation of [0, 583734605888591187546058532203790336[ with permutation
        a b
        b a

        sage: from cppyy.gbl import mpq_class # optional - gmpxxyy, pyintervalxt
        sage: l1 = mpq_class(1, 123847585) # optional - gmpxxyy, pyintervalxt
        sage: l2 = mpq_class(4928374, 3) # optional - gmpxxyy, pyintervalxt
        sage: T = pyintervalxt.IntervalExchangeTransformation((l1, l2), perm) # optional - gmpxxyy, pyintervalxt
        sage: iet_from_pyintervalxt(T) # optional - gmpxxyy, pyintervalxt
        Interval exchange transformation of [0, 610367217876793/371542755[ with permutation
        a b
        b a
        
        sage: from pyeantic.real_embedded_number_field import RealEmbeddedNumberField # optional - pyeantic, pyintervalxt
        sage: x = polygen(QQ)
        sage: K.<sqrt2> = NumberField(x^2 - 2, embedding=AA(2).sqrt())
        sage: L = RealEmbeddedNumberField(K) # optional - pyeantic, pyintervalxt
        sage: l1 = L(sqrt2).renf_elem # optional - pyeantic, pyintervalxt
        sage: l2 = L(4*sqrt2 - 3).renf_elem # optional - pyeantic, pyintervalxt
        sage: T = pyintervalxt.IntervalExchangeTransformation((l1, l2), perm) # optional - pyeantic, pyintervalxt
        sage: iet_from_pyintervalxt(T) # optional - pyeantic, pyintervalxt
        Interval exchange transformation of [0, 5*sqrt2 - 3[ with permutation
        a b
        b a
    """
    import pyintervalxt
    from cppyy.gbl import intervalxt

    if not isinstance(T, intervalxt.IntervalExchangeTransformation):
        raise TypeError("input must be an interval exchange transformation from pyintervalxt")

    top = [str(x) for x in T.top()]
    bottom = [str(x) for x in T.bottom()]
    p = iet.Permutation(top, bottom)
    if alphabet is not None:
        p.alphabet(alphabet)
    lengths = lengths_to_sage([T.lengths().get(x) for x in T.top()])

    return iet.IntervalExchangeTransformation(p, lengths)
