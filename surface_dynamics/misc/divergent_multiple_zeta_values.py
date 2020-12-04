from sage.modular.multiple_zeta import Multizetas

M = Multizetas(QQ)

# TODO: without this, it is not possible to create polynomials with Multizeta coefficients
M.is_prime_field = lambda: False

MT = M['T']  # shuffle regularization
MX = M['X']  # stuffle regularization
T = MT.gen()
X = MX.gen()

# concatenation product
def prefix(M, p, g):
    p = M.basis().keys()(p)
    ans = M.zero()
    for u, coeff in g.monomial_coefficients().items():
        ans += coeff * M.term(p + u)
    return ans

def suffix(M, f, s):
    s = M.basis().keys()(s)
    ans = M.zero()
    for u , coeff in f.monomial_coefficients().items():
        ans += coeff * M.term(u + s)
    return ans

def concatenate(M, left, right):
    r"""
    Concatenation product in the free algebra

    EXAMPLES::

        sage: concatenate(M, M((2,)) + M((3,)), M((5,)) + M((7,)))
        ζ(2,5) + ζ(2,7) + ζ(3,5) + ζ(3,7)
    """
    ans = M.zero()
    for u1, c1 in left.monomial_coefficients():
        for u2, c2 in right.monomial_coefficients():
            ans += c1 * c2 * M.term(u1 + u2)
    return ans

def stuffle_on_basis(M, w1, w2):
    r"""
    Return the stuffle product of two words

    EXAMPLES::

        sage: stuffle_on_basis(M, Word([2]), Word([5]))
        ζ(2,5) + ζ(5,2) + ζ(7)
        sage: stuffle_on_basis(M, Word([2,2]), Word([3,4]))
        ζ(2,2,3,4) + ζ(2,3,2,4) + ζ(2,3,4,2) + ζ(2,3,6) + ζ(2,5,4) + ζ(3,2,2,4) + ζ(3,2,4,2) + ζ(3,2,6) + ζ(3,4,2,2) + ζ(3,6,2) + ζ(5,2,4) + ζ(5,4,2) + ζ(5,6)
    """
    W = M.basis().keys()
    if not w1:
        return M.term(W(w2))
    if not w2:
        return M.term(W(w1))
    y1 = w1[0]
    u1 = w1[1:]
    y2 = w2[0]
    u2 = w2[1:]
    if y1 == y2:
        return prefix(M, W([y1]), stuffle_on_basis(M, u1, w2) + stuffle_on_basis(M, w1, u2)) + \
               prefix(M, W([y1+y2]), stuffle_on_basis(M, u1, u2))
    else:
        return prefix(M, W([y1]), stuffle_on_basis(M, u1, w2)) + \
               prefix(M, W([y2]), stuffle_on_basis(M, w1, u2)) + \
               prefix(M, W([y1+y2]), stuffle_on_basis(M, u1, u2))

def stuffle_product(M, left, right):
    ans = M.zero()
    for u1, c1 in left.monomial_coefficients().items():
        for u2, c2 in right.monomial_coefficients().items():
            ans += c1 * c2 * stuffle_on_basis(M, u1, u2)
    return ans

def stuffle_power(M, a, n):
    if n < 0:
        raise ValueError
    if n == 0:
        return M.one()
    apow = a
    while not (n & 1):
        apow = stuffle_product(M, apow, apow)
        n >>= 1

    res = apow
    n >>= 1
    while n:
        apow = stuffle_product(M, apow, apow)
        if n & 1:
            res = stuffle_product(M, apow, res)
        n >>= 1

    return res


def stuffle_regularization(s):
    r"""
    Given a multizeta index s returns it as an expression in MX

    EXAMPLES::

        sage: divergent_stuffle_regularization((1,))
        sage: divergent_stuffle_regularization((2,))
        sage: divergent_stuffle_regularization((2,1))
        sage: divergent_stuffle_regularization((1,1))
        sage: divergent_stuffle_regularization((2,1,1))
        sage: divergent_stuffle_regularization((1,1,1))
    """
    pass

def shuffle_regularization(s):
    r"""
    Given a multizeta index s returns it as an expression in MT
    """
    pass

def shuffle_to_harmonic(n):
    r"""
    Image of X^n under the morphism rho: M[X] -> M[T] of Ihara-Kaneko-Zagier
    """

def harmonic_to_shuffle(n):
    r"""
    Image of T^n under the morphism rho^{-1}: M[T] -> M[X] of Ihara-Kaneko-Zagier
    """
