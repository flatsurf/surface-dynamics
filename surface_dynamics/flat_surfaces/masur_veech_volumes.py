r"""
Masur-Veech volumes of Abelian strata and their connected components

.. TODO::

    Implement more formulas (Moeller-Sauvaget-Zagier for Abelian differentials,
    Chen-Moeller-Zagier for quadratic principal stratum, Elise tables, Eskin-Okounkov,
    Eskin-Okounkov-Pandharipande, etc)
"""

from sage.all import ZZ, QQ, zeta, pi
from sage.arith.misc import bernoulli, factorial

from .abelian_strata import AbelianStratum, AbelianStratumComponent
from .quadratic_strata import QuadraticStratum, QuadraticStratumComponent

# In the table below, the volume is normalized by dividing by zeta(2g)
# These values appear in
# - Eskin-Masur-Zorich "principal boundary ..."
abelian_volumes_table = {
    # dim 2
    AbelianStratum(0).hyperelliptic_component(): QQ((2,1)),
    # dim 4
    AbelianStratum(2).hyperelliptic_component(): QQ((3,4)),
    # dim 5
    AbelianStratum(1,1).hyperelliptic_component(): QQ((2,3)),
    # dim 6
    AbelianStratum(4).hyperelliptic_component(): QQ((9,64)),
    AbelianStratum(4).odd_component(): QQ((7,18)),
    # dim 7
    AbelianStratum(3,1).unique_component(): QQ((16,45)),
    AbelianStratum(2,2).hyperelliptic_component(): QQ((1,10)),
    AbelianStratum(2,2).odd_component(): QQ((7,32)),
    # dim 8
    AbelianStratum(6).hyperelliptic_component(): QQ((25, 1536)),
    AbelianStratum(6).odd_component(): QQ((1,4)),
    AbelianStratum(6).even_component(): QQ((64,405)),
    AbelianStratum(2,1,1).unique_component(): QQ((1,4)),
    # dim 9
    AbelianStratum(5,1).unique_component(): QQ((9,35)),
    AbelianStratum(4,2).odd_component(): QQ((5,42)),
    AbelianStratum(4,2).even_component(): QQ((45,512)),
    AbelianStratum(3,3).non_hyperelliptic_component(): QQ((5,27)),
    AbelianStratum(3,3).hyperelliptic_component(): QQ((1,105)),
    AbelianStratum(1,1,1,1).unique_component(): QQ((7,36)),
    # dim 10
    # AbelianStratum(8).hyperelliptic_component()
    # AbelianStratum(8).odd_component()
    # AbelianStratum(8).even_component()
    AbelianStratum(4,1,1).unique_component(): QQ((275,1728)),
    AbelianStratum(3,2,1).unique_component(): QQ((2,15)),
    AbelianStratum(2,2,2).odd_component(): QQ((155,2304)),
    AbelianStratum(2,2,2).even_component(): QQ((37,720)),
    # dim 11
    # AbelianStratum(7, 1)^c
    # AbelianStratum(6, 2)^odd
    # AbelianStratum(6, 2)^even
    # AbelianStratum(5, 3)^c
    # AbelianStratum(4^2)^hyp
    # AbelianStratum(4^2)^odd
    # AbelianStratum(4^2)^even
    AbelianStratum(3,1,1,1).unique_component(): QQ((124,1215)),
    AbelianStratum(2,2,1,1).unique_component(): QQ((131,1440))
}

def masur_veech_volume(C, rational, method):
    r"""
    Return the Masur-Veech volume of the stratum or component of stratum ``C``.

    INPUT:

    - ``rational`` (boolean) - if ``False`` (default) return the Masur-Veech volume
      and if ``True`` return the Masur-Veech volume divided by `\zeta(2g)`.

    - ``method`` - the method to use to compute the volume either

      - ``"table"`` - for a table lookup (all strata up to dimension 9 and some
        strata up to dimension 11)

      - ``"CMSZ"`` - the Chen-Möller-Sauvaget-Zagier recursion (currently only
        implemented for the principal stratum)

    TESTS::

        sage: from surface_dynamics import AbelianStratum
        sage: from surface_dynamics.flat_surfaces.masur_veech_volumes import masur_veech_volume

        sage: H4 = AbelianStratum(4)
        sage: masur_veech_volume(H4, False, 'table')
        61/108864*pi^6
        sage: masur_veech_volume(H4, False, 'CMSZ')
        61/108864*pi^6
        sage: masur_veech_volume(H4.hyperelliptic_component(), False, 'table')
        1/6720*pi^6
        sage: masur_veech_volume(H4.hyperelliptic_component(), False, 'CMSZ')
        1/6720*pi^6
        sage: masur_veech_volume(H4.odd_component(), False, 'table')
        1/2430*pi^6
        sage: masur_veech_volume(H4.odd_component(), False, 'CMSZ')
        1/2430*pi^6

        sage: H6 = AbelianStratum(6)
        sage: all(masur_veech_volume(C, True, 'table') == masur_veech_volume(C, True, 'CMSZ') for C in H6.components())
        True
    """
    if method is None:
        if (isinstance(C, AbelianStratum) and len(C.zeros()) == 1) or \
           (isinstance(C, AbelianStratumComponent) and len(C.stratum().zeros()) == 1):
            method = 'CMSZ'
        else:
            method = 'table'

    if method == 'table':
        if isinstance(C, AbelianStratumComponent):
            vol = abelian_volumes_table[C]
            S = C.stratum()
        elif isinstance(C, AbelianStratum):
            vol = sum(abelian_volumes_table[CC] for CC in C.components())
            S = C
        elif isinstance(C, QuadraticStratumComponent):
            raise NotImplementedError
        elif isinstance(C, QuadraticStratum):
            raise NotImplementedError
        else:
            raise ValueError('invalid input')

        return vol if rational else vol * zeta(2 * S.genus())

    elif method == 'CMSZ':
        if isinstance(C, AbelianStratum):
            if len(C.zeros()) != 1:
                raise NotImplementedError
            g = C.genus()
            # be careful, the output starts in genus g=1
            return minimal_strata_CMSZ(g+1, rational=rational)[g-1]
        elif isinstance(C, AbelianStratumComponent):
            S = C.stratum()
            if len(S.zeros()) != 1:
                raise NotImplementedError
            g = S.genus()
            if C._name == 'hyp':
                return minimal_strata_hyp(g, rational)

            stratum_volume = minimal_strata_CMSZ(g + 1, rational)[g-1]
            spin_diff = minimal_strata_spin_diff(g + 1,rational)[g-1]
            volume = (stratum_volume + spin_diff) / 2 if C.spin() == 0 else (stratum_volume - spin_diff) / 2
            if C.spin() == S.hyperelliptic_component().spin():
                volume -= minimal_strata_hyp(g, rational)
            return volume
        else:
            raise NotImplementedError

    else:
        raise ValueError("unknown method {!r}".format(method))


def minimal_strata_CMSZ(gmax, rational=False):
    r"""
    Return the volumes of `\cH(2g-2)` for the genus `g` going from ``1`` up to ``gmax-1``.

    The algorithm is the one from Sauvaget [Sau18]_ involving an implicit equation. As explained
    in [CheMoeSauZag20]_, one could go through Lagrange inversion. Note that they miss
    factor 2 in their theorem 4.1.

    EXAMPLES::

        sage: from surface_dynamics.flat_surfaces.masur_veech_volumes import minimal_strata_CMSZ
        sage: minimal_strata_CMSZ(6, True)
        [2, 3/4, 305/576, 87983/207360, 1019547/2867200]
        sage: minimal_strata_CMSZ(6, False)
        [1/3*pi^2,
         1/120*pi^4,
         61/108864*pi^6,
         12569/279936000*pi^8,
         12587/3311616000*pi^10]

        sage: from surface_dynamics import AbelianStratum
        sage: from surface_dynamics.flat_surfaces.masur_veech_volumes import masur_veech_volume
        sage: for rat in [True, False]:
        ....:     V0, V2, V4, V6 = minimal_strata_CMSZ(5, rational=rat)
        ....:     MV0 = masur_veech_volume(AbelianStratum(0), rat, 'table')
        ....:     assert V0 == MV0, (V0, MV0, rat)
        ....:     MV2 = masur_veech_volume(AbelianStratum(2), rat, 'table')
        ....:     assert V2 == MV2, (V2, MV2, rat)
        ....:     MV4 = masur_veech_volume(AbelianStratum(4), rat, 'table')
        ....:     assert V4 == MV4, (V4, MV4, rat)
        ....:     MV6 = masur_veech_volume(AbelianStratum(6), rat, 'table')
        ....:     assert V6 == MV6, (V6, MV6, rat)
    """
    n = 2 * gmax - 1
    R = QQ['u']
    u = R.gen()
    # B(u) = formula (15) in [CMSZ20]
    B = (2 * (u/2)._sinh_series(n+1).shift(-1)).inverse_series_trunc(n+1)
    Q = u * (sum(factorial(j-1) * B[j] * u**(j) for j in range(1,n)))._exp_series(n+1)
    # A = formula (14) in [CSMZ20]
    tA = Q.revert_series(n+1).shift(-1).inverse_series_trunc(n)
    # normalized values of volumes in [CMSZ20] are
    # v(m_1, ..., m_n) = (2m_1+1) (2m_2+1) ... (2m_n+1) Vol(m_1, m_2, ..., m_n)
    if rational:
        return [-4 * (2*g) / ZZ(2*g-1) / bernoulli(2*g) * tA[2*g] for g in range(1,gmax)]
    else:
        return [2 * (2*pi)**(2*g) * (-1)**g / ZZ(2*g-1) / factorial(2*g - 1) * tA[2*g] for g in range(1,gmax)]


def minimal_strata_hyp(g, rational=False):
    r"""
    Return the volume of the hyperelliptic component H^{hyp}(2g-2).

    The explicit formula appears in section 6.5 of [CheMoeSauZag20]_.

    EXAMPLES::

        sage: from surface_dynamics.flat_surfaces.masur_veech_volumes import minimal_strata_hyp
        sage: minimal_strata_hyp(2)
        1/120*pi^4
        sage: minimal_strata_hyp(4)
        1/580608*pi^8
        sage: minimal_strata_hyp(10)
        1/137733277917118464000*pi^20
        sage: minimal_strata_hyp(10, rational=True)
        668525/10499279483305984
    """
    if rational:
        return (-1)**(g+1) * 4 * factorial(2*g) / ((2*g-1)*2*g*(2*g+1) * 2**(4*g-2) * bernoulli(2*g) * factorial(g-1)**2)
    else:
        return 2*pi**(2*g) / ((2*g-1)*2*g*(2*g+1) * 2**(2*g-2) * factorial(g-1)**2)


def minimal_strata_spin_diff(gmax, rational=False):
    r"""
    Return the differences of volumes between even and odd components in
    H(2g-2) for the genus `g` going from ``1`` up to ``gmax-1``.

    If there are no even/odd components, the corresponding total volume is 0.
    Formulas are from [CheMoeSauZag20]_.

    EXAMPLES::

        sage: from surface_dynamics.flat_surfaces.masur_veech_volumes import minimal_strata_spin_diff
        sage: minimal_strata_spin_diff(5)
        [-1/3*pi^2, -1/120*pi^4, -143/544320*pi^6, -15697/1959552000*pi^8]
        sage: minimal_strata_spin_diff(5, rational=True)
        [-2, -3/4, -143/576, -15697/207360]
    """
    n = 2 * gmax
    R = QQ['u']
    u = R.gen()
    # B(u) = formula (15) in [CMSZ20]
    B = (2 * (u/2)._sinh_series(n).shift(-1)).inverse_series_trunc(n)
    # Pz and a in section 6.3 of [CMSZ20]
    Pz = (sum(-bernoulli(2*j) * u**(2*j) / (2*j) / 2**j for j in range(1, n // 2)))._exp_series(n)
    a = ((u*Pz.inverse_series_trunc(n)).revert_series(n).shift(-1)).inverse_series_trunc(n)
    # theorem 6.11 in [CMSZ20], normalized volume v(2g-2)=(2g-1)*Vol(2g-2), note the missing factor 2
    if rational:
        return [2* (-2) * a[2*g] * 2 * g / (2*g - 1) / bernoulli(2*g) for g in range(1, gmax)]
    else:
        return [2* (-1)**(g) * (2*pi)**(2*g) * a[2 * g] /(2*g - 1) / factorial(2*g - 1) for g in range(1, gmax)]
