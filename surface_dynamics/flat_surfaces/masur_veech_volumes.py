r"""
Masur-Veech volumes of Abelian strata

This module implements the Moeller-Sauvaget-Zagier formula that gives
a recursive computation of the Masur-Veech volumes of strata of Abelian
differentials.

TODO:

- Moeller-Sauvaget-Zagier for Abelian differentials
- Chen-Moeller-Zagier for quadratic
"""

from sage.all import QQ, zeta
from .abelian_strata import AbelianStratum, AbelianStratumComponent
from .quadratic_strata import QuadraticStratum, QuadraticStratumComponent

# In the table below, the volume is normalized by dividing by zeta(2g)
abelian_volumes_table = {
    AbelianStratum(2).hyperelliptic_component(): QQ((3,4)),
    AbelianStratum(1,1).hyperelliptic_component(): QQ((2,3)),
    AbelianStratum(4).hyperelliptic_component(): QQ((9,64)),
    AbelianStratum(4).odd_component(): QQ((7,18)),
    AbelianStratum(3,1).unique_component(): QQ((16,45)),
    AbelianStratum(2,2).hyperelliptic_component(): QQ((1,10)),
    AbelianStratum(2,2).odd_component(): QQ((7,32)),
    AbelianStratum(2,1,1).unique_component(): QQ((1,4)),
    AbelianStratum(1,1,1,1).unique_component(): QQ((7,36)),
    AbelianStratum(6).hyperelliptic_component(): QQ((25, 1536)),
    AbelianStratum(6).odd_component(): QQ((1,4)),
    AbelianStratum(6).even_component(): QQ((64,405)),
    AbelianStratum(5,1).unique_component(): QQ((9,35)),
    AbelianStratum(4,2).odd_component(): QQ((5,42)),
    AbelianStratum(4,2).even_component(): QQ((45,512)),
    AbelianStratum(3,3).non_hyperelliptic_component(): QQ((5,27)),
    AbelianStratum(3,3).hyperelliptic_component(): QQ((1,105)),
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

def masur_veech_volume(C, rational=False):
    r"""
    EXAMPLES::

        sage: from surface_dynamics import AbelianStratum
        sage: from surface_dynamics.flat_surfaces.masur_veech_volumes import masur_veech_volume
        sage: masur_veech_volume(AbelianStratum(2))
        1/120*pi^4
    """
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
        raise ValueError

    return vol if rational else vol * zeta(2 * S.genus())
