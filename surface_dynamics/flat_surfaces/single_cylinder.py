r"""
Single cylinder permutation representatives.

This file gathers the functions necessary for
:meth:`~AbelianStratum.single_cylinder_representative`,
:meth:`~AbelianStratumComponent.single_cylinder_representative`,
:meth:`~HypAbelianStratumComponent.single_cylinder_representative`,
:meth:`~EvenAbelianStratumComponent.single_cylinder_representative`, and
:meth:`~OddAbelianStratumComponent.single_cylinder_representative` in
:mod:`~surface_dynamics.flat_surfaces.abelian_strata`.

TESTS:

    sage: from surface_dynamics import *
    sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

Check for non-hyperelliptic strata::

    sage: for d in range(2, 12):
    ....:     for A in AbelianStrata(dimension=d, fake_zeros=True):
    ....:         for C in A.components():
    ....:             if C._name == 'hyp': continue
    ....:             o = C.single_cylinder_origami()
    ....:             r = o.r().cycle_tuples(singletons=True)
    ....:             u = o.u().cycle_tuples(singletons=True)
    ....:             assert len(r) == 1 and len(u) == 1
    ....:             assert o.stratum_component(True) == C

.. TODO::

     Add proper checks for hyperelliptic strata.
"""

from surface_dynamics.interval_exchanges.constructors import GeneralizedPermutation
from surface_dynamics.flat_surfaces.abelian_strata import AbelianStratum

def _cylinder_check(perm):
    r"""
    Checks for a single vertical cylinder and a single horizontal cylinder.

    INPUT:

    - ``perm`` - a permutation representative of an Abelian stratum

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: _cylinder_check(iet.Permutation('0 1 2 3 4 5', '3 2 5 4 1 0'))
        False
        sage: _cylinder_check(iet.Permutation('0 1 2 3 4 5', '2 5 4 1 3 0'))
        True
        sage: _cylinder_check(iet.Permutation('a b', 'b a'))
        True
        sage: _cylinder_check(iet.Permutation([0,3,2,1],[1,3,2,0]))
        True
        sage: _cylinder_check(iet.Permutation('1 2 3 4', '4 3 1 2'))
        False
        sage: _cylinder_check(iet.Permutation('A B C D', 'B C D A'))
        False
        sage: _cylinder_check(iet.Permutation('A C D B', 'B C D A', alphabet='ABCD'))
        True
    """
    from sage.combinat.permutation import Permutation
    from surface_dynamics.flat_surfaces.origamis.origami import Origami

    try:
        o = perm.to_origami()
    except ValueError:
        return False
    else:
        return len(o.r().cycle_tuples(singletons=True)) == 1 and len(o.u().cycle_tuples(singletons=True)) == 1

def even_zero_odd(num):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having
    a single vertical cylinder and a single horizontal cylinder in the odd
    component of the Abelian stratum with a single zero of the given order.

    Such representatives were constructed for every stratum of Abelian
    differentials by Jeffreys [Jef19]_.

    INPUT:

    - ``num`` - an even integer at least four.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = even_zero_odd(4)
        sage: perm
        0 1 2 3 4 5
        2 5 4 1 3 0
        sage: perm.stratum_component() == AbelianStratum(4).odd_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = even_zero_odd(6)
        sage: perm
        0 1 2 3 4 5 6 7
        2 5 4 7 3 1 6 0
        sage: perm.stratum_component() == AbelianStratum(6).odd_component()
        True
        sage: _cylinder_check(perm)
        True
    """
    genus = (num+2)//2
    if genus == 3:
        top_row = [0,1,2,3,4,5]
        bot_row = [2,5,4,1,3,0]
        return GeneralizedPermutation(top_row,bot_row)
    else:
        top_row = list(range(2*genus))
        bot_row = [2,5,4,7,3]
        for i in range(9,2*genus+1,2):
            bot_row += [i,i-3]
        bot_row += [1,2*genus-2,0]
        return GeneralizedPermutation(top_row,bot_row)

def no_two_odd(real_zeros):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in the odd component of an
    Abelian stratum having no zeros of order two.

    Such representatives were constructed for every stratum of Abelian
    differentials by Jeffreys [Jef19]_.

    INPUT:

    - ``real_zeros`` - a list of even positive integers none of which
      are equal to two.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = no_two_odd([6,4])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12
        2 5 4 7 3 8 6 9 12 11 1 10 0
        sage: perm.stratum_component() == AbelianStratum(6,4).odd_component()
        True
        sage: _cylinder_check(perm)
        True

    """
    perm = even_zero_odd(real_zeros[0])
    if len(real_zeros) == 1:
        return perm
    else:
        for i in range(1,len(real_zeros)):
            perm = cylinder_concatenation(perm,even_zero_odd(real_zeros[i]))
        return perm

def one_two_odd(real_zeros):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in the odd component of an
    Abelian stratum having one zero of order two.

    Such representatives were constructed for every stratum of Abelian
    differentials by Jeffreys [Jef19]_.

    INPUT:

    - ``real_zeros`` - a list of even positive integers one of which is equal to two.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = one_two_odd([4,2])
        sage: perm
        0 1 2 3 4 5 6 7 8
        2 5 8 3 6 4 1 7 0
        sage: perm.stratum_component() == AbelianStratum(4,2).odd_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = one_two_odd([8,6,2])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
        2 5 4 7 3 8 6 10 12 9 13 11 14 17 16 19 15 1 18 0
        sage: perm.stratum_component() == AbelianStratum(8,6,2).odd_component()
        True
        sage: _cylinder_check(perm)
        True

    """
    real_zeros.remove(2)
    if set(real_zeros) == {4}:
        perm = GeneralizedPermutation([0,1,2,3,4,5,6,7,8],[2,5,8,3,6,4,1,7,0])
        if len(real_zeros) == 1:
            return perm
        else:
            return cylinder_concatenation(perm,no_two_odd(real_zeros[1:]))
    else:
        perm_1 = even_zero_odd(real_zeros[0]-2)
        length_1 = len(perm_1[0])-1
        top_row_1 = perm_1[0]
        bot_row_1 = perm_1[1][:-1]

        for i in range(length_1):
            if bot_row_1[i] == 1:
                bot_row_1[i] += length_1

        top_row_2 = list(range(length_1+1,length_1+6))
        bot_row_2 = [3+length_1,5+length_1,2+length_1,1,4+length_1,0]
        top_row = top_row_1 + top_row_2
        bot_row = bot_row_1 + bot_row_2
        perm = GeneralizedPermutation(top_row,bot_row)

        if len(real_zeros) == 1:
            return perm
        else:
            return cylinder_concatenation(perm,no_two_odd(real_zeros[1:]))

def even_twos_odd(real_zeros,two_count):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in the odd component of an
    Abelian stratum having an even, at least two, number of zeros of order two.

    Such representatives were constructed for every stratum of Abelian
    differentials by Jeffreys [Jef19]_.

    INPUT:

    - ``real_zeros`` - a list of even positive integers an even number of which
      are equal to two.

    - ``two_count`` - a positive integer equal to the number of twos in ``real_zeros``.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = even_twos_odd([2,2],2)
        sage: perm
        0 1 2 3 4 5 6
        2 4 6 3 1 5 0
        sage: perm.stratum_component() == AbelianStratum(2,2).odd_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = even_twos_odd([4,2,2,2,2],4)
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17
        2 4 6 3 7 5 8 10 12 9 13 11 14 17 16 1 15 0
        sage: perm.stratum_component() == AbelianStratum({4:1,2:4}).odd_component()
        True
        sage: _cylinder_check(perm)
        True

    """
    for i in range(two_count):
        real_zeros.remove(2)

    odd_2_2 = GeneralizedPermutation([0,1,2,3,4,5,6],[2,4,6,3,1,5,0])
    twos_perm = odd_2_2

    for i in range((two_count-2)//2):
        twos_perm = cylinder_concatenation(twos_perm,odd_2_2)

    if len(real_zeros) == 0:
        return twos_perm
    else:
        return cylinder_concatenation(twos_perm,no_two_odd(real_zeros))

def odd_twos_odd(real_zeros,two_count):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in the odd component of an
    Abelian stratum having an odd, at least three, number of zeros of order two.

    Such representatives were constructed for every stratum of Abelian
    differentials by Jeffreys [Jef19]_.

    INPUT:

    - ``real_zeros`` - a list of even positive integers an odd number of which
      are equal to two.

    - ``two_count`` - a positive integer equal to the number of twos in ``real_zeros``.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = odd_twos_odd([2,2,2],3)
        sage: perm
        0 1 2 3 4 5 6 7 8 9
        2 8 6 9 4 1 3 5 7 0
        sage: perm.stratum_component() == AbelianStratum(2,2,2).odd_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = odd_twos_odd([4,2,2,2,2,2],5)
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
        2 8 6 9 4 10 3 5 7 11 13 15 12 16 14 17 20 19 1 18 0
        sage: perm.stratum_component() == AbelianStratum({4:1,2:5}).odd_component()
        True
        sage: _cylinder_check(perm)
        True

    """
    for i in range(two_count):
        real_zeros.remove(2)

    odd_2_2 = GeneralizedPermutation([0,1,2,3,4,5,6],[2,4,6,3,1,5,0])

    odd_2_2_2 = GeneralizedPermutation([0,1,2,3,4,5,6,7,8,9],[2,8,6,9,4,1,3,5,7,0])

    twos_perm = odd_2_2_2

    for i in range((two_count-3)//2):
        twos_perm = cylinder_concatenation(twos_perm,odd_2_2)

    if len(real_zeros) == 0:
        return twos_perm
    else:
        return cylinder_concatenation(twos_perm,no_two_odd(real_zeros))

def even_zero_even(num):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having
    a single vertical cylinder and a single horizontal cylinder in the even
    component of the Abelian stratum with a single zero of the given order.

    Such representatives were constructed for every stratum of Abelian
    differentials by Jeffreys [Jef19]_.

    INPUT:

    - ``num`` - an even integer at least six.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = even_zero_even(6)
        sage: perm
        0 1 2 3 4 5 6 7
        2 7 6 5 3 1 4 0
        sage: perm.stratum_component() == AbelianStratum(6).even_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = even_zero_even(8)
        sage: perm
        0 1 2 3 4 5 6 7 8 9
        2 7 6 5 3 9 4 1 8 0
        sage: perm.stratum_component() == AbelianStratum(8).even_component()
        True
        sage: _cylinder_check(perm)
        True

    """
    genus = (num+2)//2
    if genus == 4:
        top_row = [0,1,2,3,4,5,6,7]
        bot_row = [2,7,6,5,3,1,4,0]
        return GeneralizedPermutation(top_row,bot_row)
    else:
        top_row = [i for i in range(2*genus)]
        bot_row = [2,7,6,5,3,9,4]
        for i in range(11,2*genus+1,2):
            bot_row = bot_row + [i,i-3]
        bot_row = bot_row + [1,2*genus-2,0]
        return GeneralizedPermutation(top_row,bot_row)

def no_two_even(real_zeros):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in the even component of
    an Abelian stratum having no zeros of order two.

    Such representatives were constructed for every stratum of Abelian
    differentials by Jeffreys [Jef19]_.

    INPUT:

    - ``real_zeros`` - a list of even positive integers none of which
      are equal to two.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = no_two_even([4,4])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10
        2 10 7 5 8 1 9 6 4 3 0
        sage: perm.stratum_component() == AbelianStratum(4,4).even_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = no_two_even([6,4])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12
        2 7 6 5 3 8 4 9 12 11 1 10 0
        sage: perm.stratum_component() == AbelianStratum(6,4).even_component()
        True

    """
    if set(real_zeros) == {4}:
        four_count = real_zeros.count(4)
        even_4_4 = GeneralizedPermutation([0,1,2,3,4,5,6,7,8,9,10],[2,10,7,5,8,1,9,6,4,3,0])
        real_zeros.remove(4)
        real_zeros.remove(4)

        if real_zeros != []:
            odd_perm = AbelianStratum(real_zeros).odd_component().single_cylinder_representative()
            return cylinder_concatenation(even_4_4,odd_perm)
        else:
            return even_4_4
    else:
        perm = even_zero_even(real_zeros[0])
        if len(real_zeros) == 1:
            return perm
        else:
            odd_perm = AbelianStratum(real_zeros[1:]).odd_component().single_cylinder_representative()
            return cylinder_concatenation(perm,odd_perm)


def one_two_even(real_zeros):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in the even component of
    an Abelian stratum having one zero of order two.

    Such representatives were constructed for every stratum of Abelian
    differentials by Jeffreys [Jef19]_.

    INPUT:

    - ``real_zeros`` - a list of even positive integers one of which is equal to two.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = one_two_even([4,2])
        sage: perm
        0 1 2 3 4 5 6 7 8
        2 4 1 8 7 5 3 6 0
        sage: perm.stratum_component() == AbelianStratum(4,2).even_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = one_two_even([6,2])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10
        2 10 9 8 6 3 5 1 4 7 0
        sage: perm.stratum_component() == AbelianStratum(6,2).even_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = one_two_even([8,6,2])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
        2 7 6 5 3 8 4 10 12 9 13 11 14 17 16 19 15 1 18 0
        sage: perm.stratum_component() == AbelianStratum(8,6,2).even_component()
        True
        sage: _cylinder_check(perm)
        True

    """
    if real_zeros == [6,2]:
        top_row = [0,1,2,3,4,5,6,7,8,9,10]
        bot_row = [2,10,9,8,6,3,5,1,4,7,0]
        return GeneralizedPermutation(top_row,bot_row)
    elif real_zeros == [4,2]:
        top_row = [0,1,2,3,4,5,6,7,8]
        bot_row = [2,4,1,8,7,5,3,6,0]
        return GeneralizedPermutation(top_row,bot_row)
    else:
        real_zeros.remove(2)
        if set(real_zeros) == {4}:
            perm = GeneralizedPermutation([0,1,2,3,4,5,6,7,8],[2,4,1,8,7,5,3,6,0])
            odd_perm = AbelianStratum(real_zeros[1:]).odd_component().single_cylinder_representative()
            return cylinder_concatenation(perm,odd_perm)
        elif set(real_zeros) == {6} or set(real_zeros) == {4, 6}:
            perm = GeneralizedPermutation([0,1,2,3,4,5,6,7,8,9,10],[2,10,9,8,6,3,5,1,4,7,0])
            odd_perm = AbelianStratum(real_zeros[1:]).odd_component().single_cylinder_representative()
            return cylinder_concatenation(perm,odd_perm)
        else:
            perm_1 = even_zero_even(real_zeros[0]-2)
            length_1 = len(perm_1[0])-1
            top_row_1 = perm_1[0]
            bot_row_1 = perm_1[1][:-1]
            for i in range(length_1):
                if bot_row_1[i] == 1:
                    bot_row_1[i] += length_1
            top_row_2 = [i+length_1 for i in range(1,6)]
            bot_row_2 = [3+length_1,5+length_1,2+length_1,1,4+length_1,0]
            top_row = top_row_1 + top_row_2
            bot_row = bot_row_1 + bot_row_2
            perm = GeneralizedPermutation(top_row,bot_row)
            if len(real_zeros) == 1:
                return perm
            else:
                odd_perm = AbelianStratum(real_zeros[1:]).odd_component().single_cylinder_representative()
                return cylinder_concatenation(perm,odd_perm)

def two_twos_even(real_zeros):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in the even component of
    an Abelian stratum having two zeros of order two.

    Such representatives were constructed for every stratum of Abelian
    differentials by Jeffreys [Jef19]_.

    INPUT:

    - ``real_zeros`` - a list of even positive integers two of which are equal to two.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = two_twos_even([4,2,2])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11
        2 8 5 3 1 10 9 6 4 11 7 0
        sage: perm.stratum_component() == AbelianStratum(4,2,2).even_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = two_twos_even([6,2,2])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13
        2 7 6 5 3 8 4 9 11 13 10 1 12 0
        sage: perm.stratum_component() == AbelianStratum(6,2,2).even_component()
        True
        sage: _cylinder_check(perm)
        True

    """
    real_zeros.remove(2)
    real_zeros.remove(2)
    if set(real_zeros) == {4}:
        perm = GeneralizedPermutation([0,1,2,3,4,5,6,7,8,9,10,11],[2,8,5,3,1,10,9,6,4,11,7,0])
        if len(real_zeros) == 1:
            return perm
        else:
            odd_perm = AbelianStratum(real_zeros[1:]).odd_component().single_cylinder_representative()
            return cylinder_concatenation(perm,odd_perm)
    else:
        odd_2_2 = GeneralizedPermutation([0,1,2,3,4,5,6],[2,4,6,3,1,5,0])
        perm = cylinder_concatenation(even_zero_even(real_zeros[0]),odd_2_2)
        if len(real_zeros) == 1:
            return perm
        else:
            odd_perm = AbelianStratum(real_zeros[1:]).odd_component().single_cylinder_representative()
            return cylinder_concatenation(perm,odd_perm)

def even_twos_even(real_zeros,two_count):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in the even component of
    an Abelian stratum having an even, at least four, number of zeros of order two.

    Such representatives were constructed for every stratum of Abelian
    differentials by Jeffreys [Jef19]_.

    INPUT:

    - ``real_zeros`` - a list of even positive integers an even number of which
      are equal to two.

    - ``two_count`` - a positive integer equal to the number of twos in ``real_zeros``.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = even_twos_even([2,2,2,2],4)
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12
        2 5 4 1 12 3 10 7 11 9 6 8 0
        sage: perm.stratum_component() == AbelianStratum(2,2,2,2).even_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = even_twos_even([4,2,2,2,2,2,2],6)
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
        2 5 4 13 12 3 10 7 11 9 6 8 14 16 18 15 19 17 20 23 22 1 21 0
        sage: perm.stratum_component() == AbelianStratum({4:1,2:6}).even_component()
        True
        sage: _cylinder_check(perm)
        True
    """
    for i in range(two_count):
        real_zeros.remove(2)
    odd_2_2 = GeneralizedPermutation([0,1,2,3,4,5,6],[2,4,6,3,1,5,0])
    even_2_2_2_2 = GeneralizedPermutation([0,1,2,3,4,5,6,7,8,9,10,11,12],[2,5,4,1,12,3,10,7,11,9,6,8,0])
    twos_perm = even_2_2_2_2
    for i in range((two_count-4)//2):
        twos_perm = cylinder_concatenation(twos_perm,odd_2_2)
    if len(real_zeros) == 0:
        return twos_perm
    else:
        odd_perm = AbelianStratum(real_zeros).odd_component().single_cylinder_representative()
        return cylinder_concatenation(twos_perm,odd_perm)

def odd_twos_even(real_zeros,two_count):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in the even component of
    an Abelian stratum having an odd, at least three, number of zeros of order two.

    Such representatives were constructed for every stratum of Abelian
    differentials by Jeffreys [Jef19]_.

    INPUT:

    - ``real_zeros`` - a list of even positive integers an even number of which
      are equal to two.

    - ``two_count`` - a positive integer equal to the number of twos in ``real_zeros``.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = odd_twos_even([2,2,2],3)
        sage: perm
        0 1 2 3 4 5 6 7 8 9
        2 9 8 7 6 3 5 1 4 0
        sage: perm.stratum_component() == AbelianStratum(2,2,2).even_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = odd_twos_even([4,2,2,2,2,2],5)
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
        2 9 8 7 6 3 5 10 4 11 13 15 12 16 14 17 20 19 1 18 0
        sage: perm.stratum_component() == AbelianStratum({4:1,2:5}).even_component()
        True
        sage: _cylinder_check(perm)
        True

    """
    for i in range(two_count):
        real_zeros.remove(2)
    odd_2_2 = GeneralizedPermutation([0,1,2,3,4,5,6],[2,4,6,3,1,5,0])
    even_2_2_2 = GeneralizedPermutation([0,1,2,3,4,5,6,7,8,9],[2,9,8,7,6,3,5,1,4,0])
    twos_perm = even_2_2_2
    for i in range((two_count-3)//2):
        twos_perm = cylinder_concatenation(twos_perm,odd_2_2)
    if len(real_zeros) == 0:
        return twos_perm
    else:
        odd_perm = AbelianStratum(real_zeros).odd_component().single_cylinder_representative()
        return cylinder_concatenation(twos_perm,odd_perm)

def odds_right_swap(zero_pair):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in the Abelian stratum
    having a pair of zeros of the given odd orders.

    Performs a column swap on another permutation to achieve this.

    Such a method was described by Jeffreys [Jef19]_.

    INPUT:

    - ``zero_pair`` - a list of two odd positive integers at least three and
      differing by zero or two.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = odds_right_swap([3,3])
        sage: perm
        0 1 2 3 4 5 6 7 8
        2 8 6 5 7 4 1 3 0
        sage: perm.stratum_component() == AbelianStratum(3,3).non_hyperelliptic_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = odds_right_swap([5,5])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12
        2 5 4 7 3 9 6 12 8 11 1 10 0
        sage: perm.stratum_component() == AbelianStratum(5,5).non_hyperelliptic_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = odds_right_swap([5,3])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10
        2 5 4 7 3 10 6 9 1 8 0
        sage: perm.stratum_component() == AbelianStratum(5,3).unique_component()
        True
        sage: _cylinder_check(perm)
        True
    """
    if zero_pair == [3,3]:
        return GeneralizedPermutation([0,1,2,3,4,5,6,7,8],[2,8,6,5,7,4,1,3,0])
    else:
        dif = abs(zero_pair[0]-zero_pair[1])
        if dif == 0:
            j = (min(zero_pair)-3)//2
        else:
            j = (min(zero_pair)-1)//2
        perm_1 = AbelianStratum(4*j+2-dif).odd_component().single_cylinder_representative()
        perm_2 = AbelianStratum(4).odd_component().single_cylinder_representative()
        perm = cylinder_concatenation(perm_1,perm_2)
        top_row = perm[0][1:]
        bot_row = perm[1][:-1]
        top_row = top_row[:-5]+[top_row[-4],top_row[-5]]+top_row[-3:]
        bot_row = bot_row[:-5]+[bot_row[-4],bot_row[-5]]+bot_row[-3:]
        top_row = [0]+top_row
        bot_row = bot_row+[0]
        perm_3 = GeneralizedPermutation(top_row,bot_row)
        perm_3.alphabet(len(perm_3[0]))
        return perm_3

def odds_left_swap(zero_pair):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in the Abelian stratum
    having a pair of zeros of the given odd orders.

    Performs a column swap on another permutation to achieve this.

    Such a method was described by Jeffreys [Jef19]_.

    INPUT:

    - ``zero_pair`` - a list of two odd positive integers at least three and
      differing by at least four.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = odds_left_swap([7,3])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12
        2 5 4 8 3 7 9 6 12 11 1 10 0
        sage: perm.stratum_component() == AbelianStratum(7,3).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = odds_left_swap([11,5])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
        2 5 4 7 3 9 6 12 8 11 13 10 16 15 18 14 1 17 0
        sage: perm.stratum_component() == AbelianStratum(11,5).unique_component()
        True
        sage: _cylinder_check(perm)
        True
    """
    dif = abs(zero_pair[0]-zero_pair[1])
    j = (min(zero_pair)-1)//2
    perm_1 = AbelianStratum(4*j+2).odd_component().single_cylinder_representative()
    perm_2 = AbelianStratum(dif).odd_component().single_cylinder_representative()
    perm = cylinder_concatenation(perm_1,perm_2)
    swap_point = len(perm_2[0])-1
    top_row = perm[0][1:]
    bot_row = perm[1][:-1]
    top_row = top_row[:-(swap_point+1)]+[top_row[-(swap_point)],top_row[-(swap_point+1)]]+top_row[-(swap_point-1):]
    bot_row = bot_row[:-(swap_point+1)]+[bot_row[-(swap_point)],bot_row[-(swap_point+1)]]+bot_row[-(swap_point-1):]
    top_row = [0]+top_row
    bot_row = bot_row+[0]
    perm_3 = GeneralizedPermutation(top_row,bot_row)
    perm_3.alphabet(len(perm_3[0]))
    return perm_3

def no_ones_odds(odd_zeros):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in an Abelian stratum with odd
    order zeros and no zeros of order 1.

    Such representatives were constructed for every stratum of Abelian
    differentials by Jeffreys [Jef19].

    INPUT:

    - ``odd_zeros`` - an even length list of odd positive integers none of which
      are equal to one.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = no_ones_odds([7,3])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12
        2 5 4 8 3 7 9 6 12 11 1 10 0
        sage: perm.stratum_component() == AbelianStratum(7,3).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = no_ones_odds([5,3,3,3])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
        2 5 4 7 3 10 6 9 11 8 12 18 16 15 17 14 1 13 0
        sage: perm.stratum_component() == AbelianStratum(5,3,3,3).unique_component()
        True
        sage: _cylinder_check(perm)
        True
    """
    if len(odd_zeros) == 2 and abs(odd_zeros[0]-odd_zeros[1]) <= 2:
        return odds_right_swap(odd_zeros)
    elif len(odd_zeros) == 2 and abs(odd_zeros[0]-odd_zeros[1]) > 2:
        return odds_left_swap(odd_zeros)
    else:
        perm = no_ones_odds(odd_zeros[:2])
        for i in range(2,len(odd_zeros),2):
            perm = cylinder_concatenation(perm,no_ones_odds(odd_zeros[i:i+2]))
        return perm

def one_one_odds(odd_zeros):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in an Abelian stratum with odd
    order zeros and one zero of order 1.

    Such representatives were constructed for every stratum of Abelian
    differentials by Jeffreys [Jef19].

    INPUT:

    - ``odd_zeros`` - an even length list of odd positive integers one of which
      is equal to one.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = one_one_odds([3,1])
        sage: perm
        0 1 2 3 4 5 6
        2 5 1 6 4 3 0
        sage: perm.stratum_component() == AbelianStratum(3,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = one_one_odds([5,1])
        sage: perm
        0 1 2 3 4 5 6 7 8
        2 4 7 3 1 8 6 5 0
        sage: perm.stratum_component() == AbelianStratum(5,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = one_one_odds([7,3,3,1])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
        2 5 4 9 3 8 6 11 10 7 12 18 16 15 17 14 1 13 0
        sage: perm.stratum_component() == AbelianStratum(7,3,3,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
    """
    num = odd_zeros[0]
    if num == 3:
        perm = GeneralizedPermutation([0,1,2,3,4,5,6],[2,5,1,6,4,3,0])
        odd_zeros.remove(1)
        if len(odd_zeros) == 1:
            return perm
        else:
            return cylinder_concatenation(perm,no_ones_odds(odd_zeros[1:]))
    elif num == 5:
        perm = GeneralizedPermutation([0,1,2,3,4,5,6,7,8],[2,4,7,3,1,8,6,5,0])
        odd_zeros.remove(1)
        if len(odd_zeros) == 1:
            return perm
        else:
            return cylinder_concatenation(perm,no_ones_odds(odd_zeros[1:]))
    else:
        odd_zeros.remove(1)
        perm_1 = AbelianStratum(num-3).odd_component().single_cylinder_representative()
        length_1 = len(perm_1[0])-1
        top_row_1 = perm_1[0]
        bot_row_1 = perm_1[1][:-1]
        for i in range(length_1):
            if bot_row_1[i] == 1:
                bot_row_1[i] = 4+length_1
        top_row_2 = [i+length_1 for i in range(1,6)]
        bot_row_2 = [3+length_1,1+length_1,1,5+length_1,2+length_1,0]
        top_row = top_row_1 + top_row_2
        bot_row = bot_row_1 + bot_row_2
        perm = GeneralizedPermutation(top_row,bot_row)
        if len(odd_zeros) == 1:
            return perm
        else:
            return cylinder_concatenation(perm,no_ones_odds(odd_zeros[1:]))

def two_ones_odds(odd_zeros):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in an Abelian stratum with odd
    order zeros and two zeros of order 1.

    Such representatives were constructed for every stratum of Abelian
    differentials by Jeffreys [Jef19].

    INPUT:

    - ``odd_zeros`` - an even length list of odd positive integers two of which
      are equal to one.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = two_ones_odds([3,3,1,1])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12
        2 5 7 6 4 3 8 11 1 12 10 9 0
        sage: perm.stratum_component() == AbelianStratum(3,3,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = two_ones_odds([5,5,3,3,1,1])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
        2 4 7 3 9 8 6 5 10 12 15 11 17 16 14 13 18 24 22 21 23 20 1 19 0
        sage: perm.stratum_component() == AbelianStratum(5,5,3,3,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
    """
    odd_zeros.remove(1)
    odd_zeros.remove(1)
    perm = cylinder_concatenation(one_one_odds([odd_zeros[0],1]),one_one_odds([odd_zeros[1],1]))
    if len(odd_zeros) == 2:
        return perm
    else:
        return cylinder_concatenation(perm,no_ones_odds(odd_zeros[2:]))

def three_ones_odds(odd_zeros):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in an Abelian stratum with odd
    order zeros and three zeros of order 1.

    Such representatives were constructed for every stratum of Abelian
    differentials by Jeffreys [Jef19].

    INPUT:

    - ``odd_zeros`` - an even length list of odd positive integers three of which
      are equal to one.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = three_ones_odds([3,1,1,1])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10
        2 10 6 5 1 8 4 7 3 9 0
        sage: perm.stratum_component() == AbelianStratum(3,1,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = three_ones_odds([5,1,1,1])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12
        2 12 9 8 1 7 3 6 10 5 4 11 0
        sage: perm.stratum_component() == AbelianStratum(5,1,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = three_ones_odds([7,1,1,1])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
        2 14 10 9 1 4 3 6 5 12 8 11 7 13 0
        sage: perm.stratum_component() == AbelianStratum(7,1,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = three_ones_odds([9,1,1,1])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
        2 16 13 12 1 4 3 6 5 11 7 10 14 9 8 15 0
        sage: perm = three_ones_odds([3,3,3,1,1,1])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
        2 5 7 6 4 3 8 11 13 12 10 9 14 17 1 18 16 15 0
        sage: perm.stratum_component() == AbelianStratum(3,3,3,1,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
    """
    if len(odd_zeros) > 4:
        num = odd_zeros[0]
        odd_zeros.remove(1)
        odd_zeros.remove(num)
        perm = one_one_odds([num,1])
        return cylinder_concatenation(perm,two_ones_odds(odd_zeros))
    else:
        num = odd_zeros[0]
        res_4 = num % 4
        if num == 3:
            return GeneralizedPermutation([0,1,2,3,4,5,6,7,8,9,10],[2,10,6,5,1,8,4,7,3,9,0])
        elif num == 5:
            return GeneralizedPermutation([0,1,2,3,4,5,6,7,8,9,10,11,12],[2,12,9,8,1,7,3,6,10,5,4,11,0])
        elif res_4 == 3:
            top_row = [i for i in range(num+8)]
            bot_row = [2,num+7,num+3,num+2,1]
            for i in range(4,num+1,2):
                bot_row += [i,i-1]
            bot_row += [num+5,num+1,num+4,num,num+6,0]
            return GeneralizedPermutation(top_row,bot_row)
        else:
            top_row = [i for i in range(num+8)]
            bot_row = [2,num+7,num+4,num+3,1]
            for i in range(4,num-1,2):
                bot_row += [i,i-1]
            bot_row += [num+2,num-2,num+1,num+5,num,num-1,num+6,0]
            return GeneralizedPermutation(top_row,bot_row)

def even_ones_odds(odd_zeros,one_count):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in an Abelian stratum with odd
    order zeros and an even, at least four, number of zeros of order 1.

    Such representatives were constructed for every stratum of Abelian
    differentials by Jeffreys [Jef19]_.

    INPUT:

    - ``odd_zeros`` - an even length list of odd positive integers an even number of which
      are equal to one.

    - ``one_count`` - a positive integer equal to the number of ones in ``real_zeros``.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = even_ones_odds([1,1,1,1],4)
        sage: perm
        0 1 2 3 4 5 6 7 8
        2 6 5 3 1 8 4 7 0
        sage: perm.stratum_component() == AbelianStratum(1,1,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = even_ones_odds([1,1,1,1,1,1],6)
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12
        2 8 1 5 11 7 3 10 6 12 9 4 0
        sage: perm.stratum_component() == AbelianStratum(1,1,1,1,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = even_ones_odds([5,3,1,1,1,1],4)
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
        2 6 5 3 9 8 4 7 10 13 12 15 11 18 14 17 1 16 0
        sage: perm.stratum_component() == AbelianStratum(5,3,1,1,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = even_ones_odds([3,3,1,1,1,1,1,1],6)
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
        2 8 13 5 11 7 3 10 6 12 9 4 14 20 18 17 19 16 1 15 0
        sage: perm.stratum_component() == AbelianStratum(3,3,1,1,1,1,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
    """
    for i in range(one_count):
        odd_zeros.remove(1)
    four_ones = GeneralizedPermutation([0,1,2,3,4,5,6,7,8],[2,6,5,3,1,8,4,7,0])
    six_ones = GeneralizedPermutation([0,1,2,3,4,5,6,7,8,9,10,11,12],[2,8,1,5,11,7,3,10,6,12,9,4,0])
    if one_count % 4 == 0:
        perm = four_ones
        for i in range((one_count-4)//4):
            perm = cylinder_concatenation(perm,four_ones)
    else:
        perm = six_ones
        for i in range((one_count-6)//4):
            perm = cylinder_concatenation(perm,four_ones)
    if len(odd_zeros) == 0:
        return perm
    else:
        return cylinder_concatenation(perm,no_ones_odds(odd_zeros))

def odd_ones_odds(odd_zeros,one_count):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in an Abelian stratum with odd
    order zeros and an odd, at least five, number of zeros of order 1.

    Such representatives were constructed for every stratum of Abelian
    differentials by Jeffreys [Jef19]_.

    INPUT:

    - ``odd_zeros`` - an even length list of odd positive integers an odd number of which
      are equal to one.

    - ``one_count`` - a positive integer equal to the number of ones in ``real_zeros``.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = odd_ones_odds([5,1,1,1,1,1],5)
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
        2 4 7 3 9 8 6 5 10 14 13 11 1 16 12 15 0
        sage: perm.stratum_component() == AbelianStratum(5,1,1,1,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = odd_ones_odds([3,1,1,1,1,1,1,1],7)
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
        2 5 7 6 4 3 8 14 1 11 17 13 9 16 12 18 15 10 0
        sage: perm.stratum_component() == AbelianStratum(3,1,1,1,1,1,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
    """
    for i in range(one_count-1):
        odd_zeros.remove(1)
    even_ones = [1 for i in range(one_count-1)]
    return cylinder_concatenation(one_one_odds(odd_zeros),even_ones_odds(even_ones,one_count-1))

def min_on_bot(zero_pair):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in the Abelian stratum
    having a pair of zeros of the given odd orders which differ by two.

    The permutations have a particular form required for the construction
    of other representatives.

    Such representatives were constructed by Jeffreys [Jef19]_.

    INPUT:

    - ``zero_pair`` - a list of two odd positive integers at least one and
      differing by two.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = min_on_bot([3,1])
        sage: perm
        0 1 2 3 4 5 6
        2 6 5 1 4 3 0
        sage: perm.stratum_component() == AbelianStratum(3,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = min_on_bot([5,3])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10
        2 6 4 10 8 3 1 9 7 5 0
        sage: perm.stratum_component() == AbelianStratum(5,3).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = min_on_bot([7,5])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
        2 6 4 10 8 3 12 9 7 5 14 11 1 13 0
        sage: perm.stratum_component() == AbelianStratum(7,5).unique_component()
        True
        sage: _cylinder_check(perm)
        True
    """
    if zero_pair == [3,1]:
        return GeneralizedPermutation([0,1,2,3,4,5,6],[2,6,5,1,4,3,0])
    elif zero_pair == [5,3]:
        return GeneralizedPermutation([0,1,2,3,4,5,6,7,8,9,10],[2,6,4,10,8,3,1,9,7,5,0])
    else:
        num = max(zero_pair)
        top_row = [i for i in range(2*num+1)]
        bot_row = [2,6,4,10,8,3,12,9,7,5]
        for i in range(14,2*num+2,2):
            bot_row = bot_row + [i,i-3]
        bot_row = bot_row + [1,2*num-1,0]
        return GeneralizedPermutation(top_row,bot_row)

def odd_zeros_one_one(odd_zeros):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in an Abelian stratum
    having odd order zeros of the given orders.

    Such representatives were constructed for every stratum of Abelian
    differentials by Jeffreys [Jef19]_.

    INPUT:

    - ``odd_zeros`` - an even length list of odd positive integers.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = odd_zeros_one_one([5,5])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12
        2 5 4 7 3 9 6 12 8 11 1 10 0
        sage: perm.stratum_component() == AbelianStratum(5,5).non_hyperelliptic_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = odd_zeros_one_one([5,1])
        sage: perm
        0 1 2 3 4 5 6 7 8
        2 4 7 3 1 8 6 5 0
        sage: perm.stratum_component() == AbelianStratum(5,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = odd_zeros_one_one([5,3,1,1])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
        2 4 7 3 9 8 6 5 10 13 1 14 12 11 0
        sage: perm.stratum_component() == AbelianStratum(5,3,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = odd_zeros_one_one([5,1,1,1])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12
        2 12 9 8 1 7 3 6 10 5 4 11 0
        sage: perm.stratum_component() == AbelianStratum(5,1,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = odd_zeros_one_one([5,3,1,1,1,1])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
        2 6 5 3 9 8 4 7 10 13 12 15 11 18 14 17 1 16 0
        sage: perm.stratum_component() == AbelianStratum(5,3,1,1,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = odd_zeros_one_one([5,1,1,1,1,1])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
        2 4 7 3 9 8 6 5 10 14 13 11 1 16 12 15 0
        sage: perm.stratum_component() == AbelianStratum(5,1,1,1,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
    """
    one_count = odd_zeros.count(1)
    if one_count == 0:
        return no_ones_odds(odd_zeros)
    elif one_count == 1:
        return one_one_odds(odd_zeros)
    elif one_count == 2:
        return two_ones_odds(odd_zeros)
    elif one_count == 3:
        return three_ones_odds(odd_zeros)
    elif one_count >= 4 and one_count % 2 == 0:
        return even_ones_odds(odd_zeros,one_count)
    else:
        return odd_ones_odds(odd_zeros,one_count)

def only_even_2(odd_zeros):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in the Abelian stratum
    having zeros of the given odd orders and a single zero of order two.

    Such representatives were constructed for every stratum of Abelian
    differentials by Jeffreys [Jef19]_.

    INPUT:

    - ``odd_zeros`` - an even length list of odd positive integers.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = only_even_2([1,1])
        sage: perm
        0 1 2 3 4 5 6 7
        2 6 4 1 7 5 3 0
        sage: perm.stratum_component() == AbelianStratum(2,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = only_even_2([1,1,1,1,1,1])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
        2 6 4 8 7 5 3 9 13 12 10 1 15 11 14 0
        sage: perm.stratum_component() == AbelianStratum({2:1,1:6}).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = only_even_2([1,1,1,1])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11
        2 7 11 6 3 9 5 1 8 4 10 0
        sage: perm.stratum_component() == AbelianStratum(2,1,1,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = only_even_2([1,1,1,1,1,1,1,1])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
        2 7 11 6 3 9 5 12 8 4 10 13 17 16 14 1 19 15 18 0
        sage: perm.stratum_component() == AbelianStratum({2:1,1:8}).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = only_even_2([3,3,1,1])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
        2 6 4 8 7 5 3 9 15 13 12 14 11 1 10 0
        sage: perm.stratum_component() == AbelianStratum(3,3,2,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = only_even_2([3,1,1,1])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13
        2 6 4 8 7 5 3 9 12 1 13 11 10 0
        sage: perm.stratum_component() == AbelianStratum(3,2,1,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = only_even_2([5,3,3,3])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
        2 8 6 5 7 4 9 3 11 13 10 14 12 15 21 19 18 20 17 1 16 0
        sage: perm.stratum_component() == AbelianStratum(5,3,3,3,2).unique_component()
        True
        sage: _cylinder_check(perm)
        True
    """
    if odd_zeros.count(1) == len(odd_zeros) and odd_zeros.count(1) % 4 == 2:
        perm = GeneralizedPermutation([0,1,2,3,4,5,6,7],[2,6,4,1,7,5,3,0])
        if len(odd_zeros) == 2:
            return perm
        else:
            odd_zeros.remove(1)
            odd_zeros.remove(1)
            one_count = odd_zeros.count(1)
            return cylinder_concatenation(perm,even_ones_odds(odd_zeros,one_count))
    elif odd_zeros.count(1) == len(odd_zeros) and odd_zeros.count(1) % 4 == 0:
        perm = GeneralizedPermutation([0,1,2,3,4,5,6,7,8,9,10,11],[2,7,11,6,3,9,5,1,8,4,10,0])
        if len(odd_zeros) == 4:
            return perm
        else:
            odd_zeros.remove(1)
            odd_zeros.remove(1)
            odd_zeros.remove(1)
            odd_zeros.remove(1)
            one_count = odd_zeros.count(1)
            return cylinder_concatenation(perm,even_ones_odds(odd_zeros,one_count))
    elif odd_zeros.count(1) == 2 and len(odd_zeros) == 4 :
        perm = GeneralizedPermutation([0,1,2,3,4,5,6,7],[2,6,4,1,7,5,3,0])
        return cylinder_concatenation(perm,no_ones_odds(odd_zeros[:2]))
    else:
        if len(odd_zeros) == 4 and odd_zeros.count(1) == 3:
            perm = GeneralizedPermutation([0,1,2,3,4,5,6,7],[2,6,4,1,7,5,3,0])
            return cylinder_concatenation(perm,one_one_odds([odd_zeros[0],1]))
        else:
            pair_zeros = odd_zeros[:2]
            odd_zeros = odd_zeros[2:]
            dif = abs(pair_zeros[0]-pair_zeros[1])
            if 1 in pair_zeros:
                if pair_zeros == [3,1]:
                    perm = GeneralizedPermutation([0,1,2,3,4,5,6,7,8,9],[2,6,8,3,7,4,1,9,5,0])
                    if len(odd_zeros) == 0:
                        return perm
                    else:
                        one_count = odd_zeros.count(1)
                        return cylinder_concatenation(perm,even_ones_odds(odd_zeros,one_count))
                else:
                    perm_1 = one_one_odds([pair_zeros[0]-2,1])
                    length_1 = len(perm_1[0])-1
                    top_row_1 = perm_1[0]
                    bot_row_1 = perm_1[1][:-1]
                    for i in range(length_1):
                        if bot_row_1[i] == 1:
                            bot_row_1[i] += length_1
                    top_row_2 = [i+length_1 for i in range(1,6)]
                    bot_row_2 = [3+length_1,5+length_1,2+length_1,1,4+length_1,0]
                    top_row = top_row_1 + top_row_2
                    bot_row = bot_row_1 + bot_row_2
                    perm = GeneralizedPermutation(top_row,bot_row)
                    if len(odd_zeros) == 0:
                        return perm
                    else:
                        one_count = odd_zeros.count(1)
                        return cylinder_concatenation(perm,even_ones_odds(odd_zeros,one_count))
            elif dif > 0:
                pair_zeros[0] += -2
                perm_1 = no_ones_odds(pair_zeros)
                length_1 = len(perm_1[0])-1
                top_row_1 = perm_1[0]
                bot_row_1 = perm_1[1][:-1]
                for i in range(length_1):
                    if bot_row_1[i] == 1:
                        bot_row_1[i] += length_1
                top_row_2 = [i+length_1 for i in range(1,6)]
                bot_row_2 = [3+length_1,5+length_1,2+length_1,1,4+length_1,0]
                top_row = top_row_1 + top_row_2
                bot_row = bot_row_1 + bot_row_2
                perm = GeneralizedPermutation(top_row,bot_row)
                if len(odd_zeros) == 0:
                    return perm
                else:
                    perm_odd = odd_zeros_one_one(odd_zeros)
                    return cylinder_concatenation(perm,perm_odd)
            else:
                pair_zeros[1] += -2
                perm_1 = min_on_bot(pair_zeros)
                length_1 = len(perm_1[0])-1
                top_row_1 = perm_1[0]
                bot_row_1 = perm_1[1][:-1]
                for i in range(length_1):
                    if bot_row_1[i] == 1:
                        bot_row_1[i] += length_1
                top_row_2 = [i+length_1 for i in range(1,6)]
                bot_row_2 = [3+length_1,5+length_1,2+length_1,1,4+length_1,0]
                top_row = top_row_1 + top_row_2
                bot_row = bot_row_1 + bot_row_2
                perm = GeneralizedPermutation(top_row,bot_row)
                if len(odd_zeros) == 0:
                    return perm
                else:
                    odd_perm = odd_zeros_one_one(odd_zeros)
                    return cylinder_concatenation(perm,odd_perm)

def only_odds_11(even_zeros):
    r"""
    Returns a single cylinder permutation representative.

    Returns a permutation representative of a square-tiled surface having a single
    vertical cylinder and a single horizontal cylinder in the Abelian stratum
    having zeros of the given even orders and two zeros of order one.

    Such representatives were constructed for every stratum of Abelian
    differentials by Jeffreys [Jef19]_.

    INPUT:

    - ``even_zeros`` - a list of even positive integers.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

        sage: perm = only_odds_11([2])
        sage: perm
        0 1 2 3 4 5 6 7
        2 6 4 1 7 5 3 0
        sage: perm.stratum_component() == AbelianStratum(2,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = only_odds_11([2,2])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10
        2 4 9 7 3 8 5 1 10 6 0
        sage: perm.stratum_component() == AbelianStratum(2,2,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = only_odds_11([4,2])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12
        2 6 4 8 7 5 3 9 12 11 1 10 0
        sage: perm.stratum_component() == AbelianStratum(4,2,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = only_odds_11([6])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11
        2 6 4 1 8 7 10 9 11 5 3 0
        sage: perm.stratum_component() == AbelianStratum(6,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = only_odds_11([4,4])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
        2 7 4 10 9 5 8 6 3 11 14 13 1 12 0
        sage: perm.stratum_component() == AbelianStratum(4,4,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True
        sage: perm = only_odds_11([8,6])
        sage: perm
        0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
        2 7 4 14 9 8 11 10 13 5 12 6 3 15 18 17 20 16 1 19 0
        sage: perm.stratum_component() == AbelianStratum(8,6,1,1).unique_component()
        True
        sage: _cylinder_check(perm)
        True

    """
    if set(even_zeros) == {2} and len(even_zeros) % 2 == 1:
        perm = GeneralizedPermutation([0,1,2,3,4,5,6,7],[2,6,4,1,7,5,3,0])
        if len(even_zeros) == 1:
            return perm
        else:
            even_zeros.remove(2)
            even_perm = AbelianStratum(even_zeros).odd_component().single_cylinder_representative()
            return cylinder_concatenation(perm,even_perm)
    elif set(even_zeros) == {2} and len(even_zeros) % 2 == 0:
        perm = GeneralizedPermutation([0,1,2,3,4,5,6,7,8,9,10],[2,4,9,7,3,8,5,1,10,6,0])
        if len(even_zeros) == 2:
            return perm
        else:
            even_zeros.remove(2)
            even_zeros.remove(2)
            even_perm = AbelianStratum(even_zeros).odd_component().single_cylinder_representative()
            return cylinder_concatenation(perm,even_perm)
    elif even_zeros.count(2) == 1 and len(even_zeros) == 2:
        perm = GeneralizedPermutation([0,1,2,3,4,5,6,7],[2,6,4,1,7,5,3,0])
        even_perm = AbelianStratum(even_zeros[0]).odd_component().single_cylinder_representative()
        return cylinder_concatenation(perm,even_perm)
    else:
        num = even_zeros[0]
        if num % 4 == 2:
            top_row = [i for i in range(num+6)]
            bot_row = [2,6,4,1]
            for i in range(8,num+6,2):
                bot_row += [i,i-1]
            bot_row += [num+5,5,3,0]
            perm = GeneralizedPermutation(top_row,bot_row)
            if len(even_zeros) == 1:
                return perm
            else:
                even_perm = AbelianStratum(even_zeros[1:]).odd_component().single_cylinder_representative()
                return cylinder_concatenation(perm,even_perm)
        else:
            if num == 4:
                perm = GeneralizedPermutation([0,1,2,3,4,5,6,7,8,9],[2,7,4,1,9,5,8,6,3,0])
            else:
                top_row = [i for i in range(num+6)]
                bot_row = [2,7,4,1]
                for i in range(9,num+5,2):
                    bot_row += [i,i-1]
                bot_row += [num+5,5,num+4,6,3,0]
                perm = GeneralizedPermutation(top_row,bot_row)
            if len(even_zeros) == 1:
                return perm
            else:
                even_perm = AbelianStratum(even_zeros[1:]).odd_component().single_cylinder_representative()
                return cylinder_concatenation(perm,even_perm)

def cylinder_concatenation(perm_1, perm_2, alphabet=None):
    r"""
    Combines two single cylinder permutation representatives.

    Combines two single cylinder permutation representatives of connected components
    of Abelian strata to produce another single cylinder representative of a
    different stratum.

    Such a method was described by Jeffreys [Jef19].

    INPUT:

    - ``perm_1``, ``perm_2`` - two single cylinder permutation representatives.

    - ``alphabet`` - alphabet or ``None`` (defaut: ``None``):
      whether you want to specify an alphabet for your representative.

    EXAMPLES::

        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder import _cylinder_check

    We first take two single cylinder permutation representatives for the odd
    components of H_3(4)^odd and H_4(6)^odd.::

        sage: perm_1 = AbelianStratum(4).odd_component().single_cylinder_representative()
        sage: perm_1
        0 1 2 3 4 5
        2 5 4 1 3 0
        sage: perm_1.stratum_component() == AbelianStratum(4).odd_component()
        True
        sage: _cylinder_check(perm_1)
        True
        sage: perm_2 = AbelianStratum(6).odd_component().single_cylinder_representative()
        sage: perm_2
        0 1 2 3 4 5 6 7
        2 5 4 7 3 1 6 0
        sage: perm_2.stratum_component() == AbelianStratum(6).odd_component()
        True
        sage: _cylinder_check(perm_2)
        True

    We check that the cylinder_concatenation of these permutations produces
    a single cylinder permutation representative of the connected component
    H_6(6,4)^odd.::

        sage: perm_3 = cylinder_concatenation(perm_1,perm_2)
        sage: perm_3
        0 1 2 3 4 5 6 7 8 9 10 11 12
        2 5 4 6 3 7 10 9 12 8 1 11 0
        sage: perm_3.stratum_component() == AbelianStratum(6,4).odd_component()
        True
        sage: _cylinder_check(perm_3)
        True

    We now instead take the cylinder_concatenation of perm_1 with a single cylinder
    permutation representative of the connected component H_4(6)^even. We see
    that the resulting permutation is a single cylinder permutation representative
    of the connected component H_6(6,4)^even.::

        sage: perm_4 = AbelianStratum(6).even_component().single_cylinder_representative()
        sage: perm_5 = cylinder_concatenation(perm_1,perm_4,Alphabet(name='lower'))
        sage: perm_4
        0 1 2 3 4 5 6 7
        2 7 6 5 3 1 4 0
        sage: perm_4.stratum_component() == AbelianStratum(6).even_component()
        True
        sage: _cylinder_check(perm_4)
        True
        sage: perm_5
        a b c d e f g h i j k l m
        c f e g d h m l k i b j a
        sage: perm_5.stratum_component() == AbelianStratum(6,4).even_component()
        True
        sage: _cylinder_check(perm_5)
        True
    """
    from sage.combinat.words.alphabet import Alphabet
    from sage.rings.semirings.non_negative_integer_semiring import NN

    alph = Alphabet(NN)
    perm_1.alphabet(alph)
    perm_2.alphabet(alph)
    length_1 = len(perm_1[0])-1
    length_2 = len(perm_2[0])-1
    top_row = [i for i in range(length_1+length_2+1)]
    bot_row1 = perm_1[1][:-1]
    bot_row2 = perm_2[1][:-1]
    for i in range(length_1):
        if bot_row1[i] == 1:
            bot_row1[i] += length_1
    for j in range(length_2):
        if not bot_row2[j] == 1:
            bot_row2[j] += length_1
    bot_row = bot_row1 + bot_row2 + [0]
    perm = GeneralizedPermutation(top_row,bot_row)
    if not alphabet == None:
        perm.alphabet(alphabet)
    return perm
