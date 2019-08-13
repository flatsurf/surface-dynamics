def cylinder_concatenation(perm_1,perm_2):
    r"""
    Combines two single cylinder permutation representatives of connected components of Abelian strata to produce another single cylinder representative of a different stratum.
    
    INPUT:
    
    - ``perm_1``, ``perm_2`` - two single cylinder permutation representatives
    
    EXAMPLES::
    
        sage: from surface_dynamics import *
        sage: from surface_dynamics.flat_surfaces.single_cylinder_concatenation import cylinder_concatenation
        
        sage: perm_1 = AbelianStratum(4).odd_component().single_cylinder_representative()
        sage: perm_2 = AbelianStratum(6).even_component().single_cylinder_representative()
        sage: perm_3 = cylinder_concatenation(perm_1,perm_2)
        sage: perm_1
        0 1 2 3 4 5
        2 5 4 1 3 0
        sage: perm_2
        0 1 2 3 4 5 6 7
        2 7 6 5 3 1 4 0
        sage: perm_3
        0 1 2 3 4 5 6 7 8 9 10 11 12
        2 5 4 6 3 7 12 11 10 8 1 9 0
        sage: perm_3.stratum_component()
        H_6(6, 4)^even
    """
    from surface_dynamics.interval_exchanges.constructors import GeneralizedPermutation
    
    
    length_1 = len(perm_1[0])-1
    length_2 = len(perm_2[0])-1
    top_row = [i for i in range(length_1+length_2+1)]
    bot_row1 = perm_1[1][:-1]
    bot_row2 = perm_2[1][:-1]
    for i in range(length_1):
        if bot_row1[i]==1:
            bot_row1[i] += length_1
    for j in range(length_2):
        if not bot_row2[j]==1:
            bot_row2[j] += length_1
    bot_row = bot_row1 + bot_row2 + [0]
    return GeneralizedPermutation(top_row,bot_row)