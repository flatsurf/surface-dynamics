#encoding: utf-8
r"""
Cover of permutations
"""

from sage.structure.sage_object import SageObject

from sage.misc.cachefunc import cached_method
from sage.interfaces.gap import gap
from sage.rings.integer import Integer
from sage.matrix.constructor import Matrix, identity_matrix

from copy import copy
from labelled import mean_and_std_dev

class PermutationCover(SageObject):
    r"""
    Finite cover of a permutation

    Has three attributes

    - ``_base`` -- base permutation

    - ``_degree_cover`` -- degree of the cover

    - ``_permut_cover`` -- permutations of the copies over each interval
    """
    _base = None
    _degree_cover = None
    _permut_cover = None

    def __init__(self, base, degree, permut):
        self._base = base.__copy__()
        self._degree_cover = degree
        self._permut_cover = copy(permut)

        if not self._base.is_irreducible():
            raise ValueError("the base must be irreducible")

    def __len__(self):
        return(len(self._base))

    def __repr__(self):
        r"""
        A representation of the generalized permutation cover.

        INPUT:

        - ``sep`` - (default: '\n') a separator for the two intervals

        OUTPUT:

        string -- the string that represents the permutation


        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p1 = iet.Permutation('a b c', 'c b a')
            sage: p1.cover(['(1,2)', '(1,3)', '(2,3)'])
            Covering of degree 3 of the permutation:
            a b c
            c b a
        """
        s = 'Covering of degree %i of the permutation:\n'%(self._degree_cover)
        s += str(self._base)
        return s

    def __getitem__(self,i):
        return(self._base[i])

    def __copy__(self):
        q = PermutationCover(self._base.__copy__(), self._degree_cover, copy(self._permut_cover))
        return(q)

    def covering_data(self, label, list=True):
        r"""
        Returns the permutation associated to a label.

        INPUT:

        - ``label`` - label you are looking at

        - ``list`` - boolen (default: True) - wether or not retourning as a list or sage Permutation

        EXAMPLES::

            sage: from surface_dynamics.all import *

            sage: p = QuadraticStratum([1,1,-1,-1]).components()[0].permutation_representative()
            sage: pc = p.orientation_cover()
            sage: pc
            Covering of degree 2 of the permutation:
            0 1 2 3 3
            2 1 4 4 0

            sage: pc.covering_data(1)
            [1, 2]
            sage: pc.covering_data(3, list=False)
            [2, 1]

        """
        alphabet = self._base.alphabet()
        perm = self._permut_cover[alphabet.rank(label)]
        perm = [i+1 for i in perm]
        if list:
            return(perm)
        else:
            from sage.combinat.permutation import Permutation
            return Permutation(perm)

    def interval_diagram(self, glue_ends=True, sign=False):
        base_diagram = self._base.interval_diagram(glue_ends=False, sign = True)
        singularities= []

        alphabet = self._base.alphabet()
        perm = lambda label : self._permut_cover[alphabet.rank(label)]

        for orbit in base_diagram:
            init_label = orbit[0][0]
            cover_copies = set(range(self._degree_cover))
            while cover_copies:
                d_init = cover_copies.pop()
                d = d_init
                singularity = []
                while True:
                    for base_singularity in orbit:
                        label, j = base_singularity
                        singularity.append((label, j, d) if sign else (label, d))
                        d = perm(label).index(d) if j else perm(label)[d]
                        # let us remind the convention that covering data is permutation
                        # for path going through label with a 0 number in the canonical
                        # orientation
                    if d == d_init: break
                    else: cover_copies.remove(d) 

                singularities.append(singularity)

        return singularities

    def homology_matrix(self):
        r"""
        Return proejction matrix from the left from relative to absolute homology space.
        """
        singularities = self.interval_diagram(glue_ends=False, sign=True)
        nb_sing = len(singularities)
        borders = []
        for d in range(self._degree_cover):
            for a in self._base.alphabet():
                border_side = [0 for _ in xrange(nb_sing)]
                i = 0
                while i < nb_sing and not any((a,0,d) == sing for sing in singularities[i]):
                    i += 1
                border_side[i] = 1
                j = 0
                while j < nb_sing and not any((a,1,d) == sing for sing in singularities[j]):
                    j += 1
                border_side[j] = -1
                if i == j and i < nb_sing: border_side[j] = 0
                borders.append(border_side)
        return(Matrix(borders).left_kernel().matrix())
                

    def profile(self):
        r"""
        Return the profile of the surface.

        Return the list of angles of singularities in the surface divided by pi.

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p = iet.Permutation('a b', 'b a')
            sage: p.cover(['(1,2)', '(1,3)']).profile()
            [6]
            sage: p.cover(['(1,2,3)','(1,4)']).profile()
            [6, 2]
            sage: p.cover(['(1,2,3)(4,5,6)','(1,4,7)(2,5)(3,6)']).profile()
            [6, 2, 2, 2, 2]

            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: p.cover(['(1,2)', '()', '(1,2)']).profile()
            [2, 2, 2, 2]
        """
        from sage.combinat.partition import Partition
        from sage.combinat.permutation import Permutations

        s = []

        base_diagram = self._base.interval_diagram(sign=True,glue_ends=True)
        p_id = Permutations(self._degree_cover).identity()
        for orbit in base_diagram:
            flat_orbit = []
            for x in orbit:
                if isinstance(x[0], tuple):
                    flat_orbit.extend(x)
                else:
                    flat_orbit.append(x)
            p = p_id
            for lab,sign in flat_orbit:
                q = self.covering_data(lab, list=False)
                if sign: q = q.inverse()
                p = p*q
            for c in p.cycle_type():
               s.append(len(orbit)*c)

        return Partition(sorted(s,reverse=True))

    def is_orientable(self):
        r"""
        Test whether the covering has an orientable foliation.

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: from itertools import product
            sage: it = iter(product(('()', '(1,2)'), repeat=3))
            sage: it.next()
            ('()', '()', '()')

            sage: for cov in it:
            ....:     c = p.cover(cov)
            ....:     print cov, c.is_orientable()
            ('()', '()', '(1,2)') False
            ('()', '(1,2)', '()') False
            ('()', '(1,2)', '(1,2)') False
            ('(1,2)', '()', '()') False
            ('(1,2)', '()', '(1,2)') True
            ('(1,2)', '(1,2)', '()') False
            ('(1,2)', '(1,2)', '(1,2)') False
        """
        from surface_dynamics.interval_exchanges.template import PermutationIET
        if isinstance(self._base, PermutationIET):
            return True
        elif any(x%2 for x in self.profile()):
            return False
        else:
            # here we should really unfold the holonomy
            # i.e. try to orient each copy with pm 1
            signs = [+1] + [None] * (self._degree_cover - 1)
            todo = [0]
            A = self._base.alphabet()
            p0 = set(map(A.rank, self._base[0]))
            p1 = set(map(A.rank, self._base[1]))
            inv_letters = p0.symmetric_difference(p1)
            while todo:
                i = todo.pop()
                s = signs[i]
                assert s is not None
                for j in inv_letters:
                    ii = self._permut_cover[j][i]
                    if signs[ii] is None:
                        signs[ii] = -s
                        todo.append(ii)
                    elif signs[ii] != -s:
                        return False
            return True

    def stratum(self):
        r"""
        Stratum of the covering translation surface
        
        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: p.cover(['(1,2)', '()', '(1,2)']).stratum()
            H_1(0^4)
            sage: p.cover(['(1,2)', '(1,2)', '(1,2)']).stratum()
            Q_0(0^2, -1^4)
        """
        if self.is_orientable():
            from surface_dynamics.flat_surfaces.abelian_strata import AbelianStratum
            return AbelianStratum([(x-2)/2 for x in self.profile()])
        else:
            from surface_dynamics.flat_surfaces.quadratic_strata import QuadraticStratum
            return QuadraticStratum([x-2 for x in self.profile()])

    def genus(self):
        r"""
        Genus of the covering translation surface

        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: p = iet.GeneralizedPermutation('a a b', 'b c c')
            sage: p.cover(['(1,2)', '()', '(1,2)']).genus()
            1
            sage: p.cover(['(1,2)', '(1,2)', '(1,2)']).genus()
            0

        TESTS::

            sage: from surface_dynamics.all import *
            sage: o = AbelianStratum([1,2,3,4]).one_component().one_origami()
            sage: assert(o.genus() == AbelianStratum([1,2,3,4]).genus())
            sage: qc = QuadraticStratum([1,2,3,4,-1,-1]).one_component()
            sage: p = qc.permutation_representative()
            sage: assert(p.orientation_cover().genus() == qc.orientation_cover_component().genus())
        """
        p = self.profile()
        return Integer((sum(p)-2*len(p))/4+1)


    def base_stratum(self):
        return self._base.stratum()

    def base_genus(self):
        return self._base.genus()

    def lyapunov_exponents_H_plus(self, nb_vectors=None, nb_experiments=100,
                                  nb_iterations=32768, lengths=None, output_file=None, 
                                  return_speed=False, isotypic_decomposition=False, verbose=False):
        r"""
        Compute the H^+ Lyapunov exponents in  the covering locus.

        It calls the C-library lyap_exp interfaced with Cython. The computation
        might be significantly faster if ``nb_vectors=1`` (or if it is not
        provided but genus is 1).

        INPUT:
 
        - ``nb_vectors`` -- the number of exponents to compute. The number of
          vectors must not exceed the dimension of the space!

         - ``nb_experiments`` -- the number of experiments to perform. It might
           be around 100 (default value) in order that the estimation of
           confidence interval is accurate enough.

         - ``nb_iterations`` -- the number of iteration of the Rauzy-Zorich
           algorithm to perform for each experiments. The default is 2^15=32768
           which is rather small but provide a good compromise between speed and
           quality of approximation.

        - ``lengths`` -- specify some lengths for the iet you want to start with.
        None by default will imply random lengths.

        - ``output_file`` -- if provided (as a file object or a string) output
          the additional information in the given file rather than on the
          standard output.

        - ``return_speed`` -- wether or not return the lyapunov exponents list
        in a pair with the speed of the geodesic.

        - ``isotypic_decomposition`` -- whether or not decompose the space
        into isotypic subspaces for the cover automorphisms.

        - ``verbose`` -- if ``True`` provide additional informations rather than
          returning only the Lyapunov exponents (i.e. ellapsed time, confidence
          intervals, ...)


        EXAMPLES::

            sage: from surface_dynamics.all import *
            sage: QuadraticStratum([1,1,-1,-1]).one_component().lyapunov_exponents_H_plus() #abs tol .1
            [2./3]

            sage: from surface_dynamics.all import *
            sage: q = QuadraticStratum([1,1,-1,-1]).one_component()
            sage: p = q.permutation_representative(reduced=False).orientation_cover()
            sage: p.lyapunov_exponents_H_plus(isotypic_decomposition=True) # abs tol .1
            [[1., 1./3], [2./3]]
        """
        import surface_dynamics.interval_exchanges.lyapunov_exponents as lyapunov_exponents
        import time

        n = len(self)

        if nb_vectors is None:
            nb_vectors = self.genus()

        if output_file is None:
            from sys import stdout
            output_file = stdout
        elif isinstance(output_file, str):
            output_file = open(output_file, "w")

        sigma = sum(self._permut_cover, [])
        sigma = map(int, sigma)
        nb_vectors = int(nb_vectors)
        nb_experiments = int(nb_experiments)
        nb_iterations = int(nb_iterations)

        if verbose:
            output_file.write("Stratum : " + str(self.cover_stratum()))
            output_file.write("\n")

        if nb_vectors < 0 :     raise ValueError("the number of vectors must be positive")
        if nb_vectors == 0:     return []
        if nb_experiments <= 0: raise ValueError("the number of experiments must be positive")
        if nb_iterations <= 0 : raise ValueError("the number of iterations must be positive")
        if len(sigma) %n != 0 : raise ValueError("you must give a permutation for each interval")
        
        #Translate our structure to the C structure"
        k = int(len(self[0]))
        def convert((i,j)):
            return(j + i*k)

        gp, twin = range(2*n), range(2*n)

        base_twin = self._base._twin
        alphabet = self._base.alphabet()

        for i in range(2):
            for j in range(len(self[i])):
                gp[convert((i,j))] = int(alphabet.rank(self[i][j]))
                if isinstance(base_twin[i][j], int):
                       twin[convert((i,j))] = int(convert((1-i, base_twin[i][j])))
                else:
                       twin[convert((i,j))] = int(convert(base_twin[i][j]))

        if lengths != None:
            lengths = map(int, lengths)

        projections = None
        dimensions=[nb_vectors]
        if isotypic_decomposition:
            from sage.rings.all import CC
            from sage.functions.other import real

            size_of_matrix = len(self.cover_generators())
            projections = range(size_of_matrix**2 * self.n_characters())
            dimensions = range(self.n_characters())
            for i_char in xrange(self.n_characters()):
                M = self.isotypic_projection_matrix(i_char)
                H = self.homology_matrix()
                dimensions[i_char] = (H*M).rank()//2
                for i in xrange(size_of_matrix):
                    for j in xrange(size_of_matrix):
                        projections [i_char * (size_of_matrix**2) + i * size_of_matrix + j] = float(real(CC(M[j][i])))

        t0 = time.time()
        res = lyapunov_exponents.lyapunov_exponents_H_plus_cover(
            gp, k, twin, sigma, nb_experiments, nb_iterations, 
            dimensions, projections, None, verbose) #No lengths for random length
        t1 = time.time()

        res_final = []

        m,d = mean_and_std_dev(res[0])
        lexp = m

        if verbose:
            from math import log, floor, sqrt
            output_file.write("sample of %d experiments\n"%nb_experiments)
            output_file.write("%d iterations (~2^%d)\n"%(
                    nb_iterations,
                    floor(log(nb_iterations) / log(2))))
            output_file.write("ellapsed time %s\n"%time.strftime("%H:%M:%S",time.gmtime(t1-t0)))
            output_file.write("Lexp Rauzy-Zorich: %f (std. dev. = %f, conf. rad. 0.01 = %f)\n"%(
                    m,d, 2.576*d/sqrt(nb_experiments)))

        if isotypic_decomposition:
            i_0 = 1
            for i_char in xrange(self.n_characters()):
                res_int = []
                if verbose:
                    output_file.write("##### char_%d #####\n"%(i_char))
                for i in xrange(i_0, i_0 + dimensions[i_char]) :
                    m,d = mean_and_std_dev(res[i])
                    res_int.append(m)
                if verbose:
                    output_file.write("theta%d           : %f (std. dev. = %f, conf. rad. 0.01 = %f)\n"%(
                        i,m,d, 2.576*d/sqrt(nb_experiments)))
                res_final.append(res_int)
                i_0 += dimensions[i_char]

        else:
            for i in xrange(1,nb_vectors+1):
                m,d = mean_and_std_dev(res[i])
                if verbose:
                    output_file.write("theta%d           : %f (std. dev. = %f, conf. rad. 0.01 = %f)\n"%(i,m,d, 2.576*d/sqrt(nb_experiments)))
                res_final.append(m)

        if return_speed: return (lexp, res_final)
        else: return res_final


    @cached_method
    def automorphism_group(self):
        r"""
        Return the galois group of the cover 
        """
        perm = map(lambda p : [i+1 for i in p], self._permut_cover)
        perm = map(lambda p : "PermList( %s )"%(str(p)), perm)

        G = gap("Centralizer(SymmetricGroup(%s), "%self._degree_cover 
                + "Group([" + ', '.join(perm) + "]))") 
        return G

    @cached_method
    def _characters(self):
        r"""
        We want to decompose homology space with the observation that it gives a representation of the galois group.
        
        .. warning::

        Internal class! Do not use directly!
        
        OUPUT:
        
        - tuple -- composed of 5 elements which will be gien by 5 cached functions:

            - character_table: table of character, character[i][g] give the value of 
            the i-th character on g
            g is given by a number
            - character_degree: list of degree of every character
            - autmorphism_group_order : Order of the group
            - automorphism_group_permutation : table s.t. perm[g] give the permutation 
            associated to the group element g on the cover
            - n_characters : number of characters


        """
        from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField
        UCF = UniversalCyclotomicField()
        G = self.automorphism_group()
        G_order, T = gap.Order(G)._sage_(), gap.CharacterTable(G)
        irr_characters = gap.Irr(T)
        n_characters = len(irr_characters)
        character_table = [[UCF(irr_characters[i][j]._sage_()) for j in xrange(1, G_order + 1)] for i in xrange(1, n_characters + 1)]
        character_degree = [gap.Degree(irr_characters[i])._sage_() for i in xrange(1, n_characters + 1)]
        gap_size_centralizers = gap.SizesCentralizers(T)
        gap_orders = gap.OrdersClassRepresentatives(T)
        elements_group = gap.Elements(G)
        perm = [list(gap.ListPerm(elements_group[i], self._degree_cover)) 
                for i in range(1, G_order + 1)]

        #extract real characters
        real_character_table = []
        i = 0
        while i < n_characters:
            chi = character_table[i]
            if any(not x.is_real() for x in chi):
                assert chi == map(lambda z: z.conjugate(), character_table[i+1])
                real_character_table.append([chi[j] + chi[j].conjugate() for j in xrange(n_characters)])
                i += 1
            else:
                real_character_table.append(chi)
            i += 1

        n_real_characters = len(real_character_table)

        return(real_character_table, character_degree, G_order, map(lambda x: [i-1 for i in x], perm), n_real_characters)

    def character_table(self):
        r"""
        Table of characters. 

        character_table()[i][g] give the value of the i-th character on g.
        the g-th element of the group is a permutation given by the list
        automorphism_group_permutation
        """
        return self._characters()[0]

    def character_degree(self):
        r"""
        List of degrees for every character.
        """
        return self._characters()[1]

    def automorphism_group_order(self):
        r"""
        Order of the automorphism group.
        """
        return Integer(self._characters()[2])

    def automorphism_group_permutation(self):
        return self._characters()[3]

    def n_characters(self):
        return Integer(self._characters()[4])

    def cover_generators(self):
        alphabet = self._base.alphabet()
        return([(d,a) for d in range(1, self._degree_cover+1) for a in alphabet])

    def isotypic_projection_matrix(self, i_character):
        r"""
        Return the projection matrix to the isotypic space of the \chi_i
        character recall the formula of Theorem 8 in [Ser]:

        .. MATH::

            p_i_character (d, a) = \sum_{t \in G} char_i(t).conjugate (automorphism_group_permutation(t)(d), a)

        REFERENCES:
        
        .. [Ser] J.-P. Serre, "Représentation des groupes finis."
        
        """
        res = [[0 for _ in self.cover_generators()] for _ in self.cover_generators()]
        char = self.character_table()[i_character]
        coeff = self.character_degree()[i_character]/self.automorphism_group_order()
        alphabet = self._base.alphabet()
        for d, a in self.cover_generators():
            for t in range(self.automorphism_group_order()):
                i = (d - 1)*len(alphabet) + alphabet.rank(a)
                j = int(self.automorphism_group_permutation()[t][d-1])*len(alphabet) + alphabet.rank(a)
                res[i][j] += coeff*char[t]
        return Matrix(res)
