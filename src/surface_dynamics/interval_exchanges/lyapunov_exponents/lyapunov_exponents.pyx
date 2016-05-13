r"""
Python bindings for various computation of Lyapunov exponents.
"""

from libc.stdlib cimport malloc,free
from math import isnan, isinf

cdef extern from "lyapunov_exponents.h":
    ctypedef struct quad_cover:
        pass
    ctypedef struct generalized_permutation:
        pass

    # initialisation/allocation/free
    generalized_permutation * new_generalized_permutation(int *perm, int *twin, int k, int n)
    quad_cover * new_quad_cover(generalized_permutation * gp, size_t ** sigma, size_t degree, size_t nb_vectors)
    void set_lengths(quad_cover * qcc, long double	       *lengths)
    void set_random_lengths_quad_cover(quad_cover * qcc)
    void set_random_vectors(quad_cover * qcc)
    void renormalize_length_quad_cover(quad_cover *qcc)

    #int check_generalized_permutation(generalized_permutation *p)
    #int check_quad_cover(quad_cover * qcc)

    void free_generalized_permutation(generalized_permutation ** gp)
    void free_quad_cover(quad_cover ** qcc)

    # print
    #void print_generalized_permutation(generalized_permutation * p)
    void print_quad_cover(quad_cover * qcc)
    void print_vectors(quad_cover * qcc)
    void print_permutation(size_t *perm, size_t degree)

    # algorithms
    #void renormalize_length_quad_cover(quad_cover * qcc)
    #void rauzy_induction_H_plus_quad_cover(quad_cover *qcc)

    int init_GS(size_t dim)
    void free_GS()
    #void orthogonalize_GS(quad_cover * qcc, double * theta)

    void lyapunov_exponents_H_plus(quad_cover *qcc, double *theta, size_t nb_induction)
    void lyapunov_exponents_isotypic(quad_cover *qcc, double *theta, size_t nb_induction, size_t nb_char, size_t *dimensions, double *proj)
    void top_lyapunov_exponents_H_plus(quad_cover *qcc, double *theta, size_t nb_iterations)

def lyapunov_exponents_H_plus_cover(
    gp, k, twin, sigma,
    nb_experiments, nb_iterations,
    dimensions, projections, lengths, verbose):
    r"""
    Compute the Lyapunov exponents of the H^+ part of the KZ-cocycle for covering locii.

    We assume that all the inputs are clean. If not, it may cause some SEGFAULT
    which would interrupt python!

    INPUT:

    - ``gp`` -- a generalized permutation given as a list of integers

    - ``k`` -- the length of the top interval

    - ``twin`` -- the twin data of the gp

    - ``sigma`` -- covering data

    - ``nb_experiments`` -- number of experimets

    - ``nb_iterations`` -- the number of iterations of the Rauzy-Zorich
      induction to perform

    - ``dimensions`` -- number of vectors we want for each component

    - ``projections`` -- isotypic projection matrices

    - ``verbose`` -- if ``True`` print additional information concerning the
      mean and standard deviation
    """
    if not projections and len(dimensions)>1:
        print "Warning, you are computing lyapunov exponents on several familly of vectors whitout any projections"

    n = int(len(gp))/2

    cdef int *p, *t   # permutation, twin, sigma
    cdef size_t **s
    cdef size_t *tab
    cdef size_t nc = len(dimensions) if dimensions else 0
    cdef size_t degree = int(len(sigma))/n if sigma else 1
    cdef size_t i, j, nn
    cdef size_t *dim
    cdef generalized_permutation *gp_c
    cdef quad_cover *qcc
    cdef double * theta
    cdef double *proj


    # convert the data of into C values
    p = <int *> malloc(2*n * sizeof(int))
    t = <int *> malloc(2*n * sizeof(int))
    s = <size_t **> malloc(n * sizeof(size_t*))

    for i in range(2*n):
        p[i] = gp[i]
        t[i] = twin[i]
    
    for i in range(n):
        tab = <size_t *> malloc(degree * sizeof(size_t))
        for j from 0 <= j < degree:
            tab[j] = sigma[j + degree * i]
        s[i] = tab

    if dimensions is not None :
        nb_vectors = sum(dimensions)

    nn = int(nb_vectors) + 1 # We add one to keep track of the length of the geodesic

    theta = <double *> malloc(nn * sizeof(double))

    gp_c = new_generalized_permutation(p, t, k, n)
    qcc = <quad_cover *> new_quad_cover(gp_c, s, degree, nb_vectors)

    if lengths is None:
       set_random_lengths_quad_cover(qcc)
    else:
        l = <long double *> malloc(n * sizeof(long double))	
        for i from 0 <= i < n:
            l[i] = <long double> lengths[i]
        set_lengths(qcc, l)
        free(l)

    free_generalized_permutation(&(gp_c))
    free(p)
    free(t)

    res = [[] for _ in range(nn)]

    nan_or_inf, tot_isnan, tot_isinf = 0, 0, 0

    if nb_vectors == 1:
        for i in xrange(nb_experiments):
            top_lyapunov_exponents_H_plus(qcc, theta, nb_iterations)
            if any(isnan(theta[j]) or isinf(theta[j]) for j in xrange(nn)):
                print [theta[j] for j in range(nn)], "Warning: contains NaN of Inf"
                nan_or_inf  += 1
            else:
                for j in xrange(nn):
                    res[j].append(theta[j]) 
    else:
        init_GS(nb_vectors)

        if projections:
            dim = <size_t *> malloc(nc * sizeof(size_t))
            for i from 0 <= i < nc:
                dim[i] = int(dimensions[i])

            proj = <double *> malloc((n * degree)**2 * nc * sizeof(double))
            for i from 0 <= i < (n * degree)**2 * nc:
                 proj[i] = <double> projections[i]

            for i in range(nb_experiments):
                lyapunov_exponents_isotypic(qcc, theta, nb_iterations, nc, dim, proj)
                #cleaning some experiments which return NaN or inf as a lyapunov exponent
                if any(isnan(theta[j]) or isinf(theta[j]) for j in xrange(nn)):
                    print [theta[j] for j in range(nn)], "Warning: contains NaN of Inf"
                    nan_or_inf  += 1
                else:
                    for j in xrange(nn):
                        res[j].append(theta[j])
            free(proj)
            free(dim)

        else:
            for i in range(nb_experiments):
                lyapunov_exponents_H_plus(qcc, theta, nb_iterations)
                if any(isnan(theta[j]) or isinf(theta[j]) for j in xrange(nn)):
                    print [theta[j] for j in range(nn)], "Warning: contains NaN of Inf"
                    nan_or_inf  += 1
                else:
                    for j in xrange(nn):
                        res[j].append(theta[j])            

    for i in xrange(n):
        free(s[i])
    free(s)

    free_quad_cover(&qcc)
    free_GS()
    free(theta)

    if nan_or_inf > .1 * nb_experiments:
        raise RuntimeError('too many NaN or Inf: got %d among %d experiments'%(nan_or_inf, nb_experiments))

    return res
