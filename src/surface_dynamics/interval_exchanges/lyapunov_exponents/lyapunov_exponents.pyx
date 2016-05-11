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
    void lyapunov_exponents_isotopic(quad_cover *qcc, double *theta, size_t nb_induction, size_t nb_char, size_t *dimensions, double *proj)
    void top_lyapunov_exponents_H_plus(quad_cover *qcc, double *theta, size_t nb_iterations)

def lyapunov_exponents_H_plus_cover(
    gp, k, twin, sigma, degree,
    nb_vectors, nb_experiments, nb_iterations, lengths=None, nb_char=0, 
    dimensions=None, projections=None, isotopic_decomposition=False):
    r"""
    Compute the Lyapunov exponents of the H^+ part of the KZ-cocycle for covering locii.

    We assume that all the inputs are clean. If not, it may cause some SEGFAULT
    which would interrupt python!

    INPUT:

    - ``gp`` -- a generalized permutation given as a list of integers

    - ``twin`` -- the twin data of the gp

    - ``k`` -- the length of the top interval

    - ``sigma`` -- covering data

    - ``nb_vectors`` -- the number of vectors to use

    - ``nb_experiments`` -- number of experimets

    - ``nb_iterations`` -- the number of iterations of the Rauzy-Zorich
      induction to perform

    - ``verbose`` -- if ``True`` print additional information concerning the
      mean and standard deviation
    """
    cdef int *p, *t   # permutation, twin, sigma
    cdef size_t **s
    cdef size_t *tab
    cdef size_t nc
    cdef size_t *dim
    cdef generalized_permutation *gp_c
    cdef quad_cover *qcc
    cdef double * theta
    cdef double *proj

    n = len(gp)//2

    # convert the data of into C values
    p = <int *> malloc(2*n * sizeof(int))
    t = <int *> malloc(2*n * sizeof(int))
    s = <size_t **> malloc(n * sizeof(size_t*))

    for i from 0 <= i < 2*n:
        p[i] = gp[i]
        t[i] = twin[i]
    
    for i from 0 <= i < n:
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

    res = [[] for _ in xrange(nn)]

    c_isnan, c_isinf, tot_isnan, tot_isinf = 0, 0, 0, 0

    if nb_vectors == 1:
        for i in xrange(nb_experiments):
            top_lyapunov_exponents_H_plus(qcc, theta, nb_iterations)
            if isnan(theta[0]) or isnan(theta[1]):
                c_isnan  += 1
            elif isinf(theta[0]) or isinf(theta[1]):
                c_isinf  += 1
            else:
                for j in xrange(2):
                    res[j].append(theta[j])
    else:
        init_GS(nb_vectors)
        for i in xrange(nb_experiments):
            if projections and isotopic_decomposition:
                dim = <size_t *> malloc(int(nb_char) * sizeof(size_t))
                nc = nb_char
                for i from 0 <= i < nb_char:
                    dim[i] = int(dimensions[i])
                proj = <double *> malloc((n * degree)**2 * nb_char * sizeof(double))
                for i from 0 <= i < (n * degree)**2 * nb_char:
                    proj[i] = <double> projections[i]
                lyapunov_exponents_isotopic(qcc, theta, nb_iterations, nc, dim, proj)
                free(dim)
                free(proj)
            else:
                lyapunov_exponents_H_plus(qcc, theta, nb_iterations)
                
            if any(isnan(theta[i]) for i in xrange(nn)): c_isnan  += 1
            elif any(isinf(theta[i]) for i in xrange(nn)): c_isinf  += 1
            else:
                for j in xrange(nn):
                    res[j].append(theta[j])
            if (i == nb_experiments-1) and (c_isnan + c_isinf < 9*nb_experiments/10):
                i = nb_experiments - (c_isnan+c_isinf) - 1
                tot_isnan += c_isnan
                tot_isinf += c_isinf
                c_isnan, c_isinf = 0, 0

    if tot_isnan > 0:
        print("Warning, " + str(tot_isnan) + " NaN results in the experiments")
    if tot_isinf > 0:
        print("Warning, " + str(tot_isinf) + " infinity results in the experiments")
    if (c_isnan + c_isinf >= 9*nb_experiments/10):
        raise NameError('too much NaN or inf results')

    for i in xrange(n):
        free(s[i])
    free(s)

    free_quad_cover(&qcc)
    free_GS()
    free(theta)

    return res
