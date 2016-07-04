#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

typedef struct
{
	size_t degree;      /* degree (nb of squares) */
	int *pa, *pb;       /* permutations */

	size_t nb_vectors;  /* nb_vectors (at least 2 and typically genus-1) */
	double **va, **vb;  /* vectors */
}origami_data;

origami_data * new_origami_data(size_t degree, size_t nb_vectors, int *pa, int *pb);

void free_origami_data(origami_data *o);

void lyapunov_exponents(origami_data *o, size_t NB_ITERATIONS, double *ttheta);

/***** Origami with involution *****/

typedef struct
{
	size_t degree;      /* degree (nb os squares) */
	int *pa, *pb;        /* permutations */

	int *s;                /* involution such that x = sxs^-1 */
	size_t nb_vectors_p;   /* nb vectors in the + part of the involution */
	size_t nb_vectors_m;   /* nb vectors in the - part of the involution */

	double **va, **vb; /* vectors */
	                   /* the first nb_vectors_p are assumed to be in the + part while */
	                   /* the next nb_vectors_m are assumed to be in the - part        */
}origami_with_involution_data;

origami_with_involution_data * new_origami_with_involution_data(size_t degree, size_t nb_vectors_p, size_t nb_vectors_m, int *pa, int *pb, int *s);

void free_origami_with_involution_data(origami_with_involution_data *o);

void set_random_lengths_and_vectors_with_involution(origami_with_involution_data *o);

void lyapunov_exponents_with_involution(origami_with_involution_data *o, size_t NB_ITERATIONS, double * ttheta);

