#include <lyapunov_exponents.h>

/*********************************************************/
/* limits */
#define DEGREE_MAX 600
#define NB_VECTORS_MAX 60
/*********************************************************/

/*********** algorithm parameters **************************/
/* the size of theta[0] at which we add partial sum       */
/* theta in ttheta.                                       */
#define PARTIAL_THETA_SUM 0X1000  /* 4096 = 2**12: which should give number of the order of 700 */

/* the length of the interval under which we perform the  */
/* the renormalization step                               */
#define MIN_LENGTH    .000008L  /* roughly between 4 and 20 iterations of Gauss map */
#define CRITIC_LENGTH .0000000000000000000001L
/**********************************************************/

/* verbosity level */
#define _VERBOSE_ 0

/* double and long double log_e(2) */
# define LOG2      0.69314718055994530942
# define LOG2l     0.6931471805599453094172321214581766L

/*TODO: Kahan summation to preserve low bits (see Wikipedia) (... forget it) */

/* what is the expectation of the maximum partial quotient of order < N for a random initial point */

/* on intel centrino                */
/* RAND_MAX = INT_MAX  =  2^31-1    */
/* sizeof(int)         =  4*8 = 32  */
/* sizeof(long)        =  8*8 = 64  */
/* sizeof(double)      =  8*8 = 64  */
/* sizeof(long double) = 16*8 = 128 */


/* random uniform double in [0,1] */
static double drand()
{
	return (((double) random() + (double) random() / ( (double) RAND_MAX) ) / ((double) RAND_MAX+1) );
}

/* random uniform long double in [0,1] */
static long double ldrand()
{
	return ((long double) random() + (((long double) random() + (long double) random() / ((long double) RAND_MAX)) / ((long double) RAND_MAX))) / ((long double) RAND_MAX+1);
}

/* Gauss distribution in [0,1] with double precision */
static double drandGauss()
{
	return exp(LOG2 * drand()) - 1.;
}

/* Gauss distribution in [0,1] with long double precision */
static double ldrandGauss()
{
	return expl(LOG2l * ldrand()) - 1.;
}

inline size_t lcm(size_t a, size_t b)
/* the least common multiple of a and b */
{
	size_t aa = a, r=b, bb;

	while(bb = r)
	{
		r = aa%bb;
		aa = bb;
	}

	return (a/aa)*b;
}

void invariant_projection(double **v, size_t DEGREE, size_t NB_VECTORS, int *s)
/* put the vector v[i] for i=0,...,NB_VECTORS in the invariant subspace of s */
{
	size_t i,j;

#if _VERBOSE_ > 10
	printf("Projection on E+ of the %d vectors at %x ",NB_VECTORS,v);
	printf("with respect to the involution:");
	for(i=0;i<DEGREE;++i) printf("%d ",s[i]);
	printf("\n");
#endif

	for(j=0; j<DEGREE; ++j)
		if(j < s[j])
		{
#if _VERBOSE_ > 10
				printf(" project (j,s[j]) = (%d,%d)\n",j,s[j]);
#endif
			for(i=0; i<NB_VECTORS; ++i)
			{
				v[i][j] = (v[i][j] + v[i][s[j]]) / 2.0;
				v[i][s[j]] = v[i][j];
			}

		}
}

void anti_invariant_projection(double **v, size_t DEGREE, size_t NB_VECTORS, int *s)
/* put the vector v[i] for i=0,...,NB_VECTORS in the anti-invariant subspace of s */
{
	size_t i,j;

#if _VERBOSE_ > 10
	printf("Projection on E- of the %d vectors at %x ",NB_VECTORS,v);
	printf("with respect to the involution:");
	for(i=0;i<DEGREE;++i) printf("%d ",s[i]);
	printf("\n");
#endif

	for(j=0; j<DEGREE; ++j)
	{
		if(j == s[j])
			for(i=0; i<NB_VECTORS; ++i) v[i][j] = 0.;

		if(j < s[j])
		{
			for(i=0; i<NB_VECTORS; ++i)
			{
				v[i][j] = (v[i][j] - v[i][s[j]]) / 2.0;
				v[i][s[j]] = -v[i][j];
			}
		}
	}
}

void QR(double **va, double **vb, size_t DEGREE, size_t NB_VECTORS, double *theta)
/*
 * QR methods applied to the vectors given by va and vb. A vector for that function
 * is v[i] = (va[i],vb[i]) for i=0,...,NB_VECTORS. The number of components for
 * va[i] and vb[i] is DEGREE.
 *
 *  - modify va[i] and vb[i] to make them holonomy free and orthonormal
 *  - if not NULL, update theta[i] by the log of the norm of the i-th component
 *
 * WARNING: NB_VECTORS must be at least 2 (note that for NB_VECTORS=1,
 * there is no need for orthonormalization)
 */
{
	size_t i,ii,j;

	double hol_a[DEGREE_MAX],hol_b[DEGREE_MAX];
	double sscal[2*DEGREE_MAX];
	double *scal = sscal;
	double *scal_new = sscal + DEGREE_MAX;
	double norm,sqnorm,c;
	double *tmp=NULL;

	/* compute holonomy in (hol_a[i],hol_b[i]) */
	for(i=0; i<NB_VECTORS; ++i)
	{
		hol_a[i] = va[i][0];
		hol_b[i] = vb[i][0];
		for(j=1; j<DEGREE; ++j)
		{
			hol_a[i] += va[i][j];
			hol_b[i] += vb[i][j];
		}
		hol_a[i] /= DEGREE;
		hol_b[i] /= DEGREE;
	}

	/* put v in the holonomy free subspace */
	/* compute <v_0,v_1> in scal_new[0] */
	/* compute <v_0,v_0> in sqnorm */
	scal_new[0] = 0.;
	sqnorm = 0.;
	for(j=0; j < DEGREE; ++j)
	{
		for(i=0; i < NB_VECTORS; ++i)
		{
			va[i][j] -= hol_a[i];
			vb[i][j] -= hol_b[i];
		}
		scal_new[0] += va[0][j] * va[1][j] + vb[0][j] * vb[1][j];
		sqnorm  += va[0][j] * va[0][j] + vb[0][j] * vb[0][j];
	}
	/* vector by vector orhtonormalisation */
	for(i=1; i < NB_VECTORS-1; ++i)
	{
		/* sqnorm contains <v_(i-1),v_(i-1)> */
		/* v_0, v_1, ..., v_(i-2) are normalized */
		/* v_0, v_1, ..., v_(i-1) are orthogonal */
		/* scal_new contains the i vectors <v_0,v_i>, <v_1,v_i>, ..., <v_(i-1),v_i> */
		tmp = scal; scal = scal_new; scal_new = tmp;
		for(ii=0; ii<=i; ++ii)
			scal_new[ii] = 0.;

		c = scal[i-1]/sqnorm;
		if(theta != NULL) theta[i] += log(sqnorm);
		norm = sqrt(sqnorm);
		sqnorm = 0.;

		/* c = <v_i,v_(i-1)> / <v_(i-1,v_(i-1)> */
		for(j=0; j<DEGREE; ++j)
		{
			/* subtract the part of the span of v_0, v_1, ..., v_(i-2) */
			for(ii=0; ii<i-1; ++ii)
			{
				va[i][j] -= scal[ii] * va[ii][j];
				vb[i][j] -= scal[ii] * vb[ii][j];
			}
			/* subtract the part of the span of v_(i-1) */
			va[i][j] -= c * va[i-1][j];
			vb[i][j] -= c * vb[i-1][j];

			/* normalize v_(i-1) */
			va[i-1][j] /= norm;
			vb[i-1][j] /= norm;

			/* compute scalar products and norms for next loop */
			/* sqnorm = <v_i, v_i> */
			/* scal_new[ii] = <v_(i+1), v_ii>*/
			sqnorm += va[i][j]*va[i][j] + vb[i][j]*vb[i][j];
			for(ii=0; ii<=i; ++ii)
				scal_new[ii] += va[ii][j]*va[i+1][j] + vb[ii][j]*vb[i+1][j];
		}

	}

	/* here i = NB_VECTORS-1 */
	/* renormalize v_(i-1)   */
	/* put v_i in the orthogonal of the span of v_0, v_1, ..., v_(i-1) */
	c = scal_new[i-1] / sqnorm;
	if(theta != NULL) theta[i] += log(sqnorm);
	norm = sqrt(sqnorm);
	sqnorm = .0;

	for(j=0; j<DEGREE; ++j)
	{
		for(ii=0; ii<i-1; ++ii)
		{
			va[i][j] -= scal_new[ii] * va[ii][j];
			vb[i][j] -= scal_new[ii] * vb[ii][j];
		}
		va[i][j] -= c * va[i-1][j];
		vb[i][j] -= c * vb[i-1][j];

		va[i-1][j] /= norm;
		vb[i-1][j] /= norm;

		sqnorm += va[i][j]*va[i][j] + vb[i][j]*vb[i][j];
	}

	/* we renormalize v_i*/
	if(theta != NULL) theta[i+1] += log(sqnorm);
	norm = sqrt(sqnorm);
	for(j=0; j < DEGREE; ++j)
	{
		va[i][j] /= norm;
		vb[i][j] /= norm;
	}

}

/*TODO: clean allocation */
origami_data * new_origami_data(size_t degree, size_t nb_vectors, int *pa, int *pb)
{
	size_t i;
	origami_data * o = (origami_data *) malloc(sizeof(origami_data));

	o->degree = degree;
	o->nb_vectors = nb_vectors;

	if(((o->pa = (int *) malloc(degree * sizeof(int))) == NULL) ||
	   ((o->pb = (int *) malloc(degree * sizeof(int))) == NULL))
		return NULL;

	for(i=0; i<degree; ++i)
	{
		o->pa[i] = pa[i];
		o->pb[i] = pb[i];
	}

	if(((o->va = (double **) malloc(nb_vectors * sizeof(double *))) == NULL) ||
	   ((o->vb = (double **) malloc(nb_vectors * sizeof(double *))) == NULL))
		return NULL;

	for(i=0; i<nb_vectors; ++i)
		if(((o->va[i] = (double *) malloc(degree * sizeof(double))) == NULL) ||
		   ((o->vb[i] = (double *) malloc(degree * sizeof(double))) == NULL))
			return NULL;

	return o;
}

void free_origami_data(origami_data *o)
{
	size_t i;

	if(o != NULL)
	{
		if(o->pa != NULL) free(o->pa);
		if(o->pb != NULL) free(o->pb);

		if(o->va != NULL)
		{
			for(i=0; i<o->nb_vectors;++i)
				if(o->va[i] != NULL) free(o->va[i]);
			free(o->va);
		}
		if(o->vb != NULL)
		{
			for(i=0; i<o->nb_vectors;++i)
				if(o->vb[i] != NULL) free(o->vb[i]);
			free(o->vb);
		}

		free(o);
	}
}

void set_random_lengths_and_vectors(origami_data *o)
{
	size_t i,j;

	/* initialize random length */
	for(i=0; i < o->nb_vectors; ++i)
		for(j=0; j < o->degree; ++j)
		{
			o->va[i][j] = drand();
			o->vb[i][j] = drand();
		}

	QR(o->va, o->vb, o->degree, o->nb_vectors, NULL);
}


void lyapunov_exponents(origami_data *o, size_t NB_ITERATIONS, double * ttheta)
/* one trajectory of Lyapunov exponents */
/* return as many Lyapunov exponents as t contains vectors */
/* the permutations o->pa and o->pb are modified... they become the last element of the SL(2,Z)
 * orbit visited by this random walk */
{
	/* iterator */
	size_t iteration=0;
	size_t i,j,jj;

	/* local variable for cycles length */
	size_t m,mm;

	/* swap for the iet */
	long double ll;
	int *pp;
	double **vv;

	/* set data */
	set_random_lengths_and_vectors(o);
    long double la;
	long double lb;
	la = ldrandGauss();
	while(la > 0.5) la = ldrandGauss();
	lb = 1. - la;

	/* get data from the origami */
	size_t DEGREE = o->degree;
	size_t NB_VECTORS = o->nb_vectors;
	int *pa = o->pa;
	int *pb = o->pb;
	double **va = o->va;
	double **vb = o->vb;

	/* lyapunov exponents */
	double *theta = (double *) malloc((NB_VECTORS+1)*sizeof(double));

	/* data for the iet */
	size_t order;
	double vb_tot[NB_VECTORS_MAX][DEGREE_MAX];
	size_t cycles[DEGREE_MAX];

	/* initialize theta */
	for(i=0; i<NB_VECTORS+1; ++i)
		ttheta[i] = theta[i] = 0.0;

	while(iteration < NB_ITERATIONS)
	{
	    /* Rauzy induction while the length is not too small   */
		/* we assume lb > la at the begining of the loop       */
		/* one step of iteration corresponds to go from lb and */
		/* get smaller than MIN_LENGTH */
		while(lb > MIN_LENGTH)
		{
			/* compute the order of the perm pb */
			/* compute the vector vb_tot        */
			/* TODO: do not compute order but only */
			/* the length of the cycles            */
			for(j=0; j<DEGREE; cycles[j++] = DEGREE);
			order = 1;
			mm = 0;
			for(j=0; j<DEGREE; ++j)
			{
				if(cycles[j] == DEGREE)
				{
					cycles[j] = mm;
					for(i=0; i<NB_VECTORS; ++i)
						vb_tot[i][mm] = vb[i][j];
					m=1;
					jj=j;
					while((jj=pb[jj]) != j)
					{
						cycles[jj] = mm;
						for(i=0; i<NB_VECTORS; ++i)
							vb_tot[i][mm] += vb[i][jj];
						m++;
					}
					if(m != 1)
					{
						for(i=0; i<NB_VECTORS; ++i)
							vb_tot[i][mm] /= m;
						order = lcm(order,m);
					}
					mm++;
				}
			}
			/* length of the extended step */
			m = floor(lb / la);

			/* update length */
			lb -= m * la;

			/* update the integer part */
			mm = m % order;
			m = m - mm;
			for(i=0; i<NB_VECTORS; ++i)
				for(j=0; j<DEGREE; ++j)
					va[i][j] += m * vb_tot[i][cycles[pa[j]]];

			/* update the remainder part         */
			/* modify pa by precomposition by pb */
			for(m=0; m<mm; ++m)
				for(j=0; j<DEGREE; ++j)
				{
					jj = pa[j];
					for(i=0; i<NB_VECTORS; ++i)
						va[i][j] += vb[i][jj];
					pa[j] = pb[jj];
				}

			/* swap a and b */
			ll = lb; lb = la; la = ll;
			pp = pb; pb = pa; pa = pp;
			vv = vb; vb = va; va = vv;
		} /* end of induction */

		++iteration;

		/* TODO: stop the trajectory */
		if(lb < CRITIC_LENGTH) /* restart with random data */
		{
			fprintf(stderr, "Warning: reinitialize length\n");
			la = ldrandGauss();
			while(la > 0.5) la = ldrandGauss();
			lb = 1.-la;
		}

		else /* renormalization of length */
		{
			ll = la+lb;
			theta[0] -= log((double) ll);
			la /= ll;
			lb /= ll;
		}

		/* partial sum of theta in ttheta */
		if(theta[0] > PARTIAL_THETA_SUM)
		{
			for(i=0; i<NB_VECTORS+1; ++i)
			{
				ttheta[i] += theta[i];
				theta[i] = 0.0;
			}
		}

		/* orthonrmalization */
		QR(va, vb, DEGREE, NB_VECTORS, theta);

	} /* end of the iteration */

	for(i=0; i<NB_VECTORS+1; ++i)
		ttheta[i] += theta[i];

	free(theta);
}

/**********************************/
/* Origami with involution
typedef struct{
	size_t degree;
	int *pa,*pb;
	long double la, lb;

	size_t s;
	size_t nb_vectors_p;
	double **va_p, **vb_p;

	size_t nb_vectors_m;
	double **va_m, **vb_m;
}origami_with_involution_data;
*/

origami_with_involution_data * new_origami_with_involution_data(
		size_t degree, size_t nb_vectors_p, size_t nb_vectors_m,
		int *pa, int *pb, int *s)
{
	size_t i;
	size_t nb_vectors = nb_vectors_p + nb_vectors_m;

	origami_with_involution_data * o = (origami_with_involution_data *) malloc(sizeof(origami_with_involution_data));

	o->degree = degree;
	o->nb_vectors_p = nb_vectors_p;
	o->nb_vectors_m = nb_vectors_m;

	if(((o->pa = (int *) malloc(degree * sizeof(int))) == NULL) ||
	   ((o->pb = (int *) malloc(degree * sizeof(int))) == NULL) ||
	   ((o->s  = (int *) malloc(degree * sizeof(int))) == NULL))
		return NULL;

	for(i=0; i<degree; ++i)
	{
		o->pa[i] = pa[i];
		o->pb[i] = pb[i];
		o->s[i]  = s[i];
	}

	if(((o->va = (double **) malloc(nb_vectors * sizeof(double *))) == NULL) ||
	   ((o->vb = (double **) malloc(nb_vectors * sizeof(double *))) == NULL))
		return NULL;

	for(i=0; i<nb_vectors; ++i)
		if(((o->va[i] = (double *) malloc(degree * sizeof(double))) == NULL) ||
		   ((o->vb[i] = (double *) malloc(degree * sizeof(double))) == NULL))
			return NULL;

#if _VERBOSE_ > 20
	printf(" new done\n");
	printf("involution:");
	for(i=0;i<o->degree;++i) printf("%d ",o->s[i]);
	printf("\n");
#endif

	return o;
}

void free_origami_with_involution_data(origami_with_involution_data *o)
{
	size_t i;
	size_t nb_vectors = o->nb_vectors_p + o->nb_vectors_m;

#if _VERBOSE_ > 20
	printf("free...");
#endif

	if(o != NULL)
	{
		if(o->pa != NULL) free(o->pa);
		if(o->pb != NULL) free(o->pb);
		if(o->s  != NULL) free(o->s);

		if(o->va != NULL)
		{
			for(i=0; i<nb_vectors;++i)
				if(o->va[i] != NULL) free(o->va[i]);
			free(o->va);
		}
		if(o->vb != NULL)
		{
			for(i=0; i<nb_vectors;++i)
				if(o->vb[i] != NULL) free(o->vb[i]);
			free(o->vb);
		}

		free(o);
	}

#if _VERBOSE_ > 20
	printf("done\n");
#endif

}

void set_random_lengths_and_vectors_with_involution(origami_with_involution_data *o)
{
	size_t i,j;

#if _VERBOSE_ > 20
	printf("set_random_lengths_and_vectors_with_involution...\n");
	printf("o->va_p = %x, o->va_m = %x, sizeof(double *) = %x\n",o->va,o->va+o->nb_vectors_p,sizeof(double *));
#endif

	/* initialize random length */
	for(j=0; j < o->degree; ++j)
	{
		for(i=0; i < o->nb_vectors_p + o->nb_vectors_m; ++i)
		{
			o->va[i][j] = drand();
			o->vb[i][j] = drand();
		}
	}

#if _VERBOSE_ > 20
	printf("random vectors\n");
	for(i=0; i<o->nb_vectors_p; i++)
	{
		for(j=0; j<o->degree; j++)
			printf("(%f,%f) ",o->va[i][j],o->vb[i][j]);
		printf("\n");
	}
#endif

	if(o->nb_vectors_p)
	{
		invariant_projection(o->va, o->degree, o->nb_vectors_p, o->s);
		invariant_projection(o->vb, o->degree, o->nb_vectors_p, o->s);
		QR(o->va, o->vb, o->degree, o->nb_vectors_p, NULL);
	}

#if _VERBOSE_ > 20
	printf("random vectors after projection in E+ and QR\n");
	for(i=0; i<o->nb_vectors_p; i++)
	{
		for(j=0; j<o->degree; j++)
			printf("(%f,%f) ",o->va[i][j],o->vb[i][j]);
		printf("\n");
	}
#endif

#if _VERBOSE_ > 20
	printf("random vectors\n");
	for(i=o->nb_vectors_p; i<o->nb_vectors_p+o->nb_vectors_m; i++)
	{
		for(j=0; j<o->degree; j++)
			printf("(%f,%f) ",o->va[i][j], o->vb[i][j]);
		printf("\n");
	}
	printf("\n");
#endif


	if(o->nb_vectors_m)
	{
		anti_invariant_projection(o->va + o->nb_vectors_p, o->degree, o->nb_vectors_m, o->s);
		anti_invariant_projection(o->vb + o->nb_vectors_p, o->degree, o->nb_vectors_m, o->s);
		QR(o->va + o->nb_vectors_p, o->vb + o->nb_vectors_p, o->degree, o->nb_vectors_m, NULL);
	}

#if _VERBOSE_ > 20
	printf("random vectors after projection in E- and QR\n");
	for(i=o->nb_vectors_p; i<o->nb_vectors_p+o->nb_vectors_m; i++)
	{
		for(j=0; j<o->degree; j++)
			printf("(%f,%f) ",o->va[i][j], o->vb[i][j]);
		printf("\n");
	}
	printf("\n");
#endif
}


void lyapunov_exponents_with_involution(origami_with_involution_data *o, size_t NB_ITERATIONS, double * ttheta)
/* one trajectory of Lyapunov exponents */
/* return as many Lyapunov exponents as t contains vectors */
/* the permutations o->pa and o->pb are modified... they become the last element of the SL(2,Z)
 * orbit visited by this random walk */
{
	/* iterator */
	size_t iteration=0;
	size_t i,j,jj;

	/* local variable for cycles length */
	size_t m,mm;

	/* swap for the iet */
	long double ll;
	int *pp;
	double **vv;

	/* get data from the origami */
	set_random_lengths_and_vectors_with_involution(o);
	long double la = ldrandGauss();
	while(la > 0.5) la = ldrandGauss();
   	long double lb = 1.0 - la;

	size_t DEGREE = o->degree;
	size_t NB_VECTORS = o->nb_vectors_p + o->nb_vectors_m;
	int *pa = o->pa;
	int *pb = o->pb;
	double **va = o->va;
	double **vb = o->vb;

	/* lyapunov exponents */
	double *theta = (double *) malloc((NB_VECTORS+1)*sizeof(double));

	/* data for the iet */
	size_t order;
	double vb_tot[NB_VECTORS_MAX][DEGREE_MAX];
	size_t cycles[DEGREE_MAX];

#if _VERBOSE_ > 10
	printf("involution\n");
	for(i=0; i<DEGREE; i++)
		printf("%d ",(o->s)[i]);
	printf("\n");
	printf("vectors in E+\n");
	for(i=0; i<o->nb_vectors_p; i++)
	{
		for(j=0; j<DEGREE; j++)
			printf("(%f,%f) ",va[i][j],vb[i][j]);
		printf("\n");
	}
	printf("vectors in E-\n");
	for(i=o->nb_vectors_p; i<o->nb_vectors_p+o->nb_vectors_m; i++)
	{
		for(j=0; j<DEGREE; j++)
			printf("(%f,%f) ",va[i][j], vb[i][j]);
		printf("\n");
	}
#endif

	/* initialize theta */
#if _VERBOSE_ > 10
	printf("initialize ttheta and theta...");
#endif
	for(i=0; i<NB_VECTORS+1; ++i)
		ttheta[i] = theta[i] = 0.0;
#if _VERBOSE_ > 10
	printf(" done\n");
#endif

	while(iteration < NB_ITERATIONS)
	{
	    /* Rauzy induction while the length is not too small   */
		/* we assume lb > la at the begining of the loop       */
		/* one step of iteration corresponds to go from lb and */
		/* get smaller than MIN_LENGTH */
#if _VERBOSE_ > 1
		printf("iteration %d: ",iteration);
#endif
		while(lb > MIN_LENGTH)
		{
			/* compute the order of the perm pb */
			/* compute the vector vb_tot        */
			/* TODO: do not compute order but only */
			/* the length of the cycles            */
			for(j=0; j<DEGREE; cycles[j++] = DEGREE);
			order = 1;
			mm = 0;
			for(j=0; j<DEGREE; ++j)
			{
				if(cycles[j] == DEGREE)
				{
					cycles[j] = mm;
					for(i=0; i<NB_VECTORS; ++i)
						vb_tot[i][mm] = vb[i][j];
					m=1;
					jj=j;
					while((jj=pb[jj]) != j)
					{
						cycles[jj] = mm;
						for(i=0; i<NB_VECTORS; ++i)
							vb_tot[i][mm] += vb[i][jj];
						m++;
					}
					if(m != 1)
					{
						for(i=0; i<NB_VECTORS; ++i)
							vb_tot[i][mm] /= m;
						order = lcm(order,m);
					}
					mm++;
				}
			}
			/* length of the extended step */
			m = floor(lb / la);

			/* update length */
			lb -= m * la;

			/* update the integer part */
			mm = m % order;
			m = m - mm;
			for(i=0; i<NB_VECTORS; ++i)
				for(j=0; j<DEGREE; ++j)
					va[i][j] += m * vb_tot[i][cycles[pa[j]]];

			/* update the remainder part         */
			/* modify pa by precomposition by pb */
			for(m=0; m<mm; ++m)
				for(j=0; j<DEGREE; ++j)
				{
					jj = pa[j];
					for(i=0; i<NB_VECTORS; ++i)
						va[i][j] += vb[i][jj];
					pa[j] = pb[jj];
				}

			/* swap a and b */
			ll = lb; lb = la; la = ll;
			pp = pb; pb = pa; pa = pp;
			vv = vb; vb = va; va = vv;
#if _VERBOSE_ > 1
			printf("x");
#endif
		} /* end of induction */

#if _VERBOSE_ > 1
			printf("\n");
#endif

		++iteration;

		/* TODO: stop the trajectory */
		if(lb < CRITIC_LENGTH) /* restart with random data */
		{
			fprintf(stderr, "Warning: reinitialize length\n");
			la = ldrandGauss();
			while(la > 0.5) la = ldrandGauss();
			lb = 1.-la;
		}

		else /* renormalization of length */
		{
			ll = la+lb;
			theta[0] -= log((double) ll);
			la /= ll;
			lb /= ll;
		}

		/* partial sum of theta in ttheta */
		if(theta[0] > PARTIAL_THETA_SUM)
		{
			for(i=0; i<NB_VECTORS+1; ++i)
			{
				ttheta[i] += theta[i];
				theta[i] = 0.0;
			}
		}

		/* orthonormalization */
#if _VERBOSE_ > 5
		printf("vectors in E+ before QR-orthonormalization\n");
		for(i=0; i<o->nb_vectors_p; i++)
		{
			for(j=0; j<DEGREE; j++)
				printf("(%f,%f) ",va[i][j],vb[i][j]);
			printf("\n");
		}
#endif

		if(o->nb_vectors_p)
		{
			invariant_projection(va, DEGREE, o->nb_vectors_p, o->s);
			invariant_projection(vb, DEGREE, o->nb_vectors_p, o->s);
			QR(va, vb, DEGREE, o->nb_vectors_p, theta);
		}
/*		invariant_projection(va, DEGREE, o->nb_vectors_p, o->s);*/
/*		invariant_projection(vb, DEGREE, o->nb_vectors_p, o->s);*/

#if _VERBOSE_ > 5
		printf("vectors in E+ after QR-orthonormalization\n");
		for(i=0; i<o->nb_vectors_p; i++)
		{
			for(j=0; j<DEGREE; j++)
				printf("(%f,%f) ",va[i][j],vb[i][j]);
			printf("\n");
		}
		printf("vectors in E- before QR-orthonormalization\n");
		for(i=o->nb_vectors_p; i<NB_VECTORS; i++)
		{
			for(j=0; j<DEGREE; j++)
				printf("(%f,%f) ",va[i][j],vb[i][j]);
			printf("\n");
		}
#endif

		if(o->nb_vectors_m)
		{
			anti_invariant_projection(va + o->nb_vectors_p, DEGREE, o->nb_vectors_m, o->s);
			anti_invariant_projection(vb + o->nb_vectors_p, DEGREE, o->nb_vectors_m, o->s);
			QR(va+o->nb_vectors_p, vb+o->nb_vectors_p, DEGREE, o->nb_vectors_m, theta+o->nb_vectors_p);
		}
/*		anti_invariant_projection(va + o->nb_vectors_p, DEGREE, o->nb_vectors_m, o->s);	*/
/*		anti_invariant_projection(vb + o->nb_vectors_p, DEGREE, o->nb_vectors_m, o->s);*/


#if _VERBOSE_ > 5
		printf("vectors in E- after QR-orthonormalization\n");
		for(i=o->nb_vectors_p; i<NB_VECTORS; i++)
		{
			for(j=0; j<DEGREE; j++)
				printf("(%f,%f) ",va[i][j],vb[i][j]);
			printf("\n");
		}
		printf("\n");
#endif

	} /* end of the iteration */

	for(i=0; i<NB_VECTORS+1; ++i)
		ttheta[i] += theta[i];

	free(theta);
}

