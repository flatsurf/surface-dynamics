/* Some highly non optimized linear algebra routines */
/* 1) random orthogonal vectors                      */
/* 2) orthogonalization                              */

#include "lyapunov_exponents.h"

/* 0 <= i < nb_vectors   */
/* 0 <= j < degree       */
/* 0 <= k < nb_intervals */
/*  index at pos (i,j,k) = v[k + nb_intervals * (j + degree * i)] */

/* global variable used in GS orthogonalisation                                 */
/* they must be dynamic so we keep them here in order to avoid many malloc/free */
double * scal     = NULL;
double * scal_new = NULL;
size_t scal_size  = 0;

void print_vectors(quad_cover *qcc)
{
	size_t i,j,k;

	for(i=0; i<qcc->nb_vectors; ++i)
	{
		printf("v[%zu]=",i);
		for(j=0; j < qcc->nb_labels; ++j)
		{
			for(k=0; k < qcc->degree; ++k)
				printf(" %f", (qcc->labels)[j].v[i + qcc->nb_vectors * k]);
			printf(" |");
		}
		printf("\n");
	}
}

void set_random_vectors(quad_cover *qcc)
/* set random orthogonal frame */
/* warning: before calling this function call init_GS(dim) */
{
	size_t i,j;
	double *vv,norm;

	for(j=0; j<qcc->nb_labels; ++j)
	{
		vv = (qcc->labels)[j].v;
		for(i=0; i<qcc->nb_vectors*qcc->degree; i++)
			vv[i] = drand() - .5;
	}
	if(qcc->nb_vectors == 1)
	{
		for(i=0; i<qcc->nb_labels; i++)
			for(j=0; j<qcc->degree; j++)
			{
				if((qcc->labels)[i].v[j] > 0) norm += (qcc->labels)[i].v[j];
				else norm -= (qcc->labels)[i].v[j];
			}
	}
	else orthogonalize_GS(qcc, NULL);
}

int init_GS(size_t nb_vectors)
/* allocate scal and scal_new in order that it contains at least 2*nb_vectors elements */
{
	if(scal == NULL)
	{
		scal_size = 256;
		/* warning: scal and scal_new are switched in the GS orthogonalization */
		/* so it is safer to allocate them independently                       */
		scal = (double *) malloc(scal_size * sizeof(double));
		scal_new = (double *) malloc(scal_size * sizeof(double));
		//printf("info: GS allocation with size %lu at %p and %p\n",scal_size, scal, scal_new);
	}

	else if(scal_size < nb_vectors)
	{
		while(scal_size < nb_vectors) scal_size = (scal_size * 3) / 2;
		scal = (double *) realloc(scal, scal_size * sizeof(double));
		scal_new = (double *) realloc(scal_new, scal_size * sizeof(double));
		//printf("info: GS reallocation with size %lu at %p and %p\n",scal_size,scal,scal_new);
	}

	if(scal == NULL)
	{
		scal_size = 0;
		return 1;
	}

	return 0;
}

void free_GS(void)
{
	scal_size = 0;
	if(scal != NULL)
	{
		free(scal);
		free(scal_new);
		scal = NULL;
		scal_new = NULL;
	}
}

void orthogonalize_GS(quad_cover *qcc, double *theta)
/* orthogonalization using Gram-Schmidt process                                                */
/* warning: this method should be called AFTER being sure that init_GS has been called         */
/* it theta is not NULL it is updated by 2*log of the diagonal from position 1 to nb_vectors   */
/* warning: it is 2*log(entry) and not log(entry)                                              */
/* element at position (i,j,k): (qcc->labels)[j].v[i + nb_vectors*k];                          */
{
	size_t i,ii,j,k;

	double norm,sqnorm,c;
	double *tmp=NULL;
	double *vv;

	/* some check that we will remove */
	if(qcc->nb_vectors < 2)
	{
		fprintf(stderr, "calling GS with nb_vectors < 2 is not possible.\n");
		exit(EXIT_FAILURE);
	}
	if(qcc->nb_vectors > qcc->degree * qcc->nb_labels)
	{
		fprintf(stderr, "too much vectors for the dimension!\n");
		exit(EXIT_FAILURE);
	}
	if(scal == NULL || scal_new == NULL)
	{
		fprintf(stderr, "you must call init_GS before calling orthogonalize_GS.\n");
		exit(EXIT_FAILURE);
	}

	/* put v in the holonomy free subspace */
	/* compute <v_0,v_1> in scal_new[0] */
	/* compute <v_0,v_0> in sqnorm */
	scal_new[0] = 0.;
	sqnorm = 0.;
	for(j=0; j < qcc->nb_labels; ++j)
	{
		vv = (qcc->labels)[j].v;
		for(k=0; k<qcc->degree; ++k)
		{
			scal_new[0] += vv[0 + qcc->nb_vectors*k] * vv[1 + qcc->nb_vectors*k];
			sqnorm  += vv[0 + qcc->nb_vectors*k] * vv[0 + qcc->nb_vectors*k];
		}
	}


	/* vector by vector orhtonormalisation */
	for(i=1; i < qcc->nb_vectors-1; ++i)
	{
		/* sqnorm contains <v_(i-1),v_(i-1)> */
		/* v_0, v_1, ..., v_(i-2) are normalized */
		/* v_0, v_1, ..., v_(i-1) are orthogonal */
		/* scal_new contains the i scalar products <v_0,v_i>, <v_1,v_i>, ..., <v_(i-1),v_i> */
		tmp = scal; scal = scal_new; scal_new = tmp;
		for(ii=0; ii<=i; ++ii)
			scal_new[ii] = 0.;

		c = scal[i-1]/sqnorm;
		if(theta != NULL) theta[i] += log(sqnorm);
		norm = sqrt(sqnorm);
		sqnorm = 0.;

		/* c = <v_i,v_(i-1)> / <v_(i-1,v_(i-1)> */
		for(j=0; j < qcc->nb_labels; ++j)
		{
			vv = (qcc->labels)[j].v;
			for(k=0; k < qcc->degree; ++k)
			{
				/* subtract the part of the span of v_0, v_1, ..., v_(i-2) */
				for(ii=0; ii<i-1; ++ii)
					vv[i + qcc->nb_vectors * k] -= scal[ii] * vv[ii + qcc->nb_vectors * k];

				/* subtract the part of the span of v_(i-1) */
				vv[i + qcc->nb_vectors * k] -= c * vv[ii + qcc->nb_vectors*k];

				if ( i + qcc->nb_vectors * k < 0 || i + qcc->nb_vectors * k >= qcc->nb_vectors*qcc->degree)
				  printf("NO : %zu,  %zu, %zu\n", i, k, i + qcc->nb_vectors * k);

				/* normalize v_(i-1) */
				vv[(i-1) + qcc->nb_vectors*k] /= norm;

				/* compute scalar products and norms for next loop */
				/* sqnorm = <v_i, v_i> */
				/* scal_new[ii] = <v_(i+1), v_ii>*/
				sqnorm += vv[i + qcc->nb_vectors * k] * vv[i + qcc->nb_vectors *k];
				for(ii=0; ii<=i; ++ii)
					scal_new[ii] += vv[ii + qcc->nb_vectors*k] * vv[i+1 + qcc->nb_vectors*k];
			}
		}
	}

	/* here i = NB_VECTORS-1 */
	/* renormalize v_(i-1)   */
	/* put v_i in the orthogonal of the span of v_0, v_1, ..., v_(i-1) */
	c = scal_new[i-1] / sqnorm;
	if(theta != NULL) theta[i] += log(sqnorm);
	norm = sqrt(sqnorm);
	sqnorm = .0;

	for(j=0; j< qcc->nb_labels; ++j)
	{
		vv = (qcc->labels)[j].v;
		for(k=0; k<qcc->degree; ++k)
		{
			for(ii=0; ii<i-1; ++ii)
				vv[i + qcc->nb_vectors*k] -= scal_new[ii] * vv[ii + qcc->nb_vectors*k];

				if ( i + qcc->nb_vectors * k < 0 || i + qcc->nb_vectors * k >= qcc->nb_vectors*qcc->degree)
				  printf("NO : %zu,  %zu, %zu\n", i, k, i + qcc->nb_vectors * k);
			vv[i + qcc->nb_vectors*k] -= c * vv[i-1 + qcc->nb_vectors*k];

			vv[i-1 + qcc->nb_vectors * k] /= norm;

			sqnorm += vv[i + qcc->nb_vectors * k] * vv[i + qcc->nb_vectors * k];
		}
	}
	
	/* we renormalize v_i*/
	if(theta != NULL) theta[i+1] += log(sqnorm);
	norm = sqrt(sqnorm);
	for(j=0; j < qcc->nb_labels; ++j)
	{
		vv = (qcc->labels)[j].v;
		for(k=0; k < qcc->degree; ++k)
			vv[i + qcc->nb_vectors * k] /= norm;
	}
}


void orthogonalize_iso(quad_cover *qcc, double *theta, size_t nb_char, size_t* dimensions)
/* orthogonalization using Gram-Schmidt process                                                */
/* warning: this method should be called AFTER being sure that init_GS has been called         */
/* it theta is not NULL it is updated by 2*log of the diagonal from position 1 to nb_vectors   */
/* warning: it is 2*log(entry) and not log(entry)                                              */
/* element at position (i,j,k): (qcc->labels)[j].v[i + nb_vectors*k];                          */
{
  int i,ii,j,k, i_char, sum_dimensions=0;

  double norm,sqnorm,c;
  double *tmp=NULL;
  double *vv;
	
  for (i_char=0; i_char<nb_char; ++i_char) {

    /* first deal with dim 1 */
    if(dimensions[i_char] == 1) {
      sqnorm = 0.;
      for(j=0; j < qcc->nb_labels; ++j)
	{
	  vv = (qcc->labels)[j].v;
	  for(k=0; k<qcc->degree; ++k)
	    sqnorm  += vv[sum_dimensions + qcc->nb_vectors*k] * vv[sum_dimensions + qcc->nb_vectors*k];
	}

      /* we renormalize v*/
      if(theta != NULL) theta[sum_dimensions + 1] += log(sqnorm);
      norm = sqrt(sqnorm);
      for(j=0; j < qcc->nb_labels; ++j)
	{
	  vv = (qcc->labels)[j].v;
	  for(k=0; k < qcc->degree; ++k)
	    vv[sum_dimensions + qcc->nb_vectors * k] /= norm;
	}
    }

    if(dimensions[i_char]>1) {

      /* some check that we will remove */
      if(qcc->nb_vectors > qcc->degree * qcc->nb_labels)
	{
	  fprintf(stderr, "too much vectors for the dimension!\n");
	  exit(EXIT_FAILURE);
	}
      if(scal == NULL || scal_new == NULL)
	{
	  fprintf(stderr, "you must call init_GS before calling orthogonalize_GS.\n");
	  exit(EXIT_FAILURE);
	}

	/* put v in the holonomy free subspace */
	/* compute <v_0,v_1> in scal_new[0] */
	/* compute <v_0,v_0> in sqnorm */
	scal_new[0] = 0.;
	sqnorm = 0.;
	for(j=0; j < qcc->nb_labels; ++j)
	{
		vv = (qcc->labels)[j].v;
		for(k=0; k<qcc->degree; ++k)
		{
			scal_new[0] += vv[sum_dimensions + qcc->nb_vectors*k] * vv[sum_dimensions + 1 + qcc->nb_vectors*k];
			sqnorm  += vv[sum_dimensions + qcc->nb_vectors*k] * vv[sum_dimensions + qcc->nb_vectors*k];
		}
	}


	for(i=1; i < dimensions[i_char]-1; ++i)
	  {
	    /* vector by vector orhtonormalisation */
	    /* sqnorm contains <v_(i-1),v_(i-1)> */
	    /* v_0, v_1, ..., v_(i-2) are normalized */
	    /* v_0, v_1, ..., v_(i-1) are orthogonal */
	    /* scal_new contains the i scalar products <v_0,v_i>, <v_1,v_i>, ..., <v_(i-1),v_i> */
	    tmp = scal; scal = scal_new; scal_new = tmp;
	    for(ii=0; ii<=i; ++ii)
	      scal_new[ii] = 0.;
	
	    c = scal[i-1]/sqnorm;
	    if(theta != NULL) theta[sum_dimensions + i] += log(sqnorm);
	    norm = sqrt(sqnorm);
	    sqnorm = 0.;

	    /* c = <v_i,v_(i-1)> / <v_(i-1,v_(i-1)> */
	    for(j=0; j < qcc->nb_labels; ++j)
	      {
		vv = (qcc->labels)[j].v;
		for(k=0; k < qcc->degree; ++k)
		  {
		    /* subtract the part of the span of v_0, v_1, ..., v_(i-2) */
		    for(ii=0; ii<i-1; ++ii)
		      vv[sum_dimensions + i + qcc->nb_vectors * k] -= scal[ii] * vv[sum_dimensions + ii + qcc->nb_vectors * k];

		    /* subtract the part of the span of v_(i-1) */
		    vv[sum_dimensions + i + qcc->nb_vectors * k] -= c * vv[sum_dimensions + ii + qcc->nb_vectors*k];

		    /* normalize v_(i-1) */
		    vv[sum_dimensions + (i-1) + qcc->nb_vectors*k] /= norm;

		    /* compute scalar products and norms for next loop */
		    /* sqnorm = <v_i, v_i> */
		    /* scal_new[ii] = <v_(i+1), v_ii>*/
		    sqnorm += vv[sum_dimensions + i + qcc->nb_vectors * k] * vv[sum_dimensions + i + qcc->nb_vectors *k];
		    for(ii=0; ii<=i; ++ii)
		      scal_new[ii] += vv[sum_dimensions + ii + qcc->nb_vectors*k] * vv[sum_dimensions + i+1 + qcc->nb_vectors*k];
		  }
	      }
	  }
	
	/* here i = NB_VECTORS-1 */
	/* renormalize v_(i-1)   */
	/* put v_i in the orthogonal of the span of v_0, v_1, ..., v_(i-1) */
	c = scal_new[i-1] / sqnorm;
	if(theta != NULL) theta[sum_dimensions + i] += log(sqnorm);
	norm = sqrt(sqnorm);
	sqnorm = .0;
	
	for(j=0; j< qcc->nb_labels; ++j)
	  {
	    vv = (qcc->labels)[j].v;
	    for(k=0; k<qcc->degree; ++k)
	      {
		for(ii=0; ii<i-1; ++ii)
		  vv[sum_dimensions + i + qcc->nb_vectors*k] -= scal_new[ii] * vv[sum_dimensions + ii + qcc->nb_vectors*k];

		vv[sum_dimensions + i + qcc->nb_vectors*k] -= c * vv[sum_dimensions + i-1 + qcc->nb_vectors*k];
		
		vv[sum_dimensions + i-1 + qcc->nb_vectors * k] /= norm;

		sqnorm += vv[sum_dimensions + i + qcc->nb_vectors * k] * vv[sum_dimensions + i + qcc->nb_vectors * k];
	      }
	  }
	
	/* we renormalize v_i*/
	if(theta != NULL) theta[sum_dimensions + i+1] += log(sqnorm);
	norm = sqrt(sqnorm);
	for(j=0; j < qcc->nb_labels; ++j)
	  {
	    vv = (qcc->labels)[j].v;
	    for(k=0; k < qcc->degree; ++k)
	      vv[sum_dimensions + i + qcc->nb_vectors * k] /= norm;
	  }
    }
    sum_dimensions += dimensions[i_char];
  }
}
void check_orthogonality(quad_cover *qcc)
{
	double s;
	double *vv;
	size_t i1,i2,j,k;

	for(i1=0; i1<qcc->nb_vectors; ++i1)
	{
		for(i2=i1; i2<qcc->nb_vectors; ++i2)
		{
			s = 0;
			for(j = 0; j < qcc->nb_labels; ++j)
			{
				vv = (qcc->labels)[j].v;
				for(k=0; k < qcc->degree; ++k)
					s += vv[i1 + qcc->nb_vectors*k] * vv[i2 + qcc->nb_vectors*k];
			}
			printf("<v%d,v%d> = %f\t",(int) i1, (int) i2, s);
		}
		printf("\n");
	}
}


void check_orthogonality_iso(quad_cover *qcc, size_t nb_char, size_t* dimensions)
{
	double s;
	double *vv;
	size_t i1,i2,j,k, i_char;
	size_t sum_dimensions=0;

	for (i_char=0; i_char<nb_char; ++i_char)
	  sum_dimensions += dimensions[i_char];

	if (sum_dimensions != qcc->nb_vectors) {
	  fprintf(stderr, "Wrong isotypic decomposition : dimensions doesn't match\nSum of isotypic dimensions : %zu,   nb_vectors : %zu",
		  sum_dimensions, qcc->nb_vectors);
	  exit(EXIT_FAILURE);
	}
	
	sum_dimensions = 0;
	
	for (i_char=0; i_char<nb_char; ++i_char){
	  for(i1=sum_dimensions; i1<sum_dimensions + dimensions[i_char]; ++i1)
	    {
	      for(i2=i1; i2<sum_dimensions + dimensions[i_char]; ++i2)
		{
		  s = 0;
		  for(j = 0; j < qcc->nb_labels; ++j)
		    {
		      vv = (qcc->labels)[j].v;
		      for(k=0; k < qcc->degree; ++k)
			s += vv[i1 + qcc->nb_vectors*k] * vv[i2 + qcc->nb_vectors*k];
		    }

		  if (i1==i2 && (s - 1) > 0x40000000000){
		    fprintf(stderr, "Wrong normalisation\nNorm : %lf\n", s);
		    exit(EXIT_FAILURE);
		  }

		  if (i1!=i2 && s > 0x40000000000) {
		    fprintf(stderr, "Wrong orthogonalisation\nScalar product : %lf\n", s);
		    exit(EXIT_FAILURE);
		  }
		}
	    }
	  sum_dimensions += dimensions[i_char];
	}
}


inline void project_isotypic(quad_cover *qcc, size_t nb_char, size_t* dimensions, double* projections)
{
  /*
    We compute the product of the the nb_char matrices of projection, of size : size_of_matrix = nb_labels*degree
    M_(i_char)[(lab_i, deg_i), (lab_j, deg_j)] = projections [i_char*(size_of_matrix**2) + (deg_i*nb_labels + lab_i)*size_of_matrix + (deg_j*nb_labels + lab_j)]

    So the product is
    M*v(lab_i, deg_i) = sum (lab_j, deg_j)  M[(lab_i, deg_i),(lab_j, deg_j)] * v[lab_j, deg_j]

    We project the dimension[0] first vectors with M_0
    Then the dimension[1] next with M_1 and so on

    v_k[lab_i, deg_i] = labels[lab_i].v[ k + nb_vector * deg_i]

    Hence :
    M_k * v_l (lab_i, deg_i) = sum (lab_j, deg_j) projections[k*(size_of_matrix**2) + (deg_i*nb_labels + lab_i)*size_of_matrix + (deg_j*nb_labels + lab_j)]
                                                  * labels[lab_j].v[ l + nb_vector * deg_j]
   */
  size_t i_char, i_vec, lab_i, lab_j, deg_i, deg_j, sum_dimensions;

  size_t size_of_matrix = qcc->degree * qcc->nb_labels;
  size_t sq_size = size_of_matrix * size_of_matrix;

  double res;

  sum_dimensions = 0;

  for(i_char=0; i_char<nb_char; ++i_char){
    for(i_vec=0; i_vec<dimensions[i_char]; ++i_vec){
      for(lab_i=0; lab_i < qcc->nb_labels; ++lab_i)
	for(deg_i=0; deg_i < qcc->degree; ++deg_i)
	  {
	    res = 0.;
	    for(lab_j=0; lab_j<qcc->nb_labels; ++lab_j)
	      for (deg_j=0; deg_j<qcc->degree; ++deg_j)
		res += (double) projections[ i_char * (sq_size) +
					     (deg_i * qcc->nb_labels + lab_i) * size_of_matrix +
					     (deg_j * qcc->nb_labels + lab_j)]
		  * (qcc->labels)[lab_j].v[sum_dimensions + i_vec + qcc->nb_vectors * deg_j];
	    (qcc->v_buffer)[lab_i + qcc->nb_labels * deg_i] = res;
	    if (lab_i + qcc->nb_labels * deg_i >= qcc->nb_labels * qcc->degree || lab_i + qcc->nb_labels * deg_i < 0)
	      printf("ind : %zu, max : %zu\n", lab_i + qcc->nb_labels * deg_i, qcc->nb_labels * qcc->degree);
	  }

      for(lab_i=0; lab_i < qcc->nb_labels; ++lab_i)
	for(deg_i=0; deg_i < qcc->degree; ++deg_i)
	  (qcc->labels)[lab_i].v[sum_dimensions + i_vec + qcc->nb_vectors * deg_i] = (qcc->v_buffer)[lab_i + qcc->nb_labels * deg_i];
    }
    sum_dimensions += dimensions[i_char];
  }
}


inline double max_norm(quad_cover *qcc)
{
  double sqnorm = 0., max = 0.;
  size_t i_vec, lab, deg;

  for(i_vec=0; i_vec<qcc->nb_vectors; ++i_vec) {
    for(lab=0; lab < qcc->nb_labels; ++lab)
      for(deg=0; deg < qcc->degree; ++deg)
	sqnorm += (qcc->labels)[lab].v[i_vec + qcc->nb_vectors * deg] * (qcc->labels)[lab].v[i_vec + qcc->nb_vectors * deg];
    if (sqnorm > max)
      max = sqnorm;
    sqnorm = 0;
  }

  return max;
}
