/* computing Lyapunov exponents in H^+ */

#include "lyapunov_exponents.h"


int check_quad_cover(quad_cover * qcc)
{
  /*
    Check if everything is all right in our quad_cover qcc
  */
  size_t n_top=0,n_bot=0;
  long double l_top=0,l_bot=0;
  interval *i,*j;

  i = qcc->top;
  if(i == NULL)
    {
      fprintf(stderr, "no interval on the top!\n");
      return 1;
    }
  while(i != NULL)
    {
      n_top ++;
      l_top += i->lab->length;
      if((i->next != NULL) && (i->next->prev != i))
	{
	  fprintf(stderr, "Structure of the top line is broken\n");
	  return 1;
	}
      if((i->twin == i) || (i->lab != i->twin->lab))
	{
	  fprintf(stderr, "the twin is broken somwhere on the top line\n");
	  return 1;
	}
      j = i->twin;
      while(j->prev != NULL) j = j->prev;
      if((j != qcc->bot) && (j != qcc->top))
	{
	  fprintf(stderr, "oups!\n");
	  return 1;
	}
      if(i->lab->same_interval && (j != qcc->top))
	{
	  fprintf(stderr, "the field same_interval is not coherent\n");
	  return 1;
	}
      if(i->lab->length <= 0.000000000001l)
	{
	  fprintf(stderr, "critical length on the top for label=%zu length=%.16f\n",
		  i->lab - qcc->labels, (double) i->lab->length);
	  return 1;
	}
      i = i->next;
    }
  i = qcc->bot;
  if(i == NULL)
    {
      fprintf(stderr, "No interval on the bottom!\n");
      return 1;
    }
  while(i != NULL)
    {
      n_bot ++;
      l_bot += i->lab->length;
      if((i->next != NULL) && (i->next->prev != i))
	{
	  fprintf(stderr, "Structure of the bottom line is broken\n");
	  return 1;
	}
      if((i->twin == i) || (i->lab != i->twin->lab))
	{
	  fprintf(stderr, "the twin is borken on the bottom line\n");
	  return 1;
	}
      j = i->twin;
      while(j->prev != NULL) j = j->prev;
      if((j != qcc->bot) && (j != qcc->top))
	{
	  fprintf(stderr, "oups!\n");
	  return 1;
	}
      if(i->lab->same_interval && (j != qcc->bot))
	{
	  fprintf(stderr, "the field same_interval is not coherent\n");
	  return 1;
	}
      if(i->lab->length <= 0.0000000000001l)
	{
	  fprintf(stderr, "critical length on the bottom\n");
	  return 1;
	}
      i = i->next;
    }
  if(n_bot + n_top != 2 * qcc->nb_labels)
    {
      fprintf(stderr, "Number of labels do not coincide with the structure of intervals\n");
      return 1;
    }
  if((qcc->length - l_top > LENGTH_ERROR_TOLERANCE) || (l_top - qcc->length > LENGTH_ERROR_TOLERANCE))
    {
      fprintf(stderr, "The length of the top interval (=%lf) does not coincide with the length (=%lf)\n",(double) l_top,(double) qcc->length);
      return 1;
    }
  if((qcc->length - l_bot > LENGTH_ERROR_TOLERANCE) || (l_bot - qcc->length > LENGTH_ERROR_TOLERANCE))
    {
      fprintf(stderr, "The length of the bottom interval does not coincide with the length\n");
      return 1;
    }
  if((l_bot - l_top > LENGTH_ERROR_TOLERANCE) || (l_top - l_bot > LENGTH_ERROR_TOLERANCE))
    {
      fprintf(stderr, "The length of the top and bottom intervals do not match (get %lf and %lf)\n",(double) l_top, (double) l_bot);
      return 1;
    }

  return 0;
}


quad_cover * new_quad_cover(generalized_permutation *gp, size_t **sigma, size_t degree, size_t nb_vectors)

/*   Create a quad_cover of given degree over a translation surface, with the given generalized_permutation for the interval exchange, */
/*   moreover sigma is an array of permutations that gives the action of the fondamental group of the surface on the cover (we give the action */
/*   for the canonical generating familly of the group on the translation surface. */
/*   Finally, we set nb_vectors vectors, with degree*nb_intervals coordinates, to follow the monodromy afterward. */

{

  size_t i;
  quad_cover * qcc;

  qcc = (quad_cover *) malloc(sizeof(quad_cover));
  if(qcc == NULL) return NULL;

  qcc->intervals = (interval *) malloc(2 * gp->n * sizeof(interval));
  if(qcc->intervals == NULL)
    {
      free(qcc);
      return NULL;
    }

  qcc->labels = (label *) malloc(gp->n * sizeof(label));
  if(qcc->labels == NULL)
    {
      free(qcc->intervals);
      free(qcc);
      return NULL;
    }

  qcc->nb_labels = gp->n;
  qcc->degree = degree;
  qcc->nb_vectors = nb_vectors;

  for(i=0; i < gp->n; ++i) (qcc->labels)[i].v = NULL;

  for(i=0; i < gp->n; ++i)
    {
      (qcc->labels)[i].sigma = sigma[i];
      (qcc->labels)[i].v = (double *) malloc(degree * nb_vectors * sizeof(double));
      if(qcc->labels[i].v == NULL)
	{
	  free_quad_cover(&qcc); // this is safe
	  return NULL;
	}
    }

  qcc->buffer = (double *) malloc(degree * nb_vectors * sizeof(double));
  if(qcc->buffer == NULL)
    {
      free_quad_cover(&qcc);
      return NULL;
    }
  qcc->v_buffer = (double *) malloc(degree * gp->n * sizeof(double));
  if(qcc->v_buffer == NULL)
    {
      free_quad_cover(&qcc);
      return NULL;
    }

  for(i=0; i < 2*gp->n; ++i)
    {
      (qcc->labels)[(gp->perm)[i]].same_interval = ((i < gp->k) == ((gp->twin)[i] < gp->k));
      (qcc->intervals)[i].lab = qcc->labels + (gp->perm)[i];
      if((i == 0) || (i == gp->k))
	(qcc->intervals)[i].prev = NULL;
      else
	(qcc->intervals)[i].prev = qcc->intervals + i - 1;

      if((i == gp->k-1) || (i == 2*gp->n-1))
	(qcc->intervals)[i].next = NULL;
      else
	(qcc->intervals)[i].next = qcc->intervals + i + 1;

      (qcc->intervals)[i].orientation = 1;
      if(i < gp->k)
	{
	  if((gp->twin)[i] < i)
	    (qcc->intervals)[i].orientation = -1;
	}
      else
	{
	  if(((gp->twin)[i] < gp->k) || ((gp->twin)[i] > i))
	    (qcc->intervals)[i].orientation = -1;
	}

      (qcc->intervals)[i].twin = (qcc->intervals) + (gp->twin)[i];
      (qcc->intervals)[i].is_top = i < gp->k;
      (qcc->intervals)[i].give_name = ((qcc->intervals)[i].orientation == 1 && i < gp->k) || ((qcc->intervals)[i].orientation == -1 && i >= gp->k);
    }
  qcc->top = qcc->intervals + gp->k - 1;
  qcc->bot = qcc->intervals + 2*gp->n - 1;

  qcc->perm_one = (size_t *) malloc(degree * sizeof(size_t));
  qcc->perm_two = (size_t *) malloc(degree * sizeof(size_t));
  qcc->perm_buffer = (size_t *) malloc(degree * sizeof(size_t));

  return qcc;
}

inline int give_name(interval *inter)
{
  return ((inter->orientation == 1 && inter->is_top) || ((inter->orientation == -1 && !inter->is_top)));
}


void renormalize_length_quad_cover(quad_cover * qcc)
/* reset the length in the good subspace if needed */
/* warning: qcc->length should be 1                */
{
  long double ltop=0.,lbot=0.,lltop=0.,llbot=0.,length;
  interval * i;
  size_t j;
  
  i = qcc->top;
  while(i != NULL)
    {
      if(i->lab->same_interval) lltop += i->lab->length;
      ltop += i->lab->length;
      i = i->prev;
    }
  i = qcc->bot;
  while(i != NULL)
    {
      if(i->lab->same_interval) llbot += i->lab->length;
      lbot += i->lab->length;
      i = i->prev;
    }
  length = (lbot + ltop) / 2;
  
  /* rescale in order to have length 1 */
  for(j=0; j < qcc->nb_labels; ++j)
    (qcc->labels)[j].length /= length;
  qcc->length = 1;
  
  /* now if ltop and lbot are two far appart we project */
  if(((lbot - ltop)/(lbot + ltop) > EPSILON_LENGTH_PROJECTION) || ((ltop - lbot) / (ltop+lbot) > EPSILON_LENGTH_PROJECTION))
    {
      ltop /= length;
      lbot /= length;
      lltop /= length;
      llbot /= length;
      /* we modify the length */
      length = 1 + (1 - ltop) / lltop;
      if(length < 0)
	{
	  fprintf(stderr,"length error, can not renormalize\n");
	  exit(EXIT_FAILURE);
	}
      i = qcc->top;
      while(i != NULL)
	{
	  if((i->orientation == -1) && (i->lab->same_interval)) i->lab->length *= length;
	  i = i->prev;
	}
      length = 1 + (1 - lbot) / llbot;
      if(length < 0)
	{
	  fprintf(stderr, "length error, can not renormalize\n");
	  exit(EXIT_FAILURE);
	}
      i = qcc->bot;
      while(i != NULL)
	{
	  if((i->orientation == -1) && (i->lab->same_interval)) i->lab->length *= length;
	  i = i->prev;
	}
    }
}

void set_lengths(quad_cover *qcc, long double *lengths){
  size_t i;
  
  for (i=0; i< qcc->nb_labels; ++i)
    (qcc->labels)[i].length = (long double) lengths[i];
  qcc->length = 0;
  for (i=0; i< 2*qcc->nb_labels; ++i)
    qcc->length += (qcc->intervals)[i].lab->length;
  qcc->length /= 2;
}

void set_random_lengths_quad_cover(quad_cover * qcc)
{
  interval * i;
  size_t j,nb_odd_top=0,nb_odd_bot=0;
  
  qcc->length = 0.;
  for(j=0; j < qcc->nb_labels; j++)
    (qcc->labels)[j].length = 0;
  
  i = qcc->top;
  while(i != NULL)
    {
      nb_odd_top += i->lab->same_interval;
      i = i->prev;
    }
  i = qcc->bot;
  while(i != NULL)
    {
      nb_odd_bot += i->lab->same_interval;
      i = i->prev;
    }
  
  i = qcc->top;
  while(i != NULL)
    {
      if(i->lab->same_interval) i->lab->length += ((nb_odd_top + nb_odd_bot) * ldrand()) / (2 * nb_odd_top);
      else i->lab->length = ldrand();
      qcc->length += i->lab->length;
      i = i->prev;
    }
  
  i = qcc->bot;
  while(i != NULL)
    {
      if(i->lab->same_interval) i->lab->length += ((nb_odd_top + nb_odd_bot) * ldrand()) / (2 * nb_odd_bot);
      else i->lab->length = ldrand();
      qcc->length += i->lab->length;
      i = i->prev;
    }

  renormalize_length_quad_cover(qcc);
}

void randomize_length_quad_cover(quad_cover * qcc)
{
  qcc->top->lab->length += ((long double) (.5 - drand())) / 0x4000000000000; // we divide by 2^50 ~ 10^16
  qcc->bot->lab->length += ((long double) (.5 - drand())) / 0x4000000000000;
  renormalize_length_quad_cover(qcc);
}

void free_quad_cover(quad_cover ** qcc)
/* the value of *qcc is modified in order to prevent double free */
{
  size_t i;
  
  if(*qcc != NULL)
    {
      if((*qcc)->intervals != NULL)
	{
	  free((*qcc)->intervals);
			(*qcc)->intervals = NULL;
	}
		
      if((*qcc)->v_buffer != NULL)
	{
	  free((*qcc)->v_buffer);
	  (*qcc)->v_buffer = NULL;
	}
      
      if((*qcc)->buffer != NULL)
	{
	  free((*qcc)->buffer);
	  (*qcc)->buffer = NULL;
	}
      
      for(i=0; i< (*qcc)->nb_labels; ++i)
	{
	  if(((*qcc)->labels)[i].v != NULL)
	    {
	      free(((*qcc)->labels)[i].v);
	      ((*qcc)->labels)[i].v = NULL;
	    }
	}
      
      if((*qcc)->labels != NULL)
	{
	  free((*qcc)->labels);
	  (*qcc)->labels = NULL;
	}
      
      if((*qcc)->labels != NULL)
	{
	  free((*qcc)->perm_one);
	  (*qcc)->perm_one = NULL;
	}
      
      if((*qcc)->labels != NULL)
	{
	  free((*qcc)->perm_two);
	  (*qcc)->perm_two = NULL;
	}
      
      
      if((*qcc)->labels != NULL)
	{
	  free((*qcc)->perm_buffer);
	  (*qcc)->perm_buffer = NULL;
	}
      
      free(*qcc);
      *qcc = NULL;	
    } 
}

void print_quad_cover(quad_cover * qcc)
{
  interval * i;
  size_t j, n, d, k;
  int verbose = 0;
  
  i = qcc->top;
  while(i->prev != NULL)
    i = i->prev;
  while(i != NULL)
    {
      if(i->orientation == 1)
	printf(" ");
      else
	printf("-");
      printf("%zu ", i->lab - qcc->labels);
      i = i->next;
    }
  printf("\n");
  
  i = qcc->bot;
  while(i->prev != NULL)
    i = i->prev;
  while(i != NULL)
    {
      if(i->orientation == 1)
	printf(" ");
      else
	printf("-");
      printf("%zu ", i->lab - qcc->labels);
      i = i->next;
    }

  printf("\nlengths:");
  for(j=0; j < qcc->nb_labels; ++j){
    for(n=0; n < qcc->degree; ++n)
      printf(" %Lf", (qcc->labels)[j].length);
    printf("\n");
  }

  /* printf("\npermutation:"); */
  /* for(j=0; j < qcc->nb_labels; ++j){ */
  /*   for(n=0; n < qcc->degree; ++n) */
  /*     printf(" %d", (qcc->labels)[j].sigma[n]); */
  /*   printf("\n"); */
  /* } */

  printf("\n");
}


void rauzy_induction_H_plus_quad_cover(quad_cover *qcc)
/* perform one step of rauzy-zorich induction                           */
/* v shoud have dimension (nb_vectors) x (qcc->degree * qcc->nb_labels) */
{
  interval *win, *los, *win_twin;
  interval **los_ptr;
  size_t jlos,jwin;  /* index of los and win in the matrix */
  size_t i,k;
  double * tmp;
  int debug = 0, debug_ = 0;

  if(qcc->top->lab->length > qcc->bot->lab->length)
    {
      win = qcc->top;
      los = qcc->bot;
      los_ptr = &(qcc->bot);
    }
  else
    {
      win = qcc->bot;
      los = qcc->top;
      los_ptr = &(qcc->top);
    }

  /* Warning: los_ptr is intended to modify either qcc->top or qcc->bot */
  win_twin = win->twin;
  jlos = los->lab - qcc->labels;
  jwin = win->lab - qcc->labels;

  if (los->orientation == 1 && win->orientation == 1){
    permutation(0, qcc->perm_one, qcc->degree);
    inverse_permutation(win->lab->sigma, qcc->perm_two, qcc->degree);}
  if (los->orientation == 1 && win->orientation == -1){
    memcpy(qcc->perm_one, win->lab->sigma, qcc->degree * sizeof(size_t));
    memcpy(qcc->perm_two, win->lab->sigma, qcc->degree * sizeof(size_t));}
  if(los->orientation == -1 && win->orientation == 1){
    inverse_permutation(los->lab->sigma, qcc->perm_one, qcc->degree);
    permutation(0, qcc->perm_two, qcc->degree);}
  if(los->orientation == -1 && win->orientation == -1){
    inverse_permutation(los->lab->sigma, qcc->perm_buffer, qcc->degree);
    perm_product(win->lab->sigma, qcc->perm_buffer, qcc->perm_one, qcc->degree);
    permutation(0,qcc->perm_two, qcc->degree);}

  if (debug_){
    printf("win :%zu, %Lf      los :%zu, %Lf\ndiff : %Lf\n", win->lab - qcc->labels, win->lab->length, los->lab - qcc->labels, los->lab->length, win->lab->length - los->lab->length);
  }
  if (debug){
    printf("perm:\n");
    print_permutation(qcc->perm_one, qcc->degree);
    print_permutation(qcc->perm_two, qcc->degree);
  }

  /* modify the length */
  win->lab->length -= los->lab->length;
  qcc->length -= los->lab->length;

  if(win->lab->same_interval)   /* win and its twin are on the same interval */
    {
      /* move los after win_twin */
      *los_ptr = los->prev;
      los->prev->next = NULL;
      los->prev = win_twin->prev;
      los->next = win_twin;
      if(win_twin->prev != NULL) win_twin->prev->next = los;
      win_twin->prev = los;

      /* modify the position */
      los->lab->same_interval = 1 - los->lab->same_interval;
    }
  else
    {
      /* move los before win_twin */
      if(win->twin != los->next)
      {
      *los_ptr = los->prev;
      los->prev->next = NULL;
      los->next = win_twin->next;
      los->prev = win_twin;
      if(win_twin->next != NULL) win_twin->next->prev = los;
      else *los_ptr = los;
      win_twin->next = los;
      }
    }

  /* apply KZ to the vectors */
  /* we put the result of KZ in qcc->buffer and then we switch the pointers */
  for(i=0; i<qcc->nb_vectors; ++i)
    {
      for(k=0; k<qcc->degree; ++k)
	/* recall: 0 <= i < nb_vectors, 0 <= j < nb_labels, 0 <= k < degree */
	/* elt at pos (i,j,k) is (qcc->labels)[j].v[i + nb_vectors * k]     */
	{
	  (qcc->buffer)[i + qcc->nb_vectors * k] =
	    (qcc->labels)[jwin].v[i + qcc->nb_vectors * k]
	    + los->orientation * win->orientation * (qcc->labels)[jlos].v[i + qcc->nb_vectors * (qcc->perm_one[k])];
	}
      // switch qcc->labels[jlos] and qcc->buffer
    }
  tmp = qcc->buffer; qcc->buffer = (qcc->labels)[jwin].v; (qcc->labels)[jwin].v = tmp;

  for(i=0; i<qcc->nb_vectors; ++i)
    {
      for(k=0; k<qcc->degree; ++k)
	/* recall: 0 <= i < nb_vectors, 0 <= j < nb_labels, 0 <= k < degree */
	/* elt at pos (i,j,k) is (qcc->labels)[j].v[i + nb_vectors * k]     */
	{
	  (qcc->buffer)[i + qcc->nb_vectors * k] =
	    (qcc->labels)[jlos].v[i + qcc->nb_vectors * (qcc->perm_two[k])];
	}
    }
  // switch qcc->labels[jlos] and qcc->buffer
  tmp = qcc->buffer; qcc->buffer = (qcc->labels)[jlos].v; (qcc->labels)[jlos].v = tmp;


  /* modify the cover data                     */
  if (los->orientation == 1 && win->orientation == 1){
    inverse_permutation(win->lab->sigma, qcc->perm_one, qcc->degree);
    memcpy(qcc->perm_two, los->lab->sigma, qcc->degree * sizeof(size_t));
    perm_product(qcc->perm_one, qcc->perm_two, los->lab->sigma, qcc->degree);
  }
  if (los->orientation == 1 && win->orientation == -1){
    memcpy(qcc->perm_two, los->lab->sigma, qcc->degree * sizeof(size_t));
    perm_product(win->lab->sigma, qcc->perm_two, los->lab->sigma, qcc->degree);
  }
  if (los->orientation == -1 && win->orientation == 1){
    perm_product(los->lab->sigma, win->lab->sigma, los->lab->sigma, qcc->degree);
  }
  if (los->orientation == -1 && win->orientation == -1){
    inverse_permutation(win->lab->sigma, qcc->perm_two, qcc->degree);
    perm_product(los->lab->sigma, qcc->perm_two, los->lab->sigma, qcc->degree);
  }

  /*change begining of loser interval*/
  los = *los_ptr;
}


void top_lyapunov_exponents_H_plus(quad_cover *qcc, double *theta, size_t nb_iterations)
{
  size_t i,j,k,nb_ren=0;
  double norm;
#ifdef USE_KAHAN_SUMMATION
  double c0=0,c1=0;
  double y,t;
#endif

  set_random_lengths_quad_cover(qcc);
  set_random_vectors(qcc);
  
  theta[0] = theta[1] = 0;
  
  for(i=0; i < nb_iterations; ++i)
    {
      rauzy_induction_H_plus_quad_cover(qcc);
      
      if(qcc->length < 0.00001)
	{
#ifdef USE_KAHAN_SUMMATION
	  y = -logl(qcc->length) - c0;
	  t = theta[0] + y;
	  c0 = (t - theta[0]) - y;
	  theta[0] = t;
#else
	  theta[0] -= logl(qcc->length);
#endif
	  renormalize_length_quad_cover(qcc);
	  
	  norm = 0;
	  for(j=0; j<qcc->nb_labels; ++j)
	    for(k=0; k<qcc->degree; ++k)
	      {
		if((qcc->labels)[j].v[k] > 0)
		  norm += (qcc->labels)[j].v[k];
		else
		  norm -= (qcc->labels)[j].v[k];
	      }
	  for(j=0; j<qcc->nb_labels; ++j)
	    for(k=0; k<qcc->degree; ++k)
	      (qcc->labels)[j].v[k] /= norm;
#ifdef USE_KAHAN_SUMMATION
	  y = log(norm) - c1;
	  t = theta[1] + y;
	  c1 = (t - theta[1]) - y;
	  theta[1] = t;
#else
	  theta[1] += log(norm);
#endif
	  if((nb_ren+1) % 500 == 0)  // a bit of salt
	    {
	      qcc->top->lab->length += ((long double) (.5-drand())) / 0x40000000000;
	      qcc->bot->lab->length += ((long double) (.5-drand())) / 0x40000000000;
	    }
	  
	  nb_ren += 1;
	}
    }
  
  theta[1] /= theta[0];
  theta[0] /= nb_iterations;
}


void lyapunov_exponents_H_plus(quad_cover *qcc, double *theta, size_t nb_iterations)
{
  size_t i,nb_ren=0;
  
  set_random_lengths_quad_cover(qcc);
  set_random_vectors(qcc);
  
  for(i=0; i < qcc->nb_vectors+1; ++i) theta[i] = 0.;
  
  for(i=0; i < nb_iterations; ++i)
    {
      rauzy_induction_H_plus_quad_cover(qcc);

      if(qcc->length < 0.001)
	{
	  theta[0] -= logl(qcc->length);
	  renormalize_length_quad_cover(qcc);
	  orthogonalize_GS(qcc, theta);
	  
	  if((nb_ren+1) % 500 == 0)  // a bit of salt
	    {
	      qcc->top->lab->length += ((long double) (.5-drand())) / 0x40000000000;
	      qcc->bot->lab->length += ((long double) (.5-drand())) / 0x40000000000;
	    }
	  
	  nb_ren += 1;
	}
      theta[0] -= logl(qcc->length);
      renormalize_length_quad_cover(qcc);
      orthogonalize_GS(qcc, theta);
    }
  
  for(i=1; i<qcc->nb_vectors+1; i++) theta[i] /= (2*theta[0]);
  theta[0] /= nb_iterations;
}


void lyapunov_exponents_isotypic(quad_cover *qcc, double *theta, size_t nb_iterations, size_t nb_char, size_t* dimensions, double *proj)
{
  size_t i,j,nb_ren=0;
  double buffer;

  set_random_lengths_quad_cover(qcc);
  set_random_vectors(qcc);
  
  
  project_isotypic(qcc, nb_char, dimensions, proj);
  orthogonalize_iso(qcc, theta, nb_char, dimensions);		
  check_orthogonality_iso(qcc, nb_char, dimensions);
  
  for(i=0; i < qcc->nb_vectors+1; ++i) theta[i] = 0.;
  
  for(i=0; i < nb_iterations; ++i)
    {
      rauzy_induction_H_plus_quad_cover(qcc);

      if(qcc->length < 0.001)
	{
	  theta[0] -= logl(qcc->length);
	  renormalize_length_quad_cover(qcc);
	  
	  orthogonalize_iso(qcc, theta, nb_char, dimensions);
	  project_isotypic(qcc, nb_char, dimensions, proj);
	  
	  if((nb_ren+1) % 500 == 0)  // a bit of salt
	    {
	      qcc->top->lab->length += ((long double) (.5-drand())) / 0x40000000000;
	      qcc->bot->lab->length += ((long double) (.5-drand())) / 0x40000000000;
	    }
          
	  nb_ren += 1;
	}
      theta[0] -= logl(qcc->length);
      renormalize_length_quad_cover(qcc);
      orthogonalize_iso(qcc, theta, nb_char, dimensions);
      project_isotypic(qcc, nb_char, dimensions, proj);
    }

  for(i=1; i<qcc->nb_vectors+1; i++) theta[i] /= (2*theta[0]);
  theta[0] /= nb_iterations;
}
