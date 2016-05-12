#include "lyapunov_exponents.h"

void print_permutation(size_t *sigma, size_t degree)
{
  size_t i;
  for(i=0; i < degree; ++i)
    printf(" %zu", sigma[i] + 1);
  printf("\n");
}

int check_permutation(size_t *sigma, size_t degree)
{
  size_t i, j;
  size_t seen;

  for(i=0; i < degree; ++i)
    {
      seen = 0;
      for(j=0; j < degree; ++j)
	{
	  if(sigma[j] == i)
	    {
	      if(seen)
		{
		  fprintf(stderr, "%zu appears twice in the permutation", i);
		  return 1;
		}
	    }
	  else
	    seen = 1;
	  ++j;
	}
      if(!seen){
	fprintf(stderr, "%zu doesn't appear in the permutation", i);
	return 1;
      }
    }
  for(i=0; i < degree; ++i)
    {
      if(sigma[i] >= degree || sigma[i] < 0)
	{
	  fprintf(stderr, "%zu is not between 0 and %zu", i, degree);
	  return 1;
	}
    }
  return 0;
}
	

inline void inverse_permutation(size_t *sigma, size_t *perm_buffer, size_t degree)
{
  size_t i;
  for(i=0; i < degree; ++i)
    perm_buffer[sigma[i]] = i;
}

inline void permutation(int n, size_t *perm_buffer, size_t degree)
{
  size_t i;
  for(i=0; i < degree; ++i)
    perm_buffer[i] = (i + n) % degree;
}

inline void perm_name(interval *inter, size_t *perm_buffer, size_t degree)
{
  if(inter->orientation == 1)
    permutation(0, perm_buffer, degree);        //return identity
  else
    inverse_permutation(inter->lab->sigma, perm_buffer, degree);
}

inline void perm_ident_rev(interval *inter, size_t *perm_buffer, size_t degree)
{
  if(inter->orientation == 1)
    memcpy(perm_buffer, inter->lab->sigma, degree * sizeof(size_t));
  else
    inverse_permutation(inter->lab->sigma, perm_buffer, degree);
}

inline void perm_product(size_t *sigma, size_t *tau, size_t *perm_buffer, size_t degree)
{
  size_t i;
  for (i=0; i < degree; ++i)
    perm_buffer[i] = tau[sigma[i]];
}
