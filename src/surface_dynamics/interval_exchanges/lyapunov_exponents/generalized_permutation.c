#include "lyapunov_exponents.h"

generalized_permutation * new_generalized_permutation(int *p, int *t, int k, int n)
{
	generalized_permutation * gp = (generalized_permutation *) malloc(sizeof(generalized_permutation));
	gp->perm = malloc(4 * n * sizeof(int));
	gp->twin = gp->perm + 2 * n;
	memcpy(gp->perm, p, 2 * n * sizeof(int));
	memcpy(gp->twin, t, 2 * n * sizeof(int));
	gp->k = k;
	gp->n = n;
	return gp;
}

void free_generalized_permutation(generalized_permutation ** gp)
{
	if((*gp) != NULL)
	{
		if((*gp)->perm != NULL)
		{
			//fprintf(stdout, "FREE gp at %p with perm at %p\n",*gp, (*gp)->perm);
			free((*gp)->perm);
			(*gp)->perm = NULL;
			(*gp)->twin = NULL;
		}
		free(*gp);
		*gp = NULL;
	}
}

void print_generalized_permutation(generalized_permutation *gp)
{
	size_t i;
	for(i=0; i < gp->k; ++i) printf(" %2d", (gp->perm)[i]);
	printf("\n");
	for(i=0; i < gp->k; ++i) printf(" %2d", (gp->twin)[i]);
	printf("\n**********\n");
	for(i=gp->k; i < 2*gp->n; ++i) printf(" %2d", (gp->perm)[i]);
	printf("\n");
	for(i=gp->k; i < 2*gp->n; ++i) printf(" %2d", (gp->twin)[i]);
	printf("\n");
}

int check_generalized_permutation(generalized_permutation * gp)
{
	int i;
	int * multiplicities;

	multiplicities = (int *) malloc(gp->n * sizeof(int));
	if(multiplicities == NULL)
	{
		fprintf(stderr, "allocation error inside check_generalized_permutation");
		return 2;
	}
	for(i=0; i < gp->n; ++i) multiplicities[i] = 0;

	for(i=0; i < 2*gp->n; ++i)
	{
		if(((gp->perm)[i] < 0) || ((gp->perm)[i] > gp->n))
		{
			fprintf(stderr, "the range of p is wrong: p[%d] = %d\n", i, (gp->perm)[i]);
			return 1;
		}
		if(multiplicities[(gp->perm)[i]] == 2)
		{
			fprintf(stderr, "%d appears more than twice in the permutation\n", (gp->perm)[i]);
			return 1;
		}
		multiplicities[(gp->perm)[i]] += 1;

		if((i == (gp->twin)[i]) || ((gp->twin)[(gp->twin)[i]] != i))
		{
			fprintf(stderr,"Error: t is not an involution without fixed point i=%d -> t[i]=%d -> t[t[i]]=%d \n",
							i,(gp->twin)[i],(gp->twin)[(gp->twin)[i]]);
			return 1;
		}
		if((gp->perm)[i] != (gp->perm)[(gp->twin)[i]])
		{
			fprintf(stderr,"Error: p does not coincide with t at i=%d and t[i]=%d (p[i]=%d and p[t[i]]=%d)\n",
							i,(gp->twin)[i],(gp->perm)[i],(gp->perm)[(gp->twin)[i]]);
			return 1;
		}
	}
	free(multiplicities);
	return 0;
}
