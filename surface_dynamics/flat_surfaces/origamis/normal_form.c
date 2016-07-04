#include <normal_form.h>
#include <stdio.h>

/***************/
/* COMPARISONS */
/***************/

inline int origami_diff(int *x1, int *x2, unsigned int n)
{
	unsigned int i;
	int test;

	for(i=0; i<2*n; ++i)
		if(test = x1[i] - x2[i]) return test;
	return 0;
}

inline int pillowcase_cover_diff(int *g1, int *g2, unsigned int n)
{
	unsigned int i;
	int test;

	for(i=0; i<4*n; ++i)
		if(test = g1[i] - g2[i]) return test;
	return 0;
}

/********************/
/* GLOBAL VARIABLES */
/********************/

unsigned int SF_NMAX = 0;
int * SF_bloc = NULL;
int * SF_g = NULL, *SF_gg = NULL, * SF_ggg = NULL;
int * SF_xx = NULL, * SF_xxx = NULL;
int * SF_yy = NULL, * SF_yyy = NULL;
int * SF_ren = NULL, * SF_ren_inv = NULL, * SF_rren = NULL;
int * SF_seen = NULL, * SF_to_test = NULL;
int * SF_cycles, * SF_lengths, *SF_vertices;

int SF_realloc(unsigned int n)
{
#ifdef DEBUG
	FILE * f = fopen("/home/vincent/alloc.log","a");
#endif

	if(n > SF_NMAX)
	{
#ifdef DEBUG
		fprintf(f,"ask for allocation of %d squares...",n);
#endif
		SF_NMAX = 2*n;
		if(SF_NMAX < 30) SF_NMAX = 30;
#ifdef DEBUG
		fprintf(f,"alloc %d squares\n",SF_NMAX);
		fprintf(f,"old values:\n SF_xx : %u\n SF_xxx: %u\n SF_yy : %u\n SF_yyy: %u\n",SF_xx,SF_xxx,SF_yy,SF_yyy);
		fflush(f);
#endif

		SF_bloc = (int *) realloc(SF_bloc, 16 * SF_NMAX * sizeof(int));

		SF_g        = SF_bloc +  0*SF_NMAX;
		SF_gg       = SF_bloc +  4*SF_NMAX;
        SF_ren      = SF_bloc +  8*SF_NMAX;
        SF_ren_inv  = SF_bloc +  9*SF_NMAX;
        SF_rren     = SF_bloc + 10*SF_NMAX;
        SF_seen     = SF_bloc + 11*SF_NMAX;
        SF_to_test  = SF_bloc + 12*SF_NMAX;
		SF_cycles   = SF_bloc + 13*SF_NMAX;
		SF_lengths  = SF_bloc + 14*SF_NMAX;
		SF_vertices = SF_bloc + 15*SF_NMAX;

		/* for origami we prefer use (x,y) instead of g */
		SF_xx       = SF_g  +  0*SF_NMAX;
		SF_yy       = SF_g  +  1*SF_NMAX;
		SF_xxx      = SF_g  +  2*SF_NMAX;
		SF_yyy      = SF_gg +  3*SF_NMAX;

#ifdef DEBUG
		fprintf(f,"new values:\n SF_xx : %u\n SF_xxx: %u\n SF_yy : %u\n SF_yyy: %u\n",SF_xx,SF_xxx,SF_yy,SF_yyy);
		fflush(f);
#endif

	}
    
    if (SF_bloc == NULL)
	{
		SF_free();
#ifdef DEBUG
		fprintf(f,"A pointer is null, free.\n");
		fflush(f);
		return 1;
#endif
	}
	return 0;
}
	
void SF_free(void)
{
	if(SF_bloc)
	{
		free(SF_bloc);
		SF_bloc=NULL;
	}

	SF_NMAX=0;
}

/* TODO: add in argument the list of potential 0 squares */
int origami_normal_form(int *x, int *y, int *renum, unsigned int n)
{
	int i,j,k,cmax;
	int mi,m;
	int i_test = 0;
	int * tmp;

    SF_realloc(n);

	/* set SF_xxx = MAX and SF_yyy = MAX */
	/* set SF_vertices = x y x^-1 y ^-1  */
	/* set SF_seen = 1                   */
    for(i=0; i<n; ++i)
	{
		SF_xxx[i] = SF_yyy[i] = n;

		SF_vertices[y[x[i]]] = x[y[i]];
		SF_seen[i] = 1;
	}
       
	/* compute the length of the maximal cycle */
	/* set SF_cycles[i] = (cycle number of i)  */
	/* set SF_lengths[i] = (cycle length of i) */
	cmax=0;
	j=0;
	for(i=0; i<n; ++i)
	{
		if(SF_seen[i])
		{
			SF_lengths[j] = 0;
			do
			{
				SF_seen[i] = 0;
				SF_cycles[i] = j;
				++(SF_lengths[j]);
				i = SF_vertices[i];
			}while(SF_seen[i]);
			if(SF_lengths[j] > cmax)
				cmax = SF_lengths[j];
			++j;
		}
	}

	/* set SF_to_test[0..i_test] = elements of SF_vertices with maximal cycles */
    for(i=0; i<n; ++i)
        if(SF_lengths[SF_cycles[i]] == cmax)
	 		SF_to_test[i_test++] = i;

	/* test the renumerotation for all element of SF_to_test */
	/* at each step:                                         */
	/* SF_seen[i] = 0 if and only if i is renumeroted        */
    while(i_test)
	{
		for(i=0; i<n; SF_seen[i++] = 1);
		i_test--;
		i = SF_to_test[i_test];

		/* renumerote the first cycle */
		k = 0;
		do{
			SF_seen[i] = 0;
			SF_ren_inv[SF_ren[i] = k] = i;
			k++;
			i = x[i];
		}while(SF_seen[i]);

		/* renumerote the other cycles */
		for(mi=0; mi < n; ++mi)
		{
			m = y[SF_ren_inv[mi]];
			while(SF_seen[m])
			{
				do{
					SF_seen[m] = 0;
					SF_ren_inv[SF_ren[m] = k] = m;
					k++;
					m = x[m];
				}while(SF_seen[m]);

				m = y[m];
			}
		}

		/* set SF_xx = "x renumeroted" */
		/* set SF_yy = "y renumeroted" */
		for(i=0; i<n; ++i)
		{
			SF_xx[SF_ren[i]] = SF_ren[x[i]];
			SF_yy[SF_ren[i]] = SF_ren[y[i]];
		}

		/* compare (SF_xx,SF_yy) to (SF_xxx,SF_yyy) and */
		/* set to the minimum                           */
		for(i=0; i<n-1; ++i)
		{
			if(SF_xx[i] < SF_xxx[i])
			{
				tmp = SF_xx; SF_xx = SF_xxx; SF_xxx = tmp;
				tmp = SF_yy; SF_yy = SF_yyy; SF_yyy = tmp;
				tmp = SF_ren; SF_ren = SF_rren; SF_rren = tmp;
				break;
			}
			else if(SF_xx[i] > SF_xxx[i]) break;
		}
		if(i == n-1) /* here SF_xx = SF_xxx */
			for(i=0; i < n-1; ++i)
			{
				if(SF_yy[i] < SF_yyy[i])
				{
					tmp = SF_xx; SF_xx = SF_xxx; SF_xxx = tmp;
					tmp = SF_yy; SF_yy = SF_yyy; SF_yyy = tmp;
					tmp = SF_ren; SF_ren = SF_rren; SF_rren = tmp;
					break;
				}
				else if(SF_yy[i] > SF_yyy[i]) break;
			}
	}
	
	/* set x = SF_xxx */
	/* set y = SF_yyy */
	/* set renum = SF_rren */
	for(i=0; i < n; ++i)
	{
		x[i] = SF_xxx[i];
		y[i] = SF_yyy[i];
		renum[i] = SF_rren[i];
	}
}
	

int projective_normal_form(int *x, int *y, int *renum, unsigned int n)
{
	int * xx = (int *) malloc(n*sizeof(int));
	int * yy = (int *) malloc(n*sizeof(int));
	int * rrenum = (int *) malloc(n*sizeof(int));
	int i,test;

    SF_realloc(n);

	for(i=0; i<n; ++i)
	{
		xx[x[i]] = i;
		yy[y[i]] = i;
	}

	origami_normal_form(x,y,renum,n);
	origami_normal_form(xx,yy,rrenum,n);

	test = 1;
	for(i=0; i<n-1; ++i)
		if(x[i] < xx[i])
		{
			test=0;
			break;
		}
	if(test)
		for(i=0; i<n-1; ++i)
			if(y[i] < yy[i])
			{
				test=0;
				break;
			}
	if(test)
	{
		memcpy(x,xx,n*sizeof(int));
		memcpy(y,yy,n*sizeof(int));
		memcpy(renum,rrenum,n*sizeof(int));
	}

	free(xx);
	free(yy);
	free(rrenum);
}

