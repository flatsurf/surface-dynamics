#include "int_iet.h"

#include <stdlib.h>

int int_vector_first(uint64_t * x, int n, int k)
{
    int i;
    if ((k < 0) || (n < 0) || (k > n))
        return 0;
    if (k == 0) return n == 0;
    x[0] = n - k + 1;
    for (i = 1; i < k; i++) x[i] = 1;
    return 1;
}

int int_vector_next(uint64_t * x, int n, int k)
{
    int i;
    if (k == 1 || k == n) return 0;
    if (x[0] > 1)
    {
        x[0] -= 1;
        x[1] += 1;
        return 1;
    }
    for (i = 1; x[i] == 1; i++);
    if (i == k - 1)
        return 0;
    x[i+1] += 1;
    x[0] += x[i] - 2;
    x[i] = 1;
    return 1;
}

#include <stdio.h>
void int_li_vector_info(li_vector_iterator_t t)
{
    int i;
    printf("n      = %d\n", t->n);
    printf("kfree  = %d\n", t->kfree);
    printf("ktop   = %d\n", t->ktop);
    printf("kbot   = %d\n", t->kbot);
    printf("nfree  = %lu\n", t->nfree);
    printf("nother = %lu\n", t->nother);
    printf("xfree  =");
    for(i=0; i<t->kfree; i++)
        printf( " %lu", t->xfree[i]);
    printf("\n");
    printf("xtop   =");
    for(i=0; i<t->ktop; i++)
        printf(" %lu", t->xtop[i]);
    printf("\n");
    printf("xbot   =");
    for(i=0; i<t->kbot; i++)
        printf(" %lu", t->xbot[i]);
    printf("\n");
}

void int_li_vector_init(li_vector_iterator_t t, uint64_t n, int kfree, int ktop, int kbot)
{
    t->x = (uint64_t *) malloc((kfree + ktop + kbot) * sizeof(uint64_t));
    t->kfree = kfree;
    t->ktop = ktop;
    t->kbot = kbot;

    t->n = n;

    t->xfree = t->x;
    t->xtop = t->xfree + kfree;
    t->xbot = t->xtop + ktop;
}

void int_li_vector_clear(li_vector_iterator_t t)
{
    free(t->x);
}

int int_li_vector_prefirst(li_vector_iterator_t t)
{
    t->nfree = t->nother = -1;
    return 1;
}

/* we want that "free + 2*top" = "free + 2*bot" = n */
int int_li_vector_first_or_next(li_vector_iterator_t t)
{
    if (t->nfree == -1 && t->nother == -1)
    {
        if (t->kfree == 0)
            t->nfree = 0;
        else if (t->ktop > t->kbot)
            t->nfree = t->n - 2 * t->ktop;
        else
            t->nfree = t->n - 2 * t->kbot;
        if (t->nfree % 2 != t->n % 2)
            t->nfree -= 1;
        t->nother = (t->n - t->nfree) / 2;

        return int_vector_first(t->xfree, t->nfree, t->kfree) && \
               int_vector_first(t->xtop, t->nother, t->ktop) && \
               int_vector_first(t->xbot, t->nother, t->kbot);
    }

    if (int_vector_next(t->xfree, t->nfree, t->kfree)) return 1;
    int_vector_first(t->xfree, t->nfree, t->kfree);
    if (int_vector_next(t->xtop, t->nother, t->ktop)) return 1;
    int_vector_first(t->xtop, t->nother, t->ktop);
    if (int_vector_next(t->xbot, t->nother, t->kbot)) return 1;
    /* decrease nfree and increase nother */
    t->nfree -= 2;
    t->nother += 1;
    return int_vector_first(t->xfree, t->nfree, t->kfree) && \
           int_vector_first(t->xtop, t->nother, t->ktop) && \
           int_vector_first(t->xbot, t->nother, t->kbot);
}
