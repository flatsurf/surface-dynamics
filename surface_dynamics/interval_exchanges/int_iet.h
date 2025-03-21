/*****************************************************************************/
/*       int_iet.h                                                           */
/*                                                                           */
/* Interval exchange transformations with integer lengths. The functions     */
/* works for both interval exchanges and linear involutions.                 */
/*                                                                           */
/*       Copyright (C) 2015 Vincent Delecroix <vincent.delecroix@labri.fr>   */
/*                                                                           */
/*  Distributed under the terms of the GNU General Public License (GPL)      */
/*  as published by the Free Software Foundation; either version 2 of        */
/*  the License, or (at your option) any later version.                      */
/*                  https://www.gnu.org/licenses/                            */
/*****************************************************************************/

#include <stdint.h>
#include <stdio.h>
 
typedef struct {
    uint64_t length;
    uint64_t height;
    int      same_interval;  // whether or not the two belong to the same interval
} label;

typedef struct Xinterval {
    struct Xinterval *prev;    // interval on the left
    struct Xinterval *next;    // interval on the right
    struct Xinterval *twin;    // the twin interval
    label * lab;               // associated label
} interval;

typedef struct {
    unsigned int nb_labels;  // number of labels (= pairs of intervals)
    interval * top;          // first interval on top
    interval * bot;          // first interval on bot

    /* internal stuff that should not be needed beyond creation/destruction */
    label * labels;          // the labels
    interval * intervals;    // the intervals (twice the number of labels)
} int_iet;

typedef int_iet int_iet_t[1];

typedef struct{
    /* public */
    uint64_t * x;
    int n;
    int kfree;
    int ktop;
    int kbot;

    /* internal */
    uint64_t nfree;
    uint64_t nother;
    uint64_t * xfree;
    uint64_t * xtop;
    uint64_t * xbot;
} li_vector_iterator;

typedef li_vector_iterator li_vector_iterator_t[1];

uint64_t uint64_rand();

/* memory allocation */
void int_iet_init(int_iet_t t, unsigned int n);
void int_iet_clear(int_iet_t t);

/* safety check */
/* Return 0 if t is valid and 1 otherwise */
int  int_iet_check(const int_iet_t t);

/* set data */
void int_iet_set_labels_and_twin(int_iet_t t, int * labels, int * twin, int k);
void int_iet_set_lengths(int_iet_t t, uint64_t * lengths);
void int_iet_set_random_lengths(int_iet_t t, uint64_t L);

/* output */
void int_iet_fprint(FILE * stream, int_iet_t t);
void int_iet_print(int_iet_t t);

/* number of cylinders */
int int_iet_num_cylinders(uint64_t * widths, uint64_t * heights, int_iet_t t);

/* iteration through integer vectors of given sum and length */
int int_vector_first(uint64_t * x, int n, int k);
int int_vector_next(uint64_t * x, int n, int k);
void int_li_vector_init(li_vector_iterator_t t, uint64_t n, int kfree, int ktop, int kbot);
void int_li_vector_clear(li_vector_iterator_t t);
int int_li_vector_prefirst(li_vector_iterator_t t);
int int_li_vector_first_or_next(li_vector_iterator_t t);
void int_li_vector_info(li_vector_iterator_t t);
