#include <stdlib.h>
#include <string.h>


/***********************************************/
/* memory allocation for normal_form functions */
int SF_realloc(unsigned int n);
void SF_free(void);

/***********************************************/
/* Comparison functions                        */
/*   -1 if o1 < o2                             */
/*    0 if o1 == o2                            */ 
/*    1 if o1 > o2                             */
/* TODO: inline troubles*/
int origami_diff(int *o1, int *o2, int n);
/* Comparison of two origamis. */
int pillowcase_cover_diff(int *g1, int *g2, int n);
/* Comparison of pillowcase covers. */

/***********************************************/
/* Normal form                                 */
/* Modify x and y in order to make the origami (x,y) in */
/* normal form                                          */
/* renum is set to the permutation used to renumerote   */
void origami_normal_form(int *x, int *y, int *renum, int n);
/* int pillowcase_cover_normal_form(int *g, int *renum, unsigned int n); */
