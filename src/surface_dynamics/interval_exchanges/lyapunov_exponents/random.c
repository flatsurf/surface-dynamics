#include "lyapunov_exponents.h"

/* random (?) uniform double in [0,1] */
inline double drand(void)
{
  return (((double) random() + (double) random() / ( (double) RAND_MAX)) / ((double) RAND_MAX+1) );
}

/* random (?) uniform long double in [0,1] */
inline long double ldrand(void)
{
  return ((long double) random() + (((long double) random() + (long double) random() / ((long double) RAND_MAX)) / ((long double) RAND_MAX))) / ((long double) RAND_MAX+1);
}
