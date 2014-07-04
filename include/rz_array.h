/** @file rz_array.h
 *  @brief Definitions for 2D/3D arrays
 */

#ifndef _RZ_ARRAY_H_

#include "cstream.h"

typedef struct rz_array_t rz_array_t;

struct rz_array_t {
  REAL *data;
 
  int strides[3];
  int r0, z0, theta0, nr, nz;
  int ntheta;
  int dim;

  int len;

  rz_array_t *host;
};

#define R_INDX 0
#define Z_INDX 1
#define THETA_INDX 2

/* The 3D macros */
#define RZTP_OK(_A, _R, _Z, _T)						\
  ((_R) >= (_A)->r0 && (_R) < (_A)->nr + (_A)->r0 + 4			\
   && (_Z) >= (_A)->z0 && (_Z) < (_A)->nz + (_A)->z0 + 4		\
   && (_T) >= (_A)->theta0 && (_T) < (_A)->ntheta + (_A)->theta0 + 4)	\

#define __RZTP(_A, _R, _Z, _T) ((_A)->data				\
			      +	((_R) - (_A)->r0) * (_A)->strides[R_INDX] \
			      + ((_Z) - (_A)->z0) * (_A)->strides[Z_INDX] \
			      + (((_T) - (_A)->theta0) * (_A)->strides[THETA_INDX]))

#if defined (DEBUG_LEVEL) && DEBUG_LEVEL > 4
/* Look, mum! I am also doing array bound checking.
 * But, alas, since in some parts of the code we perform direct pointer
 * arithmetic, this check will not detect all possible out-of-indexes
 * (but we catch most of them).
 */
# define RZTP(_A, _R, _Z, _T) (RZTP_OK (_A, _R, _Z, _T)?		\
			       __RZTP (_A, _R, _Z, _T):			\
			       (fprintf (stderr, \
					 "%s:%d: Out of bounds ir = %d, iz = %d, itheta = %d\n", \
					 __FILE__, __LINE__, _R, _Z, _T), \
				fprintf (stderr, \
					 "->r0 = %d ->z0 = %d ->theta0 = %d ->nr = %d ->nz = %d ->ntheta = %d\n", \
					 (_A)->r0, (_A)->z0, (_A)->theta0, \
					 (_A)->nr, (_A)->nz, (_A)->ntheta), \
				exit(-1), (double*) NULL))
#else
# define RZTP(_A, _R, _Z, _T) __RZTP(_A, _R, _Z, _T)
#endif

#define RZT(_A,_R,_Z,_T) (*RZTP(_A, _R, _Z, _T))
#define RZTm(_A,_R,_Z,_T) ((fprintf (stderr, \
					 "%s:%d: CHECK ir = %d, iz = %d, itheta = %d\n", \
					 __FILE__, __LINE__, _R, _Z, _T), \
				fprintf (stderr, \
					 "->r0 = %d ->z0 = %d ->theta0 = %d ->nr = %d ->nz = %d ->ntheta = %d\n", \
					 (_A)->r0, (_A)->z0, (_A)->theta0, \
					 (_A)->nr, (_A)->nz, (_A)->ntheta), \
				fprintf (stderr, \
					 "->data = %g\n", \
					 (_A)->data), \
				exit(-1), (double*) NULL))


/* These are valid for 2D arrays.  When applied to a 3D array, they
 * return the value with theta=0.
 */
#define RZP(_A,_R,_Z) RZTP(_A, _R, _Z, (_A)->theta0)
#define RZ(_A,_R,_Z) RZT(_A, _R, _Z, (_A)->theta0)

#define BND_CND_HNEUMANN   1
#define BND_CND_HDIRICHLET -1

#define BND_INWARD -1
#define BND_OUTWARD 1

#define _RZ_ARRAY_H_
#endif
