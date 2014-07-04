/* Wrapper for the HSTCYL function. 
   Alejandro Luque Estepa - 2006
*/

#include "stdlib.h"
#include "stdio.h"
#include "../include/fishpack.h"
#include "math.h"
#include "assert.h"

const char *hstcyl_error_str[FISH_ERROR_MAX] = 
  {"No error",
   "A .LT. 0",
   "A .GE. B",
   "MBDCND .LT. 1 OR MBDCND .GT. 6",
   "C .GE. D",
   "N .LE. 2",
   "NBDCND .LT. 0 OR NBDCND .GT. 4",
   "A = 0 AND MBDCND = 1,2,3, OR 4",
   "A .GT. 0 AND MBDCND .GE. 5",
   "M .LE. 2",
   "IDIMF .LT. M",
   "LAMBDA .GT. 0",
   "A=0, MBDCND .GE. 5, ELMBDA .NE. 0"};

int number=0;

void mhstcyl_(double *r0, double *r1, int *nr, int *rbndcnd, 
	      double *bcr0, double *bcr1, 
	      double *z0, double *z1, int *nz, int *zbndcnd, 
	      double *bcz0, double *bcz1,
	      double *lambda, double *s, double *f, 
	      int *idimf, double *pertrb, int *ierror);

/* Calls hstcyl and checks for errors. If an error occurs, prints the
 corresponding message and exists. */
void 
fish_hstcyl (double r0, double r1, int nr, 
	     int rbndcnd, double *bcr0, double *bcr1,
	     double z0, double z1, int nz, 
	     int zbndcnd, double *bcz0, double *bcz1,
	     double lambda, double s, double *f, int idimf)
{
  double pertrb;
  int ierror;

  mhstcyl_(&r0, &r1, &nr, &rbndcnd, bcr0, bcr1, 
	   &z0, &z1, &nz, &zbndcnd, bcz0, bcz1,
	   &lambda, &s, f, &idimf, &pertrb, &ierror);
  number+=1;

/*  fprintf(stderr, "# calls mhstcyl       =%d\n",number); */
  if (pertrb != 0.0) {
    fprintf (stderr, "%s: ERROR: Undefined solution to the Poisson equation.",
	     __func__);
    fprintf (stderr, "\nTry to change your boundary conditions.\n");
    exit(-1);
  }

  if (0 != ierror) {
    if (ierror <= FISH_ERROR_MAX) {
      fprintf (stderr, "%s: ERROR: %s\n", __func__, hstcyl_error_str[ierror]);
      fprintf (stderr, "The call to fishpack was\n");
      fprintf (stderr, "mhstcyl_ (A = %f, B = %f, M = %d,\n",  r0, r1, nr);
      fprintf (stderr, "         MBDCND = %d, ...\n", rbndcnd);
      fprintf (stderr, "         C = %f, D = %f, N = %d,\n",  z0, z1, nz);
      fprintf (stderr, "         NBDCND = %d, ...\n", zbndcnd);
      fprintf (stderr, "         ELMBDA = %f, ...\n", lambda);
      fprintf (stderr, "         ES = %f, ...\n", s);
      fprintf (stderr, "         IDIMF = %d, ...)\n", idimf);
    } else {
      fprintf (stderr, "%s: Unknown error #%d\n", __func__, ierror);
    }

    exit(-1);
  }
}
  
