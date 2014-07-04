#ifndef _FISPHACK_H_

/* Arguments for the boundary conditions in fishpack:
   The first name is the b.c. at the lowest value of the variable,
   the second one is the b.c. at the highest value.
     DIR = Dirichlet, 
     NEU = Neumann, 
     UNS = Unspecified (when r = 0 for cylindrical calculations).
*/
#define FISH_PERIODIC 0
#define FISH_DIR_DIR  1
#define FISH_DIR_NEU  2
#define FISH_NEU_NEU  3
#define FISH_NEU_DIR  4

#define FISH_UNS_DIR  5
#define FISH_UNS_NEU  6


/* This is the maximum number of grid cells that fishpack can handle
   without accumulating larger and larger roundoff errors
   (it's only an approximation: not checked throughly so if you notice
   that you are getting strange things such as artificial lines in
   the electric field, this is the first thing you should check). */
#define FISH_MAX_GRIDPOINTS   1700

#define LOG2 0.69314718055994530941723212145818

#define FISH_WORK(M_, N_)  (13 * (M_) + 4 * (N_)		\
			    + (M_) * (int) (log(N_) / LOG2))

#define FISH_ERROR_MAX 13
const char *hstcyl_error_str[FISH_ERROR_MAX];

void 
fish_hstcyl (double r0, double r1, int nr, 
	     int rbndcnd, double *bcr0, double *bcr1,
	     double z0, double z1, int nz, 
	     int zbndcnd, double *bcz0, double *bcz1,
	     double lambda, double s, double *f, int idimf);
void 
fish_hstcrt (double r0, double r1, int nr, 
	     int rbndcnd, double *bcr0, double *bcr1,
	     double z0, double z1, int nz, 
	     int zbndcnd, double *bcz0, double *bcz1,
	     double lambda, double *f, int idimf);

#define _FISHPACK_H_
#endif /* _FISHPACK_H_ */

