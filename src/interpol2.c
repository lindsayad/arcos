/** @file interpol2.c
 *  @brief New interpolation routines
 *
 * We want to obtain a function \f$f(r, z)\f$ from another one \f$phi(r, z)\f$
 * which is known only in some points of a stencil, 
 *     \f$\phi(r0 + i dr, z0 + j dz), 0 <= {i, j} < q\f$.
 *
 * We approximate \f$f\f$ as
 *
 *    \f$f(r, z) = \sum_{nm}  (r - r0)^n (z - z0)^m c_{nm}\f$
 *
 * where \f$m + n < p\f$, \f$p\f$ being the order of the interpolating function
 * plus one.
 *
 * The coefficients \f$c_{nm}\f$ are calculated from \f$\phi\f$ as
 *
 *   \f$c_{nm} = \sum_{ij} b_{nm ij} phi(r0 + i L, z0 + j L)\f$
 *
 * Thus, to define an interpolation method, we only need the matrix
 * \f$b_{nm ij}\f$, which has \f$1/2 q^2 p(p+1)\f$ elements.
 *
 * HOWEVER, using this general interpolation matrix turns out to be too
 * slow.  Hence we define also specific routines that have the same
 * interface and that are called transparently.  Thus the calling program
 * does not have to worry about the internals of how the interpolation
 * is computed.
 */
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>

#include "grid.h"
#include "interpol2.h"
#include "parameters.h"
#include "proto.h"
#include "rz_array.h"
#include "species.h"

#define anm(N_, M_) this->stencil[(N_) * this->method->q + (M_)]

static double gen_apply (interpol_t *this, double r, double z);
void bilin_set_coeffs (interpol_t *this);
void quadlog_set_coeffs (interpol_t *this);
double quadlog_apply (interpol_t *this, double r, double z);

/** @brief Creates a new interpolator that uses a given method. */
interpol_t*
interpol_new_a (double Lr, double Lz, interpol_method_t *method)
{
  interpol_t *interpol;

  interpol = xmalloc (sizeof (interpol_t));
  
  interpol->method = method;
  interpol->Lr = Lr;
  interpol->Lz = Lz;

  /* The number of pairs {i, j} (i, j >= 0) such that i + j < p 
     is p (p + 1) / 2 */
  interpol->coeffs = xmalloc (sizeof(double) * 
			      method->p * (method->p + 1) / 2);
  interpol->stencil = xmalloc (sizeof(double) * method->q * method->q);

  return interpol;
}

/** @brief Frees the interpolator */
void
interpol_free (interpol_t *this)
{
  free (this->coeffs);
  free (this->stencil);
  free (this);
}

/** @brief Sets the stencil of the interpolation object.
 *
 * @a r0 and @a z0 are the coordinates of the bottom-left corner of the stencil
 * (the [0,0] point).
*/
void
interpol_set_stencil (interpol_t *this, double r0, double z0, ...)
{
  int i, j, indx;
  va_list ap;

  va_start (ap, z0);

  for (indx = 0, i = 0; i < this->method->q; i++) {
    for (j = 0; j < this->method->q; j++) {
      this->stencil[indx++] = va_arg (ap, double);
    }
  }

  this->r0 = r0;
  this->z0 = z0;

  va_end (ap);

  /* If the method has an optimized routine, use it. */
  if (this->method->set_coeffs != NULL) {
    this->method->set_coeffs (this);
  }
  else {
    interpol_set_coeffs (this);
  }
}

/** @brief Sets the stencil reading from an array.
 *
 * @a ir, @a iz give the index of the cell in the array corresponding to
 * the stencil origin.
 */
void
interpol_set_stencil_at (grid_t *grid, interpol_t *this, double r0, double z0, 
			 rz_array_t *ar,
			 int ir, int iz, int itheta)
{
  int indx, i, j;
  REAL *sp, *ap;

  assert (ar->ntheta > itheta);

  sp = this->stencil;
  for (indx = 0, i = 0; i < this->method->q; i++) {
     double r;
     if (this->method->masses) {
       r = cyl_r_at (ir + i - this->method->ir0, grid->level);
     } else {
       r = 1;
     }

    ap = RZTP (ar, ir + i - this->method->ir0, 
               iz - this->method->iz0, itheta);

    for (j = 0; j < this->method->q; j++, ap += ar->strides[Z_INDX]) {
    /* The case ir = -1 does not give any problem, since we set
       up the correct boundaries. */
       *sp = r * (*ap);
       sp++;
    }
  }

  this->r0 = r0;
  this->z0 = z0;

  /* If the method has an optimized routine, use it. */
  if (this->method->set_coeffs != NULL) {
      this->method->set_coeffs (this);
  }
  else {
    interpol_set_coeffs (this);
  }
}

/** @brief Calculates the coefficients of the interpolating polynomial.
 *
 * Has to be called after interpol_set_stencil. */
void 
interpol_set_coeffs (interpol_t *this)
{
  double Lrn[MAX_ORDER], Lzn[MAX_ORDER];
  double *matrix, sum;
  int n, ip, j;
  int pindx, mindx;

  matrix = this->method->matrix;

  Lrn[0] = 1.0;
  Lzn[0] = 1.0;

  for (ip = 0, pindx = 0, mindx = 0; ip < this->method->p; ip++) {
    Lrn[ip + 1] = Lrn[ip] * this->Lr;
    Lzn[ip + 1] = Lzn[ip] * this->Lz;
    for (n = 0; n <= ip; n++) {
      sum = 0;

      /* This is the coefficient of r^n z^m, with m + n = ip */
      for (j = 0; j < this->method->q * this->method->q; j++) 
	sum += 	matrix[mindx++] * this->stencil[j];

      this->coeffs[pindx++] = sum / Lrn[n] / Lzn[ip - n];
    }
  }
}

/** @brief Evaluates the interpolating polynomial. */
double
interpol_apply (interpol_t *this, double r, double z)
{
  double result;

  /* If the method has an optimized routine, call it and return */
  if (this->method->apply != NULL) {
    result = this->method->apply (this, r, z);
  } else {
  /* Else, apply the general method. */
    result = gen_apply (this, r, z);
  }

  if (this->method->masses) {
    return result / cyl_q(r);
  } else {
    return result;
  }
}

/** @brief gen_apply ?????. */
static double
gen_apply (interpol_t *this, double r, double z)
{
  double deltar, deltaz, sum;

  deltar = r - this->r0;
  deltaz = z - this->z0;

  int ipmax=this->method->p;
  if (ipmax < 3) {
     fprintf(stdout,"ipmax=%d\n",ipmax);
     fatal("interpol2: order interpolation not equal to 3\n");
  } 

  if (deltar==deltaz) {
    sum = this->coeffs[0] +
          deltar * (this->coeffs[1] +
                    this->coeffs[2] +
                    deltar * (this->coeffs[3] +
                              this->coeffs[4] +
                              this->coeffs[5] ));

  } else
  {
    sum = this->coeffs[0] +
          deltaz * (this->coeffs[1] + deltaz * this->coeffs[3]) +
          deltar * (this->coeffs[2] + 
                    deltaz * this->coeffs[4] +
                    deltar * this->coeffs[5] );
  }
  return sum;
}   

/** You are not supposed to understand these numbers.  
 *  Refer to InterpolArrays.nb
 */
double matrix_zero[] = {1.0};

interpol_method_t interpol_zero = {matrix_zero, 1, 1, 0, 0, FALSE,
				   interpol_set_coeffs, 
				   gen_apply};

interpol_method_t interpol_zero_masses = {matrix_zero, 1, 1, 0, 0, TRUE,
					  interpol_set_coeffs, 
					  gen_apply};

/* You are not supposed to understand these numbers.  
   Refer to InterpolArrays.nb */

/** Coefficients for a 9-point stencil quadratic interpolation. */
double matrix_quadratic[] = 
  {-0.1111111111111111,  0.2222222222222222, -0.1111111111111111,
    0.2222222222222222,  0.5555555555555556,  0.2222222222222222,
   -0.1111111111111111,  0.2222222222222222, -0.1111111111111111,
   -0.1666666666666667,  0                 ,  0.1666666666666667,
   -0.1666666666666667,  0                 ,  0.1666666666666667,
   -0.1666666666666667,  0                 ,  0.1666666666666667,
   -0.1666666666666667, -0.1666666666666667, -0.1666666666666667,
    0                 ,  0                 ,  0                 ,
    0.1666666666666667,  0.1666666666666667,  0.1666666666666667,
    0.1666666666666667, -0.3333333333333333,  0.1666666666666667,
    0.1666666666666667, -0.3333333333333333,  0.1666666666666667,
    0.1666666666666667, -0.3333333333333333,  0.1666666666666667,
    0.25              ,  0                 , -0.25              ,
    0.25              ,  0                 , -0.25              ,
    0                 ,  0                 ,  0                 ,
   -0.25              ,  0                 ,  0.25              ,
    0.1666666666666667,  0.1666666666666667,  0.1666666666666667,
   -0.3333333333333333, -0.3333333333333333, -0.3333333333333333,
    0.1666666666666667,  0.1666666666666667,  0.1666666666666667
  };

interpol_method_t interpol_quadratic = {matrix_quadratic, 3, 3, 1, 1, FALSE,
					interpol_set_coeffs, 
					gen_apply};

interpol_method_t interpol_quadratic_masses = {matrix_quadratic, 3, 3, 1, 1, 
					       TRUE,
					       interpol_set_coeffs, 
					       gen_apply};

/** @brief In this approach, we make sure that the interpolation error at
 * the matrix center (anm(2, 2)) is zero (this is used by J. Wackers' code) */
double matrix_wackers[] = 
  { 0                 ,  0                 ,  0                 ,
    0                 ,  1.0               ,  0                 ,
    0                 ,  0                 ,  0                 ,
   -0.1666666666666667,  0                 ,  0.1666666666666667,
   -0.1666666666666667,  0                 ,  0.1666666666666667,
   -0.1666666666666667,  0                 ,  0.1666666666666667,
   -0.1666666666666667, -0.1666666666666667, -0.1666666666666667,
    0                 ,  0                 ,  0                 ,
    0.1666666666666667,  0.1666666666666667,  0.1666666666666667,
    0.1               , -0.2               ,  0.1               ,
    0.3               , -0.6               ,  0.3               ,
    0.1               , -0.2               ,  0.1               ,
    0.25              ,  0                 , -0.25              ,
    0                 ,  0                 ,  0                 ,
   -0.25              ,  0                 ,  0.25              ,
    0.1               ,  0.3               ,  0.1               ,
   -0.2               , -0.6               , -0.2               ,
    0.1               ,  0.3               ,  0.1               };

interpol_method_t interpol_wackers = {matrix_wackers, 3, 3, 1, 1, FALSE,
				      interpol_set_coeffs, 
				      gen_apply};

interpol_method_t interpol_wackers_masses = {matrix_wackers, 3, 3, 1, 1, TRUE,
					     interpol_set_coeffs, 
					     gen_apply};


double matrix_averaged[] = 
{ -0.01880341880341880 , -0.008547008547008547, -0.01880341880341880 ,
  -0.008547008547008547,  1.109401709401709   , -0.008547008547008547,
  -0.01880341880341880 , -0.008547008547008547, -0.01880341880341880 ,
  -0.1666666666666667  ,  0                   ,  0.1666666666666667  ,
  -0.1666666666666667  ,  0                   ,  0.1666666666666667  ,
  -0.1666666666666667  ,  0                   ,  0.1666666666666667  ,
  -0.1666666666666667  , -0.1666666666666667  , -0.1666666666666667  ,
   0                   ,  0                   ,  0                   ,
   0.1666666666666667  ,  0.1666666666666667  ,  0.1666666666666667  ,
   0.1128205128205128  , -0.1987179487179487  ,  0.1128205128205128  ,
   0.3012820512820513  , -0.6564102564102564  ,  0.3012820512820513  ,
   0.1128205128205128  , -0.1987179487179487  ,  0.1128205128205128  ,
   0.25                ,  0                   , -0.25                ,
   0                   ,  0                   ,  0                   ,
  -0.25                ,  0                   ,  0.25                ,
   0.1128205128205128  ,  0.3012820512820513  ,  0.1128205128205128  ,
  -0.1987179487179487  , -0.6564102564102564  , -0.1987179487179487  ,
   0.1128205128205128  ,  0.3012820512820513  ,  0.1128205128205128   };

interpol_method_t interpol_averaged = {matrix_averaged, 3, 3, 1, 1, FALSE,
				      interpol_set_coeffs, 
				      gen_apply};

interpol_method_t interpol_averaged_masses = {matrix_averaged, 3, 3, 1, 1, TRUE,
					     interpol_set_coeffs, 
					     gen_apply};

/* 4-point bilinear interpolation. */
double matrix_bilin[] =
  { 1.0,  0.0,  0.0,  0.0, 
   -1.0,  0.0,  1.0,  0.0,
   -1.0,  1.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,
    1.0, -1.0, -1.0,  1.0,
    0.0,  0.0,  0.0,  0.0};

/** @brief For bilinear interpolations, the standard fallback routines are too slow,
   so we define optimized routines to increase performance.  The results,
   however, should be the same. */
void 
bilin_set_coeffs (interpol_t *this)
{
  double Lr, Lz;
  double *c;

  c = this->coeffs;

  Lr = this->Lr;
  Lz = this->Lz;

  c[0] = anm(0, 0);
  c[1] = (anm(0, 1) - anm(0, 0)) / Lz;
  c[2] = (anm(1, 0) - anm(0, 0)) / Lr;
  c[3] = 0;  /* Never used anyway. */
  c[4] = (anm(0, 0) - anm(1, 0) - anm(0, 1) + anm(1, 1)) / Lr / Lz;
  c[5] = 0;  /* Never used anyway. */
}  

double
bilin_apply (interpol_t *this, double r, double z)
{
  double d_r, d_z;
  double *c;

  c = this->coeffs;

  d_r = r - this->r0;
  d_z = z - this->z0;

  return c[0] + d_z * (c[1] + c[4] * d_r) + c[2] * d_r; 
}


interpol_method_t interpol_bilin = {matrix_bilin, 
				    3, 2, 
				    0, 0, 
				    FALSE,
				    bilin_set_coeffs, 
				    bilin_apply};

interpol_method_t interpol_bilin_masses = {matrix_bilin, 
					   3, 2, 
					   0, 0, TRUE,
					   bilin_set_coeffs, 
					   bilin_apply};

void
quadlog_set_coeffs (interpol_t *this)
{
  int indx, i, j;

  for (indx = 0, i = 0; i < this->method->q; i++) {
    for (j = 0; j < this->method->q; j++) {
      this->stencil[indx] = log (this->stencil[indx]);
      indx++;
    }
  }
  interpol_set_coeffs (this);
}

double
quadlog_apply (interpol_t *this, double r, double z)
{
  double rlog;

  rlog = gen_apply (this, r, z);
  return exp (rlog);
}

interpol_method_t interpol_quadlog = {matrix_quadratic, 3, 3, 1, 1,
					FALSE,
				      quadlog_set_coeffs, 
				      quadlog_apply};

interpol_method_t *interpol_methods_index[MAX_INTERPOL_METHODS] = {
  &interpol_zero_masses,
  &interpol_quadratic_masses,
  &interpol_wackers_masses,
  &interpol_averaged_masses,
  &interpol_zero,
  &interpol_quadratic,
  &interpol_wackers,
  &interpol_averaged,
  &interpol_quadlog };
