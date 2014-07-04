/** @file dft.c
 * @brief All routines related with discrete fourier transformations.
 *
 *  We make intensive use of the fftw library.\n
 *  FFTW, the Fastest Fourier Transform in the West, is a collection of fast
 *  C routines to compute the discrete Fourier transform.
 *  See http://www.fftw.org/doc/ for the manual documents of FFTW version 3.3.3.
 *
 *  Note: All calls to FFTW should be made here, so 
 *  1.  If we decide to change our FFT library, we only have to make 
 *      changes in this module.
 *  2.  We can compile a 2d version of the code that does not depend on fftw.
 *      by eliminating this module.
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <fftw3.h> 

#include "cdr.h"
#include "cstream.h"
#include "grid.h"
#include "parameters.h"
#include "proto.h"
#include "rz_array.h"
#include "species.h"

static void renormalize (cdr_grid_t *grid, rz_array_t *var);
static double rnd_gauss (double mu, double sigma);
static double ranf (void);

/** @brief Transform the discrete fourier transformations ???? */
void
dft_transform (rz_array_t *in, rz_array_t *out, int sign)
{
  int i;
  fftw_plan *plans;
  fftw_r2r_kind kind;

  debug (2, "dft_transform (..., sign = %d)\n", sign);

  if (sign > 0) kind = FFTW_R2HC;
  else kind = FFTW_HC2R;

  assert(in->dim == 3 && out->dim == 3);
  assert(in->ntheta == out->ntheta && in->nr == out->nr && in->nz == out->nz);

  /* Unfortunately, the fftw planning routines are not thread-safe. */

  /* The +4 here is because we also transform the buffer, which is
     useful to easily set boundaries, etc. */
  plans = xmalloc (sizeof (fftw_plan) * (in->nz + 4));

  for (i = 0; i < in->nz + 4; i++) {
    plans[i] = fftw_plan_many_r2r (1,                            /* Rank */
				   &in->ntheta,                  /* n[] */
				   in->nr + 4,                   /* howmany */
				   RZTP (in, in->r0, in->z0 + i, 0), 
				   /* in */
				   &in->ntheta,                  /* inembed */
				   in->strides[THETA_INDX],      /* istride */
				   in->strides[R_INDX],          /* idist */
				   RZTP (out, in->r0, out->z0 + i, 0), 
				   /* out */
				   &out->ntheta,                 /* onembed */
				   out->strides[THETA_INDX],     /* ostride */
				   out->strides[R_INDX],         /* odist */
				   &kind,                        /* kind */
				   FFTW_ESTIMATE);               /* flags */
    
  }

#pragma omp parallel
  {
#pragma omp for
    for (i = 0; i < in->nz + 4; i++) {
      fftw_execute (plans[i]);
    }
  }


  for (i = 0; i < in->nz + 4; i++) {
      fftw_destroy_plan (plans[i]);
  }

  free (plans);
  fftw_cleanup();
}


/** @brief Calculates the Fourier transform of the derivative of a function given its
   Fourier transform. */
void
dft_diff (grid_t *grid, rz_array_t *f)
{
  int ir, iz, k;

  debug (2, "dft_diff (" grid_printf_str ", ...)\n", grid_printf_args(grid));

  /* The parallelization here is not optimal if max_ntheta is not a multiple
     of 2. */
#pragma omp parallel
  {
#pragma omp for private(ir, iz)
    for (k = 0; k < grid->ntheta / 2 + 1; k++) {
      if (k < grid->ntheta / 2) {
	iter_grid_n (grid, ir, iz, 2) {
	  double re, im, re_f, im_f;
	
	  /* For k == 0, everything is real. */
	  re_f = RZT (f, ir, iz, k);
	  im_f = (k == 0? 0: RZT (f, ir, iz, grid->ntheta - k));
	  
	  re = re_f * wk[k] - (k == 0? 0: im_f * wk[grid->ntheta - k]);
	  im = (k == 0? 0: re_f * wk[grid->ntheta - k] + im_f * wk[k]);
	  
	  RZT (f, ir, iz, k) = re / r_at (ir, grid->level);
	  if (k != 0) 
	    RZT (f, ir, iz, grid->ntheta - k) = im / r_at (ir, grid->level);
	} 
      } else if ((grid->ntheta % 2) == 0) {
	iter_grid_n (grid, ir, iz, 2) {
	  /* Also for k = ntheta / 2, everything is real. */
	  RZT (f, ir, iz, k) = wk[k] * RZT (f, ir, iz, k) 
	    / r_at (ir, grid->level);
	}
      }   
    }
  }
}


/** @brief Calculates the @a weight for a given cdr grid and a variable ???
 *
 * For a given cdr grid and a variable, calculates its "weight": 
 * the integral of the variable raised to power for each Fourier mode.
 * weights[] must be able to contain at least cdr->ntheta values.
*/
void
dft_weight (cdr_grid_t *cdr, rz_array_t *var, double weights[], double power)
{
  int k, ir, iz;
  double pwr, tmp;

  debug (2, "dft_weigth(" grid_printf_str ", ...)\n", grid_printf_args(cdr));

#pragma omp parallel
  {
#pragma omp for private(ir, iz, pwr, tmp)
    for (k = 0; k < cdr->ntheta; k++) {
      pwr = 0;
      iter_grid(cdr, ir, iz) {
	/* Slow, since power is usually 1 or 2, 
	   but this is outside the main loop, so we buy generality
	   at the cost of speed. */
	tmp = pow(RZT(var, ir, iz, k), power);
	pwr += (tmp * r_at(ir, cdr->level));
      }
      weights[k] = twopi * pwr * dr[cdr->level] * dz[cdr->level];
    }
  }
}


/** @brief Outputs to a file 'weights.tsv' the results of dft_weight.
 *
 * Assumes that grid->charge is calculated but in real space then it
 * transforms it to Fourier space and leaves it like that: so do not
 * assume that grid->charge is unchanged.
 *
 *  grid->dens[electrons], on the other hand, is transformed forth and back.
 */
void
dft_out_weights (cdr_grid_t *grid, const char *prefix, double t)
{
  static FILE *fp_charge = NULL;
  static FILE *fp_electrons = NULL;
  static double *powers_charge = NULL;
  static double *powers_electrons = NULL;
  int i;

  if (NULL == fp_charge) {
    char *fname;
    asprintf (&fname, "%s/charge-weights.tsv", prefix);
    fp_charge = fopen (fname, "w");

    if (NULL == fp_charge) {
      fatal ("Unable to open %s\n", fname);
      return;
    }
    free (fname);

    asprintf (&fname, "%s/electrons-weights.tsv", prefix);
    fp_electrons = fopen (fname, "w");

    if (NULL == fp_electrons) {
      fatal ("Unable to open %s\n", fname);
      return;
    }
    free (fname);
  }

  if (NULL == powers_charge) {
    powers_charge = xmalloc (sizeof(double) * grid->ntheta);
    powers_electrons = xmalloc (sizeof(double) * grid->ntheta);
  }

  dft_transform (grid->charge, grid->charge, 1);  
  dft_weight (grid, grid->charge, powers_charge, 2.0);

  /* The electron density has to be transformed back, since we will still need
     it. */
  dft_transform (grid->dens[electrons], grid->dens[electrons], 1);  
  dft_weight (grid, grid->dens[electrons], powers_electrons, 1.0);
  dft_transform (grid->dens[electrons], grid->dens[electrons], -1);
  renormalize(grid, grid->dens[electrons]);

  fprintf (fp_charge, "%g", t);
  fprintf (fp_electrons, "%g", t);

  for (i = 0; i < grid->ntheta; i++) {
    fprintf (fp_charge, "\t%g", powers_charge[i]);
    fprintf (fp_electrons, "\t%g", powers_electrons[i]);
  }

  fprintf (fp_charge, "\n");
  fprintf (fp_electrons, "\n");

  fflush(fp_charge);
  fflush(fp_electrons);
}


/** @brief Perturbs a FFT-transformed variable.
 *
 * Assumes that the unperturbed variable is axi-symmetrical and hence
 * all modes are zero except k=0.  The perturbation is then selected
 * as a Gaussian random number epsilon_k times the zero-mode value of
 * the variable.  epsilon_k is distributed as
 * epsilon_k = N(0, perturb_epsilon).
 */
void
dft_perturb (cdr_grid_t *cdr, rz_array_t *var, double *epsilon_k)
{
  int k, ir, iz;

#pragma omp parallel
  {
#pragma omp for private(ir, iz)
    for (k = 1; k < cdr->ntheta; k++) {
      if (k > perturb_max_k && k < (cdr->ntheta - perturb_max_k))
	continue;

      iter_grid_n (cdr, ir, iz, 2) {
	RZT (var, ir, iz, k) = 
	  /* We multiply the perturbation by r to ensure that it is
	     continuous at r -> 0. */
	  epsilon_k[k - 1] * RZT (var, ir, iz, 0) * r_at (ir, cdr->level);
      }
    }
  }
}

/** @brief Takes a density in real space, transforms it to Fourier space,
 * perturbs it and then transforms back to real space. */
void 
dft_dens_perturb_r (cdr_grid_t *grid, int species, double *epsilon_k)
{
  int ir, iz, k, allocated = FALSE;
  cdr_grid_t *leaf;
  double A;

  debug (2, "dft_dens_perturb_r(" grid_printf_str ", %s)\n", 
	 grid_printf_args(grid), spec_index[species]->name);

  /* We normalize the perturbation using the maximum of the densities. */
  A =  invpi32 * seed_N / (seed_sigma_x * seed_sigma_y * seed_sigma_z);

  if (NULL == epsilon_k) {
    epsilon_k = (double*) xmalloc (sizeof(double) * (grid->ntheta - 1));
    
    for (k = 1; k < grid->ntheta; k++) {
      /* We divide by ntheta because the FFT is not normalized and we want
	 to express perturb_epsilon as a fraction of the k=0 mode. */
      epsilon_k[k - 1] = rnd_gauss (0, perturb_epsilon) / grid->ntheta / A;
    }
    
    allocated = TRUE;
  }

  iter_childs (grid, leaf) {
    dft_dens_perturb_r (leaf, species, epsilon_k);
  }

  dft_transform (grid->dens[species], grid->dens[species], 1);
  
  renormalize(grid, grid->dens[species]);

  dft_perturb (grid, grid->dens[species], epsilon_k);
  dft_transform (grid->dens[species], grid->dens[species], -1);

  if (allocated) free (epsilon_k);
}


/** @brief Renormalize @a var for some cases.
 *
 * FFTW does not normalize the Fourier transform, so after transforming
 * forth and back, the original values are multiplied by max_ntheta.
 * We use this routine to renormalize @a var in those cases.
 *
 * Note that this is taken care of when calculating the charge, so
 * in most cases you do not need to call this function.
 */
static void 
renormalize (cdr_grid_t *grid, rz_array_t *var)
{
  int ir, iz, itheta;

#pragma omp parallel
  {
#pragma omp for private(ir, iz)
    for (itheta = 0; itheta < grid->ntheta; itheta++) {
      iter_grid_n (grid, ir, iz, 2) {
	RZT(var, ir, iz, itheta) /= grid->ntheta;
      }
    }
  }
}


/** @brief Returns a random number with a gaussian distribution centered
 * around mu with width sigma.
 */
static double
rnd_gauss(double mu, double sigma)
{
  double x1, x2, w;
  static double y1, y2;
  static int has_more = FALSE;

  
  if (has_more) {
    has_more = FALSE;
    return mu + y2 * sigma;
  }

  do {
    x1 = 2.0 * ranf() - 1.0;
    x2 = 2.0 * ranf() - 1.0;
    w = x1 * x1 + x2 * x2;
  } while (w >= 1.0);
  
  w = sqrt ((-2.0 * log (w)) / w);
  y1 = x1 * w;
  y2 = x2 * w;

  has_more = TRUE;

  return mu + y1 * sigma;
}

#define AM (1.0 / RAND_MAX)

/** @brief Returns a random number uniformly distributed in [0, 1] */
static double
ranf (void)
{
  return rand() * AM;
}
