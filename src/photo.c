/** @file photo.c
 *  @brief Routines for the calculation of the photoionization source term.
 *
 * We use the Helmholtz approximation for the kernel in the photoionization
 * integral and we use the routines of the Poisson solver for the Helmholtz
 * equation (FISHPACK was modified to take a new inhomogeneous term into
 * account).
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "cdr.h"
#include "interpol2.h"
#include "mapper.h"
#include "parameters.h"
#include "photo.h"
#include "poisson.h"
#include "proto.h"
#include "rz_array.h"
#include "species.h"

void photo_copy (mapper_t *mapper, grid_t *source, grid_t *target, 
		 int ir, int iz, int itheta);
void photo_coarsen (mapper_t *mapper, grid_t *source, grid_t *target, 
		    int ir, int iz, int itheta);
int photo_interpol_set (mapper_t *mapper, grid_t *source, 
			 interpol_t *interpol,
			 int pr, int pz, int itheta);
void photo_interpol (mapper_t *mapper, grid_t *source, grid_t *target, 
		     interpol_t *interpol, 
		     int ir, int iz, int itheta);

photo_term_t *photo_terms = NULL;

mapper_t photo_mapper = {&interpol_bilin, photo_coarsen, photo_copy,
			 photo_interpol_set, photo_interpol};

mapper_t *photo_mappers[] = {&photo_mapper, NULL};

pois_problem_t *pois_photo_1;
pois_problem_t *pois_photo_2;
extern pois_problem_t *pois_electrostatic;

/** @brief Initializes photoionization QQQQ? */
void
photo_init ()
{
  /* Since it makes no sense, we use photo_extra_levels < 0 to say that
     the conditions for the photoionization are the same as for
     the electrostatic problem.  This is needed for backwards compatibility. */
  if (extra_photo_levels < 0) 
  {
    pois_photo_1 = pois_electrostatic;
    pois_photo_2 = pois_electrostatic;
    return;
  }

  /* Two different sets of refinement criteria for the two photoionization
     terms. If extra_photo_levels_2 < 0, then pois_photo_2 = pois_photo_1 */
  pois_photo_1 = xmalloc (sizeof (pois_problem_t));

  pois_photo_1->max_level = photo_max_level;
  pois_photo_1->extra_levels = extra_photo_levels;
  pois_photo_1->max_error = photo_max_error;

  pois_photo_1->bnd_right = photo_bnd_right;
  pois_photo_1->bnd_top = photo_bnd_top;
  pois_photo_1->bnd_bottom = photo_bnd_bottom;

  if (extra_photo_levels_2 < 0)
  {
     pois_photo_2 = pois_photo_1;
     return;
  }

  pois_photo_2 = xmalloc (sizeof (pois_problem_t));

  pois_photo_2->max_level = photo_max_level_2;
  pois_photo_2->extra_levels = extra_photo_levels_2;
  pois_photo_2->max_error = photo_max_error_2;

  pois_photo_2->bnd_right = photo_bnd_right_2;
  pois_photo_2->bnd_top = photo_bnd_top_2;
  pois_photo_2->bnd_bottom = photo_bnd_bottom_2;
}


/** @brief Copies the photoionization QQQQ? */
void
photo_copy (mapper_t *mapper, grid_t *source, grid_t *target, 
	    int ir, int iz, int itheta)
{
  cdr_grid_t *cdr;
  pois_grid_t *pois;

  cdr = (cdr_grid_t*) target;
  pois = (pois_grid_t*) source;

  RZT (cdr->photo, ir, iz, itheta) = RZ (pois->phi, ir, iz);
}

/** @brief Coarsens the photoionization QQQQ? */
void 
photo_coarsen (mapper_t *mapper, grid_t *source, grid_t *target, 
	       int ir, int iz, int itheta)
{
  cdr_grid_t *cdr;
  pois_grid_t *pois;
  int level_diff, z, r;

  cdr = (cdr_grid_t*) target;
  pois = (pois_grid_t*) source;

  level_diff = pois->level - cdr->level;

  z = (iz << level_diff) + (1 << (level_diff - 1));
  r = (ir << level_diff) + (1 << (level_diff - 1));

  if (grid_contains (source, r, z, GRID_INSIDE) && 
      grid_contains (source, r - 1, z, GRID_INSIDE) &&
      grid_contains (source, r, z - 1, GRID_INSIDE) &&
      grid_contains (source, r - 1, z - 1, GRID_INSIDE)) {
    
    RZT (cdr->photo, ir, iz, itheta)  = 0.25 *
      (RZ (pois->phi, r, z) 
       + RZ (pois->phi, r - 1, z) 
       + RZ (pois->phi, r, z - 1) 
       + RZ (pois->phi, r - 1, z - 1));
  }

}

/** @brief photo_interpol_set QQQQ? */
int 
photo_interpol_set (mapper_t *mapper, grid_t *source, interpol_t *interpol,
		      int pr, int pz, int itheta)
{
  pois_grid_t *pois;

  pois = (pois_grid_t*) source;

  interpol_set_stencil (interpol,
			r_at (pr, pois->level),
			z_at (pz, pois->level),
			RZ (pois->phi, pr, pz), 
			RZ (pois->phi, pr, pz + 1),
			RZ (pois->phi, pr + 1, pz), 
			RZ (pois->phi, pr + 1, pz + 1));

  return TRUE;
}

/** @brief photo_interpol QQQQ? */
void
photo_interpol (mapper_t *mapper, grid_t *source, grid_t *target, 
		interpol_t *interpol, int ir, int iz, int itheta)
{
  double r, z;
  cdr_grid_t *cdr;

  cdr = (cdr_grid_t *) target;

  r = r_at (ir, cdr->level);
  z = z_at (iz, cdr->level);

  RZT (cdr->photo, ir, iz, itheta) = interpol_apply (interpol, r, z);
}


/** @brief Registers a photoionization term with given @a A and @a lambda. */
void
photo_register (double A, double lambda)
{
  photo_term_t *t;

  t = (photo_term_t *) xmalloc (sizeof (photo_term_t));

  t->A = A;
  t->lambda = lambda;
  t->next = photo_terms;

  photo_terms = t;
}

/** @brief Copies a list of photo terms into @a *dest. */
void
photo_copy_list (photo_term_t *src, photo_term_t **dest)
{
  photo_term_t *ptr_term;
  photo_term_t *p = NULL, *p1;

  *dest = NULL;

  while (src) {
    p1 = (photo_term_t *) xmalloc (sizeof (photo_term_t));
    if (p) {
      p->next = p1;
    } else {
      *dest = p1;
    }

    p = p1;
    p->A = src->A;
    p->lambda = src->lambda;
    src = src->next;
  }

  if (p) p->next = NULL;
}


/** @brief Unregisters all the photoionization terms and
 *  frees the allocated space.
 */
void
photo_unregister_all (void)
{
  photo_term_t *t, *next;

  for (t = photo_terms; t; t = next) {
    next = t->next;
    free (t);
  }

  photo_terms = NULL;
}

/** @brief Transforms back the photoionization calculation into real space. */
void
photo_dft_r (cdr_grid_t *grid, int sign)
{
  cdr_grid_t *leaf;

  debug (3, "photo_dft_r(" grid_printf_str ", %d)\n", 
	 grid_printf_args(grid), sign);

  iter_childs (grid, leaf) {
    photo_dft_r (leaf, sign);
  }

  dft_transform (grid->photo, grid->photo, sign);
}

/** @brief Copies the derivative of the ion density into cdr->charge
 *
 * The ion density, which is supposed to be at this point, the impact
 * ionization) into cdr->charge, which will be used as the source for
 * the Poisson/Helmholtz solver. */
void
photo_copy_source (cdr_grid_t *grid)
{
  int ir, iz, itheta;

   debug(3, "photo_copy_source (" grid_printf_str ")\n",
	 grid_printf_args(grid));

#pragma omp parallel
  {
#pragma omp for private (ir, iz)
    iter_grid_3d_n (grid, ir, iz, itheta, 2) {
      RZT (grid->charge, ir, iz, itheta) = 
	RZT (grid->d_dens[ions], ir, iz, itheta) / grid->ntheta;
    }
  }
}

/** @brief Recursive version of @a photo_copy_source. */
mk_recursive (photo_copy_source, cdr_grid_t)

/** @brief Once a photoionization term is computed, we add it to d_dens.
 *
 * Note that we already copied the contents of d_dens[ions] into
 * charge, so when we solve again the Helholtz equation, the source
 * is still the same.
 */
void
photo_add_term (photo_term_t *term, cdr_grid_t *cdr)
{
  int s, ir, iz, itheta;
  int updated[2] = {electrons, photo_ions};
  
  debug (3, "photo_add_term (" photo_printf_str ", " grid_printf_str ")\n",
	 photo_printf_args(term), grid_printf_args(cdr));
    
#pragma omp parallel
    {
#pragma omp for private (ir, iz, s)
      iter_grid_3d_n (cdr, ir, iz, itheta, 2) {
	for (s = 0; s < 2; s++) {
	  RZT (cdr->d_dens[updated[s]], ir, iz, itheta) += 
	    term->A * RZT (cdr->photo, ir, iz, itheta);
	}
      }
    }  
}

/** @brief ...and the recursive version of @a photo_add_term */
void
photo_add_term_r (photo_term_t *term, cdr_grid_t *cdr)
{
  cdr_grid_t *child;

  photo_add_term (term, cdr);

  iter_childs (cdr, child) {
    photo_add_term_r (term, child);
  }
}

/** @brief photo_calc_term QQQQ */
pois_grid_t **
photo_calc_term (photo_term_t *term, cdr_grid_t *cdr, int i)
{
  /* Call to the generic Poisson/Helmholtz solver */
  if (i == 1) {
	return pois_gen_solve_a (cdr, pois_photo_2, photo_mappers, term->lambda);
  } else {
	return pois_gen_solve_a (cdr, pois_photo_1, photo_mappers, term->lambda);
  }
}

/** @brief Calculates the photoionization and adds it to the derivatives
 *  of the species densities. */
void
photo_calc (photo_term_t *terms, cdr_grid_t *cdr)
{
  pois_grid_t **pois_modes;
  photo_term_t *term;
  // photo_term_t *ptr_term;

  photo_copy_source_r (cdr);
  
  if (cdr->ntheta != 1)
    cdr_dft_charge_r (cdr, 1);
  
  // for (ptr_term = terms; ptr_term != NULL; ptr_term = ptr_term->next) {
  //   printf("terms: A=%g lambda=%g ptr=%d\n", ptr_term->A, ptr_term->lambda, (int)ptr_term->next);
  // }

  int j = 0;
  for (term = terms; term; term = term->next ) {
    int i;
    pois_modes = photo_calc_term (term, cdr, j);
    j++;
    if (cdr->ntheta != 1)
      photo_dft_r (cdr, -1);

    photo_add_term_r (term, cdr);
    debug (3, "photo_calc (" photo_printf_str ", " grid_printf_str ")\n",
        	 photo_printf_args(term), grid_printf_args(cdr));

    /* Free the allocated memory. */
    for (i = 0; i < max_ntheta; i++) {      
      pois_free_r (pois_modes[i]);
    }
    free (pois_modes);
  }
}

/** @brief Loads a photoionization file, consisting in a series of rows with
   @a A and @a lambda. */
void
photo_load_file (char *fname)
{
  FILE *fp;
  double A, lambda;
  int c;
  int i;

  printf("\n");
  printf ("Loading photoionization data from `%s'...\n", fname);
  printf("\n");
  fp = fopen (fname, "r");

  if (NULL == fp) {
    warning ("Unable to open photoionization file `%s'\n", fname);
    exit (-1);
    return;
  }

  
  i=0;
  do {
    c = fscanf (fp, "%lf %lf", &A, &lambda);
    if (c != 2) {
      break;
    }

    printf ("Registring photoionization term A = %.5g, lambda = %.4g\n",
	    A, lambda);

    photo_register (A, lambda);
    i++;
  } while (TRUE);
  printf("\n");

  fclose (fp);
}
