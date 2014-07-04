/** @file cdr.c
*   @brief Routines for handling the CDR (convection-diffusion-reaction) equations
*
* Default interpolation methods.\n
* May be changed elsewhere.\n
* \n
* interpol_method_t *cdr_interpol_inside = \&interpol_quadratic;\n
* interpol_method_t *cdr_interpol_bnd = \&interpol_bilin;\n
* \n
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#include "cdr.h"
#include "cstream.h"
#include "grid.h"
#include "interpol2.h"
#include "mapper.h"
#include "parameters.h"
#include "poisson.h"
#include "proto.h"
#include "react_table.h"
#include "rz_array.h"
#include "species.h"
#include "tree.h"

/** Declaration of the functions to interpolate the densities. */
decl_mapper_funcs(dens);
mapper_t **dens_mappers, **dens_bnd_mappers;

int cdr_ext_bnd_cond[4];
double z_cutoff = 10000000.0;
double max_err = 0.0;
react_table mu;
react_table diff;

extern pois_problem_t *pois_electrostatic;

static void set_all_bnd (cdr_grid_t *grid,
			 int sign, int ir, int iz, int itheta,
			 int dim0, int inout,
			 int dim1, int dim1_from, int dim1_to,
			 int dim2, int dim2_from, int dim2_to);
static void set_axis_bnd (cdr_grid_t *grid);
static void max_update (double *x, double cmp);
static double f_ad (rz_array_t *dens, double efield, double dc, double mc,
		    int ir, int iz, int itheta, int dim);

static double psi (double theta);
static void prepare_grid (cdr_grid_t *grid);
static double curv_at (cdr_grid_t *grid, rz_array_t *ar,
		       int ir, int iz, int itheta, double (*f) (double));
static int needs_refinement (cdr_grid_t *grid, int ir, int iz, int itheta,
			     int *in_edge);
static int any_needs_refinement (cdr_grid_t *grid, int ir, int iz,
				 int *in_edge);
static int brick_needs_refinement (cdr_grid_t *grid, int r0, int z0,
				   int r1, int z1, int *in_edge);
static void refine_in (cdr_grid_t *grid, int cr0, int cz0, int cr1, int cz1,
		       int contains_edge);

static void restrict_from (cdr_grid_t *parent, cdr_grid_t *child);
static int match_grids (cdr_grid_t *fro, cdr_grid_t *to);

double gauss2_xyz (double x, double y, double z);
int curr_seed = 0;
react_table mu;
react_table diff;

void aux_dump_frames_r (cdr_grid_t *grid, FILE *fp);

/** @brief Some initialization is needed for this module which has to be
 * inside this function.
 */
void
cdr_init (void)
{
  cdr_ext_bnd_cond[BND_BOTTOM] = cdr_bnd_bottom;
  cdr_ext_bnd_cond[BND_TOP] = cdr_bnd_top;
  cdr_ext_bnd_cond[BND_RIGHT] = cdr_bnd_right;

  dens_mappers = cdr_mappers_a (interpol_methods_index[cdr_interp_in]);
  dens_bnd_mappers = cdr_mappers_a (interpol_methods_index[cdr_interp_bnd]);

  react_table_read("input/default.mobility",&mu);
  react_table_read("input/default.diffusion",&diff);
}

/** @brief And the cleaning up has to be here. */
void
cdr_end (void)
{
  cdr_free_mappers (dens_mappers);
  cdr_free_mappers (dens_bnd_mappers);
}

/** @brief Creates a new 3D cdr grid */
cdr_grid_t*
cdr_new_3d_a (int r0, int z0, int r1, int z1, int ntheta)
{
  rz_array_t *ez, *er, *etheta, *eabs, *charge, **dens, **d_dens, *photo;
  cdr_grid_t *grid;
  int i;
  REAL *max_dens;

  debug (2, "cdr_new_3d_a (%d, %d, %d, %d)\n", r0, z0, r1, z1);

  ez = rz_new_3d_a (r0, z0, r1, z1, ntheta);
  er = rz_new_3d_a (r0, z0, r1, z1, ntheta);
  etheta = rz_new_3d_a (r0, z0, r1, z1, ntheta);
  eabs = rz_new_3d_a (r0, z0, r1, z1, ntheta);
  charge = rz_new_3d_a (r0, z0, r1, z1, ntheta);
  photo = has_photoionization? rz_new_3d_a (r0, z0, r1, z1, ntheta): NULL;

  /* Apart from the real species, (electron, ions...) we can create
   * a virtual species that enter into the possible reactions.
   * For example, the absolute value of the electric field, eabs.
   */
  dens = (rz_array_t **) xmalloc (sizeof(rz_array_t*)
				  * (no_species + N_VIRTUAL_SPECIES));
  d_dens = (rz_array_t **) xmalloc (sizeof(rz_array_t*) * no_species);

  max_dens = (REAL *) xmalloc (sizeof(REAL) * no_species);

  for (i = 0; i < no_species; i++) {
    dens[i] = rz_new_3d_a (r0, z0, r1, z1, ntheta);
    d_dens[i] = rz_new_3d_a (r0, z0, r1, z1, ntheta);
  }

  dens[no_species] = eabs;

  grid = (cdr_grid_t *) xmalloc (sizeof(cdr_grid_t));

  grid->r0 = r0;
  grid->r1 = r1;
  grid->z0 = z0;
  grid->z1 = z1;
  grid->charge = charge;
  grid->er = er;
  grid->ez = ez;
  grid->etheta = etheta;
  grid->eabs = eabs;
  grid->photo = photo;

  grid->dens = dens;
  grid->d_dens = d_dens;
  grid->max_dens = max_dens;

  grid->ntheta = ntheta;

  grid->contains_edge = FALSE;

  init_leaf (grid);

  return grid;
}

/** @brief Creates a guest cdr grid:
 *
 *  one that does not have memory of its own but
 *  points to the data of another grid (the @a host grid).  This is useful
 *  to easily set boundary conditions.  Probably, this function and
 *  cdr_new_3d_a should be somehow merged.
 */
cdr_grid_t*
cdr_guest (cdr_grid_t *host, int r0, int z0, int r1, int z1)
{
  rz_array_t *ez, *er, *etheta, *eabs, *charge, **dens, **d_dens, *photo;
  cdr_grid_t *grid;
  int i;
  REAL *max_dens;

  debug (2, "cdr_guest (" grid_printf_str ",%d, %d, %d, %d)\n",
	 grid_printf_args(host), r0, z0, r1, z1);

  ez = rz_guest (host->ez, r0, z0, r1, z1);
  er = rz_guest (host->er, r0, z0, r1, z1);
  etheta = rz_guest (host->etheta, r0, z0, r1, z1);
  eabs = rz_guest (host->eabs, r0, z0, r1, z1);
  charge = rz_guest (host->charge, r0, z0, r1, z1);
  photo = (has_photoionization?
	   rz_guest (host->photo, r0, z0, r1, z1): NULL);

  dens = (rz_array_t **) xmalloc (sizeof(rz_array_t*)
				  * (no_species + N_VIRTUAL_SPECIES));
  d_dens = (rz_array_t **) xmalloc (sizeof(rz_array_t*) * no_species);

  max_dens = (REAL *) xmalloc (sizeof(REAL) * no_species);

  for (i = 0; i < no_species; i++) {
    dens[i] = rz_guest (host->dens[i], r0, z0, r1, z1);
    d_dens[i] = rz_guest (host->d_dens[i], r0, z0, r1, z1);
  }

  dens[no_species] = eabs;

  grid = (cdr_grid_t *) xmalloc (sizeof(cdr_grid_t));

  grid->r0 = r0;
  grid->r1 = r1;
  grid->z0 = z0;
  grid->z1 = z1;

  grid->charge = charge;
  grid->er = er;
  grid->ez = ez;
  grid->etheta = etheta;
  grid->eabs = eabs;
  grid->photo = photo;

  grid->dens = dens;
  grid->d_dens = d_dens;
  grid->max_dens = max_dens;

  grid->ntheta = host->ntheta;

  grid->contains_edge = FALSE;

  set_leaf (grid, host->parent, NULL, NULL, host->level);

  return grid;
}

/** @brief Returns a grid with the same dimensions as grid and at the same
 * level.
 *
 * Realize that the data is _not_ copied into the new grid.
 */
cdr_grid_t*
cdr_like_a (cdr_grid_t *grid)
{
  cdr_grid_t *n;

  debug (2, "cdr_like_a(" grid_printf_str ")\n",
	 grid_printf_args(grid));

  n = cdr_new_3d_a (grid->r0, grid->z0, grid->r1, grid->z1, grid->ntheta);
  n->level = grid->level;
  n->ext_bound = grid->ext_bound;
  n->contains_edge = grid->contains_edge;

  return n;
}

/** @brief Returns a grid exactly like te one received, with the densities
 * copied.
 *
 * But the descendants are ignored.
 */
cdr_grid_t*
cdr_clone_a (cdr_grid_t *grid)
{
  cdr_grid_t *n;
  int itheta, s, nr, nz, nspecs;

  n = cdr_like_a (grid);

  nr = grid->r1 - grid->r0 + 4;
  nz = grid->z1 - grid->z0 + 4;

  nspecs = no_species;

  /* When we use the efield as a refinement criterium, we also have to copy
   * eabs (which is a @i virtual density).
   */
  if (ref_level_eabs >= 0 && ref_threshold_eabs >= 0.0) nspecs++;

#pragma omp parallel private(s)
  {
#pragma omp for
    iter_grid_theta_n (grid, itheta, 2) {
      for (s = 0; s < nspecs; s++) {
	rz_copy_modes (grid->dens[s], grid->r0 - 2, grid->z0 - 2,
		       n->dens[s], n->r0 - 2, n->z0 - 2,
		       nr, nz, itheta, itheta);
      }
    }
  }
  return n;
}

/** @brief Impose periodic boundary conditions in \f$\theta\f$ on a grid. */
void
cdr_set_periodic (cdr_grid_t *grid)
{
  int i;

  rz_set_periodic (grid->etheta);

  for (i = 0; i < no_species; i++) {
    /* For negative r, we have the values
     * of the densities on an opposite itheta.
     */

    rz_set_periodic (grid->dens[i]);
    rz_set_periodic (grid->d_dens[i]);

  }

  /* I am not completely sure that this is actually neccesary */
  rz_set_periodic (grid->er);
  rz_set_periodic (grid->ez);
}

/** @brief Recursive version of cdr_set_periodic. */
mk_recursive (cdr_set_periodic, cdr_grid_t)

/** @brief Frees a cdr grid. */
void
cdr_free (cdr_grid_t *grid)
{
  int i;

  debug (2, "cdr_free\n");

  rz_free (grid->charge);
  rz_free (grid->er);
  rz_free (grid->ez);
  rz_free (grid->etheta);
  rz_free (grid->eabs);
  if (has_photoionization) rz_free (grid->photo);

  for (i = 0; i < no_species; i++) {
    rz_free (grid->dens[i]);
    rz_free (grid->d_dens[i]);
  }

  free (grid->dens);
  free (grid->d_dens);
  free (grid->max_dens);

  free (grid);
}

/** @brief Recursive version of cdr_free. */
void
cdr_free_r (cdr_grid_t *grid)
{
  cdr_grid_t *leaf;

  debug (2, "cdr_free_r (" grid_printf_str ")\n", grid_printf_args (grid));

  free_childs (grid, leaf, cdr_free_r);

  cdr_free (grid);
}

/** @brief Calculates the charge in a CDR grid.
 *
 *  Important Note:  The charge is here divided by ntheta because
 *  our FFT functions are unnormalized and hence when they are applied
 *  forth and back they return the initial value multiplied
 *  by ntheta.  I believe that this is the best point to normalize the
 *  data.
 */
void
cdr_calc_charge (cdr_grid_t *grid)
{
  int s, r, z, t;
  double s_charge;

  debug (3, "cdr_calc_charge\n");

  rz_set_zero (grid->charge);

  for (s = 0; s < no_species; s++) {
    s_charge = spec_index[s]->charge;

#pragma omp parallel private(r, z)
    {
#pragma omp for
      iter_grid_theta (grid, t) {
	for (r = grid->r0 - 2; r < grid->r1 + 2; r++) {
	  for (z = grid->z0 - 2; z < grid->z1 + 2; z++) {
	    RZT (grid->charge, r, z, t) += (RZT (grid->dens[s], r, z, t)
					    * s_charge) / grid->ntheta;
	  }
	}
      }
    }
  }
}

/** @brief Recursive version of cdr_calc_charge. */
mk_recursive (cdr_calc_charge, cdr_grid_t)

/** @brief Fourier-transforms the charge in grid and all its descendants. */
void
cdr_dft_charge_r (cdr_grid_t *grid, int sign)
{
  cdr_grid_t *leaf;

  debug (2, "cdr_dft_charge_r(" grid_printf_str ")\n",
	 grid_printf_args(grid));

  iter_childs (grid, leaf) {
    cdr_dft_charge_r (leaf, sign);
  }

  dft_transform (grid->charge, grid->charge, sign);
}

/** @brief Fourier-transforms the electric field in grid and all its descendants. */
void
cdr_dft_field_r (cdr_grid_t *grid, int sign)
{
  cdr_grid_t *leaf;

  debug (2, "cdr_dft_field_r(" grid_printf_str ")\n",
	 grid_printf_args(grid));

  iter_childs (grid, leaf) {
    cdr_dft_field_r (leaf, sign);
  }

  dft_transform (grid->er, grid->er, sign);
  dft_transform (grid->ez, grid->ez, sign);

  dft_diff ((grid_t *) grid, grid->etheta);
  dft_transform (grid->etheta, grid->etheta, sign);
}

/** @brief Creates a grid one level coarser than GRID that covers the same area
 *    or somewhat more.
 *
 *  Note:  The new grid contains information only about the charge,
 *  not the other functions.
 *  That's why we don't use cdr_restrict.
 */
cdr_grid_t*
cdr_create_coarser_a (cdr_grid_t *grid)
{
  int r0_half, r1_half, z0_half, z1_half, itheta;
  int r, z;

  cdr_grid_t *new_grid;

  debug (2, "cdr_create_coarser_a\n");

  r0_half = grid->r0 >> 1;
  r1_half = (grid->r1 + 1) >> 1;

  z0_half = grid->z0 >> 1;
  z1_half = (grid->z1 + 1) >> 1;

  new_grid = cdr_new_3d_a (r0_half, z0_half, r1_half, z1_half, grid->ntheta);

#pragma omp parallel private (r, z)
  {
    /** Note that here we may be going one cell out of the aparent boundaries
       of the input grids but we allocated enough memory for that. */
#pragma omp for
    iter_grid_3d_n (new_grid, r, z, itheta, 1) {
      /* Start cylindrical */
      /* Here we are interpolating the masses, but in many other places
       * we interpolate the densities. Why?
       */
      RZT (new_grid->charge, r, z, itheta) =
	0.25
	* (cyl_q (r + 0.25) * (RZT (grid->charge, 2 * r, 2 * z, itheta)
			      + RZT (grid->charge, 2 * r, 2 * z + 1, itheta))
	   + cyl_q (r + 0.75) * (RZT (grid->charge, 2 * r + 1, 2 * z, itheta)
				+ RZT (grid->charge, 2 * r + 1, 2 * z + 1,
				       itheta)))
	/ cyl_q (r + 0.5);
      /* End cylindrical */
    }
  }
  return new_grid;
}

/** @brief Takes a tree of input nodes and adds n coarser grids, returning a
 *  pointer to the root leaf.
 */
cdr_grid_t*
cdr_add_coarser_grids_a (cdr_grid_t *prev_root, int n)
{
  int i;
  cdr_grid_t *grid = NULL;

  debug (2, "cdr_add_coarser_grids_a\n");

  assert (n > 0);

  for (i = 0; i < n; i++) {
    grid = cdr_create_coarser_a (prev_root);
    set_leaf (grid, NULL, NULL, prev_root, prev_root->level - 1);

    prev_root->parent = grid;

    prev_root = grid;
  }
  return grid;
}

/** @brief Frees the memory allocated by cdr_add_coarser_grids_a.
 *
 * Note: be careful to pass here only a tree returned by
 * cdr_add_coarser_grids_a, since some grids may lose their references and
 * stay allocated forever otherwise.
 * (see assert below).
 */
void
cdr_free_coarser_grids (cdr_grid_t *prev_root, int n)
{
  int i;
  cdr_grid_t *grid;

  debug (2, "cdr_free_coarser_grids_a\n");

  for (i = 0; i < n; i++) {
    assert (NULL == prev_root->next);
    grid = prev_root->first_child;
    cdr_free (prev_root);
    prev_root = grid;
  }

  prev_root->parent = NULL;
}

/** @brief The interface between the cdr and the poisson parts
 *
 *  If return_pois is true, returns the Poisson trees and they have
 *  to be de-allocated by the calling program.  If not, takes care himself
 *  of deallocation and returns NULL.
 */
pois_grid_t**
cdr_calc_field_r (cdr_grid_t *cdr, int return_pois)
{
  pois_grid_t **pois_modes;

  debug (2, "cdr_calc_field_r (" grid_printf_str ", %d)\n",
	 grid_printf_args(cdr), return_pois);

  if (cdr->ntheta != 1)
    cdr_dft_charge_r (cdr, 1);

  pois_modes = pois_solve_a (cdr, pois_electrostatic);

  if (cdr->ntheta != 1)
    cdr_dft_field_r (cdr, -1);

  cdr_add_ext_field_r (cdr);

  /* Formerly, I added the inhomogeneous field here.  But this is problematic:
   *  see pois_add_inhom_phi_r in poisson.c for details.
   *
   *  To return back to adding electric fields and not potentials, uncomment
   *  these lines:

   *  if (pois_inhom) {
   *     q_factor = pois_inhom_q_factor (pois_modes[0]);
   *     debug (1, "q_factor = %g\n", q_factor);
   *
   *     cdr_add_inhom_field_r (cdr, q_factor);
   *  }
   */

  /* Calculates the absolute value of the field. */
  cdr_calc_eabs_r (cdr);

  /* Impose periodic boundary conditions in theta. */
  if (cdr->ntheta != 1)
    cdr_set_periodic_r (cdr);

  if (!return_pois) {
    int i;

    for (i = 0; i < max_ntheta; i++) {
      pois_free_r (pois_modes[i]);
    }
    free (pois_modes);
    pois_modes = NULL;
  }
  return pois_modes;
}

/** @brief Adds an external electric field to a cdr @a grid.
 *
 * Really, it would be more efficient to include that field in some other
 * routine and we will avoid one more grid sweep. However, for the sake
 * of code clarity, extendability and maintenability I prefer to add it here.
 *
 * ext_e_r, ext_e_z and ext_e_theta are functions of three real numbers
 * (r, z, \f$\theta\f$) that return each of the components of the electric field.
 *
 *  Note that it is the responsibility of e_theta to behave correctly
 *  (i.e. periodically) when it receives thetas outside \f$[0, 2 \pi]\f$.
 */
void
cdr_add_ext_field (cdr_grid_t *grid)
{
  int ir, iz, itheta;

  debug (2, "cdr_add_ext_field (" grid_printf_str ")\n",
	 grid_printf_args(grid));

#pragma omp parallel private(ir, iz)
  {
    #pragma omp for
    iter_grid_3d_n (grid, ir, iz, itheta, 2) {
      RZT(grid->er, ir, iz, itheta) +=
	ext_e_r (er_r_at (ir, grid->level),
		 er_z_at (iz, grid->level),
		 theta_at (itheta));

      RZT(grid->ez, ir, iz, itheta) +=
	ext_e_z (ez_r_at (ir, grid->level),
		 ez_z_at (iz, grid->level),
		 theta_at (itheta));

      RZT(grid->etheta, ir, iz, itheta) +=
	ext_e_theta (r_at (ir, grid->level),
		     z_at (iz, grid->level),
		     etheta_theta_at (itheta));
    }
  }
}

/** @brief Recursive version of cdr_add_ext_field. */
mk_recursive (cdr_add_ext_field, cdr_grid_t)


/** @brief Adds to a cdr grid the inhomogeneous field created by a charge
 * q located at (r = 0, z = pois_inhom_z).
 *
 * THIS FUNCTION IS KEPT FOR REFERENCE ONLY.  It was superseded by
 * pois_add_inhom_phi_r (see poisson.c)
 */
void
cdr_add_inhom_field_r (cdr_grid_t *grid, double q)
{
  int ir, iz, itheta;
  cdr_grid_t *child;

  debug (2, "cdr_add_inhom_field_r (" grid_printf_str ", q = %f)\n",
	 grid_printf_args(grid), q);

  assert (0 != pois_inhom);

#pragma omp parallel private(ir, iz)
  {
    #pragma omp for
    iter_grid_3d_n (grid, ir, iz, itheta, 2) {
      RZT(grid->er, ir, iz, itheta) +=
	q * pois_inhom_er (er_r_at (ir, grid->level),
			   er_z_at (iz, grid->level));

      RZT(grid->ez, ir, iz, itheta) +=
	q * pois_inhom_ez (ez_r_at (ir, grid->level),
			   ez_z_at (iz, grid->level));
    }
  }
  iter_childs (grid, child) {
    cdr_add_inhom_field_r (child, q);
  }
}

/** @brief Calculates the absolute value of the electric field by linear
 *  interpolation of its components.
 *
 *  Of course, it requires that all the components are already present. */
void
cdr_calc_eabs (cdr_grid_t *grid)
{
  int ir, iz, itheta;
  double er, ez, etheta;

#pragma omp parallel
  {
#pragma omp for private(ir, iz, er, ez, etheta)
    iter_grid_3d (grid, ir, iz, itheta) {
      er = 0.5 * (RZT (grid->er, ir, iz, itheta)
		  + RZT (grid->er, ir - 1, iz, itheta));
      ez = 0.5 * (RZT (grid->ez, ir, iz, itheta)
		  + RZT (grid->ez, ir, iz - 1, itheta));
      if (grid->ntheta != 1) {
	etheta = 0.5 * (RZT (grid->etheta, ir, iz, itheta)
			+ RZT (grid->etheta, ir, iz, itheta - 1));
      } else {
	etheta = 0;
      }
      RZT (grid->eabs, ir, iz, itheta) = sqrt (er * er + ez * ez
					     + etheta * etheta);

    }
  }
}

/** @brief Recursive version of cdr_calc_eabs. */
mk_recursive (cdr_calc_eabs, cdr_grid_t)

/********************
 * Time integration *
*********************/

/** @brief Makes sure that all species of one grid are non-negative. */
void
cdr_nonegative (cdr_grid_t *grid)
{
  int s, ir, iz, itheta;

#pragma omp parallel
  {
#pragma omp for private(ir, iz, s)
    iter_grid_theta_n (grid, itheta, 2) {
      iter_grid_n(grid, ir, iz, 2) {
	   for (s = 0; s < no_species; s++)
	        *RZTP (grid->dens[s], ir, iz, itheta) =
	        MYMAX(0, *RZTP (grid->dens[s], ir, iz, itheta));
         }
       }
     }
}

/** @brief Recursive version of cdr_nonegative. */
mk_recursive (cdr_nonegative, cdr_grid_t)

/** @brief Sets the boundary values for the root grid and those grids that
   have a boundary coinciding with any of the boundaries of the root grid. */
void
cdr_set_ext_bnd (cdr_grid_t *grid)
{
  int r0, z0, r1, z1, ntheta;

  debug (2, "cdr_set_ext_bnd (" grid_printf_str " [grid->ext_bound = 0x%x])\n",
	 grid_printf_args(grid), grid->ext_bound);

  r0 = grid->r0;
  z0 = grid->z0;
  r1 = grid->r1;
  z1 = grid->z1;
  ntheta = grid->ntheta;

  /* Matching conditions (reduced to Hom. Neumann for ntheta = 1 at
     r = r0 */
  if (grid->ext_bound & BND_MASK (BND_LEFT)) {
    set_axis_bnd (grid);
  }

  /* At r = r1. */
  if (grid->ext_bound & BND_MASK (BND_RIGHT))
    set_all_bnd (grid, cdr_ext_bnd_cond[BND_RIGHT], r1 - 1, z0 - 2, 0,
		 R_INDX, 1,
		 Z_INDX, 0, z1 - z0 + 4,
		 THETA_INDX, 0, ntheta);

  /* At z = z0. */
  if (grid->ext_bound & BND_MASK (BND_BOTTOM))
    set_all_bnd (grid, cdr_ext_bnd_cond[BND_BOTTOM], r0 - 2, z0, 0,
		 Z_INDX, -1,
		 R_INDX, 0, r1 - r0 + 4,
		 THETA_INDX, 0, ntheta);

  /* At z = z1. */
  if (grid->ext_bound & BND_MASK (BND_TOP))
    set_all_bnd (grid, cdr_ext_bnd_cond[BND_TOP], r0 - 2, z1 - 1, 0,
		 Z_INDX, 1,
		 R_INDX, 0, r1 - r0 + 4,
		 THETA_INDX, 0, ntheta);
}

/** @brief Recursive version of cdr_set_ext_bnd. */
mk_recursive (cdr_set_ext_bnd, cdr_grid_t)

/** @brief Sets some boundary conditions for all the relevant variables in a
 * grid, to understand the parameters, see rz_array.c:rz_set_bnd(..)
 */
static void
set_all_bnd (cdr_grid_t *grid,
	     int sign, int ir, int iz, int itheta, int dim0, int inout,
	     int dim1, int dim1_from, int dim1_to,
	     int dim2, int dim2_from, int dim2_to)
{
  int s;

  for (s = 0; s < no_species; s++) {
    rz_set_bnd (grid->dens[s], sign, RZTP (grid->dens[s], ir, iz, itheta),
		dim0, inout,
		dim1, dim1_from, dim1_to,
		dim2, dim2_from, dim2_to);
  }
}

/** @brief Sets the proper boundary conditions for a grid that has its
 * left boundary at r = 0 (e.g. the root grid).
 *
 * We don't set Hom. Neumann at r = r0 since r = r0 is
 * not a real boundary: the proper condition here is matching the
 * opposite grid.  But note that this condition reduces to Hom.
 * Neumann in the case ntheta = 1 (because itheta2 = itheta = 0).
 */
static void
set_axis_bnd (cdr_grid_t *grid)
{
  int s, itheta, itheta2, ntheta;

  debug (2, "set_axis_bnd (" grid_printf_str ")\n",
	 grid_printf_args(grid));

  assert (0 == grid->r0);

  ntheta = grid->ntheta;

  for (s = 0; s < no_species; s++) {
    iter_grid_theta (grid, itheta) {
      itheta2 = (itheta + ntheta / 2) % ntheta;

      rz_copy_bnd (grid->dens[s], grid->dens[s],
		   /* sign  = */ 1,
		   /* start_from =  */ RZTP (grid->dens[s],
					     grid->r0, grid->z0 - 2, itheta2),
		   /* start_to = */    RZTP (grid->dens[s],
					     grid->r0, grid->z0 - 2, itheta),
		   R_INDX, BND_OUTWARD, BND_INWARD,
		   Z_INDX, 0, grid->z1 - grid->z0 + 4,
		   THETA_INDX, 0, 1 /* <- relative to start_xxx */);
    }
  }
}

/** @brief Calculates the derivatives of the densities for each point of the
 *  grid.
 */
void
cdr_calc_d_dens (cdr_grid_t *grid)
{
  react_apply_all (grid);

  cdr_advect_diffu_r (grid);
}

/* See "An Adaptive Grid Refinement Strategy for the Simulation
 * of Negative Streamers", C. Montijn, W. Hundsdorfer and U. Ebert
 * for the computations and notation used in this function.
 */

/** @brief This thing looks really ugly, but what does it mean?
 * For example, for stride = dens[electrons]->strides[R_INDX],
 *  it is the \f$p_{ij}\f$ of the aforementioned paper:
 *  \f$p_{i,j} = (\sigma_{i,j} - \sigma_{i-1,j}) / (\sigma_{i+1,j} - \sigma_{i,j})\f$
 */
#define pij(p_, stride_)					\
  (((*(p_)) - *((p_) - (stride_))) / (*((p_) + (stride_)) - *(p_)))

/** @brief This is the F^a + F^d in the paper (times \f$\delta\f${r, z, \f$\theta\f$}). */
static double
f_ad (rz_array_t *dens, double efield, double dc, double mc,
      int ir, int iz, int itheta, int dim)
{
  double sigma, psi_p, sigmasigma;
  int ishift;
  double *p;

  /* a p is the unmodified density-value in the cell (ir,iz,itheta).
   * Normally this should be n_e, as only electrons are mobile.
   */
  p = RZTP (dens, ir, iz, itheta);

  ishift = (efield > 0? 0: 1);

  /* sigma is just
   * sigma_{i    , j} for efield < 0,
   * sigma_{i + 1, j} for efield > 0
   */
  sigma = *(p + ishift * dens->strides[dim]);

  /* psi_p has to be
   *  psi(  p_{i  , j}) for efield < 0
   * -psi(1/p_{i+1, j}) for efield > 0
   */
  psi_p = pij (p + ishift * dens->strides[dim], dens->strides[dim]);

  /* E > 0            ---> psi_p = -psi( 1 / p_ij )
   * E < 0, p_ij == 0 ---> psi_p = 1
   * E < 0, p_ij != 0 ---> psi_p = psi( 1 / p_ij )
   */
  psi_p = (efield > 0)? -psi (psi_p) : psi_p == 0? 1.0: psi (1.0 / psi_p);

  /* sigmasigma is (sigma_{i+1, j} - sigma_{i, j}) */
  sigmasigma = (*p) - *(p + dens->strides[dim]);

  /* The first term corresponds to the advection term, the second to the
   * diffusion term.
   *
   * See eq. 3.1.5~3.1.8 in the thesis of C. Montijn. Division by cell-size
   * is done in the function calling this one.
   */
  return efield * mc * (sigma + psi_p * sigmasigma) + dc * sigmasigma;
}

/** @brief Adds the advection and diffusion terms to the derivatives of the
 * densities.
 */
void
cdr_advect_diffu (cdr_grid_t *grid)
{
  int ir, iz, itheta, s;
  int ii;
  double gdr_inv, gdz_inv, d_iso0, d_iso;

  double field;
  int dimd, t, sg;
  REAL *d_dens_p;
  double r_inv, er_r;
  rz_array_t *er,*ez,*etheta;

  debug (2, "cdr_advect_diffu (" grid_printf_str ")\n",
	 grid_printf_args(grid));

  gdr_inv = 1.0 / dr[grid->level];
  gdz_inv = 1.0 / dz[grid->level];
  er     = grid->er;
  ez     = grid->ez;
  etheta = grid->etheta;

  for (s = 0; s < no_species; s++) {
    /* The charge / mass ratio of a species. */
    double cm0, cm;
    int *strides;
    int ir, rmax, dim2;
    rz_array_t *dens, *d_dens, **d_dens_copy;

    if (spec_index[s]->mass <= 0) continue;

    dens = grid->dens[s];
    d_dens = grid->d_dens[s];
    strides = dens->strides;
    rmax = strides[1]; dim2 = strides[2]/strides[1];

    /* Make a copy of d_dens to be able to compare both versions
     * of arcos_cdr_advect
    int nr = grid->r1 - grid->r0 + 4;
    int nz = grid->z1 - grid->z0 + 4;
    d_dens_copy = (rz_array_t **) xmalloc (sizeof(rz_array_t*) * no_species);
    rz_copy(d_dens,grid->r0,grid->z0,*d_dens_copy,grid->r0,grid->z0,nr,nz);

    strides = dens->strides;
    */

    arcos_cdr_advect_vec(spec_index[s]->mass,spec_index[s]->charge,
                     dens->data,d_dens->data,er->data,er->len,
		     ez->data,ez->len,diffusion_coeff,
                     dr[grid->level],dz[grid->level],
                     sprite_module,grid->r0,grid->r1,grid->z0,grid->z1);

    /* rz_free(*d_dens_copy); */
  }
}

mk_recursive (cdr_advect_diffu, cdr_grid_t)

/** @brief Limiter function, psi(\f$\theta\f$) = max (0, min(1, 1/3 + \f$\theta\f$/6, \f$\theta\f$))
 */
static double 
psi (double theta) 
{
  if (theta < 0) return 0;
  if (theta < 0.4) return theta;
  if (theta < 4) return (1 + 0.5 * theta) / 3.0;
  else return 1.0;
}

/** @brief Returns the minimum tau that satisfies the Courant-Friedrichs-Lewy
   restriction (CFL). */
double
cdr_courant (cdr_grid_t *grid)
{
  int s, ir, iz, itheta, depth, is_3d;

  double dr_min, dz_min, rdtheta_min, emax_r, emax_z, emax_theta, dmax, ddmax;
  double tau_a, tau_d, tau_rt, tau_f, max_diff;

  debug (2, "cdr_courant (" grid_printf_str ")\n", grid_printf_args(grid));

  is_3d = (grid->ntheta == 1? 0: 1);

  /* Should we better keep track of the level of the finest grid?  */
  depth = grid_max_depth_r ((grid_t *) grid);
  dr_min = dr[depth];
  dz_min = dz[depth];
  rdtheta_min = dtheta * grid_rmin_r ((grid_t *) grid);

  /* This part is not parallel, but if needed it can be parallelized
   * easily by creating vectors xmax[0...ntheta-1], where we would store
   * the maximum for each ntheta.
   */
  emax_r = -1;
  emax_z = -1;
  emax_theta = -1;
  dmax = -1;
  ddmax = 0;

  iter_grid_3d (grid, ir, iz, itheta) {
    double mu;
    if (sprite_module) {
	    mu = 1.0 / spr_density_at (z_at (iz, grid->level));
    } else {
	    mu = 1.0;
    }
    //mu = sprite_module? 1.0 / spr_density_at (z_at (iz, grid->level)): 1.0;

    max_update (&emax_r, mu * fabs (RZT (grid->er, ir, iz, itheta)));
    max_update (&emax_z, mu * fabs (RZT (grid->ez, ir, iz, itheta)));
    max_update (&emax_theta, mu * fabs (RZT(grid->etheta, ir, iz, itheta)));

    for (s = 0; s < no_species; s++) {
      if (spec_index[s]->mass > 0) {
        max_update (&dmax, mu * RZT(grid->dens[s], ir, iz, itheta));
        max_update (&ddmax, -(RZT(grid->d_dens[s], ir, iz, itheta) / 
    			  RZT(grid->dens[s], ir, iz, itheta)));
      }
    }
  }

  /* This assert fails when the grid has become completely nan.  In that
   * case, let us not waste cpu cycles any longer.
   */
  assert (emax_r >= 0 && emax_z >= 0 && emax_theta >= 0 && dmax >= 0);

  /* Advection.  This one usually dominates. */
  tau_a = nu_a / (emax_r / dr_min + emax_z / dz_min 
		  + is_3d * emax_theta / rdtheta_min);

  /* Diffusion. */
  /* <ISOTROPIC DIFFUSION> */
  if (sprite_module) {
    double dens_z1, dens_z0;
    /* Since we do not know where the density is higher, we check at both
     * extremes of the domain.  We are assuming that the density is monotonic
     * but nothing else (if it is not an exp. it works still).
     */
    dens_z0 = spr_density_at (z_at (grid->z0, grid->level));
    dens_z1 = spr_density_at (z_at (grid->z1, grid->level));
    max_diff = diffusion_coeff / MYMIN(dens_z0, dens_z1);
  } else {
    max_diff = 0.1; /* diffusion_coeff */
  }

  tau_d = nu_d / (max_diff / dr_min / dr_min 
		  + max_diff / dz_min / dz_min
		  + is_3d * max_diff / rdtheta_min / rdtheta_min);

  /* </ISOTROPIC DIFFUSION> */

  /* Mainly reaction growth / decay time */
  /* tau_f = nu_f / ddmax; */

  /* If sigma is zero somewhere, tau_f turns inf even if nu_f is very large.
   * We correct this here. */
  /*if (tau_f == 0)*/ tau_f = tau_d;  /* Just ignore tau_f please. */

  /* Dielectric relaxation time */
  tau_rt = nu_rt / dmax;

  return MYMIN (tau_a, MYMIN(tau_d, MYMIN(tau_f, tau_rt)));
}

static
void max_update (double *x, double cmp)
{
  if (*x < cmp) *x = cmp;
}

/** @brief Writes into the densities of the grid dest the densities of grid orig
 * plus its derivatives times h.
 *
 * If orig = dest, this implements an Euler step.
 * But it can also be used as a tool to implement a RK */
void
cdr_update (cdr_grid_t *orig, cdr_grid_t *dest, double h)
{
  int ir, iz, itheta, s;

  debug (2, "cdr_update (orig = " grid_printf_str ", dest = " 
	 grid_printf_str ", %f)\n",
	 grid_printf_args(orig), grid_printf_args(dest), h);

  /* We assume that orig and dest have the same shape. */
  assert (orig->r0 == dest->r0 && orig->z0 == dest->z0 && 
	  orig->r1 == dest->r1 && orig->z1 == dest->z1 &&
	  orig->ntheta == dest->ntheta);

#pragma omp parallel
  {
#pragma omp for private(ir, iz, s)
    iter_grid_3d_n (orig, ir, iz, itheta, 2) {
      for (s = 0; s < no_species; s++) {
	RZT (dest->dens[s], ir, iz, itheta) = 
	  RZT (orig->dens[s], ir, iz, itheta) 
	  + h * RZT (orig->d_dens[s], ir, iz, itheta);
      }
    }
  }
  cdr_calc_charge (dest);
}

/** @brief Sets dest->dens = dens_0->dens + h/2 * d_dens_1->d_dens +\n
 *                                          h/2 * d_dens_2->d_dens
 *
 * Useful to implement 2nd order RK time-stepping. */
void
cdr_rk2_update (cdr_grid_t *dens_0, cdr_grid_t *d_dens_1, 
		cdr_grid_t *d_dens_2, cdr_grid_t *dest, 
		double h)
{
  int ir, iz, itheta, s;

  debug (2, "cdr_rk2_update (dens_0 = " grid_printf_str 
	 ",\n\t\td_dens_1 = " grid_printf_str
	 ",\n\t\td_dens_2 = " grid_printf_str
	 ",\n\t\tdest = " grid_printf_str ")\n",
	 grid_printf_args(dens_0), 
	 grid_printf_args(d_dens_1), 
	 grid_printf_args(d_dens_2), 
	 grid_printf_args(dest));

#pragma omp parallel
  {
#pragma omp for private(ir, iz, s)
    iter_grid_3d_n (dest, ir, iz, itheta, 2) {
      for (s = 0; s < no_species; s++) {
	RZT (dest->dens[s], ir, iz, itheta) = 
	  RZT (dens_0->dens[s], ir, iz, itheta) 
	  + 0.5 * h * RZT (d_dens_1->d_dens[s], ir, iz, itheta)
	  + 0.5 * h * RZT (d_dens_2->d_dens[s], ir, iz, itheta);
      }
    }
  }
  cdr_calc_charge (dest);
}

/** @brief Recursive version of cdr_rk2_update
 *
 * Note: we assume here that all the trees are congruent: this is,
 * they have the same structure with nodes of the same size.
 */
void
cdr_rk2_update_r (cdr_grid_t *dens_0, cdr_grid_t *d_dens_1, 
		  cdr_grid_t *d_dens_2, cdr_grid_t *dest, 
		  double h)
{
  cdr_grid_t *c_dens_0, *c_d_dens_1, *c_d_dens_2, *c_dest;

  cdr_rk2_update (dens_0, d_dens_1, d_dens_2, dest, h);

  c_dens_0 = dens_0->first_child;
  c_d_dens_1 = d_dens_1->first_child;
  c_d_dens_2 = d_dens_2->first_child;

  iter_childs (dest, c_dest) {
    cdr_rk2_update_r (c_dens_0, c_d_dens_1, c_d_dens_2, c_dest, h);

    c_dens_0 = c_dens_0->next;
    c_d_dens_1 = c_d_dens_1->next;
    c_d_dens_2 = c_d_dens_2->next;
  }
}

/** @brief Recursive version of the former QQQQ, used to implement Euler
 * integration
 */
void
cdr_self_update_r (cdr_grid_t *grid, double h)
{
  cdr_grid_t *child;
  cdr_update (grid, grid, h);

  iter_childs (grid, child) {
    cdr_self_update_r (child, h);
  }
}

/** @brief Yet another recursive version of cdr_update, used to implement 
 *  higher order Runge-Kutta integrations.
 *  
 * It creates a new tree of cdr grids.
 */
void 
cdr_like_update_ar (cdr_grid_t *grid, cdr_grid_t *new_grid, double h)
{
  cdr_grid_t *child, *new_child, *prev = NULL;

  cdr_update (grid, new_grid, h);

  iter_childs (grid, child) {
    new_child = cdr_like_a (child);

    /* We can't use add_child here because we read the children of the grid
     * in inverse order as we inserted them.
     */
    if (NULL == prev) {
      new_grid->first_child = new_child;
    } else {
      prev->next = new_child;
    }
    set_leaf (new_child, new_grid, NULL, NULL, new_grid->level + 1);
    grid_inherit_ext_bound ((grid_t *) new_child);

    prev = new_child;

    cdr_like_update_ar (child, new_child, h);
  }
}

///** @brief Computes the grid error. QQQQ */
//double
//cdr_grid_error (cdr_grid_t *grid, cdr_grid_t *inter)
//{
//  int ir, iz, itheta;
//  int spec_nr;
//  double vol, err, val;
//
//  spec_nr = find_species_by_name("electrons");
//
//  assert(spec_nr > -1);
//
//  vol = dz[grid->level] * dr[grid->level];
//
//  iter_grid_3d(grid, ir, iz, itheta)
//    {
//
//      val = RZT(grid->dens[spec_nr], ir, iz, itheta);
//      err = val > 1e-12 ? (val - RZT(inter->dens[spec_nr], ir, iz, itheta)) / val : 0;
//     err = err > 0 ? err : -err;
//      if (max_err < err) 
//	{
//	  max_err = err;
//	}
//    }
//    return err;
//}
//
///** @brief Recursive version of cdr_grid_error. QQQQ */
//void
//cdr_grid_error_r (cdr_grid_t *grid, cdr_grid_t *inter)
//{
//  cdr_grid_t *g_child, *i_child;
//
//  cdr_grid_error(grid, inter);
//
//  i_child = inter->first_child;
//
//  for (g_child = grid->first_child; g_child; g_child = g_child->next)
//    {
//      cdr_grid_error_r(g_child, i_child);
//
//      i_child = i_child->next;
//    }
//}

/** @brief Makes a full second order Runge-Kutta step.
 *
 * Returns the timestep actually performed, which can be smaller that h
 * if it does not satisfy the Courant criterium.
 */
double
cdr_rk2 (cdr_grid_t *grid, double h, double t)
{
  /* intermediate step. */
  cdr_grid_t *inter, *child;
  double courant_h;

  /* After this number of warnings about small timesteps, the code stops. */
  static int count_min_timestep = 5;

  debug (2, "cdr_rk2 (" grid_printf_str ", %f)\n",
	 grid_printf_args(grid), h);

  cdr_nonegative_r (grid);

  /* We want to be sure that the boundaries are properly set. */
  prepare_grid (grid);

  /* Experimental: Moving z_max boundary -- DISABLED */

  cdr_calc_field_r (grid, /*return_pois = */ FALSE);

  /* If we are using the sprite module, we set here the proper magnitudes
     corresponding to the head altitude here. */
  if (sprite_module) spr_hook (grid);

  cdr_calc_d_dens (grid);

  /* The actual timestep is limited by the Courant stability criterium.
   *
   * Note that cdr_courant has to be called AFTER the computation of the fields.
   */
  courant_h = cdr_courant (grid);
  h = MYMIN (h, courant_h);

  if (h < warn_min_timestep) {
    warning ("Time step [h = %g] is rather small.  If this is ok, reduce "
	     "warn_time_step to remove this warning\n", h);

    if (0 == count_min_timestep--) {
      fatal ("Too many warnings about small timesteps.\n");
    }
  }

  inter = cdr_like_a (grid);

  /* We create a new tree and put there y + dy/dt * h at each level. */
  cdr_like_update_ar (grid, inter, h);

  /* Set the proper boundary conditions also in the intermediate step. */
  prepare_grid (inter);

  cdr_calc_field_r (inter, /* return_pois = */ FALSE);
  cdr_calc_d_dens (inter);

  /* This is something like grid = grid + h / 2 * F(grid) + h / 2 F(inter) */
  cdr_rk2_update_r (grid, grid, inter, grid, h);

  /* Also, for the refinement functions, the boundary conditions have to be
   * appropriately set.
   */
  prepare_grid (grid);

  cdr_free_r (inter);

  return h;
}

/** @brief Prepares a grid family :
 * a) Set the boundaries for all subgrids,\n
 * b) Restrict the coarser grids from the fine ones.\n
 */
static void
prepare_grid (cdr_grid_t *grid)
{
  /* Restrict the densities from finer to coarser grids. */
  cdr_restrict_r (grid);

  /* Sets the external boundaries. */
  cdr_set_ext_bnd_r (grid);

  /* Sets the boundaries of children grids interpolating from their parents. */
  cdr_set_bnd_r (grid);

  /* If two grids share a boundary, use the values of one to set the 
     boundaries of the other. */
  cdr_match_r (grid, NULL);

  /* Sets periodic boundary conditions in theta. */
  cdr_set_periodic_r (grid);
}  

/*****************************************************************
 * Refinement and mapping. 
 *****************************************************************/

/** @brief Takes a pointer to a family of cdr grids, creates a new family
 * which is stored in that pointer and takes rid of the initial family. 
 *
 * This function is the entry point of the refinement code.
*/
void
cdr_update_refined (cdr_grid_t **ptree)
{
  cdr_grid_t *old, *new;

  old = *ptree;

  new = cdr_clone_a (old);
  cdr_calc_charge (new);
  cdr_refine_r (new, old);

  *ptree = new;

  cdr_free_r (old);
}

/** @brief Updates the values of max_dens and max_charge.
 *
 * This values are needed for the refinement criterium.
 */
void
cdr_calc_maxs (cdr_grid_t *grid)
{
  int start = TRUE;
  int s, ir, iz, itheta;

  debug (2, "cdr_calc_maxs (" grid_printf_str ")\n",
	 grid_printf_args(grid));

  /* If we have a parent, we inherit the maxs from him. */
  if (grid->parent) {
    grid->max_charge = grid->parent->max_charge;
    for (s = 0; s < no_species; s++) {
      grid->max_dens[s] = grid->parent->max_dens[s];
    }
    return;  
  }

  iter_grid_3d (grid, ir, iz, itheta) {
    if (start) {
      start = FALSE;
      grid->max_charge = fabs(RZT (grid->charge, ir, iz, itheta));
      for (s = 0; s < no_species; s++) {
	grid->max_dens[s] = RZT (grid->dens[s], ir, iz, itheta);
      }
    } else {
      max_update (&grid->max_charge, fabs(RZT (grid->charge, ir, iz, itheta)));
      for (s = 0; s < no_species; s++) {
	max_update (&grid->max_dens[s], RZT (grid->dens[s], ir, iz, itheta));
      }
    }
  }

  debug (2, "->max_charge = %g, ->max_dens[0] = %g, ->max_dens[1] = %g\n",
	 grid->max_charge, grid->max_dens[0], grid->max_dens[1]);
}

/** @brief Determines the curvature of ar, which has to be a component of 
 *  grid at @a ir, @a iz.
 *
 *  If @a f is nonnull, determines the curvature of \f$f(ar)\f$. \n
 *  Note that we only consider the curvatures in @a r and @a z (not in
 *  \f$\theta\f$).
 *  This is because \f$\theta\f$ is never refined. (can that be changed in some 
 *  future?)
 */
static double
curv_at (cdr_grid_t *grid, rz_array_t *ar, int ir, int iz, int itheta,
	 double (*f) (double))
{
  int i, level;
  REAL x[3];
  double dur, duz;

  for (i = 0; i < 3; i++) {
    x[i] = RZT (ar, ir - 1 + i, iz, itheta);
    x[i] = f? f(x[i]): x[i];
  }

  level = grid->level;
  /* <CYLINDRICAL> */
  dur = (cyl_er_r_at (ir, level) * (x[2] - x[1]) 
	 - cyl_er_r_at (ir - 1, level) * (x[1] - x[0])) / cyl_r_at (ir, level);
  /* </CYLINDRICAL> */

  for (i = 0; i < 3; i++) {
    x[i] = RZT (ar, ir, iz - 1 + i, itheta);
    x[i] = f? f(x[i]): x[i];
  }

  duz = (x[2] - 2 * x[1] + x[0]);

  return fabs(dur) + fabs(duz);
}

/** @brief Does the point in grid at @a ir, @a iz, @a itheta require to be
 * further refined? 
 */
static int
needs_refinement (cdr_grid_t *grid, int ir, int iz, int itheta, int *in_edge)
{
  double curv;
  int s;
  /* For sprites, the eabs refinement criterium may depend on the location. */
  int s_ref_level_eabs;
  double s_ref_threshold_eabs;

  /* If the maximum charge is below this number, we ignore the refinement
   * criterium.
   */
  const double epsilon_charge = 1e-8;

  if (!sprite_module || ref_level_eabs < 0) {
    s_ref_level_eabs = ref_level_eabs;
    s_ref_threshold_eabs = ref_threshold_eabs;
  } else {
    double z, back_dens;
    z = z_at (iz, grid->level);
    back_dens = spr_density_at (z);

    s_ref_threshold_eabs = ref_threshold_eabs * back_dens;

    /* log2 appears because whenever we increase the background density by a
     * factor 2, we should refine up to one more level.
     */
    s_ref_level_eabs = ref_level_eabs + (int) floor (log2 (back_dens));
  }

  /* Electric field refinement. If eabs > ref_threshold_eabs, refine up to
   * some predefined level.
   */

  if (grid->level < s_ref_level_eabs && 
      RZT (grid->eabs, ir, iz, itheta) > s_ref_threshold_eabs) {
    debug (4, "Refine grid " grid_printf_str 
	   " at ir = %d, iz = %d, itheta = %d\n"
	   "\t because too high EABS [eabs = %g]\n",
	   grid_printf_args(grid), ir, iz, itheta, 
	   RZT (grid->eabs, ir, iz, itheta));

    *in_edge = FALSE;
    return TRUE;
  }

  curv = curv_at (grid, grid->charge, ir, iz, itheta, NULL);

  if (grid->max_charge > epsilon_charge 
      && curv / grid->max_charge > ref_threshold_charge) {
    debug (4, "Refine grid " grid_printf_str 
	   " at ir = %d, iz = %d, itheta = %d\n"
	   "\t because too high CHARGE curvature [curv / max_charge = %g]\n",
	   grid_printf_args(grid), ir, iz, itheta, curv / grid->max_charge);

    *in_edge = FALSE;
    return TRUE;
  }

  for (s = 0; s < no_species; s++) {
    /* We ignore immobile species: */
    if (spec_index[s]->mass <= 0) continue;

    curv = curv_at (grid, grid->dens[s], ir, iz, itheta, NULL);
    if (curv / grid->max_dens[s] > ref_threshold_dens) {
      debug (4, "Refine grid " grid_printf_str 
	     " at ir = %d, iz = %d, itheta = %d\n"
	     "\t because too high DENS[%s] curvature [curv / max_dens = %g]\n",
	     grid_printf_args(grid), ir, iz, itheta, spec_index[s]->name, 
	     curv / grid->max_dens[s]);

      *in_edge = FALSE;
      return TRUE;
    }

    /* If the grid contains the leading edge, the criterium of
     * refinement also takes into account the density.
     */
    if (*in_edge && RZT (grid->dens[s], ir, iz, itheta) > ref_threshold_edge) {
      debug (4, "Refine grid " grid_printf_str 
	     " at ir = %d, iz = %d, itheta = %d\n"
	     "\t because too high DENS[%s] in leading edge [dens = %g]\n",
	     grid_printf_args(grid), ir, iz, itheta, spec_index[s]->name, 
	     RZT (grid->dens[s], ir, iz, itheta));
      return TRUE;
    }
  }    

  /* If we are close to the border, we also check the border itself. 
   *
   * NOTE: Does this have any sense?  We are anyway checking for
   * the boundaries, aren't we?
   */
  if (ir == grid->r0 + 1 && needs_refinement (grid, ir - 1, iz, itheta, 
					      in_edge)) 
    return TRUE;

  if (iz == grid->z0 + 1 && needs_refinement (grid, ir, iz - 1, itheta,
					      in_edge)) 
    return TRUE;

  if (ir == grid->r1 - 2 && needs_refinement (grid, ir + 1, iz, itheta,
					      in_edge)) 
    return TRUE;

  if (iz == grid->z1 - 2 && needs_refinement (grid, ir, iz + 1, itheta,
					      in_edge)) 
    return TRUE;

  return FALSE;
} 

/** @brief Calls the needs_refinement at all \f$\theta\f$ and returns true if
 * the grids needs to be refined in _any_ of the angles.
 *
 * This is the safest thing to do, but some other criterium could be imagined
 * as well.  The search is stopped as soon as a cell that requires refinement
 * is found, but not if it requires refinement due to the leading edge density
 * criterium.  This is because we need to find where we have to stop applying
 *  that criterium.
 */
static int
any_needs_refinement (cdr_grid_t *grid, int ir, int iz, int *in_edge)
{
  int itheta;
  int ret = FALSE;
  iter_grid_theta (grid, itheta) {
    if (needs_refinement (grid, ir, iz, itheta, in_edge)) {
      debug (6, grid_printf_str 
	     " needs refinement at ir = %d, iz = %d itheta = %d\n",
	     grid_printf_args(grid), ir, iz, itheta);
      ret = TRUE;
    }
    if (ret && ! *in_edge) return TRUE;
  }
  return ret;
}

/** @brief Determines whether a @a brick needs refinement.
 *
 * Now we look if we need refinement in a @a brick \n
 * i.e. a rectangle in @a r, @a z space
 * (usually, @a height = @a cdr_brick_dz,
 *           @a width  = @a cdr_brick_dr),
 * but can be smaller if we reach a boundary. 
 *
 * See above, any_needs_refinement to see when do we stop the search.
 */
static int
brick_needs_refinement (cdr_grid_t *grid, int r0, int z0, int r1, int z1,
			int *in_edge)
{
  int ir, iz;
  int ret = FALSE;

  for (ir = r0; ir < r1; ir ++) {
    for (iz = z1 - 1; iz >= z0; iz--) {
      ret = (any_needs_refinement (grid, ir, iz, in_edge) || ret);
      if (ret && ! *in_edge) return TRUE;
    }
  }
  return ret;
}

/** @brief Uses the previous routines to refine (if needed) a cdr grid */
void
cdr_refine (cdr_grid_t *grid)
{
  int ir;

  debug (2, "cdr_refine (" grid_printf_str ")\n",
	 grid_printf_args(grid));

  assert (grid->level <= cdr_max_level);

  /* Checks if the maximum refinement level has been reached. */
  if (grid->level == cdr_max_level) return;

  cdr_calc_maxs (grid);

  for (ir = grid->r0; ir < grid->r1; ir += cdr_brick_dr) {
    int irmax = MYMIN (grid->r1, ir + cdr_brick_dr);
    int tainted, building, z0 = -1, z1 = -1;
    int iz;
    int in_edge;
    int first_edge;

    first_edge = in_edge = grid->contains_edge;

    /* The leading edge extends downwards. */
    for (iz = grid->z0, building = FALSE; 
	 iz < grid->z1; iz += cdr_brick_dz) {
      int izmax = MYMIN (grid->z1, iz + cdr_brick_dz);

      tainted = brick_needs_refinement (grid, ir, iz, irmax, izmax, 
					&first_edge);
      if (tainted) {
	z1 = izmax;
	if (!building) {
	  z0 = iz;
	  building = TRUE;
	}
      } else /* ! tainted */ {
	if (building) {
	  assert (z0 >= 0 && z1 > 0);
	  /* If a grid does not satisfy the curvature criterium anywhere
	   * and is here only because of the density threshold (edge) 
	   * criterium, he does not deserve to be refined.
	   */
	  building = FALSE;
	  if (!first_edge) {
	    refine_in (grid, ir, z0, irmax, z1, in_edge);
	    in_edge = FALSE;
	  }
	}
      } 
    }
    if (building) {
      assert (z0 >= 0 && z1 > 0);
      if (!first_edge) refine_in (grid, ir, z0, irmax, z1, in_edge);
    }
  }
}

/** @brief Adds a child to a cdr grid.
 *
 * The coordinates are in the parent's units
 */
static void
refine_in (cdr_grid_t *grid, int cr0, int cz0, int cr1, int cz1, 
	   int contains_edge)
{
  cdr_grid_t *child;
  int nr0, nr1, nz0, nz1;

  nr0 = cr0 << 1;
  nz0 = cz0 << 1;
  nr1 = cr1 << 1;
  nz1 = cz1 << 1;

  child = cdr_new_3d_a (nr0, nz0, nr1, nz1, grid->ntheta);

  add_child (grid, child);
  grid_inherit_ext_bound ((grid_t*) child);

  child->contains_edge = contains_edge;

  debug (3, "new child created {r0 = %d, z0 = %d, r1 = %d, z1 = %d, "
	 "level = %d}\n", child->r0, child->z0, child->r1, child->z1, 
	 child->level);
}

/** @brief Recursively refines a cdr grid and constructs a tree based on it.
 *
 *  if source == NULL\n
 *  then uses the initial densities\n
 *  else sets the densities of the newly created grids by interpolating and
 *  copying from the tree source. 
 */
void
cdr_refine_r (cdr_grid_t *grid, cdr_grid_t *source)
{
  cdr_grid_t *child;
  int itheta;

  cdr_calc_charge (grid);

  cdr_refine (grid);

  iter_childs (grid, child) {
    if (NULL != source) {
#pragma omp parallel
      {    
#pragma omp for
	iter_grid_theta (child, itheta) {
	  map_grid_r (dens_mappers, (grid_t *) source, (grid_t *) child, 
		      itheta, 
		      /* copy = */ TRUE, 
		      /* interpol = */ TRUE, 
		      /* coarsen = */ FALSE, 
		      /* s_buff = */ 0, /* t_buff = */ 2);
	}
      }
    } else {
      cdr_init_dens (child);
    }
    cdr_refine_r (child, source);
  }
}

/* Matching of the boundaries:
 *
 * The boundary conditions of a grid that has a common border with another
 * one at the same level have to be matched. This is the purpose of this set
 * of functions.
 */

/** @brief Examines pairs of grids and matches them and all their descendants.
 * 
 * If grid2 == NULL, matches only the descendats of grid1. Hence, this
 *  routine is usually called as cdr_match_r (root, NULL).
 */
void
cdr_match_r (cdr_grid_t *grid1, cdr_grid_t *grid2)
{
  cdr_grid_t *child1, *child2;
  int match = TRUE;

  if (grid2) {
    debug (2, "cdr_match_r (grid1 = " grid_printf_str 
	   ", grid2 = " grid_printf_str ")\n",
	   grid_printf_args (grid1), grid_printf_args (grid2));
  } else {
    debug (2, "cdr_match_r (grid1 = " grid_printf_str 
	   ", grid2 = NULL)\n", grid_printf_args (grid1));
  }

  if (NULL != grid2) {
    match = match_grids (grid1, grid2);
  }

  if (!match) return;

  iter_childs (grid1, child1) {
    debug (4, "\tchild1 = " grid_printf_str "\n",
	   grid_printf_args (child1));
    iter_childs (grid1, child2) {
      debug (4, "\t\tchild2 = " grid_printf_str "\n",
	     grid_printf_args (child2));

      cdr_match_r (child1, (child1 != child2)? child2: NULL);
    }

    if (NULL == grid2)
      continue;

    iter_childs (grid2, child2) {
      cdr_match_r (child1, child2);
    }
  }
}

/** @brief Matches the boundaries of grid to reading from grid @a fro.
 *
 * Returns TRUE if there was a matching, FALSE otherwise.
 * Maybe this function should be rewritten to make it more compact and clear.
 */
static int
match_grids (cdr_grid_t *fro, cdr_grid_t *to)
{
  int s, x0, x1;

  debug (2, "match_grids (to = " grid_printf_str 
	 ", fro = " grid_printf_str ")\n",
	 grid_printf_args(to), grid_printf_args(fro));

  assert (to->level == fro->level);

  if (to->r1 == fro->r0) {
    x0 = MYMAX (to->z0 - 2, fro->z0 - 2);
    x1 = MYMIN (to->z1 + 2, fro->z1 + 2);
    if (x1 <= x0) return FALSE;

    for (s = 0; s < no_species; s++) {
      rz_copy_bnd (fro->dens[s], to->dens[s], 1,
		   RZTP (fro->dens[s], fro->r0, x0, 0),
		   RZTP (to->dens[s], to->r1 - 1, x0, 0),
		   R_INDX, BND_OUTWARD, BND_OUTWARD,
		   Z_INDX, 0, x1 - x0,
		   THETA_INDX, 0, fro->ntheta);
    }
    return TRUE;
  }

  if (to->r0 == fro->r1) {
    x0 = MYMAX (to->z0 - 2, fro->z0 - 2);
    x1 = MYMIN (to->z1 + 2, fro->z1 + 2);
    if (x1 <= x0) return FALSE;

    for (s = 0; s < no_species; s++) {
      rz_copy_bnd (fro->dens[s], to->dens[s], 1,
		   RZTP (fro->dens[s], fro->r1 - 1, x0, 0),
		   RZTP (to->dens[s], to->r0, x0, 0),
		   R_INDX, BND_INWARD, BND_INWARD,
		   Z_INDX, 0, x1 - x0,
		   THETA_INDX, 0, fro->ntheta);
    }
    return TRUE;
  }

  if (to->z1 == fro->z0) {
    x0 = MYMAX (to->r0 - 2, fro->r0 - 2);
    x1 = MYMIN (to->r1 + 2, fro->r1 + 2);
    if (x1 <= x0) return FALSE;

    for (s = 0; s < no_species; s++) {
      rz_copy_bnd (fro->dens[s], to->dens[s], 1,
		   RZTP (fro->dens[s], x0, fro->z0, 0),
		   RZTP (to->dens[s], x0, to->z1 - 1, 0),
		   Z_INDX, BND_OUTWARD, BND_OUTWARD,
		   R_INDX, 0, x1 - x0,
		   THETA_INDX, 0, fro->ntheta);
    }
    return TRUE;
  }

  if (to->z0 == fro->z1) {
    x0 = MYMAX (to->r0 - 2, fro->r0 - 2);
    x1 = MYMIN (to->r1 + 2, fro->r1 + 2);
    if (x1 <= x0) return FALSE;

    for (s = 0; s < no_species; s++) {
      rz_copy_bnd (fro->dens[s], to->dens[s], 1,
		   RZTP (fro->dens[s], x0, fro->z1 - 1, 0),
		   RZTP (to->dens[s], x0, to->z0, 0),
		   Z_INDX, BND_INWARD, BND_INWARD,
		   R_INDX, 0, x1 - x0,
		   THETA_INDX, 0, fro->ntheta);
    }
    return TRUE;
  }
  return FALSE;
}

/** @brief Sets the boundary conditions of all children of grid by
 * interpolating from grid itself.
 */
void
cdr_set_bnd (cdr_grid_t *grid)
{
  cdr_grid_t *top, *bottom, *left, *right, *parent;
  int itheta;

  debug (2, "cdr_set_bnd (" grid_printf_str ")\n",
	 grid_printf_args(grid));

  parent = grid->parent;
  if (NULL == parent) {
    assert (0 == grid->level);
    return;
  }

  /* Since our mapping functions are easier to implement when we map
   * between rectangular regions (grids) we map the boundary conditions 
   * by creating from grid four "guest" grids that share memory with him
   * and that contains each of the four boundaries.
   */
  top = cdr_guest (grid, grid->r0 - 2, grid->z1, 
		   grid->r1 + 2, grid->z1 + 2);
  bottom = cdr_guest (grid, grid->r0 - 2, grid->z0 - 2, 
		      grid->r1 + 2, grid->z0);
  left = cdr_guest (grid, grid->r0 - 2, grid->z0, 
		    grid->r0, grid->z1);
  right = cdr_guest (grid, grid->r1, grid->z0, 
		       grid->r1 + 2, grid->z1);

#pragma omp parallel
    {    
#pragma omp for
      iter_grid_theta (grid, itheta) {
	if (0 == (grid->ext_bound & BND_MASK (BND_TOP)))
	  map_grid (dens_bnd_mappers, (grid_t *) parent, (grid_t *) top, itheta, 
		    FALSE, TRUE, FALSE, 1, 0);
	
	if (0 == (grid->ext_bound & BND_MASK (BND_BOTTOM)))
	  map_grid (dens_bnd_mappers, (grid_t *) parent, (grid_t *) bottom, itheta,
		    FALSE, TRUE, FALSE, 1, 0);

	if (0 == (grid->ext_bound & BND_MASK (BND_LEFT))) {
	  assert (grid->r0 != 0);
	  map_grid (dens_bnd_mappers, (grid_t *) parent, (grid_t *) left, itheta, 
		  FALSE, TRUE, FALSE, 1, 0);
	  
	}

	if (0 == (grid->ext_bound & BND_MASK (BND_RIGHT)))
	  map_grid (dens_bnd_mappers, (grid_t *) parent, (grid_t *) right, itheta, 
		    FALSE, TRUE, FALSE, 1, 0);
      }
    }

    cdr_free (top);
    cdr_free (bottom);
    cdr_free (left);
    cdr_free (right);
}

/** @brief Recursive version of cdr_set_bnd. */
mk_recursive (cdr_set_bnd, cdr_grid_t)

/** @brief Mapping of the densities. */
mapper_t**
cdr_mappers_a (interpol_method_t *interp_method)
{
  mapper_t *mapper, **mappers;
  int s, nmappers;

  nmappers = no_species;

  /* mapper->extra represent the species whose density we will map.
   * There are no_species species, but no_species + 1 represents
   * the abs value of the electric field, which is mapped only to use
   * as an electric field-based refinement criterium.
   */
  if (ref_level_eabs >= 0 && ref_threshold_eabs >= 0.0) nmappers++;

  mappers = (mapper_t **) xmalloc (sizeof(mapper_t*) * (nmappers + 1));

  for (s = 0; s < nmappers; s++) {
    mapper = (mapper_t *) xmalloc (sizeof(mapper_t));
    mapper->extra = s;
    mapper->interpol_method = interp_method;
    mapper->copy = dens_copy;
    mapper->interpol_set = dens_interpol_set;
    mapper->interpol = dens_interpol;
    mapper->coarsen = NULL;
    mapper->shift_r = 0;
    mapper->shift_z = 0;
    mappers[s] = mapper;
  }

  mappers[nmappers] = NULL;

  return mappers;
}

/** @brief Frees the memory of all mappers QQQQ */
void
cdr_free_mappers (mapper_t **mappers)
{
  int s;

  for (s = 0; s < no_species; s++) {
    free (mappers[s]);
  }
  free (mappers);
}  

/** @brief dens_copy Copies QQQQ */
void
dens_copy (mapper_t *mapper, grid_t *source, grid_t *target, 
	   int ir, int iz, int itheta)
{
  cdr_grid_t *cdr_s, *cdr_d;

  cdr_s = (cdr_grid_t*) source;
  cdr_d = (cdr_grid_t*) target;

  RZT (cdr_d->dens[mapper->extra], ir, iz, itheta) = 
    RZT (cdr_s->dens[mapper->extra], ir, iz, itheta);
}

/** @brief Interpolates a set QQQQ */
int
dens_interpol_set (mapper_t *mapper, grid_t *source, interpol_t *interpol,
		   int pr, int pz, int itheta)
{
  cdr_grid_t *cdr;

  cdr = (cdr_grid_t*) source;

  interpol_set_stencil_at (source, interpol, 
			   r_at (pr, cdr->level),
			   z_at (pz, cdr->level),
			   cdr->dens[mapper->extra], pr, pz, itheta);
  return TRUE;
}

/** @brief dens_interpol Interpolates QQQQ */
void
dens_interpol (mapper_t *mapper, grid_t *source, grid_t *target, 
	       interpol_t *interpol, 
	       int ir, int iz, int itheta)
{
  double r, z;
  cdr_grid_t *cdr;

  cdr = (cdr_grid_t *) target;

  r = r_at (ir, cdr->level);
  z = z_at (iz, cdr->level);

  RZT (cdr->dens[mapper->extra], ir, iz, itheta) = 
    interpol_apply (interpol, r, z);
}

/** @brief Restricts the values of the densities from the child grid into the 
 * parent.
 */
static void
restrict_from (cdr_grid_t *parent, cdr_grid_t *child)
{
  int ir, iz, cr, cz, itheta, s;
  REAL sum;

  debug (2, "restrict_from (parent =" grid_printf_str 
	 ", child = " grid_printf_str ")\n",
	 grid_printf_args(parent), grid_printf_args(child));

  assert (parent->level == child->level - 1);

#pragma omp parallel private(sum, ir, iz, cr, cz, s)
  {
#pragma omp for
    iter_grid_theta (parent, itheta) {
      iter_grid_parent (child, ir, iz) {
	for (s = 0; s < no_species; s++) {
	  sum = 0;
	  for (cr = (ir << 1); cr <= (ir << 1) + 1; cr++)
	    for (cz = (iz << 1); cz <= (iz << 1) + 1; cz++)
	      sum += cyl_r_at (cr, child->level) 
		* RZT (child->dens[s], cr, cz, itheta);
	  RZT (parent->dens[s], ir, iz, itheta) = 
	    0.25 * sum / cyl_r_at (ir, parent->level);
	}
      }
    }
  }
}

/** @brief Restricts a grid with all its children.
 */
void
cdr_restrict (cdr_grid_t *grid)
{
  cdr_grid_t *child;

  iter_childs(grid, child) {
    restrict_from (grid, child);
  }
}

/** @brief And its recursive version.
 *
 * Note that this has to be tail recursive, sice we first have to restrict
 * the childs.
 */
mk_tail_recursive (cdr_restrict, cdr_grid_t)


/** @brief Now the routines to initialize a grid with given densities. 
 * gauss2 initializes one or two gaussian seeds
 */
double
gauss2_xyz (double x, double y, double z)
{
  double q;

  q =  invpi32 * 1 / (seed_index[curr_seed]->sigma_x * seed_index[curr_seed]->sigma_y *
		      seed_index[curr_seed]->sigma_z)
    * exp (- SQ(x - seed_index[curr_seed]->x0) / SQ(seed_index[curr_seed]->sigma_x)
	   - SQ(y - seed_index[curr_seed]->y0) / SQ(seed_index[curr_seed]->sigma_y)
	   - SQ(z - seed_index[curr_seed]->z0) / SQ(seed_index[curr_seed]->sigma_z));
  return q;
}

/** @brief Sets the densities of species species according to the function
 *  \f$f (x, y, z)\f$. 
 *
 *  mode can be one of the following:\n
 *  SET_DENS_OVERWRITE:  Set the dens to f(x, y, x), ignoring what was there 
 *                       before.\n
 *  SET_DENS_ADD:        Adds f to the former density.\n
 *  SET_DENS_SUB:        Substracts f from the former density.\n
 */
void 
cdr_set_dens (cdr_grid_t *cdr, int species, int mode, double factor,
	       double (*f) (double, double, double))
{
  int ir, iz, itheta;
  double ctheta, stheta, x, y;
  debug (2, "cdr_set_dens (" grid_printf_str ", %s, ...)\n",
	 grid_printf_args(cdr), spec_index[species]->name);

#pragma omp parallel private(ir, iz, ctheta, stheta, x, y)
  {
#pragma omp for
    iter_grid_theta_n (cdr, itheta, 2) {
      //sincos (theta_at(itheta), &stheta, &ctheta);
      stheta=sin(theta_at(itheta));
      ctheta=cos(theta_at(itheta));
      iter_grid_r_n(cdr, ir, 2) {
	x = r_at(ir, cdr->level) * ctheta;
	y = r_at(ir, cdr->level) * stheta;
	iter_grid_z_n(cdr, iz, 2) {
	  switch (mode) {
	  case SET_DENS_OVERWRITE:
	    *RZTP(cdr->dens[species], ir, iz, itheta) = 
	      factor * f (x, y, z_at(iz, cdr->level));
	    break;
	  case SET_DENS_ADD:
	    *RZTP(cdr->dens[species], ir, iz, itheta) += 
	      factor * f (x, y, z_at(iz, cdr->level));
	    break;
	  case SET_DENS_SUB:
	    *RZTP(cdr->dens[species], ir, iz, itheta) -= 
	      factor * f (x, y, z_at(iz, cdr->level));
	    break;
	  }
	}
      }
    }
  }
}

/** @brief Returns 1.0 */
double f_one(double x, double y, double z)
{
  return 1.0;
}

/** @brief Inits the densities of all species. */
void
cdr_init_dens (cdr_grid_t *cdr)
{
  int cnt;

  debug (2, "cdr_init_dens (" grid_printf_str ")\n",
	 grid_printf_args(cdr));

  for (cnt = 0; cnt < no_species; cnt++)
  {
    cdr_set_dens(cdr, cnt, SET_DENS_OVERWRITE, 0.0, f_one);
  }

  for (cnt = 0; cnt < no_seed; cnt++)
  {
    curr_seed = cnt;
    if (seed_index[cnt]->type == 0)
      cdr_set_dens(cdr, seed_index[cnt]->species, SET_DENS_ADD, seed_index[cnt]->value, gauss2_xyz);
    if (seed_index[cnt]->type == 1)
      cdr_set_dens(cdr, seed_index[cnt]->species, SET_DENS_ADD, seed_index[cnt]->value, f_one);
  }
}

/** @brief Initializes the densities allocating a new grid family.*/
cdr_grid_t*
cdr_scratch_init (void)
{
  cdr_grid_t *cdr;

  cdr = cdr_new_3d_a (0, 0, gridpoints_r, gridpoints_z, max_ntheta);

  cdr->level = 0;
  cdr->ext_bound = BND_MASK_ALL;
  cdr->contains_edge = TRUE;

  cdr_init_dens (cdr);
  cdr_refine_r (cdr, NULL);

  if (perturb_epsilon > 0.0 && max_ntheta > 1) {
    dft_dens_perturb_r (cdr, electrons, NULL);

    /* Maybe the perturbation has messed up the boundary conditions, so we
     * have to repair them.
     */
    cdr_set_ext_bnd_r (cdr);
    cdr_set_periodic_r (cdr);
  }
  return cdr;
}

/** @brief Dumps all the contents of the given grid into filenames given by
 * prefix and name
 */
void
cdr_dump (cdr_grid_t *grid, const char *prefix, const char *name)
{
  char *fname;
  int s, nt;

  int m = cdr_output_margin;

  /* To make it easier to produce plots in theta, we save one extra
   *  theta that has to give the same data as the first one.  Except if
   *  we are working in 2D, where anyway we have only one possible theta (0)
   */
  nt = grid->ntheta == 1? 1: grid->ntheta + 1;


#pragma omp parallel sections private(fname)
  {
#pragma omp section
    {
      asprintf (&fname, "%s/r.%s.tsv", prefix, name);
      rz_axis_dump (fname, grid->r0 - m, grid->r1 + m, dr[grid->level]);
      free (fname);
    }
#pragma omp section
    {
      asprintf (&fname, "%s/z.%s.tsv", prefix, name);
      rz_axis_dump (fname, grid->z0 - m , grid->z1 + m, dz[grid->level]);
      free (fname);
    }
#pragma omp section
    {
      /* In all this function we use ntheta + 1 to make matlab's life easier.
       * Note that some of these variables have not been made periodic
       * and hence the value for ntheta may be undefined. 
       */
      asprintf (&fname, "%s/theta.%s.tsv", prefix, name);
      rz_axis_dump (fname, 0, nt, dtheta);
      free (fname);
    }
#pragma omp section
    {
      asprintf (&fname, "%s/charge.%s.tsv", prefix, name);
      rz_dump_3d (grid->charge, fname, "w", grid->r0 - m, grid->z0 - m, 
		  grid->r1 + m, grid->z1 + m, nt);
      free (fname);
    }
#pragma omp section
    {
      asprintf (&fname, "%s/er.%s.tsv", prefix, name);
      rz_dump_3d (grid->er, fname, "w", grid->r0 - m, grid->z0 - m, 
		  grid->r1 + m, grid->z1 + m, nt);
      free (fname);
    }
#pragma omp section
    {
      asprintf (&fname, "%s/ez.%s.tsv", prefix, name);
      rz_dump_3d (grid->ez, fname, "w", grid->r0 - m, grid->z0 - m, 
		  grid->r1 + m, grid->z1 + m, nt);
      free (fname);
    }
#pragma omp section
    {
      asprintf (&fname, "%s/etheta.%s.tsv", prefix, name);
      rz_dump_3d (grid->etheta, fname, "w", grid->r0 - m, grid->z0 - m, 
		  grid->r1 + m, grid->z1 + m, nt);
      free (fname);
    }

#pragma omp section
    {
      asprintf (&fname, "%s/eabs.%s.tsv", prefix, name);
      rz_dump_3d (grid->eabs, fname, "w", grid->r0 - m, grid->z0 - m, 
		  grid->r1 + m, grid->z1 + m, nt);
      free (fname);
    }

#pragma omp section
    {
      if (has_photoionization) {
	asprintf (&fname, "%s/photo.%s.tsv", prefix, name);
	rz_dump_3d (grid->photo, fname, "w", grid->r0 - m, grid->z0 - m, 
		    grid->r1 + m, grid->z1 + m, nt);
	free (fname);
      }
    }

#pragma omp section
    {
      for (s = 0; s < no_species; s++) {
	/* The densities */
	asprintf (&fname, "%s/%s.%s.tsv", prefix, spec_index[s]->name, name);
	rz_dump_3d (grid->dens[s], fname, "w", grid->r0 - m, grid->z0 - m, 
		    grid->r1 + m, grid->z1 + m, nt);
	free (fname);

	/* ...and their derivatives. */
	asprintf (&fname, "%s/d_%s.%s.tsv", prefix, spec_index[s]->name, name);
	rz_dump_3d (grid->d_dens[s], fname, "w", grid->r0 - m, grid->z0 - m, 
		    grid->r1 + m, grid->z1 + m, nt);
	free (fname);
      }
    }
  }
}

/** @brief Dumps all the contents of the given grid into filenames given by
 * @a prefix and @a name.
 *
 * It also writes in tree.NAME.dat the tree structure of grid, in a format
 * appropiate for cdr_load_tree to understand.

   @a infp has to be NULL when called from outside.

   The format of a tree.NAME.dat is as follows.  Each line is composed
   by a "command" and some parameters. the "commands" can be

   * time
     parameter: sim_time
     Writes the current simulation time, as received in sim_time

   * open 
     parameters: gridname r0 z0 r1 z1 level ext_bound margin

     Reads the densities from the files given by gridname, with
     the corresponding values of the grid parameters.  Leaves the grid
     open to add childrens to it.  Those children will be all the grids
     read until the grid is closed.

   * close
     parameters: gridname

     Closes the latest opened grid.  All subsequent grids will be regarded
     as brothers of gridname.
 */
void
cdr_dump_r (cdr_grid_t *grid, const char *prefix, const char *name, 
	    FILE *infp, double sim_time)
{
  cdr_grid_t *child;
  char *cname;
  int i, nchilds;
  char codes[] = "abcdefghijklmnopqrstuvwxyz";
  FILE *fp;

  if (NULL == infp) {
    /* Root call */
    asprintf (&cname, "%s/tree.%s.dat", prefix, name);
    fp = fopen (cname, "w");
    if (NULL == fp) {
      fatal ("Could not open file %s/tree.%s.dat to write\n", 
	     prefix, name);
    }
    free (cname);
    fprintf (fp, "time %g\n", sim_time);
  } else {
    fp = infp;
  }

  /* We want to be sure that the charge that we plot is actually the
   * charge and not its Fourier transform.
   */
  cdr_calc_charge_r (grid);

  fprintf (fp, "open %s %d %d %d %d %d %d %d\n", name, grid->r0, grid->z0, 
	   grid->r1, grid->z1, grid->level, grid->ext_bound,
	   cdr_output_margin);

  cdr_dump (grid, prefix, name);

  nchilds = grid_howmany_children ((grid_t *) grid);

  for (i = 0; i < nchilds; i++) {
    child = (cdr_grid_t *) grid_get_child ((grid_t*) grid, nchilds - i - 1);

    assert (NULL != child);

    if (i > sizeof(codes) - 2)
      asprintf (&cname, "%s{%d}", name, i);
    else
      asprintf (&cname, "%s%c", name, codes[i]);

    cdr_dump_r (child, prefix, cname, fp, sim_time);
    free (cname);
  }

  fprintf (fp, "close %s\n", name);

  if (NULL == infp) {
    /* Root call */
    fclose (fp);
  }
}

/** @brief Reads the data format produced by cdr_dump_r and creates a complete
 * tree from it.
 *
 * This function is really quick-and-dirty: one should check the fscanfs
 * and do everything cleanlier. 
 */
cdr_grid_t *
cdr_load_tree_r (const char *prefix, const char *name, FILE *infp)
{
  char command[16], gridname[32], *fname;
  int r0, z0, r1, z1, level, ext_bound, margin, nt, s;
  int open_close;
  cdr_grid_t *grid = NULL, *leaf;
  FILE *fp;

  debug (2, "cdr_load_tree_r (prefix = \"%s\", name = \"%s\", ...)\n",
	 prefix, name);

  nt = max_ntheta == 1? 1: max_ntheta + 1;

  if (NULL == infp) {
    asprintf (&fname, "%s/tree.%s.dat", prefix, name);
    fp = fopen (fname, "r");
    if (NULL == fp) {
      fatal ("Could not open file %s/tree.%s.dat to read\n", 
	     prefix, name);
    }
    free (fname);
  } else {
    fp = infp;
  }

  do {
    fscanf (fp, "%15s %15s", command, gridname);
    open_close = TRUE;

    /* We analyze some commands that do not create/close grids, now
     * restricted to set time.
     */
    if (0 == strcmp (command, "time")) {
      double new_time;
      sscanf (gridname, "%lf\n", &new_time);
      warning ("Starting time [%f] was read from %s/tree.%s.dat\n",
	       new_time, prefix, name);

      start_t = new_time;
      open_close = FALSE;
    }
  } while (!open_close);

  if (0 == strcmp (command, "open")) {
    debug (3, "opening %s\n", gridname);

    fscanf (fp, "%d %d %d %d %d %d %d\n", 
	    &r0, &z0, &r1, &z1, &level, &ext_bound, &margin);

    grid = cdr_new_3d_a (r0, z0, r1, z1, max_ntheta);
    grid->ext_bound = ext_bound;
    grid->level = level;

    for (s = 0; s < no_species; s++) {
      /* The densities */
      asprintf (&fname, "%s/%s.%s.tsv", prefix, spec_index[s]->name, gridname);
      debug (3, "Loading %s\n", fname);
      rz_dump_3d (grid->dens[s], fname, "r", 
		  grid->r0 - margin, grid->z0 - margin, 
		  grid->r1 + margin, grid->z1 + margin, nt);
      free (fname);
    }

    do {
      leaf = cdr_load_tree_r (prefix, name, fp);
      if (NULL != leaf) {
	assert (leaf->level == grid->level + 1);
	add_child (grid, leaf);
      }

    } while (NULL != leaf);
  } else if (0 == strcmp (command, "close")) {
    debug (3, "closing %s\n", gridname);
    grid = NULL;
  }

  if (NULL == infp) {
    fclose (fp);
  }
  return grid;
}

/** @brief Given a grid family, writes in the file named @a fname the
 * coordinates of the frames that define the grid and his descendants.
 */
void
cdr_dump_frames (cdr_grid_t *grid, const char *prefix, const char *name)
{
  FILE *fp;
  char *fname;

  asprintf (&fname, "%s/frames.%s.tsv", prefix, name);
  fp = fopen (fname, "w");
  free (fname);

  if (fp == NULL) {
    warning ("Unable to open %s\n", fname);
    return;
  }
  aux_dump_frames_r (grid, fp);

  fclose (fp);
}

/** @brief aux_dump_frames_r QQQQ */
void
aux_dump_frames_r (cdr_grid_t *grid, FILE *fp)
{
  cdr_grid_t *child;
  int level;

  level = grid->level;

  fprintf (fp, "%g %g %g %g %d\n", 
	   er_r_at (grid->r0 - 1, level), ez_z_at (grid->z0 - 1, level),
	   er_r_at (grid->r1 - 1, level), ez_z_at (grid->z1 - 1, level), 
	   level);

  iter_childs(grid, child) {
    aux_dump_frames_r (child, fp);
  }
}
