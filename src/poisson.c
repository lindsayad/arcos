/** @file poisson.c
 *  @brief Poisson/Helmholtz solver, including the routines for
 *  inhomogeneous (point-plane) configurations.
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "cdr.h"
#include "cstream.h"
#include "fishpack.h"
#include "grid.h"
#include "interpol2.h"
#include "mapper.h"
#include "parameters.h"
#include "poisson.h"
#include "proto.h"
#include "rz_array.h"
#include "species.h"

static void get_bound_cond (pois_grid_t *grid, pois_problem_t *prob,
			    int *r_bound_cond, int *z_bound_cond, 
			    double lambda);
static int needs_refinement (pois_grid_t *grid, int ir, int iz, 
			     double threshold);
static int refine_in (pois_grid_t *grid, int cr0, int cz0, int cr1, int cz1);

decl_mapper_funcs(er);
decl_mapper_funcs(ez);
decl_mapper_funcs(etheta);
decl_mapper_funcs(charge);
/*!< Here we set three mapper_t objects that specify how we should do
 * the mapping between the potential (once it is calculated) and the three
 * components of the electric field.
 */

/* Global variables for the computations of inhomogeneous fields:
 *  initialized in pois_inhom_init. */

int pois_gridpoints_z, pois_gridpoints_r;
/*!< Number of gridpoints at level 0 of the poisson solver
 *
 * if !pois_inhom, this has to be equal to the parameter gridpoins_z.
 * We also allow a different extension in r, albeit this is not yet used.
 */

double pois_inhom_L;    /*!< Total length of the poisson domain. */
double pois_inhom_z;    /*!< Location of the floating charge. */
double pois_inhom_dphi;
/*!< Potential difference between the electrodes created by a unit charge 
 *  located at (0, pois_inhom_z)
 */

double pois_inhom_fixed_q_t;
/*!< Time-dependent fixed charge for the inhomogeneous poisson code
 *  It is updated in cstream_set_field_at_time acording to the parameter
 *  rise_time.
 */

extern double E_x, E_y, E_z;

static double inhom_phi_term (double r, double z, double b);
static double inhom_e_term (double r, double z, double b);
static double inhom_er_term (double r, double z, double b);
static double inhom_ez_term (double r, double z, double b);
/*!< To calculate the potential and electric fields created by a unit charge
 *  between two planar electrodes, we use the "mirror charges method"
 *  the contribution of each mirror charge is given by these functions.
 */

mapper_t er_mapper = mk_mapper_staggered (er, &interpol_bilin, 0, -1);
mapper_t ez_mapper = mk_mapper_staggered (ez, &interpol_bilin, -1, 0);
mapper_t etheta_mapper = mk_mapper (etheta, &interpol_quadratic);
mapper_t *e_mappers[] = {&er_mapper, &ez_mapper, &etheta_mapper, NULL};
/*!< These are the mapper structures to translate between cdr and poisson
 *  grids.  See mapper.c and mapper.h for more info.
 */

mapper_t charge_mapper = mk_mapper_down(charge, &interpol_wackers);
mapper_t *charge_mappers[] = {&charge_mapper, NULL};

pois_problem_t *pois_electrostatic;

/** @brief  All initialization of the Poisson solver has to be here. */
void
pois_init (void)
{
  if (pois_inhom) {
    pois_inhom_init ();
  }
  else {
    pois_gridpoints_z = gridpoints_z;
  }
  pois_gridpoints_r = gridpoints_r;

  pois_electrostatic = xmalloc (sizeof (pois_problem_t));

  pois_electrostatic->max_level = pois_max_level;
  pois_electrostatic->extra_levels = extra_pois_levels;
  pois_electrostatic->max_error = pois_max_error;

  pois_electrostatic->bnd_right = pois_bnd_right;
  pois_electrostatic->bnd_top = pois_bnd_top;
  pois_electrostatic->bnd_bottom = pois_bnd_bottom;
}

/** @brief Creates a 2D Poisson grid. */
pois_grid_t*
pois_new_a (int r0, int z0, int r1, int z1)
{
  pois_grid_t *grid;

  grid = pois_new_3d_a (r0, z0, r1, z1, 1);

  return grid;
}

/** @brief Creates a new 3D Poisson grid.
 *
 * @a r0, @a z0, @a r1, @a z1 are in GRID units
 */
pois_grid_t*
pois_new_3d_a (int r0, int z0, int r1, int z1, int ntheta)
{
  rz_array_t *phi, *charge, *error;
  pois_grid_t *grid;

  debug (2, "pois_new_a (r0 = %d, z0 = %d, r1 = %d, z1 = %d)\n", 
	 r0, z0, r1, z1);

  phi = rz_new_3d_a (r0, z0, r1, z1, ntheta);
  error = rz_new_3d_a (r0, z0, r1, z1, ntheta);

#ifdef F_PHI_SAME_AS_CHARGE
  charge = phi;
#else
  charge = rz_new_3d_a (r0, z0, r1, z1, ntheta);
#endif /* F_PHI_SAME_AS_CHARGE */

  grid = (pois_grid_t *) xmalloc (sizeof (pois_grid_t));

  grid->r0 = r0;
  grid->r1 = r1;
  grid->z0 = z0;
  grid->z1 = z1;
  grid->ntheta = ntheta;

  grid->phi = phi;
  grid->charge = charge;
  grid->error = error;

  /* By default, the array does not have any external boundary. */
  grid->ext_bound = BND_NONE;

  /* Initially, a grid does not have any relatives. */
  init_leaf (grid);

  return grid;
}

/** @brief  Frees the memory allocated by pois_new_a */
void
pois_free (pois_grid_t *grid)
{
  debug (3, "pois_free (...)\n");

  /* The decission of whether to use the same memory for the charge and
   * the potential can be made elsewhere.
   */
  if (grid->phi != grid->charge)
    rz_free (grid->charge);
  rz_free (grid->phi);
  rz_free (grid->error);

  free (grid);
}

/** @brief Recursively frees a tree of Poisson grids. */
void
pois_free_r (pois_grid_t *grid)
{
  pois_grid_t *leaf;

  debug (3, "pois_free_r (" grid_printf_str ")\n", grid_printf_args (grid));

  free_childs (grid, leaf, pois_free_r);

  pois_free (grid);
}

/** @brief Creates a new Poisson grid.
 *
 * @a r0, @a r1... are in global units (i.e.  independent of the grid level). 
 */
pois_grid_t*
pois_new_glob_a (int r0, int z0, int r1, int z1, int level)
{
  debug (3, "pois_new_glob_a (%d, %d, %d, %d)\n", r0, r1, z0, z1);
  pois_grid_t *grid;

  if (level >= 0) {
    grid = pois_new_a (r0 << level, z0 << level, 
		       r1 << level, z1 << level);
  } else {
    level = -level;
    grid = pois_new_a (r0 >> level, z0 >> level, 
		       r1 >> level, z1 >> level);
  }
  return grid;
}

/** @brief Starts the tree of Poisson grids with the two coarsest ones.
 *
 * Receives the grid dimensions at level 0.
 */
pois_grid_t*
pois_init_tree_a (int r0, int z0, int r1, int z1)
{
  pois_grid_t *coarsest, *sub_coarsest;
  int level;

  debug (3, "pois_init_tree_a (%d, %d, %d, %d)\n", r0, z0, r1, z1);

  level = -extra_pois_levels;

  coarsest = pois_new_glob_a (r0, z0, r1, z1, level);
  coarsest->ext_bound = BND_MASK_ALL;
  coarsest->level = level;

  sub_coarsest = pois_new_glob_a (r0, z0, r1, z1, level + 1);

  add_child (coarsest, sub_coarsest);
  grid_inherit_ext_bound ((grid_t *) sub_coarsest);

  return coarsest;
}

/** @brief Is it possible to use a single function to handle all four
 *  boundaries?
 *
 * Maybe this should be rewritten to use the interpol2 module? 
 */
REAL*
pois_boundary_a (pois_grid_t *grid, int boundary)
{
  pois_grid_t *parent;

  /* Forget the 0s, they are there to shut up the compiler. */
  int stride = 0, ortho_stride = 0, icells, invert, n;
  int rb, zb;
  double interp_array[3][3] = {{-4.0, 1.0, 18.0},
			       {-2.0, 78.0, 14.0},
  			       {-6.0, -7.0, 4.0}};

  const double prefactor = 1 / 96.0;
  REAL *bnd, *p, *lstart;

  debug (3, "pois_boundary_a (..., %d)\n", boundary);

  parent = grid->parent;

  /* First we check if the boundary is in r or in z. */
  if (boundary & BND_AT_Z) {
    icells = grid->r1 - grid->r0;
    zb = 1; rb = 0;
    if (NULL != parent) {
      stride = parent->phi->strides[R_INDX];
      ortho_stride = parent->phi->strides[Z_INDX];
    }
  } else {
    zb = 0; rb = 1;
    icells = grid->z1 - grid->z0;
    if (NULL != parent) {
      stride = parent->phi->strides[Z_INDX];
      ortho_stride = parent->phi->strides[R_INDX];
    }
  }

  debug (3, "   grid->ext_bound = %d\n", grid->ext_bound);

  if (grid->ext_bound & BND_MASK (boundary)) {
    /* If the boundary is external, we apply homogeneous Dirichlet/Neumann. */
    bnd = (REAL *) xcalloc (sizeof (REAL), icells);
    debug (3, "pois_boundary_a (...) -> [0.0] * %d\n", icells);

    assert (grid->r1 != 0);
    return bnd;
  } else {
    bnd = (REAL *) xmalloc (sizeof (REAL) * icells);
  }
  /* If the grid is root it should have returned in the previous if */
  assert (NULL != parent);

  /* Now we check if the boundary is at {r,z} minimum or {r,z} maximum */
  if (boundary & BND_MAX) {
    lstart = RZP (parent->phi, 
		(grid->r1 >> 1) - 1 + rb, 
		(grid->z1 >> 1) - 1 + zb);
    stride = -stride;
    ortho_stride = -ortho_stride;
    invert = TRUE;
  } else {
    lstart = RZP (parent->phi, (grid->r0 >> 1) - rb, (grid->z0 >> 1) - zb);
    invert = FALSE;
  }

  for (p = lstart, n = 0; n < icells; p += stride, n += 2) {
    /* p points to the cell in the parent grid and p + stride points to 
       an adjacent cell. */
    double s = 0.0, s2 = 0.0;
    int i, j;
    for (i = -1; i <= 1; i++) {
      for (j = -1; j <= 1; j++) {
	s += *(p + i * stride + j * ortho_stride) * interp_array[i + 1][j + 1];
	s2 += *(p - i * stride + j * ortho_stride) * interp_array[i + 1][j + 1];
      }
    }
    if (!invert) {
      bnd[n] = prefactor * s;
      if (n < icells - 1) bnd[n + 1] = prefactor * s2;
    } else {
      bnd[icells - n - 1] = prefactor * s;
      if (n < icells - 1) bnd[icells - n - 2] = prefactor * s2;
    }
  }
  return bnd;
}

/** @brief Solves the Poisson equation, calling the FISHPACK routine.
 *
 * Assumes that grid->charge is already set.
 */
void
pois_solve_grid (pois_grid_t *grid, pois_problem_t *prob, 
		 double lambda, double s) 
{
  REAL *boundaries[4];
  int i, z_bound_cond, r_bound_cond;
  double rmin, rmax, zmin, zmax;

  debug (3, "pois_solve_grid (..., lambda = %f, s = %f)\n", lambda, s);

  for (i = 0; i < 4; i++) {
    boundaries[i] = pois_boundary_a (grid, i);
    debug (3, "boundaries[%d] = %p\n", i, boundaries[i]);
  }

#ifndef F_PHI_SAME_AS_CHARGE
  rz_copy (grid->charge, grid->r0, grid->z0,
	   grid->phi, grid->r0, grid->z0,
	   grid->r1 - grid->r0, grid->z1 - grid->z0);
#endif

  rmin = grid->r0 * dr[grid->level];
  rmax = grid->r1 * dr[grid->level];
  zmin = grid->z0 * dz[grid->level];
  zmax = grid->z1 * dz[grid->level];

  debug (3, "HSTCYL/HSTCRT in a grid %d x %d\n", grid->r1 - grid->r0,
	 grid->z1 - grid->z0);

  get_bound_cond (grid, prob, &r_bound_cond, &z_bound_cond, lambda);

#ifndef TRUE2D
  fish_hstcyl (rmin, rmax, grid->r1 - grid->r0, 
	       r_bound_cond, boundaries[BND_LEFT], boundaries[BND_RIGHT],
	       zmin, zmax, grid->z1 - grid->z0,
	       z_bound_cond, boundaries[BND_BOTTOM], boundaries[BND_TOP],
	       lambda, s,
	       RZP (grid->phi, grid->r0, grid->z0), 
	       grid->phi->strides[Z_INDX]);
#else
  fish_hstcrt (rmin, rmax, grid->r1 - grid->r0, 
	       r_bound_cond, boundaries[BND_LEFT], boundaries[BND_RIGHT],
	       zmin, zmax, grid->z1 - grid->z0,
	       z_bound_cond, boundaries[BND_BOTTOM], boundaries[BND_TOP],
	       lambda,
	       RZP (grid->phi, grid->r0, grid->z0), 
	       grid->phi->strides[Z_INDX]);
#endif
  debug (3, "Finished HSTCYL/HSTCRT in a grid %d x %d\n", grid->r1 - grid->r0,
	 grid->z1 - grid->z0);

  pois_set_phi_boundaries (grid, boundaries, 
			   (grid->ext_bound & BND_MASK (BND_LEFT))
			    && 0. == lambda,
			   
			   prob->bnd_right == BND_CND_HNEUMANN
			   && (grid->ext_bound & BND_MASK (BND_RIGHT)),
			   
			   prob->bnd_bottom == BND_CND_HNEUMANN
			   && (grid->ext_bound & BND_MASK (BND_BOTTOM)),
			   
			   prob->bnd_top == BND_CND_HNEUMANN
			   && (grid->ext_bound & BND_MASK (BND_TOP)));

  for (i = 0; i < 4; i++) {
    debug (3, "boundaries[%d] = %p\n", i, boundaries[i]);
    free (boundaries[i]);
  }
  debug (3, "  <- pos_solve_grid (...)\n");
}

/** @brief Gets the boundary conditions that we pass to FISHPACK.
 *
 * Note that we use the same conditions for the Poisson equation
 * and for the Helmholtz equation used to find the photoionization
 * source.  From my point of view, it doesn't make sense to use different
 * conditions, so that will simply add an unneccesary complication.
 */
static void
get_bound_cond (pois_grid_t *grid, pois_problem_t *prob,
		int *r_bound_cond, int *z_bound_cond, double lambda)
{
  static int right_inside[] = {FISH_UNS_DIR, FISH_NEU_DIR, FISH_DIR_DIR};
  static int right_ext_neu[] = {FISH_UNS_NEU, FISH_NEU_NEU, FISH_DIR_NEU};
  static int *right_ext_dir = right_inside;
  static int fish_order[4] = { FISH_DIR_DIR, FISH_DIR_NEU, 
			       FISH_NEU_DIR, FISH_NEU_NEU};

  int bottom, top;
  int *right_sel;

  if ((grid->ext_bound & BND_MASK (BND_TOP))
      && prob->bnd_top == BND_CND_HNEUMANN)
    top = 1;
  else top = 0;

  if ((grid->ext_bound & BND_MASK (BND_BOTTOM))
      && prob->bnd_bottom == BND_CND_HNEUMANN)
    bottom = 1;
  else bottom = 0;

  *z_bound_cond = fish_order[(bottom << 1) + top];

  if ((grid->ext_bound & BND_MASK (BND_RIGHT))) {
    if (prob->bnd_right == BND_CND_HNEUMANN)
      right_sel = right_ext_neu;
    else 
      right_sel = right_ext_dir;
  } else {
    right_sel = right_inside;
  }

#ifndef TRUE2D
  *r_bound_cond = ((grid->ext_bound & BND_MASK (BND_LEFT)) && 0. == lambda)? 
		   right_sel[0]: right_sel[2];
#else
  *r_bound_cond = ((grid->ext_bound & BND_MASK (BND_LEFT))? 
		   right_sel[1]: right_sel[2]);
#endif

}

/** @brief Uses the @a boundaries[] vectors to set the boundaries by
 *  interpolation for the phi array of a grid.
 *
 *  (in the extra space we allocated when we created the array).
 *  xxx_neu tells us if the xxx boundary has Neumann b.c.
 *  (or unspecified, as FISHPACK calls them if they are applied in the axis)
 *
 *  TODO: Please rewrite this mess.
 */
void
pois_set_phi_boundaries (pois_grid_t *grid, REAL *boundaries[], 
			 int left_neu, int right_neu, int bottom_neu, 
			 int top_neu)
{
  int i;

  debug (3, "pois_set_phi_boundaries (...)\n");

  if ((grid->ext_bound & BND_MASK(BND_BOTTOM)) && bottom_neu) {
    for (i = grid->r0; i < grid->r1; i++) {
      RZ (grid->phi, i, grid->z0 - 1) = RZ (grid->phi, i, grid->z0);
    }
  } else {
    for (i = grid->r0; i < grid->r1; i++) {
      RZ (grid->phi, i, grid->z0 - 1) = 
	2 * boundaries[BND_BOTTOM][i - grid->r0] - RZ (grid->phi, i, grid->z0);
    }
  }

  if ((grid->ext_bound & BND_MASK(BND_TOP)) && top_neu) {
    for (i = grid->r0; i < grid->r1; i++) {
      RZ (grid->phi, i, grid->z1) = RZ (grid->phi, i, grid->z1 - 1);
    }
  } else {
    for (i = grid->r0; i < grid->r1; i++) {
      RZ (grid->phi, i, grid->z1) = 2 * boundaries[BND_TOP][i - grid->r0] 
	- RZ(grid->phi, i, grid->z1 - 1);
    }
  }

  if (left_neu) {
    for (i = grid->z0; i < grid->z1; i++) {
      RZ (grid->phi, grid->r0 - 1, i) = RZ(grid->phi, grid->r0, i);
    }
  } else {
    for (i = grid->z0; i < grid->z1; i++) {
      RZ(grid->phi, grid->r0 - 1, i) = 2 * boundaries[BND_LEFT][i - grid->z0] 
	- RZ(grid->phi, grid->r0, i);
    }
  }

  if (right_neu) {
    for (i = grid->z0; i < grid->z1; i++) {
      RZ (grid->phi, grid->r1, i) = RZ(grid->phi, grid->r1 - 1, i);
    }
  } else {
    for (i = grid->z0; i < grid->z1; i++) {      
      RZ(grid->phi, grid->r1, i) = 2 * boundaries[BND_RIGHT][i - grid->z0] 
	- RZ(grid->phi, grid->r1 - 1, i);
    }
  }

  /*  This is quite nightmarish.  Think about rewritting. */
  if ((grid->ext_bound & BND_MASK(BND_RIGHT)) && right_neu) {
    RZ(grid->phi, grid->r1, grid->z0 - 1) = RZ(grid->phi, grid->r1 - 1, 
						   grid->z0 - 1);
    RZ(grid->phi, grid->r1, grid->z1) = RZ(grid->phi, grid->r1 - 1, 
					       grid->z1);
  } else {
    RZ(grid->phi, grid->r1, grid->z0 - 1) =
      RZ(grid->phi, grid->r1 - 1, grid->z0 - 1) 
      + RZ(grid->phi, grid->r1, grid->z0)
      - 0.5 * (RZ(grid->phi, grid->r1 - 2, grid->z0 - 1) 
	       + RZ(grid->phi, grid->r1, grid->z0 + 1));

    RZ(grid->phi, grid->r1, grid->z1) =
      RZ(grid->phi, grid->r1 - 1, grid->z1) 
      + RZ(grid->phi, grid->r1, grid->z1 - 1)
      - 0.5 * (RZ(grid->phi, grid->r1 - 2, grid->z1) 
	       + RZ(grid->phi, grid->r1, grid->z1 - 2));
  }

  if ((grid->ext_bound & BND_MASK(BND_LEFT)) && left_neu) {
    RZ(grid->phi, grid->r0 - 1, grid->z0 - 1) = RZ(grid->phi, grid->r0, 
						   grid->z0 - 1);
    RZ(grid->phi, grid->r0 - 1, grid->z1) = RZ(grid->phi, grid->r0, 
					       grid->z1);
  } else {
    RZ(grid->phi, grid->r0 - 1, grid->z0 - 1) =  
      RZ(grid->phi, grid->r0, grid->z0 - 1)
      + RZ(grid->phi, grid->r0 - 1, grid->z0)
      - 0.5 * (RZ(grid->phi, grid->r0 + 1, grid->z0 - 1) 
	       + RZ(grid->phi, grid->r0 - 1, grid->z0 + 1));

    RZ(grid->phi, grid->r0 - 1, grid->z1) =  
      RZ(grid->phi, grid->r0, grid->z1)
      + RZ(grid->phi, grid->r0 - 1, grid->z1 - 1)
      - 0.5 * (RZ(grid->phi, grid->r0 + 1, grid->z1) 
	       + RZ(grid->phi, grid->r0 - 1, grid->z1 - 2));
  }
}

/** @brief  Estimates the error on this grid by comparing with the calculations
 *  for his parent and performing a Richardson extrapolation.
 */
void
pois_set_error (pois_grid_t *grid)
{
  pois_grid_t *parent;
  interpol_t *interpol;
  int pr, pz;

  debug (3, "pois_set_error(...)\n");

  parent = grid->parent;
  assert (NULL != parent);

  interpol = interpol_new_a (2.0, 2.0, &interpol_wackers);

  iter_grid_parent (grid, pr, pz) {
    int i, j;

    interpol_set_stencil_at ((grid_t *) grid, interpol, 
			     pr * 2.0 + 1, pz * 2.0 + 1,
			     parent->phi, pr, pz, 0);
    for (i = 0; i < 2; i++)
      for (j = 0; j < 2; j++) {
	int cr, cz;
	double v;

	cr = (pr << 1) + i;
	cz = (pz << 1) + j;
	
	v = interpol_apply (interpol, (double) cr + 0.5, (double) cz + 0.5);

	/* The error will only be used as its abs value.  Therefore
	   we calculate it here and forever. */
	RZ (grid->error, cr, cz) = fabs (RZ (grid->phi, cr, cz) - v);
	
      }
  }
  interpol_free (interpol);
}

/** @brief The refinement routine.
 *
 * Threshold is initially pois_max_error, but if we get too large refinement
 * it can be increased (after warning the user that we are not reaching
 * the accuracy that he asked for).
 */
int
pois_refine (pois_grid_t *grid, double threshold)
{
  int ir, iz;
  int zmin, zmax, rmin, rmax;
  /* Be careful: cz1 and cr1 are inclusive: when the child is actually
     created, one has to increse them by 1. */
  int cr0, cr1, cz0, cz1;
  int all_ok = TRUE;

  int tainted_line = FALSE;

  debug (3, "pois_refine(...)\n");

  /* Initially we set values out of the grid, so the first comparisons
   * always change the values.
   */
  cr0 = grid->r1; cr1 = grid->r0 - 1;
  cz0 = grid->z1; cz1 = grid->z0 - 1;
 
  /* The boundary can only coincide with the parent's boundary if this
   * is an external boundary.
   */
  rmin = (grid->ext_bound & BND_MASK(BND_LEFT))? grid->r0: (grid->r0 + 1);
  rmax = (grid->ext_bound & BND_MASK(BND_RIGHT))? grid->r1: (grid->r1 - 1);
  zmin = (grid->ext_bound & BND_MASK(BND_BOTTOM))? grid->z0: (grid->z0 + 1);
  zmax = (grid->ext_bound & BND_MASK(BND_TOP))? grid->z1: (grid->z1 - 1);

  for (iz = zmin; iz < zmax && all_ok; iz++) {
    tainted_line = FALSE;
    for (ir = rmin; ir < rmax; ir++) {
      if (needs_refinement (grid, ir, iz, threshold)) {
	tainted_line = TRUE;
	if (ir < cr0) cr0 = ir;
	if (ir > cr1) cr1 = ir;
      }
    }
    if (tainted_line) {
      if (iz < cz0) cz0 = iz;
      if (iz > cz1) cz1 = iz;
    } else if (cr1 >= cr0 && cz1 >= cz0) {
      all_ok = all_ok && refine_in (grid, cr0, cz0, cr1, cz1);
      cr0 = grid->r1; cr1 = grid->r0 - 1;
      cz0 = grid->z1; cz1 = grid->z0 - 1;
    }
  }

  if (cr1 >= cr0 && cz1 >= cz0) {
    all_ok = all_ok && refine_in (grid, cr0, cz0, cr1, cz1);
  }

  /* If there was some error with some refinement, we delete the created
   * sub-grids and go back to the caller, which has to increase the refinement
   * threshold.
   */
  if (!all_ok) {
    pois_grid_t *leaf;
    free_childs (grid, leaf, pois_free_r);
    set_childless (grid);
  }
  return all_ok;
}

#define REFINE_IN_MAX_WARN_CNT 20;

/** @brief Refines and creates a new grid according to the coordinates given.
 *
 * (if the FISHPACK limit allows it). Returns TRUE if the refinement
 * succeeded, FALSE otherwise (if the FISHPACK limit was exceeded).
*/
static int
refine_in (pois_grid_t *grid, int cr0, int cz0, int cr1, int cz1)
{
  pois_grid_t *child;
  int nr0, nr1, nz0, nz1;
  static int warn_cnt = REFINE_IN_MAX_WARN_CNT;

  nr0 = cr0 << 1;
  nz0 =  cz0 << 1;
  nr1 = (++cr1) << 1;
  nz1 = (++cz1) << 1;

  if ((nr1 - nr0) > FISH_MAX_GRIDPOINTS ||  
      (nz1 - nz0) > FISH_MAX_GRIDPOINTS) {
    if (warn_cnt > 0) {
      warning ("FISHPACK limit exceeded.  Not refining grid "
	       grid_printf_str ".\n", nr0, nz0, nr1, nz1, grid->level + 1);
      warn_cnt --;
    }
    if (warn_cnt == 0) {
      warning ("No more warnings about FISHPACK limit.\n");
      warn_cnt --;
    }
    return FALSE;
  }

  child = pois_new_a (nr0, nz0, nr1, nz1);

  add_child (grid, child);
  grid_inherit_ext_bound ((grid_t *) child);
  debug (3, "new child created {r0 = %d, z0 = %d, r1 = %d, z1 = %d, "
	 "level = %d}\n", child->r0, child->z0, child->r1, child->z1, 
	 child->level);

  return TRUE;
}

/** @brief Calculates the error measure of this grid.
 *
 * For debug/testing purposes
 */
void
pois_error_measures (pois_grid_t *grid, double *L1, double *L2, double *Lmax)
{
  int ir, iz, i;
  double err;
  *L2 = 0.0;
  *L1 = 0.0;
  *Lmax = 0.0;
  i = 0;

  iter_grid (grid, ir, iz) {
    i++;
    err = RZ (grid->error, ir, iz);
    if (err > *Lmax) *Lmax = err;
    *L1 += err;
    *L2 += err * err;
  }

  *L1 /= i;
  *L2 = sqrt (*L2 / i);

}

/** @brief Writes some error measures of this grid and its descendants
 *  into @a fp
 */
void
pois_write_error_r (pois_grid_t *grid, FILE *fp)
{
  double L1, L2, Lmax;
  pois_grid_t *child;

  pois_error_measures (grid, &L1, &L2, &Lmax);

  fprintf (fp, "%d %g %g %g %g %g\n", grid->level, dr[grid->level],
	   dz[grid->level], L1, L2, Lmax);

  iter_childs (grid, child) {
    pois_write_error_r (child, fp);
  }
}

/** @brief Takes a @a cdr tree and solves the Poisson equation for it.
 *
 * Returns a poisson tree.
 */
pois_grid_t **
pois_solve_a (cdr_grid_t *cdr, pois_problem_t *prob)
{
  /* Call to the generic Poisson/Helmholtz solver */
  return pois_gen_solve_a (cdr, prob, e_mappers, 0.0);
}

/** @brief Solves a Poisson/Helmholtz equation with FISHPACK and maps
 *  the result using a given set of mapper.
 *
 *  From this method downwards (in the sense of calling order, not file
 *  organization), all the code is shared by the Poisson solver and the
 *  Photoionization (Helmholtz) solver.
 */
pois_grid_t **
pois_gen_solve_a (cdr_grid_t *cdr, pois_problem_t *prob,
		  mapper_t **mappers, double es)
{
  cdr_grid_t *tree;
  pois_grid_t **pois_modes;
  int mode;

  pois_modes = (pois_grid_t**) xmalloc (sizeof(pois_grid_t*) * max_ntheta);
  assert (0 == cdr->level);
  tree = cdr_add_coarser_grids_a (cdr, prob->extra_levels);

  /* This is the high-level parallelizable loop! */
#pragma omp parallel
  {
#pragma omp for schedule(dynamic)
    for (mode = 0; mode < max_ntheta; mode++) {
      pois_modes[mode] = pois_solve_mode (tree, cdr, prob, mode, es);

      if (pois_inhom && es == 0.0 && mode == 0) {
	/* Add the inhomogeneous laplacian potential to phi.
	 * It has only to be added to the zero mode.
	 */
	double q_factor;
	q_factor = pois_inhom_q_factor (pois_modes[0]);
	
	/* Was:
	   pois_add_inhom_phi_r (pois_modes[0], q_factor);
	*/
	pois_add_inhom_phi_r (pois_modes[0], -q_factor);
      }

      /* For the interpolation, we require an extra buffer cell
       * in some boundaries and we check for overflows inside the
       * mapper methods.
       */
      map_trees_r (mappers, (grid_t*) pois_modes[mode], (grid_t*) cdr, 
		   mode, FALSE, TRUE, FALSE, 
		   /*s_buf = */ 2, 
		   /*t_buf = */ 4);
      map_trees_r (mappers, (grid_t*) pois_modes[mode], (grid_t*) cdr, 
		   mode, TRUE, FALSE, TRUE, 
		   /*s_buf = */ 1, 
		   /*t_buf = */ 2);
    }
  }
  cdr_free_coarser_grids (tree, extra_pois_levels);

  return pois_modes;
}

/** @brief Solves a single Fourier mode.
 *
 * @a tree is the root of the @a cdr tree with extra_pois_levels added.
 * @a cdr is the cdr grid at level 0, contained in tree:
 *
 *   cdr = tree->first_child->firts_child...->first_child
 *            \_________ extra_pois_levels __________/
 */
pois_grid_t *
pois_solve_mode (cdr_grid_t *tree, cdr_grid_t *cdr, pois_problem_t *prob,
		 int mode, double es)
{
  pois_grid_t *pois;

  pois = pois_init_tree_a (0, 0, pois_gridpoints_r, pois_gridpoints_z);

  /* Solves the coarsest Poisson grid. */
  /* Changed Wed Mar 19 13:50:23 2008: was s_buf = 1. */
  map_trees_r (charge_mappers, (grid_t *) tree, (grid_t *) pois, 
	       mode, TRUE, TRUE, FALSE, 0, 2);
  pois_solve_grid (pois, prob, -w2k[mode], es);

  /* Wed May 14 11:00:30 2008:
   * I divide pois_max_error by abs(wk[mode]) because in the calculation
   * of e_theta, the errors get multiplied by wk[mode].  Hence if we
   * want to get similar errors in all modes, we have to remove it from the
   * threshold.  The 1 is because wk[0] = 0, but we do not want to have
   * a \inf threshold.
   */
  pois_solve_r (pois->first_child, tree, prob, mode, es, 
		prob->max_error / (1 + abs(wk[mode])));

  return pois;
}

#define POIS_SOLVE_R_MAX_WARN 50

/** @brief  Recursively solves one Poisson grid. */
void
pois_solve_r (pois_grid_t *pois, cdr_grid_t *cdr, pois_problem_t *prob,
	      int mode, double es, double threshold)
{
  pois_grid_t *child;
  int succeeded = FALSE;
  static int warn_cnt = POIS_SOLVE_R_MAX_WARN;

  /* Changed Wed Mar 19 13:51:00 2008: was s_buf = 1 */
  map_grid_r (charge_mappers, (grid_t *) cdr, (grid_t *) pois, 
	      mode, TRUE, TRUE, FALSE, 0, 2);
  debug(3, "w2k[mode = %d] = %f\n", mode, w2k[mode]);

  pois_solve_grid (pois, prob, -w2k[mode], es);

  if (pois->level >= prob->max_level) return;

  pois_set_error (pois);

  while (!succeeded) {
    succeeded = pois_refine (pois, threshold);

    if (!succeeded) {
      threshold *= 2;
      if (warn_cnt > 0) {
	warning ("Poisson error threshold was too small: trying %g.\
 Consider increasing pois_max_error\n",
		 threshold);
	if (warn_cnt >= 0) warn_cnt --;
      }
    }
  }
  iter_childs (pois, child) {
    pois_solve_r (child, cdr, prob, mode, es, threshold);
  }
}

/** @brief Returns 1 if the given cell has to be refined.
 *
 * which means that it satisfies at least one of the following:
 * a. The error in the cell is larger than the threshold.
 * b. The cell has a neightbour with an error larger than the threshold
 */
static int
needs_refinement (pois_grid_t *grid, int ir, int iz, double threshold)
{
  int i, j;
  int r = FALSE;
  int irmin, irmax, izmin, izmax;

  irmin = ir > grid->r0 + 3? ir - 1: grid->r0;
  irmax = ir < grid->r1 - 4? ir + 1: grid->r1 - 1;
  izmin = iz > grid->z0 + 3? iz - 1: grid->z0;
  izmax = iz < grid->z1 - 4? iz + 1: grid->z1 - 1;

  for (i = irmin; i <= irmax && !r; i++)
    for (j = izmin; j <= izmax; j++)
      if (RZ(grid->error, i, j) > threshold) {
	r = TRUE;
	break;
      }

  return r;
}

/*****************************************************************************
 * Functions for the implemetation of a needle-plate discharge.              *
 * The method used here is that of the "floating charge":                    *
 * A charge is located at pois_inhom_z, which is used to keep the potential  *
 * at L_z constant.                                                          *
 ****************************************************************************/

/** @brief pois_inhom_init QQQQ */
void
pois_inhom_init (void)
{
  /* The domain of solution of the Poisson equation will be different
   * from that of the CDR equation.  The difference between them is
   * given by the parameter needle_length, but has to be rounded
   * to a multiple of the grid size at level 0.
   */
  pois_gridpoints_z = (gridpoints_z
		       + (int) nearbyint (needle_length / dz[0]));
  pois_inhom_L = pois_gridpoints_z * dz[0];
		       
  pois_inhom_z = L_z + needle_radius;
  pois_inhom_dphi = pois_inhom_phi (0, pois_inhom_L) - pois_inhom_phi (0, L_z);

}

/* We use the following functions to calculate the fields and potential
 * created by a unit charge located at (r = 0, z = b).
 * These are created by the method of mirror charges, with a total of
 * pois_inhom_reflections charges.
 */

/* Single term of phi. */
static double
inhom_phi_term (double r, double z, double b)
{
  double insqrt;

  insqrt = (SQ(z - b) + SQ(r));

  return invfourpi / sqrt (insqrt);
}

/** @brief Common decay factor for the electric fields. */
static double
inhom_e_term (double r, double z, double b)
{
  double insqrt;

  insqrt = (SQ(z - b) + SQ(r));
  insqrt = insqrt * insqrt * insqrt;

  return invfourpi / sqrt (insqrt);
}

/** @brief Single term of E_r. */
static double
inhom_er_term (double r, double z, double b)
{
  return r * inhom_e_term (r, z, b);
}

/** @brief Single term of E_z */
static double
inhom_ez_term (double r, double z, double b)
{
  return (z - b) * inhom_e_term (r, z, b);
}

/* Since we have to sum over the same locations of the mirror charges
 * but with different "field functions", we use the following macro.
 *
 * Note that using a single function for inhom_phi, inhom_er, and inhom_ez
 * would not be a good idea since they are evaluated at different points.
 * Passing that function a pointer to another function would have been somewhat
 * slower and, I think, not so readable.
 */
#define sum_reflections(total_, field_, r_, z_)					\
       do {									\
	 int i__;								\
	 total_ = (field_ (r_, z_, pois_inhom_z) -				\
		   field_ (r_, z_, -pois_inhom_z));				\
										\
	 for (i__ = 2; i__ <= pois_inhom_reflections; i__ += 2) {		\
	   total_ += (field_ (r_, z_, pois_inhom_z + i__ * pois_inhom_L) +	\
		      field_ (r_, z_, pois_inhom_z - i__ * pois_inhom_L));	\
	 									\
	   total_ -= (field_ (r_, z_, -pois_inhom_z + i__ * pois_inhom_L) +	\
		      field_ (r_, z_, -pois_inhom_z - i__ * pois_inhom_L));	\
	 }								\
       } while (0)

/** @brief Returns the potential created by a unit charge located at pois_inhom_z
 *  between an electrode at 0 and a second one at pois_inhom_L.
 */
double
pois_inhom_phi (double r, double z)
{
  double res;

  sum_reflections (res, inhom_phi_term, r, z);

  return res;
}

/** @brief Radial component of the field. */
double
pois_inhom_er (double r, double z)
{
  double res;

  sum_reflections (res, inhom_er_term, r, z);

  return res;
}

/** @brief Axial component of the field. */
double
pois_inhom_ez (double r, double z)
{
  double res;

  sum_reflections (res, inhom_ez_term, r, z);

  return res;
}

/** @brief Once we have the electrostatic potential created by the space charges,
 * we use this function to compute the multiplying factor of the floating
 * charge.
 *
 * One can also forget about keeping the potential fixed somewhere and impose
 * a constant (non-floating, but may depend on time) charge pois_inhom_fixed_q,
 * which is used if nonzero. The use of that is to simulate a charged cloud
 * close to the earth, which acts as an electrode: usually you would also combine
 * that with the sprite module, that allows varing neutral densities.
 */
double
pois_inhom_q_factor (pois_grid_t *pois)
{
  double u;
  /* Note that the signs preceding the numerical potentials are inverted since
   * in our notation, we use the opposite of the physical one, in order
   * to use the charge as the source of the Poisson equation.
   */

  if (pois_inhom_fixed_q != 0.0) {
    return pois_inhom_fixed_q_t;
  }

  /* Was: 
  u = - E_z * (pois_inhom_L - L_z) + pois_phi_at (pois, 0.0, L_z, 0.0);
  */

  u = E_z * (L_z - pois_inhom_L) + pois_phi_at (pois, 0.0, L_z, 0.0);

  return  -u / pois_inhom_dphi;
}

/** @brief Adds the potential created by the needle to a Poisson grid and their
 *  descendants.
 *
 * Note that \f$\phi\f$ here is in Fourier space and the Laplacian
 * potential is always axi-symmetric, so you only have to call this function
 * with the zero-mode of the Poisson grids.
 *
 * The reason that it is better to add the inhomogeneous field to the potential
 * and not later to the electric fields is that if we do that, we will
 * have to substract two large and very curved functions that give something
 * close to zero: the error we make then will be comparable to the result.
 */
void
pois_add_inhom_phi_r (pois_grid_t *grid, double q)
{
  int ir, iz;
  pois_grid_t *child;
  double res;

  debug (2, "pois_add_inhom_phi_r (" grid_printf_str ", q = %f)\n",
	 grid_printf_args(grid), q);

  /* Note that there is no factor grid->ntheta coming from thr FFT 
     un-normalization.  The reason is that actually q is sqrt(N) smaller
     than it has to be and later, in the inverse FFT we will multiply by
     sqrt(N).  So everything is correct, I believe. */
  iter_grid_n (grid, ir, iz, 2) {

    sum_reflections (res, inhom_phi_term, r_at(ir, grid->level), z_at(iz, grid->level));


    RZ(grid->phi, ir, iz) += q * res;
  }

  iter_childs (grid, child) {
    pois_add_inhom_phi_r (child, q);
  }
}

/***************************************************************
 *  Mapper functions for the components of the electric field. *
 *  See mapper.h for an explanation of each of these methods.  *
 ***************************************************************/

/** @brief er_copy QQQQ */
void
er_copy (mapper_t *mapper, grid_t *source, grid_t *target, 
	 int ir, int iz, int itheta)
{
  cdr_grid_t *cdr;
  pois_grid_t *pois;

  cdr = (cdr_grid_t*) target;
  pois = (pois_grid_t*) source;

  if (ir < pois->r1)
    RZT (cdr->er, ir, iz, itheta) = ER_RZ (pois, ir, iz);
}

/** @brief ez_copy QQQQ */
void
ez_copy (mapper_t *mapper, grid_t *source, grid_t *target, 
	 int ir, int iz, int itheta)
{
  cdr_grid_t *cdr;
  pois_grid_t *pois;

  cdr = (cdr_grid_t*) target;
  pois = (pois_grid_t*) source;

  if (iz < pois->z1)
    RZT (cdr->ez, ir, iz, itheta) = EZ_RZ (pois, ir, iz);
}

/** @brief etheta_copy QQQQ */
void
etheta_copy (mapper_t *mapper, grid_t *source, grid_t *target, 
	     int ir, int iz, int itheta)
{
  cdr_grid_t *cdr;
  pois_grid_t *pois;

  cdr = (cdr_grid_t*) target;
  pois = (pois_grid_t*) source;

  RZT (cdr->etheta, ir, iz, itheta) = RZ (pois->phi, ir, iz);
}

/** @brief er_coarsen QQQQ */
void 
er_coarsen (mapper_t *mapper, grid_t *source, grid_t *target, 
	    int ir, int iz, int itheta)
{
  cdr_grid_t *cdr;
  pois_grid_t *pois;
  int level_diff, er_ir, er_iz;

  cdr = (cdr_grid_t*) target;
  pois = (pois_grid_t*) source;

  level_diff = pois->level - cdr->level;

  er_iz = (iz << level_diff) + (1 << (level_diff - 1));
  er_ir = ((ir + 1) << level_diff) - 1;

  if (grid_contains (source, er_ir, er_iz - 1, GRID_INSIDE) && 
      grid_contains (source, er_ir, er_iz, GRID_INSIDE)) {

    RZT (cdr->er, ir, iz, itheta) = 0.5 * (ER_RZ(pois, er_ir, er_iz) +
					   ER_RZ(pois, er_ir, er_iz - 1));
  }
}

/** @brief ez_coarsen QQQQ */
void 
ez_coarsen (mapper_t *mapper, grid_t *source, grid_t *target, 
	    int ir, int iz, int itheta)
{
  cdr_grid_t *cdr;
  pois_grid_t *pois;
  int level_diff, ez_ir, ez_iz;

  cdr = (cdr_grid_t*) target;
  pois = (pois_grid_t*) source;

  level_diff = pois->level - cdr->level;

  ez_iz = ((iz + 1) << level_diff) - 1;
  ez_ir = (ir << level_diff) + (1 << (level_diff - 1));

  if (grid_contains (source, ez_ir - 1, ez_iz, GRID_INSIDE) && 
      grid_contains (source, ez_ir, ez_iz, GRID_INSIDE)) {

    RZT(cdr->ez, ir, iz, itheta) = 0.5 * (EZ_RZ(pois, ez_ir, ez_iz) +
					  EZ_RZ(pois, ez_ir - 1, ez_iz));
  }
}

/** @brief etheta_coarsen QQQQ */
void 
etheta_coarsen (mapper_t *mapper, grid_t *source, grid_t *target, 
		int ir, int iz, int itheta)
{
  cdr_grid_t *cdr;
  pois_grid_t *pois;
  int level_diff, ez_ir, er_iz;

  cdr = (cdr_grid_t*) target;
  pois = (pois_grid_t*) source;

  level_diff = pois->level - cdr->level;

  er_iz = (iz << level_diff) + (1 << (level_diff - 1));
  ez_ir = (ir << level_diff) + (1 << (level_diff - 1));

  if (grid_contains (source, ez_ir, er_iz, GRID_INSIDE) && 
      grid_contains (source, ez_ir - 1, er_iz, GRID_INSIDE) &&
      grid_contains (source, ez_ir, er_iz - 1, GRID_INSIDE) &&
      grid_contains (source, ez_ir - 1, er_iz - 1, GRID_INSIDE)) {

    RZT (cdr->etheta, ir, iz, itheta)  = 0.25 *
      (RZ (pois->phi, ez_ir, er_iz) 
       + RZ (pois->phi, ez_ir - 1, er_iz) 
       + RZ (pois->phi, ez_ir, er_iz - 1) 
       + RZ (pois->phi, ez_ir - 1, er_iz - 1));
  }
}

/** @brief er_interpol QQQQ */
int 
er_interpol_set (mapper_t *mapper, grid_t *source, interpol_t *interpol,
		      int pr, int pz, int itheta)
{
  pois_grid_t *pois;

  pois = (pois_grid_t*) source;

  /* When we set Neumann b.c. in one boundary and Dirichlet in the opposite
   * a larger-than-usual error appears in the boundaries because we are
   * matching potetntials in one side and fields in the other.  The only way
   * I see of avoiding this is to forget about the extrapolated values of
   * the grids.  But we still need those values on the external boundaries.
   */
  if (pr < pois->r0 || pz < pois->z0
      || pr > pois->r1 - 1 || pz > pois->z1)
    return FALSE;

  interpol_set_stencil (interpol, 
			er_r_at (pr - 1, pois->level), 
			er_z_at (pz - 1, pois->level),
			ER_RZ (pois, pr - 1, pz - 1), 
			ER_RZ (pois, pr - 1, pz),
			ER_RZ (pois, pr, pz - 1), 
			ER_RZ (pois, pr, pz));

  return TRUE;
}

/** @brief ez_interpol_set QQQQ */
int
ez_interpol_set (mapper_t *mapper, grid_t *source, interpol_t *interpol,
		      int pr, int pz, int itheta)
{
  pois_grid_t *pois;

  pois = (pois_grid_t*) source;

  if (pr < pois->r0 || pz < pois->z0
      || pr > pois->r1 || pz > pois->z1 - 1)
    return FALSE;

  interpol_set_stencil (interpol, 
			ez_r_at (pr - 1, pois->level), 
			ez_z_at (pz - 1, pois->level),
			EZ_RZ (pois, pr - 1, pz - 1), 
			EZ_RZ (pois, pr - 1, pz),
			EZ_RZ (pois, pr, pz - 1), 
			EZ_RZ (pois, pr, pz));
  return TRUE;
}

/** @brief etheta_interpol_set QQQQ */
int 
etheta_interpol_set (mapper_t *mapper, grid_t *source, interpol_t *interpol,
		      int pr, int pz, int itheta)
{
  pois_grid_t *pois;

  pois = (pois_grid_t*) source;

  if (pr < pois->r0 || pz < pois->z0
      || pr > pois->r1 - 1|| pz > pois->z1 - 1)
    return FALSE;

  interpol_set_stencil_at (source, interpol, 
			   r_at (pr, pois->level),
			   z_at (pz, pois->level),
			   pois->phi, pr, pz, 0);
  return TRUE;
}

/** @brief er_interpol QQQQ */
void
er_interpol (mapper_t *mapper, grid_t *source, grid_t *target, 
	     interpol_t *interpol, int ir, int iz, int itheta)
{
  double er_r, er_z;
  cdr_grid_t *cdr;
  pois_grid_t *pois;
  int level_diff;

  cdr = (cdr_grid_t *) target;
  pois = (pois_grid_t *) source;

  level_diff = cdr->level - pois->level;

  /* Note that when pr = pois->r1, the electric field is calculated
   * from the extrapolation of phi and it is hence less accurate.
   * The same applies below to pz and pois->z1.
   */
  if ((ir >> level_diff) >= pois->r1) return;

  er_r = er_r_at (ir, cdr->level);
  er_z = er_z_at (iz, cdr->level);

  RZT (cdr->er, ir, iz, itheta) = interpol_apply (interpol, er_r, er_z);
}

/** @brief ez_interpol QQQQ */
void
ez_interpol (mapper_t *mapper, grid_t *source, grid_t *target, 
	     interpol_t *interpol, int ir, int iz, int itheta)
{
  double ez_r, ez_z;
  cdr_grid_t *cdr;
  pois_grid_t *pois;
  int level_diff;

  cdr = (cdr_grid_t *) target;
  pois = (pois_grid_t *) source;

  level_diff = cdr->level - pois->level;

  if ((iz >> level_diff) >= pois->z1) return;

  ez_r = ez_r_at (ir, cdr->level);
  ez_z = ez_z_at (iz, cdr->level);

  RZT (cdr->ez, ir, iz, itheta) = interpol_apply (interpol, ez_r, ez_z);
}

/** @brief etheta_interpol QQQQ */
void
etheta_interpol (mapper_t *mapper, grid_t *source, grid_t *target, 
		 interpol_t *interpol, int ir, int iz, int itheta)
{
  double eth_r, eth_z;
  cdr_grid_t *cdr;

  cdr = (cdr_grid_t *) target;

  eth_r = r_at (ir, cdr->level);
  eth_z = z_at (iz, cdr->level);

  RZT (cdr->etheta, ir, iz, itheta) = 
    interpol_apply (interpol, eth_r, eth_z);
}

/** @brief Mapping of the charge. */
void
charge_copy (mapper_t *mapper, grid_t *source, grid_t *target, 
	     int ir, int iz, int itheta)
{
  cdr_grid_t *cdr;
  pois_grid_t *pois;

  pois = (pois_grid_t*) target;
  cdr = (cdr_grid_t*) source;

  RZ (pois->charge, ir, iz) = RZT (cdr->charge, ir, iz, itheta);
}

/** @brief charge_interpol_set QQQQ. */
int
charge_interpol_set (mapper_t *mapper, grid_t *source, interpol_t *interpol,
		     int pr, int pz, int itheta)
{
  cdr_grid_t *cdr;

  cdr = (cdr_grid_t*) source;

  if (pr < cdr->r0 || pz < cdr->z0
      || pr > cdr->r1 - 1|| pz > cdr->z1 - 1)
    return FALSE;

  interpol_set_stencil_at (source, interpol, 
			   r_at (pr, cdr->level),
			   z_at (pz, cdr->level),
			   cdr->charge, pr, pz, itheta);
  return TRUE;
}

/** @brief charge_interpol QQQQ. */
void
charge_interpol (mapper_t *mapper, grid_t *source, grid_t *target, 
		 interpol_t *interpol, 
		 int ir, int iz, int itheta)
{
  double r, z;
  pois_grid_t *pois;

  pois = (pois_grid_t *) target;

  r = r_at (ir, pois->level);
  z = z_at (iz, pois->level);

  RZT (pois->charge, ir, iz, 0) = interpol_apply (interpol, r, z);
}

/** @brief Finds the best-possible approximation for the potential at a given
 * point (@a r, @a z, @a theta) by interpolating from the finest grid that covers
 * that point.
 */
double
pois_phi_at (pois_grid_t *grid, double r, double z, double theta)
{
  int ir, iz, itheta, it, itn;
  REAL phi_it[2];
  interpol_t *phi_at_interpol = NULL;

  debug (3, "pois_phi_at (" grid_printf_str ", r = %g, z = %g, theta = %g)\n",
	 grid_printf_args (grid), r, z, theta);

  grid = (pois_grid_t *) grid_finest_containing_r ((grid_t*) grid, r, z);

  if (NULL == grid) {
    fatal ("In pois_phi_at: potential outside the grid");
  }

  phi_at_interpol = interpol_new_a (dr[grid->level], dz[grid->level],
				    &interpol_quadratic);

  ir = (int) (r / dr[grid->level]);
  iz = (int) (z / dz[grid->level]);

  /* If r == 0, it doesn't matter which theta we use, so we use 0.
   * also if we are in a 2D simulation. 
   */
  if (grid->ntheta > 1 && r > 0.0) {
    itheta = (int) (theta / dtheta);
    itn = 2;
  } else {
    itheta = 0;
    itn = 1;
  }

  for (it = 0; it < itn; it++) {
    interpol_set_stencil_at ((grid_t *) grid, phi_at_interpol,
			     r_at (ir, grid->level),
			     z_at (iz, grid->level), 
			     grid->phi, ir, iz, itheta + it);

    phi_it[it] = interpol_apply (phi_at_interpol, r, z);
  }

  interpol_free (phi_at_interpol);

  if (grid->ntheta > 1 && r > 0.0) {
    return (phi_it[0] * (theta - theta_at (itheta)) +
	    phi_it[1] * (theta_at (itheta + 1) - theta)) / dtheta;
  } else {
    return phi_it[0];
  }
}

/** @brief Dumps all the contents of the given grid into filenames given by
 * prefix and name
 */
void
pois_dump (pois_grid_t *grid, const char *prefix, const char *name)
{
  char *fname;
  int m = pois_output_margin;

  asprintf (&fname, "%s/r.%s.tsv", prefix, name);
  rz_axis_dump (fname, grid->r0 - m, grid->r1 + m, dr[grid->level]);
  free (fname);

  asprintf (&fname, "%s/z.%s.tsv", prefix, name);
  rz_axis_dump (fname, grid->z0 - m, grid->z1 + m, dz[grid->level]);
  free (fname);

  asprintf (&fname, "%s/phi.%s.tsv", prefix, name);
  rz_dump (grid->phi, fname, "w", grid->r0 - m, grid->z0 - m, 
	   grid->r1 + m, grid->z1 + m);
  free (fname);

  asprintf (&fname, "%s/charge.%s.tsv", prefix, name);
  rz_dump (grid->charge, fname, "w", grid->r0 - m, grid->z0 - m, 
	   grid->r1 + m, grid->z1 + m);
  free (fname);

  asprintf (&fname, "%s/error.%s.tsv", prefix, name);
  rz_dump (grid->error, fname, "w", grid->r0 - m, grid->z0 - m, 
	   grid->r1 + m, grid->z1 + m);
  free (fname);
}

/** @brief Dumps all the contents of the given grid into filenames given by
 * prefix and name
 */
void
pois_dump_r (pois_grid_t *grid, const char *prefix, const char *name)
{
  pois_grid_t *child;
  char *cname;
  int i = 0, nchilds;
  char codes[] = "abcdefghijklmnopqrstuvwxyz";

  pois_dump (grid, prefix, name);

  nchilds = grid_howmany_children ((grid_t *) grid);

  for (i = 0; i < nchilds; i++) {
    child = (pois_grid_t *) grid_get_child ((grid_t*) grid, nchilds - i - 1);

    assert (NULL != child);

    if (i > sizeof(codes) - 2)
      asprintf (&cname, "%s{%d}", name, i);
    else
      asprintf (&cname, "%s%c", name, codes[i]);

    pois_dump_r (child, prefix, cname);
    free (cname);
  }
}
