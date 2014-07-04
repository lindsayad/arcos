/** @file mapper.c
 *  @brief All functions related to the mapping of whatever variable
 *  (densities, electric field...) from one tree into another.
 *
 *  Here the relevant data structure is defined by mapper_t, which tells
 *  us what functions should we call to perform copying, interpolation
 *  and coarsening from one grid into another (not neccesarily ofthe same type). */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "cstream.h"
#include "grid.h"
#include "interpol2.h"
#include "mapper.h"
#include "proto.h"
#include "species.h"

static void map_same (mapper_t **mappers, 
		      grid_t *source, grid_t *target, int ntheta,
		      int r0, int z0, int r1, int z1);
static void map_interpol (mapper_t **mappers, grid_t *source, grid_t *target, 
			  int ntheta, int r0, int z0, int r1, int z1, 
			  int level_diff);
static void interpol_inner (mapper_t **mappers, grid_t *source, grid_t *target,
			    int pr, int pz, int ntheta,	
			    interpol_t **interpols, 
			    int *extend_r, int *extend_z);
static void map_coarsen (mapper_t **mappers, grid_t *source, grid_t *target, 
			 int ntheta, int r0, int z0, int r1, int z1, 
			 int level_diff);


/** @brief Maps two grids.
 *
 * The arguments @a copy, @a interpol and @a coarsen are
   booleans that specify whether we should apply these methods
   in the mapping. */
void
map_grid (mapper_t **mappers, grid_t *source, grid_t *target, int ntheta, 
	  int copy, int interpol, int coarsen, int s_buff, int t_buff)
{
  int r0, z0, r1, z1, level_diff;
  int br0, bz0, br1, bz1;
  mapper_t **m;
  int overlap;

  debug (3, "map_grid (source = " grid_printf_str 
	 ",\n\t dest = " grid_printf_str 
	 ",\n\t s_buff = %d, t_buff = %d)\n",
	 grid_printf_args (source), grid_printf_args (target),
	 s_buff, t_buff);

  overlap = grid_overlap (source, target, s_buff, t_buff, 
			  &r0, &z0, &r1, &z1, &level_diff);

  if (!overlap) {
    debug (3, " -> no overlap\n");
    return;
  }
  
  if (level_diff == 0 && copy) {
    map_same (mappers, source, target, ntheta, r0, z0, r1, z1);
  } else if (level_diff < 0 && interpol) {
    /* When we interpolate electric fields, the stencil used in the coarse
       grid does not map exactly to cells in the finer grid by the standard
       rules.  Hence may have to extend the range of fine-grid cells that
       can be affected by a coarse grid. 
       But later we must be careful, because not all fine-grid cells can
       be interpolated.
    */
    for (m = mappers; *m; m++) {
      int max_shift;
     
      max_shift =  MYMAX(abs((*m)->shift_r) << (-level_diff - 1),
			 abs((*m)->shift_z) << (-level_diff - 1));

      debug (3, "max_shift = %d, shift_r = %d, shift_z = %d\n", 
	     max_shift, (*m)->shift_r, (*m)->shift_z);

      grid_overlap_with_shifts (source, target, 
				s_buff, 
				MYMAX (t_buff, max_shift),
				&br0, &bz0, &br1, &bz1, &level_diff,
				(*m)->shift_r, (*m)->shift_z);
      if (br0 < r0) r0 = br0;
      if (bz0 < z0) z0 = bz0;
      if (br1 > r1) r1 = br1;
      if (bz1 > z1) z1 = bz1;
    }
    map_interpol (mappers, source, target, ntheta, r0, z0, r1, z1, 
		  -level_diff);
  } else if (level_diff > 0 && coarsen) {
    map_coarsen (mappers, source, target, ntheta, r0, z0, r1, z1, level_diff);
  }
}


/** @brief Recursive version of map_grid.
 *
 * Recurses over the source tree depth-first
 */
void
map_grid_r (mapper_t **mappers, grid_t *source, grid_t *target, int ntheta,
	    int copy, int interpol, int coarsen, int s_buff, int t_buff)
{
  grid_t *child;

  debug (3, "map_grid_r (...)\n");

  map_grid (mappers, source, target, ntheta, copy, interpol, coarsen,
	    s_buff, t_buff);

  iter_childs (source, child) {
    map_grid_r (mappers, child, target, ntheta, copy, interpol, coarsen,
		s_buff, t_buff);
  }
}

/** @brief Double-recursive version:
 * 
 * Also recurses over the target tree.
 */
void
map_trees_r (mapper_t **mappers, grid_t *source, grid_t *target, int ntheta,
	     int copy, int interpol, int coarsen, int s_buff, int t_buff)
{
  grid_t *child;

  debug (3, "map_trees_r (...)\n");

  iter_childs (target, child) {
    map_trees_r (mappers, source, child, ntheta, copy, interpol, coarsen,
		s_buff, t_buff);
  }

  map_grid_r (mappers, source, target, ntheta, copy, interpol, coarsen,
	      s_buff, t_buff);
}

/** @brief The level of the source and target grids is the same.
 *
 * This is the easiest case.
 */
static void
map_same (mapper_t **mappers, 
	  grid_t *source, grid_t *target, int ntheta,
	  int r0, int z0, int r1, int z1)
{
  int ir, iz;
  mapper_t **m;

  debug (3, "map_same (source = " grid_printf_str 
	 ",\n\t dest = " grid_printf_str 
	 ",\n\t r0 = %d, z0 = %d, r1 = %d, z1 = %d)\n",
	 grid_printf_args (source), grid_printf_args (target),
	 r0, z0, r1, z1);

  for (ir = r0; ir < r1; ir++) {
    for (iz = z0; iz < z1; iz++) {
      for (m = mappers; *m; m++) {
	(*m)->copy (*m, source, target, ir, iz, ntheta);
      }
    }
  }
  debug (3, " < map_same (...)\n");
}

/** @brief With this, we limit the number of mappers.
 *
 * But it is extremely improbable that someone would like to use more
 * than 32 mappers at the same time. */
#define MAX_MAPPERS  32

/** @brief The source grid is coarser than the target grid. */
static void
map_interpol (mapper_t **mappers, grid_t *source, grid_t *target, int ntheta,
			 int r0, int z0, int r1, int z1, int level_diff)
{
  int pr, pz;
  int i, extend_r, extend_z;
  interpol_t *interpols[MAX_MAPPERS];
  mapper_t **m;

  debug (3, "map_interpol (source = " grid_printf_str 
	 ",\n\t dest = " grid_printf_str 
	 ",\n\t r0 = %d, z0 = %d, r1 = %d, z1 = %d, level_diff = %d)\n",
	 grid_printf_args (source), grid_printf_args (target),
	 r0, z0, r1, z1, level_diff);

  for (m = mappers, i = 0; *m; m++, i++) 
    interpols[i] = interpol_new_a (dr[source->level], dz[source->level], 
				   (*m)->interpol_method);

  /* The +- 1 are there because even when we are including the
     boundaries to calculate the overlap, we have to let one buffer cell
     for the 9-point interpolation. */
  r0 = MYMAX (r0 >> level_diff, source->r0 - 1);

  /* Explanation of the formula used: r1 is the index of the first cell not
     in the overlap 
     -> r1 - 1 is the last cell in the overlap 
     -> (r1 - 1) >> level_diff indexes the last cell in the coarser grid
        that contains part of the overlap.
     -> since the loop has to include the former index, we add one.
  */
  r1 = MYMIN (((r1 - 1) >> level_diff) + 1, source->r1 + 1);
  z0 = MYMAX (z0 >> level_diff, source->z0 - 1);
  z1 = MYMIN (((z1 - 1) >> level_diff) + 1, source->z1 + 1);


  debug (4, "r0 = %d, z0 = %d, r1 = %d, z1 = %d\n",
	 r0, z0, r1, z1);

  /* assert (r0 < r1 && z0 < z1); */

  for (pr = r0, extend_r = TRUE; pr < r1; pr++) {
    /* The 0 initialization is to shut up the compiler. */
    int extend_r_in = 0;
    for (pz = z0, extend_z = TRUE; pz < z1; pz++) {
      debug (6, "pr = %d, pz = %d, extend_r = %d, extend_z = %d\n", 
	     pr, pz, extend_r, extend_z);
      extend_r_in = extend_r;

      interpol_inner (mappers, source, target, pr, pz, ntheta,
		      interpols, &extend_r_in, &extend_z);
    }     
    extend_r = extend_r_in;
  }

  for (m = mappers, i= 0; *m; m++, i++) 
    interpol_free (interpols[i]);

  debug (3, "< map_interpol (...)\n");
}

/** @brief Calculates all the points in a target grid, finer than a source
 *  grid that are interpolated from a stencil centered at @a pr, @a pz in
 *  the source grid.
 *
 *  @a interpol_xxx are interpolation objects that can thus be reused.
 *  @a extend_r and @a extend_z are here because one call to @a interpol_inner
 *  can actually set (1 << level_diff + 1) fine grid points, but the boundaries
 *  are shared with the next/previous call.
 *
 *  In the first call in a row/column, we set all the points, since the
 *  right/lower boundary is not shared.
 */
static void
interpol_inner (mapper_t **mappers, grid_t *source, grid_t *target, 
		int pr, int pz, int ntheta,
		interpol_t **interpols, int *extend_r, int *extend_z) 
{
  int level_diff;
  int i, ir, iz;
  int rmin, rmax, zmin, zmax;
  int inside;

  mapper_t **m;

  level_diff = target->level - source->level;

  debug (6, "interpol_inner (source = " grid_printf_str ",\n\t\t target = "
	 grid_printf_str ",\n\t\t pr = %d, pz = %d, ntheta = %d)\n",
	 grid_printf_args (source), grid_printf_args (target), pr, pz, ntheta);

  for (m = mappers, i = 0; *m; m++, i++) {
    int fine_shift_r, fine_shift_z;

    inside = (*m)->interpol_set (*m, source, interpols[i], pr, pz, ntheta);
  
    if (!inside) continue;

    fine_shift_r = (*m)->shift_r << (level_diff - 1);
    fine_shift_z = (*m)->shift_z << (level_diff - 1);

    debug (6, "shift_r = %d, shift_z = %d\n", (*m)->shift_r, (*m)->shift_z);

    rmin = (pr << level_diff) + fine_shift_r;
    /* (*m)->shift_z tells us if the grid is staggered.  If it is not,
       it doesn't make sense to extend it. */
    if (*extend_r && (*m)->shift_z != 0) {
      *extend_r = FALSE;
      rmin--;
    }

    zmin = (pz << level_diff) + fine_shift_z;
    if (*extend_z && (*m)->shift_r != 0) {
      *extend_z = FALSE;
      zmin--;
    }

    /* The 2 are there because there is no problem in writing in the
       buffer margin of the target grid. */
    rmin = MYMAX (target->r0 - 2, rmin);
    rmax = MYMIN (target->r1 + 2, ((pr + 1) << level_diff) + fine_shift_r);
    zmin = MYMAX (target->z0 - 2, zmin);
    zmax = MYMIN (target->z1 + 2, ((pz + 1) << level_diff) + fine_shift_z);
    

    debug(6, "\trmin = %d, rmax = %d, zmin = %d, zmax = %d\n",
	  rmin, rmax, zmin, zmax);

    for (ir = rmin; ir < rmax; ir++) {
      for (iz = zmin; iz < zmax; iz++) {
	(*m)->interpol (*m, source, target, interpols[i], 
			ir, iz, ntheta);
      }   
    }
  }
}

/** @brief The target grid is coarser than the source grid. */
static void
map_coarsen (mapper_t **mappers, grid_t *source, grid_t *target, int ntheta,
	     int r0, int z0, int r1, int z1, int level_diff)
{
  int i, ir, iz;
  int rmin, rmax, zmin, zmax;
  mapper_t **m;

  debug (3, "map_coarsen (source = " grid_printf_str 
	 ",\n\t dest = " grid_printf_str 
	 ",\n\t r0 = %d, z0 = %d, r1 = %d, z1 = %d, level_diff = %d)\n",
	 grid_printf_args (source), grid_printf_args (target),
	 r0, z0, r1, z1, level_diff);

  rmin = MYMAX (r0 >> level_diff, target->r0 - 2);

  /* See above for an explanation for this formula. */
  rmax = MYMIN (((r1 - 1) >> level_diff) + 1, target->r1 + 2);
  zmin = MYMAX (z0 >> level_diff, target->z0 - 2);
  zmax = MYMIN (((z1 - 1) >> level_diff) + 1, target->z1 + 2);
  
  for (iz = zmin; iz < zmax; iz++) {
    for (ir = rmin; ir < rmax; ir++) {
      for (m = mappers, i = 0; *m; m++, i++) {
	(*m)->coarsen (*m, source, target, ir, iz, ntheta);
      }
    } 
  }
  debug (3, "< map_coarsen (...)\n");
}

