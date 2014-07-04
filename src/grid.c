/** @file grid.c
 *  @brief General grid functions common for cdr and poisson grids.
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "grid.h"
#include "parameters.h"
#include "proto.h"
#include "species.h"

/** @brief Gets the maximum depth of a tree (the level of its deepest node) */
int
grid_max_depth_r (grid_t *grid)
{
  int l, lc;
  grid_t *child;

  l = grid->level;
  iter_childs (grid, child) {
    lc = grid_max_depth_r (child);
    l = (l > lc)? l: lc;
  }
  return l;
}

/** @brief Gets the minimum r of a grid and all its descendants */
double
grid_rmin_r (grid_t *grid)
{
  int lev;
  double minr, cminr;
  grid_t *child;

  lev = grid->level;
  minr = r_at (grid->r0, lev);

  iter_childs (grid, child) {
    cminr = grid_rmin_r (child);
    minr = (minr > cminr)? cminr: minr;
  }
  return minr;
}

/** @brief Recursively prints one grid and all its descendants. */
void
grid_print_r (grid_t *grid, int indent)
{
  int i;
  grid_t *child;

  for (i = 0; i < indent; i++) fputs ("  ", stdout);
  printf ("grid (" grid_printf_str") {\n", grid_printf_args(grid));

  iter_childs (grid, child) {
    grid_print_r (child, indent + 1);
  }
  for (i = 0; i < indent; i++) fputs ("  ", stdout);
  fputs ("}\n", stdout);
}

/** @brief Checks whether a point is in the interior or on a given boundary
 * of a grid.
 *
 * Does the grid contain the point (i, j) (in its own coordinates)?\n
 * @a check has the following meaning:\n
 *  GRID_INSIDE : Checks if the point is in the interior of the grid\n
 *  BND_MASK_*  : Checks if it is on the given boundary (BND_MASK_ALL checks for
 *               any boundary).\n
 */
int
grid_contains (grid_t *grid, int i, int j, int check)
{
  int r0, r1, z0, z1;
  int r;

  r0 = grid->r0;
  z0 = grid->z0;
  r1 = grid->r1 - 1;
  z1 = grid->z1 - 1;

  r = FALSE;

  if (check & BND_MASK_LEFT) {
    r0 -= 1;
  }

  if (check & BND_MASK_RIGHT) {
    r1 += 1;
  }

  if (check & BND_MASK_BOTTOM) {
    z0 -= 1;
  }

  if (check & BND_MASK_TOP) {
    z1 += 1;
  }

  if (r0 < i && i < r1 && z0 < j && j < z1) 
    r = TRUE;

  if (check & GRID_INSIDE) {
    return r;
  } else {
    return (r && !grid_contains (grid, i, j, GRID_INSIDE));
  }
}

/** @brief Checks whether two grids overlap and return the rectangle that they
 * in common.
 *
 * Checks whether two grids overlap and return the rectangle that they 
 * have in common in the coordinates of the finest grid in 
 * (*left, *bottom) - (*right, *top)\n
 *
 * Note that *left, *bottom... can be modified even if the two grids do
 * not overlap.
 * @a buff1 and @a buff2 are integers that specify the numbers of cells
 * that are taken as buffer in the boundaries for @a grid1 and @a grid2.
 */
int
grid_overlap (grid_t *grid1, grid_t *grid2, int buff1, int buff2, 
	      int *left, int *bottom, int *right, int *top, int *level_diff)
{
  return grid_overlap_with_shifts (grid1, grid2, buff1, buff2, 
				  left, bottom, right, top, level_diff,
				  0, 0);
}

/** @brief Checks the overlap between two grids.
 *
 * Checks the overlap between two grids before it shifts the
 * coarser grid by (shift_r << (level_diff - 1), shift_z << (level_diff - 1)
 * (in finer grid coordinates).
 *
 *  Q: WTF! WHY!!!? \n
 *  A: We use this to implement the mapping between electric fields.  Since
 *     The electric field is evaluated not in the cell centers but in the
 *     borders, when we interpolate we have to take a shift into account.
 *     See mapper.c and mapper.h for more info.
 */
int
grid_overlap_with_shifts (grid_t *grid1, grid_t *grid2, int buff1, int buff2, 
			 int *left, int *bottom, int *right, int *top, 
			 int *level_diff, int shift_r, int shift_z)
{
  int sign;

  debug (3, "grid_overlap(" grid_printf_str", " grid_printf_str 
	 ", buff1 = %d, buff2 = %d)\n",
	 grid_printf_args(grid1), grid_printf_args(grid2),
	 buff1, buff2);

  sign = 1;

  /* If level_diff < 0, grid2 is finer */
  *level_diff = grid1->level - grid2->level;

  if (*level_diff < 0) {
    sign = -1;
    *level_diff = - *level_diff;
    XCHG (grid1, grid2);
    XCHG (buff1, buff2);
  }

  /* At this point, level_diff >= 0 and grid1 is finer 
     or equivalent to grid2 */

  *left = MAX_AT_LEVEL_WITH_SHIFT (grid1->r0 - buff1, grid2->r0 - buff2, 
				   *level_diff, shift_r);
  *right = MIN_AT_LEVEL_WITH_SHIFT (grid1->r1 + buff1, grid2->r1 + buff2, 
				    *level_diff, shift_r);

  /* Note: *right is not contained in the grid.  Hence, if *left == *right
     the two grids do not overlap.  Same for *bottom and *top, below. */
  debug (3, "*left = %d, *right = %d\n", *left, *right);
  if (*left >= *right) {
    *level_diff *= sign;
    return FALSE;
  }

  *bottom = MAX_AT_LEVEL_WITH_SHIFT (grid1->z0 - buff1, grid2->z0 - buff2, 
				     *level_diff, shift_z);
  *top = MIN_AT_LEVEL_WITH_SHIFT (grid1->z1 + buff1, grid2->z1 + buff2, 
				  *level_diff, shift_z);

  debug (3, "*bottom = %d, *top = %d\n", *bottom, *top);
  if (*bottom >= *top) {
    *level_diff *= sign;
    return FALSE;
  }
  debug (3, "pois_overlap = TRUE\n");

  *level_diff *= sign;
  return TRUE;
}

/** @brief Makes a grid inherit the ext_bound field of its parent,
 *
 * with the appropiate modifications.
 */
void
grid_inherit_ext_bound (grid_t *grid)
{
  grid_t *parent = grid->parent;

  debug (3, "grid_inherit_ext_bound (...)\n");

  grid->ext_bound = BND_NONE;

  /* If the parent does not have any external boundary, the child will also
     not have */
  if (BND_NONE == parent->ext_bound) {
    return;
  }

  if (parent->ext_bound & BND_MASK_LEFT) {
    if ((parent->r0 << 1) == grid->r0) grid->ext_bound |= BND_MASK_LEFT;
  }

  if (parent->ext_bound & BND_MASK_RIGHT) {
    if ((parent->r1 << 1) == grid->r1) grid->ext_bound |= BND_MASK_RIGHT;
  }

  if (parent->ext_bound & BND_MASK_BOTTOM) {
    if ((parent->z0 << 1) == grid->z0) grid->ext_bound |= BND_MASK_BOTTOM;
  }

  if (parent->ext_bound & BND_MASK_TOP) {
    if ((parent->z1 << 1) == grid->z1) grid->ext_bound |= BND_MASK_TOP;
  }

}


/** @brief Finds the finest descendant of grid that contains the point given by
 * the (@a r, @a z) coordinates.
 *
 * Or NULL if the point is not contained in grid.
 * We assume that only ONE sub-grid can contain a given point.
 */
grid_t *
grid_finest_containing_r (grid_t *grid, double r, double z)
{
  grid_t *leaf, *ret;
  
  if (r < er_r_at(grid->r0 - 1, grid->level) || 
      r > er_r_at(grid->r1 - 1, grid->level))
    return NULL;

  if (z < ez_z_at(grid->z0 - 1, grid->level) || 
      z > ez_z_at(grid->z1 - 1, grid->level))
    return NULL;

  iter_childs (grid, leaf) {
    ret = grid_finest_containing_r (leaf, r, z);
    if (NULL != ret)
      return ret;
  }
  return grid;
}

/** @brief How many children does one grid have? */
int
grid_howmany_children (grid_t *grid)
{
  grid_t *leaf;
  int count;

  count = 0;
  iter_childs (grid, leaf) {
    count ++;
  }
  return count;
}

/** @brief Finds child number @a n of the given grid.
 *
 * If it doesn't have enough children just returns NULL
 */
grid_t*
grid_get_child (grid_t *grid, int n)
{
  grid_t *leaf;
  int count;

  count = 0;
  iter_childs (grid, leaf) {
    if (count == n) return leaf;
    count ++;
  }
  return NULL;
}  
