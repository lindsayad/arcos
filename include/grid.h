/** @file grid.h
 *  @brief Definitions for working with general grids.
 */
#ifndef _GRID_H_

#ifndef _TREE_H_
# include "tree.h"
#endif

#define RECT_COORDS  int r0, r1, z0, z1, ntheta
/*!< Coordinate names for the grids */

#define BND_AT_R        0
#define BND_AT_Z        2
#define BND_MIN         0
#define BND_MAX         1

#define BND_LEFT        (BND_AT_R | BND_MIN)
#define BND_RIGHT       (BND_AT_R | BND_MAX)
#define BND_BOTTOM      (BND_AT_Z | BND_MIN)
#define BND_TOP         (BND_AT_Z | BND_MAX)

#define BND_MASK(B_)    (1 << (B_))
#define BND_NONE             0
#define BND_MASK_LEFT        BND_MASK(BND_LEFT)
#define BND_MASK_RIGHT       BND_MASK(BND_RIGHT)
#define BND_MASK_BOTTOM      BND_MASK(BND_BOTTOM)
#define BND_MASK_TOP         BND_MASK(BND_TOP)
#define BND_MASK_ALL         (BND_MASK_LEFT | BND_MASK_RIGHT\
			      | BND_MASK_BOTTOM | BND_MASK_TOP)

#define GRID_INSIDE    (1 << (BND_TOP + 1))

/* Locations of the cell centers. */
#define r_at(r_, level_) (((double) (r_) + 0.5) * dr[level_])
#define z_at(z_, level_) (((double) (z_) + 0.5) * dz[level_])
#define theta_at(theta_) ((double) (theta_) * dtheta)

/* Locations of the electric field. */
#define er_r_at(i_, level_) ((i_ + 1) * dr[level_])
#define er_z_at(j_, level_) (((double) (j_) + 0.5) * dz[level_])
#define ez_r_at(i_, level_) (((double) (i_) + 0.5) * dr[level_])
#define ez_z_at(j_, level_) ((j_ + 1) * dz[level_])
#define etheta_theta_at(i_) (((double) (i_) + 0.5) * dtheta)


/** We can produce easily a true 2d code by replacing r when it appears
 *  due to the cylindrical coordinates by this macros.  This has a rather
 *  large performance penalty, so it is better to do that in compile time
 */

/* Locations of the electric field. */
#ifdef TRUE2D
# define cyl_q(X_) 1
# define cyl_er_r_at(i_, level_) 1
# define cyl_r_at(r_, level_) 1
#else
# define cyl_q(X_) (X_)
# define cyl_er_r_at(i_, level_) er_r_at(i_, level_) 
# define cyl_r_at(r_, level_) r_at(r_, level_)
#endif

/* These are shortcuts to make the iteration over grid cells easier.
 *  Note that it is important that the inner loop has the smallest
 *  stride to minimize the cache faults.
 */
#define iter_grid_z(grid_, iz_) \
  for(iz_ = grid_->z0; iz_ < grid_->z1; iz_++) 

#define iter_grid_r(grid_, ir_) \
  for(ir_ = grid_->r0; ir_ < grid_->r1; ir_++) 

#define iter_grid_theta(grid_, it_) \
  for(it_ = 0; it_ < grid_->ntheta; it_++) 

#define iter_grid(grid_, ir_, iz_) \
  iter_grid_z(grid_, iz_)	   \
       iter_grid_r(grid_, ir_)

#define iter_grid_3d(grid_, ir_, iz_, it_)		\
  iter_grid_theta(grid_, it_)				\
    iter_grid_z(grid_, iz_)				\
       iter_grid_r(grid_, ir_)


/* Functions to iterate over the grids _including_ 1 extra boundary */
#define iter_grid_z_1(grid_, iz_)			\
  for(iz_ = grid_->z0 - 1; iz_ < grid_->z1 + 1; iz_++) 

#define iter_grid_r_1(grid_, ir_) \
  for(ir_ = grid_->r0 - 1; ir_ < grid_->r1 + 1; ir_++) 

#define iter_grid_1(grid_, ir_, iz_) \
  iter_grid_z_1(grid_, iz_)	   \
       iter_grid_r(grid_, ir_)

#define iter_grid_3d_1(grid_, ir_, iz_, it_)		\
  iter_grid_theta(grid_, it_)				\
    iter_grid_z_1(grid_, iz_)				\
       iter_grid_r_1(grid_, ir_)


/* These include n extra cells in the boundary. We avoid that when 
 * ntheta == 1 because we are in the 2D case and there are no boundaries
 * in theta*/
#define iter_grid_theta_n(grid_, it_, n_)			\
  for(it_ = (grid_->ntheta == 1? 0: -n_); it_ < grid_->ntheta	\
	+ (grid_->ntheta == 1? 0: n_); it_++) 

#define iter_grid_z_n(grid_, iz_, n_)			\
  for(iz_ = grid_->z0 - n_; iz_ < grid_->z1 + n_; iz_++) 

#define iter_grid_r_n(grid_, ir_, n_)			\
  for(ir_ = grid_->r0 - n_; ir_ < grid_->r1 + n_; ir_++) 

#define iter_grid_n(grid_, ir_, iz_, n_)		\
  iter_grid_z_n(grid_, iz_, n_)				\
       iter_grid_r_n(grid_, ir_, n_)

#define iter_grid_3d_n(grid_, ir_, iz_, it_, n_)	\
  iter_grid_theta_n(grid_, it_, n_)			\
     iter_grid_z_n(grid_, iz_, n_)			\
       iter_grid_r_n(grid_, ir_, n_)


/* Iter through the parent indices that fill one grid. */
#define iter_grid_parent(grid_, ir_, iz_) \
  for(ir_ = (grid_->r0 >> 1) + (grid_->r0 % 2);				\
      ir_ < (grid_->r1 >> 1); ir_++)					\
    for(iz_ = (grid_->z0 >> 1) + (grid_->z0 % 2);			\
	iz_ < (grid_->z1 >> 1); iz_++)


/** @brief The cdr grids and poisson grids types are cdr_grid_t and
 *  pois_grid_t, defined in cdr.h and poisson.h.
 *
 *  This structure contains the common parts of both structures and
 *  thus can be used in general routines that accept any kind of grid.
 *  The C standard promises that we can du this.
 */
typedef struct grid_t grid_t;
struct grid_t {
  RECT_COORDS;
  LEAF_FIELDS(grid_t);
  /* ext_bound indicates which boundaries are external boundaries.
     e.g. ext_bound & BND_LEFT_TOP tells us whether the left-top boundary
     is an external one. */
  int ext_bound;
};

/** @brief To ease typing, you can write e.g.
 * printf ("Solving grid (" grid_printf_str ")\n", grid_printf_args(grid));
 */
#define grid_printf_str "{r0 = %d, z0 = %d, r1 = %d, z1 = %d, level = %d}"
#define grid_printf_args(G_) (G_)->r0, (G_)->z0, (G_)->r1, (G_)->z1,	\
    (G_)->level

#define _GRID_H_
#endif
