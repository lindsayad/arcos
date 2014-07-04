/** @file poisson.h
 *  @brief Structures and functions for the Poisson solver
 */

#ifndef _POISSON_H_

#include "cstream.h"
#include "rz_array.h"

#ifndef _GRID_H_
#include "grid.h"
#endif

typedef struct pois_grid_t pois_grid_t;

struct pois_grid_t {
  RECT_COORDS;
  LEAF_FIELDS(pois_grid_t);
  int ext_bound;

  rz_array_t *phi;
  rz_array_t *charge;
  rz_array_t *error;
};

/*!< With this structure we define a "problem" for the Poisson/Helmholtz solver.
 * Thus we differentiate between the electrostatic (Poisson) or the
 * photo-ionization (Helmholtz) uses.
 */
typedef struct pois_problem_t pois_problem_t;

struct pois_problem_t {
  int max_level;
  int extra_levels;
  double max_error;

  int bnd_right;
  int bnd_top;
  int bnd_bottom;
};

typedef struct pois_boundaries_t {
  REAL *left, *right, *top, *bottom;
  int r, z;
} pois_boundaries_t;

/*!< Note that the electric fields are computed here as the derivatives
 * of \f$\phi\f$, and not _minus_ the derivatives.
 *
 * This is because our \f$\phi\f$ is not actually the electrostatic potential
 * in its standard definition but its opposite. This simplifies the computations
 * since we can use the charge as the source of the Poisson equation.
 */
#define UNCHECK_ER_RZ(grid_, ir_, iz_)					\
  (((RZ(grid_->phi, (ir_) + 1, iz_) - RZ(grid_->phi, ir_, iz_))		\
    / dr[grid_->level]))

#define UNCHECK_EZ_RZ(grid_, ir_, iz_)					\
  (((RZ(grid_->phi, ir_, (iz_) + 1) - RZ(grid_->phi, ir_, iz_))	\
    / dz[grid_->level]))

/*!< The easiest way (and not too performance costly) is to check the
 *  boundaries where the field can be calculated.
 *
 *  Note that, if everything is OK, the 0 there that we set should never be
 *  used again (and hence any number would be ok (use NaN for debugging) */
#define ER_RZ(grid_, ir_, iz_)			\
  (((ir_) < ((grid_)->r1 + 1) && ((ir_) > ((grid_)->r0 - 2)))? \
   UNCHECK_ER_RZ(grid_, ir_, iz_): 0)

#define EZ_RZ(grid_, ir_, iz_)			\
  (((iz_) < ((grid_)->z1 + 1) && ((iz_) > ((grid_)->z0 - 2)))? \
   UNCHECK_EZ_RZ(grid_, ir_, iz_): 0)

#define _POISSON_H_
#endif
