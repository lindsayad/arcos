/** @file cdr.h
 *  @brief Header file for cdr grids. */

#include "cstream.h"
#include "tree.h"

#ifndef _GRID_H_
#include "grid.h"
#endif

#ifndef  _RZ_ARRAY_H_
#include "rz_array.h"
#endif

#ifndef  _TREE_H_
#include "tree.h"
#endif

#ifndef _CDR_H_
typedef struct cdr_grid_t cdr_grid_t;

struct cdr_grid_t {
  RECT_COORDS;
  LEAF_FIELDS(cdr_grid_t);
  int ext_bound;

  /* Pointer to each of the species. And their time derivatives */
  rz_array_t **dens;
  rz_array_t **d_dens;

  /* Components and magnitude of the electric field */
  rz_array_t *er, *ez, *etheta, *eabs;
  rz_array_t *charge;

  rz_array_t *photo;

  REAL *max_dens;
  REAL max_charge;
  REAL max_eabs;

  int contains_edge;
};

#define SET_DENS_OVERWRITE 1
#define SET_DENS_ADD 2
#define SET_DENS_SUB 3

#define _CDR_H_
#endif
