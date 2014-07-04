/** @file mapper.h
 *  @brief Definitions for mapper objects.
 *
 * A mapper objects defines how and what to map between two grid families.
 * For example, a mapper object can be used to map the charge from
 * a @a cdr grid to a @a Poisson grid, or to map the @a potential of a Poisson
 * grid into one of the components of the @a electric field in a @a cdr grid.
 * For that, a mapper object has to implement some methods that tell us
 * how to map things in the cases\n
 * a) the source grid is coarser than the target grid (interpolation),\n
 * b) the source and the target are at the same level (copy),\n
 * c) the target grid is coarser (coarsening).
 */

#ifndef _MAPPER_H_

typedef struct mapper_t mapper_t;

#ifndef _INTERPOL2_H_
# include "interpol2.h"
#endif

#define _coarsen_params_ (mapper_t *mapper, grid_t *source, grid_t *target, \
			  int ir, int iz, int itheta)
#define _copy_params_ (mapper_t *mapper, grid_t *source, grid_t *target, \
		       int ir, int iz, int itheta)
#define _interpol_set_params_ (mapper_t *mapper, grid_t *source, \
			       interpol_t *interpol, int ir, \
			       int iz, int itheta)
#define _interpol_params_ (mapper_t *mapper, grid_t *source, grid_t *target, \
			   interpol_t *interpol, int ir, int iz, int itheta)


struct mapper_t {
  interpol_method_t *interpol_method;

  void (*coarsen) _coarsen_params_;
  void (*copy) _copy_params_;
  /* may return false if we are outside the source grid. */
  int (*interpol_set) _interpol_set_params_;

  void (*interpol) _interpol_params_;

  /* This is used for the interpolation of electric fields and indicates
   * the staggering of the cells.  When we call the interpolator to set
   * the stencil at coordinates pr, pz, the fine-grid cells whose value will
   * later be calculated ir, iz with
   *   pr << level_diff + shift_r << (level_diff - 1) <= ir 
   * < pr << level_diff + shift_r << (level_diff - 1)
   * 
   * (and the equivalent for iz). 
   * Hence note that for interpolation in not-staggered grids, 
   * (i.e. charge, densities, ...)
   * shift_r = shift_z = 0.
   */
  int shift_r, shift_z;

  /* We provide this extra member to allow the easy creation of different
   * mappers that share the same functions but have a slightly different
   * behaviour.  In particular this is used in cdr.c to define
   * different interpolators for each species.  In that case, extra is the
   * species index.
   */
  int extra;
};


/** Usually, we will define functions called xxxx_coarsen, xxxx_copy, etc
   to implement the above interface.  We define these macros to facilitate
   the declaration of such functions. */
#define decl_mapper_funcs(_VAR)				\
  void _VAR ## _coarsen _coarsen_params_;		\
  void _VAR ## _copy _copy_params_;			\
  int _VAR ## _interpol_set _interpol_set_params_;	\
  void _VAR ## _interpol _interpol_params_

/** Useful to init a staggered mapper object */
#define mk_mapper_staggered(_C, _INTERPOL_METHOD, _SHIFT_R, _SHIFT_Z)	\
  {_INTERPOL_METHOD, _C ## _coarsen, _C ## _copy,			\
      _C ## _interpol_set, _C ##_interpol, _SHIFT_R, _SHIFT_Z, 0}

/** Useful to init a mapper object */
#define mk_mapper(_C, _INTERPOL_METHOD)					\
  mk_mapper_staggered (_C, _INTERPOL_METHOD, 0, 0)

/** Mappers that only copy or interpolate. */
#define mk_mapper_down(_C, _INTERPOL_METHOD)				\
  {_INTERPOL_METHOD, NULL, _C ## _copy,					\
      _C ## _interpol_set, _C ##_interpol, 0, 0, 0}

#define _MAPPER_H_
#endif
