/** @file interpol2.h
 * @brief Header file for interpol2.c.
 */
#ifndef _INTERPOL2_H_

/* We will certainly never use a 16th order scheme,
 * but in case we do, remember to change it here
 */
#define MAX_ORDER 16

#define MAX_INTERPOL_METHODS 9

typedef struct interpol_t interpol_t;
typedef struct interpol_method_t interpol_method_t;

struct interpol_t
{
  double *coeffs;
  double *stencil;

  /* The separation between stencil points. */
  double Lr, Lz;

  /* The interpolation method to use. */
  interpol_method_t *method;

  /* Origin of the interpolation polynomial, which is thus written in
     powers of (r - r0) and (z - z0). */
  double r0, z0;
};

struct interpol_method_t
{
  double *matrix;

  /* p is the order of the polynomial plus one.
   * q is the number of stencil points in one direction (the total
   * number is, of course, q^2. */
  int p, q;

  /* Location of the origin with respect to the stencil in units of L. */
  int ir0, iz0;

  int masses;

  void (*set_coeffs) (interpol_t *this);
  double (*apply) (interpol_t *this, double r, double z);
};

extern interpol_method_t interpol_zero;
extern interpol_method_t interpol_wackers;
extern interpol_method_t interpol_luque;
extern interpol_method_t interpol_quadratic;
extern interpol_method_t interpol_bilin;
extern interpol_method_t interpol_quadlog;
extern interpol_method_t interpol_zero_masses;
extern interpol_method_t interpol_wackers_masses;
extern interpol_method_t interpol_luque_masses;
extern interpol_method_t interpol_quadratic_masses;
extern interpol_method_t interpol_bilin_masses;

extern interpol_method_t *interpol_methods_index[MAX_INTERPOL_METHODS];

#define _INTERPOL2_H_
#endif
