/** @file sprites.c
 *  @brief Routines for the "sprites module."
 *
 * The idea here is to cope with strongly varying densities along
 * the streamer propagation.
 *
 * For photoionization we use a "local maximum approximation", i.e. the
 * quenching and absorption lengths for all the volume are taken as those
 * corresponding to the streamer head.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "cdr.h"
#include "parameters.h"
#include "photo.h"
#include "proto.h"
#include "rz_array.h"
#include "species.h"

static int var_max_axis (cdr_grid_t *grid, rz_array_t *var, int sign);

double spr_nu_a, spr_nu_d, spr_nu_rt;
photo_term_t *spr_photo_terms;
extern photo_term_t *photo_terms;

/** @brief The density at a given height, assuming dens = 1.0 at height = 0.0.
 *
 * We are assuming an exponential profile here but in principle any other
 * profile would also do.
 */
double
spr_density_at (double altitude)
{
  return exp (-altitude / dens_decay_len);
}


/** @brief Initialize the sprite module.
 *
 * We have to remember the values of the global variables that we will adapt
 * to altitude.
 */
void
spr_init ()
{
  /* We make the Courant numbers depend on the altitude: this is because
     when we apply the Courant-Levy criterium we assume v = -E, which is not
     true if the mobility now is != 1.0 */
  spr_nu_a = nu_a;
  spr_nu_d = nu_d;
  spr_nu_rt = nu_rt;

  photo_copy_list (photo_terms, &spr_photo_terms);
}

/** @brief Hook to update the variables that depend on the head altitude.
 *
 * Has to be called AFTER the charge is calculated.
 */
void
spr_hook (cdr_grid_t *grid)
{
  double altitude;
  
  altitude = spr_head_altitude (grid, sprite_sign);

  spr_update (altitude);
}

/** @brief Updates the magnitudes that we are tracking according to a
 * given altitude. 
 */
void
spr_update (double altitude)
{
  double back_dens = spr_density_at (altitude);
  photo_term_t *p0, *p;

  nu_a = spr_nu_a;
  nu_d = spr_nu_d * back_dens;
  nu_rt = spr_nu_rt;

  for (p0 = spr_photo_terms, p = photo_terms; p0; p0 = p0->next, p = p->next) {
    p->lambda = p0->lambda * back_dens * back_dens;
    p->A = p0->A * (sprite_dens_0 + sprite_dens_q) 
      / (back_dens + sprite_dens_q);
  }
  
}

/** @brief Finds the position of the streamer head by locating the maximum
 * (or -maximum) of the charge along the streamer axis.
 *
 * @a sign has to be +1 for positive streamers and -1 for negative streamers.
 */
double
spr_head_altitude (cdr_grid_t *grid, int sign)
{
  int ih;
  ih = var_max_axis (grid, grid->charge, sign);
  if (0 == ih) ih = var_max_axis (grid, grid->dens[electrons], 1);
  return z_at (ih, grid->level);
}

/** @brief var_max_axis QQQQ */
static int
var_max_axis (cdr_grid_t *grid, rz_array_t *var, int sign)
{
  int iz, iz_max = 0;
  double var_max = 0.0;

  iter_grid_z (grid, iz) {
    if (sign * RZT (var, 0, iz, 0) > var_max) {
      iz_max = iz;
      var_max = RZT (var, 0, iz, 0);
    }
  }
  
  return iz_max;
}
