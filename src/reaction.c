/** @file reaction.c
 *  @brief Functions to handle the "reaction" part of
 *  convection-diffusion-reaction equation.
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "cdr.h"
#include "grid.h"
#include "parameters.h"
#include "photo.h"
#include "proto.h"
#include "react_table.h"
#include "species.h"
#include "reaction.h"

static void fill_react_gaps ();

reaction_t *reactions_list = NULL;
extern double z_cutoff;
double species_written[20];

/** @brief Returns the position in the species-array of a given species.
 *
 * Returns -1 if not found. */
int
find_species_by_name(const char *spec_name)
{
  int res = -1;
  int cnt;

  for (cnt = 0; cnt < no_species; cnt++)
  {
    if (strcmp(spec_index[cnt]->name, spec_name) == 0) 
      res = cnt;
  }

  if (res == -1) {
      printf("Species-lookup-failure: %s!\n",spec_name);
      assert(1 == 0);
  }
  return res;
}

/** @brief Adds a given reaction to the reaction list */
void
react_add (reaction_t *react)
{
  react_table *rt;

  debug (3, "react_add (...)\n");

  printf("Adding reaction with table: %s\n",react->tablefile);

  rt = (react_table *) xmalloc (sizeof(react_table));

  if (react->tablefile != NULL) {
    react_table_read(react->tablefile, rt);
    react->rt = rt;
  }

  react->next = reactions_list;
  reactions_list = react;
}

/** @brief Applies a reaction to the given @a grid */
void 
react_apply (reaction_t *react, cdr_grid_t *grid, int overwrite)
{
  int i, ir, iz, itheta;
  rz_array_t *grid_in[REACTION_MAX_IN], *grid_out[REACTION_MAX_IN + REACTION_MAX_OUT],
    *grid_eabs;
  double *in = NULL, *out = NULL, eabs;
  double test;
  double rate;
  int pos;
  double e1, e2, val1, val2, log_e, res, r_mod;
  int cnt, curr_species;

  debug (3, "react_apply (..., " grid_printf_str ")\n",
         grid_printf_args (grid));

  grid_eabs = grid->dens[no_species];

  for (i = 0; i < react->nin; i++) {
    grid_in[i] = grid->dens[react->input[i]];
    grid_out[i] = grid->d_dens[react->input[i]];
  }

  for (i = react->nin; i < react->nin + react->nout; i++) {
    grid_out[i] = grid->d_dens[react->output[i - react->nin]];
  }

#pragma omp parallel private(ir, iz, i, in, out)
  {
    /* malloc(0) is legal, but I do not want to play with fire. */
    if (react->nin > 0) {
      in = (double *) xmalloc (sizeof(double) * react->nin);
    }

    /* Do not know what use nout == 0 may have, (perhaps some debugging?)
       but we leave it here as theoretically possible. */
    if (react->nout > 0) {
      out = (double *) xmalloc (sizeof(double) * (react->nout + react->nin));
    }


    #pragma omp for
    iter_grid_theta(grid, itheta) { //ITER3
      iter_grid_z(grid, iz) { //ITER2
        double back_dens;
        if (sprite_module) {
          back_dens = spr_density_at (z_at (iz, grid->level));
        } else {
          back_dens = 1.0;
        }


        iter_grid_r(grid, ir) { //ITER1
          if ( z_at (iz, grid->level) < z_cutoff ) { //IF1
            eabs = fabs (RZT(grid_eabs, ir, iz, itheta));
            for (i = 0; i < react->nin; i++) {
              in[i] = fabs (RZT(grid_in[i], ir, iz, itheta));
	    }

            log_e = log10(eabs);

            /* If the supplied fieldstrength falls outside the boundaries of the table,
              return predetermined under-/overflow values */
            if (log_e < react->rt->e_min) {
              rate = react->rt->underflow;
            } else if (log_e > react->rt->e_min + react->rt->e_step * react->rt->steps) {
              rate = react->rt->overflow;
            } else {
              pos = floor((log_e - react->rt->e_min) / react->rt->e_step);
              val1 = react->rt->values[pos];
              val2 = react->rt->values[pos+1];
              e1 = pow(10, react->rt->e_min + react->rt->e_step * (double) pos);
              e2 = e1 * pow(10, react->rt->e_step);

              rate = val1 + ((val2 - val1) / (e2 - e1)) * (eabs - e1);
            }

            for (i = 0; i < react->nin; i++){ rate *=  MYMAX(0, in[i]); }

            for (i = 0; i < react->nin + react->nout; i++) { //FOR1
              if (i < react->nin) {
                curr_species = react->input[i]; r_mod = -rate;
              } else { curr_species = react->output[i - react->nin]; r_mod = rate; }
    
              if (spec_index[curr_species]->charge != 0.0) {
                 RZT(grid_out[i], ir, iz, itheta) += r_mod;
              }
            } //FOR1
          } //IF1
        } //ITER1
      } //ITER2
    } //ITER3
    free (in);
    free (out);
  }
}

/** @brief Recursive version of @a react_apply  */
void
react_apply_r (reaction_t *react, cdr_grid_t *grid, int overwrite)
{
  cdr_grid_t *child;

  react_apply (react, grid, overwrite);

  iter_childs (grid, child) {
    react_apply_r (react, child, overwrite);
  }
}

/** @brief Sets the d_dens field of a grid to zero */
void 
zero_fill (cdr_grid_t* grid)
{
  int ir, iz, itheta, i;

  iter_grid_3d (grid, ir, iz, itheta)
  for (i = 0; i < no_species; i++) {
    RZT(grid->d_dens[i], ir, iz, itheta) = 0.0;
  }
}

/** @brief Applies all reactions to the given @a grid and his descendants */
void 
react_apply_all (cdr_grid_t *grid)
{
  reaction_t *react;
  int overwrite;
  int last = -1;
  int cnt;

  zero_fill(grid);

  overwrite = TRUE;

  for (react = reactions_list; react; react = react->next) {
    if (react->is_photo) {
      photo_calc (photo_terms, grid);
    } else {
      react_apply_r (react, grid, overwrite);
    }
    overwrite = FALSE;
  }
}

/** @brief Fill in the gaps that we left in the definitions */
static void
fill_react_gaps ()
{
}

/** @brief Initializes the list of reactions. */
void 
react_init ()
{
  /* Note that the reactions are applied in inverse order than listed here. */
  /* Rest of the kinetic model: */
  kinetic_init ();
 
}

/** Below this electric field, we do not waste time calculating anything.
   Besides, this avoid NaNs for eabs == 0. */
#define EPS_EABS 1e-6
