/** @file parameters.h
 *  @brief  The declarations of all global parameters, i.e. those
 *          that the user should be able to set.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libconfig.h>
#include <stdbool.h>

#ifndef _CSTREAM_H_
# include "cstream.h"
#endif

#ifndef _RZ_ARRAY_H_
# include "rz_array.h"
#endif

#ifndef _PARAMETERS_H_

/** @brief These are the global parameters. */
char* prog_id;
char* output_dir;
char* kin_input;
char* load_file;
char* photoionization_file;

int   cdr_bnd_bottom,cdr_bnd_right,cdr_bnd_top,
      cdr_brick_dr,cdr_brick_dz,cdr_max_level,cdr_interp_bnd,cdr_interp_in,
      extra_photo_levels,extra_pois_levels,extra_photo_levels_2,
      gridpoints_r,gridpoints_z,
      has_photoionization,
      max_disk_space,max_disk_space_mb,max_levels,max_ntheta,
      perturb_max_k,photo_bnd_bottom,photo_bnd_bottom_2,
      photo_bnd_right,photo_bnd_right_2,
      photo_bnd_top,photo_bnd_top_2,photo_bnd_top_2_st,
      photo_max_level,photo_max_level_2,
      pois_bnd_bottom,pois_bnd_top,pois_bnd_right,
      pois_inhom,pois_inhom_reflections,pois_max_level,
      cdr_output_margin,pois_output_margin,
      restart,ref_level_eabs,spec_total,sprite_module,sprite_sign;

double attachment_rate,attachment_E0,
       attempt_dt,diffusion_coeff,
       background_ionization,background_increase_length,
       constant_source,dens_decay_len,
       end_t,E0_x,E0_y,E0_z,L_r,L_z,
       needle_length,needle_radius,nu_a,nu_d,nu_f,nu_rt,
       off_time,output_dt,
       perturb_epsilon,
       photo_max_error,photo_max_error_2,
       pois_inhom_fixed_q,pois_max_error,
       ref_threshold_charge,ref_threshold_dens,ref_threshold_eabs,
       ref_threshold_edge,
       rise_time,start_t,
       seed_sigma_x,seed_sigma_y,seed_sigma_z,seed_N,
       sprite_dens_0,sprite_dens_q,
       warn_min_timestep;
int   pois_output;

#define _PARAMETERS_H_
#endif
