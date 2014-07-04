/** @file configuration.c
 *  @brief Module for input/output of parameters
 */
/* ----------------------------------------------------------------------------
   libconfig - A library for processing structured configuration files
   Copyright (C) 2005-2010  Mark A Lindner

   This file is part of libconfig.

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public License
   as published by the Free Software Foundation; either version 2.1 of
   the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Library General Public
   License along with this library; if not, see
   <http://www.gnu.org/licenses/>.
   ----------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <libconfig.h>

#include <cstream.h>
#include <parameters.h>
#include <configuration.h>
#include <reaction.h>
#include <species.h>

/** @brief Looks in a libconfig file *.cfg for a parameter of type 'int', with 
 *  name 'sstring' and add its value to the global variable with the same name.
 *  If already initialized, its value will be overwritten.
 *
 *  @param[in] name		name of the parameter (see e.g., input file default.cfg)
 *  @param[in] type		type of parameter     (see e.g., input file default.cfg)
 *  @param[in] comment		comment to be added   (see e.g., input file default.cfg)
 *  @param[in] sstring  	name to be found
 *  @param[out] *par		global 'int' variable named 'name';
 *                              on return *par will have value 'value'
 *  @param[in] value		of type 'int'
 *  @param[in] count		number of elements in array 'parameter_names' to examine
 *  @param[in] *change		if TRUE then value in setting_default will be changed
 *  @param[in] setting_default	related to default configuration file
 *  @param[out] i		position in global array 'parameter_names' 
 *                              defined in configuration.h
 *  
 */
bool
change_cfg_integer(const char* name,
		   const char *type,
		   const char *comment,
                   const char *sstring,
		   int *par,
		   int value,
		   int count,
                   bool *change,
		   config_setting_t *setting_default,
		   int i)
{
  int j;
  config_setting_t *elem,*setting_value,*setting_comment;

  if(!(strncmp(name,sstring,COMPARE_LIMIT))) {
    if(!(strncmp(type,"int",COMPARE_LIMIT))) {
      *par = value;
      if(change) {
        /* Find the corresponding index in the file default.cfg */
        for(j=0; j<count; j++) {
          if(!(strncmp(name,parameter_names[j].name,COMPARE_LIMIT))) {
            i=j;
            break;
          }
        }
        /* Get access to the right entry in file default.cfg */
        elem            = config_setting_get_elem(setting_default,i);
        setting_value   = config_setting_get_member(elem,"value");
        setting_comment = config_setting_get_member(elem,"comment");
        /* Adapt the value in the configuration file default.cfg */
        if(setting_value && setting_comment) {
          config_setting_set_int   (setting_value,value);
          config_setting_set_string(setting_comment,comment);
        }
      } else
        parameter_names[i].name   = name;
      return(TRUE);
    } else {
      printf("Warning: wrong type for %s; type =%s\n",name,type);
      printf("Warning: type must be INT for %s\n",sstring);
    }
  } 
  return(FALSE);
}

/** @brief Looks in a libconfig file *.cfg for a parameter of type 'double', with 
 *  name 'sstring' and add its value to the global variable with the same name.
 *  If already initialized, its value will be overwritten.
 *
 *  @param[in] name		name of the parameter (see e.g., input file default.cfg)
 *  @param[in] type		type of parameter     (see e.g., input file default.cfg)
 *  @param[in] comment		comment to be added to name (see e.g., input file default.cfg)
 *  @param[in] sstring  	name to be found
 *  @param[out] *par		global 'double' variable named 'name';
 *                              on return *par will have value 'value'
 *  @param[in] value		of type 'double'
 *  @param[in] count		number of elements in array 'parameter_names' to examine
 *  @param[in] change		if TRUE then value in setting_default will be changed
 *  @param[in] setting_default	related to default configuration file
 *  @param[out] i		position in global array 'parameter_names' 
 *                              defined in configuration.h
 *  
 */
bool
change_cfg_double(const char* name,
		  const char *type,
		  const char *comment,
                  const char *sstring,
		  double *par,
		  double value,
		  int count,
                  bool *change,
		  config_setting_t *setting_default,
		  int i)
{
  int j;
  config_setting_t *elem,*setting_value,*setting_comment;

  if(!(strncmp(name,sstring,COMPARE_LIMIT))) {
    if(!(strncmp(type,"double",COMPARE_LIMIT))) {
      *par = value;
      if(change) {
        /* Find the corresponding index in the file default.cfg */
        for(j=0; j<count; j++) {
          if(!(strncmp(name,parameter_names[j].name,COMPARE_LIMIT))) {
            i=j;
            break;
          }
        }
        /* Get access to the right entry in file default.cfg */
        elem            = config_setting_get_elem(setting_default,i);
        setting_value   = config_setting_get_member(elem,"value");
        setting_comment = config_setting_get_member(elem,"comment");

        /* Adapt the value in the configuration file default.cfg */
        if(setting_value && setting_comment) {
          config_setting_set_float(setting_value,value);
          config_setting_set_string(setting_comment,comment);
        }
      } else
        parameter_names[i].name   = name;
      return(TRUE);
    } else {
      printf("Warning: wrong type for %s; type =%s\n",name,type);
      printf("Warning: type must be DOUBLE for %s\n",sstring);
    }
  } 
  return(FALSE);
}
/** @brief
 *
 * Looks in a libconfig file  *.cfg for a parameter of type 'string', with 
 *  name 'sstring' and add its value to the global variable with the same name.
 *  If already initialized, its value will be overwritten.
 *
 *  @param[in] name		name of the parameter (see e.g., input file default.cfg)
 *  @param[in] type		type of parameter     (see e.g., input file default.cfg)
 *  @param[in] comment		comment to be added to name
 *                             (see e.g., input file default.cfg)
 *  @param[in] sstring  	name to be found
 *  @param[out] par		global variable named 'name'; on return par will have value 'value'
 *  @param[in] value		of type 'string'
 *  @param[in] count		number of elements in array 'parameter_names' to examine
 *  @param[in] change		if TRUE then value in setting_default will be changed
 *  @param[in] setting_default	related to default configuration file
 *  @param[out] i		position in global array 'parameter_names' 
 *                              defined in configuration.h
 *
 */  
bool
change_cfg_string(const char *name,
	  	  const char *type,
		  const char *comment,
                  const char *sstring,
		  char** par,
		  char* value,
		  int count,
                  bool *change,
		  config_setting_t *setting_default,
		  int i)
{
  int j;
  config_setting_t *elem,*setting_value,*setting_comment;

  if(!(strncmp(name,sstring,COMPARE_LIMIT))) {
    if(!(strncmp(type,"string",COMPARE_LIMIT))) {
      *par = value;
      if(change) {
/* Find the corresponding index in the file default.cfg */
        for(j=0; j<count; j++) {
          if(!(strncmp(name,parameter_names[j].name,COMPARE_LIMIT))) {
            i=j;
            break;
          }
        }
        /* Get access to the right entry in file default.cfg */
        elem            = config_setting_get_elem(setting_default,i);
        setting_value   = config_setting_get_member(elem,"value");
        setting_comment = config_setting_get_member(elem,"comment");

        /* Adapt the value in the configuration file default.cfg */
        if(setting_value && setting_comment) {
          config_setting_set_string(setting_value,value);
          config_setting_set_string(setting_comment,comment);
        }
      } else
        parameter_names[i].name   = name;
      return(TRUE);
    } else {
      printf("Warning: wrong type for %s; type =%s\n",name,type);
      printf("Warning: type must be STRING for %s\n",sstring);
    }
  } 
  return(FALSE);
}

/** @brief File that searches for all parameters in configuration file connected
 *  to setting.
 *
 *  Calls for each parameter the right change_cfg_... function
 */
 /*
  * !!!! DOXYGEN ???
 *  @param[in] name				name of the parameter (see e.g., input file default.cfg)
 *  @param[in] type				type of parameter     (see e.g., input file default.cfg)
 *  @param[in] comment				comment to be added to name (see e.g., input file default.cfg)
 *  @param[out] ii				position in global array 'parameter_names' 
 *  @param[out] ivalue,dvalue,bbool,astring	name to be found, may be a 'int','double','bool'
 *                                              or 'char*'
 *  @param[in] config_setting_t	*setting	represents s configuration setting
 *  @param[in] count				number of elements in array 'parameter_names' to examine
 *  @param[in] *change				if TRUE then value in setting_default will be changed
 *  
 */
void
change_cfg_parameters (const char *name,
		       const char *type,
		       const char *comment,
                       int ii,
		       int ivalue,
		       double dvalue,
		       bool bbool,
                       const char* astring,
		       config_setting_t *setting,
                       int count,
		       bool *change)
{
  /* An identification name for this run, prog_id */
  if(change_cfg_string(name,type,comment,"prog_id",&prog_id,(char*)astring,count,change,setting,ii))
    return;

  /* Output directory, output_dir */
  if(change_cfg_string(name,type,comment,"output_dir",&output_dir,(char*)astring,count,change,setting,ii))
    return;

  /* Kinetics input file, kin_input */
  if(change_cfg_string(name,type,comment,"kin_input",&kin_input,(char*)astring,count,change,setting,ii))
    return;

  /* Restart from a previous file */
  if(change_cfg_integer(name,type,comment,"restart",&restart,ivalue,count,change,setting,ii))
    return;

  /* The simulation will continue from file */
  if(change_cfg_string(name,type,comment,"load_file",&load_file,(char*)astring,count,change,setting,ii))
    return;

  /* The name of a file with photoionization parameters */
  if(change_cfg_string(name,type,comment,"photoionization_file",&photoionization_file,
                (char*)astring,count,change,setting,ii))
    return;

  /* Time interval for output to be written to disk */
  if(change_cfg_double(name,type,comment,"output_dt",&output_dt,dvalue,count,change,setting,ii))
    return;

  /* Output of the Poisson grids, including the potential */
  if(change_cfg_integer(name,type,comment,"pois_output",&pois_output,ivalue,count,change,setting,ii))
    return;

  /* Margin outside the grids in the output of the cdr equation */
  if(change_cfg_integer(name,type,comment,"cdr_output_margin",&cdr_output_margin,ivalue,count,change,setting,ii))
    return;

  /* Output of the Poisson grids, including the potential */
  if(change_cfg_integer(name,type,comment,"photo_bnd_right",&photo_bnd_right,ivalue,count,change,setting,ii))
    return;

  /* Margin outside the grids in the output of the Poisson equation */
  if(change_cfg_integer(name,type,comment,"pois_output_margin",&pois_output_margin,ivalue,count,change,setting,ii))
    return;

  /* Sometimes, we can issue warnings about strange things happening in the 
   * program.  Most of the time something is going wrong and you should
   * revise your parameters.
   */
  if(change_cfg_double(name,type,comment,"warn_min_timestep",&warn_min_timestep,dvalue,count,change,setting,ii))
    return;

  /* Imposing a limit for the disk space is hard but it is the law.
   *
   * The default is 1 Tb -> You should use something smaller. 
   */
  if(change_cfg_integer (name,type,comment,"max_disk_space_mb",&max_disk_space_mb,ivalue,count,change,setting,ii))
    return;

  /************************ 
   * Numerical parameters *
   ************************/

  /* Number of gridpoints in r direction at level 0. */
  if(change_cfg_integer (name,type,comment,"gridpoints_r",&gridpoints_r,ivalue,count,change,setting,ii))
    return;

  /* Number of gridpoints in z direction at level 0. */
  if(change_cfg_integer (name,type,comment,"gridpoints_z",&gridpoints_z,ivalue,count,change,setting,ii))
    return;

  /* Number of azimuthal modes. */
  if(change_cfg_integer (name,type,comment,"max_ntheta",&max_ntheta,ivalue,count,change,setting,ii))
    return;

  /* Initial time may be away from 0 */
  if(change_cfg_double(name,type,comment,"start_t",&start_t,dvalue,count,change,setting,ii))
    return;

  /* End time of simulation, originally 500.0 */
  if(change_cfg_double(name,type,comment,"end_t",&end_t,dvalue,count,change,setting,ii))
    return;

  /* The actual timestep may be smaller that this one, if it is needed
   * to satisfy the Courant constraint.
   */
  if(change_cfg_double(name,type,comment,"attempt_dt",&attempt_dt,dvalue,count,change,setting,ii))
    return;

  /* Number of levels that are used for the Poisson equation
   * but NOT for the cdr integrator.
   */
  if(change_cfg_integer (name,type,comment,"extra_pois_levels",&extra_pois_levels,ivalue,count,change,setting,ii))
    return;

  /* Maximum nesting depth */
  if(change_cfg_integer (name,type,comment,"max_levels",&max_levels,ivalue,count,change,setting,ii))
    return;

  /* Total number of species */
  if(change_cfg_integer (name,type,comment,"no_species",&no_species,ivalue,count,change,setting,ii))
    return;

  /* Maximum error allowed in the Poisson solver.
   * Cells with larger errors will be further refined.
   */
  if(change_cfg_double(name,type,comment,"pois_max_error",&pois_max_error,dvalue,count,change,setting,ii))
    return;

  /* Maximum refinement level for the Poisson solver. */
  if(change_cfg_integer(name,type,comment,"pois_max_level",&pois_max_level,ivalue,count,change,setting,ii))
    return;

  /* These are photo-ionization parameters equivalent to the Poisson ones.
   *  If extra_photo_levels < 0, then these parameters are ignored and the
   *  poisson-parameters are used. If extra_photo_levels_2 < 0 then the
   *  parameters for the second term are ignored and the parameters for the
   *  first term are used for both terms.
   */

  /* First photo-ionization term */
  if(change_cfg_integer (name,type,comment,"extra_photo_levels",&extra_photo_levels,ivalue,count,change,setting,ii))
    return;

  /* Maximum level of refinement in the photo-ionization solver.*/
  if(change_cfg_integer (name,type,comment,"photo_max_level",&photo_max_level,ivalue,count,change,setting,ii))
    return;

  /* Error threshold that leads to refinement in the photo-ionization code. */
  if(change_cfg_double(name,type,comment,"photo_max_error",&photo_max_error,dvalue,count,change,setting,ii))
    return;

  /* Photo-ionization boundary condition at r = L_r
       *  1 means homogeneous Neumann boundary conditions
       * -1 means homogeneous Dirichlet boundary conditions.
       */
  if(change_cfg_integer(name,type,comment,"photo_bnd_right" ,&photo_bnd_right ,ivalue,count,change,setting,ii))
    return;

  if(change_cfg_integer(name,type,comment,"photo_bnd_bottom",&photo_bnd_bottom,ivalue,count,change,setting,ii))
    return;

  if(change_cfg_integer(name,type,comment,"photo_bnd_top"   ,&photo_bnd_top   ,ivalue,count,change,setting,ii))
    return;


  /* Second photo-ionization term: Extra levels for the photo-ionization solver. */
  if(change_cfg_integer(name,type,comment,"extra_photo_levels_2",&extra_photo_levels_2,ivalue,count,change,setting,ii))
    return;

  /* Maximum level of refinement in the photo-ionization solver.*/
  if(change_cfg_integer(name,type,comment,"photo_max_level_2",&photo_max_level_2,ivalue,count,change,setting,ii))
    return;


  /* Error threshold that leads to refinement in the photo-ionization code.*/
  if(change_cfg_double(name,type,comment,"photo_max_error_2",&photo_max_error_2,dvalue,count,change,setting,ii))
    return;

  /* Photo-ionization boundary condition at r = L_r.
   *  1 for homogeneous Neumann,
   * -1 for homogeneous Dirichlet
   */
  if(change_cfg_integer(name,type,comment,"photo_bnd_right_2",&photo_bnd_right_2,ivalue,count,change,setting,ii))
    return;

  /* Photo-ionization boundary condition at  z = 0.
   *  1 for homogeneous Neumann,
   * -1 for homogeneous Dirichlet
   */
  if(change_cfg_integer(name,type,comment,"photo_bnd_bottom_2",&photo_bnd_bottom_2,ivalue,count,change,setting,ii))
    return;

  /* Photo-ionization boundary condition at  z = L_z.
   *  1 for homogeneous Neumann,
   * -1 for homogeneous Dirichlet
   */
  if(change_cfg_integer(name,type,comment,"photo_bnd_top_2",&photo_bnd_top_2,ivalue,count,change,setting,ii))
    return;

  /* Boundary conditions for the CDR system. */
  if(change_cfg_integer(name,type,comment,"cdr_bnd_bottom",&cdr_bnd_bottom,ivalue,count,change,setting,ii))
    return;

  if(change_cfg_integer(name,type,comment,"cdr_bnd_top",&cdr_bnd_top,ivalue,count,change,setting,ii))
    return;

  if(change_cfg_integer(name,type,comment,"cdr_bnd_right",&cdr_bnd_right,ivalue,count,change,setting,ii))
    return;

  if(change_cfg_integer(name,type,comment,"pois_bnd_right",&pois_bnd_right,ivalue,count,change,setting,ii))
    return;

  if(change_cfg_integer(name,type,comment,"pois_bnd_bottom",&pois_bnd_bottom,ivalue,count,change,setting,ii))
    return;

  if(change_cfg_integer(name,type,comment,"pois_bnd_top",&pois_bnd_top,ivalue,count,change,setting,ii))
    return;

  /* Courant numbers. */
  if(change_cfg_double(name,type,comment,"nu_a",&nu_a,dvalue,count,change,setting,ii))
    return;

  if(change_cfg_double(name,type,comment,"nu_d",&nu_d,dvalue,count,change,setting,ii))
    return;

  if(change_cfg_double(name,type,comment,"nu_rt",&nu_rt,dvalue,count,change,setting,ii))
    return;

  if(change_cfg_double(name,type,comment,"nu_f",&nu_f,dvalue,count,change,setting,ii))
    return;

  /* Refinement criteria for the CDR equation. */
  if(change_cfg_double(name,type,comment,"ref_threshold_eabs",&ref_threshold_eabs,dvalue,count,change,setting,ii))
    return;

  /* Refinement criteria for the CDR equation. */
  if(change_cfg_integer(name,type,comment,"ref_level_eabs",&ref_level_eabs,ivalue,count,change,setting,ii))
    return;

  if(change_cfg_double(name,type,comment,"ref_threshold_charge",&ref_threshold_charge,dvalue,count,change,setting,ii))
    return;

  if(change_cfg_double(name,type,comment,"ref_threshold_dens",&ref_threshold_dens,dvalue,count,change,setting,ii))
    return;

  if(change_cfg_double(name,type,comment,"ref_threshold_edge",&ref_threshold_edge,dvalue,count,change,setting,ii))
    return;

  /* r- and z-length of the minimal refinement area in the cdr equation*/
  if(change_cfg_integer(name,type,comment,"cdr_brick_dr",&cdr_brick_dr,ivalue,count,change,setting,ii))
    return;

  if(change_cfg_integer(name,type,comment,"cdr_brick_dz",&cdr_brick_dz,ivalue,count,change,setting,ii))
    return;

  /* Maximum refinement level for the cdr solver. */
  if(change_cfg_integer(name,type,comment,"cdr_max_level",&cdr_max_level,ivalue,count,change,setting,ii))
    return;

  /* Interpolation method. */
  if(change_cfg_integer(name,type,comment,"cdr_interp_in",&cdr_interp_in,ivalue,count,change,setting,ii))
    return;

  if(change_cfg_integer(name,type,comment,"cdr_interp_bnd",&cdr_interp_bnd,ivalue,count,change,setting,ii))
    return;

  /*********************** 
   * Physical parameters *
   ***********************/
  if(change_cfg_double(name,type,comment,"L_r",&L_r,dvalue,count,change,setting,ii))
    return;

  if(change_cfg_double(name,type,comment,"L_z",&L_z,dvalue,count,change,setting,ii))
    return;

  if(change_cfg_double(name,type,comment,"diffusion_coeff",&diffusion_coeff,dvalue,count,change,setting,ii))
    return;

  if(change_cfg_double(name,type,comment,"attachment_rate",&attachment_rate,dvalue,count,change,setting,ii))
    return;

  if(change_cfg_double(name,type,comment,"attachment_E0",&attachment_E0,dvalue,count,change,setting,ii))
    return;

  /* Constant external electric field. */
  if(change_cfg_double(name,type,comment,"E0_x",&E0_x,dvalue,count,change,setting,ii))
    return;

  if(change_cfg_double(name,type,comment,"E0_y",&E0_y,dvalue,count,change,setting,ii))
    return;

  if(change_cfg_double(name,type,comment,"E0_z",&E0_z,dvalue,count,change,setting,ii))
    return;

  /* Constant external electric field. */
  if(change_cfg_double(name,type,comment,"rise_time",&rise_time,dvalue,count,change,setting,ii))
    return;

  if(change_cfg_double(name,type,comment,"off_time",&off_time,dvalue,count,change,setting,ii))
    return;

  if(change_cfg_integer(name,type,comment,"has_photoionization",&has_photoionization,ivalue,count,change,setting,ii))
    return;

  /* Width of the initial seed in x-, y- and z-direction */
  if(change_cfg_double(name,type,comment,"seed_sigma_x",&seed_sigma_x,dvalue,count,change,setting,ii))
    return;

  if(change_cfg_double(name,type,comment,"seed_sigma_y",&seed_sigma_y,dvalue,count,change,setting,ii))
    return;

  if(change_cfg_double(name,type,comment,"seed_sigma_z",&seed_sigma_z,dvalue,count,change,setting,ii))
    return;

  if(change_cfg_double(name,type,comment,"seed_N",&seed_N,dvalue,count,change,setting,ii))
    return;

  /* Initial background ionization */
  if(change_cfg_double(name,type,comment,"background_ionization",&background_ionization,dvalue,count,change,setting,ii))
    return;

  /* Length of exponential increase of the pre-ionization for atmospherical
   * models
   */
  if(change_cfg_double(name,type,comment,"background_increase_length",&background_increase_length,
                    dvalue,count,change,setting,ii))
    return;

  /* Use the point-plane geometry? */
  if(change_cfg_integer(name,type,comment,"pois_inhom",&pois_inhom,ivalue,count,change,setting,ii))
    return;

  /* Number of mirror charges to use*/
  if(change_cfg_integer(name,type,comment,"pois_inhom_reflections",&pois_inhom_reflections,ivalue,count,change,setting,ii))
    return;

  /* Length of the needle */
  if(change_cfg_double(name,type,comment,"needle_length",&needle_length,dvalue,count,change,setting,ii))
    return;

  /* Radius of the needle */
  if(change_cfg_double(name,type,comment,"needle_radius",&needle_radius,dvalue,count,change,setting,ii))
    return;

  /* If nonzero, the charge is fixed, not floating
       *  Simulation of charged clouds close to the earth surface.
       */
  if(change_cfg_double(name,type,comment,"pois_inhom_fixed_q",&pois_inhom_fixed_q,dvalue,count,change,setting,ii))
    return;

  /* Constant ionization rate. */
  if(change_cfg_double(name,type,comment,"constant_source",&constant_source,dvalue,count,change,setting,ii))
    return;

  /* Random perturbations for stability analysis. */
  if(change_cfg_double(name,type,comment,"perturb_epsilon",&perturb_epsilon,dvalue,count,change,setting,ii))
    return;

  /* Perturb only modes up to perturb_max_k, i.e. large number to perturb all */
  if(change_cfg_integer(name,type,comment,"perturb_max_k",&perturb_max_k,ivalue,count,change,setting,ii))
    return;

  /****************** 
   * Sprites module *
   ******************/
  if(change_cfg_integer(name,type,comment,"sprite_module",&sprite_module,ivalue,count,change,setting,ii))
    return;

  if(change_cfg_double(name,type,comment,"dens_decay_len",&dens_decay_len,dvalue,count,change,setting,ii))
    return;

  if(change_cfg_double(name,type,comment,"sprite_dens_0",&sprite_dens_0,dvalue,count,change,setting,ii))
    return;

  if(change_cfg_double(name,type,comment,"sprite_dens_q",&sprite_dens_q,dvalue,count,change,setting,ii))
    return;

  if(change_cfg_integer(name,type,comment,"sprite_sign",&sprite_sign,ivalue,count,change,setting,ii))
    return;

}

/* @brief Read just one parameter from configuration file connected
 *        to setting2.
 *
 *  In case setting1 and setting2 are different then
 *  the value in setting1 will be overwritten by the value of setting2.
 *
 *  In case setting1 and setting2 are the same then
 *  the value in setting1 will be added.
 *
 *  @param[in] config_setting_t	*setting1	represents s configuration setting
 *  @param[in] config_setting_t	*setting2	represents s configuration setting
 *  @param[in] ii				position in array 'parameter_names' 
 *  @param[in] count				number of elements in parameter array to examine
 *  @param[in] *change				if TRUE then value in setting_default will be changed
 *  
 */
bool
read_parameter(config_setting_t *setting1,
	       config_setting_t *setting2,
               int ii,
	       int count,
	       bool *change) 
{
  const char *type,*name,*comment;
  double      dvalue;
  int         ivalue;
  int         bbool;
  const char* astring;
  config_setting_t *elem = config_setting_get_elem(setting2,ii);

  /* Only output the record if all of the expected fields are present. */
  if(config_setting_lookup_string(elem, "type",   &type) &&
     config_setting_lookup_string(elem, "name",   &name) &&
     config_setting_lookup_string(elem, "comment",&comment) &&
     (config_setting_lookup_float (elem,"value",  &dvalue) ||
      config_setting_lookup_int   (elem,"value",  &ivalue) ||
      config_setting_lookup_bool  (elem,"value",  &bbool ) ||
      config_setting_lookup_string(elem,"value",  &astring))
  ) {
    printf("\n# %s\n",comment);
    if(!(strncmp(type,"string",COMPARE_LIMIT)))
      printf("(%s) %s=%s\n",type,name,astring);
    if(!(strncmp(type,"double",COMPARE_LIMIT)))
      printf("(%s) %s=%g\n",type,name,dvalue);
    if(!(strncmp(type,"int",COMPARE_LIMIT)))
      printf("(%s) %s=%i\n",type,name,ivalue);
    if(!(strncmp(type,"bool",COMPARE_LIMIT)))
      printf("(%s) %s=%i\n",type,name,bbool);
    } else {
      printf("read_parameter: type is no string, double, int or bool but something else\n");
      return(FALSE);
    }

  change_cfg_parameters(name,type,comment,ii,ivalue,dvalue,bbool,astring,
                        setting1,count,change);
  return(TRUE);
}

/* @brief Reads the specifications of specie number 'ii' in the configuration
 *        file related to *setting. Initializes the fields of *temp_s.
 *
 *  @param[in] *setting		Of type config_setting_t. It represents a configuration
 *                              setting
 *  @param[in] ii		position in global 'parameter_names' array
 *  @param[out] *temp_s		Of type 'species_t'. Its fields will be initialized
 *                              with values "name","mass" and "charge" read
 *                              from setting connected to configuration file
 */
bool
read_specie(config_setting_t *setting,
	    int ii,
	    species_t *temp_s) 
{
  config_setting_t *elem = config_setting_get_elem(setting,ii);

  /* Only output the record if all of the expected fields are present. */
  if (!(config_setting_lookup_string(elem,"name",&temp_s->name) &&
        config_setting_lookup_float (elem,"mass",&temp_s->mass) &&
        config_setting_lookup_float (elem,"charge",&temp_s->charge)))
  {
	  printf("wrong types in kinetic file for species\n");
	  return(FALSE);
  }
  return(TRUE);
}

/* @brief Read the specifications of seed number 'ii' in the configuration
 *        file related to *setting. Initialize  *temp_se.
 *
 *  @param[in] *setting		of type config_setting_t. It represents a configuration
 *                              setting
 *  @param[in] ii		position in global array '*seed_index' 
 *  @param[out] *temp_se	of type 'seed_t' describing a seed
 */
bool
read_seed(config_setting_t *setting,
	  int ii,
	  seed_t *temp_se) 
{
  config_setting_t *elem = config_setting_get_elem(setting,ii);

  /* Only output the record if all of the expected fields are present. */
  if (!(config_setting_lookup_string(elem,"species",&temp_se->kind_species) &&
        config_setting_lookup_float (elem,"value",&temp_se->value) &&
        config_setting_lookup_string(elem,"type",&temp_se->kind_type)))
  {
    printf("wrong types in kinetic file for seed: species,value or type\n");
    return(FALSE);
  } 

  /* Find the position in the species-array of a given species */
  temp_se->species = find_species_by_name(temp_se->kind_species);

  /* Translate string 'kind_type' into integer value,
   * temp_se->type = -1 denotes unknown */ 
  if (strcmp(temp_se->kind_type,"gaussian") == 0)
    temp_se->type = 0;
  else if (strcmp(temp_se->kind_type,"constant") == 0)
    temp_se->type = 1;
  else
    temp_se->type = -1;

  if (!(config_setting_lookup_float(elem,"x0",&temp_se->x0)))
  { temp_se->x0=0.0; }
  if (!(config_setting_lookup_float(elem,"y0",&temp_se->y0)))
  { temp_se->y0=0.0; }
  if (!(config_setting_lookup_float(elem,"z0",&temp_se->z0)))
  { temp_se->z0=0.0; }

  if (!(config_setting_lookup_float (elem,"sigma_x",&temp_se->sigma_x)))
  { temp_se->sigma_x=0.0; }
  if (!(config_setting_lookup_float(elem,"sigma_y",&temp_se->sigma_y)))
  { temp_se->sigma_y=0.0; }
  if (!(config_setting_lookup_float(elem,"sigma_z",&temp_se->sigma_z)))
  { temp_se->sigma_z=0.0; }

  return(TRUE);
}

/* @brief Read the specifications of reaction number 'ii' in the configuration
 *        file related to *setting. Initialize *temp_r. 
 *
 *  @param[in] *setting		of type config_setting_t. It represents a configuration
 *                              setting
 *  @param[in] ii		position in global array '*reaction_index' 
 *  @param[out] *temp_r		of type 'reaction_t' describing a reaction
 */
bool
read_reaction(config_setting_t *setting,
	      int ii,
	      reaction_t *temp_r) 
{
  config_setting_t *elem = config_setting_get_elem(setting,ii);
  const char *table;
  const char *error;
  int cnt;

  /* Only output the record if fields nin and nout are present. */
  if (!(config_setting_lookup_string(elem,"reacttable",&table)))
  {
    printf("wrong types in kinetic file for reaction.reacttable\n");
    return(FALSE);
  } else
    temp_r->tablefile = (char *) table;

  /* Only output the record if fields nin and nout are present. */
  if (!(config_setting_lookup_int(elem,"nin",&temp_r->nin) &&
        config_setting_lookup_int(elem,"nout",&temp_r->nout)))
  {
    printf("wrong types in kinetic file for reactions\n");
    return(FALSE);
  }

  error="not_init";
  for(cnt = 0; cnt < REACTION_MAX_IN; ++cnt) {
    temp_r->input[cnt] = -1;
    temp_r->inname[cnt]=error;
  }

  config_setting_lookup_string(elem,"specin0",&temp_r->inname[0]);
  config_setting_lookup_string(elem,"specin1",&temp_r->inname[1]);
  config_setting_lookup_string(elem,"specin2",&temp_r->inname[2]);
  config_setting_lookup_string(elem,"specin3",&temp_r->inname[3]);

  /* Find the position in the species-array of a given species */
  for(cnt = 0; cnt < temp_r->nin; cnt++) {
    if(strcmp(temp_r->inname[cnt],error)==0)
    {
      printf("NOT all reaction.specin initialized for reaction %d and specie %d\n",ii,cnt);
      exit(1);
    }
    else
      temp_r->input[cnt] = find_species_by_name(temp_r->inname[cnt]);
  }

  for(cnt = 0; cnt < REACTION_MAX_OUT; ++cnt) {
    temp_r->output[cnt] = -1;
    temp_r->outname[cnt]=error;
  }

  config_setting_lookup_string(elem,"specout0",&temp_r->outname[0]);
  config_setting_lookup_string(elem,"specout1",&temp_r->outname[1]);
  config_setting_lookup_string(elem,"specout2",&temp_r->outname[2]);
  config_setting_lookup_string(elem,"specout3",&temp_r->outname[3]);
  config_setting_lookup_string(elem,"specout4",&temp_r->outname[4]);
  config_setting_lookup_string(elem,"specout5",&temp_r->outname[5]);

  /* Find the position in the species-array of a given species */
  for(cnt = 0; cnt < temp_r->nout; cnt++) {
    if(strcmp(temp_r->outname[cnt],error)==0)
    {
      debug(1,"temp_r->outname[%d]=%s; error=%s\n",cnt,temp_r->outname[cnt],error);
      printf("NOT all reaction.specout initialized for reaction %d and specie %d\n",ii,cnt);
      exit(1);
    }
    else
      temp_r->output[cnt] = find_species_by_name(temp_r->outname[cnt]);
  }

  return(TRUE);
}
