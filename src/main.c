/**i @file main.c
 *  @brief Code regarding input of parameters + starting of the code. */

#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
//#include <getopt.h>
#include <math.h>

#include "cdr.h"
#include "cstream.h"
#include "grid.h"
#include "parameters.h"
#include "proto.h"
#include "species.h"

int main (int argc, char *argv[]);
void start_process (cdr_grid_t *cdr);
void output_tree (cdr_grid_t *cdr, int step, double t);
void check_params (void);

/** @brief The main program */
int
main (int argc, char *argv[])
{
  // param_t **p;
  // char *s;
  // struct option *opt;
  // struct option *long_options;
  // int c;
  cdr_grid_t *cdr;
  int count,count_user;
  static const char *output_file;
  config_t cfg,cfg_user;
  config_setting_t *setting,*setting_user;
  const char *prog_name,*f_defname,*f_username;

  f_defname   = "input/default.cfg";
  f_username  = "input/user_init.cfg";
  output_file = "input/user_continue.cfg";
  
  /* Open configuration cfg */
  config_init(&cfg);

  /* Read the default file.
   * The default file might be a result of a previous run, containing values
   * for all global variables.
   * If there is an error, report it and exit. */
  if(! config_read_file(&cfg, f_defname))
  {
    fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
            config_error_line(&cfg), config_error_text(&cfg));
    config_destroy(&cfg);
    return(EXIT_FAILURE);
  }

  /* Get the program name. */
  if(config_lookup_string(&cfg, "name", &prog_name))
    printf("Example : %s\n\n", prog_name);
  else
    fprintf(stderr, "No ' program name' setting in configuration file.\n");

  /* Output a list of all parameters in the file. */
  setting = config_lookup(&cfg, "param");
  if(setting != NULL)
  {
    int i;
    count = config_setting_length(setting);

    fprintf(stdout, "# params: count=%i\n",count);

    for(i = 0; i < count; ++i)
    {
      if (read_parameter(setting,setting,i,count,FALSE))
        continue;
    }

    /* Write out the updated configuration. */
    if(! config_write_file(&cfg, output_file))
    {
      fprintf(stderr, "Error while writing file.\n");
      config_destroy(&cfg);
      return(EXIT_FAILURE);
    }
    putchar('\n');
  }
  fprintf(stdout, "Default configuration successfully written to: %s\n",
          output_file);

/* ========================================================================== */
 
  /* Open configuration cfg_user */
  config_init(&cfg_user);

  /* Read the file made by the user, with adapted values.
   * This file does not have to contain entries for variables which remain
   * unchanged.
   * If there is an error, report it and exit. */
  if(! config_read_file(&cfg_user, f_username))
  {
    fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg_user),
            config_error_line(&cfg_user), config_error_text(&cfg_user));
    config_destroy(&cfg_user);
    return(EXIT_FAILURE);
  }

  /* Output a list of all parameters in the file. */
  setting_user = config_lookup(&cfg_user, "param");
  if(setting_user) 
  {
    int i;
    count_user = config_setting_length(setting_user);

    fprintf(stdout, "# params: count_user=%i\n",count_user);

    for(i = 0; i < count_user; ++i)
    {
      if (read_parameter(setting,setting_user,i,count,TRUE))
        continue;
    }

    /* Write out the updated configuration. */
    if(! config_write_file(&cfg, output_file))
    {
      fprintf(stderr, "Error while writing file.\n");
      config_destroy(&cfg);
      return(EXIT_FAILURE);
    } else {
      fprintf(stdout, "Updated configuration successfully written to: %s\n",
              output_file);
    }
    putchar('\n');
  }

  check_params ();

  cstream_init ();

  if (restart) {
     if (0 == strcmp (load_file, "")) {
       printf("restart file is empty\n");
       return 1;
     } else
       cdr = cdr_load_tree_r (output_dir, load_file, NULL);
  } else 
    cdr = cdr_scratch_init ();

  start_process(cdr);

  /* Destroy both configurations */
  config_destroy(&cfg_user);
  config_destroy(&cfg);

  return 0;
}

/** @brief Starts a given run on a given @a cdr tree. */
void 
start_process (cdr_grid_t *cdr)
{ 
  int i, i0;
  double t, dt;
  double t0, t1;
  struct timeval tv;
  FILE *fp;

  i0 = (int) (start_t / output_dt);

  debug(1,"before loop :start_t= %e; t = %e; dt = %e\n",start_t,t, dt);
  /* This is the main loop of the code. */
  for (i = i0, t = start_t; t < end_t; t += dt)
    {
      debug(1,"inside loop: i, t, dt : %d %e %e\n",i,t,dt);
    /* Set the field and the fixed charge when we have time-dependent 
       external fields (controled by rise_time) */
      cstream_set_field_at_time (t);
      dt = cdr_rk2 (cdr, attempt_dt, t);
      cdr_update_refined (&cdr);
      printf("t = %lf\tdt = %lf\tnext output at: %lf\n",t,dt,i * output_dt);
      if (output_dt <= 0.0 || t >= output_dt * i) 
	{
	  output_tree (cdr, i, t);
	  if (cdr->ntheta > 1) 
	    {
	      dft_out_weights (cdr, output_dir, t);
	    }
	  i++;
	}
    }
  cdr_dump (cdr, output_dir, "C");
  cdr_free_r (cdr);
  cstream_end ();
}

/** @brief Outputs the tree contained in cdr.
 *
 * If pois_output (global parameter) is true, it also dumps the Poisson tree.
 */
void
output_tree (cdr_grid_t *cdr, int step, double t)
{
  char *gname;
  pois_grid_t **pois_modes;
  FILE *fp;
  char logfile[100];

  asprintf (&gname, "C%.3d", step);
  cdr_set_ext_bnd_r (cdr);

  /* !!!!!!!!!!!! */
  if (1 || t > start_t || 0 == strcmp (load_file, "")) {
    cdr_dump_r (cdr, output_dir, gname, NULL, t);
    cdr_dump_frames (cdr, output_dir, gname);

    printf ("CDR grid %s saved\n", gname);
  }

  free (gname);

  if (pois_output) {
    int mode;
    pois_modes = cdr_calc_field_r (cdr, TRUE);

    for (mode = 0; mode < max_ntheta; mode++) {
      asprintf (&gname, "P%.3d%.2d", step, mode);
      pois_dump_r (pois_modes[mode], output_dir, gname);
      pois_free_r (pois_modes[mode]);
      printf ("Poisson grid %s saved\n", gname);
      free (gname);
    }

    free (pois_modes);
  }
}

/** @brief Checks for some errors in the parameters that will produce
 *  nonsense output.
 */
void
check_params ()
{
  if ((gridpoints_z % 4) != 0 && pois_inhom == 1) {
    fatal ("To use inhomogeneous fields, gridpoints_z must be a"
	   " multiple of 4\n");
  }
}
