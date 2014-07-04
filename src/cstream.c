/** @file cstream.c
 *  @brief General initialization/termination functions.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#define ALLOC_PARAMS

#include "cstream.h"
#include "parameters.h"
#include "proto.h"
#include "species.h"

double E_x, E_y, E_z;
/**< For time-dependent electric fields E_x, E_y, E_z there is a difference between
   E0 (the highest field) and E (the field at a given time). */

extern double pois_inhom_fixed_q_t;
/**< The same is true for time-dependent fixed charges (useful to model
   the charging of a cloud in modeling sprites. */

double *dr, *dz, dtheta;
/**< dr[i] and dz[i] represent the grid size at level @a i.  Note that
   @a i can go from -extra_pois_levels...max_levels. */
double *dr_start, *dz_start;

double *w2k, *wk;
/**< w2k[k] = 2 / dtheta^2 (1 - cos(k dtheta)). This numbers appear due to the
   finite-differences derivatives of a Fourier transform.  If we would make
   a continuous fourier transform, they will approach k^2.
*/


char *invok_name = "cstream";
/**< Name of this program. */

const double twopi = 6.283185307179586477L;
/**< \f$ 2 \times \pi \f$. */

const double invfourpi = 0.079577471545947667884441882L;
/**< \f$ 1 / (4 \pi )\f$ */

const double invpi32 = 0.179587122125166561689081L;
/**< \f$ 1 / (\pi ^ {3/2} )\f$ */

/** @brief Initializes the parameters.
 *
 *  The global parameters will get their default values.
 *  These values correspond with the values in input/defaults.cfg
 */
void
init_parameters(void)
{
    debug (1, "entry init_parameters\n");
    /* An identification name for this run */
    prog_id = "example";

    /* Output directory.  */
    output_dir = "output/";

    /* Kinetics input file */
    kin_input = "input/kinetic_example.cfg";

    /* If restart is TRUE, the simulation will continue with data from a previous run */
    restart = 0;

    /* If restart is TRUE, the name of the file with data from previous run, otherwise empty */
    load_file = "";

    /* Time interval for output to be written to disk */
    output_dt = 0.100;

    /* Output of the Poisson grids, including the potential? */
    pois_output = 0;

    /* Margin outside the grids in the output of the cdr equation */
    cdr_output_margin =  0;

    /* Margin outside the grids in the output of the poisson equation */
    pois_output_margin = 1;

    /* If the time steps are smaller than this number, the program issues a warning */
    warn_min_timestep = 1e-06;

    /* Maximum disk space, in Mb, to use */
    max_disk_space_mb = 1048576;

    /* Number of R and Z gridpoints at level 0 */
    gridpoints_r = 600;
    gridpoints_z = 600;

    /* Number of azimuthal gridcells and modes */
    max_ntheta = 1;

    /* Initial and end time */
    start_t = 0.0;
    end_t = 0.12;

    /* Attempted timestep.  The actual timestep may be larger */
    attempt_dt = 50.0;

    /* Extra levels for the Poisson solver */
    extra_pois_levels = 2;

    /* Maximum level of refinement. Use a big number here */
    max_levels = 64;

    /* Error threshold that leads to refinement in the Poisson code. */
    pois_max_error = 0.001;

    /* Maximum level of refinement in the Poisson equation. */
    pois_max_level = 3;

    /* Extra levels for the photo-ionization solver */
    extra_photo_levels = -1;

    /* Maximum level of refinement in the photo-ionization solver. */
    photo_max_level = 4;

    /* Error threshold that leads to refinement in the photo-ionization code. */
    photo_max_error = 0.01;

    /* Photo-ionization boundary condition at r = L_r, z = 0, z = L_z.  1 for Hom. Neumann, -1 for Hom. Dirichlet */
    photo_bnd_right  = -1;
    photo_bnd_bottom = -1;
    photo_bnd_top    = -1;

    /* Extra levels for the photo-ionization solver */
    extra_photo_levels_2 = -1;

    /* Maximum level of refinement in the photo-ionization solver. */
    photo_max_level_2 = 4;

    /* Error threshold that leads to refinement in the photo-ionization code. */
    photo_max_error_2 = 0.01;

    /* Photo-ionization boundary condition at r = L_r, z = 0, z = L_z.  1 for Hom. Neumann, -1 for Hom. Dirichlet */
    photo_bnd_right_2  = -1;
    photo_bnd_bottom_2 = -1;
    photo_bnd_top_2    = -1;

    /* Particles boundary condition at z = 0, z = L_z, r = L_r.  1 for Hom. Neumann, -1 for Hom. Dirichlet */
    cdr_bnd_bottom = 1;
    cdr_bnd_top    = 1;
    cdr_bnd_right  = 1;

    /* Potential boundary condition at r = L_r, z = 0,  z = L_z.  1 for Hom. Neumann, -1 for Hom. Dirichlet */
    pois_bnd_right  = -1;
    pois_bnd_bottom = -1;
    pois_bnd_top    = -1;

    /* Maximum advection and diffusion Courant number */
    nu_a = 0.2;
    nu_d = 0.2;

    /* Maximum ratio of dt/relaxation time */
    nu_rt = 0.2;

    /* Maximum ratio of change of the densities (set to a very large number to ignore) */
    nu_f = 1e+20;

    /* Refinement threshold for the electric field */
    ref_threshold_eabs = 0.2;

    /* Maximum refinement level reached through ref_threshold_eabs */
    ref_level_eabs = 4;

    /* Refinement threshold for the curvature of the charge, densities */
    ref_threshold_charge = 0.004;
    ref_threshold_dens = 0.004;

    /* Refinement threshold for the densities in the leading edge */
    ref_threshold_edge = 10000.0;

    /* r-length and  z-length  of the minimal refinement area in the cdr equation */
    cdr_brick_dr = 8;
    cdr_brick_dz = 8;

    /* Maximum level of refinement in the Fluid equation. */
    cdr_max_level = 3;

    /* Interpolation method for the grid interior, and grid boundaries (0=zero_masses, 1=quadratic_masses [default], 2=wackers_masses, 3=quadlog */
    cdr_interp_in = 1;
    cdr_interp_bnd = 1;

    /* Length in r and z of the complete domain */
    L_r = 13044.0;
    L_z = 13044.0;

    /* Isotropic difussion coefficient */
    diffusion_coeff = 0.1;

    /* Whether the code includes photoionization or not */
    has_photoionization = 1;

    /* The name of a file from which we can read the photoionization parameters */
    photoionization_file = "input/air760torr.photo";

    /* Rate of dissociative attachment */
    attachment_rate = 0.0;

    /* E0 in the exp(-E0/E) factor in the attachment expression. */
    attachment_E0 = 0.0;

    /* x-, y- and z-component of the external electric field */
    E0_x = 0.0;
    E0_y = 0.0;
    E0_z = -0.06;

    /* Rise time of the electric field (0 for instantaneous rise) */
    rise_time = 0.0;

    /* Time to switch off the electric field (0.0 means never) */
    off_time = 0.0;

    /* x-, y- and z-width of the initial seed */
    seed_sigma_x = 0.0;
    seed_sigma_y = 0.0;
    seed_sigma_z = 0.0;

    /* Number of electrons in the initial seed */
    seed_N = 0.0;

    /* Initial at z=0 densities of electrons and ions */
    background_ionization = 0.0;

    /* Length of exponential increase of the pre-ionization (for atmospherical models) */
    background_increase_length = 0.0;

    /* Use the point-plane geometry? */
    pois_inhom = 1;

    /* Number of mirror charges to use */
    pois_inhom_reflections = 4;

    /* Length and radius of the needle */
    needle_length = 2500.0;
    needle_radius = 400.0;

    /* If nonzero, the charge is fixed, not floating (simulation of charged clouds close to the earth surface) */
    pois_inhom_fixed_q = 0.0;

    /* Constant ionization rate */
    constant_source = 0.0;

    /* Initial perturbation to the axisymmetric configuration */
    perturb_epsilon = 0.0;

    /* Perturb only modes up to perturb_max_k (large number to perturb all) */
    perturb_max_k = 1024;

    /* 1 if the sprite module is activated, 0 otherwise */
    sprite_module = 0;

    /* Lenght of exponential decay of the density w/r to altitude */
    dens_decay_len = 0.0;

    /* Density at z = 0 */
    sprite_dens_0 = 0.0;

    /* Quenching density */
    sprite_dens_q = 0.0;

    /* Sign of the sprite head that we are following (the other will not be reliable */
    sprite_sign = -1;
    debug (1, "exit init_parameters\n");
}

  int n, i;
/** @brief Initializes the grid sizes.  The parameters dr_root and dz_root
 *  specify the grid size at level 0 (dz_root = dz[0], dr_root[0] = dr[0])
 */
static void
init_gridsizes_a (void)
{
  int n, i;
  double root_dr, root_dz;

  root_dr = L_r / gridpoints_r;
  root_dz = L_z / gridpoints_z;

  n = max_levels + extra_pois_levels + 1;

  dr_start = (double *) xmalloc(sizeof(double) * n);
  dz_start = (double *) xmalloc(sizeof(double) * n);

  dr = dr_start + extra_pois_levels;
  dz = dz_start + extra_pois_levels;

  dr_start[0] = root_dr * (1 << extra_pois_levels);
  dz_start[0] = root_dz * (1 << extra_pois_levels);

  for (i = 1; i < n; i++) {
    dr_start[i] = dr_start[i - 1] / 2.0;
    dz_start[i] = dz_start[i - 1] / 2.0;

    debug (3, "dr[%d] = %e\n", i - extra_pois_levels,
	   dr[i - extra_pois_levels]);
    debug (3, "dz[%d] = %e\n", i - extra_pois_levels,
	   dz[i - extra_pois_levels]);
  }

  dtheta = twopi / max_ntheta;
}

/** @brief Initializes the vector of wk's for the Helmholtz equation. */
static void
init_wk_a (void)
{
  int k;
  double twobydtheta2;

  debug (2, "init_w2k_a()\n");

  w2k = (double*) xmalloc (sizeof(double) * max_ntheta);
  wk  = (double*) xmalloc (sizeof(double) * max_ntheta);

  twobydtheta2 = 2. / (dtheta * dtheta);

  for (k = 0; k < max_ntheta / 2 + (max_ntheta % 2); k++) {
    double w2k_ = twobydtheta2 * (1 - cos (k * dtheta));
    double re_wk, im_wk;

    re_wk = (cos (k * dtheta) - 1) / dtheta;
    im_wk = sin (k * dtheta) / dtheta;

    wk[k] = re_wk;
    w2k[k]= w2k_;

    if (k != 0) {
      wk[max_ntheta - k] = im_wk;
      w2k[max_ntheta - k] = w2k_;
    }
  }

  if ((max_ntheta % 2) == 0) {
    w2k[max_ntheta / 2] = twobydtheta2 * (1 - cos (k * dtheta));
    wk[max_ntheta / 2] = (cos (k * dtheta) - 1) / dtheta;
  }

  assert (wk[0] == 0.);
}

/** @brief Frees the space allocated for dr and dz. */
static void
free_gridsizes (void)
{
  free (dr_start);
  free (dz_start);
}

/** @brief Frees the allocated space for \f$ \lambda \f$. */
static void
free_wk (void)
{
  free (wk);
  free (w2k);
}

/** @brief  Here all the initialization calls. */
void
cstream_init (void)
{
  /* For the perturbations, it is better to start with different seeds at
     each run. */
  srand (time (0));

  init_gridsizes_a ();
  init_wk_a ();
  react_init ();
  cdr_init ();
  pois_init ();
  photo_init ();
  if (has_photoionization) {
    printf("Photoionization\n");
    photo_load_file (photoionization_file);
  } else {
    printf("No photoionization\n");
  }

  if (sprite_module) spr_init ();
}

/** @brief Here all the cleaning and finishing calls. */
void
cstream_end (void)
{
  photo_unregister_all ();
  free_gridsizes ();
  free_wk ();
  cdr_end ();
}


/** @brief When we use a time-varying external field, we update the field
 *  components at each timesteps. */
void
cstream_set_field_at_time (double t)
{
  static int rise_reached = FALSE, off_reached = FALSE;
  double factor;

  if (rise_reached && off_reached)
    return;

  if (t >= rise_time) {
    rise_reached = TRUE;
  }

  if (t >= off_time || off_time == 0.0) {
    off_reached = TRUE;
  }

  /* We make sure that we never use a field _larger_ than E0 and that
     if rise_time is 0 we do not calculate t / 0.0.  */
  factor = t >= rise_time? 1.0: (t / rise_time);

  factor = (off_time > 0.0 && t >= off_time)? 0.0: factor;

  E_x = factor * E0_x;
  E_y = factor * E0_y;
  E_z = factor * E0_z;

  if (pois_inhom_fixed_q != 0.0) {
    pois_inhom_fixed_q_t = factor * pois_inhom_fixed_q;
  }

  return;
}

/** @brief Functions for a constant external electric field given in (x, y, z) components. */
double
e0_r (double r, double z, double theta)
{
  return E_x * cos (theta) + E_y * sin (theta);
}

/** @brief Functions for a constant external electric field given in (x, y, z) components. */
double
e0_z (double r, double z, double theta)
{
  return E_z;
}

/** @brief Functions for a constant external electric field given in (x, y, z) components. */
double
e0_theta (double r, double z, double theta)
{
  return -E_x * cos (theta) + E_y * sin (theta);
}

/** Initializes the component of the external electric field in \f$r\f$-direction. */
decl_field_comp(r) = e0_r;
/** Initializes the component of the external electric field in \f$z\f$-direction. */
decl_field_comp(z) = e0_z;
/** Initializes the component of the external electric field in \f$\theta\f$-direction. */
decl_field_comp(theta) = e0_theta;
