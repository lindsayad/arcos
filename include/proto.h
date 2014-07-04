/** @file proto.h
 *  @brief Function prototypes for exported functions.
 */

#ifndef _RZ_ARRAY_H_
# include "rz_array.h"
#endif

/* rz_array.c */
rz_array_t *rz_new_3d_a (int r0, int z0, int r1, int z1, int ntheta);
rz_array_t *rz_guest (rz_array_t *host, int r0, int z0, int r1, int z1);
void rz_set_zero (rz_array_t *array);
void rz_set_periodic (rz_array_t *array);
void rz_copy_bnd (rz_array_t *from, rz_array_t *to, 
		  int sign, REAL *start_from, REAL *start_to, 
		  int dim0, int inout_from, int inout_to,
		  int dim1, int dim1_0, int dim1_1, 
		  int dim2, int dim2_0, int dim2_1);
void rz_set_bnd (rz_array_t *array, int sign, REAL *start, int dim0, int inout,
		 int dim1, int dim1_from, int dim1_to, 
		 int dim2, int dim2_from, int dim2_to);
void rz_free (rz_array_t *array);
void rz_copy (rz_array_t *fro, int rfro, int zfro, 
	      rz_array_t *to, int rto, int zto,
	      int rn, int zn);
void rz_copy_modes (rz_array_t *fro, int rfro, int zfro, 
		   rz_array_t *to, int rto, int zto,
		    int rn, int zn, int nmode_fro, int nmode_to);
void rz_dump (rz_array_t *rz_array, const char *fname, const char *mode,
	      int r0, int z0, int r1, int z1);
void rz_dump_3d (rz_array_t *rz_array, const char *fname, const char *mode,
		 int r0, int z0, int r1, int z1, int ntheta);
void rz_axis_dump (const char *fname, int x0, int x1, double delta);


#ifndef _GRID_H_
# include "grid.h"
#endif

/* grid.c */
int grid_max_depth_r (grid_t *grid);
double grid_rmin_r (grid_t *grid);
void grid_print_r (grid_t *grid, int indent);
int grid_contains (grid_t *grid, int i, int j, int check);
int grid_overlap (grid_t *grid1, grid_t *grid2, int buff1, int buff2,
		  int *left, int *bottom, int *right, int *top, 
		  int *level_diff);
int grid_overlap_with_shifts (grid_t *grid1, grid_t *grid2, 
			      int buff1, int buff2, 
			      int *left, int *bottom, int *right, int *top, 
			      int *level_diff, int shift_r, int shift_z);
void grid_inherit_ext_bound (grid_t *grid);
grid_t *grid_finest_containing_r (grid_t *grid, double r, double z);
int grid_howmany_children (grid_t *grid);
grid_t *grid_get_child (grid_t *grid, int n);


#ifndef _INTERPOL2_H_
# include "interpol2.h"
#endif

interpol_t *interpol_new_a (double Lr, double Lz, interpol_method_t *method);
void interpol_free (interpol_t *this);
void interpol_set_stencil (interpol_t *this, double r0, double z0, ...);
void interpol_set_stencil_at (grid_t *grid,
			      interpol_t *this, double r0, double z0, 
			      rz_array_t *ar,
			      int ir, int iz, int itheta);
void interpol_set_coeffs (interpol_t *this);
double interpol_apply (interpol_t *this, double r, double z);

#ifndef _MAPPER_H_
# include "mapper.h"
#endif

void map_grid (mapper_t **mappers, grid_t *source, grid_t *target, int ntheta,
	       int copy, int interpol, int coarsen, int s_buff, int t_buff);
void map_grid_r (mapper_t **mappers, grid_t *source, grid_t *target, 
		 int ntheta, int copy, int interpol, int coarsen, 
		 int s_buff, int t_buff);
void map_trees_r (mapper_t **mappers, grid_t *source, grid_t *target, 
		  int ntheta, int copy, int interpol, int coarsen,
		  int s_buff, int t_buff);

#ifndef _CDR_H_
# include "cdr.h"
#endif

#ifndef _POISSON_H_
# include "poisson.h"
#endif

/* cdr.c */
void cdr_init (void);
void cdr_end (void);

cdr_grid_t *cdr_new_3d_a (int r0, int z0, int r1, int z1, int ntheta);
cdr_grid_t *cdr_like_a (cdr_grid_t *grid);
cdr_grid_t *cdr_clone_a (cdr_grid_t *grid);
void cdr_set_periodic (cdr_grid_t *grid);
void cdr_free (cdr_grid_t *grid);
void cdr_free_r (cdr_grid_t *grid);
void cdr_calc_charge (cdr_grid_t *grid);
void cdr_calc_charge_r (cdr_grid_t *root);
void cdr_dft_charge_r (cdr_grid_t *grid, int sign);
cdr_grid_t *cdr_create_coarser_a (cdr_grid_t *grid);
cdr_grid_t *cdr_add_coarser_grids_a (cdr_grid_t *prev_root, int n);
void cdr_free_coarser_grids (cdr_grid_t *prev_root, int n);
void cdr_add_ext_field (cdr_grid_t *grid);
void cdr_add_ext_field_r (cdr_grid_t *grid);
void cdr_add_inhom_field_r (cdr_grid_t *cdr, double q);
void cdr_calc_eabs (cdr_grid_t *grid);
void cdr_calc_eabs_r (cdr_grid_t *grid);
void cdr_nonegative (cdr_grid_t *grid);
void cdr_nonegative_r (cdr_grid_t *grid);
pois_grid_t **cdr_calc_field_r (cdr_grid_t *grid, int return_pois);
void cdr_set_ext_bnd (cdr_grid_t *grid);
void cdr_set_ext_bnd_r (cdr_grid_t *grid);
void cdr_calc_d_dens (cdr_grid_t *grid);
void cdr_calc_d_dens_r (cdr_grid_t *grid);
void cdr_advect_diffu (cdr_grid_t *grid);
void cdr_advect_diffu_r (cdr_grid_t *grid);

double cdr_courant (cdr_grid_t *grid);
void cdr_update (cdr_grid_t *orig, cdr_grid_t *dest, double h);
void cdr_rk2_update (cdr_grid_t *dens_0, cdr_grid_t *d_dens_1, 
		     cdr_grid_t *d_dens_2, cdr_grid_t *dest, 
		     double h);
void cdr_rk2_update_r (cdr_grid_t *dens_0, cdr_grid_t *d_dens_1, 
		       cdr_grid_t *d_dens_2, cdr_grid_t *dest, 
		       double h);
void cdr_self_update_r (cdr_grid_t *grid, double h);
void cdr_like_update_ar (cdr_grid_t *grid, cdr_grid_t *new_grid, double h);
double cdr_rk2 (cdr_grid_t *grid, double h, double t);

void cdr_update_refined (cdr_grid_t **ptree);
void cdr_calc_maxs (cdr_grid_t *grid);
void cdr_refine (cdr_grid_t *grid);
void cdr_refine_r (cdr_grid_t *grid, cdr_grid_t *source);
void cdr_match_r (cdr_grid_t *grid1, cdr_grid_t *grid2);
void cdr_set_bnd (cdr_grid_t *grid);
void cdr_set_bnd_r (cdr_grid_t *grid);

void cdr_restrict (cdr_grid_t *grid);
void cdr_restrict_r (cdr_grid_t *grid);

void cdr_set_dens (cdr_grid_t *cdr, int species, int mode, double factor,
		   double (*f) (double, double, double));
void cdr_init_dens (cdr_grid_t *cdr);
cdr_grid_t *cdr_scratch_init (void);

mapper_t **cdr_mappers_a (interpol_method_t *interp_method);
void cdr_free_mappers (mapper_t **mappers);

void cdr_dump (cdr_grid_t *grid, const char *prefix, const char *name);
void cdr_dump_r (cdr_grid_t *grid, const char *prefix, const char *name,
		 FILE *infp, double sim_time);
cdr_grid_t *cdr_load_tree_r (const char *prefix, const char *name, FILE *infp);
void cdr_dump_frames (cdr_grid_t *grid, const char *prefix, const char *name);

/* poisson.c */
void pois_init (void);
pois_grid_t *pois_new_a (int r0, int z0, int r1, int z1);
pois_grid_t *pois_new_3d_a (int r0, int z0, int r1, int z1, int ntheta);
void pois_free (pois_grid_t *grid);
void pois_free_r (pois_grid_t *grid);
pois_grid_t *pois_new_glob_a (int r0, int z0, int r1, int z1, int level);
pois_grid_t *pois_init_tree_a (int r0, int z0, int r1, int z1);
REAL *pois_boundary_a (pois_grid_t *grid, int boundary);
void pois_solve_grid (pois_grid_t *grid, pois_problem_t *prob, 
		      double lambda, double s);
void pois_set_phi_boundaries (pois_grid_t *grid, REAL *boundaries[], 
			      int left_uns, int right_neu, int bottom_neu,
			      int top_neu);
int pois_map_charge (pois_grid_t *pois, cdr_grid_t *cdr, int mode);
void pois_set_error (pois_grid_t *grid);
int pois_refine (pois_grid_t *grid, double threshold);
pois_grid_t **pois_solve_a (cdr_grid_t *cdr, pois_problem_t *prob);
pois_grid_t **pois_gen_solve_a (cdr_grid_t *cdr, pois_problem_t *prob, 
				mapper_t **mappers, double es);
pois_grid_t *pois_solve_mode (cdr_grid_t *tree, cdr_grid_t *cdr,
			      pois_problem_t *prob, int mode, double es);
void pois_solve_r (pois_grid_t *pois, cdr_grid_t *cdr, pois_problem_t *prob, 
		   int mode, double es, double threshold);
double pois_phi_at (pois_grid_t *grid, double r, double z, double theta);
void pois_dump (pois_grid_t *grid, const char *prefix, const char *name);
void pois_dump_r (pois_grid_t *grid, const char *prefix, const char *name);
void pois_write_error_r (pois_grid_t *grid, FILE *fp);
void pois_error_measures (pois_grid_t *grid, double *L1, double *L2, 
			  double *Lmax);
void pois_inhom_init (void);
double pois_inhom_phi (double r, double z);
double pois_inhom_er (double r, double z);
double pois_inhom_ez (double r, double z);
double pois_inhom_q_factor (pois_grid_t *pois);
void pois_add_inhom_phi_r (pois_grid_t *grid, double q);

#ifndef _PHOTO_H_
# include "photo.h"
#endif

/* photo.c */
void photo_init ();
void photo_register (double A, double lambda);
void photo_copy_list (photo_term_t *src, photo_term_t **dest);
void photo_unregister_all (void);
void photo_dft_r (cdr_grid_t *grid, int sign);
void photo_copy_source (cdr_grid_t *grid);
void photo_add_term (photo_term_t *term, cdr_grid_t *cdr);
void photo_add_term_r (photo_term_t *term, cdr_grid_t *cdr);
pois_grid_t **photo_calc_term (photo_term_t *term, cdr_grid_t *cdr, int i);
void photo_calc (photo_term_t *terms, cdr_grid_t *cdr);
void photo_load_file (char *fname);

/* cstream.c */
void init_parameters(void);
void cstream_init (void);
void cstream_end (void);
void cstream_set_field_at_time (double t);

/* dft.c */
void dft_transform (rz_array_t *in, rz_array_t *out, int sign);
void dft_diff (grid_t *grid, rz_array_t *f);
void dft_weight (cdr_grid_t *cdr, rz_array_t *var, double weights[], 
		 double power);
void dft_out_weights (cdr_grid_t *grid, const char *prefix, double t);
void dft_perturb (cdr_grid_t *cdr, rz_array_t *var, double *epsilon_k);
void dft_dens_perturb_r (cdr_grid_t *grid, int species, double *epsilon_k);


#ifndef _REACTION_H_
# include "reaction.h"
#endif

/* reaction.c */
int find_species_by_name(const char *spec_name);
void react_add (reaction_t *react);
void react_init ();
void react_apply_r (reaction_t *react, cdr_grid_t *grid, int overwrite);
void react_apply_all (cdr_grid_t *grid);

#ifndef _REACT_TABLE_H_
# include "react_table.h"
#endif

/* kinetic module (usually minimal.c) */
void kinetic_init (void);

/* sprites module */
double spr_density_at (double altitude);
void spr_init ();
void spr_hook (cdr_grid_t *grid);
void spr_update (double altitude);
double spr_head_altitude (cdr_grid_t *grid, int sign);

/* rt.c */
void kinetic_init (void);
void read_input_file(const char *f_kinetic_name, const char *filename);





