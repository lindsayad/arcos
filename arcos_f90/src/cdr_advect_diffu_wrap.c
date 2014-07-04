/** @file cdr_advect_diffu_wrap.c
*   @brief Wrapper function for the cdr_advect_diffu function. 
*
*   This function makes it possible to call the Fortran routine cdr_advect_vec
*   a c function. The c routine is named arcos_cdr_advect_vec.
*
*   @author Margreet Nool - 2009
*/

#include "stdlib.h"
#include "stdio.h"
#include "math.h"

int number_cdr=0;

/* cdr_advect: vectorized fortran90 routine; see f90/src/cdr_advect_vec.f90 */
void cdr_advect_vec_(double *mass,double *charge,double *dens_array, 
                 double *d_dens_array,double *er_array,double *ez_array,
		 double *diffusion_coeff,double *dr,double *dz,
		 int *sprite_module,int *r0,int *r1,int *z0,int *z1);

/* Calls hstcrt and checks for errors. If an error occurs, prints the
 corresponding message and exists. */
void 
arcos_cdr_advect_vec(double mass,double charge,double *dens_array, 
                 double *d_dens_array,double *er_array,int er_len,
		 double *ez_array,int ez_len,double diffusion_coeff,
                 double dr,double dz,int sprite_module,
                 int r0,int r1,int z0,int z1)
{
  if ( number_cdr % 1000 == 0 ) {

     fprintf(stdout, "wrap:# calls cdr_advec_vec       =%d\n",number_cdr);
     
  }

  cdr_advect_vec_(&mass,&charge,dens_array,d_dens_array,er_array,ez_array,
              &diffusion_coeff,&dr,&dz,&sprite_module,&r0,&r1,&z0,&z1);
  number_cdr += 1;

}

