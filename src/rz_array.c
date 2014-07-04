/** @file rz_array.c
 *  @brief General functions to work with FORTRAN-compatible 2d/3d arrays.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "parameters.h"
#include "proto.h"
#include "rz_array.h"
#include "species.h"

long int used_disk_space = 0;
/*!< This is the total disk space used so far for the array dumps.
   It is used to limit the disk space used by a simulation and thus
   avoid that one simulation gets crazy and stops other simulations
   that are sharing the same disk space. */

static void check_disk_space(void);

/** @brief Creates a new 3D array structure. */
rz_array_t *
rz_new_3d_a (int r0, int z0, int r1, int z1, int ntheta)
{
  int rmax, zmax;
  rz_array_t *array;

  debug (3, "rz_new_3d_a (%d, %d, %d, %d, %d)\n", r0, z0, r1, z1, ntheta);

  /* We allow a margin of two cells to set the boundary conditions (see below)*/
  r0 -= 2; z0 -= 2; r1 += 2; z1 += 2;

  rmax = r1 - r0;
  zmax = z1 - z0;

  assert (rmax > 0 && zmax > 0 && ntheta > 0);

  debug (3, "r0 = %d, z0 = %d, r1 = %d, z1 = %d\n", r0, z0, r1, z1);

  array = (rz_array_t *) xmalloc (sizeof(rz_array_t));

  if (ntheta != 1) {
    array->len = rmax * zmax * (ntheta + 4);
    array->theta0 = -2;
  } else {
    array->len = rmax * zmax;
    array->theta0 = 0;
  }

  /*  Initially all arrays are zeroed. */
  array->data = (REAL *) xcalloc (array->len, sizeof(REAL));


  /* These strides are to make the array compatible with array(ir, iz) 
     in FORTRAN. */
  array->strides[R_INDX] = 1;
  array->strides[Z_INDX] = rmax;
  array->strides[THETA_INDX] = rmax * zmax;

  array->dim = ntheta == 1? 2: 3;

  debug (3, "rmax = %d\n", rmax);
  debug (3, "strides = {%d, %d}\n", array->strides[0], array->strides[1]);

  array->r0 = r0;
  array->z0 = z0;
  /* Note that these are NOT the r0, z0 we received.  The idea here is that
     RZ(array, r0, z0) will produce the correct result, where the allocated
     array looks like (in the picture, r0, z0 ARE the received ones).
     Note: in the picture there is only one buffer cell, while we actually
     allocate two.

       +--------------+--------------+...+--------------+--------------+
       |r0 - 1, z0 - 1|              |   |              |              |
       +--------------+--------------+...+--------------+--------------+
       |              |   r0, z0     |   |              |              |
       +--------------+--------------+...+--------------+--------------+
                     ...             |...|             ...
       +--------------+--------------+...+--------------+--------------+
       |              |              |   |r1 - 1, z1 - 1|              |
       +--------------+--------------+...+--------------+--------------+
       |              |              |   |              |    r1, z1    |
       +--------------+--------------+...+--------------+--------------+

     but the rows z0 +/- 1 and the columns r0 +/- 1 do not belong to the
     physical space and are used only to specify boundary conditions.
  */

  array->nr = r1 - r0 - 4;
  array->nz = z1 - z0 - 4;
  array->ntheta = ntheta;

  array->host = NULL;

  return array;
}

/** @brief Creates guest array: This is an array that contains
 *  part of another one.
 *
 *  We use the memory of the original array and do not allocate
 *  new memory.
 */
rz_array_t *
rz_guest (rz_array_t *host, int r0, int z0, int r1, int z1)
{
  int ntheta, rmax, zmax;
  rz_array_t *array;

  ntheta = host->ntheta;
  array = (rz_array_t *) xmalloc (sizeof(rz_array_t));

  rmax = r1 - r0;
  zmax = z1 - z0;

  if (ntheta != 1) {
    array->len = rmax * zmax * (ntheta + 4);
  } else {
    array->len = rmax * zmax;
  }
  
  array->r0 = r0;
  array->z0 = z0;

  array->theta0 = host->theta0;

  array->strides[R_INDX] = host->strides[R_INDX];
  array->strides[Z_INDX] = host->strides[Z_INDX];
  array->strides[THETA_INDX] = host->strides[THETA_INDX];

  array->dim = host->dim;

  array->nr = r1 - r0 - 4;
  array->nz = z1 - z0 - 4;
  array->ntheta = host->ntheta;

  array->data = RZTP (host, r0, z0, array->theta0);
  array->host = host;

  return array;
}

/** @brief Sets the data of an array to zero */
void
rz_set_zero (rz_array_t *array)
{
  debug (3, "rz_set_zero\n");
  memset (array->data, 0, array->len * sizeof(REAL));
}

/** @brief Sets periodic boundary conditions for theta in an array.*/
void
rz_set_periodic (rz_array_t *array)
{
  /* If the array is 2D, we do nothing. */
  if (array->ntheta == 1) return;

  rz_copy_modes (array, array->r0, array->z0, 
		 array, array->r0, array->z0, 
		 array->nr + 4, array->nz + 4, array->ntheta - 1, -1);

  rz_copy_modes (array, array->r0, array->z0, 
		 array, array->r0, array->z0, 
		 array->nr + 4, array->nz + 4, array->ntheta - 2, -2);

  rz_copy_modes (array, array->r0, array->z0, 
		 array, array->r0, array->z0, 
		 array->nr + 4, array->nz + 4, 0, array->ntheta);

  rz_copy_modes (array, array->r0, array->z0, 
		 array, array->r0, array->z0, 
		 array->nr + 4, array->nz + 4, 1, array->ntheta + 1);

}


/** @brief Sets boundaries for one array reading data from other array
 * (but both can be the same: see rz_set_bnd).
 *
 * Starting at *start_xxx, it sweeps a 2-dimensional subspace of the grid
 * in dimensions dim1 and dim2 (and perpendicular to dim0). There
 * it sets the boudary conditions.
 *
 * The value of from is multiplied by sign.
 * Thus, if the boundary is itself, is -1 for homogeneous dirichlet
 * and 1 for homogeneous Neumann (though you should better use
 * BND_CND_HNEUMANN and BND_CND_HDIRICHLET)
 *
 * @a inout is -1 if we look/set the values on smaller values
 * of dim0 1 if we look at larger values.
 * (though it is better to use BND_INWARD=-1 and BND_OUTWARD=1).
 *
 * Note: it is better to call this function _before_ rz_set_periodic.
 */
void
rz_copy_bnd (rz_array_t *from, rz_array_t *to, 
	     int sign, REAL *start_from, REAL *start_to, 
	     int dim0, int inout_from, int inout_to,
	     int dim1, int dim1_0, int dim1_1, 
	     int dim2, int dim2_0, int dim2_1)
{
  int i, j;
  REAL *pfrom, *pto;

  for (i = dim1_0; i < dim1_1; i++) {
    for (j = dim2_0; j < dim2_1; j++) {
      pfrom = start_from + i * from->strides[dim1] + j * from->strides[dim2];
      pto = start_to + i * to->strides[dim1] + j * to->strides[dim2];

      *(pto + inout_to * to->strides[dim0]) = sign * (*pfrom);

      *(pto + 2 * inout_to * to->strides[dim0]) =  
	sign * (*(pfrom + inout_from * from->strides[dim0]));      
    }
  }
}


/** @brief Sets boundary conditions (Neumann/Dirichlet) reading from
 * the array itself. 
 *
 * For example, to set Neumann conditions at r = 0 (dim0 = R_INDX),
 * we call it as 
 *
 * rz_set_bnd (array, 1, RZTP (array, 0, 0, 0), R_INDX, 1, 
 *             Z_INDX, z0, z1, THETA_INDX, 0, ntheta);
 */
void
rz_set_bnd (rz_array_t *array, int sign, REAL *astart, int dim0, int inout,
	    int dim1, int dim1_0, int dim1_1, 
	    int dim2, int dim2_0, int dim2_1)
{
  rz_copy_bnd (array, array, sign, astart, astart, dim0, -inout, inout, 
	       dim1, dim1_0, dim1_1, dim2, dim2_0, dim2_1); 
}


/** @brief Frees an array */
void
rz_free (rz_array_t *array)
{
  debug (3, "rz_free\n");

  if (NULL == array->host) free (array->data);
  free (array);
}

/** @brief Copies a sub-array of @a rn x @a zn from @a fro
 * (starting at (@a rfro, z@a fro)) to @a to (starting at
 * (@a rto, @a zto)). */
void 
rz_copy (rz_array_t *fro, int rfro, int zfro, 
	     rz_array_t *to, int rto, int zto,
	     int rn, int zn)
{
  int i;
  REAL *pfro, *pto;
  
  debug (3, "rz_copy\n");

  pfro = RZP (fro, rfro, zfro);
  pto = RZP (to, rto, zto);
  
  for (i = 0; i < zn; i++) {
    memcpy (pto, pfro, sizeof(REAL) * rn);

    pfro += fro->strides[Z_INDX];
    pto += to->strides[Z_INDX];
  }
}

/** @brief Copies a sub-array of the mode nmode of @a rn x @a zn from
 * @a fro (starting at (@a rfro, @a zfro)) to @a to (starting at
 * (@a rto, @a zto)). 
*/
void 
rz_copy_modes (rz_array_t *fro, int rfro, int zfro, 
	       rz_array_t *to, int rto, int zto,
	       int rn, int zn, int nmode_fro, int nmode_to)
{
  int i, j;
  REAL *pfro, *pto;
  
  debug (3, "rz_copy\n");

  pfro = RZTP (fro, rfro, zfro, nmode_fro);
  pto = RZTP (to, rto, zto, nmode_to);
  
  for (i = 0; i < zn; i++) {
    /* memcpy (pto, pfro, sizeof(REAL) * rn); Maybe not thread safe? */
    for (j = 0; j < rn; j++ ) *(pto + j) = *(pfro + j);

    pfro += fro->strides[Z_INDX];
    pto += to->strides[Z_INDX];
  }
}


/** @brief Reads or writes a @a rz_array, depending on the string mode,
 * "r" for reading, "w" for writing.
 *
 * Note however that we must already have the dimensions of the array in it.
 */
void
rz_dump (rz_array_t *rz_array, const char *fname, const char *mode,
	 int r0, int z0, int r1, int z1)
{
  FILE *fp;
  int ir, iz;
  int writing;
  double f;

  debug (3, "rz_dump(\"%s\", \"%s\", %d, %d, %d, %d)\n", fname, mode, 
	 r0, z0, r1, z1);
 
  if (0 == strcmp (mode, "w")) {
    writing = TRUE;
  } else if (0 == strcmp (mode, "r")) {
    writing = FALSE;
  } else {
    fatal ("Unknown mode for rz_dump\n");
  }

  fp = fopen (fname, mode);

  if (fp == NULL) {
    warning("Unable to open %s for %s\n", fname, 
	    writing? "writing": "reading");
    return;
  }

  for (ir = r0; ir < r1; ir++)
    for (iz = z0; iz < z1; iz++) {
      if (writing)
	fprintf (fp, "%15.5e\n", RZ(rz_array, ir, iz));
      else {
	if (fscanf (fp, "%lf", &f) != 1) {
	  warning ("Error reading file %s, at line %d, (ir, iz) = (%d, %d)\n",
		   fname, iz + ir * (z1 - z0), ir, iz);
	  f = 0.0;
	}
	*RZP (rz_array, ir, iz) = f;
      }
    }
  
  used_disk_space += ftell (fp);
  fclose(fp);

  check_disk_space ();
}

/** @brief rz_dump_3d QQQQ */
void
rz_dump_3d (rz_array_t *rz_array, const char *fname, const char *mode,
	    int r0, int z0, int r1, int z1, int ntheta)
{
  FILE *fp;
  int ir, iz, itheta;
  int writing;
  double f;

  debug (3, "rz_dump_3d(\"%s\", \"%s\", %d, %d, %d, %d, %d)\n", 
	 fname, mode, r0, z0, r1, z1, ntheta);
 
  if (0 == strcmp (mode, "w")) {
    writing = TRUE;
  } else if (0 == strcmp (mode, "r")) {
    writing = FALSE;
  } else {
    fatal ("Unknown mode for rz_dump_3d\n");
  }

  fp = fopen (fname, mode);

  if (fp == NULL) {
    warning ("Unable to open %s for %s\n", fname, 
	     writing? "writing": "reading");
    return;
  }

  for (itheta = 0; itheta < ntheta; itheta++) {
    for (ir = r0; ir < r1; ir++) {
      for (iz = z0; iz < z1; iz++) {
	if (writing)
	  fprintf (fp, "%15.5e\n", RZT (rz_array, ir, iz, itheta));
	else {
	  if (fscanf (fp, "%lf", &f) != 1) {
	    warning ("Error reading file %s, at line %d, "
		     "(ir, iz, itheta) = (%d, %d, %d)\n",
		     fname, iz + (ir + itheta * (r1 - r0)) * (z1 - z0), 
		     ir, iz, itheta);
	    f = 0.0;
	  }
	  *RZTP (rz_array, ir, iz, itheta) = f;
	}
      }
    }
  }
  
  /* Add the disk space used to the total disk space that we have used so far
   */
  
  used_disk_space += ftell (fp);
  fclose(fp);

  check_disk_space ();
}
 
/** @brief rz_axis_dump QQQQ */
void
 rz_axis_dump (const char *fname, int x0, int x1, double delta)
{
  FILE *fp;
  int i;

  debug (3, "rz_axis_dump(\"%s\", %d, %d, %f)\n", fname, x0, x1, delta);
 
  fp = fopen (fname, "w");

  if (fp == NULL) {
    warning ("Unable to open %s\n", fname);
    return;
  }

  for (i = x0; i < x1; i++)
    fprintf (fp, "%15.5e\n", ((double) i + 0.5) * delta);


  used_disk_space += ftell (fp);
  fclose(fp);

  check_disk_space ();
}

/** @brief check_disk_space QQQQ */
static void
check_disk_space (void)
{
  if (used_disk_space > (((long int) max_disk_space_mb) << 20)) {
    fatal ("The disk space limit has been surpassed. "
	   "Increase max_disk_space_mb if you really need more space\n");
  }
}

