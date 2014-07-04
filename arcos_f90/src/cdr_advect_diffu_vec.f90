!> @file cdr_advect_diffu_vec.f90
!! @brief This code replaces the original c code of cdr_advect_diffu in file
!!  cdr.c
!!
!! @author Margreet Nool
!!
!!     The arrays er_array and ez_array are in C stored as
!!     1D arrays. In the wrapper routine cdr_advect_diffu_wrap.c it is
!!     checked of their sizes, er_len, ez_len are equal to rmax-by-zmax.
!!
!!     Note that these are NOT the r0, z0 we received.  The idea here is that
!!     RZ(array, r0, z0) will produce the correct result, where the allocated
!!     array looks like (in the picture, r0, z0 ARE the received ones).
!!     Note: in the picture there is only one buffer cell, while we actually
!!     allocate two.
!!
!!      +--------------+--------------+...+--------------+--------------+
!!      |r0 - 1, z0 - 1|              |   |              |              |
!!      +--------------+--------------+...+--------------+--------------+
!!      |              |   r0, z0     |   |              |              |
!!      +--------------+--------------+...+--------------+--------------+
!!                    ...             |...|             ...
!!      +--------------+--------------+...+--------------+--------------+
!!      |              |              |   |r1 - 1, z1 - 1|              |
!!      +--------------+--------------+...+--------------+--------------+
!!      |              |              |   |              |    r1, z1    |
!!      +--------------+--------------+...+--------------+--------------+
!!
!!     but the rows z0 +/- 1 and the columns r0 +/- 1 do not belong to the
!!     physical space and are used only to specify boundary conditions.
!!
!!     In the C routines r0, and z0 are zero, and the calculations start 
!!     at the third element (r0-2, r0-1, r0, ...., r1, r1+1), and
!!     (z0-2, z0-1, r0, ...., z1, z1+1), respectively. 
!!     rmax = r1 - r0 + 4 and zmax = z1 - z0 + 4.
!!-----------------------------------------------

! subroutine cdr_advect_vec
! subroutine mat_efield_r
! subroutine mat_efield_z
! subroutine mat_cdr_f_ad_r
! subroutine mat_cdr_f_ad_z
! subroutine make_vecs_er

   subroutine cdr_advect_vec(mass,charge,dens_array,d_dens_array, &
                         er_array,ez_array,diffusion_coeff,dr,dz, &
                         sprite_module,r0,r1,z0,z1)
      implicit none
!     ..
!     .. Parameters ..
      double precision,parameter :: ONE  = 1.0D0
!     ..
!     .. Dummy Variables ..
      logical,intent(in)          :: sprite_module
      integer,intent(in)          :: r0,r1,z0,z1
      double precision,intent(in) :: charge,dr,dz,mass,diffusion_coeff
!     ..
!     .. Dummy Arrays ..
      double precision,dimension(r0-2:r1+1,z0-2:z1+1),intent(in)  :: &
                       dens_array,er_array,ez_array
      double precision,dimension(r0-2:r1+1,z0-2:z1+1),intent(out) :: &
                       d_dens_array
!     ..
!     .. Local Variables ..
      integer          :: ir,ishift
      double precision :: cm0,cm,gdz_inv,d_iso0,d_iso
!     ..
!     .. Local arrays ..
      double precision,dimension(r0:r1-1)             :: er_min
      double precision,dimension(r0+1:r1)             :: er
      double precision,dimension(z1-z0)               :: d_iso_r
      double precision,dimension(z1-z0+1)             :: d_iso_z
      double precision,dimension(r1-r0+1,z1-z0  )     :: mat_r
      double precision,dimension(r1-r0  ,z1-z0+1)     :: mat_z
      double precision,dimension(r1-r0+1,z1-z0  )     :: efield_r
      double precision,dimension(r0:r1-1,z0-1:z1-1)   :: efield_z
!     ..
!     .. External subroutines ..
      external :: mat_efield_r,mat_efield_z,mat_cdr_f_ad_r,mat_cdr_f_ad_z
!     ..

      if (mass <= 0) return
  
      ! <ISOTROPIC DIFFUSSION>
      d_iso0 = diffusion_coeff / mass
      ! </ISOTROPIC DIFFUSION>

      gdz_inv = ONE / dz
  
      cm0 = charge / mass
      
      if (.not.sprite_module) then
        cm = cm0
        d_iso = d_iso0
      end if
  
!-----------------------------------------------
!     Increase the values r0, r1, z0 and z1 with 3:
!     because of 2 buffer cells at the left side and to compensate the C 0
!-----------------------------------------------
      ishift = 3

      ! Compute the electric field in radial direction
      call mat_efield_r(er_array(r0-1:r1-1,z0:z1-1), &
                        efield_r,d_iso_r,cm0,d_iso0, &
                        r0,r1,dz,z0,z1,ishift,dr,sprite_module)
      ! Calculate the radial advection from left and right 
      call mat_cdr_f_ad_r(dens_array(r0-2:r1+1,z0:z1-1), &
                          efield_r(1:r1-r0+1,1:z1-z0),d_iso_r, &
                          mat_r,r0,r1,z0,z1)
      ! Add the radial advection 
      call  make_vecs_er (er_min,er,r0,r1,dr)
      do ir = r0,r1-1
         d_dens_array(ir,z0:z1-1) = d_dens_array(ir,z0:z1-1) + &
                                    er_min(ir) * mat_r(ir-r0+1,1:z1-z0) - &
                                    er(ir+1)   * mat_r(ir-r0+2,1:z1-z0)
      end do

      ! Compute the electric field in axial direction
      call mat_efield_z(ez_array(r0  :r1-1,z0-1:z1-1), &
                        efield_z(r0:r1-1,z0-1:z1-1),d_iso_z,cm0,d_iso0, &
                        r0,r1,dz,z0,z1,ishift,dz,sprite_module)
      ! Calculate the axial advection from the north and south
      call mat_cdr_f_ad_z(dens_array(r0  :r1-1,z0-2:z1+1), &
                          efield_z(r0:r1-1,z0-1:z1-1),d_iso_z, &
                          mat_z,r0,r1,z0,z1)
      ! Add the axial advection 
      d_dens_array(r0:r1-1,z0:z1-1) = d_dens_array(r0:r1-1,z0:z1-1) + &
                                      gdz_inv * (mat_z(1:,1:z1-z0) - &
                                                 mat_z(1:,2:z1-z0+1))

   end subroutine cdr_advect_vec

   subroutine mat_efield_r(er_array,efield,d_iso,cm0,d_iso0, &
                           r0,r1,dz,z0,z1,z_shift,dx,sprite_module)
!     ..
      implicit none
!     ..
!     .. Parameters ..
      double precision,parameter :: HALF = 0.5D0,ONE  = 1.0D0
!     ..
!     .. Dummy Variables ..
      logical,intent(in)          :: sprite_module
      integer,intent(in)          :: r0,r1,z0,z1,z_shift
      double precision,intent(in) :: cm0,d_iso0,dz,dx
!     ..
!     .. Dummy Arrays ..
      double precision,dimension(r1-r0+1,z1-z0),intent(in)  :: er_array
      double precision,dimension(r1-r0+1,z1-z0),intent(out) :: efield
      double precision,dimension(z1-z0),intent(out)         :: d_iso
!     ..
!     .. Local Variables ..
      integer                           :: iz
!     ..
!     .. Local arrays ..
      double precision,dimension(z1-z0) :: back_dens_inv,cm
!-----------------------------------------------

      ! If the densities are varying...
      if (sprite_module) then
         do iz = z0,z1
            ! spr_density_at(iz) = (iz-z_shift + half) * dz
            ! back_dens_inv = ONE / spr_density_at
            back_dens_inv(iz-z0+1) = ONE / ((iz-z_shift + half) * dz)
         end do
         cm    = cm0 * back_dens_inv
         d_iso = d_iso0 * back_dens_inv
      else
         back_dens_inv = ONE
         cm            = cm0
         d_iso         = d_iso0
      end if

      d_iso = d_iso / dx

      ! efield = cm * er_array
      do iz = 1,z1-z0
         efield(:,iz) = cm(iz) * er_array(:,iz)
      end do

   end subroutine mat_efield_r

   subroutine mat_efield_z (ez_array,efield,d_iso,cm0,d_iso0, &
                            r0,r1,dz,z0,z1,z_shift,dx,sprite_module)
!     ..
      implicit none
!     ..
!     .. Parameters ..
      double precision,parameter :: HALF = 0.5D0,ONE  = 1.0D0
!     ..
!     .. Dummy Variables ..
      logical,intent(in)          :: sprite_module
      integer,intent(in)          :: r0,r1,z0,z1,z_shift
      double precision,intent(in) :: cm0,d_iso0,dz,dx
!     ..
!     .. Dummy Arrays ..
      double precision,dimension(r0:r1-1,z0-1:z1-1),intent(in)  :: ez_array
      double precision,dimension(r0:r1-1,z0-1:z1-1),intent(out) :: efield
      double precision,dimension(z0-1:z1-1),intent(out)     :: d_iso
!     ..
!     .. Local Variables ..
      integer                               :: iz
!     ..
!     .. Local arrays ..
      double precision,dimension(z0-1:z1-1) :: back_dens_inv,cm
!     ..
!-----------------------------------------------

      ! If the densities are varying...
      if (sprite_module) then
         do iz = z0-1,z1+1
            ! spr_density_at(iz) = (iz-z_shift + half) * dz
            ! back_dens_inv = ONE / spr_density_at
            back_dens_inv(iz) = ONE / ((iz-z_shift - half) * dz)
         end do
         cm    = cm0 * back_dens_inv
         d_iso = d_iso0 * back_dens_inv
      else
         back_dens_inv = ONE
         cm            = cm0
         d_iso         = d_iso0
      end if

      d_iso = d_iso / dx

      ! efield = cm * ez_array
      do iz = z0-1,z1-1
         efield(:,iz) = cm(iz) * ez_array(:,iz)
      end do

   end subroutine mat_efield_z

   subroutine mat_cdr_f_ad_r(data_array,efield,d_iso,mat_r,r0,r1,z0,z1)
!     ..
      implicit none
!     ..
!     .. Parameters ..
      double precision,parameter :: ZERO = 0.0D0,HALF = 0.5D0, &
                                    ONE  = 1.0D0,THREE = 3.0D0, &
                                    verysmall = 1.0e-20
!     ..
!     .. Dummy Variables ..
      integer,intent(in) :: r0,r1,z0,z1
!     ..
!     .. Dummy Arrays ..
      double precision,dimension(z1-z0),intent(in)          :: d_iso
      double precision,dimension(r1-r0+1,z1-z0),intent(in)  :: efield
      double precision,dimension(r1-r0+4,z1-z0),intent(in)  :: data_array
      double precision,dimension(r1-r0+1,z1-z0),intent(out) :: mat_r
!     ..
!     .. Local Variables ..
      integer            :: iz,length
!     ..
!     .. Local arrays ..
      double precision,dimension(r1-r0+1,z1-z0) :: aux,psi_p
      double precision,dimension(r1-r0+3,z1-z0) :: mat_diff
!     ..
!-----------------------------------------------
      length = r1 - r0 + 1
      mat_diff(1:length+2,:) = data_array(2:length+3,:) - &
                               data_array(1:length+2,:)

      where (abs(mat_diff(2:length+1,:)) > verysmall)
          where (efield>=ZERO)
            ! p_min  =data_array(ir-k1,iz-k2)
            ! p      =data_array(ir,iz)
            ! p_plus =data_array(ir+k1,iz+k2)
   
             psi_p = mat_diff(1:length,:) / mat_diff(2:length+1,:)
         elsewhere (efield<ZERO)
            ! p       =data_array(ir,iz)
            ! p_plus  =data_array(ir+  k1,iz+  k2)
            ! p_2plus =data_array(ir+2*k1,iz+2*k2)
   
             psi_p = mat_diff(3:length+2,:) / mat_diff(2:length+1,:) 
         end where
         where (psi_p <= ZERO)
            ! psi_p = psi_MN (psi_p)
            psi_p = ZERO
         elsewhere (psi_p >= 4)
            psi_p = ONE
         elsewhere (psi_p >= 0.4)
            psi_p = (ONE + HALF * psi_p) / THREE
         end where
      elsewhere  ! abs(mat_diff(2:length+1,:)) < verysmall
         psi_p = -ONE
      end where

      ! sigmasigma is (sigma_{i+1,j} - sigma_{i,j})*/
      ! sigmasigma = p - p_plus 
      aux = psi_p * mat_diff(2:length+1,:)
      ! aux = p + psi_p * sigmasigma
      where (efield>=ZERO)
         aux = data_array(2:length+1,:) + aux
      elsewhere
         aux = data_array(3:length+2,:) - aux
      end where

      ! aux = efield * (p + psi_p * sigmasigma)
      aux = aux * efield
      do iz = 1,z1-z0
         ! mat_r = d_iso * sigmasigma
         mat_r(1:length,iz) = d_iso(iz) * mat_diff(2:length+1,iz)
      end do
      ! mat_r =  efield * (p + psi_p * sigmasigma) + d_iso * sigmasigma
      mat_r = aux - mat_r

   end subroutine mat_cdr_f_ad_r

   subroutine mat_cdr_f_ad_z(data_array,efield,d_iso,mat_z,r0,r1,z0,z1)
!     ..
      implicit none
!     ..
!     .. Parameters ..
      double precision,parameter :: ZERO = 0.0D0,HALF = 0.5D0, &
                                    ONE  = 1.0D0,THREE = 3.0D0, &
                                    verysmall = 1.0e-20
!     ..
!     .. Dummy Variables ..
      integer,intent(in)          :: r0,r1,z0,z1
!     ..
!     .. Dummy Arrays ..
      double precision,dimension(z1-z0+1),intent(in)           :: d_iso
      double precision,dimension(r0:r1-1,z0-1:z1-1),intent(in) :: efield
      double precision,dimension(r1-r0,z1-z0+4),intent(in)     :: data_array
      double precision,dimension(r1-r0,z1-z0+1),intent(out)    :: mat_z
!     ..
!     .. Local Variables ..
      integer          :: iz,length
!     ..
!     .. Local Arrays ..
      double precision,dimension(r1-r0,z1-z0+1) :: aux,psi_p,pij
      double precision,dimension(r1-r0,z1-z0+3) :: mat_diff
!-----------------------------------------------
      length = z1 - z0 + 1
      mat_diff(:,1:length+2) = data_array(:,2:length+3) - &
                               data_array(:,1:length+2)

      where (abs(mat_diff(:,2:length+1)) > verysmall)
         ! psi_p = (p - p_min) / (p_plus- p)
         where (efield>ZERO)
            ! p_min  =data_array(ir-k1,iz-k2)
            ! p      =data_array(ir,iz)
            ! p_plus =data_array(ir+k1,iz+k2)
   
            psi_p = mat_diff(:,1:length) / mat_diff(:,2:length+1)
         elsewhere ! efield=<ZERO
            ! p       =data_array(ir,iz)
            ! p_plus  =data_array(ir+  k1,iz+  k2)
            ! p_2plus =data_array(ir+2*k1,iz+2*k2)

            psi_p = mat_diff(:,3:length+2) / mat_diff(:,2:length+1) 
         end where
         ! psi_p = psi_MN (psi_p)
         pij   = psi_p
         where ((psi_p > 0.4) .and. (psi_p < 4))
            pij =  (ONE + HALF * psi_p) / THREE
         elsewhere (psi_p >= 4)
            pij = ONE
         else where (psi_p <= ZERO)
            pij = ZERO
         end where
      elsewhere ! abs(mat_diff(:,2:length+1)) < verysmall
         pij = ONE
      end where

      ! sigmasigma is (sigma_{i+1,j} - sigma_{i,j})*/
      ! sigmasigma = p - p_plus 
      aux = pij * mat_diff(:,2:length+1)
      ! aux = p + pij * sigmasigma
      where (efield>=0)
         aux = data_array(:,2:length+1) + aux
      elsewhere
         aux = data_array(:,3:length+2) - aux
      end where

      ! aux = efield * (p + pij * sigmasigma)
      aux = aux * efield
      do iz = 1,length
         ! mat_z = d_iso * sigmasigma
         mat_z(:,iz) = d_iso(iz) * mat_diff(:,iz+1)
      end do
      ! mat_z =  efield * (p + pij * sigmasigma) + d_iso * sigmasigma
      mat_z = aux - mat_z

   end subroutine mat_cdr_f_ad_z

   subroutine make_vecs_er (er_min,er,r0,r1,dr)
!     ..
      implicit none
!
!     .. Parameters ..
      double precision,parameter :: HALF = 0.5D0,ONE  = 1.0D0
!     ..
!     .. Dummy Variables ..
      integer,intent(in)                              :: r0,r1
      double precision,intent(in)                     :: dr
!     ..
!     .. Dummy Arrays ..
      double precision,dimension(r0:r1-1),intent(out) :: er_min
      double precision,dimension(r0+1:r1),intent(out) :: er
!     ..
!     .. Local Variables ..
      integer          :: ir
      double precision :: r_inv
!     ..
! .............................................................
      do ir = r0,r1-1
         r_inv = ONE / ((ir+ HALF) * dr)
         er_min(ir) = ir * r_inv
         er(ir+1)   = (ir+1) * r_inv
      end do

   end subroutine make_vecs_er
