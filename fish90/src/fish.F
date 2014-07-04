!
!     file fish.f
!
!  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!  .                                                             .
!  .                  copyright (c) 2004 by UCAR                 .
!  .                                                             .
!  .       UNIVERSITY CORPORATION for ATMOSPHERIC RESEARCH       .
!  .                                                             .
!  .                      all rights reserved                    .
!  .                                                             .
!  .                                                             .
!  .                      FISHPACK version 5.0                   .
!  .                                                             .
!  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                        F I S H P A C K                        *
!     *                                                               *
!     *                                                               *
!     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
!     *                                                               *
!     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
!     *                                                               *
!     *                  (Version 5.0 , JUNE 2004)                    *
!     *                                                               *
!     *                             BY                                *
!     *                                                               *
!     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
!     *                                                               *
!     *                             OF                                *
!     *                                                               *
!     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
!     *                                                               *
!     *                BOULDER, COLORADO  (80307)  U.S.A.             *
!     *                                                               *
!     *                   WHICH IS SPONSORED BY                       *
!     *                                                               *
!     *              THE NATIONAL SCIENCE FOUNDATION                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!

!     this module is used by all fishpack solvers to allocate
!     real and complex work space
      MODULE fish

         TYPE fishworkspace
            DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: rew
            COMPLEX,DIMENSION(:),ALLOCATABLE :: cxw
         END TYPE fishworkspace

      CONTAINS
         SUBROUTINE allocatfish(irwk,icwk,wsave,ierror)
         IMPLICIT NONE
         TYPE (fishworkspace) :: wsave
!        irwk is the required real work space length
!        icwk is the required integer work space length
         INTEGER, INTENT(IN) :: irwk,icwk
!        ierror is set to 20 if the dynamic allocation is unsuccessful
!        (e.g., this would happen if m,n are too large for the computers memory
         INTEGER, INTENT(INOUT) :: ierror
         INTEGER :: istatus
!        first deallocate to avoid memory leakage
         write(*,*) 'fish.F: allocatfish'
         if(allocated(wsave%rew)) then
            write(*,*) 'allocated wsave%rew, size= ',size(wsave%rew)
            DEALLOCATE(wsave%rew,STAT=istatus)
            write(*,*) 'deallocate wsave%rew, istatus=',istatus
         end if
!        allocate irwk words of real work space
         if (irwk > 0) then
            ALLOCATE(wsave%rew(irwk),STAT = istatus)
            write(*,*) 'fish.F: allocatfish, irwk:',irwk
         end if
!
         if(allocated(wsave%cxw)) then
            write(*,*) 'allocated wsave%cxw, size=',size(wsave%cxw)
            DEALLOCATE(wsave%cxw,STAT=istatus)
            write(*,*) 'deallocate wsave%cxw, istatus=',istatus
         end if
!        allocate icwk words of complex work space
         if (icwk > 0) then
            ALLOCATE(wsave%cxw(icwk),STAT = istatus)
            write(*,*) 'fish.F: allocatfish, icrk:',icwk
         end if
         ierror = 0
!        flag fatal error if allocation fails
         if (istatus .ne. 0 ) then
            ierror = 20
         END IF
         RETURN
        END SUBROUTINE allocatfish

        SUBROUTINE BLK_space(N,M,irwk,icwk)
!       this subroutine computes the real and complex work space
!       requirements (generous estimate) of blktri for N,M values
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: N,M
        INTEGER,INTENT(OUT) :: irwk,icwk
        INTEGER :: L,log2n
!       compute nearest integer greater than or equal to
!       log base 2 of n+1, i.e., log2n is smallest integer
!       such that 2**log2n >= n+1
        log2n = 1
        do
           log2n = log2n+1
           if (n+1 <= 2**log2n) EXIT
        end do
        L = 2**(log2n+1)
        irwk = (log2n-2)*L+5+MAX0(2*N,6*M)+log2n+2*n
        icwk = ((log2n-2)*L+5+log2n)/2+3*M+N
        RETURN
        END SUBROUTINE BLK_space

        SUBROUTINE GEN_space(N,M,irwk)
!       this subroutine computes the real work space
!       requirement (generously) of genbun for the current N,M
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: N,M
        INTEGER,INTENT(OUT) :: irwk
        INTEGER :: log2n
!       compute nearest integer greater than or equal to
!       log base 2 of n+1, i.e., log2n is smallest integer
!       such that 2**log2n >= n+1
        log2n = 1
        do
           log2n = log2n+1
           if (n+1 <= 2**log2n) EXIT
        end do
        irwk = 4*N + (10 + log2n)*M
        RETURN
        END SUBROUTINE GEN_space

        SUBROUTINE fishfin(wsave,ierror)
!       this subroutine releases allocated work space
!       fishfin should be called after a fishpack solver has finished
!       TYPE (fishworkspace) variable wsave.
        IMPLICIT NONE
        TYPE (fishworkspace) :: wsave
        INTEGER,intent(out)  :: ierror
        INTEGER              :: istatus
!
        write(*,*) 'fish.F: fishfin'
        if(allocated(wsave%rew)) then
           write(*,*) 'deallocated wsave%rew, size= ',size(wsave%rew)
           DEALLOCATE(wsave%rew,STAT=istatus)
        end if
        if(allocated(wsave%cxw)) then
           write(*,*) 'deallocated wsave%cxw, size=',size(wsave%cxw)
           DEALLOCATE(wsave%cxw,STAT=istatus)
        end if
! #ifndef G95
         ! if(associated(wsave%rew))DEALLOCATE(wsave%rew)
         ! if(associated(wsave%cxw))DEALLOCATE(wsave%cxw)
! #endif
        if (istatus .ne. 0 ) then
           ierror = 20
        END IF

        END SUBROUTINE fishfin

      END MODULE fish
