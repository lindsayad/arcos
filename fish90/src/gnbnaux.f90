!
!     file gnbnaux.f
!
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
! PACKAGE GNBNAUX
!
! LATEST REVISION        June 2004
!
! PURPOSE                TO PROVIDE AUXILIARY ROUTINES FOR FISHPACK
!                        ENTRIES GENBUN AND POISTG.
!
! USAGE                  THERE ARE NO USER ENTRIES IN THIS PACKAGE.
!                        THE ROUTINES IN THIS PACKAGE ARE NOT INTENDED
!                        TO BE CALLED BY USERS, BUT RATHER BY ROUTINES
!                        IN PACKAGES GENBUN AND POISTG.
!
! SPECIAL CONDITIONS     NONE
!
! I/O                    NONE
!
! PRECISION              SINGLE
!
!
! LANGUAGE               FORTRAN 90
!
! HISTORY                WRITTEN IN 1979 BY ROLAND SWEET OF NCAR'S
!                        SCIENTIFIC COMPUTING DIVISION.  MADE AVAILABLE
!                        ON NCAR'S PUBLIC LIBRARIES IN JANUARY, 1980.
!                        Revised by John Adams in June 2004 incorporating
!                        Fortran 90 features
!
! PORTABILITY            FORTRAN 90
! ********************************************************************
      SUBROUTINE COSGEN(N, IJUMP, FNUM, FDEN, vecA, IA)
      implicit none
      
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, &
                                     ONE  = 1.0D0, TWO  = 2.0D0, &
                                     FOUR = 4.0D0

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)                         :: N, IJUMP, IA
      DOUBLE PRECISION, INTENT(IN)                :: FNUM, FDEN
      DOUBLE PRECISION, DIMENSION(IA),INTENT(OUT) :: vecA
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER          :: K3, K4, K, K1, K5, I, K2, NP1
      DOUBLE PRECISION :: PI, PIBYN, X, Y
!-----------------------------------------------
!
!
!     THIS SUBROUTINE COMPUTES REQUIRED COSINE VALUES IN ASCENDING
!     ORDER.  WHEN IJUMP .GT. 1 THE ROUTINE COMPUTES VALUES
!
!        2*COS(J*PI/L) , J=1,2,...,L AND J .NE. 0(MOD N/IJUMP+1)
!
!     WHERE L = IJUMP*(N/IJUMP+1).
!
!
!     WHEN IJUMP = 1 IT COMPUTES
!
!            2*COS((J-FNUM)*PI/(N+FDEN)) ,  J=1, 2, ... ,N
!
!     WHERE
!        FNUM = 0.5, FDEN = 0.0,  FOR REGULAR REDUCTION VALUES
!        FNUM = 0.0, FDEN = 1.0, FOR B-R AND C-R WHEN ISTAG = 1
!        FNUM = 0.0, FDEN = 0.5, FOR B-R AND C-R WHEN ISTAG = 2
!        FNUM = 0.5, FDEN = 0.5, FOR B-R AND C-R WHEN ISTAG = 2
!                                IN POISN2 ONLY.
!
!
      PI = 4.0*ATAN(1.0)
      IF (N /= 0) THEN
         IF (IJUMP /= 1) THEN
            K3 = N/IJUMP + 1
            K4 = K3 - 1
            PIBYN = PI/FLOAT(N + IJUMP)
            DO K = 1, IJUMP
               K1 = (K - 1)*K3
               K5 = (K - 1)*K4
               DO I = 1, K4
                  X = K1 + I
                  K2 = K5 + I
                  vecA(K2) = -2.*COS(X*PIBYN)
               END DO
            END DO
         ELSE
            NP1 = N + 1
            Y = PI/(FLOAT(N) + FDEN)
            DO I = 1, N
               X = FLOAT(NP1 - I) - FNUM
               vecA(I) = 2.*COS(X*Y)
            END DO
         ENDIF
      ENDIF
!
      END SUBROUTINE COSGEN

      SUBROUTINE MERGE(TCOS, I1, M1, I2, M2, I3, itcos)
      implicit none

      
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, &
                                     ONE  = 1.0D0, TWO  = 2.0D0, &
                                     FOUR = 4.0D0

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: I1, M1, I2, M2, I3, ITCOS
      DOUBLE PRECISION, DIMENSION(ITCOS), INTENT(INOUT) :: TCOS
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER          :: J11, J3, J1, J2, J, L, K, M
      DOUBLE PRECISION :: X, Y
!-----------------------------------------------
!
!     THIS SUBROUTINE MERGES TWO ASCENDING STRINGS OF NUMBERS IN THE
!     ARRAY TCOS.  THE FIRST STRING IS OF LENGTH M1 AND STARTS AT
!     TCOS(I1+1).  THE SECOND STRING IS OF LENGTH M2 AND STARTS AT
!     TCOS(I2+1).  THE MERGED STRING GOES INTO TCOS(I3+1).
!
!
      J1 = 1
      J2 = 1
      J = I3
      IF (M1 == 0) GO TO 107
      IF (M2 == 0) GO TO 104
  101 CONTINUE
      J11 = J1
      J3 = MAX(M1,J11)
      DO J1 = J11, J3
         J = J + 1
         L = J1 + I1
         X = TCOS(L)
         L = J2 + I2
         Y = TCOS(L)
         IF (X - Y > 0.) GO TO 103
         TCOS(J) = X
      END DO
      GO TO 106
  103 CONTINUE
      TCOS(J) = Y
      J2 = J2 + 1
      IF (J2 <= M2) GO TO 101
      IF (J1 > M1) GO TO 109
  104 CONTINUE
      K = J - J1 + 1
      DO J = J1, M1
         M = K + J
         L = J + I1
         TCOS(M) = TCOS(L)
      END DO
      GO TO 109
  106 CONTINUE
      IF (J2 > M2) GO TO 109
  107 CONTINUE
      K = J - J2 + 1
      DO J = J2, M2
         M = K + J
         L = J + I2
         TCOS(M) = TCOS(L)
      END DO
  109 CONTINUE
!
      END SUBROUTINE MERGE


      SUBROUTINE TRIX(IDEGBR, IDEGCR, M, vecA, vecB, vecC, vecY, TCOS, ITCOS,  &
                      vecD, vecW)
      implicit none
      
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, &
                                     ONE  = 1.0D0, TWO  = 2.0D0, &
                                     FOUR = 4.0D0

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: IDEGBR, IDEGCR, M, ITCOS
      DOUBLE PRECISION, DIMENSION(M), INTENT(IN)     :: vecA, vecB, vecC
      DOUBLE PRECISION, DIMENSION(ITCOS), INTENT(IN) :: TCOS
      DOUBLE PRECISION, DIMENSION(M), INTENT(INOUT)  :: vecY, vecD, vecW
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER          :: MM1, IFB, IFC, L, LINT, K, I, IP
      DOUBLE PRECISION :: X, XX, Z
!-----------------------------------------------
!
!     SUBROUTINE TO SOLVE A SYSTEM OF LINEAR EQUATIONS WHERE THE
!     COEFFICIENT MATRIX IS A RATIONAL FUNCTION IN THE MATRIX GIVEN BY
!     TRIDIAGONAL  ( . . . , vecA(I), vecB(I), vecC(I), . . . ).
!
      MM1 = M - 1
      IFB = IDEGBR + 1
      IFC = IDEGCR + 1
      L = IFB/IFC
      LINT = 1
      DO K = 1, IDEGBR
         X = TCOS(K)
         IF (K == L) THEN
            I = IDEGBR + LINT
            XX = X - TCOS(I)
            vecW(:M) = vecY(:M)
            vecY(:M) = XX*vecY(:M)
         ENDIF
         Z = 1./(vecB(1)-X)
         vecD(1) = vecC(1)*Z
         vecY(1) = vecY(1)*Z
         DO I = 2, MM1
            Z = 1./(vecB(I)-X-vecA(I)*vecD(I-1))
            vecD(I) = vecC(I)*Z
            vecY(I) = (vecY(I)-vecA(I)*vecY(I-1))*Z
         END DO
         Z = vecB(M) - X - vecA(M)*vecD(MM1)
         IF (Z == 0.) THEN
            vecY(M) = 0.
         ELSE
            vecY(M) = (vecY(M)-vecA(M)*vecY(MM1))/Z
         ENDIF
         DO IP = 1, MM1
            vecY(M-IP) = vecY(M-IP) - vecD(M-IP)*vecY(M+1-IP)
         END DO
         IF (K /= L) CYCLE 
         vecY(:M) = vecY(:M) + vecW(:M)
         LINT = LINT + 1
         L = (LINT*IFB)/IFC
      END DO
!
      END SUBROUTINE TRIX


      SUBROUTINE TRI3(M, vecA, vecB, vecC, ivecK, vecY1, vecY2, vecY3, TCOS, &
                      ITCOS, vecD, vecW1, vecW2, vecW3)
      implicit none
      
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, &
                                     ONE  = 1.0D0, TWO  = 2.0D0, &
                                     FOUR = 4.0D0

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: M,ITCOS
      INTEGER, DIMENSION(4),INTENT(IN) :: ivecK
      DOUBLE PRECISION, DIMENSION(M), INTENT(IN)     :: vecA, vecB, vecC
      DOUBLE PRECISION, DIMENSION(ITCOS), INTENT(IN) :: TCOS
      DOUBLE PRECISION, DIMENSION(M), INTENT(INOUT) :: vecY1, vecY2, vecY3, &
                                                       vecD, vecW1, vecW2, vecW3
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MM1, K1, K2, K3, K4, IF1, IF2, IF3, IF4, K2K3K4, &
                 L1, L2, L3, LINT1, LINT2, LINT3, KINT1, KINT2, KINT3, &
                 N, I, IP
      DOUBLE PRECISION :: X, Z, XX
!-----------------------------------------------
!
!     SUBROUTINE TO SOLVE THREE LINEAR SYSTEMS WHOSE COMMON COEFFICIENT
!     MATRIX IS A RATIONAL FUNCTION IN THE MATRIX GIVEN BY
!
!                  TRIDIAGONAL (...,vecA(I),vecB(I),vecC(I),...)
!
      MM1 = M - 1
      K1 = ivecK(1)
      K2 = ivecK(2)
      K3 = ivecK(3)
      K4 = ivecK(4)
      IF1 = K1 + 1
      IF2 = K2 + 1
      IF3 = K3 + 1
      IF4 = K4 + 1
      K2K3K4 = K2 + K3 + K4
      IF (K2K3K4 /= 0) THEN
         L1 = IF1/IF2
         L2 = IF1/IF3
         L3 = IF1/IF4
         LINT1 = 1
         LINT2 = 1
         LINT3 = 1
         KINT1 = K1
         KINT2 = KINT1 + K2
         KINT3 = KINT2 + K3
      ELSE
         write(*,*) 'warning tri3: l1,l2,l3,kint1,kint2,kint3 uninitialized'
         stop 'stop in tri3: l1,l2,l3,kint1,kint2,kint3 uninitialized'
      ENDIF
      DO N = 1, K1
         X = TCOS(N)
         IF (K2K3K4 /= 0) THEN
            IF (N == L1) THEN
               vecW1(:M) = vecY1(:M)
            ENDIF
            IF (N == L2) THEN
               vecW2(:M) = vecY2(:M)
            ENDIF
            IF (N == L3) THEN
               vecW3(:M) = vecY3(:M)
            ENDIF
         ENDIF
         Z = 1./(vecB(1)-X)
         vecD(1) = vecC(1)*Z
         vecY1(1) = vecY1(1)*Z
         vecY2(1) = vecY2(1)*Z
         vecY3(1) = vecY3(1)*Z
         DO I = 2, M
            Z = 1./(vecB(I)-X-vecA(I)*vecD(I-1))
            vecD(I) = vecC(I)*Z
            vecY1(I) = (vecY1(I)-vecA(I)*vecY1(I-1))*Z
            vecY2(I) = (vecY2(I)-vecA(I)*vecY2(I-1))*Z
            vecY3(I) = (vecY3(I)-vecA(I)*vecY3(I-1))*Z
         END DO
         DO IP = 1, MM1
            vecY1(M-IP) = vecY1(M-IP) - vecD(M-IP)*vecY1(M+1-IP)
            vecY2(M-IP) = vecY2(M-IP) - vecD(M-IP)*vecY2(M+1-IP)
            vecY3(M-IP) = vecY3(M-IP) - vecD(M-IP)*vecY3(M+1-IP)
         END DO
         IF (K2K3K4 == 0) CYCLE 
         IF (N == L1) THEN
            I = LINT1 + KINT1
            XX = X - TCOS(I)
            vecY1(:M) = XX*vecY1(:M) + vecW1(:M)
            LINT1 = LINT1 + 1
            L1 = (LINT1*IF1)/IF2
         ENDIF
         IF (N == L2) THEN
            I = LINT2 + KINT2
            XX = X - TCOS(I)
            vecY2(:M) = XX*vecY2(:M) + vecW2(:M)
            LINT2 = LINT2 + 1
            L2 = (LINT2*IF1)/IF3
         ENDIF
         IF (N /= L3) CYCLE 
         I = LINT3 + KINT3
         XX = X - TCOS(I)
         vecY3(:M) = XX*vecY3(:M) + vecW3(:M)
         LINT3 = LINT3 + 1
         L3 = (LINT3*IF1)/IF4
      END DO
      RETURN 
!
! REVISION HISTORY---
!
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! OCTOBER   1980    CHANGED SEVERAL DIVIDES OF FLOATING INTEGERS
!                   TO INTEGER DIVIDES TO ACCOMODATE CRAY-1 ARITHMETIC.
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
!-----------------------------------------------------------------------
      END SUBROUTINE TRI3
