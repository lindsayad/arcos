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
         write(*,*) 'warning', &
                    'tri3: l1,l2,l3,kint1,kint2,kint3 uninitialized'
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
