!
!     file hstcrt.f
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
!     SUBROUTINE HSTCRT (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
!    +                   ELMBDA,F,IDIMF,PERTRB,IERROR)
!
! DIMENSION OF           BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N)
! ARGUMENTS
!
! LATEST REVISION        June 2004
!
! PURPOSE                 SOLVES THE STANDARD FIVE-POINT FINITE
!                         DIFFERENCE APPROXIMATION TO THE HELMHOLTZ
!                         EQUATION
!                           (D/DX)(DU/DX) + (D/DY)(DU/DY) + LAMBDA*U
!                           = F(X,Y)
!                         ON A STAGGERED GRID IN CARTESIAN COORDINATES.
!
! USAGE                   CALL HSTCRT (A,B,M,MBDCND,BDA,BDB,C,D
!                                      N,NBDCND,BDC,BDD,ELMBDA,
!                                      F,IDIMF,PERTRB,IERROR)
!
! ARGUMENTS
! ON INPUT
!
!                        A,B
!                          THE RANGE OF X, I.E. A .LE. X .LE. B.
!                          A MUST BE LESS THAN B.
!
!                        M
!                          THE NUMBER OF GRID POINTS IN THE
!                          INTERVAL (A,B).  THE GRID POINTS
!                          IN THE X-DIRECTION ARE GIVEN BY
!                          X(I) = A + (I-0.5)DX FOR I=1,2,...,M
!                          WHERE DX =(B-A)/M.  M MUST BE GREATER
!                          THAN 2.
!
!                        MBDCND
!                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
!                          AT X = A AND X = B.
!
!                          = 0  IF THE SOLUTION IS PERIODIC IN X,
!                               U(M+I,J) = U(I,J).
!
!                          = 1  IF THE SOLUTION IS SPECIFIED AT
!                               X = A AND X = B.
!
!                          = 2  IF THE SOLUTION IS SPECIFIED AT
!                               X = A AND THE DERIVATIVE
!                               OF THE SOLUTION WITH RESPECT TO X
!                               IS SPECIFIED AT X = B.
!
!                          = 3  IF THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO X IS SPECIFIED
!                               AT X = A  AND X = B.
!
!                          = 4  IF THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO X IS SPECIFIED
!                               AT X = A  AND THE SOLUTION IS
!                               SPECIFIED AT X = B.
!
!                        BDA
!                          A ONE-DIMENSIONAL ARRAY OF LENGTH N
!                          THAT SPECIFIES THE BOUNDARY VALUES
!                          (IF ANY) OF THE SOLUTION AT X = A.
!
!                          WHEN MBDCND = 1 OR 2,
!                            BDA(J) = U(A,Y(J)) ,         J=1,2,...,N.
!
!                          WHEN MBDCND = 3 OR 4,
!                            BDA(J) = (D/DX)U(A,Y(J)) ,   J=1,2,...,N.
!
!                        BDB
!                          A ONE-DIMENSIONAL ARRAY OF LENGTH N
!                          THAT SPECIFIES THE BOUNDARY VALUES
!                          OF THE SOLUTION AT X = B.
!
!                          WHEN MBDCND = 1 OR 4
!                            BDB(J) = U(B,Y(J)) ,        J=1,2,...,N.
!
!                          WHEN MBDCND = 2 OR 3
!                            BDB(J) = (D/DX)U(B,Y(J)) ,  J=1,2,...,N.
!
!                        C,D
!                          THE RANGE OF Y, I.E. C .LE. Y .LE. D.
!                          C MUST BE LESS THAN D.
!
!
!                        N
!                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
!                          (C,D).  THE UNKNOWNS IN THE Y-DIRECTION
!                          ARE GIVEN BY Y(J) = C + (J-0.5)DY,
!                          J=1,2,...,N, WHERE DY = (D-C)/N.
!                          N MUST BE GREATER THAN 2.
!
!                        NBDCND
!                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
!                          AT Y = C   AND Y = D.
!
!
!                          = 0  IF THE SOLUTION IS PERIODIC IN Y, I.E.
!                               U(I,J) = U(I,N+J).
!
!                          = 1  IF THE SOLUTION IS SPECIFIED AT Y = C
!                               AND Y = D.
!
!                          = 2  IF THE SOLUTION IS SPECIFIED AT Y = C
!                               AND THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO Y IS SPECIFIED AT
!                               Y = D.
!
!                          = 3  IF THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO Y IS SPECIFIED AT
!                               Y = C AND Y = D.
!
!                          = 4  IF THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO Y IS SPECIFIED AT
!                               Y = C AND THE SOLUTION IS SPECIFIED
!                               AT Y = D.
!
!                        BDC
!                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
!                          SPECIFIES THE BOUNDARY VALUES OF THE
!                          SOLUTION AT Y = C.
!
!                          WHEN NBDCND = 1 OR 2,
!                            BDC(I) = U(X(I),C) ,        I=1,2,...,M.
!
!                          WHEN NBDCND = 3 OR 4,
!                            BDC(I) = (D/DY)U(X(I),C),   I=1,2,...,M.
!
!                          WHEN NBDCND = 0, BDC IS A DUMMY VARIABLE.
!
!                        BDD
!                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
!                          SPECIFIES THE BOUNDARY VALUES OF THE
!                          SOLUTION AT Y = D.
!
!                          WHEN NBDCND = 1 OR 4,
!                            BDD(I) = U(X(I),D) ,        I=1,2,...,M.
!
!                          WHEN NBDCND = 2 OR 3,
!                            BDD(I) = (D/DY)U(X(I),D) ,  I=1,2,...,M.
!
!                          WHEN NBDCND = 0, BDD IS A DUMMY VARIABLE.
!
!                        ELMBDA
!                          THE CONSTANT LAMBDA IN THE HELMHOLTZ
!                          EQUATION.  IF LAMBDA IS GREATER THAN 0,
!                          A SOLUTION MAY NOT EXIST. HOWEVER,
!                          HSTCRT WILL  ATTEMPT TO FIND A SOLUTION.
!
!                        F
!                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
!                          THE VALUES OF THE RIGHT SIDE OF THE
!                          HELMHOLTZ EQUATION.  FOR I=1,2,...,M
!                          AND J=1,2,...,N
!
!                            F(I,J) = F(X(I),Y(J)) .
!
!                          F MUST BE DIMENSIONED AT LEAST M X N.
!
!                        IDIMF
!                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
!                          F AS IT APPEARS IN THE PROGRAM CALLING
!                          HSTCRT.  THIS PARAMETER IS USED TO SPECIFY
!                          THE VARIABLE DIMENSION OF F.
!                          IDIMF MUST BE AT LEAST M.
!
!
! ON OUTPUT              F
!                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
!                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
!                          (X(I),Y(J)) FOR  I=1,2,...,M, J=1,2,...,N.
!
!                        PERTRB
!                          IF A COMBINATION OF PERIODIC OR DERIVATIVE
!                          BOUNDARY CONDITIONS IS SPECIFIED FOR A
!                          POISSON EQUATION (LAMBDA = 0), A SOLUTION
!                          MAY NOT EXIST.  PERTRB IS A CONSTANT,
!                          CALCULATED AND SUBTRACTED FROM F, WHICH
!                          ENSURES THAT A SOLUTION EXISTS.  HSTCRT
!                          THEN COMPUTES THIS SOLUTION, WHICH IS A
!                          LEAST SQUARES SOLUTION TO THE ORIGINAL
!                          APPROXIMATION.  THIS SOLUTION PLUS ANY
!                          CONSTANT IS ALSO A SOLUTION; HENCE, THE
!                          SOLUTION IS NOT UNIQUE.  THE VALUE OF
!                          PERTRB SHOULD BE SMALL COMPARED TO THE
!                          RIGHT SIDE F.  OTHERWISE, A SOLUTION IS
!                          OBTAINED TO AN ESSENTIALLY DIFFERENT PROBLEM.
!                          THIS COMPARISON SHOULD ALWAYS BE MADE TO
!                          INSURE THAT A MEANINGFUL SOLUTION HAS BEEN
!                          OBTAINED.
!
!                        IERROR
!                          AN ERROR FLAG THAT INDICATES INVALID INPUT
!                          PARAMETERS.  EXCEPT TO NUMBERS 0 AND  6,
!                          A SOLUTION IS NOT ATTEMPTED.
!
!                          =  0  NO ERROR
!
!                          =  1  A .GE. B
!
!                          =  2  MBDCND .LT. 0 OR MBDCND .GT. 4
!
!                          =  3  C .GE. D
!
!                          =  4  N .LE. 2
!
!                         =  5  NBDCND .LT. 0 OR NBDCND .GT. 4
!
!                         =  6  LAMBDA .GT. 0
!
!                         =  7  IDIMF .LT. M
!
!                         =  8  M .LE. 2
!
!                         SINCE THIS IS THE ONLY MEANS OF INDICATING
!                         A POSSIBLY INCORRECT CALL TO HSTCRT, THE
!                         USER SHOULD TEST IERROR AFTER THE CALL.
!
!                          = 20 If the dynamic allocation of real and
!                               complex work space required for solution
!                               fails (for example if N,M are too large
!                               for your computer)
!
!
! I/O                    NONE
!
! PRECISION              SINGLE
!
! REQUIRED LIBRARY       fish.f,comf.f,genbun.f,gnbnaux.f,poistg.f
! FILES
!
! LANGUAGE               FORTRAN 90
!
! HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN 1977.
!                        RELEASED ON NCAR'S PUBLIC SOFTWARE LIBRARIES
!                        IN JANUARY 1980.
!                        Revised in June 2004 by John Adams using
!                        Fortran 90 dynamically allocated work space.
!
! PORTABILITY            FORTRAN 90
!
! ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE-DIFFERENCE
!                        EQUATIONS, INCORPORATES BOUNDARY DATA, ADJUSTS
!                        THE RIGHT SIDE WHEN THE SYSTEM IS SINGULAR
!                        AND CALLS EITHER POISTG OR GENBUN WHICH SOLVES
!                        THE LINEAR SYSTEM OF EQUATIONS.
!
! TIMING                 FOR LARGE M AND N, THE OPERATION COUNT
!                        IS ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
!
! ACCURACY               THE SOLUTION PROCESS EMPLOYED RESULTS IN A
!                        LOSS OF NO MORE THAN FOUR SIGNIFICANT DIGITS
!                        FOR N AND M AS LARGE AS 64.  MORE DETAILED
!                        INFORMATION ABOUT ACCURACY CAN BE FOUND IN
!                        THE DOCUMENTATION FOR PACKAGE POISTG WHICH
!                        SOLVES THE FINITE DIFFERENCE EQUATIONS.
!
! REFERENCES             U. SCHUMANN AND R. SWEET,"A DIRECT METHOD
!                        FOR THE SOLUTION OF POISSON'S EQUATION WITH
!                        BOUNDARY CONDITIONS ON A STAGGERED GRID OF
!                        ARBITRARY SIZE," J. COMP. PHYS. 20(1976),
!                        PP. 171-182.
!***********************************************************************
      SUBROUTINE HSTCRT(A, B, M, MBDCND, vecBDA, vecBDB, C, D, N, NBDCND, &
                        vecBDC, vecBDD, ELMBDA, matF, IDIMF, PERTRB, IERROR)
      USE fish
      implicit none
      
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, &
                                     ONE  = 1.0D0, TWO  = 2.0D0, &
                                     FOUR = 4.0D0

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,INTENT(IN)  :: M,MBDCND,N,NBDCND,IDIMF
      INTEGER,INTENT(OUT) :: IERROR
      DOUBLE PRECISION,INTENT(IN)                       :: A,B,C,D,ELMBDA
      DOUBLE PRECISION,INTENT(OUT)                      :: PERTRB
      DOUBLE PRECISION,DIMENSION(N),INTENT(IN)          :: vecBDA,vecBDB
      DOUBLE PRECISION,DIMENSION(M),INTENT(IN)          :: vecBDC,vecBDD
      DOUBLE PRECISION,DIMENSION(IDIMF,N),INTENT(INOUT) :: matF
!-----------------------------------------------
!   Allocatable arrays
!-----------------------------------------------
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE         :: work
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: irwk, istatus

      ! JANNIS: add interface
      interface
         SUBROUTINE HSTCRTT(A, B, M, MBDCND, vecBDA, vecBDB, C, D, N, NBDCND, &
              vecBDC, vecBDD, ELMBDA, matF, IDIMF, PERTRB, IERROR, vecW, IW)
           INTEGER,INTENT(IN)  :: M,MBDCND,N,NBDCND,IDIMF,IW
           INTEGER,INTENT(OUT) :: IERROR
           DOUBLE PRECISION,INTENT(IN)                       :: A,B,C,D,ELMBDA
           DOUBLE PRECISION,INTENT(OUT)                      :: PERTRB
           DOUBLE PRECISION,DIMENSION(N),INTENT(IN)          :: vecBDA,vecBDB
           DOUBLE PRECISION,DIMENSION(M),INTENT(IN)          :: vecBDC,vecBDD
           DOUBLE PRECISION,DIMENSION(IDIMF,N),INTENT(INOUT) :: matF
           DOUBLE PRECISION,DIMENSION(IW),INTENT(INOUT)      :: vecW
         end subroutine HSTCRTT
      end interface
!-----------------------------------------------
!
!     CHECK FOR INVALID PARAMETERS.
!
      IERROR = 0

      IF (A >= B) IERROR = 1
      IF (MBDCND<0 .OR. MBDCND>4) IERROR = 2
      IF (C >= D) IERROR = 3
      IF (N <= 2) IERROR = 4
      IF (NBDCND<0 .OR. NBDCND>4) IERROR = 5
      IF (IDIMF < M) IERROR = 7
      IF (M <= 2) IERROR = 8
      IF (IERROR /= 0) RETURN 
!C!      write(*,*) 'hstcrt: vecBDC:',vecBDC(1:size(vecBDC))
!C!      write(*,*) 'hstcrt: BDD:',BDD(:)
!C!      write(*,*) 'hstcrt: ELMBDA:',ELMBDA
!C!      write(*,*) 'hstcrt: matF:',matF(1:IDIMF,:)
!C!      write(*,*) 'hstcrt: IDIMF:',IDIMF
!     compute and allocate required real work space
      CALL GEN_SPACE (N, M, IRWK)
      irwk = irwk + 3*M
      allocate(work(irwK),STAT=istatus)
!@!      write(*,*) 'HSTCRT: allocate work(irwk); irwk=',irwk
!     return if allocation failed (e.g., if n,m are too large)
      IF (istatus > 0)  THEN
         write(*,*) 'HSTCRT: error allocate work(irwk); irwk=',irwk
         RETURN 
      END IF

!     check that allocation was successful
      call hstcrtt(a,b,m,mbdcnd,vecbda,vecbdb,c,d,n,nbdcnd, &
                   vecbdc,vecbdd,elmbda,matf,idimf,pertrb,ierror, &
                   work,size(work))
!
!@!      write(*,*) 'HSTCRT: deallocate work(irwk); irwk=',irwk
      deallocate(work,STAT=istatus)
      IF (istatus > 0)  THEN
         write(*,*) 'HSTCRT: error deallocate work'
         RETURN 
      END IF
!
      END SUBROUTINE HSTCRT

      SUBROUTINE HSTCRTT(A, B, M, MBDCND, vecBDA, vecBDB, C, D, N, NBDCND, &
                         vecBDC, vecBDD, ELMBDA, matF, IDIMF, PERTRB, &
                         IERROR, vecW, IW)

!      USE genbunal
!      USE poisson

      implicit none
      
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, &
                                     ONE  = 1.0D0, TWO  = 2.0D0, &
                                     FOUR = 4.0D0

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,INTENT(IN)  :: M,MBDCND,N,NBDCND,IDIMF,IW
      INTEGER,INTENT(OUT) :: IERROR
      DOUBLE PRECISION,INTENT(IN)                       :: A,B,C,D,ELMBDA
      DOUBLE PRECISION,INTENT(OUT)                      :: PERTRB
      DOUBLE PRECISION,DIMENSION(N),INTENT(IN)          :: vecBDA,vecBDB
      DOUBLE PRECISION,DIMENSION(M),INTENT(IN)          :: vecBDC,vecBDD
      DOUBLE PRECISION,DIMENSION(IDIMF,N),INTENT(INOUT) :: matF
      DOUBLE PRECISION,DIMENSION(IW),INTENT(INOUT)      :: vecW
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER          :: NPEROD, MPEROD, NP, MP, ID2, ID3, ID4, &
                          J, IERR1, IW2
      DOUBLE PRECISION :: DELTAX,TWDELX,DELXSQ,DELTAY,TWDELY,DELYSQ,TWDYSQ,S,ST2
!-----------------------------------------------
 
      NPEROD = NBDCND
      MPEROD = 0
      IF (MBDCND > 0) MPEROD = 1
      DELTAX = (B - A)/DBLE(M)
      TWDELX = ONE/DELTAX
      DELXSQ = TWO/DELTAX**2
      DELTAY = (D - C)/DBLE(N)
      TWDELY = ONE/DELTAY
      DELYSQ = DELTAY**2
      TWDYSQ = TWO/DELYSQ
      NP = NBDCND + 1
      MP = MBDCND + 1
!
!     DEFINE THE A,B,C COEFFICIENTS IN W-ARRAY.
!
      ID2 = M
      ID3 = ID2 + M
      ID4 = ID3 + M
      S = (DELTAY/DELTAX)**2
      ST2 = TWO*S
      vecW(:M) = S
      vecW(ID2+1:M+ID2) = (-ST2) + ELMBDA*DELYSQ
      vecW(ID3+1:M+ID3) = S
!
!     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
!
      GO TO (111,102,102,104,104) MP
  102 CONTINUE
      matF(1,:N) = matF(1,:N) - vecBDA(:N)*DELXSQ
      vecW(ID2+1) = vecW(ID2+1) - vecW(1)
      GO TO 106
  104 CONTINUE
      matF(1,:N) = matF(1,:N) + vecBDA(:N)*TWDELX
      vecW(ID2+1) = vecW(ID2+1) + vecW(1)
  106 CONTINUE
      GO TO (111,107,109,109,107) MP
  107 CONTINUE
      matF(M,:N) = matF(M,:N) - vecBDB(:N)*DELXSQ
      vecW(ID3) = vecW(ID3) - vecW(1)
      GO TO 111
  109 CONTINUE
      matF(M,:N) = matF(M,:N) - vecBDB(:N)*TWDELX
      vecW(ID3) = vecW(ID3) + vecW(1)
  111 CONTINUE
      GO TO (121,112,112,114,114) NP
  112 CONTINUE
      matF(:M,1) = matF(:M,1) - vecBDC(:M)*TWDYSQ
      GO TO 116
  114 CONTINUE
      matF(:M,1) = matF(:M,1) + vecBDC(:M)*TWDELY
  116 CONTINUE
      GO TO (121,117,119,119,117) NP
  117 CONTINUE
      matF(:M,N) = matF(:M,N) - vecBDD(:M)*TWDYSQ
      GO TO 121
  119 CONTINUE
      matF(:M,N) = matF(:M,N) - vecBDD(:M)*TWDELY
  121 CONTINUE
      matF(:M,:N) = matF(:M,:N)*DELYSQ
      IF (MPEROD /= 0) THEN
         vecW(1) = ZERO
         vecW(ID4) = ZERO
      ENDIF
      PERTRB = ZERO
      IF (ELMBDA >= ZERO) THEN
         IF (ELMBDA /= ZERO) THEN
            IERROR = 6
         ELSE
            GO TO (127,133,133,127,133) MP
  127       CONTINUE
            GO TO (128,133,133,128,133) NP
!
!     FOR SINGULAR PROBLEMS MUST ADJUST DATA TO INSURE THAT A SOLUTION
!     WILL EXIST.
!
  128       CONTINUE
            S = ZERO
            DO J = 1, N
               S = S + SUM(matF(:M,J))
            END DO
            PERTRB = S/DBLE(M*N)
            matF(:M,:N) = matF(:M,:N) - PERTRB
            PERTRB = PERTRB/DELYSQ
!
!     SOLVE THE EQUATION.
!
         ENDIF
      ENDIF
  133 CONTINUE
      IERR1 = 0
!     Workarray splitted into parts of length M and IW
      IW2 = size(vecW)-ID4
      IF (NPEROD /= 0) THEN
         CALL POISTGG (NPEROD, N, MPEROD, M, vecW(1:M), vecW(ID2+1:ID2+M), &
                       vecW(ID3+1:ID3+M),IDIMF, matF, IERR1, vecW(ID4+1:),IW2)
      ELSE
         CALL GENBUNN (NPEROD, N, MPEROD, M, vecW(1:M), vecW(ID2+1:ID2+M), &
                       vecW(ID3+1:ID3+M), IDIMF, matF, IERR1, vecW(ID4+1:),IW2)
      ENDIF
!
! REVISION HISTORY---
!
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
! June      2004    Version 5.0, Fortran 90 changes
!-----------------------------------------------------------------------
      END SUBROUTINE HSTCRTT


!     file hstcyl.f
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
!     SUBROUTINE HSTCYL (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
!    +                   ELMBDA,F,IDIMF,PERTRB,IERROR)
!
! DIMENSION OF           BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N)
! ARGUMENTS
!
! LATEST REVISION        June 2004
!
! PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
!                        DIFFERENCE APPROXIMATION ON A STAGGERED
!                        GRID TO THE MODIFIED HELMHOLTZ EQUATION
!                        IN CYLINDRICAL COORDINATES. THIS EQUATION
!
!                          (1/R)(D/DR)(R(DU/DR)) + (D/DZ)(DU/DZ)
!
!                          + LAMBDA*(1/R**2)*U = F(R,Z)
!
!                        IS A TWO-DIMENSIONAL MODIFIED HELMHOLTZ
!                        EQUATION RESULTING FROM THE FOURIER TRANSFORM
!                        OF A THREE-DIMENSIONAL POISSON EQUATION.
!
! USAGE                  CALL HSTCYL (A,B,M,MBDCND,BDA,BDB,C,D,N,
!                                     NBDCND,BDC,BDD,ELMBDA,F,IDIMF,
!                                     PERTRB,IERROR)
!
! ARGUMENTS
! ON INPUT               A,B
!
!                          THE RANGE OF R, I.E. A .LE. R .LE. B.
!                          A MUST BE LESS THAN B AND A MUST BE
!                          BE NON-NEGATIVE.
!
!                        M
!                          THE NUMBER OF GRID POINTS IN THE INTERVAL
!                          (A,B).  THE GRID POINTS IN THE R-DIRECTION
!                          R-DIRECTION ARE GIVEN BY
!                          R(I) = A + (I-0.5)DR FOR I=1,2,...,M
!                          WHERE DR =(B-A)/M.
!                          M MUST BE GREATER THAN 2.
!
!                        MBDCND
!                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
!                          AT R = A AND R = B.
!
!                          = 1  IF THE SOLUTION IS SPECIFIED AT R = A
!                               (SEE NOTE BELOW) AND R = B.
!
!                          = 2  IF THE SOLUTION IS SPECIFIED AT R = A
!                               (SEE NOTE BELOW) AND THE DERIVATIVE
!                               OF THE SOLUTION WITH RESPECT TO R IS
!                               SPECIFIED AT R = B.
!
!                          = 3  IF THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO R IS SPECIFIED AT
!                               R = A (SEE NOTE BELOW) AND R = B.
!
!                          = 4  IF THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO R IS SPECIFIED AT
!                               R = A (SEE NOTE BELOW) AND THE
!                               SOLUTION IS SPECIFIED AT R = B.
!
!                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
!                               R = A = 0 AND THE SOLUTION IS
!                               SPECIFIED AT R = B.
!
!                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
!                               R = A = 0 AND THE DERIVATIVE OF THE
!                               SOLUTION WITH RESPECT TO R IS SPECIFIED
!                               AT R = B.
!
!                          NOTE:
!                          IF A = 0, DO NOT USE MBDCND = 1,2,3, OR 4,
!                          BUT INSTEAD USE MBDCND = 5 OR 6.
!                          THE RESULTING APPROXIMATION GIVES THE ONLY
!                          MEANINGFUL BOUNDARY CONDITION,
!                          I.E. DU/DR = 0.
!                          (SEE D. GREENSPAN, 'INTRODUCTORY NUMERICAL
!                          ANALYSIS OF ELLIPTIC BOUNDARY VALUE
!                          PROBLEMS,' HARPER AND ROW, 1965, CHAPTER 5.)
!
!                        BDA
!                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
!                          SPECIFIES THE BOUNDARY VALUES (IF ANY)
!                          OF THE SOLUTION AT R = A.
!
!                          WHEN MBDCND = 1 OR 2,
!                            BDA(J) = U(A,Z(J)) ,       J=1,2,...,N.
!
!                          WHEN MBDCND = 3 OR 4,
!                            BDA(J) = (D/DR)U(A,Z(J)) ,   J=1,2,...,N.
!
!                          WHEN MBDCND = 5 OR 6, BDA IS A DUMMY
!                          VARIABLE.
!
!                        BDB
!                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
!                          SPECIFIES THE BOUNDARY VALUES OF THE
!                          SOLUTION AT R = B.
!
!                          WHEN MBDCND = 1,4,OR 5,
!                            BDB(J) = U(B,Z(J)) ,        J=1,2,...,N.
!
!                          WHEN MBDCND = 2,3, OR 6,
!                            BDB(J) = (D/DR)U(B,Z(J)) ,   J=1,2,...,N.
!
!                        C,D
!                          THE RANGE OF Z, I.E. C .LE. Z .LE. D.
!                          C MUST BE LESS THAN D.
!
!                        N
!                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
!                          (C,D).  THE UNKNOWNS IN THE Z-DIRECTION
!                          ARE GIVEN BY Z(J) = C + (J-0.5)DZ,
!                          J=1,2,...,N, WHERE DZ = (D-C)/N.
!                          N MUST BE GREATER THAN 2.
!
!                        NBDCND
!                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
!                          AT Z = C  AND Z = D.
!
!                          = 0  IF THE SOLUTION IS PERIODIC IN Z, I.E.
!                               U(I,J) = U(I,N+J).
!
!                          = 1  IF THE SOLUTION IS SPECIFIED AT Z = C
!                               AND Z = D.
!
!                          = 2  IF THE SOLUTION IS SPECIFIED AT Z = C
!                               AND THE DERIVATIVE OF THE SOLUTION WITH
!                               RESPECT TO Z IS SPECIFIED AT Z = D.
!
!                          = 3  IF THE DERIVATIVE OF THE SOLUTION WITH
!                               RESPECT TO Z IS SPECIFIED AT Z = C
!                               AND Z = D.
!
!                          = 4  IF THE DERIVATIVE OF THE SOLUTION WITH
!                               RESPECT TO Z IS SPECIFIED AT Z = C AND
!                               THE SOLUTION IS SPECIFIED AT Z = D.
!
!                        BDC
!                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
!                          SPECIFIES THE BOUNDARY VALUES OF THE
!                          SOLUTION AT Z = C.
!
!                          WHEN NBDCND = 1 OR 2,
!                            BDC(I) = U(R(I),C) ,        I=1,2,...,M.
!
!                          WHEN NBDCND = 3 OR 4,
!                            BDC(I) = (D/DZ)U(R(I),C),    I=1,2,...,M.
!
!                          WHEN NBDCND = 0, BDC IS A DUMMY VARIABLE.
!
!                        BDD
!                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
!                          SPECIFIES THE BOUNDARY VALUES OF THE
!                          SOLUTION AT Z = D.
!
!                          WHEN NBDCND = 1 OR 4,
!                            BDD(I) = U(R(I),D) ,       I=1,2,...,M.
!
!                          WHEN NBDCND = 2 OR 3,
!                            BDD(I) = (D/DZ)U(R(I),D) ,   I=1,2,...,M.
!
!                          WHEN NBDCND = 0, BDD IS A DUMMY VARIABLE.
!
!                        ELMBDA
!                          THE CONSTANT LAMBDA IN THE MODIFIED
!                          HELMHOLTZ EQUATION.  IF LAMBDA IS GREATER
!                          THAN 0, A SOLUTION MAY NOT EXIST.
!                          HOWEVER, HSTCYL WILL ATTEMPT TO FIND A
!                          SOLUTION.  LAMBDA MUST BE ZERO WHEN
!                          MBDCND = 5 OR 6.
!
!                        F
!                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
!                          THE VALUES OF THE RIGHT SIDE OF THE
!                          MODIFIED HELMHOLTZ EQUATION.
!                          FOR I=1,2,...,M   AND J=1,2,...,N
!                            F(I,J) = F(R(I),Z(J)) .
!                          F MUST BE DIMENSIONED AT LEAST M X N.
!
!                        IDIMF
!                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
!                          F AS IT APPEARS IN THE PROGRAM CALLING
!                          HSTCYL.  THIS PARAMETER IS USED TO SPECIFY
!                          THE VARIABLE DIMENSION OF F.  IDIMF MUST
!                          BE AT LEAST M.
!
! ON OUTPUT
!
!                        F
!                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
!                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
!                          (R(I),Z(J)) FOR  I=1,2,...,M, J=1,2,...,N.
!
!                        PERTRB
!                          IF A COMBINATION OF PERIODIC, DERIVATIVE,
!                          OR UNSPECIFIED BOUNDARY CONDITIONS IS
!                          SPECIFIED FOR A POISSON EQUATION
!                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
!                          PERTRB IS A CONSTANT, CALCULATED AND
!                          SUBTRACTED FROM F, WHICH ENSURES THAT A
!                          SOLUTION EXISTS.  HSTCYL THEN COMPUTES
!                          THIS SOLUTION, WHICH IS A LEAST SQUARES
!                          SOLUTION TO THE ORIGINAL APPROXIMATION.
!                          THIS SOLUTION PLUS ANY CONSTANT IS ALSO
!                          A SOLUTION; HENCE, THE SOLUTION IS NOT
!                          UNIQUE.  THE VALUE OF PERTRB SHOULD BE
!                          SMALL COMPARED TO THE RIGHT SIDE F.
!                          OTHERWISE, A SOLUTION IS OBTAINED TO AN
!                          ESSENTIALLY DIFFERENT PROBLEM.
!                          THIS COMPARISON SHOULD ALWAYS BE MADE TO
!                          INSURE THAT A MEANINGFUL SOLUTION HAS BEEN
!                          OBTAINED.
!
!                        IERROR
!                          AN ERROR FLAG THAT INDICATES INVALID INPUT
!                          PARAMETERS. EXCEPT TO NUMBERS 0 AND 11,
!                          A SOLUTION IS NOT ATTEMPTED.
!
!                          =  0  NO ERROR
!
!                          =  1  A .LT. 0
!
!                          =  2  A .GE. B
!
!                          =  3  MBDCND .LT. 1 OR MBDCND .GT. 6
!
!                          =  4  C .GE. D
!
!                          =  5  N .LE. 2
!
!                          =  6  NBDCND .LT. 0 OR NBDCND .GT. 4
!
!                          =  7  A = 0 AND MBDCND = 1,2,3, OR 4
!
!                          =  8  A .GT. 0 AND MBDCND .GE. 5
!
!                          =  9  M .LE. 2
!
!                          = 10  IDIMF .LT. M
!
!                          = 11  LAMBDA .GT. 0
!
!                          = 12  A=0, MBDCND .GE. 5, ELMBDA .NE. 0
!
!                          SINCE THIS IS THE ONLY MEANS OF INDICATING
!                          A POSSIBLY INCORRECT CALL TO HSTCYL, THE
!                          USER SHOULD TEST IERROR AFTER THE CALL.
!
!                          = 20 If the dynamic allocation of real and
!                               complex work space required for solution
!                               fails (for example if N,M are too large
!                               for your computer)
! I/O                    NONE
!
! PRECISION              SINGLE
!
! REQUIRED LIBRARY       fish.f,comf.f,genbun.f,gnbnaux.f,poistg.f
! FILES
!
! LANGUAGE               FORTRAN 90
!
! HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN 1977.
!                        RELEASED ON NCAR'S PUBLIC SOFTWARE LIBRARIES
!                        IN JANUARY 1980.
!                        Revised in June 2004 by John Adams using
!                        Fortran 90 dynamically allocated work space.
!
! PORTABILITY            FORTRAN 90
!
! ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE-DIFFERENCE
!                        EQUATIONS, INCORPORATES BOUNDARY DATA, ADJUSTS
!                        THE RIGHT SIDE WHEN THE SYSTEM IS SINGULAR AND
!                        CALLS EITHER POISTG OR GENBUN WHICH SOLVES THE
!                        LINEAR SYSTEM OF EQUATIONS.
!
! TIMING                 FOR LARGE M AND N, THE OPERATION COUNT
!                        IS ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
!
! ACCURACY               THE SOLUTION PROCESS RESULTS IN A LOSS
!                        OF NO MORE THAN FOUR SIGNIFICANT DIGITS
!                        FOR N AND M AS LARGE AS 64.
!                        MORE DETAILED INFORMATION ABOUT ACCURACY
!                        CAN BE FOUND IN THE DOCUMENTATION FOR
!                        SUBROUTINE POISTG WHICH IS THE ROUTINE THAT
!                        ACTUALLY SOLVES THE FINITE DIFFERENCE
!                        EQUATIONS.
!
! REFERENCES             U. SCHUMANN AND R. SWEET, "A DIRECT METHOD FOR
!                        THE SOLUTION OF POISSON'S EQUATION WITH NEUMANN
!                        BOUNDARY CONDITIONS ON A STAGGERED GRID OF
!                        ARBITRARY SIZE," J. COMP. PHYS. 20(1976),
!                        PP. 171-182.
!***********************************************************************
      SUBROUTINE HSTCYL(A, B, M, MBDCND, vecBDA, vecBDB, C, D, N, &
                        NBDCND, vecBDC, vecBDD, ELMBDA, matF, IDIMF, &
                        PERTRB, IERROR)
      USE fish
      implicit none
      
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, &
                                     ONE  = 1.0D0, TWO  = 2.0D0, &
                                     FOUR = 4.0D0

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)   :: M, MBDCND, N, NBDCND, IDIMF
      INTEGER, INTENT(OUT)  :: IERROR
      DOUBLE PRECISION, INTENT(IN)  :: A, B, C, D, ELMBDA
      DOUBLE PRECISION, INTENT(OUT) :: PERTRB
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: vecBDA, vecBDB, vecBDC, &
                                                    vecBDD 
      DOUBLE PRECISION, DIMENSION(IDIMF,N), INTENT(INOUT) :: matF
!-----------------------------------------------
!   Allocatable arrays
!-----------------------------------------------
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE         :: work
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: irwk, istatus
!-----------------------------------------------
      ! JANNIS: add interface
      interface
         SUBROUTINE HSTCYLL(A, B, M, MBDCND, vecBDA, vecBDB, C, D, N, NBDCND, &
              vecBDC, vecBDD, ELMBDA, matF, IDIMF, PERTRB, IERROR, W)
           INTEGER, INTENT(IN)  :: M, MBDCND, N, NBDCND, IDIMF
           INTEGER, INTENT(OUT) :: IERROR
           DOUBLE PRECISION, INTENT(IN) :: A, B, C, D, ELMBDA
           DOUBLE PRECISION, INTENT(OUT) :: PERTRB
           DOUBLE PRECISION, DIMENSION(M), INTENT(IN) :: vecBDA, vecBDB,vecBDC, &
                vecBDD  
           DOUBLE PRECISION, DIMENSION(IDIMF,N), INTENT(INOUT) :: matF
           DOUBLE PRECISION, DIMENSION(:), INTENT(OUT)         :: W
         end subroutine hstcyll
      end interface
      
      IERROR = 0
      IF (A < ZERO) IERROR = 1
      IF (A >= B) IERROR = 2
      IF (MBDCND<=0 .OR. MBDCND>=7) IERROR = 3
      IF (C >= D) IERROR = 4
      IF (N <= 2) IERROR = 5
      IF (NBDCND<0 .OR. NBDCND>=5) IERROR = 6
      IF (A==ZERO .AND. MBDCND/=5 .AND. MBDCND/=6) IERROR = 7
      IF (A>ZERO .AND. MBDCND>=5) IERROR = 8
      IF (IDIMF < M) IERROR = 10
      IF (M <= 2) IERROR = 9
      IF (A==ZERO .AND. MBDCND>=5 .AND. ELMBDA/=ZERO) IERROR = 12
      IF (IERROR /= 0) RETURN 
!     allocate real work space
!     compute and allocate required real work space
      CALL GEN_SPACE (N, M, IRWK)
      IRWK = IRWK + 3*M
      allocate(work(irwK),STAT=istatus)
!@!      write(*,*) 'HSTCYL: allocate work(irwk); irwk=',irwk
!     return if allocation failed (e.g., if n,m are too large)
      IF (istatus > 0)  THEN
         write(*,*) 'HSTCYL: error allocate work(irwk); irwk=',irwk
         RETURN 
      END IF
!     check that allocation was successful
      call HSTCYLL(A, B, M, MBDCND, vecBDA, vecBDB, C, D, N, NBDCND,  &
                   vecBDC, vecBDD, ELMBDA, matF, IDIMF, PERTRB, IERROR, work)
!     release allocated work space
!@!      write(*,*) 'HSTCYL: deallocate work(irwk); irwk=',irwk
      deallocate(work,STAT=istatus)
      IF (istatus > 0)  THEN
         write(*,*) 'HSTCYL: error deallocate work'
         RETURN 
      END IF
!
!      RETURN 
      END SUBROUTINE HSTCYL
 
      SUBROUTINE HSTCYLL(A, B, M, MBDCND, vecBDA, vecBDB, C, D, N, NBDCND, &
                         vecBDC, vecBDD, ELMBDA, matF, IDIMF, PERTRB, IERROR, W)

!      USE poisson
!      USE genbunal

      implicit none
      
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, &
                                     ONE  = 1.0D0, TWO  = 2.0D0, &
                                     FOUR = 4.0D0

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: M, MBDCND, N, NBDCND, IDIMF
      INTEGER, INTENT(OUT) :: IERROR
      DOUBLE PRECISION, INTENT(IN) :: A, B, C, D, ELMBDA
      DOUBLE PRECISION, INTENT(OUT) :: PERTRB
      DOUBLE PRECISION, DIMENSION(M), INTENT(IN) :: vecBDA, vecBDB,vecBDC, &
                                                    vecBDD  
      DOUBLE PRECISION, DIMENSION(IDIMF,N), INTENT(INOUT) :: matF
      DOUBLE PRECISION, DIMENSION(:), INTENT(OUT)         :: W
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NP, IWB, IWC, IWR, I, J, K, LP, IERR1
      DOUBLE PRECISION :: DELTAR, DLRSQ, DELTHT, DLTHSQ, A1
!-----------------------------------------------
      DELTAR = (B - A)/DBLE(M)
      DLRSQ = DELTAR**2
      DELTHT = (D - C)/DBLE(N)
      DLTHSQ = DELTHT**2
      NP = NBDCND + 1
!
!     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
!
      IWB = M
      IWC = IWB + M
      IWR = IWC + M
      DO I = 1, M
         J = IWR + I
         W(J) = A + (DBLE(I) - HALF)*DELTAR
         W(I) = (A + DBLE(I - 1)*DELTAR)/(DLRSQ*W(J))
         K = IWC + I
         W(K) = (A + DBLE(I)*DELTAR)/(DLRSQ*W(J))
         K = IWB + I
         W(K) = ELMBDA/W(J)**2 - TWO/DLRSQ
      END DO
!
!     ENTER BOUNDARY DATA FOR R-BOUNDARIES.
!
      GO TO (102,102,104,104,106,106) MBDCND
  102 CONTINUE
      A1 = TWO*W(1)
      W(IWB+1) = W(IWB+1) - W(1)
      matF(1,:N) = matF(1,:N) - A1*vecBDA(:N)
      GO TO 106
  104 CONTINUE
      A1 = DELTAR*W(1)
      W(IWB+1) = W(IWB+1) + W(1)
      matF(1,:N) = matF(1,:N) + A1*vecBDA(:N)
  106 CONTINUE
      GO TO (107,109,109,107,107,109) MBDCND
  107 CONTINUE
      W(IWC) = W(IWC) - W(IWR)
      A1 = TWO*W(IWR)
      matF(M,:N) = matF(M,:N) - A1*vecBDB(:N)
      GO TO 111
  109 CONTINUE
      W(IWC) = W(IWC) + W(IWR)
      A1 = DELTAR*W(IWR)
      matF(M,:N) = matF(M,:N) - A1*vecBDB(:N)
!
!     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
!
  111 CONTINUE
      A1 = TWO/DLTHSQ
      GO TO (121,112,112,114,114) NP
  112 CONTINUE
      matF(:M,1) = matF(:M,1) - A1*vecBDC(:M)
      GO TO 116
  114 CONTINUE
      A1 = ONE/DELTHT
      matF(:M,1) = matF(:M,1) + A1*vecBDC(:M)
  116 CONTINUE
      A1 = TWO/DLTHSQ
      GO TO (121,117,119,119,117) NP
  117 CONTINUE
      matF(:M,N) = matF(:M,N) - A1*vecBDD(:M)
      GO TO 121
  119 CONTINUE
      A1 = ONE/DELTHT
      matF(:M,N) = matF(:M,N) - A1*vecBDD(:M)
  121 CONTINUE
      PERTRB = ZERO
      IF (ELMBDA >= ZERO) THEN
         IF (ELMBDA /= ZERO) THEN
            IERROR = 11
         ELSE
            GO TO (130,130,124,130,130,124) MBDCND
  124       CONTINUE
            GO TO (125,130,130,125,130) NP
  125       CONTINUE
            DO I = 1, M
               A1 = ZERO
               A1 = SUM(matF(I,:N))
               J = IWR + I
               PERTRB = PERTRB + A1*W(J)
            END DO
            PERTRB = PERTRB/(DBLE(M*N)*HALF*(A + B))
            matF(:M,:N) = matF(:M,:N) - PERTRB
         ENDIF
      ENDIF
  130 CONTINUE
      W(:M) = W(:M)*DLTHSQ
      W(IWC+1:M+IWC) = W(IWC+1:M+IWC)*DLTHSQ
      W(IWB+1:M+IWB) = W(IWB+1:M+IWB)*DLTHSQ
      matF(:M,:N) = matF(:M,:N)*DLTHSQ
      LP = NBDCND
      W(1) = ZERO
      W(IWR) = ZERO
!
!     SOLVE THE SYSTEM OF EQUATIONS.
!
      IERR1 = 0
      IF (NBDCND /= 0) THEN
         CALL POISTGG (LP, N, 1, M, W(1:M), W(IWB+1:IWB+M), W(IWC+1:IWC+M), &
                       IDIMF, matF, IERR1, W(IWR+1:), SIZE(W)-IWR)
      ELSE
         CALL GENBUNN (LP, N, 1, M, W(1:M), W(IWB+1:IWB+M), W(IWC+1:IWC+M), &
                       IDIMF, matF, IERR1, W(IWR+1:), SIZE(W)-IWR)
      ENDIF
!      RETURN 
!
! REVISION HISTORY---
!
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
! June      2004    Version 5.0, Fortran 90 changes
!-----------------------------------------------------------------------
      END SUBROUTINE HSTCYLL
!
!     file mhstcyl.f
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
!     SUBROUTINE MHSTCYL (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD,
!    +                   ELMBDA,ES,F,IDIMF,PERTRB,IERROR)
!
! DIMENSION OF           BDA(N),BDB(N),BDC(M),BDD(M),F(IDIMF,N)
! ARGUMENTS
!
! LATEST REVISION        June 2004
!
! PURPOSE                SOLVES THE STANDARD FIVE-POINT FINITE
!                        DIFFERENCE APPROXIMATION ON A STAGGERED
!                        GRID TO THE MODIFIED HELMHOLTZ EQUATION
!                        IN CYLINDRICAL COORDINATES. THIS EQUATION
!
!                          (D2U/DR2) + (D/DZ)(DU/DZ)
!
!                          + LAMBDA*(1/R**2)*U + ES * U  = F(R,Z)
!
!                        IS A TWO-DIMENSIONAL MODIFIED HELMHOLTZ
!                        EQUATION RESULTING FROM THE FOURIER TRANSFORM
!                        OF A THREE-DIMENSIONAL POISSON EQUATION.
!
! USAGE                  CALL MHSTCYL (A,B,M,MBDCND,BDA,BDB,C,D,N,
!                                     NBDCND,BDC,BDD,ELMBDA,ES,F,IDIMF,
!                                     PERTRB,IERROR)
!
! ARGUMENTS
! ON INPUT               A,B
!
!                          THE RANGE OF R, I.E. A .LE. R .LE. B.
!                          A MUST BE LESS THAN B AND A MUST BE
!                          BE NON-NEGATIVE.
!
!                        M
!                          THE NUMBER OF GRID POINTS IN THE INTERVAL
!                          (A,B).  THE GRID POINTS IN THE R-DIRECTION
!                          R-DIRECTION ARE GIVEN BY
!                          R(I) = A + (I-0.5)DR FOR I=1,2,...,M
!                          WHERE DR =(B-A)/M.
!                          M MUST BE GREATER THAN 2.
!
!                        MBDCND
!                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
!                          AT R = A AND R = B.
!
!                          = 1  IF THE SOLUTION IS SPECIFIED AT R = A
!                               (SEE NOTE BELOW) AND R = B.
!
!                          = 2  IF THE SOLUTION IS SPECIFIED AT R = A
!                               (SEE NOTE BELOW) AND THE DERIVATIVE
!                               OF THE SOLUTION WITH RESPECT TO R IS
!                               SPECIFIED AT R = B.
!
!                          = 3  IF THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO R IS SPECIFIED AT
!                               R = A (SEE NOTE BELOW) AND R = B.
!
!                          = 4  IF THE DERIVATIVE OF THE SOLUTION
!                               WITH RESPECT TO R IS SPECIFIED AT
!                               R = A (SEE NOTE BELOW) AND THE
!                               SOLUTION IS SPECIFIED AT R = B.
!
!                          = 5  IF THE SOLUTION IS UNSPECIFIED AT
!                               R = A = 0 AND THE SOLUTION IS
!                               SPECIFIED AT R = B.
!
!                          = 6  IF THE SOLUTION IS UNSPECIFIED AT
!                               R = A = 0 AND THE DERIVATIVE OF THE
!                               SOLUTION WITH RESPECT TO R IS SPECIFIED
!                               AT R = B.
!
!                          NOTE:
!                          IF A = 0, DO NOT USE MBDCND = 1,2,3, OR 4,
!                          BUT INSTEAD USE MBDCND = 5 OR 6.
!                          THE RESULTING APPROXIMATION GIVES THE ONLY
!                          MEANINGFUL BOUNDARY CONDITION,
!                          I.E. DU/DR = 0.
!                          (SEE D. GREENSPAN, 'INTRODUCTORY NUMERICAL
!                          ANALYSIS OF ELLIPTIC BOUNDARY VALUE
!                          PROBLEMS,' HARPER AND ROW, 1965, CHAPTER 5.)
!
!                        BDA
!                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
!                          SPECIFIES THE BOUNDARY VALUES (IF ANY)
!                          OF THE SOLUTION AT R = A.
!
!                          WHEN MBDCND = 1 OR 2,
!                            BDA(J) = U(A,Z(J)) ,       J=1,2,...,N.
!
!                          WHEN MBDCND = 3 OR 4,
!                            BDA(J) = (D/DR)U(A,Z(J)) ,   J=1,2,...,N.
!
!                          WHEN MBDCND = 5 OR 6, BDA IS A DUMMY
!                          VARIABLE.
!
!                        BDB
!                          A ONE-DIMENSIONAL ARRAY OF LENGTH N THAT
!                          SPECIFIES THE BOUNDARY VALUES OF THE
!                          SOLUTION AT R = B.
!
!                          WHEN MBDCND = 1,4,OR 5,
!                            BDB(J) = U(B,Z(J)) ,        J=1,2,...,N.
!
!                          WHEN MBDCND = 2,3, OR 6,
!                            BDB(J) = (D/DR)U(B,Z(J)) ,   J=1,2,...,N.
!
!                        C,D
!                          THE RANGE OF Z, I.E. C .LE. Z .LE. D.
!                          C MUST BE LESS THAN D.
!
!                        N
!                          THE NUMBER OF UNKNOWNS IN THE INTERVAL
!                          (C,D).  THE UNKNOWNS IN THE Z-DIRECTION
!                          ARE GIVEN BY Z(J) = C + (J-0.5)DZ,
!                          J=1,2,...,N, WHERE DZ = (D-C)/N.
!                          N MUST BE GREATER THAN 2.
!
!                        NBDCND
!                          INDICATES THE TYPE OF BOUNDARY CONDITIONS
!                          AT Z = C  AND Z = D.
!
!                          = 0  IF THE SOLUTION IS PERIODIC IN Z, I.E.
!                               U(I,J) = U(I,N+J).
!
!                          = 1  IF THE SOLUTION IS SPECIFIED AT Z = C
!                               AND Z = D.
!
!                          = 2  IF THE SOLUTION IS SPECIFIED AT Z = C
!                               AND THE DERIVATIVE OF THE SOLUTION WITH
!                               RESPECT TO Z IS SPECIFIED AT Z = D.
!
!                          = 3  IF THE DERIVATIVE OF THE SOLUTION WITH
!                               RESPECT TO Z IS SPECIFIED AT Z = C
!                               AND Z = D.
!
!                          = 4  IF THE DERIVATIVE OF THE SOLUTION WITH
!                               RESPECT TO Z IS SPECIFIED AT Z = C AND
!                               THE SOLUTION IS SPECIFIED AT Z = D.
!
!                        BDC
!                          A ONE DIMENSIONAL ARRAY OF LENGTH M THAT
!                          SPECIFIES THE BOUNDARY VALUES OF THE
!                          SOLUTION AT Z = C.
!
!                          WHEN NBDCND = 1 OR 2,
!                            BDC(I) = U(R(I),C) ,        I=1,2,...,M.
!
!                          WHEN NBDCND = 3 OR 4,
!                            BDC(I) = (D/DZ)U(R(I),C),    I=1,2,...,M.
!
!                          WHEN NBDCND = 0, BDC IS A DUMMY VARIABLE.
!
!                        BDD
!                          A ONE-DIMENSIONAL ARRAY OF LENGTH M THAT
!                          SPECIFIES THE BOUNDARY VALUES OF THE
!                          SOLUTION AT Z = D.
!
!                          WHEN NBDCND = 1 OR 4,
!                            BDD(I) = U(R(I),D) ,       I=1,2,...,M.
!
!                          WHEN NBDCND = 2 OR 3,
!                            BDD(I) = (D/DZ)U(R(I),D) ,   I=1,2,...,M.
!
!                          WHEN NBDCND = 0, BDD IS A DUMMY VARIABLE.
!
!                        ELMBDA
!                          THE CONSTANT LAMBDA IN THE MODIFIED
!                          HELMHOLTZ EQUATION.  IF LAMBDA IS GREATER
!                          THAN 0, A SOLUTION MAY NOT EXIST.
!                          HOWEVER, MHSTCYL WILL ATTEMPT TO FIND A
!                          SOLUTION.  LAMBDA MUST BE ZERO WHEN
!                          MBDCND = 5 OR 6.
!
!                        F
!                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES
!                          THE VALUES OF THE RIGHT SIDE OF THE
!                          MODIFIED HELMHOLTZ EQUATION.
!                          FOR I=1,2,...,M   AND J=1,2,...,N
!                            F(I,J) = F(R(I),Z(J)) .
!                          F MUST BE DIMENSIONED AT LEAST M X N.
!
!                        IDIMF
!                          THE ROW (OR FIRST) DIMENSION OF THE ARRAY
!                          F AS IT APPEARS IN THE PROGRAM CALLING
!                          MHSTCYL.  THIS PARAMETER IS USED TO SPECIFY
!                          THE VARIABLE DIMENSION OF F.  IDIMF MUST
!                          BE AT LEAST M.
!
! ON OUTPUT
!
!                        F
!                          CONTAINS THE SOLUTION U(I,J) OF THE FINITE
!                          DIFFERENCE APPROXIMATION FOR THE GRID POINT
!                          (R(I),Z(J)) FOR  I=1,2,...,M, J=1,2,...,N.
!
!                        PERTRB
!                          IF A COMBINATION OF PERIODIC, DERIVATIVE,
!                          OR UNSPECIFIED BOUNDARY CONDITIONS IS
!                          SPECIFIED FOR A POISSON EQUATION
!                          (LAMBDA = 0), A SOLUTION MAY NOT EXIST.
!                          PERTRB IS A CONSTANT, CALCULATED AND
!                          SUBTRACTED FROM F, WHICH ENSURES THAT A
!                          SOLUTION EXISTS.  MHSTCYL THEN COMPUTES
!                          THIS SOLUTION, WHICH IS A LEAST SQUARES
!                          SOLUTION TO THE ORIGINAL APPROXIMATION.
!                          THIS SOLUTION PLUS ANY CONSTANT IS ALSO
!                          A SOLUTION; HENCE, THE SOLUTION IS NOT
!                          UNIQUE.  THE VALUE OF PERTRB SHOULD BE
!                          SMALL COMPARED TO THE RIGHT SIDE F.
!                          OTHERWISE, A SOLUTION IS OBTAINED TO AN
!                          ESSENTIALLY DIFFERENT PROBLEM.
!                          THIS COMPARISON SHOULD ALWAYS BE MADE TO
!                          INSURE THAT A MEANINGFUL SOLUTION HAS BEEN
!                          OBTAINED.
!
!                        IERROR
!                          AN ERROR FLAG THAT INDICATES INVALID INPUT
!                          PARAMETERS. EXCEPT TO NUMBERS 0 AND 11,
!                          A SOLUTION IS NOT ATTEMPTED.
!
!                          =  0  NO ERROR
!
!                          =  1  A .LT. 0
!
!                          =  2  A .GE. B
!
!                          =  3  MBDCND .LT. 1 OR MBDCND .GT. 6
!
!                          =  4  C .GE. D
!
!                          =  5  N .LE. 2
!
!                          =  6  NBDCND .LT. 0 OR NBDCND .GT. 4
!
!                          =  7  A = 0 AND MBDCND = 1,2,3, OR 4
!
!                          =  8  A .GT. 0 AND MBDCND .GE. 5
!
!                          =  9  M .LE. 2
!
!                          = 10  IDIMF .LT. M
!
!                          = 11  LAMBDA .GT. 0
!
!                          = 12  A=0, MBDCND .GE. 5, ELMBDA .NE. 0
!
!                          SINCE THIS IS THE ONLY MEANS OF INDICATING
!                          A POSSIBLY INCORRECT CALL TO MHSTCYL, THE
!                          USER SHOULD TEST IERROR AFTER THE CALL.
!
!                          = 20 If the dynamic allocation of real and
!                               complex work space required for solution
!                               fails (for example if N,M are too large
!                               for your computer)
! I/O                    NONE
!
! PRECISION              SINGLE
!
! REQUIRED LIBRARY       fish.f,comf.f,genbun.f,gnbnaux.f,poistg.f
! FILES
!
! LANGUAGE               FORTRAN 90
!
! HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN 1977.
!                        RELEASED ON NCAR'S PUBLIC SOFTWARE LIBRARIES
!                        IN JANUARY 1980.
!                        Revised in June 2004 by John Adams using
!                        Fortran 90 dynamically allocated work space.
!
! PORTABILITY            FORTRAN 90
!
! ALGORITHM              THIS SUBROUTINE DEFINES THE FINITE-DIFFERENCE
!                        EQUATIONS, INCORPORATES BOUNDARY DATA, ADJUSTS
!                        THE RIGHT SIDE WHEN THE SYSTEM IS SINGULAR AND
!                        CALLS EITHER POISTG OR GENBUN WHICH SOLVES THE
!                        LINEAR SYSTEM OF EQUATIONS.
!
! TIMING                 FOR LARGE M AND N, THE OPERATION COUNT
!                        IS ROUGHLY PROPORTIONAL TO M*N*LOG2(N).
!
! ACCURACY               THE SOLUTION PROCESS RESULTS IN A LOSS
!                        OF NO MORE THAN FOUR SIGNIFICANT DIGITS
!                        FOR N AND M AS LARGE AS 64.
!                        MORE DETAILED INFORMATION ABOUT ACCURACY
!                        CAN BE FOUND IN THE DOCUMENTATION FOR
!                        SUBROUTINE POISTG WHICH IS THE ROUTINE THAT
!                        ACTUALLY SOLVES THE FINITE DIFFERENCE
!                        EQUATIONS.
!
! REFERENCES             U. SCHUMANN AND R. SWEET, "A DIRECT METHOD FOR
!                        THE SOLUTION OF POISSON'S EQUATION WITH NEUMANN
!                        BOUNDARY CONDITIONS ON A STAGGERED GRID OF
!                        ARBITRARY SIZE," J. COMP. PHYS. 20(1976),
!                        PP. 171-182.
!***********************************************************************
      SUBROUTINE MHSTCYL(A, B, M, MBDCND, vecBDA, vecBDB, C, D, &
                         N, NBDCND, vecBDC, vecBDD, ELMBDA, ES, &
                         matF, IDIMF, PERTRB, IERROR)
      USE fish
      implicit none
      
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, &
                                     ONE  = 1.0D0, TWO  = 2.0D0, &
                                     FOUR = 4.0D0

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: M, MBDCND, N, NBDCND, IDIMF 
      INTEGER, INTENT(OUT) :: IERROR
      DOUBLE PRECISION, INTENT(IN)  :: A, B, C, D, ELMBDA, ES
      DOUBLE PRECISION, INTENT(OUT) :: PERTRB
      DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: vecBDA, vecBDB
      DOUBLE PRECISION, DIMENSION(M), INTENT(IN) :: vecBDC, vecBDD
      DOUBLE PRECISION, DIMENSION(IDIMF,N), INTENT(OUT) :: matF
!-----------------------------------------------
!   Allocatable arrays
!-----------------------------------------------
      DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE         :: work
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: irwk, istatus
!-----------------------------------------------
      IERROR = 0
      IF (A < ZERO) IERROR = 1
      IF (A >= B) IERROR = 2
      IF (MBDCND<=0 .OR. MBDCND>=7) IERROR = 3
      IF (C >= D) IERROR = 4
      IF (N <= 2) IERROR = 5
      IF (NBDCND<0 .OR. NBDCND>=5) IERROR = 6
!Mod      IF (A==ZERO .AND. MBDCND/=5 .AND. MBDCND/=6) IERROR = 7
      IF (A>ZERO .AND. MBDCND>=5) IERROR = 8
      IF (IDIMF < M) IERROR = 10
      IF (M <= 2) IERROR = 9
!Mod      IF (A==ZERO .AND. MBDCND>=5 .AND. ELMBDA/=ZERO) IERROR = 12
      IF (IERROR /= 0) RETURN 
!     allocate real work space
!     compute and allocate required real work space
      CALL GEN_SPACE (N, M, IRWK)
      IRWK = IRWK + 3*M
!@!      write(*,*) 'MHSTCYL: allocate work(irwk); irwk=',irwk
      allocate(work(irwk),STAT=istatus)
!     return if allocation failed (e.g., if n,m are too large)
      IF (istatus > 0)  THEN
         write(*,*) 'HSTCRT: error allocate work(irwk); irwk=',irwk
         RETURN 
      END IF
      call MHSTCYLL(A, B, M, MBDCND, vecBDA, vecBDB, C, D, N, NBDCND, &
                    vecBDC, vecBDD, ELMBDA, ES, matF, IDIMF, PERTRB, &
                    IERROR, work, IRWK )
!     release allocated work space
!@!      write(*,*) 'MHSTCYL: deallocate work(irwk); irwk=',irwk
      deallocate(work,STAT=istatus)
      IF (istatus > 0)  THEN
         write(*,*) 'MHSTCYL: error deallocate work'
         RETURN 
      END IF
!
      END SUBROUTINE MHSTCYL
 
      SUBROUTINE MHSTCYLL(A, B, M, MBDCND, vecBDA, vecBDB, C, D, N, NBDCND, &
                          vecBDC, vecBDD, ELMBDA, ES, matF, IDIMF, &
                          PERTRB, IERROR, W, IW)
!      USE genbunal
!      USE poisson

      implicit none
      
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, &
                                     ONE  = 1.0D0, TWO  = 2.0D0, &
                                     FOUR = 4.0D0

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: M,  MBDCND, N, NBDCND, IDIMF, IW
      INTEGER, INTENT(OUT) :: IERROR
      DOUBLE PRECISION, INTENT(IN) :: A, B, C, D, ELMBDA, ES
      DOUBLE PRECISION, INTENT(OUT) :: PERTRB
      DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: vecBDA, vecBDB
      DOUBLE PRECISION, DIMENSION(M), INTENT(IN) :: vecBDC, vecBDD
      DOUBLE PRECISION, DIMENSION(IDIMF,N), INTENT(OUT) :: matF
      DOUBLE PRECISION, DIMENSION(IW), INTENT(INOUT)   :: W
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NP, IWB, IWC, IWR, I, J, K, LP, IERR1
      DOUBLE PRECISION :: DELTAR, DLRSQ, DELTHT, DLTHSQ, A1
!-----------------------------------------------
      DELTAR = (B - A)/DBLE(M)
      DLRSQ = DELTAR**2
      DELTHT = (D - C)/DBLE(N)
      DLTHSQ = DELTHT**2
      NP = NBDCND + 1
!
!     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
!
      IWB = M
      IWC = IWB + M
      IWR = IWC + M
      DO I = 1, M
         J = IWR + I
         W(J) = A + (DBLE(I) - HALF)*DELTAR
         W(I) = (W(J) - HALF * DELTAR)/(DLRSQ*W(J))
         K = IWC + I
         W(K) = (W(J) + HALF * DELTAR)/(DLRSQ*W(J))
         K = IWB + I
         W(K) = ELMBDA/W(J)**2 + ES - TWO/DLRSQ
      END DO
!
!     ENTER BOUNDARY DATA FOR R-BOUNDARIES.
!
      GO TO (102,102,104,104,106,106) MBDCND
  102 CONTINUE
      A1 = TWO*W(1)
      W(IWB+1) = W(IWB+1) - W(1)
      matF(1,:N) = matF(1,:N) - A1*vecBDA(:N)
      GO TO 106
  104 CONTINUE
      A1 = DELTAR*W(1)
      W(IWB+1) = W(IWB+1) + W(1)
      matF(1,:N) = matF(1,:N) + A1*vecBDA(:N)
  106 CONTINUE
      GO TO (107,109,109,107,107,109) MBDCND
  107 CONTINUE
      W(IWC) = W(IWC) - W(IWR)
      A1 = TWO*W(IWR)
      matF(M,:N) = matF(M,:N) - A1*vecBDB(:N)
      GO TO 111
  109 CONTINUE
      W(IWC) = W(IWC) + W(IWR)
      A1 = DELTAR*W(IWR)
      matF(M,:N) = matF(M,:N) - A1*vecBDB(:N)
!
!     ENTER BOUNDARY DATA FOR THETA-BOUNDARIES.
!
  111 CONTINUE
      A1 = TWO/DLTHSQ
      GO TO (121,112,112,114,114) NP
  112 CONTINUE
      matF(:M,1) = matF(:M,1) - A1*vecBDC(:M)
      GO TO 116
  114 CONTINUE
      A1 = ONE/DELTHT
      matF(:M,1) = matF(:M,1) + A1*vecBDC(:M)
  116 CONTINUE
      A1 = TWO/DLTHSQ
      GO TO (121,117,119,119,117) NP
  117 CONTINUE
      matF(:M,N) = matF(:M,N) - A1*vecBDD(:M)
      GO TO 121
  119 CONTINUE
      A1 = ONE/DELTHT
      matF(:M,N) = matF(:M,N) - A1*vecBDD(:M)
  121 CONTINUE
      PERTRB = ZERO
      IF (ELMBDA >= ZERO) THEN
         IF (ELMBDA /= ZERO) THEN
            IERROR = 11
         ELSE
            GO TO (130,130,124,130,130,124) MBDCND
  124       CONTINUE
            GO TO (125,130,130,125,130) NP
  125       CONTINUE
            DO I = 1, M
               A1 = ZERO
               A1 = SUM(matF(I,:N))
               J = IWR + I
               PERTRB = PERTRB + A1*W(J)
            END DO
            PERTRB = PERTRB/(DBLE(M*N)*HALF*(A + B))
            matF(:M,:N) = matF(:M,:N) - PERTRB
         ENDIF
      ENDIF
  130 CONTINUE
      W(:M) = W(:M)*DLTHSQ
      W(IWC+1:M+IWC) = W(IWC+1:M+IWC)*DLTHSQ
      W(IWB+1:M+IWB) = W(IWB+1:M+IWB)*DLTHSQ
      matF(:M,:N) = matF(:M,:N)*DLTHSQ
      LP = NBDCND
      W(1) = ZERO
      W(IWR) = ZERO
!
!     SOLVE THE SYSTEM OF EQUATIONS.
!
      IERR1 = 0
      IF (NBDCND /= 0) THEN
         CALL POISTGG (LP, N, 1, M, W, W(IWB+1:IWB+M), W(IWC+1:IWC+M), &
                       IDIMF, matF,IERR1, W(IWR+1:), size(W)-IWR)
      ELSE
         CALL GENBUNN (LP, N, 1, M, W, W(IWB+1:IWB+M), W(IWC+1:IWC+M), &
                       IDIMF, matF,IERR1, W(IWR+1:), size(W)-IWR)
      ENDIF
!
! REVISION HISTORY---
!
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
! June      2004    Version 5.0, Fortran 90 changes
! February  2009    Modified version of HSTCYL
!-----------------------------------------------------------------------
      END SUBROUTINE MHSTCYLL
