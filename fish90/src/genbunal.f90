!
!     file genbun.f
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
!     SUBROUTINE GENBUN (NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR)
!
!
! DIMENSION OF           A(M),B(M),C(M),Y(IDIMY,N)
! ARGUMENTS
!
! LATEST REVISION        JUNE 2004
!
! PURPOSE                THE NAME OF THIS PACKAGE IS A MNEMONIC FOR THE
!                        GENERALIZED BUNEMAN ALGORITHM.
!
!                        IT SOLVES THE REAL LINEAR SYSTEM OF EQUATIONS
!
!                        A(I)*X(I-1,J) + B(I)*X(I,J) + C(I)*X(I+1,J)
!                        + X(I,J-1) - 2.*X(I,J) + X(I,J+1) = Y(I,J)
!
!                        FOR I = 1,2,...,M  AND  J = 1,2,...,N.
!
!                        INDICES I+1 AND I-1 ARE EVALUATED MODULO M,
!                        I.E., X(0,J) = X(M,J) AND X(M+1,J) = X(1,J),
!                        AND X(I,0) MAY EQUAL 0, X(I,2), OR X(I,N),
!                        AND X(I,N+1) MAY EQUAL 0, X(I,N-1), OR X(I,1)
!                        DEPENDING ON AN INPUT PARAMETER.
!
! USAGE                  CALL GENBUN (NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,
!                                     IERROR)
!
! ARGUMENTS
!
! ON INPUT               NPEROD
!
!                          INDICATES THE VALUES THAT X(I,0) AND
!                          X(I,N+1) ARE ASSUMED TO HAVE.
!
!                          = 0  IF X(I,0) = X(I,N) AND X(I,N+1) =
!                               X(I,1).
!                          = 1  IF X(I,0) = X(I,N+1) = 0  .
!                          = 2  IF X(I,0) = 0 AND X(I,N+1) = X(I,N-1).
!                          = 3  IF X(I,0) = X(I,2) AND X(I,N+1) =
!                               X(I,N-1).
!                          = 4  IF X(I,0) = X(I,2) AND X(I,N+1) = 0.
!
!                        N
!                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.
!                          N MUST BE GREATER THAN 2.
!
!                        MPEROD
!                          = 0 IF A(1) AND C(M) ARE NOT ZERO
!                          = 1 IF A(1) = C(M) = 0
!
!                        M
!                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.
!                          N MUST BE GREATER THAN 2.
!
!                        A,B,C
!                          ONE-DIMENSIONAL ARRAYS OF LENGTH M THAT
!                          SPECIFY THE COEFFICIENTS IN THE LINEAR
!                          EQUATIONS GIVEN ABOVE.  IF MPEROD = 0
!                          THE ARRAY ELEMENTS MUST NOT DEPEND UPON
!                          THE INDEX I, BUT MUST BE CONSTANT.
!                          SPECIFICALLY, THE SUBROUTINE CHECKS THE
!                          FOLLOWING CONDITION .
!
!                            A(I) = C(1)
!                            C(I) = C(1)
!                            B(I) = B(1)
!
!                          FOR I=1,2,...,M.
!
!                        IDIMY
!                          THE ROW (OR FIRST) DIMENSION OF THE
!                          TWO-DIMENSIONAL ARRAY Y AS IT APPEARS
!                          IN THE PROGRAM CALLING GENBUN.
!                          THIS PARAMETER IS USED TO SPECIFY THE
!                          VARIABLE DIMENSION OF Y.
!                          IDIMY MUST BE AT LEAST M.
!
!                        Y
!                          A TWO-DIMENSIONAL COMPLEX ARRAY THAT
!                          SPECIFIES THE VALUES OF THE RIGHT SIDE
!                          OF THE LINEAR SYSTEM OF EQUATIONS GIVEN
!                          ABOVE.
!                          Y MUST BE DIMENSIONED AT LEAST M*N.
!
!
!  ON OUTPUT             Y
!
!                          CONTAINS THE SOLUTION X.
!
!                        IERROR
!                          AN ERROR FLAG WHICH INDICATES INVALID
!                          INPUT PARAMETERS  EXCEPT FOR NUMBER
!                          ZERO, A SOLUTION IS NOT ATTEMPTED.
!
!                          = 0  NO ERROR.
!                          = 1  M .LE. 2  .
!                          = 2  N .LE. 2
!                          = 3  IDIMY .LT. M
!                          = 4  NPEROD .LT. 0 OR NPEROD .GT. 4
!                          = 5  MPEROD .LT. 0 OR MPEROD .GT. 1
!                          = 6  A(I) .NE. C(1) OR C(I) .NE. C(1) OR
!                               B(I) .NE. B(1) FOR
!                               SOME I=1,2,...,M.
!                          = 7  A(1) .NE. 0 OR C(M) .NE. 0 AND
!                                 MPEROD = 1
!                          = 20 If the dynamic allocation of real and
!                               complex work space required for solution
!                               fails (for example if N,M are too large
!                               for your computer)
!
!
! SPECIAL CONDITONS      NONE
!
! I/O                    NONE
!
! PRECISION              SINGLE
!
! REQUIRED FILES         comf.f,gnbnaux.f,fish.f
! FILES
!
! LANGUAGE               FORTRAN 90
!
! HISTORY                WRITTEN IN 1979 BY ROLAND SWEET OF NCAR'S
!                        SCIENTIFIC COMPUTING DIVISION.  MADE AVAILABLE
!                        ON NCAR'S PUBLIC LIBRARIES IN JANUARY, 1980.
!                        Revised in June 2004 by John Adams using
!                        Fortran 90 dynamically allocated work space.
!
! ALGORITHM              THE LINEAR SYSTEM IS SOLVED BY A CYCLIC
!                        REDUCTION ALGORITHM DESCRIBED IN THE
!                        REFERENCE.
!
! PORTABILITY            FORTRAN 90 --
!                        THE MACHINE DEPENDENT CONSTANT PI IS
!                        DEFINED IN FUNCTION PIMACH.
!
! REFERENCES             SWEET, R., "A CYCLIC REDUCTION ALGORITHM FOR
!                        SOLVING BLOCK TRIDIAGONAL SYSTEMS OF ARBITRARY
!                        DIMENSIONS," SIAM J. ON NUMER. ANAL., 14 (1977)
!                        PP. 706-720.
!
! ACCURACY               THIS TEST WAS PERFORMED ON a platform with
!                        64 bit floating point arithmetic.
!                        A UNIFORM RANDOM NUMBER GENERATOR WAS USED
!                        TO CREATE A SOLUTION ARRAY X FOR THE SYSTEM
!                        GIVEN IN THE 'PURPOSE' DESCRIPTION ABOVE
!                        WITH
!                          A(I) = C(I) = -0.5*B(I) = 1, I=1,2,...,M
!
!                        AND, WHEN MPEROD = 1
!
!                          A(1) = C(M) = 0
!                          A(M) = C(1) = 2.
!
!                        THE SOLUTION X WAS SUBSTITUTED INTO THE
!                        GIVEN SYSTEM  AND, USING DOUBLE PRECISION
!                        A RIGHT SIDE Y WAS COMPUTED.
!                        USING THIS ARRAY Y, SUBROUTINE GENBUN
!                        WAS CALLED TO PRODUCE APPROXIMATE
!                        SOLUTION Z.  THEN RELATIVE ERROR
!                          E = MAX(ABS(Z(I,J)-X(I,J)))/
!                              MAX(ABS(X(I,J)))
!                        WAS COMPUTED, WHERE THE TWO MAXIMA ARE TAKEN
!                        OVER I=1,2,...,M AND J=1,...,N.
!
!                        THE VALUE OF E IS GIVEN IN THE TABLE
!                        BELOW FOR SOME TYPICAL VALUES OF M AND N.
!
!                   M (=N)    MPEROD    NPEROD        E
!                   ------    ------    ------      ------
!
!                     31        0         0         6.E-14
!                     31        1         1         4.E-13
!                     31        1         3         3.E-13
!                     32        0         0         9.E-14
!                     32        1         1         3.E-13
!                     32        1         3         1.E-13
!                     33        0         0         9.E-14
!                     33        1         1         4.E-13
!                     33        1         3         1.E-13
!                     63        0         0         1.E-13
!                     63        1         1         1.E-12
!                     63        1         3         2.E-13
!                     64        0         0         1.E-13
!                     64        1         1         1.E-12
!                     64        1         3         6.E-13
!                     65        0         0         2.E-13
!                     65        1         1         1.E-12
!                     65        1         3         4.E-13
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE GENBUN(NPEROD, N, MPEROD, M, A, B, C, IDIMY, Y, IERROR)
      USE fish
      implicit none
      
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, &
                                     ONE  = 1.0D0, TWO  = 2.0D0, &
                                     FOUR = 4.0D0

      TYPE(fishworkspace) :: w
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: NPEROD, N, MPEROD, M, IDIMY
      INTEGER, INTENT(OUT) :: IERROR
      DOUBLE PRECISION, DIMENSION(:)  :: A, B, C
      DOUBLE PRECISION, INTENT(INOUT) :: Y(IDIMY,*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IRWK, ICWK
!-----------------------------------------------

      write(*,*) 'fish90: genbun'
      IERROR = 0
!     check input arguments
      IF (M <= 2) then
         ierror = 1
         return
      end if
      IF (N <= 2) then
         ierror = 2
         return
      end if
      IF (IDIMY < M) then
         ierror = 3
         return
      end if
      IF (NPEROD<0 .OR. NPEROD>4) then
         ierror = 4
         return
      end if
      IF (MPEROD<0 .OR. MPEROD>1) then
         ierror = 5
         return
      end if
!     compute and allocate real work space for genbun
      CALL GEN_SPACE (N, M, IRWK)
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     return if allocation failed (e.g., if n,m are too large)
      IF (IERROR == 20)  THEN
         write(*,*) 'error call ALLOCATFISH'
         RETURN 
      END IF
      call genbunn(NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR,w%rew,IRWK)
!     release allocated work space
      CALL FISHFIN (W,IERROR)
      IF (IERROR == 20)  THEN
         write(*,*) 'error call FISHFIN'
         RETURN 
      END IF
!
      END SUBROUTINE GENBUN

 
      SUBROUTINE GENBUNN(NPEROD,N,MPEROD,M,vecA,vecB,vecC,IDIMY,matY, &
                         IERROR,vecW,IW)
      implicit none
      
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, &
                                     ONE  = 1.0D0, TWO  = 2.0D0, &
                                     FOUR = 4.0D0

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)    :: NPEROD, N, MPEROD, M, IDIMY,IW
      INTEGER, INTENT(INOUT) :: IERROR
      DOUBLE PRECISION, DIMENSION(M), INTENT(IN)   :: vecA, vecB, vecC
      DOUBLE PRECISION, DIMENSION(IW), INTENT(OUT) :: vecW
      DOUBLE PRECISION, DIMENSION(IDIMY,N),INTENT(INOUT) :: matY
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, MP1, IWBA, IWBB, IWBC, IWB2, IWB3, IWW1, IWW2, IWW3, &
                 IWD, IWTCOS, IWP, IW2,  K, J, MP, NP, &
                 IPSTOR, IREV, MH, MHM1, MODD, NBY2, MSKIP
      DOUBLE PRECISION :: A1
!-----------------------------------------------
      IF (MPEROD /= 1) THEN
         DO I = 2, M
            IF (vecA(I) /= vecC(1)) GO TO 103
            IF (vecC(I) /= vecC(1)) GO TO 103
            IF (vecB(I) /= vecB(1)) GO TO 103
         END DO
         GO TO 104
      ENDIF
      IF (vecA(1)/=ZERO .OR. vecC(M)/=ZERO) IERROR = 7
      GO TO 104
  103 CONTINUE
      IERROR = 6
  104 CONTINUE
      IF (IERROR /= 0) RETURN 
      MP1 = M + 1
      IWBA = MP1
      IWBB = IWBA + M
      IWBC = IWBB + M
      IWB2 = IWBC + M
      IWB3 = IWB2 + M
      IWW1 = IWB3 + M
      IWW2 = IWW1 + M
      IWW3 = IWW2 + M
      IWD = IWW3 + M
      IWTCOS = IWD + M
      IWP = IWTCOS + 4*N
      vecW(IWBA:M-1+IWBA) = -vecA(:M)
      vecW(IWBC:M-1+IWBC) = -vecC(:M)
      vecW(IWBB:M-1+IWBB) = TWO - vecB(:M)
      matY(:M,:N) = -matY(:M,:N)
      MP = MPEROD + 1
      NP = NPEROD + 1
      GO TO (114,107) MP
  107 CONTINUE
      GO TO (108,109,110,111,123) NP
  108 CONTINUE
      IW2 = IW - IWP + 1
      CALL POISP2 (M, N, vecW(IWBA:IWBA+M-1), vecW(IWBB:IWBB+M-1), &
                   vecW(IWBC:IWBC+M-1), matY, IDIMY, &
                   vecW(1:M), vecW(IWB2:IWB2+M-1), vecW(IWB3:IWB3+M-1), &
                   vecW(IWW1:IWW1+M-1), vecW(IWW2:IWW2+M-1), &
                   vecW(IWW3:IWW3+M-1), vecW(IWD:IWD+M-1), &
                   vecW(IWTCOS:IWTCOS+4*N-1), vecW(IWP:),IW2)
      GO TO 112
  109 CONTINUE
      IW2 = IW - IWP + 1
      CALL POISD2 (M, N, 1, vecW(IWBA:IWBA+M-1), vecW(IWBB:IWBB+M-1),  &
                   vecW(IWBC::IWBC+M-1), matY, IDIMY, &
                   vecW(1:M), vecW(IWW1:IWW1+M-1), vecW(IWD:IWD+M-1), &
                   vecW(IWTCOS:IWTCOS+4*N-1), vecW(IWP:), IW2)
      GO TO 112
  110 CONTINUE
      IW2 = IW - IWP + 1
      CALL POISN2 (M, N, 1, 2, vecW(IWBA:IWBA+M-1), vecW(IWBB:IWBB+M-1), &
                   vecW(IWBC:IWBC+M-1), matY, IDIMY, vecW(1:M),  &
                   vecW(IWB2:IWB2+M-1), vecW(IWB3:IWB3+M-1), &
                   vecW(IWW1:IWW1+M-1), vecW(IWW2:IWW2+M-1), &
                   vecW(IWW3::IWW3+M-1), vecW(IWD:IWD+M-1), &
                   vecW(IWTCOS:IWTCOS+4*N-1),vecW(IWP:),IW2)
      GO TO 112
  111 CONTINUE
      IW2 = IW - IWP + 1
      CALL POISN2 (M, N, 1, 1, vecW(IWBA:IWBA+M-1), vecW(IWBB:IWBB+M-1), &
                   vecW(IWBC::IWBC+M-1), matY, IDIMY, vecW(1:M), &
                   vecW(IWB2:IWW1+M-1), vecW(IWB3:IWB3+M-1), &
                   vecW(IWW1:IWW1+M-1), vecW(IWW2:IWW2+M-1), &
                   vecW(IWW3::IWW3+M-1), vecW(IWD:IWD+M-1), &
                   vecW(IWTCOS:IWTCOS+4*N-1),vecW(IWP:),IW2)
  112 CONTINUE
      IPSTOR = vecW(IWW1)
      IREV = 2
      IF (NPEROD == 4) GO TO 124
  113 CONTINUE
      GO TO (127,133) MP
  114 CONTINUE
      MH = (M + 1)/2
      MHM1 = MH - 1
      MODD = 1
      IF (MH*2 == M) MODD = 2
      DO J = 1, N
         vecW(:MHM1) = matY(MH-1:MH-MHM1:(-1),J) - matY(MH+1:MHM1+MH,J)
         vecW(MH+1:MHM1+MH) = matY(MH-1:MH-MHM1:(-1),J) + matY(MH+1:MHM1+MH,J)
         vecW(MH) = TWO*matY(MH,J)
         GO TO (117,116) MODD
  116    CONTINUE
         vecW(M) = TWO*matY(M,J)
  117    CONTINUE
         matY(:M,J) = vecW(:M)
      END DO
      K = IWBC + MHM1 - 1
      I = IWBA + MHM1
      vecW(K) = ZERO
      vecW(I) = ZERO
      vecW(K+1) = TWO*vecW(K+1)
      SELECT CASE (MODD) 
      CASE DEFAULT
         K = IWBB + MHM1 - 1
         vecW(K) = vecW(K) - vecW(I-1)
         vecW(IWBC-1) = vecW(IWBC-1) + vecW(IWBB-1)
      CASE (2) 
         vecW(IWBB-1) = vecW(K+1)
      END SELECT
      GO TO 107
!
!     REVERSE COLUMNS WHEN NPEROD = 4.
!
  123 CONTINUE
      IREV = 1
      NBY2 = N/2
  124 CONTINUE
      DO J = 1, NBY2
         MSKIP = N + 1 - J
         DO I = 1, M
            A1 = matY(I,J)
            matY(I,J) = matY(I,MSKIP)
            matY(I,MSKIP) = A1
         END DO
      END DO
      GO TO (110,113) IREV
  127 CONTINUE
      DO J = 1, N
         vecW(MH-1:MH-MHM1:(-1)) = HALF*(matY(MH+1:MHM1+MH,J)+matY(:MHM1,J))
         vecW(MH+1:MHM1+MH) = HALF*(matY(MH+1:MHM1+MH,J)-matY(:MHM1,J))
         vecW(MH) = HALF*matY(MH,J)
         GO TO (130,129) MODD
  129    CONTINUE
         vecW(M) = HALF*matY(M,J)
  130    CONTINUE
         matY(:M,J) = vecW(:M)
      END DO
  133 CONTINUE
      vecW(1) = IPSTOR + IWP - 1
      RETURN 
      END SUBROUTINE GENBUNN


      SUBROUTINE POISD2(MR,NR,ISTAG,vecBA,vecBB,vecBC,matY,IDIMY, &
                        vecB,vecW,vecD,TCOS,vecP,IDIMP)
      implicit none
      
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, &
                                     ONE  = 1.0D0, TWO  = 2.0D0, &
                                     FOUR = 4.0D0

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: MR, NR, ISTAG, IDIMY,IDIMP
      DOUBLE PRECISION, DIMENSION(MR), INTENT(IN)         :: vecBA, vecBB, &
                                                             vecBC
      DOUBLE PRECISION, DIMENSION(IDIMY,NR),INTENT(INOUT) :: matY
      DOUBLE PRECISION, DIMENSION(MR), INTENT(INOUT)      :: vecB, vecD, vecW
      DOUBLE PRECISION, DIMENSION(4*NR), INTENT(INOUT)    :: TCOS
      DOUBLE PRECISION, DIMENSION(IDIMP), INTENT(INOUT)   :: vecP
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: M, N, JSH, IP, IPSTOR, KR, IRREG, JSTSAV, I, LR, NUN, &
                 JST, JSP, L, NODD, J, JM1, JP1, JM2, JP2, JM3, JP3, &
                 NODDPR, KRPI, IDEG, JDEG
      DOUBLE PRECISION :: FI, T
!-----------------------------------------------
!
!     SUBROUTINE TO SOLVE POISSON'S EQUATION FOR DIRICHLET BOUNDARY
!     CONDITIONS.
!
!     ISTAG = 1 IF THE LAST DIAGONAL BLOCK IS THE MATRIX A.
!     ISTAG = 2 IF THE LAST DIAGONAL BLOCK IS THE MATRIX A+I.
!
      M = MR
      N = NR
      JSH = 0
      FI = ONE/FLOAT(ISTAG)
      IP = -M
      IPSTOR = 0
      SELECT CASE (ISTAG) 
      CASE DEFAULT
         KR = 0
         IRREG = 1
         IF (N > 1) GO TO 106
         TCOS(1) = ZERO
      CASE (2) 
         KR = 1
         JSTSAV = 1
         IRREG = 2
         IF (N > 1) GO TO 106
         TCOS(1) = -ONE
      END SELECT
      vecB(:M) = matY(:M,1)
      CALL TRIX (1, 0, M, vecBA, vecBB, vecBC, vecB, TCOS, 4*NR, vecD, vecW)
      matY(:M,1) = vecB(:M)
      GO TO 183
  106 CONTINUE
      LR = 0
      vecP(:M) = ZERO
      NUN = N
      JST = 1
      JSP = N
!
!     IRREG = 1 WHEN NO IRREGULARITIES HAVE OCCURRED, OTHERWISE IT IS 2.
!
  108 CONTINUE
      L = 2*JST
      NODD = 2 - 2*((NUN + 1)/2) + NUN
!
!     NODD = 1 WHEN NUN IS ODD, OTHERWISE IT IS 2.
!
      SELECT CASE (NODD) 
      CASE DEFAULT
         JSP = JSP - L
      CASE (1) 
         JSP = JSP - JST
         IF (IRREG /= 1) JSP = JSP - L
      END SELECT
      CALL COSGEN (JST, 1, HALF, ZERO, TCOS, 4*NR)
      IF (L <= JSP) THEN
         DO J = L, JSP, L
            JM1 = J - JSH
            JP1 = J + JSH
            JM2 = J - JST
            JP2 = J + JST
            JM3 = JM2 - JSH
            JP3 = JP2 + JSH
            IF (JST == 1) THEN
               vecB(:M) = TWO*matY(:M,J)
               matY(:M,J) = matY(:M,JM2) + matY(:M,JP2)
            ELSE
               DO I = 1, M
                  T = matY(I,J) - matY(I,JM1) - matY(I,JP1) + matY(I,JM2) + &
                      matY(I,JP2)
                  vecB(I) = T + matY(I,J) - matY(I,JM3) - matY(I,JP3)
                  matY(I,J) = T
               END DO
            ENDIF
            CALL TRIX (JST, 0, M, vecBA, vecBB, vecBC, vecB, TCOS, 4*NR, vecD, vecW)
            matY(:M,J) = matY(:M,J) + vecB(:M)
         END DO
      ENDIF
!
!     REDUCTION FOR LAST UNKNOWN
!
      SELECT CASE (NODD) 
      CASE DEFAULT
         GO TO (152,120) IRREG
!
!     ODD NUMBER OF UNKNOWNS
!
  120    CONTINUE
         JSP = JSP + L
         J = JSP
         JM1 = J - JSH
         JP1 = J + JSH
         JM2 = J - JST
         JP2 = J + JST
         JM3 = JM2 - JSH
         GO TO (123,121) ISTAG
  121    CONTINUE
         IF (JST /= 1) GO TO 123
         vecB(:M) = matY(:M,J)
         matY(:M,J) = ZERO
         GO TO 130
  123    CONTINUE
         SELECT CASE (NODDPR) 
         CASE DEFAULT
            vecB(:M) = HALF*(matY(:M,JM2)-matY(:M,JM1)-matY(:M,JM3)) + &
                       vecP(IP+1:M+IP)+ matY(:M,J)
         CASE (2) 
            vecB(:M) = HALF*(matY(:M,JM2)-matY(:M,JM1)-matY(:M,JM3)) +  &
                       matY(:M,JP2) - matY(:M,JP1) + matY(:M,J)
         END SELECT
         matY(:M,J) = HALF*(matY(:M,J)-matY(:M,JM1)-matY(:M,JP1))
  130    CONTINUE
         CALL TRIX (JST, 0, M, vecBA, vecBB, vecBC, vecB, TCOS, 4*NR, vecD, vecW)
         IP = IP + M
         IPSTOR = MAX0(IPSTOR,IP + M)
         vecP(IP+1:M+IP) = matY(:M,J) + vecB(:M)
         vecB(:M) = matY(:M,JP2) + vecP(IP+1:M+IP)
         IF (LR == 0) THEN
            DO I = 1, JST
               KRPI = KR + I
               TCOS(KRPI) = TCOS(I)
            END DO
         ELSE
            CALL COSGEN (LR, JSTSAV, ZERO, FI, TCOS(JST+1), 4*NR-JST)
            CALL MERGE (TCOS, 0, JST, JST, LR, KR, 4*NR)
         ENDIF
         CALL COSGEN (KR, JSTSAV, ZERO, FI, TCOS, 4*NR)
         CALL TRIX (KR, KR, M, vecBA, vecBB, vecBC, vecB, TCOS, 4*NR, vecD, vecW)
         matY(:M,J) = matY(:M,JM2) + vecB(:M) + vecP(IP+1:M+IP)
         LR = KR
         KR = KR + L
!
!     EVEN NUMBER OF UNKNOWNS
!
      CASE (2) 
         JSP = JSP + L
         J = JSP
         JM1 = J - JSH
         JP1 = J + JSH
         JM2 = J - JST
         JP2 = J + JST
         JM3 = JM2 - JSH
         SELECT CASE (IRREG) 
         CASE DEFAULT
            JSTSAV = JST
            IDEG = JST
            KR = L
         CASE (2) 
            CALL COSGEN (KR, JSTSAV, ZERO, FI, TCOS, 4*NR)
            CALL COSGEN (LR, JSTSAV, ZERO, FI, TCOS(KR+1), 4*NR-KR)
            IDEG = KR
            KR = KR + JST
         END SELECT
         IF (JST == 1) THEN
            IRREG = 2
            vecB(:M) = matY(:M,J)
            matY(:M,J) = matY(:M,JM2)
         ELSE
            vecB(:M) = matY(:M,J) + HALF*(matY(:M,JM2)-matY(:M,JM1)- &
                       matY(:M,JM3))
            SELECT CASE (IRREG) 
            CASE DEFAULT
               matY(:M,J) = matY(:M,JM2) + HALF*(matY(:M,J)- &
                            matY(:M,JM1)-matY(:M,JP1))
               IRREG = 2
            CASE (2) 
               SELECT CASE (NODDPR) 
               CASE DEFAULT
                  matY(:M,J) = matY(:M,JM2) + vecP(IP+1:M+IP)
                  IP = IP - M
               CASE (2) 
                  matY(:M,J) = matY(:M,JM2) + matY(:M,J) - matY(:M,JM1)
               END SELECT
            END SELECT
         ENDIF
         CALL TRIX (IDEG, LR, M, vecBA, vecBB, vecBC, vecB, TCOS, 4*NR, vecD, vecW)
         matY(:M,J) = matY(:M,J) + vecB(:M)
      END SELECT
  152 CONTINUE
      NUN = NUN/2
      NODDPR = NODD
      JSH = JST
      JST = 2*JST
      IF (NUN >= 2) GO TO 108
!
!     START SOLUTION.
!
      J = JSP
      vecB(:M) = matY(:M,J)
      SELECT CASE (IRREG) 
      CASE DEFAULT
         CALL COSGEN (JST, 1, HALF, ZERO, TCOS, 4*NR)
         IDEG = JST
      CASE (2) 
         KR = LR + JST
         CALL COSGEN (KR, JSTSAV, ZERO, FI, TCOS, 4*NR)
         CALL COSGEN (LR, JSTSAV, ZERO, FI, TCOS(KR+1), 4*NR-KR)
         IDEG = KR
      END SELECT
      CALL TRIX (IDEG, LR, M, vecBA, vecBB, vecBC, vecB, TCOS, 4*NR, vecD, vecW)
      JM1 = J - JSH
      JP1 = J + JSH
      SELECT CASE (IRREG) 
      CASE DEFAULT
         matY(:M,J) = HALF*(matY(:M,J)-matY(:M,JM1)-matY(:M,JP1)) + vecB(:M)
      CASE (2) 
         SELECT CASE (NODDPR) 
         CASE DEFAULT
            matY(:M,J) = vecP(IP+1:M+IP) + vecB(:M)
            IP = IP - M
         CASE (2) 
            matY(:M,J) = matY(:M,J) - matY(:M,JM1) + vecB(:M)
         END SELECT
      END SELECT
  164 CONTINUE
      JST = JST/2
      JSH = JST/2
      NUN = 2*NUN
      IF (NUN > N) GO TO 183
      DO J = JST, N, L
         JM1 = J - JSH
         JP1 = J + JSH
         JM2 = J - JST
         JP2 = J + JST
         IF (J <= JST) THEN
            vecB(:M) = matY(:M,J) + matY(:M,JP2)
         ELSE
            IF (JP2 <= N) GO TO 168
            vecB(:M) = matY(:M,J) + matY(:M,JM2)
            IF (JST < JSTSAV) IRREG = 1
            GO TO (170,171) IRREG
  168       CONTINUE
            vecB(:M) = matY(:M,J) + matY(:M,JM2) + matY(:M,JP2)
         ENDIF
  170    CONTINUE
         CALL COSGEN (JST, 1, HALF, ZERO, TCOS, 4*NR)
         IDEG = JST
         JDEG = 0
         GO TO 172
  171    CONTINUE
         IF (J + L > N) LR = LR - JST
         KR = JST + LR
         CALL COSGEN (KR, JSTSAV, ZERO, FI, TCOS, 4*NR)
         CALL COSGEN (LR, JSTSAV, ZERO, FI, TCOS(KR+1), 4*NR-KR)
         IDEG = KR
         JDEG = LR
  172    CONTINUE
         CALL TRIX (IDEG, JDEG, M, vecBA, vecBB, vecBC, vecB, TCOS, 4*NR, vecD, vecW)
         IF (JST <= 1) THEN
            matY(:M,J) = vecB(:M)
         ELSE
            IF (JP2 > N) GO TO 177
  175       CONTINUE
            matY(:M,J) = HALF*(matY(:M,J)-matY(:M,JM1)-matY(:M,JP1)) + vecB(:M)
            CYCLE 
  177       CONTINUE
            GO TO (175,178) IRREG
  178       CONTINUE
            IF (J + JSH <= N) THEN
               matY(:M,J) = vecB(:M) + vecP(IP+1:M+IP)
               IP = IP - M
            ELSE
               matY(:M,J) = vecB(:M) + matY(:M,J) - matY(:M,JM1)
            ENDIF
         ENDIF
      END DO
      L = L/2
      GO TO 164
  183 CONTINUE
      vecW(1) = IPSTOR
      RETURN 
      END SUBROUTINE POISD2


      SUBROUTINE POISN2(M, N, ISTAG, MIXBND, vecA, vecBB, vecC, &
                        matQ, IDIMQ, vecB, vecB2,vecB3, &
                        vecW, vecW2, vecW3, vecD, TCOS, vecP, IDIMP)
      implicit none
      
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, &
                                     ONE  = 1.0D0, TWO  = 2.0D0, &
                                     FOUR = 4.0D0

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: M, N, ISTAG, MIXBND, IDIMQ,IDIMP
      DOUBLE PRECISION, DIMENSION(M), INTENT(IN)   :: vecA, vecBB, vecC
      DOUBLE PRECISION, DIMENSION(IDIMQ,N), INTENT(INOUT) :: matQ
      DOUBLE PRECISION, DIMENSION(M),INTENT(INOUT) :: vecB, vecB2, vecB3,  &
                                                      vecD, vecW, vecW2, vecW3
      DOUBLE PRECISION, DIMENSION(IDIMP),INTENT(INOUT) :: vecP
      DOUBLE PRECISION, DIMENSION(4*N),INTENT(INOUT)   :: TCOS
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(4) :: K
      INTEGER :: K1, K2, K3, K4, MR, IP, IPSTOR, I2R, JR, NR, NLAST, &
                 KR, LR, I, NROD, JSTART, JSTOP, I2RBY2, &
                 J, JP1, JP2, JP3, JM1,JM2, JM3, NRODPR, II, I1, I2, &
                 JR2, NLASTP, JSTEP
      DOUBLE PRECISION :: FISTAG, FNUM, FDEN, FI, T
!-----------------------------------------------
!
!     SUBROUTINE TO SOLVE POISSON'S EQUATION WITH NEUMANN BOUNDARY
!     CONDITIONS.
!
!     ISTAG = 1 IF THE LAST DIAGONAL BLOCK IS A.
!     ISTAG = 2 IF THE LAST DIAGONAL BLOCK IS A-I.
!     MIXBND = 1 IF HAVE NEUMANN BOUNDARY CONDITIONS AT BOTH BOUNDARIES.
!     MIXBND = 2 IF HAVE NEUMANN BOUNDARY CONDITIONS AT BOTTOM AND
!     DIRICHLET CONDITION AT TOP.  (FOR THIS CASE, MUST HAVE ISTAG = 1.)
!
      EQUIVALENCE (K(1), K1), (K(2), K2), (K(3), K3), (K(4), K4)
      FISTAG = 3 - ISTAG
      FNUM = ONE/FLOAT(ISTAG)
      FDEN = HALF*FLOAT(ISTAG - 1)
      MR = M
      IP = -MR
      IPSTOR = 0
      I2R = 1
      JR = 2
      NR = N
      NLAST = N
      KR = 1
      LR = 0
      GO TO (101,103) ISTAG
  101 CONTINUE
      matQ(:MR,N) = HALF*matQ(:MR,N)
      GO TO (103,104) MIXBND
  103 CONTINUE
      IF (N <= 3) GO TO 155
  104 CONTINUE
      JR = 2*I2R
      NROD = 1
      IF ((NR/2)*2 == NR) NROD = 0
      SELECT CASE (MIXBND) 
      CASE DEFAULT
         JSTART = 1
      CASE (2) 
         JSTART = JR
         NROD = 1 - NROD
      END SELECT
      JSTOP = NLAST - JR
      IF (NROD == 0) JSTOP = JSTOP - I2R
      CALL COSGEN (I2R, 1, HALF, ZERO, TCOS, 4*NR)
      I2RBY2 = I2R/2
      IF (JSTOP < JSTART) THEN
         J = JR
      ELSE
         DO J = JSTART, JSTOP, JR
            JP1 = J + I2RBY2
            JP2 = J + I2R
            JP3 = JP2 + I2RBY2
            JM1 = J - I2RBY2
            JM2 = J - I2R
            JM3 = JM2 - I2RBY2
            IF (J == 1) THEN
               JM1 = JP1
               JM2 = JP2
               JM3 = JP3
            ENDIF
            IF (I2R == 1) THEN
               IF (J == 1) JM2 = JP2
               vecB(:MR) = TWO*matQ(:MR,J)
               matQ(:MR,J) = matQ(:MR,JM2) + matQ(:MR,JP2)
            ELSE
               DO I = 1, MR
                  FI = matQ(I,J)
                  matQ(I,J)=matQ(I,J)-matQ(I,JM1)-matQ(I,JP1)+ &
                            matQ(I,JM2)+matQ(I,JP2)
                  vecB(I) = FI + matQ(I,J) - matQ(I,JM3) - matQ(I,JP3)
               END DO
            ENDIF
            CALL TRIX (I2R, 0, MR, vecA, vecBB, vecC, vecB, TCOS, 4*NR, vecD, vecW)
            matQ(:MR,J) = matQ(:MR,J) + vecB(:MR)
!
!     END OF REDUCTION FOR REGULAR UNKNOWNS.
!
         END DO
!
!     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN.
!
         J = JSTOP + JR
      ENDIF
      NLAST = J
      JM1 = J - I2RBY2
      JM2 = J - I2R
      JM3 = JM2 - I2RBY2
      IF (NROD /= 0) THEN
!
!     ODD NUMBER OF UNKNOWNS
!
         IF (I2R == 1) THEN
            vecB(:MR) = FISTAG*matQ(:MR,J)
            matQ(:MR,J) = matQ(:MR,JM2)
         ELSE
            vecB(:MR) = matQ(:MR,J) + HALF*(matQ(:MR,JM2)- &
                        matQ(:MR,JM1)-matQ(:MR,JM3))
            IF (NRODPR == 0) THEN
               matQ(:MR,J) = matQ(:MR,JM2) + vecP(IP+1:MR+IP)
               IP = IP - MR
            ELSE
               matQ(:MR,J) = matQ(:MR,J) - matQ(:MR,JM1) + matQ(:MR,JM2)
            ENDIF
            IF (LR /= 0) THEN
               CALL COSGEN (LR, 1, HALF, FDEN, TCOS(KR+1), 4*NR-KR)
            ELSE
               vecB(:MR) = FISTAG*vecB(:MR)
            ENDIF
         ENDIF
         CALL COSGEN (KR, 1, HALF, FDEN, TCOS, 4*NR)
         CALL TRIX (KR, LR, MR, vecA, vecBB, vecC, vecB, TCOS, 4*NR, vecD, vecW)
         matQ(:MR,J) = matQ(:MR,J) + vecB(:MR)
         KR = KR + I2R
      ELSE
         JP1 = J + I2RBY2
         JP2 = J + I2R
         IF (I2R == 1) THEN
            vecB(:MR) = matQ(:MR,J)
            CALL TRIX (1, 0, MR, vecA, vecBB, vecC, vecB, TCOS, 4*NR, vecD, vecW)
            IP = 0
            IPSTOR = MR
            SELECT CASE (ISTAG) 
            CASE DEFAULT
               vecP(:MR) = vecB(:MR)
               vecB(:MR) = vecB(:MR) + matQ(:MR,N)
               TCOS(1) = ONE
               TCOS(2) = ZERO
               CALL TRIX (1, 1, MR, vecA, vecBB, vecC, vecB, TCOS, 4*NR, vecD, vecW)
               matQ(:MR,J) = matQ(:MR,JM2) + vecP(:MR) + vecB(:MR)
               GO TO 150
            CASE (1) 
               vecP(:MR) = vecB(:MR)
               matQ(:MR,J) = matQ(:MR,JM2) + TWO*matQ(:MR,JP2) + 3.*vecB(:MR)
               GO TO 150
            END SELECT
         ENDIF
         vecB(:MR) = matQ(:MR,J) + HALF*(matQ(:MR,JM2)- &
                     matQ(:MR,JM1)-matQ(:MR,JM3))
         IF (NRODPR == 0) THEN
            vecB(:MR) = vecB(:MR) + vecP(IP+1:MR+IP)
         ELSE
            vecB(:MR) = vecB(:MR) + matQ(:MR,JP2) - matQ(:MR,JP1)
         ENDIF
         CALL TRIX (I2R, 0, MR, vecA, vecBB, vecC, vecB, TCOS, 4*NR, vecD, vecW)
         IP = IP + MR
         IPSTOR = MAX0(IPSTOR,IP + MR)
         vecP(IP+1:MR+IP) = vecB(:MR) + HALF*(matQ(:MR,J)- &
                            matQ(:MR,JM1)-matQ(:MR,JP1))
         vecB(:MR) = vecP(IP+1:MR+IP) + matQ(:MR,JP2)
         IF (LR /= 0) THEN
            CALL COSGEN (LR, 1, HALF, FDEN, TCOS(I2R+1), 4*NR-I2R)
            CALL MERGE (TCOS, 0, I2R, I2R, LR, KR, 4*NR)
         ELSE
            DO I = 1, I2R
               II = KR + I
               TCOS(II) = TCOS(I)
            END DO
         ENDIF
         CALL COSGEN (KR, 1, HALF, FDEN, TCOS, 4*NR)
         IF (LR == 0) THEN
            GO TO (146,145) ISTAG
         ENDIF
  145    CONTINUE
         CALL TRIX (KR, KR, MR, vecA, vecBB, vecC, vecB, TCOS, 4*NR, vecD, vecW)
         GO TO 148
  146    CONTINUE
         vecB(:MR) = FISTAG*vecB(:MR)
  148    CONTINUE
         matQ(:MR,J) = matQ(:MR,JM2) + vecP(IP+1:MR+IP) + vecB(:MR)
  150    CONTINUE
         LR = KR
         KR = KR + JR
      ENDIF
      SELECT CASE (MIXBND) 
      CASE DEFAULT
         NR = (NLAST - 1)/JR + 1
         IF (NR <= 3) GO TO 155
      CASE (2) 
         NR = NLAST/JR
         IF (NR <= 1) GO TO 192
      END SELECT
      I2R = JR
      NRODPR = NROD
      GO TO 104
  155 CONTINUE
      J = 1 + JR
      JM1 = J - I2R
      JP1 = J + I2R
      JM2 = NLAST - I2R
      IF (NR /= 2) THEN
         IF (LR /= 0) GO TO 170
         IF (N == 3) THEN
!
!     CASE N = 3.
!
            GO TO (156,168) ISTAG
  156       CONTINUE
            vecB(:MR) = matQ(:MR,2)
            TCOS(1) = ZERO
            CALL TRIX (1, 0, MR, vecA, vecBB, vecC, vecB, TCOS, 4*NR, vecD, vecW)
            matQ(:MR,2) = vecB(:MR)
            vecB(:MR) = 4.*vecB(:MR) + matQ(:MR,1) + TWO*matQ(:MR,3)
            TCOS(1) = -TWO
            TCOS(2) = TWO
            I1 = 2
            I2 = 0
            CALL TRIX (I1, I2, MR, vecA, vecBB, vecC, vecB, TCOS, 4*NR, vecD, vecW)
            matQ(:MR,2) = matQ(:MR,2) + vecB(:MR)
            vecB(:MR) = matQ(:MR,1) + TWO*matQ(:MR,2)
            TCOS(1) = ZERO
            CALL TRIX (1, 0, MR, vecA, vecBB, vecC, vecB, TCOS, 4*NR, vecD, vecW)
            matQ(:MR,1) = vecB(:MR)
            JR = 1
            I2R = 0
            GO TO 194
         ENDIF
!
!     CASE N = 2**P+1
!
         GO TO (162,170) ISTAG
  162    CONTINUE
         vecB(:MR) = matQ(:MR,J) + HALF*matQ(:MR,1) - &
                     matQ(:MR,JM1) + matQ(:MR,NLAST) - matQ(:MR,JM2)
         CALL COSGEN (JR, 1, HALF, ZERO, TCOS, 4*NR)
         CALL TRIX (JR, 0, MR, vecA, vecBB, vecC, vecB, TCOS, 4*NR, vecD, vecW)
         matQ(:MR,J) = HALF*(matQ(:MR,J)-matQ(:MR,JM1)-matQ(:MR,JP1)) + vecB(:MR)
         vecB(:MR) = matQ(:MR,1) + TWO*matQ(:MR,NLAST) + 4.*matQ(:MR,J)
         JR2 = 2*JR
         CALL COSGEN (JR, 1, ZERO, ZERO, TCOS, 4*NR)
         TCOS(JR+1:JR*2) = -TCOS(JR:1:(-1))
         CALL TRIX (JR2, 0, MR, vecA, vecBB, vecC, vecB, TCOS, 4*NR, vecD, vecW)
         matQ(:MR,J) = matQ(:MR,J) + vecB(:MR)
         vecB(:MR) = matQ(:MR,1) + TWO*matQ(:MR,J)
         CALL COSGEN (JR, 1, HALF, ZERO, TCOS, 4*NR)
         CALL TRIX (JR, 0, MR, vecA, vecBB, vecC, vecB, TCOS, 4*NR, vecD, vecW)
         matQ(:MR,1) = HALF*matQ(:MR,1) - matQ(:MR,JM1) + vecB(:MR)
         GO TO 194
!
!     CASE OF GENERAL N WITH NR = 3 .
!
  168    CONTINUE
         vecB(:MR) = matQ(:MR,2)
         matQ(:MR,2) = ZERO
         vecB2(:MR) = matQ(:MR,3)
         vecB3(:MR) = matQ(:MR,1)
         JR = 1
         I2R = 0
         J = 2
         GO TO 177
  170    CONTINUE
         vecB(:MR) = HALF*matQ(:MR,1) - matQ(:MR,JM1) + matQ(:MR,J)
         IF (NROD == 0) THEN
            vecB(:MR) = vecB(:MR) + vecP(IP+1:MR+IP)
         ELSE
            vecB(:MR) = vecB(:MR) + matQ(:MR,NLAST) - matQ(:MR,JM2)
         ENDIF
         DO I = 1, MR
            T = HALF*(matQ(I,J)-matQ(I,JM1)-matQ(I,JP1))
            matQ(I,J) = T
            vecB2(I) = matQ(I,NLAST) + T
            vecB3(I) = matQ(I,1) + TWO*T
         END DO
  177    CONTINUE
         K1 = KR + 2*JR - 1
         K2 = KR + JR
         TCOS(K1+1) = -TWO
         K4 = K1 + 3 - ISTAG
         CALL COSGEN (K2 + ISTAG - 2, 1, ZERO, FNUM, TCOS(K4), 4*NR-K4+1)
         K4 = K1 + K2 + 1
         CALL COSGEN (JR - 1, 1, ZERO, ONE, TCOS(K4), 4*NR-K4+1)
         CALL MERGE (TCOS, K1, K2, K1 + K2, JR - 1, 0, 4*NR)
         K3 = K1 + K2 + LR
         CALL COSGEN (JR, 1, HALF, ZERO, TCOS(K3+1), 4*NR-K3)
         K4 = K3 + JR + 1
         CALL COSGEN (KR, 1, HALF, FDEN, TCOS(K4), 4*NR-K4+1)
         CALL MERGE (TCOS, K3, JR, K3 + JR, KR, K1, 4*NR)
         IF (LR /= 0) THEN
            CALL COSGEN (LR, 1, HALF, FDEN, TCOS(K4), 4*NR-K4+1)
            CALL MERGE (TCOS, K3, JR, K3 + JR, LR, K3 - LR, 4*NR)
            CALL COSGEN (KR, 1, HALF, FDEN, TCOS(K4), 4*NR-K4+1)
         ENDIF
         K3 = KR
         K4 = KR
         CALL TRI3 (MR, vecA, vecBB, vecC, K, vecB, vecB2, vecB3, &
                    TCOS, 4*NR, vecD, vecW, vecW2, vecW3)
         vecB(:MR) = vecB(:MR) + vecB2(:MR) + vecB3(:MR)
         TCOS(1) = TWO
         CALL TRIX (1, 0, MR, vecA, vecBB, vecC, vecB, TCOS, 4*NR, vecD, vecW)
         matQ(:MR,J) = matQ(:MR,J) + vecB(:MR)
         vecB(:MR) = matQ(:MR,1) + TWO*matQ(:MR,J)
         CALL COSGEN (JR, 1, HALF, ZERO, TCOS, 4*NR)
         CALL TRIX (JR, 0, MR, vecA, vecBB, vecC, vecB, TCOS, 4*NR, vecD, vecW)
         IF (JR == 1) THEN
            matQ(:MR,1) = vecB(:MR)
            GO TO 194
         ENDIF
         matQ(:MR,1) = HALF*matQ(:MR,1) - matQ(:MR,JM1) + vecB(:MR)
         GO TO 194
      ENDIF
      IF (N == 2) THEN
!
!     CASE  N = 2
!
         vecB(:MR) = matQ(:MR,1)
         TCOS(1) = ZERO
         CALL TRIX (1, 0, MR, vecA, vecBB, vecC, vecB, TCOS, 4*NR, vecD, vecW)
         matQ(:MR,1) = vecB(:MR)
         vecB(:MR) = TWO*(matQ(:MR,2)+vecB(:MR))*FISTAG
         TCOS(1) = -FISTAG
         TCOS(2) = TWO
         CALL TRIX (2, 0, MR, vecA, vecBB, vecC, vecB, TCOS, 4*NR, vecD, vecW)
         matQ(:MR,1) = matQ(:MR,1) + vecB(:MR)
         JR = 1
         I2R = 0
         GO TO 194
      ENDIF
      vecB3(:MR) = ZERO
      vecB(:MR) = matQ(:MR,1) + TWO*vecP(IP+1:MR+IP)
      matQ(:MR,1) = HALF*matQ(:MR,1) - matQ(:MR,JM1)
      vecB2(:MR) = TWO*(matQ(:MR,1)+matQ(:MR,NLAST))
      K1 = KR + JR - 1
      TCOS(K1+1) = -TWO
      K4 = K1 + 3 - ISTAG
      CALL COSGEN (KR + ISTAG - 2, 1, ZERO, FNUM, TCOS(K4), 4*NR-K4+1)
      K4 = K1 + KR + 1
      CALL COSGEN (JR - 1, 1, ZERO, ONE, TCOS(K4), 4*NR-K4+1)
      CALL MERGE (TCOS, K1, KR, K1 + KR, JR - 1, 0, 4*NR)
      CALL COSGEN (KR, 1, HALF, FDEN, TCOS(K1+1), 4*NR-K1)
      K2 = KR
      K4 = K1 + K2 + 1
      CALL COSGEN (LR, 1, HALF, FDEN, TCOS(K4), 4*NR-K4+1)
      K3 = LR
      K4 = 0
      CALL TRI3 (MR, vecA, vecBB, vecC, K, vecB, vecB2, vecB3, TCOS, &
                 4*NR, vecD, vecW, vecW2, vecW3)
      vecB(:MR) = vecB(:MR) + vecB2(:MR)
      TCOS(1) = TWO
      CALL TRIX (1, 0, MR, vecA, vecBB, vecC, vecB, TCOS, 4*NR, vecD, vecW)
      matQ(:MR,1) = matQ(:MR,1) + vecB(:MR)
      GO TO 194
  192 CONTINUE
      vecB(:MR) = matQ(:MR,NLAST)
      GO TO 196
  194 CONTINUE
      J = NLAST - JR
      vecB(:MR) = matQ(:MR,NLAST) + matQ(:MR,J)
  196 CONTINUE
      JM2 = NLAST - I2R
      IF (JR == 1) THEN
         matQ(:MR,NLAST) = ZERO
      ELSE
         IF (NROD == 0) THEN
            matQ(:MR,NLAST) = vecP(IP+1:MR+IP)
            IP = IP - MR
         ELSE
            matQ(:MR,NLAST) = matQ(:MR,NLAST) - matQ(:MR,JM2)
         ENDIF
      ENDIF
      CALL COSGEN (KR, 1, HALF, FDEN, TCOS, 4*NR)
      CALL COSGEN (LR, 1, HALF, FDEN, TCOS(KR+1), 4*NR-KR)
      IF (LR == 0) THEN
         vecB(:MR) = FISTAG*vecB(:MR)
      ENDIF
      CALL TRIX (KR, LR, MR, vecA, vecBB, vecC, vecB, TCOS, 4*NR, vecD, vecW)
      matQ(:MR,NLAST) = matQ(:MR,NLAST) + vecB(:MR)
      NLASTP = NLAST
  206 CONTINUE
      JSTEP = JR
      JR = I2R
      I2R = I2R/2
      IF (JR == 0) GO TO 222
      SELECT CASE (MIXBND) 
      CASE DEFAULT
         JSTART = 1 + JR
      CASE (2) 
         JSTART = JR
      END SELECT
      KR = KR - JR
      IF (NLAST + JR <= N) THEN
         KR = KR - JR
         NLAST = NLAST + JR
         JSTOP = NLAST - JSTEP
      ELSE
         JSTOP = NLAST - JR
      ENDIF
      LR = KR - JR
      CALL COSGEN (JR, 1, HALF, ZERO, TCOS, 4*NR)
      DO J = JSTART, JSTOP, JSTEP
         JM2 = J - JR
         JP2 = J + JR
         IF (J == JR) THEN
            vecB(:MR) = matQ(:MR,J) + matQ(:MR,JP2)
         ELSE
            vecB(:MR) = matQ(:MR,J) + matQ(:MR,JM2) + matQ(:MR,JP2)
         ENDIF
         IF (JR == 1) THEN
            matQ(:MR,J) = ZERO
         ELSE
            JM1 = J - I2R
            JP1 = J + I2R
            matQ(:MR,J) = HALF*(matQ(:MR,J)-matQ(:MR,JM1)-matQ(:MR,JP1))
         ENDIF
         CALL TRIX (JR, 0, MR, vecA, vecBB, vecC, vecB, TCOS, 4*NR, vecD, vecW)
         matQ(:MR,J) = matQ(:MR,J) + vecB(:MR)
      END DO
      NROD = 1
      IF (NLAST + I2R <= N) NROD = 0
      IF (NLASTP /= NLAST) GO TO 194
      GO TO 206
  222 CONTINUE
      vecW(1) = IPSTOR
      RETURN 
      END SUBROUTINE POISN2


      SUBROUTINE POISP2(M, N, vecA, vecBB, vecC, matQ, IDIMQ, &
                        vecB, vecB2, vecB3, vecW, vecW2, vecW3, &
                        vecD, TCOS, vecP, IP)
      implicit none
      
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, &
                                     ONE  = 1.0D0, TWO  = 2.0D0, &
                                     FOUR = 4.0D0

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: M, N, IDIMQ, IP
      DOUBLE PRECISION, DIMENSION(M), INTENT(IN) :: vecA, vecBB, vecC
      DOUBLE PRECISION, DIMENSION(4*N), INTENT(INOUT) :: TCOS
      DOUBLE PRECISION, DIMENSION(IDIMQ,N), INTENT(INOUT) :: matQ
      DOUBLE PRECISION, DIMENSION(IP),INTENT(INOUT) :: vecP
      DOUBLE PRECISION, DIMENSION(M),INTENT(INOUT) :: vecB, vecB2, vecB3, &
                                                      vecD, vecW, vecW2, vecW3
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER          :: MR, NR, NRM1, J, NRMJ, NRPJ, I, IPSTOR, LH
      DOUBLE PRECISION :: S, T
!-----------------------------------------------
!
!     SUBROUTINE TO SOLVE POISSON EQUATION WITH PERIODIC BOUNDARY
!     CONDITIONS.
!
      MR = M
      NR = (N + 1)/2
      NRM1 = NR - 1
      IF (2*NR == N) THEN
!
!     EVEN NUMBER OF UNKNOWNS
!
         DO J = 1, NRM1
            NRMJ = NR - J
            NRPJ = NR + J
            DO I = 1, MR
               S = matQ(I,NRMJ) - matQ(I,NRPJ)
               T = matQ(I,NRMJ) + matQ(I,NRPJ)
               matQ(I,NRMJ) = S
               matQ(I,NRPJ) = T
            END DO
         END DO
         matQ(:MR,NR) = TWO*matQ(:MR,NR)
         matQ(:MR,N)  = TWO*matQ(:MR,N)
         CALL POISD2 (MR, NRM1, 1, vecA, vecBB, vecC, matQ, IDIMQ, &
                      vecB, vecW, vecD, TCOS, vecP, IP)
         IPSTOR = vecW(1)
         CALL POISN2 (MR, NR + 1, 1, 1, vecA, vecBB, vecC, matQ(:,NR:),  &
                      IDIMQ, vecB, vecB2, vecB3, vecW, vecW2, vecW3, vecD, &
                      TCOS, vecP, IP)
         IPSTOR = MAX0(IPSTOR,INT(vecW(1)))
         DO J = 1, NRM1
            NRMJ = NR - J
            NRPJ = NR + J
            DO I = 1, MR
               S = HALF*(matQ(I,NRPJ)+matQ(I,NRMJ))
               T = HALF*(matQ(I,NRPJ)-matQ(I,NRMJ))
               matQ(I,NRMJ) = S
               matQ(I,NRPJ) = T
            END DO
         END DO
         matQ(:MR,NR) = HALF*matQ(:MR,NR)
         matQ(:MR,N)  = HALF*matQ(:MR,N)
      ELSE
         DO J = 1, NRM1
            NRPJ = N + 1 - J
            DO I = 1, MR
               S = matQ(I,J) - matQ(I,NRPJ)
               T = matQ(I,J) + matQ(I,NRPJ)
               matQ(I,J) = S
               matQ(I,NRPJ) = T
            END DO
         END DO
         matQ(:MR,NR) = TWO*matQ(:MR,NR)
         LH = NRM1/2
         DO J = 1, LH
            NRMJ = NR - J
            DO I = 1, MR
               S = matQ(I,J)
               matQ(I,J) = matQ(I,NRMJ)
               matQ(I,NRMJ) = S
            END DO
         END DO
         CALL POISD2 (MR, NRM1, 2, vecA, vecBB, vecC, matQ, IDIMQ,  &
                      vecB, vecW, vecD, TCOS, vecP, IP)
         IPSTOR = vecW(1)
         CALL POISN2 (MR, NR, 2, 1, vecA, vecBB, vecC, matQ(:,NR:), IDIMQ, &
                      vecB, vecB2, vecB3, vecW, vecW2, vecW3, vecD, TCOS, &
                      vecP, IP)
         IPSTOR = MAX0(IPSTOR,INT(vecW(1)))
         DO J = 1, NRM1
            NRPJ = NR + J
            DO I = 1, MR
               S = HALF*(matQ(I,NRPJ)+matQ(I,J))
               T = HALF*(matQ(I,NRPJ)-matQ(I,J))
               matQ(I,NRPJ) = T
               matQ(I,J) = S
            END DO
         END DO
         matQ(:MR,NR) = HALF*matQ(:MR,NR)
         DO J = 1, LH
            NRMJ = NR - J
            DO I = 1, MR
               S = matQ(I,J)
               matQ(I,J)    = matQ(I,NRMJ)
               matQ(I,NRMJ) = S
            END DO
         END DO
      ENDIF
      vecW(1) = IPSTOR
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
      END SUBROUTINE POISP2
