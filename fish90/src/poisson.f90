      SUBROUTINE POISTGG(NPEROD,N,MPEROD,M,vecA,vecB,vecC,IDIMY,matY, &
                         IERROR,W,IW)
      implicit none
      
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, &
                                     ONE  = 1.0D0, TWO  = 2.0D0, &
                                     FOUR = 4.0D0

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)  :: NPEROD, N, MPEROD, M, IDIMY,IW
      INTEGER, INTENT(OUT) :: IERROR
      DOUBLE PRECISION, DIMENSION(M), INTENT(IN)         :: vecA,vecB,vecC
      DOUBLE PRECISION, DIMENSION(IDIMY,N),INTENT(INOUT) :: matY
      DOUBLE PRECISION, DIMENSION(IW)                    :: W
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IWBA, IWBB, IWBC, IWB2, IWB3, IWW1, IWW2, IWW3, &
                 IWD,IWTCOS, IWP, I, K, J, NP, MP, IPSTOR, &
                 IREV, MH, MHM1, MODD, NBY2, MSKIP
      DOUBLE PRECISION :: A1
!-----------------------------------------------
      IERROR = 0

      IWBA = M + 1
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
      DO I = 1, M
         K = IWBA + I - 1
         ! write(*,*) '0: i,k:',i,k
         ! write(*,*) '0: vecA(I):',vecA(I)
         W(K) = -vecA(I)
         ! write(*,*) 'A: i,k,w(k):',i,k,w(k)
         K = IWBC + I - 1
         W(K) = -vecC(I)
         ! write(*,*) 'B: i,k,w(k):',i,k,w(k)
         K = IWBB + I - 1
         W(K) = TWO - vecB(I)
         ! write(*,*) 'C: i,k,w(k):',i,k,w(k)
         matY(I,:N) = -matY(I,:N)
         ! write(*,*) 'D: Y(i,:N):',Y(I,:N)
      END DO
      NP = NPEROD
      MP = MPEROD + 1
      GO TO (110,107) MP
  107 CONTINUE
      GO TO (108,108,108,119) NPEROD
  108 CONTINUE
      CALL POSTG2 (NP, N, M, W(IWBA:IWBA+M-1), W(IWBB:IWBB+M-1), &
                   W(IWBC:IWBC+M-1), IDIMY, matY, W(1:M), W(IWB2:IWB2+M-1),  &
                   W(IWB3:IWB3+M-1), W(IWW1:IWW1+M-1), W(IWW2:IWW2+M-1), &
                   W(IWW3:IWW3+M-1), W(IWD:IWD+M-1), W(IWTCOS:IWTCOS+4*N-1), &
                   W(IWP:),IW-IWP)
      IPSTOR = W(IWW1)
      IREV = 2
      IF (NPEROD == 4) GO TO 120
  109 CONTINUE
      GO TO (123,129) MP
  110 CONTINUE
      MH = (M + 1)/2
      MHM1 = MH - 1
      MODD = 1
      IF (MH*2 == M) MODD = 2
      DO J = 1, N
         DO I = 1, MHM1
            W(I) = matY(MH-I,J) - matY(I+MH,J)
            W(I+MH) = matY(MH-I,J) + matY(I+MH,J)
         END DO
         W(MH) = TWO*matY(MH,J)
         GO TO (113,112) MODD
  112    CONTINUE
         W(M) = TWO*matY(M,J)
  113    CONTINUE
         matY(:M,J) = W(:M)
      END DO
      K = IWBC + MHM1 - 1
      I = IWBA + MHM1
      W(K) = 0.
      W(I) = 0.
      W(K+1) = TWO*W(K+1)
      SELECT CASE (MODD) 
      CASE DEFAULT
         K = IWBB + MHM1 - 1
         W(K) = W(K) - W(I-1)
         W(IWBC-1) = W(IWBC-1) + W(IWBB-1)
      CASE (2) 
         W(IWBB-1) = W(K+1)
      END SELECT
      GO TO 107
  119 CONTINUE
      IREV = 1
      NBY2 = N/2
      NP = 2
  120 CONTINUE
      DO J = 1, NBY2
         MSKIP = N + 1 - J
         DO I = 1, M
            A1 = matY(I,J)
            matY(I,J) = matY(I,MSKIP)
            matY(I,MSKIP) = A1
         END DO
      END DO
      GO TO (108,109) IREV
  123 CONTINUE
      DO J = 1, N
         W(MH-1:MH-MHM1:(-1)) = HALF*(matY(MH+1:MHM1+MH,J)+matY(:MHM1,J))
         W(MH+1:MHM1+MH) = HALF*(matY(MH+1:MHM1+MH,J)-matY(:MHM1,J))
         W(MH) = HALF*matY(MH,J)
         GO TO (126,125) MODD
  125    CONTINUE
         W(M) = HALF*matY(M,J)
  126    CONTINUE
         matY(:M,J) = W(:M)
      END DO
  129 CONTINUE
      W(1) = IPSTOR + IWP - 1
!
      END SUBROUTINE POISTGG

      SUBROUTINE POSTG2(NPEROD, N, M, vecA, vecBB, vecC, IDIMQ, matQ, &
                        vecB, vecB2, vecB3, vecW,vecW2, vecW3, vecD,  &
                        TCOS, vecP,IvecP)
      implicit none
      
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, HALF = 0.5D0, &
                                     ONE  = 1.0D0, TWO  = 2.0D0, &
                                     FOUR = 4.0D0

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NPEROD,N,M,IDIMQ,IvecP
      DOUBLE PRECISION, DIMENSION(M),INTENT(IN)        :: vecA,vecBB,vecC
      DOUBLE PRECISION, DIMENSION(M),INTENT(OUT)       :: vecB
      DOUBLE PRECISION, DIMENSION(IDIMQ,N),INTENT(INOUT) :: matQ
      DOUBLE PRECISION, DIMENSION(M),INTENT(INOUT)     :: vecB2,vecB3,vecD, &
                                                          vecW2,vecW3, vecW
      DOUBLE PRECISION, DIMENSION(IvecP),INTENT(INOUT) :: vecP
      DOUBLE PRECISION, DIMENSION(4*N),INTENT(INOUT)   :: TCOS
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(4) :: K
      INTEGER :: K1, K2, K3, K4, NP, MR, IP, IPSTOR, I2R, JR, NR, NLAST, &
                 KR, LR, NROD, JSTART, JSTOP, I2RBY2, &
                 J, IJUMP, JP1, JP2, JP3, JM1, JM2, JM3, I, NRODPR, II,  &
                 NLASTP, JSTEP
      DOUBLE PRECISION :: FNUM, FNUM2, FI, T
!-----------------------------------------------
!
!     SUBROUTINE TO SOLVE POISSON'S EQUATION ON A STAGGERED GRID.
!
      EQUIVALENCE (K(1), K1), (K(2), K2), (K(3), K3), (K(4), K4)
      NP = NPEROD
      FNUM = HALF*DBLE(NP/3)
      FNUM2 = HALF*DBLE(NP/2)
      MR = M
      IP = -MR
      IPSTOR = 0
      I2R = 1
      JR = 2
      NR = N
      NLAST = N
      KR = 1
      LR = 0
      IF (NR > 3) THEN
  101    CONTINUE
         JR = 2*I2R
         NROD = 1
         IF ((NR/2)*2 == NR) NROD = 0
         JSTART = 1
         JSTOP = NLAST - JR
         IF (NROD == 0) JSTOP = JSTOP - I2R
         I2RBY2 = I2R/2
         IF (JSTOP < JSTART) THEN
            J = JR
         ELSE
            IJUMP = 1
            DO J = JSTART, JSTOP, JR
               JP1 = J + I2RBY2
               JP2 = J + I2R
               JP3 = JP2 + I2RBY2
               JM1 = J - I2RBY2
               JM2 = J - I2R
               JM3 = JM2 - I2RBY2
               IF (J == 1) THEN
                  CALL COSGEN (I2R, 1, FNUM, HALF, TCOS, size(TCOS) )
                  IF (I2R == 1) THEN
                     vecB(:MR)   = matQ(:MR,1)
                     matQ(:MR,1) = matQ(:MR,2)
                     GO TO 112
                  ENDIF
                  vecB(:MR)   = matQ(:MR,1) + &
                                HALF*(matQ(:MR,JP2)-matQ(:MR,JP1)-matQ(:MR,JP3))
                  matQ(:MR,1) = matQ(:MR,JP2) + matQ(:MR,1) - matQ(:MR,JP1)
                  GO TO 112
               ENDIF
               GO TO (107,108) IJUMP
  107          CONTINUE
               IJUMP = 2
               CALL COSGEN (I2R, 1, HALF, ZERO, TCOS, size(TCOS) )
  108          CONTINUE
               IF (I2R == 1) THEN
                  vecB(:MR)   = TWO*matQ(:MR,J)
                  matQ(:MR,J) = matQ(:MR,JM2) + matQ(:MR,JP2)
               ELSE
                  DO I = 1, MR
                     FI = matQ(I,J)
                     matQ(I,J)=matQ(I,J)-matQ(I,JM1)-matQ(I,JP1)+ &
                               matQ(I,JM2)+matQ(I,JP2)
                     vecB(I)  =FI + matQ(I,J) - matQ(I,JM3) - matQ(I,JP3)
                  END DO
               ENDIF
  112          CONTINUE
               CALL TRIX (I2R, 0, MR, vecA, vecBB, vecC, vecB, TCOS, &
                          size(TCOS), vecD, vecW)
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
               vecB(:MR)   = matQ(:MR,J)
               matQ(:MR,J) = matQ(:MR,JM2)
            ELSE
               vecB(:MR)=matQ(:MR,J)+ &
                         HALF*(matQ(:MR,JM2)-matQ(:MR,JM1)-matQ(:MR,JM3))
               IF (NRODPR == 0) THEN
                  matQ(:MR,J) = matQ(:MR,JM2) + vecP(IP+1:MR+IP)
                  IP = IP - MR
               ELSE
                  matQ(:MR,J) = matQ(:MR,J) - matQ(:MR,JM1) + matQ(:MR,JM2)
               ENDIF
               IF (LR /= 0) CALL COSGEN (LR, 1, FNUM2, HALF, TCOS(KR+1), &
                                         size(tcos)-KR)
            ENDIF
            CALL COSGEN (KR, 1, FNUM2, HALF, TCOS, size(tcos))
            CALL TRIX (KR, LR, MR, vecA, vecBB, vecC, vecB, TCOS,  &
                       size(TCOS), vecD, vecW)
            matQ(:MR,J) = matQ(:MR,J) + vecB(:MR)
            KR = KR + I2R
         ELSE
            JP1 = J + I2RBY2
            JP2 = J + I2R
            IF (I2R == 1) THEN
               vecB(:MR) = matQ(:MR,J)
               TCOS(1) = 0.
               CALL TRIX (1, 0, MR, vecA, vecBB, vecC, vecB, TCOS,  &
                          size(TCOS), vecD, vecW)
               IP = 0
               IPSTOR = MR
               vecP(:MR) = vecB(:MR)
               vecB(:MR) = vecB(:MR) + matQ(:MR,N)
               TCOS(1) = -ONE + TWO*DBLE(NP/2)
               TCOS(2) = 0.
               CALL TRIX (1, 1, MR, vecA, vecBB, vecC, vecB, TCOS,  &
                          size(TCOS), vecD, vecW)
               matQ(:MR,J) = matQ(:MR,JM2) + vecP(:MR) + vecB(:MR)
            ELSE
               vecB(:MR)=matQ(:MR,J)+ &
                         HALF*(matQ(:MR,JM2)-matQ(:MR,JM1)-matQ(:MR,JM3))
               IF (NRODPR == 0) THEN
                  vecB(:MR) = vecB(:MR) + vecP(IP+1:MR+IP)
               ELSE
                  vecB(:MR) = vecB(:MR) + matQ(:MR,JP2) - matQ(:MR,JP1)
               ENDIF
               CALL COSGEN (I2R, 1, HALF, ZERO, TCOS, size(tcos))
               CALL TRIX (I2R, 0, MR, vecA, vecBB, vecC, vecB, TCOS,  &
                          size(TCOS), vecD, vecW)
               IP = IP + MR
               IPSTOR = MAX0(IPSTOR,IP + MR)
               vecP(IP+1:MR+IP) = vecB(:MR) +  &
                                  HALF*(matQ(:MR,J)-matQ(:MR,JM1)-matQ(:MR,JP1))
               vecB(:MR) = vecP(IP+1:MR+IP) + matQ(:MR,JP2)
               IF (LR /= 0) THEN
                  CALL COSGEN (LR, 1, FNUM2, HALF, TCOS(I2R+1), size(tcos)-I2R)
                  CALL MERGE (TCOS, 0, I2R, I2R, LR, KR, size(TCOS))
               ELSE
                  DO I = 1, I2R
                     II = KR + I
                     TCOS(II) = TCOS(I)
                  END DO
               ENDIF
               CALL COSGEN (KR, 1, FNUM2, HALF, TCOS, size(TCOS))
               CALL TRIX (KR, KR, MR, vecA, vecBB, vecC, vecB, TCOS,  &
                          size(TCOS), vecD, vecW)
               matQ(:MR,J) = matQ(:MR,JM2) + vecP(IP+1:MR+IP) + vecB(:MR)
            ENDIF
            LR = KR
            KR = KR + JR
         ENDIF
         NR = (NLAST - 1)/JR + 1
         IF (NR <= 3) GO TO 142
         I2R = JR
         NRODPR = NROD
         GO TO 101
      ENDIF ! IF (NR > 3)
  142 CONTINUE
      J = 1 + JR
      JM1 = J - I2R
      JP1 = J + I2R
      JM2 = NLAST - I2R
      IF (NR /= 2) THEN
         IF (LR == 0) THEN
            IF (N == 3) THEN
!
!     CASE N = 3.
!
               GO TO (143,148,143) NP
  143          CONTINUE
               vecB(:MR) = matQ(:MR,2)
               vecB2(:MR) = matQ(:MR,1) + matQ(:MR,3)
               vecB3(:MR) = 0.
               SELECT CASE (NP) 
               CASE DEFAULT
                  TCOS(1) = -ONE
                  TCOS(2) = ONE
                  K1 = 1
               CASE (1:2) 
                  TCOS(1) = -TWO
                  TCOS(2) = ONE
                  TCOS(3) = -ONE
                  K1 = 2
               END SELECT
               K2 = 1
               K3 = 0
               K4 = 0
               GO TO 150
  148          CONTINUE
               vecB(:MR) = matQ(:MR,2)
               vecB2(:MR) = matQ(:MR,3)
               vecB3(:MR) = matQ(:MR,1)
               CALL COSGEN (3, 1, HALF, ZERO, TCOS, size(TCOS))
               TCOS(4) = -ONE
               TCOS(5) = ONE
               TCOS(6) = -ONE
               TCOS(7) = ONE
               K1 = 3
               K2 = 2
               K3 = 1
               K4 = 1
  150          CONTINUE
               CALL TRI3(MR,vecA,vecBB,vecC,K,vecB,vecB2,vecB3,TCOS, &
                         size(TCOS),vecD,vecW,vecW2,vecW3)
               vecB(:MR) = vecB(:MR) + vecB2(:MR) + vecB3(:MR)
               GO TO (153,153,152) NP
  152          CONTINUE
               TCOS(1) = TWO
               CALL TRIX (1, 0, MR, vecA, vecBB, vecC, vecB, TCOS,  &
                          size(TCOS), vecD, vecW)
  153          CONTINUE
               matQ(:MR,2) = vecB(:MR)
               vecB(:MR) = matQ(:MR,1) + vecB(:MR)
               TCOS(1) = -ONE + FOUR*FNUM
               CALL TRIX (1, 0, MR, vecA, vecBB, vecC, vecB, TCOS,  &
                          size(TCOS), vecD, vecW)
               matQ(:MR,1) = vecB(:MR)
               JR = 1
               I2R = 0
               GO TO 188
            ENDIF
!
!     CASE N = 2**P+1
!
            vecB(:MR)=matQ(:MR,J)+matQ(:MR,1)-matQ(:MR,JM1)+ &
                      matQ(:MR,NLAST)-matQ(:MR,JM2)
            GO TO (158,160,158) NP
  158       CONTINUE
            vecB2(:MR) = matQ(:MR,1) + matQ(:MR,NLAST) + &
                         matQ(:MR,J) - matQ(:MR,JM1) -matQ(:MR,JP1)
            vecB3(:MR) = 0.
            K1 = NLAST - 1
            K2 = NLAST + JR - 1
            CALL COSGEN (JR - 1, 1, ZERO, ONE, TCOS(NLAST), size(TCOS)-NLAST+1)
            TCOS(K2) = TWO*DBLE(NP - 2)
            CALL COSGEN (JR, 1, HALF - FNUM, HALF, TCOS(K2+1),size(TCOS)-K2)
            K3 = (3 - NP)/2
            CALL MERGE (TCOS, K1, JR - K3, K2 - K3, JR + K3, 0, size(TCOS))
            K1 = K1 - 1 + K3
            CALL COSGEN (JR, 1, FNUM, HALF, TCOS(K1+1),size(TCOS)-K1)
            K2 = JR
            K3 = 0
            K4 = 0
            GO TO 162
  160       CONTINUE
            DO I = 1, MR
               FI = HALF*(matQ(I,J)-matQ(I,JM1)-matQ(I,JP1))
               vecB2(I) = matQ(I,1) + FI
               vecB3(I) = matQ(I,NLAST) + FI
            END DO
            K1 = NLAST + JR - 1
            K2 = K1 + JR - 1
            CALL COSGEN (JR - 1, 1, ZERO, ONE, TCOS(K1+1), size(TCOS)-K1)
            CALL COSGEN (NLAST, 1, HALF, ZERO, TCOS(K2+1), size(TCOS)-K2)
            CALL MERGE (TCOS, K1, JR - 1, K2, NLAST, 0, size(TCOS))
            K3 = K1 + NLAST - 1
            K4 = K3 + JR
            CALL COSGEN (JR, 1, HALF, HALF, TCOS(K3+1),SIZE(TCOS)-K3)
            CALL COSGEN (JR, 1, ZERO, HALF, TCOS(K4+1),SIZE(TCOS)-K4)
            CALL MERGE (TCOS, K3, JR, K4, JR, K1, size(TCOS) )
            K2 = NLAST - 1
            K3 = JR
            K4 = JR
  162       CONTINUE
            CALL TRI3 (MR, vecA, vecBB, vecC, K, vecB, vecB2, vecB3,  &
                       TCOS, size(TCOS), vecD, vecW, vecW2, vecW3)
            vecB(:MR) = vecB(:MR) + vecB2(:MR) + vecB3(:MR)
            IF (NP == 3) THEN
               TCOS(1) = TWO
               CALL TRIX (1, 0, MR, vecA, vecBB, vecC, vecB, TCOS,  &
                          size(TCOS), vecD, vecW)
            ENDIF
            matQ(:MR,J) = vecB(:MR) +  &
                          HALF*(matQ(:MR,J)-matQ(:MR,JM1)-matQ(:MR,JP1))
            vecB(:MR) = matQ(:MR,J) + matQ(:MR,1)
            CALL COSGEN (JR, 1, FNUM, HALF, TCOS,size(TCOS))
            CALL TRIX (JR, 0, MR, vecA, vecBB, vecC, vecB, TCOS,  &
                          size(TCOS), vecD, vecW)
            matQ(:MR,1) = matQ(:MR,1) - matQ(:MR,JM1) + vecB(:MR)
            GO TO 188
         ENDIF
!
!     CASE OF GENERAL N WITH NR = 3 .
!
         vecB(:MR) = matQ(:MR,1) - matQ(:MR,JM1) + matQ(:MR,J)
         IF (NROD == 0) THEN
            vecB(:MR) = vecB(:MR) + vecP(IP+1:MR+IP)
         ELSE
            vecB(:MR) = vecB(:MR) + matQ(:MR,NLAST) - matQ(:MR,JM2)
         ENDIF
         DO I = 1, MR
            T = HALF*(matQ(I,J)-matQ(I,JM1)-matQ(I,JP1))
            matQ(I,J) = T
            vecB2(I) = matQ(I,NLAST) + T
            vecB3(I) = matQ(I,1) + T
         END DO
         K1 = KR + 2*JR
         CALL COSGEN (JR - 1, 1, ZERO, ONE, TCOS(K1+1),size(TCOS)-K1)
         K2 = K1 + JR
         TCOS(K2) = TWO*DBLE(NP - 2)
         K4 = (NP - 1)*(3 - NP)
         K3 = K2 + 1 - K4
         CALL COSGEN(KR+JR+K4,1,HALF*DBLE(K4),ONE-DBLE(K4),TCOS(K3), &
                     size(TCOS)-K3+1 )
         K4 = 1 - NP/3
         CALL MERGE (TCOS, K1, JR - K4, K2 - K4, KR + JR + K4, 0, size(TCOS) )
         IF (NP == 3) K1 = K1 - 1
         K2 = KR + JR
         K4 = K1 + K2
         CALL COSGEN (KR, 1, FNUM2, HALF, TCOS(K4+1),size(TCOS)-K4)
         K3 = K4 + KR
         CALL COSGEN (JR, 1, FNUM, HALF, TCOS(K3+1),size(TCOS)-K3)
         CALL MERGE (TCOS, K4, KR, K3, JR, K1, size(TCOS) )
         K4 = K3 + JR
         CALL COSGEN (LR, 1, FNUM2, HALF, TCOS(K4+1),size(TCOS)-K4)
         CALL MERGE (TCOS, K3, JR, K4, LR, K1 + K2, size(TCOS) )
         CALL COSGEN (KR, 1, FNUM2, HALF, TCOS(K3+1),size(TCOS)-K3)
         K3 = KR
         K4 = KR
         CALL TRI3 (MR, vecA, vecBB, vecC, K, vecB, vecB2, vecB3, &
                    TCOS, SIZE(TCOS), vecD, vecW, vecW2, vecW3)
         vecB(:MR) = vecB(:MR) + vecB2(:MR) + vecB3(:MR)
         IF (NP == 3) THEN
            TCOS(1) = TWO
            CALL TRIX (1, 0, MR, vecA, vecBB, vecC, vecB, TCOS,  &
                          size(TCOS), vecD, vecW)
         ENDIF
         matQ(:MR,J) = matQ(:MR,J) + vecB(:MR)
         vecB(:MR) = matQ(:MR,1) + matQ(:MR,J)
         CALL COSGEN (JR, 1, FNUM, HALF, TCOS,size(TCOS))
         CALL TRIX (JR, 0, MR, vecA, vecBB, vecC, vecB, TCOS,  &
                          size(TCOS), vecD, vecW)
         IF (JR == 1) THEN
            matQ(:MR,1) = vecB(:MR)
            GO TO 188
         ENDIF
         matQ(:MR,1) = matQ(:MR,1) - matQ(:MR,JM1) + vecB(:MR)
         GO TO 188
      ENDIF
      vecB3(:MR) = 0.
      vecB (:MR) = matQ(:MR,1) + vecP(IP+1:MR+IP)
      matQ (:MR,1) = matQ(:MR,1) - matQ(:MR,JM1)
      vecB2(:MR) = matQ(:MR,1) + matQ(:MR,NLAST)
      K1 = KR + JR
      K2 = K1 + JR
      CALL COSGEN (JR - 1, 1, ZERO, ONE, TCOS(K1+1),SIZE(TCOS)-K1)
      GO TO (182,183,182) NP
  182 CONTINUE
      TCOS(K2) = TWO*DBLE(NP - 2)
      CALL COSGEN (KR, 1, ZERO, ONE, TCOS(K2+1),SIZE(TCOS)-K2)
      GO TO 184
  183 CONTINUE
      CALL COSGEN (KR + 1, 1, HALF, ZERO, TCOS(K2),SIZE(TCOS)-K2)
  184 CONTINUE
      K4 = 1 - NP/3
      CALL MERGE (TCOS, K1, JR - K4, K2 - K4, KR + K4, 0, size(TCOS) )
      IF (NP == 3) K1 = K1 - 1
      K2 = KR
      CALL COSGEN (KR, 1, FNUM2, HALF, TCOS(K1+1),SIZE(TCOS)-K1)
      K4 = K1 + KR
      CALL COSGEN (LR, 1, FNUM2, HALF, TCOS(K4+1),SIZE(TCOS)-K4)
      K3 = LR
      K4 = 0
      CALL TRI3 (MR, vecA, vecBB, vecC, K, vecB, vecB2, vecB3, &
                 TCOS, SIZE(TCOS), vecD, vecW, vecW2, vecW3)
      vecB(:MR) = vecB(:MR) + vecB2(:MR)
      IF (NP == 3) THEN
         TCOS(1) = TWO
         CALL TRIX (1, 0, MR, vecA, vecBB, vecC, vecB, TCOS,  &
                    size(TCOS), vecD, vecW)
      ENDIF
      matQ(:MR,1) = matQ(:MR,1) + vecB(:MR)
  188 CONTINUE
      J = NLAST - JR
      vecB(:MR) = matQ(:MR,NLAST) + matQ(:MR,J)
      JM2 = NLAST - I2R
      IF (JR == 1) THEN
         matQ(:MR,NLAST) = 0.
      ELSE
         IF (NROD == 0) THEN
            matQ(:MR,NLAST) = vecP(IP+1:MR+IP)
            IP = IP - MR
         ELSE
            matQ(:MR,NLAST) = matQ(:MR,NLAST) - matQ(:MR,JM2)
         ENDIF
      ENDIF
      CALL COSGEN (KR, 1, FNUM2, HALF, TCOS,SIZE(TCOS))
      CALL COSGEN (LR, 1, FNUM2, HALF, TCOS(KR+1),SIZE(TCOS)-KR)
      CALL TRIX (KR, LR, MR, vecA, vecBB, vecC, vecB, TCOS,  &
                 size(TCOS), vecD, vecW)
      matQ(:MR,NLAST) = matQ(:MR,NLAST) + vecB(:MR)
      NLASTP = NLAST
  197 CONTINUE
      JSTEP = JR
      JR = I2R
      I2R = I2R/2
      IF (JR == 0) GO TO 210
      JSTART = 1 + JR
      KR = KR - JR
      IF (NLAST + JR <= N) THEN
         KR = KR - JR
         NLAST = NLAST + JR
         JSTOP = NLAST - JSTEP
      ELSE
         JSTOP = NLAST - JR
      ENDIF
      LR = KR - JR
      CALL COSGEN (JR, 1, HALF, ZERO, TCOS,SIZE(TCOS))
      DO J = JSTART, JSTOP, JSTEP
         JM2 = J - JR
         JP2 = J + JR
         IF (J == JR) THEN
            vecB(:MR) = matQ(:MR,J) + matQ(:MR,JP2)
         ELSE
            vecB(:MR) = matQ(:MR,J) + matQ(:MR,JM2) + matQ(:MR,JP2)
         ENDIF
         IF (JR == 1) THEN
            matQ(:MR,J) = 0.
         ELSE
            JM1 = J - I2R
            JP1 = J + I2R
            matQ(:MR,J) = HALF*(matQ(:MR,J)-matQ(:MR,JM1)-matQ(:MR,JP1))
         ENDIF
         CALL TRIX (JR, 0, MR, vecA, vecBB, vecC, vecB, TCOS,  &
                    size(TCOS), vecD, vecW)
         matQ(:MR,J) = matQ(:MR,J) + vecB(:MR)
      END DO
      NROD = 1
      IF (NLAST + I2R <= N) NROD = 0
      IF (NLASTP /= NLAST) GO TO 188
      GO TO 197
  210 CONTINUE
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
! June      2004    Version 5.0, Fortran 90 Changes
!-----------------------------------------------------------------------
      END SUBROUTINE POSTG2
