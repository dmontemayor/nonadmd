ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!>\brief
!!     PRAXIS RETURNS THE MINIMUM OF THE FUNCTION F(X,N) OF N VARIABLES
!!     USING THE PRINCIPAL AXIS METHOD.  THE GRADIENT OF THE FUNCTION IS
!!     NOT REQUIRED.
!!\details
!!     FOR A DESCRIPTION OF THE ALGORITHM, SEE CHAPTER SEVEN OF
!!     ^ALGORITHMS FOR FINDING ZEROS AND EXTREMA OF FUNCTIONS WITHOUT
!!     CALCULATING DERIVATIVES^ BY RICHARD BRENT; PUBLISHED BY PRENTICE-
!!     HALL, NEW JERSEY, 1973.  THIS IS A FORTRAN TRANSLATION OF
!!     PROCEDURE PRAXIS GIVEN IN THE ABOVE BOOK.
!! 
!!  
!!     THE PARAMETERS ARE:
!!\param     T0       IS A TOLERANCE.  PRAXIS ATTEMPTS TO RETURN PRAXIS=F(X)
!!              SUCH THAT IF X0 IS THE TRUE LOCAL MINIMUM NEAR X, THEN
!!              NORM(X-X0) < T0 + SQUAREROOT(MACHEP)*NORM(X).
!!              A REASONABLE VALUE FOR MOST PROBLEMS IS 1D-5
!!\param     MACHEP   IS THE MACHINE PRECISION, THE SMALLEST NUMBER SUCH THAT
!!              1 + MACHEP > 1.  MACHEP SHOULD BE 16.**-13 (ABOUT
!!              2.22D-16) FOR REAL*8 ARITHMETIC ON THE IBM 360,
!!              OR 2**-59 FOR DOUBLE PRECISION ARITHMETIC ON THE U1108.
!!\param     H0       IS THE MAXIMUM STEP SIZE.  H0 SHOULD BE SET TO ABOUT THE
!!              MAXIMUM DISTANCE FROM THE INITIAL GUESS TO THE MINIMUM.
!!              (IF H0 IS SET TOO LARGE OR TOO SMALL, THE INITIAL RATE OF
!!              CONVERGENCE MAY BE SLOW.)
!!\param     N        (AT LEAST TWO) IS THE NUMBER OF VARIABLES UPON WHICH
!!              THE FUNCTION DEPENDS.  IF N > 20 DIMENSION STATEMENTS
!!              BELOW MUST BE ALTERED (SEE FURTHER COMMENTS BELOW).
!!\param     PRIN     CONTROLS THE PRINTING OF INTERMEDIATE RESULTS.
!!              IF PRIN=0, ONLY FATAL ERROR MESSAGES ARE PRINTED.
!!              IF PRIN=1, F IS PRINTED AFTER EVERY N+1 OR N+2 LINEAR
!!              MINIMIZATIONS.  FINAL X IS PRINTED, BUT INTERMEDIATE X IS
!!              PRINTED ONLY IF N IS AT MOST 5.
!!              IF PRIN=2, SCALE FACTORS (IF ANY) AND PRINCIPAL VALUES
!!              OF THE APPROXIMATING QUADRATIC FORM ARE ALSO PRINTED.
!!              IF PRIN=3, X IS ALSO PRINTED AFTER EVERY FEW LINEAR
!!              MINIMIZATIONS.
!!              IF PRIN=4, THE PRINCIPAL VECTORS OF THE APPROXIMATING
!!              QUADRATIC FORM, SEARCH DIRECTIONS, AND SECOND DIFFERENCES
!!              ARE ALSO PRINTED.
!!              PRINTING IS ON UNIT 6  (SO ON 360/50 ADD A CARD
!!              //G.FT06F001 DD SYSOUT=A IF NECESSARY).
!!\param     X        IS AN ARRAY CONTAINING ON ENTRY A GUESS OF THE POINT OF
!!              MINIMUM, ON RETURN THE ESTIMATED POINT OF MINIMUM.
!!\param     F(X,N)   IS THE FUNCTION TO BE MINIMIZED.  F SHOULD BE A REAL*8
!!              FUNCTION DECLARED EXTERNAL IN THE CALLING PROGRAM.
!!
!!     THE APPROXIMATING QUADRATIC FORM IS
!!              Q(X') = F(X,N) + (1/2) * (X'-X)-TRANSPOSE * A * (X'-X)
!!     WHERE X IS THE BEST ESTIMATE OF THE MINIMUM AND A IS
!!              INVERSE(V-TRANSPOSE) * D * INVERSE(V)
!!     (V(*,*) IS THE MATRIX OF SEARCH DIRECTIONS; D(*) IS THE ARRAY
!!     OF SECOND DIFFERENCES). IF F HAS CONTINUOUS SECOND DERIVATIVES
!!     NEAR X0, A WILL TEND TO THE HESSIAN OF F AT X0 AS X APPROACHES X0.
!!
!!     IT IS ASSUMED THAT ON FLOATING-POINT UNDERFLOW THE RESULT IS SET
!!     TO ZERO.
!!     THE USER SHOULD OBSERVE THE COMMENT ON HEURISTIC NUMBERS AFTER
!!     THE INITIALIZATION OF MACHINE DEPENDENT NUMBERS.
!! \authors
!! Adapted from Coker Lab Codes for NonAdMD by Daniel Montemayor
!! \date
!! Nov 2012
!<
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DOUBLE PRECISION FUNCTION PRAXIS(T0,MACHEP,H0,N,PRIN,X,F,idum)
      IMPLICIT REAL*8 (A-H, O-Z)
      REAL*8 T0,MACHEP,H0,X(N),F
      integer idum
      EXTERNAL F,min
      INTEGER PRIN
!!
      LOGICAL ILLC
      INTEGER NL,NF,KL,KT,KTM
      REAL*8 S,SL,DN,DMIN,FX,F1,LDS,LDT,T,H,SF,DF,QF1,QD0,QD1,QA,QB,QC
      REAL*8 M2,M4,SMALL,VSMALL,LARGE,VLARGE,SCBD,LDFAC,T2,DNI,VALUE
!!      REAL*8 ran0
!!
!!.....IF N>20 (OR N<20 AND YOU NEED MORE SPACE) CHANGE '20' TO THE
!!     LARGEST VALUE OF N IN THE NEXT CARD, IN THE CARD 'IDIM=20', AND
!!     IN THE DIMENSION STATEMENTS IN SUBROUTINES MINFIT,MIN,FLIN,QUAD.
!!
      REAL*8 D(30),Y(30),Z(30),Q0(30),Q1(30),V(30,30)
      COMMON /GLOBAL/ FX,LDT,DMIN,NF,NL
     .       /Q/ V,Q0,Q1,QA,QB,QC,QD0,QD1,QF1
      external ran0
      IDIM=30
      IF (N.LE.IDIM) GO TO 5
      WRITE (6, 3) N
3     FORMAT (' N =', I4, ' IS TOO LARGE - CHANGE DIMENSION STATEMENTS')
      STOP
!!
!!.....INITIALIZATION.....
!!     MACHINE DEPENDENT NUMBERS:
!!
5     SMALL=MACHEP*MACHEP
      VSMALL=SMALL*SMALL
      LARGE=1.D0/SMALL
      VLARGE=1.D0/VSMALL
      M2=DSQRT(MACHEP)
      M4=DSQRT(M2)
!!
!!     HEURISTIC NUMBERS:
!!     IF THE AXES MAY BE BADLY SCALED (WHICH IS TO BE AVOIDED IF
!!     POSSIBLE), THEN SET SCBD=10.  OTHERWISE SET SCBD=1.
!!     IF THE PROBLEM IS KNOWN TO BE ILL-CONDITIONED, SET ILLC=TRUE.
!!     OTHERWISE SET ILLC=FALSE.
!!     KTM IS THE NUMBER OF ITERATIONS WITHOUT IMPROVEMENT BEFORE THE
!!     ALGORITHM TERMINATES.  KTM=4 IS VERY CAUTIOUS; USUALLY KTM=1
!!     IS SATISFACTORY.
!!
      SCBD=1.D0
      ILLC=.FALSE.
      KTM=1
!!
      LDFAC=0.01D0
      IF (ILLC) LDFAC=0.1D0
      KT=0
      NL=0
      NF=1
      FX=F(X,N)
      QF1=FX
      T=SMALL+ABS(T0)
      T2=T
      DMIN=SMALL
      H=H0
      IF (H.LT.100D0*T) H = 100D0*T
      LDT=H
!!.....THE FIRST SET OF SEARCH DIRECTIONS V IS THE IDENTITY MATRIX.....
      DO 20 I=1,N
           DO 10 J=1,N
10              V(I,J)=0.0D0
20         V(I,I)=1.D0
      D(1)=0.D0
      QD0=0.D0
      DO 30 I=1,N
           Q0(I)=X(I)
30         Q1(I)=X(I)
      IF (PRIN.GT.0) CALL PRINT(N,X,PRIN)
!!
!!.....THE MAIN LOOP STARTS HERE.....
40    SF=D(1)
      D(1)=0.D0
      S=0.D0
!!
!!.....MINIMIZE ALONG THE FIRST DIRECTION V(*,1).
!!     FX MUST BE PASSED TO MIN BY VALUE.
      VALUE=FX
      CALL MINP(N,1,2,D(1),S,VALUE,.FALSE.,F,X,T,MACHEP,H)
      IF (S.GT.0.D0) GO TO 50
           DO 45 I=1,N
45              V(I,1)=-V(I,1)
50    IF (SF.GT.0.9D0*D(1).AND.0.9D0*SF.LT.D(1)) GO TO 70
           DO 60 I=2,N
60              D(I)=0.D0
!!
!!.....THE INNER LOOP STARTS HERE.....
70    DO 170 K=2,N
           DO 75 I=1,N
75              Y(I)=X(I)
           SF=FX
           IF (KT.GT.0) ILLC=.TRUE.
80         KL=K
           DF=0.D0
!!
!!.....A RANDOM STEP FOLLOWS (TO AVOID RESOLUTION VALLEYS).
!!     PRAXIS ASSUMES THAT DRAND RETURNS A RANDOM NUMBER UNIFORMLY
!!     DISTRIBUTED IN (0,1).
!!
           IF(.NOT.ILLC) GO TO 95
                DO 90 I=1,N
                     S=(0.1D0*LDT+T2*(10**KT))*(ran0(idum)-0.5D0)
                     Z(I)=S
                     DO 85 J=1,N
85                        X(J)=X(J)+S*V(J,I)
90              CONTINUE
                FX=F(X,N)
                NF=NF+1
!!
!!.....MINIMIZE ALONG THE ^NON-CONJUGATE^ DIRECTIONS V(*,K),...,V(*,N)
!!
95         DO 105 K2=K,N
                SL=FX
                S=0.D0
                VALUE=FX
                CALL MINP(N,K2,2,D(K2),S,VALUE,.FALSE.,F,X,T,MACHEP,H)
                IF (ILLC) GO TO 97
                     S=SL-FX
                     GO TO 99
97              S=D(K2)*((S+Z(K2))**2)
99              IF (DF.GT.S) GO TO 105
                     DF=S
                     KL=K2
105        CONTINUE
           IF (ILLC.OR.(DF.GE.ABS((100.D0*MACHEP)*FX))) GO TO 110
!!
!!.....IF THERE WAS NOT MUCH IMPROVEMENT ON THE FIRST TRY, SET
!!     ILLC=TRUE AND START THE INNER LOOP AGAIN.....
!!
           ILLC=.TRUE.
           GO TO 80
110        IF (K.EQ.2.AND.PRIN.GT.3) CALL VCPRNT(1,D,N)
!!
!!.....MINIMIZE ALONG THE ^CONJUGATE^ DIRECTIONS V(*,1),...,V(*,K-1)
!!
           KM1=K-1
           DO 120 K2=1,KM1
           S=0
           VALUE=FX
           CALL MINP(N,K2,2,D(K2),S,VALUE,.FALSE.,F,X,T,MACHEP,H)
120        CONTINUE
           F1=FX
           FX=SF
           LDS=0
           DO 130 I=1,N
                SL=X(I)
                X(I)=Y(I)
                SL=SL-Y(I)
                Y(I)=SL
130             LDS=LDS+SL*SL
           LDS=DSQRT(LDS)
           IF (LDS.LE.SMALL) GO TO 160
!!
!!.....DISCARD DIRECTION V(*,KL).
!!     IF NO RANDOM STEP WAS TAKEN, V(*,KL) IS THE ^NON-CONJUGATE^
!!     DIRECTION ALONG WHICH THE GREATEST IMPROVEMENT WAS MADE.....
!!
           KLMK=KL-K
           IF (KLMK.LT.1) GO TO 141
           DO 140 II=1,KLMK
                I=KL-II
                DO 135 J=1,N
135                  V(J,I+1)=V(J,I)
140             D(I+1)=D(I)
141        D(K)=0
           DO 145 I=1,N
145             V(I,K)=Y(I)/LDS
!!
!!.....MINIMIZE ALONG THE NEW ^CONJUGATE^ DIRECTION V(*,K), WHICH IS
!!     THE NORMALIZED VECTOR:  (NEW X) - (OLD X).....
!!
           VALUE=F1
           CALL MINP(N,K,4,D(K),LDS,VALUE,.TRUE.,F,X,T,MACHEP,H)
           IF (LDS.GT.0.D0) GO TO 160
                LDS=-LDS
                DO 150 I=1,N
150                  V(I,K)=-V(I,K)
160        LDT=LDFAC*LDT
           IF (LDT.LT.LDS) LDT=LDS
           IF (PRIN.GT.0) CALL PRINT(N,X,PRIN)
           T2=0.D0
           DO 165 I=1,N
165             T2=T2+X(I)**2
           T2=M2*DSQRT(T2)+T
!!
!!.....SEE WHETHER THE LENGTH OF THE STEP TAKEN SINCE STARTING THE
!!     INNER LOOP EXCEEDS HALF THE TOLERANCE.....
!!
           IF (LDT.GT.(0.5D0*T2)) KT = -1
           KT=KT+1
           IF (KT.GT.KTM) GO TO 400
170   CONTINUE
!!.....THE INNER LOOP ENDS HERE.
!!
!!     TRY QUADRATIC EXTRAPOLATION IN CASE WE ARE IN A CURVED VALLEY.
!!
171   CALL QUAD(N,F,X,T,MACHEP,H)
      DN=0.D0
      DO 175 I=1,N
           D(I)=1.D0/DSQRT(D(I))
           IF (DN.LT.D(I)) DN=D(I)
175   CONTINUE
      IF (PRIN.GT.3) CALL MAPRNT(1,V,IDIM,N)
      DO 180 J=1,N
           S=D(J)/DN
           DO 180 I=1,N
180             V(I,J)=S*V(I,J)
!!
!!.....SCALE THE AXES TO TRY TO REDUCE THE CONDITION NUMBER.....
!!
      IF (SCBD.LE.1.D0) GO TO 200
           S=VLARGE
           DO 185 I=1,N
                SL=0.D0
                DO 182 J=1,N
182                  SL=SL+V(I,J)*V(I,J)
                Z(I)=DSQRT(SL)
                IF (Z(I).LT.M4) Z(I)=M4
                IF (S.GT.Z(I)) S=Z(I)
185        CONTINUE
           DO 195 I=1,N
                SL=S/Z(I)
                Z(I)=1.D0/SL
                IF (Z(I).LE.SCBD) GO TO 189
                     SL=1.D0/SCBD
                     Z(I)=SCBD
189             DO 190 J=1,N
190                  V(I,J)=SL*V(I,J)
195        CONTINUE
!!
!!.....CALCULATE A NEW SET OF ORTHOGONAL DIRECTIONS BEFORE REPEATING
!!     THE MAIN LOOP.
!!     FIRST TRANSPOSE V FOR MINFIT:
!!
200   DO 220 I=2,N
           IM1=I-1
           DO 210 J=1,IM1
                S=V(I,J)
                V(I,J)=V(J,I)
210             V(J,I)=S
220   CONTINUE
!!
!!.....CALL MINFIT TO FIND THE SINGULAR VALUE DECOMPOSITION OF V.
!!     THIS GIVES THE PRINCIPAL VALUES AND PRINCIPAL DIRECTIONS OF THE
!!     APPROXIMATING QUADRATIC FORM WITHOUT SQUARING THE CONDITION
!!     NUMBER.....
!!
      CALL MINFIT(IDIM,N,MACHEP,VSMALL,V,D)
!!
!!.....UNSCALE THE AXES.....
!!
      IF (SCBD.LE.1.D0) GO TO 250
           DO 230 I=1,N
                S=Z(I)
                DO 225 J=1,N
225                  V(I,J)=S*V(I,J)
230        CONTINUE
           DO 245 I=1,N
                S=0.D0
                DO 235 J=1,N
235                  S=S+V(J,I)**2
                S=DSQRT(S)
                D(I)=S*D(I)
                S=1/S
                DO 240 J=1,N
240                  V(J,I)=S*V(J,I)
245        CONTINUE
!!
250   DO 270 I=1,N
           DNI=DN*D(I)
           IF (DNI.GT.LARGE) GO TO 265
                IF (DNI.LT.SMALL) GO TO 260
                     D(I)=1/(DNI*DNI)
                     GO TO 270
260             D(I)=VLARGE
                GO TO 270
265        D(I)=VSMALL
270   CONTINUE
!!
!!.....SORT THE EIGENVALUES AND EIGENVECTORS.....
!!
      CALL SORT(IDIM,N,D,V)
      DMIN=D(N)
      IF (DMIN.LT.SMALL) DMIN=SMALL
      ILLC=.FALSE.
      IF (M2*D(1).GT.DMIN) ILLC=.TRUE.
      IF (PRIN.GT.1.AND.SCBD.GT.1.D0) CALL VCPRNT(2,Z,N)
      IF (PRIN.GT.1) CALL VCPRNT(3,D,N)
      IF (PRIN.GT.3) CALL MAPRNT(2,V,IDIM,N)
!!.....THE MAIN LOOP ENDS HERE.....
!!
      GO TO 40
!!
!!.....RETURN.....
!!
400   IF ((PRIN.GT.0).AND.(PRIN.LE.2).AND.(N.GT.5)) CALL VCPRNT(4,X,N)
      PRAXIS=FX
      RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> \brief QUAD LOOKS FOR THE MINIMUM OF F ALONG A CURVE DEFINED BY Q0,Q1,X
!<ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE QUAD(N,F,X,T,MACHEP,H)
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL F,MIN
      REAL*8 X(N),MACHEP,LDT,L
      DIMENSION V(30,30),Q0(30),Q1(30)
      COMMON /GLOBAL/ FX,LDT,DMIN,NF,NL
     .       /Q/ V,Q0,Q1,QA,QB,QC,QD0,QD1,QF1
      S = FX
      FX = QF1
      QF1 = S
      QD1 = 0.D0
      DO 1 I=1,N
         S = X(I)
         L = Q1(I)
         X(I) = L
         Q1(I) = S
1        QD1 = QD1 + (S-L)**2
      QD1 = DSQRT(QD1)
      L = QD1
      S = 0.D0
      IF (QD0. LE. 0.D0 .OR. QD1 .LE. 0.D0 .OR. NL .LT. 3*N*N) GO TO 2
      VALUE=QF1
      CALL MINP(N,0,2,S,L,VALUE,.TRUE.,F,X,T,MACHEP,H)
      QA = (L*(L-QD1))/(QD0*(QD0+QD1))
      QB = ((L+QD0)*(QD1-L))/(QD0*QD1)
      QC = (L*(L+QD0))/(QD1*(QD0+QD1))
      GO TO 3
2     FX = QF1
      QA = 0.D0
      QB = QA
      QC = 1.D0
3     QD0 = QD1
      DO 4 I=1,N
         S = Q0(I)
         Q0(I) = X(I)
4        X(I) = (QA*S + QB*X(I)) + QC*Q1(I)
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> \brief
!!   AN IMPROVED VERSION OF MINFIT (SEE GOLUB AND REINSCH, 1969)
!!   RESTRICTED TO SQUARE AB.
!! \details
!!   THE SINGULAR VALUES OF THE ARRAY AB ARE RETURNED IN Q AND AB IS
!!   OVERWRITTEN WITH THE ORTHOGONAL MATRIX V SUCH THAT U.DIAG(Q) = AB.V,
!!   WHERE U IS ANOTHER ORTHOGONAL MATRIX.
!<cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE MINFIT(M,N,MACHEP,TOL,AB,Q)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 MACHEP
      DIMENSION AB(M,M),Q(M)
      DIMENSION E(30)
!!...HOUSEHOLDER'S REDUCTION TO BIDIAGONAL FORM...
      IF (N.EQ.1) GO TO 200
      EPS = MACHEP
      G = 0.D0
      X = 0.D0
      DO 11 I=1,N
         E(I) = G
         S = 0.D0
         L = I + 1
         DO 1 J=I,N
1           S = S + AB(J,I)**2
         G = 0.D0
         IF (S.LT.TOL) GO TO 4
            F = AB(I,I)
            G = DSQRT(S)
            IF (F.GE.0.D0) G = -G
            H = F*G - S
            AB(I,I)=F-G
            IF (L.GT.N) GO TO 4
            DO 3 J=L,N
               F = 0.D0
               DO 2 K=I,N
2                 F = F + AB(K,I)*AB(K,J)
               F = F/H
               DO 3 K=I,N
3                 AB(K,J) = AB(K,J) + F*AB(K,I)
4        Q(I) = G
         S = 0.D0
         IF (I.EQ.N) GO TO 6
         DO 5 J=L,N
5           S = S + AB(I,J)*AB(I,J)
6        G = 0.D0
         IF (S.LT.TOL) GO TO 10
            IF (I.EQ.N) GO TO 16
            F = AB(I,I+1)
16          G = DSQRT(S)
            IF (F.GE.0.D0) G = -G
            H = F*G - S
            IF (I.EQ.N) GO TO 10
            AB(I,I+1) = F - G
            DO 7 J=L,N
7              E(J) = AB(I,J)/H
            DO 9 J=L,N
               S = 0.D0
               DO 8 K=L,N
8                 S = S + AB(J,K)*AB(I,K)
               DO 9 K=L,N
9                 AB(J,K) = AB(J,K) + S*E(K)
10       Y = DABS(Q(I)) + DABS(E(I))
11       IF (Y.GT.X) X = Y
!!...ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS...
      AB(N,N) = 1.D0
      G = E(N)
      L = N
      DO 25 II=2,N
         I = N - II + 1
         IF (G.EQ.0.D0) GO TO 23
         H = AB(I,I+1)*G
         DO 20 J=L,N
20          AB(J,I) = AB(I,J)/H
         DO 22 J=L,N
            S = 0.D0
            DO 21 K=L,N
21             S = S + AB(I,K)*AB(K,J)
            DO 22 K=L,N
22             AB(K,J) = AB(K,J) + S*AB(K,I)
23       DO 24 J=L,N
            AB(I,J) = 0.D0
24          AB(J,I) = 0.D0
         AB(I,I) = 1.D0
         G = E(I)
25       L = I
!!...DIAGONALIZATION OF THE BIDIAGONAL FORM...
100   EPS = EPS*X
      DO 150 KK=1,N
         K = N - KK + 1
         KT = 0
101      KT = KT + 1
         IF (KT.LE.30) GO TO 102
            E(K) = 0.D0
      CALL MAPRNT (3, AB, M, N)
102      DO 103 LL2=1,K
            L2 = K - LL2 + 1
            L = L2
            IF (DABS(E(L)).LE.EPS) GO TO 120
            IF (L.EQ.1) GO TO 103
            IF (DABS(Q(L-1)).LE.EPS) GO TO 110
103         CONTINUE
!!...CANCELLATION OF E(L) IF L>1...
110      C = 0.D0
         S = 1.D0
         DO 116 I=L,K
            F = S*E(I)
            E(I) = C*E(I)
            IF (DABS(F).LE.EPS) GO TO 120
            G = Q(I)
!!...Q(I) = H = DSQRT(G*G + F*F)...
            IF (DABS(F).LT.DABS(G)) GO TO 113
            IF (F) 112,111,112
111         H = 0.D0
            GO TO 114
112         H = DABS(F)*DSQRT(1.D0 + (G/F)**2)
            GO TO 114
113         H = DABS(G)*DSQRT(1.D0 + (F/G)**2)
114         Q(I) = H
            IF (H.NE.0.D0) GO TO 115
               G = 1.D0
               H = 1.D0
115         C = G/H
116         S = -F/H
!!...TEST FOR CONVERGENCE...
120      Z = Q(K)
         IF (L.EQ.K) GO TO 140
!!...SHIFT FROM BOTTOM 2*2 MINOR...
         X = Q(L)
         Y = Q(K-1)
         G = E(K-1)
         H = E(K)
         F = ((Y - Z)*(Y + Z) + (G - H)*(G + H))/(2*H*Y)
         G = DSQRT(F*F + 1.0D0)
         TEMP = F - G
         IF (F.GE.0.D0) TEMP = F + G
         F = ((X - Z)*(X + Z) + H*(Y/TEMP - H))/X
!!...NEXT QR TRANSFORMATION...
         C = 1.D0
         S = 1.D0
         LP1 = L + 1
         IF (LP1.GT.K) GO TO 133
         DO 132 I=LP1,K
            G = E(I)
            Y = Q(I)
            H = S*G
            G = G*C
            IF (DABS(F).LT.DABS(H)) GO TO 123
            IF (F) 122,121,122
121         Z = 0.D0
            GO TO 124
122         Z = DABS(F)*DSQRT(1.D0 + (H/F)**2)
            GO TO 124
123         Z = DABS(H)*DSQRT(1.D0 + (F/H)**2)
124         E(I-1) = Z
            IF (Z.NE.0.D0) GO TO 125
               F = 1.D0
               Z = 1.D0
125         C = F/Z
            S = H/Z
            F = X*C + G*S
            G = -X*S + G*C
            H = Y*S
            Y = Y*C
            DO 126 J=1,N
               X = AB(J,I-1)
               Z = AB(J,I)
               AB(J,I-1) = X*C + Z*S
126            AB(J,I) = -X*S + Z*C
            IF (DABS(F).LT.DABS(H)) GO TO 129
            IF (F) 128,127,128
127         Z = 0.D0
            GO TO 130
128         Z = DABS(F)*DSQRT(1.D0 + (H/F)**2)
            GO TO 130
129         Z = DABS(H)*DSQRT(1.D0 + (F/H)**2)
130         Q(I-1) = Z
            IF (Z.NE.0.D0) GO TO 131
               F = 1.D0
               Z = 1.D0
131         C = F/Z
            S = H/Z
            F = C*G + S*Y
132         X = -S*G + C*Y
133      E(L) = 0.D0
         E(K) = F
         Q(K) = X
         GO TO 101
!!...CONVERGENCE:  Q(K) IS MADE NON-NEGATIVE...
140      IF (Z.GE.0.D0) GO TO 150
         Q(K) = -Z
         DO 141 J=1,N
141         AB(J,K) = -AB(J,K)
150      CONTINUE
      RETURN
200   Q(1) = AB(1,1)
      AB(1,1) = 1.D0
      RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!> \brief
!!   FLIN IS THE FUNCTION OF ONE REAL VARIABLE L THAT IS MINIMIZED
!!   BY THE SUBROUTINE MIN
!<ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DOUBLE PRECISION FUNCTION FLIN (N,J,L,F,X,NF)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 L,X(N)
      EXTERNAL F
      DIMENSION V(30,30),Q0(30),Q1(30)
      COMMON /Q/ V,Q0,Q1,QA,QB,QC,QD0,QD1,QF1
      DIMENSION T(30)
      IF (J .EQ. 0) GO TO 2
!!...THE SEARCH IS LINEAR...
      DO 1 I=1,N
1        T(I) = X(I) + L*V(I,J)
      GO TO 4
!!...THE SEARCH IS ALONG A PARABOLIC SPACE CURVE...
2     QA = (L*(L - QD1))/(QD0*(QD0 + QD1))
      QB = ((L + QD0)*(QD1 - L))/(QD0*QD1)
      QC = (L*(L + QD0))/(QD1*(QD0 + QD1))
      DO 3 I=1,N
3        T(I) = (QA*Q0(I) + QB*X(I)) + QC*Q1(I)
!!...THE FUNCTION EVALUATION COUNTER NF IS INCREMENTED...
4     NF = NF + 1
      FLIN = F(T,N)
      RETURN
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!>\details
!!   THE SUBROUTINE MIN MINIMIZES F FROM X IN THE DIRECTION V(*,J) UNLESS
!!   J IS LESS THAN 1, WHEN A QUADRATIC SEARCH IS MADE IN THE PLANE
!!   DEFINED BY Q0,Q1,X.
!!   D2 IS EITHER ZERO OR AN APPROXIMATION TO HALF F^.
!!   ON ENTRY, X1 IS AN ESTIMATE OF THE DISTANCE FROM X TO THE MINIMUM
!!   ALONG V(*,J) (OR, IF J=0, A CURVE). ON RETURN, X1 IS THE DISTANCE
!!   FOUND.
!!   IF FK=.TRUE., THEN F1 IS FLIN(X1).  OTHERWISE X1 AND F1 ARE IGNORED
!!   ON ENTRY UNLESS FINAL FX IS GREATER THAN F1.
!!   NITS CONTROLS THE NUMBER OF TIMES AN ATTEMPT WILL BE MADE TO HALVE
!!   THE INTERVAL.
!<
      SUBROUTINE MINP(N,J,NITS,D2,X1,F1,FK,F,X,T,MACHEP,H)
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL F
      LOGICAL FK
      REAL*8 MACHEP,X(N),LDT
      DIMENSION V(30,30),Q0(30),Q1(30)
      COMMON /GLOBAL/ FX,LDT,DMIN,NF,NL
     .       /Q/ V,Q0,Q1,QA,QB,QC,QD0,QD1,QF1
!!
      LOGICAL DZ
      REAL*8 M2,M4
      SMALL = MACHEP**2
      M2 = DSQRT(MACHEP)
      M4 = DSQRT(M2)
      SF1 = F1
      SX1 = X1
      K = 0
      XM = 0.D0
      FM = FX
      F0 = FX
      DZ = D2.LT.MACHEP
!!...FIND THE STEP SIZE...
      S = 0.D0
      DO 1 I=1,N
1        S = S + X(I)**2
      S = DSQRT(S)
      TEMP = D2
      IF (DZ) TEMP = DMIN
      T2 = M4*DSQRT(DABS(FX)/TEMP + S*LDT) + M2*LDT
      S = M4*S + T
      IF (DZ.AND.T2.GT.S) T2 = S
      T2 = DMAX1(T2,SMALL)
      T2 = DMIN1(T2,.01D0*H)
      IF (.NOT.FK.OR.F1.GT.FM) GO TO 2
      XM = X1
      FM = F1
2     IF (FK.AND.DABS(X1).GE.T2) GO TO 3
      TEMP=1.D0
      IF (X1.LT.0.D0) TEMP=-1.D0
      X1=TEMP*T2
      F1 = FLIN(N,J,X1,F,X,NF)
3     IF (F1.GT.FM) GO TO 4
      XM = X1
      FM = F1
4     IF (.NOT.DZ) GO TO 6
!!...EVALUATE FLIN AT ANOTHER POINT AND ESTIMATE THE SECOND DERIVATIVE...
      X2 = -X1
      IF (F0.GE.F1) X2 = 2.D0*X1
      F2 = FLIN(N,J,X2,F,X,NF)
      IF (F2.GT.FM) GO TO 5
         XM = X2
         FM = F2
5     D2 = (X2*(F1 - F0)-X1*(F2 - F0))/((X1*X2)*(X1 - X2))
!!...ESTIMATE THE FIRST DERIVATIVE AT 0...
6     D1 = (F1 - F0)/X1 - X1*D2
      DZ = .TRUE.
!!...PREDICT THE MINIMUM...
      IF (D2.GT.SMALL) GO TO 7
         X2 = H
         IF (D1.GE.0.D0) X2 = -X2
         GO TO 8
7        X2 = (-.5D0*D1)/D2
8     IF (DABS(X2).LE.H) GO TO 11
         IF (X2) 9,9,10
9        X2 = -H
         GO TO 11
10       X2 = H
!!...EVALUATE F AT THE PREDICTED MINIMUM...
11    F2 = FLIN(N,J,X2,F,X,NF)
      IF (K.GE.NITS.OR.F2.LE.F0) GO TO 12
!!...NO SUCCESS, SO TRY AGAIN...
         K = K + 1
         IF (F0.LT.F1.AND.(X1*X2).GT.0.D0) GO TO 4
         X2 = 0.5D0*X2
          GO TO 11
!!...INCREMENT THE ONE-DIMENSIONAL SEARCH COUNTER...
12    NL = NL + 1
      IF (F2.LE.FM) GO TO 13
      X2 = XM
      GO TO 14
13    FM = F2
!!...GET A NEW ESTIMATE OF THE SECOND DERIVATIVE...
14    IF (DABS(X2*(X2 - X1)).LE.SMALL) GO TO 15
         D2 = (X2*(F1-F0) - X1*(FM-F0))/((X1*X2)*(X1 - X2))
         GO TO 16
15       IF (K.GT.0) D2 = 0.D0
16    IF (D2.LE.SMALL) D2 = SMALL
      X1 = X2
      FX = FM
      IF (SF1.GE.FX) GO TO 17
         FX = SF1
         X1 = SX1
!!...UPDATE X FOR LINEAR BUT NOT PARABOLIC SEARCH...
17    IF (J.EQ.0) RETURN
      DO 18 I=1,N
18       X(I) = X(I) + X1*V(I,J)
      RETURN
      END
!!-----------------------------------------------------------------------
!!
       subroutine print(n,x,prin)
       implicit real*8(a-h,o-z)
       integer prin
       real*8 x(n)
       write(6,*)'                         *** Values of Parameters ***'
       do 1 i=1,n,4
       write(6,10) (j,x(j),j=i,min(i+3,n))
  1    continue
  10   format(4(' x(',i2,')=',d13.6))
       return
       end
!!-----------------------------------------------------------------------
!!
       subroutine vcprnt(iflag,x,n)
       implicit real*8(a-h,o-z)
       real*8 x(n)
       if(iflag.eq.1) write(6,100)
       if(iflag.eq.2) write(6,200)
       if(iflag.eq.3) write(6,300)
       if(iflag.eq.4) write(6,400)
  100  format(20x,'    *** The Second Difference Array D is ***')
  200  format(20x,'         *** The Scale Factors are ***')
  300  format(20x,' *** The Approx. Quad. Form has Princ. Vals. ***')
  400  format(20x,'                   *** X is ***')
       do 1 i=1,n,4
       if(iflag.eq.1) write(6,10) (j,x(j),j=i,min(i+3,n))
       if(iflag.eq.2) write(6,20) (j,x(j),j=i,min(i+3,n))
       if(iflag.eq.3) write(6,30) (j,x(j),j=i,min(i+3,n))
       if(iflag.eq.4) write(6,40) (j,x(j),j=i,min(i+3,n))
  1    continue
  10   format(4(' d(',i2,')=',d13.6))
  20   format(4(' z(',i2,')=',d13.6))
  30   format(4(' d(',i2,')=',d13.6))
  40   format(4(' x(',i2,')=',d13.6))
       return
       end
!!-----------------------------------------------------------------------
!!
       subroutine maprnt(iflag,v,idim,n)
       implicit real*8(a-h,o-z)
       real*8 v(idim,idim)
       if(iflag.eq.1) write(6,100)
       if(iflag.eq.2) write(6,200)
       if(iflag.eq.3) write(6,300)
  100  format(20x,'     *** The New Directions are ***')
  200  format(20x,'     *** The Principle Axes are ***')
  300  format(20x,'        *** The AB Matrix is ***')
       do 5 k=1,n
       write(6,*)' Col. No.',k
       do 1 i=1,n,3
       if(iflag.eq.1) write(6,10) (j,k,v(j,k),j=i,min(i+2,n))
       if(iflag.eq.2) write(6,20) (j,k,v(j,k),j=i,min(i+2,n))
       if(iflag.eq.3) write(6,30) (j,k,v(j,k),j=i,min(i+2,n))
  1    continue
  5    continue
  10   format(3(' v(',i2,',',i2,')=',d13.6))
  20   format(3(' v(',i2,',',i2,')=',d13.6))
  30   format(3(' a(',i2,',',i2,')=',d13.6))
       return
       end
!!---------------------------------------------------------------------
!!
       subroutine sort(idim,n,d,v)
       implicit real*8 (a-h,o-z)
       parameter(id=30)
       real*8 d(idim),v(idim,idim),vhold(id)
       nleft=n
  6    continue
       imin=1
       dmin=d(1)
       do 1 i=2,nleft
       if(d(i).gt.dmin) goto 1
       imin=i
       dmin=d(i)
  1    continue
       do 2 i=1,n
       vhold(i)=v(i,imin)
  2    continue
       do 3 i=imin+1,nleft
       d(i-1)=d(i)
          do 4 j=1,n
          v(j,i-1)=v(j,i)
  4       continue
  3    continue
       d(nleft)=dmin
       do 5 j=1,n
       v(j,nleft)=vhold(j)
  5    continue
       nleft=nleft-1
       if(nleft.gt.1) goto 6
       return
       end
!!ccccccccccccccccccccccccccccccccccccccccc
!> \brief
!! Praxis' random number generator
!! \remarks
!! Has same name as uniform random number generator in rand_class
!! \todo
!! Use the standard rand_class to avoid conflicts
!<ccccccccccccccccccccccccccccccccccccccccc
      DOUBLE PRECISION FUNCTION  ran0(idum)
	implicit none
	INTEGER idum,IA,IM,IQ,IR,MASK
!!	real*8 ran0,AM
        real*8  AM
	PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     1    MASK=123459876)
	INTEGER k
	idum=ieor(idum,MASK)
	k=idum/IQ
	idum=IA*(idum-k*IQ)-IR*k
	if (idum.lt.0) idum=idum+IM
	ran0=AM*idum
	idum=ieor(idum,MASK)
	return
	END
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
