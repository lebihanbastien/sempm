C**********************************************************************
C UTILITY ROUTINES FOR COMPUTING THE UNSTABLE, STABLE AND CENTRAL
C PART ASSOCIATED WITH A PRODUCT OF MATRICES.
C
C   DIRH  Computation of directions, dilatations and small test.
C  VEPRO  Computation of directions.
C DIAHIP  Computation of invar. subspace associated with central part.
C POWERS  Power method for a product of matrices.
C
C AUXILIARY ROUTINES: ORTO, MATINV, FACTOR, SUBST, VNORM, GET, PUT,
C                     INVERS
C**********************************************************************
       	 subroutine fortfunc(ii,ff)
         integer ii
         real*4 ff

         write(6,100) ii, ff
100 	 format ('ii=', i2, ' ff =', f6.3)

	     return
	     end

       SUBROUTINE DIRHC (DAT,VECS,DILAT,NIT,MAXIN)
C**********************************************************************
C  IT COMPUTES THE UNSTABLE, STABLE AND "CENTRAL" DIRECTIONS OF THE
C  6*6 MATRIX A=DAT(*,*,NIT)*DAT(*,*,NIT-1)*...*DAT(*,*,1).
C  IT ALSO COMPUTES THE "DILATATIONS" OF THESE DIRECTIONS WHEN THEY
C  "CROSS" THE MATRICES.
C  ALSO, A TEST ABOUT THE "CONTINUITY" OF THE DIRECTIONS CAN BE DONE
C  DECOMENTING FEW LINES.
C
C   INPUT:
C
C     DAT(i,j,k)   i=1,6  j=1,6  k=1,NIT  CONTAIN THE M 6*6 MATRICES
C           NIT    NUMBER OF MATRICES STORED IN DAT.
C         MAXIN    DIMENSION OF THE LAST COMPONENT OF SOME MATRICES.
C
C   OUTPUT:
C
C    VECS(i,j,k)  i=1,6  j=1,6  k=0,NIT.  DIRECTIONS (MODULUS 1) AT
C                 THE BEGINING OF EACH INTERVAL DEFINED BY DAT(*,*,k).
C    DILAT(j,k)   j=1,6  k=1,NIT. DILATATIONS OF THE ABOVE VECTORS
C                 WHEN "CROSSING" THE MATRIX DAT(*,*,k). THIS IS
C                  DILAT(j,k)*VECS(*,j,k)=DAT(*,*,k)*VECS(*,j,k-1)
C                 FOR j=1,6 AND MATRIX k=1,NIT.
C                 THE DIRECTIONS j=1 and 2 REFER TO THE UNSTABLE AND
C                 STABLE RESPECTIVELY, AND j=3,6 TO THE CENTRAL ONES.
C**********************************************************************
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION DAT(6,6,MAXIN),VECS(6,6,0:MAXIN),A(6,6),BUX(6)
       DIMENSION CUX(6),ACM(6),DILAT(6,MAXIN)
C COMPUTING THE DIRECTIONS AND THEN TESTING THE CONTINUITY OF THE
C RESULTS AT THE SAME TIME THAT MATRIX DILAT IS FILLED.
       CALL VEPRO(DAT,NIT,VECS,MAXIN)
       DO 2 N=1,6
C       WRITE (*,*) 'DIRHC. VECTOR NUMBER: ',N
       ACUM=1.D0
       IF (N.LE.2) ACUM=0.D0
       DO 10 IK=1,NIT
       DO 3 I=1,6
3      BUX(I)=VECS(I,N,IK-1)
C       WRITE (6,*) 'DIRHC. VECTOR ',N,' STEP: ',IK
       VNE=VNORM(BUX,6)
       CALL GET(DAT,A,MAXIN,6,IK)
       CALL MATVEC(6,A,BUX,CUX,0)
       VNS=VNORM(CUX,6)
C       WRITE (6,*) 'DIRHC. INPUT NORM AND OUTPUT NORM: ',VNE,VNS
       DILAT(N,IK)=VNS/VNE
       DO 5 I=1,6
       BUX(I)=VECS(I,N,IK)-CUX(I)/VNS
5      CONTINUE
       IF (N.LE.2) ACUM=ACUM+DLOG10(VNS)
       IF (N.GE.3) ACUM=ACUM*VNS
       VND=VNORM(BUX,6)
C       WRITE (*,*) 'DIRHC. ERROR BETWEEN IN-OUT VECTORS: ',VND
C       WRITE (*,*) '        ACCUMULATED DILATATION: ',ACUM
10     CONTINUE
       ACM(N)=ACUM
2      CONTINUE
       WRITE (6,*) 'DIRHC. ACCUMULATED DILATATIONS: '
       WRITE (6,100) ACM(1),ACM(2)
       WRITE (6,105) (ACM(N),N=3,6)
100    FORMAT (1X,'   COMPONENTS 1 & 2 (LOG10):',6F13.5)
105    FORMAT (1X,'   COMPONENTS 3 to 6:       ',4D13.5)
       END


      SUBROUTINE VEPRO(DAT,M,VECS,NMX)
C**********************************************************************
C  IT COMPUTES THE UNSTABLE, STABLE AND "CENTRAL" DIRECTIONS OF THE
C  6*6 MATRIX A=DAT(*,*,M)*DAT(*,*,M-1)*...*DAT(*,*,1).
C
C   INPUT:
C
C     DAT(i,j,k)   i=1,6  j=1,6  k=1,M  CONTAIN THE M 6*6 MATRICES
C             M    NUMBER OF MATRICES STORED IN DAT.
C           NMX    DIMENSION OF THE LAST COMPONENT OF DAT.
C
C   OUTPUT:
C
C    VECS(i,j,k)  i=1,6  j=1,6  k=0,M  DIRECTIONS (MODULUS 1) AT THE
C                 BEGINING OF EACH INTERVAL DEFINED BY DAT(*,*,k).
C                 THIS IS VECS(*,j,k)=DAT(*,*,k)*VECS(*,j,k-1) EXCEPT
C                 FOR A MULTIPLICATIVE CONSTANT FOR EACH DIRECTION
C                 j=1,6 AND MATRIX k=1,M.
C                 THE DIRECTIONS j=1 and 2 REFER TO THE UNSTABLE AND
C                 STABLE RESPECTIVELY, AND j=3,6 TO THE CENTRAL ONES.
C
C    SUBROUTINES USED: INVERS,POWERS,GET,MATVEC,ORTO.
C    FUNCTIONS USED: VNORM.
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MMAX=1001)
      DIMENSION DAT(6,6,NMX),DATINV(6,6,MMAX),AUX(6,6),BUX1(6),IPIUX(6),
     &  A(6,6),VEP1(6),VEP2(6),CUX1(6),DELTA(2,MMAX),VECS(6,6,0:NMX),
     &  BUX2(6),CUX2(6),ACM(4,4)
      IF (M.GT.MMAX) THEN
      WRITE (6,*) 'VEPRO. PARAMETER NMAX MUST BE AT LEAST: ',M
      STOP
      ENDIF
      N=6
C
C  COMPUTING THE INVERSE MATRICES OF Ai STORED IN DAT. THE
C  INVERSE MATRICES OF Ai ARE STORED IN DATINV IN AN
C  ANALOGOUS WAY.
C
      CALL INVERS(DAT,DATINV,M,N,AUX,A,BUX1,IPIUX)
      DO 2 I=1,N
      VEP1(I)=0.D0
      VEP2(I)=0.D0
2     CONTINUE
      VEP1(1)=1.D0
      VEP2(2)=1.D0
C
C  COMPUTING THE UNSTABLE AND STABLE EIGENVECTORS USING THE
C  POWER METHOD. THEY ARE STORED IN VEP1 AND VEP2 RESPECT.
C  AFTER THIS COMPUTATION, THE SAME DIRECTIONS ARE COMPUTED
C  IN THE INTERMEDIATE STEPS AND STORED IN VECS.
C
      CALL POWERS(DAT,M,N,VEP1,VAP1,AUX,BUX1,CUX1,1)
      CALL POWERS(DATINV,M,N,VEP2,VAP2,AUX,BUX1,CUX1,-1)
C---- Per si es vol aprofitar la simetria WU WS ----------
C      VEP1(1)=VEP2(1)
C      VEP1(2)=-VEP2(2)
C      VEP1(3)=VEP2(3)
C      VEP1(4)=-VEP2(4)
C      VEP1(5)=VEP2(5)
C      VEP1(6)=-VEP2(6)
C---------------------------------------------------------
*      WRITE (*,*) 'VEPRO. LOG10 DELS VAPS HIPERBOLICS: ',VAP1,-VAP2
*C     WRITE(*,*)
*      DO 1 I=1,N
*      WRITE(*,*)'VEP1',VEP1(I),'   VEP2',VEP2(I)
*1     CONTINUE

      DO 3 I=1,N
      VECS(I,1,0)=VEP1(I)
      VECS(I,2,M)=VEP2(I)
3     CONTINUE
      DO 4 K=1,M
      DO 5 I=1,N
      BUX1(I)=VECS(I,1,K-1)
      BUX2(I)=VECS(I,2,M-K+1)
5     CONTINUE
      CALL GET (DAT,AUX,M,N,K)
      CALL MATVEC (N,AUX,BUX1,CUX1,0)
      CALL GET (DATINV,AUX,M,N,M-K+1)
      CALL MATVEC (N,AUX,BUX2,CUX2,0)
      VN1=VNORM (CUX1,N)
      VN2=VNORM (CUX2,N)
      DELTA(1,K)=VN1
      DELTA(2,M-K+1)=VN2
      DO 6 I=1,N
      VECS(I,1,K)=CUX1(I)/VN1
      VECS(I,2,M-K)=CUX2(I)/VN2
6     CONTINUE
4     CONTINUE
C
C  COMPUTING A BASIS FOR THE REMAINING CENTRAL PART.
C
      DO 10 IJ=3,N
*      WRITE(*,*)'VEPRO. CENTRAL COMPONENT, IJ=',IJ
      DO 11 I=1,N
      VECS(I,IJ,0)=0.D0
      BUX1(I)=0.D0
11    CONTINUE
      VECS(IJ,IJ,0)=1.D0
      BUX1(IJ)=1.D0
      DO 12 K=1,M
      CALL GET (DAT,AUX,M,N,K)
      CALL MATVEC (N,AUX,BUX1,CUX1,0)
      GAMMA=0.D0
      DO 13 I=1,N
      GAMMA=GAMMA+CUX1(I)*VECS(I,1,K)
13    CONTINUE
      DO 14 I=1,N
      BUX1(I)=CUX1(I)-GAMMA*VECS(I,1,K)
      VECS(I,IJ,K)=BUX1(I)
14    CONTINUE
      DEL=1.D0
      DO 15 IK=K-1,0,-1
      DEL=DEL*DELTA(1,IK+1)
      DO 16 I=1,N
      VECS(I,IJ,IK)=VECS(I,IJ,IK)-(GAMMA/DEL)*VECS(I,1,IK)
16    CONTINUE
15    CONTINUE
12    CONTINUE
10    CONTINUE


C La seguent rutina posa la part central DEL PRODUCTE de les
C matrius DAT com invariant. La matriu restringida al subespai
C invariant central ve donada per ACM. (Aquest pas segons per
C a que es pot saltar).
      CALL DIAHIP(VECS,DELTA,ACM,M)

      DO 20 IK=0,M
      DO 18 IJ=3,N
      VN1=0.D0
      DO 17 I=1,N
      VN1=VN1+VECS(I,IJ,IK)*VECS(I,IJ,IK)
17    CONTINUE
      VN1=DSQRT(VN1)
      DO 19 I=1,N
      VECS(I,IJ,IK)=VECS(I,IJ,IK)/VN1
19    CONTINUE
18    CONTINUE
C
C  IF WANTED, THE OBTAINED VECTORS CAN BE ORTONORMALIZED
C  BUT THEN THE CONTINUITY IS LOST. REMOVE THE C OF THE
C  NEXT LINE IF THIS IS WANTED.
C      CALL ORTO (VECS,IK,M)
20    CONTINUE

C  Print on screen.

*      WRITE (*,*) 'VEPRO. VECS[0] AFTER DIAHIP (1:2): '
*      DO 21 I=1,N
*      WRITE (*,*)  (VECS(I, J, 0), J = 1, 2)
*21    CONTINUE
*      WRITE (*,*) 'VEPRO. VECS[0] AFTER DIAHIP (3:4): '
*      DO 22 I=1,N
*      WRITE (*,*)  (VECS(I, J, 0), J = 3, 4)
*22    CONTINUE
*      WRITE (*,*) 'VEPRO. VECS[0] AFTER DIAHIP (5:6): '
*      DO 23 I=1,N
*      WRITE (*,*)  (VECS(I, J, 0), J = 5, 6)
*23    CONTINUE
*
*      WRITE (*,*) 'VEPRO. VECS[M] AFTER DIAHIP (1:2): '
*      DO 31 I=1,N
*      WRITE (*,*)  (VECS(I, J, M), J = 1, 2)
*31    CONTINUE
*      WRITE (*,*) 'VEPRO. VECS[M] AFTER DIAHIP (3:4): '
*      DO 32 I=1,N
*      WRITE (*,*)  (VECS(I, J, M), J = 3, 4)
*32    CONTINUE
*      WRITE (*,*) 'VEPRO. VECS[M] AFTER DIAHIP (5:6): '
*      DO 33 I=1,N
*      WRITE (*,*)  (VECS(I, J, M), J = 5, 6)
*33    CONTINUE


      RETURN
      END


      SUBROUTINE DIAHIP(VECS,DELTA,ACM,M)
C************************************************************************
C GIVEN M+1 SETS OF BASIS OF 6 VECTORS STORED BY COLUMNS IN VECS(6,6,0:M)
C ASSUMING THAT VECS(*,j,k+1)=A_k*VECS(*,j,k) WHERE A_k is a CERTAIN
C 6*6 MATRIX, AND VECS(*,j,M)=Di*VECS(*,j,0) FOR j=1,2, WHERE
C        D1=DELTA(1,1)*DELTA(1,2)*..*DELTA(1,M) and
C        D2=1.D0/(DELTA(2,1)*DELTA(2,2)*...*DELTA(2,M))
C (THIS IS VECS(*,1,0) and VECS(*,2,0) ARE EIGENVECTORS OF THE PRODUCT
C OF THE A_k MATRICES, AND DELTA MUST BE GIVEN AT THE INPUT).
C THIS ROUTINE MODIFIES THE VECTORS VECS(*,j,k) j=3,6 k=0,M, KEEPING
C ALWAYS VECS(*,j,k+1)=A_k*VECS(*,j,k), IN SUCH A WAY THAT THE NEW
C VECTORS VECS(*,j,M) j=3,6 SPAN THE SAME 4th DIMENSIONAL SPACE THAT
C THE VECTORS VECS(*,j,0) j=3,6.
C THEN RESTRICTED TO THE SUBSPACE j=3,6 VECS(*,j,M)=ACM*VECS(*,j,0)
C WHERE ACM IS A 4*4 MATRIX WHOSE jth COLUMN CONTAIN THE DESCOMPOSITION
C OF VECS(*,j,M) in terms of VECS(*,i,0) i=3,6. ACM IS GIVEN AS AN
C OUTPUT VALUE.
C************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VECS(6,6,0:M),DELTA(2,M),DESC(6,6),AUX(6,6),V(6),VB(4)
      DIMENSION X(6),IPIV(6),ALF(6),BET(6),A4M(4,4),B4M(4,4),VA(4)
      DIMENSION ACM(4,4)
      DO 3 J=1,6
      DO 5 I=1,6
      AUX(I,J)=VECS(I,J,0)
5     CONTINUE
3     CONTINUE
      CALL FACTOR(AUX,IPIV,V,6,IFLAG)
      IF (IFLAG.EQ.2) THEN
      WRITE (*,*) 'DIAHIP. INITIAL BASIS SEEMS DEPENDENT'
      STOP
      ENDIF
      DO 10 J=1,6
      DO 15 I=1,6
      V(I)=VECS(I,J,M)
15    CONTINUE
      CALL SUBST(AUX,V,X,IPIV,6)
      DO 17 I=1,6
      DESC(I,J)=X(I)
17    CONTINUE
10    CONTINUE

*      WRITE (*,*) 'DIAHIP. MATRIX IN INITIAL BASIS:'
*      WRITE (*,*) 'DIAHIP. DECS (1:2): '
*      DO 51 I=1,6
*      WRITE (*,*)  (DESC(I, J), J = 1, 2)
*51    CONTINUE
*      WRITE (*,*) 'DIAHIP. DECS (3:4): '
*      DO 52 I=1,6
*      WRITE (*,*)  (DESC(I, J), J = 3, 4)
*52    CONTINUE
*      WRITE (*,*) 'DIAHIP. DECS (5:6): '
*      DO 53 I=1,6
*      WRITE (*,*)  (DESC(I, J), J = 5, 6)
*53    CONTINUE

      AQD=1.D0
      BQD=1.D0
      DO 23 K=1,M
      AQD=AQD/DELTA(1,K)
      BQD=BQD/DELTA(2,K)
23    CONTINUE

*      WRITE(*,*)'DIAHIP. AQD = ',AQD,', BQD = ', BQD


      DO 25 J=1,4
      DO 26 I=1,4
      A4M(I,J)=DESC(J+2,I+2)*AQD
      B4M(I,J)=DESC(J+2,I+2)
C Puc carregar ACM abans del canvi ja que no varia pel canvi.
      ACM(I,J)=DESC(I+2,J+2)
26    CONTINUE
      A4M(J,J)=A4M(J,J)-1.D0
      B4M(J,J)=B4M(J,J)-BQD
      VA(J)=DESC(1,J+2)
      VB(J)=DESC(2,J+2)
25    CONTINUE

*      WRITE (*,*) 'DIAHIP. A4M (1:2): '
*      DO 31 I=1,4
*      WRITE (*,*)  (A4M(I, J), J = 1, 2)
*31    CONTINUE
*      WRITE (*,*) 'DIAHIP. A4M (3:4): '
*      DO 32 I=1,4
*      WRITE (*,*)  (A4M(I, J), J = 3, 4)
*32    CONTINUE
*
*      WRITE (*,*) 'DIAHIP. B4M (1:2): '
*      DO 41 I=1,4
*      WRITE (*,*)  (B4M(I, J), J = 1, 2)
*41    CONTINUE
*      WRITE (*,*) 'DIAHIP. B4M (3:4): '
*      DO 42 I=1,4
*      WRITE (*,*)  (B4M(I, J), J = 3, 4)
*42    CONTINUE


C NOTA: Si alguna de les matrius seguents (A4M o B4M) surt singular
C caldria entrar a la rutina uns altres vectors relacionats amb la
C part central ja que aquests semblen dependents. O be mirar les
C magnituds dels elements de la matriu i la tolerancia posada a FACTOR.
      CALL FACTOR(A4M,IPIV,V,4,IFLAG)
      IF (IFLAG.EQ.2) THEN
      WRITE (*,*) 'DIAHIP. MATRIX A4M SEEMS SINGULAR'
      STOP
      ENDIF
      CALL SUBST(A4M,VA,ALF(3),IPIV,4)
      CALL FACTOR(B4M,IPIV,V,4,IFLAG)
      IF (IFLAG.EQ.2) THEN
      WRITE (*,*) 'DIAHIP. MATRIX B4M SEEMS SINGULAR'
      STOP
      ENDIF
      CALL SUBST(B4M,VB,BET(3),IPIV,4)

*      WRITE (*,*) 'DIAHIP. ALF: '
*      DO 101 I=1,6
*      WRITE (*,*)  ALF(I)
*101    CONTINUE
*
*      WRITE (*,*) 'DIAHIP. BET: '
*      DO 102 I=1,6
*      WRITE (*,*)  BET(I)
*102    CONTINUE


      DO 35 J=3,6
      DO 50 K=M,0,-1
      DO 55 I=1,6
      VECS(I,J,K)=VECS(I,J,K)+ALF(J)*VECS(I,1,K)
55    CONTINUE
      IF (K.NE.0) ALF(J)=ALF(J)/DELTA(1,K)
50    CONTINUE
      DO 57 K=0,M
      DO 59 I=1,6
      VECS(I,J,K)=VECS(I,J,K)+BET(J)*VECS(I,2,K)
59    CONTINUE
      IF (K.NE.M) BET(J)=BET(J)/DELTA(2,K+1)
57    CONTINUE
35    CONTINUE
      RETURN
      END


      SUBROUTINE POWERS(DAT,M,N,VEP,VAP,AUX,BUX,CUX,IS)
C**********************************************************************
C  POWER METHOD IN A PRODUCT OF M MATRICES IN DIRECT OR INVERSE SENSE.
C  GIVEN MATRICES N*N STORED IN DAT(*,*,k), k=1,NIT, IT COMPUTES THE
C  LOG10 OF THE DOMINANT EIGENVALUE, VAP, AND THE ASSOCIATED
C  EIGENVECTOR, VEP, OF THE MATRIX A WHERE A IS:
C   A=DAT(*,*,M)*DAT(*,*,M-1)*...*DAT(*,*,1) IF IS=1, or,
C   A=DAT(*,*,1)*DAT(*,*,2)*...*DAT(*,*,M)   IF IS=-1.
C  USING THE POWER METHOD.
C  AUX, BUX AND CUX ARE AUXILIARY WORKING SPACE.
C
C  ROUTINES USED: GET,MATVEC.
C  FUNCTIONS USED: VNORM.
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C PRECISION AND MAX NUMBER OF ITERATIONS
      PARAMETER (PREC=1.D-14,ITM=10)
      DIMENSION DAT(N,N,M),VEP(N),AUX(N,N),BUX(N),CUX(N)
      KI=1
      KF=M
      INC=1
      IF (IS.EQ.-1) THEN
      KI=M
      KF=1
      INC=-1
      ENDIF
      IT=0
6     VAP=0.D0
      DO 1 K=KI,KF,INC
      CALL GET (DAT,AUX,M,N,K)
      CALL MATVEC (N,AUX,VEP,BUX,0)
      VN=VNORM(BUX,N)
      VAP=VAP+DLOG10(VN)
      DO 2 I=1,N
      VEP(I)=BUX(I)/VN
2     CONTINUE
1     CONTINUE
      IF (IT.EQ.0) THEN
      VAPV=VAP
      DO 3 I=1,N
      CUX(I)=VEP(I)
3     CONTINUE
      ELSE
      DVAP=DABS((VAP-VAPV)/VAP)
      DO 4 I=1,N
      CUX(I)=VEP(I)-CUX(I)
4     CONTINUE
      VN=VNORM(CUX,N)
      IF (IT.GT.ITM-5) THEN
      WRITE (*,*) 'POWERS. REL. ERR IN VAP: ',DVAP
      WRITE (*,*) 'POWERS. ABS. ERR IN VECTOR: ',VN
      ENDIF
      IF (DVAP.LT.PREC.AND.VN.LT.PREC) RETURN
      VAPV=VAP
      DO 5 I=1,N
      CUX(I)=VEP(I)
5     CONTINUE
      ENDIF
      IT=IT+1
      IF (IT.GT.ITM) THEN
      WRITE (*,*) 'POWERS. POWER METHOD SEEMS NOT CONVERGENT'
      WRITE (*,*) '        VALUE OF IS: ',IS
      STOP
      ENDIF
      GOTO 6
      END


         SUBROUTINE ORTO (A,IK,M)
C**********************************************************************
C  AUXILIARY ROUTINE. IT ORTONORMALIZES SOME VECTORS F A SET OF IK
C  BASIS STORED IN MATRIX A.
C**********************************************************************
         IMPLICIT REAL*8 (A-H,O-Z)
         DIMENSION A(6,6,0:M),B(6)
         DO 1 I=3,6
         DO 2 J=1,6
2        B(J)=A(J,I,IK)
         DO 3 K=2,I-1
         C=0.D0
         DO 4 L=1,6
4        C=C+B(L)*A(L,K,IK)
         DO 5 L=1,6
5        B(L)=B(L)-C*A(L,K,IK)
3        CONTINUE
         C=0.D0
         DO 6 L=1,6
6        C=C+B(L)*B(L)
         C=DSQRT(C)
         DO 7 L=1,6
7        B(L)=B(L)/C
         DO 8 L=1,6
8        A(L,I,IK)=B(L)
1        CONTINUE
         RETURN
         END


      SUBROUTINE MATINV (A,AINV,B,N,IPIVOT)
C**********************************************************************
C  AUXILIARY ROUTINE. GIVEN A N*N MATRIX A, IT COMPUTES ITS INVERSE
C  AND IT IS STORED IN THE N*N MATRIX AINV. B AND IPIVOT ARE TWO
C  AUXILYARY VECTORS OF LENGTH N NEEDED AS WORKING VECTORS.
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N*N),AINV(N*N),B(N),IPIVOT(N)
      CALL FACTOR(A,IPIVOT,B,N,IFLAG)
      IF(IFLAG.EQ.2)THEN
      WRITE (*,*) 'MATINV. SINGULAR MATRIX'
      STOP
      ENDIF
      DO 1 I=1,N
      B(I)=0
1     CONTINUE
      IBEG=1
      DO 2 J=1,N
      B(J)=1
      CALL SUBST (A,B,AINV(IBEG),IPIVOT,N)
      B(J)=0
      IBEG=IBEG+N
2     CONTINUE
      RETURN
      END


      SUBROUTINE FACTOR (A,IPIVOT,D,N,IFLAG)
C**********************************************************************
C  AUXILIARY ROUTINE. IT COMPUTES THE LU DESCOMPOSITION OF A N*N MATRIX
C  A. THE OUTPUT IS WRITTEN IN THE SAME MATRIX A, AND THE VECTOR IPIVOT
C  CONTAINS THE PERMUTATION USED IN THE DESCOMPOSITION.
C  IFLAG IS AN ERROR FLAG. IF ALL HAS BEEN O.K. IFLAG=1. IF THE
C  MATRIX A IS SINGULAR (LESS THAN PARAMETER TOLP) THEN IFLAG=2.
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (TOLP=1.D-8)
      DIMENSION A(N,N),IPIVOT(N),D(N)
      IFLAG=1
      DO 1 I=1,N
      IPIVOT(I)=I
      ROWMAX=0.D0
      DO 2 J=1,N
      ROWMAX=DMAX1(ROWMAX,DABS(A(I,J)))
2     CONTINUE
      IF (ROWMAX.LE.TOLP) THEN
      IFLAG=2
      RETURN
      ENDIF
      D(I)=ROWMAX
1     CONTINUE
      NM1=N-1
      IF(NM1.EQ.0)RETURN
      DO 3 K=1,NM1
      J=K
      KP1=K+1
      IP=IPIVOT(K)
      COLMAX=DABS(A(IP,K))/D(IP)
      DO 4 I=KP1,N
      IP=IPIVOT(I)
      AWIKOV=DABS(A(IP,K))/D(IP)
      IF (AWIKOV.LE.COLMAX) GOTO 4
      COLMAX=AWIKOV
      J=I
4     CONTINUE
      IF(COLMAX.LE.0)THEN
      IFLAG=2
      RETURN
      ENDIF
      IPK=IPIVOT(J)
      IPIVOT(J)=IPIVOT(K)
      IPIVOT(K)=IPK
      DO 5 I=KP1,N
      IP=IPIVOT(I)
      A(IP,K)=A(IP,K)/A(IPK,K)
      RATIO=-A(IP,K)
      DO 6 J=KP1,N
      A(IP,J)=RATIO*A(IPK,J)+A(IP,J)
6     CONTINUE
5     CONTINUE
3     CONTINUE
      RETURN
      END


      SUBROUTINE SUBST(A,B,X,IPIVOT,N)
C**********************************************************************
C  AUXILIARY ROUTINE. GIVEN A MATRIX A N*N FACTORIZED LU, AND THE
C  PERMUTATION VECTOR IPIVOT, IT SOLVES THE LINEAR SYSTEM AX=B.
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N,N),B(N),X(N),IPIVOT(N)
      IF(N.GT.1)GO TO 1
      X(1)=B(1)/A(1,1)
      RETURN
1     IP=IPIVOT(1)
      X(1)=B(IP)
      DO 2 K=2,N
      IP=IPIVOT(K)
      KM1=K-1
      SUM=0
      DO 3 J=1,KM1
      SUM=A(IP,J)*X(J)+SUM
3     CONTINUE
      X(K)=B(IP)-SUM
2     CONTINUE
      X(N)=X(N)/A(IP,N)
      K=N
      DO 4 NP1MK=2,N
      KP1=K
      K=K-1
      IP=IPIVOT(K)
      SUM=0
      DO 5 J=KP1,N
      SUM=A(IP,J)*X(J)+SUM
5     CONTINUE
      X(K)=(X(K)-SUM)/A(IP,K)
4     CONTINUE
      RETURN
      END


      REAL*8 FUNCTION VNORM(A,N)
C**********************************************************************
C  AUXILIARY FUNCTION. IT COMPUTES THE L_2 NORM OF A VECTOR A(N)
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N)
      VNORM=0
      DO 1 I=1,N
      VNORM=VNORM+A(I)*A(I)
1     CONTINUE
      VNORM=DSQRT(VNORM)
      RETURN
      END


      SUBROUTINE GET (DAT,A,M,N,K)
C**********************************************************************
C AUXILIARY ROUTINE. IT TAKES THE MATRIX A, N*N, FROM THE K-th PLACE
C OF THE MATRIX DAT.
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DAT(N,N,M),A(N,N)
      DO 1 I=1,N
      DO 2 J=1,N
      A(I,J)=DAT(I,J,K)
2     CONTINUE
1     CONTINUE
      RETURN
      END


      SUBROUTINE PUT (DAT,A,M,N,K)
C**********************************************************************
C AUXILIARY ROUTINE. IT STORES GIVEN A MATRIX A, N*N, IN THE K-th
C PLACE OF THE MATRIX DAT.
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DAT(N,N,M),A(N,N)
      DO 1 I=1,N
      DO 2 J=1,N
      DAT(I,J,K)=A(I,J)
2     CONTINUE
1     CONTINUE
      RETURN
      END


      SUBROUTINE INVERS (DAT,DATINV,M,N,AUX,AUXIN,BUX,IPIUX)
C**********************************************************************
C  AUXILIARY ROUTINE. IT COMPUTES THE INVERS MATRIX OF THE MATRICES
C  STORED IN DAT. RESULTS ARE STORED IN DATINV.
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DAT(N,N,M),DATINV(N,N,M),AUX(N,N),BUX(N),IPIUX(N),
     $ AUXIN(N,N)
      DO 1 K=1,M
      CALL GET (DAT,AUX,M,N,K)
      CALL MATINV (AUX,AUXIN,BUX,N,IPIUX)
      CALL PUT (DATINV,AUXIN,M,N,K)
1     CONTINUE
      RETURN
      END


      subroutine saxpy (n, alpha, x, y)
      integer n
      real*8 alpha, x(*), y(*)
c
c Saxpy: Compute y := alpha*x + y,
c where x and y are vectors of length n (at least).
c
c Local variables
      integer i
c
      do 10 i = 1, n
         y(i) = alpha*x(i) + y(i)
   10 continue
c
      return
      end


      subroutine MATVEC (n, A, x, y, lda)
      integer n, lda
      real*8 x(*), y(*), A(n,*)
c
c Compute y = A*x, where A is n by n and stored in an array
c with leading dimension n.
c
c Local variables
      integer i, j

c Initialize y
      do 10 i = 1, n
         y(i) = 0.0
   10 continue

c Matrix-vector product by saxpy on columns in A.
c Notice that the length of each column of A is m, not n!
      do 20 j = 1, n
         call saxpy (n, x(j), A(1,j), y)
   20 continue

      return
      end


