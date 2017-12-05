        SUBROUTINE ENDENR(XI,PERI,X1,X2,XMO1,XMO2,CAMP)
C**********************************************************************
C Donada una condicio inicial XI i un temps PERI per a un cert camp
C CAMP, aquesta rutina integra des de T=0 a T=PERI/2 endevant i
C guarda la pos+vel final a X1 i la matriu variacional a XMO1.
C Despres integra la mateixa c.i. enrera des de T=0 a T=-PERI/2 i
C guarda la pos+vel final a X2 i la matriu variacional a XMO2.
C (es pot usar per buscar orbites periodiques molt inestables)
C Esta preparada per fer arxius de dibuix pel canal 8.
C**********************************************************************
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION XI(42),X(42),X1(6),X2(6),XMO1(6,6),XMO2(6,6)
        DIMENSION B(42),F(42),R(13,42),XMO(6,6)
        EQUIVALENCE (X(7),XMO(1,1))
        SPERI=0.5D0*PERI
        DO 11 I=1,6
        X(I)=XI(I)
        DO 12 J=1,6
        XMO(J,I)=0.D0
12      CONTINUE
        XMO(I,I)=1.D0
11      CONTINUE
        TOL=1.D-13
        HMIN=1.D-4
        HMAX=1.D0
        H=1.D-2
        N=42
        T=0.D0
        open(8,file='op.dib')
        write (8,100) t,(x(i),I=1,6)
100     format(1x,7e24.16)
13      CALL RK78(T,X,N,H,HMIN,HMAX,TOL,R,B,F,CAMP)
        write (8,100) t,(x(i),I=1,6)
        IF (T.LT.SPERI) GOTO 13
        H=SPERI-T
        CALL RK78(T,X,N,H,HMIN,HMAX,TOL,R,B,F,CAMP)
        WRITE (NCS,*) 'ENDENR. CHECK INTEGRATION TIME: ',T-SPERI
        close(8)
        DO 20 I=1,6
        X1(I)=X(I)
        DO 21 J=1,6
        XMO1(J,I)=XMO(J,I)
21      CONTINUE
20      CONTINUE
C INTEGRACIO ENRERA
        DO 31 I=1,6
        X(I)=XI(I)
        DO 32 J=1,6
        XMO(J,I)=0.D0
32      CONTINUE
        XMO(I,I)=1.D0
31      CONTINUE
        H=-1.D-2
        T=0.D0
        open(8,file='op2.dib')
        write (8,100) t,(x(i),I=1,6)
33      CALL RK78(T,X,N,H,HMIN,HMAX,TOL,R,B,F,CAMP)
        write (8,100) t,(x(i),I=1,6)
        IF (T.GT.-SPERI) GOTO 33
        H=-SPERI-T
        CALL RK78(T,X,N,H,HMIN,HMAX,TOL,R,B,F,CAMP)
        WRITE (NCS,*) 'ENDENR. CHECK INTEGRATION TIME: ',T+SPERI
        close(8)
        DO 40 I=1,6
        X2(I)=X(I)
        DO 41 J=1,6
        XMO2(J,I)=XMO(J,I)
41      CONTINUE
40      CONTINUE
        RETURN
        END
