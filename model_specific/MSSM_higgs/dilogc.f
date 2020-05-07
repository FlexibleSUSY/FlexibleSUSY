CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     SPENCE-FUNKTION KOMPLEX, FREI NACH HOLLIK                     C
C---------------------------------------------------------------------C
C     20.07.83    LAST CHANGED 10.05.89        ANSGAR DENNER        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION pCSPEN(ZZ)
      
      COMPLEX*16 pCSPEN,W,SUM,ZZ,Z,U
      DOUBLE PRECISION RZ,AZ,A1
      DOUBLE PRECISION B(9)
C     BEACHTE:                 B(N)=B2N
C     B(1)=1./6.
C     B(2)=-1./30.
C     B(3)=1./42.
C     B(4)=-1./30.
C     B(5)=5./66.
C     B(6)=-691./2730.
C     B(7)=7./6.
C     B(8)=-3617./510.
C     B(9)=43867./798.
C     PI=3.1415926535897932384
C     PI*PI/6.=1.6449..., PI*PI/3=3.28986...
      
      B(1)=0.1666666666666666666666666667d0
      B(2)=-0.0333333333333333333333333333d0
      B(3)=0.0238095238095238095238095238d0
      B(4)=-0.0333333333333333333333333333d0
      B(5)=0.0757575757575757575757575758d0
      B(6)=-0.2531135531135531135531135531d0
      B(7)=1.1666666666666666666666666667d0
      B(8)=-7.09215686274509804d0
      B(9)=54.97117794486215539d0
      
c     write(*,*) 'z:',z
      Z =ZZ*DCMPLX(1D0)
      RZ=DREAL(Z)
      AZ=CDABS(Z)
      A1=CDABS(1D0-Z)
c     write(*,*)'z, rz, az, a1:',z,rz,az,a1
C     IF((SNGL(RZ) .EQ. 0.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN
C     ---> CHANGED  10.5.89
      IF(AZ .LT. 1D-20) THEN
         pCSPEN=-CDLOG(1D0-Z)
c     write(*,*) 'cspen:', cspen
         RETURN
      END IF
      IF((SNGL(RZ) .EQ. 1.0) .AND. (SNGL(DIMAG(Z)) .EQ. 0.0)) THEN
         pCSPEN=1.64493406684822643D0
c     write(*,*) 'cspen:', cspen
         RETURN
      END IF
      IF(RZ.GT.5D-1) GOTO 20
      IF(AZ.GT.1D0) GOTO 10
      W=-CDLOG(1D0-Z)
      SUM=W-0.25D0*W*W
      U=W
      IF(CDABS(U).LT.1D-10) GOTO 2
c     write(*,*) 'u:',u
c     write(*,*) 'sum:',sum
      DO 1 K=1,9
         U=U*W*W/DFLOAT(2*K*(2*K+1))
         IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 2
         SUM=SUM+U*B(K)
 1    CONTINUE
 2    pCSPEN=SUM
c     write(*,*) 'cspen:', cspen
      RETURN
 10   W=-CDLOG(1D0-1D0/Z)
      SUM=W-0.25D0*W*W
      U=W
      IF(CDABS(U).LT.1D-10) GOTO 12
      
      DO 11 K=1,9
         U=U*W*W/DFLOAT(2*K*(2*K+1))
         IF(CDABS(B(K)*U/SUM).LT.1D-20) GOTO 12
         SUM=SUM+U*B(K)
 11   CONTINUE
 12   pCSPEN=-SUM-1.64493406684822643D0-.5D0*CDLOG(-Z)**2
c     write(*,*) 'cspen:', cspen
      RETURN
 20   IF(A1.GT.1D0) GOTO 30
      W=-CDLOG(Z)
      SUM=W-0.25D0*W*W
      U=W
      IF(CDABS(U).LT.1D-10) GOTO 22
      DO 21 K=1,9
         U=U*W*W/DFLOAT(2*K*(2*K+1))
         IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 22
         SUM=SUM+U*B(K)
 21   CONTINUE
 22   pCSPEN=-SUM+1.64493406684822643D0-CDLOG(Z)*CDLOG(1D0-Z)
c     write(*,*) 'cspen:', cspen
      RETURN
 30   W=CDLOG(1D0-1D0/Z)
      SUM=W-0.25D0*W*W
      U=W
      IF(CDABS(U).LT.1D-10) GOTO 32
      DO 31 K=1,9
         U=U*W*W/DFLOAT(2*K*(2*K+1))
         IF(CDABS(U*B(K)/SUM).LT.1D-20) GOTO 32
         SUM=SUM+U*B(K)
 31   CONTINUE
 32   pCSPEN=SUM+3.28986813369645287D0
     *     +.5D0*CDLOG(Z-1D0)**2-CDLOG(Z)*CDLOG(1D0-Z)
c     write(*,*) 'cspen:', cspen
      END



!>    @brief Complex dilogarithm.
!!    Introduced by Ben Allanach to provide a wrapper to C class
!!    @param a real part of argument
!!    @param b imag part of argument
!!    @param c real part of output
!!    @param d imag part of output
      subroutine dilogc(a, b, c, d)
     $ bind(C, name="dilogc_")
      
      implicit none
      
      double precision, intent(in) :: a, b
      double precision, intent(out) :: c, d
      complex*16 x, l, pCSPEN
      
      x = dcmplx(a, b)
      l = pCSPEN(x)
      
      c = dble(l)
      d = dimag(l)
      
      return
      end
