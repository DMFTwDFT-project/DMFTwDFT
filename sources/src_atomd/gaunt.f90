
SUBROUTINE SlaterF(Fn, U, J)
  IMPLICIT NONE
  REAL*8, intent(in)  :: U, J
  REAL*8, intent(out) :: Fn(0:3,0:3)

  Fn=0
  ! F0 for s-electrons
  Fn(0,0) = U
  ! F2 for p-electrons
  Fn(0,1) = U
  Fn(1,1) = 5*J
  ! F2 and F4 for d-electrons
  Fn(0,2) = U
  Fn(1,2) = 14./(1+0.625)*J
  Fn(2,2) = 0.625*Fn(1,2)
  ! F2, F4 and F6 for f-electrons
  Fn(0,3) = U
  Fn(1,3) = 6435./(286+195*0.668+250*0.494)*J
  Fn(2,3) = 0.668*Fn(1,3)
  Fn(3,3) = 0.494*Fn(1,3)
END SUBROUTINE SlaterF

SUBROUTINE cmp_all_Gaunt(gck)
  IMPLICIT NONE
  REAL*8, intent(out) :: gck(0:3,-3:3,-3:3,0:3)
  ! External function
  REAL*8 :: Gaunt
  ! Temporaries
  REAL*8, parameter :: pi = 3.14159266
  INTEGER :: m1, m2, k, l
  REAL*8  :: c
  DO l=0,3
     DO m1=-l,l
        DO m2=-l,l
           DO k=0,2*l,2
              c = Gaunt(l,m1,k,m1-m2,l,m2)*sqrt(4*pi/(2*k+1.))
              if (abs(c)<1e-10) c=0
              gck(l,m1,m2,k/2) = c
           ENDDO
        ENDDO
     ENDDO
  ENDDO
END SUBROUTINE cmp_all_Gaunt

REAL*8 function iFactorial(j)
  IMPLICIT NONE
  INTEGER, intent(in) :: j
  INTEGER :: i
  REAL*8 :: x
  if (j<0) print *, "iFactorial defined only for non-negative numbers!"
  x=1
  iFactorial = x
  if (j.eq.1) return
  DO i=2,j
     x = x*i
  END DO
  iFactorial = x
  return
end function iFactorial

REAL*8 function dFactorial(x)
  IMPLICIT NONE
  REAL*8, intent(in) :: x
  REAL*8, PARAMETER :: spi2 = 0.8862269254527579
  REAL*8 :: y, r
  r=1
  y=x
  DO WHILE(y.gt.1.0)
     r = r * y
     y = y -1.
  ENDDO
  IF (abs(y-0.5).LT.1e-10) r = r*spi2
  dFactorial = r
  return
END function dFactorial

REAL*8 function mone(i)
  INTEGER, intent(in) :: i
  mone = 1 - 2*MOD(abs(i),2)
  return
end function mone

REAL*8 function Delta(j1, j2, j)
  IMPLICIT NONE
  REAL*8, intent(in) :: j1, j2, j
  ! function calls
  REAL*8 :: dFactorial
  Delta = sqrt(dFactorial(j1+j2-j)*dFactorial(j1-j2+j)*dFactorial(-j1+j2+j)/dFactorial(j1+j2+j+1))
  return
END function Delta

REAL*8 function f3j(j1, m1, j2, m2, j3, m3)
  IMPLICIT NONE
  REAL*8, intent(in) :: j1, j2, j3, m1, m2, m3
  INTEGER            :: tmin, tmax, t
  REAL*8             :: sum, v1, v2, dn
  ! function calls
  REAL*8             :: dFactorial
  REAL*8             :: iFactorial
  REAL*8             :: Delta
  REAL*8             :: mone
  f3j=0
  IF (abs(m1+m2+m3) .GT. 1e-10) return
  IF (abs(j1-j2)-1e-14 .GT. j3 .OR. j3 .GT. j1+j2+1e-14) return
  if (abs(m1) .GT. j1 .OR. abs(m2) .GT. j2 .OR. abs(m3) .GT. j3) return
  tmin = INT(max(max(0.0,j2-j3-m1),j1-j3+m2)+1e-14)
  tmax = INT(min(min(j1+j2-j3,j1-m1),j2+m2)+1e-14)
  sum=0
  DO t=tmin, tmax
     v1 = dFactorial(j3-j2+m1+t)*dFactorial(j3-j1-m2+t)
     v2 = dFactorial(j1+j2-j3-t)*dFactorial(j1-m1-t)*dFactorial(j2+m2-t)
     sum = sum + mone(t)/(iFactorial(t)*v1*v2)
  END DO
  dn = dFactorial(j1+m1)*dFactorial(j1-m1)*dFactorial(j2+m2)*dFactorial(j2-m2)*dFactorial(j3+m3)*dFactorial(j3-m3)
  f3j = mone(INT(j1-j2-m3))*Delta(j1,j2,j3)*sqrt(dn)*sum
  return
END function f3j

REAL*8 function Gaunt(l1, m1, l2, m2, l3, m3)
  IMPLICIT NONE
  INTEGER, intent(in) :: l1, m1, l2, m2, l3, m3
  REAL*8, PARAMETER   :: pi = 3.14159265358979d0
  REAL*8 :: l1_, l2_, l3_, mm1_, m2_, m3_, zero
  ! function calls
  REAL*8             :: f3j
  REAL*8             :: mone
  l1_ = l1;   l2_ = l2;   l3_ = l3
  mm1_ = -m1; m2_ = m2; m3_ = m3
  zero = 0
  ! Calculates <Y_{l1m1}|Y_{l2m2}|Y_{l3m3}>
  if (l1.LT.0 .OR. l2.LT.0 .OR. l3.LT.0) print *, "Quantum number l must be non-negative!"
  Gaunt = mone(m1)*sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/(4*pi))*f3j(l1_,zero,l2_,zero,l3_,zero)*f3j(l1_,mm1_,l2_,m2_,l3_,m3_)
  return
END function Gaunt

REAL*8 function ClebschG(j,m,j1,m1,j2,m2)
  IMPLICIT NONE
  REAL*8, intent(in) :: j,m,j1,m1,j2,m2
  INTEGER            :: tmin, tmax, t
  REAL*8             :: sum, v1, v2
  ! function calls
  REAL*8             :: iFactorial
  REAL*8             :: dFactorial
  REAL*8             :: mone
  REAL*8             :: Delta

  ClebschG = 0
  IF (m1+m2 .NE. m) return
  tmin = INT(max(max(0.0,j2-j-m1),j1-j+m2)+1e-14)
  tmax = INT(min(min(j1+j2-j,j1-m1),j2+m2)+1e-14)
  sum=0;
  DO t=tmin, tmax
     v1 = sqrt((2*j+1)*dFactorial(j1+m1)*dFactorial(j1-m1)*dFactorial(j2+m2)*dFactorial(j2-m2)*dFactorial(j+m)*dFactorial(j-m))
     v2 = iFactorial(t)*dFactorial(j1+j2-j-t)*dFactorial(j1-m1-t)*dFactorial(j2+m2-t)*dFactorial(j-j2+m1+t)*dFactorial(j-j1-m2+t)
     sum = sum + mone(t)*v1/v2
  END DO
  ClebschG = sum*Delta(j1,j2,j)
  return
END function ClebschG
