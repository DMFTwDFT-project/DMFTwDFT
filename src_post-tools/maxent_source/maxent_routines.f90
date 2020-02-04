! @Copyright 2007 Kristjan Haule
! 

!-------------------------------------
!   Classic Maximum Entropy Method
!-------------------------------------
SUBROUTINE MAXENT(a,rfac,alpha,temp,ker,sxt,axt,mm,dw,nt,nw,steps,iseed)
  !use IFPORT    ! This should be used for intel compiler!
  IMPLICIT NONE
  REAL*8, intent(inout) :: a(-nw:nw)
  REAL*8, intent(inout) :: rfac, alpha, temp
  REAL*8, intent(in)    :: ker(-nw:nw,nt), sxt(nt), axt(nt), mm(-nw:nw), dw(-nw:nw)
  INTEGER, intent(in)   :: nt, nw, steps, iseed
  ! locals
  REAL*8 :: da(-nw:nw),xt1(nt),xt2(nt)
  INTEGER:: i,j1,j2,jj,acc,try
  REAL*8 :: am,eps,x
  REAL*8 :: dj1,dj2,x1,x2,wgt
  REAL*8 :: aa1,aa2,de,p,mm1,mm2,arat
  
  am = maxval(a)
  da(:) = (a(:)/10.+0.01)*rfac
  
  acc=0
  try=0
  eps=1.E-12

  DO i=1,steps
     ! current G(tau) = A(x)*K(x,tau)
     xt1 = matmul(a , ker)
     ! computing current chi2
     x1 = sum(sxt*(axt-xt1)**2)

     ! number of Metroolis steps equal to the number of frequency points
     DO jj=1,2*nw+1
        ! finds two points (j1,j2) and the move for the spectra at these two points
        ! a(j1) -> a(j1)+dj1
        ! a(j2) -> a(j2)+dj2
        do     
           CALL random_number(x)
           j1=int(x*(2*nw+1))-nw
           CALL random_number(x)
           j2=int(x*(2*nw))-nw
           if (j2.GE.j1) j2=j2+1   ! j2!=j1
           
           do  ! need positive a(j1)+dj1
              CALL random_number(x)
              dj1=da(j1)*(x-0.5)              ! dj1 -> -0.5*da(j1) ... 0.5*da(j1)
              if (a(j1)+dj1.GE.0.) exit
           enddo
           CALL random_number(x)
           dj2=-dj1+0.05*da(j2)*(x-0.5)    ! dj2 -> -dj1 + -0.025*da(j2) ... 0.025*da(j2)
           IF (a(j2)+dj2.GE.0.) exit
        enddo
        
        IF (a(j1).GT.0.1*am) try=try+1  ! trial steps
        ! computing chi2 for trial step
        xt2(:)=xt1(:)+dj1*ker(j1,:)+dj2*ker(j2,:)
        x2 = sum(sxt*(axt-xt2)**2) ! chi2
     
        aa1=a(j1)+dj1
        aa2=a(j2)+dj2
        mm1=mm(j1)
        mm2=mm(j2)
        de=0.
        IF (aa1.GT.eps)   de=de-dw(j1)*aa1*LOG(aa1/mm1)
        IF (a(j1).GT.eps) de=de+dw(j1)*a(j1)*LOG(a(j1)/mm1)
        IF (aa2.GT.eps)   de=de-dw(j2)*aa2*LOG(aa2/mm2)
        IF (a(j2).GT.eps) de=de+dw(j2)*a(j2)*LOG(a(j2)/mm2)            
        wgt=((x1-x2)/2.+alpha*de)/temp ! difference in chi2+alpha*S between trial and current configuration

        ! Metropolis probability
        IF (wgt.GT.0.d0) THEN
           p=1.d0
        ELSEIF (wgt.GT.-100.d0) THEN
           p=EXP(wgt)
        ELSE
           p=0.
        ENDIF
        
        ! accepting with Metropolis probability
        CALL random_number(x)
        IF (x.LT.p) THEN 
           IF (a(j1).GT.0.1*am) acc=acc+1
           a(j1)=aa1
           a(j2)=aa2
           xt1(:) = xt2(:)
           x1=x2
        ENDIF
     ENDDO
     
     IF (MOD(i,100).EQ.0) THEN
        arat=real(acc)/real(try)
        IF (arat.GT.0.1) THEN ! if accepting probability is large, reduce the size of future steps
           IF (rfac.LT.0.01) rfac=rfac*1.5
        ELSE                  ! if accepting probability is small, increase the size of future steps
           IF (rfac.GT.0.001) rfac=rfac/1.5
        ENDIF
        
        am = maxval(a)
        da(:) = (0.1*a(:)+0.01*am)*rfac
        
        IF (MOD(i,steps/10).EQ.0) THEN
           WRITE(*,'(A,I7,2x,A,f12.5,2x,A,I7,2x,A,I9,2x,A,f7.4)') 'ann-step=', i,'chi2=', x1, 'accepted=', acc, &
                'tried=', try, 'new stepsize=', rfac
        ENDIF

        temp=temp/1.5
        acc=0
        try=0
     ENDIF
  ENDDO
  
  RETURN
END SUBROUTINE MAXENT

!---------------------------------------------------------
! Computes Entropy term as required for maxent algorithm
!---------------------------------------------------------
REAL*8 FUNCTION ENTROPY(a,mm,dw,nw)
  IMPLICIT NONE
  REAL*8, intent(in)  :: a(-nw:nw), mm(-nw:nw), dw(-nw:nw)
  INTEGER, intent(in) :: nw
  ! locals
  INTEGER :: i
  REAL*8 :: ent, eps
  ent=0.
  eps=1.E-12
  DO i=-nw,nw  ! 10
     IF (a(i).GT.eps.AND.mm(i).GT.eps) THEN
        ent=ent-a(i)*LOG(a(i)/mm(i))*dw(i) + (a(i)-mm(i))*dw(i)
     ENDIF
  ENDDO
  ENTROPY=ent
  RETURN
END FUNCTION ENTROPY

!------------------------------------------------------------------------------
! computes : Tr( lam * (lam+alpha)^-1 )
! where      lam = sqrt(a(w1))*ker(w1,:)*sxt(:)*ker(w2,:)*sqrt(a(w2))/dw
!------------------------------------------------------------------------------
REAL*8 FUNCTION LAMBDAC(alpha,a,w,dw,dlda,nw)
  IMPLICIT NONE
  !
  ! computes : Tr( lam * (lam+alpha)^-1 )
  ! where      lam = sqrt(a(w1))*ker(w1,:)*sxt(:)*ker(w2,:)*sqrt(a(w2))/dw
  !            
  REAL*8, intent(in):: alpha
  REAL*8, intent(in) :: a(-nw:nw), w(-nw:nw), dw(-nw:nw), dlda(2*nw+1,2*nw+1)
  INTEGER, intent(in):: nw
  ! locals
  REAL*8  :: tr
  INTEGER :: j, i, iw, jw
  REAL*8, allocatable  :: lam(:,:),alam(:,:),alami(:,:)
  
  allocate( lam(2*nw+1,2*nw+1), alam(2*nw+1,2*nw+1), alami(2*nw+1,2*nw+1) )
  
  lambdac=0.0

  DO j=1,2*nw+1
     DO i=1,2*nw+1
        iw=i-nw-1
        jw=j-nw-1
        lam(i,j) = SQRT(a(iw))*dlda(i,j)*SQRT(a(jw))*SQRT(dw(iw)*dw(jw))
        alam(i,j)=lam(i,j)
     ENDDO
     alam(j,j)=lam(j,j)+alpha
  ENDDO
  
  CALL INVERSE_MATRIX(alami, alam, 2*nw+1, 2*nw+1)
  
  deallocate( alam )
  
  tr=0.
  DO i=1,2*nw+1
     DO j=1,2*nw+1
        tr=tr+lam(i,j)*alami(j,i)
     ENDDO
  ENDDO
  LAMBDAC = tr

  deallocate( lam, alami )
  RETURN
END FUNCTION LAMBDAC

!------------------------------------------------------------------------------
! computes dlda(w1,w2)=ker(w1,:)*sxt(:)*ker(w2,:)/dw**2
!    which is needed for computation of lambda
!------------------------------------------------------------------------------
SUBROUTINE INITDLDA(dlda, w, dw, ker, sxt, nt, nw)
  IMPLICIT NONE
  REAL*8, intent(out) :: dlda(2*nw+1,2*nw+1)
  REAL*8, intent(in)  :: w(-nw:nw), dw(-nw:nw), ker(-nw:nw,nt), sxt(nt)
  INTEGER, intent(in) :: nw, nt
  ! locals
  INTEGER:: i,j,iw,jw,it
  REAL*8 :: dsum
  !
  ! computing the following quantity:
  !    dlda(w1,w2)=ker(w1,:)*sxt(:)*ker(w2,:)/dw**2
  !
  DO j=1,2*nw+1
     DO i=1,2*nw+1
        iw=i-nw-1
        jw=j-nw-1
        dsum=0.0
        DO it=1,nt
           dsum = dsum + ker(iw,it)*sxt(it)*ker(jw,it)
        ENDDO
        dlda(i,j)=dsum/(dw(iw)*dw(jw))
     ENDDO
  ENDDO
  call random_seed()
  RETURN
END SUBROUTINE INITDLDA

!--------------------------------------
! Needed for smoothing the maxent data
!--------------------------------------
SUBROUTINE SMOOTH(a,ns,nw)
  IMPLICIT NONE
  REAL*8, intent(inout) :: a(-nw:nw)
  INTEGER, intent(in)   :: ns, nw
  ! locals
  REAL*8 :: a1(-nw:nw)
  INTEGER i,j,i1,i2,ii
  DO i=-nw,nw
     i1=i-ns
     IF (i1.LT.-nw) i1=-nw
     i2=i+ns
     IF (i2.GT.nw) i2=nw
     a1(i)=0.
     ii=0
     DO j=i1,i2
        ii=ii+1
        a1(i)=a1(i)+a(j)
     ENDDO
     a1(i)=a1(i)/real(ii)
  ENDDO
  a(:) = a1(:)
  RETURN
END SUBROUTINE SMOOTH


SUBROUTINE INVERSE_MATRIX(Ainv, A, N, LDA)
  IMPLICIT NONE
  REAL*8, intent(out):: Ainv(LDA,LDA)
  REAL*8, intent(in) :: A(LDA,LDA)
  INTEGER, intent(in):: N, LDA
  !
  INTEGER, ALLOCATABLE :: IPIV(:)
  REAL*8, ALLOCATABLE  :: WORK(:)
  INTEGER :: LWORK, INFO
  !
  external DGETRF
  external DGETRI
  !
  LWORK = N*N
  ALLOCATE (WORK(LWORK), IPIV(N))
  
  ! DGETRF computes an LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges.
  Ainv(:,:) = A(:,:)

  CALL DGETRF( N, N, Ainv, LDA, IPIV, INFO )
  IF(info.NE.0) THEN
     print *, 'ERROR: LU decomposition failed. INFO=', info
  ENDIF
  
  !  DGETRI computes the inverse of a matrix using the LU factorization computed by DGETRF.
  CALL DGETRI(N, Ainv, LDA, IPIV, WORK, LWORK, INFO)
  
  IF (info.NE.0) THEN
     print *,  'ERROR: Matrix inversion failed. INFO=', info
  ENDIF
  DEALLOCATE (WORK,IPIV)
END SUBROUTINE INVERSE_MATRIX


!c------------------------------------------------------------------------------
SUBROUTINE FLATMODEL(mm,nw)
  IMPLICIT NONE
  REAL*8,  intent(out):: mm(-nw:nw)
  INTEGER, intent(in) :: nw
  !locals
  REAL*8  :: omega
  INTEGER :: i
  print*,'Reading model model.dat'
  OPEN(UNIT=99,FILE='model.dat',STATUS='old')
  DO i=-nw,nw
     READ(10,*) omega, mm(i)
  ENDDO 
  RETURN
END SUBROUTINE FLATMODEL
!c------------------------------------------------------------------------------
SUBROUTINE FLATMODEL0(mm,nw)
  IMPLICIT NONE
  REAL*8,  intent(out):: mm(-nw:nw)
  INTEGER, intent(in) :: nw
  !locals
  INTEGER :: i
  print*,'using FLATMODEL'
  DO i=-nw,nw 
     mm(i)=1.
  ENDDO
  RETURN
END SUBROUTINE FLATMODEL0
!c------------------------------------------------------------------------------
SUBROUTINE INITF0(f0,f1,f2,w,nw)
  IMPLICIT NONE
  REAL*8, intent(out):: f0(-nw:nw), f1(-nw:nw), f2(-nw:nw)
  REAL*8, intent(in) :: w(-nw:nw)
  INTEGER, intent(in):: nw
  ! locals
  INTEGER:: iw
  REAL*8 :: dw
  dw=w(2)-w(1)
  DO iw=-nw,nw   ! 10
     f0(iw)=dw
     f1(iw)=dw*w(iw)
     f2(iw)=dw*w(iw)**2
  ENDDO
  RETURN
END SUBROUTINE INITF0
!c------------------------------------------------------------------------------
SUBROUTINE INITKER_FERMION(ker,w,dw,beta,t,nt,nw)
  IMPLICIT NONE
  REAL*8, intent(out) :: ker(-nw:nw,nt)
  REAL*8, intent(in)  :: w(-nw:nw), dw(-nw:nw), t(nt), beta
  INTEGER, intent(in) :: nw, nt
  ! locals
  INTEGER :: it,iw
  REAL*8  :: tau, omega
  integer, volatile :: iCopy
  DO iw=-nw,nw
     DO it=1,nt
        iCopy = it  ! Intel ifort has a bug, which optimizes out this loop, and code fails
        tau = t(it)
        omega = w(iw)
        if (omega*beta<-50.)then
           ker(iw,it) = dw(iw) * EXP( (beta-tau)*omega )
        else 
           ker(iw,it) = dw(iw) * EXP(-tau*omega)/(1.+EXP(-beta*omega))
        endif
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE INITKER_FERMION

!c------------------------------------------------------------------------------
SUBROUTINE INITKER_BOSON_SYMM(ker,w,dw,beta,t,nt,nw)
  !!! Symmetric version for the case of G(tau)=G(beta-tau) hence
  !!! chi(x)=chi(-x)
  IMPLICIT NONE
  REAL*8, intent(out) :: ker(-nw:nw,nt)
  REAL*8, intent(in)  :: w(-nw:nw), dw(-nw:nw), t(nt), beta
  INTEGER, intent(in) :: nw, nt
  ! locals
  INTEGER :: it,iw
  REAL*8  :: tw, bw, tau, omega
  DO iw=-nw,nw
     DO it=1,nt
        tau = t(it)
        omega = w(iw)
        if (omega*beta<-50.)then
           ker(iw,it) = -0.5 * dw(iw) * omega * ( EXP((beta-tau)*omega) + EXP(tau*omega) )
        else if(abs(omega)<1e-6) then
           ker(iw,it) = dw(iw)/beta
        else
           tw = EXP(-tau*omega) + EXP(-(beta-tau)*omega)
           bw = 1.-EXP(-beta*omega)
           ker(iw,it) = 0.5 * dw(iw) * omega * tw/bw
        endif
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE INITKER_BOSON_SYMM

!c------------------------------------------------------------------------------
SUBROUTINE INITKER_BOSON(ker,w,dw,beta,t,nt,nw)
  IMPLICIT NONE
  REAL*8, intent(out) :: ker(-nw:nw,nt)
  REAL*8, intent(in)  :: w(-nw:nw), dw(-nw:nw), t(nt), beta
  INTEGER, intent(in) :: nw, nt
  ! locals
  INTEGER :: it,iw
  REAL*8  :: tw, bw, tau, omega
  DO iw=-nw,nw
     DO it=1,nt
        tau = t(it)
        omega = w(iw)
        if (omega*beta<-50.)then
           ker(iw,it) = -dw(iw) * omega * EXP((beta-tau)*omega) 
        else if(abs(omega)<1e-6) then
           ker(iw,it) = dw(iw)/beta
        else
           if(abs(omega*beta)>100.) then
              bw = 1.
           else
              bw = 1.-EXP(-beta*omega)
           endif
           tw = EXP(-tau*omega)
           ker(iw,it) = dw(iw) * omega * tw/bw
        endif
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE INITKER_BOSON

!-----------------------------
! Inverse Fourier inside part
!-----------------------------
SUBROUTINE FourPart(Gtau,t,Gm,om,ah,beta,nom)
  IMPLICIT NONE
  REAL*8, intent(out)    :: Gtau
  COMPLEX*16, intent(in) :: Gm(nom)
  REAL*8, intent(in)     :: om(nom), t, ah, beta
  INTEGER, intent(in)    :: nom
  ! locals
  INTEGER:: im
  REAL*8 :: dsum
  dsum=0.
  do im=1,nom
     dsum = dsum + cos(om(im)*t)*real(Gm(im)) + sin(om(im)*t)*(aimag(Gm(im))+ah/om(im))
  enddo
  Gtau = 2*dsum/beta-0.5*ah
END SUBROUTINE FourPart
SUBROUTINE FourPartB(Gtau,t,Gm,om,beta,nom)
  IMPLICIT NONE
  REAL*8, intent(out)    :: Gtau
  COMPLEX*16, intent(in) :: Gm(nom)
  REAL*8, intent(in)     :: om(nom), t, beta
  INTEGER, intent(in)    :: nom
  ! locals
  INTEGER:: im
  REAL*8 :: dsum
  dsum=0.
  do im=2,nom
     dsum = dsum + cos(om(im)*t)*real(Gm(im))
  enddo
  dsum = dsum + 0.5*real(Gm(1))
  Gtau = 2*dsum/beta
END SUBROUTINE FourPartB

!  Recursion for Pade coefficient (J.Serene)
!****************************
SUBROUTINE Padecof(Pt,gn,zn,nn)
!****************************
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: Pt(nn)
  COMPLEX*16, intent(in)  :: gn(nn),zn(nn)
  INTEGER, intent(in)     :: nn
  ! local
  COMPLEX*16  :: p(nn,nn)
  INTEGER :: i, j
  p=0.0
  do j=1,nn
     p(1,j) = gn(j)
  enddo
  do j=2,nn
     do i=2,j
        p(i,j)=(p(i-1,i-1)-p(i-1,j))/(zn(j)-zn(i-1))/p(i-1,j)
     enddo
  enddo
  do j=1,nn
     Pt(j)=p(j,j)
  enddo
  return
end SUBROUTINE Padecof

!  Calculation of a Green's function for a given pade-coeff-p(i,j)
!  on the real axis e=e+i0
!*****************************
SUBROUTINE PadeG(Gx,x,zn,Pt,nn)
!*****************************
  IMPLICIT NONE
  COMPLEX*16, intent(out) :: Gx
  COMPLEX*16, intent(in)  :: x, zn(nn),Pt(nn)
  INTEGER, intent(in)     :: nn
  ! locals
  COMPLEX*16 :: aw(0:nn),bw(0:nn)
  INTEGER :: i
  aw(0)=(0.d0,0.d0)
  aw(1)=Pt(1)
  bw(0)=(1.d0,0.d0)
  bw(1)=(1.d0,0.d0)
  do i=1,nn-1
     aw(i+1)=aw(i)+(x-zn(i))*Pt(i+1)*aw(i-1)
     bw(i+1)=bw(i)+(x-zn(i))*Pt(i+1)*bw(i-1)
  enddo
  Gx=aw(nn)/bw(nn)
  return
end SUBROUTINE PadeG

