! @Copyright 2007 Kristjan Haule
! 

program maxentropy
  IMPLICIT NONE
  INTEGER :: nt, nw, seed
  REAL*8  :: beta
  REAL*8, allocatable :: t(:), axt(:), sxt(:), xt1(:)
  REAL*8, allocatable :: w(:), a(:)
  INTEGER :: it,nr,iw,steps,iflat,idg
  REAL*8  :: deltag, alpha0,dw,PI
  REAL*8  ::  normalization

  PI=Dacos(-1.D0)
  
  READ(*,*) nt
  WRITE(*,*) 'Number of tau points: ', nt
  READ(*,*) idg,deltag
  WRITE(*,*) 'idg=', idg, 'deltag=', deltag
  READ(*,*) beta
  WRITE(*,*) 'beta = ', beta
  READ(*,*)nw
  WRITE(*,*)'Number of frequencies=', nw
  READ(*,*)dw
  WRITE(*,*)'Delta frequency=', dw
  READ(*,*)steps
  WRITE(*,*)'Number of annealing steps=', steps
  READ(*,*)alpha0
  WRITE(*,*)'Starting Alpha=', alpha0
  READ(*,*)nr
  WRITE(*,*)'Number of smoothed runs=', nr
  seed=53718551
  READ(*,*)seed
  WRITE(*,*)'Random number seed=', seed
  READ(*,*)iflat
  WRITE(*,*)'iflat=', iflat
  READ(*,*) normalization
  WRITE(*,*) 'Normalization=', normalization


  ALLOCATE( t(nt), axt(nt), sxt(nt), xt1(nt) )
  ALLOCATE( w(-nw:nw), a(-nw:nw) )
  
  OPEN(UNIT=10,FILE='Gtau1.dat',STATUS='old')
  DO it=1,nt
     READ(10,*) t(it),axt(it)
     !print *, t(it), axt(it)
     if(idg.eq.0) then
        sxt(it) = axt(it)*deltag 
     else
        sxt(it)=deltag 
     endif
     axt(it)=-axt(it)                      !  This makes G(tau) positive ! THINK
     IF (ABS(sxt(it)).LT.0.00001) sxt(it)=0.00001
     sxt(it)=1./sxt(it)**2
  ENDDO
  
  ! frequency mesh
  DO iw=-nw,nw
     w(iw)=iw*dw
  ENDDO

  OPEN(UNIT=10,FILE='dos',STATUS='unknown')
  OPEN(UNIT=44,FILE='inst_dos',STATUS='unknown')
  
  call srand(seed)
  CALL main_maxent(a,w,axt,sxt,t,iflat,steps,beta,alpha0,nr,normalization,nt,nw)
  
  DO iw=-nw,nw
     WRITE(10,'(2F15.8)')w(iw),a(iw)
  ENDDO
  
  CLOSE(10)
  CLOSE(44)
      
  DEALLOCATE( t, axt, sxt, xt1 )
  DEALLOCATE( w, a )
  
END program maxentropy

!c------------------------------------------------------------------------------
SUBROUTINE main_maxent(a,w,axt,sxt,t,iflat,steps,beta,alpha0,nr,n0,nt,nw)
  IMPLICIT NONE
  REAL*8, intent(inout):: a(-nw:nw), w(-nw:nw)
  REAL*8, intent(in)   :: axt(nt), sxt(nt), t(nt)
  INTEGER, intent(in)  :: iflat, steps, nr, nt, nw
  REAL*8, intent(in)   :: alpha0, n0, beta
  ! locals
  REAL*8  :: ker(-nw:nw,nt)
  REAL*8  :: dlda(2*nw+1,2*nw+1)
  REAL*8  :: f0(-nw:nw)
  REAL*8  :: f1(-nw:nw)
  REAL*8  :: f2(-nw:nw)
  REAL*8  :: mm(-nw:nw)
  REAL*8  :: rfac,temp,alpha,ent,tr,ls !,adiv,s0,s1,s2,
  INTEGER :: st, i, iw
  ! externals
  REAL*8  :: ENTROPY, LAMBDA

  CALL INITKER(ker,w,beta,t,nt,nw)
  CALL INITF0(f0,f1,f2,w,nw)
  CALL INITDLDA(dlda,w,ker,sxt,nt,nw)
  
  if(iflat.eq.0) then
     CALL FLATMODEL0(mm,nw)
  else
     CALL FLATMODEL(mm,nw)
  endif
  mm(:) = mm(:) * (n0/dot_product(mm,f0))  ! normalize model to n0
  
  DO i=-nw,nw
     a(i)=rand()
  ENDDO
  a(:) = a(:) * (n0/dot_product(a,f0))  ! normalize a to n0
  
  temp=10.
  rfac=1.
  alpha=alpha0

  ! Will iterate unti ls/tr ~ 1.
  ! ls = -2*S*alpha
  ! tr = Tr( lam * (lam+alpha)^-1 )
  !       where lam = sqrt(a(w1))*ker(w1,:)*sxt(:)*ker(w2,:)*sqrt(a(w2))/dw
  DO
     WRITE(*,'(A,f12.5,2x,A,f12.5)') 'Restarting maxent with rfac=', rfac, 'alpha=', alpha
     CALL MAXENT(a,rfac,alpha,temp,ker,sxt,axt,mm,f0,nt,nw,steps)
     ent=ENTROPY(a,mm,f0,nw)
     tr=LAMBDA(alpha,a,w,dlda,nw)
     ls=-2.*ent*alpha
     WRITE(*,'(A,f12.5,2x,A,f10.5,2x,A,f10.5)') 'Finished maxent with alpha=',alpha,'-2*alpha*S=',ls,'Trace=',tr
     WRITE(*,*)'ratio= ', ls/tr
     temp=0.001
     rfac=0.05
     
     ! write the dos after every iteration
     DO iw=-nw,nw
        WRITE(44,'(3F15.8)') w(iw),a(iw)
     ENDDO
     WRITE(44,*)
     
     IF (ABS(ls/tr-1.).LT.0.1) THEN
        EXIT
     ELSE IF (tr/ls.LT.0.05) THEN
        alpha=alpha*0.05
     ELSE
        alpha=alpha*(tr/ls)*(1.+0.001*(rand()-0.5))
     ENDIF
     
  ENDDO
      
  DO st=1,nr
     WRITE(*,*)'Starting smoothing ',st         
     CALL SMOOTH(a,3,nw)
     a(:) = a(:) * (n0/dot_product(a,f0)) ! Normalize a to n0
     temp=0.005
     rfac=0.005
     CALL MAXENT(a,rfac,alpha,temp,ker,sxt,axt,mm,f0,nt,nw,steps)
  ENDDO
  RETURN
END SUBROUTINE main_maxent

!c-------------------------------------
!c   Classic Maximum Entropy Method
!c-------------------------------------
SUBROUTINE MAXENT(a,rfac,alpha,temp,ker,sxt,axt,mm,f0,nt,nw,steps)
  IMPLICIT NONE
  REAL*8, intent(inout) :: a(-nw:nw)
  REAL*8, intent(inout) :: rfac, alpha, temp
  REAL*8, intent(in)    :: ker(-nw:nw,nt), sxt(nt), axt(nt), mm(-nw:nw), f0(-nw:nw)
  INTEGER, intent(in)   :: nt, nw, steps
  ! locals
  REAL*8 :: da(-nw:nw),xt1(nt),xt2(nt)
  INTEGER:: i,j1,j2,jj,acc,try
  REAL*8 :: am,eps
  REAL*8 :: dj1,dj2,x1,x2,wgt
  REAL*8 :: aa1,aa2,de,p,mm1,mm2,arat

  am = maxval(a)
  da(:) = (a(:)/10.+0.01)*rfac
  
  acc=0
  try=0
  eps=1.E-12
  DO i=1,steps

     ! computing current chi2
     xt1 = matmul(a , ker)
     x1 = sum(sxt*(axt-xt1)**2)
     
     ! number of Metroolis steps equal to the number of frequency points
     DO jj=1,2*nw+1

        do     ! need positive a(j2)+dj2
           j1=MIN(INT(rand()*(2*nw+1)),2*nw)-nw ! random number between -nw and nw
           do  ! need positive a(j1)+dj1
              dj1=da(j1)*(rand()-0.5)              ! dj1 -> -0.5*da(j1) ... 0.5*da(j1)
              if (a(j1)+dj1.GE.0.) exit
           enddo
           do  ! need j1!=j2
              j2=MIN(INT(rand()*(2*nw+1)),2*nw)-nw ! random number between -nw and nw
              IF (j1.NE.j2) exit
           enddo
           dj2=-dj1+0.05*da(j2)*(rand()-0.5)    ! dj2 -> -dj1 + -0.025*da(j2) ... 0.025*da(j2)
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
        IF (aa1.GT.eps)   de=de-f0(j1)*aa1*LOG(aa1/mm1)
        IF (a(j1).GT.eps) de=de+f0(j1)*a(j1)*LOG(a(j1)/mm1)
        IF (aa2.GT.eps)   de=de-f0(j2)*aa2*LOG(aa2/mm2)
        IF (a(j2).GT.eps) de=de+f0(j2)*a(j2)*LOG(a(j2)/mm2)            
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
        IF (rand().LT.p) THEN 
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
           WRITE(*,'(A,I7,2x,A,f12.5,2x,A,I10)') 'ann-step=', i,'chi2=', x1, 'accepted=', acc
        ENDIF

        temp=temp/1.5
        acc=0
        try=0
     ENDIF
  ENDDO
  
  RETURN
END SUBROUTINE MAXENT

!c------------------------------------------------------------------------------
REAL*8 FUNCTION ENTROPY(a,mm,f0,nw)
  IMPLICIT NONE
  REAL*8, intent(in)  :: a(-nw:nw), mm(-nw:nw), f0(-nw:nw)
  INTEGER, intent(in) :: nw
  ! locals
  INTEGER :: i
  REAL*8 :: ent, eps
  ent=0.
  eps=1.E-12
  DO i=-nw,nw  ! 10
     IF (a(i).GT.eps.AND.mm(i).GT.eps) THEN
        ent=ent-a(i)*LOG(a(i)/mm(i))*f0(i)
     ENDIF
  ENDDO
  ENTROPY=ent
  RETURN
END FUNCTION ENTROPY

!c------------------------------------------------------------------------------
REAL*8 FUNCTION LAMBDA(alpha,a,w,dlda,nw)
  IMPLICIT NONE
  !
  ! computes : Tr( lam * (lam+alpha)^-1 )
  ! where      lam = sqrt(a(w1))*ker(w1,:)*sxt(:)*ker(w2,:)*sqrt(a(w2))/dw
  !            
  REAL*8, intent(in):: alpha
  REAL*8, intent(in) :: a(-nw:nw), w(-nw:nw), dlda(2*nw+1,2*nw+1)
  INTEGER, intent(in):: nw
  ! locals
  REAL*8  :: dw, tr
  INTEGER :: j, i, iw, jw
  REAL*8  :: lam(2*nw+1,2*nw+1),alam(2*nw+1,2*nw+1),alami(2*nw+1,2*nw+1)
  
  dw=w(2)-w(1)
  DO j=1,2*nw+1
     DO i=1,2*nw+1
        iw=i-nw-1
        jw=j-nw-1
        lam(i,j) = SQRT(a(iw))*dlda(i,j)*SQRT(a(jw))*dw
        alam(i,j)=lam(i,j)
     ENDDO
     alam(j,j)=lam(j,j)+alpha
  ENDDO
  
  CALL INVERSE_MATRIX(alami, alam, 2*nw+1, 2*nw+1)
    
  tr=0.
  DO i=1,2*nw+1
     DO j=1,2*nw+1
        tr=tr+lam(i,j)*alami(j,i)
     ENDDO
  ENDDO
  LAMBDA = tr
  RETURN
END FUNCTION LAMBDA

!c------------------------------------------------------------------------------
SUBROUTINE INITDLDA(dlda, w, ker, sxt, nt, nw)
  IMPLICIT NONE
  REAL*8, intent(out) :: dlda(2*nw+1,2*nw+1)
  REAL*8, intent(in)  :: w(-nw:nw), ker(-nw:nw,nt), sxt(nt)
  INTEGER, intent(in) :: nw, nt
  ! locals
  INTEGER:: i,j,iw,jw,it
  REAL*8 :: dw, dsum
  !
  ! computing the following quantity:
  !    dlda(w1,w2)=ker(w1,:)*sxt(:)*ker(w2,:)/dw**2
  !
  dw=w(2)-w(1)
  DO j=1,2*nw+1
     DO i=1,2*nw+1
        iw=i-nw-1
        jw=j-nw-1
        dsum=0.0
        DO it=1,nt
           dsum = dsum + ker(iw,it)*sxt(it)*ker(jw,it)
        ENDDO
        dlda(i,j)=dsum/dw**2
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE INITDLDA

!c------------------------------------------------------------------------------
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
SUBROUTINE INITKER(ker,w,beta,t,nt,nw)
  IMPLICIT NONE
  REAL*8, intent(out) :: ker(-nw:nw,nt)
  REAL*8, intent(in)  :: w(-nw:nw), t(nt), beta
  INTEGER, intent(in) :: nw, nt
  ! locals
  INTEGER :: it,iw
  REAL*8  :: bw,tw,dw
  dw=w(2)-w(1)
  DO iw=-nw,nw
     DO it=1,nt
        tw=EXP(-t(it)*w(iw))
        if (w(iw)*beta<-50)then
           ker(iw,it)=dw*EXP(-t(it)*w(iw)+beta*w(iw))
        else 
           bw=1.+EXP(-beta*w(iw))
           ker(iw,it)=dw*tw/bw
        endif
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE INITKER
