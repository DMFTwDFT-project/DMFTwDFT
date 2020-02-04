      !program ksum_dmft
      include 'mpif.h'
      !implicit none
      REAL*8 tot_n,total_n,low_mu,high_mu,mix_Sig,Sig0,Sig1!,tot_n2,total_n2
      real*8 rdotk,kx,ky,kz,pi,mu,Ekin,Epot_st,Epot_dy,T,tot_Ekin,Nd
      integer qx,qy,qz,c1,c2,c3,nr,r,i,j,l,nk,ncor_orb,norb,nom,im
      integer,allocatable :: cor_idx(:),tran(:,:),IWORK(:),IPIV(:)
      integer x,y,z,rx,ry,rz,INFO,noms,iter,mu_iter,Niter,ispin,nspin
      complex*16,allocatable :: id(:,:),eval(:)
      complex*16 Gs,Gs0,DMFT_UU,DMFT_U0,ctmp
      complex*16,allocatable :: UMATRIX(:,:),hr(:,:,:),hamr(:,:,:)
      real*8,allocatable :: U(:),Up(:),UC(:,:,:),Jh(:),Sig(:),dm_ud(:)
      real*8,allocatable :: Sig_st(:),t_eval0(:,:,:),mom(:),tot_mom(:)
      real*8,allocatable :: new_Sig(:)
      real*8,allocatable :: write_Gloc(:,:,:,:)
      real*8,allocatable :: SigMdc(:,:,:),Mdc(:),eigvals(:,:)
      real*8 fermi, broaden, Gsr,n_tot,mix_mu,mu_conv
      real*8,allocatable :: dm(:),eval0(:),tot_dm(:),UM(:,:),Dval(:)
      complex*16,allocatable :: hk(:,:), evalue(:), t_eval(:,:,:,:)
      complex*16,allocatable :: Gsum(:,:)
      complex*16,allocatable :: Ginv(:,:),Gk(:,:,:)
      complex*16,allocatable :: scaler(:)
      complex*16,allocatable :: evec0(:,:),evl(:,:)
      complex*16,allocatable :: evr(:,:), tot_Gk(:,:,:)
      real*8,allocatable:: kweight(:),kpt_array(:), t_kpt_array(:)
      logical LOWB, HIGHB
      integer rank, size, ierr
      integer pr_proc,numk,ikp
      character(len=7) :: Ucfile    
      character(len=10) :: write_format   
      ! Wannier
      LOGICAL have_disentangled, iffile, lmatsubara
      character(len=33) :: header
      character(len=20) :: checkpoint
      integer :: num_bands, num_exclude_bands, num_kpts, nntot, num_wann
      integer :: num_band_max, nbmin, nbmax, nkpts
      real*8 :: real_latt(3,3), recip_latt(3,3),sweight, omega_invariant
      real*8 qq,qqx,qqy,qqz 
      real*8,allocatable:: qlist(:,:), kpt_latt(:,:)
      real*8,allocatable:: DMFT_w(:),t_DMFT_w(:)
      integer :: mp_grid(3), nb, nbb, IP, IPP, JJ, KK, nkp,num_tot_bands
      logical,allocatable:: lwindow(:,:), lexclude_band(:)
      integer,allocatable:: ndimwin(:), excl_bands(:), band_win(:,:)
      complex*16,allocatable:: u_matrix_opt(:,:,:), u_matrix(:,:,:)
      complex*16,allocatable:: UOPTU(:,:), UOPTU2(:,:)
      complex*16,allocatable:: DMFT_M(:,:)
      complex*16,allocatable:: DMFT_U(:),t_DMFT_U(:)

      call MPI_init(ierr)
      call MPI_Comm_size(MPI_Comm_World, size, ierr)
      call MPI_Comm_rank(MPI_Comm_World, rank, ierr)
      print *, 'Hello, World! I am process ',rank,' of ',size,'.'
      inquire(file='ksum.input',exist=iffile)
      if (iffile.eqv. .false.) then
        write(*,*) 'input is needed for ksum!'
        STOP
      else
        open(65,file="ksum.input")
        read(65,*) nkpts 
        read(65,*) nom 
        read(65,*) nspin
        read(65,*) ncor_orb
        read(65,*) norb
        !allocate(cor_idx(ncor_orb/norb))
        !read(65,*) (cor_idx(i),i=1,ncor_orb/norb)
        read(65,*) T
        read(65,*) n_tot
        read(65,*) mu_iter
        !allocate(U(ncor_orb/norb))
        !allocate(Up(ncor_orb/norb))
        !allocate(Jh(ncor_orb/norb))
        !read(65,*) (U(i),i=1,ncor_orb/norb)   
        !read(65,*) (Up(i),i=1,ncor_orb/norb)   
        !read(65,*) (Jh(i),i=1,ncor_orb/norb)   
        close(65)
      endif

      inquire(file='klist.dat',exist=iffile)
      if (iffile.eqv. .false.) then
        write(*,*) 'input is needed for ksum!'
        STOP
      else
        allocate(qlist(3,nkpts))
        open(65,file="klist.dat")
        do i=1,nkpts
          read(65,*) qq, qqx, qqy, qqz
          qlist(1,i)=qqx
          qlist(2,i)=qqy
          qlist(3,i)=qqz
        enddo
        close(65)
      endif
      !write(*,*) transpose(qlist)
      Niter=1
      inquire(file='DMFT_mu.out',exist=iffile)
      if (iffile.eqv. .false.) then
        write(*,*) 'input is needed for new mu!'
        STOP
      else
        open(65,file="DMFT_mu.out")
        read(65,*) mu
        close(65)
      endif
!      allocate(UC(ncor_orb/norb,2*norb,2*norb))
!      do i=1,ncor_orb/norb
! 112    FORMAT('UC',I1,'.dat')
!        WRITE(Ucfile,112) cor_idx(i)
!        inquire(file=Ucfile,exist=iffile)
!        if (iffile .eqv. .true.) then
!           open(unit=20,file=Ucfile,form='formatted')
!           do j=1,2*norb
!             read(20,*) (UC(i,j,l),l=1,2*norb)
!           enddo
!           close(20)
!        else
!           write(*,*) 'UC file must exist!!!'
!           STOP
!        endif
!      enddo
      sweight=1.0/nspin 
      allocate(SigMdc(2*ncor_orb+1,nom,nspin))
      inquire(file='SigMoo_real.out',exist=iffile)
      if (iffile .eqv. .true.) then
         open(unit=20,file='SigMoo_real.out',form='formatted')
         do j=1,nom
            read(20,*) (SigMdc(i,j,1),i=1,2*ncor_orb+1)
         enddo
         close(20)
      else
         write(*,*) 'SigMoo_real.out file must exist!!!'
         STOP
      endif
      if (nspin .eq. 2) then
         inquire(file='SigMoo_dn_real.out',exist=iffile)
         if (iffile .eqv. .true.) then
            open(unit=20,file='SigMoo_dn_real.out',form='formatted')
            do j=1,nom
               read(20,*) (SigMdc(i,j,2),i=1,2*ncor_orb+1)
            enddo
            close(20)
         else
            write(*,*) 'SigMoo_dn_real.out file must exist!!!'
            STOP
         endif
      endif

      allocate(Mdc(ncor_orb))
      inquire(file='SigMdc.out',exist=iffile)
      if (iffile .eqv. .true.) then
         open(unit=20,file='SigMdc.out',form='formatted')
         read(20,*) (Mdc(i),i=1,ncor_orb)
         close(20)
      else
         write(*,*) 'SigMdc.out file must exist!!!'
         STOP
      endif
      do i=1,nom
        do j=1,ncor_orb
          SigMdc(2*j,i,1)=SigMdc(2*j,i,1)+Mdc(j)
          !if (SigMdc(2*j+1,i).gt.0) SigMdc(2*j+1,i)=0.0
        enddo
      enddo
      if (nspin .eq. 2) then
         inquire(file='SigMdc_dn.out',exist=iffile)
         if (iffile .eqv. .true.) then
            open(unit=20,file='SigMdc_dn.out',form='formatted')
            read(20,*) (Mdc(i),i=1,ncor_orb)
            close(20)
         else
            write(*,*) 'SigMdc_dn.out file must exist!!!'
            STOP
         endif
         do i=1,nom
           do j=1,ncor_orb
             SigMdc(2*j,i,2)=SigMdc(2*j,i,2)+Mdc(j)
             !if (SigMdc(2*j+1,i).gt.0) SigMdc(2*j+1,i)=0.0
           enddo
         enddo
      endif
      deallocate(Mdc)

      !!!!!!!!!!!!!! READ WANNIER DATA !!!!!!!!!!!!!!!!!!!!!!!!!


      inquire(file='wannier90.chk',exist=iffile)
      if (iffile.eqv. .false.)then
         write(*,*) 'wannier90.chk must be present!!'
         STOP
      else
         open(unit=20,file='wannier90.chk',status='old',
     1           form='unformatted')
         read(20) header
         read(20) num_bands
         read(20) num_exclude_bands
         allocate(excl_bands(num_exclude_bands))
         read(20) (excl_bands(i),i=1,num_exclude_bands)
         read(20) ((real_latt(i,j),i=1,3),j=1,3)  ! Real lattice
         read(20) ((recip_latt(i,j),i=1,3),j=1,3)  ! Reciprocal lattice
         read(20) num_kpts                ! K-points
         read(20) (mp_grid(i),i=1,3)         ! M-P grid
         allocate(kpt_latt(3,num_kpts))
         read(20) ((kpt_latt(i,nkp),i=1,3),nkp=1,num_kpts)
         read(20) nntot                ! nntot
         read(20) num_wann                ! num_wann
         !if 
         read(20) checkpoint             ! checkpoint
         read(20) have_disentangled      ! whether a disentanglement has been performed
         if (have_disentangled) then
            read(20) omega_invariant     ! omega invariant
         ! lwindow
            allocate(lwindow(num_bands,num_kpts))
            read(20) ((lwindow(i,nkp),i=1,num_bands),nkp=1,num_kpts)
         ! ndimwin
            allocate(ndimwin(num_kpts))
            read(20) (ndimwin(nkp),nkp=1,num_kpts)
         ! U_matrix_opt
            allocate(u_matrix_opt(num_bands,num_wann,num_kpts))
            read(20) (((u_matrix_opt(i,j,nkp),i=1,num_bands),
     1               j=1,num_wann),nkp=1,num_kpts)
         else
            write(*,*) 'No U_matrix_opt ? Probably set identity'
            STOP
         endif
         ! U_matrix
         allocate(u_matrix(num_wann,num_wann,num_kpts))
         read(20) (((u_matrix(i,j,k),i=1,num_wann),j=1,num_wann)
     1            ,k=1,num_kpts)
         close(20)
      endif

      num_tot_bands=num_exclude_bands+num_bands
      allocate(lexclude_band(num_tot_bands))
      lexclude_band=.FALSE.
      DO nb=1,num_exclude_bands
         lexclude_band(excl_bands(nb))=.TRUE.
      ENDDO

      allocate(band_win(2,num_kpts))
      num_band_max=1
      DO nkp=1,num_kpts
        IP=0;IPP=0;nbmin=num_tot_bands;nbmax=1;
        DO nb=1,num_tot_bands
          IF (lexclude_band(nb)) CYCLE
          IP=IP+1
          IF (.NOT.lwindow(IP,nkp)) CYCLE
          IPP=IPP+1
          IF (nb<nbmin) nbmin=nb
          IF (nb>nbmax) nbmax=nb
        ENDDO
        band_win(1,nkp)=nbmin;band_win(2,nkp)=nbmax
        IF (num_band_max<nbmax-nbmin+1) num_band_max=nbmax-nbmin+1
      ENDDO

      allocate(eigvals(num_tot_bands,num_kpts))
      inquire(file='wannier90.eig',exist=iffile)
      if (iffile.eqv. .false.)then
         write(*,*) 'wannier90.eig must be present!!'
         STOP
      else
         open(unit=20,file='wannier90.eig',status='old',
     1           form='formatted')
      DO nkp=1,num_kpts
        DO nb=1,num_tot_bands
          read(20,*) x,y,eigvals(nb,nkp)
        ENDDO
      ENDDO
      endif
      !write(*,*) eigvals

!!!!!!!! K parallel !!!!!!!!!!!!!!!!
      pr_proc  = floor(num_kpts/DBLE(size)+0.999)
      numk = pr_proc
      if ((rank+1)*pr_proc .gt. num_kpts) then
         if (num_kpts-rank*pr_proc.gt.0) then
            numk = num_kpts-rank*pr_proc
         else
            numk = 0
         endif
      endif



      pi=3.14159265359;broaden=0.03
      !mu_iter=100
      !write(*,*) real_latt
      !STOP
      rx=mp_grid(1)/2;ry=mp_grid(2)/2;rz=mp_grid(3)/2
      nr=(2*rx+1)*(2*ry+1)*(2*rz+1)
      allocate(tran(3,nr))
      allocate(hk(num_wann,num_wann))
      allocate(UMATRIX(num_bands,num_wann))
      allocate(hr(nr,num_wann,num_wann))
      allocate(hamr(nr,num_wann,num_wann))
      hr=0
      hamr=0
      tran=0
      r=0
      do x=-rx,rx
        do y=-ry,ry
          do z=-rz,rz
            r=r+1
            tran(1,r)=x;tran(2,r)=y;tran(3,r)=z
          enddo
        enddo
      enddo
      do nkp=1,num_kpts
!      do c1=1,qx
!        do c2=1,qy
!          do c3=qz/2+1,qz
            !ikp = c3-qz/2+(c2-qy/2-1)*(qz/2)+(c1-qx/2-1)*(qy/2)*(qz/2)
!            ikp = c3-qz/2+(c2-1)*(qz/2)+(c1-1)*(qy)*(qz/2)
!     1             -rank*pr_proc
        ikp = nkp-rank*pr_proc
        if (ikp.GT.numk) EXIT
        if (ikp.gt.0) then
!          if (rank.eq.0) then
!            WRITE(199,*) kpt_latt(1,nkp), kpt_latt(2,nkp), kpt_latt(3,nkp)
!          kx=2*pi*(2*c1-qx-1)/(2*qx)
!          ky=2*pi*(2*c2-qy-1)/(2*qy)
!          kz=2*pi*(2*c3-qz-1)/(2*qz)
          ! create H(r)...
          !write(*,*) cmplx(0.0,rdotk)
          !write(*,*) kpt_latt(:,nkp) 
          nbmin=band_win(1,nkp)
          nbmax=band_win(2,nkp)
          num_band_max=nbmax-nbmin+1
          hk=0
          UMATRIX=0
          UMATRIX=MATMUL(u_matrix_opt(:,:,nkp),u_matrix(:,:,nkp))
          do x=1,num_wann
            do y=1,num_wann
              do z=1,num_band_max
                !write(*,*) UMATRIX(z,x)
                !write(*,*) eigvals(nbmin+z-1,nkp)
                hk(x,y)=hk(x,y)+conjg(UMATRIX(z,x))*UMATRIX(z,y)
     1                  *eigvals(nbmin+z-1,nkp)
              enddo
            enddo
          enddo
          do r=1,nr
            !write(*,*) dfloat(tran(:,r))
            rdotk=-2*pi*dot_product(kpt_latt(:,nkp),dfloat(tran(:,r)))
            hr(r,:,:)=hr(r,:,:)+hk*exp(dcmplx(0.0,rdotk))
            !write(*,*) exp(dcmplx(0.0,rdotk))
          enddo
          !! create H(k)...
          !hk=0
          !do r=1,nr
          !  hk=hk+ham(r,:,:)*exp(
     1    !  (0.0,1.0)*(kx*tran(r,1)+ky*tran(r,2)+kz*tran(r,3)))
          !enddo
        endif 
      enddo
      CALL MPI_ALLREDUCE(hr,hamr,nr*num_wann*num_wann,MPI_DOUBLE_COMPLEX
     1      ,MPI_SUM,MPI_Comm_World,ierr)
      hamr=hamr/num_kpts
      !if (rank.eq.0) write(*,*) hamr(63,1,:)
      !write(*,*) hamr
      deallocate(eigvals)
      deallocate(hr)
      deallocate(UMATRIX)
      deallocate(hk)
      deallocate(u_matrix)
      deallocate(u_matrix_opt)
      deallocate(kpt_latt)

      pr_proc  = floor(nkpts/DBLE(size)+0.999)
      numk = pr_proc
      if ((rank+1)*pr_proc .gt. nkpts) then
         if (nkpts-rank*pr_proc.gt.0) then
            numk = nkpts-rank*pr_proc
         else
            numk = 0
         endif
      endif
                  



      allocate(hk(num_wann,num_wann))
      allocate(Ginv(num_wann,num_wann))
      allocate(Gk(nom,nkpts,nspin))
      allocate(tot_Gk(nom,nkpts,nspin))
      allocate(IPIV(num_wann))
      allocate(id(num_wann,num_wann))

      Gk=0;IPIV=0
      do nkp=1,nkpts
        ikp = nkp-rank*pr_proc
        if (ikp.GT.numk) EXIT
        if (ikp.gt.0) then
          hk=0
          do r=1,nr
            rdotk=2*pi*dot_product(qlist(:,nkp),dfloat(tran(:,r)))
            hk=hk+hamr(r,:,:)*exp(dcmplx(0.0,rdotk))
          enddo
          do ispin=1,nspin
            do im=1,nom
              id=0; Ginv=0;
              do j=1,num_wann; id(j,j)=1.0; enddo
              Ginv=hk*(-1.0)
              do j=1,num_wann
                Ginv(j,j)=Ginv(j,j)+mu+dcmplx(SigMdc(1,im,ispin)
     1          ,broaden)
                !Ginv(j,j)=Ginv(j,j)+mu+ommesh(i)+(0.0,1.0)*broaden
              enddo
              do j=1,ncor_orb
                Ginv(j,j)=Ginv(j,j)-dcmplx(SigMdc(2*j,im,ispin)
     1          ,SigMdc(2*j+1,im,ispin))
              enddo
              call ZGESV(num_wann,num_wann,Ginv,num_wann,IPIV,id,
     1             num_wann,INFO)
              do i=1,num_wann
                Gk(im,nkp,ispin)=Gk(im,nkp,ispin)+id(i,i)
              enddo
            enddo
          enddo
        endif    
      enddo
      tot_Gk=0;
      CALL MPI_REDUCE(Gk,tot_Gk,nom*nkpts*nspin,
     1      MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_Comm_World,ierr)
      deallocate(hk)
      deallocate(Ginv)
      deallocate(Gk)
      deallocate(IPIV)
      deallocate(id)
      if (rank.eq.0) then
        !write(*,*) 2*num_wann+1
        !allocate(write_Gloc(2*num_wann+1,nom,nkpts,nspin))
        allocate(write_Gloc(3,nom,nkpts,nspin))
        do ispin=1,nspin
          do ik=1,nkpts
            do im=1,nom
              write_Gloc(1,im,ik,ispin)=SigMdc(1,im,ispin)
              !do i=1,num_wann !ncor_orb
              write_Gloc(2,im,ik,ispin)=real(tot_Gk(im,ik,ispin))
              write_Gloc(3,im,ik,ispin)=aimag(tot_Gk(im,ik,ispin
     1          ))
              !enddo
            enddo
          enddo
        enddo
        OPEN(UNIT=99,FILE='Gk.out',FORM='FORMATTED')
        !write(*,*) write_format
        DO ik=1,nkpts
          WRITE(99,*) 'k=',qlist(1,ik),qlist(2,ik),qlist(3,ik) 
        DO im=1,nom
          WRITE(99,'(3F18.12)') (write_Gloc(i,im,ik,1),i=1,3)
        ENDDO
        ENDDO
        CLOSE(99)
        if (nspin.eq.2) then
          OPEN(UNIT=99,FILE='Gk_dn.out',FORM='FORMATTED')
          !write(*,*) write_format
          DO ik=1,nkpts
          WRITE(99,*) 'k=',qlist(1,ik),qlist(2,ik),qlist(3,ik) 
          DO im=1,nom
            WRITE(99,'(3F18.12)') (write_Gloc(i,im,ik,2),i=1,3)
          ENDDO
          ENDDO
          CLOSE(99)
        endif
        deallocate(write_Gloc)

      endif

      deallocate(qlist)
      deallocate(tot_Gk)
      deallocate(hamr)
      deallocate(tran)
      deallocate(SigMdc)
      !STOP



      CALL MPI_Finalize(ierr)

      end             

      SUBROUTINE EIGVALH(ndim,evector,evalue)

      INTEGER ndim, i, j, IWORK(3+5*ndim), INFO
      REAL*8 evalue(ndim), RWORK(1+5*ndim+2*ndim**2)
      COMPLEX*16 evector(ndim,ndim), WORK(2*ndim+ndim**2)

      call ZHEEVD('N','U',ndim,evector,ndim,evalue,WORK,2*ndim+ndim**2,
     1 RWORK,1+5*ndim+2*ndim**2,IWORK,3+5*ndim,INFO)
      
      END SUBROUTINE

      SUBROUTINE EIGVAL(ndim,ComplexM,evalue,evl,evr)

      INTEGER ndim, i, j, IWORK(3+5*ndim), INFO
      REAL*8 RWORK(1+5*ndim+2*ndim**2)
      COMPLEX*16 evalue(ndim), ComplexM(ndim,ndim)
      COMPLEX*16 evl(ndim,ndim), evr(ndim,ndim), WORK(2*ndim+ndim**2)

      call ZGEEV('N','N',ndim,ComplexM,ndim,evalue,evl,ndim,evr,ndim,
     1 WORK,2*ndim+ndim**2,RWORK,INFO)

      END SUBROUTINE

      SUBROUTINE EIGH(ndim,evector,evalue)

      INTEGER ndim, i, j, IWORK(3+5*ndim), INFO
      REAL*8 ctmp, smalleps, evalue(ndim), RWORK(1+5*ndim+2*ndim**2)
      COMPLEX*16 scaler(ndim), evector(ndim,ndim), WORK(2*ndim+ndim**2)
      complex*16 scalprod

      smalleps= 1e-5
      call ZHEEVD('V','U',ndim,evector,ndim,evalue,WORK,2*ndim+ndim**2,
     1 RWORK,1+5*ndim+2*ndim**2,IWORK,3+5*ndim,INFO)

      !========== Make degenerate eigenvectors orthogonal
      DO i=1,ndim
        DO j=1,i-1
          IF (abs(evalue(j)-evalue(i)).LT.smalleps .AND.
     1 abs(scalprod(conjg(evector(:,j)),evector(:,i),ndim)) 
     1 .GT.smalleps) THEN
            evector(:,i) = evector(:,i) -
     1      scalprod(conjg(evector(:,j)),evector(:,i),ndim)/
     1      scalprod(conjg(evector(:,j)),evector(:,j),ndim)*evector(:,j)
          ENDIF
        ENDDO
      ENDDO
      !========= Normalize eigenvectors
      DO i = 1,ndim
         ctmp = 0.0
         DO j = 1,ndim
            ctmp = ctmp+conjg(evector(j,i))*evector(j,i)
         ENDDO
         scaler(i) = SQRT(ctmp)
      ENDDO
      DO i = 1,ndim
         evector(:,i) = evector(:,i)/scaler(i)
      ENDDO

      END SUBROUTINE


      SUBROUTINE EIG(ndim,ComplexM,evalue,evl,evr)

      INTEGER ndim, i, j, IWORK(3+5*ndim), INFO
      REAL*8 smalleps, RWORK(1+5*ndim+2*ndim**2)
      COMPLEX*16 ctmp, scaler(ndim), evalue(ndim), ComplexM(ndim,ndim),
     1 evl(ndim,ndim), evr(ndim,ndim), WORK(2*ndim+ndim**2)
      complex*16 scalprod

      smalleps= 1e-5
      call ZGEEV('V','V',ndim,ComplexM,ndim,evalue,evl,ndim,evr,ndim,
     1 WORK,2*ndim+ndim**2,RWORK,INFO)
      evl = conjg(TRANSPOSE(evl))

      !========== Step 5, Make degenerate eigenvectors orthogonal
      DO i=1,ndim
        DO j=1,i-1
          IF (abs(evalue(j)-evalue(i)).LT.smalleps .AND.
     1       abs(scalprod(evl(j,:),evr(:,i),ndim)).GT.smalleps) THEN
            evr(:,i) = evr(:,i) -
     1      scalprod(evl(j,:),evr(:,i),ndim)/
     1      scalprod(evl(j,:),evr(:,j),ndim)*evr(:,j)
          ENDIF
        ENDDO
        DO j=1,i-1
          IF (abs(evalue(j)-evalue(i)).LT.smalleps .AND.
     1        abs(scalprod(evl(i,:),evr(:,j),ndim)).GT.smalleps) THEN
            evl(i,:) = evl(i,:) -
     1      scalprod(evl(i,:),evr(:,j),ndim)/
     1      scalprod(evl(j,:),evr(:,j),ndim)*evl(j,:)
          ENDIF
        ENDDO
      ENDDO
        !========= Step 6, Normalize eigenvectors
      DO i = 1,ndim
         ctmp = 0.d0
         DO j = 1,ndim
            ctmp = ctmp+evl(i,j)*evr(j,i)
         ENDDO
         scaler(i) = SQRT(ctmp)
      ENDDO
      DO i = 1,ndim
         evl(i,:) = evl(i,:)/scaler(i)
         evr(:,i) = evr(:,i)/scaler(i)
      ENDDO

      END SUBROUTINE

      COMPLEX*16 FUNCTION scalprod(a,b,norb)
      IMPLICIT NONE
      COMPLEX*16 a(norb), b(norb)
      INTEGER    norb
      INTEGER    i
      scalprod = 0.0
      DO i=1,norb
         scalprod = scalprod + a(i)*b(i)
      ENDDO
      END FUNCTION scalprod


      real*8 function fermi(eps)
      implicit none
      real*8 eps
      
      if (eps>100.0) then 
        fermi=0.0
        return
      else if (eps<-100) then
        fermi=1.0
        return
      else 
        fermi=1.0/(1.0+exp(eps))
        return
      endif
      end function fermi
