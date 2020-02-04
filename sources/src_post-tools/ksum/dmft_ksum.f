      !program ksum_dmft
      include 'mpif.h'
      !implicit none
      REAL*8 tot_n,total_n,low_mu,high_mu,mix_Sig,Sig0,Sig1!,tot_n2,total_n2
      real*8 rdotk,kx,ky,kz,pi,mu,Ekin,Epot_st,Epot_dy,T,tot_Ekin,Nd
      integer qx,qy,qz,c1,c2,c3,nr,r,i,j,l,nk,ncor_orb,norb,nom,im
      integer,allocatable :: cor_idx(:),tran(:,:),IWORK(:),IPIV(:)
      integer x,y,z,rx,ry,rz,INFO,noms,iter,mu_iter,Niter
      complex*16,allocatable :: id(:,:),eval(:)
      complex*16 Gs,Gs0,DMFT_UU,DMFT_U0,ctmp
      complex*16,allocatable :: UMATRIX(:,:),hr(:,:,:),hamr(:,:,:)
      real*8,allocatable :: U(:),Up(:),UC(:,:,:),Jh(:),Sig(:),old_Sig(:)
      real*8,allocatable :: Sig_st(:)
      real*8,allocatable :: Ed(:),tot_Ed(:),new_Sig(:),write_Gloc(:,:)
      real*8,allocatable :: SigMdc(:,:),Mdc(:),eigvals(:,:),t_eval0(:,:)
      real*8 fermi, broaden, Gsr,n_tot,mix_mu,mu_conv
      real*8,allocatable :: dm(:),eval0(:),tot_dm(:),UM(:,:),Dval(:)
      real*8,allocatable :: dm_full(:,:), tot_dm_full(:,:)
      complex*16,allocatable :: hk(:,:), evalue(:), t_eval(:,:,:)
      complex*16,allocatable :: Gsum(:,:)
      complex*16,allocatable :: Ginv(:,:),Gloc(:,:)
      complex*16,allocatable :: scaler(:)
      complex*16,allocatable :: evec0(:,:),evl(:,:)
      complex*16,allocatable :: evr(:,:), tot_Gloc(:,:)
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
      integer :: num_band_max, nbmin, nbmax
      real*8 :: real_latt(3,3), recip_latt(3,3), omega_invariant
      real*8,allocatable:: kpt_latt(:,:)
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
        read(65,*) Niter 
        read(65,*) qx, qy, qz 
        read(65,*) noms, nom 
        read(65,*) ncor_orb
        read(65,*) norb
        allocate(cor_idx(ncor_orb/norb))
        read(65,*) (cor_idx(i),i=1,ncor_orb/norb)
        read(65,*) T
        read(65,*) n_tot
        read(65,*) mix_Sig
        allocate(U(ncor_orb/norb))
        allocate(Up(ncor_orb/norb))
        allocate(Jh(ncor_orb/norb))
        read(65,*) (U(i),i=1,ncor_orb/norb)   
        read(65,*) (Up(i),i=1,ncor_orb/norb)   
        read(65,*) (Jh(i),i=1,ncor_orb/norb)   
        close(65)
      endif
      inquire(file='DMFT_mu.out',exist=iffile)
      if (iffile.eqv. .false.) then
        write(*,*) 'input is needed for new mu!'
        STOP
      else
        open(65,file="DMFT_mu.out")
        read(65,*) mu
        close(65)
      endif
      allocate(UC(ncor_orb/norb,2*norb,2*norb))
      do i=1,ncor_orb/norb
 112    FORMAT('UC',I1,'.dat')
        WRITE(Ucfile,112) cor_idx(i)
        inquire(file=Ucfile,exist=iffile)
        if (iffile .eqv. .true.) then
           open(unit=20,file=Ucfile,form='formatted')
           do j=1,2*norb
             read(20,*) (UC(i,j,l),l=1,2*norb)
           enddo
           close(20)
        else
           write(*,*) 'UC file must exist!!!'
           STOP
        endif
      enddo
      
      allocate(SigMdc(2*ncor_orb+1,nom))
      inquire(file='SigMoo.out',exist=iffile)
      if (iffile .eqv. .true.) then
         open(unit=20,file='SigMoo.out',form='formatted')
         do j=1,nom
            read(20,*) (SigMdc(i,j),i=1,2*ncor_orb+1)
         enddo
         close(20)
      else
         write(*,*) 'SigMoo.out file must exist!!!'
         STOP
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
          SigMdc(2*j,i)=SigMdc(2*j,i)+Mdc(j)
          !if (SigMdc(2*j+1,i).gt.0) SigMdc(2*j+1,i)=0.0
        enddo
      enddo
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



      pi=3.14159265359!;smalleps = 1e-5;broaden=0.0
      mu_iter=100
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

      num_kpts=qx*qx*(qz/2+1)
      allocate(kpt_latt(3,num_kpts))
      allocate(kweight(num_kpts))
      nkp=0
      if (MOD(qz,2) .eq. 0) then
        do c1=1,qx
          do c2=1,qy
            do c3=1,qz/2+1
              nkp=nkp+1
              kpt_latt(1,nkp)=dfloat(c1-qx/2-1)/dfloat(qx)
              kpt_latt(2,nkp)=dfloat(c2-qy/2-1)/dfloat(qy)
              kpt_latt(3,nkp)=dfloat(c3-qz/2-1)/dfloat(qz)
              if (c3.eq.1 .OR. c3.eq.(qz/2+1)) then
                kweight(nkp)=1.0/dfloat(qx*qy*qz)
              else
                kweight(nkp)=2.0/dfloat(qx*qy*qz)
              endif
            enddo
          enddo
        enddo
      else
        do c1=1,qx
          do c2=1,qy
            do c3=1,qz/2+1
              nkp=nkp+1
              kpt_latt(1,nkp)=dfloat(c1-qx/2-1)/dfloat(qx)
              kpt_latt(2,nkp)=dfloat(c2-qy/2-1)/dfloat(qy)
              kpt_latt(3,nkp)=dfloat(c3-qz/2-1)/dfloat(qz)
              if (c3.eq.(qz/2+1)) then
                kweight(nkp)=1.0/dfloat(qx*qy*qz)
              else
                kweight(nkp)=2.0/dfloat(qx*qy*qz)
              endif
            enddo
          enddo
        enddo
      endif
      !write(*,*) kpt_latt
      !write(*,*) kweight
      pr_proc  = floor(num_kpts/DBLE(size)+0.999)
      numk = pr_proc
      if ((rank+1)*pr_proc .gt. num_kpts) then
         if (num_kpts-rank*pr_proc.gt.0) then
            numk = num_kpts-rank*pr_proc
         else
            numk = 0
         endif
      endif
                  

      allocate(Sig_st(ncor_orb))
      allocate(tot_Ed(ncor_orb))
      allocate(tot_dm(num_wann))
      allocate(tot_dm_full(num_wann,num_wann))
      allocate(tot_Gloc(num_wann,nom))
      do iter=1,Niter 
        allocate(hk(num_wann,num_wann))
        allocate(evec0(num_wann,num_wann))
        allocate(eval0(num_wann))
        allocate(t_eval0(num_wann,numk))
        allocate(t_eval(num_wann,noms,numk))
        allocate(evalue(num_wann))
        allocate(evl(num_wann,num_wann))
        allocate(evr(num_wann,num_wann))
        do nkp=1,num_kpts
          ikp = nkp-rank*pr_proc
          if (ikp.GT.numk) EXIT
          if (ikp.gt.0) then
            hk=0
            do r=1,nr
              !write(*,*) dfloat(tran(:,r))
              !write(*,*) 'hi'
              rdotk=2*pi*dot_product(kpt_latt(:,nkp),dfloat(tran(:,r)))
              hk=hk+hamr(r,:,:)*exp(dcmplx(0.0,rdotk))
            enddo    
            !!!!!!!! Compute the infinity value !!!!
            evec0=0;eval0=0
            evec0=hk*(1.0)
            do i=1,ncor_orb
              evec0(i,i)=evec0(i,i)+SigMdc(2*i,nom)
            enddo
            CALL EIGVALH(num_wann,evec0,eval0) 
            do nb=1,num_wann
              t_eval0(nb,ikp)=eval0(nb)
            enddo
            do im=1,noms
              evec0=0;evalue=0
              evec0=hk*(1.0) 
              do i=1,ncor_orb
                evec0(i,i)=evec0(i,i)
     1          +dcmplx(SigMdc(2*i,im),SigMdc(2*i+1,im))
              enddo
              CALL EIGVAL(num_wann,evec0,evalue,evl,evr)
              do nb=1,num_wann
                t_eval(nb,im,ikp)=evalue(nb)
              enddo
            enddo
          endif
        enddo
        deallocate(hk)
        deallocate(evec0)
        deallocate(eval0)
        deallocate(evalue)
        deallocate(evl)
        deallocate(evr)
        low_mu=mu;high_mu=mu
        LOWB=.FALSE.;HIGHB=.FALSE.
        do l=1,mu_iter 
          tot_n=0
          mu=(low_mu+high_mu)/2
          do nkp=1,num_kpts
            ikp = nkp-rank*pr_proc
            if (ikp.GT.numk) EXIT
            if (ikp.gt.0) then
              do nb=1,num_wann
                tot_n=tot_n+2*fermi((t_eval0(nb,ikp)-mu)/T)*kweight(nkp)
                do im=1,noms
                  Gs0=1.0/(mu+dcmplx(0.0,SigMdc(1,im))-t_eval0(nb,ikp))
                  Gs=1.0/(mu+dcmplx(0.0,SigMdc(1,im))-t_eval(nb,im,ikp))
                  tot_n=tot_n+4*T*real(Gs-Gs0)*kweight(nkp)
                enddo
              enddo
            endif
          enddo
          total_n=0
          CALL MPI_ALLREDUCE(tot_n,total_n,1,MPI_DOUBLE_PRECISION,
     1          MPI_SUM,MPI_Comm_World,ierr)
          if (abs(total_n-n_tot)<1E-10) EXIT
          if (total_n<n_tot) then
            if (HIGHB.eqv..FALSE.) then
              high_mu=mu+(n_tot-total_n)*2/(0.1*n_tot)
            endif
            low_mu=mu
            LOWB=.TRUE.
          else
            if (LOWB.eqv..FALSE.) then
              low_mu=mu+(n_tot-total_n)*2/(0.1*n_tot)
            endif
            high_mu=mu
            HIGHB=.TRUE.
          endif
          !if (rank.eq.0) write(*,*) mu, total_n, n_tot
          if (l.eq.mu_iter) then
            write(*,*) 'Fail to find mu, increase mu_iter'
            !STOP
          endif
        enddo

        deallocate(t_eval0)
        deallocate(t_eval)

        allocate(Ed(ncor_orb))
        allocate(dm(num_wann))
        allocate(dm_full(num_wann,num_wann))
        allocate(hk(num_wann,num_wann))
        allocate(Gsum(num_wann,num_wann))
        allocate(Ginv(num_wann,num_wann))
        allocate(Gloc(num_wann,nom))
        allocate(IPIV(num_wann))
        allocate(id(num_wann,num_wann))
        allocate(evec0(num_wann,num_wann))
        allocate(eval0(num_wann))
        allocate(evalue(num_wann))
        allocate(evl(num_wann,num_wann))
        allocate(evr(num_wann,num_wann))

        Ed=0;Gloc=0;Ekin=0;dm=0;dm_full=0;IPIV=0
        do nkp=1,num_kpts
          ikp = nkp-rank*pr_proc
          if (ikp.GT.numk) EXIT
          if (ikp.gt.0) then
            hk=0
            do r=1,nr
              rdotk=2*pi*dot_product(kpt_latt(:,nkp),dfloat(tran(:,r)))
              hk=hk+hamr(r,:,:)*exp(dcmplx(0.0,rdotk))
            enddo
            do i=1,ncor_orb
              Ed(i)=Ed(i)+real(hk(i,i))*kweight(nkp)
            enddo
            !!!!!!!! Compute the infinity value !!!!
            evec0=0;eval0=0
            evec0=hk*(1.0)
            do i=1,ncor_orb
              evec0(i,i)=evec0(i,i)+SigMdc(2*i,nom)
            enddo
            CALL EIGH(num_wann,evec0,eval0)
            Gsum=0
            do im=1,noms
              Ginv=0;evalue=0
              Ginv=hk*(1.0)
              do i=1,ncor_orb
                Ginv(i,i)=Ginv(i,i)
     1          +dcmplx(SigMdc(2*i,im),SigMdc(2*i+1,im))
              enddo
              CALL EIG(num_wann,Ginv,evalue,evl,evr)
              do nb=1,num_wann
                Gs0=1.0/(mu+dcmplx(0.0,SigMdc(1,im))-eval0(nb))
                Gs=1.0/(mu+dcmplx(0.0,SigMdc(1,im))-evalue(nb))
                do i=1,num_wann
                  do j=1,num_wann
                    DMFT_UU=evr(i,nb)*evl(nb,j)*Gs
                    if (i.eq.j) Gloc(i,im)=Gloc(i,im)+DMFT_UU
     1                          *kweight(nkp)
                    DMFT_U0=evec0(i,nb)*conjg(evec0(j,nb))*Gs0
                    Gsum(i,j)=Gsum(i,j)+T*(DMFT_UU-DMFT_U0)*kweight(nkp)
                    Gsum(j,i)=Gsum(j,i)+T*conjg(DMFT_UU-DMFT_U0)*
     1                        kweight(nkp)
                  enddo
                enddo
              enddo
            enddo
            do nb=1,num_wann
              Gsr=fermi((eval0(nb)-mu)/T)
              do i=1,num_wann
                do j=1,num_wann
                  Gsum(i,j)=Gsum(i,j)+Gsr*evec0(i,nb)*conjg(evec0(j,nb))
     1                      *kweight(nkp)
                enddo
              enddo
            enddo
            do i=1,num_wann
              dm(i)=dm(i)+real(Gsum(i,i))
              do j=1,num_wann
                !dm_full(i,j)=dm_full(i,j)+real(Gsum(i,j))
                dm_full(i,j)=dm_full(i,j)+aimag(Gsum(i,j))
                Ekin=Ekin+real(hk(i,j)*Gsum(j,i))
              enddo
            enddo

            do im=noms+1,nom
              id=0; Ginv=0;
              do j=1,num_wann; id(j,j)=1.0; enddo
              Ginv=hk*(-1.0)
              do j=1,num_wann
                Ginv(j,j)=Ginv(j,j)+mu+dcmplx(0.0,SigMdc(1,im))
                !Ginv(j,j)=Ginv(j,j)+mu+ommesh(i)+(0.0,1.0)*broaden
              enddo
              do j=1,ncor_orb
                Ginv(j,j)=Ginv(j,j)-dcmplx(SigMdc(2*j,im),SigMdc(2*j+1,
     1          im))
              enddo
              call ZGESV(num_wann,num_wann,Ginv,num_wann,IPIV,id,
     1             num_wann,INFO)
              do i=1,num_wann
                Gloc(i,im)=Gloc(i,im)+id(i,i)*kweight(nkp)
              enddo
            enddo
          endif    
        enddo
        tot_Gloc=0;tot_dm=0;tot_Ekin=0;tot_Ed=0
        CALL MPI_REDUCE(Gloc,tot_Gloc,num_wann*nom,MPI_DOUBLE_COMPLEX,
     1        MPI_SUM,0,MPI_Comm_World,ierr)
        CALL MPI_ALLREDUCE(Ekin,tot_Ekin,1,MPI_DOUBLE_PRECISION,
     1        MPI_SUM,MPI_Comm_World,ierr)
        CALL MPI_ALLREDUCE(Ed,tot_Ed,ncor_orb,MPI_DOUBLE_PRECISION,
     1        MPI_SUM,MPI_Comm_World,ierr)
        CALL MPI_ALLREDUCE(dm,tot_dm,num_wann,MPI_DOUBLE_PRECISION,
     1        MPI_SUM,MPI_Comm_World,ierr)
        CALL MPI_ALLREDUCE(dm_full,tot_dm_full,num_wann**2,
     1        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_Comm_World,ierr)
        allocate(Sig(ncor_orb))
        allocate(old_Sig(ncor_orb))
        allocate(new_Sig(ncor_orb))
        Epot_st=0;Epot_dy=0;Sig_st=0;Sig=0;old_Sig=0;new_Sig=0;
        do i=1,ncor_orb/norb
          Nd=0
          do j=1,2*norb
            !UC(i,j,j)=UC(i,j,j)-U(i)
            Nd=Nd+tot_dm((i-1)*norb+MOD(j-1,norb)+1)
          enddo
          !write(*,*) Nd
          do j=1,norb
            do l=1,2*norb
              if (j.eq.l) then
                Sig((i-1)*norb+j)=Sig((i-1)*norb+j)+UC(i,j,l)
     1          *tot_dm((i-1)*norb+MOD(l-1,norb)+1)
              else
                Sig((i-1)*norb+j)=Sig((i-1)*norb+j)+(UC(i,j,l)+U(i))
     1          *tot_dm((i-1)*norb+MOD(l-1,norb)+1)
              endif
            enddo
            Epot_st=Epot_st+Sig((i-1)*norb+j)*tot_dm((i-1)*norb+j)
            Sig((i-1)*norb+j)=Sig((i-1)*norb+j)
     1                     -(Up(i)*(Nd-0.5)-Jh(i)*(Nd-1)/2)
!            Sig_st((i-1)*norb+j)=Sig_st((i-1)*norb+j)
!     1                     -(Up(i)*(Nd-0.5)-Jh(i)*(Nd-1)/2)
          enddo
          Epot_st=Epot_st-Up(i)*Nd*(Nd-1.0)/2.0+Jh(i)*Nd*(Nd-2.0)/4.0
        enddo
        do i=1,ncor_orb
          old_Sig(i)=SigMdc(2*i,nom)
        enddo
        new_Sig=old_Sig+mix_Sig*(Sig-old_Sig)
        !if (rank.eq.0) write(*,*) new_Sig(1:5) 
        do i=1,ncor_orb
          new_Sig(i)=new_Sig(i)-old_Sig(i) 
          !write(*,*) new_Sig(i)
          do im=1,nom
            SigMdc(2*i,im)=SigMdc(2*i,im)+new_Sig(i)
          enddo
        enddo
        tot_dm=tot_dm*2;tot_Ekin=tot_Ekin*2
        !if (rank.eq.0) write(*,*) tot_dm(1:10)!, tot_Ekin
        

        deallocate(Sig)
        deallocate(old_Sig)
        deallocate(new_Sig)
        deallocate(Ed)
        deallocate(dm)
        deallocate(dm_full)
        deallocate(hk)
        deallocate(Gsum)
        deallocate(Ginv)
        deallocate(Gloc)
        deallocate(IPIV)
        deallocate(id)
        deallocate(evec0)
        deallocate(eval0)
        deallocate(evalue)
        deallocate(evl)
        deallocate(evr)
      enddo
      if (rank.eq.0) then
        do i=1,ncor_orb
          Sig0=SigMdc(2*i,noms);Sig1=SigMdc(2*i+1,noms)*SigMdc(1,noms)
          do j=1,nom
            Epot_dy=Epot_dy+2.0*T*real(tot_Gloc(i,j)*(dcmplx(SigMdc(2*i
     1      ,j)-Sig0,SigMdc(2*i+1,j)))-Sig1/(SigMdc(1,j)**2))
          enddo
          Epot_dy=Epot_dy+Sig1/(4*T)
        enddo
        allocate(write_Gloc(2*ncor_orb+1,nom))
        do im=1,nom
          write_Gloc(1,im)=SigMdc(1,im)
          do i=1,ncor_orb
            write_Gloc(2*i,im)=real(tot_Gloc(i,im))
            write_Gloc(2*i+1,im)=aimag(tot_Gloc(i,im))
          enddo
        enddo
        OPEN(UNIT=99,FILE='G_loc.out',FORM='FORMATTED')
 110    FORMAT('(',I2,'F18.12)')
        WRITE(write_format,110) 2*ncor_orb+1
        !write(*,*) write_format
        DO im=1,nom
          WRITE(99,write_format) (write_Gloc(i,im),i=1,2*ncor_orb+1)
        ENDDO
        CLOSE(99)
        deallocate(write_Gloc)
        OPEN(UNIT=99,FILE='SigMdc.out',FORM='FORMATTED')
 120    FORMAT('(',I2,'F18.12)')
        WRITE(write_format,120) ncor_orb
        !write(*,*) write_format
        WRITE(99,write_format) (SigMdc(2*i,nom),i=1,ncor_orb)
        CLOSE(99)
        OPEN(UNIT=99,FILE='Sig_st.out',FORM='FORMATTED')
 122    FORMAT('(',I2,'F18.12)')
        WRITE(write_format,122) ncor_orb
        !write(*,*) write_format
        WRITE(99,write_format) (Sig_st(i),i=1,ncor_orb)
        CLOSE(99)
        OPEN(UNIT=99,FILE='Ed.out',FORM='FORMATTED')
 130    FORMAT('(',I2,'F18.12)')
        WRITE(write_format,130) ncor_orb
        WRITE(99,write_format) (tot_Ed(i),i=1,ncor_orb)
        CLOSE(99)
        OPEN(UNIT=99,FILE='DMINFO',STATUS='old',
     1   FORM='FORMATTED',ACCESS='APPEND')
 132    FORMAT('(',I2,'F18.12)')
        WRITE(write_format,132) ncor_orb
        WRITE(99,write_format) (tot_dm(i),i=1,ncor_orb)
        CLOSE(99)
        OPEN(UNIT=99,FILE='FULLDMINFO_IMAG',
     1   FORM='FORMATTED')
 134    FORMAT('(',I2,'F18.12)')
        WRITE(write_format,134) ncor_orb
        do i=1,ncor_orb
           WRITE(99,write_format) (tot_dm_full(i,j),j=1,ncor_orb)
        enddo
        CLOSE(99)
        OPEN(UNIT=99,FILE='DMFT_mu.out',FORM='FORMATTED')
        WRITE(99,'(F18.12)') mu 
        CLOSE(99)
        OPEN(UNIT=99,FILE='ITERINFO',STATUS='old',
     1   FORM='FORMATTED',ACCESS='APPEND')
        WRITE(99,'(9F12.6)') mu, total_n, SUM(tot_dm(1:norb)),
     1SUM(tot_dm(ncor_orb-norb+1:ncor_orb)), tot_Ekin, Epot_st+Epot_dy,
     1   SigMdc(2,nom), SigMdc(2*ncor_orb,nom)  
        CLOSE(99)
      endif

      deallocate(Sig_st)
      deallocate(tot_Gloc)
      deallocate(tot_Ed)
      deallocate(tot_dm)
      deallocate(tot_dm_full)
      deallocate(cor_idx)
      deallocate(U)
      deallocate(Uc)
      deallocate(Up)
      deallocate(Jh)
      deallocate(hamr)
      deallocate(kpt_latt)
      deallocate(kweight)
      deallocate(tran)
      deallocate(SigMdc)
      !STOP



      CALL MPI_Finalize(ierr)

      !deallocate(IPIV)
      !deallocate(id)
      !deallocate(eval)
      !deallocate(Dval)
      !deallocate(ham)
      !deallocate(WORK)
      !deallocate(ommesh)
      !deallocate(RWORK)
      !deallocate(dm)
      !deallocate(eval0)
      !deallocate(tot_dm)
      !deallocate(sigi)
      !deallocate(hk)
      !deallocate(Gsum)
      !deallocate(Ginv)
      !deallocate(Gloc)
      !deallocate(scaler)
      !deallocate(evec0)
      !deallocate(evl)
      !deallocate(evr)
      !deallocate(tot_Gloc)
      !deallocate(UM)
      !deallocate(ndimwin)
      !deallocate(lwindow)
      !deallocate(excl_bands)
      !deallocate(lexclude_band)
      !deallocate(UOPTU)
      !deallocate(UOPTU2)
      !deallocate(DMFT_M)
      !deallocate(band_win)
      !deallocate(kpt_array)
      !deallocate(DMFT_w)
      !deallocate(DMFT_U)
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
