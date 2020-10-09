
module dmft_ksum 

  use constants, only: dp
  use io, only: stdout, maxlen
  use read_inputs 

  implicit none

  real(kind=dp), allocatable, save :: DMFT_eval_loc(:,:)
  complex(kind=dp), allocatable, save :: DMFT_evec_loc(:,:,:)
  complex(kind=dp), allocatable, save :: DMFT_EDM_loc(:,:,:)

contains

  subroutine compute_tetra0()
    use constants
    use io
    use utility
    use tetrahedron
    use generate_ham, only: tran
    use comms, only: on_root,comms_barrier,comms_reduce,comms_allreduce,my_node_id, num_nodes, comms_array_split, comms_bcast
   
    implicit none

    logical :: iffile
    character(len=10) :: write_format
    logical :: LOWB,HIGHB
    real(kind=dp) :: low_mu, high_mu, tot_n
    complex(kind=dp) :: Gs, Gs0, new_eval
    integer :: x,y,z,i,j,l,ik,r,ib,im,ispin,loc_i,ierr
    integer :: nbmin,nbmax,num_band_max,ION,IDIR 
    integer :: numk_loc,ntet_loc!,mu_iter 
    ! Needed to split an array on different nodes
    integer, dimension(0:num_nodes - 1) :: counts
    integer, dimension(0:num_nodes - 1) :: displs
    integer, dimension(0:num_nodes - 1) :: counts_tet
    integer, dimension(0:num_nodes - 1) :: displs_tet
    real(kind=dp) :: Ekin, Nd, sweight, dEkin(nions,ndir)
    real(kind=dp), allocatable :: eval0(:,:)
    real(kind=dp) :: rdotk
    real(kind=dp) :: tdos(1)
    real(kind=dp) :: cdos(1)
    real(kind=dp), allocatable :: kweight(:,:)
    real(kind=dp), allocatable :: occ_wann(:)
    complex(kind=dp), allocatable :: evec0(:,:,:)
    complex(kind=dp), allocatable :: Umat(:,:)
    complex(kind=dp), allocatable :: Umat_loc(:,:,:,:)
    complex(kind=dp), allocatable :: dUmat_loc(:,:,:,:)

    call comms_array_split(num_new_kpts, counts, displs)    
    numk_loc=counts(my_node_id)

    if (.not. allocated(eval0)) then
      allocate (eval0(num_new_kpts,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval0 in compute_DMFT_mu')
    endif
    if (.not. allocated(evec0)) then
      allocate (evec0(numk_loc,num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval0 in compute_DMFT_mu')
    endif

    call comms_array_split(num_new_kpts*6, counts_tet, displs_tet)
    ntet_loc=counts_tet(my_node_id)

!    eval0=0.0_dp
!    do ik=1,num_new_kpts
!      nbmin=band_win(1,ik)
!      do j=1,num_wann
!        eval0(ik,j)=eigvals(nbmin+j-1,ik)
!      enddo
!!      print *, eval0(ik,:)
!    enddo
   
    eval0=0.0_dp
    evec0=cmplx_0
    loc_i=1
    do ik=displs(my_node_id)+1,displs(my_node_id)+numk_loc
      nbmin=band_win(1,ik)
      nbmax=band_win(2,ik)
      num_band_max=nbmax-nbmin+1
      do x=1,num_wann
        do y=1,num_wann
          do z=1,num_band_max
            evec0(loc_i,x,y)=evec0(loc_i,x,y)+dconjg(UMatrix(z,x,ik))*UMatrix(z,y,ik)*eigvals(nbmin+z-1,ik)
          enddo
        enddo
      enddo
      CALL EIGENH(evec0(loc_i,:,:),num_wann,eval0(ik,:))
      loc_i=loc_i+1
    enddo
    call comms_allreduce(eval0,num_new_kpts,num_wann,'SUM')

    low_mu=mu_DFT;high_mu=mu_DFT
    LOWB=.FALSE.;HIGHB=.FALSE.
    do l=1,mu_iter
      tot_n=0.0_dp
      mu_DFT=(low_mu+high_mu)/2
      call tetra_totdos(0,num_wann,num_new_kpts,num_new_kpts*6,displs_tet(my_node_id)+1,displs_tet(my_node_id)+ntet_loc,kibz,tetptr,tet_idx,eval0(:,:),mu_DFT,0.0_dp,tdos,cdos)
      tot_n=tot_n+2.0_dp*cdos(1)
      call comms_allreduce(tot_n,1,'SUM')
      !if (on_root) write(*,*) mu, tot_n, n_elec
      if (abs(tot_n-n_elec)<1E-10) EXIT
      if (tot_n<n_elec) then
        if (HIGHB.eqv..FALSE.) then
          high_mu=mu_DFT+(n_elec-tot_n)*2/(0.1*n_elec)
        endif
        low_mu=mu_DFT
        LOWB=.TRUE.
      else
        if (LOWB.eqv..FALSE.) then
          low_mu=mu_DFT+(n_elec-tot_n)*2/(0.1*n_elec)
        endif
        high_mu=mu_DFT
        HIGHB=.TRUE.
      endif
      !if (on_root) write(*,*) mu, tot_n, n_elec
      !STOP
      if (l.eq.mu_iter) then
        write(*,*) 'Fail to find mu, increase mu_iter'
          !STOP
      endif
    enddo

    if (on_root) then   
      OPEN(UNIT=99,FILE='DFT_mu.out',FORM='FORMATTED')
      WRITE(99,'(F20.15)') mu_DFT
      CLOSE(99)
    endif

    if (.not. allocated(kweight)) then
      allocate (kweight(num_new_kpts,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval0 in compute_DMFT_mu')
    endif

    call tetra_poccup(num_wann,num_new_kpts,num_new_kpts*6,displs_tet(my_node_id)+1,displs_tet(my_node_id)+ntet_loc,kibz,tetptr,tet_idx,eval0(:,:),mu_DFT,kweight(:,:))
    call comms_allreduce(kweight,num_new_kpts,num_wann,'SUM')
    if (allocated(eval0)) deallocate (eval0)

    if (.not. allocated(occ_wann)) then
      allocate (occ_wann(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Gloc_sum in compute_G_loc')
    endif

    sweight=2.0_dp/dfloat(num_new_kpts)
    occ_wann=0.0_dp
    tot_n=0.0_dp;Nd=0.0_dp;Ekin=0.0_dp
    dEkin=0.0_dp
    !print *, lforce
    loc_i=1
    do ik=displs(my_node_id)+1,displs(my_node_id)+numk_loc
      do ib=1,num_wann
        do i=1,num_wann
          occ_wann(i)=occ_wann(i)+real(evec0(loc_i,i,ib)*dconjg(evec0(loc_i,i,ib)))*kweight(ik,ib)*sweight
        enddo
      enddo
      nbmin=band_win(1,ik)
      nbmax=band_win(2,ik)
      num_band_max=nbmax-nbmin+1
      if (.not. allocated(Umat)) then
        allocate (Umat(num_band_max,num_wann), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating Gloc_sum in compute_G_loc')
      endif
      if (lforce .eq. .true.) then 
        if (.not. allocated(Umat_loc)) then
          allocate (Umat_loc(num_band_max,num_wann,ndir,nions), stat=ierr)
          if (ierr /= 0) call io_error('Error allocating Gloc_sum in compute_G_loc')
        endif
        if (.not. allocated(dUmat_loc)) then
          allocate (dUmat_loc(num_band_max,num_wann,ndir,nions), stat=ierr)
          if (ierr /= 0) call io_error('Error allocating Gloc_sum in compute_G_loc')
        endif
        Umat_loc=cmplx_0
        dUmat_loc=cmplx_0
        do ION=1,nions
          do IDIR=1,ndir   
            Umat_loc(:,:,IDIR,ION)=matmul(UMatrix_loc(:,:,ik,IDIR,ION),evec0(loc_i,:,:))
            dUmat_loc(:,:,IDIR,ION)=matmul(dUMatrix(:,:,ik,IDIR,ION),evec0(loc_i,:,:))
          enddo
        enddo
      endif
      Umat=cmplx_0
      Umat=matmul(UMatrix(:,:,ik),evec0(loc_i,:,:))
      do ib=1,num_band_max
        do i=1,num_wann
        !do i=1,5
          Ekin=Ekin+eigvals(nbmin+ib-1,ik)*kweight(ik,i)*sweight*real(Umat(ib,i)*dconjg(Umat(ib,i)))
          if (lforce.eq..true.) then
            do ION=1,nions
              do IDIR=1,ndir 
                !dEkin(ION,IDIR)=dEkin(ION,IDIR)+eigvals(nbmin+ib-1,ik)*kweight(ik,i)*sweight*real(Umat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION)))
                dEkin(ION,IDIR)=dEkin(ION,IDIR)+deig(nbmin+ib-1,ik,IDIR,ION)*kweight(ik,i)*sweight*real(Umat(ib,i)*dconjg(Umat(ib,i)))
                !dEkin(ION,IDIR)=dEkin(ION,IDIR)+eigvals(nbmin+ib-1,ik)*kweight(ik,i)*sweight*real(dUmat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION))+Umat_loc(ib,i,IDIR,ION)*dconjg(dUmat_loc(ib,i,IDIR,ION)))
              enddo
            enddo
          endif
          !if (lforce.eq..true.) dEkin=dEkin+eigvals(nbmin+ib-1,ik)*kweight(ik,i)*sweight*real(dUmat_loc(ib,i)*dconjg(Umat_loc(ib,i))+Umat_loc(ib,i)*dconjg(dUmat_loc(ib,i)))
        enddo
      enddo
      if (allocated(Umat)) deallocate (Umat,stat=ierr)
      if (allocated(Umat_loc)) deallocate (Umat_loc,stat=ierr)
      if (allocated(dUmat_loc)) deallocate (dUmat_loc,stat=ierr)
      loc_i=loc_i+1
    enddo
    do ib=1,num_wann
      tot_n=tot_n+occ_wann(ib)
    enddo
    do ib=1,5
      Nd=Nd+occ_wann(ib)
    enddo
    !print *, dEkin
    call comms_allreduce(occ_wann(1),num_wann,'SUM')
    call comms_allreduce(Ekin,1,'SUM')
    call comms_allreduce(dEkin,nions,ndir,'SUM')
    call comms_allreduce(tot_n,1,'SUM')
    call comms_allreduce(Nd,1,'SUM')
     
    if (on_root) then   
!
      OPEN(UNIT=99,FILE='INFO_KSUM',STATUS='old',FORM='FORMATTED',ACCESS='APPEND')
      WRITE(99,'(11F12.6)') mu_DFT, tot_n, SUM(occ_wann(1:n_orbs)), &
      SUM(occ_wann(n_orbs*(n_atoms-1)+1:n_orbs*n_atoms)), Ekin, 0.0, dEkin(1,1) 
      CLOSE(99)

      OPEN(UNIT=99,FILE='INFO_DM',STATUS='old',FORM='FORMATTED',ACCESS='APPEND')
90    FORMAT('(',I2,'F20.15)')
      WRITE(write_format,90) n_atoms*n_orbs
      WRITE(99,write_format,advance='no') (occ_wann(ib), ib=1,n_atoms*n_orbs)
      WRITE(99,write_format) (0.0, ib=1,n_atoms*n_orbs)
      CLOSE(99)
      if (lforce.eq..true.) then
        OPEN(UNIT=99,FILE='FORCE0',FORM='FORMATTED')
        do ION=1,nions
          WRITE(99,'(3F12.8)') dEkin(ION,1), dEkin(ION,2), dEkin(ION,3)
        enddo
        CLOSE(99)
      endif
    endif
    if (allocated(occ_wann)) deallocate (occ_wann)

  end subroutine compute_tetra0

  subroutine compute_tetra0_ibz(n_kpts_loc,kpt_dft,wght_dft)
    use constants
    use io
    use utility
    use tetrahedron
    use generate_ham, only: tran
    use comms, only: on_root,comms_barrier,comms_reduce,comms_allreduce,my_node_id, num_nodes, comms_array_split, comms_bcast
   
    implicit none

    integer, intent(in) :: n_kpts_loc
    real(kind=dp), intent(in) :: kpt_dft(3,n_kpts_loc)
    real(kind=dp), intent(in) :: wght_dft(n_kpts_loc)
    logical :: iffile
    logical :: LOWB,HIGHB
    real(kind=dp) :: low_mu, high_mu, tot_n, sweight
    complex(kind=dp) :: Gs, Gs0, new_eval, DMFT_UU
    integer :: x,y,z,i,j,l,ik,r,ib,im,ispin,loc_i,ierr
    integer :: nbmin,nbmax,num_band_max,kx_loc,ky_loc,kz_loc 
    integer :: numk_loc,ntet_loc,nibz_loc!,mu_iter 
    ! Needed to split an array on different nodes
    integer, dimension(0:num_nodes - 1) :: counts
    integer, dimension(0:num_nodes - 1) :: displs
    integer, dimension(0:num_nodes - 1) :: counts_tet
    integer, dimension(0:num_nodes - 1) :: displs_tet
    integer, dimension(0:num_nodes - 1) :: counts_ibz
    integer, dimension(0:num_nodes - 1) :: displs_ibz
    real(kind=dp), allocatable :: kweight(:,:)
    real(kind=dp), allocatable :: occ_wann(:)
    real(kind=dp), allocatable :: eval0(:,:)
    real(kind=dp) :: rdotk
    real(kind=dp) :: tdos(1)
    real(kind=dp) :: cdos(1)
    complex(kind=dp), allocatable :: Hk(:,:)
    complex(kind=dp), allocatable :: evec0(:,:,:)

    call comms_array_split(num_new_kpts, counts, displs)    
    numk_loc=counts(my_node_id)
    call comms_array_split(num_new_kpts*6, counts_tet, displs_tet)
    ntet_loc=counts_tet(my_node_id)
    call comms_array_split(n_kpts_loc, counts_ibz, displs_ibz)    
    nibz_loc=counts_ibz(my_node_id)
    inquire(file='DFT_mu.out',exist=iffile)
    if (iffile.eqv. .false.) then
      write(*,*) 'input is needed for new mu!'
      STOP
    else
      open(65,file="DFT_mu.out")
      read(65,*) mu
      close(65)
    endif

    if (.not. allocated(eval0)) then
      allocate (eval0(num_new_kpts,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval0 in compute_DMFT_mu')
    endif
    if (.not. allocated(Hk)) then
      allocate (Hk(num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Hk in compute_DMFT_mu')
    endif
    if (.not. allocated(evec0)) then
      allocate (evec0(num_new_kpts,num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Hk in compute_DMFT_mu')
    endif
    if (.not. allocated(kweight)) then
      allocate (kweight(num_new_kpts,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Hk in compute_DMFT_mu')
    endif
    if (.not. allocated(DMFT_evec_loc)) then
      allocate (DMFT_evec_loc(num_bands,num_wann,n_kpts_loc), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating DMFT_evec_loc in compute_G_loc')
    endif


    eval0=0.0_dp
    evec0=cmplx_0
    DMFT_evec_loc=cmplx_0
    do ik=displs(my_node_id)+1,displs(my_node_id)+numk_loc
      nbmin=band_win(1,ik)
      nbmax=band_win(2,ik)
      num_band_max=nbmax-nbmin+1
      do x=1,num_wann
        do y=1,num_wann
          do z=1,num_band_max
            evec0(ik,x,y)=evec0(ik,x,y)+dconjg(UMatrix(z,x,ik))*UMatrix(z,y,ik)*eigvals(nbmin+z-1,ik)
          enddo
        enddo
      enddo
      CALL EIGENH(evec0(ik,:,:),num_wann,eval0(ik,:))
    enddo
    call comms_allreduce(eval0,num_new_kpts,num_wann,'SUM')
    call comms_allreduce(evec0,num_new_kpts,num_wann,num_wann,'SUM')

    low_mu=mu;high_mu=mu
    LOWB=.FALSE.;HIGHB=.FALSE.
    do l=1,mu_iter
      tot_n=0.0_dp
      mu=(low_mu+high_mu)/2
      kweight=0.0_dp
      call tetra_poccup(num_wann,num_new_kpts,num_new_kpts*6,displs_tet(my_node_id)+1,displs_tet(my_node_id)+ntet_loc,kibz,tetptr,tet_idx,eval0(:,:),mu,kweight(:,:))
      call comms_allreduce(kweight,num_new_kpts,num_wann,'SUM')
      do ik=displs_ibz(my_node_id)+1,displs_ibz(my_node_id)+nibz_loc
        kx_loc=MODULO(NINT(kpt_dft(1,ik)*mp_grid(1)-0.25),mp_grid(1))
        ky_loc=MODULO(NINT(kpt_dft(2,ik)*mp_grid(2)-0.25),mp_grid(2))
        kz_loc=MODULO(NINT(kpt_dft(3,ik)*mp_grid(3)-0.25),mp_grid(3))
        do i=1,num_new_kpts
          if (kx_loc.eq.ikpt(1,i) .and. ky_loc.eq.ikpt(2,i) .and. kz_loc.eq.ikpt(3,i)) then
            loc_i=i
            exit
          endif
        enddo
        do ib=1,num_wann
          tot_n=tot_n+2.0_dp*kweight(loc_i,ib)*wght_dft(ik)
        enddo
      enddo
      call comms_allreduce(tot_n,1,'SUM')
!      if (on_root) write(*,*) mu, tot_n, n_elec
      if (abs(tot_n-n_elec)<1E-10) EXIT
      if (tot_n<n_elec) then
        if (HIGHB.eqv..FALSE.) then
          high_mu=mu+(n_elec-tot_n)*2/(0.1*n_elec)
        endif
        low_mu=mu
        LOWB=.TRUE.
      else
        if (LOWB.eqv..FALSE.) then
          low_mu=mu+(n_elec-tot_n)*2/(0.1*n_elec)
        endif
        high_mu=mu
        HIGHB=.TRUE.
      endif
      !if (on_root) write(*,*) mu, tot_n, n_elec
      !STOP
      if (l.eq.mu_iter) then
        write(*,*) 'Fail to find mu, increase mu_iter'
          !STOP
      endif
    enddo

    DMFT_evec_loc=cmplx_0
    do ik=displs_ibz(my_node_id)+1,displs_ibz(my_node_id)+nibz_loc
      kx_loc=MODULO(NINT(kpt_dft(1,ik)*mp_grid(1)-0.25),mp_grid(1))
      ky_loc=MODULO(NINT(kpt_dft(2,ik)*mp_grid(2)-0.25),mp_grid(2))
      kz_loc=MODULO(NINT(kpt_dft(3,ik)*mp_grid(3)-0.25),mp_grid(3))
      do i=1,num_new_kpts
        if (kx_loc.eq.ikpt(1,i) .and. ky_loc.eq.ikpt(2,i) .and. kz_loc.eq.ikpt(3,i)) then
          loc_i=i
          exit
        endif
      enddo 
      do ib=1,num_wann
        do i=1,num_wann
          do j=1,num_wann
            DMFT_UU=evec0(loc_i,i,ib)*dconjg(evec0(loc_i,j,ib))*kweight(loc_i,ib)
            DMFT_evec_loc(i,j,ik)=DMFT_evec_loc(i,j,ik)-DMFT_UU
          enddo
        enddo
      enddo
    enddo

    if (allocated(evec0)) deallocate (evec0)
    if (allocated(eval0)) deallocate (eval0)
    if (allocated(kweight)) deallocate (kweight)

    if (on_root) then   
      OPEN(UNIT=99,FILE='DFT_mu.out',FORM='FORMATTED')
      WRITE(99,'(F20.15)') mu
      CLOSE(99)
    endif

  end subroutine compute_tetra0_ibz

  subroutine compute_tetra_ibz_mu(n_kpts_loc,kpt_dft,wght_dft)
    use constants
    use io
    use utility
    use tetrahedron
    use generate_ham, only: tran
    use comms, only: on_root,comms_barrier,comms_reduce,comms_allreduce,my_node_id, num_nodes, comms_array_split, comms_bcast
   
    implicit none

    integer, intent(in) :: n_kpts_loc
    real(kind=dp), intent(in) :: kpt_dft(3,n_kpts_loc)
    real(kind=dp), intent(in) :: wght_dft(n_kpts_loc)
    logical :: iffile
    logical :: LOWB,HIGHB
    real(kind=dp) :: low_mu, high_mu, tot_n, sweight
    complex(kind=dp) :: Gs, Gs0, new_eval
    integer :: x,y,z,i,j,l,ik,r,ib,im,ispin,loc_i,ierr
    integer :: nbmin,nbmax,num_band_max,kx_loc,ky_loc,kz_loc 
    integer :: numk_loc,ntet_loc,nibz_loc!,mu_iter 
    ! Needed to split an array on different nodes
    integer, dimension(0:num_nodes - 1) :: counts
    integer, dimension(0:num_nodes - 1) :: displs
    integer, dimension(0:num_nodes - 1) :: counts_tet
    integer, dimension(0:num_nodes - 1) :: displs_tet
    integer, dimension(0:num_nodes - 1) :: counts_ibz
    integer, dimension(0:num_nodes - 1) :: displs_ibz
    real(kind=dp), allocatable :: kweight(:,:)
    real(kind=dp), allocatable :: occ_wann(:)
    real(kind=dp), allocatable :: eval0(:,:,:)
    real(kind=dp) :: rdotk
    real(kind=dp) :: tdos(1)
    real(kind=dp) :: cdos(1)
    complex(kind=dp), allocatable :: Hk(:,:)
    complex(kind=dp), allocatable :: evec0(:,:)

    call comms_array_split(num_new_kpts, counts, displs)    
    numk_loc=counts(my_node_id)
    call comms_array_split(num_new_kpts*6, counts_tet, displs_tet)
    ntet_loc=counts_tet(my_node_id)
    call comms_array_split(n_kpts_loc, counts_ibz, displs_ibz)    
    nibz_loc=counts_ibz(my_node_id)

    if (.not. allocated(eval0)) then
      allocate (eval0(nspin,num_new_kpts,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval0 in compute_DMFT_mu')
    endif
    if (.not. allocated(Hk)) then
      allocate (Hk(num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Hk in compute_DMFT_mu')
    endif
    if (.not. allocated(evec0)) then
      allocate (evec0(num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Hk in compute_DMFT_mu')
    endif
    if (.not. allocated(kweight)) then
      allocate (kweight(num_new_kpts,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Hk in compute_DMFT_mu')
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


!    eval0=0.0_dp
!    do ik=1,num_new_kpts
!      nbmin=band_win(1,ik)
!      do j=1,num_wann
!        eval0(1,ik,j)=eigvals(nbmin+j-1,ik)
!      enddo
!    enddo
!!    print *, (eval0(1,1,:))
    eval0=0.0_dp
    do ik=displs(my_node_id)+1,displs(my_node_id)+numk_loc
      Hk=cmplx_0
      nbmin=band_win(1,ik)
      nbmax=band_win(2,ik)
      num_band_max=nbmax-nbmin+1
      do x=1,num_wann
        do y=1,num_wann
          do z=1,num_band_max
            Hk(x,y)=Hk(x,y)+dconjg(UMatrix(z,x,ik))*UMatrix(z,y,ik)*eigvals(nbmin+z-1,ik)
          enddo
        enddo
      enddo
      do ispin=1,nspin
        evec0=cmplx_0
        call zcopy(num_wann**2,Hk(:,:),1,evec0(:,:),1)
        do i=1,n_atoms
          do j=1,n_orbs
            if (sym_idx(ispin,i,j)>0) evec0((i-1)*n_orbs+j,(i-1)*n_orbs+j)=evec0((i-1)*n_orbs+j,(i-1)*n_orbs+j)+Sigoo(sym_idx(ispin,i,j))
          enddo
        enddo
        CALL EIGENVALH(evec0(:,:),num_wann,eval0(ispin,ik,:))
      enddo
    enddo
    if (allocated(Hk)) deallocate (Hk)
    if (allocated(evec0)) deallocate (evec0)
    call comms_allreduce(eval0,nspin,num_new_kpts,num_wann,'SUM')

    low_mu=mu;high_mu=mu
    LOWB=.FALSE.;HIGHB=.FALSE.
    do l=1,mu_iter
      tot_n=0.0_dp
      mu=(low_mu+high_mu)/2
      do ispin=1,nspin
        kweight=0.0_dp
        call tetra_poccup(num_wann,num_new_kpts,num_new_kpts*6,displs_tet(my_node_id)+1,displs_tet(my_node_id)+ntet_loc,kibz,tetptr,tet_idx,eval0(ispin,:,:),mu,kweight(:,:))
        call comms_allreduce(kweight,num_new_kpts,num_wann,'SUM')
        do ik=displs_ibz(my_node_id)+1,displs_ibz(my_node_id)+nibz_loc
          kx_loc=MODULO(NINT(kpt_dft(1,ik)*mp_grid(1)-0.25),mp_grid(1))
          ky_loc=MODULO(NINT(kpt_dft(2,ik)*mp_grid(2)-0.25),mp_grid(2))
          kz_loc=MODULO(NINT(kpt_dft(3,ik)*mp_grid(3)-0.25),mp_grid(3))
          do i=1,num_new_kpts
            if (kx_loc.eq.ikpt(1,i) .and. ky_loc.eq.ikpt(2,i) .and. kz_loc.eq.ikpt(3,i)) then
              loc_i=i
              exit
            endif
          enddo
          do ib=1,num_wann
            tot_n=tot_n+2.0_dp*kweight(loc_i,ib)*wght_dft(ik)/dfloat(nspin)
          enddo
        enddo
      enddo
      call comms_allreduce(tot_n,1,'SUM')
!      if (on_root) write(*,*) mu, tot_n, n_elec
      if (abs(tot_n-n_elec)<1E-10) EXIT
      if (tot_n<n_elec) then
        if (HIGHB.eqv..FALSE.) then
          high_mu=mu+(n_elec-tot_n)*2/(0.1*n_elec)
        endif
        low_mu=mu
        LOWB=.TRUE.
      else
        if (LOWB.eqv..FALSE.) then
          low_mu=mu+(n_elec-tot_n)*2/(0.1*n_elec)
        endif
        high_mu=mu
        HIGHB=.TRUE.
      endif
      !if (on_root) write(*,*) mu, tot_n, n_elec
      !STOP
      if (l.eq.mu_iter) then
        write(*,*) 'Fail to find mu, increase mu_iter'
          !STOP
      endif
    enddo

    if (allocated(eval0)) deallocate (eval0)
    if (allocated(kweight)) deallocate (kweight)

    if (on_root) then   
      OPEN(UNIT=99,FILE='DMFT_mu.out',FORM='FORMATTED')
      WRITE(99,'(F20.15)') mu
      CLOSE(99)
    endif

  end subroutine compute_tetra_ibz_mu

  subroutine compute_tetra_mu()
    use constants
    use io
    use utility
    use tetrahedron
    use generate_ham, only: tran
    use comms, only: on_root,comms_barrier,comms_reduce,comms_allreduce,my_node_id, num_nodes, comms_array_split, comms_bcast
   
    implicit none

    logical :: iffile
    logical :: LOWB,HIGHB
    real(kind=dp) :: low_mu, high_mu, tot_n, sweight
    complex(kind=dp) :: Gs, Gs0, new_eval
    integer :: x,y,z,i,j,l,ik,r,ib,im,ispin,loc_i,ierr
    integer :: nbmin,nbmax,num_band_max 
    integer :: numk_loc,ntet_loc!,mu_iter 
    ! Needed to split an array on different nodes
    integer, dimension(0:num_nodes - 1) :: counts
    integer, dimension(0:num_nodes - 1) :: displs
    integer, dimension(0:num_nodes - 1) :: counts_tet
    integer, dimension(0:num_nodes - 1) :: displs_tet
    real(kind=dp), allocatable :: kweight(:,:)
    real(kind=dp), allocatable :: occ_wann(:)
    real(kind=dp), allocatable :: eval0(:,:,:)
    real(kind=dp) :: rdotk
    real(kind=dp) :: tdos(1)
    real(kind=dp) :: cdos(1)
    complex(kind=dp), allocatable :: Hk(:,:)
    complex(kind=dp), allocatable :: evec0(:,:)

    call comms_array_split(num_new_kpts, counts, displs)    
    numk_loc=counts(my_node_id)
    inquire(file='DMFT_mu.out',exist=iffile)
    if (iffile.eqv. .false.) then
      write(*,*) 'input is needed for new mu!'
      STOP
    else
      open(65,file="DMFT_mu.out")
      read(65,*) mu
      close(65)
    endif

    if (.not. allocated(eval0)) then
      allocate (eval0(nspin,num_new_kpts,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval0 in compute_DMFT_mu')
    endif
    if (.not. allocated(Hk)) then
      allocate (Hk(num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Hk in compute_DMFT_mu')
    endif
    if (.not. allocated(evec0)) then
      allocate (evec0(num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Hk in compute_DMFT_mu')
    endif

    call comms_array_split(num_new_kpts*6, counts_tet, displs_tet)
    ntet_loc=counts_tet(my_node_id)

!    eval0=0.0_dp
!    do ik=1,num_new_kpts
!      nbmin=band_win(1,ik)
!      do j=1,num_wann
!        eval0(1,ik,j)=eigvals(nbmin+j-1,ik)
!      enddo
!    enddo
!!    print *, (eval0(1,1,:))
    eval0=0.0_dp
    do ik=displs(my_node_id)+1,displs(my_node_id)+numk_loc
      Hk=cmplx_0
      nbmin=band_win(1,ik)
      nbmax=band_win(2,ik)
      num_band_max=nbmax-nbmin+1
      do x=1,num_wann
        do y=1,num_wann
          do z=1,num_band_max
            Hk(x,y)=Hk(x,y)+dconjg(UMatrix(z,x,ik))*UMatrix(z,y,ik)*eigvals(nbmin+z-1,ik)
          enddo
        enddo
      enddo
      do ispin=1,nspin
        evec0=cmplx_0
        call zcopy(num_wann**2,Hk(:,:),1,evec0(:,:),1)
        do i=1,n_atoms
          do j=1,n_orbs
            if (sym_idx(ispin,i,j)>0) evec0((i-1)*n_orbs+j,(i-1)*n_orbs+j)=evec0((i-1)*n_orbs+j,(i-1)*n_orbs+j)+Sigoo(sym_idx(ispin,i,j))
          enddo
        enddo
        CALL EIGENVALH(evec0(:,:),num_wann,eval0(ispin,ik,:))
      enddo
    enddo
    if (allocated(Hk)) deallocate (Hk)
    call comms_allreduce(eval0,nspin,num_new_kpts,num_wann,'SUM')

    low_mu=mu;high_mu=mu
    LOWB=.FALSE.;HIGHB=.FALSE.
    sweight=2.0_dp/nspin
    do l=1,mu_iter
      tot_n=0.0_dp
      mu=(low_mu+high_mu)/2
      do ispin=1,nspin
        call tetra_totdos(0,num_wann,num_new_kpts,num_new_kpts*6,displs_tet(my_node_id)+1,displs_tet(my_node_id)+ntet_loc,kibz,tetptr,tet_idx,eval0(ispin,:,:),mu,0.0_dp,tdos,cdos)
        tot_n=tot_n+sweight*cdos(1)
      enddo
      call comms_allreduce(tot_n,1,'SUM')
      if (on_root) write(*,*) mu, tot_n, n_elec
      if (abs(tot_n-n_elec)<1E-10) EXIT
      if (tot_n<n_elec) then
        if (HIGHB.eqv..FALSE.) then
          high_mu=mu+(n_elec-tot_n)*2/(0.1*n_elec)
        endif
        low_mu=mu
        LOWB=.TRUE.
      else
        if (LOWB.eqv..FALSE.) then
          low_mu=mu+(n_elec-tot_n)*2/(0.1*n_elec)
        endif
        high_mu=mu
        HIGHB=.TRUE.
      endif
      !if (on_root) write(*,*) mu, tot_n, n_elec
      !STOP
      if (l.eq.mu_iter) then
        write(*,*) 'Fail to find mu, increase mu_iter'
          !STOP
      endif
    enddo

    if (allocated(eval0)) deallocate (eval0)

    if (on_root) then   
      OPEN(UNIT=99,FILE='DMFT_mu.out',FORM='FORMATTED')
      WRITE(99,'(F20.15)') mu
      CLOSE(99)
    endif

  end subroutine compute_tetra_mu

  subroutine compute_tetra()
    use constants
    use io
    use utility
    use tetrahedron
    use generate_ham, only: tran
    use comms, only: on_root,comms_allreduce,my_node_id, num_nodes, comms_array_split,comms_reduce, comms_barrier
   
    implicit none

    character(len=10) :: write_format
    real(kind=dp) :: tot_n
    complex(kind=dp) :: Gs, Gs0
    integer :: i,j,l,ik,r,ib,im,ispin,loc_i,ierr,x,y,z
    integer :: nbmin,nbmax,num_band_max, iter,Niter 
    integer :: numk_loc,nfine_tot, ntet_loc, nbin
    !integer,intent(in) :: n_kpts_loc,n_wann
    ! Needed to split an array on different nodes
    integer, dimension(0:num_nodes - 1) :: counts
    integer, dimension(0:num_nodes - 1) :: displs
    integer, dimension(0:num_nodes - 1) :: counts_tet
    integer, dimension(0:num_nodes - 1) :: displs_tet
    real(kind=dp) :: rdotk, fermi0, tot, IDIR, ION
    real(kind=dp) :: eps,deps, dosef, Ekin, Nd, sweight,dm_loc, dEkin(nions,ndir), dN(nspin,n_atoms*n_orbs,nions,ndir), dmu(nions,ndir)
    real(kind=dp) :: dSig(nspin,n_atoms,nions,ndir), dN2(num_bands,nions,ndir)
    real(kind=dp), allocatable :: eval0(:,:,:)
    real(kind=dp), allocatable :: Ed(:)
    real(kind=dp), allocatable :: kweight(:,:,:)
    real(kind=dp), allocatable :: occ_wann(:)
    real(kind=dp), allocatable :: mom(:)
    real(kind=dp), allocatable :: cpdos(:,:,:,:)
    real(kind=dp), allocatable :: pdos(:,:,:,:)
    complex(kind=dp), allocatable :: Umat(:,:)
    complex(kind=dp), allocatable :: Umat_loc(:,:,:,:)
    complex(kind=dp), allocatable :: dUmat_loc(:,:,:,:)
    complex(kind=dp), allocatable :: Hk(:,:)
    complex(kind=dp), allocatable :: evec0(:,:,:,:)
!    complex(kind=dp), allocatable :: (:,:,:)
    complex(kind=dp), allocatable :: Delta(:,:)
    complex(kind=dp), allocatable :: Gloc(:,:,:)
    complex(kind=dp), allocatable :: sym_Gloc(:,:)
    complex(kind=dp) :: DMFT_UU,DMFT_U0 

    call comms_array_split(num_new_kpts, counts, displs)
    numk_loc=counts(my_node_id)
    call comms_array_split(num_new_kpts*6, counts_tet, displs_tet)
    ntet_loc=counts_tet(my_node_id)

    if (.not. allocated(eval0)) then
      allocate (eval0(nspin,num_new_kpts,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval0 in compute_DMFT_mu')
    endif
    if (.not. allocated(evec0)) then
      allocate (evec0(nspin,numk_loc,num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval0 in compute_DMFT_mu')
    endif
    if (.not. allocated(Hk)) then
      allocate (Hk(num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Hk in compute_G_loc')
    endif
    if (.not. allocated(kweight)) then
      allocate (kweight(nspin,num_new_kpts,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval0 in compute_DMFT_mu')
    endif


    eval0=0.0_dp
    evec0=cmplx_0
    loc_i=1
    do ik=displs(my_node_id)+1,displs(my_node_id)+numk_loc
      Hk=cmplx_0
      nbmin=band_win(1,ik)
      nbmax=band_win(2,ik)
      num_band_max=nbmax-nbmin+1
      do x=1,num_wann
        do y=1,num_wann
          do z=1,num_band_max
            Hk(x,y)=Hk(x,y)+dconjg(UMatrix(z,x,ik))*UMatrix(z,y,ik)*eigvals(nbmin+z-1,ik)
          enddo
        enddo
      enddo
      do ispin=1,nspin
        call zcopy(num_wann**2,Hk(:,:),1,evec0(ispin,loc_i,:,:),1)
        do i=1,n_atoms
          do j=1,n_orbs
            if (sym_idx(ispin,i,j)>0) evec0(ispin,loc_i,(i-1)*n_orbs+j,(i-1)*n_orbs+j)=evec0(ispin,loc_i,(i-1)*n_orbs+j,(i-1)*n_orbs+j)+Sigoo(sym_idx(ispin,i,j))
          enddo
        enddo
        CALL EIGENH(evec0(ispin,loc_i,:,:),num_wann,eval0(ispin,ik,:))
      enddo
      loc_i=loc_i+1
    enddo
    if (allocated(Hk)) deallocate (Hk)
    call comms_allreduce(eval0,nspin,num_new_kpts,num_wann,'SUM')

    nbin=0
    if (.not. allocated(pdos)) then
      allocate (pdos(nspin,nbin+1,num_new_kpts,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval0 in compute_DMFT_mu')
    endif
    if (.not. allocated(cpdos)) then
      allocate (cpdos(nspin,nbin+1,num_new_kpts,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval0 in compute_DMFT_mu')
    endif
    cpdos=0.0_dp
    pdos=0.0_dp
    
    do ispin=1,nspin
      call tetra_poccup(num_wann,num_new_kpts,num_new_kpts*6,displs_tet(my_node_id)+1,displs_tet(my_node_id)+ntet_loc,kibz,tetptr,tet_idx,eval0(ispin,:,:),mu,kweight(ispin,:,:))
      call tetra_pdos(nbin,num_wann,num_new_kpts,num_new_kpts*6,displs_tet(my_node_id)+1,displs_tet(my_node_id)+ntet_loc,kibz,tetptr,tet_idx,eval0(ispin,:,:),mu,mu,pdos(ispin,:,:,:),cpdos(ispin,:,:,:))
      !do ik=1,num_new_kpts
      !enddo
    enddo
    call comms_allreduce(kweight,nspin,num_new_kpts,num_wann,'SUM')
    call comms_allreduce(pdos,nspin,nbin+1,num_new_kpts,num_wann,'SUM')
    call comms_allreduce(cpdos,nspin,nbin+1,num_new_kpts,num_wann,'SUM')
!    print *, sum(kweight(1,1,:))

    dosef=0.0_dp
    sweight=2.0_dp/dfloat(nspin)/dfloat(num_new_kpts)
    do ik=displs(my_node_id)+1,displs(my_node_id)+numk_loc
      do ispin=1,nspin
        do i=1,num_wann
          dosef=dosef+pdos(ispin,1,ik,i)*sweight
        enddo
      enddo
    enddo
    call comms_allreduce(dosef,1,'SUM')

    if (.not. allocated(occ_wann)) then
      allocate (occ_wann(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Gloc_sum in compute_G_loc')
    endif
    if (.not. allocated(mom)) then
      allocate (mom(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Gloc_sum in compute_G_loc')
    endif


    occ_wann=0.0_dp
    tot_n=0.0_dp;Nd=0.0_dp;Ekin=0.0_dp;mom=0.0_dp
    dEkin=0.0_dp
    dmu=0.0_dp
    !Gloc_sum=cmplx_0
    loc_i=1
    do ik=displs(my_node_id)+1,displs(my_node_id)+numk_loc
      do ispin=1,nspin
        do ib=1,num_wann
          do i=1,num_wann
            dm_loc=real(evec0(ispin,loc_i,i,ib)*dconjg(evec0(ispin,loc_i,i,ib)))*kweight(ispin,ik,ib)*sweight
            occ_wann(i)=occ_wann(i)+dm_loc
            if (ispin.eq.1 .AND. nspin.eq.2) then
              mom(i)=mom(i)+dm_loc
            endif
            if (ispin.eq.2) mom(i)=mom(i)-dm_loc
          enddo
        enddo
      enddo
      nbmin=band_win(1,ik)
      nbmax=band_win(2,ik)
      num_band_max=nbmax-nbmin+1
      do ispin=1,nspin
        if (.not. allocated(Umat)) then
          allocate (Umat(num_band_max,num_wann), stat=ierr)
          if (ierr /= 0) call io_error('Error allocating Gloc_sum in compute_G_loc')
        endif
        Umat=cmplx_0
        Umat=matmul(UMatrix(:,:,ik),evec0(ispin,loc_i,:,:))
        if (lforce .eq. .true.) then 
          if (.not. allocated(Umat_loc)) then
            allocate (Umat_loc(num_band_max,num_wann,ndir,nions), stat=ierr)
            if (ierr /= 0) call io_error('Error allocating Gloc_sum in compute_G_loc')
          endif
          if (.not. allocated(dUmat_loc)) then
            allocate (dUmat_loc(num_band_max,num_wann,ndir,nions), stat=ierr)
            if (ierr /= 0) call io_error('Error allocating Gloc_sum in compute_G_loc')
          endif
          Umat_loc=cmplx_0
          dUmat_loc=cmplx_0
          do ION=1,nions
            do IDIR=1,ndir   
              Umat_loc(:,:,IDIR,ION)=matmul(UMatrix_loc(:,:,ik,IDIR,ION),evec0(ispin,loc_i,:,:))
              dUmat_loc(:,:,IDIR,ION)=matmul(dUMatrix(:,:,ik,IDIR,ION),evec0(ispin,loc_i,:,:))
              !dUmat_loc(:,:,IDIR,ION)=matmul(dconjg(transpose(dUMatrix(:,:,ik,IDIR,ION))),Umat_loc(:,:,IDIR,ION))
            enddo
          enddo
        endif
        do ib=1,num_band_max
          do i=1,num_wann
            Ekin=Ekin+eigvals(nbmin+ib-1,ik)*kweight(ispin,ik,i)*sweight*real(Umat(ib,i)*dconjg(Umat(ib,i)))
          enddo
        enddo
        if (lforce.eq..true.) then
          do ION=1,nions
            do IDIR=1,ndir 
              do i=1,num_wann
                deps=0.0_dp
                eps=0.0_dp
                do ib=1,num_band_max
                  !deps=deps+eigvals(nbmin+ib-1,ik)*real(dUmat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION))+Umat_loc(ib,i,IDIR,ION)*dconjg(dUmat_loc(ib,i,IDIR,ION)))+deig(nbmin+ib-1,ik,IDIR,ION)*real(Umat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION)))
                  deps=deps+deig(nbmin+ib-1,ik,IDIR,ION)*real(Umat(ib,i)*dconjg(Umat(ib,i)))
                  eps=eps+eigvals(nbmin+ib-1,ik)*real(Umat(ib,i)*dconjg(Umat(ib,i)))
                enddo
                dmu(ION,IDIR)=dmu(ION,IDIR)+deps*pdos(ispin,1,ik,i)*sweight/dosef
                !dmu(ION,IDIR)=dmu(ION,IDIR)+eigvals(nbmin+ib-1,ik)*pdos(ispin,1,ik,i)*sweight*real(dUmat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION))+Umat_loc(ib,i,IDIR,ION)*dconjg(dUmat_loc(ib,i,IDIR,ION)))/dosef
              !dmu(ION,IDIR)=dmu(ION,IDIR)+pdos(ispin,1,ik,i)*sweight*real(Umat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION)))
              !dEkin(ION,IDIR)=dEkin(ION,IDIR)+pdos(ispin,1,ik,i)*sweight*real(Umat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION)))
              !dmu(ION,IDIR)=dmu(ION,IDIR)+deig(nbmin+ib-1,ik,IDIR,ION)*pdos(ispin,1,ik,i)*sweight*real(Umat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION)))!/dosef
              !dmu(ION,IDIR)=dmu(ION,IDIR)+eigvals(nbmin+ib-1,ik)*pdos(ispin,1,ik,i)*sweight*real(dUmat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION))+Umat_loc(ib,i,IDIR,ION)*dconjg(dUmat_loc(ib,i,IDIR,ION)))!/dosef
              !dEkin(ION,IDIR)=dEkin(ION,IDIR)+(eval0(ispin,ik,i)-mu)*deig(nbmin+ib-1,ik,IDIR,ION)*pdos(ispin,1,ik,i)*sweight!*real(Umat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION)))
              !dEkin(ION,IDIR)=dEkin(ION,IDIR)+eigvals(nbmin+ib-1,ik)*kweight(ispin,ik,i)*sweight*real(Umat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION)))
              !dEkin(ION,IDIR)=dEkin(ION,IDIR)+deig(nbmin+ib-1,ik,IDIR,ION)*kweight(ispin,ik,i)*sweight*real(Umat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION)))
              !dEkin(ION,IDIR)=dEkin(ION,IDIR)+eigvals(nbmin+ib-1,ik)*kweight(ispin,ik,i)*sweight*real(dUmat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION))+Umat_loc(ib,i,IDIR,ION)*dconjg(dUmat_loc(ib,i,IDIR,ION)))
!              dEkin(ION,IDIR)=dEkin(ION,IDIR)+sig(nbmin+ib-1,ik)*kweight(ispin,ik,i)*sweight*real(dUmat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION))+Umat_loc(ib,i,IDIR,ION)*dconjg(dUmat_loc(ib,i,IDIR,ION)))
                !dEkin(ION,IDIR)=dEkin(ION,IDIR)+deps*kweight(ispin,ik,i)*sweight
              enddo
            enddo
          enddo
        endif
            !if (lforce.eq..true.) dEkin=dEkin+eigvals(nbmin+ib-1,ik)*kweight(ik,i)*sweight*real(dUmat_loc(ib,i)*dconjg(Umat_loc(ib,i))+Umat_loc(ib,i)*dconjg(dUmat_loc(ib,i)))
        !if (lforce.eq..true.) then
        !  do ION=1,nions
        !    do IDIR=1,ndir 
        !      do i=1,n_atoms
        !        do j=1,n_orbs
        !          do ib=1,num_wann
        !            if (sym_idx(ispin,i,j)>0) dEkin(ION,IDIR)=dEkin(ION,IDIR)+Sigoo(sym_idx(ispin,i,j))*kweight(ispin,ik,ib)*sweight*real(dUmat_loc((i-1)*n_orbs+j,ib,IDIR,ION)*dconjg(Umat_loc((i-1)*n_orbs+j,ib,IDIR,ION))+Umat_loc((i-1)*n_orbs+j,ib,IDIR,ION)*dconjg(dUmat_loc((i-1)*n_orbs+j,ib,IDIR,ION)))
        !          enddo
        !        enddo
        !      enddo
        !    enddo
        !  enddo
        !endif
        if (allocated(Umat)) deallocate (Umat_loc,stat=ierr)
        if (allocated(Umat_loc)) deallocate (Umat_loc,stat=ierr)
        if (allocated(dUmat_loc)) deallocate (dUmat_loc,stat=ierr)
      enddo
      loc_i=loc_i+1
    enddo
    do ib=1,num_wann
      tot_n=tot_n+occ_wann(ib)
    enddo
    do ib=1,5
      Nd=Nd+occ_wann(ib)
    enddo
    
!
    if (lforce.eq..true.) call comms_allreduce(dmu,nions,ndir,'SUM')
    call comms_allreduce(occ_wann(1),num_wann,'SUM')
    call comms_allreduce(mom(1),num_wann,'SUM')
    call comms_allreduce(Ekin,1,'SUM')
    call comms_allreduce(tot_n,1,'SUM')
    call comms_allreduce(Nd,1,'SUM')


    dSig=0.0_dp
    Niter=10 
    if (lforce.eq..true.) then
      do iter=1,Niter
        dN=0.0_dp
        dN2=0.0_dp
        loc_i=1
        do ik=displs(my_node_id)+1,displs(my_node_id)+numk_loc
          nbmin=band_win(1,ik)
          nbmax=band_win(2,ik)
          num_band_max=nbmax-nbmin+1
          do ispin=1,nspin
            if (.not. allocated(Umat)) then
              allocate (Umat(num_band_max,num_wann), stat=ierr)
              if (ierr /= 0) call io_error('Error allocating Gloc_sum in compute_G_loc')
            endif
            Umat=cmplx_0
            Umat=matmul(UMatrix(:,:,ik),evec0(ispin,loc_i,:,:))
            if (.not. allocated(Umat_loc)) then
              allocate (Umat_loc(num_band_max,num_wann,ndir,nions), stat=ierr)
              if (ierr /= 0) call io_error('Error allocating Gloc_sum in compute_G_loc')
            endif
            if (.not. allocated(dUmat_loc)) then
              allocate (dUmat_loc(num_band_max,num_wann,ndir,nions), stat=ierr)
              if (ierr /= 0) call io_error('Error allocating Gloc_sum in compute_G_loc')
            endif
            Umat_loc=cmplx_0
            dUmat_loc=cmplx_0
            do ION=1,nions
              do IDIR=1,ndir   
                Umat_loc(:,:,IDIR,ION)=matmul(UMatrix_loc(:,:,ik,IDIR,ION),evec0(ispin,loc_i,:,:))
                dUmat_loc(:,:,IDIR,ION)=matmul(dUMatrix(:,:,ik,IDIR,ION),evec0(ispin,loc_i,:,:))
                !dUmat_loc(:,:,IDIR,ION)=matmul(dconjg(transpose(dUMatrix(:,:,ik,IDIR,ION))),Umat_loc(:,:,IDIR,ION))
              enddo
            enddo
            do ION=1,nions
              do IDIR=1,ndir 
                do i=1,num_wann
                  deps=0.0_dp
                  eps=0.0_dp
                  do ib=1,num_band_max
                    !deps=deps+real(Umat(ib,i)*dconjg(Umat(ib,i)))
                    !deps=deps+real(dUmat_loc(ib,i,IDIR,ION)*dconjg(dUmat_loc(ib,i,IDIR,ION)))
                    !deps=deps+real(dUmat_loc(ib,i,IDIR,ION)*dconjg(dUmat_loc(ib,i,IDIR,ION)))
                    !deps=deps+eigvals(nbmin+ib-1,ik)*real(dUmat_loc(ib,i,IDIR,ION)*dconjg(dUmat_loc(ib,i,IDIR,ION)))
                    !deps=deps+eigvals(nbmin+ib-1,ik)*real(dUmat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION))+Umat_loc(ib,i,IDIR,ION)*dconjg(dUmat_loc(ib,i,IDIR,ION)))+deig(nbmin+ib-1,ik,IDIR,ION)*real(Umat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION)))
                    !deps=deps+eigvals(nbmin+ib-1,ik)*real(dUmat_loc(ib,i,IDIR,ION)*dconjg(Umat(ib,i,IDIR,ION))+Umat(ib,i,IDIR,ION)*dconjg(dUmat_loc(ib,i,IDIR,ION)))+deig(nbmin+ib-1,ik,IDIR,ION)*real(Umat(ib,i,IDIR,ION)*dconjg(Umat(ib,i,IDIR,ION)))
                    !deps=deps+eigvals(nbmin+ib-1,ik)*real(dUmat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION))+Umat_loc(ib,i,IDIR,ION)*dconjg(dUmat_loc(ib,i,IDIR,ION)))
                    !deps=deps+eigvals(nbmin+ib-1,ik)*real(Umat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION)))!-dUmat_loc(ib,i,IDIR,ION)*dconjg(dUmat_loc(ib,i,IDIR,ION)))/0.0002_dp
                    !deps=deps+eigvals(nbmin+ib-1,ik)*real(Umat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION)))!-dUmat_loc(ib,i,IDIR,ION)*dconjg(dUmat_loc(ib,i,IDIR,ION)))
                    deps=deps+deig(nbmin+ib-1,ik,IDIR,ION)*real(Umat(ib,i)*dconjg(Umat(ib,i)))
                    !eps=eps+eigvals(nbmin+ib-1,ik)*real(Umat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION)))
                    !dmu(ION,IDIR)=dmu(ION,IDIR)+pdos(ispin,1,ik,i)*sweight*real(Umat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION)))
                    !dEkin(ION,IDIR)=dEkin(ION,IDIR)+pdos(ispin,1,ik,i)*sweight*real(Umat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION)))
                    !dmu(ION,IDIR)=dmu(ION,IDIR)+deig(nbmin+ib-1,ik,IDIR,ION)*pdos(ispin,1,ik,i)*sweight*real(Umat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION)))/dosef
                    !dmu(ION,IDIR)=dmu(ION,IDIR)+eigvals(nbmin+ib-1,ik)*pdos(ispin,1,ik,i)*sweight*real(dUmat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION))+Umat_loc(ib,i,IDIR,ION)*dconjg(dUmat_loc(ib,i,IDIR,ION)))/dosef
                  enddo
                  do ib=1,n_atoms
                    do j=1,n_orbs
                      deps=deps+dSig(ispin,ib,ION,IDIR)*real(evec0(ispin,loc_i,(ib-1)*n_orbs+j,i)*dconjg(evec0(ispin,loc_i,(ib-1)*n_orbs+j,i)))
                    enddo 
                  enddo
                  !dEkin(ION,IDIR)=dEkin(ION,IDIR)+deps*kweight(ispin,ik,i)*sweight
                  !do ib=1,num_band_max
                  !  dEkin(ION,IDIR)=dEkin(ION,IDIR)+eigvals(nbmin+ib-1,ik)*kweight(ispin,ik,i)*sweight*real(dUmat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION))+Umat_loc(ib,i,IDIR,ION)*dconjg(dUmat_loc(ib,i,IDIR,ION)))
                  !  !dEkin(ION,IDIR)=dEkin(ION,IDIR)+eigvals(nbmin+ib-1,ik)*kweight(ispin,ik,i)*sweight*real(Umat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION)))
                  !  dEkin(ION,IDIR)=dEkin(ION,IDIR)+deig(nbmin+ib-1,ik,IDIR,ION)*kweight(ispin,ik,i)*sweight*real(Umat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION)))
                  !enddo
                  !do ib=1,n_atoms
                  !  do j=1,n_orbs
                  !    !if (sym_idx(ispin,ib,j)>0) dEkin(ION,IDIR)=dEkin(ION,IDIR)+(deps-dmu(ION,IDIR))*pdos(ispin,1,ik,i)*sweight*real(Umat_loc((ib-1)*n_orbs+j,i,IDIR,ION)*dconjg(Umat_loc((ib-1)*n_orbs+j,i,IDIR,ION)))*Sigoo(sym_idx(ispin,ib,j))
                  !    if (sym_idx(ispin,ib,j)>0) dEkin(ION,IDIR)=dEkin(ION,IDIR)+kweight(ispin,ik,i)*sweight*real(dUmat_loc((ib-1)*n_orbs+j,i,IDIR,ION)*dconjg(Umat((ib-1)*n_orbs+j,i))+Umat((ib-1)*n_orbs+j,i)*dconjg(dUmat_loc((ib-1)*n_orbs+j,i,IDIR,ION)))*Sigoo(sym_idx(ispin,ib,j))
                  !  enddo
                  !enddo
                  do ib=1,n_atoms*n_orbs
                    dN(ispin,ib,ION,IDIR)=dN(ispin,ib,ION,IDIR)+(deps-dmu(ION,IDIR))*pdos(ispin,1,ik,i)*sweight*real(evec0(ispin,loc_i,ib,i)*dconjg(evec0(ispin,loc_i,ib,i)))!*Sigoo(sym_idx(ispin,ib,j)) 
                  enddo
                  !do ib=1,num_band_max
                  !  dN2(ib,ION,IDIR)=dN(ib,ION,IDIR)+(deps-dmu(ION,IDIR))*pdos(ispin,1,ik,i)*sweight*real(Umat(ib,i)*dconjg(Umat(ib,i)))!*Sigoo(sym_idx(ispin,ib,j)) 
                     !dEkin(ION,IDIR)=dEkin(ION,IDIR)+(deps-dmu(ION,IDIR))*pdos(ispin,1,ik,i)*sweight*real(evec0(ispin,loc_i,ib,i)*dconjg(evec0(ispin,loc_i,ib,i)))!*Sigoo(sym_idx(ispin,ib,j)) 
                    !dEkin(ION,IDIR)=dEkin(ION,IDIR)+eps*pdos(ispin,1,ik,i)*sweight*real(dUmat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION))+Umat_loc(ib,i,IDIR,ION)*dconjg(dUmat_loc(ib,i,IDIR,ION)))
                    !dEkin(ION,IDIR)=dEkin(ION,IDIR)+eigvals(nbmin+ib-1,ik)*pdos(ispin,1,ik,i)*sweight*real(dUmat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION))+Umat_loc(ib,i,IDIR,ION)*dconjg(dUmat_loc(ib,i,IDIR,ION)))
                    !dEkin(ION,IDIR)=dEkin(ION,IDIR)+eigvals(nbmin+ib-1,ik)*pdos(ispin,1,ik,i)*sweight*real(dUmat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION))+Umat_loc(ib,i,IDIR,ION)*dconjg(dUmat_loc(ib,i,IDIR,ION)))
                    !dEkin(ION,IDIR)=dEkin(ION,IDIR)+deig(nbmin+ib-1,ik,IDIR,ION)*pdos(ispin,1,ik,i)*sweight*real(Umat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION)))
                    !if (ib.eq.i) dEkin(ION,IDIR)=dEkin(ION,IDIR)-dmu(ION,IDIR)*pdos(ispin,1,ik,i)*sweight!*real(Umat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION)))
                    !dEkin(ION,IDIR)=dEkin(ION,IDIR)+(eval0(ispin,ik,i)-mu)*deig(nbmin+ib-1,ik,IDIR,ION)*pdos(ispin,1,ik,i)*sweight!*real(Umat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION)))
!                    dEkin(ION,IDIR)=dEkin(ION,IDIR)+sig(nbmin+ib-1,ik)*kweight(ispin,ik,i)*sweight*real(dUmat_loc(ib,i,IDIR,ION)*dconjg(Umat_loc(ib,i,IDIR,ION))+Umat_loc(ib,i,IDIR,ION)*dconjg(dUmat_loc(ib,i,IDIR,ION)))
                !  enddo
                enddo
                !if (lforce.eq..true.) dEkin=dEkin+eigvals(nbmin+ib-1,ik)*kweight(ik,i)*sweight*real(dUmat_loc(ib,i)*dconjg(Umat_loc(ib,i))+Umat_loc(ib,i)*dconjg(dUmat_loc(ib,i)))
              enddo
            enddo
          enddo
          loc_i=loc_i+1
        enddo
        call comms_allreduce(dN,nspin,n_atoms*n_orbs,nions,ndir,'SUM')
        call comms_allreduce(dN2,num_bands,nions,ndir,'SUM')
        if (on_root) print *, iter, dmu(5,3), dN(2,1,5,3)
        dSig=0.0_dp
        do ispin=1,nspin
          do ION=1,nions
            do IDIR=1,ndir 
              do ib=1,n_atoms
                do j=1,n_orbs
                  dSig(ispin,ib,ION,IDIR)=dSig(ispin,ib,ION,IDIR)+5.0*dN(ispin,(ib-1)*n_orbs+j,ION,IDIR)
                enddo 
              enddo
            enddo
          enddo
        enddo
        loc_i=1
        dmu=0.0_dp
        do ik=displs(my_node_id)+1,displs(my_node_id)+numk_loc
          do ispin=1,nspin
            Umat=cmplx_0
            Umat=matmul(UMatrix(:,:,ik),evec0(ispin,loc_i,:,:))
            do ION=1,nions
              do IDIR=1,ndir 
                do i=1,num_wann
                  deps=0.0_dp
                  eps=0.0_dp
                  do ib=1,num_band_max
                    deps=deps+deig(nbmin+ib-1,ik,IDIR,ION)*real(Umat(ib,i)*dconjg(Umat(ib,i)))
                    eps=eps+eigvals(nbmin+ib-1,ik)*real(Umat(ib,i)*dconjg(Umat(ib,i)))
                  enddo
                  do ib=1,n_atoms
                    do j=1,n_orbs
                      deps=deps+dSig(ispin,ib,ION,IDIR)*real(evec0(ispin,loc_i,(ib-1)*n_orbs+j,i)*dconjg(evec0(ispin,loc_i,(ib-1)*n_orbs+j,i)))
                    enddo 
                  enddo
                  dmu(ION,IDIR)=dmu(ION,IDIR)+deps*pdos(ispin,1,ik,i)*sweight/dosef
                enddo
              enddo
            enddo
          enddo
          loc_i=loc_i+1
        enddo
        call comms_allreduce(dmu,nions,ndir,'SUM')
      enddo
      dEkin=0.0_dp
      do ION=1,nions
        do IDIR=1,ndir
          do ispin=1,nspin
            do ib=1,n_atoms
              do j=1,n_orbs
                dEkin(ION,IDIR)=dEkin(ION,IDIR)+Sigoo(sym_idx(ispin,ib,j))*dN(ispin,(ib-1)*n_orbs+j,ION,IDIR)
              enddo
            enddo
          enddo
        enddo
      enddo
    endif
    if (allocated(Umat)) deallocate (Umat,stat=ierr)
    if (allocated(Umat_loc)) deallocate (Umat_loc,stat=ierr)
    if (allocated(dUmat_loc)) deallocate (dUmat_loc,stat=ierr)



    if (allocated(eval0)) deallocate (eval0,stat=ierr)
    if (allocated(evec0)) deallocate (evec0,stat=ierr)
    if (allocated(kweight)) deallocate (kweight)
    if (allocated(pdos)) deallocate (pdos)
    if (allocated(cpdos)) deallocate (cpdos)


    if (on_root) then   
      OPEN(UNIT=99,FILE='INFO_KSUM',STATUS='old',FORM='FORMATTED',ACCESS='APPEND')
      WRITE(99,'(11F12.6)') mu, tot_n, SUM(occ_wann(1:n_orbs)), &
      SUM(occ_wann(n_orbs*(n_atoms-1)+1:n_orbs*n_atoms)), Ekin, Sigoo(1), dEkin(1,1)
      CLOSE(99)

      OPEN(UNIT=99,FILE='INFO_DM',STATUS='old',FORM='FORMATTED',ACCESS='APPEND')
90    FORMAT('(',I2,'F20.15)')
      WRITE(write_format,90) n_atoms*n_orbs
      WRITE(99,write_format,advance='no') (occ_wann(ib), ib=1,n_atoms*n_orbs)
      WRITE(99,write_format) (mom(ib), ib=1,n_atoms*n_orbs)
      CLOSE(99)
      if (lforce.eq..true.) then
        OPEN(UNIT=99,FILE='FORCE',FORM='FORMATTED')
        do ION=1,nions
          WRITE(99,'(3F12.8)') dN(1,1,ION,1), dN(1,1,ION,2), dN(1,1,ION,3)
          !WRITE(99,'(3F12.8)') dSig(1,1,ION,1), dSig(1,1,ION,2), dSig(1,1,ION,3)
          !WRITE(99,'(3F12.8)') dEkin(ION,1), dEkin(ION,2), dEkin(ION,3)
          !WRITE(99,'(3F12.8)') dmu(ION,1), dmu(ION,2), dmu(ION,3)
        enddo
        CLOSE(99)
      endif
    endif
    if (allocated(occ_wann)) deallocate (occ_wann)
    if (allocated(mom)) deallocate (mom)
  end subroutine compute_tetra

  subroutine compute_tetra_EDM(n_kpts_loc,kpt_dft,wght_dft)
    use constants
    use io
    use utility
    use tetrahedron
    use generate_ham, only: tran
    use comms, only: on_root,comms_allreduce,my_node_id, num_nodes, comms_array_split,comms_reduce, comms_barrier
   
    implicit none

    integer, intent(in) :: n_kpts_loc
    real(kind=dp), intent(in) :: kpt_dft(3,n_kpts_loc)
    real(kind=dp), intent(in) :: wght_dft(n_kpts_loc)
    character(len=10) :: write_format
    real(kind=dp) :: tot_n
    complex(kind=dp) :: Gs, Gs0
    integer :: i,j,l,ik,r,ib,im,ispin,loc_i,ierr,x,y,z,kx_loc,ky_loc,kz_loc
    integer :: nbmin,nbmax,num_band_max 
    integer :: numk_loc,nfine_tot, ntet_loc, nibz_loc
    !integer,intent(in) :: n_kpts_loc,n_wann
    ! Needed to split an array on different nodes
    integer, dimension(0:num_nodes - 1) :: counts
    integer, dimension(0:num_nodes - 1) :: displs
    integer, dimension(0:num_nodes - 1) :: counts_tet
    integer, dimension(0:num_nodes - 1) :: displs_tet
    integer, dimension(0:num_nodes - 1) :: counts_ibz
    integer, dimension(0:num_nodes - 1) :: displs_ibz
    real(kind=dp) :: rdotk, fermi0, tot
    real(kind=dp) :: Ekin, Nd, sweight,dm_loc
    real(kind=dp), allocatable :: eval0(:,:)
    real(kind=dp), allocatable :: Ed(:)
    real(kind=dp), allocatable :: kweight(:,:)
    real(kind=dp), allocatable :: occ_wann(:)
    real(kind=dp), allocatable :: mom(:)
    complex(kind=dp), allocatable :: Umat_loc(:,:)
    complex(kind=dp), allocatable :: Hk(:,:)
    complex(kind=dp), allocatable :: evec0(:,:,:)
    complex(kind=dp), allocatable :: Umat(:,:)
    complex(kind=dp) :: DMFT_UU,DMFT_U0 

    call comms_array_split(num_new_kpts, counts, displs)
    numk_loc=counts(my_node_id)
    call comms_array_split(num_new_kpts*6, counts_tet, displs_tet)
    ntet_loc=counts_tet(my_node_id)
    call comms_array_split(n_kpts_loc, counts_ibz, displs_ibz)    
    nibz_loc=counts_ibz(my_node_id)

    if (.not. allocated(eval0)) then
      allocate (eval0(num_new_kpts,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval0 in compute_DMFT_mu')
    endif
    if (.not. allocated(evec0)) then
      allocate (evec0(num_new_kpts,num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval0 in compute_DMFT_mu')
    endif
    if (.not. allocated(kweight)) then
      allocate (kweight(num_new_kpts,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval0 in compute_DMFT_mu')
    endif


    if (.not. allocated(occ_wann)) then
      allocate (occ_wann(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Gloc_sum in compute_G_loc')
    endif
    if (.not. allocated(DMFT_EDM_loc)) then
      allocate (DMFT_EDM_loc(num_bands,num_bands,n_kpts_loc), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating DMFT_eval_loc in compute_G_loc')
    endif
!    if (.not. allocated(DMFT_evec_loc)) then
!      allocate (DMFT_evec_loc(num_bands,num_wann,n_kpts_loc), stat=ierr)
!      if (ierr /= 0) call io_error('Error allocating DMFT_evec_loc in compute_G_loc')
!    endif

!    DMFT_eval_loc=0.0_dp; 
    DMFT_EDM_loc=cmplx_0

    sweight=2.0_dp/dfloat(nspin)/dfloat(num_new_kpts)
    tot_n=0.0_dp;Nd=0.0_dp;Ekin=0.0_dp
    do ispin=1,nspin
      occ_wann=0.0_dp
      eval0=0.0_dp
      evec0=cmplx_0
      kweight=0.0_dp
      do ik=displs(my_node_id)+1,displs(my_node_id)+numk_loc
        nbmin=band_win(1,ik)
        nbmax=band_win(2,ik)
        num_band_max=nbmax-nbmin+1
        do x=1,num_wann
          do y=1,num_wann
            do z=1,num_band_max
              evec0(ik,x,y)=evec0(ik,x,y)+dconjg(UMatrix(z,x,ik))*UMatrix(z,y,ik)*eigvals(nbmin+z-1,ik)
            enddo
          enddo
        enddo
        do i=1,n_atoms
          do j=1,n_orbs
            if (sym_idx(ispin,i,j)>0) evec0(ik,(i-1)*n_orbs+j,(i-1)*n_orbs+j)=evec0(ik,(i-1)*n_orbs+j,(i-1)*n_orbs+j)+Sigoo(sym_idx(ispin,i,j))
          enddo
        enddo
        CALL EIGENH(evec0(ik,:,:),num_wann,eval0(ik,:))
      enddo

      call comms_allreduce(eval0,num_new_kpts,num_wann,'SUM')
      call comms_allreduce(evec0,num_new_kpts,num_wann,num_wann,'SUM')

      call tetra_poccup(num_wann,num_new_kpts,num_new_kpts*6,displs_tet(my_node_id)+1,displs_tet(my_node_id)+ntet_loc,kibz,tetptr,tet_idx,eval0(:,:),mu,kweight(:,:))
      call comms_allreduce(kweight,num_new_kpts,num_wann,'SUM')


      do ik=displs(my_node_id)+1,displs(my_node_id)+numk_loc
        do ib=1,num_wann
          do i=1,num_wann
            dm_loc=real(evec0(ik,i,ib)*dconjg(evec0(ik,i,ib)))*kweight(ik,ib)*sweight
            occ_wann(i)=occ_wann(i)+dm_loc
          enddo
        enddo
        nbmin=band_win(1,ik)
        nbmax=band_win(2,ik)
        num_band_max=nbmax-nbmin+1
        if (.not. allocated(Umat_loc)) then
          allocate (Umat_loc(num_band_max,num_wann), stat=ierr)
          if (ierr /= 0) call io_error('Error allocating Gloc_sum in compute_G_loc')
        endif
        Umat_loc=cmplx_0
        Umat_loc=matmul(UMatrix(:,:,ik),evec0(ik,:,:))
        do ib=1,num_band_max
          do i=1,num_wann
            Ekin=Ekin+eigvals(nbmin+ib-1,ik)*kweight(ik,i)*sweight*real(Umat_loc(ib,i)*dconjg(Umat_loc(ib,i)))
          enddo
        enddo
        if (allocated(Umat_loc)) deallocate (Umat_loc,stat=ierr)
      enddo
!      do ib=1,num_wann
!        tot_n=tot_n+occ_wann(ib)
!      enddo
!      do ib=1,5
!        Nd=Nd+occ_wann(ib)
!      enddo

      do ik=displs_ibz(my_node_id)+1,displs_ibz(my_node_id)+nibz_loc
        kx_loc=MODULO(NINT(kpt_dft(1,ik)*mp_grid(1)-0.25),mp_grid(1))
        ky_loc=MODULO(NINT(kpt_dft(2,ik)*mp_grid(2)-0.25),mp_grid(2))
        kz_loc=MODULO(NINT(kpt_dft(3,ik)*mp_grid(3)-0.25),mp_grid(3))
        do i=1,num_new_kpts
          if (kx_loc.eq.ikpt(1,i) .and. ky_loc.eq.ikpt(2,i) .and. kz_loc.eq.ikpt(3,i)) then
            loc_i=i
            exit
          endif
        enddo 
        nbmin=band_win(1,loc_i)
        if (.not. allocated(Umat)) then
          allocate (Umat(num_band_max,num_wann), stat=ierr)
          if (ierr /= 0) call io_error('Error allocating Gloc_sum in compute_G_loc')
        endif
        Umat=cmplx_0
        Umat=matmul(UMatrix(:,:,loc_i),evec0(loc_i,:,:))
        do ib=1,num_wann
          do i=1,num_band_max
            DMFT_EDM_loc(i,i,ik)=DMFT_EDM_loc(i,i,ik)-eigvals(nbmin+i-1,loc_i)*Umat(i,ib)*dconjg(Umat(i,ib))*kweight(loc_i,ib)/dfloat(nspin)
            tot_n=tot_n+2.0_dp*Umat(i,ib)*dconjg(Umat(i,ib))*kweight(loc_i,ib)*wght_dft(ik)/dfloat(nspin)
            if (i .le. 5) Nd=Nd+2.0_dp*evec0(loc_i,i,ib)*dconjg(evec0(loc_i,i,ib))*kweight(loc_i,ib)*wght_dft(ik)/dfloat(nspin)
            do j=1,num_band_max
              DMFT_EDM_loc(i,j,ik)=DMFT_EDM_loc(i,j,ik)+Umat(i,ib)*dconjg(Umat(j,ib))*eval0(loc_i,ib)*kweight(loc_i,ib)/dfloat(nspin)
            enddo
          enddo
        enddo
      enddo
    enddo

    if (allocated(eval0)) deallocate (eval0,stat=ierr)
    if (allocated(evec0)) deallocate (evec0,stat=ierr)
    if (allocated(kweight)) deallocate (kweight)
    if (allocated(occ_wann)) deallocate (occ_wann)
!
    call comms_allreduce(Ekin,1,'SUM')
    call comms_allreduce(tot_n,1,'SUM')
    call comms_allreduce(Nd,1,'SUM')


    call comms_allreduce(DMFT_EDM_loc,num_bands,num_bands,n_kpts_loc,'SUM')
    
    if (allocated(UMatrix)) deallocate (UMatrix,stat=ierr)


   if (on_root) then
      OPEN(UNIT=99,FILE='INFO_DFT_loop',STATUS='old',FORM='FORMATTED',ACCESS='APPEND')
      WRITE(99,'(4F12.6)') mu, tot_n, Nd, Ekin 
      CLOSE(99)
    endif
  end subroutine compute_tetra_EDM

  subroutine compute_tetra_DM(n_kpts_loc,kpt_dft,wght_dft)
    use constants
    use io
    use utility
    use tetrahedron
    use generate_ham, only: tran
    use comms, only: on_root,comms_allreduce,my_node_id, num_nodes, comms_array_split,comms_reduce, comms_barrier
   
    implicit none

    integer, intent(in) :: n_kpts_loc
    real(kind=dp), intent(in) :: kpt_dft(3,n_kpts_loc)
    real(kind=dp), intent(in) :: wght_dft(n_kpts_loc)
    character(len=10) :: write_format
    real(kind=dp) :: tot_n
    complex(kind=dp) :: Gs, Gs0
    integer :: i,j,l,ik,r,ib,im,ispin,loc_i,ierr,x,y,z,kx_loc,ky_loc,kz_loc
    integer :: nbmin,nbmax,num_band_max 
    integer :: numk_loc,nfine_tot, ntet_loc, nibz_loc
    !integer,intent(in) :: n_kpts_loc,n_wann
    ! Needed to split an array on different nodes
    integer, dimension(0:num_nodes - 1) :: counts
    integer, dimension(0:num_nodes - 1) :: displs
    integer, dimension(0:num_nodes - 1) :: counts_tet
    integer, dimension(0:num_nodes - 1) :: displs_tet
    integer, dimension(0:num_nodes - 1) :: counts_ibz
    integer, dimension(0:num_nodes - 1) :: displs_ibz
    real(kind=dp) :: rdotk, fermi0, tot
    real(kind=dp) :: Ekin, Nd, sweight,dm_loc
    real(kind=dp), allocatable :: eval0(:,:)
    real(kind=dp), allocatable :: Ed(:)
    real(kind=dp), allocatable :: kweight(:,:)
    real(kind=dp), allocatable :: occ_wann(:)
    real(kind=dp), allocatable :: mom(:)
    complex(kind=dp), allocatable :: Umat_loc(:,:)
    complex(kind=dp), allocatable :: Hk(:,:)
    complex(kind=dp), allocatable :: evec0(:,:,:)
    complex(kind=dp) :: DMFT_UU,DMFT_U0 

    call comms_array_split(num_new_kpts, counts, displs)
    numk_loc=counts(my_node_id)
    call comms_array_split(num_new_kpts*6, counts_tet, displs_tet)
    ntet_loc=counts_tet(my_node_id)
    call comms_array_split(n_kpts_loc, counts_ibz, displs_ibz)    
    nibz_loc=counts_ibz(my_node_id)

    if (.not. allocated(eval0)) then
      allocate (eval0(num_new_kpts,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval0 in compute_DMFT_mu')
    endif
    if (.not. allocated(evec0)) then
      allocate (evec0(num_new_kpts,num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval0 in compute_DMFT_mu')
    endif
!    if (.not. allocated(Hk)) then
!      allocate (Hk(num_wann,num_wann), stat=ierr)
!      if (ierr /= 0) call io_error('Error allocating Hk in compute_G_loc')
!    endif
    if (.not. allocated(kweight)) then
      allocate (kweight(num_new_kpts,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval0 in compute_DMFT_mu')
    endif


    if (.not. allocated(occ_wann)) then
      allocate (occ_wann(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Gloc_sum in compute_G_loc')
    endif
    if (.not. allocated(DMFT_eval_loc)) then
      allocate (DMFT_eval_loc(num_wann,n_kpts_loc), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating DMFT_eval_loc in compute_G_loc')
    endif
!    if (.not. allocated(DMFT_evec_loc)) then
!      allocate (DMFT_evec_loc(num_bands,num_wann,n_kpts_loc), stat=ierr)
!      if (ierr /= 0) call io_error('Error allocating DMFT_evec_loc in compute_G_loc')
!    endif

    DMFT_eval_loc=0.0_dp; 
!    DMFT_evec_loc=cmplx_0

    sweight=2.0_dp/dfloat(nspin)/dfloat(num_new_kpts)
    tot_n=0.0_dp;Nd=0.0_dp;Ekin=0.0_dp
    do ispin=1,nspin
      occ_wann=0.0_dp
      eval0=0.0_dp
      evec0=cmplx_0
      kweight=0.0_dp
      do ik=displs(my_node_id)+1,displs(my_node_id)+numk_loc
        nbmin=band_win(1,ik)
        nbmax=band_win(2,ik)
        num_band_max=nbmax-nbmin+1
        do x=1,num_wann
          do y=1,num_wann
            do z=1,num_band_max
              evec0(ik,x,y)=evec0(ik,x,y)+dconjg(UMatrix(z,x,ik))*UMatrix(z,y,ik)*eigvals(nbmin+z-1,ik)
            enddo
          enddo
        enddo
        do i=1,n_atoms
          do j=1,n_orbs
            if (sym_idx(ispin,i,j)>0) evec0(ik,(i-1)*n_orbs+j,(i-1)*n_orbs+j)=evec0(ik,(i-1)*n_orbs+j,(i-1)*n_orbs+j)+Sigoo(sym_idx(ispin,i,j))
          enddo
        enddo
        CALL EIGENH(evec0(ik,:,:),num_wann,eval0(ik,:))
      enddo

      call comms_allreduce(eval0,num_new_kpts,num_wann,'SUM')
      call comms_allreduce(evec0,num_new_kpts,num_wann,num_wann,'SUM')

      call tetra_poccup(num_wann,num_new_kpts,num_new_kpts*6,displs_tet(my_node_id)+1,displs_tet(my_node_id)+ntet_loc,kibz,tetptr,tet_idx,eval0(:,:),mu,kweight(:,:))
      call comms_allreduce(kweight,num_new_kpts,num_wann,'SUM')


      do ik=displs(my_node_id)+1,displs(my_node_id)+numk_loc
        do ib=1,num_wann
          do i=1,num_wann
            dm_loc=real(evec0(ik,i,ib)*dconjg(evec0(ik,i,ib)))*kweight(ik,ib)*sweight
            occ_wann(i)=occ_wann(i)+dm_loc
          enddo
        enddo
        nbmin=band_win(1,ik)
        nbmax=band_win(2,ik)
        num_band_max=nbmax-nbmin+1
        if (.not. allocated(Umat_loc)) then
          allocate (Umat_loc(num_band_max,num_wann), stat=ierr)
          if (ierr /= 0) call io_error('Error allocating Gloc_sum in compute_G_loc')
        endif
        Umat_loc=cmplx_0
        Umat_loc=matmul(UMatrix(:,:,ik),evec0(ik,:,:))
        do ib=1,num_band_max
          do i=1,num_wann
            Ekin=Ekin+eigvals(nbmin+ib-1,ik)*kweight(ik,i)*sweight*real(Umat_loc(ib,i)*dconjg(Umat_loc(ib,i)))
          enddo
        enddo
        if (allocated(Umat_loc)) deallocate (Umat_loc,stat=ierr)
      enddo
!      do ib=1,num_wann
!        tot_n=tot_n+occ_wann(ib)
!      enddo
!      do ib=1,5
!        Nd=Nd+occ_wann(ib)
!      enddo

      do ik=displs_ibz(my_node_id)+1,displs_ibz(my_node_id)+nibz_loc
        kx_loc=MODULO(NINT(kpt_dft(1,ik)*mp_grid(1)-0.25),mp_grid(1))
        ky_loc=MODULO(NINT(kpt_dft(2,ik)*mp_grid(2)-0.25),mp_grid(2))
        kz_loc=MODULO(NINT(kpt_dft(3,ik)*mp_grid(3)-0.25),mp_grid(3))
        do i=1,num_new_kpts
          if (kx_loc.eq.ikpt(1,i) .and. ky_loc.eq.ikpt(2,i) .and. kz_loc.eq.ikpt(3,i)) then
            loc_i=i
            exit
          endif
        enddo 
        do ib=1,num_wann
          do i=1,num_wann
            tot_n=tot_n+2.0_dp*evec0(loc_i,i,ib)*dconjg(evec0(loc_i,i,ib))*kweight(loc_i,ib)*wght_dft(ik)/dfloat(nspin)
            if (i .le. 5) Nd=Nd+2.0_dp*evec0(loc_i,i,ib)*dconjg(evec0(loc_i,i,ib))*kweight(loc_i,ib)*wght_dft(ik)/dfloat(nspin)
            do j=1,num_wann
              DMFT_UU=evec0(loc_i,i,ib)*dconjg(evec0(loc_i,j,ib))*kweight(loc_i,ib)/dfloat(nspin)
              DMFT_evec_loc(i,j,ik)=DMFT_evec_loc(i,j,ik)+DMFT_UU
            enddo
          enddo
        enddo
      enddo
    enddo
    if (allocated(eval0)) deallocate (eval0,stat=ierr)
    if (allocated(evec0)) deallocate (evec0,stat=ierr)
    if (allocated(kweight)) deallocate (kweight)
    if (allocated(occ_wann)) deallocate (occ_wann)
!
    call comms_allreduce(Ekin,1,'SUM')
    call comms_allreduce(tot_n,1,'SUM')
    call comms_allreduce(Nd,1,'SUM')


    do ik=displs_ibz(my_node_id)+1,displs_ibz(my_node_id)+nibz_loc
      call EIGENH(DMFT_evec_loc(1:num_wann,:,ik),num_wann,DMFT_eval_loc(:,ik))
      DMFT_evec_loc(:,:,ik)=MATMUL(UMatrix(:,:,ik),DMFT_evec_loc(1:num_wann,:,ik))
    enddo
    call comms_allreduce(DMFT_evec_loc,num_bands,num_wann,n_kpts_loc,'SUM')
    call comms_allreduce(DMFT_eval_loc,num_wann,n_kpts_loc,'SUM')
    
    if (allocated(UMatrix)) deallocate (UMatrix,stat=ierr)


   if (on_root) then
      OPEN(UNIT=99,FILE='INFO_DFT_loop',STATUS='old',FORM='FORMATTED',ACCESS='APPEND')
      WRITE(99,'(4F12.6)') mu, tot_n, Nd, Ekin 
      CLOSE(99)
    endif
  end subroutine compute_tetra_DM

  subroutine compute_DMFT_mu()
    use constants
    use io
    use utility
    use generate_kpts, only: kpts, weight, num_new_kpts
    use generate_ham, only: tran
    use comms, only: on_root,comms_barrier,comms_reduce,comms_allreduce,my_node_id, num_nodes, comms_array_split, comms_bcast
   
    implicit none

    logical :: iffile
    logical :: LOWB,HIGHB
    real(kind=dp) :: low_mu, high_mu, tot_n
    complex(kind=dp) :: Gs, Gs0, new_eval
    !real(kind=dp) :: Os_real,Oinf_real,Os,Oinf
    !real(kind=dp) :: coeff_a,coeff_b,Os,Oinf
    !real(kind=dp) :: coeff_real,coeff_imag
    !real(kind=dp) :: at, at0
    integer :: i,j,l,ik,r,ib,im,ispin,loc_i,ierr
    integer :: nbmin,nbmax,num_band_max 
    integer :: numk_loc!,mu_iter 
    ! Needed to split an array on different nodes
    integer, dimension(0:num_nodes - 1) :: counts
    integer, dimension(0:num_nodes - 1) :: displs
    real(kind=dp), allocatable :: eval0(:,:,:)
    complex(kind=dp), allocatable :: eval(:,:,:,:)
    complex(kind=dp), allocatable :: Hk_loc(:,:)
    !real(kind=dp), allocatable :: eval0_loc(:)
    !complex(kind=dp), allocatable :: eval_loc(:)
    real(kind=dp) :: rdotk
    complex(kind=dp), allocatable :: Hk(:,:)

!    noms=400

!    write(*,*) num_new_kpts
    !write(*,*) 'hi1'
    !ierr=0 
    call comms_array_split(num_new_kpts, counts, displs)    
    !write(*,*) my_node_id, counts(my_node_id), displs(my_node_id)
    numk_loc=counts(my_node_id)
    !mu_iter=100
    if (.not. allocated(eval0)) then
      allocate (eval0(num_wann,nspin,numk_loc), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval0 in compute_DMFT_mu')
    endif
    if (.not. allocated(eval)) then
      allocate (eval(num_wann,nom,nspin,numk_loc), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval in compute_DMFT_mu')
    endif
    if (.not. allocated(Hk)) then
      allocate (Hk(num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Hk in compute_DMFT_mu')
    endif
    if (.not. allocated(Hk_loc)) then
      allocate (Hk_loc(num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Hk_loc in compute_DMFT_mu')
    endif
    !if (.not. allocated(eval0_loc)) then
    !  allocate (eval0_loc(num_wann), stat=ierr)
    !  if (ierr /= 0) call io_error('Error allocating eval0_loc in compute_DMFT_mu')
    !endif
    !if (.not. allocated(eval_loc)) then
    !  allocate (eval_loc(num_wann), stat=ierr)
    !  if (ierr /= 0) call io_error('Error allocating eval_loc in compute_DMFT_mu')
    !endif

    eval=cmplx_0
    eval0=0.0_dp
    loc_i=0
    do ik=displs(my_node_id)+1,displs(my_node_id)+numk_loc
      loc_i=loc_i+1
      Hk=cmplx_0
      !write(*,*) kpts(1,1),tran(1,1),HamR(1,1,1)
      do r=1,nR
        rdotk=twopi*dot_product(kpts(:,ik),dfloat(tran(:,r)))
        Hk(:,:)=Hk(:,:)+HamR(r,:,:)*exp(cmplx_i*rdotk)
      enddo
      do ispin=1,nspin
        !!!!!!!! Compute the w=infinity value !!!!
        !eval0_loc=0.0_dp
        !Hk_loc=Hk*1.0_dp
        Hk_loc=cmplx_0
        call zcopy(num_wann**2,Hk,1,Hk_loc,1)
        do i=1,n_atoms
          do j=1,n_orbs
            if (sym_idx(ispin,i,j)>0) Hk_loc((i-1)*n_orbs+j,(i-1)*n_orbs+j)=Hk_loc((i-1)*n_orbs+j,(i-1)*n_orbs+j)+Sigoo(sym_idx(ispin,i,j))
          enddo
        enddo
        !do i=1,num_wann
        !  if (sym_idx(ispin,i)>0) Hk_loc(i,i)=Hk_loc(i,i)+Sigoo(sym_idx(ispin,i))
        !!do i=1,ncor_orb
        !!  Hk_loc(i,i)=Hk_loc(i,i)+Sigoo(i)
        !enddo
        CALL EIGENVALH(Hk_loc,num_wann,eval0(:,ispin,loc_i))
        !do ib=1,num_wann
        !!  eval0(ib,ispin,loc_i)=eval0_loc(ib)
        !  write(*,*) loc_i,ib,eval0(ib,ispin,loc_i)
        !enddo
        do im=1,nom
          !eval_loc=cmplx_0
          !Hk_loc=Hk*1.0_dp
          Hk_loc=cmplx_0
          call zcopy(num_wann**2,Hk,1,Hk_loc,1)
          do i=1,n_atoms
            do j=1,n_orbs
              if (sym_idx(ispin,i,j)>0) Hk_loc((i-1)*n_orbs+j,(i-1)*n_orbs+j)=Hk_loc((i-1)*n_orbs+j,(i-1)*n_orbs+j)+Sigma(sym_idx(ispin,i,j),im)
            enddo
          enddo
          !do j=1,num_wann
          !  if (sym_idx(ispin,j)>0) Hk_loc(j,j)=Hk_loc(j,j)+Sigma(sym_idx(ispin,j),im)
          !enddo
          CALL EIGENVAL(Hk_loc,num_wann,eval(:,im,ispin,loc_i))
        enddo
        !Hk_loc=Hk*1.0_dp
        !do j=1,ncor_orb
        !  Hk_loc(j,j)=Hk_loc(j,j)+Sigma(j,nom)-real(Sigma(j,1200))
        !enddo
        !CALL EIGENVAL(Hk_loc,num_wann,eval(:,nom,ispin,loc_i)) 
        !call quicksort_cmplx(eval(:,noms,ispin,loc_i))
        !call quicksort_cmplx(eval(:,nom,ispin,loc_i))
        !!if (on_root) then
        !!  do ib=1,num_wann
        !!    write(*,*) ib, eval(ib,noms,ispin,loc_i)
        !!    write(*,*) ib, eval0(ib,ispin,loc_i)
        !!    write(*,*) ib, eval(ib,nom,ispin,loc_i)
        !!  enddo
        !!endif
        !!STOP

        !do ib=1,num_wann
        !  Os=aimag(eval(ib,noms,ispin,loc_i))
        !  Oinf=aimag(eval(ib,nom,ispin,loc_i))
        !  Os_real=real(eval(ib,noms,ispin,loc_i))
        !  Oinf_real=real(eval(ib,nom,ispin,loc_i))
        !  !if (abs(Os-Oinf).gt.eps5) then
        !  !  coeff_a=(Oinf*om(nom)-Os*om(noms))/(Os-Oinf)
        !  !  coeff_b=(Oinf*Os)/(Os-Oinf)*(om(nom)-om(noms))
        !  !else
        !  !  coeff_a=0.0_dp
        !  !  coeff_b=Os*om(noms)
        !  !endif
        !  do im=noms+1,nom-1
        !    !eval(:,im,ispin,loc_i)=eval(:,noms,ispin,loc_i)
        !    !eval(ib,im,ispin,loc_i)=cmplx(real(eval(ib,noms,ispin,loc_i)),coeff_b/(om(im)+coeff_a))
        !    eval(ib,im,ispin,loc_i)=cmplx(Oinf_real+(Os_real-Oinf_real)*om(noms)/om(im),Oinf+(Os-Oinf)*om(noms)/om(im))
        !    !eval(ib,im,ispin,loc_i)=cmplx(real(eval(ib,noms,ispin,loc_i)),Oinf+(Os-Oinf)*om(noms)/om(im))
        !    !eval(ib,im,ispin,loc_i)=cmplx(Oinf,aimag(eval(ib,noms,ispin,loc_i)))
        !    !eval(ib,im,ispin,loc_i)=cmplx(Oinf+(Os-Oinf)*om(noms)**2/om(im)**2,aimag(eval(ib,noms,ispin,loc_i)))
        !    !write(*,*) ib,im,eval(ib,im,ispin,loc_i)
        !  enddo
        !enddo
        !    !coeff_real=(real(eval(ib,noms,ispin,loc_i))-eval0(ib,ispin,loc_i))*om(noms)**2
        !    !write(*,*) ib,eval(ib,noms,ispin,loc_i),eval(ib,noms+1,ispin,loc_i)
        !    !coeff_imag=aimag(eval(ib,noms,ispin,loc_i))*om(noms)
      enddo
    enddo
    !if (allocated(tran)) deallocate (tran)
    if (allocated(Hk)) deallocate (Hk)
    !if (allocated(eval0_loc)) deallocate (eval0_loc)
    !if (allocated(eval_loc)) deallocate (eval_loc)
    if (allocated(Hk_loc)) deallocate (Hk_loc)

!    if (on_root) then 
!     endif
!     call comms_bcast(mu,1) 
!    n_elec=48.0_dp
!    T=0.03_dp

    !write(*,*) my_node_id, 'hi3'
    low_mu=mu;high_mu=mu
    LOWB=.FALSE.;HIGHB=.FALSE.
    do l=1,mu_iter
      !call comms_barrier
      !write(*,*) l,'hi33'
      tot_n=0.0_dp
      loc_i=0
      mu=(low_mu+high_mu)/2
      do ik=displs(my_node_id)+1,displs(my_node_id)+numk_loc
        loc_i=loc_i+1
        do ispin=1,nspin
          do ib=1,num_wann
            tot_n=tot_n+2.0_dp*fermi((eval0(ib,ispin,loc_i)-mu)/T)*weight(ik)
            do im=1,nom
              Gs0=1.0_dp/(mu+cmplx_i*om(im)-eval0(ib,ispin,loc_i))
              Gs=1.0_dp/(mu+cmplx_i*om(im)-eval(ib,im,ispin,loc_i))
              tot_n=tot_n+4.0_dp*T*real(Gs-Gs0)*weight(ik)
            enddo
              !do im=noms+1,nom
              !  new_eval=cmplx(real(eval(ib,noms,ispin,loc_i)),coeff_b/(om(im)+coeff_a))
              !  Gs0=1.0_dp/(mu+cmplx_i*om(im)-eval0(ib,ispin,loc_i))
              !  Gs=1.0_dp/(mu+cmplx_i*om(im)-new_eval)
              !  tot_n=tot_n+4.0_dp*T*real(Gs-Gs0)*weight(ik)
              !enddo
            !at0=atan(om(noms)/(mu-eval0(ib,ispin,loc_i)))-atan(om(nom)/(mu-eval0(ib,ispin,loc_i)))
            !at=atan((om(noms)-aimag(eval(ib,noms,ispin,loc_i)))/(mu-real(eval(ib,noms,ispin,loc_i))))-atan((om(nom)-aimag(eval(ib,noms,ispin,loc_i)))/(mu-real(eval(ib,noms,ispin,loc_i))))
            !write(*,*) at, at0
            !tot_n=tot_n+2.0_dp*(at0-at)/pi
            !write(*,*) ik, ib, tot_n!*weight(ik)
          enddo
        enddo
        !tot_n=tot_n*weight(ik)
        !write(*,*) weight(ik)
        !STOP
      enddo
      call comms_allreduce(tot_n,1,'SUM')
!      if (on_root) write(*,*) mu, tot_n, n_elec
      if (abs(tot_n-n_elec)<1E-10) EXIT
      if (tot_n<n_elec) then
        if (HIGHB.eqv..FALSE.) then
          high_mu=mu+(n_elec-tot_n)*2/(0.1*n_elec)
        endif
        low_mu=mu
        LOWB=.TRUE.
      else
        if (LOWB.eqv..FALSE.) then
          low_mu=mu+(n_elec-tot_n)*2/(0.1*n_elec)
        endif
        high_mu=mu
        HIGHB=.TRUE.
      endif
      !if (on_root) write(*,*) mu, tot_n, n_elec
      !STOP
      if (l.eq.mu_iter) then
        write(*,*) 'Fail to find mu, increase mu_iter'
          !STOP
      endif
    enddo
    !call comms_barrier
    !write(*,*) 'hi4'

    if (allocated(eval0)) deallocate (eval0)
    if (allocated(eval)) deallocate (eval)

  end subroutine compute_DMFT_mu

  subroutine compute_dos()
    use constants
    use io
    use utility
    use generate_kpts, only: kpts, weight, num_new_kpts
    use generate_ham, only: tran
    use comms, only: on_root,comms_allreduce,my_node_id, num_nodes, comms_array_split,comms_reduce, comms_barrier
   
    implicit none

    character(len=10) :: write_format
    real(kind=dp) :: tot_n
    complex(kind=dp) :: Gs, Gs0
    integer :: i,j,l,ik,r,ib,im,ispin,loc_i,ierr
    integer :: nbmin,nbmax,num_band_max 
    integer :: numk_loc
    ! Needed to split an array on different nodes
    integer, dimension(0:num_nodes - 1) :: counts
    integer, dimension(0:num_nodes - 1) :: displs
    real(kind=dp) :: rdotk, fermi0

    real(kind=dp) :: Ekin, dEkin, mul_fac(ncor_orb)
    real(kind=dp), allocatable :: dm(:),mom(:)
    real(kind=dp), allocatable :: eval0(:)
    real(kind=dp), allocatable :: Ed(:)
    real(kind=dp), allocatable :: sym_Ed(:)
    complex(kind=dp), allocatable :: eval(:)
    complex(kind=dp), allocatable :: dHk(:,:)
    complex(kind=dp), allocatable :: Hk(:,:)
    complex(kind=dp), allocatable :: evec0(:,:)
    complex(kind=dp), allocatable :: evec(:,:)
    complex(kind=dp), allocatable :: Gloc_sum(:,:)
    complex(kind=dp), allocatable :: Delta(:,:)
    complex(kind=dp), allocatable :: Gloc(:,:,:)
    complex(kind=dp), allocatable :: sym_Gloc(:,:)
    complex(kind=dp) :: DMFT_UU,DMFT_U0 

    call comms_array_split(num_new_kpts, counts, displs)    
    numk_loc=counts(my_node_id)
    if (.not. allocated(Gloc)) then
      allocate (Gloc(num_wann,nom,nspin), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Gloc in compute_G_loc')
    endif
    if (.not. allocated(Hk)) then
      allocate (Hk(num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Hk in compute_G_loc')
    endif
    if (.not. allocated(evec)) then
      allocate (evec(num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating evec in compute_G_loc')
    endif
    if (.not. allocated(Ed)) then
      allocate (Ed(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating evec in compute_G_loc')
    endif

    loc_i=0; Ed=0.0_dp
    Gloc=cmplx_0;
    do ik=displs(my_node_id)+1,displs(my_node_id)+numk_loc
      loc_i=loc_i+1
      Hk=cmplx_0
      do r=1,nR
        rdotk=twopi*dot_product(kpts(:,ik),dfloat(tran(:,r)))
        Hk(:,:)=Hk(:,:)+HamR(r,:,:)*exp(cmplx_i*rdotk)
      enddo
      do ib=1,num_wann
        Ed(ib)=Ed(ib)+Hk(ib,ib)*weight(ik)*nspin
      enddo
      do ispin=1,nspin
        !!!!!!!! Compute the w=infinity value !!!!
        do im=1,nom
          evec=cmplx_0
          call zcopy(num_wann**2,(-1.0_dp)*Hk,1,evec,1)
          do j=1,num_wann
            evec(j,j)=evec(j,j)+mu+om(im)+cmplx_i*broaden
          enddo
          do i=1,n_atoms
            do j=1,n_orbs
              if (sym_idx(ispin,i,j)>0) evec((i-1)*n_orbs+j,(i-1)*n_orbs+j)=evec((i-1)*n_orbs+j,(i-1)*n_orbs+j)-Sigma(sym_idx(ispin,i,j),im)
            enddo
          enddo
          !do j=1,num_wann
          !  if (sym_idx(ispin,j)>0) evec(j,j)=evec(j,j)-Sigma(sym_idx(ispin,j),im)
          !!do j=1,ncor_orb
          !!  evec(j,j)=evec(j,j)-Sigma(j,im)+real(Sigma(j,1200))
          !enddo
          call MAT_INV(evec,num_wann)
          do i=1,num_wann
            Gloc(i,im,ispin)=Gloc(i,im,ispin)+evec(i,i)*weight(ik)*nspin
          enddo
        enddo
      enddo
    enddo
    if (allocated(HamR)) deallocate (HamR,stat=ierr)
    if (allocated(Hk)) deallocate (Hk,stat=ierr)
    if (allocated(evec)) deallocate (evec,stat=ierr)

    call comms_reduce(Gloc,num_wann,nom,nspin,'SUM')
    call comms_allreduce(Ed(1),num_wann,'SUM')

    if (on_root) then

      call write_complex_function('G_loc.out',Gloc(:,:,1),om)

      if (.not. allocated(sym_Gloc)) then
        allocate (sym_Gloc(ncor_orb,nom), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating sym_Gloc in compute_G_loc')
      endif
      sym_Gloc=cmplx_0
      mul_fac=0.0_dp
      do ispin=1,nspin
        !do ib=1,num_wann
        do i=1,n_atoms
          do j=1,n_orbs
            if (sym_idx(ispin,i,j)>0) then
              sym_Gloc(sym_idx(ispin,i,j),:)=sym_Gloc(sym_idx(ispin,i,j),:)+Gloc((i-1)*n_orbs+j,:,ispin)
              mul_fac(sym_idx(ispin,i,j))=mul_fac(sym_idx(ispin,i,j))+1.0_dp
            endif
          enddo
        enddo
      enddo
      do ib=1,ncor_orb
        sym_Gloc(ib,:)=sym_Gloc(ib,:)/mul_fac(ib)
      enddo
      
      !call write_complex_function('G_loc.out',sym_Gloc,om)

      if (.not. allocated(sym_Ed)) then
        allocate (sym_Ed(ncor_orb), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating sym_Ed in compute_G_loc')
      endif
      sym_Ed=0.0_dp
      mul_fac=0.0_dp
      !write(*,*) Ed
      !do ib=1,num_wann
      do ispin=1,nspin
        do i=1,n_atoms
          do j=1,n_orbs
            if (sym_idx(ispin,i,j)>0) then
              sym_Ed(sym_idx(ispin,i,j))=sym_Ed(sym_idx(ispin,i,j))+Ed((i-1)*n_orbs+j)!real(HamR((nR+1)/2,ib,ib))
              mul_fac(sym_idx(ispin,i,j))=mul_fac(sym_idx(ispin,i,j))+1.0_dp
            endif
          enddo
        enddo
      enddo
      do ib=1,ncor_orb
        sym_Ed(ib)=sym_Ed(ib)/mul_fac(ib)
      enddo
       
      if (.not. allocated(Delta)) then
        allocate (Delta(ncor_orb,nom), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating Delta in compute_G_loc')
      endif
      Delta=cmplx_0
      do ib=1,ncor_orb
        do im=1,nom
          Delta(ib,im)=mu+om(im)+cmplx_i*broaden-sym_Ed(ib)-Sigma(ib,im)-1.0_dp/sym_Gloc(ib,im)
        enddo
      enddo
      call write_complex_function('Delta.out',Delta,om)
      call write_float('Ed.out',sym_Ed)
      if (allocated(Delta)) deallocate (Delta, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating Delta in compute_G_loc')
      if (allocated(sym_Gloc)) deallocate (sym_Gloc, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating sym_Gloc in compute_G_loc')
      if (allocated(sym_Ed)) deallocate (sym_Ed, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating sym_Ed in compute_G_loc')
    endif

    if (allocated(UMatrix)) deallocate (UMatrix)
    if (allocated(kpts)) deallocate (kpts)
    if (allocated(weight)) deallocate (weight)
    if (allocated(Ed)) deallocate (Ed)
    if (allocated(Gloc)) deallocate (Gloc)

  end subroutine compute_dos

  subroutine compute_G_loc()
    use constants
    use io
    use utility
    use generate_kpts, only: kpts, weight, num_new_kpts
    use generate_ham, only: tran
    use comms, only: on_root,comms_allreduce,my_node_id, num_nodes, comms_array_split,comms_reduce, comms_barrier
   
    implicit none

    character(len=10) :: write_format
    real(kind=dp) :: tot_n
    complex(kind=dp) :: Gs, Gs0
    integer :: i,j,l,ik,r,ib,im,ispin,loc_i,ierr
    integer :: nbmin,nbmax,num_band_max 
    integer :: numk_loc
    ! Needed to split an array on different nodes
    integer, dimension(0:num_nodes - 1) :: counts
    integer, dimension(0:num_nodes - 1) :: displs
    real(kind=dp) :: rdotk, fermi0

    real(kind=dp) :: Ekin, dEkin, mul_fac(ncor_orb)
    real(kind=dp), allocatable :: dm(:),mom(:)
    real(kind=dp), allocatable :: eval0(:)
    real(kind=dp), allocatable :: Ed(:)
    real(kind=dp), allocatable :: sym_Ed(:)
    complex(kind=dp), allocatable :: eval(:)
    complex(kind=dp), allocatable :: dHk(:,:)
    complex(kind=dp), allocatable :: Hk(:,:)
    complex(kind=dp), allocatable :: evec0(:,:)
    complex(kind=dp), allocatable :: evec(:,:)
    complex(kind=dp), allocatable :: Gloc_sum(:,:)
    complex(kind=dp), allocatable :: Delta(:,:)
    complex(kind=dp), allocatable :: Gloc(:,:,:)
    complex(kind=dp), allocatable :: sym_Gloc(:,:)
    complex(kind=dp) :: DMFT_UU,DMFT_U0 

    !write(*,*) num_new_kpts
    call comms_array_split(num_new_kpts, counts, displs)    
    !write(*,*) my_node_id, counts(my_node_id), displs(my_node_id)
    numk_loc=counts(my_node_id)
!    mu_iter=100
    if (.not. allocated(Gloc)) then
      allocate (Gloc(num_wann,nom,nspin), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Gloc in compute_G_loc')
    endif
    if (.not. allocated(Gloc_sum)) then
      allocate (Gloc_sum(num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Gloc_sum in compute_G_loc')
    endif
    if (.not. allocated(mom)) then
      allocate (mom(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating mom in compute_G_loc')
    endif
    if (.not. allocated(dm)) then
      allocate (dm(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating dm in compute_G_loc')
    endif
    !write(*,*) my_node_id, allocated(eval0)
    if (.not. allocated(eval0)) then
      allocate (eval0(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval0 in compute_G_loc')
    endif
    !if (.not. allocated(eval)) then
    !  allocate (eval(num_wann), stat=ierr)
    !  if (ierr /= 0) call io_error('Error allocating eval in compute_G_loc')
    !endif
    if (.not. allocated(Hk)) then
      allocate (Hk(num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Hk in compute_G_loc')
    endif
    if ((.not. allocated(dHk)).and. (lforce.eq..true.)) then
      allocate (dHk(num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating dHk in compute_G_loc')
    endif

    if (.not. allocated(evec0)) then
      allocate (evec0(num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating evec0 in compute_G_loc')
    endif
    if (.not. allocated(evec)) then
      allocate (evec(num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating evec in compute_G_loc')
    endif
    if (.not. allocated(Ed)) then
      allocate (Ed(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating evec in compute_G_loc')
    endif
    !if (.not. allocated(evl)) then
    !  allocate (evl(num_wann,num_wann), stat=ierr)
    !  if (ierr /= 0) call io_error('Error allocating evl in compute_G_loc')
    !endif
    !if (.not. allocated(evr)) then
    !  allocate (evr(num_wann,num_wann), stat=ierr)
    !  if (ierr /= 0) call io_error('Error allocating evr in compute_G_loc')
    !endif

    !write(*,*) my_node_id, 'hi0'
    loc_i=0; Ed=0.0_dp
    Ekin=0.0_dp;dm=0.0_dp;mom=0.0_dp;Gloc=cmplx_0;
    dEkin=0.0_dp
    do ik=displs(my_node_id)+1,displs(my_node_id)+numk_loc
      loc_i=loc_i+1
      Hk=cmplx_0
      if (lforce.eq..true.) dHk=cmplx_0
      !write(*,*) kpts(1,1),tran(1,1),HamR(1,1,1)
      do r=1,nR
        rdotk=twopi*dot_product(kpts(:,ik),dfloat(tran(:,r)))
        Hk(:,:)=Hk(:,:)+HamR(r,:,:)*exp(cmplx_i*rdotk)
        if (lforce.eq..true.) dHk(:,:)=dHk(:,:)+dHamR(r,:,:)*exp(cmplx_i*rdotk)
      enddo
      do ib=1,num_wann
        Ed(ib)=Ed(ib)+Hk(ib,ib)*weight(ik)*nspin
      enddo
      do ispin=1,nspin
        !!!!!!!! Compute the w=infinity value !!!!
        eval0=0.0_dp
        evec0=cmplx_0
        call zcopy(num_wann**2,Hk,1,evec0,1)
        Gloc_sum=cmplx_0
        do i=1,n_atoms
          do j=1,n_orbs
            if (sym_idx(ispin,i,j)>0) evec0((i-1)*n_orbs+j,(i-1)*n_orbs+j)=evec0((i-1)*n_orbs+j,(i-1)*n_orbs+j)+Sigoo(sym_idx(ispin,i,j))
          enddo
        enddo
        !do i=1,num_wann
        !  if (sym_idx(ispin,i)>0) evec0(i,i)=evec0(i,i)+Sigoo(sym_idx(ispin,i))
        !!do i=1,ncor_orb
        !!  evec0(i,i)=evec0(i,i)+Sigoo(i)
        !enddo
        CALL EIGENH(evec0,num_wann,eval0)
        do im=1,nom
          evec=cmplx_0
          !evec=Hk*(-1.0_dp)
          call zcopy(num_wann**2,(-1.0_dp)*Hk,1,evec,1)
          do j=1,num_wann
            evec(j,j)=evec(j,j)+mu+cmplx_i*om(im)
          enddo
          do i=1,n_atoms
            do j=1,n_orbs
              if (sym_idx(ispin,i,j)>0) evec((i-1)*n_orbs+j,(i-1)*n_orbs+j)=evec((i-1)*n_orbs+j,(i-1)*n_orbs+j)-Sigma(sym_idx(ispin,i,j),im)
            enddo
          enddo
          !do j=1,num_wann
          !  if (sym_idx(ispin,j)>0) evec(j,j)=evec(j,j)-Sigma(sym_idx(ispin,j),im)
          !!do j=1,ncor_orb
          !!  evec(j,j)=evec(j,j)-Sigma(j,im)+real(Sigma(j,1200))
          !enddo
          call MAT_INV(evec,num_wann)
          do i=1,num_wann
            Gloc(i,im,ispin)=Gloc(i,im,ispin)+evec(i,i)*weight(ik)*nspin
          enddo
          Gloc_sum(:,:)=Gloc_sum(:,:)+T*(evec(:,:)+transpose(conjg(evec(:,:))))*weight(ik)
          !eval=cmplx_0; evec=cmplx_0; evl=cmplx_0; evr=cmplx_0; 
          !evec=Hk*1.0_dp
          !do j=1,ncor_orb
          !  evec(j,j)=evec(j,j)+Sigma(j,im)-real(Sigma(j,1200))
          !enddo
          !CALL EIGEN(evec,num_wann,eval,evl,evr)
          !do ib=1,num_wann
          !  Gs0=1.0_dp/(mu+cmplx_i*om(im)-eval0(ib))
          !  !Gs=1.0_dp/(mu+cmplx_i*om(im)-eval(ib))
!         !       Ekin=Ekin-2*T*(log(abs(Gs))-log(abs(Gs0)))
!         !                 *kweight(nkp)*sweight
          !  do i=1,num_wann
          !    do j=1,num_wann
          !      !DMFT_UU=evr(i,ib)*evl(ib,j)*Gs*weight(ik)
          !      !if (i.eq.j) Gloc(i,im,ispin)=Gloc(i,im,ispin)+DMFT_UU*nspin
          !      DMFT_U0=evec0(i,ib)*conjg(evec0(j,ib))*Gs0*weight(ik)
          !      Gloc_sum(i,j)=Gloc_sum(i,j)-T*DMFT_U0
          !      Gloc_sum(j,i)=Gloc_sum(j,i)-T*conjg(DMFT_U0)
          !    enddo
          !  enddo
          !enddo
          evec=cmplx_0
          call zcopy(num_wann**2,(-1.0_dp)*Hk,1,evec,1)
          !evec=Hk*(-1.0_dp)
          do j=1,num_wann
            evec(j,j)=evec(j,j)+mu+cmplx_i*om(im)
          enddo
          do i=1,n_atoms
            do j=1,n_orbs
              if (sym_idx(ispin,i,j)>0) evec((i-1)*n_orbs+j,(i-1)*n_orbs+j)=evec((i-1)*n_orbs+j,(i-1)*n_orbs+j)-Sigoo(sym_idx(ispin,i,j))
            enddo
          enddo
         ! do j=1,num_wann
         !   if (sym_idx(ispin,j)>0) evec(j,j)=evec(j,j)-Sigoo(sym_idx(ispin,j))
         ! enddo
          call MAT_INV(evec,num_wann)
          Gloc_sum(:,:)=Gloc_sum(:,:)-T*(evec(:,:)+transpose(conjg(evec(:,:))))*weight(ik)
        enddo
        !do im=noms+1,nom
        !enddo

        do ib=1,num_wann
          !Gsr=fermi((eval0(ib)-mu)/T)
!              Ekin=Ekin-T*logfermi((eval0(nb)-mu)/T)
!                    *kweight(nkp)*sweight
          !write(*,*) mu,eval0(ib),T,fermi((eval0(ib)-mu)/T),weight(ik)
          fermi0=fermi((eval0(ib)-mu)/T)
          do i=1,num_wann
            do j=1,num_wann
              Gloc_sum(i,j)=Gloc_sum(i,j)+fermi0* &
              evec0(i,ib)*conjg(evec0(j,ib))*weight(ik)
            enddo
          enddo
        enddo



        !write(*,*) Gloc_sum(1,1)
        do i=1,num_wann
          dm(i)=dm(i)+real(Gloc_sum(i,i))*2
          if (ispin.eq.1 .AND. nspin.eq.2) then
            mom(i)=mom(i)+real(Gloc_sum(i,i))*2
          endif
          if (ispin.eq.2) mom(i)=mom(i)-real(Gloc_sum(i,i))*2
          do j=1,num_wann
            Ekin=Ekin+real(Hk(i,j)*Gloc_sum(j,i))*2
            if (lforce.eq..true.) dEkin=dEkin+real(dHk(i,j)*Gloc_sum(j,i))*2
          enddo
        enddo


      enddo
    enddo
!    write(*,*) lforce, dEkin
    if (allocated(HamR)) deallocate (HamR,stat=ierr)
    if (allocated(dHamR)) deallocate (dHamR,stat=ierr)
    if (allocated(Hk)) deallocate (Hk,stat=ierr)
    if (allocated(dHk)) deallocate (dHk,stat=ierr)
    if (allocated(eval0)) deallocate (eval0,stat=ierr)
    if (allocated(evec0)) deallocate (evec0,stat=ierr)
    if (allocated(evec)) deallocate (evec,stat=ierr)
    if (allocated(Gloc_sum)) deallocate (Gloc_sum,stat=ierr)

    call comms_reduce(Gloc,num_wann,nom,nspin,'SUM')
    call comms_allreduce(Ekin,1,'SUM')
    call comms_allreduce(dEkin,1,'SUM')
    call comms_allreduce(Ed(1),num_wann,'SUM')
    call comms_allreduce(dm(1),num_wann,'SUM')
    call comms_allreduce(mom(1),num_wann,'SUM')
    tot_n=0.0_dp
    do i=1,num_wann
      tot_n=tot_n+dm(i)
    enddo
!    call comms_barrier

    !if (on_root) write(*,*) tot_n,Gloc(1,1,1), Ekin, Ed, dm, mom

    if (on_root) then

      if (.not. allocated(sym_Gloc)) then
        allocate (sym_Gloc(ncor_orb,nom), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating sym_Gloc in compute_G_loc')
      endif
      sym_Gloc=cmplx_0
      mul_fac=0.0_dp
      do ispin=1,nspin
        !do ib=1,num_wann
        do i=1,n_atoms
          do j=1,n_orbs
            if (sym_idx(ispin,i,j)>0) then
              sym_Gloc(sym_idx(ispin,i,j),:)=sym_Gloc(sym_idx(ispin,i,j),:)+Gloc((i-1)*n_orbs+j,:,ispin)
              mul_fac(sym_idx(ispin,i,j))=mul_fac(sym_idx(ispin,i,j))+1.0_dp
            endif
          enddo
        enddo
      enddo
      do ib=1,ncor_orb
        sym_Gloc(ib,:)=sym_Gloc(ib,:)/mul_fac(ib)
      enddo
      
!      write(*,*) 'hi1'
      call write_complex_function('G_loc.out',sym_Gloc,om)

      if (.not. allocated(sym_Ed)) then
        allocate (sym_Ed(ncor_orb), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating sym_Ed in compute_G_loc')
      endif
      sym_Ed=0.0_dp
      mul_fac=0.0_dp
      !write(*,*) Ed
      !do ib=1,num_wann
      do ispin=1,nspin
        do i=1,n_atoms
          do j=1,n_orbs
            if (sym_idx(ispin,i,j)>0) then
              sym_Ed(sym_idx(ispin,i,j))=sym_Ed(sym_idx(ispin,i,j))+Ed((i-1)*n_orbs+j)!real(HamR((nR+1)/2,ib,ib))
              mul_fac(sym_idx(ispin,i,j))=mul_fac(sym_idx(ispin,i,j))+1.0_dp
            endif
          enddo
        enddo
      enddo
      do ib=1,ncor_orb
        sym_Ed(ib)=sym_Ed(ib)/mul_fac(ib)
      enddo
      !write(*,*) real(HamR((nR+1)/2,1,1)),real(HamR((nR+1)/2,2,2)),real(HamR((nR+1)/2,3,3)),real(HamR((nR+1)/2,4,4)),real(HamR((nR+1)/2,5,5)),real(HamR((nR+1)/2,6,6))
      !write(*,*) mul_fac
       
      if (.not. allocated(Delta)) then
        allocate (Delta(ncor_orb,nom), stat=ierr)
        if (ierr /= 0) call io_error('Error allocating Delta in compute_G_loc')
      endif
      Delta=cmplx_0
      do ib=1,ncor_orb
        do im=1,nom
          Delta(ib,im)=mu+cmplx_i*om(im)-sym_Ed(ib)-Sigma(ib,im)-1.0_dp/sym_Gloc(ib,im)
        enddo
      enddo
!      write(*,*) Delta(1,1) 
!      write(*,*) 'hi2'
      call write_complex_function('Delta.out',Delta,om)
      call write_float('Ed.out',sym_Ed)
!      call write_float('DMFT_mu.out',mu)

!      write(*,*) 'hi3'
!      write(*,*) sym_Ed 
!      print *, "Hello world", my_node_id, num_nodes
!!      OPEN(UNIT=99,FILE='Ed.out',FORM='FORMATTED')
!! 110  FORMAT('(',I2,'F20.15)')
!!      WRITE(write_format,110) ncor_orb
!!      write(*,*) write_format
!!      WRITE(99,write_format) (sym_Ed(ib), ib=1,ncor_orb)
!!      CLOSE(99)
!      write(*,*) 'hi4'
      OPEN(UNIT=99,FILE='DMFT_mu.out',FORM='FORMATTED')
      WRITE(99,'(F20.15)') mu
      CLOSE(99)
!      write(*,*) 'hi5'
      OPEN(UNIT=99,FILE='INFO_KSUM',STATUS='old',FORM='FORMATTED',ACCESS='APPEND')
      WRITE(99,'(11F12.6)') mu, tot_n, SUM(dm(1:n_orbs)), &
      !SUM(dm(ncor_orb-norb+1:ncor_orb)), tot_Ekin,!+mu*total_n,
      SUM(dm(n_orbs*(n_atoms-1)+1:n_orbs*n_atoms)), Ekin, Sigoo(1), dEkin
      CLOSE(99)
!      write(*,*) 'hi6'
      OPEN(UNIT=99,FILE='INFO_DM',STATUS='old',FORM='FORMATTED',ACCESS='APPEND')
100   FORMAT('(',I2,'F20.15)')
      WRITE(write_format,100) n_atoms*n_orbs
      WRITE(99,write_format,advance='no') (dm(ib), ib=1,n_atoms*n_orbs)
      WRITE(99,write_format) (mom(ib), ib=1,n_atoms*n_orbs)
      CLOSE(99)

      !write(*,*) allocated(Delta),allocated(sym_Gloc),allocated(sym_Ed)
      if (allocated(Delta)) deallocate (Delta, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating Delta in compute_G_loc')
      if (allocated(sym_Gloc)) deallocate (sym_Gloc, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating sym_Gloc in compute_G_loc')
      if (allocated(sym_Ed)) deallocate (sym_Ed, stat=ierr)
      if (ierr /= 0) call io_error('Error deallocating sym_Ed in compute_G_loc')
      !write(*,*) 'hi2'
    endif
    !call comms_barrier


    if (allocated(UMatrix)) deallocate (UMatrix)
    if (allocated(kpts)) deallocate (kpts)
    if (allocated(weight)) deallocate (weight)
    if (allocated(Ed)) deallocate (Ed)
    if (allocated(Gloc)) deallocate (Gloc)
    if (allocated(dm)) deallocate (dm)
    if (allocated(mom)) deallocate (mom)

  end subroutine compute_G_loc

  subroutine compute_DMFT_DM()
    use constants
    use io
    use utility
    use generate_kpts
    use generate_ham, only: tran
    use comms, only: on_root,comms_allreduce,my_node_id, num_nodes, comms_array_split,comms_reduce, comms_barrier
   
    implicit none

    character(len=10) :: write_format
    real(kind=dp) :: tot_n
    complex(kind=dp) :: Gs, Gs0
    integer :: i,j,l,ik,r,ib,im,ispin,loc_i,ierr
    integer :: nbmin,nbmax,num_band_max 
    integer :: numk_loc,nfine_tot
    !integer,intent(in) :: n_kpts_loc,n_wann
    ! Needed to split an array on different nodes
    integer, dimension(0:num_nodes - 1) :: counts
    integer, dimension(0:num_nodes - 1) :: displs
    real(kind=dp) :: rdotk, fermi0, tot

    real(kind=dp) :: Ekin, Nd, sweight
    real(kind=dp), allocatable :: eval0(:)
    real(kind=dp), allocatable :: Ed(:)
    complex(kind=dp), allocatable :: eval(:)
    complex(kind=dp), allocatable :: Hk(:,:)
    complex(kind=dp), allocatable :: evec0(:,:)
    complex(kind=dp), allocatable :: evec(:,:)
    complex(kind=dp), allocatable :: Gloc_sum(:,:,:)
    complex(kind=dp), allocatable :: Delta(:,:)
    complex(kind=dp), allocatable :: Gloc(:,:,:)
    complex(kind=dp), allocatable :: sym_Gloc(:,:)
    complex(kind=dp) :: DMFT_UU,DMFT_U0 

    call comms_array_split(num_new_kpts, counts, displs)    
    numk_loc=counts(my_node_id)
    if (.not. allocated(Gloc_sum)) then
      allocate (Gloc_sum(num_wann,num_wann,num_new_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Gloc_sum in compute_G_loc')
    endif
    if (.not. allocated(eval0)) then
      allocate (eval0(num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating eval0 in compute_G_loc')
    endif
    if (.not. allocated(Hk)) then
      allocate (Hk(num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating Hk in compute_G_loc')
    endif

    if (.not. allocated(evec0)) then
      allocate (evec0(num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating evec0 in compute_G_loc')
    endif
    if (.not. allocated(evec)) then
      allocate (evec(num_wann,num_wann), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating evec in compute_G_loc')
    endif

    sweight=1.0_dp/dfloat(nspin)
    loc_i=0
    tot_n=0.0_dp;Nd=0.0_dp;Ekin=0.0_dp
    Gloc_sum=cmplx_0
    do ik=displs(my_node_id)+1,displs(my_node_id)+numk_loc
      loc_i=loc_i+1
      Hk=cmplx_0
      do r=1,nR
        rdotk=twopi*dot_product(kpts(:,ik),dfloat(tran(:,r)))
        Hk(:,:)=Hk(:,:)+HamR(r,:,:)*exp(cmplx_i*rdotk)
      enddo
      do ispin=1,nspin
        !!!!!!!! Compute the w=infinity value !!!!
        eval0=0.0_dp
        evec0=cmplx_0
        call zcopy(num_wann**2,Hk,1,evec0,1)
        !evec0=Hk*1.0_dp
        do i=1,n_atoms
          do j=1,n_orbs
            if (sym_idx(ispin,i,j)>0) evec0((i-1)*n_orbs+j,(i-1)*n_orbs+j)=evec0((i-1)*n_orbs+j,(i-1)*n_orbs+j)+Sigoo(sym_idx(ispin,i,j))
          enddo
        enddo
        !do i=1,num_wann
        !  if (sym_idx(ispin,i)>0) evec0(i,i)=evec0(i,i)+Sigoo(sym_idx(ispin,i))
        !!do i=1,ncor_orb
        !!  evec0(i,i)=evec0(i,i)+Sigoo(i)
        !enddo
        CALL EIGENH(evec0,num_wann,eval0)
        do im=1,nom
          evec=cmplx_0
          !evec=Hk*(-1.0_dp)
          call zcopy(num_wann**2,(-1.0_dp)*Hk,1,evec,1)
          do j=1,num_wann
            evec(j,j)=evec(j,j)+mu+cmplx_i*om(im)
          enddo
          do i=1,n_atoms
            do j=1,n_orbs
              if (sym_idx(ispin,i,j)>0) evec((i-1)*n_orbs+j,(i-1)*n_orbs+j)=evec((i-1)*n_orbs+j,(i-1)*n_orbs+j)-Sigma(sym_idx(ispin,i,j),im)
            enddo
          enddo
          !do j=1,num_wann
          !  if (sym_idx(ispin,j)>0) evec(j,j)=evec(j,j)-Sigma(sym_idx(ispin,j),im)
          !enddo
          call MAT_INV(evec,num_wann)
          Gloc_sum(:,:,ik)=Gloc_sum(:,:,ik)+T*(evec(:,:)+transpose(conjg(evec(:,:))))*sweight

          evec=cmplx_0
          !evec=Hk*(-1.0_dp)
          call zcopy(num_wann**2,(-1.0_dp)*Hk,1,evec,1)
          do j=1,num_wann
            evec(j,j)=evec(j,j)+mu+cmplx_i*om(im)
          enddo
          do i=1,n_atoms
            do j=1,n_orbs
              if (sym_idx(ispin,i,j)>0) evec((i-1)*n_orbs+j,(i-1)*n_orbs+j)=evec((i-1)*n_orbs+j,(i-1)*n_orbs+j)-Sigoo(sym_idx(ispin,i,j))
            enddo
          enddo
          !do j=1,num_wann
          !  if (sym_idx(ispin,j)>0) evec(j,j)=evec(j,j)-Sigoo(sym_idx(ispin,j))
          !enddo
          call MAT_INV(evec,num_wann)
          Gloc_sum(:,:,ik)=Gloc_sum(:,:,ik)-T*(evec(:,:)+transpose(conjg(evec(:,:))))*sweight
        enddo
        !do im=noms+1,nom
        !enddo

        do ib=1,num_wann
          fermi0=fermi((eval0(ib)-mu)/T)
          do i=1,num_wann
            do j=1,num_wann
              Gloc_sum(i,j,ik)=Gloc_sum(i,j,ik)+fermi0* &
              evec0(i,ib)*conjg(evec0(j,ib))*sweight
            enddo
          enddo
        enddo
        do i=1,num_wann
          tot_n=tot_n+real(Gloc_sum(i,i,ik))*2*weight(ik)
          do j=1,num_wann
            Ekin=Ekin+real(Hk(i,j)*Gloc_sum(j,i,ik))*2*weight(ik)
          enddo
        enddo
        do i=1,5
          Nd=Nd+real(Gloc_sum(i,i,ik))*2*weight(ik)
        enddo
      enddo
    enddo
    if (allocated(HamR)) deallocate (HamR,stat=ierr)
    if (allocated(Hk)) deallocate (Hk,stat=ierr)
    if (allocated(eval0)) deallocate (eval0,stat=ierr)
    if (allocated(evec0)) deallocate (evec0,stat=ierr)
    if (allocated(kpts)) deallocate (kpts)
    if (allocated(weight)) deallocate (weight)

    call comms_allreduce(Gloc_sum,num_wann,num_wann,num_new_kpts,'SUM')
    call comms_allreduce(Ekin,1,'SUM')
    call comms_allreduce(tot_n,1,'SUM')
    call comms_allreduce(Nd,1,'SUM')


    call comms_array_split(n_kpts, counts, displs)    
    numk_loc=counts(my_node_id)

    if (.not. allocated(DMFT_eval_loc)) then
      allocate (DMFT_eval_loc(num_wann,n_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating DMFT_eval_loc in compute_G_loc')
    endif
    if (.not. allocated(DMFT_evec_loc)) then
      allocate (DMFT_evec_loc(num_wann,num_wann,n_kpts), stat=ierr)
      if (ierr /= 0) call io_error('Error allocating DMFT_evec_loc in compute_G_loc')
    endif

    nfine_tot=nfine(1)*nfine(2)*nfine(3)
    DMFT_eval_loc=0.0_dp; DMFT_evec_loc=cmplx_0
    do ik=displs(my_node_id)+1,displs(my_node_id)+numk_loc
      evec=cmplx_0
      do i=1,nfine_tot
        evec(:,:)=evec(:,:)+Gloc_sum(:,:,i+(ik-1)*nfine_tot)/dfloat(nfine_tot)
      enddo
      call EIGENH(evec,num_wann,DMFT_eval_loc(:,ik))
      !DMFT_evec_loc(:,:,ik)=evec(:,:)*1.0_dp
      !do i=1,num_wann
      !  tot=0.0_dp
      !  do j=1,num_wann
      !    tot=tot+UMatrix(band_win(1,ik)+i-1,j,ik)*conjg(UMatrix(band_win(1,ik)+i-1,j,ik))
      !  enddo
      !  write(*,*) band_win(1,ik),num_wann,i,tot
      !enddo
      DMFT_evec_loc(:,:,ik)=MATMUL(UMatrix(1:num_wann,:,ik),evec(:,:))
    enddo
    call comms_allreduce(DMFT_evec_loc,num_wann,num_wann,n_kpts,'SUM')
    call comms_allreduce(DMFT_eval_loc,num_wann,n_kpts,'SUM')
    
    if (allocated(evec)) deallocate (evec,stat=ierr)
    if (allocated(Gloc_sum)) deallocate (Gloc_sum,stat=ierr)
    if (allocated(UMatrix)) deallocate (UMatrix)

    !write(*,*) 'a', my_node_id, DMFT_eval_loc(:,1)
    if (on_root) then

      !write(*,*) tot_n, Nd, Ekin
      OPEN(UNIT=99,FILE='DMFT_mu.out',FORM='FORMATTED')
      WRITE(99,'(F20.15)') mu
      CLOSE(99)

      OPEN(UNIT=99,FILE='INFO_DFT_loop',STATUS='old',FORM='FORMATTED',ACCESS='APPEND')
      WRITE(99,'(4F12.6)') mu, tot_n, Nd, Ekin 
      CLOSE(99)

    endif
  end subroutine compute_DMFT_DM

!  subroutine compute_Delta()
!    use constants
!    use io
!    use utility
!    use generate_kpts, only: kpts, weight, num_new_kpts
!    use generate_ham, only: tran
!    use comms, only: on_root,comms_allreduce,my_node_id, num_nodes, comms_array_split,comms_reduce
!   
!    implicit none
!
!    logical :: iffile
!    logical :: LOWB,HIGHB
!    real(kind=dp) :: low_mu, high_mu, tot_n
!    complex(kind=dp) :: Gs, Gs0
!    integer :: i,j,l,ik,r,ib,im,ispin,loc_i,ierr
!    integer :: nbmin,nbmax,num_band_max 
!    integer :: numk_loc,mu_iter 
!    ! Needed to split an array on different nodes
!    integer, dimension(0:num_nodes - 1) :: counts
!    integer, dimension(0:num_nodes - 1) :: displs
!    real(kind=dp) :: rdotk
!
!    real(kind=dp) :: Ekin
!    real(kind=dp), allocatable :: dm(:),mom(:)
!    real(kind=dp), allocatable :: eval0(:)
!    complex(kind=dp), allocatable :: eval(:)
!    complex(kind=dp), allocatable :: Hk(:,:)
!    complex(kind=dp), allocatable :: evec0(:,:)
!    complex(kind=dp), allocatable :: evec(:,:)
!    complex(kind=dp), allocatable :: evl(:,:)
!    complex(kind=dp), allocatable :: evr(:,:)
!    complex(kind=dp), allocatable :: Gloc_sum(:,:)
!    complex(kind=dp) :: DMFT_UU,DMFT_U0 
!
!    !write(*,*) num_new_kpts
!    call comms_array_split(num_new_kpts, counts, displs)    
!    !write(*,*) my_node_id, counts(my_node_id), displs(my_node_id)
!    numk_loc=counts(my_node_id)
!    mu_iter=100
!    if (.not. allocated(Gloc)) then
!      allocate (Gloc(num_wann,nom,nspin), stat=ierr)
!      if (ierr /= 0) call io_error('Error allocating Gloc in compute_G_loc')
!    endif
!    if (.not. allocated(Gloc_sum)) then
!      allocate (Gloc_sum(num_wann,num_wann), stat=ierr)
!      if (ierr /= 0) call io_error('Error allocating Gloc_sum in compute_G_loc')
!    endif

end module dmft_ksum 
 
