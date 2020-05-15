
module dmft_ksum 

  use constants, only: dp
  use io, only: stdout, maxlen
  use read_inputs 

  implicit none

  real(kind=dp), allocatable, save :: DMFT_eval_loc(:,:)
  complex(kind=dp), allocatable, save :: DMFT_evec_loc(:,:,:)

contains
 
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
      !if (on_root) write(*,*) mu, tot_n, n_elec
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
      do i=1,n_atoms
        do j=1,n_orbs
          if (sym_idx(1,i,j)>0) then
            sym_Ed(sym_idx(1,i,j))=sym_Ed(sym_idx(1,i,j))+Ed((i-1)*n_orbs+j)!real(HamR((nR+1)/2,ib,ib))
            mul_fac(sym_idx(1,i,j))=mul_fac(sym_idx(1,i,j))+1.0_dp
          endif
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
    call comms_barrier

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
      do i=1,n_atoms
        do j=1,n_orbs
          if (sym_idx(1,i,j)>0) then
            sym_Ed(sym_idx(1,i,j))=sym_Ed(sym_idx(1,i,j))+Ed((i-1)*n_orbs+j)!real(HamR((nR+1)/2,ib,ib))
            mul_fac(sym_idx(1,i,j))=mul_fac(sym_idx(1,i,j))+1.0_dp
          endif
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
 
