!------------------------------------------------------------!
!                                                            !
!                       DMFTwDFT                             !
!                                                            !
!          DMFT self-consistent calculation using            !
!                  Wannier functions                         !
!                                                            !
! Please cite XXXXX                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Compute_DMFT(n_kpts_loc,n_wann,kpt_dft,wght_dft,band_win_loc,DMFT_eval,DMFT_evec)
  !!! main program
  USE constants
  USE comms, only: my_node_id,num_nodes,comms_setup_vars
  USE read_inputs
  USE generate_ham
  USE generate_kpts
  USE dmft_ksum

  implicit none

  integer, intent(in) :: n_kpts_loc,n_wann!,mpi_comm_world_loc,mpi_sum_loc
  !integer :: my_node_id,ierr,num_nodes,b!,mpi_integer_loc
  real(kind=dp), intent(in) :: kpt_dft(3,n_kpts_loc)
  real(kind=dp), intent(in) :: wght_dft(n_kpts_loc)
  integer,intent(out) :: band_win_loc(2,n_kpts_loc)
  real(kind=dp), intent(out) :: DMFT_eval(n_wann,n_kpts_loc)
  complex(kind=dp), intent(out) :: DMFT_evec(n_wann,n_wann,n_kpts_loc)
  integer :: i,j
  real(kind=dp) :: tot
!  integer :: id,nodes

!  write(*,*) kpt_dft
!!
!!  mpi_comm_world=mpi_comm_world_loc
!  call mpi_comm_rank(mpi_comm_world_loc, my_node_id, ierr)
!  call mpi_comm_size(mpi_comm_world_loc, num_nodes, ierr)
!  call mpi_allreduce(my_node_id,b,1,mpi_integer_loc,mpi_sum_loc,mpi_comm_world_loc,ierr)
  call comms_setup_vars
!  print *, "Hello world", my_node_id, num_nodes
  !call comms_end
!  call test(mpi_comm_world)
  call Read_wan_chk()
  call Read_wan_eig()
  call Read_sig_inp()
  call Read_dmft_params()
  call generate_hamr()
  !call generate_uniform_kmesh()
  call generate_dense_kmesh(n_kpts_loc,kpt_dft,wght_dft)
  call compute_DMFT_mu()
  !call compute_G_loc()

  call compute_DMFT_DM()
  !write(*,*) band_win(1,1)
  band_win_loc=0
  do i=1,n_kpts_loc
    band_win_loc(:,i)=band_win(:,i)
  enddo
  !write(*,*) band_win_loc(1,1)
  DMFT_eval=0.0_dp;DMFT_evec=cmplx_0
  call dcopy(n_wann*n_kpts,DMFT_eval_loc,1,DMFT_eval,1)
  call zcopy(n_wann**2*n_kpts,DMFT_evec_loc,1,DMFT_evec,1)

end subroutine Compute_DMFT

subroutine Compute_DMFT_from_amn(n_kpts_loc,n_wann,kpt_dft,wght_dft,band_win_loc,DMFT_eval,DMFT_evec)
  !!! main program
  USE constants
  USE comms, only: my_node_id,num_nodes,comms_setup_vars
  USE read_inputs
  USE generate_ham
  USE generate_kpts
  USE dmft_ksum

  implicit none

  integer, intent(in) :: n_kpts_loc,n_wann!,mpi_comm_world_loc,mpi_sum_loc
  !integer :: my_node_id,ierr,num_nodes,b!,mpi_integer_loc
  real(kind=dp), intent(in) :: kpt_dft(3,n_kpts_loc)
  real(kind=dp), intent(in) :: wght_dft(n_kpts_loc)
  integer,intent(out) :: band_win_loc(2,n_kpts_loc)
  real(kind=dp), intent(out) :: DMFT_eval(n_wann,n_kpts_loc)
  complex(kind=dp), intent(out) :: DMFT_evec(n_wann,n_wann,n_kpts_loc)
  integer :: i
!  integer :: id,nodes

!  write(*,*) kpt_dft
!!
!!  mpi_comm_world=mpi_comm_world_loc
!  call mpi_comm_rank(mpi_comm_world_loc, my_node_id, ierr)
!  call mpi_comm_size(mpi_comm_world_loc, num_nodes, ierr)
!  call mpi_allreduce(my_node_id,b,1,mpi_integer_loc,mpi_sum_loc,mpi_comm_world_loc,ierr)
  call comms_setup_vars
!  print *, "Hello world", my_node_id, num_nodes
  !call comms_end
!  call test(mpi_comm_world)
  call Read_wan_chk()
  call Read_wan_eig()
  call Read_wan_amn()
  call Compute_UNI_from_amn()
  call Read_sig_inp()
  call Read_dmft_params()
  call generate_hamr()
  !STOP
  !call generate_uniform_kmesh()
  call generate_dense_kmesh(n_kpts_loc,kpt_dft,wght_dft)
  call compute_DMFT_mu()
  !call compute_G_loc()

  call compute_DMFT_DM()
  band_win_loc=0
  do i=1,n_kpts_loc
    band_win_loc(:,i)=band_win(:,i)
  enddo
  DMFT_eval=0.0_dp;DMFT_evec=cmplx_0
  call dcopy(n_wann*n_kpts,DMFT_eval_loc,1,DMFT_eval,1)
  call zcopy(n_wann**2*n_kpts,DMFT_evec_loc,1,DMFT_evec,1)

end subroutine Compute_DMFT_from_amn
