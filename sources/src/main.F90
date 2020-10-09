!------------------------------------------------------------!
!                                                            !
!                       DMFTwDFT                             !
!                                                            !
!          DMFT self-consistent calculation using            !
!                  Wannier functions                         !
!                                                            !
! Please cite XXXXX                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program DMFTwDFT 
  !!! main program
  USE constants 
  USE comms
  USE read_inputs
  USE generate_ham
  USE generate_kpts
  USE dmft_ksum 

  implicit none

  !integer :: i

  call comms_setup
  !print *, "Hello world", my_node_id, num_nodes
  call Read_wan_chk() 
  call Read_wan_eig() 
  !call Read_wan_amn() 
  !call Compute_UNI_from_amn()
  if (on_root) then
    call Read_wan_amn() 
    !call Check_Unitarity() 
    call Print_overlap() 
  endif
  !STOP
  call Read_sig_inp()
  call Read_dmft_params()
  call generate_hamr()
  call generate_uniform_kmesh()
  call compute_DMFT_mu()
  call compute_G_loc()
  call comms_end

end program DMFTwDFT
