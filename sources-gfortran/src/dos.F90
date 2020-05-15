!------------------------------------------------------------!
!                                                            !
!                       DMFTwDFT                             !
!                                                            !
!          DMFT self-consistent calculation using            !
!                  Wannier functions                         !
!                                                            !
! Please cite XXXXX                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program DMFTwDFT_dos 
  !!! main program
  USE constants 
  USE comms
  USE read_inputs
  USE generate_ham
  USE generate_kpts
  USE dmft_ksum 

  implicit none

  call comms_setup
  call Read_wan_chk() 
  call Read_wan_eig() 
  !write(*,*) 'hi'
  call Read_sig_inp_real()
  !write(*,*) 'hi2'
  call Read_dmft_params()
  !write(*,*) 'hi3'
  call generate_hamr()
  call generate_uniform_kmesh()
  call compute_dos()
  call comms_end

end program DMFTwDFT_dos
