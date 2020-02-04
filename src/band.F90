!------------------------------------------------------------!
!                                                            !
!                       DMFTwDFT                             !
!                                                            !
!          DMFT self-consistent calculation using            !
!                  Wannier functions                         !
!                                                            !
! Please cite XXXXX                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program DMFTwDFT_band
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
  call Read_sig_inp_real()
  call Read_dmft_params()
  call generate_hamr()
  call generate_klist()
  call compute_dos()
  call compute_Gk()
  call comms_end

end program DMFTwDFT_band
