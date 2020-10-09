!-*- mode: F90 -*-!
!------------------------------------------------------------!
! This file is distributed as part of the Wannier90 code and !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the Wannier90      !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the Wannier90 code is www.wannier.org       !
!                                                            !
! The Wannier90 code is hosted on GitHub:                    !
!                                                            !
! https://github.com/wannier-developers/wannier90            !
!------------------------------------------------------------!

module constants

  !! This module contains the definitions of constants
  !! used in Wannier90 - both numerical constants such as pi
  !! and numerical convergence tolerances, but also physical
  !! constant such as the speed of light
  !!
  !! Values of the fundamental constants are taken from
  !! http://physics.nist.gov/cuu/Constants/index.html
  !! By default CODATA2010 is used (CODATA2006 can be selected
  !! using an appropriate compile-time flag (see INSTALL guide)

  implicit none

  private

  !~~ GENERIC CONSTANTS ~~!
!aam_2012-04-11; fix to run on MacBook Air
  integer, parameter, public          :: dp = kind(1.0d0)
  !! double precision
!  integer, parameter, public          :: dp = selected_real_kind(14,200)
!  integer, parameter, public          :: dp = selected_real_kind(15,300)
  real(kind=dp), parameter, public    :: pi = 3.141592653589793238462643383279_dp
  !! $$\pi$$
  real(kind=dp), parameter, public    :: twopi = 2*pi
  !! $$2\pi$$
  complex(kind=dp), parameter, public :: cmplx_i = (0.0_dp, 1.0_dp)
  !! i as a complex variable
  complex(kind=dp), parameter, public :: cmplx_0 = (0.0_dp, 0.0_dp)
  !! 0 as a complex variable
  complex(kind=dp), parameter, public :: cmplx_1 = (1.0_dp, 0.0_dp)
  !! 1 as a complex variable

  !~~ NUMERICAL CONVERGENCE CONSTANTS ~~!
  real(kind=dp), parameter, public    :: eps2 = 1.0e-2_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: eps5 = 1.0e-5_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: eps6 = 1.0e-6_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: eps7 = 1.0e-7_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: eps8 = 1.0e-8_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: eps10 = 1.0e-10_dp
  !! numerical convergence constant
  real(kind=dp), parameter, public    :: smearing_cutoff = 10._dp
  !! Cutoff for the smearing functions
  real(kind=dp), parameter, public    :: min_smearing_binwidth_ratio = 2._dp
  !! Don't smear but simply add the contribution to the
  !! relevant bin if the smearing/binwidth ratio is smaller than this value

end module constants
