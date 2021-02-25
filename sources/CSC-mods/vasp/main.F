#include "symbol.inc"
#define usgrid
!****************** PROGRAM VASP  Version 5.0 (f90)********************
! RCS:  $Id: main.F,v 1.18 2003/06/27 13:22:18 kresse Exp kresse $
! Vienna Ab initio total energy and Molecular-dynamics Program
!            written  by   Kresse Georg
!                     and  Juergen Furthmueller
! Georg Kresse                       email: Georg.Kresse@univie.ac.at
! Juergen Furthmueller               email: furth@ifto.physik.uni-jena.de
! Institut fuer Materialphysik         voice: +43-1-4277-51402
! Uni Wien, Sensengasse 8/12           fax:   +43-1-4277-9514 (or 9513)
! A-1090 Wien, AUSTRIA                 http://cms.mpi.univie.ac.at/kresse
!
! This program comes without any waranty.
! No part of this program must be distributed, modified, or supplied
! to any other person for any reason whatsoever
! without prior written permission of the Institut of Materials Science
! University Vienna
!
! This program performs total energy calculations using
! a selfconsistency cylce (i.e. mixing + iterative matrix diagonal.)
! or a direct optimisation of the one electron wavefunctions
! most of the algorithms implemented are described in
! G. Kresse and J. Furthmueller
!  Efficiency of ab--initio total energy calculations for
!   metals and semiconductors using a plane--wave basis set
!  Comput. Mat. Sci. 6,  15-50 (1996)
! G. Kresse and J. Furthmueller
!  Efficient iterative schemes for ab--initio total energy
!   calculations using a plane--wave basis set
!   Phys. Rev. B 54, 11169 (1996)
!
! The iterative matrix diagonalization is based
! a) on the conjugated gradient eigenvalue minimisation proposed by
!  D.M. Bylander, L. Kleinmann, S. Lee, Phys Rev. B 42, 1394 (1990)
! and is a variant of an algorithm proposed by
!  M.P. Teter, M.C. Payne and D.C. Allan, Phys. Rev. B 40,12255 (1989)
!  T.A. Arias, M.C. Payne, J.D. Joannopoulos, Phys Rev. B 45,1538(1992)
! b) or the residual vector minimization method (RMM-DIIS) proposed by
!  P. Pulay,  Chem. Phys. Lett. 73, 393 (1980).
!  D. M. Wood and A. Zunger, J. Phys. A, 1343 (1985)
! For the mixing a Broyden/Pulay like method is used (see for instance):
!  D. Johnson, Phys. Rev. B 38, 12807 (1988)
!
! The program can use normconserving PP,
! generalised ultrasoft-PP (Vanderbilt-PP Vanderbilt Phys Rev B 40,
! 12255 (1989)) and PAW (P.E. Bloechl, Phys. Rev. B{\bf 50}, 17953 (1994))
! datasets. Partial core corrections can be handled
! Spin and GGA and exact exchange functionals are implemented
!
! The units used in the programs are electron-volts and angstroms.
! The unit cell is arbitrary, and arbitrary species of ions are handled.
! A full featured symmetry-code is included, and calculation of
! Monkhorst-Pack special-points is possible (tetrahedron method can be
! used as well). This part was written by J. Furthmueller.
!
! The original version was written by  M.C. Payne
! at Professor J. Joannopoulos research  group at the MIT
! (3000 lines, excluding FFT, July 1989)
! The program was completely rewritten and vasply extended by
! Kresse Georg (gK) and Juergen Furthmueller. Currently the
! code has about 200000 source lines
!
!** The following parts have been taken from other programs
! - Tetrahedron method (original author unknown)
!
! please refer to the README file to learn about new features
! notes on singe-precision:
! USAGE NOT RECOMMENDED DUE TO FINITE DIFFERENCES IN FEW SPECIFIC
! FORCE-SUBROUTINE
! (except for native 64-bit-REAL machines like CRAY style machines)
!**********************************************************************

      PROGRAM VAMP
      USE prec

      USE charge
      USE pseudo
      USE lattice
      USE steep
      USE us
      USE pawm
      USE pot
      USE force
      USE fileio
      USE nonl_high
      USE rmm_diis
      USE ini
      USE ebs
      USE wave_high
      USE choleski
      USE mwavpre
      USE mwavpre_noio
      USE msphpro
      USE broyden
      USE msymmetry
      USE subrot
      USE melf
      USE base
      USE mpimy
      USE mgrid
      USE mkpoints
      USE constant
      USE setexm
      USE poscar
      USE wave
      USE hamil
      USE main_mpi
      USE chain
      USE pardens
      USE finite_differences
      USE LDAPLUSU_MODULE
      USE cl
      USE Constrained_M_modular
      USE writer
      USE sym_prec
      USE elpol
      USE mdipol
      USE wannier
      USE vaspxml
      USE full_kpoints
      USE kpoints_change
      USE fock
      USE compat_gga
      USE mlr_main
      USE mlrf_main
      USE mlr_optic
      USE pwkli
      USE gridq
      USE twoelectron4o
      USE dfast
      USE aedens
      USE xi
      USE subrotscf
      USE pead
      USE egrad
      USE hamil_high
      USE morbitalmag
      USE relativistic
      USE rhfatm
      USE meta
      USE mkproj
      USE classicfields
      USE rpa_force
! Thomas Bucko's code contributions
#ifdef tbdyn
      USE random_seeded
      USE dynconstr
#endif
      USE vdwforcefield
      USE dimer_heyden
      USE dvvtrajectory

      USE mlwf
#ifdef VASP2WANNIER90
      USE dmft
      USE crpa
#endif
      USE chgfit
      USE stockholder
      USE mlr_main_nmr
      USE hyperfine
      USE wannier_interpolation
      USE auger
      USE dmatrix

      USE lcao
      USE wnpr
! solvation__
      USE solvation
! solvation__
      USE locproj
#ifdef PROFILING
      USE profiling
#endif
! elphon_
      USE elphon
! elphon
! bexternal__
      USE bexternal
! bexternal__
      USE openmp
#ifdef CUDA_GPU
      USE cuda_interface
      USE main_gpu
#endif
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

!=======================================================================
!  a small set of parameters might be set here
!  but this is really rarely necessary :->
!=======================================================================
!----I/O-related things (adapt on installation or for special purposes)
!     IU6    overall output ('console protocol'/'OUTCAR' I/O-unit)
!     IU0    very important output ('standard [error] output I/O-unit')
!     IU5    input-parameters ('standard input'/INCAR I/O-unit)
!     ICMPLX size of complex items (in bytes/complex item)
!     MRECL  maximum record length for direct access files
!            (if no restictions set 0 or very large number)
      INTEGER,PARAMETER :: ICMPLX=16,MRECL=10000000

!=======================================================================
!  structures
!=======================================================================
      TYPE (potcar),ALLOCATABLE :: P(:)
      TYPE (wavedes)     WDES
      TYPE (nonlr_struct), TARGET :: NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (wavespin)    W          ! wavefunction
      TYPE (wavespin)    W_F        ! wavefunction for all bands simultaneous
      TYPE (wavespin)    W_G        ! same as above
      TYPE (wavefun)     WUP
      TYPE (wavefun)     WDW
      TYPE (wavefun)     WTMP       ! temporary
      TYPE (latt)        LATT_CUR
      TYPE (latt)        LATT_INI
      TYPE (type_info)   T_INFO
      TYPE (dynamics)    DYN
      TYPE (info_struct) INFO
      TYPE (in_struct)   IO
      TYPE (mixing)      MIX
      TYPE (kpoints_struct) KPOINTS
      TYPE (symmetry)    SYMM
      TYPE (grid_3d)     GRID       ! grid for wavefunctions
      TYPE (grid_3d)     GRID_SOFT  ! grid for soft chargedensity
      TYPE (grid_3d)     GRIDC      ! grid for potentials/charge
      TYPE (grid_3d)     GRIDUS     ! very find grid temporarily used in us.F
      TYPE (grid_3d)     GRIDB      ! Broyden grid
      TYPE (transit)     B_TO_C     ! index table between GRIDB and GRIDC
      TYPE (transit)     SOFT_TO_C  ! index table between GRID_SOFT and GRIDC
      TYPE (transit)     C_TO_US    ! index table between GRID_SOFT and GRIDC
      TYPE( prediction)  PRED
      TYPE (smear_struct) SMEAR_LOOP
      TYPE (paco_struct) PACO
      TYPE (energy)      E
      TYPE (ham_handle)  HAMILTONIAN
      TYPE (tau_handle)  KINEDEN

#ifdef tbdyn
       TYPE (gadget_io)   g_io
       INTEGER :: SEED(3),SEED_INIT(3)
       INTEGER :: DTVALUES(8)
       INTEGER :: IDUMLONG
       INTEGER :: CLOCK
       INTEGER, PARAMETER :: SEED1_MAX=900000000, K_SEED=3
       INTEGER, PARAMETER :: SEED2_MAX=1000000
#endif

      INTEGER :: NGX,NGY,NGZ,NGXC,NGYC,NGZC
      INTEGER :: NRPLWV,LDIM,LMDIM,LDIM2,LMYDIM
      INTEGER :: IRMAX,IRDMAX,ISPIND
      INTEGER :: NPLWV,MPLWV,NPLWVC,MPLWVC,NTYPD,NIOND,NIONPD,NTYPPD
      INTEGER :: NEDOS
      LOGICAL :: LNBANDS=.TRUE.
      INTEGER :: TIU6, TIU0
      INTEGER :: ISPECIAL=0         ! allows to select special undocumented features
      INTEGER :: MDALGO=0           ! dublicates MDALGO in tbdyn
!=======================================================================
!  begin array dimensions ...
!=======================================================================
#if defined(CUDA_GPU) && defined(USE_PINNED_MEMORY)
      COMPLEX(q),POINTER    :: CHTOT(:,:)    ! charge-density in real / reciprocal space
      TYPE(C_PTR)           :: CHTOT_PTR
      RGRID ,POINTER:: SV(:,:)               ! soft part of local potential
      TYPE(C_PTR) :: SV_PTR
      RGRID :: fake
      COMPLEX(q) :: fakec
#else
      COMPLEX(q),ALLOCATABLE:: CHTOT(:,:)    ! charge-density in real / reciprocal space
      RGRID  ,ALLOCATABLE:: SV(:,:)          ! soft part of local potential
#endif
      COMPLEX(q),ALLOCATABLE:: CHTOTL(:,:)   ! old charge-density
      RGRID     ,ALLOCATABLE:: DENCOR(:)     ! partial core
      COMPLEX(q),ALLOCATABLE:: CVTOT(:,:)    ! local potential
      COMPLEX(q),ALLOCATABLE:: CSTRF(:,:)    ! structure-factor
!-----non-local pseudopotential parameters
      OVERLAP,ALLOCATABLE:: CDIJ(:,:,:,:) ! strength of PP
      OVERLAP,ALLOCATABLE:: CQIJ(:,:,:,:) ! overlap of PP
      OVERLAP,ALLOCATABLE:: CRHODE(:,:,:,:) ! augmentation occupancies
!-----elements required for mixing in PAW method
      REAL(q)   ,ALLOCATABLE::   RHOLM(:,:),RHOLM_LAST(:,:)
!-----charge-density and potential on small grid
      COMPLEX(q),ALLOCATABLE:: CHDEN(:,:)    ! pseudo charge density
!-----description how to go from one grid to the second grid
!-----density of states
      REAL(q)   ,ALLOCATABLE::  DOS(:,:),DOSI(:,:)
      REAL(q)   ,ALLOCATABLE::  DDOS(:,:),DDOSI(:,:)
!-----local l-projected wavefunction characters
      REAL(q)   ,ALLOCATABLE::   PAR(:,:,:,:,:),DOSPAR(:,:,:,:)
!  all-band-simultaneous-update arrays
      GDEF   ,POINTER::   CHF(:,:,:,:),CHAM(:,:,:,:)
!  optics stuff
      GDEF   ,ALLOCATABLE::   NABIJ(:,:)
!
      LOGICAL :: LVCADER

!-----remaining mainly work arrays
      COMPLEX(q), ALLOCATABLE,TARGET :: CWORK1(:),CWORK2(:),CWORK(:,:)
      TYPE (wavefun1)    W1            ! current wavefunction
      TYPE (wavedes1)    WDES1         ! descriptor for one k-point

      GDEF, ALLOCATABLE  ::  CPROTM(:),CMAT(:,:)
!=======================================================================
!  a few fixed size (or small) arrays
!=======================================================================
!-----Forces and stresses
      REAL(q)   VTMP(3), XCSIF(3,3), EWSIF(3,3), TSIF(3,3), D2SIF(3,3)
!-----forces on ions
      REAL(q)   ,ALLOCATABLE::  EWIFOR(:,:), TIFOR(:,:)
!-----data for STM simulation (Bardeen)
      REAL(q)  STM(7)
!-----Temporary data for tutorial messages ...
      INTEGER,PARAMETER :: NTUTOR=1000
      REAL(q)     RTUT(NTUTOR),RDUM
      INTEGER  ITUT(NTUTOR),IDUM
      COMPLEX(q)  CDUM  ; LOGICAL  LDUM
!=======================================================================
!  end array dimensions ...
!=======================================================================
      INTEGER NTYP_PP      ! number of types on POTCAR file

      INTEGER I,J,N,NT,K
!---- used for creation of param.inc
      REAL(q)    WFACT,PSRMX,PSDMX
      REAL(q)    XCUTOF,YCUTOF,ZCUTOF

!---- timing information
      INTEGER IERR

      INTEGER NORDER   !   order of smearing
!---- a few logical and string variables
      LOGICAL    LTMP,LSTOP2
      LOGICAL    LPAW           ! paw is used
      LOGICAL    LPARD          ! partial band decomposed charge density
      LOGICAL    LREALLOCATE    ! reallocation of proj operators required
      LOGICAL    L_NO_US        ! no ultrasoft PP
      LOGICAL    LADDGRID       ! additional support grid


      LOGICAL    LBERRY         ! calculate electronic polarisation

#ifdef libbeef
      LOGICAL    LBEEFENS       ! calculate BEEF-vdW ensemble XC energies
      LOGICAL    LBEEFBAS       ! only print basis energies
#endif

      CHARACTER (LEN=40)  SZ
      CHARACTER (LEN=1)   CHARAC
      CHARACTER (LEN=5)   IDENTIFY
!-----parameters for sphpro.f
      INTEGER :: LDIMP,LMDIMP,LTRUNC=3
!=======================================================================
! All COMMON blocks
!=======================================================================
      INTEGER IXMIN,IXMAX,IYMIN,IYMAX,IZMIN,IZMAX
      COMMON /WAVCUT/ IXMIN,IXMAX,IYMIN,IYMAX,IZMIN,IZMAX

      REAL(q)  RHOTOT(4)
      INTEGER(8) IL,I1,I2_0,I3,I4
#ifdef gammareal
      CHARACTER (LEN=80),PARAMETER :: VASP = &
        'vasp.5.4.4.18Apr17-6-g9f103f2a35' // ' ' // &
        '(build ' // __DATE__// ' ' //__TIME__// ') ' // &
        'gamma-only'
#else
      CHARACTER (LEN=80),PARAMETER :: VASP = &
        'vasp.5.4.4.18Apr17-6-g9f103f2a35' // ' ' // &
        '(build ' // __DATE__// ' ' //__TIME__// ') ' // &
        'complex'
#endif

      INTEGER  ISYMOP,NROT,IGRPOP,NROTK,INVMAP,NPCELL
      REAL(q)  GTRANS,AP
      LOGICAL, EXTERNAL :: USE_OEP_IN_GW

      COMMON /SYMM/   ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

#ifdef libbeef
      LOGICAL  LBEEFCALCBASIS
      COMMON /BEEFENS/ LBEEFCALCBASIS,LBEEFBAS
#endif

      !integer, parameter :: dp = kind(1.0d0)
      !real(kind=dp),allocatable :: DMFT_eval_(:,:)
      !complex(kind=dp),allocatable :: DMFT_evec_(:,:,:)
      INTEGER NB, NK, NP, NBG
      REAL(q) tot_n
      COMPLEX(q),allocatable :: WaveW(:,:,:),WaveW_tot(:,:,:)
      COMPLEX(q),allocatable :: WProj(:,:,:),WProj_tot(:,:,:)
      REAL(q),allocatable :: FERWE(:,:)

      integer eig_i,eig_k,num_wannier

      INTERFACE
        subroutine Compute_DMFT(n_kpts_loc,n_wann,kpt_dft,wght_dft,band_win_loc,DMFT_eval,DMFT_evec)

        implicit none
        !integer, intent(in) :: mpi_comm_world,mpi_sum,mpi_integer
        integer, parameter :: dp = kind(1.0d0)
        integer, intent(in) :: n_kpts_loc,n_wann
        real(kind=dp), intent(in) :: kpt_dft(3,n_kpts_loc)
        real(kind=dp), intent(in) :: wght_dft(n_kpts_loc)
        integer, intent(out) :: band_win_loc(2,n_kpts_loc)
        real(kind=dp), intent(out) :: DMFT_eval(n_wann,n_kpts_loc)
        complex(kind=dp), intent(out) :: DMFT_evec(n_wann,n_wann,n_kpts_loc)
        end subroutine
        subroutine Compute_DMFT_from_amn(n_kpts_loc,n_wann,kpt_dft,wght_dft,band_win_loc,DMFT_eval,DMFT_evec)

        implicit none
        !integer, intent(in) :: mpi_comm_world,mpi_sum,mpi_integer
        integer, parameter :: dp = kind(1.0d0)
        integer, intent(in) :: n_kpts_loc,n_wann
        real(kind=dp), intent(in) :: kpt_dft(3,n_kpts_loc)
        real(kind=dp), intent(in) :: wght_dft(n_kpts_loc)
        integer, intent(out) :: band_win_loc(2,n_kpts_loc)
        real(kind=dp), intent(out) :: DMFT_eval(n_wann,n_kpts_loc)
        complex(kind=dp), intent(out) :: DMFT_evec(n_wann,n_wann,n_kpts_loc)
        end subroutine
      END INTERFACE

#ifdef PROFILING
!=======================================================================
!  initialise profiling
!=======================================================================
      CALL INIT_PROFILING
#endif
!=======================================================================
!  initialise / set constants and parameters ...
!=======================================================================

#ifdef libbeef
      LBEEFCALCBASIS = .FALSE.
#endif

      IO%LOPEN =.TRUE.  ! open all files with file names
      IO%IU0   = 6
      IO%IU6   = 8
      IO%IU5   = 7

      IO%ICMPLX=ICMPLX
      IO%MRECL =MRECL
      PRED%ICMPLX=ICMPLX

#ifdef tbdyn
      g_io%REPORT=66
      g_io%REFCOORD=67
      g_io%CONSTRAINTS=69
      g_io%STRUCTINPUT=533
      g_io%PENALTY=534
#endif
!
! with all the dynamic libraries (blas, lapack, scalapack etc.
! VASP has a size of at least 30 Mbyte
!
      CALL INIT_FINAL_TIMING
      CALL REGISTER_ALLOCATE(30000000._q, "base")

! switch off kill
!     CALL sigtrp()

      NPAR=1
      IUXML_SET=20

!$!-----------------------------------------------------------------------
!$!  initialise openMP
!$!-----------------------------------------------------------------------
!$! prevent nesting in openMP
!$! read thread related stuff from INCAR

#if defined(MPI) || defined(MPI_CHAIN)
      CALL INIT_MPI(NPAR,IO)
      NODE_ME= COMM%NODE_ME
      IONODE = COMM%IONODE
      IF (NODE_ME/=IONODE) THEN
         IUXML_SET=-1
      ENDIF
      IF (KIMAGES>0) IUXML_SET=-1
#endif

#ifdef CUDA_GPU
      CALL GPU_BANNER(IO)
#endif

      TIU6 = IO%IU6
      TIU0 = IO%IU0
      CALL START_XML( IUXML_SET, "vasprun.xml" )

!-----------------------------------------------------------------------
!  open Files
!-----------------------------------------------------------------------
      IF (IO%IU0>=0) WRITE(TIU0,*) VASP
#if defined(makeparam)
      IF (IO%IU6/=6 .AND. IO%IU6>0) &
      OPEN(UNIT=IO%IU6,FILE=DIR_APP(1:DIR_LEN)//'OUTPAR',STATUS='UNKNOWN')
#elif defined(makekpoints)
      IF (IO%IU6/=6 .AND. IO%IU6>0) &
      OPEN(UNIT=IO%IU6,FILE=DIR_APP(1:DIR_LEN)//'OUTKPT',STATUS='UNKNOWN')
#else
      IF (IO%IU6/=6 .AND. IO%IU6>0) &
      OPEN(UNIT=IO%IU6,FILE=DIR_APP(1:DIR_LEN)//'OUTCAR',STATUS='UNKNOWN')
#endif
      OPEN(UNIT=18,FILE=DIR_APP(1:DIR_LEN)//'CHGCAR',STATUS='UNKNOWN')
#ifdef logflow
      OPEN(UNIT=19,FILE=DIR_APP(1:DIR_LEN)//'JOBFLOW',ACCESS='APPEND',STATUS='UNKNOWN')
      WRITE(19,'(A)') 'SECTION = "----------------------------------------------------------------------------------------------------"'
#endif

#ifndef makekpoints
      CALL OPENWAV(IO, COMM)

#ifndef makeparam
      io_begin
      IF (KIMAGES==0) THEN
      OPEN(UNIT=22,FILE=DIR_APP(1:DIR_LEN)//'EIGENVAL',STATUS='UNKNOWN')
      OPEN(UNIT=13,FILE=DIR_APP(1:DIR_LEN)//'CONTCAR',STATUS='UNKNOWN')
      OPEN(UNIT=16,FILE=DIR_APP(1:DIR_LEN)//'DOSCAR',STATUS='UNKNOWN')
      OPEN(UNIT=17,FILE=DIR_APP(1:DIR_LEN)//'OSZICAR',STATUS='UNKNOWN')
      OPEN(UNIT=60,FILE=DIR_APP(1:DIR_LEN)//'PCDAT',STATUS='UNKNOWN')
      OPEN(UNIT=61,FILE=DIR_APP(1:DIR_LEN)//'XDATCAR',STATUS='UNKNOWN')
      OPEN(UNIT=70,FILE=DIR_APP(1:DIR_LEN)//'CHG',STATUS='UNKNOWN')
#ifdef tbdyn
      OPEN(UNIT=g_io%REPORT,FILE=DIR_APP(1:DIR_LEN)//'REPORT',STATUS='UNKNOWN')
#endif
      ENDIF
      io_end
#endif
      IF (IO%IU6>=0) WRITE(IO%IU6,*) VASP
      CALL XML_GENERATOR
#if defined(MPI)
      CALL PARSE_GENERATOR_XML(VASP//" parallel")
#else
      CALL PARSE_GENERATOR_XML(VASP//" serial")
#endif
      CALL MY_DATE_AND_TIME(IO%IU6)
      CALL XML_CLOSE_TAG

      CALL WRT_DISTR(IO%IU6)
#endif

! unit for extrapolation of wavefunction
      PRED%IUDIR =21
! unit for broyden mixing
      MIX%IUBROY=23
! unit for total potential
      IO%IUVTOT=62

 130  FORMAT (5X, //, &
     &'----------------------------------------------------', &
     &'----------------------------------------------------'//)

 140  FORMAT (5X, //, &
     &'----------------------------------------- Iteration ', &
     &I4,'(',I4,')  ---------------------------------------'//)
!-----------------------------------------------------------------------
! read header of POSCAR file to get NTYPD, NTYPDD, NIOND and NIONPD
!-----------------------------------------------------------------------
      CALL RD_POSCAR_HEAD(LATT_CUR, T_INFO, &
     &           NIOND,NIONPD, NTYPD,NTYPPD, IO%IU0, IO%IU6)

      ALLOCATE(T_INFO%ATOMOM(3*NIOND),T_INFO%RWIGS(NTYPPD),T_INFO%ROPT(NTYPD),T_INFO%POMASS(NTYPD), &
               T_INFO%DARWIN_V(NTYPD), T_INFO%DARWIN_R(NTYPD),T_INFO%VCA(NTYPD))

      IF (IO%IU6>=0) THEN
         WRITE(TIU6,130)
         WRITE(TIU6,*)'INCAR:'
      ENDIF
!  first scan of POTCAR to get LDIM, LMDIM, LDIM2 ...
      LDIM =16
      LDIM2=(LDIM*(LDIM+1))/2
      LMDIM=64

      ALLOCATE(P(NTYPD))
      T_INFO%POMASS=0
      T_INFO%RWIGS=0

      INFO%NLSPLINE=.FALSE.
!-----------------------------------------------------------------------
! read pseudopotentials
!-----------------------------------------------------------------------
      CALL RD_PSEUDO(INFO,P, &
     &           NTYP_PP,NTYPD,LDIM,LDIM2,LMDIM, &
     &           T_INFO%POMASS,T_INFO%RWIGS,T_INFO%TYPE,T_INFO%VCA, &
     &           IO%IU0,IO%IU6,-1,LPAW)

!-----------------------------------------------------------------------
! read INCAR
!-----------------------------------------------------------------------
      CALL XML_TAG("incar")

      CALL READER( &
          IO%IU5,IO%IU0,IO%INTERACTIVE,INFO%SZNAM1,INFO%ISTART,INFO%IALGO,MIX%IMIX,MIX%MAXMIX,MIX%MREMOVE, &
          MIX%AMIX,MIX%BMIX,MIX%AMIX_MAG,MIX%BMIX_MAG,MIX%AMIN, &
          MIX%WC,MIX%INIMIX,MIX%MIXPRE,MIX%MIXFIRST,IO%LFOUND,INFO%LDIAG,INFO%LSUBROT,INFO%LREAL,IO%LREALD,IO%LPDENS, &
          DYN%IBRION,INFO%ICHARG,INFO%INIWAV,INFO%NELM,INFO%NELMALL,INFO%NELMIN,INFO%NELMDL,INFO%EDIFF,DYN%EDIFFG, &
          DYN%NSW,DYN%ISIF,PRED%IWAVPR,SYMM%ISYM,DYN%NBLOCK,DYN%KBLOCK,INFO%ENMAX,DYN%POTIM,DYN%TEBEG, &
          DYN%TEEND,DYN%NFREE, &
          PACO%NPACO,PACO%APACO,T_INFO%NTYP,NTYPD,DYN%SMASS,SCALEE,T_INFO%POMASS, &
          T_INFO%DARWIN_V,T_INFO%DARWIN_R,T_INFO%VCA,LVCADER, &
          T_INFO%RWIGS,INFO%NELECT,INFO%NUP_DOWN,INFO%TIME, &
          KPOINTS%EMIN,KPOINTS%EMAX,KPOINTS%EFERMI,KPOINTS%ISMEAR,KPOINTS%SPACING,KPOINTS%LGAMMA, &
          DYN%PSTRESS,INFO%NDAV, &
          KPOINTS%SIGMA,KPOINTS%LTET,INFO%WEIMIN,INFO%EBREAK,INFO%DEPER,IO%NWRITE,INFO%LCORR, &
          IO%IDIOT,T_INFO%NIONS,T_INFO%NTYPP,IO%LMUSIC,IO%LOPTICS,STM, &
          INFO%ISPIN,T_INFO%ATOMOM,NIOND,IO%LWAVE,IO%LDOWNSAMPLE,IO%LCHARG,IO%LVTOT,IO%LVHAR,INFO%SZPREC, &
          INFO%ENAUG,IO%LORBIT,IO%LELF,T_INFO%ROPT,INFO%ENINI, &
          NGX,NGY,NGZ,NGXC,NGYC,NGZC,NBANDS,NEDOS,NBLK,LATT_CUR, &
          LPLANE_WISE,LCOMPAT,LMAX_CALC,SET_LMAX_MIX_TO,WDES%NSIM,LPARD,LPAW,LADDGRID, &
          WDES%LNONCOLLINEAR,WDES%LSORBIT,WDES%SAXIS,INFO%LMETAGGA, &
          WDES%LSPIRAL,WDES%LZEROZ,WDES%QSPIRAL,WDES%LORBITALREAL, &
          INFO%LASPH,INFO%TURBO,INFO%IRESTART,INFO%NREBOOT,INFO%NMIN,INFO%EREF, &
          INFO%NLSPLINE,ISPECIAL,MDALGO &
#ifdef libbeef
         ,LBEEFENS,LBEEFBAS &
#endif
         )
#ifdef tbdyn
      SEED=0
      !c user provided SEED for random number generator
      CALL RDATAB(IO%LOPEN,INCAR,IO%IU5,'RANDOM_SEED','=','#',';','I', &
            SEED,RDUM,CDUM,LDUM,CHARAC,N,K_SEED,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
                          ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
            WRITE(IO%IU0,*)'Error reading item ''RANDOM_SEED'' from file INCAR.'
         STOP
      ELSEIF (IERR==3) THEN
         CALL DATE_AND_TIME(VALUES=DTVALUES)

         !c with this choice we'r sure not to exceed
         !c SEED1_MAX (maximal value we can get here is 893581000)
         IDUMLONG=  DTVALUES(3)*24*60*60*300 &
                +DTVALUES(5)*60*60*1000 &
                +DTVALUES(6)*60*1000 &
                +DTVALUES(7)*1000 &
                +DTVALUES(8)
         SEED(1)=MOD(IDUMLONG,SEED1_MAX)
      ELSEIF ((SEED(1)<0).OR.(SEED(1)>SEED1_MAX).OR. &
               (SEED(2)<0).OR.(SEED(2)>SEED2_MAX).OR. &
                (SEED(3)<0)) THEN
         CALL VTUTOR('E','RANDOMSEED',RTUT,1, &
     &            SEED,3,CDUM,1,LDUM,1,IO%IU0,IO%IDIOT)
         CALL VTUTOR('S','RANDOMSEED',RTUT,1, &
     &            SEED,3,CDUM,1,LDUM,1,IO%IU6,IO%IDIOT)
      ELSE
         CALL XML_INCAR_V('RANDOM_SEED','I',SEED,RDUM,CDUM,LDUM,CHARAC,N)
      ENDIF
      CALLMPI( M_bcast_i(COMM_WORLD,SEED,K_SEED))
      CALL RANE_ION(RDUM,PUT=SEED(1:K_SEED))
      SEED_INIT=SEED
#endif

      KPOINTS%NKPX=MAX(1,CEILING(LATT_CUR%BNORM(1)*PI*2/KPOINTS%SPACING))
      KPOINTS%NKPY=MAX(1,CEILING(LATT_CUR%BNORM(2)*PI*2/KPOINTS%SPACING))
      KPOINTS%NKPZ=MAX(1,CEILING(LATT_CUR%BNORM(3)*PI*2/KPOINTS%SPACING))

      CALL GGA_COMPAT_MODE(IO%IU5, IO%IU0, LCOMPAT)

      IF (WDES%LNONCOLLINEAR) THEN
         INFO%ISPIN = 1
      ENDIF
! METAGGA not implemented for non collinear magnetism
!      IF (WDES%LNONCOLLINEAR .AND. INFO%LMETAGGA) THEN
!         WRITE(*,*) 'METAGGA for non collinear magnetism not supported.'
!         WRITE(*,*) 'exiting VASP; sorry for the inconveniences.'
!         STOP
!      ENDIF
!-MM- Spin spirals require LNONCOLLINEAR=.TRUE.
      IF (.NOT.WDES%LNONCOLLINEAR .AND. WDES%LSPIRAL) THEN
         WRITE(*,*) 'Spin spirals require LNONCOLLINEAR=.TRUE. '
         WRITE(*,*) 'exiting VASP; sorry dude!'
         STOP
      ENDIF
!-MM- end of addition

      IF (LCOMPAT) THEN
              CALL VTUTOR('W','VASP.4.4',RTUT,1, &
     &                 ITUT,1,CDUM,1,LDUM,1,IO%IU6,IO%IDIOT)
              CALL VTUTOR('W','VASP.4.4',RTUT,1, &
     &                 ITUT,1,CDUM,1,LDUM,1,IO%IU0,IO%IDIOT)
      ENDIF
! WRITE out an advice if some force dependent ionic algorithm and METAGGA
! or ASPH
!      IF ((INFO%LMETAGGA) .AND. &
!     &      (DYN%IBRION>0 .OR. (DYN%IBRION==0 .AND. DYN%SMASS/=-2))) THEN
!         CALL VTUTOR('A','METAGGA and forces',RTUT,1, &
!     &                 ITUT,1,CDUM,1,(/INFO%LASPH, INFO%LMETAGGA /),2, &
!     &                 IO%IU0,IO%IDIOT)
!      ENDIF
! The meaning of LVTOT has changed w.r.t. previous VASP version,
! therefore we write a warning
      IF (IO%LVTOT) THEN
         CALL VTUTOR('W','LVTOT',RTUT,1, &
     &                 ITUT,1,CDUM,1,LDUM,1,IO%IU0,IO%IDIOT)
      ENDIF

#if defined(MPI) || defined(MPI_CHAIN)
      IF ( REAL(COMM%NCPU,q)/COMM_INB%NCPU/COMM_INB%NCPU>4 .OR. &
           REAL(COMM%NCPU,q)/COMM_INB%NCPU/COMM_INB%NCPU<0.25_q) THEN
           CALL VTUTOR('W','NPAR efficiency',RTUT,1, &
           ITUT,1,CDUM,1,LDUM,1,IO%IU6,IO%IDIOT)
           CALL VTUTOR('W','NPAR efficiency',RTUT,1, &
           ITUT,1,CDUM,1,LDUM,1,IO%IU0,IO%IDIOT)
      ENDIF
      IF (KPAR.GT.1) THEN
              CALL VTUTOR('W','KPAR',RTUT,1, &
     &                 ITUT,1,CDUM,1,LDUM,1,IO%IU6,IO%IDIOT)
              CALL VTUTOR('W','KPAR',RTUT,1, &
     &                 ITUT,1,CDUM,1,LDUM,1,IO%IU0,IO%IDIOT)
      END IF
#endif
!-----------------------------------------------------------------------
! core level shift related items (parses INCAR)
!-----------------------------------------------------------------------
      CALL INIT_CL_SHIFT(IO%IU5,IO%IU0, T_INFO%NIONS, T_INFO%NTYP )
! Berry phase read INCAR
      CALL READER_ADD_ON(IO%IU5,IO%IU0,LBERRY,IGPAR,NPPSTR, &
            INFO%ICHARG,KPOINTS%ISMEAR,KPOINTS%SIGMA)
! Do we want to write the AE-densities?
      CALL INIT_AEDENS(IO%IU0,IO%IU5)

      ISPIND=INFO%ISPIN

      DYN%TEMP =DYN%TEBEG
      INFO%RSPIN=3-INFO%ISPIN

      CALL RESPONSE_READER(IO%IU5, IO%IU6, IO%IU0)
      CALL CRPA_READER(IO%IU5, IO%IU6, IO%IU0, NBANDS, L2E4W, L2E4W_ALL ,IO%LDOWNSAMPLE)
      CALL PEAD_READER(IO%IU5, IO%IU6, IO%IU0)
      CALL DMATRIX_READER(IO%IU5, IO%IU6, IO%IU0)
      CALL XC_FOCK_READER(IO%IU5, IO%IU0, IO%IU6, INFO%SZPREC, DYN%ISIF, SYMM%ISYM, INFO%IALGO, &
         LATT_CUR%OMEGA, T_INFO%NTYP, T_INFO%NIONS, MIX%IMIX, MIX%AMIX, MIX%AMIX_MAG, MIX%BMIX, MIX%BMIX_MAG, IO%LVTOT)
      CALL EGRAD_READER(IO%IU5, IO%IU6, IO%IU0, T_INFO%NTYP)
      CALL CLASSICFIELDS_READER(IO%IU5, IO%IU6, IO%IU0)
      CALL HYPERFINE_READER(IO%IU5, IO%IU6, IO%IU0, T_INFO%NTYP)
!-----------------------------------------------------------------------
! loop over different smearing parameters
!-----------------------------------------------------------------------
#ifndef makekpoints
      SMEAR_LOOP%ISMCNT=0
      IF (KPOINTS%ISMEAR==-3) THEN
        IF(IO%IU6>=0)   WRITE(TIU6,7219)
 7219   FORMAT('   Loop over smearing-parameters in INCAR')
        CALL RDATAB(IO%LOPEN,INCAR,IO%IU5,'SMEARINGS','=','#',';','F', &
     &            IDUM,SMEAR_LOOP%SMEARS(1),CDUM,LDUM,CHARAC,N,200,IERR)
        IF ((IERR/=0).OR.((IERR==0).AND. &
     &          ((N<2).OR.(N>200).OR.(MOD(N,2)/=0)))) THEN
           IF (IO%IU0>=0) &
           WRITE(TIU0,*)'Error reading item ''SMEARINGS'' from file INCAR.'
           STOP
        ENDIF
        SMEAR_LOOP%ISMCNT=N/2
        DYN%NSW   =SMEAR_LOOP%ISMCNT+1
        DYN%KBLOCK=DYN%NSW
        KPOINTS%LTET  =.TRUE.
        DYN%IBRION=-1
        KPOINTS%ISMEAR=-5
      ENDIF
!=======================================================================
!  now read in Pseudopotential
!  modify the potential if required (core level shifts)
!=======================================================================
      LMDIM=0
      LDIM=0
      DO NT=1,NTYP_PP
        LMDIM=MAX(LMDIM,P(NT)%LMMAX)
        LDIM =MAX(LDIM ,P(NT)%LMAX)
      END DO
      CALL DEALLOC_PP(P,NTYP_PP)

      LDIM2=(LDIM*(LDIM+1))/2
      LMYDIM=9
! second scan with correct setting
      CALL RD_PSEUDO(INFO,P, &
     &           NTYP_PP,NTYPD,LDIM,LDIM2,LMDIM, &
     &           T_INFO%POMASS,T_INFO%RWIGS,T_INFO%TYPE,T_INFO%VCA, &
     &           IO%IU0,IO%IU6,IO%NWRITE,LPAW)
      CALL CL_MODIFY_PP( NTYP_PP, P, INFO%ENAUG )
! now check everything
      CALL POST_PSEUDO(NTYPD,NTYP_PP,T_INFO%NTYP,T_INFO%NIONS,T_INFO%NITYP,T_INFO%VCA,P,INFO, &
     &        IO%LREALD,T_INFO%ROPT, IO%IDIOT,IO%IU6,IO%IU0,LMAX_CALC,L_NO_US,WDES%LSPIRAL)
      CALL LDIM_PSEUDO(IO%LORBIT, NTYPD, P, LDIMP, LMDIMP)
! check the difference between ENINI and ENMAX for spin spiral calculations
      IF (WDES%LSPIRAL) CALL CHECK_SPIRAL_ENCUT(WDES,INFO,LATT_CUR,IO)
!-----------------------------------------------------------------------
! LDA+U initialisation (parses INCAR)
!-----------------------------------------------------------------------
      CALL LDAU_READER(T_INFO%NTYP,IO%IU5,IO%IU0)
      IF ( (LGW0 .OR. LGW .OR. LscQPGW) .AND. USELDApU()) THEN
         CALL VTUTOR('E','Double counting',RTUT,1, &
     &               ITUT,1,CDUM,1,LDUM,1,IO%IU0,IO%IDIOT)
         CALL VTUTOR('E','Double counting',RTUT,1, &
     &               ITUT,1,CDUM,1,LDUM,1,IO%IU6,IO%IDIOT)
      ENDIF

      IF (USELDApU().OR.LCALC_ORBITAL_MOMENT()) &
     &   CALL INITIALIZE_LDAU(T_INFO%NIONS,T_INFO%NTYP,P,WDES%LNONCOLLINEAR,IO%IU0,IO%IDIOT)

      CALL SET_PAW_AUG(T_INFO%NTYP, P, IO%IU6, LMAX_CALC, LCOMPAT)
!-----------------------------------------------------------------------
! optics initialisation (parses INCAR)
!-----------------------------------------------------------------------
      IF (IO%LOPTICS) CALL SET_NABIJ_AUG(P,T_INFO%NTYP)
#endif
!-----------------------------------------------------------------------
!  Read UNIT=15: POSCAR Startjob and Continuation-job
!-----------------------------------------------------------------------
      CALL RD_POSCAR(LATT_CUR, T_INFO, DYN, &
     &           NIOND,NIONPD, NTYPD,NTYPPD, &
     &           IO%IU0,IO%IU6)

!-----------------------------------------------------------------------
! diverse INCAR readers
!-----------------------------------------------------------------------
      CALL RESPONSE_SET_ENCUT(INFO%ENMAX)
      CALL CONSTRAINED_M_READER(T_INFO,WDES,IO%IU0,IO%IU5)
      CALL WRITER_READER(IO%IU0,IO%IU5)
      CALL WANNIER_READER(IO%IU0,IO%IU5,IO%IU6,IO%IDIOT)
      CALL FIELD_READER(T_INFO,P,LATT_CUR,INFO%NELECT,IO%IU0,IO%IU5,IO%IU6)
      CALL LR_READER(INFO%EDIFF,IO%IU0,IO%IU5,IO%IU6)
      CALL ORBITALMAG_READER(IO%IU0, IO%IU5, IO%IU6, T_INFO%NIONS)
      CALL RHFATM_READER(IO)
      CALL XC_META_READER(IO,T_INFO%NTYP)
      CALL CHGFIT_READER(IO,T_INFO%NTYP)
      CALL STOCKHOLDER_READER(IO)
      CALL LJ_READER(IO)
      CALL MLWF_READER(IO%IU5,IO%IU6,IO%IU0,IO%LDOWNSAMPLE)
      CALL WNPR_READER(T_INFO%NIONS,IO%IU0,IO%IU5)
! locproj_
      CALL LPRJ_READER(T_INFO,P,IO%IU0)
! locproj_
! elphon_
      CALL ELPH_READER(IO%IU0,IO%IU5)
! elphon_
      CALL AUGER_READER(IO)
      CALL WANNIER_INTERPOL_READER(IO)
! solvation__
      CALL SOL_READER(T_INFO%NIONS,INFO%EDIFF,IO)
! solvation__
! bexternal__
      CALL BEXT_READER(IO%IU0,IO%IU5)
! bexternal__
      CALL CLASSICFIELDS_WRITE(IO%IU6)
!-----------------------------------------------------------------------
! exchange correlation table
!-----------------------------------------------------------------------
      CALL PUSH_XC_TYPE_FOR_GW ! switch now to AEXX=1.0 ; ALDAC = 0.0
      IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
         CALL SETUP_LDA_XC(2,IO%IU6,IO%IU0,IO%IDIOT)
      ELSE
         CALL SETUP_LDA_XC(1,IO%IU6,IO%IU0,IO%IDIOT)
      ENDIF
!-----------------------------------------------------------------------
! init all chains (INCAR reader)
!-----------------------------------------------------------------------
      CALL chain_init( T_INFO, IO)
!-----------------------------------------------------------------------
!xml finish copying parameters from INCAR to xml file
! no INCAR reading from here
      CALL XML_CLOSE_TAG("incar")
!-----------------------------------------------------------------------
#if CUDA_GPU
      CALL GPU_TEST(INFO,NPAR,IO)
#endif

      CALL COUNT_DEGREES_OF_FREEDOM( T_INFO, NDEGREES_OF_FREEDOM, &
          IO%IU6, IO%IU0, DYN%IBRION)

!-----for static calculations or relaxation jobs DYN%VEL is meaningless
      IF (DYN%INIT == -1) THEN
        CALL INITIO(T_INFO%NIONS,T_INFO%LSDYN,NDEGREES_OF_FREEDOM, &
               T_INFO%NTYP,T_INFO%ITYP,DYN%TEMP, &
               T_INFO%POMASS,DYN%POTIM, &
               DYN%POSION,DYN%VEL,T_INFO%LSFOR,LATT_CUR%A,LATT_CUR%B,DYN%INIT, COMM, IO%IU6)
#ifdef tbdyn
        ! TB dynamics will initialize the velocities for MDALGO>0
        IF (MDALGO==0) DYN%INIT=0
#else
        DYN%INIT=0
#endif
      ENDIF
      IF (DYN%IBRION/=0 .AND. DYN%IBRION/=40 .AND. DYN%IBRION/=44) THEN
          DYN%VEL=0._q
      ENDIF
      IF (IO%IU6>=0) THEN
         WRITE(TIU6,*)
         WRITE(TIU6,130)
      ENDIF


#ifdef tbdyn
      !c some MD methods (e.q. Langevin dynamics) do not conserve total momentum!
      IF ( T_INFO%LSDYN ) THEN
         CALL SET_SELECTED_VEL_ZERO(T_INFO, DYN%VEL,LATT_CUR)
      ELSE IF (MDALGO/=3 .AND. MDALGO/=31 .AND. DYN%IBRION/=44 .AND. DYN%IBRION/=40) THEN
         CALL SYMVEL_WARNING( T_INFO%NIONS, T_INFO%NTYP, T_INFO%ITYP, &
         T_INFO%POMASS, DYN%VEL, IO%IU6, IO%IU0 )
      ENDIF
#else
      IF ( T_INFO%LSDYN ) THEN
         CALL SET_SELECTED_VEL_ZERO(T_INFO, DYN%VEL,LATT_CUR)
      ELSE
         !tb beg
         !c do not remove CoM motion if iDM or DVV is to be used
         IF (DYN%IBRION/=44 .AND. DYN%IBRION/=40) THEN
           CALL SYMVEL_WARNING( T_INFO%NIONS, T_INFO%NTYP, T_INFO%ITYP, &
             T_INFO%POMASS, DYN%VEL, IO%IU6, IO%IU0 )
         ENDIF
         !tb end
      ENDIF
#endif

      CALL NEAREST_NEIGHBOUR(IO%IU6, IO%IU0, T_INFO, LATT_CUR, P%RWIGS)
!-----------------------------------------------------------------------
!  initialize the symmetry stuff
!-----------------------------------------------------------------------
      ALLOCATE(SYMM%ROTMAP(NIOND,1,1), &
               SYMM%TAU(NIOND,3), &
     &         SYMM%TAUROT(NIOND,3),SYMM%WRKROT(3*(NIOND+2)), &
     &         SYMM%PTRANS(NIOND+2,3),SYMM%INDROT(NIOND+2))
      IF (INFO%ISPIN==2) THEN
         ALLOCATE(SYMM%MAGROT(48,NIOND))
      ELSE
         ALLOCATE(SYMM%MAGROT(1,1))
      ENDIF
      ! break symmetry parallel to IGPAR
      IF (LBERRY) THEN
         LATT_CUR%A(:,IGPAR)=LATT_CUR%A(:,IGPAR)*(1+TINY*10)
         CALL LATTIC(LATT_CUR)
      ENDIF
! Rotate the initial magnetic moments to counter the clockwise
! rotation brought on by the spin spiral
      IF (WDES%LSPIRAL) CALL COUNTER_SPIRAL(WDES%QSPIRAL,T_INFO%NIONS,T_INFO%POSION,T_INFO%ATOMOM)

      IF (SYMM%ISYM>0) THEN
! Finite temperature allows no symmetry by definition ...
         NCDIJ=INFO%ISPIN
         IF (WDES%LNONCOLLINEAR) NCDIJ=4
         CALL INISYM(LATT_CUR%A,DYN%POSION,DYN%VEL,T_INFO%LSFOR, &
                     T_INFO%LSDYN,T_INFO%NTYP,T_INFO%NITYP,NIOND, &
                     SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP, &
                     SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
                     SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,NCDIJ,IO%IU6)
         IF (IO%NWRITE>=3) CALL WRTSYM(T_INFO%NIONS,NIOND,SYMM%PTRANS,SYMM%ROTMAP,SYMM%MAGROT,NCDIJ,IO%IU6)
      ELSE
! ... so take nosymm!
         CALL NOSYMM(LATT_CUR%A,T_INFO%NTYP,T_INFO%NITYP,NIOND,SYMM%PTRANS, &
        &   SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%MAGROT,INFO%ISPIN,IO%IU6)
      END IF
      CALL SETUP_PRIMITIVE_CELL(LATT_CUR%A, DYN%POSION, IO%IU0)
!=======================================================================
!  Read UNIT=14: KPOINTS
!  number of k-points and k-points in reciprocal lattice
!=======================================================================
      IF(IO%IU6>=0)  WRITE(TIU6,*)
      ! use full k-point grid if finite differences are used or
      ! linear response is applied
      IF (DYN%IBRION==5 .OR. DYN%IBRION==6 .OR. DYN%IBRION==7.OR. DYN%IBRION==8 &
          .OR. LEPSILON .OR. LVEL .OR. KINTER/=0 .OR. LMAGBLOCH  &
          .OR. LCHIMAG .OR. LTIME_EVOLUTION) CALL USE_FULL_KPOINTS
      ! apply preferentially time inversion symmetry to generate orbitals at -k
      IF (WDES%LORBITALREAL) CALL USE_TIME_INVERSION

      IF (LBERRY) THEN
         CALL RD_KPOINTS_BERRY(KPOINTS,NPPSTR,IGPAR, &
        &   LATT_CUR, &
        &   SYMM%ISYM>=0.AND..NOT.WDES%LSORBIT.AND..NOT.WDES%LSPIRAL, &
        &   IO%IU6,IO%IU0)
          IF (LBERRY) THEN
            LATT_CUR%A(:,IGPAR)=LATT_CUR%A(:,IGPAR)/(1+TINY*10)
            CALL LATTIC(LATT_CUR)
         ENDIF
      ELSE
#ifdef oldsym
         CALL SETUP_KPOINTS(KPOINTS,LATT_CUR, &
            SYMM%ISYM>=0.AND. &
            .NOT.WDES%LSORBIT.AND. &
            .NOT.WDES%LSPIRAL, &
            SYMM%ISYM<0,IO%IU6,IO%IU0)

         CALL SETUP_FULL_KPOINTS(KPOINTS,LATT_CUR,T_INFO%NIOND, &
            SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM, &
            SYMM%ISYM>=0.AND. &
            .NOT.WDES%LSORBIT.AND. &
            .NOT.WDES%LSPIRAL, &
            IO%IU6,IO%IU0,LSYMGRAD)
#else
         CALL SETUP_KPOINTS(KPOINTS,LATT_CUR, &
            SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
            SYMM%ISYM<0,IO%IU6,IO%IU0)

         CALL SETUP_FULL_KPOINTS(KPOINTS,LATT_CUR,T_INFO%NIOND, &
            SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM, &
            SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
            IO%IU6,IO%IU0, LSYMGRAD)
#endif
      ENDIF
      CALL SETUP_ORIG_KPOINTS

!=======================================================================
!  at this point we have enough information to
!  create a param.inc file
!=======================================================================
      XCUTOF =SQRT(INFO%ENMAX /RYTOEV)/(2*PI/(LATT_CUR%ANORM(1)/AUTOA))
      YCUTOF =SQRT(INFO%ENMAX /RYTOEV)/(2*PI/(LATT_CUR%ANORM(2)/AUTOA))
      ZCUTOF =SQRT(INFO%ENMAX /RYTOEV)/(2*PI/(LATT_CUR%ANORM(3)/AUTOA))
!
!  setup NGX, NGY, NGZ if required
!
! high precision do not allow for wrap around
      IF (INFO%SZPREC(1:1)=='h' .OR. INFO%SZPREC(1:1)=='a') THEN
        WFACT=4
      ELSE
! medium-low precision allow for small wrap around
        WFACT=3
      ENDIF
      GRID%NGPTAR(1)=XCUTOF*WFACT+0.5_q
      GRID%NGPTAR(2)=YCUTOF*WFACT+0.5_q
      GRID%NGPTAR(3)=ZCUTOF*WFACT+0.5_q
      IF (NGX /= -1)   GRID%NGPTAR(1)=  NGX
      IF (NGY /= -1)   GRID%NGPTAR(2)=  NGY
      IF (NGZ /= -1)   GRID%NGPTAR(3)=  NGZ
      CALL FFTCHK(GRID%NGPTAR)
!
!  setup NGXC, NGYC, NGZC if required
!
      IF (INFO%LOVERL) THEN
        IF (INFO%ENAUG==0) INFO%ENAUG=INFO%ENMAX*1.5_q
        IF (INFO%SZPREC(1:1)=='h') THEN
          WFACT=16._q/3._q
        ELSE IF (INFO%SZPREC(1:1)=='l') THEN
          WFACT=3
        ELSE
          WFACT=4
        ENDIF
        XCUTOF =SQRT(INFO%ENAUG /RYTOEV)/(2*PI/(LATT_CUR%ANORM(1)/AUTOA))
        YCUTOF =SQRT(INFO%ENAUG /RYTOEV)/(2*PI/(LATT_CUR%ANORM(2)/AUTOA))
        ZCUTOF =SQRT(INFO%ENAUG /RYTOEV)/(2*PI/(LATT_CUR%ANORM(3)/AUTOA))
        GRIDC%NGPTAR(1)=XCUTOF*WFACT
        GRIDC%NGPTAR(2)=YCUTOF*WFACT
        GRIDC%NGPTAR(3)=ZCUTOF*WFACT
        ! prec Single no double grid technique
        IF (INFO%SZPREC(1:1)=='s') THEN
           GRIDC%NGPTAR(1)=GRID%NGPTAR(1)
           GRIDC%NGPTAR(2)=GRID%NGPTAR(2)
           GRIDC%NGPTAR(3)=GRID%NGPTAR(3)
        ELSE IF (INFO%SZPREC(1:1)=='a' .OR. INFO%SZPREC(1:1)=='n') THEN
           GRIDC%NGPTAR(1)=GRID%NGPTAR(1)*2
           GRIDC%NGPTAR(2)=GRID%NGPTAR(2)*2
           GRIDC%NGPTAR(3)=GRID%NGPTAR(3)*2
        ENDIF
        IF (NGXC /= -1)  GRIDC%NGPTAR(1)=NGXC
        IF (NGYC /= -1)  GRIDC%NGPTAR(2)=NGYC
        IF (NGZC /= -1)  GRIDC%NGPTAR(3)=NGZC
        CALL FFTCHK(GRIDC%NGPTAR)
      ELSE
        GRIDC%NGPTAR(1)= 1
        GRIDC%NGPTAR(2)= 1
        GRIDC%NGPTAR(3)= 1
      ENDIF

      GRIDC%NGPTAR(1)=MAX(GRIDC%NGPTAR(1),GRID%NGPTAR(1))
      GRIDC%NGPTAR(2)=MAX(GRIDC%NGPTAR(2),GRID%NGPTAR(2))
      GRIDC%NGPTAR(3)=MAX(GRIDC%NGPTAR(3),GRID%NGPTAR(3))
      GRIDUS%NGPTAR=GRIDC%NGPTAR
      IF (LADDGRID) GRIDUS%NGPTAR=GRIDC%NGPTAR*2

      NGX = GRID %NGPTAR(1); NGY = GRID %NGPTAR(2); NGZ = GRID %NGPTAR(3)
      NGXC= GRIDC%NGPTAR(1); NGYC= GRIDC%NGPTAR(2); NGZC= GRIDC%NGPTAR(3)

      IF (NBANDS == -1) THEN
         LNBANDS=.FALSE.  ! NBANDS was not read from INCAR
         IF (WDES%LNONCOLLINEAR)  THEN
             NMAG=MAX(SUM(T_INFO%ATOMOM(1:T_INFO%NIONS*3-2:3)), &
                      SUM(T_INFO%ATOMOM(2:T_INFO%NIONS*3-1:3)), &
                      SUM(T_INFO%ATOMOM(3:T_INFO%NIONS*3:3)))
         ELSE IF (INFO%ISPIN > 1) THEN
             NMAG=SUM(T_INFO%ATOMOM(1:T_INFO%NIONS))
         ELSE
             NMAG=0
         ENDIF
         NMAG = (NMAG+1)/2
         NBANDS=MAX(NINT(INFO%NELECT+2)/2+MAX(T_INFO%NIONS/2,3),INT(0.6*INFO%NELECT))+NMAG
         IF (WDES%LNONCOLLINEAR) NBANDS = NBANDS*2
         NBANDS=((NBANDS+NPAR-1)/NPAR)*NPAR
      ENDIF

      IF (NBANDS/=((NBANDS+NPAR-1)/NPAR)*NPAR) THEN
         ITUT(1)=NBANDS
         ITUT(2)=((NBANDS+NPAR-1)/NPAR)*NPAR
         CALL VTUTOR('W','NBANDS changed',RTUT,1, &
     &        ITUT,2,CDUM,1,LDUM,1,IO%IU0,IO%IDIOT)
         CALL VTUTOR('W','NBANDS changed',RTUT,1, &
     &        ITUT,2,CDUM,1,LDUM,1,IO%IU6,IO%IDIOT)

      ENDIF

      NBANDS=((NBANDS+NPAR-1)/NPAR)*NPAR

      IF (INFO%EBREAK == -1) INFO%EBREAK=0.25_q*MIN(INFO%EDIFF,ABS(DYN%EDIFFG)/10)/NBANDS

      IF(INFO%TURBO==0)THEN
         IF (((.NOT. WDES%LNONCOLLINEAR).AND.  INFO%NELECT>REAL(NBANDS*2,KIND=q)).OR. &
                 ((WDES%LNONCOLLINEAR).AND.((INFO%NELECT*2)>REAL(NBANDS*2,KIND=q)))) THEN
            ITUT(1)=INFO%NELECT ; ITUT(2)=NBANDS
            CALL VTUTOR('E','Number of electrons',RTUT,1, &
     &                  ITUT,2,CDUM,1,LDUM,1,IO%IU0,IO%IDIOT)
            CALL VTUTOR('S','Number of electrons',RTUT,1, &
     &                  ITUT,2,CDUM,1,LDUM,1,IO%IU6,IO%IDIOT)
         ENDIF
      ELSE
         IF (( (.NOT. WDES%LNONCOLLINEAR).AND.  INFO%NELECT>REAL(NBANDS*2,KIND=q)).OR. &
                    ((WDES%LNONCOLLINEAR).AND.((INFO%NELECT*2)>REAL(NBANDS*2,KIND=q)))) THEN
            IF(KPOINTS%EFERMI==0) THEN
               ITUT(1)=INFO%NELECT ; ITUT(2)=NBANDS
               CALL VTUTOR('E','Number of electrons',RTUT,1, &
     &              ITUT,2,CDUM,1,LDUM,1,IO%IU0,IO%IDIOT)
               CALL VTUTOR('S','Number of electrons',RTUT,1, &
     &              ITUT,2,CDUM,1,LDUM,1,IO%IU6,IO%IDIOT)
            ENDIF
         ENDIF
      ENDIF

      NRPLWV=4*PI*SQRT(INFO%ENMAX /RYTOEV)**3/3* &
     &     LATT_CUR%OMEGA/AUTOA**3/(2*PI)**3*1.1_q+50
#ifdef gammareal
      NRPLWV=NRPLWV/2
#endif
      WDES%NRPLWV=NRPLWV
      PSRMX=0
      PSDMX=0
      DO NT=1,T_INFO%NTYP
        PSRMX=MAX(PSRMX,P(NT)%PSRMAX)
        PSDMX=MAX(PSDMX,P(NT)%PSDMAX)
      ENDDO
      IF (INFO%LREAL) THEN
       IRMAX=4*PI*PSRMX**3/3/(LATT_CUR%OMEGA/ &
     &        (GRID%NGPTAR(1)*GRID%NGPTAR(2)*GRID%NGPTAR(3)))+50
      ELSE
       IRMAX=1
      ENDIF
      IRDMAX=1
      IF (INFO%LOVERL) THEN
       IRDMAX=4*PI*PSDMX**3/3/(LATT_CUR%OMEGA/ &
     &        (GRIDC%NGPTAR(1)*GRIDC%NGPTAR(2)*GRIDC%NGPTAR(3)))+200
      ENDIF
#ifdef usgrid
       IRDMAX=4*PI*PSDMX**3/3/(LATT_CUR%OMEGA/ &
     &        (GRIDUS%NGPTAR(1)*GRIDUS%NGPTAR(2)*GRIDUS%NGPTAR(3)))+200
#endif

      NPLWV =NGX *NGY *NGZ;
      MPLWV =NGX *NGY *NGZ
      NPLWVC=NGXC*NGYC*NGZC;
      MPLWVC=xm(NGXC)*NGYC*zm(NGZC)
!=======================================================================
!  set the basic quantities in WDES
!  and set the grids
!=======================================================================

      WDES%ENMAX =INFO%ENMAX

      WDES%NB_PAR=NPAR
      WDES%NB_TOT=NBANDS
      WDES%NBANDS=NBANDS/NPAR
#ifdef MPI
      WDES%NB_LOW=COMM_INTER%NODE_ME
#else
      WDES%NB_LOW=1
#endif
      CALL INIT_KPOINT_WDES(WDES, KPOINTS )
      WDES%ISPIN =INFO%ISPIN
      WDES%COMM  =>COMM
      WDES%COMM_INB    =>COMM_INB
      WDES%COMM_INTER  =>COMM_INTER
      WDES%COMM_KIN    =>COMM_KIN
      WDES%COMM_KINTER =>COMM_KINTER
      NULLIFY( WDES%COMM_SHMEM )
#ifdef use_shmem
      WDES%COMM_SHMEM  =>COMM_SHMEM
      WDES%COMM_intra_node =>COMM_intra_node
      WDES%COMM_inter_node =>COMM_inter_node
#endif

      CALL  SET_FULL_KPOINTS(WDES%NKPTS_FOR_GEN_LAYOUT,WDES%VKPT)
      CALL  SET_FOCK_KPOINTS(WDES%NKDIM)

      IF (WDES%LNONCOLLINEAR) then
         WDES%NRSPINORS = 2
         INFO%RSPIN = 1
      ELSE
         WDES%NRSPINORS = 1
      ENDIF
      WDES%RSPIN = INFO%RSPIN

      CALL WDES_SET_NPRO(WDES,T_INFO,P,INFO%LOVERL)
!
! set up the descriptor for the initial wavefunctions
! (read from file)
      LATT_INI=LATT_CUR
! get header from WAVECAR file (LATT_INI is important)
! also set INFO%ISTART
      IF (INFO%ISTART > 0) THEN
        CALL INWAV_HEAD(WDES, LATT_INI, LATT_CUR, ENMAXI,INFO%ISTART, IO%IU0)
        IF (INFO%ISTART == 0 .AND. INFO%ICHARG == 0) INFO%ICHARG=2
      ENDIF

      CALL INIT_SCALAAWARE( WDES%NB_TOT, NRPLWV, WDES%COMM_KIN )

!=======================================================================
!  Write all important information
!=======================================================================
      IF (DYN%IBRION==5 .OR. DYN%IBRION==6 ) THEN
         DYN%NSW=4*(3*T_INFO%NIONS+9)+1
         IF (DYN%NFREE /= 1 .AND. DYN%NFREE /= 2 .AND. DYN%NFREE /= 4)  DYN%NFREE =2
         DYN%KBLOCK=DYN%NSW
      ENDIF

      IF (IO%IU6>=0) THEN

      WRITE(TIU6,130)
      WRITE(TIU6,7205) KPOINTS%NKPTS,WDES%NKDIM,WDES%NB_TOT,NEDOS, &
     &              T_INFO%NIONS,LDIM,LMDIM, &
     &              NPLWV,IRMAX,IRDMAX, &
     &              NGX,NGY,NGZ, &
     &              NGXC,NGYC,NGZC,GRIDUS%NGPTAR,T_INFO%NITYP

      XAU= (NGX*PI/(LATT_CUR%ANORM(1)/AUTOA))
      YAU= (NGY*PI/(LATT_CUR%ANORM(2)/AUTOA))
      ZAU= (NGZ*PI/(LATT_CUR%ANORM(3)/AUTOA))
      WRITE(TIU6,7211) XAU,YAU,ZAU
      XAU= (NGXC*PI/(LATT_CUR%ANORM(1)/AUTOA))
      YAU= (NGYC*PI/(LATT_CUR%ANORM(2)/AUTOA))
      ZAU= (NGZC*PI/(LATT_CUR%ANORM(3)/AUTOA))
      WRITE(TIU6,7212) XAU,YAU,ZAU

      ENDIF

 7211 FORMAT('   NGX,Y,Z   is equivalent  to a cutoff of ', &
     &           F6.2,',',F6.2,',',F6.2,' a.u.')
 7212 FORMAT('   NGXF,Y,Z  is equivalent  to a cutoff of ', &
     &           F6.2,',',F6.2,',',F6.2,' a.u.'/)

      XCUTOF =SQRT(INFO%ENMAX /RYTOEV)/(2*PI/(LATT_CUR%ANORM(1)/AUTOA))
      YCUTOF =SQRT(INFO%ENMAX /RYTOEV)/(2*PI/(LATT_CUR%ANORM(2)/AUTOA))
      ZCUTOF =SQRT(INFO%ENMAX /RYTOEV)/(2*PI/(LATT_CUR%ANORM(3)/AUTOA))
! high precision do not allow for wrap around
      IF (INFO%SZPREC(1:1)=='h'.OR.INFO%SZPREC(1:1)=='a') THEN
        WFACT=4
      ELSE
! medium-low precision allow for small wrap around
        WFACT=3
      ENDIF
      ITUT(1)=XCUTOF*WFACT+0.5_q
      ITUT(2)=YCUTOF*WFACT+0.5_q
      ITUT(3)=ZCUTOF*WFACT+0.5_q
      CALL FFTCHK(ITUT(1:3))
      IF (IO%IU6>=0 .AND. (ITUT(1)/=NGX .OR. ITUT(2)/=NGY .OR. ITUT(3)/=NGZ)) &
        WRITE(TIU6,72111) ITUT(1),ITUT(2),ITUT(3)

72111 FORMAT(' Based on PREC, I would recommend the following setting:'/ &
     &       '   dimension x,y,z NGX = ',I5,' NGY =',I5,' NGZ =',I5 /)

      IF (NGX<ITUT(1) .OR. NGY<ITUT(2) .OR. NGZ<ITUT(3)) THEN
               CALL VTUTOR('W','FFT-GRID IS NOT SUFFICIENT',RTUT,1, &
     &                  ITUT,3,CDUM,1,LDUM,1,IO%IU6,IO%IDIOT)
               CALL VTUTOR('W','FFT-GRID IS NOT SUFFICIENT',RTUT,1, &
     &                  ITUT,3,CDUM,1,LDUM,1,IO%IU0,IO%IDIOT)
      ENDIF


      AOMEGA=LATT_CUR%OMEGA/T_INFO%NIONS
      QF=(3._q*PI*PI*INFO%NELECT/(LATT_CUR%OMEGA))**(1._q/3._q)*AUTOA

! chose the mass so that the typical Nose frequency is 40 timesteps
!-----just believe this (or look out  for all this factors in  STEP)
      IF (DYN%SMASS==0)  DYN%SMASS= &
         ((40._q*DYN%POTIM*1E-15_q/2._q/PI/LATT_CUR%ANORM(1))**2)* &
         2.E20_q*BOLKEV*EVTOJ/AMTOKG*NDEGREES_OF_FREEDOM*MAX(DYN%TEBEG,DYN%TEEND)
!      IF (DYN%SMASS<0)  DYN%SMASS= &
!         ((ABS(DYN%SMASS)*DYN%POTIM*1E-15_q/2._q/PI/LATT_CUR%ANORM(1))**2)* &
!         2.E20_q*BOLKEV*EVTOJ/AMTOKG*NDEGREES_OF_FREEDOM*MAX(DYN%TEBEG,DYN%TEEND)

      SQQ=  DYN%SMASS*(AMTOKG/EVTOJ)*(1E-10_q*LATT_CUR%ANORM(1))**2
      SQQAU=SQQ/RYTOEV
      IF (DYN%SMASS>0) THEN
        WOSZI= SQRT(2*BOLKEV*DYN%TEMP*NDEGREES_OF_FREEDOM/SQQ)
      ELSE
        WOSZI=1E-30_q
      ENDIF
!-----initial temperature
      CALL EKINC(EKIN,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,DYN%POTIM,LATT_CUR%A,DYN%VEL)
      TEIN = 2*EKIN/BOLKEV/NDEGREES_OF_FREEDOM
!-----be carefull about division by 0
      DYN%NBLOCK=MAX(1,DYN%NBLOCK)
      DYN%KBLOCK=MAX(1,DYN%KBLOCK)
      IF (DYN%NSW<DYN%KBLOCK*DYN%NBLOCK) DYN%KBLOCK=1
      IF (DYN%NSW<DYN%KBLOCK*DYN%NBLOCK) DYN%NBLOCK=MAX(DYN%NSW,1)

      DYN%NSW=INT(DYN%NSW/DYN%NBLOCK/DYN%KBLOCK)*DYN%NBLOCK*DYN%KBLOCK
      IF (IO%IU6>=0) THEN

      WRITE(TIU6,7210) INFO%SZNAM1,T_INFO%SZNAM2

      WRITE(TIU6,7206) IO%NWRITE,INFO%SZPREC,INFO%ISTART,INFO%ICHARG,WDES%ISPIN,WDES%LNONCOLLINEAR, &
     &      WDES%LSORBIT, INFO%INIWAV, &
     &      INFO%LASPH,INFO%LMETAGGA, &
     &      INFO%ENMAX,INFO%ENMAX/RYTOEV,SQRT(INFO%ENMAX/RYTOEV), &
     &      XCUTOF,YCUTOF,ZCUTOF,INFO%ENINI, &
     &      INFO%ENAUG, &
     &      INFO%NELM,INFO%NELMIN,INFO%NELMDL,INFO%EDIFF,INFO%LREAL,INFO%NLSPLINE,LCOMPAT,GGA_COMPAT, &
     &      LMAX_CALC,SET_LMAX_MIX_TO,LFCI, &
     &      T_INFO%ROPT
      WRITE(TIU6,7204) &
     &      DYN%EDIFFG,DYN%NSW,DYN%NBLOCK,DYN%KBLOCK, &
     &      DYN%IBRION,DYN%NFREE,DYN%ISIF,PRED%IWAVPR,SYMM%ISYM,INFO%LCORR

      TMP=0
      IF (DYN%POTIM>0) TMP=1/(WOSZI*(DYN%POTIM*1E-15_q)/2._q/PI)

      WRITE(TIU6,7207) &
     &      DYN%POTIM,TEIN,DYN%TEBEG,DYN%TEEND, &
     &      DYN%SMASS,WOSZI,TMP,SQQAU,SCALEE, &
     &      PACO%NPACO,PACO%APACO,DYN%PSTRESS

      WRITE(TIU6,7215) (T_INFO%POMASS(NI),NI=1,T_INFO%NTYP)
      RTUT(1:T_INFO%NTYP)=P(1:T_INFO%NTYP)%ZVALF ! work around IBM bug
      WRITE(TIU6,7216) (RTUT(NI),NI=1,T_INFO%NTYP)
      WRITE(TIU6,7203) (T_INFO%RWIGS(NI),NI=1,T_INFO%NTYP)
      WRITE(TIU6,72031) (T_INFO%VCA(NI),NI=1,T_INFO%NTYP)

      WRITE(TIU6,7208) &
     &      INFO%NELECT,INFO%NUP_DOWN, &
     &      KPOINTS%EMIN,KPOINTS%EMAX,KPOINTS%EFERMI,KPOINTS%ISMEAR,KPOINTS%SIGMA

      WRITE(TIU6,7209) &
     &      INFO%IALGO,INFO%LDIAG,INFO%LSUBROT, &
     &      INFO%TURBO,INFO%IRESTART,INFO%NREBOOT,INFO%NMIN,INFO%EREF,&
     &      MIX%IMIX,MIX%AMIX,MIX%BMIX,MIX%AMIX_MAG,MIX%BMIX_MAG,MIX%AMIN, &
     &      MIX%WC,MIX%INIMIX,MIX%MIXPRE,MIX%MAXMIX, &
     &      INFO%WEIMIN,INFO%EBREAK,INFO%DEPER,INFO%TIME, &
     &      AOMEGA,AOMEGA/(AUTOA)**3, &
     &      QF,QF/AUTOA,QF**2*RYTOEV,QF**2,SQRT(QF/AUTOA/AUTOA*4/PI)
      WRITE(TIU6,*)
      WRITE(TIU6,7224) IO%LWAVE,IO%LDOWNSAMPLE,IO%LCHARG,IO%LVTOT,IO%LVHAR,IO%LELF,IO%LORBIT

      CALL WRITE_EFIELD(TIU6)
      ENDIF

 7210 FORMAT( &
     &       ' SYSTEM =  ',A40/ &
     &       ' POSCAR =  ',A40/)

 7205 FORMAT(//' Dimension of arrays:'/ &
     &       '   k-points           NKPTS = ',I6, &
     &       '   k-points in BZ     NKDIM = ',I6, &
     &       '   number of bands    NBANDS= ',I6/ &
     &       '   number of dos      NEDOS = ',I6, &
     &       '   number of ions     NIONS = ',I6/ &
     &       '   non local maximal  LDIM  = ',I6, &
     &       '   non local SUM 2l+1 LMDIM = ',I6/ &
     &       '   total plane-waves  NPLWV = ',I6/ &
     &       '   max r-space proj   IRMAX = ',I6, &
     &       '   max aug-charges    IRDMAX= ',I6/ &
     &       '   dimension x,y,z NGX = ',I5,' NGY =',I5,' NGZ =',I5/ &
     &       '   dimension x,y,z NGXF= ',I5,' NGYF=',I5,' NGZF=',I5/ &
     &       '   support grid    NGXF= ',I5,' NGYF=',I5,' NGZF=',I5/ &
     &       '   ions per type =            ',10I4/)

 7206 FORMAT(' Startparameter for this run:'/ &
     &       '   NWRITE = ',I6,  '    write-flag & timer' / &
     &       '   PREC   = ',A6,  '    normal or accurate (medium, high low for compatibility)'/ &
     &       '   ISTART = ',I6,  '    job   : 0-new  1-cont  2-samecut'/ &
     &       '   ICHARG = ',I6,  '    charge: 1-file 2-atom 10-const'/ &
     &       '   ISPIN  = ',I6,  '    spin polarized calculation?'/ &
     &       '   LNONCOLLINEAR = ',L6, ' non collinear calculations'/ &
     &       '   LSORBIT = ',L6, '    spin-orbit coupling'/ &
     &       '   INIWAV = ',I6,  '    electr: 0-lowe 1-rand  2-diag'/ &
     &       '   LASPH  = ',L6,  '    aspherical Exc in radial PAW'/ &
     &       '   METAGGA= ',L6,  '    non-selfconsistent MetaGGA calc.'// &
     &       ' Electronic Relaxation 1'/ &
     &       '   ENCUT  = ', &
     &              F6.1,' eV ',F6.2,' Ry  ',F6.2,' a.u. ', &
     &              3F6.2,'*2*pi/ulx,y,z'/ &
     &       '   ENINI  = ',F6.1,'     initial cutoff'/ &
     &       '   ENAUG  = ',F6.1,' eV  augmentation charge cutoff'/ &
     &       '   NELM   = ',I6,  ';   NELMIN=',I3,'; NELMDL=',I3, &
     &         '     # of ELM steps '    / &
     &       '   EDIFF  = ',E7.1,'   stopping-criterion for ELM'/ &
     &       '   LREAL  = ',L6,  '    real-space projection'     / &
     &       '   NLSPLINE    = ',L1,'    spline interpolate recip. space projectors'/ &
     &       '   LCOMPAT= ',L6,  '    compatible to vasp.4.4'/&
     &       '   GGA_COMPAT  = ',L1,'    GGA compatible to vasp.4.4-vasp.4.6'/&
     &       '   LMAXPAW     = ',I4,' max onsite density'/&
     &       '   LMAXMIX     = ',I4,' max onsite mixed and CHGCAR'/&
     &       '   VOSKOWN= ',I6,  '    Vosko Wilk Nusair interpolation'/&
     &      ('   ROPT   = ',4F10.5))
 7204 FORMAT( &
     &       ' Ionic relaxation'/ &
     &       '   EDIFFG = ',E7.1,'   stopping-criterion for IOM'/ &
     &       '   NSW    = ',I6,  '    number of steps for IOM' / &
     &       '   NBLOCK = ',I6,  ';   KBLOCK = ',I6, &
     &         '    inner block; outer block '/ &
     &       '   IBRION = ',I6, &
     &         '    ionic relax: 0-MD 1-quasi-New 2-CG'/ &
     &       '   NFREE  = ',I6,  &
     &         '    steps in history (QN), initial steepest desc. (CG)'/ &
     &       '   ISIF   = ',I6,  '    stress and relaxation' / &
     &       '   IWAVPR = ',I6, &
     &         '    prediction:  0-non 1-charg 2-wave 3-comb' / &
     &       '   ISYM   = ',I6, &
     &         '    0-nonsym 1-usesym 2-fastsym' / &
     &       '   LCORR  = ',L6, &
     &         '    Harris-Foulkes like correction to forces' /)

 7207 FORMAT( &
     &       '   POTIM  =' ,F7.4,'    time-step for ionic-motion'/ &
     &       '   TEIN   = ',F6.1,'    initial temperature'       / &
     &       '   TEBEG  = ',F6.1,';   TEEND  =',F6.1, &
     &               ' temperature during run'/ &
     &       '   SMASS  = ',F6.2,'    Nose mass-parameter (am)'/ &
     &       '   estimated Nose-frequenzy (Omega)   = ',E9.2, &
     &           ' period in steps =',F6.2,' mass=',E12.3,'a.u.'/ &
     &       '   SCALEE = ',F6.4,'    scale energy and forces'       / &
     &       '   NPACO  = ',I6,  ';   APACO  = ',F4.1, &
     &       '  distance and # of slots for P.C.'  / &
     &       '   PSTRESS= ',F6.1,' pullay stress'/)

!    &       '   damping for Cell-Motion     SIDAMP = ',F6.2/
!    &       '   mass for Cell-Motion        SIMASS = ',F6.2/

 7215 FORMAT('  Mass of Ions in am'/ &
     &       ('   POMASS = ',8F6.2))
 7216 FORMAT('  Ionic Valenz'/ &
     &       ('   ZVAL   = ',8F6.2))
 7203 FORMAT('  Atomic Wigner-Seitz radii'/ &
     &       ('   RWIGS  = ',8F6.2))
72031 FORMAT('  virtual crystal weights '/ &
     &       ('   VCA    = ',8F6.2))

 7208 FORMAT( &
     &       '   NELECT = ',F12.4,  '    total number of electrons'/ &
     &       '   NUPDOWN= ',F12.4,  '    fix difference up-down'// &
     &       ' DOS related values:'/ &
     &       '   EMIN   = ',F6.2,';   EMAX   =',F6.2, &
     &       '  energy-range for DOS'/ &
     &       '   EFERMI = ',F6.2,/ &
     &       '   ISMEAR =',I6,';   SIGMA  = ',F6.2, &
     &       '  broadening in eV -4-tet -1-fermi 0-gaus'/)

 7209 FORMAT( &
     &       ' Electronic relaxation 2 (details)'/  &
     &       '   IALGO  = ',I6,  '    algorithm'            / &
     &       '   LDIAG  = ',L6,  '    sub-space diagonalisation (order eigenvalues)' / &
     &       '   LSUBROT= ',L6,  '    optimize rotation matrix (better conditioning)' / &
     &       '   TURBO    = ',I6,  '    0=normal 1=particle mesh'/ &
     &       '   IRESTART = ',I6,  '    0=no restart 2=restart with 2 vectors'/ &
     &       '   NREBOOT  = ',I6,  '    no. of reboots'/ &
     &       '   NMIN     = ',I6,  '    reboot dimension'/ &
     &       '   EREF     = ',F6.2,  '    reference energy to select bands'/ &
     &       '   IMIX   = ',I6,  '    mixing-type and parameters'/ &
     &       '     AMIX     = ',F6.2,';   BMIX     =',F6.2/ &
     &       '     AMIX_MAG = ',F6.2,';   BMIX_MAG =',F6.2/ &
     &       '     AMIN     = ',F6.2/ &
     &       '     WC   = ',F6.0,';   INIMIX=',I4,';  MIXPRE=',I4,';  MAXMIX=',I4// &
     &       ' Intra band minimization:'/ &
     &       '   WEIMIN = ',F6.4,'     energy-eigenvalue tresh-hold'/ &
     &       '   EBREAK = ',E9.2,'  absolut break condition' / &
     &       '   DEPER  = ',F6.2,'     relativ break condition  ' // &
     &       '   TIME   = ',F6.2,'     timestep for ELM'          // &
     &       '  volume/ion in A,a.u.               = ',F10.2,3X,F10.2/ &
     &       '  Fermi-wavevector in a.u.,A,eV,Ry     = ',4F10.6/ &
     &       '  Thomas-Fermi vector in A             = ',2F10.6/)


 7224 FORMAT( &
     &       ' Write flags'/  &
     &       '   LWAVE        = ',L6,  '    write WAVECAR' / &
     &       '   LDOWNSAMPLE  = ',L6,  '    k-point downsampling of WAVECAR' / &
     &       '   LCHARG       = ',L6,  '    write CHGCAR' / &
     &       '   LVTOT        = ',L6,  '    write LOCPOT, total local potential' / &
     &       '   LVHAR        = ',L6,  '    write LOCPOT, Hartree potential only' / &
     &       '   LELF         = ',L6,  '    write electronic localiz. function (ELF)'/&
     &       '   LORBIT       = ',I6,  '    0 simple, 1 ext, 2 COOP (PROOUT), +10 PAW based schemes'//)

       CALL WRITE_CL_SHIFT(IO%IU6)
       IF (USELDApU()) CALL WRITE_LDApU(IO%IU6)

       CALL WRITE_FOCK(IO%IU6)

       CALL WRITE_BERRY_PARA(IO%IU6,LBERRY,IGPAR,NPPSTR)
       CALL LR_WRITER(IO%IU6)
       CALL WRITE_ORBITALMAG(IO%IU6)
       CALL WRITE_RESPONSE(IO%IU6)
       CALL PEAD_WRITER(IO%IU6)
       CALL CHGFIT_WRITER(IO)
! solvation__
       CALL SOL_WRITER(IO)
! solvation__

#ifdef logflow
       CLOSE(19)
#endif
       CALL XML_TAG("parameters")
       CALL XML_WRITER( &
          NPAR, &
          INFO%SZNAM1,INFO%ISTART,INFO%IALGO,MIX%IMIX,MIX%MAXMIX,MIX%MREMOVE, &
          MIX%AMIX,MIX%BMIX,MIX%AMIX_MAG,MIX%BMIX_MAG,MIX%AMIN, &
          MIX%WC,MIX%INIMIX,MIX%MIXPRE,MIX%MIXFIRST,IO%LFOUND,INFO%LDIAG,INFO%LSUBROT,INFO%LREAL,IO%LREALD,IO%LPDENS, &
          DYN%IBRION,INFO%ICHARG,INFO%INIWAV,INFO%NELM,INFO%NELMIN,INFO%NELMDL,INFO%EDIFF,DYN%EDIFFG, &
          DYN%NSW,DYN%ISIF,PRED%IWAVPR,SYMM%ISYM,DYN%NBLOCK,DYN%KBLOCK,INFO%ENMAX,DYN%POTIM,DYN%TEBEG, &
          DYN%TEEND,DYN%NFREE, &
          PACO%NPACO,PACO%APACO,T_INFO%NTYP,NTYPD,DYN%SMASS,SCALEE,T_INFO%POMASS,T_INFO%DARWIN_V,T_INFO%DARWIN_R,  &
          T_INFO%RWIGS,INFO%NELECT,INFO%NUP_DOWN,INFO%TIME,KPOINTS%EMIN,KPOINTS%EMAX,KPOINTS%EFERMI, &
          KPOINTS%ISMEAR,KPOINTS%SPACING,KPOINTS%LGAMMA,DYN%PSTRESS,INFO%NDAV, &
          KPOINTS%SIGMA,KPOINTS%LTET,INFO%WEIMIN,INFO%EBREAK,INFO%DEPER,IO%NWRITE,INFO%LCORR, &
          IO%IDIOT,T_INFO%NIONS,T_INFO%NTYPP,IO%LMUSIC,IO%LOPTICS,STM, &
          INFO%ISPIN,T_INFO%ATOMOM,NIOND,IO%LWAVE,IO%LDOWNSAMPLE,IO%LCHARG,IO%LVTOT,IO%LVHAR,INFO%SZPREC, &
          INFO%ENAUG,IO%LORBIT,IO%LELF,T_INFO%ROPT,INFO%ENINI, &
          NGX,NGY,NGZ,NGXC,NGYC,NGZC,NBANDS,NEDOS,NBLK,LATT_CUR, &
          LPLANE_WISE,LCOMPAT,LMAX_CALC,SET_LMAX_MIX_TO,WDES%NSIM,LPARD,LPAW,LADDGRID, &
          WDES%LNONCOLLINEAR,WDES%LSORBIT,WDES%SAXIS,INFO%LMETAGGA, &
          WDES%LSPIRAL,WDES%LZEROZ,WDES%QSPIRAL, &
          INFO%LASPH,WDES%LORBITALREAL, &
          INFO%TURBO,INFO%IRESTART,INFO%NREBOOT,INFO%NMIN,INFO%EREF, &
          INFO%NLSPLINE,ISPECIAL,MDALGO)

       CALL  XML_WRITE_GGA_COMPAT_MODE
       CALL  XML_WRITE_BERRY(LBERRY, IGPAR, NPPSTR)
       CALL  XML_WRITE_CL_SHIFT
       CALL  XML_WRITE_LDAU
       CALL  XML_WRITE_CONSTRAINED_M(T_INFO%NIONS)
       CALL  XML_WRITE_XC_FOCK
       CALL  XML_WRITE_LR
       CALL  XML_WRITE_ORBITALMAG
       CALL  XML_WRITE_RESPONSE
       CALL  XML_WRITE_CLASSICFIELDS
! solvation__
       CALL XML_WRITE_SOL
! solvation__

       CALL XML_CLOSE_TAG("parameters")
!=======================================================================
!  set some important flags and write out text information
!  DYN%IBRION        selects dynamic
!  INFO%LCORR =.TRUE. calculate Harris corrections to forces
!=======================================================================
       IF (MIX%AMIN>=0.1_q .AND. MAXVAL(LATT_CUR%ANORM(:))>50) THEN
          CALL VTUTOR('W','long lattice',RTUT,1, &
     &         ITUT,1,CDUM,1,LDUM,1,IO%IU6,IO%IDIOT)
          CALL VTUTOR('W','long lattice',RTUT,1, &
     &         ITUT,1,CDUM,1,LDUM,1,IO%IU0,IO%IDIOT)
       ENDIF
!---- relaxation related information
      IF (DYN%IBRION==10) THEN
         INFO%NELMDL=ABS(INFO%NELM)
         INFO%LCORR=.TRUE.
         IF (DYN%POTIM <= 0.0001_q ) DYN%POTIM=1E-20_q
      ENDIF

      IF (IO%IU6>=0) THEN

      WRITE(TIU6,130)

      IF (DYN%IBRION == -1) THEN
        WRITE(TIU6,*)'Static calculation'
      ELSE IF (DYN%IBRION==0) THEN
        WRITE(TIU6,*)'molecular dynamics for ions'
        IF (DYN%SMASS>0) THEN
          WRITE(TIU6,*)'  using nose mass (canonical ensemble)'
        ELSE IF (DYN%SMASS==-3) THEN
          WRITE(TIU6,*)'  using a microcanonical ensemble'
        ELSE IF (DYN%SMASS==-1) THEN
          WRITE(TIU6,*)'  scaling velocities every NBLOCK steps'
        ELSE IF (DYN%SMASS==-2) THEN
          WRITE(TIU6,*)'  keeping initial velocities unchanged'
        ENDIF
      ELSE IF (DYN%IBRION==1) THEN
           WRITE(TIU6,*)'quasi-Newton-method for relaxation of ions'
      ELSE IF (DYN%IBRION==2) THEN
           WRITE(TIU6,*)'conjugate gradient relaxation of ions'
      ELSE IF (DYN%IBRION==3) THEN
              WRITE(TIU6,*)'quickmin algorithm: (dynamic with friction)'
      ELSE IF (DYN%IBRION==5) THEN
              WRITE(TIU6,*)'finite differences'
      ELSE IF (DYN%IBRION==6) THEN
              WRITE(TIU6,*)'finite differences with symmetry'
      ELSE IF (DYN%IBRION==7) THEN
              WRITE(TIU6,*)'linear response'
      ELSE IF (DYN%IBRION==8) THEN
              WRITE(TIU6,*)'linear response using symmetry'
      ELSE IF (DYN%IBRION==10) THEN
           WRITE(TIU6,*)'relaxation of ions and charge simultaneously'
      ELSE IF (DYN%IBRION==11) THEN
           WRITE(TIU6,*)'interactive mode, write forces and read positions'
      !tb beg
      ELSE IF (DYN%IBRION==44) THEN
          WRITE(TIU6,*)'improved dimer method for transition state relaxation'
      !tb end
      ENDIF

      IF (DYN%IBRION/=-1 .AND. T_INFO%LSDYN) THEN
        WRITE(TIU6,*)'  using selective dynamics as specified on POSCAR'
        IF (.NOT.T_INFO%LDIRCO) THEN
          WRITE(TIU6,*)'  WARNING: If single coordinates had been '// &
     &                'selected the selection of coordinates'
          WRITE(TIU6,*)'           is made according to the '// &
     &                'corresponding   d i r e c t   coordinates!'
          WRITE(TIU6,*)'           Don''t support selection of '// &
     &                'single cartesian coordinates -- sorry ... !'
        ENDIF
      ENDIF
      ENDIF

      IF (INFO%ICHARG>=10) THEN
        INFO%LCHCON=.TRUE.
        IF(IO%IU6>=0)  WRITE(TIU6,*)'charge density and potential remain constant during run'
        MIX%IMIX=0
      ELSE
        INFO%LCHCON=.FALSE.
        IF(IO%IU6>=0)  WRITE(TIU6,*)'charge density and potential will be updated during run'
      ENDIF

      IF ((WDES%ISPIN==2).AND.INFO%LCHCON.AND.(DYN%IBRION/=-1)) THEN
         IF (IO%IU0>=0)  &
         WRITE(TIU0,*) &
          'Spin polarized Harris functional dynamics is a good joke ...'
         IF (IO%IU6>=0) &
         WRITE(TIU6,*) &
          'Spin polarized Harris functional dynamics is a good joke ...'
         STOP
      ENDIF
      IF (IO%IU6>=0) THEN
        IF (WDES%ISPIN==1 .AND. .NOT. WDES%LNONCOLLINEAR ) THEN
          WRITE(TIU6,*)'non-spin polarized calculation'
        ELSE IF ( WDES%LNONCOLLINEAR ) THEN
          WRITE(TIU6,*)'non collinear spin polarized calculation'
        ELSE
          WRITE(TIU6,*)'spin polarized calculation'
        ENDIF
      ENDIF

! paritial dos
      JOBPAR=1
!     IF (DYN%IBRION>=0) JOBPAR=0
      DO NT=1,T_INFO%NTYP
         IF (T_INFO%RWIGS(NT)<=0._q) JOBPAR=0
      ENDDO

!  INFO%LCDIAG  call EDDIAG after  eigenvalue optimization
!  INFO%LPDIAG  call EDDIAG before eigenvalue optimization
!  INFO%LDIAG   perform sub space rotation (when calling EDDIAG)
!  INFO%LORTHO  orthogonalization of wavefcuntions within optimization
!                     no Gramm-Schmidt required
!  INFO%LRMM    use RMM-DIIS minimization
!  INFO%LDAVID  use blocked Davidson
!  INFO%LCHCON  charge constant during run
!  INFO%LCHCOS  charge constant during band minimisation
!  INFO%LONESW  all band simultaneous

      INFO%LCHCOS=.TRUE.
      INFO%LONESW=.FALSE.
      INFO%LONESW_AUTO=.FALSE.
      INFO%LDAVID=.FALSE.
      INFO%LRMM  =.FALSE.
      INFO%LORTHO=.TRUE.
      INFO%LCDIAG=.FALSE.
      INFO%LPDIAG=.TRUE.
      INFO%LPRECONDH=.FALSE.
      INFO%IHARMONIC=0
      INFO%LEXACT_DIAG=.FALSE.

!  all bands conjugate gradient (CG) or damped orbital optimization
!  with fall back to DIIS when convergence is save
      IF (INFO%IALGO>=100) THEN
        IF (INFO%LDIAG) THEN
           IF(IO%IU6>=0)  WRITE(TIU6,*) 'Conjugate gradient for all bands (Freysoldt, et al. PRB 79, 241103 (2009))'
        ELSE
           IF(IO%IU6>=0)  WRITE(TIU6,*) 'Conjugate gradient for all bands (Kresse, et al. variant)'
        ENDIF
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Fall back to RMM-DIIS when convergence is save'
        INFO%IALGO=MOD(INFO%IALGO,10)
        INFO%LCHCOS=INFO%LCHCON
        INFO%LONESW=.TRUE.
        INFO%LONESW_AUTO=.TRUE.
        INFO%LRMM  =.TRUE.    ! this is tricky, set some usefull defaults for calling
        INFO%LORTHO=.FALSE.   !  routines in electron.F  (fall back is DIIS)
!  exact diagonalization
      ELSE IF (INFO%IALGO>=90) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Exact diagonalization'
        INFO%IALGO=MOD(INFO%IALGO,10)
        INFO%LEXACT_DIAG=.TRUE.
        INFO%LORTHO=.TRUE.
        INFO%LCDIAG=.FALSE.
        INFO%LPDIAG=.FALSE.
!  routines implemented in david_inner (Harmonic Jacobi Davidson, and Davidson)
      ELSE IF (INFO%IALGO>=80 .OR. (INFO%IALGO>=70 .AND. INFO%EREF==0)) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Davidson algorithm suitable for deep iteration (NRMM>8)'
        INFO%IALGO=MOD(INFO%IALGO,10)
        INFO%IHARMONIC=2
        INFO%LORTHO=.FALSE.
        INFO%LCDIAG=.TRUE.
        INFO%LPDIAG=.FALSE.
      ELSE IF (INFO%IALGO>=70) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Jacobi Davidson Harmonic (JDH) for inner eigenvalue problems'
        INFO%IALGO=INFO%IALGO-70
        INFO%IHARMONIC=1
        INFO%LORTHO=.FALSE.
        INFO%LCDIAG=.TRUE.
        INFO%LPDIAG=.FALSE.
!  RMM-DIIS + Davidson
      ELSE IF (INFO%IALGO>=60) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'RMM-DIIS sequential band-by-band and'
        IF(IO%IU6>=0)  WRITE(TIU6,*) ' variant of blocked Davidson during initial phase'
        INFO%IALGO=INFO%IALGO-60
        INFO%LRMM   =.TRUE.
        INFO%LDAVID =.TRUE.
        INFO%LORTHO =.FALSE.
        INFO%LDIAG  =.TRUE.        ! subspace rotation is always selected
!  all bands conjugate gradient (CG) or damped orbital optimization
      ELSE IF (INFO%IALGO>=50) THEN
        IF (INFO%LDIAG) THEN
           IF(IO%IU6>=0)  WRITE(TIU6,*) 'Conjugate gradient for all bands (Freysoldt, et al. PRB 79, 241103 (2009))'
        ELSE
           IF(IO%IU6>=0)  WRITE(TIU6,*) 'Conjugate gradient for all bands (Kresse, et al. variant)'
        ENDIF
        INFO%IALGO=MOD(INFO%IALGO,10)
        INFO%LCHCOS=INFO%LCHCON
        INFO%LONESW=.TRUE.
!  RMM-DIIS
      ELSE IF (INFO%IALGO>=40) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'RMM-DIIS sequential band-by-band'
        INFO%IALGO=MOD(INFO%IALGO,10)
        INFO%LRMM  =.TRUE.
        INFO%LORTHO=.FALSE.
!  blocked Davidson (Liu)
      ELSE IF (INFO%IALGO>=30) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Variant of blocked Davidson'
        INFO%IALGO=INFO%IALGO-30
        IF (INFO%LDIAG) THEN    ! if LDIAG is set
           IF(IO%IU6>=0)  WRITE(TIU6,*) 'Davidson routine will perform the subspace rotation'
           INFO%LCDIAG=.FALSE.  ! routine does the diagonalisation itself
           INFO%LPDIAG=.FALSE.  ! hence LPDIAG and LCDIAG are set to .FALSE.
        ENDIF
        INFO%LDAVID=.TRUE.
      ELSE IF (INFO%IALGO>=20) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Conjugate gradient sequential band-by-band (Teter, Alan, Payne)'
        INFO%IALGO  =INFO%IALGO-20
        INFO%LORTHO=.FALSE.
         CALL VTUTOR('S','IALGO8',RTUT,1,ITUT,1,CDUM,1,LDUM,1,IO%IU0,IO%IDIOT)
      ELSE IF (INFO%IALGO>=10) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Compatibility mode'
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Conjugate gradient sequential band-by-band (Teter, Alan, Payne)'
        INFO%IALGO=INFO%IALGO-10
        INFO%LCDIAG=.TRUE.
        INFO%LPDIAG=.FALSE.
         CALL VTUTOR('S','IALGO8',RTUT,1,ITUT,1,CDUM,1,LDUM,1,IO%IU0,IO%IDIOT)
      ELSE IF (INFO%IALGO>=5 .OR. INFO%IALGO==0) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Conjugate gradient sequential band-by-band (Teter, Alan, Payne)'
         CALL VTUTOR('S','IALGO8',RTUT,1,ITUT,1,CDUM,1,LDUM,1,IO%IU0,IO%IDIOT)
      ELSE IF (INFO%IALGO<0) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Performance tests'
      ELSE IF (INFO%IALGO <1) THEN
        IF (IO%IU0>=0) &
        WRITE(TIU0,*) 'Algorithms no longer implemented'
        STOP
      ELSE IF (INFO%IALGO==2) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'None: do nothing, only one-electron occupancies are recalculated'
        INFO%LDIAG =.FALSE.
      ELSE IF (INFO%IALGO==3) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Eigenval: update one-electron energies, occupancies fixed'
        INFO%LDIAG =.FALSE.
      ELSE IF (INFO%IALGO==4) THEN
        IF(IO%IU6>=0)  WRITE(TIU6,*) 'Subrot: only subspace diagonalization (rotation)'
      ENDIF

      SZ=' '
      IF (INFO%LCHCOS) THEN
         SZ='   charged. constant during bandupdate'
      ELSE
         INFO%LCORR=.FALSE.
!        IMIX used to be set to 0, removed
!         MIX%IMIX=0
      ENDIF

      IF (IO%IU6>=0) THEN

      IF (.NOT. INFO%LRMM .AND. .NOT. INFO%LDAVID) THEN
        IF (INFO%IALGO==5) THEN
          WRITE(TIU6,*)'steepest descent',SZ
        ELSEIF (INFO%IALGO==6) THEN
          WRITE(TIU6,*)'conjugated gradient',SZ
        ELSEIF (INFO%IALGO==7) THEN
          WRITE(TIU6,*)'preconditioned steepest descent',SZ
        ELSEIF (INFO%IALGO==8) THEN
          WRITE(TIU6,*)'preconditioned conjugated gradient',SZ
        ELSEIF (INFO%IALGO==0) THEN
          WRITE(TIU6,*)'preconditioned conjugated gradient (Jacobi prec)',SZ
        ENDIF
        IF (.NOT.INFO%LONESW) THEN
          WRITE(TIU6,*)'   band-by band algorithm'
        ENDIF
      ENDIF

      IF (INFO%LDIAG) THEN
        WRITE(TIU6,*)'perform sub-space diagonalisation'
      ELSE
        WRITE(TIU6,*)'perform Loewdin sub-space diagonalisation'
        WRITE(TIU6,*)'   ordering is kept fixed'
      ENDIF

      IF (INFO%LPDIAG) THEN
        WRITE(TIU6,*)'   before iterative eigenvector-optimisation'
      ELSE
        WRITE(TIU6,*)'   after iterative eigenvector-optimisation'
      ENDIF

      IF (MIX%IMIX==1 .OR. MIX%IMIX==2 .OR. MIX%IMIX==3) THEN
        WRITE(TIU6,*)'Kerker-like  mixing scheme'
      ELSE IF (MIX%IMIX==4) THEN
       WRITE(TIU6,'(A,F10.1)')' modified Broyden-mixing scheme, WC = ',MIX%WC
       IF (MIX%INIMIX==1) THEN
         WRITE(TIU6,'(A,F8.4,A,F12.4)') &
     &     ' initial mixing is a Kerker type mixing with AMIX =',MIX%AMIX, &
     &     ' and BMIX =',MIX%BMIX
       ELSE IF (MIX%INIMIX==2) THEN
         WRITE(TIU6,*)'initial mixing equals unity matrix (no mixing!)'
       ELSE
         WRITE(TIU6,'(A,F8.4)') &
     &     ' initial mixing is a simple linear mixing with ALPHA =',MIX%AMIX
       ENDIF
       IF (MIX%MIXPRE==1) THEN
         WRITE(TIU6,*)'Hartree-type preconditioning will be used'
       ELSE IF (MIX%MIXPRE==2) THEN
         WRITE(TIU6,'(A,A,F12.4)') &
     &     ' (inverse) Kerker-type preconditioning will be used', &
     &     ' corresponding to BMIX =',MIX%BMIX
       ELSE
         WRITE(TIU6,*)'no preconditioning will be used'
       ENDIF
      ELSE
        WRITE(TIU6,*)'no mixing'
      ENDIF
      IF (WDES%NB_TOT*2==NINT(INFO%NELECT)) THEN
        WRITE(TIU6,*)'2*number of bands equal to number of electrons'
        IF (MIX%IMIX/=0 .AND..NOT.INFO%LCHCOS) THEN
          WRITE(TIU6,*) &
     &      'WARNING: mixing without additional bands will not converge'
        ELSE IF (MIX%IMIX/=0) THEN
          WRITE(TIU6,*) 'WARNING: mixing has no effect'
        ENDIF

      ELSE
        WRITE(TIU6,*)'using additional bands ',INT(WDES%NB_TOT-INFO%NELECT/2)
        IF (KPOINTS%SIGMA<=0) THEN
          WRITE(TIU6,*) &
     &  'WARNING: no broadening specified (might cause bad convergence)'
        ENDIF
      ENDIF

      IF (INFO%LREAL) THEN
        WRITE(TIU6,*)'real space projection scheme for non local part'
      ELSE
        WRITE(TIU6,*)'reciprocal scheme for non local part'
      ENDIF

      IF (INFO%LCORE) THEN
        WRITE(TIU6,*)'use partial core corrections'
      ENDIF

      IF (INFO%LCORR) THEN
        WRITE(TIU6,*)'calculate Harris-corrections to forces ', &
     &              '  (improved forces if not selfconsistent)'
      ELSE
        WRITE(TIU6,*)'no Harris-corrections to forces '
      ENDIF

      IF (ISGGA()) THEN
        WRITE(TIU6,*)'use gradient corrections '
        IF (INFO%LCHCON) THEN
           IF (IO%IU0>=0) &
           WRITE(TIU0,*)'WARNING: stress and forces are not correct'
           WRITE(TIU6,*)'WARNING: stress and forces are not correct'
           WRITE(TIU6,*)' (second derivative of E(xc) not defined)'
        ENDIF
      ENDIF

      IF (INFO%LOVERL) THEN
         WRITE(TIU6,*)'use of overlap-Matrix (Vanderbilt PP)'
      ENDIF
      IF (KPOINTS%ISMEAR==-1) THEN
        WRITE(TIU6,7213) KPOINTS%SIGMA
 7213 FORMAT(' Fermi-smearing in eV        SIGMA  = ',F6.2)

      ELSE IF (KPOINTS%ISMEAR==-2) THEN
        WRITE(TIU6,7214)
 7214 FORMAT(' partial occupancies read from INCAR or WAVECAR (fixed during run)')
      ELSE IF (KPOINTS%ISMEAR==-4) THEN
        WRITE(TIU6,7222)
 7222 FORMAT(' Fermi weights with tetrahedron method witout', &
     &       ' Bloechl corrections')
      ELSE IF (KPOINTS%ISMEAR==-5) THEN
        WRITE(TIU6,7223)
 7223 FORMAT(' Fermi weights with tetrahedron method with', &
     &       ' Bloechl corrections')

      ELSE IF (KPOINTS%ISMEAR>0) THEN
        WRITE(TIU6,7217) KPOINTS%ISMEAR,KPOINTS%SIGMA
 7217 FORMAT(' Methfessel and Paxton  Order N=',I2, &
     &       ' SIGMA  = ',F6.2)
      ELSE
        WRITE(TIU6,7218) KPOINTS%SIGMA
 7218 FORMAT(' Gauss-broadening in eV      SIGMA  = ',F6.2)
      ENDIF

      WRITE(TIU6,130)
!=======================================================================
!  write out the lattice parameters
!=======================================================================
      WRITE(TIU6,7220) INFO%ENMAX,LATT_CUR%OMEGA, &
     &    ((LATT_CUR%A(I,J),I=1,3),(LATT_CUR%B(I,J),I=1,3),J=1,3), &
     &    (LATT_CUR%ANORM(I),I=1,3),(LATT_CUR%BNORM(I),I=1,3)

      WRITE(TIU6,*)

      IF (INFO%ISTART==1 .OR.INFO%ISTART==2) THEN

      WRITE(TIU6,*)'old parameters found on file WAVECAR:'
      WRITE(TIU6,7220) ENMAXI,LATT_INI%OMEGA, &
     &    ((LATT_INI%A(I,J),I=1,3),(LATT_INI%B(I,J),I=1,3),J=1,3)


      WRITE(TIU6,*)
 7220 FORMAT('  energy-cutoff  :  ',F10.2/ &
     &       '  volume of cell :  ',F10.2/ &
     &       '      direct lattice vectors',17X,'reciprocal lattice vectors'/ &
     &       3(2(3X,3F13.9)/) / &
     &       '  length of vectors'/ &
     &        (2(3X,3F13.9)/) /)

      ENDIF
!=======================================================================
!  write out k-points,weights,size & positions
!=======================================================================

 7104 FORMAT(' k-points in units of 2pi/SCALE and weight: ',A40)
 7105 FORMAT(' k-points in reciprocal lattice and weights: ',A40)
 7016 FORMAT(' position of ions in fractional coordinates (direct lattice) ')
 7017 FORMAT(' position of ions in cartesian coordinates  (Angst):')
 7009 FORMAT(1X,3F12.8,F12.3)
 7007 FORMAT(1X,3F12.8)

      WRITE(TIU6,7104) KPOINTS%SZNAMK

      DO NKP=1,KPOINTS%NKPTS
        VTMP(1)=WDES%VKPT(1,NKP)
        VTMP(2)=WDES%VKPT(2,NKP)
        VTMP(3)=WDES%VKPT(3,NKP)
        CALL DIRKAR(1,VTMP,LATT_CUR%B)
        WRITE(TIU6,7009) VTMP(1)*LATT_CUR%SCALE,VTMP(2)*LATT_CUR%SCALE, &
                  VTMP(3)*LATT_CUR%SCALE,KPOINTS%WTKPT(NKP)
      ENDDO

      WRITE(TIU6,*)
      WRITE(TIU6,7105) KPOINTS%SZNAMK
      DO NKP=1,KPOINTS%NKPTS
        WRITE(TIU6,7009) WDES%VKPT(1,NKP),WDES%VKPT(2,NKP),WDES%VKPT(3,NKP),KPOINTS%WTKPT(NKP)
      ENDDO
      WRITE(TIU6,*)

      WRITE(TIU6,7016)
      WRITE(TIU6,7007) ((DYN%POSION(I,J),I=1,3),J=1,T_INFO%NIONS)
      WRITE(TIU6,*)
      WRITE(TIU6,7017)

      DO J=1,T_INFO%NIONS
        VTMP(1)=DYN%POSION(1,J)
        VTMP(2)=DYN%POSION(2,J)
        VTMP(3)=DYN%POSION(3,J)
        CALL  DIRKAR(1,VTMP,LATT_CUR%A)
        WRITE(TIU6,7007) (VTMP(I),I=1,3)
      ENDDO
      WRITE(TIU6,*)
      CALL  WRITE_EULER(IO%IU6, WDES%LNONCOLLINEAR, WDES%SAXIS)

      WRITE(TIU6,130)
      ENDIF

      io_begin
      IF (KIMAGES==0) THEN
!=======================================================================
!  write out initial header for PCDAT, XDATCAR
!=======================================================================
      CALL PCDAT_HEAD(60,T_INFO, LATT_CUR, DYN, PACO, INFO%SZNAM1)
      IF (.NOT. (DYN%ISIF==3 .OR. DYN%ISIF>=7 )) THEN
        ! for DYN%ISIF=3 and >=7, the header is written in each step
        CALL XDAT_HEAD(61, T_INFO, LATT_CUR, DYN, INFO%SZNAM1)
      ENDIF
!=======================================================================
!  write out initial header for DOS
!=======================================================================
      JOBPAR_=JOBPAR
      IF (IO%LORBIT>=10 ) JOBPAR_=1

      WRITE(16,'(4I4)') T_INFO%NIONP,T_INFO%NIONS,JOBPAR_,WDES%NCDIJ
      WRITE(16,'(5E15.7)')AOMEGA,((LATT_CUR%ANORM(I)*1E-10),I=1,3),DYN%POTIM*1E-15
      WRITE(16,*) DYN%TEMP
      WRITE(16,*) ' CAR '
      WRITE(16,*) INFO%SZNAM1

!=======================================================================
!  write out initial header for EIGENVALUES
!=======================================================================
      WRITE(22,'(4I5)') T_INFO%NIONS,T_INFO%NIONS,DYN%NBLOCK*DYN%KBLOCK,WDES%ISPIN
      WRITE(22,'(5E15.7)') &
     &         AOMEGA,((LATT_CUR%ANORM(I)*1E-10_q),I=1,3),DYN%POTIM*1E-15_q
      WRITE(22,*) DYN%TEMP
      WRITE(22,*) ' CAR '
      WRITE(22,*) INFO%SZNAM1
      WRITE(22,'(3I7)') NINT(INFO%NELECT),KPOINTS%NKPTS,WDES%NB_TOT
      ENDIF
      io_end

      IF (IO%IU0>=0) &
      WRITE(TIU0,*)'POSCAR, INCAR and KPOINTS ok, starting setup'
!=======================================================================
! initialize the required grid structures
!=======================================================================
      CALL INILGRD(NGX,NGY,NGZ,GRID)
      CALL INILGRD(NGX,NGY,NGZ,GRID_SOFT)
      CALL INILGRD(NGXC,NGYC,NGZC,GRIDC)
      CALL INILGRD(GRIDUS%NGPTAR(1),GRIDUS%NGPTAR(2),GRIDUS%NGPTAR(3),GRIDUS)
#ifdef MPI
! only wavefunction grid uses local communication
      GRID%COMM     =>COMM_INB
! all other grids use world wide communication at the moment set their
! communication boards to the world wide communicator

!PK Charge density grids are replicated per set of distributed k-points
      GRID_SOFT%COMM=>COMM_KIN
      GRIDC%COMM    =>COMM_KIN
      GRIDUS%COMM   =>COMM_KIN
      GRIDB%COMM    =>COMM_KIN

      GRIDC%COMM_KIN    =>COMM_KIN
      GRIDC%COMM_KINTER =>COMM_KINTER
#endif
#ifdef usgrid
      CALL GEN_RC_GRID(GRIDUS)
      CALL GEN_RC_SUB_GRID(GRIDC,GRIDUS, C_TO_US, .TRUE.,.TRUE.)
#else
      CALL GEN_RC_GRID(GRIDC)
#endif
      CALL GEN_RC_SUB_GRID(GRID_SOFT, GRIDC, SOFT_TO_C, .TRUE.,.TRUE.)
      CALL GEN_RC_GRID(GRIDUS)
!=======================================================================
!  allocate work arrays
!=======================================================================
!
! GEN_LAYOUT determines the data layout (distribution) of the columns on parallel
! computers and allocates all required arrays of the WDES descriptor
      IF (INFO%ISTART==1) THEN
         CALL GEN_LAYOUT(GRID,WDES, LATT_CUR%B,LATT_CUR%B,IO%IU6,.TRUE.)
         CALL GEN_INDEX(GRID,WDES, LATT_CUR%B,LATT_CUR%B,IO%IU6,IO%IU0,.TRUE.)
! all other cases use LATT_INI for setup of GENSP
      ELSE
         DWRITE0 'call to genlay'
         CALL GEN_LAYOUT(GRID,WDES, LATT_CUR%B,LATT_INI%B,IO%IU6,.TRUE.)
         DWRITE0 'call to genind'
         CALL GEN_INDEX(GRID,WDES, LATT_CUR%B,LATT_INI%B,IO%IU6,IO%IU0,.TRUE.)
      ENDIF
#ifdef CUDA_GPU
      CALL GPU_INIT(WDES)
#endif
      CALL SET_NBLK_NSTRIP( WDES)
!
! non local projection operators
!
      CALL NONL_ALLOC(NONL_S,T_INFO,P,WDES,INFO%LREAL)
      CALL NONLR_SETUP(NONLR_S,T_INFO,P,INFO%LREAL,WDES%LSPIRAL)
!  optimize grid for real space representation and calculate IRMAX, IRALLOC
      NONLR_S%IRMAX=0 ; NONLR_S%IRALLOC=0
      CALL REAL_OPTLAY(GRID,LATT_CUR,NONLR_S,LPLANE_WISE,LREALLOCATE, IO%IU6, IO%IU0)
! allign GRID_SOFT with GRID in real space
      CALL SET_RL_GRID(GRID_SOFT,GRID)
! allocate real space projectors
      CALL NONLR_ALLOC(NONLR_S)
!  init FFT
      CALL FFTINI(WDES%NINDPW(1,1),WDES%NGVECTOR(1),KPOINTS%NKPTS,WDES%NGDIM,GRID)
#ifdef MPI
      CALL MAPSET(GRID)   ! generate the communication maps (patterns) for FFT
      IF (IO%IU6 >=0) THEN
         IF (GRID%RL%NFAST==1) THEN
            WRITE(TIU6,"(/' serial   3D FFT for wavefunctions')")
         ELSE IF (GRID%IN_RL%LOCAL  ) THEN
            WRITE(TIU6,"(/' parallel 3D FFT for wavefunctions:'/'    minimum data exchange during FFTs selected (reduces bandwidth)')")
         ENDIF
      ENDIF
      DWRITE0  'mapset aug done'
      CALL MAPSET(GRID_SOFT)
      DWRITE0  'mapset soft done'
      CALL MAPSET(GRIDC)
      IF (GRIDC%IN_RL%LOCAL .AND. IO%IU6 >=0) THEN
         WRITE(TIU6,"(' parallel 3D FFT for charge:'/'    minimum data exchange during FFTs selected (reduces bandwidth)'/)")
      ENDIF
      DWRITE0  'mapset wave done'
#ifdef usgrid
      CALL MAPSET(GRIDUS)
#endif
#endif

! allocate all other arrays
!
      ISP     =WDES%ISPIN
      NCDIJ   =WDES%NCDIJ
      NTYP    =T_INFO%NTYP
      NIOND   =T_INFO%NIOND
      NIOND_LOC=WDES%NIONS
      LMDIM   =P(1)%LMDIM

#if defined(CUDA_GPU) && defined(USE_PINNED_MEMORY)
      CALL nvpinnedmalloc(SV_PTR,DIMREAL(GRID%MPLWV)*NCDIJ*int(c_sizeof(fake),c_size_t))
      CALL c_f_pointer(SV_PTR,SV,(/DIMREAL(GRID%MPLWV),NCDIJ/))
! charges, potentials and so on
      CALL nvpinnedmalloc(CHTOT_PTR, GRIDC%MPLWV * NCDIJ*int(c_sizeof(fakec),c_size_t))
      CALL c_f_pointer(CHTOT_PTR, CHTOT, (/GRIDC%MPLWV,NCDIJ/))
#else
      ALLOCATE(CHTOT(GRIDC%MPLWV,NCDIJ),SV(DIMREAL(GRID%MPLWV),NCDIJ))
#endif
      ALLOCATE(CHTOTL(GRIDC%MPLWV,NCDIJ),DENCOR(GRIDC%RL%NP), &
               CVTOT(GRIDC%MPLWV,NCDIJ),CSTRF(GRIDC%MPLWV,NTYP), &
! small grid quantities
               CHDEN(GRID_SOFT%MPLWV,NCDIJ), &
! non local things
               CDIJ(LMDIM,LMDIM,NIOND_LOC,NCDIJ), &
               CQIJ(LMDIM,LMDIM,NIOND_LOC,NCDIJ), &
               CRHODE(LMDIM,LMDIM,NIOND_LOC,NCDIJ), &
! forces (depend on NIOND)
               EWIFOR(3,NIOND),TIFOR(3,NIOND), &
! dos
               DOS(NEDOS,NCDIJ),DOSI(NEDOS,NCDIJ), &
               DDOS(NEDOS,NCDIJ),DDOSI(NEDOS,NCDIJ), &
               PAR(1,1,1,1,NCDIJ),DOSPAR(1,1,1,NCDIJ), &
! paco
               PACO%SIPACO(0:PACO%NPACO))

      CALL REGISTER_ALLOCATE(16._q*(SIZE(CHTOT)+SIZE(CHTOTL)+SIZE(DENCOR)+SIZE(CVTOT)+SIZE(CSTRF)+SIZE(CHDEN)+SIZE(SV)), "grid")
      IF (WDES%LGAMMA) THEN
         CALL REGISTER_ALLOCATE(8._q*(SIZE(CDIJ)+SIZE(CQIJ)+SIZE(CRHODE)), "one-center")
      ELSE
         CALL REGISTER_ALLOCATE(16._q*(SIZE(CDIJ)+SIZE(CQIJ)+SIZE(CRHODE)), "one-center")
      ENDIF

      CALL ALLOCATE_AVEC(HAMILTONIAN%AVEC, HAMILTONIAN%AVTOT, GRID, GRIDC)

      CALL ALLOCATE_MU(HAMILTONIAN%MU, HAMILTONIAN%MUTOT, GRID, GRIDC, WDES)
! test_
!     CALL GENERATE_TAU_HANDLE(KINEDEN, GRIDC, WDES%NCDIJ)
      CALL ALLOCATE_TAU(KINEDEN%TAU,KINEDEN%TAUC,KINEDEN%TAUL,GRIDC,WDES%NCDIJ)
! test_
      CALL CREATE_CMBJ_AUX(GRIDC,T_INFO,LATT_CUR)

      DWRITE0 'allocation done'
!
      ALLOCATE(CWORK1(GRID%MPLWV))
      IF (IO%IU0>=0) WRITE(TIU0,*)'FFT: planning ...'
      CALL INIDAT(GRID%RC%NP,CWORK1)
      CALL FFTMAKEPLAN(CWORK1(1),GRID)
      DEALLOCATE(CWORK1)


      MPLMAX=MAX(GRIDC%MPLWV,GRID_SOFT%MPLWV,GRID%MPLWV)
#ifdef MPI
! give T3D opportunity to allocate all required shmem workspace
      CALL SHM_MAX(WDES, MPLMAX, MALLOC)
      CALL SHM_ALLOC(MALLOC)
#endif
      MIX%NEIG=0
! calculate required numbers elements which must be mixed in PAW
! set table for Clebsch-Gordan coefficients, maximum L is 2*3 (f states)
      LMAX_TABLE=6;  CALL YLM3ST_(LMAX_TABLE)

      N_MIX_PAW=0
      CALL SET_RHO_PAW_ELEMENTS(WDES, P , T_INFO, INFO%LOVERL, N_MIX_PAW )
      ALLOCATE( RHOLM(N_MIX_PAW,WDES%NCDIJ), RHOLM_LAST(N_MIX_PAW,WDES%NCDIJ))
! solve pseudo atomic problem
      IF (WANPROJ().OR.LPRJ_LCAO()) CALL LCAO_INIT(P,INFO,IO%IU0,IO%IU6)
! setup fock (requires YLM)
      CALL SETUP_FOCK(T_INFO, P, WDES, GRID, LATT_CUR, LMDIM,  INFO%SZPREC, IO%IU6, IO%IU0 )
! setup atomic PAW terms
      IF (LRHFATM()) THEN
         CALL SET_PAW_ATOM_POT_RHF(P,T_INFO,INFO%LOVERL,INFO%EALLAT,IO)
      ELSE
        CALL RHFATM_CROP_PSEUDO(P,T_INFO,INFO%LOVERL,IO)
        CALL SET_PAW_ATOM_POT(P , T_INFO, INFO%LOVERL,  &
             LMDIM, INFO%EALLAT, INFO%LMETAGGA, IO%IU6  )
      ENDIF
! possibly setup of the scf determination of the subspace rotation
      CALL SETUP_SUBROT_SCF(INFO,WDES,LATT_CUR,GRID,GRIDC,GRID_SOFT,SOFT_TO_C,IO%IU0,IO%IU5,IO%IU6)
! setup pead
      CALL PEAD_SETUP(WDES,GRID,GRIDC,GRIDUS,C_TO_US,KPOINTS,LATT_CUR,LATT_INI,T_INFO,P,LMDIM,INFO%LOVERL,IRDMAX,IO)
!=======================================================================
! allocate wavefunctions
!=======================================================================
      CALL ALLOCW(WDES,W,WUP,WDW)
      IF (INFO%LONESW) THEN
        CALL ALLOCW(WDES,W_F,WTMP,WTMP)
        CALL ALLOCW(WDES,W_G,WTMP,WTMP)
        ALLOCATE(CHAM(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN), &
                 CHF (WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN))
      ELSE
        ALLOCATE(CHAM(1,1,1,1), &
                 CHF (1,1,1,1))
      ENDIF
      CALL DUMP_ALLOCATE(IO%IU6)
!=======================================================================
! now read in wavefunctions
!=======================================================================

      W%CELTOT=0
      CALL START_TIMING("G")
      IF (INFO%ISTART>0) THEN
         IF (IO%IU0>=0) WRITE(TIU0,*)'reading WAVECAR'
         IF (IO%LDOWNSAMPLE ) THEN
            CALL INWAV_DOWNFOLD(IO, WDES, W, SYMM, KPOINTS, GRID, LATT_CUR, LATT_INI, INFO%ISTART,  EFERMI )
         ELSE
            CALL INWAV_FAST(IO, WDES, W, GRID, LATT_CUR, LATT_INI, INFO%ISTART,  EFERMI )
         ENDIF
      ELSE
         IF (IO%IU0>=0) WRITE(TIU0,*)'WAVECAR not read'
      ENDIF
      CALL CLOSEWAV

      IF (WDES%LSPIRAL.AND.(INFO%ISTART>0)) CALL CLEANWAV(WDES,W,INFO%ENINI)
      CALL STOP_TIMING("G",IO%IU6,'INWAV')


      IF (INFO%ISTART/=2) LATT_INI=LATT_CUR
!=======================================================================
! At this very point everything has been read in
! and we are ready to write all important information
! to the xml file
!=======================================================================
      CALL XML_ATOMTYPES(T_INFO%NIONS, T_INFO%NTYP, T_INFO%NITYP, T_INFO%ITYP, P%ELEMENT, P%POMASS, P%ZVALF, P%SZNAMP )

      CALL XML_TAG("structure","initialpos")
      CALL XML_CRYSTAL(LATT_CUR%A, LATT_CUR%B, LATT_CUR%OMEGA)
      CALL XML_POSITIONS(T_INFO%NIONS, DYN%POSION)
      IF (T_INFO%LSDYN) CALL XML_LSDYN(T_INFO%NIONS,T_INFO%LSFOR(1,1))
      IF (DYN%IBRION<=0 .AND. DYN%NSW>0 ) CALL XML_VEL(T_INFO%NIONS, DYN%VEL)
      IF (T_INFO%LSDYN) CALL XML_NOSE(DYN%SMASS)
      CALL XML_CLOSE_TAG("structure")
!=======================================================================
! initialize index tables for broyden mixing
!=======================================================================
      IF (((MIX%IMIX==4).AND.(.NOT.INFO%LCHCON)).OR.DYN%IBRION==10) THEN
! Use a reduced mesh but only if using preconditioning ... :
         IF (EXXOEP>=1) THEN
            CALL BRGRID(GRIDC,GRIDB,MAX(INFO%ENMAX*2,ENCUTGW_IN_CHI()),IO%IU6,LATT_CUR%B)
         ELSE
            CALL BRGRID(GRIDC,GRIDB,INFO%ENMAX,IO%IU6,LATT_CUR%B)
         ENDIF
         CALL INILGRD(GRIDB%NGX,GRIDB%NGY,GRIDB%NGZ,GRIDB)
         CALL GEN_RC_SUB_GRID(GRIDB,GRIDC, B_TO_C, .FALSE.,.TRUE.)
      ENDIF
! calculate the structure-factor, initialize some arrays
      IF (INFO%TURBO==0) CALL STUFAK(GRIDC,T_INFO,CSTRF)
      CHTOT=0 ; CHDEN=0; CVTOT=0
!=======================================================================
! construct initial  charge density:  a bit of heuristic is used
!  to get sensible defaults if the user specifies stupid values in the
!  INCAR files
! for the initial charge density there are several possibilties
! if INFO%ICHARG= 1 read in charge-density from a File
! if INFO%ICHARG= 2-3 construct atomic charge-densities of overlapping atoms
! if INFO%ICHARG= 4 read potential from file
! if INFO%ICHARG= 5 read GAMMA file and update orbitals accordingly
! if INFO%ICHARG >=10 keep chargedensity constant
!
!=======================================================================
      ! subtract 10 from ICHARG (10 means fixed charge density)
      IF (INFO%ICHARG>10) THEN
        INFO%INICHG= INFO%ICHARG-10
      ELSE
        INFO%INICHG= INFO%ICHARG
      ENDIF

      IF (INFO%INICHG==5) THEN
         ! no convergence correction for forces (reread GAMMA in last electronic step)
         IF (.NOT. MIX%MIXFIRST) INFO%LCORR=.FALSE.
         CALL CREATE_VASP_LOCK(COMM)
         ! and reset INICHG to 1 so that CHGCAR is read as well
         INFO%INICHG=1
      ENDIF

      ! then initialize CRHODE and than RHOLM (PAW related occupancies)
      CALL DEPATO(WDES, LMDIM, CRHODE, INFO%LOVERL, P, T_INFO)
      CALL SET_RHO_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
           CRHODE, RHOLM)

! MM: Stuff for Janos Angyan: write the AE-charge density for the
      IF (LWRT_AECHG()) THEN
      ! overlapping atomic charges to AECCAR1
         ! add overlapping atomic charges on dense regular grid
         CALL RHOATO_WORK(.TRUE.,.FALSE.,GRIDC,T_INFO,LATT_CUR%B,P,CSTRF,CHTOT)
         ! add (n_ae - n_ps - n_comp) as defined on radial grid
         CALL AUGCHG(WDES,GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
        &             LATT_CUR,P,T_INFO,SYMM, INFO%LOVERL, SOFT_TO_C,&
        &              LMDIM,CRHODE, CHTOT,CHDEN, IRDMAX,.FALSE.,.FALSE.)

#ifdef MPI
         IF (WDES%COMM_KINTER%NODE_ME.EQ.1) THEN
#endif
         io_begin
         ! write AECCAR1
         OPEN(UNIT=99,FILE=DIR_APP(1:DIR_LEN)//'AECCAR1',STATUS='UNKNOWN')
         ! write header
         CALL OUTPOS(99,.FALSE.,INFO%SZNAM1,T_INFO,LATT_CUR%SCALE,LATT_CUR%A, &
        &             .FALSE., DYN%POSION)
         io_end
         ! write AE charge density
         CALL OUTCHG(GRIDC,99,.TRUE.,CHTOT)
         io_begin
         CLOSE(99)
         io_end
         ! reset charge densities to zero again
         CHTOT=0; CHDEN=0
      ! AE core density to AECCAR0
         ! add partial core density on dense regular grid
         IF (INFO%LCORE) THEN
            CALL RHOATO_WORK(.TRUE.,.TRUE.,GRIDC,T_INFO,LATT_CUR%B,P,CSTRF,CHTOT)
         ENDIF
         ! add core densities (nc_ae-nc_ps) as defined on radial grid
         CALL AUGCHG(WDES,GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
        &             LATT_CUR,P,T_INFO,SYMM, INFO%LOVERL, SOFT_TO_C,&
        &              LMDIM,CRHODE, CHTOT,CHDEN, IRDMAX,.FALSE.,.TRUE.)
         io_begin
         ! write AECCAR0
         OPEN(UNIT=99,FILE=DIR_APP(1:DIR_LEN)//'AECCAR0',STATUS='UNKNOWN')
         ! write header
         CALL OUTPOS(99,.FALSE.,INFO%SZNAM1,T_INFO,LATT_CUR%SCALE,LATT_CUR%A, &
        &             .FALSE., DYN%POSION)
         io_end
         ! write AE core charge density
         CALL OUTCHG(GRIDC,99,.TRUE.,CHTOT)
         io_begin
         CLOSE(99)
         io_end
         ! and reset charge densities to zero again

#ifdef MPI
         END IF
#endif
         CHTOT=0; CHDEN=0
      ENDIF

     ! initial set of wavefunctions from diagonalization of Hamiltonian
     ! set INICHG to 2

      IF (INFO%INICHG==0 .AND. INFO%INIWAV==2) THEN
        INFO%INICHG=2
        IF (IO%IU6>=0) &
        WRITE(TIU6,*)'WARNING: no initial charge-density supplied,', &
                       ' atomic charge-density will be used'
      ENDIF

      IF (INFO%INICHG==1 .OR.INFO%INICHG==2 .OR.INFO%INICHG==3) THEN
         IF (IO%IU6>=0) WRITE(TIU6,*)'initial charge density was supplied:'
      ENDIF

      IF (INFO%INICHG==1) THEN

!PK Only first k-point group reads charge density
!PK Broadcast results outside of READCH due to awkward error logic in READCH
#ifdef MPI
         IF (COMM_KINTER%NODE_ME.EQ.1 ) THEN
#endif
           CALL READCH(GRIDC, INFO%LOVERL, T_INFO, CHTOT, RHOLM, INFO%INICHG, WDES%NCDIJ, &
             LATT_CUR, P, CSTRF(1,1), 18, IO%IU0)
#ifdef MPI
         END IF
#endif
        CALLMPI( M_bcast_i( COMM_KINTER, INFO%INICHG, 1))
        IF (INFO%INICHG.NE.0) THEN
           CALLMPI( M_bcast_z( COMM_KINTER, CHTOT, GRIDC%MPLWV*WDES%NCDIJ))
           CALLMPI( M_bcast_d( COMM_KINTER, RHOLM, SIZE(RHOLM)) )
        END IF

        IF (INFO%ICHARG>10 .AND. INFO%INICHG==0) THEN
           WRITE(*,*)'ERROR: charge density could not be read from file CHGCAR', &
               ' for ICHARG>10'
           STOP
        ENDIF
        ! error on reading CHGCAR, set INFO%INICHG to 2
        IF (INFO%INICHG==0)  INFO%INICHG=2
        ! no magnetization density set it according to MAGMOM
        IF (INFO%INICHG==-1) THEN
           IF (WDES%NCDIJ>1) &
                CALL MRHOATO(.FALSE.,GRIDC,T_INFO,LATT_CUR%B,P,CHTOT(1,2),WDES%NCDIJ-1)
           INFO%INICHG=1
           IF( TIU0 >=0) WRITE(TIU0,*)'magnetization density of overlapping atoms calculated'
        ENDIF

      ENDIF

      IF (INFO%INICHG==1) THEN
      ELSE IF (INFO%INICHG==2 .OR.INFO%INICHG==3) THEN
         IF(INFO%TURBO==0)THEN
            CALL RHOATO_WORK(.FALSE.,.FALSE.,GRIDC,T_INFO,LATT_CUR%B,P,CSTRF,CHTOT)
         ELSE
            CALL RHOATO_PARTICLE_MESH(.FALSE.,.FALSE.,GRIDC,LATT_CUR,T_INFO,INFO,P,CHTOT,IO%IU6)
         ENDIF
         IF (WDES%NCDIJ>1) &
              CALL MRHOATO(.FALSE.,GRIDC,T_INFO,LATT_CUR%B,P,CHTOT(1,2),WDES%NCDIJ-1)

         IF (IO%IU6>=0) WRITE(TIU6,*)'charge density of overlapping atoms calculated'
      ELSE IF (INFO%INICHG==4) THEN
         IF (IO%IU6>=0) WRITE(TIU6,*)'potential read from file POT'
         INFO%INICHG=4
      ELSE
         INFO%INICHG=0
      ENDIF

      IF (INFO%INICHG==1 .OR.INFO%INICHG==2 .OR.INFO%INICHG==3) THEN
         DO I=1,WDES%NCDIJ
            RHOTOT(I) =RHO0(GRIDC, CHTOT(1,I))
         ENDDO
         IF(IO%IU6>=0)  WRITE(TIU6,200) RHOTOT(1:WDES%NCDIJ)
 200     FORMAT(' number of electron ',F15.7,' magnetization ',3F15.7)
      ENDIF

      ! set the partial core density
      DENCOR=0
      IF (INFO%LCORE) CALL RHOPAR(GRIDC,T_INFO,INFO,LATT_CUR,P,CSTRF,DENCOR,IO%IU6)
      IF (LDO_METAGGA()) CALL TAUPAR(GRIDC,T_INFO,LATT_CUR%B,LATT_CUR%OMEGA,P,CSTRF,KINEDEN%TAUC)

      IF (INFO%INIWAV==2) THEN
         IF (IO%IU0>=0) &
         WRITE(TIU0,*) 'ERROR: this version does not support INIWAV=2'
         STOP
      ENDIF

      IF (IO%IU6>=0) THEN
        IF (INFO%INICHG==0 .OR. INFO%INICHG==4 .OR. (.NOT.INFO%LCHCOS.AND. INFO%NELMDL==0) ) THEN
          WRITE(TIU6,*)'charge density for first step will be calculated', &
           ' from the start-wavefunctions'
        ELSE
          WRITE(TIU6,*)'keeping initial charge density in first step'
        ENDIF
        WRITE(TIU6,130)
      ENDIF

      ! add Gaussian "charge-transfer" charges, if required
      CALL RHOADD_GAUSSIANS(T_INFO,LATT_CUR,P,GRIDC,NCDIJ,CHTOT,CSTRF)
      CALL RHOADD_GAUSSIANS_LIST(LATT_CUR,GRIDC,NCDIJ,CHTOT)

      DWRITE0 'atomic charge done'
!========================subroutine SPHER ==============================
! RSPHER calculates the real space projection operators
!    (VOLUME)^.5 Y(L,M)  VNLR(L) EXP(-i r k)
! subroutine SPHER calculates the nonlocal pseudopotential
! multiplied by the spherical harmonics and (1/VOLUME)^.5:
!    1/(VOLUME)^.5 Y(L,M)  VNL(L)
! (routine must be called if the size of the unit cell is changed)
!=======================================================================
      IF (INFO%LREAL) THEN
         CALL RSPHER(GRID,NONLR_S,LATT_CUR)

         INDMAX=0
         DO NI=1,T_INFO%NIONS
            INDMAX=MAX(INDMAX,NONLR_S%NLIMAX(NI))
         ENDDO
         IF (IO%IU6>=0) &
         WRITE(TIU6,*)'Maximum index for non-local projection operator ',INDMAX
      ELSE
         CALL SPHER(GRID, NONL_S, P, WDES, LATT_CUR,  1)
         CALL PHASE(WDES,NONL_S,0)
      ENDIF

      DWRITE0 'non local setup done'
!=======================================================================
! set the coefficients of the plane wave basis states to zero before
! initialising the wavefunctions by calling WFINIT (usually random)
!=======================================================================
      IF (INFO%ISTART<=0) THEN
        W%CPTWFP=0
        CALL WFINIT(WDES, W, INFO%ENINI)
      IF (INFO%INIWAV==1 .AND. INFO%NELMDL==0 .AND. INFO%INICHG/=0) THEN
        IF (IO%IU0>=0) &
        WRITE(TIU0,*) 'WARNING: random wavefunctions but no delay for ', &
                        'mixing, default for NELMDL'
        INFO%NELMDL=-5
        IF (INFO%LRMM .AND. .NOT. INFO%LDAVID) INFO%NELMDL=-12
      ENDIF

      ! initialize the occupancies, in the spin polarized case
      ! we try to set up and down occupancies individually
      ELEKTR=INFO%NELECT
      CALL DENINI(W%FERTOT(1,1,1),WDES%NB_TOT,KPOINTS%NKPTS,ELEKTR,WDES%LNONCOLLINEAR)

      IF (WDES%ISPIN==2) THEN
        RHOMAG=RHO0(GRIDC,CHTOT(1,2))
        ELEKTR=INFO%NELECT+RHOMAG
        CALL DENINI(W%FERTOT(1,1,1),WDES%NB_TOT,KPOINTS%NKPTS,ELEKTR,WDES%LNONCOLLINEAR)
        ELEKTR=INFO%NELECT-RHOMAG
        CALL DENINI(W%FERTOT(1,1,2),WDES%NB_TOT,KPOINTS%NKPTS,ELEKTR,WDES%LNONCOLLINEAR)
      ENDIF

      ENDIF
      DWRITE0 'wavefunctions initialized'
!=======================================================================
! INFO%LONESW initialize W_F%CELEN fermi-weights and augmentation charge
!=======================================================================
      E%EENTROPY=0

      IF ((INFO%LONESW .AND. KPOINTS%ISMEAR/=-2) .OR. (INFO%IALGO==3 .AND. KPOINTS%ISMEAR/=-2)) THEN
         E%EENTROPY=0
         IF (INFO%LONESW) W_F%CELTOT = W%CELTOT
         CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
               INFO%NUP_DOWN,  E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
               NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
      ENDIF
!=======================================================================
!  read Fermi-weigths from INCAR if supplied
!=======================================================================
      OPEN(UNIT=IO%IU5,FILE=INCAR,STATUS='OLD')
      CALL RDATAB(.FALSE.,INCAR,IO%IU5,'FERWE','=','#',';','F', &
     &     IDUM,W%FERTOT(1,1,1),CDUM,LDUM,CHARAC,N,KPOINTS%NKPTS*WDES%NB_TOT,IERR)

      IF ( ((IERR/=0) .AND. (IERR/=3)) .OR. &
     &     ((IERR==0).AND.(N<(KPOINTS%NKPTS*WDES%NB_TOT)))) THEN
         IF (IO%IU0>=0) &
         WRITE(TIU0,*)'Error reading item ''FERWE'' from file INCAR.'
         STOP
      ENDIF
! attention this feature is not supported by the xml writer
!      CALL XML_INCAR_V('FERWE','F',IDUM,W%FERTOT(1,1,1),CDUM,LDUM,CHARAC,N)

      IF (WDES%ISPIN==2) THEN
         CALL RDATAB(.FALSE.,INCAR,IO%IU5,'FERDO','=','#',';','F', &
     &        IDUM,W%FERTOT(1,1,INFO%ISPIN),CDUM,LDUM,CHARAC,N,KPOINTS%NKPTS*WDES%NB_TOT,IERR)
         IF ( ((IERR/=0) .AND. (IERR/=3)) .OR. &
     &        ((IERR==0).AND.(N<(KPOINTS%NKPTS*WDES%NB_TOT)))) THEN
            IF (IO%IU0>=0) &
            WRITE(TIU0,*)'Error reading item ''FERDO'' from file INCAR.'
            STOP
         ENDIF
! attention this feature is not supported by the xml writer
!         CALL XML_INCAR_V('FERDO','F',IDUM,W%FERTOT(1,1,INFO%ISPIN),CDUM,LDUM,CHARAC,N)
      ENDIF
      CLOSE(IO%IU5)

! if ISMEAR == -2 occupancies will be kept fixed
      IF (KPOINTS%ISMEAR==-2) THEN
         KPOINTS%SIGMA=-ABS(KPOINTS%SIGMA)
      ENDIF

!=======================================================================
! write out STM file if wavefunctions exist
!=======================================================================
       IF (STM(1) < STM(2) .AND. STM(3) /= 0 .AND. STM(4) < 0 .AND. INFO%ISTART == 1) THEN
        IF (STM(7) == 0) THEN
          CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
               INFO%NUP_DOWN,  E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
               NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
        ELSEIF (STM(7) /= 0) THEN
          EFERMI = STM(7)
          do_io WRITE(*,*)'-------------  VASP_IETS  --------------'
          do_io write(*,*) 'Efermi FROM INPUT= ',EFERMI
          do_io WRITE(*,*)'STM(7)=',STM(7)
          do_io WRITE(*,*)'-------------  VASP_IETS  --------------'

        ELSE
          do_io WRITE(*,*)'-------------  VASP_IETS  --------------'
          do_io WRITE(*,*)'ERROR WITH STM(7) STOP NOW'
          do_io WRITE(*,*)'STM(7)=',STM(7)
          STOP
        ENDIF
         do_io WRITE(*,*) 'Writing STM file'
         CALL  WRT_STM_FILE(GRID, WDES, WUP, WDW, EFERMI, LATT_CUR, &
               STM, T_INFO)
         do_io WRITE(*,*) "STM File written, ..."
         STOP
       ENDIF

!=======================================================================
! calculate the projections of the wavefunctions onto the projection
! operators using real-space projection scheme or reciprocal scheme
! then perform an orthogonalisation of the wavefunctions
!=======================================================================
      ! first call SETDIJ to set the array CQIJ
      CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
                  LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

      CALL CALC_JPAW_HAMIL(T_INFO, P)
      ! set vector potential
      CALL VECTORPOT(GRID, GRIDC, GRID_SOFT, SOFT_TO_C,  WDES%COMM_INTER, &
                 LATT_CUR, T_INFO%POSION, HAMILTONIAN%AVEC, HAMILTONIAN%AVTOT)
      ! first call to SETDIJ_AVEC only sets phase twisted projectors
      CALL SETDIJ_AVEC(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
                 LMDIM, CDIJ, HAMILTONIAN%AVTOT, NONLR_S, NONL_S, IRDMAX )

      DWRITE0 'setdij done'

      IF (IO%IU6>=0) &
      WRITE(TIU6,*)'Maximum index for augmentation-charges ', &
                     IRDMAA,'(set IRDMAX)'
      CALL WVREAL_PRECISE(W)

      CALL PROALL (GRID,LATT_CUR,NONLR_S,NONL_S,W)
      DWRITE0 'proall done'

      CALL ORTHCH(WDES,W, INFO%LOVERL, LMDIM,CQIJ)
      CALL REDIS_PW_OVER_BANDS(WDES, W)
      DWRITE0 'orthch done'

      IF (IO%IU6>=0) WRITE(TIU6,130)

      ! strictly required for HF, since SET_CHARGE fails for ISYM=3
      CALL KPAR_SYNC_ALL(WDES,W)
!=======================================================================
! partial band decomposed chargedensities PARDENS
!=======================================================================
      IF (LPARD) THEN
           CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
                INFO%NUP_DOWN,  E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
                NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
           CALL PARCHG(W,WUP,WDW,WDES,CHDEN,CHTOT,CRHODE,INFO,GRID, &
                GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
                LATT_CUR,P,T_INFO,SOFT_TO_C,SYMM,IO, &
                DYN,EFERMI,LMDIM,IRDMAX,NIOND)
      ENDIF
!=======================================================================
! calculate initial chargedensity if
! ) we do not have any chargedensity
! ) we use a selfconsistent minimization scheme and do not have a delay
!=======================================================================
      IF (INFO%INICHG==0 .OR. INFO%INICHG==4 .OR. (.NOT.INFO%LCHCOS .AND. INFO%NELMDL==0)  ) THEN
         IF (IO%IU6>=0) &
              WRITE(TIU6,*)'initial charge from wavefunction'

         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'initial charge from wavefunction'

         !IF (IO%IU6>=0) THEN
         !  write(*,*) KPOINTS_FULL%NKPX,KPOINTS_FULL%NKPY,KPOINTS_FULL%NKPZ
         !  write(*,*) WDES%NKPTS
         !  DO NK=1,WDES%NKPTS
         !     write(*,*) WDES%VKPT(1,NK), WDES%VKPT(2,NK), WDES%VKPT(3,NK), WDES%WTKPT(NK)
         !  ENDDO
         !ENDIF
         !write(*,*) COMM%MPI_COMM, MPI_sum
         inquire(file='dmft_params.dat',exist=LDMFT)
         if (LDMFT.eq..TRUE.) then
            IF (IO%IU6>=0) &
                WRITE(TIU6,*)'DFT+DMFT run with DMFTwDFT'

            IF (IO%IU0>=0) &
                WRITE(TIU0,*)'DFT+DMFT run with DMFTwDFT'


           CALL ParseWannier90(num_wannier)


           CALL CHECK_FULL_KPOINTS
           !num_wann_=14
           num_wann_ = num_wannier
           !allocate(DMFT_eval_(28,KPOINTS_FULL%NKPTS))
           !allocate(DMFT_evec_(28,28,KPOINTS_FULL%NKPTS))
           !allocate(band_win_(2,KPOINTS_FULL%NKPTS))
           allocate(DMFT_eval_(num_wannier,WDES%NKPTS))
           allocate(DMFT_evec_(num_wannier,num_wannier,WDES%NKPTS))
           allocate(band_win_(2,WDES%NKPTS))
           !CALL MLWF_WANNIER90_AMN(WDES,W,P,CQIJ,T_INFO,LATT_CUR,INFO,IO)
           !CALL MPI_barrier( WDES%COMM%MPI_COMM, ierr )
           !write(*,*) 'hi'
           !CALL Compute_DMFT(KPOINTS_FULL%NKPTS,num_wann_,KPOINTS_FULL%VKPT,KPOINTS_FULL%WTKPT,band_win_,DMFT_eval_,DMFT_evec_)
           !write(*,*) 'hi2'
           !write(*,*) "Computing DMFT..."
           CALL Compute_DMFT(WDES%NKPTS,num_wannier,WDES%VKPT,WDES%WTKPT,band_win_,DMFT_eval_,DMFT_evec_)

            ! Print the eigen values for the wannier manifold for each
            ! k-point

!            write(*,*) "RSPIN: "
!            write(*,*) W%WDES%RSPIN
            open(1, file = 'dmft-nkij_eig_vasp.dat',status="REPLACE")
            write(1,*) "k, wk, band, n_kij_eig"
            do eig_k=1,WDES%NKPTS
                do eig_i=1,num_wannier
                    ! Normalization
                    DMFT_eval_norm = W%WDES%RSPIN*WDES%WTKPT(eig_k)*DMFT_eval_(eig_i,eig_k)
                    write(1,"(i8,f10.4,i8,f10.4)") eig_k,WDES%WTKPT(eig_k),eig_i,DMFT_eval_norm
                enddo
            enddo
            close(1)




           !CALL Compute_DMFT_from_amn(WDES%NKPTS,num_wann_,WDES%VKPT,WDES%WTKPT,band_win_,DMFT_eval_,DMFT_evec_)
           !IF (IO%IU6>=0) then
           !  DO IK=1,WDES%NKPTS
           !    write(*,*) IK,DMFT_eval_(:,IK),DMFT_evec_(:,1,IK)
           !  ENDDO
           !ENDIF
           !STOP
         endif

!!!!!!!!!  !!!!!!!! Rotate WAVE function !!!!!!!!!!!!!!!!!!

!           ALLOCATE(FERWE(WDES%NBANDS,WDES%NKPTS))
!           ALLOCATE(WaveW(WDES%NRPLWV,WDES%NB_TOT,WDES%NKPTS))
!           ALLOCATE(WaveW_tot(WDES%NRPLWV,WDES%NB_TOT,WDES%NKPTS))
!           ALLOCATE(WProj(WDES%NPRO_TOT,WDES%NB_TOT,WDES%NKPTS))
!           ALLOCATE(WProj_tot(WDES%NPRO_TOT,WDES%NB_TOT,WDES%NKPTS))
!           FERWE=0;WaveW=0;WaveW_tot=0;WProj=0;WProj_tot=0
!           DO NK=1,WDES%NKPTS
!#ifdef MPI
!           IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
!#endif
!           DO NB=1,WDES%NBANDS
!             FERWE(NB,NK)=W%FERWE(NB,NK,1)
!             DO NP=1,WDES%NRPLWV
!                WaveW(NP,WDES%NB_LOW+(NB-1)*WDES%NB_PAR,NK)=W%CW(NP,NB,NK,1)!/WDES%COMM_KINTER%NCPU
!             ENDDO
!             DO NP=1,WDES%NPRO_TOT
!                WProj(NP,WDES%NB_LOW+(NB-1)*WDES%NB_PAR,NK)=W%CPROJ(NP,NB,NK,1)!/WDES%COMM_KINTER%NCPU
!             ENDDO
!           ENDDO
!           ENDDO
!
!!           write(*,*) 'hi'
!           CALL MPI_ALLREDUCE(WaveW,WaveW_tot,WDES%NRPLWV*WDES%NB_TOT*WDES%NKPTS,MPI_DOUBLE_COMPLEX,MPI_SUM,WDES%COMM%MPI_COMM,ierr)
!           CALL MPI_ALLREDUCE(WProj,WProj_tot,WDES%NPRO_TOT*WDES%NB_TOT*WDES%NKPTS,MPI_DOUBLE_COMPLEX,MPI_SUM,WDES%COMM%MPI_COMM,ierr)
!!           write(*,*) 'hi2'
!           DO NK=1,WDES%NKPTS
!#ifdef MPI
!             IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
!#endif
!             nbmin=band_win_(1,NK); nbmax=band_win_(2,NK)
!             num_band_max=nbmax-nbmin+1
!             DO NB=1,WDES%NBANDS
!               NBG=WDES%NB_LOW+(NB-1)*WDES%NB_PAR
!               if (NBG .GE. nbmin .AND. NBG .LE. nbmax) then
!                 W%FERWE(NB,NK,1)=0
!                 if (NBG-nbmin+1 .LE. num_wann_) W%FERWE(NB,NK,1)=DMFT_EVAL_(NBG-nbmin+1,NK)
!                 W%CW(:,NB,NK,1)=0
!                 W%CPROJ(:,NB,NK,1)=0
!                 DO NP=1,num_band_max
!                   W%CW(:,NB,NK,1)=W%CW(:,NB,NK,1)+WaveW_tot(:,nbmin+NP-1,NK)*DMFT_EVEC_(NP,NBG-nbmin+1,NK)
!                   W%CPROJ(:,NB,NK,1)=W%CPROJ(:,NB,NK,1)+WProj_tot(:,nbmin+NP-1,NK)*DMFT_EVEC_(NP,NBG-nbmin+1,NK)
!                 ENDDO
!               endif
!             ENDDO
!           ENDDO
!           DEALLOCATE(WaveW_tot)
!           DEALLOCATE(WProj_tot)
!         endif
         !tot_n=0
         !DO NK=1,WDES%NKPTS
         !  DO NB=1,WDES%NBANDS
         !    tot_n=tot_n+W%FERWE(NB,NK,1)*WDES%WTKPT(NK)*2
         !  ENDDO
         !ENDDO
         !write(*,*) WDES%NB_LOW, WDES%COMM_KINTER%NODE_ME,tot_n
         !STOP
         !DO NK=1,WDES%NKPTS
         !  DO NB=1,num_wann_
         !    tot_n=0
         !    DO NP=1,num_wann_
         !      tot_n=tot_n+DMFT_EVEC_(NB,NP,NK)*CONJG(DMFT_EVEC_(NB,NP,NK))
         !    ENDDO
         !    !write(*,*) NK,NB,tot_n
         !  ENDDO
         !ENDDO


         CALL SET_CHARGE(W, WDES, INFO%LOVERL, &
              GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
              LATT_CUR, P, SYMM, T_INFO, &
              CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

         CALL SET_KINEDEN(GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM, &
              T_INFO%NIONS,W,WDES,KINEDEN)

!         if (LDMFT.eq..TRUE.) then
!           DO NK=1,WDES%NKPTS
!#ifdef MPI
!             IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
!#endif
!             DO NB=1,WDES%NBANDS
!                W%FERWE(NB,NK,1)=FERWE(NB,NK)
!                NBG=WDES%NB_LOW+(NB-1)*WDES%NB_PAR
!                if (NBG .GE. nbmin .AND. NBG .LE. nbmax) then
!                  W%CW(:,NB,NK,1)=WaveW(:,NBG,NK)
!                  W%CPROJ(:,NB,NK,1)=WProj(:,NBG,NK)
!                endif
!             ENDDO
!           ENDDO
!           DEALLOCATE(FERWE)
!           DEALLOCATE(WaveW)
!           DEALLOCATE(WProj)
!         endif
      ENDIF

!=======================================================================
! initialise the predictor for the wavefunction with
! the first available ionic positions
! PRED%INIPRE=2 continue with existing file
! PRED%INIPRE=1 initialise
!=======================================================================
      IF (INFO%ISTART==3) THEN
        PRED%INIPRE=2
      ELSE
        PRED%INIPRE=1
      ENDIF

      IF (DYN%IBRION/=-1 .AND. PRED%IWAVPR >= 12 ) THEN
        CALL WAVPRE_NOIO(GRIDC,P,PRED,T_INFO,W,WDES,LATT_CUR,IO%LOPEN, &
           CHTOT,RHOLM,N_MIX_PAW, CSTRF, KINEDEN, LMDIM,CQIJ,INFO%LOVERL,IO%IU0)
      ELSE IF (DYN%IBRION/=-1 .AND. PRED%IWAVPR >= 2 .AND. PRED%IWAVPR <10 )  THEN
        CALL WAVPRE(GRIDC,P,PRED,T_INFO,W,WDES,LATT_CUR,IO%LOPEN, &
           CHTOT,RHOLM,N_MIX_PAW, CSTRF, LMDIM,CQIJ,INFO%LOVERL,IO%IU0)
      ENDIF
      DWRITE0 "wavpre is ok"

!=======================================================================
! if INFO%IALGO < 0 make some performance tests
!=======================================================================
      IF (INFO%IALGO<0) THEN
         CALL PERFORMANCE_TEST
         GOTO 5100
      ENDIF

!======================== SUBROUTINE FEWALD ============================
! calculate ewald energy forces, and stress
!=======================================================================
      CALL START_TIMING("G")
      IF(INFO%TURBO==0) THEN
         CALL FEWALD(DYN%POSION,EWIFOR,LATT_CUR%A,LATT_CUR%B,LATT_CUR%ANORM,LATT_CUR%BNORM, &
              LATT_CUR%OMEGA,EWSIF,E%TEWEN,T_INFO%NTYP,P%ZVALF,T_INFO%VCA, &
              T_INFO%NIONS,NIOND,T_INFO%ITYP,T_INFO%NITYP,IO%IU6,.TRUE.)
      ELSE
         ALLOCATE(CWORK1(GRIDC%MPLWV))
         CALL POTION_PARTICLE_MESH(GRIDC,P,LATT_CUR,T_INFO,CWORK1,E%PSCENC,E%TEWEN,EWIFOR)
         DEALLOCATE(CWORK1)
      ENDIF
      CALL STOP_TIMING("G",TIU6,'FEWALD')
!=======================================================================
! Set INFO%LPOTOK to false: this requires  a recalculation of the
! total local potential
!=======================================================================
      INFO%LPOTOK=.FALSE.
      TOTENG=0
      INFO%LSTOP=.FALSE.
      INFO%LSOFT=.FALSE.

      IF (INFO%INICHG==4) THEN

         IF (IO%LOPEN) OPEN(IO%IUVTOT,FILE='POT',STATUS='UNKNOWN')
         REWIND IO%IUVTOT

         CALL READPOT(GRIDC, T_INFO, CVTOT, INFO%INICHG, WDES, &
              LATT_CUR, P, LMDIM, CDIJ, N_MIX_PAW, IO%IUVTOT, IO%IU0)

        IF (INFO%INICHG==4) THEN
           DO ISP=1,WDES%NCDIJ
              CALL FFT_RC_SCALE(CVTOT(1,ISP),CVTOT(1,ISP),GRIDC)
           ENDDO

           CALL SET_SV( GRID, GRIDC, GRID_SOFT, WDES%COMM_INTER, SOFT_TO_C, WDES%NCDIJ, SV, CVTOT)
           INFO%LPOTOK=.TRUE.
        ELSE
           IF (IO%IU0>=0) THEN
              WRITE(IO%IU0,*)'reading POT failed can not continue'
           ENDIF

           STOP
        ENDIF
      ENDIF
!=======================================================================
! some special calculation routines such as GW, MP2, MP2NO, ...
!=======================================================================
      IF (LMP2 .OR. LMP2KPAR .OR. LMP2NO .OR. LFCIDUMP  .OR. LRPAX .OR.  LCCSD .OR. LBRACKETST) THEN
         CALL PSI_CORR_MAIN(P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI,  &
          & T_INFO, DYN, IO,KPOINTS,SYMM,GRID,INFO,MIX%AMIX,MIX%BMIX)
      ENDIF

      IF(FOURORBIT==1) THEN
         CALL TWOELECTRON4O_MAIN( &
              P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
              T_INFO,DYN,INFO,IO,KPOINTS,SYMM,GRID, LMDIM)
         CALL STOP_TIMING("G",IO%IU6,'TWOE4O')
      ENDIF
      IF (LMP2 .OR. LMP2KPAR .OR. LMP2NO .OR. LFCIDUMP  .OR. LRPAX .OR.  LCCSD .OR. &
        & LBRACKETST .OR. FOURORBIT==1 ) THEN
         STOP
      ENDIF

#ifdef VASP2WANNIER90
      IF (L2E4W) THEN
         CALL DMFT_MATRIX_ELEMENTS( &
        &   WDES,W,KPOINTS,GRID,T_INFO,INFO,P, &
        &   NONL_S,NONLR_S,SYMM,LATT_CUR,CQIJ,LMDIM,IO)
         ! stop execution of VASP (wave functions have been messed up so let's quit)
         GOTO 5100
      ENDIF
#endif
      IF (LCHI) THEN
         IF (.NOT. LUSEPEAD() .AND. .NOT. LNABLA .AND. .NOT. EXXOEP>=1 .AND. .NOT. USE_OEP_IN_GW()) THEN
         ! optical properties do not work for conventional GW
         ! because the potential is non-local
         ! OEP are fine however
            IO%LOPTICS=.FALSE.
         ENDIF
         IF (.NOT. LRPAFORCE) THEN
            CALL DO_RPA
            INFO%LSTOP=.TRUE.   ! stop ionic loop immediately
         ENDIF
      ENDIF
!***********************************************************************
!***********************************************************************
!
! ++++++++++++ do 1000 n=1,required number of timesteps ++++++++++++++++
!
! this is the main loop of the program during which one complete step of
! the electron dynamics and ionic movements is performed.
! NSTEP           loop counter for ionic movement
! N               loop counter for self-consistent loop
! INFO%NELM            number of electronic movement loops
!***********************************************************************
!***********************************************************************

      IF (IO%IU0>=0) WRITE(TIU0,*)'entering main loop'

      NSTEP = 0
      CALL START_TIMING("LOOP+")

      CALL CONSTRAINED_M_INIT(T_INFO,GRIDC,LATT_CUR)
      CALL INIT_WRITER(P,T_INFO,WDES)
!=======================================================================
      ion: DO
      IF (INFO%LSTOP) EXIT ion
!=======================================================================

!  reset broyden mixing
      MIX%LRESET=.TRUE.
io_begin
      IF (IO%LOPEN) CALL WFORCE(IO%IU6)            ! write of OUTCAR
      IF (IO%LOPEN .AND. NSTEP/=0) CALL WFORCE(17) ! write of OSZICAR
io_end
!=======================================================================
! initialize pair-correlation funtion to zero
! also set TMEAN und SMEAN to zero
! TMEAN/SMEAN  is the mean temperature
!=======================================================================
      IF (MOD(NSTEP,DYN%NBLOCK*DYN%KBLOCK)==0) THEN
        SMEANP=0
        PACO%SIPACO=0
        TMEAN=0
        TMEAN0=0
        SMEAN=0
        SMEAN0=0
        DDOSI=0._q
        DDOS =0._q
      ENDIF

      NSTEP = NSTEP + 1

    IF (INFO%NELM==0 .AND. IO%IU6>=0) WRITE(TIU6,140) NSTEP,0

!***********************************************************************
!***********************************************************************
! this part performs the total energy minimisation
! INFO%NELM loops are made, if the difference between two steps
! is less then INFO%EDIFF the loop is aborted
!***********************************************************************
!***********************************************************************
    IF (.NOT. LJ_ONLY .AND. .NOT. LCHI) THEN
      CALL XML_TAG("calculation")

      CALL ELECTRONIC_OPTIMIZATION

      ! possibly orbitals have been updated
      ! force them to be real again be calling WVREAL_PRECISE
      ! costs some extra work, but it is ultimately safer to do it here
      IF (W%WDES%LORBITALREAL .AND. (INFO%LDIAG .OR. INFO%LEXACT_DIAG)) THEN
        CALL WVREAL_PRECISE(W)
        CALL PROALL (GRID,LATT_CUR,NONLR_S,NONL_S,W)
        CALL STOP_TIMING("G",IO%IU6,"WVREAL",XMLTAG="wvreal")
        CALL ORTHCH(WDES,W, INFO%LOVERL, LMDIM,CQIJ)
        CALL STOP_TIMING("G",IO%IU6,"ORTHCH",XMLTAG="orth")
      ENDIF

      IF (NSTEP==1) THEN
         IF (.NOT.INFO%LONESW.AND.(LPEAD_CALC_EPS().OR.LPEAD_NONZERO_EFIELD())) THEN
            CALL ALLOCW(WDES,W_F,WTMP,WTMP)
            CALL ALLOCW(WDES,W_G,WTMP,WTMP)
            DEALLOCATE(CHAM, CHF)
            ALLOCATE(CHAM(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN), &
              CHF (WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN))
!           ! setup for scf subspace rotation
!           INFO%LONESW=.TRUE.
!           CALL SETUP_SUBROT_SCF(INFO,WDES,LATT_CUR,GRID,GRIDC,GRID_SOFT,SOFT_TO_C,IO%IU0,IO%IU5,IO%IU6)
!           INFO%LONESW=.FALSE.
         ENDIF

         CALL PEAD_ELMIN( &
          HAMILTONIAN,KINEDEN, &
          P,WDES,NONLR_S,NONL_S,W,W_F,W_G,LATT_CUR,LATT_INI, &
          T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI,CHF,CHAM,ECONV, &
          NSTEP,LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI,LDIMP,LMDIMP)

         IF (.NOT.INFO%LONESW.AND.(LPEAD_CALC_EPS().OR.LPEAD_NONZERO_EFIELD())) THEN
            CALL DEALLOCW(W_F)
            CALL DEALLOCW(W_G)
            DEALLOCATE(CHAM, CHF)
            ALLOCATE(CHAM(1,1,1,1),CHF (1,1,1,1))
         ENDIF
      ENDIF

#ifdef VASP2WANNIER90
      IF ((.NOT.(KINTER/=0)).AND.(.NOT.L2E4W) .AND. WANNIER90()) THEN
         CALL KPAR_SYNC_ALL(WDES,W)
         !switch off symmetry if UWAN files are written
         IF ((SYMM%ISYM>=0.AND. .NOT.WDES%LGAMMA) .AND. LWRITE_WANPROJ) THEN
            ! switch of symmetry
            CALL NOSYMM(LATT_CUR%A,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS,SYMM%PTRANS, &
                 SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%MAGROT,WDES%ISPIN,IO%IU6)
            ! reread k-points with LINVERSION=.FALSE. to generate full mesh
            CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR,.FALSE., &
                 T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IO%IU6,IO%IU0)
            CALL KPAR_SYNC_ALL(WDES,W)
            CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_CUR,-1, IO%IU0)
            CALL PEAD_RESETUP_WDES(WDES, GRID, KPOINTS, LATT_CUR, LATT_CUR, IO)
            CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
         ENDIF
         !LWANNIER_RUN=.TRUE.
         !WRITE_MMN_AMN=.FALSE.
         CALL MLWF_WANNIER90(WDES,W,P,CQIJ,T_INFO,LATT_CUR,INFO,IO)
         IF ((SYMM%ISYM>=0.AND. .NOT.WDES%LGAMMA) .AND. LWRITE_WANPROJ) THEN
            CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
                 T_INFO%LSDYN,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
                 SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
                 SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,-1)
            CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
                 SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
                 T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,-1,IO%IU0)
            CALL KPAR_SYNC_ALL(WDES,W)
            CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
            CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
         ENDIF
      ENDIF
#endif
      CALL WNPR_PROJECT(W,WDES,T_INFO,INFO,P,CQIJ,LATT_CUR,IO)
    ELSE IF (LCHI) THEN
       CALL DO_RPA
    END IF ! LJ_ONLY

!  'soft stop': stop after the next ionic step finished
!  in order to do this create file STOPCAR and set an entry INFO%LSTOP=.TRUE.
      LTMP=.FALSE.
      CALL RDATAB(IO%LOPEN,'STOPCAR',99,'LSTOP','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LTMP,CHARAC,NCOUNT,1,IERR)
      ITMP=0 ; IF (LTMP) ITMP=1 ; CALLMPI ( M_sum_i(W%WDES%COMM, ITMP, 1))
      IF (ITMP>0) THEN
        IF (IO%IU0>=0) &
             WRITE(TIU0,*) 'soft stop encountered!  aborting job ...'
        IF (IO%IU6>=0) &
        WRITE(TIU6,*) 'soft stop encountered!  aborting job ...'

        INFO%LSOFT=.TRUE.
      ENDIF

    IF (.NOT. LJ_ONLY) THEN
!=======================================================================
! do some check for occupation numbers (do we have enough bands]:
!=======================================================================
      IOCCUP=0
      IOCCVS=0
      SUMNEL=0._q
      DO ISP=1,WDES%ISPIN
      DO NN=1,KPOINTS%NKPTS
         IF (ABS(W%FERTOT(WDES%NB_TOT,NN,ISP))>1.E-2_q) THEN
            IOCCUP=IOCCUP+1
! total occupancy ('number of electrons in this band'):
            SUMNEL=SUMNEL+WDES%RSPIN*W%FERTOT(WDES%NB_TOT,NN,ISP)*KPOINTS%WTKPT(NN)
            IF (IOCCUP<NTUTOR) THEN
               ITUT(IOCCUP)=NN
               RTUT(IOCCUP)=W%FERTOT(WDES%NB_TOT,NN,ISP)
            ENDIF
         ENDIF
! count seperately 'seriously large occupations' ...
         IF (ABS(W%FERTOT(WDES%NB_TOT,NN,ISP))>2.E-1_q) IOCCVS=IOCCVS+1
      ENDDO
      ENDDO

      IF (((IOCCUP/=0).AND.(SUMNEL>1E-3_q)).OR.(IOCCVS/=0)) THEN
         IOCC=MIN(IOCCUP+1,NTUTOR)
         ITUT(IOCC)=WDES%NB_TOT
         RTUT(IOCC)=SUMNEL
         CALL VTUTOR('U','HIGHEST BANDS OCCUPIED', &
     &               RTUT,IOCC,ITUT,IOCC,CDUM,1,LDUM,1,IO%IU6,IO%IDIOT)
         CALL VTUTOR('U','HIGHEST BANDS OCCUPIED', &
     &               RTUT,1,ITUT,IOCC,CDUM,1,LDUM,1,IO%IU0,IO%IDIOT)
! no matter how the IO%IDIOT-flag is set: for seriously large occupations
! give always a message to the user because this is very important!!
         IF ((IO%IDIOT<=0).AND.((IOCCVS/=0).OR.(SUMNEL>1E-2_q))) THEN
            CALL VTUTOR('U','HIGHEST BANDS OCCUPIED', &
     &                  RTUT,IOCC,ITUT,IOCC,CDUM,1,LDUM,1,IO%IU6,1)
            CALL VTUTOR('U','HIGHEST BANDS OCCUPIED', &
     &                  RTUT,1,ITUT,IOCC,CDUM,1,LDUM,1,IO%IU0,1)
         ENDIF
      ELSE IF (IOCCUP/=0) THEN
! for less serious occupancies give just some 'good advice' ...
         IOCC=MIN(IOCCUP+1,NTUTOR)
         ITUT(IOCC)=WDES%NB_TOT
         RTUT(IOCC)=SUMNEL
         CALL VTUTOR('A','HIGHEST BANDS OCCUPIED', &
     &               RTUT,IOCC,ITUT,IOCC,CDUM,1,LDUM,1,IO%IU6,IO%IDIOT)
         CALL VTUTOR('A','HIGHEST BANDS OCCUPIED', &
     &               RTUT,1,ITUT,IOCC,CDUM,1,LDUM,1,IO%IU0,IO%IDIOT)
      ENDIF

! Okay when we arrive here first then we got it for the first ionic
! configuration and if more steps follow we might not need any further
! 'delay' switching on the selfconsistency? (if NELMDL<0 switch off!).
      IF (INFO%NELMDL<0) INFO%NELMDL=0

! calculation partial DOS
      IF (IO%LORBIT>=10) THEN
         CALL SPHPRO_FAST( &
          GRID,LATT_CUR, P,T_INFO,W, WDES, 71,IO%IU6,&
          INFO%LOVERL,LMDIM,CQIJ, LDIMP, LDIMP,LMDIMP,.FALSE., IO%LORBIT,PAR, &
          EFERMI, KPOINTS%EMIN, KPOINTS%EMAX)
      ENDIF
!***********************************************************************
!***********************************************************************
! Now perform the ion movements:
!***********************************************************************
!***********************************************************************

!====================== FORCES+STRESS ==================================
!
!=======================================================================
      NORDER=0
      IF (KPOINTS%ISMEAR>=0) NORDER=KPOINTS%ISMEAR

 7261 FORMAT(/ &
     &        A/ &
     &        '  ---------------------------------------------------'/ &
     &        '  free  energy   TOTEN  = ',F18.8,' eV'// &
     &        '  energy  without entropy=',F18.8, &
     &        '  energy(sigma->0) =',F18.8)

      ! no forces for OEP methods (EXXOEP/=0)
      ! no forces for IALGO=1-4
      ! no forces if potential was read in INFO%INICHG==4
      ! no forces if inner eigenvalue solver are used
      ! no forces for RPA (has already been done)
      IF ( EXXOEP==0 .AND.  (INFO%LONESW .OR. INFO%LDIAG .OR. INFO%IALGO>4) &
           .AND. INFO%INICHG/=4 .AND. INFO%IHARMONIC==0 .AND. (.NOT. LCHI)) THEN
        CALL FORCE_AND_STRESS( &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR, &
          T_INFO,T_INFO,DYN,INFO,IO,MIX,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV, &
          LMDIM, IRDMAX, (NSTEP==1 .OR.NSTEP==DYN%NSW).OR.IO%NWRITE>=1, &
          DYN%ISIF/=0, DYN%ISIF/=0,.TRUE.,  &
          XCSIF, EWSIF, TSIF, EWIFOR, TIFOR, PRESS, TOTEN, KPOINTS )
      ENDIF

!     trigger calculation of ensemble energies in next call to POTLOK
#ifdef libbeef
      IF(LBEEFENS) LBEEFCALCBASIS = .TRUE.
#endif

      IF (IO%IU6>=0) THEN
         WRITE(TIU6,130)
         WRITE(TIU6,7261) '  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)', &
            TOTEN,TOTEN-E%EENTROPY,TOTEN-E%EENTROPY/(2+NORDER)
      ENDIF

      IF (DYN%PSTRESS/=0) THEN
         TOTEN=TOTEN+DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA

         IF (IO%IU6>=0) &
              WRITE(TIU6,7264) TOTEN,DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA
 7264    FORMAT ('  enthalpy is  TOTEN    = ',F18.8,' eV   P V=',F18.8/)

      ELSE
         IF (IO%IU6>=0) WRITE(TIU6,*)
      ENDIF

      CALL XML_TAG("structure")
      CALL XML_CRYSTAL(LATT_CUR%A, LATT_CUR%B, LATT_CUR%OMEGA)
      CALL XML_POSITIONS(T_INFO%NIONS, DYN%POSION)
      CALL XML_CLOSE_TAG("structure")
      CALL XML_FORCES(T_INFO%NIONS, TIFOR)
      IF (DYN%ISIF/=0) THEN
         CALL XML_STRESS(TSIF*EVTOJ*1E22_q/LATT_CUR%OMEGA)
      ENDIF


! check the consistency of forces and total energy
      CALL CHECK(T_INFO%NIONS,DYN%POSION,TIFOR,EWIFOR,TOTEN,E%TEWEN,LATT_CUR%A,IO%IU6)
      IF (IO%IU6>=0) WRITE(IO%IU6,130)

!-----------------------------------------------------------------------
! we have maybe to update CVTOT here (which might have been destroyed
!  by the force routines)
!-----------------------------------------------------------------------
      IF (LVCADER) THEN
         ! derivative with respect to VCA parameters
         CALL VCA_DER( &
              HAMILTONIAN,KINEDEN, &
              P,WDES,NONLR_S,NONL_S,W,LATT_CUR, &
              T_INFO,INFO,IO,KPOINTS,GRID,GRID_SOFT, &
              GRIDC,GRIDUS,C_TO_US,SOFT_TO_C,SYMM, &
              CHTOT,DENCOR,CVTOT,CSTRF, &
              CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM, &
              CHDEN,SV,DOS,DOSI,CHAM, &
              LMDIM,IRDMAX,NEDOS)
          CALL STOP_TIMING("G",TIU6,'VCADER')
      ENDIF

      IF (.NOT.INFO%LPOTOK) THEN


      CALL POTLOK(GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES,  &
                  INFO,P,T_INFO,E,LATT_CUR, &
                  CHTOT,CSTRF,CVTOT,DENCOR,SV, SOFT_TO_C,XCSIF)

      CALL POTLOK_METAGGA(KINEDEN, &
                  GRID,GRIDC,GRID_SOFT,WDES%COMM_INTER,WDES,INFO,P,T_INFO,E,LATT_CUR, &
                  CHDEN,CHTOT,DENCOR,CVTOT,SV,HAMILTONIAN%MUTOT,HAMILTONIAN%MU,SOFT_TO_C,XCSIF)

      CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
                  LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

      CALL SETDIJ_AVEC(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
                  LMDIM,CDIJ,HAMILTONIAN%AVTOT, NONLR_S, NONL_S, IRDMAX)

      CALL SET_DD_MAGATOM(WDES, T_INFO, P, LMDIM, CDIJ)

      CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
         WDES%NCDIJ, LMDIM, CDIJ, RHOLM, CRHODE, &
         E, LMETA=.FALSE., LASPH =INFO%LASPH, LCOREL=.FALSE.)

      CALL UPDATE_CMBJ(GRIDC,T_INFO,LATT_CUR,IO%IU6)

      CALL STOP_TIMING("G",TIU6,'POTLOK')

      ENDIF

      CALL EGRAD_EFG_PW_HAR_ONLY(T_INFO,LATT_CUR,WDES,GRIDC,P,CSTRF,CHTOT(:,1))
      CALL EGRAD_WRITE_EFG(T_INFO,WDES,IO)

      IF (WDES%NCDIJ>1) THEN
         CALL HYPERFINE_PW(T_INFO,LATT_CUR,GRIDC,CHTOT(:,2))
         CALL HYPERFINE_WRITE_TENSORS(T_INFO,WDES,IO)
      ENDIF

      IF (WDES%ISPIN==2) THEN
         CALL DMATRIX_CALC(WDES,W,T_INFO,LATT_CUR,P,NONL_S,NONLR_S,GRID,GRIDC,CHTOT,LMDIM,CRHODE,IO%IU6,IO%IU0)
      ENDIF

!=======================================================================
!
!  electronic polarisation
!
!=======================================================================
      IF (LBERRY) THEN
      io_begin
         IF (IO%IU6>=0) THEN
            WRITE(TIU6,7230) IGPAR,NPPSTR
         ENDIF
      io_end
         CALL BERRY(WDES,W,GRID,GRIDC,GRIDUS,C_TO_US,LATT_CUR,KPOINTS, &
           IGPAR,NPPSTR,P,T_INFO,LMDIM,CQIJ,IRDMAX,INFO%LOVERL,IO%IU6, &
           DIP%POSCEN)
      ENDIF

 7230 FORMAT(/' Berry-Phase calculation of electronic polarization'// &
     &        '   IGPAR = ',I1,'   NPPSTR = ',I4/)
      IF (IO%IU6>=0) WRITE(TIU6,130)

      IF (LWANNIER()) THEN
         CALL LOCALIZE(W,WDES,GRID,GRIDC,GRIDUS,C_TO_US,LATT_CUR,T_INFO,P, &
        &     LMDIM,INFO%LOVERL,IRDMAX,IO%IU0,IO%IU6)
      ENDIF
      IF ( (DYN%IBRION /= -1 .AND. DYN%IBRION /=8 .AND. DYN%IBRION /=7) .AND. LEPSILON ) THEN
         ! if LEPSILON is set LR_SKELETON is called every NBLOCK*KBLOCK steps
         ! to calculate autocorrelation functions of the polarization
         IF (MOD(NSTEP-1,DYN%NBLOCK*DYN%KBLOCK)==0) THEN
            CALL LR_SKELETON( &
              KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
              T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
              GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
              CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
              CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
              CHDEN,SV,DOS,DOSI,CHAM, &
              DYN%IBRION,LMDIM,IRDMAX,NEDOS, &
              TOTEN,EFERMI, TIFOR, LFAST=.TRUE. )
         ENDIF

         CALL STORE_DIPOLE_AND_AUTOCORRELATE( T_INFO, LATT_CUR, DYN%POTIM, DYN%NSW, NSTEP, IO)
      ENDIF
!=======================================================================
! add chain forces and constrain forces
!=======================================================================
      CALL VCA_FORCE(T_INFO, TIFOR)

      IF (DYN%IBRION /=5 .AND. DYN%IBRION /=6  .AND. DYN%IBRION /=7 ) &
      CALL SET_SELECTED_FORCES_ZERO(T_INFO,DYN%VEL,TIFOR,LATT_CUR)

      ! scale forces and total energy by scaling argument
      ! this allows  integration from ideal gas to liquid state
      ! alternatively if the DYNMATFULL file exists add forces from
      ! second order Hessian matrix (allows thermodynamic integration from
      ! harmonic model)
      ! maybe this should go over to the force.F file
      ! of course it is unique since it  changes radically energy
      CALL DYNMATFULL_ENERGY_FORCE(SCALEE, T_INFO%NIONS, DYN%POSION, TOTEN, TIFOR, LATT_CUR%A , IO%IU0)
      IF (SCALEE/=1) THEN
         EENTROPY=E%EENTROPY*SCALEE
         IF (IO%IU6>=0) WRITE(TIU6,7261) '  SCALED FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)', &
            TOTEN,TOTEN-EENTROPY,TOTEN-EENTROPY/(2+NORDER)
      ENDIF

      CALL CHAIN_FORCE(T_INFO%NIONS,DYN%POSION,TOTEN,TIFOR, &
           LATT_CUR%A,LATT_CUR%B,IO%IU6)

      CALL PARALLEL_TEMPERING(NSTEP,T_INFO%NIONS,DYN%POSION,DYN%VEL,TOTEN,TIFOR,DYN%TEBEG,DYN%TEEND, &
           LATT_CUR%A,LATT_CUR%B,IO%IU6)

      IF (DYN%IBRION /=5 .AND. DYN%IBRION /=6 .AND. DYN%IBRION /=7 ) &
      CALL SET_SELECTED_FORCES_ZERO(T_INFO,DYN%VEL,TIFOR,LATT_CUR)

      EKIN=0
      EKIN_LAT=0
      TEIN=0
      ES =0
      EPS=0

!     CYCLE ion ! if uncommented VASP iterates forever without ionic upd

      DO I=1,WDES%NCDIJ
         RHOTOT(I)=RHO0(GRIDC, CHTOT(1,I))
      END DO
    END IF ! LJ_ONLY

      IF (LJ_IS_ACTIVE()) THEN
        !FIXME: TOTEN should be declared explicitly
        TOTEN = 0.0_q
        CALL FORCE_AND_STRESS_LJ(GRIDC,IO,WDES,P,LATT_CUR,DYN,T_INFO,TSIF,TIFOR,TOTEN)

        IF (DYN%IBRION /=5 .AND. DYN%IBRION /=6 .AND. DYN%IBRION /=7 ) &
          CALL SET_SELECTED_FORCES_ZERO(T_INFO,DYN%VEL,TIFOR,LATT_CUR)
      END IF

      DYN%AC=  LATT_CUR%A

    IF (.NOT. LJ_ONLY) THEN
!=======================================================================
!
! write statement for static calculations with possibly extensive
! post processing (GW, dielectric response, ...)
!
!=======================================================================
      IF (DYN%IBRION/=0) THEN
         CALL XML_TAG("energy")
         EENTROPY=E%EENTROPY*SCALEE
         CALL XML_ENERGY(TOTEN,TOTEN-EENTROPY/(2+NORDER),EENTROPY)
         CALL XML_CLOSE_TAG
      ENDIF

      IF (DYN%IBRION==-1 .OR. DYN%IBRION==5 .OR. DYN%IBRION==6 .OR. &
           DYN%IBRION==7.OR.DYN%IBRION==8) THEN
         EENTROPY=E%EENTROPY*SCALEE
         io_begin
         WRITE(17 ,7281,ADVANCE='NO') NSTEP,TOTEN, &
              TOTEN-EENTROPY/(2+NORDER),EENTROPY
         IF (IO%IU0>=0) &
              WRITE(TIU0,7281,ADVANCE='NO') NSTEP,TOTEN, &
              TOTEN-EENTROPY/(2+NORDER),EENTROPY

         IF ( WDES%NCDIJ>=2 ) THEN
           WRITE(17,77281) RHOTOT(2:WDES%NCDIJ)
           IF (IO%IU0>=0) WRITE(TIU0,77281) RHOTOT(2:WDES%NCDIJ)
         ELSE
           WRITE(17,*)
           IF (IO%IU0>=0) WRITE(TIU0,*)
         ENDIF

         IF (M_CONSTRAINED()) CALL WRITE_CONSTRAINED_M(17,.TRUE.)

         io_end
      ENDIF
    END IF ! LJ_ONLY
!=======================================================================
! IBRION = 0  molecular dynamics
! ------------------------------
! ) calculate the accelerations in reduced units
!   the scaling is brain-damaging, so here it is in a more native form
!    convert dE/dr from EV/Angst to J/m   EVTOJ/1E-10
!    divide by mass                       1/T_INFO%POMASS*AMTOKG
!    multiply by timestep*2               (DYN%POTIM*1E-15)**2
! ) transform to direct mesh KARDIR
! ) integrate the equations of motion for the ions using nose dynamic
!   and a predictor-corrector scheme
!=======================================================================
ibrion: IF (DYN%IBRION==0) THEN

        FACT=(DYN%POTIM**2)*EVTOJ/AMTOKG *1E-10_q
        NI=1
        DO NT=1,T_INFO%NTYP
        DO NI=NI,T_INFO%NITYP(NT)+NI-1
          DYN%D2C(1,NI)=TIFOR(1,NI)*FACT/2/T_INFO%POMASS(NT)
          DYN%D2C(2,NI)=TIFOR(2,NI)*FACT/2/T_INFO%POMASS(NT)
          DYN%D2C(3,NI)=TIFOR(3,NI)*FACT/2/T_INFO%POMASS(NT)
        ENDDO
        ENDDO

        CALL KARDIR(T_INFO%NIONS,DYN%D2C,LATT_CUR%B)

        ISCALE=0
        IF (DYN%SMASS==-1 .AND. MOD(NSTEP-1,DYN%NBLOCK)==0 ) THEN
          ISCALE=1
        ENDIF

        DYN%TEMP=DYN%TEBEG+(DYN%TEEND-DYN%TEBEG)*NSTEP/ABS(DYN%NSW)

#ifdef tbdyn
        IF (NSTEP==1) THEN
           IF (IO%IU6>0) THEN
             WRITE(IO%IU6,ADVANCE='NO',FMT='(3X,A22)') '        RANDOM_SEED = '
             DO i=1,K_SEED
               WRITE(IO%IU6,ADVANCE='NO',FMT='(X,I16)') SEED_INIT(i)
             ENDDO
             WRITE(IO%IU6,ADVANCE='NO',FMT='(/)')

             WRITE(g_io%REPORT,ADVANCE='NO',FMT='(3X,A22)') '        RANDOM_SEED = '
             DO i=1,K_SEED
               WRITE(g_io%REPORT,ADVANCE='NO',FMT='(X,I16)') SEED_INIT(i)
             ENDDO
             WRITE(g_io%REPORT,ADVANCE='NO',FMT='(/)')
           ENDIF
        ENDIF

        CALL STEP_tb(DYN,T_INFO,INFO,LATT_CUR,EKIN,EKIN_LAT,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM,&
                          IO,TOTEN,g_io,WDES,TSIF,ISCALE,TEIN,TIFOR)


        CALL RANE_ION(RDUM,GET=SEED(1:K_SEED))

        IF (IO%IU6>0) THEN
           WRITE(IO%IU6,ADVANCE='NO',FMT='(3X,A22)') '        RANDOM_SEED = '
           DO i=1,K_SEED
              WRITE(IO%IU6,ADVANCE='NO',FMT='(X,I16)') SEED(i)
           ENDDO
           WRITE(IO%IU6,ADVANCE='NO',FMT='(/)')

           WRITE(g_io%REPORT,ADVANCE='NO',FMT='(/,3X,A22)')  '        RANDOM_SEED = '
           DO i=1,K_SEED
             WRITE(g_io%REPORT,ADVANCE='NO',FMT='(X,I16)') SEED(i)
           ENDDO
           WRITE(g_io%REPORT,ADVANCE='NO',FMT='(/)')
        ENDIF
#else
        IF (.NOT. T_INFO%LSDYN) THEN
           CALL SYMVEL (T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS, &
           DYN%POSION,DYN%D2C,LATT_CUR%A,LATT_CUR%B)
        ENDIF

        CALL STEP(DYN%INIT,ISCALE,T_INFO%NIONS,LATT_CUR%A,LATT_CUR%ANORM,DYN%D2C,DYN%SMASS,DYN%POSION,DYN%POSIOC, &
             DYN%POTIM,T_INFO%POMASS,T_INFO%NTYP,T_INFO%ITYP,DYN%TEMP,DYN%VEL,DYN%D2,DYN%D3,DYN%SNOSE, &
             EKIN,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM, IO%IU6)

        TEIN = 2*EKIN/BOLKEV/NDEGREES_OF_FREEDOM
#endif
        CALL STOP_TIMING("G",TIU6,'IONSTEP')

! sum energy of images along chain

        CALLMPI_C( sum_chain( TOTEN ))
        CALLMPI_C( sum_chain( EKIN ))
        CALLMPI_C( sum_chain( EKIN_LAT ))
        CALLMPI_C( sum_chain( ES ))
        CALLMPI_C( sum_chain( EPS ))

        ETOTAL=TOTEN+EKIN+ES+EPS

!  report  energy  of electrons + kinetic energy + nose-energy
        io_begin

        WRITE(TIU6,7260) TOTEN,EKIN-EKIN_LAT,EKIN_LAT,TEIN,ES,EPS,ETOTAL

        CALL OFIELD_WRITE(TIU6)
        IF (LJ_IS_ACTIVE()) THEN
          WRITE(TIU6, 7262) TEIN * BOLKEV / LJ_EPSILON
    7262  FORMAT('  reduced temperature T*= ',F16.6, ' eps/DOF'/)
        END IF


        CALL XML_TAG("energy")
        CALL XML_ENERGY(TOTEN,TOTEN-E%EENTROPY/(2+NORDER),E%EENTROPY)
        CALL XML_TAG_REAL("kinetic",EKIN-EKIN_LAT)
        CALL XML_TAG_REAL("lattice kinetic",EKIN_LAT)
        CALL XML_TAG_REAL("nosepot",ES)
        CALL XML_TAG_REAL("nosekinetic",EPS)
        CALL XML_TAG_REAL("total",ETOTAL)
        CALL XML_CLOSE_TAG

        WRITE(17,7280,ADVANCE='NO') NSTEP,TEIN,ETOTAL,TOTEN, &
             TOTEN-E%EENTROPY/(2+NORDER),EKIN,ES,EPS
        IF (IO%IU0>=0) &
             WRITE(TIU0,7280,ADVANCE='NO')NSTEP,TEIN,ETOTAL,TOTEN, &
             TOTEN-E%EENTROPY/(2+NORDER),EKIN,ES,EPS
        IF (WDES%NCDIJ>=2) THEN
           WRITE(17,77280) RHOTOT(2:WDES%NCDIJ)
           IF (IO%IU0>=0) WRITE(TIU0,77280) RHOTOT(2:WDES%NCDIJ)
        ELSE
           WRITE(17,*)
           IF (IO%IU0>=0) WRITE(TIU0,*)
        ENDIF

    IF (.NOT. LJ_ONLY) THEN
        WRITE(TIU6,7270) DISMAX
    END IF
        io_end

7260  FORMAT(/ &
     &        '  ENERGY OF THE ELECTRON-ION-THERMOSTAT SYSTEM (eV)'/ &
     &        '  ---------------------------------------------------'/ &
     &        '% ion-electron   TOTEN  = ',F16.6,'  see above'/ &
     &        '  kinetic energy EKIN   = ',F16.6/ &
     &        '  kin. lattice  EKIN_LAT= ',F16.6, &
     &        '  (temperature',F8.2,' K)'/ &
     &        '  nose potential ES     = ',F16.6/ &
     &        '  nose kinetic   EPS    = ',F16.6/ &
     &        '  ---------------------------------------------------'/ &
     &        '  total energy   ETOTAL = ',F16.6,' eV'/)

 7270 FORMAT( '  maximum distance moved by ions :',E14.2/)

 7280 FORMAT(I5,' T= ',F6.0,' E= ',E14.8, &
     &   ' F= ',E14.8,' E0= ',E14.8,1X,' EK= ',E11.5, &
     &   ' SP= ',E8.2,' SK= ',E8.2)
77280 FORMAT(' mag=',3F11.3)

!=======================================================================
! DYN%IBRION =
! 1  quasi-Newton algorithm
! 2  conjugate gradient
! 3  quickmin
! 4  not supported yet
! 5  finite differences
! 6  finite differences with symmetry
! 7  linear response
!  IBRION ==5 finite differences
!=======================================================================
      ELSE IF (DYN%IBRION==5) THEN ibrion
       IF (.NOT. LJ_ONLY) THEN
        DYN%POSIOC=DYN%POSION
        CALL FINITE_DIFF( INFO%LSTOP, DYN%POTIM, T_INFO%NIONS, T_INFO%NTYP, &
             T_INFO%NITYP, T_INFO%POMASS, DYN%POSION, TIFOR, DYN%NFREE, &
             T_INFO%LSDYN,T_INFO%LSFOR, LATT_CUR%A, LATT_CUR%B,  &
             IO%IU6, IO%IU0 , IO%NWRITE)
        CALL LATTIC(LATT_CUR)
        ! we need to reinitialise the symmetry code at this point
        ! the number of k-points is changed on the fly
        IF (SYMM%ISYM>0) THEN
          CALL INISYM(LATT_CUR%A,DYN%POSION,DYN%VEL,T_INFO%LSFOR, &
             T_INFO%LSDYN,T_INFO%NTYP,T_INFO%NITYP,NIOND, &
             SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
             SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,NCDIJ,IO%IU6)
#ifdef oldsym
          CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
               SYMM%ISYM>=0.AND..NOT.WDES%LSORBIT.AND..NOT.WDES%LSPIRAL, &
               T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IO%IU6,IO%IU0)
#else
          CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
               SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
               T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IO%IU6,IO%IU0)
#endif
          CALL KPAR_SYNC_ALL(WDES,W)
          CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
          CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
          IF (INFO%LONESW) THEN
             CALL DEALLOCW(W_F)
             CALL DEALLOCW(W_G)
             DEALLOCATE(CHAM, CHF)
             CALL ALLOCW(WDES,W_F,WTMP,WTMP)
             CALL ALLOCW(WDES,W_G,WTMP,WTMP)
             ALLOCATE(CHAM(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN), &
               CHF (WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN))
          ENDIF
        ENDIF
       END IF !LJ_ONLY
      ELSE IF (DYN%IBRION==6) THEN ibrion
       IF (.NOT. LJ_ONLY) THEN
        DYN%POSIOC=DYN%POSION
        CALL FINITE_DIFF_ID( INFO%LSTOP, DYN%POTIM, T_INFO,  &
             DYN%POSION, TOTEN, TIFOR, TSIF, DYN%NFREE, &
             DYN%ISIF>=3, LATT_CUR%A, LATT_CUR%B, WDES%SAXIS, SYMM, NCDIJ, &
             ISPECIAL, DYN%TEBEG, IO%IU6, IO%IU0 , IO%NWRITE)
        CALL LATTIC(LATT_CUR)
        ! we need to reinitialise the symmetry code at this point
        ! the number of k-points is changed on the fly
        IF (SYMM%ISYM>0) THEN
          CALL INISYM(LATT_CUR%A,DYN%POSION,DYN%VEL,T_INFO%LSFOR, &
             T_INFO%LSDYN,T_INFO%NTYP,T_INFO%NITYP,NIOND, &
             SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
             SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,NCDIJ,IO%IU6)
#ifdef oldsym
          CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
               SYMM%ISYM>=0.AND..NOT.WDES%LSORBIT.AND..NOT.WDES%LSPIRAL, &
               T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IO%IU6,IO%IU0)
#else
          CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
               SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
               T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IO%IU6,IO%IU0)
#endif
          CALL KPAR_SYNC_ALL(WDES,W)
          CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
          CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
          IF (INFO%LONESW) THEN
             CALL DEALLOCW(W_F)
             CALL DEALLOCW(W_G)
             DEALLOCATE(CHAM, CHF)
             CALL ALLOCW(WDES,W_F,WTMP,WTMP)
             CALL ALLOCW(WDES,W_G,WTMP,WTMP)
             ALLOCATE(CHAM(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN), &
             CHF (WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN))
          ENDIF
        ENDIF
       END IF !LJ_ONLY

      ELSE IF (DYN%IBRION==7.OR.DYN%IBRION==8) THEN ibrion

       IF (.NOT. LJ_ONLY) THEN
         CALL LR_SKELETON( &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
          T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI,CHAM, &
          DYN%IBRION,LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI, TIFOR)

       END IF ! LJ_ONLY
!=======================================================================
! meaning of DYN%ISIF :
!  DYN%ISIF  calculate                           relax
!        force     stress                    ions      lattice
!   0     X                                   X
!   1     X        uniform                    X
!   2     X          X                        X
!   3     X          X                        X          X
!   4     X          X                        X          X **
!   5     X          X                                   X **
!   6     X          X                                   X
!   7     X          X                                  uniform
!
!   **  for DYN%ISIF=4 & DYN%ISIF=5 isotropic pressure will be subtracted
!       (-> cell volume constant, optimize only the cell shape)
!
!=======================================================================
      ELSE IF (DYN%IBRION>0) THEN ibrion
       IF (.NOT. LJ_ONLY) THEN
! sum energy of images along chain
        EENTROPY=E%EENTROPY*SCALEE

        IF (IMAGES==0) THEN
           CALLMPI_C( sum_chain( TOTEN ))
           CALLMPI_C( sum_chain( EENTROPY ))
        ENDIF

        io_begin
        WRITE(17,7281,ADVANCE='NO') NSTEP,TOTEN, &
             TOTEN-EENTROPY/(2+NORDER),TOTEN-TOTENG
        IF (IO%IU0>=0) &
             WRITE(TIU0,7281,ADVANCE='NO')NSTEP,TOTEN, &
             TOTEN-EENTROPY/(2+NORDER),TOTEN-TOTENG

        IF ( WDES%NCDIJ>=2 ) THEN
           WRITE(17,77281) RHOTOT(2:WDES%NCDIJ)
           IF (IO%IU0>=0) WRITE(TIU0,77281) RHOTOT(2:WDES%NCDIJ)
        ELSE
           WRITE(17,*)
           IF (IO%IU0>=0) WRITE(TIU0,*)
        ENDIF

        io_end
 7281 FORMAT(I4,' F= ',E14.8,' E0= ',E14.8,1X,' d E =',E12.6)
77281 FORMAT('  mag=',3F11.4)
!-----------------------------------------------------------------------
!  set DYN%D2C to forces in cartesian coordinates multiplied by FACT
!  FACT is determined from timestep in a way, that a stable timestep
!   gives a good trial step
!-----------------------------------------------------------------------
        FACT=0
        IF (DYN%ISIF<5) FACT=10*DYN%POTIM*EVTOJ/AMTOKG *1E-10_q
        LSTOP2=.TRUE.

        NI=1
        DO NT=1,T_INFO%NTYP
        DO NI=NI,T_INFO%NITYP(NT)+NI-1
           DYN%D2C(1,NI)=TIFOR(1,NI)*FACT
           DYN%D2C(2,NI)=TIFOR(2,NI)*FACT
           DYN%D2C(3,NI)=TIFOR(3,NI)*FACT
           IF (SQRT(TIFOR(1,NI)**2+TIFOR(2,NI)**2+TIFOR(3,NI)**2) &
                &       >ABS(DYN%EDIFFG)) LSTOP2=.FALSE.
        ENDDO
        ENDDO
! for all DYN%ISIF greater or equal 3 cell shape optimisations will be done
        FACTSI = 0
        IF (DYN%ISIF>=3) FACTSI=10*DYN%POTIM*EVTOJ/AMTOKG/T_INFO%NIONS *1E-10_q

        DO I=1,3
        DO K=1,3
           D2SIF(I,K)=TSIF(I,K)*FACTSI
        ENDDO
        D2SIF(I,I)=D2SIF(I,I)-DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA*FACTSI
        ENDDO
! For DYN%ISIF =4 or =5 we take only pure shear stresses: subtract pressure
        IF ((DYN%ISIF==4).OR.(DYN%ISIF==5)) THEN
           DO I=1,3
              D2SIF(I,I)=D2SIF(I,I)-PRESS*FACTSI
           ENDDO
        ENDIF
! For DYN%ISIF =7 take only pressure (volume relaxation)
        IF (DYN%ISIF==7) THEN
           DO I=1,3
           DO J=1,3
              D2SIF(J,I)=0
           ENDDO
           D2SIF(I,I)=PRESS*FACTSI
           ENDDO
        ENDIF

        CALL CONSTR_CELL_RELAX(D2SIF)

        IF (FACTSI/=0) THEN
           DO I=1,3
           DO J=1,3
              IF (FACTSI/=0) THEN
                 IF (ABS(D2SIF(J,I))/FACTSI/T_INFO%NIONS>ABS(DYN%EDIFFG)) LSTOP2=.FALSE.
              ENDIF
           ENDDO
          ENDDO
        ENDIF
!-----------------------------------------------------------------------
!  do relaxations using diverse algorithms
!-----------------------------------------------------------------------
        ! change of the energy between two ionic step used as stopping criterion
        INFO%LSTOP=(ABS(TOTEN-TOTENG)<DYN%EDIFFG)

        CALLMPI_C ( and_chain( LSTOP2 ))
        CALLMPI_C ( and_chain( INFO%LSTOP ))

        ! IFLAG=0 means no reinit of wavefunction prediction
        IFLAG=0
        IF (DYN%IBRION==1) THEN
           CALL BRIONS(T_INFO%NIONS,DYN%POSION,DYN%POSIOC,DYN%D2C,LATT_CUR%A,LATT_CUR%B,D2SIF, &
                MAX(DYN%NSW+1,DYN%NFREE+1),DYN%NFREE,IO%IU6,IO%IU0,FACT,FACTSI,E1TEST)
! Sometimes there is the danger that the optimisation scheme (especially
! the Broyden scheme) fools itself by performing a 'too small' step - to
! avoid this use a second break condition ('trial step energy change'):
           INFO%LSTOP=INFO%LSTOP .AND. (ABS(E1TEST) < DYN%EDIFFG)
! if we have very small forces (small trial energy change) we can stop
           INFO%LSTOP=INFO%LSTOP .OR. (ABS(E1TEST) < 0.1_q*DYN%EDIFFG)
           TOTENG=TOTEN
!-----------------------------------------------------------------------
        ELSE IF (DYN%IBRION==3 .AND. DYN%SMASS<=0) THEN
           CALL ION_VEL_QUENCH(T_INFO%NIONS,LATT_CUR%A,LATT_CUR%B,IO%IU6,IO%IU0, &
                T_INFO%LSDYN, &
                DYN%POSION,DYN%POSIOC,FACT,DYN%D2C,FACTSI,D2SIF,DYN%D2,E1TEST)
           IF (IFLAG==1) INFO%LSTOP=INFO%LSTOP .OR. (ABS(E1TEST) < 0.1_q*DYN%EDIFFG)
           TOTENG=TOTEN
!-----------------------------------------------------------------------
        ELSE IF (DYN%IBRION==3) THEN
           CALL IONDAMPED(T_INFO%NIONS,LATT_CUR%A,LATT_CUR%B,IO%IU6,IO%IU0, &
                T_INFO%LSDYN, &
                DYN%POSION,DYN%POSIOC,FACT,DYN%D2C,FACTSI,D2SIF,DYN%D2,E1TEST,DYN%SMASS)
           IF (IFLAG==1) INFO%LSTOP=INFO%LSTOP .OR. (ABS(E1TEST) < 0.1_q*DYN%EDIFFG)
           TOTENG=TOTEN
!-----------------------------------------------------------------------
        ELSE IF (DYN%IBRION==2) THEN
           IFLAG=1
           IF (NSTEP==1) IFLAG=0
!  set accuracy of energy (determines whether cubic interpolation is used)
           EACC=MAX(ABS(INFO%EDIFF),ABS(ECONV))

           CALLMPI_C( sum_chain( EACC))
           IF ( LHYPER_NUDGE() ) EACC=1E10    ! energy not very accurate, use only force

           CALL IONCGR(IFLAG,T_INFO%NIONS,TOTEN,LATT_CUR%A,LATT_CUR%B,DYN%NFREE,DYN%POSION,DYN%POSIOC, &
                FACT,DYN%D2C,FACTSI,D2SIF,DYN%D2,DYN%D3,DISMAX,IO%IU6,IO%IU0, &
                EACC,DYN%EDIFFG,E1TEST,LSTOP2)
!    if IFLAG=1 new trial step -> reinit of waveprediction
           INFO%LSTOP=.FALSE.
           IF (IFLAG==1) THEN
              INFO%LSTOP=(ABS(TOTEN-TOTENG)<DYN%EDIFFG)
              TOTENG=TOTEN
           ENDIF
           IF (IFLAG==2) INFO%LSTOP=.TRUE.
! if we have very small forces (small trial energy change) we can stop
           IF (IFLAG==1) INFO%LSTOP=INFO%LSTOP .OR. (ABS(E1TEST) < 0.1_q*DYN%EDIFFG)
!tb start
        ! dimer method
        ELSE IF (DYN%IBRION==44) THEN
           CALL dimer(DYN,T_INFO,INFO,LATT_CUR,IO,TOTEN,FACT,LSTOP2)
           CALLMPI_C ( and_chain( LSTOP2 ))
           CALLMPI_C ( and_chain( INFO%LSTOP ))
           TOTENG=TOTEN

        ! damped velocity Verler to calculate IRC (not optimization!!!)
        ELSE IF (DYN%IBRION==40) THEN
           FACT=(DYN%POTIM**2)*EVTOJ/AMTOKG *1E-10_q
           NI=1
           DO NT=1,T_INFO%NTYP
             DO NI=NI,T_INFO%NITYP(NT)+NI-1
               DYN%D2C(1,NI)=TIFOR(1,NI)*FACT/2/T_INFO%POMASS(NT)
               DYN%D2C(2,NI)=TIFOR(2,NI)*FACT/2/T_INFO%POMASS(NT)
               DYN%D2C(3,NI)=TIFOR(3,NI)*FACT/2/T_INFO%POMASS(NT)
             ENDDO
           ENDDO
           CALL KARDIR(T_INFO%NIONS,DYN%D2C,LATT_CUR%B)

           CALL dvv(DYN,T_INFO,INFO,LATT_CUR,IO,TOTEN,FACT)
           FACT=10*DYN%POTIM*EVTOJ/AMTOKG *1E-10_q
           NI=1
           DO NT=1,T_INFO%NTYP
             DO NI=NI,T_INFO%NITYP(NT)+NI-1
               DYN%D2C(1,NI)=TIFOR(1,NI)*FACT
               DYN%D2C(2,NI)=TIFOR(2,NI)*FACT
               DYN%D2C(3,NI)=TIFOR(3,NI)*FACT
             ENDDO
           ENDDO
           TOTENG=TOTEN
!tb end

!-----------------------------------------------------------------------
! interactive mode
!-----------------------------------------------------------------------
        ELSE IF (DYN%IBRION==11) THEN
           CALL INPOS(LATT_CUR, T_INFO, DYN, IO%IU6, IO%IU0, INFO%LSTOP, WDES%COMM)
        ENDIF

! restrict volume for constant volume relaxation
        IF (DYN%ISIF==4 .OR. DYN%ISIF==5) THEN
           OMEGA_OLD=LATT_CUR%OMEGA
           CALL LATTIC(LATT_CUR)
           SCALEQ=(ABS(OMEGA_OLD) / ABS(LATT_CUR%OMEGA))**(1._q/3._q)
           DO I=1,3
              LATT_CUR%A(1,I)=SCALEQ*LATT_CUR%A(1,I)
              LATT_CUR%A(2,I)=SCALEQ*LATT_CUR%A(2,I)
              LATT_CUR%A(3,I)=SCALEQ*LATT_CUR%A(3,I)
           ENDDO
        ENDIF
        CALL LATTIC(LATT_CUR)
!-----------------------------------------------------------------------
!  reinitialize the prediction algorithm for the wavefunction if needed
!-----------------------------------------------------------------------
        PRED%INIPRE=3
        IF ( PRED%IWAVPR >=12 .AND. &
             &     (ABS(TOTEN-TOTENG)/T_INFO%NIONS>1.0_q .OR. IFLAG==1)) THEN
           CALL WAVPRE_NOIO(GRIDC,P,PRED,T_INFO,W,WDES,LATT_CUR,IO%LOPEN, &
                CHTOT,RHOLM,N_MIX_PAW, CSTRF, KINEDEN, LMDIM,CQIJ,INFO%LOVERL,IO%IU0)

        ELSE IF ( PRED%IWAVPR >=2 .AND. PRED%IWAVPR <10   .AND. &
             &     (ABS(TOTEN-TOTENG)/T_INFO%NIONS>1.0_q .OR. IFLAG==1)) THEN
           CALL WAVPRE(GRIDC,P,PRED,T_INFO,W,WDES,LATT_CUR,IO%LOPEN, &
                CHTOT,RHOLM,N_MIX_PAW, CSTRF, LMDIM,CQIJ,INFO%LOVERL,IO%IU0)
        ENDIF

        ! use forces as stopping criterion if EDIFFG<0
        IF (DYN%EDIFFG<0) INFO%LSTOP=LSTOP2
        io_begin
        WRITE(TIU6,130)

        IF (INFO%LSTOP) THEN
           IF (IO%IU0>=0) &
                WRITE(TIU0,*) 'reached required accuracy - stopping ', &
                'structural energy minimisation'
           WRITE(TIU6,*) ' '
           WRITE(TIU6,*) 'reached required accuracy - stopping ', &
                'structural energy minimisation'
        ENDIF
        io_end
       END IF !LJ_ONLY
      ELSE IF ( DYN%IBRION == -1 .AND. (LEPSILON  .OR. KINTER/=0 .OR. LMAGBLOCH)) THEN ibrion
       IF (.NOT. LJ_ONLY) THEN
         CALL LR_SKELETON( &
              KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
              T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
              GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
              CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
              CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
              CHDEN,SV,DOS,DOSI,CHAM, &
              DYN%IBRION,LMDIM,IRDMAX,NEDOS, &
              TOTEN,EFERMI, TIFOR)
       END IF !LJ_ONLY
      ELSE IF ( DYN%IBRION == -1 .AND. LAUGER) THEN ibrion
       IF (.NOT. LJ_ONLY) THEN
         CALL CALCULATE_AUGER(W,WDES,LATT_CUR,SYMM,T_INFO,GRID,P,NONL_S,KPOINTS,IO)
       END IF !LJ_ONLY
      ELSE IF ( LCHIMAG ) THEN ibrion
       IF (.NOT. LJ_ONLY) THEN
         CALL MLR_SKELETON( &
              HAMILTONIAN,W,WDES,GRID,GRID_SOFT,GRIDC,GRIDUS,C_TO_US,SOFT_TO_C,KPOINTS,LATT_CUR,LATT_INI, &
              T_INFO,DYN,SYMM,P,NONL_S,NONLR_S,LMDIM,CDIJ,CQIJ,SV,E,INFO,IRDMAX,IO)
       END IF !LJ_ONLY
      ENDIF ibrion
!=======================================================================
!  update of ionic positions performed
!  in any case POSION should now hold the new positions and
!  POSIOC the old one
!=======================================================================

!   reached required number of time steps
      IF (NSTEP>=DYN%NSW) INFO%LSTOP=.TRUE.
#ifdef tbdyn
      IF (INFO%LSTOP) THEN
        CLOSE(g_io%REPORT)
      END IF
#endif
!   soft stop or hard stop
      IF (INFO%LSOFT) INFO%LSTOP=.TRUE.

!   if we need to pull the brake, then POSION is reset to POSIOC
!   except for molecular dynamics, where the next electronic
!   step should indeed correspond to POSIOC
      IF ( INFO%LSTOP .AND. DYN%IBRION>0 ) THEN
         LATT_CUR%A=DYN%AC
         DYN%POSION=DYN%POSIOC
      ENDIF

      CALLMPI ( M_bcast_d(WDES%COMM, DYN%POSION , T_INFO%NIONS*3))
      CALLMPI ( M_bcast_d(WDES%COMM, LATT_CUR%A(1,1) , 9))
      CALL LATTIC(LATT_CUR)

!=======================================================================
!  update mean temperature mean energy
!=======================================================================
      SMEAN =SMEAN +1._q/DYN%SNOSE(1)
      SMEAN0=SMEAN0+DYN%SNOSE(1)
      TMEAN =TMEAN +TEIN/DYN%SNOSE(1)
      TMEAN0=TMEAN0+TEIN

    IF (.NOT. LJ_ONLY) THEN
!=======================================================================
!  SMEAR_LOOP%ISMCNT != 0 Loop over several KPOINTS%SIGMA-values
!  set new smearing parameters and continue main loop
!=======================================================================
      IF (SMEAR_LOOP%ISMCNT/=0) THEN
         KPOINTS%ISMEAR=NINT(SMEAR_LOOP%SMEARS(2*SMEAR_LOOP%ISMCNT-1))
         KPOINTS%SIGMA=SMEAR_LOOP%SMEARS(2*SMEAR_LOOP%ISMCNT)

         SMEAR_LOOP%ISMCNT=SMEAR_LOOP%ISMCNT-1

         io_begin
         IF (IO%IU0>=0) &
              WRITE(TIU0,7283) KPOINTS%ISMEAR,KPOINTS%SIGMA
         WRITE(17,7283) KPOINTS%ISMEAR,KPOINTS%SIGMA
         io_end

7283     FORMAT('ISMEAR = ',I4,' SIGMA = ',F10.3)

         KPOINTS%LTET=((KPOINTS%ISMEAR<=-4).OR.(KPOINTS%ISMEAR>=30))
         IF (KPOINTS%ISMEAR==-6) KPOINTS%ISMEAR=-1
         IF (KPOINTS%ISMEAR>=0)  KPOINTS%ISMEAR=MOD(KPOINTS%ISMEAR,30)
         SIGMA=ABS(KPOINTS%SIGMA)

         CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
              INFO%NUP_DOWN, E%EENTROPY, EFERMI, SIGMA, .FALSE.,  &
              NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
      ENDIF
!=======================================================================
! tasks which are done all DYN%NBLOCK   steps:
!=======================================================================
!-----------------------------------------------------------------------
!  write out the position of the IONS to XDATCAR
!-----------------------------------------------------------------------
      IF (INFO%INICHG==3) THEN
        IF (MOD(NSTEP,DYN%NBLOCK)==1) THEN
           IF (IO%IU0>=0) &
           WRITE(TIU0,*) 'non selfconsistent'
           INFO%LCHCON=.FALSE.
        ELSE
           IF (IO%IU0>=0) &
           WRITE(TIU0,*) 'selfconsistent'
           INFO%LCHCON=.TRUE.
        ENDIF
      ENDIF
    END IF !LJ_ONLY

   nblock: IF (MOD(NSTEP,DYN%NBLOCK)==0) THEN
        io_begin
        ! tb: write lattice vectors to XDATCAR
        IF (DYN%ISIF>=3) THEN
           CALL XDAT_HEAD_AC(61, T_INFO, LATT_CUR, DYN, INFO%SZNAM1)
        ENDIF
        WRITE(61,'(A,I6)') 'Direct configuration=',NSTEP
        WRITE(61,7007) ((DYN%POSIOC(I,J),I=1,3),J=1,T_INFO%NIONS)
        IF (IO%LOPEN) CALL WFORCE(61)
        io_end

    IF (.NOT. LJ_ONLY) THEN
!-----------------------------------------------------------------------
! acummulate dos
!-----------------------------------------------------------------------
        DO ISP=1,WDES%NCDIJ
        DO I=1,NEDOS
          DDOSI(I,ISP)=DDOSI(I,ISP)+DOSI(I,ISP)
          DDOS (I,ISP)=DDOS (I,ISP)+DOS (I,ISP)
        ENDDO
        ENDDO
    END IF !LJ_ONLY

!-----------------------------------------------------------------------
! evaluate the pair-correlation function  using the exact places
! also sum up mean temperatur
!-----------------------------------------------------------------------
        SMEANP =SMEANP +1._q
        CALL SPACO(T_INFO%NIONS,1._q,DYN%POSIOC,DYN%AC,LATT_CUR%BNORM,PACO%SIPACO(0),PACO%NPACO,PACO%APACO)
      ENDIF nblock

!=======================================================================
! tasks which are done all DYN%NBLOCK*DYN%KBLOCK   steps
!=======================================================================
!-----------------------------------------------------------------------
! write  pair-correlation and density of states
! quantities are initialized to 0 at the beginning of the main-loop
!-----------------------------------------------------------------------
    wrtpair: IF (MOD(NSTEP,DYN%NBLOCK*DYN%KBLOCK)==0 .AND. DYN%IBRION ==0) THEN
       io_begin

       PCFAK = 1.5_q/PI/T_INFO%NIONS**2*LATT_CUR%OMEGA*(PACO%NPACO/PACO%APACO)**3
       WRITE(60,'(3E15.7)') TMEAN0/(DYN%NBLOCK*DYN%KBLOCK),TMEAN/SMEAN
       WRITE(TIU6,8022) SMEAN0/(DYN%NBLOCK*DYN%KBLOCK),TMEAN0/(DYN%NBLOCK*DYN%KBLOCK), &
            &                TMEAN/SMEAN
8022   FORMAT(/' mean value of Nose-termostat <S>:',F10.3, &
            &        ' mean value of <T> :',F10.3/ &
            &        ' mean temperature <T/S>/<1/S>  :',F10.3/)

       DO  I=0,PACO%NPACO-1
          WRITE(60,'(F7.3)') &
               PCFAK*PACO%SIPACO(I)/ REAL( 3*I*(I+1)+1 ,KIND=q) /SMEANP
       ENDDO
       IF (IO%LOPEN) CALL WFORCE(60)

    IF (.NOT. LJ_ONLY) THEN
       DELTAE=(KPOINTS%EMAX-KPOINTS%EMIN)/(NEDOS-1)

       WRITE(16,'(2F16.8,I5,2F16.8)') KPOINTS%EMAX,KPOINTS%EMIN,NEDOS,EFERMI,1.0
       DO I=1,NEDOS
          EN=KPOINTS%EMIN+DELTAE*(I-1)
          WRITE(16,7062) EN,(DDOS(I,ISP)/DYN%KBLOCK,ISP=1,WDES%ISPIN),(DDOSI(I,ISP)/DYN%KBLOCK,ISP=1,WDES%ISPIN)
       ENDDO
       IF (IO%LOPEN) CALL WFORCE(16)
7062   FORMAT(3X,F8.3,8E12.4)
       CALL XML_DOS(EFERMI, KPOINTS%EMIN, KPOINTS%EMAX, .FALSE., &
            DDOS, DDOSI, DOSPAR, NEDOS, LPAR, T_INFO%NIONP, WDES%NCDIJ)
    END IF !LJ_ONLY
       io_end
    ENDIF wrtpair
!-----------------------------------------------------------------------
!  update file CONTCAR
!-----------------------------------------------------------------------

!-----write out positions (only done on IONODE)
#ifdef MPI
    IF (WDES%COMM_KINTER%NODE_ME.EQ.1) THEN
#endif
     io_begin

      CALL OUTPOS(13,.TRUE.,T_INFO%SZNAM2,T_INFO,LATT_CUR%SCALE,DYN%AC,T_INFO%LSDYN,DYN%POSIOC)
      CALL OUTPOS_TRAIL(13,IO%LOPEN, LATT_CUR, T_INFO, DYN)

     io_end

    IF (.NOT. LJ_ONLY) THEN
!=======================================================================
!  append new chargedensity to file CHG
!=======================================================================
      IF (IO%LCHARG .AND. MOD(NSTEP,10)==1) THEN

      io_begin
      CALL OUTPOS(70,.FALSE.,INFO%SZNAM1,T_INFO,LATT_CUR%SCALE,LATT_CUR%A,.FALSE.,DYN%POSIOC)
      io_end

      DO ISP=1,WDES%NCDIJ
         CALL OUTCHG(GRIDC,70,.FALSE.,CHTOT(:,ISP))
      ENDDO
      io_begin
      IF (IO%LOPEN) CALL WFORCE(70)
      io_end
      ENDIF
    END IF !LJ_ONLY

#ifdef MPI
    ENDIF
#endif

    IF (.NOT. LJ_ONLY) THEN
!=======================================================================
! if ions were moved recalculate some quantities
!=======================================================================
!=======================================================================
! WAVPRE prediction of the new wavefunctions and charge-density
! if charge-density constant during ELM recalculate the charge-density
! according to overlapping atoms
! for relaxation jobs do not predict in the last step
!=======================================================================
      INFO%LPOTOK=.FALSE.
  prepare_next_step: &
    & IF ( .NOT. INFO%LSTOP .OR. DYN%IBRION==0 ) THEN

      ! extrapolate charge using  atomic charges
      IF (PRED%IWAVPR==1 .OR.PRED%IWAVPR==11) PRED%INIPRE=5
      ! extrapolate wavefunctions and charge
      IF (PRED%IWAVPR==2 .OR.PRED%IWAVPR==12) PRED%INIPRE=0
      ! mixed mode
      IF (PRED%IWAVPR==3 .OR.PRED%IWAVPR==13) PRED%INIPRE=4
      PRED%IPRE=0

      IF (PRED%IWAVPR >=11) THEN
        CALL WAVPRE_NOIO(GRIDC,P,PRED,T_INFO,W,WDES,LATT_CUR,IO%LOPEN, &
           CHTOT,RHOLM,N_MIX_PAW, CSTRF, KINEDEN, LMDIM,CQIJ,INFO%LOVERL,IO%IU0)
      ELSE IF (PRED%IWAVPR >=1 .AND. PRED%IWAVPR<10 ) THEN
        CALL WAVPRE(GRIDC,P,PRED,T_INFO,W,WDES,LATT_CUR,IO%LOPEN, &
           CHTOT,RHOLM,N_MIX_PAW, CSTRF, LMDIM,CQIJ,INFO%LOVERL,IO%IU0)
      ENDIF

      IF (PRED%IPRE<0) THEN
        io_begin
        IF (IO%IU0>=0) &
        WRITE(TIU0,*)'bond charge predicted'
        io_end
        PRED%IPRE=ABS(PRED%IPRE)
      ELSE
      ! PRED%IPRE < 0 then WAVPRE calculated new structure factor
      !     in all other cases we have to recalculate s.f.
         IF (INFO%TURBO==0) CALL STUFAK(GRIDC,T_INFO,CSTRF)
      ENDIF

      IF (INFO%LCHCON.AND.INFO%INICHG==2) THEN
        IF (IO%IU0>=0)  WRITE(TIU0,*)'charge from overlapping atoms'
        ! then initialize CRHODE and than RHOLM (PAW related occupancies)
        CALL DEPATO(WDES, LMDIM, CRHODE, INFO%LOVERL, P, T_INFO)
        CALL SET_RHO_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
           CRHODE, RHOLM)

        CALL RHOATO_WORK(.FALSE.,.FALSE.,GRIDC,T_INFO,LATT_CUR%B,P,CSTRF,CHTOT)
        ! set magnetization to 0
        DO ISP=2,WDES%NCDIJ
           CALL RC_ADD(CHTOT(1,ISP),0.0_q,CHTOT(1,ISP),0.0_q,CHTOT(1,ISP),GRIDC)
        ENDDO
        ! add Gaussian "charge-transfer" charges, if required
        CALL RHOADD_GAUSSIANS(T_INFO,LATT_CUR,P,GRIDC,NCDIJ,CHTOT,CSTRF)
        CALL RHOADD_GAUSSIANS_LIST(LATT_CUR,GRIDC,NCDIJ,CHTOT)
      ENDIF

      IF (INFO%LCORE) CALL RHOPAR(GRIDC,T_INFO,INFO,LATT_CUR,P,CSTRF,DENCOR,IO%IU6)
      IF (LDO_METAGGA()) CALL TAUPAR(GRIDC,T_INFO,LATT_CUR%B,LATT_CUR%OMEGA,P,CSTRF,KINEDEN%TAUC)

      CALL STOP_TIMING("G",TIU6,'WAVPRE')
!-----------------------------------------------------------------------
! call the ewald program to get the energy of the new ionic
! configuration
!-----------------------------------------------------------------------
      IF(INFO%TURBO==0) THEN
         CALL FEWALD(DYN%POSION,EWIFOR,LATT_CUR%A,LATT_CUR%B,LATT_CUR%ANORM,LATT_CUR%BNORM, &
              &     LATT_CUR%OMEGA,EWSIF,E%TEWEN,T_INFO%NTYP,P%ZVALF,T_INFO%VCA, &
              &     T_INFO%NIONS,NIOND,T_INFO%ITYP,T_INFO%NITYP,IO%IU6,.TRUE.)
      ELSE
         ALLOCATE(CWORK1(GRIDC%MPLWV))
         CALL POTION_PARTICLE_MESH(GRIDC,P,LATT_CUR,T_INFO,CWORK1,E%PSCENC,E%TEWEN,EWIFOR)
         DEALLOCATE(CWORK1)
      ENDIF

      CALL STOP_TIMING("G",TIU6,'FEWALD')

! volume might have changed restet IRDMAX
      IRDMAX=4*PI*PSDMX**3/3/(LATT_CUR%OMEGA/ &
     &     (GRIDC%NGPTAR(1)*GRIDC%NGPTAR(2)*GRIDC%NGPTAR(3)))+200
#ifdef usgrid
       IRDMAX=4*PI*PSDMX**3/3/(LATT_CUR%OMEGA/ &
     &        (GRIDUS%NGPTAR(1)*GRIDUS%NGPTAR(2)*GRIDUS%NGPTAR(3)))+200
#endif
!-----------------------------------------------------------------------
! if basis cell changed recalculate kinetic-energy array and tables
!-----------------------------------------------------------------------
      IF (DYN%ISIF>=3) THEN
        CALL GEN_INDEX(GRID,WDES, LATT_CUR%B,LATT_INI%B,-1,-1,.TRUE.)
        CALL STOP_TIMING("G",TIU6,'GENKIN')
        IF (WDES%LSPIRAL) CALL CLEANWAV(WDES,W,INFO%ENINI)
      ENDIF
!-----------------------------------------------------------------------
!  recalculate the real-space projection operators
!  if volume changed also recalculate reciprocal projection operators
!  and reset the cache for the phase-factor
!-----------------------------------------------------------------------
      IF (INFO%LREAL) THEN
! reset IRMAX, IRALLOC if required (no redistribution of GRIDS allowed)

        CALL REAL_OPTLAY(GRID,LATT_CUR,NONLR_S,.TRUE.,LREALLOCATE, IO%IU6,IO%IU0)
        IF (LREALLOCATE) THEN
           ! reallocate real space projectors
           CALL NONLR_DEALLOC(NONLR_S)
           CALL NONLR_ALLOC(NONLR_S)
        END IF
        CALL RSPHER(GRID,NONLR_S,LATT_CUR)

      ELSE
        IF (DYN%ISIF>=3) THEN
          CALL SPHER(GRID,NONL_S,P,WDES,LATT_CUR, 1)
        ENDIF
        CALL PHASE(WDES,NONL_S,0)
      ENDIF

      CALL RESETUP_FOCK( WDES, LATT_CUR)
!-----------------------------------------------------------------------
! recalculate projections and perform Gramm-Schmidt orthogonalization
!-----------------------------------------------------------------------
      CALL WVREAL(WDES,GRID,W) ! only for gamma some action
      CALL START_TIMING("G")
      CALL PROALL (GRID,LATT_CUR,NONLR_S,NONL_S,W)
      CALL ORTHCH(WDES,W, INFO%LOVERL, LMDIM,CQIJ)
      CALL REDIS_PW_OVER_BANDS(WDES, W)

      CALL KPAR_SYNC_ALL(WDES,W)

      CALL STOP_TIMING("G",TIU6,'ORTHCH')
!-----------------------------------------------------------------------
! set  INFO%LPOTOK to .F. this requires a recalculation of the local pot.
!-----------------------------------------------------------------------
      INFO%LPOTOK=.FALSE.
!-----------------------------------------------------------------------
! if prediction of wavefunctions was performed and
! diagonalization of sub-space-matrix is selected then
! )  POTLOK: calculate potential according to predicted charge-density
! )  SETDIJ: recalculate the energy of the augmentation charges
! )  then perform a sub-space-diagonal. and generate new Fermi-weights
! )  set INFO%LPOTOK to true because potential is OK
! )  recalculate total energy
! if charge density not constant during band minimization
! )  calculate charge-density according to new wavefunctions
!     and set LPOTOK to false (requires recalculation of loc. potential)
! in all other cases the predicted charge density is used in the next
!   step of ELM
!-----------------------------------------------------------------------
  pre_done: IF (PRED%IPRE>1) THEN
      io_begin
      IF (IO%IU0>=0) &
      WRITE(TIU0,*)'prediction of wavefunctions'

      WRITE(TIU6,2450) PRED%ALPHA,PRED%BETA
      io_end
 2450 FORMAT(' Prediction of Wavefunctions ALPHA=',F6.3,' BETA=',F6.3)

!   wavefunctions are not diagonal, so if they are written to the file
!   or if we do not diagonalize before the optimization
!   rotate them now
  pre_subrot: IF (.NOT.INFO%LDIAG .OR. INFO%LONESW .OR. &
     &     INFO%LSTOP .OR. (MOD(NSTEP,10)==0 &
     &     .AND. PRED%IWAVPR >= 2 .AND.  PRED%IWAVPR < 10) ) THEN

      CALL START_TIMING("G")
      IF (IO%IU0>=0) &
      WRITE(TIU0,*)'wavefunctions rotated'
      CALL POTLOK(GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES, &
                  INFO,P,T_INFO,E,LATT_CUR, &
                  CHTOT,CSTRF,CVTOT,DENCOR,SV, SOFT_TO_C,XCSIF)

      CALL POTLOK_METAGGA(KINEDEN, &
                  GRID,GRIDC,GRID_SOFT,WDES%COMM_INTER,WDES,INFO,P,T_INFO,E,LATT_CUR, &
                  CHDEN,CHTOT,DENCOR,CVTOT,SV,HAMILTONIAN%MUTOT,HAMILTONIAN%MU,SOFT_TO_C,XCSIF)

      CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
                  LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

      CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
         WDES%NCDIJ, LMDIM, CDIJ,  RHOLM, CRHODE, &
          E, LMETA =  .FALSE., LASPH =INFO%LASPH , LCOREL=.FALSE.)

      CALL UPDATE_CMBJ(GRIDC,T_INFO,LATT_CUR,IO%IU6)

      CALL STOP_TIMING("G",TIU6,'POTLOK')
      INFO%LPOTOK=.TRUE.

      IFLAG=3
      CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
          LMDIM,CDIJ,CQIJ, IFLAG,SV,T_INFO,P,IO%IU0,E%EXHF)

      CALL REDIS_PW_OVER_BANDS(WDES, W)

      SIGMA=ABS(KPOINTS%SIGMA)
      CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, &
               INFO%NELECT, INFO%NUP_DOWN, E%EENTROPY, EFERMI, SIGMA, .FALSE.,  &
               NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)

     ! for the selfconsistent update set W_F%CELTOT and TOTEN now
      IF (INFO%LONESW) W_F%CELTOT=W%CELTOT
      E%EBANDSTR= BANDSTRUCTURE_ENERGY(WDES, W)
      TOTEN=E%EBANDSTR+E%DENC+E%XCENC+E%TEWEN+E%PSCENC+E%EENTROPY+Ediel_sol

      CALL STOP_TIMING("G",IO%IU6,'EDDIAG')

  ENDIF pre_subrot

  ENDIF pre_done

    IF (INFO%LONESW) THEN

      CALL START_TIMING("G")

      CALL SET_CHARGE(W, WDES, INFO%LOVERL, &
           GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
           LATT_CUR, P, SYMM, T_INFO, &
           CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

      CALL SET_KINEDEN(GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM, &
           T_INFO%NIONS,W,WDES,KINEDEN)

      INFO%LPOTOK=.FALSE.

      CALL STOP_TIMING("G",IO%IU6,'CHARGE')

    ENDIF

  ELSE prepare_next_step
! in any case we have to call RSPHER at this point even if the ions do not move
! since the force routine overwrites the required arrays
      IF (INFO%LREAL) THEN
         CALL RSPHER(GRID,NONLR_S,LATT_CUR)
      ENDIF
  ENDIF prepare_next_step
!=======================================================================
!  update file WAVECAR if INFO%LSTOP = .TRUE.
!  or if wavefunctions on file TMPCAR are rotated
!=======================================================================
! after 10 steps rotate the wavefunctions on the file
      LTMP= MOD(NSTEP,10) == 0 .AND. PRED%IWAVPR >= 2 .AND. PRED%IWAVPR < 10
  wrtwave: IF ( IO%LWAVE .AND. ( INFO%LSTOP .OR. LTMP .OR. LHFCALC)  ) THEN

      CALL OUTWAV(IO, WDES, W, LATT_INI,  EFERMI)
!-----------------------------------------------------------------------
! rotate wavefunctions on file (gives a better prediction)
!------------------------------------------------------------------------
      IF (LTMP) THEN
        PRED%INIPRE=10

        CALL WAVPRE(GRIDC,P,PRED,T_INFO,W,WDES,LATT_CUR,IO%LOPEN, &
           CHTOT,RHOLM,N_MIX_PAW, CSTRF, LMDIM,CQIJ,INFO%LOVERL,IO%IU0)

        IF (IO%IU0>=0) &
             WRITE(TIU0,*)'wavefunctions on file TMPCAR rotated'
! and read in wavefunctions  (destroyed by WAVPRE)
        CALL OPENWAV(IO, COMM)
        CALL INWAV_FAST(IO, WDES, W, GRID, LATT_CUR, LATT_INI, INFO%ISTART,  EFERMI )
        CALL CLOSEWAV

        CALL PROALL (GRID,LATT_CUR,NONLR_S,NONL_S,W)
        CALL ORTHCH(WDES,W, INFO%LOVERL, LMDIM,CQIJ)
        CALL REDIS_PW_OVER_BANDS(WDES, W)
      ENDIF
   ENDIF wrtwave

!elphon_
   IF (ELPH_LWSWQ()) CALL ELPH_OVERLAP_W_WP(W,LATT_CUR,LMDIM,CQIJ,NONLR_S,NONL_S,IO)
!elphon_
!=======================================================================
! next electronic energy minimisation
      CALL STOP_TIMING("LOOP+",IO%IU6,XMLTAG='totalsc')

      IF (.NOT. INFO%LSTOP )  THEN
         CALL XML_CLOSE_TAG("calculation")
         CALL XML_FLUSH
      ENDIF
    END IF !LJ_ONLY

      ENDDO ion
!=======================================================================
! here we are at the end of the required number of timesteps
!=======================================================================
      io_begin
      IF (IO%LOPEN) CALL WFORCE(IO%IU6)
      io_end

    IF (.NOT. LJ_ONLY) THEN
      IF (LWRT_CHGFIT()) THEN
         ! subtract the overlapping atomic charge density from CHTOT
         CALL ATOMIC_CHARGES(T_INFO,LATT_CUR,P,GRIDC,CSTRF,1.0_q,-1.0_q,CHTOT)
         ! write CHGFIT
         CALL WRITE_CHARGE_RC(INFO,T_INFO,LATT_CUR,GRIDC,CHTOT,IO,99)
         ! restore CHTOT
         CALL ATOMIC_CHARGES(T_INFO,LATT_CUR,P,GRIDC,CSTRF,1.0_q,1.0_q,CHTOT)
      ENDIF
      ! calculate stockholder charge partitioning
      CALL STOCKHOLDER_ANALYSIS(INFO,T_INFO,LATT_CUR,P,GRIDC,NCDIJ,CHTOT,DENCOR,IO)
!=======================================================================
!  write out some additional information
!  create the File CHGCAR
!=======================================================================
      IF (IO%LCHARG) THEN
#ifdef MPI
         IF (WDES%COMM_KINTER%NODE_ME.EQ.1) THEN
#endif
         io_begin
         REWIND 18
         CALL OUTPOS(18,.FALSE.,INFO%SZNAM1,T_INFO,LATT_CUR%SCALE,LATT_CUR%A,.FALSE.,DYN%POSION)
         io_end
! if you uncomment the following lines the pseudo core charge density
! is added to the pseudo charge density
!         CALL FFT3D(CHTOT(1,1),GRIDC,1)
!         CALL RL_ADD(CHTOT(1,1),1._q/GRIDC%NPLWV,DENCOR(1),1._q/GRIDC%NPLWV,CHTOT(1,1),GRIDC)
!         CALL FFT3D(CHTOT(1,1),GRIDC,-1)
         CALL OUTCHG(GRIDC,18,.TRUE.,CHTOT)
         CALL WRT_RHO_PAW(P, T_INFO, INFO%LOVERL, RHOLM(:,1), GRIDC%COMM, 18 )
         DO ISP=2,WDES%NCDIJ
            do_io WRITE(18,'(5E20.12)') (T_INFO%ATOMOM(NI),NI=1,T_INFO%NIONS)
            CALL OUTCHG(GRIDC,18,.TRUE.,CHTOT(:,ISP))
            CALL WRT_RHO_PAW(P, T_INFO, INFO%LOVERL, RHOLM(:,ISP), GRIDC%COMM, 18 )
         ENDDO
         IF (IO%LOPEN) THEN
            do_io CALL REOPEN(18)
         ELSE
            do_io REWIND 18
         ENDIF
#ifdef MPI
         END IF
#endif
      ENDIF
!-----if we are interested in the total (local) potential write it here:
      IF (IO%LVTOT.OR.IO%LVHAR) THEN
#ifdef MPI
         IF (WDES%COMM_KINTER%NODE_ME.EQ.1) THEN
#endif
         io_begin
         IF (IO%LOPEN) OPEN(IO%IUVTOT,FILE='LOCPOT',STATUS='UNKNOWN')
         REWIND IO%IUVTOT
         CALL OUTPOS(IO%IUVTOT,.FALSE.,INFO%SZNAM1,T_INFO,LATT_CUR%SCALE,LATT_CUR%A,.FALSE.,DYN%POSION)
         io_end
! comment out the following line to add  exchange correlation
         IF (IO%LVHAR) CALL SET_LEXCH(-1)
         CALL POTLOK(GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES, &
                  INFO,P,T_INFO,E,LATT_CUR, &
                  CHTOT,CSTRF,CVTOT,DENCOR,SV, SOFT_TO_C,XCSIF)

         IF (.NOT. IO%LVHAR) THEN
           CALL POTLOK_METAGGA(KINEDEN, &
                 GRID,GRIDC,GRID_SOFT,WDES%COMM_INTER,WDES,INFO,P,T_INFO,E,LATT_CUR, &
                 CHDEN,CHTOT,DENCOR,CVTOT,SV,HAMILTONIAN%MUTOT,HAMILTONIAN%MU,SOFT_TO_C,XCSIF)
         ENDIF

         ! call the dipol routine without changing the potential
         IF ( DIP%IDIPCO >0 ) THEN
           DIP%LCOR_DIP=.FALSE.
           CALL CDIPOL_CHTOT_REC(GRIDC, LATT_CUR,P,T_INFO, &
               CHTOT,CSTRF,CVTOT, WDES%NCDIJ, INFO%NELECT, E%PSCENC )

           CALL WRITE_VACUUM_LEVEL(IO%IU6)
         ENDIF

!        CALL OUTPOT(GRIDC, IO%IUVTOT,.TRUE.,HAMILTONIAN%MUTOT)
         CALL OUTPOT(GRIDC, IO%IUVTOT,.TRUE.,CVTOT)
         DO ISP=2,WDES%NCDIJ
            do_io WRITE(IO%IUVTOT,'(5E20.12)') (T_INFO%ATOMOM(NI),NI=1,T_INFO%NIONS)
            CALL OUTPOT(GRIDC, IO%IUVTOT,.TRUE.,CVTOT(1,ISP))
         ENDDO
         IF (IO%LOPEN) THEN
            do_io CALL REOPEN(IO%IUVTOT)
         ELSE
            do_io REWIND IO%IUVTOT
         ENDIF
#ifdef MPI
         ENDIF
#endif
      ENDIF
!=======================================================================
!  Write out the Eigenvalues
!=======================================================================

#ifdef MPI
      IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
         CALL KPAR_SYNC_CELTOT(WDES,W)
      END IF
#endif

      io_begin
      DO NK=1,KPOINTS%NKPTS
        WRITE(22,*)
        WRITE(22,'(4E15.7)') WDES%VKPT(1,NK),WDES%VKPT(2,NK),WDES%VKPT(3,NK),KPOINTS%WTKPT(NK)
        DO N=1,WDES%NB_TOT
          IF (INFO%ISPIN==1) WRITE(22,852) N,REAL( W%CELTOT(N,NK,1) ,KIND=q), W%FERTOT(N,NK,1)
          IF (INFO%ISPIN==2) &
            WRITE(22,8852) N,REAL( W%CELTOT(N,NK,1) ,KIND=q) ,REAL( W%CELTOT(N,NK,INFO%ISPIN) ,KIND=q), W%FERTOT(N,NK,1), W%FERTOT(N,NK,INFO%ISPIN)
        ENDDO
      ENDDO
      IF (IO%LOPEN) CALL WFORCE(22)
      io_end
      CALL XML_EIGENVALUES(W%CELTOT, W%FERTOT, WDES%NB_TOT, KPOINTS%NKPTS, INFO%ISPIN)

  852 FORMAT(1X,I4,4X,F12.6,2X,F9.6)
 8852 FORMAT(1X,I4,4X,F12.6,2X,F12.6,2X,F9.6,2X,F9.6)
!=======================================================================
!  get current density and write output
!=======================================================================
      CALL CURRENT( W, GRID_SOFT, GRIDC, GRIDUS, C_TO_US, SOFT_TO_C, P, LATT_CUR, &
          HAMILTONIAN%AVEC, HAMILTONIAN%AVTOT, CHTOT, NONLR_S, NONL_S, &
          T_INFO, LMDIM, CRHODE, CDIJ, CQIJ, IRDMAX, IO%IU6, IO%IU0,  IO%NWRITE)

      CALL WRITE_ORBITALMAGOUT(IO%IU6)
      CALL XML_WRITE_ORBITALMAGOUT
!=======================================================================
!  calculate optical matrix elements
!=======================================================================

      IF (IO%LOPTICS.AND.KINTER==0) THEN
        CALL START_TIMING("G")
        ! VASP onboard optics
        CALL LR_OPTIC( &
             P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
             T_INFO,INFO,IO,KPOINTS,SYMM,GRID,GRID_SOFT, &
             GRIDC,GRIDUS,C_TO_US,SOFT_TO_C, &
             CHTOT,DENCOR,CVTOT,CSTRF, &
             CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM, &
             CHDEN,SV,LMDIM,IRDMAX,EFERMI,NEDOS, &
             LPOT= EXXOEP==0 .AND. INFO%INICHG/=4 .AND.  .NOT. USE_OEP_IN_GW() .AND. .NOT. LDO_METAGGA())
         ! potential must not be updated for OEP methods
         ! if potential was read from a file INFO%INICHG/=4
         ! or if a meta-GGA is used
        IF (NPAR ==1 .AND. KPAR==1 .AND. LPAW) THEN
           ! offboard optics by Juergen Furthmueller
           ALLOCATE(NABIJ(WDES%NB_TOT,WDES%NB_TOT))

           CALL CALC_NABIJ(NABIJ,W,WDES,P,KPOINTS,GRID_SOFT,LATT_CUR, &
                IO,INFO,T_INFO,COMM,IU0,55)
           DEALLOCATE(NABIJ)
        ENDIF
        CALL STOP_TIMING("G",IO%IU6,'OPTICS')
      ENDIF
!=======================================================================
!  calculate four orbital integrals
!=======================================================================
      CALL START_TIMING("G")

      CALL TWOELECTRON4O_MAIN( &
      P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
      T_INFO,DYN,INFO,IO,KPOINTS,SYMM,GRID,LMDIM )

      IF (FOURORBIT/=0) THEN
         ! Set all DFT exchange contributions to zero
         CALL PUSH_XC_TYPE(LEXCH,1._q,1._q,0._q,1._q,LDASCREEN)
         ! Initialize xc tables
         CALL SETUP_LDA_XC(1,IO%IU6,IO%IU0,IO%IDIOT)
         ! calculate correlation contributions on PW grid
         CALL POTLOK(GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES,  &
                     INFO,P,T_INFO,E,LATT_CUR, &
                     CHTOT,CSTRF,CVTOT,DENCOR,SV, SOFT_TO_C,XCSIF)
         CALL POTLOK_METAGGA(KINEDEN, &
                     GRID,GRIDC,GRID_SOFT,WDES%COMM_INTER,WDES,INFO,P,T_INFO,E,LATT_CUR, &
                     CHDEN,CHTOT,DENCOR,CVTOT,SV,HAMILTONIAN%MUTOT,HAMILTONIAN%MU,SOFT_TO_C,XCSIF)
         io_begin
            WRITE(*,'(A,F14.7)') 'LDA correlation energy (PW grid):',E%EXCG
         io_end
         ! Restore the original situation
         CALL POP_XC_TYPE
         ! reset the xc tables
         CALL SETUP_LDA_XC(1,IO%IU6,IO%IU0,IO%IDIOT)
         ! recalculate plane wave contributions
         CALL POTLOK(GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES,  &
                     INFO,P,T_INFO,E,LATT_CUR, &
                     CHTOT,CSTRF,CVTOT,DENCOR,SV, SOFT_TO_C,XCSIF)
         CALL POTLOK_METAGGA(KINEDEN, &
                     GRID,GRIDC,GRID_SOFT,WDES%COMM_INTER,WDES,INFO,P,T_INFO,E,LATT_CUR, &
                     CHDEN,CHTOT,DENCOR,CVTOT,SV,HAMILTONIAN%MUTOT,HAMILTONIAN%MU,SOFT_TO_C,XCSIF)
      ENDIF
      CALL STOP_TIMING("G",IO%IU6,'4ORBIT')

!=======================================================================
! MM: Stuff for Janos Angyan: write the AE-charge density to AECCAR2
!=======================================================================
      IF (LWRT_AECHG()) THEN
         CALL AUGCHG(WDES,GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
        &             LATT_CUR,P,T_INFO,SYMM, INFO%LOVERL, SOFT_TO_C,&
        &              LMDIM,CRHODE, CHTOT,CHDEN, IRDMAX,.TRUE.,.FALSE.)

#ifdef MPI
         IF (WDES%COMM_KINTER%NODE_ME.EQ.1) THEN
#endif
         OPEN(UNIT=99,FILE=DIR_APP(1:DIR_LEN)//'AECCAR2',STATUS='UNKNOWN')
         io_begin
         ! write header
         CALL OUTPOS(99,.FALSE.,INFO%SZNAM1,T_INFO,LATT_CUR%SCALE,LATT_CUR%A, &
        &             .FALSE., DYN%POSION)
         io_end
         ! write AE charge density
         CALL OUTCHG(GRIDC,99,.TRUE.,CHTOT)
         CLOSE(99)
         ! write Fourier transform of AE charge density
         IF (LWRTSTRF()) THEN
            OPEN(UNIT=99,FILE=DIR_APP(1:DIR_LEN)//'STRFAC',STATUS='UNKNOWN')
            CALL WRTSTRF(GRIDC,LATT_CUR,CHTOT,99)
            CLOSE(99)
         ENDIF
!        OPEN(UNIT=99,FILE=DIR_APP(1:DIR_LEN)//'RADCHG',STATUS='UNKNOWN')
!        CALL WRT_RHO_RAD(WDES,P,T_INFO,INFO%LOVERL,LMDIM,CRHODE,99)
!        CLOSE(99)

         END IF

         ! and restore the total charge density
         CALL DEPLE(WDES,GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
                  LATT_CUR,P,T_INFO,SYMM, INFO%LOVERL, SOFT_TO_C,&
                  LMDIM,CRHODE, CHTOT,CHDEN, IRDMAX)
#ifdef MPI
      ENDIF
#endif
! MM: end of addition

!=======================================================================
!  calculate ELF
!=======================================================================
      IF (IO%LELF) THEN
      ALLOCATE(CWORK(GRID_SOFT%MPLWV,WDES%NCDIJ))

      CALL ELF(GRID,GRID_SOFT,LATT_CUR,SYMM,NIOND, W,WDES,  &
               CHDEN,CWORK)

#ifdef MPI
     IF (WDES%COMM_KINTER%NODE_ME.EQ.1) THEN
#endif
! write ELF to file ELFCAR
      io_begin
      OPEN(UNIT=53,FILE='ELFCAR',STATUS='UNKNOWN')
      CALL OUTPOS(53,.FALSE.,INFO%SZNAM1,T_INFO,LATT_CUR%SCALE,LATT_CUR%A,.FALSE.,DYN%POSION)
      io_end

      DO ISP=1,WDES%NCDIJ
         CALL OUTCHG(GRID_SOFT,53,.FALSE.,CWORK(1,ISP))
      ENDDO

      DEALLOCATE(CWORK)

      do_io CLOSE(53)

      END IF
#ifdef MPI
      ENDIF
#endif
#ifdef MPI
      IF (WDES%COMM_KINTER%NODE_ME.EQ.1) THEN
#endif
      IF (WRITE_MOMENTS()) CALL WR_MOMENTS(GRID,LATT_CUR,P,T_INFO,W,WDES,.TRUE.)
      IF (WRITE_DENSITY()) CALL WR_PROJ_CHARG(GRID,P,LATT_CUR,T_INFO,WDES)
#ifdef MPI
      END IF
#endif
      IF (LCALC_ORBITAL_MOMENT().AND.WDES%LNONCOLLINEAR) CALL WRITE_ORBITAL_MOMENT(WDES,T_INFO%NIONS,IO%IU6)
      IF (WDES%LSORBIT) CALL WRITE_SPINORB_MATRIX_ELEMENTS(WDES,T_INFO,IO)

!=======================================================================
!-WAH Write augmentation charges and projectors for differential waves
!=======================================================================

      IF (INFO%LOVERL) THEN
       IF (STM(4) > 0) THEN
#ifdef MPI
         IF (WDES%COMM_KINTER%NODE_ME.EQ.1) THEN
#endif
          do_io WRITE(*,*) "writing IETS projectors"
          CALL WRT_IETS(LMDIM, WDES%NIONS, WDES%NRSPINORS, CQIJ, WDES, W)
          do_io WRITE(*,*) "IETS projectors  written, exiting"
         END IF
#ifdef MPI
       END IF
#endif
      END IF
!=======================================================================
!  possibly decompose into Bloch states of a given primitive cell
!=======================================================================
!     CALL  KPROJ(IO%IU5, IO%IU0, IO%IU6, GRID, NONL_S, T_INFO, SYMM, P, LATT_CUR, KPOINTS, W)
      CALL  KPROJ(IO%IU5, IO%IU0, IO%IU6, GRID, LATT_CUR, W, SYMM, CQIJ)

!=======================================================================
!  total DOS, calculate ion and lm decomposed occupancies and dos
!=======================================================================
!     IF (JOBPAR/=0  .AND. IO%LORBIT<10 .AND. NPAR /=1) THEN
!        CALL VTUTOR('W','partial DOS',RTUT,1, &
!             &                  ITUT,1,CDUM,1,LDUM,1,IO%IU6,IO%IDIOT)
!        CALL VTUTOR('W','partial DOS',RTUT,1, &
!             &                  ITUT,1,CDUM,1,LDUM,1,IO%IU0,IO%IDIOT)
!     ELSE IF (JOBPAR/=0 .OR. IO%LORBIT>=10  ) THEN
      IF (JOBPAR/=0 .OR. IO%LORBIT>=10  ) THEN
         DEALLOCATE(PAR,DOSPAR)

         IF (IO%LORBIT==1 .OR. IO%LORBIT==2 .OR. (IO%LORBIT>=11 .AND. IO%LORBIT<=14)) THEN
            LPAR=LMDIMP
         ELSE
            LPAR=LDIMP
         ENDIF

         ALLOCATE(PAR(WDES%NB_TOT,WDES%NKPTS,LPAR,T_INFO%NIONP,WDES%NCDIJ))

         IF (IO%LORBIT>=10) THEN
            CALL SPHPRO_FAST( &
                 GRID,LATT_CUR, P,T_INFO,W, WDES, 71,IO%IU6,&
                 INFO%LOVERL,LMDIM,CQIJ, LPAR, LDIMP, LMDIMP, .TRUE., IO%LORBIT,PAR, &
                 EFERMI, KPOINTS%EMIN, KPOINTS%EMAX)
! locproj_
                 CALL LPRJ_PROALL(W,WDES,GRID,P,CQIJ,LATT_CUR,T_INFO,INFO,IO%IU6,IO%IU0)
! locproj_
         ELSE
            CALL SPHPRO( &
                 GRID,LATT_CUR, P,T_INFO,W, WDES, 71,IO%IU6,&
                 INFO%LOVERL,LMDIM,CQIJ, LPAR, LDIMP, LMDIMP, LTRUNC, IO%LORBIT,PAR)
         ENDIF

         CALL CHGLOC(WDES%NB_TOT,KPOINTS%NKPTS,LPAR,T_INFO%NIONS,WDES%ISPIN,PAR,W%FERWE)

         !  get and write partial / projected DOS ...

         !  some compilers require to remove this statment
         DEALLOCATE(W%CW)         ! make space free so that DOSPAR can take this space
         ALLOCATE (DOSPAR(NEDOS,LPAR,T_INFO%NIONP,WDES%NCDIJ))
         SIGMA=ABS(KPOINTS%SIGMA)
         CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
              INFO%NUP_DOWN, E%EENTROPY, EFERMI, SIGMA, .FALSE.,  &
              NEDOS, LPAR, T_INFO%NIONP, DOS, DOSI, PAR, DOSPAR)
         io_begin
         DELTAE=(KPOINTS%EMAX-KPOINTS%EMIN)/(NEDOS-1)

         WRITE(16,'(2F16.8,I5,2F16.8)') KPOINTS%EMAX,KPOINTS%EMIN,NEDOS,EFERMI,1.0
         DO I=1,NEDOS
            EN=KPOINTS%EMIN+DELTAE*(I-1)
            WRITE(16,7062) EN,(DOS(I,ISP)/DYN%KBLOCK,ISP=1,WDES%ISPIN),(DOSI(I,ISP)/DYN%KBLOCK,ISP=1,WDES%ISPIN)
         ENDDO

         DO NI=1,T_INFO%NIONP
            WRITE(16,'(2F16.8,I5,2F16.8)') KPOINTS%EMAX,KPOINTS%EMIN,NEDOS,EFERMI,1.0
            DO I=1,NEDOS
               EN=KPOINTS%EMIN+DELTAE*(I-1)
               WRITE(16,'(3X,F8.3,36E12.4)') &
                    &            EN,((DOSPAR(I,LPRO,NI,ISP),ISP=1,WDES%NCDIJ),LPRO=1,LPAR)
            ENDDO
         ENDDO

         CALL XML_DOS(EFERMI, KPOINTS%EMIN, KPOINTS%EMAX, .TRUE., &
              DOS, DOSI, DOSPAR, NEDOS, LPAR, T_INFO%NIONP, WDES%NCDIJ)
         CALL XML_PROCAR(PAR, W%CELTOT, W%FERTOT, WDES%NB_TOT, WDES%NKPTS, LPAR ,T_INFO%NIONP,WDES%NCDIJ)
         io_end
      ELSE IF (DYN%IBRION/=1) THEN
         CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
              INFO%NUP_DOWN, E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE.,  &
              NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)

         io_begin
         DELTAE=(KPOINTS%EMAX-KPOINTS%EMIN)/(NEDOS-1)

         WRITE(16,'(2F16.8,I5,2F16.8)') KPOINTS%EMAX,KPOINTS%EMIN,NEDOS,EFERMI,1.0
         DO I=1,NEDOS
            EN=KPOINTS%EMIN+DELTAE*(I-1)
            WRITE(16,7062) EN,(DOS(I,ISP)/DYN%KBLOCK,ISP=1,WDES%ISPIN),(DOSI(I,ISP)/DYN%KBLOCK,ISP=1,WDES%ISPIN)
         ENDDO
         CALL XML_DOS(EFERMI, KPOINTS%EMIN, KPOINTS%EMAX, .FALSE., &
              DOS, DOSI, DOSPAR, NEDOS, LPAR, T_INFO%NIONP, WDES%NCDIJ)
         io_end
! locproj_
         CALL LPRJ_PROALL(W,WDES,GRID,P,CQIJ,LATT_CUR,T_INFO,INFO,IO%IU6,IO%IU0)
! locproj_

      ELSE
! locproj_
         CALL LPRJ_PROALL(W,WDES,GRID,P,CQIJ,LATT_CUR,T_INFO,INFO,IO%IU6,IO%IU0)
! locproj_
      ENDIF

      CALL XML_CLOSE_TAG("calculation")
    END IF !LJ_ONLY
! locproj_
!=======================================================================
!  project the wave functions onto local functions and write to
!  PROJCAR and LOCPROJ
!=======================================================================
      CALL LPRJ_WRITE(IO%IU6, IO%IU0, W)
      CALL LPRJ_DEALLOC_COVL
! locproj_


      CALL XML_TAG("structure","finalpos")
      CALL XML_CRYSTAL(LATT_CUR%A, LATT_CUR%B, LATT_CUR%OMEGA)
      CALL XML_POSITIONS(T_INFO%NIONS, DYN%POSION)
      IF (T_INFO%LSDYN) CALL XML_LSDYN(T_INFO%NIONS,T_INFO%LSFOR(1,1))
      IF (DYN%IBRION<=0 .AND. DYN%NSW>0 ) CALL XML_VEL(T_INFO%NIONS, DYN%VEL)
      IF (T_INFO%LSDYN) CALL XML_NOSE(DYN%SMASS)
      CALL XML_CLOSE_TAG("structure")
!=======================================================================
! breath a sigh of relief - you have finished
! this jump is just a jump to the END statement
!=======================================================================
 5100 CONTINUE

      IF (MIX%IMIX==4 .AND. INFO%IALGO.NE.-1) THEN
        CALL CLBROYD(MIX%IUBROY)
      ENDIF

      IF (INFO%LSOFT) THEN
         io_begin
         IF (IO%IU0>0) &
         WRITE(TIU0,*) 'deleting file STOPCAR'
         IF (IO%LOPEN) OPEN(99,FILE='STOPCAR',ERR=5111)
         CLOSE(99,STATUS='DELETE',ERR=5111)
 5111    CONTINUE
         io_end
      ENDIF
#ifdef CUDA_GPU
#ifdef USE_PINNED_MEMORY
      CALL nvpinnedfree(CHTOT_PTR)
      CALL nvpinnedfree(SV_PTR)
#endif
      CALL GPU_FINALIZE()
#endif
      CALL DUMP_ALLOCATE(IO%IU6)
      CALL DUMP_FINAL_TIMING(IO%IU6)
      CALL STOP_XML
#ifdef PROFILING
      ! started in INIT_PROFILING
      PROFILING_STOP('total_time')
      IF (IO%IU0>=0) THEN
        POUT=0
        PTHRESHOLD=0.0_q
        CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'PTHRESHOLD','=','#',';','F', &
        &            IDUM,PTHRESHOLD,CDUM,LDUM,CHARAC,N,1,IERR)
        IF ((IERR/=0).AND.(IERR/=3)) THEN
          WRITE(IU0,*)'Error reading item ''PTHRESHOLD'' from file INCAR.'
          STOP
        ENDIF
        IF (PTHRESHOLD/=0._q) POUT=POUT+1
        PROUTINE=''
        CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'PROUTINE','=','#',';','S', &
        &            IDUM,RDUM,CDUM,LDUM,PROUTINE,N,40,IERR)
        IF ((IERR/=0).AND.(IERR/=3)) THEN
         WRITE(IU0,*)'Error reading item ''PROUTINE'' from file INCAR.'
         STOP
        ENDIF
        CALL STRIP(PROUTINE,N,'A')
        CALL LOWER(PROUTINE)
        IF (PROUTINE/='') POUT=POUT+1
        PLEVEL=0
        CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'PLEVEL','=','#',';','I', &
        &            PLEVEL,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
        IF ((IERR/=0).AND.(IERR/=3)) THEN
          WRITE(IU0,*)'Error reading item ''PLEVEL'' from file INCAR.'
          STOP
        ENDIF
        IF (PLEVEL/=0) POUT=POUT+1
        PFLAT=.FALSE.
        CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'PFLAT','=','#',';','L', &
        &            IDUM,RDUM,CDUM,PFLAT,CHARAC,N,1,IERR)
        IF ((IERR/=0).AND.(IERR/=3)) THEN
          WRITE(IU0,*)'Error reading item ''PFLAT'' from file INCAR.'
          STOP
        ENDIF

        IF (PTHRESHOLD/=0.0_q) CALL PRINT_FULL_PROFILE(IO%IU6,PTHRESHOLD=PTHRESHOLD,PFLAT=PFLAT)
        IF (PROUTINE/='') CALL PRINT_ROUTINE_PROFILE(IO%IU6,PROUTINE=PROUTINE,PFLAT=PFLAT)
        IF (PLEVEL/=0) CALL PRINT_LEVEL_PROFILE(IO%IU6,PLEVEL=PLEVEL,PFLAT=PFLAT)
        IF (POUT==0)  CALL PRINT_FULL_PROFILE(IO%IU6,PFLAT=PFLAT)
      ENDIF
#endif
      CALLMPI_C(M_exit())

      CONTAINS

!**********************************************************************
!
!  internal subroutine to perform optimization of electronic
!  degrees of freedom
!
!**********************************************************************

      SUBROUTINE ELECTRONIC_OPTIMIZATION

      IF ( EXXOEP==1 ) THEN
         CALL ELMIN_LHF( &
          HAMILTONIAN,KINEDEN, &
          P,WDES,NONLR_S,NONL_S,W,W_F,W_G,LATT_CUR,LATT_INI, &
          T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI,CHF,CHAM,ECONV,XCSIF, &
          NSTEP,LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI,LDIMP,LMDIMP)
      ELSEIF ( EXXOEP>1) THEN
         CALL ELMIN_OEP( &
          HAMILTONIAN,KINEDEN, &
          P,WDES,NONLR_S,NONL_S,W,W_F,W_G,LATT_CUR,LATT_INI, &
          T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI,CHF,CHAM,ECONV,XCSIF, &
          NSTEP,LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI,LDIMP,LMDIMP)
      ELSEIF ( INFO%LONESW ) THEN
         CALL ELMIN_ALL( &
          HAMILTONIAN,KINEDEN, &
          P,WDES,NONLR_S,NONL_S,W,W_F,W_G,LATT_CUR,LATT_INI, &
          T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI,CHF,CHAM,ECONV,XCSIF, &
          NSTEP,LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI,LDIMP,LMDIMP)

          ! convergence reached switch back to DIIS
          IF (.NOT. INFO%LABORT .AND. INFO%LONESW_AUTO) THEN
             INFO%LONESW=.FALSE.
          ENDIF

      ELSE
         CALL ELMIN( &
          HAMILTONIAN,KINEDEN, &
          P,WDES,NONLR_S,NONL_S,W,W_F,W_G,LATT_CUR,LATT_INI, &
          T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI,CHF,CHAM,ECONV,XCSIF, &
          NSTEP,LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI,LDIMP,LMDIMP)


         IF (INFO%LABORT .AND. INFO%LONESW_AUTO) THEN
          INFO%LONESW=.TRUE.
          W_F%CELTOT=W%CELTOT ! copy current weights to W_F%CELTOT
          MIX%HARD_RESET=.TRUE.
         ! no convergence
         ! try to switch to CG for all bands
          CALL ELMIN_ALL( &
          HAMILTONIAN,KINEDEN, &
          P,WDES,NONLR_S,NONL_S,W,W_F,W_G,LATT_CUR,LATT_INI, &
          T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI,CHF,CHAM,ECONV,XCSIF, &
          NSTEP,LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI,LDIMP,LMDIMP)

         ENDIF
      ENDIF

      END SUBROUTINE ELECTRONIC_OPTIMIZATION

!**********************************************************************
!
!  internal subroutine to perform RPA calculations
!  TODO: there are still logical flaws in this routine
!  - LRPAFORCE will not work if LNBANDS is .TRUE.
!    since STORE_WDES_GROUNDSTATE is not called
!
!**********************************************************************

      SUBROUTINE DO_RPA
        CALL START_TIMING("LOOPX")

            CALL CALCULATE_XI( &
                KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
                GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C,E, &
                CHTOT, CHTOTL, DENCOR, CVTOT, CSTRF, IRDMAX, &
                T_INFO, DYN, INFO, IO, KPOINTS, SYMM, MIX, &
                LMDIM, CQIJ, CDIJ, CRHODE, N_MIX_PAW, RHOLM, RHOLM_LAST, CHDEN, SV, &
                EFERMI, NEDOS, DOS, DOSI)
         CALL PEAD_RESETUP_WDES(WDES,GRID,KPOINTS,LATT_CUR,LATT_INI,IO)

#ifdef VASP2WANNIER90
         CALL KPAR_SYNC_ALL(WDES,W)
         CALL MLWF_WANNIER90(WDES,W,P,CQIJ,T_INFO,LATT_CUR,INFO,IO)
#endif
         CALL SEPERATOR_TIMING(IO%IU6)
         CALL STOP_TIMING("LOOPX",IO%IU6,'GWTOTAL')

         IF (LRPAFORCE) THEN
            CALL LR_SKELETON( &
                 KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
                 T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
                 GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
                 CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
                 CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
                 CHDEN,SV,DOS,DOSI,CHAM, &
                 8,LMDIM,IRDMAX,NEDOS, &
                 TOTEN,EFERMI, TIFOR)

            CALL WRITE_RPA_FORCE( TIFOR, T_INFO, LATT_CUR, SYMM, DYN%POSIOC, DYN%IBRION, IO)
            CALL SEPERATOR_TIMING(IO%IU6)
            CALL STOP_TIMING("LOOPX",IO%IU6,'LINEAR')
         ENDIF

         ! now we restore the HF xc-correlation functional
         CALL PUSH_XC_TYPE_FOR_GW
         IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
            CALL SETUP_LDA_XC(2,-1,-1,IO%IDIOT)
         ELSE
            CALL SETUP_LDA_XC(1,-1,-1,IO%IDIOT)
         ENDIF
         ! now update the PAW one-center terms to current functional
         CALL SET_PAW_ATOM_POT( P , T_INFO, WDES%LOVERL, LMDIM, INFO%EALLAT, INFO%LMETAGGA, IO%IU6 )
         ! and restore the convergence corrections
         DO NK=1,WDES%NKPTS
            FSG_STORE(NK)=SET_FSG(GRIDHF, LATT_CUR, NK)
         ENDDO
    END SUBROUTINE DO_RPA


!**********************************************************************
!
!  internal subroutine to perform a few performance tests
!  Output is written to IUT
!
!**********************************************************************

    SUBROUTINE PERFORMANCE_TEST
      io_begin
      IUT=IO%IU0

      IF (IUT>0) WRITE(IUT,5001)
 5001 FORMAT(/ &
     & ' All results refer to a run over all bands and one k-point'/ &
     & ' VNLACC   non local part of H'/ &
     & ' PROJ     calculate projection of all bands (contains FFTWAV)'/ &
     & ' RACC     non local part of H in real space (contains FFTEXT)'/ &
     & ' RPRO     calculate projection of all bands in real space '/ &
     & '          both calls contain on FFT (to be subtracted)'/ &
     & ' FFTWAV   FFT of a wavefunction to real space'/ &
     & ' FFTEXT   FFT to real space'/ &
     & ' ECCP     internal information only (subtract FFTWAV)'/ &
     & ' POTLOK   update of local potential (including one FFT)'/ &
     & ' SETDIJ   calculate stregth of US PP'/ &
     & ' ORTHCH   gramm-schmidt orth.  applying Choleski decomp.'/ &
     & ' LINCOM   unitary transformation of wavefunctions'/ &
     & ' LINUP    upper triangle transformation of wavefunctions'/ &
     & ' ORTHON   orthogonalisation of one band to all others')
       io_end


! set the wavefunction descriptor
      ISP=1
      NK=1
      CALL SETWDES(WDES,WDES1,NK)

      INFO%ISPIN=1
      INFO%RSPIN=2

      NPL=WDES%NPLWKP(NK)
      ALLOCATE(CWORK1(GRID%MPLWV),CWORK2(GRID%MPLWV),CPROTM(LMDIM*NIOND))
      W1%CR=>CWORK1

      CALLMPI( MPI_barrier( WDES%COMM%MPI_COMM, ierror ))
#ifdef MPI
      IF (WDES%COMM_INTER%NCPU/=1) THEN
      CALL START_TIMING("G")
      DO I=1,10
         CALL REDIS_PW(WDES1, WDES%NBANDS, W%CW   (1,1,NK,1))
      ENDDO
      CALL STOP_TIMING("G",IUT,'10xRED')
      ENDIF
#endif

      IF(INFO%TURBO==0)THEN
         CALL START_TIMING("G")
         CALL STUFAK(GRIDC,T_INFO,CSTRF)
         CALL STOP_TIMING("G",IUT,'STUFAK')
      ENDIF

      CALL POTLOK(GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES, &
                  INFO,P,T_INFO,E,LATT_CUR, &
                  CHTOT,CSTRF,CVTOT,DENCOR,SV, SOFT_TO_C,XCSIF)
      CALL POTLOK_METAGGA(KINEDEN, &
                  GRID,GRIDC,GRID_SOFT,WDES%COMM_INTER,WDES,INFO,P,T_INFO,E,LATT_CUR, &
                  CHDEN,CHTOT,DENCOR,CVTOT,SV,HAMILTONIAN%MUTOT,HAMILTONIAN%MU,SOFT_TO_C,XCSIF)
      CALL STOP_TIMING("G",IUT,'POTLOK ')

      CALL FORLOC(GRIDC,P,T_INFO,LATT_CUR, CHTOT,TIFOR)
      CALL STOP_TIMING("G",IUT,'FORLOC')

      CALL SOFT_CHARGE(GRID,GRID_SOFT,W,WDES, CHDEN(1,1))
      CALL STOP_TIMING("G",IUT,'CHSP  ')


      CALL DEPLE(WDES,GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
               LATT_CUR,P,T_INFO,SYMM, INFO%LOVERL, SOFT_TO_C,&
               LMDIM,CRHODE, CHTOT,CHDEN, IRDMAX)

      CALL SET_RHO_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
           CRHODE, RHOLM)
      CALL STOP_TIMING("G",IUT,'DEPLE ')

      IF (INFO%LREAL) THEN

        CALL RSPHER(GRID,NONLR_S,LATT_CUR)
        CALL STOP_TIMING("G",IUT,'RSPHER')
        CWORK2=0
        CALL START_TIMING("G")
        CALL RACCT(NONLR_S,WDES,W,GRID,CDIJ,CQIJ, ISP, LMDIM, NK)
        CALL STOP_TIMING("G",IUT,'RACC')

      ELSE

        CALL PHASE(WDES,NONL_S,NK)
        NPL=WDES%NPLWKP(NK)
        CALL START_TIMING("G")
        DO N=1,WDES%NBANDS
          EVALUE=W%CELEN(N,1,1)
          CALL SETWAV(W,W1,WDES1,N,1)  ! allocation for W1%CR done above
          CALL VNLACC(NONL_S,W1,CDIJ,CQIJ, ISP, EVALUE,  CWORK2)
        ENDDO
        CALL STOP_TIMING("G",IUT,'VNLACC')
      ENDIF


      CALL START_TIMING("G")
      IF (INFO%LREAL) THEN
        CALL START_TIMING("G")
        CALL RPRO(NONLR_S,WDES,W,GRID,NK)
        CALL STOP_TIMING("G",IUT,'RPRO  ')
      ELSE
        CALL PROJ(NONL_S,WDES,W,NK)
        CALL STOP_TIMING("G",IUT,'PROJ  ')
      ENDIF

      CALL START_TIMING("G")
      DO  N=1,WDES%NBANDS
        CALL FFTWAV(NPL,WDES%NINDPW(1,NK),CWORK1,W%CPTWFP(1,N,NK,1),GRID)
      ENDDO
      CALL STOP_TIMING("G",IUT,'FFTWAV')

      DO N=1,WDES%NBANDS
        CALL INIDAT(GRID%RC%NP,CWORK1)
      ENDDO
      CALL STOP_TIMING("G",IUT,'FFTINI')


      DO N=1,WDES%NBANDS
        CALL INIDAT(GRID%RC%NP,CWORK1)
        CALL FFT3D(CWORK1,GRID,1)
      ENDDO
      CALL STOP_TIMING("G",IUT,'FFT3DF')

      DO N=1,WDES%NBANDS
        CALL INIDAT(GRID%RC%NP,CWORK1)
        CALL FFT3D(CWORK1,GRID,1)
        CALL FFT3D(CWORK1,GRID,-1)
      ENDDO
      CALL STOP_TIMING("G",IUT,'FFTFB ')

      DO N=1,WDES%NBANDS
        CALL INIDAT(GRID%RL%NP,CWORK1)
        CALL FFTEXT(NPL,WDES%NINDPW(1,NK),CWORK1,CWORK2,GRID,.FALSE.)
      ENDDO

      CALL STOP_TIMING("G",IUT,'FFTEXT ')

      DO N=1,WDES%NBANDS
          CALL FFTWAV(NPL,WDES%NINDPW(1,NK),CWORK1,W%CPTWFP(1,N,NK,1),GRID)
          CALL SETWAV(W,W1,WDES1,N,1)  ! allocation for W1%CR done above
          IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
             CALL ECCP_TAU(WDES1,W1,W1,LMDIM,CDIJ,GRID,SV,LATT_CUR,HAMILTONIAN%MU,W1%CELEN)
          ELSE
             CALL ECCP(WDES1,W1,W1,LMDIM,CDIJ,GRID,SV, W1%CELEN)
          ENDIF
      ENDDO

      CALL STOP_TIMING("G",IUT,'ECCP ')

      CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
                  LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

      CALL STOP_TIMING("G",IUT,'SETDIJ ')

      CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
         INFO%ISPIN, LMDIM, CDIJ,  RHOLM, CRHODE, &
          E,  LMETA =  .FALSE., LASPH =INFO%LASPH , LCOREL=.FALSE.)

      CALL UPDATE_CMBJ(GRIDC,T_INFO,LATT_CUR,IO%IU6)

      CALL STOP_TIMING("G",IUT,'SETPAW ')

      CALL ORTHCH(WDES,W, INFO%LOVERL, LMDIM,CQIJ)
      CALL REDIS_PW_OVER_BANDS(WDES, W)
      CALL STOP_TIMING("G",IUT,'ORTHCH')

      IF (INFO%LDIAG) THEN
        IFLAG=3
      ELSE
        IFLAG=4
      ENDIF
      CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
          LMDIM,CDIJ,CQIJ, IFLAG,SV,T_INFO,P,IO%IU0,E%EXHF)

      CALL REDIS_PW_OVER_BANDS(WDES, W)
      CALL STOP_TIMING("G",IUT,'EDDIAG')

! avoid that MATMUL is too clever
      ALLOCATE(CMAT(WDES%NB_TOT,WDES%NB_TOT))
      DO N1=1,WDES%NB_TOT
      DO N2=1,WDES%NB_TOT
        IF (N1==N2)  THEN
         CMAT(N1,N2)=0.99999_q
        ELSE
         CMAT(N1,N2)=EXP((0.7_q,0.5_q)/100)
       ENDIF
      ENDDO
      ENDDO

      NPRO= WDES%NPRO
      CALLMPI( SET_NPL_NPRO(WDES1, NPL, NPRO))
#ifdef MPI
      NCPU=WDES%COMM_INTER%NCPU ! number of procs involved in band dis.
#else
      NCPU=1
#endif
      NRPLWV_RED=WDES%NRPLWV/NCPU
      NPROD_RED =WDES%NPROD /NCPU

      CALL START_TIMING("G")
      CALL LINCOM('F',W%CPTWFP(:,:,NK,1),W%CPROJ(:,:,NK,1),CMAT(1,1), &
       WDES%NB_TOT,WDES%NB_TOT,NPL,0,NRPLWV_RED,NPROD_RED,WDES%NB_TOT, &
       W%CPTWFP(:,:,NK,1),W%CPROJ(:,:,NK,1))
      CALL STOP_TIMING("G",IUT,'LINCOM')

      CALL LINCOM('F',W%CPTWFP(:,:,NK,1),W%CPROJ(:,:,NK,1),CMAT(1,1), &
       WDES%NB_TOT,WDES%NB_TOT,NPL,NPRO,NRPLWV_RED,NPROD_RED,WDES%NB_TOT, &
       W%CPTWFP(:,:,NK,1),W%CPROJ(:,:,NK,1))
      CALL STOP_TIMING("G",IUT,'LINCOM')

      CALL LINCOM('U',W%CPTWFP(:,:,NK,1),W%CPROJ(:,:,NK,1),CMAT(1,1), &
       WDES%NB_TOT,WDES%NB_TOT,NPL,NPRO,NRPLWV_RED,NPROD_RED,WDES%NB_TOT, &
       W%CPTWFP(:,:,NK,1),W%CPROJ(:,:,NK,1))
      CALL STOP_TIMING("G",IUT,'LINUP ')
    END SUBROUTINE PERFORMANCE_TEST

  END PROGRAM
