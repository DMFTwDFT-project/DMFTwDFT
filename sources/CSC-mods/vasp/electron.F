#include "symbol.inc"
!**********************************************************************
! RCS:  $Id: electron.F,v 1.12 2003/06/27 13:22:15 kresse Exp kresse $
!
! subroutine for performing electronic minimization in VASP
!
!**********************************************************************

      SUBROUTINE ELMIN( &
          HAMILTONIAN,KINEDEN, &
          P,WDES,NONLR_S,NONL_S,W,W_F,W_G,LATT_CUR,LATT_INI, &
          T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI,CHF,CHAM,ECONV,XCSIF, &
          NSTEP,LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI,LDIMP,LMDIMP)

      USE mlwf
      USE prec
      USE hamil_high
      USE morbitalmag
      USE pseudo
      USE lattice
      USE steep
      USE us
      USE pot
      USE force
      USE fileio
      USE nonl_high
      USE rmm_diis
      USE david
      USE david_inner
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
      USE pawm
      USE cl
      USE vaspxml
      USE mdipol
      USE pawfock
      USE Constrained_M_modular
      USE ini
      USE LDAPLUSU_MODULE
      USE core_rel
      USE pp_data
      USE gw_model
      USE meta
      USE locproj
! solvation__
      USE solvation
! solvation__
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)
!=======================================================================
!  structures
!=======================================================================
      TYPE (ham_handle)  HAMILTONIAN
      TYPE (tau_handle)  KINEDEN
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (wavedes)     WDES
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (wavespin)    W          ! wavefunction
      TYPE (wavespin)    W_F        ! wavefunction for all bands simultaneous
      TYPE (wavespin)    W_G        ! same as above
      TYPE (latt)        LATT_CUR
      TYPE (dynamics)    DYN
      TYPE (info_struct) INFO
      TYPE (in_struct)   IO
      TYPE (mixing)      MIX
      TYPE (kpoints_struct) KPOINTS
      TYPE (symmetry)    SYMM
      TYPE (grid_3d)     GRID       ! grid for wavefunctions
      TYPE (grid_3d)     GRID_SOFT  ! grid for soft chargedensity
      TYPE (grid_3d)     GRIDC      ! grid for potentials/charge
      TYPE (grid_3d)     GRIDUS     ! temporary grid in us.F
      TYPE (grid_3d)     GRIDB      ! Broyden grid
      TYPE (transit)     B_TO_C     ! index table between GRIDB and GRIDC
      TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
      TYPE (transit)     SOFT_TO_C  ! index table between GRID_SOFT and GRIDC
      TYPE (energy)      E
      TYPE (latt)        LATT_INI

      INTEGER NSTEP,LMDIM,IRDMAX,NEDOS
      REAL(q) :: TOTEN,EFERMI

      COMPLEX(q)  CHTOT(GRIDC%MPLWV,WDES%NCDIJ) ! charge-density in real / reciprocal space
      COMPLEX(q)  CHTOTL(GRIDC%MPLWV,WDES%NCDIJ)! old charge-density
      RGRID       DENCOR(GRIDC%RL%NP)           ! partial core
      COMPLEX(q)  CVTOT(GRIDC%MPLWV,WDES%NCDIJ) ! local potential
      COMPLEX(q)  CSTRF(GRIDC%MPLWV,T_INFO%NTYP)! structure factor

!   augmentation related quantities
      OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
               CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
               CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
!  paw sphere charge density
      INTEGER N_MIX_PAW
      REAL(q)  RHOLM(N_MIX_PAW,WDES%NCDIJ),RHOLM_LAST(N_MIX_PAW,WDES%NCDIJ)
!  charge-density and potential on soft grid
      COMPLEX(q)  CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
      RGRID       SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ)
!  density of states
      REAL(q)    DOS(NEDOS,WDES%ISPIN),DOSI(NEDOS,WDES%ISPIN)
!  Hamiltonian
      GDEF       CHF(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN), &
                 CHAM(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)
      REAL(q) :: XCSIF(3,3)

! local
      REAL(q) :: TOTENL=0
      REAL(q) :: DESUM1,DESUM(INFO%NELM)
      INTEGER :: IONODE, NODE_ME
!  needed temporary for aspherical GGA calculation
      OVERLAP,ALLOCATABLE ::  CDIJ_TMP(:,:,:,:)
! local l-projected wavefunction characters (not really used here)
      REAL(q)    PAR(1,1,1,1,WDES%NCDIJ),DOSPAR(1,1,1,WDES%NCDIJ)

      REAL(q), EXTERNAL :: RHO0
      INTEGER N,ISP,ICONJU,IROT,ICEL,I,II,IRDMAA, &
              IERR,IDUM,IFLAG,ICOUEV,ICOUEV2,NN,NORDER,IERRBR,L,LP, &
              NCOUNT
      REAL(q) BTRIAL,RDUM,RMS,ORT,TOTEN2,RMS2,RMST, &
              WEIGHT,BETATO,DESUM2,RMSC,RMSP
      REAL(q) RHOAUG(WDES%NCDIJ),RHOTOT(WDES%NCDIJ)
      COMPLEX(q) CDUM
      CHARACTER (LEN=1) CHARAC
      LOGICAL LDELAY, LABORT_WITHOUT_CONV
! parameters for FAST_SPHPRO
      INTEGER :: LDIMP,LMDIMP
      REAL(q) :: TIFOR(3,T_INFO%NIONS),tot_n
      INTEGER NB, NK, NP, NBG
      COMPLEX(q),allocatable :: WaveW(:,:,:),WaveW_tot(:,:,:)
      COMPLEX(q),allocatable :: WProj(:,:,:),WProj_tot(:,:,:)
      REAL(q),allocatable :: FERWE(:,:)

      IONODE=0
      NODE_ME=0
#ifdef MPI
      IONODE  = WDES%COMM%IONODE
      NODE_ME = WDES%COMM%NODE_ME
#endif

      NELM=INFO%NELM
      NORDER=0 ; IF (KPOINTS%ISMEAR>=0) NORDER=KPOINTS%ISMEAR
      ! to make timing more sensefull syncronize now
      CALLMPI( MPI_barrier( WDES%COMM%MPI_COMM, ierror ))
      CALL START_TIMING("LOOP")
      CALL START_TIMING("G")

     io_begin
      IF (IO%IU0>=0) WRITE(IO%IU0,142)
      WRITE(17,142)
  142 FORMAT('       N       E                     dE             ' &
            ,'d eps       ncg     rms          rms(c)')
     io_end

      DESUM1=0
      INFO%LMIX=.FALSE.

 130  FORMAT (5X, //, &
     &'----------------------------------------------------', &
     &'----------------------------------------------------'//)

 140  FORMAT (5X, //, &
     &'--------------------------------------- Iteration ', &
     &I6,'(',I4,')  ---------------------------------------'//)
     DWRITE0 'electron entered'

      CALL DIPOL_RESET()

      CALL SPAWN_PP(T_INFO,SYMM,WDES,P,IO)
      CALL INIT_CORE_REL(WDES,CRHODE,IO%IU0,IO%IU5,IO%IU6)

!-----------------------------------------------------------------------

      IF (MIX%IMIX/=0 .AND. .NOT. INFO%LCHCON  .AND. MIX%MIXFIRST) THEN
         INFO%LPOTOK=.FALSE.
         ! copy current density to last density
         DO ISP=1,WDES%NCDIJ
            CALL RC_ADD(CHTOT(1,ISP),1.0_q,CHTOT(1,ISP),0.0_q,CHTOTL(1,ISP),GRIDC)
         ENDDO

         IF (LDO_METAGGA().AND.LMIX_TAU()) THEN
            DO ISP=1,WDES%NCDIJ
               CALL RC_ADD(KINEDEN%TAU(1,ISP),1.0_q,KINEDEN%TAU(1,ISP),0.0_q,KINEDEN%TAUL(1,ISP),GRIDC)
            ENDDO
         ENDIF

         RHOLM_LAST=RHOLM

         INFO%LPOTOK=.FALSE.

         ! add the available GAMMA file
         IF (MOD(INFO%ICHARG,10)==5) THEN
            CALL ADD_GAMMA_FROM_FILE( WDES, W, KPOINTS, INFO%NELECT, INFO%NUP_DOWN, .TRUE., IO )
         ENDIF

         ! recalculate charge density and  kinetic energy density
         CALL SET_CHARGE(W, WDES, INFO%LOVERL, &
           GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
           LATT_CUR, P, SYMM, T_INFO, &
           CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

         CALL SET_KINEDEN(GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM, &
           T_INFO%NIONS,W,WDES,KINEDEN)

         CALL STOP_TIMING("G",IO%IU6,"CHARGE")

         INFO%LMIX=.TRUE.
! update one center charge density
         IF (LCORREL()) THEN
            CALL SET_RHO_ONE_CENTRE(CRHODE,RHO_ONE_CENTRE)
         ENDIF

         IF (MIX%IMIX==4) THEN
!  broyden mixing ... :
            IF (LCORREL()) THEN
               CALL BRMIX(KINEDEN,GRIDB,GRIDC,IO,MIX,B_TO_C, &
                    (2*GRIDC%MPLWV),CHTOT,CHTOTL,WDES%NCDIJ,LATT_CUR%B, &
                    LATT_CUR%OMEGA, N_MIX_PAW, RHOLM, RHOLM_LAST, &
                    RMST,RMSC,RMSP,WEIGHT,.TRUE.,IERRBR, &
                    NMIX_ONE_CENTRE,RHO_ONE_CENTRE,RHO_ONE_CENTRE_LAST)
               MIX%LRESET=.FALSE.
            ELSE
               CALL BRMIX(KINEDEN,GRIDB,GRIDC,IO,MIX,B_TO_C, &
                    (2*GRIDC%MPLWV),CHTOT,CHTOTL,WDES%NCDIJ,LATT_CUR%B, &
                    LATT_CUR%OMEGA, N_MIX_PAW, RHOLM, RHOLM_LAST, &
                    RMST,RMSC,RMSP,WEIGHT,.TRUE.,IERRBR)
               MIX%LRESET=.FALSE.
            ENDIF
         ELSE
!  simple mixing ... :
            RMST=0
            CALL MIX_SIMPLE(GRIDC,MIX,WDES%NCDIJ, CHTOT,CHTOTL, &
                 N_MIX_PAW, RHOLM, RHOLM_LAST, LATT_CUR%B, LATT_CUR%OMEGA, RMST)
! test
            IF (LCORREL()) THEN
               RHO_ONE_CENTRE(:)=RHO_ONE_CENTRE_LAST(:)+0.2*(RHO_ONE_CENTRE(:)-RHO_ONE_CENTRE_LAST(:))
               RHO_ONE_CENTRE_LAST(:)=RHO_ONE_CENTRE(:)
            ENDIF
! test
         ENDIF
         CALL STOP_TIMING("G",IO%IU6,"MIXING")
      ENDIF
!=======================================================================
      electron: DO N=1,NELM

      CALL XML_TAG("scstep")

!======================================================================
      io_begin
      WRITE(IO%IU6,140) NSTEP,N
      io_end
!=======================================================================
! if recalculation of total lokal potential is necessary (INFO%LPOTOK=.F.)
! call POTLOK: the subroutine calculates
! ) the hartree potential from the electronic  charge density
! ) the exchange correlation potential
! ) and the total lokal potential
!  in addition all double counting correction and forces are calculated
! &
! call SETDIJ
! calculates the Integral of the depletion charges * local potential
! and sets CDIJ
!=======================================================================

      CALL WVREAL(WDES,GRID,W) ! only for gamma some action

      IF (.NOT. INFO%LPOTOK) THEN

      ! core relaxation and repseudization
      IF (LCORREL()) THEN
!        CALL PW_TO_RADIAL(WDES,GRID_SOFT,CHDEN(:,1),LATT_CUR,T_INFO)
         CALL GET_AVERAGEPOT_PW( GRIDC, LATT_CUR, IRDMAX,  &
        &   T_INFO, P, WDES%NCDIJ, CVTOT, MAX(INFO%ENAUG,INFO%ENMAX), IO%IU6)
         CALL CORREL(RHO_ONE_CENTRE)
      ENDIF
#define no_update_potential
#ifdef no_update_potential
      CALL POTLOK(GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES, &
                  INFO,P,T_INFO,E,LATT_CUR, &
                  CHTOT,CSTRF,CVTOT,DENCOR,SV, SOFT_TO_C,XCSIF)

      CALL POTLOK_METAGGA(KINEDEN, &
                  GRID,GRIDC,GRID_SOFT,WDES%COMM_INTER,WDES,INFO,P,T_INFO,E,LATT_CUR, &
                  CHDEN,CHTOT,DENCOR,CVTOT,SV,HAMILTONIAN%MUTOT,HAMILTONIAN%MU,SOFT_TO_C,XCSIF)

      CALL VECTORPOT(GRID, GRIDC, GRID_SOFT, SOFT_TO_C,  WDES%COMM_INTER, &
                 LATT_CUR, T_INFO%POSION, HAMILTONIAN%AVEC, HAMILTONIAN%AVTOT)

      CALL STOP_TIMING("G",IO%IU6,"POTLOK")
      DWRITE0 'potlok is ok'


      CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
                  LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

      CALL SETDIJ_AVEC(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
                  LMDIM,CDIJ,HAMILTONIAN%AVTOT, NONLR_S, NONL_S, IRDMAX)

      CALL SET_DD_MAGATOM(WDES, T_INFO, P, LMDIM, CDIJ)

      CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
         WDES%NCDIJ, LMDIM, CDIJ(1,1,1,1),  RHOLM, CRHODE(1,1,1,1), &
         E,  LMETA=.FALSE., LASPH=INFO%LASPH, LCOREL= .FALSE.  )

      CALL UPDATE_CMBJ(GRIDC,T_INFO,LATT_CUR,IO%IU6)

      CALL STOP_TIMING("G",IO%IU6,"SETDIJ")
      DWRITE0 'setdij is ok'
#else
      CALL  UPDATE_POTENTIAL(  &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,LATT_CUR, &
          T_INFO,INFO,IO, &
          GRID,GRID_SOFT,GRIDC,GRIDUS,C_TO_US,SOFT_TO_C, &
          CHTOT,DENCOR,CVTOT,CSTRF, &
          LMDIM,IRDMAX,CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,CHDEN,SV)
#endif

      IF (USELDApU()) CALL LDAPLUSU_PRINTOCC(WDES,T_INFO%NIONS,T_INFO%ITYP,IO%IU6)
!remove
      INFO%LPOTOK=.TRUE.
      ENDIF

      IF (LCORREL()) THEN
         CALL ORTHCH(WDES,W, INFO%LOVERL, LMDIM,CQIJ)
         DWRITE0 'orthch is ok'
      ENDIF



!======================== SUBROUTINE EDDSPX ============================
!
! these subroutines improve the electronic degrees of freedom
! using band by band schemes
! the Harris functional is used for the calculation
! of the total (free) energy so
! E  =  Tr[ H rho ] - d.c. (from input potential)
!
!=======================================================================
      DESUM1=0
      RMS   =0
      ICOUEV=0

      LDELAY=.FALSE.
      ! if Davidson and RMM are selected, use Davidsons algorithm during
      ! delay phase
      IF (INFO%LRMM  .AND. INFO%LDAVID .AND. (N <= ABS(INFO%NELMDL) .OR. N==1)) LDELAY=.TRUE.
      ! if LDELAY is set, subspace rotation and orthogonalisations can be bypassed
      ! since they are done by the Davidson algorithm

!
! sub space rotation before eigenvalue optimization
!
      IF (INFO%LPDIAG .AND. .NOT. LDELAY ) THEN

        IF (INFO%LDIAG) THEN
           IFLAG=3    ! exact diagonalization
        ELSE
           IFLAG=4    ! using Loewdin perturbation theory
        ENDIF
        IF (INFO%IALGO==3) THEN
           IFLAG=0
        ENDIF
        IF (N < ABS(INFO%NELMDL)) IFLAG=13

        IF (INFO%IALGO/=2) THEN
           CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
                LMDIM,CDIJ,CQIJ, IFLAG,SV,T_INFO,P,IO%IU0,E%EXHF,EXHF_ACFDT=E%EXHF_ACFDT)
        ENDIF

        CALL STOP_TIMING("G",IO%IU6,"EDDIAG",XMLTAG="diag")
        DWRITE0 "eddiag is ok"
      ELSEIF (MOD(INFO%ICHARG,10)==5) THEN
         ! for ICHARG=5, the states have been rotated according to the supplied GAMMA
         ! diagonalize back to our current one-electron Hamiltonian
         CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
              LMDIM,CDIJ,CQIJ, IFLAG,SV,T_INFO,P,IO%IU0,E%EXHF,EXHF_ACFDT=E%EXHF_ACFDT)
         CALL STOP_TIMING("G",IO%IU6,"EDDIAG",XMLTAG="diag")
      ENDIF
!-----------------------------------------------------------------------
      select_algo: IF (INFO%LEXACT_DIAG) THEN
        CALL EDDIAG_EXACT(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
                LMDIM,CDIJ,CQIJ, IFLAG,SV,T_INFO,P,IO%IU0,IO%IU6,E%EXHF,E%EXHF_ACFDT)

        CALL STOP_TIMING("G",IO%IU6,"EDDIAG",XMLTAG="diag")
        DWRITE0 "EDDIAG_EXACT is ok"

      ELSE IF (INFO%LRMM .AND. .NOT. LDELAY) THEN
!
! RMM-DIIS algorithm
!
        CALL EDDRMM(HAMILTONIAN,GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,WDES, &
             LMDIM,CDIJ,CQIJ, RMS,DESUM1,ICOUEV, SV, IO%IU6,IO%IU0, &
             N < ABS(INFO%NELMDL)-ABS(INFO%NELMDL)/4)
        ! previous line selects  special algorithm during delay

        CALL STOP_TIMING("G",IO%IU6,"RMM-DIIS",XMLTAG="diis")
        DWRITE0 "eddrmm is ok"

      ELSE IF (INFO%LDAVID) THEN
!
! blocked Davidson algorithm,
!
        NSIM=WDES%NSIM*2
#ifdef MPI
        NSIM=((WDES%NSIM*2+WDES%COMM_INTER%NCPU-1)/WDES%COMM_INTER%NCPU)*WDES%COMM_INTER%NCPU
#endif
        CALL EDDAV(HAMILTONIAN, P, GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,WDES, NSIM, &
             LMDIM,CDIJ,CQIJ, RMS,DESUM1,ICOUEV, SV, E%EXHF, IO%IU6,IO%IU0, &
             LDELAY=.FALSE.,  LSUBROTI=INFO%LDIAG, LEMPTY=.FALSE., LHF= N>=ABS(INFO%NELMDL), &
             EXHF_ACFDT= E%EXHF_ACFDT)
        CALL STOP_TIMING("G",IO%IU6,"EDDAV",XMLTAG="dav")
        DWRITE0 "edddav is ok"

      ELSE IF (INFO%IHARMONIC>0) THEN
         NSIM=WDES%NSIM
#ifdef MPI
         NSIM=((WDES%NSIM+WDES%COMM_INTER%NCPU-1)/WDES%COMM_INTER%NCPU)*WDES%COMM_INTER%NCPU
#endif
         IF (INFO%EREF==0) THEN
            CALL EDDAV_INNER(HAMILTONIAN, P, GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,WDES, NSIM, &
                 LMDIM,CDIJ,CQIJ, RMS,DESUM1,ICOUEV, SV, E%EXHF, IO%IU6,IO%IU0, &
                 LEMPTY=.FALSE. )
         ELSE
            CALL EDDAV_INNER(HAMILTONIAN, P, GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,WDES, NSIM, &
                 LMDIM,CDIJ,CQIJ, RMS,DESUM1,ICOUEV, SV, E%EXHF, IO%IU6,IO%IU0, &
                 LEMPTY=.FALSE., EREF=INFO%EREF)
         ENDIF
         CALL STOP_TIMING("G",IO%IU6,"EINNER",XMLTAG="dav")
         ! since we iterate deep it is safer to recalculate projeciton operators
         ! and reorthogonalize
         CALL PROALL (GRID,LATT_CUR,NONLR_S,NONL_S,W)

      ELSE IF (INFO%IALGO==5 .OR.INFO%IALGO==6 .OR. &
     &         INFO%IALGO==7 .OR.INFO%IALGO==8 .OR. INFO%IALGO==0) THEN select_algo

!
! CG (Teter, Alan, Payne) potential is fixed !!
!
        CALL EDSTEP(GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,WDES, &
             LMDIM,CDIJ,CQIJ, RMS,DESUM1,ICOUEV, SV,IO%IU6, IO%IU0)

        CALL STOP_TIMING("G",IO%IU6,"EDSTEP",XMLTAG="cg")
        DWRITE0 "edstep is ok"

      ENDIF select_algo
!-----------------------------------------------------------------------
! orthogonalise all bands (necessary only for residuum-minimizer
! or inner eigenvalue problems, since they iterate very deep)
!
      IF (.NOT.INFO%LORTHO .AND. .NOT. LDELAY) THEN

        CALL ORTHCH(WDES,W, INFO%LOVERL, LMDIM,CQIJ)

        CALL STOP_TIMING("G",IO%IU6,"ORTHCH",XMLTAG="orth")
        DWRITE0 "ortch is ok"
      ENDIF
!
! sub space rotation after eigen value optimization
!
      IF (INFO%LCDIAG .AND. .NOT. LDELAY) THEN

        IF (INFO%LDIAG) THEN
           IFLAG=3
        ELSE
           IFLAG=4
        ENDIF

        CALL REDIS_PW_OVER_BANDS(WDES, W)
        CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
             LMDIM,CDIJ,CQIJ, IFLAG,SV,T_INFO,P,IO%IU0,E%EXHF,EXHF_ACFDT=E%EXHF_ACFDT)

        CALL STOP_TIMING("G",IO%IU6,"EDDIAG",XMLTAG="diag")
      ENDIF
!=======================================================================
! recalculate the broadened density of states and fermi-weights
! recalculate depletion charge size
!=======================================================================
      CALL MRG_CEL(WDES,W)
      IF (INFO%IALGO/=3) THEN
         E%EENTROPY=0
         DOS=0
         DOSI=0
         CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
            INFO%NUP_DOWN, E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE.,  &
            NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
      ENDIF
      DWRITE0 "densta is ok"
!=======================================================================
! calculate free-energy and bandstructur-energy
! EBANDSTR = sum of the energy eigenvalues of the electronic states
!         weighted by the relative weight of the special k point
! TOTEN = total free energy of the system
!=======================================================================
      E%EBANDSTR=BANDSTRUCTURE_ENERGY(WDES, W)
      TOTEN=E%EBANDSTR+E%DENC+E%XCENC+E%TEWEN+E%PSCENC+E%EENTROPY+E%PAWPS+E%PAWAE+INFO%EALLAT+E%EXHF+ECORE()+ Ediel_sol
!-MM- Added to accomodate constrained moment calculations
      IF (M_CONSTRAINED()) TOTEN=TOTEN+E_CONSTRAINT()
      io_begin
      CALL WRITE_CONSTRAINED_M(17,.FALSE.)
      io_end
!-MM- end of additions
!---- write total energy to OSZICAR file and stdout
      DESUM(N)=TOTEN-TOTENL
      ECONV=DESUM(N)
      io_begin
  305 FORMAT('CG : ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5, &
     &       I6,'  ',E10.3)
  302 FORMAT('NONE ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5, &
     &       I6,'  ',E10.3)
  303 FORMAT('EIG: ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5, &
     &       I6,'  ',E10.3)
  304 FORMAT('DIA: ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5, &
     &       I6,'  ',E10.3)
 1303 FORMAT('RMM: ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5, &
     &       I6,'  ',E10.3)
10303 FORMAT('DAV: ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5, &
     &       I6,'  ',E10.3)
20303 FORMAT('JDH: ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5, &
     &       I6,'  ',E10.3)
30303 FORMAT('DAVI:',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5, &
     &       I6,'  ',E10.3)

      IF (INFO%LEXACT_DIAG) THEN
        WRITE(17,  304,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
        IF (IO%IU0>=0) &
        WRITE(IO%IU0, 304,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
      ELSE IF (INFO%LRMM .AND. .NOT. LDELAY) THEN
        WRITE(17,  1303,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
        IF (IO%IU0>=0) &
        WRITE(IO%IU0, 1303,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
      ELSE IF (INFO%IHARMONIC==1) THEN
        WRITE(17, 20303,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
        IF (IO%IU0>=0) &
        WRITE(IO%IU0,20303,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
      ELSE IF (INFO%IHARMONIC==2) THEN
        WRITE(17, 30303,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
        IF (IO%IU0>=0) &
        WRITE(IO%IU0,30303,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
      ELSE IF (INFO%LDAVID) THEN
        WRITE(17, 10303,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
        IF (IO%IU0>=0) &
        WRITE(IO%IU0,10303,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
      ELSE IF (INFO%IALGO==4) THEN
        WRITE(17,   304,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
        IF (IO%IU0>=0) &
        WRITE(IO%IU0,  304,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
      ELSE IF (INFO%IALGO==3) THEN
        WRITE(17,   303,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
        IF (IO%IU0>=0) &
        WRITE(IO%IU0,  303,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
      ELSE IF (INFO%IALGO==2) THEN
        WRITE(17,   302,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
        IF (IO%IU0>=0) &
        WRITE(IO%IU0,  302,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
      ELSE
        WRITE(17,   305,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
        IF (IO%IU0>=0) &
        WRITE(IO%IU0,  305,ADVANCE='NO')  N,TOTEN,DESUM(N),DESUM1,ICOUEV,RMS
      ENDIF

      CALL STOP_TIMING("G",IO%IU6,"DOS")
      io_end
!=======================================================================
!  Test for Break condition
!=======================================================================


      INFO%LABORT=.FALSE.
      LABORT_WITHOUT_CONV=.FALSE.

!-----eigenvalues and energy must be converged
      IF(ABS(DESUM(N))<INFO%EDIFF.AND.ABS(DESUM1)<INFO%EDIFF) INFO%LABORT=.TRUE.
!-----charge-density not constant and in last cycle no change of charge
      IF (.NOT. INFO%LMIX .AND. .NOT. INFO%LCHCON .AND. MIX%IMIX/=0) INFO%LABORT=.FALSE.
!-----do not stop during the non-selfconsistent startup phase
      IF (N <= ABS(INFO%NELMDL)) INFO%LABORT=.FALSE.
!-----do not stop before minimum number of iterations is reached
      IF (N < ABS(INFO%NELMIN)) INFO%LABORT=.FALSE.
!-----but stop after INFO%NELM steps no matter where we are now
      IF (N>=INFO%NELM) THEN
         IF (.NOT.INFO%LABORT) LABORT_WITHOUT_CONV=.TRUE.
         INFO%LABORT=.TRUE.
      ENDIF

      IF ((IO%LORBIT>=10).AND.(MOD(N,5)==0).AND.WDES%LNONCOLLINEAR) THEN
         CALL SPHPRO_FAST( &
          GRID,LATT_CUR, P,T_INFO,W, WDES, 71,IO%IU6,&
          INFO%LOVERL,LMDIM,CQIJ, LDIMP, LDIMP,LMDIMP,.FALSE., IO%LORBIT, PAR, &
          EFERMI, KPOINTS%EMIN, KPOINTS%EMAX)
      ENDIF
! ======================================================================
! If the end of the electronic loop is reached
! calculate accurate initial state core level shifts
! if required
! ======================================================================
      IF (INFO%LABORT .AND. ACCURATE_CORE_LEVEL_SHIFTS()) THEN

         ALLOCATE(CDIJ_TMP(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ))
         CDIJ_TMP=CDIJ

         CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
              WDES%NCDIJ, LMDIM, CDIJ_TMP(1,1,1,1), RHOLM, CRHODE, &
              E, LMETA= .FALSE. , LASPH=INFO%LASPH , LCOREL=.TRUE. )
         DEALLOCATE(CDIJ_TMP)
      ENDIF
!========================= subroutine CHSP  ============================
! if charge density is updated
!  ) first copy current charge to CHTOTL
!  ) set  INFO%LPOTOK to .F. this requires a recalculation of the local pot.
!  ) set INFO%LMIX to .T.
!  ) call subroutine SET_CHARGE to generate the new charge density
!  ) then performe mixing
! MIND:
! ) if delay is selected  do not update
! ) if convergence corrections to forces are calculated do not update charge
!   in last iteration
!=======================================================================
      IF (MOD(INFO%ICHARG,10)==5) THEN
         IF (IO%LORBIT==14) CALL SPHPRO_FAST( &
              GRID,LATT_CUR, P,T_INFO,W, WDES, 71,IO%IU6,&
              INFO%LOVERL,LMDIM,CQIJ, LDIMP, LDIMP,LMDIMP,.FALSE., IO%LORBIT,PAR, &
              EFERMI, KPOINTS%EMIN, KPOINTS%EMAX)
         CALL LPRJ_PROALL(W,WDES,GRID,P,CQIJ,LATT_CUR,T_INFO,INFO,IO%IU6,IO%IU0)
         CALL LPRJ_WRITE(IO%IU6, IO%IU0, W)
         CALL LPRJ_LDApU(IO, W)   ! write a LDA+U GAMMA file
         CALL LPRJ_DEALLOC_COVL
      ENDIF

      INFO%LMIX=.FALSE.
      MIX%NEIG=0

      IF (.NOT. INFO%LCHCON .AND. .NOT. (INFO%LABORT .AND. INFO%LCORR ) &
     &    .AND. N >= ABS(INFO%NELMDL)  ) THEN

      DO ISP=1,WDES%NCDIJ
         CALL RC_ADD(CHTOT(1,ISP),1.0_q,CHTOT(1,ISP),0.0_q,CHTOTL(1,ISP),GRIDC)
      ENDDO

      IF (LDO_METAGGA().AND.LMIX_TAU()) THEN
         DO ISP=1,WDES%NCDIJ
            CALL RC_ADD(KINEDEN%TAU(1,ISP),1.0_q,KINEDEN%TAU(1,ISP),0.0_q,KINEDEN%TAUL(1,ISP),GRIDC)
         ENDDO
      ENDIF

      RHOLM_LAST=RHOLM

      INFO%LPOTOK=.FALSE.

      IF (MOD(INFO%ICHARG,10)==5) THEN
         CALL REMOVE_VASP_LOCK(WDES%COMM) ! remove vasp.lock file
         CALL WAIT_VASP_LOCK(WDES%COMM)   ! wait that it is re-created by other process
         CALL ADD_GAMMA_FROM_FILE( WDES, W, KPOINTS, INFO%NELECT, INFO%NUP_DOWN, INFO%LABORT, IO  )
      ENDIF

      !LWANNIER90_RUN=.TRUE.
      !WRITE_MMN_AMN=.FALSE.
      !if (LDMFT.eq..TRUE..AND. MOD(N,4).eq.0) then
      if (LDMFT.eq..TRUE.) then
        CALL MLWF_WANNIER90_AMN(WDES,W,P,CQIJ,T_INFO,LATT_CUR,INFO,IO)
        !CALL MLWF_WANNIER90(WDES,W,P,CQIJ,T_INFO,LATT_CUR,INFO,IO)
!      else
      !CALL MLWF_WANNIER90_EIG(WDES,W,P,CQIJ,T_INFO,LATT_CUR,INFO,IO)
!      endif
        CALL MPI_barrier( WDES%COMM%MPI_COMM, ierr )
        !STOP
        !CALL Compute_DMFT_from_amn(KPOINTS_FULL%NKPTS,num_wann_,KPOINTS_FULL%VKPT,KPOINTS_FULL%WTKPT,band_win_,DMFT_eval_,DMFT_evec_)
        !write(*,*) "Computing DMFT from AMN..."
        !write(*,'(A,I4)') "num_wann_:", num_wann_
        CALL Compute_DMFT_from_amn(WDES%NKPTS,num_wann_,WDES%VKPT,WDES%WTKPT,band_win_,DMFT_eval_,DMFT_evec_)
        !CALL Compute_DMFT(WDES%NKPTS,num_wann_,WDES%VKPT,WDES%WTKPT,band_win_,DMFT_eval_,DMFT_evec_)
        !IF (IO%IU6>=0) then
        !  DO IK=1,WDES%NKPTS
        !    write(*,*) IK,DMFT_eval_(:,IK),DMFT_evec_(:,1,IK)
        !  ENDDO
        !ENDIF
        !STOP
      endif

!!!!!  !!!!!!!!!!! Rotate WAVE function !!!!!!!!!!!!!!!!!!

!        ALLOCATE(FERWE(WDES%NBANDS,WDES%NKPTS))
!        ALLOCATE(WaveW(WDES%NRPLWV,WDES%NB_TOT,WDES%NKPTS))
!        ALLOCATE(WaveW_tot(WDES%NRPLWV,WDES%NB_TOT,WDES%NKPTS))
!        ALLOCATE(WProj(WDES%NPRO_TOT,WDES%NB_TOT,WDES%NKPTS))
!        ALLOCATE(WProj_tot(WDES%NPRO_TOT,WDES%NB_TOT,WDES%NKPTS))
!        FERWE=0;WaveW=0;WaveW_tot=0;WProj=0;WProj_tot=0
!        DO NK=1,WDES%NKPTS
!#ifdef MPI
!        IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
!#endif
!        DO NB=1,WDES%NBANDS
!          FERWE(NB,NK)=W%FERWE(NB,NK,1)
!          DO NP=1,WDES%NRPLWV
!             WaveW(NP,WDES%NB_LOW+(NB-1)*WDES%NB_PAR,NK)=W%CW(NP,NB,NK,1)!/WDES%COMM_KINTER%NCPU
!          ENDDO
!          DO NP=1,WDES%NPRO_TOT
!             WProj(NP,WDES%NB_LOW+(NB-1)*WDES%NB_PAR,NK)=W%CPROJ(NP,NB,NK,1)!/WDES%COMM_KINTER%NCPU
!          ENDDO
!        ENDDO
!        ENDDO
!
!        CALL MPI_ALLREDUCE(WaveW,WaveW_tot,WDES%NRPLWV*WDES%NB_TOT*WDES%NKPTS,MPI_DOUBLE_COMPLEX,MPI_SUM,WDES%COMM%MPI_COMM,ierr)
!        CALL MPI_ALLREDUCE(WProj,WProj_tot,WDES%NPRO_TOT*WDES%NB_TOT*WDES%NKPTS,MPI_DOUBLE_COMPLEX,MPI_SUM,WDES%COMM%MPI_COMM,ierr)
!        DO NK=1,WDES%NKPTS
!#ifdef MPI
!          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
!#endif
!          nbmin=band_win_(1,NK); nbmax=band_win_(2,NK)
!          num_band_max=nbmax-nbmin+1
!          DO NB=1,WDES%NBANDS
!            NBG=WDES%NB_LOW+(NB-1)*WDES%NB_PAR
!            if (NBG .GE. nbmin .AND. NBG .LE. nbmax) then
!              W%FERWE(NB,NK,1)=0
!              if (NBG-nbmin+1 .LE. num_wann_) W%FERWE(NB,NK,1)=DMFT_EVAL_(NBG-nbmin+1,NK)
!              W%CW(:,NB,NK,1)=0
!              W%CPROJ(:,NB,NK,1)=0
!              DO NP=1,num_band_max
!                W%CW(:,NB,NK,1)=W%CW(:,NB,NK,1)+WaveW_tot(:,nbmin+NP-1,NK)*DMFT_EVEC_(NP,NBG-nbmin+1,NK)
!                W%CPROJ(:,NB,NK,1)=W%CPROJ(:,NB,NK,1)+WProj_tot(:,nbmin+NP-1,NK)*DMFT_EVEC_(NP,NBG-nbmin+1,NK)
!              ENDDO
!            endif
!          ENDDO
!        ENDDO
!        DEALLOCATE(WaveW_tot)
!        DEALLOCATE(WProj_tot)
        !tot_n=0
        !DO NK=1,WDES%NKPTS
        !  DO NB=1,WDES%NBANDS
        !    tot_n=tot_n+W%FERWE(NB,NK,1)*WDES%WTKPT(NK)*2
        !  ENDDO
        !ENDDO
        !write(*,*) tot_n
        !DO NK=1,WDES%NKPTS
        !  DO NB=1,num_wann_
        !    tot_n=0
        !    DO NP=1,num_wann_
        !      tot_n=tot_n+DMFT_EVEC_(NB,NP,NK)*CONJG(DMFT_EVEC_(NB,NP,NK))
        !    ENDDO
        !    write(*,*) NK,NB,tot_n
        !  ENDDO
        !ENDDO
!      endif

      CALL SET_CHARGE(W, WDES, INFO%LOVERL, &
           GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
           LATT_CUR, P, SYMM, T_INFO, &
           CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

      CALL SET_KINEDEN(GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM, &
           T_INFO%NIONS,W,WDES,KINEDEN)


!      if (LDMFT.eq..TRUE.) then
!      !if (LDMFT.eq..TRUE..AND. MOD(N,4).eq.0) then
!        DO NK=1,WDES%NKPTS
!#ifdef MPI
!          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
!#endif
!          DO NB=1,WDES%NBANDS
!             W%FERWE(NB,NK,1)=FERWE(NB,NK)
!             NBG=WDES%NB_LOW+(NB-1)*WDES%NB_PAR
!             if (NBG .GE. nbmin .AND. NBG .LE. nbmax) then
!               W%CW(:,NB,NK,1)=WaveW(:,NBG,NK)
!               W%CPROJ(:,NB,NK,1)=WProj(:,NBG,NK)
!             endif
!          ENDDO
!        ENDDO
!        DEALLOCATE(FERWE)
!        DEALLOCATE(WaveW)
!        DEALLOCATE(WProj)
!      endif

      CALL STOP_TIMING("G",IO%IU6,"CHARGE")
!-----------------------------------------------------------------------

      IF (MIX%IMIX/=0) THEN
      INFO%LMIX=.TRUE.
! update one center charge density
      IF (LCORREL()) THEN
         CALL SET_RHO_ONE_CENTRE(CRHODE,RHO_ONE_CENTRE)
      ENDIF

      IF (MIX%IMIX==4) THEN
!  broyden mixing ... :
        IF (LCORREL()) THEN
           CALL BRMIX(KINEDEN,GRIDB,GRIDC,IO,MIX,B_TO_C, &
          &   (2*GRIDC%MPLWV),CHTOT,CHTOTL,WDES%NCDIJ,LATT_CUR%B, &
          &   LATT_CUR%OMEGA, N_MIX_PAW, RHOLM, RHOLM_LAST, &
          &   RMST,RMSC,RMSP,WEIGHT,.TRUE.,IERRBR, &
          &   NMIX_ONE_CENTRE,RHO_ONE_CENTRE,RHO_ONE_CENTRE_LAST)
           MIX%LRESET=.FALSE.
        ELSE
           CALL BRMIX(KINEDEN,GRIDB,GRIDC,IO,MIX,B_TO_C, &
          &   (2*GRIDC%MPLWV),CHTOT,CHTOTL,WDES%NCDIJ,LATT_CUR%B, &
          &   LATT_CUR%OMEGA, N_MIX_PAW, RHOLM, RHOLM_LAST, &
          &   RMST,RMSC,RMSP,WEIGHT,.TRUE.,IERRBR)
           MIX%LRESET=.FALSE.
        ENDIF
      ELSE
!  simple mixing ... :
        RMST=0
        CALL MIX_SIMPLE(GRIDC,MIX,WDES%NCDIJ, CHTOT,CHTOTL, &
             N_MIX_PAW, RHOLM, RHOLM_LAST, LATT_CUR%B, LATT_CUR%OMEGA, RMST)
! test
        IF (LCORREL()) THEN
           RHO_ONE_CENTRE(:)=RHO_ONE_CENTRE_LAST(:)+ &
           &   0.2*(RHO_ONE_CENTRE(:)-RHO_ONE_CENTRE_LAST(:))
           RHO_ONE_CENTRE_LAST(:)=RHO_ONE_CENTRE(:)
        ENDIF
! test
      ENDIF
      CALL STOP_TIMING("G",IO%IU6,"MIXING")
      DWRITE0 "mixing is ok"

!---- ENDIF (MIX%IMIX/=0)     end of mixing
      ENDIF
!-----ENDIF (.NOT.INFO%LCHCON)   end of charge update
      ENDIF

      IF (W%OVER_BAND) THEN
         CALL REDIS_PW_OVER_BANDS(WDES, W)
         CALL STOP_TIMING("G",IO%IU6,"REDIS")
      ENDIF
!=======================================================================
! total time used for this step
!=======================================================================
     CALL SEPERATOR_TIMING(IO%IU6)
     CALL STOP_TIMING("LOOP",IO%IU6,XMLTAG='total')

!=======================================================================
!  important write statements
!=======================================================================

 2440 FORMAT(/' eigenvalue-minimisations  :',I6,/ &
     &       ' total energy-change (2. order) :',E14.7,'  (',E14.7,')')
 2441 FORMAT(/ &
     &       ' Broyden mixing:'/ &
     &       '  rms(total) =',E12.5,'    rms(broyden)=',E12.5,/ &
     &       '  rms(prec ) =',E12.5/ &
     &       '  weight for this iteration ',F10.2)

 2442 FORMAT(/' eigenvalues of (default mixing * dielectric matrix)' / &
             '  average eigenvalue GAMMA= ',F8.4,/ (10F8.4))

 200  FORMAT(' number of electron ',F15.7,' magnetization ',3F15.7)
 201  FORMAT(' augmentation part  ',F15.7,' magnetization ',3F15.7)


      DO I=1,WDES%NCDIJ
         RHOTOT(I)=RHO0(GRIDC, CHTOT(1,I))
         RHOAUG(I)=RHOTOT(I)-RHO0(GRID_SOFT, CHDEN(1,I))
      END DO

      io_begin

    ! iteration counts
      WRITE(IO%IU6,2440) ICOUEV,DESUM(N),DESUM1

    ! charge density
      WRITE(IO%IU6,200) RHOTOT
      IF (INFO%LOVERL) THEN
        WRITE(IO%IU6,201) RHOAUG
      ENDIF
    ! dipol moment
      IF (DIP%LCOR_DIP) CALL WRITE_DIP(IO%IU6)

    ! mixing
      IF ( INFO%LMIX .AND. MIX%IMIX==4 ) THEN
        IF (IERRBR/=0) THEN
          IF (IO%IU0>=0) &
          WRITE(IO%IU0,*) 'ERROR: Broyden mixing failed, tried ''simple '// &
                      'mixing'' now and reset mixing at next step!'
          IF (IO%IU6>=0) &
          WRITE(IO%IU6,*) 'ERROR: Broyden mixing failed, tried ''simple '// &
                       'mixing'' now and reset mixing at next step!'
        ENDIF

        IF (IO%NWRITE>=2 .OR. NSTEP==1) THEN
          WRITE(IO%IU6,2441) RMST,RMSC,RMSP,WEIGHT
          IF (ABS(RMST-RMSC)/RMST> 0.1_q) THEN
            WRITE(IO%IU6,*) ' WARNING: grid for Broyden might be to small'
          ENDIF
        ENDIF
        IF (IO%IU0>=0) WRITE(IO%IU0,308) RMST
        WRITE(17,308) RMST
   308  FORMAT('   ',E10.3)
        IF (MIX%NEIG > 0) THEN
           WRITE(IO%IU6,2442) MIX%AMEAN,MIX%EIGENVAL(1:MIX%NEIG)
        ENDIF
      ELSE IF (INFO%LMIX) THEN
        IF (IO%IU0>=0) &
        WRITE(IO%IU0,308) RMST
        WRITE(17,308) RMST
      ELSE
        IF (IO%IU0>=0) &
        WRITE(IO%IU0,*)
        WRITE(17 ,*)
      ENDIF
 io1: IF (IO%NWRITE>=2 .OR. (NSTEP==1)) THEN
    ! energy
      IF (LCORREL()) THEN
         WRITE(IO%IU6,7241) E%PSCENC,E%TEWEN,E%DENC,E%EXHF,E%XCENC,E%PAWPS,E%PAWAE, &
                         E%EENTROPY,E%EBANDSTR,INFO%EALLAT+ECORE(),Ediel_sol,TOTEN, &
                         TOTEN-E%EENTROPY,TOTEN-E%EENTROPY/(2+NORDER)
      ELSE
         WRITE(IO%IU6,7240) E%PSCENC,E%TEWEN,E%DENC,E%EXHF,E%XCENC,E%PAWPS,E%PAWAE, &
                         E%EENTROPY,E%EBANDSTR,INFO%EALLAT,Ediel_sol,TOTEN, &
                         TOTEN-E%EENTROPY,TOTEN-E%EENTROPY/(2+NORDER)
      ENDIF

      IF (LHFCALC) THEN
         WRITE(IO%IU6,'( "  exchange ACFDT corr.  = ",F18.8,"  see jH, gK, PRB 81, 115126")') E%EXHF_ACFDT
      ENDIF

7240  FORMAT(/ &
              ' Free energy of the ion-electron system (eV)'/ &
     &        '  ---------------------------------------------------'/ &
     &        '  alpha Z        PSCENC = ',F18.8/ &
     &        '  Ewald energy   TEWEN  = ',F18.8/ &
     &        '  -Hartree energ DENC   = ',F18.8/ &
     &        '  -exchange      EXHF   = ',F18.8/ &
     &        '  -V(xc)+E(xc)   XCENC  = ',F18.8/ &
     &        '  PAW double counting   = ',2F18.8/ &
     &        '  entropy T*S    EENTRO = ',F18.8/ &
     &        '  eigenvalues    EBANDS = ',F18.8/ &
     &        '  atomic energy  EATOM  = ',F18.8/ &
     &        '  Solvation  Ediel_sol  = ',F18.8/ &
     &        '  ---------------------------------------------------'/ &
     &        '  free energy    TOTEN  = ',F18.8,' eV'// &
     &        '  energy without entropy =',F18.8, &
     &        '  energy(sigma->0) =',F18.8)
7241  FORMAT(/ &
              ' Free energy of the ion-electron system (eV)'/ &
     &        '  ---------------------------------------------------'/ &
     &        '  alpha Z        PSCENC = ',F18.8/ &
     &        '  Ewald energy   TEWEN  = ',F18.8/ &
     &        '  -Hartree energ DENC   = ',F18.8/ &
     &        '  -exchange      EXHF   = ',F18.8/ &
     &        '  -V(xc)+E(xc)   XCENC  = ',F18.8/ &
     &        '  PAW double counting   = ',2F18.8/ &
     &        '  entropy T*S    EENTRO = ',F18.8/ &
     &        '  eigenvalues    EBANDS = ',F18.8/ &
     &        '  core contrib.  ECORE  = ',F18.8/ &
     &        '  Solvation  Ediel_sol  = ',F18.8/ &
     &        '  ---------------------------------------------------'/ &
     &        '  free energy    TOTEN  = ',F18.8,' eV'// &
     &        '  energy without entropy =',F18.8, &
     &        '  energy(sigma->0) =',F18.8)
72612 FORMAT(//&
     &        '  METAGGA EXCHANGE AND CORRELATION (eV)'/ &
     &        '  ---------------------------------------------------'/ &
     &        '  LDA+GGA E(xc)  EXCG   = ',F18.6/ &
     &        '  LDA+GGA PAW    PS : AE= ',2F18.6/ &
     &        '  core xc             AE= ',1F18.6/ &
     &        '  metaGGA E(xc)  EXCM   = ',F18.6/ &
     &        '  metaGGA PAW    PS : AE= ',2F18.6/ &
     &        '  metaGGA core xc     AE= ',1F18.6/ &
     &        '  ---------------------------------------------------'/ &
     &        '  METAGGA result:'/ &
     &        '  free  energy   TOTEN  = ',F18.6,' eV'// &
     &        '  energy  without entropy=',F18.6, &
     &        '  energy(sigma->0) =',F16.6)
      ELSE io1
      WRITE(IO%IU6,7242) TOTEN,TOTEN-E%EENTROPY
 7242 FORMAT(/'  free energy = ',E20.12, &
     &        '  energy without entropy= ',E20.12)

      ENDIF io1
!     too slow on many servers nowadays
!     IF (IO%LOPEN) CALL WFORCE(IO%IU6)
!     IF (IO%LOPEN) CALL WFORCE(17)
      WRITE(IO%IU6,130)
      io_end
!=======================================================================
!  perform some additional write statments if required
!=======================================================================
!-----Eigenvalues and weights
      IF (((NSTEP==1 .OR. NSTEP==DYN%NSW).AND.INFO%LABORT).OR. &
     &     (IO%NWRITE>=1 .AND.INFO%LABORT).OR.IO%NWRITE>=3) THEN

        ! calculate and print the core level shifts
      IF (INFO%LOVERL) THEN
        CALL CL_SHIFT_PW( GRIDC, LATT_CUR, IRDMAX,  &
           T_INFO, P, WDES%NCDIJ, CVTOT, MAX(INFO%ENAUG,INFO%ENMAX), IO%IU6)
      ELSE
        IF (IO%IU0>=0) WRITE(IO%IU0,*) "WARNING: NC-PP core level shifts not calculated"
      ENDIF
      ENDIF

      IF (((NSTEP==1 .OR. NSTEP==DYN%NSW).AND.INFO%LABORT).OR. &
     &     (IO%NWRITE>=1 .AND.INFO%LABORT).OR.IO%NWRITE>=3) THEN
        CALL KPAR_SYNC_CELTOT( WDES, W )

      io_begin
        CALL RHOAT0(P,T_INFO, BETATO,LATT_CUR%OMEGA)

        WRITE(IO%IU6,2202) EFERMI,REAL( E%CVZERO ,KIND=q) ,E%PSCENC/INFO%NELECT+BETATO
 2202   FORMAT(' E-fermi : ', F8.4,'     XC(G=0): ',F8.4, &
     &         '     alpha+bet :',F8.4/)

        IF(INFO%IHARMONIC==1)THEN
           CALL WRITE_EIGENVAL_RESIDUAL( WDES, W, IO%IU6)
        ELSE
           CALL WRITE_EIGENVAL( WDES, W, IO%IU6)
        END IF
      io_end
      ENDIF

      IF (((NSTEP==1 .OR. NSTEP==DYN%NSW).AND.INFO%LABORT).OR. &
     &     (IO%NWRITE>=1 .AND.INFO%LABORT).OR.IO%NWRITE>=3) THEN
      io_begin
!-----Charge-density along one line
        WRITE(IO%IU6,130)
        DO I=1,WDES%NCDIJ
           WRITE(IO%IU6,*)'soft charge-density along one line, spin component',I
           WRITE(IO%IU6,'(10(6X,I4))') (II,II=0,9)
           CALL WRT_RC_LINE(IO%IU6,GRID_SOFT, CHDEN(1,I))
           IF (INFO%LOVERL) THEN
              WRITE(IO%IU6,*)'total charge-density along one line'
              CALL WRT_RC_LINE(IO%IU6,GRIDC, CHTOT(1,I))
           ENDIF
           WRITE(IO%IU6,*)
        ENDDO
!-----pseudopotential strength and augmentation charge
        DO NI=1,1
        DO I=1,WDES%NCDIJ
           WRITE(IO%IU6,*) 'pseudopotential strength for first ion, spin component:',I
           DO LP=1,P(1)%LMMAX
              WRITE(IO%IU6,'(16(F7.3,1X))') &
     &             (CDIJ(L,LP,NI,I),L=1,MIN(8,P(1)%LMMAX))
!     &             (REAL(CDIJ(L,LP,NI,I),q),AIMAG(CDIJ(L,LP,1,I))*1000,L=1,MIN(16,P(1)%LMMAX))
           ENDDO
        ENDDO
        ENDDO

        IF (INFO%LOVERL) THEN
        DO NI=1,1
        DO I=1,WDES%NCDIJ
           WRITE(IO%IU6,*) 'total augmentation occupancy for first ion, spin component:',I
           DO LP=1,P(1)%LMMAX
              WRITE(IO%IU6,'(16(F7.3,1X))') &
     &             (REAL(CRHODE(L,LP,NI,I),q),L=1,MIN(16,P(1)%LMMAX))
           ENDDO
!           DO LP=1,P(1)%LMMAX
!              WRITE(IO%IU6,'(16(F7.3,1X))') &
!     &             (AIMAG(CRHODE(L,LP,1,I))*1E6,L=1,MIN(16,P(1)%LMMAX))
!           ENDDO
        ENDDO
        ENDDO
        ENDIF
      io_end

      ENDIF
!=======================================================================
!  xml related output
!=======================================================================
      CALL XML_TAG("energy")
      IF (INFO%LABORT .OR. N==1) THEN
         CALL XML_TAG_REAL("alphaZ",E%PSCENC)
         CALL XML_TAG_REAL("ewald", E%TEWEN)
         CALL XML_TAG_REAL("hartreedc",E%DENC)
         CALL XML_TAG_REAL("XCdc",E%XCENC)
         CALL XML_TAG_REAL("pawpsdc",E%PAWPS)
         CALL XML_TAG_REAL("pawaedc",E%PAWAE)
         CALL XML_TAG_REAL("eentropy",E%EENTROPY)
         CALL XML_TAG_REAL("bandstr",E%EBANDSTR)
         CALL XML_TAG_REAL("atom",INFO%EALLAT)
         CALL XML_ENERGY(TOTEN, TOTEN-E%EENTROPY, TOTEN-E%EENTROPY/(2+NORDER))
      ELSE
         CALL XML_ENERGY(TOTEN, TOTEN-E%EENTROPY, TOTEN-E%EENTROPY/(2+NORDER))
      ENDIF
      CALL XML_CLOSE_TAG
      CALL XML_CLOSE_TAG("scstep")
!=======================================================================
! relaxed core related output
!=======================================================================
      IF (INFO%LABORT.AND.LCORREL()) THEN
         CALL REPORT(.TRUE.)
! uncomment this line for specific postprocessing of relaxed core stuff
!        CALL RCPOSTPROC
!
      ENDIF
!======================== end of loop ENDLSC ===========================
! This is the end of the selfconsistent calculation loop
!=======================================================================
      IF (INFO%LABORT) THEN
         io_begin
        WRITE(IO%IU6,131)
 131    FORMAT (5X, //, &
     &  '------------------------ aborting loop because EDIFF', &
     &  ' is reached ----------------------------------------'//)
         io_end
        EXIT electron
      ENDIF
      INFO%LSOFT=.FALSE.
#ifndef noSTOPCAR
#ifndef F90_T3D
      CALL RDATAB(IO%LOPEN,'STOPCAR',99,'LABORT','=','#',';','L', &
     &            IDUM,RDUM,CDUM,INFO%LSOFT,CHARAC,NCOUNT,1,IERR)
      ITMP=0 ; IF (INFO%LSOFT) ITMP=1 ; CALLMPI ( M_sum_i(W%WDES%COMM, ITMP, 1))
      IF (ITMP>0) INFO%LSOFT=.TRUE.
#endif
#endif
      IF (INFO%LSOFT) THEN
        io_begin
        IF (IO%IU0>=0) &
        WRITE(IO%IU0,*) 'hard stop encountered!  aborting job ...'
        WRITE(IO%IU6,13131)
13131   FORMAT (5X, //, &
     &  '------------------------ aborting loop because hard', &
     &  ' stop was set ---------------------------------------'//)
        io_end
        EXIT electron
      ENDIF
      TOTENL=TOTEN

      ENDDO electron

!
! calculate dipol corrections now
!
      IF ( DIP%IDIPCO >0 ) THEN
         IF (.NOT. DIP%LCOR_DIP) THEN
            CALL CDIPOL_CHTOT_REC(GRIDC, LATT_CUR,P,T_INFO, &
                 CHTOT,CSTRF,CVTOT, WDES%NCDIJ, INFO%NELECT )

            CALL WRITE_DIP(IO%IU6)
            IF (IO%IU6>0) THEN
               WRITE(IO%IU6,*)
               WRITE(IO%IU6,*) &
                    " *************** adding dipol energy to TOTEN NOW **************** "
            ENDIF
            TOTEN=TOTEN+ DIP%ECORR
         ENDIF
      ENDIF

      ! notify calling routine whether convergence has been reached
      INFO%LABORT=LABORT_WITHOUT_CONV

      DWRITE0 'electron left'

      RETURN
    END SUBROUTINE ELMIN


!**********************************************************************
!
! write eigenvalues
!
!**********************************************************************

    SUBROUTINE WRITE_EIGENVAL(  WDES, W, IU6)
      USE wave
      USE prec
      TYPE (wavedes)     WDES
      TYPE (wavespin)    W          ! wavefunction
      INTEGER IU6
    ! local
      INTEGER ISP, NN, I
      INTEGER NB

      NB=WDES%NB_TOT

      IF (IU6>=0) THEN

        DO ISP=1,WDES%ISPIN
        IF (WDES%ISPIN==2) WRITE(IU6,'(/A,I1)') ' spin component ',ISP
        DO NN=1,WDES%NKPTS
        WRITE(IU6,2201)NN,WDES%VKPT(1,NN),WDES%VKPT(2,NN),WDES%VKPT(3,NN), &
     &      (I,REAL( W%CELTOT(I,NN,ISP) ,KIND=q) ,W%FERTOT(I,NN,ISP)*WDES%RSPIN,I=1,NB)
        ENDDO
        ENDDO

 2201   FORMAT(/' k-point ',I5,' :',3X,3F10.4/ &
     &         '  band No.  band energies     occupation '/ &
     &           (1X,I6,3X,F10.4,3X,F10.5))
     ENDIF

   END SUBROUTINE WRITE_EIGENVAL


   SUBROUTINE WRITE_EIGENVAL_RESIDUAL(  WDES, W, IU6)
     USE wave
     USE prec
     TYPE (wavedes)     WDES
     TYPE (wavespin)    W          ! wavefunction
     INTEGER IU6
     ! local
     INTEGER ISP, NN, I
     INTEGER NB

     NB=WDES%NB_TOT

     IF (IU6>=0) THEN

        DO ISP=1,WDES%ISPIN
           IF (WDES%ISPIN==2) WRITE(IU6,'(/A,I1)') ' spin component ',ISP
           DO NN=1,WDES%NKPTS
              WRITE(IU6,2201)NN,WDES%VKPT(1,NN),WDES%VKPT(2,NN),WDES%VKPT(3,NN), &
                   &      (I,REAL( W%CELTOT(I,NN,ISP) ,KIND=q) ,W%FERTOT(I,NN,ISP)*WDES%RSPIN,W%AUXTOT(I,NN,ISP),I=1,NB)
           ENDDO
        ENDDO

2201    FORMAT(/' k-point ',I5,' :',3X,3F10.4/ &
             &         '  band No.  band energies  occupation  residual'/ &
             &           (1X,I6,3X,F10.4,5X,F10.5,3X,E10.3))
     ENDIF

   END SUBROUTINE WRITE_EIGENVAL_RESIDUAL


   SUBROUTINE WRITE_EIGENVAL_NBANDS(  WDES, W, IU6, NBANDS)
     USE wave
     USE prec
     TYPE (wavedes)     WDES
     TYPE (wavespin)    W          ! wavefunction
     INTEGER IU6
     INTEGER NBANDS
     ! local
     INTEGER ISP, NN, I
     INTEGER NB

     NB=NBANDS

     IF (IU6>=0) THEN

        DO ISP=1,WDES%ISPIN
           IF (WDES%ISPIN==2) WRITE(IU6,'(/A,I1)') ' spin component ',ISP
           DO NN=1,WDES%NKPTS
              WRITE(IU6,2201)NN,WDES%VKPT(1,NN),WDES%VKPT(2,NN),WDES%VKPT(3,NN), &
                   &      (I,REAL( W%CELTOT(I,NN,ISP) ,KIND=q) ,W%FERTOT(I,NN,ISP)*WDES%RSPIN,I=1,NB)
           ENDDO
        ENDDO

2201    FORMAT(/' k-point ',I5,' :',3X,3F10.4/ &
             &         '  band No.  band energies     occupation '/ &
             &           (1X,I6,3X,F10.4,3X,F10.5))
     ENDIF

   END SUBROUTINE WRITE_EIGENVAL_NBANDS
