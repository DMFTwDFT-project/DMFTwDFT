!#define debug
#include "symbol.inc"
      MODULE mlwf
      USE prec
      USE pead

      ! write the wave functions to UNK files
      LOGICAL, PRIVATE, SAVE :: WRITE_UNK
      ! write the mmn and amn files when WANNIER90 runs in lib mode
      LOGICAL, SAVE :: WRITE_MMN_AMN

      ! read amn file instead of computing it
      LOGICAL, PRIVATE, SAVE :: READ_AMN

      ! wannier90_run variables: output
      COMPLEX(q), ALLOCATABLE, SAVE :: U_matrix(:,:,:,:)
      COMPLEX(q), ALLOCATABLE, SAVE :: U_matrix_opt(:,:,:,:)
      LOGICAL, ALLOCATABLE, SAVE    :: lwindow(:,:,:)
      LOGICAL, ALLOCATABLE, SAVE    :: lexclude_band(:)
      REAL(q), ALLOCATABLE, PRIVATE, SAVE :: wann_centres(:,:,:)
      REAL(q), ALLOCATABLE, PRIVATE, SAVE :: wann_spreads(:,:)
      REAL(q), ALLOCATABLE, PRIVATE, SAVE :: spread(:,:)
      integer, parameter :: dp = kind(1.0d0)
      integer, save:: num_band_max, nbmin, nbmax, num_wann_
      integer,allocatable, save :: band_win_(:,:)
      real(kind=dp),allocatable, save :: DMFT_eval_(:,:)
      complex(kind=dp),allocatable, save :: DMFT_evec_(:,:,:)
      logical, save :: LDMFT

      ! we keep this information because it comes in handy sometimes
      ! when we need to use U_matrix(_opt) outside of this module
      REAL(q), ALLOCATABLE :: kpt_latt(:,:)
      INTEGER :: num_kpts

      ! will be set to num_wann if WANNIER90 has run in lib mode
      INTEGER, SAVE :: MLWF_num_wann=-1

      ! prepare the input for wannier90
      LOGICAL, SAVE :: LWANNIER90=.FALSE.
      ! setup and run wannier90 in library mode
      LOGICAL, SAVE :: LWANNIER90_RUN=.FALSE.

      !filename of Wannier projection matrices
      CHARACTER(LEN=7),PUBLIC,PARAMETER :: WANPROJFILE='WANPROJ'
      ! write U_nm(k) matrix
      LOGICAL, PUBLIC, SAVE :: LWRITE_WANPROJ

#ifdef VASP2WANNIER90
      INTERFACE
        subroutine wannier_setup( &
          seed_name_, mp_grid_, num_kpts_, real_lattice_, recip_lattice_, &
          kpt_latt_, num_bands_tot_, num_atoms_, atom_symbols_, atoms_cart_, &
          gamma_only_, spinors_, nntot_, nnlist_, nncell_, num_bands_, &
          num_wann_, proj_site_, proj_l_, proj_m_, proj_radial_, proj_z_, &
          proj_x_, proj_zona_, exclude_bands_ &
#ifdef VASP2WANNIER90v2
          , proj_s_, proj_s_qaxis_ &
#endif
          )
          implicit none
          integer, parameter :: dp = kind(1.0d0) ! statt 'USE w90_constants'
          integer, parameter :: num_nnmax=12     ! statt 'USE w90_parameters'
          character(len=*), intent(in) :: seed_name_
          integer, dimension(3), intent(in) :: mp_grid_
          integer, intent(in) :: num_kpts_
          real(kind=dp), dimension(3,3), intent(in) :: real_lattice_
          real(kind=dp), dimension(3,3), intent(in) :: recip_lattice_
          real(kind=dp), dimension(3,num_kpts_), intent(in) :: kpt_latt_
          integer, intent(in) :: num_bands_tot_
          integer, intent(in) :: num_atoms_
          character(len=*), dimension(num_atoms_), intent(in) :: atom_symbols_
          real(kind=dp), dimension(3,num_atoms_), intent(in) :: atoms_cart_
          logical, intent(in) :: gamma_only_
          logical, intent(in) :: spinors_
          integer, intent(out) :: nntot_
          integer, dimension(num_kpts_,num_nnmax), intent(out) :: nnlist_
          integer,dimension(3,num_kpts_,num_nnmax), intent(out) :: nncell_
          integer, intent(out) :: num_bands_
          integer, intent(out) :: num_wann_
          real(kind=dp), dimension(3,num_bands_tot_), intent(out) :: proj_site_
          integer, dimension(num_bands_tot_), intent(out) :: proj_l_
          integer, dimension(num_bands_tot_), intent(out) :: proj_m_
          integer, dimension(num_bands_tot_), intent(out) :: proj_radial_
          real(kind=dp), dimension(3,num_bands_tot_), intent(out) :: proj_z_
          real(kind=dp), dimension(3,num_bands_tot_), intent(out) :: proj_x_
          real(kind=dp), dimension(num_bands_tot_), intent(out) :: proj_zona_
          integer, dimension(num_bands_tot_), intent(out) :: exclude_bands_
#ifdef VASP2WANNIER90v2
          integer, dimension(num_bands_tot_), optional, intent(out) :: proj_s_
          real(kind=dp), dimension(3,num_bands_tot_), optional, intent(out) :: proj_s_qaxis_
#endif
        end subroutine wannier_setup
        !
        subroutine wannier_run( &
          seed_name_, mp_grid_, num_kpts_, real_lattice_, recip_lattice_, &
          kpt_latt_, num_bands_, num_wann_, nntot_, num_atoms_, atom_symbols_, &
          atoms_cart_, gamma_only_, M_matrix_, A_matrix_, eigenvalues_, &
          U_matrix_, U_matrix_opt_, lwindow_, wann_centres_, wann_spreads_, spread_)
          implicit none
          integer, parameter :: dp = kind(1.0d0) ! statt 'USE w90_constants'
          character(len=*), intent(in) :: seed_name_
          integer, dimension(3), intent(in) :: mp_grid_
          integer, intent(in) :: num_kpts_
          real(kind=dp), dimension(3,3), intent(in) :: real_lattice_
          real(kind=dp), dimension(3,3), intent(in) :: recip_lattice_
          real(kind=dp), dimension(3,num_kpts_), intent(in) :: kpt_latt_
          integer, intent(in) :: num_bands_
          integer, intent(in) :: num_wann_
          integer, intent(in) :: nntot_
          integer, intent(in) :: num_atoms_
          character(len=*), dimension(num_atoms_), intent(in) :: atom_symbols_
          real(kind=dp), dimension(3,num_atoms_), intent(in) :: atoms_cart_
          logical, intent(in) :: gamma_only_
          complex(kind=dp), dimension(num_bands_,num_bands_,nntot_,num_kpts_), intent(in) :: M_matrix_
          complex(kind=dp), dimension(num_bands_,num_wann_,num_kpts_), intent(in) :: A_matrix_
          real(kind=dp), dimension(num_bands_,num_kpts_), intent(in) :: eigenvalues_
          complex(kind=dp), dimension(num_wann_,num_wann_,num_kpts_), intent(out) :: U_matrix_
          complex(kind=dp), dimension(num_bands_,num_wann_,num_kpts_), optional, intent(out) :: U_matrix_opt_
          logical, dimension(num_bands_,num_kpts_), optional, intent(out) :: lwindow_
          real(kind=dp), dimension(3,num_wann_), optional, intent(out) :: wann_centres_
          real(kind=dp), dimension(num_wann_), optional, intent(out) :: wann_spreads_
          real(kind=dp), dimension(3), optional, intent(out) :: spread_
        end subroutine wannier_run
      END INTERFACE
#endif

      CONTAINS

!******************** SUBROUTINE MLWF_READER ***************************
!
!***********************************************************************
      SUBROUTINE MLWF_READER(IU5,IU6,IU0,LDOWNSAMPLE)
      USE base
      USE vaspxml
      USE full_kpoints
      IMPLICIT NONE
      INTEGER IU5,IU6,IU0
      LOGICAL LDOWNSAMPLE
      ! local variables
      INTEGER IDUM, N, IERR
      REAL(q) RDUM
      COMPLEX(q) CDUM
      LOGICAL LOPEN,LDUM,LRPA
      CHARACTER (1) :: CHARAC

      LOPEN=.FALSE.
      OPEN(UNIT=IU5,FILE=INCAR,STATUS='OLD')

      CALL RDATAB(LOPEN,INCAR,IU5,'LWANNIER90','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LWANNIER90,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
            WRITE(IU0,*)'Error reading item ''LWANNIER90'' from file INCAR.'
         LWANNIER90=.FALSE.
      ENDIF
      CALL XML_INCAR('LWANNIER90','L',IDUM,RDUM,CDUM,LWANNIER90,CHARAC,N)

      CALL RDATAB(LOPEN,INCAR,IU5,'LWANNIER90_RUN','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LWANNIER90_RUN,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
            WRITE(IU0,*)'Error reading item ''LWANNIER90_RUN'' from file INCAR.'
         LWANNIER90_RUN=.FALSE.
      ENDIF
      CALL XML_INCAR('LWANNIER90_RUN','L',IDUM,RDUM,CDUM,LWANNIER90_RUN,CHARAC,N)

      ! Do we want to write UWAN files?
      LWRITE_WANPROJ=.FALSE.
      CALL RDATAB(.FALSE.,INCAR,IU5,'LWRITE_WANPROJ','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LWRITE_WANPROJ,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
     &                    ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) &
            WRITE(IU0,*)'Error reading item ''LWRITE_WANPROJ'' from file INCAR.'
         LWRITE_WANPROJ=.FALSE.
      ENDIF
      IF ( LWRITE_WANPROJ .AND. .NOT.LWANNIER90_RUN) THEN
         IF(IU0>=0)WRITE(IU0,*)'LWANNIER90_RUN forced by LWRITE_WANPROJ'
         LWANNIER90_RUN=.TRUE.
      ENDIF

      !WANPROJ needs downsampling
      IF ( LWRITE_WANPROJ .AND. .NOT.LDOWNSAMPLE) THEN
         IF(IU0>=0)WRITE(IU0,*)'LDOWNSAMPLE forced by LWRITE_WANPROJ'
         LDOWNSAMPLE=.TRUE.
      ENDIF

      IF (WANNIER90()) THEN
         ! Switch on the PEAD routines
         CALL PEAD_REQUEST

         ! Do we want to write UNK files?
         WRITE_UNK=.FALSE.
         CALL RDATAB(.FALSE.,INCAR,IU5,'LWRITE_UNK','=','#',';','L', &
        &            IDUM,RDUM,CDUM,WRITE_UNK,CHARAC,N,1,IERR)
         IF (((IERR/=0).AND.(IERR/=3)).OR. &
        &                    ((IERR==0).AND.(N<1))) THEN
            IF (IU0>=0) &
               WRITE(IU0,*)'Error reading item ''LWRITE_UNK'' from file INCAR.'
            WRITE_UNK=.FALSE.
         ENDIF
         ! Do we want to write the mmn and amn files
         ! even though WANNIER90 runs in library mode?
         WRITE_MMN_AMN=.NOT.LWANNIER90_RUN
         CALL RDATAB(.FALSE.,INCAR,IU5,'LWRITE_MMN_AMN','=','#',';','L', &
        &            IDUM,RDUM,CDUM,WRITE_MMN_AMN,CHARAC,N,1,IERR)
         IF (((IERR/=0).AND.(IERR/=3)).OR. &
        &                    ((IERR==0).AND.(N<1))) THEN
            IF (IU0>=0) &
               WRITE(IU0,*)'Error reading item ''LWRITE_MMN_AMN'' from file INCAR.'
            WRITE_MMN_AMN=.NOT.LWANNIER90_RUN
         ENDIF
         ! Do we want to read the amn file instead of
         ! computing it anew?
         READ_AMN=.FALSE.
         CALL RDATAB(.FALSE.,INCAR,IU5,'LREAD_AMN','=','#',';','L', &
        &            IDUM,RDUM,CDUM,READ_AMN,CHARAC,N,1,IERR)
         IF (((IERR/=0).AND.(IERR/=3)).OR. &
        &                    ((IERR==0).AND.(N<1))) THEN
            IF (IU0>=0) &
               WRITE(IU0,*)'Error reading item ''LREAD_AMN'' from file INCAR.'
            READ_AMN=.FALSE.
         ENDIF
         ! reading the amn file does not make sense if
         ! we do not want to run WANNIER90 in library mode.
         READ_AMN=(READ_AMN.AND.LWANNIER90_RUN)

         !
         CALL USE_FULL_KPOINTS
      ENDIF

      CLOSE(IU5)

#ifndef VASP2WANNIER90
      IF (WANNIER90()) THEN
         IF (IU0>=0) WRITE(IU0,'(A)') &
        &    'MLWF_READER: ERROR: VASP was compiled without wannier90 library, exiting now ...'
         STOP
      ENDIF
#endif
      RETURN
      END SUBROUTINE MLWF_READER


!***********************************************************************
!
! WANNIER90
!
!***********************************************************************
      FUNCTION WANNIER90()
      IMPLICIT NONE
      LOGICAL WANNIER90
      WANNIER90=LWANNIER90.OR.LWANNIER90_RUN
      END FUNCTION WANNIER90


      FUNCTION WANNIER90RUN()
      IMPLICIT NONE
      LOGICAL WANNIER90RUN
      WANNIER90RUN=LWANNIER90_RUN
      END FUNCTION WANNIER90RUN

      SUBROUTINE WANNIER90_RUN_REQUEST
      IMPLICIT NONE
      LWANNIER90_RUN=.TRUE.
      RETURN
      END SUBROUTINE WANNIER90_RUN_REQUEST

!******************** SUBROUTINE MLWF_WANNIER90 ************************
!
!***********************************************************************
      SUBROUTINE MLWF_WANNIER90_AMN(&
     & WDES,W,P,CQIJ,T_INFO,LATT_CUR,INFO,IO &
     & )
      USE ini
      USE base
      USE constant
      USE pseudo
      USE poscar
      USE lattice
      USE wave_high
      USE radial
      USE full_kpoints

      IMPLICIT NONE

      TYPE(wavedes) WDES
      TYPE(wavespin) W
      TYPE(potcar) P(:)
      TYPE(type_info) T_INFO
      TYPE(latt) LATT_CUR
      TYPE(info_struct) INFO
      TYPE(in_struct) IO

      OVERLAP CQIJ(:,:,:,:)

      INTEGER, ALLOCATABLE :: band_win(:,:)
      INTEGER :: i,ratio,remainder
      integer, dimension(0:WDES%COMM_KINTER%NCPU - 1) :: counts
      integer, dimension(0:WDES%COMM_KINTER%NCPU - 1) :: displs
      COMPLEX(q), ALLOCATABLE :: A_matrix_loc(:,:,:,:)
      ! wannier90_setup variables: input
      CHARACTER(len=9) :: seed_name='wannier90'
      INTEGER :: mp_grid(3)
      REAL(q) :: real_lattice(3,3)
      REAL(q) :: recip_lattice(3,3)
      INTEGER :: num_bands_tot
      INTEGER :: num_atoms
      CHARACTER(len=20), ALLOCATABLE :: atom_symbols(:)
      REAL(q), ALLOCATABLE :: atoms_cart(:,:)
#ifdef gammareal
      LOGICAL :: gamma_only=.TRUE.
#else
      LOGICAL :: gamma_only=.FALSE.
#endif
      LOGICAL :: spinors
      INTEGER, PARAMETER :: num_nnmax=12
      ! wannier90_setup variables: output
      INTEGER :: nntot
      INTEGER, ALLOCATABLE :: nnlist(:,:)
      INTEGER, ALLOCATABLE :: nncell(:,:,:)
      INTEGER :: num_bands
      INTEGER :: num_wann
      REAL(q), ALLOCATABLE :: proj_site(:,:)
      INTEGER, ALLOCATABLE :: proj_l(:)
      INTEGER, ALLOCATABLE :: proj_m(:)
      INTEGER, ALLOCATABLE :: proj_radial(:)
      REAL(q), ALLOCATABLE :: proj_z(:,:)
      REAL(q), ALLOCATABLE :: proj_x(:,:)
      REAL(q), ALLOCATABLE :: proj_zona(:)
      INTEGER, ALLOCATABLE :: exclude_bands(:)
      ! wannier90_run variables: input
      COMPLEX(q), ALLOCATABLE :: M_matrix(:,:,:,:,:)
      COMPLEX(q), ALLOCATABLE :: A_matrix(:,:,:,:)
      REAL(q), ALLOCATABLE :: eigenvalues(:,:,:)

      ! local variables
      INTEGER NI,MI,NK,NKP,ICNTR
      INTEGER NKI,NKJ,ISP,ISPINOR,L,M,N,NP
      INTEGER NEXCLB,NPROJ
      REAL(q) POS(3)
      REAL(q) KI(3),KJ(3)
      INTEGER IDUM,IERR
      REAL(q) RDUM
      COMPLEX(q) CDUM,CPROJ
      LOGICAL LDUM,LWIN_FOUND,LAMN_FOUND
      CHARACTER (1) :: CHARAC

      TYPE(rgrid) R
      REAL(q) :: RSTART=0.0025_q,REND=10._q,H=0.025_q
      REAL(q), ALLOCATABLE :: FTMP(:),FR(:,:),FGTMP(:,:),FG(:,:,:)
      INTEGER, PARAMETER :: NMAX=1000

      COMPLEX(q),ALLOCATABLE ::  S(:,:)
      GDEF, ALLOCATABLE :: A(:,:),AP(:,:)
      REAL(q), ALLOCATABLE :: ROTYLM(:,:)
      REAL(q), ALLOCATABLE :: HYBRID_ORBITAL(:)

      LOGICAL, ALLOCATABLE :: EXCLUDE_BAND(:)
      LOGICAL LUSE_BLOCH_PHASES

      INTEGER, PARAMETER :: LMAX=3
      CHARACTER(LEN=11):: UNKFILE
      CHARACTER(LEN=2) :: SP(2)=(/"up","dn"/)

      INTEGER num_bands_on_file,num_kpts_on_file,NPROJ_on_file

      IF (.NOT.WANNIER90()) RETURN

#ifdef MPI
!     IF (WDES%COMM_KINTER%NCPU.NE.1) THEN
!        CALLMPI( M_stop('MLWF_WANNIER90: KPAR>1 not implemented, sorry.') )
!        STOP
!     END IF
#endif
      ! WANNIER90 cannot plot spinors
      IF (WRITE_UNK.AND.WDES%LNONCOLLINEAR) THEN
         IF (IO%IU0>=0) WRITE(IO%IU0,*) &
        & 'MLWF_WANNIER90: ERROR: will not write spinors to UNK files, sorry ...'
         WRITE_UNK=.FALSE.
      ENDIF
!
! allocation
      ALLOCATE(S(WDES%NB_TOT,WDES%NB_TOT))

      CALL CHECK_FULL_KPOINTS

      mp_grid(1)=KPOINTS_FULL%NKPX
      mp_grid(2)=KPOINTS_FULL%NKPY
      mp_grid(3)=KPOINTS_FULL%NKPZ

      num_kpts=mp_grid(1)*mp_grid(2)*mp_grid(3)

      !ALLOCATE(band_win(2,num_kpts))

      IF (ALLOCATED(kpt_latt)) DEALLOCATE(kpt_latt)
      ALLOCATE(kpt_latt(3,num_kpts))
      kpt_latt=KPOINTS_FULL%VKPT

      real_lattice(1,:)=LATT_CUR%A(:,1)
      real_lattice(2,:)=LATT_CUR%A(:,2)
      real_lattice(3,:)=LATT_CUR%A(:,3)

      recip_lattice(1,:)=LATT_CUR%B(:,1)
      recip_lattice(2,:)=LATT_CUR%B(:,2)
      recip_lattice(3,:)=LATT_CUR%B(:,3)

      recip_lattice=recip_lattice*TPI

      num_bands_tot=WDES%NB_TOT

      num_atoms=T_INFO%NIONS

      ALLOCATE(atom_symbols(num_atoms),atoms_cart(3,num_atoms))
      DO NI=1,T_INFO%NIONS
         atom_symbols(NI)=T_INFO%TYPE(T_INFO%ITYP(NI))
!        WRITE(atom_symbols(NI),'(A2,18X)') T_INFO%TYPE(T_INFO%ITYP(NI))

         POS(1)=LATT_CUR%A(1,1)*T_INFO%POSION(1,NI)+ &
        &       LATT_CUR%A(1,2)*T_INFO%POSION(2,NI)+ &
        &       LATT_CUR%A(1,3)*T_INFO%POSION(3,NI)
         POS(2)=LATT_CUR%A(2,1)*T_INFO%POSION(1,NI)+ &
        &       LATT_CUR%A(2,2)*T_INFO%POSION(2,NI)+ &
        &       LATT_CUR%A(2,3)*T_INFO%POSION(3,NI)
         POS(3)=LATT_CUR%A(3,1)*T_INFO%POSION(1,NI)+ &
        &       LATT_CUR%A(3,2)*T_INFO%POSION(2,NI)+ &
        &       LATT_CUR%A(3,3)*T_INFO%POSION(3,NI)
         atoms_cart(:,NI)=POS(:)
      ENDDO

      spinors=WDES%LNONCOLLINEAR

      ! A minimal wannier90.win file must exist; it must at least
      ! contain the keyword "num_wann"
      IF (IO%IU6>=0) THEN
         INQUIRE(FILE=seed_name//'.win',EXIST=LWIN_FOUND)
         IF (LWIN_FOUND) THEN
            OPEN(UNIT=99,FILE=seed_name//'.win',STATUS='OLD')
            CALL RDATAB(.FALSE.,seed_name//'.win',99,'num_wann','=','#',';','I', &
           &            IDUM,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
            IF (((IERR/=0).AND.(IERR/=3)).OR. &
           &                    ((IERR==0).AND.(N<1))) THEN
               IF (IO%IU0>=0) &
                  WRITE(IO%IU0,*)'Error reading item ''num_wann'' from file '//seed_name//'.win'
               STOP
            ENDIF
            IF (IERR==3) WRITE(99,'(A,I6,2X,A)') ' num_wann =',WDES%NB_TOT,'! set to NBANDS by VASP'
         ELSE
            OPEN(UNIT=99,FILE=seed_name//'.win',STATUS='REPLACE')
            WRITE(99,'(A,I6,2X,A)') ' num_wann =',WDES%NB_TOT,'! set to NBANDS by VASP'
         ENDIF
         CLOSE(99)
      ENDIF

      ! Initialize everything to zero
      nntot=0; num_bands=0; num_wann=0

      ALLOCATE(nnlist(num_kpts,num_nnmax),nncell(3,num_kpts,num_nnmax))
      nnlist=0; nncell=0

      ALLOCATE(proj_site(3,num_bands_tot),proj_l(num_bands_tot),proj_m(num_bands_tot), &
     &   proj_radial(num_bands_tot),proj_z(3,num_bands_tot),proj_x(3,num_bands_tot), &
     &   proj_zona(num_bands_tot),exclude_bands(num_bands_tot))
      proj_site=0; proj_l=0; proj_m=0; proj_radial=0; proj_z=0; proj_x=0; proj_zona=0; exclude_bands=0

      ! Only one node will do the actual work,
      ! otherwise all will write to wannier90.wout
#ifdef VASP2WANNIER90
      IF (IO%IU6>=0) &
      CALL wannier_setup(seed_name,mp_grid,num_kpts,real_lattice,recip_lattice, &
     &                   kpt_latt,num_bands_tot,num_atoms,atom_symbols,atoms_cart, &
     &                   gamma_only,spinors,nntot,nnlist,nncell,num_bands,num_wann, &
     &                   proj_site,proj_l,proj_m,proj_radial,proj_z,proj_x,proj_zona, &
     &                   exclude_bands)
#endif
      ! Now communicate the output to the other nodes
      CALLMPI( M_sum_i(WDES%COMM,nntot,1) )
      CALLMPI( M_sum_i(WDES%COMM,nnlist,num_kpts*num_nnmax) )
      CALLMPI( M_sum_i(WDES%COMM,nncell,3*num_kpts*num_nnmax) )
      CALLMPI( M_sum_i(WDES%COMM,num_bands,1) )
      CALLMPI( M_sum_i(WDES%COMM,num_wann,1) )
      CALLMPI( M_sum_d(WDES%COMM,proj_site,3*num_bands_tot) )
      CALLMPI( M_sum_i(WDES%COMM,proj_l,num_bands_tot) )
      CALLMPI( M_sum_i(WDES%COMM,proj_m,num_bands_tot) )
      CALLMPI( M_sum_i(WDES%COMM,proj_radial,num_bands_tot) )
      CALLMPI( M_sum_d(WDES%COMM,proj_z,3*num_bands_tot) )
      CALLMPI( M_sum_d(WDES%COMM,proj_x,3*num_bands_tot) )
      CALLMPI( M_sum_d(WDES%COMM,proj_zona,num_bands_tot) )
      CALLMPI( M_sum_i(WDES%COMM,exclude_bands,num_bands_tot) )

      ALLOCATE(EXCLUDE_BAND(num_bands_tot))
      EXCLUDE_BAND=.FALSE. ; NEXCLB=0
      DO N=1,num_bands_tot
         M=exclude_bands(N)
         IF (M>num_bands_tot.OR.M<0) THEN
            WRITE(*,*) 'MLWF_WANNIER90: ERROR: exclude_bands seems corrupt',N,M
         ELSEIF (M==0) THEN
            EXIT
         ELSE
            EXCLUDE_BAND(M)=.TRUE.
            NEXCLB=NEXCLB+1
         ENDIF
      ENDDO

      IF ((num_bands_tot-NEXCLB)/=num_bands) THEN
         IF (IO%IU0>=0) WRITE(*,*) 'MLWF_WANNIER90: ERROR: num_bands seems corrupt', &
        &   num_bands_tot-NEXCLB,num_bands
         STOP
      ENDIF

      IF (num_bands<num_wann) THEN
         IF (IO%IU0>=0) WRITE(*,*) 'MLWF_WANNIER90: ERROR: num_bands < num_wann', &
        &   num_bands,num_wann
         STOP
      ENDIF

#ifdef debug
      IF (IO%IU6>=0) THEN
!        WRITE(*,*) num_bands_tot,num_bands,num_wann
!        WRITE(*,*) EXCLUDE_BAND
         DO NK=1,num_kpts
            WRITE(*,'(I4,3F14.7)') NK,kpt_latt(:,NK)
            DO NKP=1,nntot
            WRITE(*,'(I4,3F14.7,3I4)') nnlist(NK,NKP),kpt_latt(:,nnlist(NK,NKP)),nncell(:,NK,NKP)
            ENDDO
            WRITE(*,*)
         ENDDO
      ENDIF
#endif

      ! Here the work starts in earnest
      ! We need to calculate M(m,n,j,i) = < u_m,k_i | u_n,k_j >
      ! where k_i = kpt_latt(:,i), and k_j = kpt_latt(:,nnlist(i,j))+nncell(:,i,j)
      ! U(m,n,j,i) is then written to the file seed_name.mmn

      ALLOCATE(A_matrix_loc(num_bands,num_wann,num_kpts,WDES%ISPIN))
      A_matrix_loc=0
      ALLOCATE(A_matrix(num_bands,num_wann,num_kpts,WDES%ISPIN))
      ALLOCATE(eigenvalues(num_bands,num_kpts,WDES%ISPIN))

      ! just a precaution, in case the compiler does not
      ! nullify pointers by default
      R%R=>NULL(); R%SI=>NULL()

      ! loop over spin
      spin: DO ISP=1,WDES%ISPIN

!      IF (IO%IU6>=0.AND.WRITE_MMN_AMN) THEN
!         IF (WDES%ISPIN==1) THEN
!            OPEN(UNIT=99,FILE=seed_name//'.mmn',STATUS='REPLACE')
!         ELSE
!            OPEN(UNIT=99,FILE=seed_name//'.'//SP(ISP)//'.mmn',STATUS='REPLACE')
!         ENDIF
!         WRITE(99,'(A)') 'File generated by VASP: '//INFO%SZNAM1
!         WRITE(99,'(3I12)') num_bands,num_kpts,nntot
!      ENDIF
!      ! runs over all k-points
!      ki_loop: DO NKI=1,num_kpts
!         KI(:)=kpt_latt(:,NKI)
!         ! runs over the number of nearest-neighbours of k_i
!         kj_loop: DO NKJ=1,nntot
!            KJ(:)=kpt_latt(:,nnlist(NKI,NKJ))+REAL(nncell(:,NKI,NKJ),q)
!#ifdef debug
!            IF (IO%IU6>=0) WRITE(*,'(2I4,3F14.7,2X,3F14.7)') NKI,NKJ,KI,KJ
!#endif
!            CALL PEAD_CALC_OVERLAP(W,KI,KJ,ISP,P,CQIJ,LATT_CUR,T_INFO,S,LQIJB=.TRUE.)
!
!            IF (IO%IU6>=0.AND.WRITE_MMN_AMN) WRITE(99,'(5I5)') NKI,nnlist(NKI,NKJ),nncell(:,NKI,NKJ)
!
!            NI=0
!            n_loop: DO N=1,WDES%NB_TOT
!               IF (EXCLUDE_BAND(N)) CYCLE n_loop
!               NI=NI+1; MI=0
!               m_loop: DO M=1,WDES%NB_TOT
!                  IF (EXCLUDE_BAND(M)) CYCLE m_loop
!                  MI=MI+1
!                  M_matrix(MI,NI,NKJ,NKI,ISP)=S(M,N)
!                  IF (IO%IU6>=0.AND.WRITE_MMN_AMN) WRITE(99,'(2F18.12)') S(M,N)
!               ENDDO m_loop
!            ENDDO n_loop
!            ! check consistency
!            IF (NI/=MI.OR.(NI/=num_bands)) THEN
!               IF (IO%IU0>=0) WRITE(*,*) 'MLWF_WANNIER90: ERROR: num_bands seems corrupt (2)',num_bands,NI
!               STOP
!            ENDIF
!
!         ENDDO kj_loop
!      ENDDO ki_loop
!
!      IF (IO%IU6>=0.AND.WRITE_MMN_AMN) CLOSE(99)


      ! and we might want to read an existing amn file
      IF (READ_AMN) THEN
      ! compute NPROJ
      NPROJ=0
      DO ISPINOR=1,WDES%NRSPINORS
      DO ICNTR=1,num_bands_tot
         IF (proj_l(ICNTR)==0.AND.proj_m(ICNTR)==0.AND.proj_radial(ICNTR)==0) CYCLE
         NPROJ=NPROJ+1
      ENDDO
      ENDDO
      ! check if amn file exists
      IF (WDES%ISPIN==1) THEN
         INQUIRE(FILE=seed_name//'.amn',EXIST=LAMN_FOUND)
      ELSE
         INQUIRE(FILE=seed_name//'.'//SP(ISP)//'.amn',EXIST=LAMN_FOUND)
      ENDIF
      IF (LAMN_FOUND) THEN
         ! open amn file
         IF (WDES%ISPIN==1) THEN
            OPEN(UNIT=99,FILE=seed_name//'.amn',STATUS='OLD')
         ELSE
            OPEN(UNIT=99,FILE=seed_name//'.'//SP(ISP)//'.amn',STATUS='OLD')
         ENDIF
         ! read header
         READ(99,*)
         READ(99,*) num_bands_on_file,num_kpts_on_file,NPROJ_on_file
         ! check dimensions
         IF (num_bands_on_file/=num_bands.OR.num_kpts_on_file/=num_kpts.OR.NPROJ_on_file/=NPROJ) THEN
            IF (IO%IU0>=0) WRITE(*,'(A,6I5,/A)') &
           &   'MLWF_WANNIER90: ERROR: amn file does not apply to present calculation:', &
           &   num_bands_on_file,num_bands,num_kpts_on_file,num_kpts,NPROJ_on_file,NPROJ, &
           &   'defaulting to LREAD_AMN=.FALSE.'
            READ_AMN=.FALSE.
         ELSE
            A_matrix=0
            DO NKI=1,num_kpts_on_file
            DO N=1,NPROJ_on_file
            DO M=1,num_bands_on_file
               READ(99,*) IDUM,IDUM,IDUM,A_matrix(M,N,NKI,ISP)
            ENDDO
            ENDDO
            ENDDO
         ENDIF
         CLOSE(99)
      ELSE
         IF (IO%IU0>=0) WRITE(*,*)'MLWF_WANNIER90: ERROR: amn file not found, defaulting to LREAD_AMN=.FALSE.'
         READ_AMN=.FALSE.
      ENDIF
      ENDIF ! finished reading A_matrix

      IF (.NOT.READ_AMN) THEN
      ! we need to calculate the projection of the wave functions
      ! onto a set of trial functions {g}: A(m,n,i) = < \psi_m,k_i | g_n >
      ! A(m,n,i) is then written to the file seed_name.amn
#ifdef debug
      IF (IO%IU6>=0) THEN
         DO ICNTR=1,num_bands_tot
            WRITE(*,'(I4,3F14.7,3I4)')  &
           &   ICNTR,proj_site(1:3,ICNTR),proj_l(ICNTR),proj_m(ICNTR),proj_radial(ICNTR)
         ENDDO
      ENDIF
#endif
      ALLOCATE(A(WDES%NB_TOT,(LMAX+1)**2),AP(WDES%NB_TOT,(LMAX+1)**2))
      ALLOCATE(ROTYLM((LMAX+1)**2,(LMAX+1)**2),HYBRID_ORBITAL((LMAX+1)**2))

      NPROJ=0
!     CALL START_TIMING("T1")
      ! runs over all projector sites
      spinor: DO ISPINOR=1,WDES%NRSPINORS
      sites : DO ICNTR=1,num_bands_tot
         IF (proj_l(ICNTR)==0.AND.proj_m(ICNTR)==0.AND.proj_radial(ICNTR)==0) CYCLE sites
         NPROJ=NPROJ+1
         ! setup the Ylm rotation matrix in accordance with proj_z and proj_x
         CALL SETROTYLM(proj_x(:,ICNTR),proj_z(:,ICNTR),LMAX,ROTYLM)
#ifdef debug
         IF (IO%IU6>=0) THEN
            WRITE(*,'(I4,3F14.7)') ICNTR,proj_site(1:3,ICNTR)
            DO N=10,MIN((LMAX+1)**2,16)
               WRITE(*,'(9F10.5)') (ROTYLM(N,M),M=10,MIN((LMAX+1)**2,16))
            ENDDO
         ENDIF
#endif
         ! translate between VASP and the orbital definition of wannier90
         CALL WANNIER90_ORBITAL_DEFINITIONS(proj_l(ICNTR),proj_m(ICNTR),HYBRID_ORBITAL)
         ! setup the radial functions (real space)
         CALL SETRGRID(RSTART/proj_zona(ICNTR),REND,H,R)
         ALLOCATE(FTMP(R%NMAX),FR(R%NMAX,LMAX+1))
         ALLOCATE(FGTMP(NMAX,5),FG(NMAX,5,LMAX+1))
         ! For now we use the same radial function for all L
         DO L=0,LMAX
            CALL RADIAL_FUNCTION(proj_radial(ICNTR),R,proj_zona(ICNTR),FTMP)
            CALL BESSEL_TRANSFORM_RADIAL_FUNCTION(L,R,FTMP,REAL(SQRT(2._q*INFO%ENMAX/HSQDTM)/NMAX,KIND=q),FGTMP)
            FR(:,L+1)=FTMP; FG(:,:,L+1)=FGTMP
         ENDDO


         !write(*,*) WDES%COMM_KINTER%NCPU, WDES%COMM_KINTER%NODE_ME, WDES%NB_LOW,WDES%NB_PAR,WDES%NBANDS
         ratio = num_kpts/WDES%COMM_KINTER%NCPU
         remainder = MOD(num_kpts, WDES%COMM_KINTER%NCPU)

         do i = 0, WDES%COMM_KINTER%NCPU - 1
           if (i < remainder) then
             counts(i) = ratio + 1
             displs(i) = i*(ratio + 1)
           else
             counts(i) = ratio
             displs(i) = remainder*(ratio + 1) + (i - remainder)*ratio
           end if
         end do
         !write(*,*) WDES%COMM_KINTER%NODE_ME,counts(WDES%COMM_KINTER%NODE_ME-1),displs(WDES%COMM_KINTER%NODE_ME-1)
         !STOP

         CALL START_TIMING("T2")
         ! runs over all k-points
         !kpoints : DO NKI=1,num_kpts
         kpoints : DO NKI=displs(WDES%COMM_KINTER%NODE_ME-1)+1,displs(WDES%COMM_KINTER%NODE_ME-1)+counts(WDES%COMM_KINTER%NODE_ME-1)
           ! CALL CALC_OVERLAP_GN( &
           !&   LMAX,FG,proj_site(:,ICNTR),W,kpt_latt(:,NKI),ISP,ISPINOR,P,CQIJ,LATT_CUR,T_INFO,A)
            CALL CALC_OVERLAP_GN_band( &
           &   band_win_(1,NKI),band_win_(2,NKI),LMAX,FG,proj_site(:,ICNTR),W,kpt_latt(:,NKI),ISP,ISPINOR,P,CQIJ,LATT_CUR,T_INFO,A)
            ! and rotate them in accordance with proj_z and proj_x
            AP=0
            DO M=1,WDES%NB_TOT
            DO N=1,(LMAX+1)**2
               DO NP=1,(LMAX+1)**2
                  AP(M,N)=AP(M,N)+A(M,NP)*ROTYLM(N,NP)
               ENDDO
            ENDDO
            ENDDO
#ifdef debug
            IF (IO%IU6>=0) THEN
               WRITE(*,'(I4,3F14.7)') NKI,kpt_latt(1:3,NKI)
               DO M=1,WDES%NB_TOT
                  WRITE(*,'(4(2F14.7,2X))') (A(M,N),N=1,4)
               ENDDO
               WRITE(*,*)
               DO M=1,WDES%NB_TOT
                  WRITE(*,'(4(2F14.7,2X))') (AP(M,N),N=1,4)
               ENDDO
               WRITE(*,*)
            ENDIF
#endif
            ! Make the desired linear combinations
            MI=0
            DO M=1,WDES%NB_TOT
               IF (EXCLUDE_BAND(M)) CYCLE
               CPROJ=0
               DO N=1,(LMAX+1)**2
                  CPROJ=CPROJ+AP(M,N)*HYBRID_ORBITAL(N)
               ENDDO
               MI=MI+1
               A_matrix_loc(MI,NPROJ,NKI,ISP)=CPROJ!/WDES%NB_PAR
            ENDDO

         ENDDO kpoints

         DEALLOCATE(FTMP,FR,FGTMP,FG)

         !IF (IO%IU0>=0) WRITE(IO%IU0,'(X,A,I3,A)') 'Projection ',NPROJ,' done.'
!        CALL STOP_TIMING("T2",IO%IU6,'  OVL')
!        CALL STOP_TIMING("T1",IO%IU6,'ICNTR')
      ENDDO sites
      ENDDO spinor

      DEALLOCATE(A,AP,ROTYLM,HYBRID_ORBITAL)

      ! check consistency
      !IF (IO%IU0>=0) write(*,*) EXCLUDE_BAND
      !IF (IO%IU0>=0) write(*,*) MI, WDES%NB_TOT, num_bands
      IF (MI/=num_bands) THEN
         IF (IO%IU0>=0) WRITE(*,*) 'MLWF_WANNIER90: ERROR: num_bands seems corrupt (3)',num_bands,MI
         STOP
      ENDIF
      IF (NPROJ/=0.AND.NPROJ/=num_wann) THEN
         IF (IO%IU0>=0) WRITE(*,*) 'MLWF_WANNIER90: ERROR: number of projections not equal to num_wann',NPROJ,num_wann
         STOP
      ENDIF

      ENDIF ! finished computing A_matrix

      A_matrix=0
      CALL MPI_ALLREDUCE(A_matrix_loc,A_matrix,num_bands*num_wann*num_kpts*WDES%ISPIN,MPI_DOUBLE_COMPLEX,MPI_SUM,WDES%COMM%MPI_COMM,ierr)

      ! write A_matrix to file
      LUSE_BLOCH_PHASES=.FALSE.
      IF (NPROJ/=0) THEN
         IF (IO%IU6>=0.AND.WRITE_MMN_AMN.AND..NOT.READ_AMN) THEN
            IF (WDES%ISPIN==1) THEN
               OPEN(UNIT=99,FILE=seed_name//'.amn',STATUS='REPLACE')
            ELSE
               OPEN(UNIT=99,FILE=seed_name//'.'//SP(ISP)//'.amn',STATUS='REPLACE')
            ENDIF
            WRITE(99,'(A)') 'File generated by VASP: '//INFO%SZNAM1
            WRITE(99,'(3I12)') num_bands,num_kpts,NPROJ
            DO NKI=1,num_kpts
            DO N=1,NPROJ
            DO M=1,num_bands
               WRITE(99,'(3I5,2F18.12)') M,N,NKI,A_matrix(M,N,NKI,ISP)
            ENDDO
            ENDDO
            ENDDO
            CLOSE(99)
         ENDIF
      ELSE
         ! read use_bloch_phases from the .win file
         OPEN(UNIT=99,FILE=seed_name//'.win',STATUS='OLD')
         CALL RDATAB(.FALSE.,seed_name//'.win',99,'use_bloch_phases','=','#',';','L', &
        &            IDUM,RDUM,CDUM,LUSE_BLOCH_PHASES,CHARAC,N,1,IERR)
         IF (((IERR/=0).AND.(IERR/=3)).OR. &
        &                    ((IERR==0).AND.(N<1))) THEN
            IF (IO%IU0>=0) &
               WRITE(IO%IU0,*)'Error reading item ''use_bloch_phases'' from file '//seed_name//'.win'
            STOP
         ENDIF
         IF (IERR==3) THEN
            ! if it is not found, then set it
            IF (IO%IU6>=0) WRITE(99,'(/A)') 'use_bloch_phases = T'
            LUSE_BLOCH_PHASES=.TRUE.
         ENDIF
         CLOSE(99)
      ENDIF
      ! if we do not have starting projections and
      ! use_bloch_phases=.false. in the .win file
      ! we will have a problem running WANNIER90
      ! in library mode

      ! and ...
      ! we should write the eigenvalues onto seed_name.eig
      !IF (IO%IU6>=0.AND.WRITE_MMN_AMN) THEN
      IF (IO%IU6>=0) THEN
         IF (WDES%ISPIN==1) THEN
            OPEN(UNIT=99,FILE=seed_name//'.eig',STATUS='REPLACE')
         ELSE
            OPEN(UNIT=99,FILE=seed_name//'.'//SP(ISP)//'.eig',STATUS='REPLACE')
         ENDIF
      ENDIF

      DO NKI=1,num_kpts
         NK=KPOINT_IN_FULL_GRID(kpt_latt(:,NKI),KPOINTS_FULL)
         NK=KPOINTS_FULL%NEQUIV(NK)
         MI=0
         DO M=1,WDES%NB_TOT
            IF (EXCLUDE_BAND(M)) CYCLE
            MI=MI+1
            eigenvalues(MI,NKI,ISP)=REAL(W%CELTOT(M,NK,ISP))
            !IF (IO%IU6>0.AND.WRITE_MMN_AMN) WRITE(99,'(2I12,F22.12)') MI,NKI,REAL(W%CELTOT(M,NK,ISP))
            IF (IO%IU6>0) WRITE(99,'(2I12,F22.12)') MI,NKI,REAL(W%CELTOT(M,NK,ISP))
         ENDDO
      ENDDO
      !IF (IO%IU6>=0.AND.WRITE_MMN_AMN) CLOSE(99)
      IF (IO%IU6>=0) CLOSE(99)

      ENDDO spin

      ! Here we can try to run WANNIER90 in library mode
      IF (WANNIER90RUN()) THEN
         IF (ALLOCATED(U_matrix)) DEALLOCATE( &
        &   U_matrix,U_matrix_opt,lwindow,lexclude_band,wann_centres,wann_spreads,spread)

         ALLOCATE(U_matrix(num_wann,num_wann,num_kpts,WDES%ISPIN), &
        &   U_matrix_opt(num_bands,num_wann,num_kpts,WDES%ISPIN), &
        &   lwindow(num_bands,num_kpts,WDES%ISPIN), &
        &   lexclude_band(num_bands_tot), &
        &   wann_centres(3,num_wann,WDES%ISPIN), &
        &   wann_spreads(num_wann,WDES%ISPIN), &
        &   spread(3,WDES%ISPIN))

         U_matrix=0._q; U_matrix_opt=0._q; lwindow=.FALSE.; lexclude_band=.FALSE.
         wann_centres=0._q; wann_spreads=0._q; spread=0._q

         ! Execute wannier_run only on the masternode
#ifdef VASP2WANNIER90
         IF (IO%IU6>=0) THEN
            IF (WDES%ISPIN==1) THEN
               CALL wannier_run(seed_name,mp_grid,num_kpts,real_lattice,recip_lattice, &
              &                 kpt_latt,num_bands,num_wann,nntot,num_atoms,atom_symbols, &
              &                 atoms_cart,gamma_only,&
              &                 M_matrix(:,:,:,:,1),A_matrix(:,:,:,1),eigenvalues(:,:,1), &
              &                 U_matrix(:,:,:,1),U_matrix_opt(:,:,:,1),lwindow(:,:,1), &
              &                 wann_centres(:,:,1),wann_spreads(:,1),spread(:,1))
           ELSE
               DO ISP=1,WDES%ISPIN
               ! copy wannier90.win to wannier90.(up/dn).win
               CALL COPYFILE(seed_name//'.win',seed_name//'.'//SP(ISP)//'.win')

               ! run wannier90 for spin component ISP
               CALL wannier_run(seed_name//'.'//SP(ISP), &
              &                 mp_grid,num_kpts,real_lattice,recip_lattice, &
              &                 kpt_latt,num_bands,num_wann,nntot,num_atoms,atom_symbols, &
              &                 atoms_cart,gamma_only,&
              &                 M_matrix(:,:,:,:,ISP),A_matrix(:,:,:,ISP),eigenvalues(:,:,ISP), &
              &                 U_matrix(:,:,:,ISP),U_matrix_opt(:,:,:,ISP),lwindow(:,:,ISP), &
              &                 wann_centres(:,:,ISP),wann_spreads(:,ISP),spread(:,ISP))

               ENDDO
            ENDIF
         ENDIF
#endif
         MLWF_num_wann=num_wann
         ! Communicate the results to the other nodes
         CALLMPI( M_sum_z(WDES%COMM,U_matrix,num_wann*num_wann*num_kpts*WDES%ISPIN) )
         CALLMPI( M_sum_z(WDES%COMM,U_matrix_opt,num_bands*num_wann*num_kpts*WDES%ISPIN) )
         CALLMPI( M_bcast_l(WDES%COMM,lwindow,num_bands*num_kpts*WDES%ISPIN) )
         CALLMPI( M_sum_d(WDES%COMM,wann_centres,3*num_wann*WDES%ISPIN) )
         CALLMPI( M_sum_d(WDES%COMM,wann_spreads,num_wann*WDES%ISPIN) )
         CALLMPI( M_sum_d(WDES%COMM,spread,3*WDES%ISPIN) )
         ! and store a copy of EXCLUDE_BAND
         lexclude_band(:)=EXCLUDE_BAND(:)
         DEALLOCATE(wann_centres,wann_spreads,spread)

         !possibly write UWAN files
         IF (LWRITE_WANPROJ .AND. IO%IU6>=0) THEN
            WRITE(IO%IU0,*)' Writing '//WANPROJFILE//' file'
            CALL WRITE_WANPROJ_FILE(WDES)
         ENDIF
      ENDIF

      DEALLOCATE(atom_symbols,atoms_cart)
      DEALLOCATE(nnlist,nncell)
      DEALLOCATE(proj_site,proj_l,proj_m,proj_radial,proj_z,proj_x,proj_zona,exclude_bands)
      DEALLOCATE(EXCLUDE_BAND)
      DEALLOCATE(A_matrix,A_matrix_loc,eigenvalues)
      DEALLOCATE(S)

      IF (ASSOCIATED(R%R)) THEN
         DEALLOCATE(R%R); NULLIFY(R%R)
      ENDIF
      IF (ASSOCIATED(R%SI)) THEN
         DEALLOCATE(R%SI); NULLIFY(R%SI)
      ENDIF

      RETURN
      END SUBROUTINE MLWF_WANNIER90_AMN

      SUBROUTINE MLWF_WANNIER90_EIG(&
     & WDES,W,P,CQIJ,T_INFO,LATT_CUR,INFO,IO &
     & )
      USE ini
      USE base
      USE constant
      USE pseudo
      USE poscar
      USE lattice
      USE wave_high
      USE radial
      USE full_kpoints

      IMPLICIT NONE

      TYPE(wavedes) WDES
      TYPE(wavespin) W
      TYPE(potcar) P(:)
      TYPE(type_info) T_INFO
      TYPE(latt) LATT_CUR
      TYPE(info_struct) INFO
      TYPE(in_struct) IO

      OVERLAP CQIJ(:,:,:,:)

      INTEGER, ALLOCATABLE :: band_win(:,:)
      INTEGER :: i,ratio,remainder
      integer, dimension(0:WDES%COMM_KINTER%NCPU - 1) :: counts
      integer, dimension(0:WDES%COMM_KINTER%NCPU - 1) :: displs
      COMPLEX(q), ALLOCATABLE :: A_matrix_loc(:,:,:,:)
      ! wannier90_setup variables: input
      CHARACTER(len=9) :: seed_name='wannier90'
      INTEGER :: mp_grid(3)
      REAL(q) :: real_lattice(3,3)
      REAL(q) :: recip_lattice(3,3)
      INTEGER :: num_bands_tot
      INTEGER :: num_atoms
      CHARACTER(len=20), ALLOCATABLE :: atom_symbols(:)
      REAL(q), ALLOCATABLE :: atoms_cart(:,:)
#ifdef gammareal
      LOGICAL :: gamma_only=.TRUE.
#else
      LOGICAL :: gamma_only=.FALSE.
#endif
      LOGICAL :: spinors
      INTEGER, PARAMETER :: num_nnmax=12
      ! wannier90_setup variables: output
      INTEGER :: nntot
      INTEGER, ALLOCATABLE :: nnlist(:,:)
      INTEGER, ALLOCATABLE :: nncell(:,:,:)
      INTEGER :: num_bands
      INTEGER :: num_wann
      REAL(q), ALLOCATABLE :: proj_site(:,:)
      INTEGER, ALLOCATABLE :: proj_l(:)
      INTEGER, ALLOCATABLE :: proj_m(:)
      INTEGER, ALLOCATABLE :: proj_radial(:)
      REAL(q), ALLOCATABLE :: proj_z(:,:)
      REAL(q), ALLOCATABLE :: proj_x(:,:)
      REAL(q), ALLOCATABLE :: proj_zona(:)
      INTEGER, ALLOCATABLE :: exclude_bands(:)
      ! wannier90_run variables: input
      COMPLEX(q), ALLOCATABLE :: M_matrix(:,:,:,:,:)
      COMPLEX(q), ALLOCATABLE :: A_matrix(:,:,:,:)
      REAL(q), ALLOCATABLE :: eigenvalues(:,:,:)

      ! local variables
      INTEGER NI,MI,NK,NKP,ICNTR
      INTEGER NKI,NKJ,ISP,ISPINOR,L,M,N,NP
      INTEGER NEXCLB,NPROJ
      REAL(q) POS(3)
      REAL(q) KI(3),KJ(3)
      INTEGER IDUM,IERR
      REAL(q) RDUM
      COMPLEX(q) CDUM,CPROJ
      LOGICAL LDUM,LWIN_FOUND,LAMN_FOUND
      CHARACTER (1) :: CHARAC

      TYPE(rgrid) R
      REAL(q) :: RSTART=0.0025_q,REND=10._q,H=0.025_q
      REAL(q), ALLOCATABLE :: FTMP(:),FR(:,:),FGTMP(:,:),FG(:,:,:)
      INTEGER, PARAMETER :: NMAX=1000

      COMPLEX(q),ALLOCATABLE ::  S(:,:)
      GDEF, ALLOCATABLE :: A(:,:),AP(:,:)
      REAL(q), ALLOCATABLE :: ROTYLM(:,:)
      REAL(q), ALLOCATABLE :: HYBRID_ORBITAL(:)

      LOGICAL, ALLOCATABLE :: EXCLUDE_BAND(:)
      LOGICAL LUSE_BLOCH_PHASES

      INTEGER, PARAMETER :: LMAX=3
      CHARACTER(LEN=11):: UNKFILE
      CHARACTER(LEN=2) :: SP(2)=(/"up","dn"/)

      INTEGER num_bands_on_file,num_kpts_on_file,NPROJ_on_file

      IF (.NOT.WANNIER90()) RETURN

      CALL CHECK_FULL_KPOINTS

      mp_grid(1)=KPOINTS_FULL%NKPX
      mp_grid(2)=KPOINTS_FULL%NKPY
      mp_grid(3)=KPOINTS_FULL%NKPZ
      num_kpts=mp_grid(1)*mp_grid(2)*mp_grid(3)

      IF (ALLOCATED(kpt_latt)) DEALLOCATE(kpt_latt)
      ALLOCATE(kpt_latt(3,num_kpts))
      kpt_latt=KPOINTS_FULL%VKPT

      IF (IO%IU6>=0) THEN
         IF (WDES%ISPIN==1) THEN
            OPEN(UNIT=99,FILE='wannier90.eig',STATUS='REPLACE')
         ENDIF
      ENDIF


      DO NKI=1,num_kpts
         NK=KPOINT_IN_FULL_GRID(kpt_latt(:,NKI),KPOINTS_FULL)
         NK=KPOINTS_FULL%NEQUIV(NK)
         MI=0
         DO M=1,WDES%NB_TOT
            MI=MI+1
            IF (IO%IU6>0) WRITE(99,'(2I12,F22.12)') MI,NKI,REAL(W%CELTOT(M,NK,1))
         ENDDO
      ENDDO
      IF (IO%IU6>=0) CLOSE(99)
      DEALLOCATE(kpt_latt)
      END SUBROUTINE MLWF_WANNIER90_EIG

      SUBROUTINE MLWF_WANNIER90(&
     & WDES,W,P,CQIJ,T_INFO,LATT_CUR,INFO,IO &
     & )
      USE ini
      USE base
      USE constant
      USE pseudo
      USE poscar
      USE lattice
      USE wave_high
      USE radial
      USE full_kpoints

      IMPLICIT NONE

      TYPE(wavedes) WDES
      TYPE(wavespin) W
      TYPE(potcar) P(:)
      TYPE(type_info) T_INFO
      TYPE(latt) LATT_CUR
      TYPE(info_struct) INFO
      TYPE(in_struct) IO

      OVERLAP CQIJ(:,:,:,:)

      INTEGER, ALLOCATABLE :: band_win(:,:)
      INTEGER :: i,ratio,remainder
      integer, dimension(0:WDES%COMM_KINTER%NCPU - 1) :: counts
      integer, dimension(0:WDES%COMM_KINTER%NCPU - 1) :: displs
      COMPLEX(q), ALLOCATABLE :: A_matrix_loc(:,:,:,:)
      ! wannier90_setup variables: input
      CHARACTER(len=9) :: seed_name='wannier90'
      INTEGER :: mp_grid(3)
      REAL(q) :: real_lattice(3,3)
      REAL(q) :: recip_lattice(3,3)
      INTEGER :: num_bands_tot
      INTEGER :: num_atoms
      CHARACTER(len=20), ALLOCATABLE :: atom_symbols(:)
      REAL(q), ALLOCATABLE :: atoms_cart(:,:)
#ifdef gammareal
      LOGICAL :: gamma_only=.TRUE.
#else
      LOGICAL :: gamma_only=.FALSE.
#endif
      LOGICAL :: spinors
      INTEGER, PARAMETER :: num_nnmax=12
      ! wannier90_setup variables: output
      INTEGER :: nntot
      INTEGER, ALLOCATABLE :: nnlist(:,:)
      INTEGER, ALLOCATABLE :: nncell(:,:,:)
      INTEGER :: num_bands
      INTEGER :: num_wann
      REAL(q), ALLOCATABLE :: proj_site(:,:)
      INTEGER, ALLOCATABLE :: proj_l(:)
      INTEGER, ALLOCATABLE :: proj_m(:)
      INTEGER, ALLOCATABLE :: proj_radial(:)
      REAL(q), ALLOCATABLE :: proj_z(:,:)
      REAL(q), ALLOCATABLE :: proj_x(:,:)
      REAL(q), ALLOCATABLE :: proj_zona(:)
      INTEGER, ALLOCATABLE :: exclude_bands(:)
      ! wannier90_run variables: input
      COMPLEX(q), ALLOCATABLE :: M_matrix(:,:,:,:,:)
      COMPLEX(q), ALLOCATABLE :: A_matrix(:,:,:,:)
      REAL(q), ALLOCATABLE :: eigenvalues(:,:,:)

      ! local variables
      INTEGER NI,MI,NK,NKP,ICNTR
      INTEGER NKI,NKJ,ISP,ISPINOR,L,M,N,NP
      INTEGER NEXCLB,NPROJ
      REAL(q) POS(3)
      REAL(q) KI(3),KJ(3)
      INTEGER IDUM,IERR
      REAL(q) RDUM
      COMPLEX(q) CDUM,CPROJ
      LOGICAL LDUM,LWIN_FOUND,LAMN_FOUND
      CHARACTER (1) :: CHARAC

      TYPE(rgrid) R
      REAL(q) :: RSTART=0.0025_q,REND=10._q,H=0.025_q
      REAL(q), ALLOCATABLE :: FTMP(:),FR(:,:),FGTMP(:,:),FG(:,:,:)
      INTEGER, PARAMETER :: NMAX=1000

      COMPLEX(q),ALLOCATABLE ::  S(:,:)
      GDEF, ALLOCATABLE :: A(:,:),AP(:,:)
      REAL(q), ALLOCATABLE :: ROTYLM(:,:)
      REAL(q), ALLOCATABLE :: HYBRID_ORBITAL(:)

      LOGICAL, ALLOCATABLE :: EXCLUDE_BAND(:)
      LOGICAL LUSE_BLOCH_PHASES

      INTEGER, PARAMETER :: LMAX=3
      CHARACTER(LEN=11):: UNKFILE
      CHARACTER(LEN=2) :: SP(2)=(/"up","dn"/)

      INTEGER num_bands_on_file,num_kpts_on_file,NPROJ_on_file

      IF (.NOT.WANNIER90()) RETURN

#ifdef MPI
!     IF (WDES%COMM_KINTER%NCPU.NE.1) THEN
!        CALLMPI( M_stop('MLWF_WANNIER90: KPAR>1 not implemented, sorry.') )
!        STOP
!     END IF
#endif
      ! WANNIER90 cannot plot spinors
      IF (WRITE_UNK.AND.WDES%LNONCOLLINEAR) THEN
         IF (IO%IU0>=0) WRITE(IO%IU0,*) &
        & 'MLWF_WANNIER90: ERROR: will not write spinors to UNK files, sorry ...'
         WRITE_UNK=.FALSE.
      ENDIF
!
! allocation
      ALLOCATE(S(WDES%NB_TOT,WDES%NB_TOT))

      CALL CHECK_FULL_KPOINTS

      mp_grid(1)=KPOINTS_FULL%NKPX
      mp_grid(2)=KPOINTS_FULL%NKPY
      mp_grid(3)=KPOINTS_FULL%NKPZ

      num_kpts=mp_grid(1)*mp_grid(2)*mp_grid(3)

      IF (ALLOCATED(kpt_latt)) DEALLOCATE(kpt_latt)
      ALLOCATE(kpt_latt(3,num_kpts))
      kpt_latt=KPOINTS_FULL%VKPT

      real_lattice(1,:)=LATT_CUR%A(:,1)
      real_lattice(2,:)=LATT_CUR%A(:,2)
      real_lattice(3,:)=LATT_CUR%A(:,3)

      recip_lattice(1,:)=LATT_CUR%B(:,1)
      recip_lattice(2,:)=LATT_CUR%B(:,2)
      recip_lattice(3,:)=LATT_CUR%B(:,3)

      recip_lattice=recip_lattice*TPI

      num_bands_tot=WDES%NB_TOT

      num_atoms=T_INFO%NIONS

      ALLOCATE(atom_symbols(num_atoms),atoms_cart(3,num_atoms))
      DO NI=1,T_INFO%NIONS
         atom_symbols(NI)=T_INFO%TYPE(T_INFO%ITYP(NI))
!        WRITE(atom_symbols(NI),'(A2,18X)') T_INFO%TYPE(T_INFO%ITYP(NI))

         POS(1)=LATT_CUR%A(1,1)*T_INFO%POSION(1,NI)+ &
        &       LATT_CUR%A(1,2)*T_INFO%POSION(2,NI)+ &
        &       LATT_CUR%A(1,3)*T_INFO%POSION(3,NI)
         POS(2)=LATT_CUR%A(2,1)*T_INFO%POSION(1,NI)+ &
        &       LATT_CUR%A(2,2)*T_INFO%POSION(2,NI)+ &
        &       LATT_CUR%A(2,3)*T_INFO%POSION(3,NI)
         POS(3)=LATT_CUR%A(3,1)*T_INFO%POSION(1,NI)+ &
        &       LATT_CUR%A(3,2)*T_INFO%POSION(2,NI)+ &
        &       LATT_CUR%A(3,3)*T_INFO%POSION(3,NI)
         atoms_cart(:,NI)=POS(:)
      ENDDO

      spinors=WDES%LNONCOLLINEAR

      ! A minimal wannier90.win file must exist; it must at least
      ! contain the keyword "num_wann"
      IF (IO%IU6>=0) THEN
         INQUIRE(FILE=seed_name//'.win',EXIST=LWIN_FOUND)
         IF (LWIN_FOUND) THEN
            OPEN(UNIT=99,FILE=seed_name//'.win',STATUS='OLD')
            CALL RDATAB(.FALSE.,seed_name//'.win',99,'num_wann','=','#',';','I', &
           &            IDUM,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
            IF (((IERR/=0).AND.(IERR/=3)).OR. &
           &                    ((IERR==0).AND.(N<1))) THEN
               IF (IO%IU0>=0) &
                  WRITE(IO%IU0,*)'Error reading item ''num_wann'' from file '//seed_name//'.win'
               STOP
            ENDIF
            IF (IERR==3) WRITE(99,'(A,I6,2X,A)') ' num_wann =',WDES%NB_TOT,'! set to NBANDS by VASP'
         ELSE
            OPEN(UNIT=99,FILE=seed_name//'.win',STATUS='REPLACE')
            WRITE(99,'(A,I6,2X,A)') ' num_wann =',WDES%NB_TOT,'! set to NBANDS by VASP'
         ENDIF
         CLOSE(99)
      ENDIF

      ! Initialize everything to zero
      nntot=0; num_bands=0; num_wann=0

      ALLOCATE(nnlist(num_kpts,num_nnmax),nncell(3,num_kpts,num_nnmax))
      nnlist=0; nncell=0

      ALLOCATE(proj_site(3,num_bands_tot),proj_l(num_bands_tot),proj_m(num_bands_tot), &
     &   proj_radial(num_bands_tot),proj_z(3,num_bands_tot),proj_x(3,num_bands_tot), &
     &   proj_zona(num_bands_tot),exclude_bands(num_bands_tot))
      proj_site=0; proj_l=0; proj_m=0; proj_radial=0; proj_z=0; proj_x=0; proj_zona=0; exclude_bands=0

      ! Only one node will do the actual work,
      ! otherwise all will write to wannier90.wout
#ifdef VASP2WANNIER90
      IF (IO%IU6>=0) &
      CALL wannier_setup(seed_name,mp_grid,num_kpts,real_lattice,recip_lattice, &
     &                   kpt_latt,num_bands_tot,num_atoms,atom_symbols,atoms_cart, &
     &                   gamma_only,spinors,nntot,nnlist,nncell,num_bands,num_wann, &
     &                   proj_site,proj_l,proj_m,proj_radial,proj_z,proj_x,proj_zona, &
     &                   exclude_bands)
#endif
      ! Now communicate the output to the other nodes
      CALLMPI( M_sum_i(WDES%COMM,nntot,1) )
      CALLMPI( M_sum_i(WDES%COMM,nnlist,num_kpts*num_nnmax) )
      CALLMPI( M_sum_i(WDES%COMM,nncell,3*num_kpts*num_nnmax) )
      CALLMPI( M_sum_i(WDES%COMM,num_bands,1) )
      CALLMPI( M_sum_i(WDES%COMM,num_wann,1) )
      CALLMPI( M_sum_d(WDES%COMM,proj_site,3*num_bands_tot) )
      CALLMPI( M_sum_i(WDES%COMM,proj_l,num_bands_tot) )
      CALLMPI( M_sum_i(WDES%COMM,proj_m,num_bands_tot) )
      CALLMPI( M_sum_i(WDES%COMM,proj_radial,num_bands_tot) )
      CALLMPI( M_sum_d(WDES%COMM,proj_z,3*num_bands_tot) )
      CALLMPI( M_sum_d(WDES%COMM,proj_x,3*num_bands_tot) )
      CALLMPI( M_sum_d(WDES%COMM,proj_zona,num_bands_tot) )
      CALLMPI( M_sum_i(WDES%COMM,exclude_bands,num_bands_tot) )

      ALLOCATE(EXCLUDE_BAND(num_bands_tot))
      EXCLUDE_BAND=.FALSE. ; NEXCLB=0
      DO N=1,num_bands_tot
         M=exclude_bands(N)
         IF (M>num_bands_tot.OR.M<0) THEN
            WRITE(*,*) 'MLWF_WANNIER90: ERROR: exclude_bands seems corrupt',N,M
         ELSEIF (M==0) THEN
            EXIT
         ELSE
            EXCLUDE_BAND(M)=.TRUE.
            NEXCLB=NEXCLB+1
         ENDIF
      ENDDO

      IF ((num_bands_tot-NEXCLB)/=num_bands) THEN
         IF (IO%IU0>=0) WRITE(*,*) 'MLWF_WANNIER90: ERROR: num_bands seems corrupt', &
        &   num_bands_tot-NEXCLB,num_bands
         STOP
      ENDIF

      IF (num_bands<num_wann) THEN
         IF (IO%IU0>=0) WRITE(*,*) 'MLWF_WANNIER90: ERROR: num_bands < num_wann', &
        &   num_bands,num_wann
         STOP
      ENDIF

#ifdef debug
      IF (IO%IU6>=0) THEN
!        WRITE(*,*) num_bands_tot,num_bands,num_wann
!        WRITE(*,*) EXCLUDE_BAND
         DO NK=1,num_kpts
            WRITE(*,'(I4,3F14.7)') NK,kpt_latt(:,NK)
            DO NKP=1,nntot
            WRITE(*,'(I4,3F14.7,3I4)') nnlist(NK,NKP),kpt_latt(:,nnlist(NK,NKP)),nncell(:,NK,NKP)
            ENDDO
            WRITE(*,*)
         ENDDO
      ENDIF
#endif

      ! Here the work starts in earnest
      ! We need to calculate M(m,n,j,i) = < u_m,k_i | u_n,k_j >
      ! where k_i = kpt_latt(:,i), and k_j = kpt_latt(:,nnlist(i,j))+nncell(:,i,j)
      ! U(m,n,j,i) is then written to the file seed_name.mmn

      ALLOCATE(A_matrix_loc(num_bands,num_wann,num_kpts,WDES%ISPIN))
      A_matrix_loc=0
      ALLOCATE(M_matrix(num_bands,num_bands,nntot,num_kpts,WDES%ISPIN))
      ALLOCATE(A_matrix(num_bands,num_wann,num_kpts,WDES%ISPIN))
      ALLOCATE(eigenvalues(num_bands,num_kpts,WDES%ISPIN))

      ! just a precaution, in case the compiler does not
      ! nullify pointers by default
      R%R=>NULL(); R%SI=>NULL()

      ! loop over spin
      spin: DO ISP=1,WDES%ISPIN

      IF (IO%IU6>=0.AND.WRITE_MMN_AMN) THEN
         IF (WDES%ISPIN==1) THEN
            OPEN(UNIT=99,FILE=seed_name//'.mmn',STATUS='REPLACE')
         ELSE
            OPEN(UNIT=99,FILE=seed_name//'.'//SP(ISP)//'.mmn',STATUS='REPLACE')
         ENDIF
         WRITE(99,'(A)') 'File generated by VASP: '//INFO%SZNAM1
         WRITE(99,'(3I12)') num_bands,num_kpts,nntot
      ENDIF
      ! runs over all k-points
      ki_loop: DO NKI=1,num_kpts
         KI(:)=kpt_latt(:,NKI)
         ! runs over the number of nearest-neighbours of k_i
         kj_loop: DO NKJ=1,nntot
            KJ(:)=kpt_latt(:,nnlist(NKI,NKJ))+REAL(nncell(:,NKI,NKJ),q)
#ifdef debug
            IF (IO%IU6>=0) WRITE(*,'(2I4,3F14.7,2X,3F14.7)') NKI,NKJ,KI,KJ
#endif
            CALL PEAD_CALC_OVERLAP(W,KI,KJ,ISP,P,CQIJ,LATT_CUR,T_INFO,S,LQIJB=.TRUE.)

            IF (IO%IU6>=0.AND.WRITE_MMN_AMN) WRITE(99,'(5I5)') NKI,nnlist(NKI,NKJ),nncell(:,NKI,NKJ)

            NI=0
            n_loop: DO N=1,WDES%NB_TOT
               IF (EXCLUDE_BAND(N)) CYCLE n_loop
               NI=NI+1; MI=0
               m_loop: DO M=1,WDES%NB_TOT
                  IF (EXCLUDE_BAND(M)) CYCLE m_loop
                  MI=MI+1
                  M_matrix(MI,NI,NKJ,NKI,ISP)=S(M,N)
                  IF (IO%IU6>=0.AND.WRITE_MMN_AMN) WRITE(99,'(2F18.12)') S(M,N)
               ENDDO m_loop
            ENDDO n_loop
            ! check consistency
            IF (NI/=MI.OR.(NI/=num_bands)) THEN
               IF (IO%IU0>=0) WRITE(*,*) 'MLWF_WANNIER90: ERROR: num_bands seems corrupt (2)',num_bands,NI
               STOP
            ENDIF

         ENDDO kj_loop
      ENDDO ki_loop

      IF (IO%IU6>=0.AND.WRITE_MMN_AMN) CLOSE(99)

      ! and we might want to read an existing amn file
      IF (READ_AMN) THEN
      ! compute NPROJ
      NPROJ=0
      DO ISPINOR=1,WDES%NRSPINORS
      DO ICNTR=1,num_bands_tot
         IF (proj_l(ICNTR)==0.AND.proj_m(ICNTR)==0.AND.proj_radial(ICNTR)==0) CYCLE
         NPROJ=NPROJ+1
      ENDDO
      ENDDO
      ! check if amn file exists
      IF (WDES%ISPIN==1) THEN
         INQUIRE(FILE=seed_name//'.amn',EXIST=LAMN_FOUND)
      ELSE
         INQUIRE(FILE=seed_name//'.'//SP(ISP)//'.amn',EXIST=LAMN_FOUND)
      ENDIF
      IF (LAMN_FOUND) THEN
         ! open amn file
         IF (WDES%ISPIN==1) THEN
            OPEN(UNIT=99,FILE=seed_name//'.amn',STATUS='OLD')
         ELSE
            OPEN(UNIT=99,FILE=seed_name//'.'//SP(ISP)//'.amn',STATUS='OLD')
         ENDIF
         ! read header
         READ(99,*)
         READ(99,*) num_bands_on_file,num_kpts_on_file,NPROJ_on_file
         ! check dimensions
         IF (num_bands_on_file/=num_bands.OR.num_kpts_on_file/=num_kpts.OR.NPROJ_on_file/=NPROJ) THEN
            IF (IO%IU0>=0) WRITE(*,'(A,6I5,/A)') &
           &   'MLWF_WANNIER90: ERROR: amn file does not apply to present calculation:', &
           &   num_bands_on_file,num_bands,num_kpts_on_file,num_kpts,NPROJ_on_file,NPROJ, &
           &   'defaulting to LREAD_AMN=.FALSE.'
            READ_AMN=.FALSE.
         ELSE
            A_matrix=0
            DO NKI=1,num_kpts_on_file
            DO N=1,NPROJ_on_file
            DO M=1,num_bands_on_file
               READ(99,*) IDUM,IDUM,IDUM,A_matrix(M,N,NKI,ISP)
            ENDDO
            ENDDO
            ENDDO
         ENDIF
         CLOSE(99)
      ELSE
         IF (IO%IU0>=0) WRITE(*,*)'MLWF_WANNIER90: ERROR: amn file not found, defaulting to LREAD_AMN=.FALSE.'
         READ_AMN=.FALSE.
      ENDIF
      ENDIF ! finished reading A_matrix

      IF (.NOT.READ_AMN) THEN
      ! we need to calculate the projection of the wave functions
      ! onto a set of trial functions {g}: A(m,n,i) = < \psi_m,k_i | g_n >
      ! A(m,n,i) is then written to the file seed_name.amn
#ifdef debug
      IF (IO%IU6>=0) THEN
         DO ICNTR=1,num_bands_tot
            WRITE(*,'(I4,3F14.7,3I4)')  &
           &   ICNTR,proj_site(1:3,ICNTR),proj_l(ICNTR),proj_m(ICNTR),proj_radial(ICNTR)
         ENDDO
      ENDIF
#endif
      ALLOCATE(A(WDES%NB_TOT,(LMAX+1)**2),AP(WDES%NB_TOT,(LMAX+1)**2))
      ALLOCATE(ROTYLM((LMAX+1)**2,(LMAX+1)**2),HYBRID_ORBITAL((LMAX+1)**2))

      NPROJ=0
!     CALL START_TIMING("T1")
      ! runs over all projector sites
      spinor: DO ISPINOR=1,WDES%NRSPINORS
      sites : DO ICNTR=1,num_bands_tot
         IF (proj_l(ICNTR)==0.AND.proj_m(ICNTR)==0.AND.proj_radial(ICNTR)==0) CYCLE sites
         NPROJ=NPROJ+1
         ! setup the Ylm rotation matrix in accordance with proj_z and proj_x
         CALL SETROTYLM(proj_x(:,ICNTR),proj_z(:,ICNTR),LMAX,ROTYLM)
#ifdef debug
         IF (IO%IU6>=0) THEN
            WRITE(*,'(I4,3F14.7)') ICNTR,proj_site(1:3,ICNTR)
            DO N=10,MIN((LMAX+1)**2,16)
               WRITE(*,'(9F10.5)') (ROTYLM(N,M),M=10,MIN((LMAX+1)**2,16))
            ENDDO
         ENDIF
#endif
         ! translate between VASP and the orbital definition of wannier90
         CALL WANNIER90_ORBITAL_DEFINITIONS(proj_l(ICNTR),proj_m(ICNTR),HYBRID_ORBITAL)
         ! setup the radial functions (real space)
         CALL SETRGRID(RSTART/proj_zona(ICNTR),REND,H,R)
         ALLOCATE(FTMP(R%NMAX),FR(R%NMAX,LMAX+1))
         ALLOCATE(FGTMP(NMAX,5),FG(NMAX,5,LMAX+1))
         ! For now we use the same radial function for all L
         DO L=0,LMAX
            CALL RADIAL_FUNCTION(proj_radial(ICNTR),R,proj_zona(ICNTR),FTMP)
            CALL BESSEL_TRANSFORM_RADIAL_FUNCTION(L,R,FTMP,REAL(SQRT(2._q*INFO%ENMAX/HSQDTM)/NMAX,KIND=q),FGTMP)
            FR(:,L+1)=FTMP; FG(:,:,L+1)=FGTMP
         ENDDO

         ratio = num_kpts/WDES%COMM_KINTER%NCPU
         remainder = MOD(num_kpts, WDES%COMM_KINTER%NCPU)

         do i = 0, WDES%COMM_KINTER%NCPU - 1
           if (i < remainder) then
             counts(i) = ratio + 1
             displs(i) = i*(ratio + 1)
           else
             counts(i) = ratio
             displs(i) = remainder*(ratio + 1) + (i - remainder)*ratio
           end if
         end do
         !write(*,*) WDES%COMM_KINTER%NODE_ME,counts(WDES%COMM_KINTER%NODE_ME-1),displs(WDES%COMM_KINTER%NODE_ME-1)
         !STOP

         CALL START_TIMING("T2")
         ! runs over all k-points
         !kpoints : DO NKI=1,num_kpts
         kpoints : DO NKI=displs(WDES%COMM_KINTER%NODE_ME-1)+1,displs(WDES%COMM_KINTER%NODE_ME-1)+counts(WDES%COMM_KINTER%NODE_ME-1)

            !CALL CALC_OVERLAP_GN_band( &
           !&   15,42,LMAX,FG,proj_site(:,ICNTR),W,kpt_latt(:,NKI),ISP,ISPINOR,P,CQIJ,LATT_CUR,T_INFO,A)
            CALL CALC_OVERLAP_GN( &
           &   LMAX,FG,proj_site(:,ICNTR),W,kpt_latt(:,NKI),ISP,ISPINOR,P,CQIJ,LATT_CUR,T_INFO,A)
            ! and rotate them in accordance with proj_z and proj_x
            AP=0
            DO M=1,WDES%NB_TOT
            DO N=1,(LMAX+1)**2
               DO NP=1,(LMAX+1)**2
                  AP(M,N)=AP(M,N)+A(M,NP)*ROTYLM(N,NP)
               ENDDO
            ENDDO
            ENDDO
#ifdef debug
            IF (IO%IU6>=0) THEN
               WRITE(*,'(I4,3F14.7)') NKI,kpt_latt(1:3,NKI)
               DO M=1,WDES%NB_TOT
                  WRITE(*,'(4(2F14.7,2X))') (A(M,N),N=1,4)
               ENDDO
               WRITE(*,*)
               DO M=1,WDES%NB_TOT
                  WRITE(*,'(4(2F14.7,2X))') (AP(M,N),N=1,4)
               ENDDO
               WRITE(*,*)
            ENDIF
#endif
            ! Make the desired linear combinations
            MI=0
            DO M=1,WDES%NB_TOT
               IF (EXCLUDE_BAND(M)) CYCLE
               CPROJ=0
               DO N=1,(LMAX+1)**2
                  CPROJ=CPROJ+AP(M,N)*HYBRID_ORBITAL(N)
               ENDDO
               MI=MI+1
               A_matrix_loc(MI,NPROJ,NKI,ISP)=CPROJ
            ENDDO

         ENDDO kpoints

         DEALLOCATE(FTMP,FR,FGTMP,FG)

         IF (IO%IU0>=0) WRITE(IO%IU0,'(X,A,I3,A)') 'Projection ',NPROJ,' done.'
!        CALL STOP_TIMING("T2",IO%IU6,'  OVL')
!        CALL STOP_TIMING("T1",IO%IU6,'ICNTR')
      ENDDO sites
      ENDDO spinor

      DEALLOCATE(A,AP,ROTYLM,HYBRID_ORBITAL)

      ! check consistency
      IF (MI/=num_bands) THEN
         IF (IO%IU0>=0) WRITE(*,*) 'MLWF_WANNIER90: ERROR: num_bands seems corrupt (3)',num_bands,MI
         STOP
      ENDIF
      IF (NPROJ/=0.AND.NPROJ/=num_wann) THEN
         IF (IO%IU0>=0) WRITE(*,*) 'MLWF_WANNIER90: ERROR: number of projections not equal to num_wann',NPROJ,num_wann
         STOP
      ENDIF

      ENDIF ! finished computing A_matrix

      A_matrix=0
      CALL MPI_ALLREDUCE(A_matrix_loc,A_matrix,num_bands*num_wann*num_kpts*WDES%ISPIN,MPI_DOUBLE_COMPLEX,MPI_SUM,WDES%COMM%MPI_COMM,ierr)

      ! write A_matrix to file
      LUSE_BLOCH_PHASES=.FALSE.
      IF (NPROJ/=0) THEN
         IF (IO%IU6>=0.AND.WRITE_MMN_AMN.AND..NOT.READ_AMN) THEN
            IF (WDES%ISPIN==1) THEN
               OPEN(UNIT=99,FILE=seed_name//'.amn',STATUS='REPLACE')
            ELSE
               OPEN(UNIT=99,FILE=seed_name//'.'//SP(ISP)//'.amn',STATUS='REPLACE')
            ENDIF
            WRITE(99,'(A)') 'File generated by VASP: '//INFO%SZNAM1
            WRITE(99,'(3I12)') num_bands,num_kpts,NPROJ
            DO NKI=1,num_kpts
            DO N=1,NPROJ
            DO M=1,num_bands
               WRITE(99,'(3I5,2F18.12)') M,N,NKI,A_matrix(M,N,NKI,ISP)
            ENDDO
            ENDDO
            ENDDO
            CLOSE(99)
         ENDIF
      ELSE
         ! read use_bloch_phases from the .win file
         OPEN(UNIT=99,FILE=seed_name//'.win',STATUS='OLD')
         CALL RDATAB(.FALSE.,seed_name//'.win',99,'use_bloch_phases','=','#',';','L', &
        &            IDUM,RDUM,CDUM,LUSE_BLOCH_PHASES,CHARAC,N,1,IERR)
         IF (((IERR/=0).AND.(IERR/=3)).OR. &
        &                    ((IERR==0).AND.(N<1))) THEN
            IF (IO%IU0>=0) &
               WRITE(IO%IU0,*)'Error reading item ''use_bloch_phases'' from file '//seed_name//'.win'
            STOP
         ENDIF
         IF (IERR==3) THEN
            ! if it is not found, then set it
            IF (IO%IU6>=0) WRITE(99,'(/A)') 'use_bloch_phases = T'
            LUSE_BLOCH_PHASES=.TRUE.
         ENDIF
         CLOSE(99)
      ENDIF
      ! if we do not have starting projections and
      ! use_bloch_phases=.false. in the .win file
      ! we will have a problem running WANNIER90
      ! in library mode
      IF (NPROJ==0.AND.(.NOT.LUSE_BLOCH_PHASES).AND.WANNIER90RUN()) THEN
!
      ENDIF

      ! and ...
      ! we should write the eigenvalues onto seed_name.eig
      !IF (IO%IU6>=0.AND.WRITE_MMN_AMN) THEN
      IF (IO%IU6>=0) THEN
         IF (WDES%ISPIN==1) THEN
            OPEN(UNIT=99,FILE=seed_name//'.eig',STATUS='REPLACE')
         ELSE
            OPEN(UNIT=99,FILE=seed_name//'.'//SP(ISP)//'.eig',STATUS='REPLACE')
         ENDIF
      ENDIF

      DO NKI=1,num_kpts
         NK=KPOINT_IN_FULL_GRID(kpt_latt(:,NKI),KPOINTS_FULL)
         NK=KPOINTS_FULL%NEQUIV(NK)
         MI=0
         DO M=1,WDES%NB_TOT
            IF (EXCLUDE_BAND(M)) CYCLE
            MI=MI+1
            eigenvalues(MI,NKI,ISP)=REAL(W%CELTOT(M,NK,ISP))
            !IF (IO%IU6>0.AND.WRITE_MMN_AMN) WRITE(99,'(2I12,F22.12)') MI,NKI,REAL(W%CELTOT(M,NK,ISP))
            IF (IO%IU6>0) WRITE(99,'(2I12,F22.12)') MI,NKI,REAL(W%CELTOT(M,NK,ISP))
         ENDDO
      ENDDO
      !IF (IO%IU6>=0.AND.WRITE_MMN_AMN) CLOSE(99)
      IF (IO%IU6>=0) CLOSE(99)

      ! and last but not least we may need to write UNK files
      ! which is a bit tricky since the bands at a certain k-point
      ! more often than not do not reside on the same node and
      ! to make matters worse, WANNIER90 works with the full k-mesh
      ! whereas VASP uses the symmetry reduced one.
 100  FORMAT('UNK',I5.5,'.',I1)
      IF (WRITE_UNK) THEN
         DO NKI=1,num_kpts
            IF (IO%IU6>=0) THEN
               WRITE(UNKFILE,100) NKI,ISP
               OPEN(UNIT=99,FILE=UNKFILE,FORM='UNFORMATTED',STATUS='REPLACE')
               WRITE(99) W%WDES%GRID%NGX,W%WDES%GRID%NGY,W%WDES%GRID%NGZ,NKI,num_bands
            ENDIF
            CALL WRITE_WAVE_FUNCTIONS(W,kpt_latt(:,NKI),ISP,EXCLUDE_BAND,P,LATT_CUR,99)
            IF (IO%IU6>=0) CLOSE(99)
         ENDDO
      ENDIF

      ENDDO spin

      ! Write information on seed_name.win
      IF (IO%IU6>=0) THEN
         OPEN(UNIT=99,FILE=seed_name//'.win',STATUS='OLD')

         IF (spinors) THEN
            CALL OCCURS_IN_FILE(99,'spinors',N)
            IF (N==0) THEN
               WRITE(99,'(/A)') 'spinors = .true.'
            ELSE
               WRITE(*,'(A,/A/)') 'MLWF_WANNIER90: WARNING: '//seed_name//'.win seems to contains', &
              &                  'the SPINORS tag, I hope it applies to the current setup ...'
            ENDIF
         ENDIF

         CALL OCCURS_IN_FILE(99,'unit_cell_cart',N)
         IF (N==0) THEN
            WRITE(99,'(/A)') 'begin unit_cell_cart'
            WRITE(99,'(3F14.7)') real_lattice(1,:)
            WRITE(99,'(3F14.7)') real_lattice(2,:)
            WRITE(99,'(3F14.7)') real_lattice(3,:)
            WRITE(99,'(A)')  'end unit_cell_cart'
         ELSE
            WRITE(*,'(A,/A/)') 'MLWF_WANNIER90: WARNING: '//seed_name//'.win seems to contains', &
           &                  'a UNIT_CELL_CART block already, I hope it applies to the current setup ...'
         ENDIF

         CALL OCCURS_IN_FILE(99,'atoms_frac',N)
         IF (N==0) THEN
            CALL OCCURS_IN_FILE(99,'atoms_cart',N)
            IF (N==0) THEN
               WRITE(99,'(/A)') 'begin atoms_cart'
               DO NI=1,num_atoms
                  WRITE(99,'(A2,2X,3F14.7)') atom_symbols(NI),atoms_cart(:,NI)
               ENDDO
               WRITE(99,'(A)')  'end atoms_cart'
            ELSE
               WRITE(*,'(A,/A/)') 'MLWF_WANNIER90: WARNING: '//seed_name//'.win seems to contains', &
              &                  'a ATOMS_CART block already, I hope it applies to the current setup ...'
            ENDIF
         ELSE
            WRITE(*,'(A,/A/)') 'MLWF_WANNIER90: WARNING: '//seed_name//'.win seems to contains', &
           &                  'a ATOMS_FRAC block already, I hope it applies to the current setup ...'
         ENDIF

         CALL OCCURS_IN_FILE(99,'mp_grid',N)
         IF (N==0) THEN
            WRITE(99,'(/A,3I6)') 'mp_grid =',mp_grid
         ELSE
            WRITE(*,'(A,/A/)') 'MLWF_WANNIER90: WARNING: '//seed_name//'.win seems to contains', &
           &                  'the MP_GRID tag already, I hope it applies to the current setup ...'
         ENDIF

         CALL OCCURS_IN_FILE(99,'kpoints',N)
         IF (N==0) THEN
            WRITE(99,'(/A)') 'begin kpoints'
            DO NK=1,num_kpts
               WRITE(99,'(3F20.12)') kpt_latt(:,NK)
            ENDDO
            WRITE(99,'(A)')  'end kpoints'
         ELSE
            WRITE(*,'(A,/A/)') 'MLWF_WANNIER90: WARNING: '//seed_name//'.win seems to contains', &
           &                  'a KPOINTS block already, I hope it applies to the current setup ...'
         ENDIF

         CLOSE(99)
      ENDIF

      ! Here we can try to run WANNIER90 in library mode
      IF (WANNIER90RUN()) THEN
         IF (ALLOCATED(U_matrix)) DEALLOCATE( &
        &   U_matrix,U_matrix_opt,lwindow,lexclude_band,wann_centres,wann_spreads,spread)

         ALLOCATE(U_matrix(num_wann,num_wann,num_kpts,WDES%ISPIN), &
        &   U_matrix_opt(num_bands,num_wann,num_kpts,WDES%ISPIN), &
        &   lwindow(num_bands,num_kpts,WDES%ISPIN), &
        &   lexclude_band(num_bands_tot), &
        &   wann_centres(3,num_wann,WDES%ISPIN), &
        &   wann_spreads(num_wann,WDES%ISPIN), &
        &   spread(3,WDES%ISPIN))

         U_matrix=0._q; U_matrix_opt=0._q; lwindow=.FALSE.; lexclude_band=.FALSE.
         wann_centres=0._q; wann_spreads=0._q; spread=0._q

         ! Execute wannier_run only on the masternode
#ifdef VASP2WANNIER90
         IF (IO%IU6>=0) THEN
            IF (WDES%ISPIN==1) THEN
               CALL wannier_run(seed_name,mp_grid,num_kpts,real_lattice,recip_lattice, &
              &                 kpt_latt,num_bands,num_wann,nntot,num_atoms,atom_symbols, &
              &                 atoms_cart,gamma_only,&
              &                 M_matrix(:,:,:,:,1),A_matrix(:,:,:,1),eigenvalues(:,:,1), &
              &                 U_matrix(:,:,:,1),U_matrix_opt(:,:,:,1),lwindow(:,:,1), &
              &                 wann_centres(:,:,1),wann_spreads(:,1),spread(:,1))
           ELSE
               DO ISP=1,WDES%ISPIN
               ! copy wannier90.win to wannier90.(up/dn).win
               CALL COPYFILE(seed_name//'.win',seed_name//'.'//SP(ISP)//'.win')

               ! run wannier90 for spin component ISP
               CALL wannier_run(seed_name//'.'//SP(ISP), &
              &                 mp_grid,num_kpts,real_lattice,recip_lattice, &
              &                 kpt_latt,num_bands,num_wann,nntot,num_atoms,atom_symbols, &
              &                 atoms_cart,gamma_only,&
              &                 M_matrix(:,:,:,:,ISP),A_matrix(:,:,:,ISP),eigenvalues(:,:,ISP), &
              &                 U_matrix(:,:,:,ISP),U_matrix_opt(:,:,:,ISP),lwindow(:,:,ISP), &
              &                 wann_centres(:,:,ISP),wann_spreads(:,ISP),spread(:,ISP))

               ENDDO
            ENDIF
         ENDIF
#endif
         MLWF_num_wann=num_wann
         ! Communicate the results to the other nodes
         CALLMPI( M_sum_z(WDES%COMM,U_matrix,num_wann*num_wann*num_kpts*WDES%ISPIN) )
         CALLMPI( M_sum_z(WDES%COMM,U_matrix_opt,num_bands*num_wann*num_kpts*WDES%ISPIN) )
         CALLMPI( M_bcast_l(WDES%COMM,lwindow,num_bands*num_kpts*WDES%ISPIN) )
         CALLMPI( M_sum_d(WDES%COMM,wann_centres,3*num_wann*WDES%ISPIN) )
         CALLMPI( M_sum_d(WDES%COMM,wann_spreads,num_wann*WDES%ISPIN) )
         CALLMPI( M_sum_d(WDES%COMM,spread,3*WDES%ISPIN) )
         ! and store a copy of EXCLUDE_BAND
         lexclude_band(:)=EXCLUDE_BAND(:)
         DEALLOCATE(wann_centres,wann_spreads,spread)

         !possibly write UWAN files
         IF (LWRITE_WANPROJ .AND. IO%IU6>=0) THEN
            WRITE(IO%IU0,*)' Writing '//WANPROJFILE//' file'
            CALL WRITE_WANPROJ_FILE(WDES)
         ENDIF
      ENDIF

      DEALLOCATE(atom_symbols,atoms_cart)
      DEALLOCATE(nnlist,nncell)
      DEALLOCATE(proj_site,proj_l,proj_m,proj_radial,proj_z,proj_x,proj_zona,exclude_bands)
      DEALLOCATE(EXCLUDE_BAND)
      DEALLOCATE(M_matrix,A_matrix,A_matrix_loc,eigenvalues)
      DEALLOCATE(S)

      IF (ASSOCIATED(R%R)) THEN
         DEALLOCATE(R%R); NULLIFY(R%R)
      ENDIF
      IF (ASSOCIATED(R%SI)) THEN
         DEALLOCATE(R%SI); NULLIFY(R%SI)
      ENDIF

      RETURN
      END SUBROUTINE MLWF_WANNIER90


!******************** SUBROUTINE MLWF_ROTATE_ORBITALS ******************
!
!***********************************************************************
      SUBROUTINE MLWF_ROTATE_ORBITALS( &
     &   WDES,W,KPOINTS,GRID,T_INFO,P,NONL_S,SYMM,LATT_CUR,IO)
      USE base
      USE lattice
      USE pseudo
      USE poscar
      USE mgrid
      USE msymmetry
      USE nonl_high
      USE wave_high
      USE full_kpoints
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      TYPE(kpoints_struct) KPOINTS
      TYPE(grid_3d) GRID
      TYPE(type_info) T_INFO
      TYPE(potcar) P(T_INFO%NTYP)
      TYPE(nonl_struct) NONL_S
      TYPE(symmetry) SYMM
      TYPE(latt) LATT_CUR
      TYPE(in_struct) IO
      ! local variables

      IF (SYMM%ISYM>=0.AND. .NOT.WDES%LGAMMA) THEN
         ! switch of symmetry
         CALL NOSYMM(LATT_CUR%A,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS,SYMM%PTRANS, &
        &   SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%MAGROT,WDES%ISPIN,IO%IU6)
         ! reread k-points with LINVERSION=.FALSE. to generate full mesh
         CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR,.FALSE., &
        &   T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IO%IU6,IO%IU0)
         CALL KPAR_SYNC_ALL(WDES,W)
         CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_CUR,-1, IO%IU0)
         CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
      ENDIF

      CALL MLWF_ROTATE_ORBITALS_FULLK(WDES,W)

      RETURN
      END SUBROUTINE MLWF_ROTATE_ORBITALS


!******************** SUBROUTINE MLWF_ROTATE_ORBITALS_FULLK ************
!
!***********************************************************************
      SUBROUTINE MLWF_ROTATE_ORBITALS_FULLK(WDES,W)
      USE dfast
      USE wave_high
      USE full_kpoints
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      ! local variables
      TYPE(wavedes1) WDES1
      TYPE(wavefuna) WA
      INTEGER I,J,K,IP,IPP
      INTEGER ISP,NK,NKP,NW,NB,NBEXCL,NBwin
      COMPLEX(q) CTMP
      GDEF, ALLOCATABLE :: U(:,:)

      CALL CHECK_FULL_KPOINTS

      IF ((.NOT.ALLOCATED(U_matrix)).OR.(.NOT.ALLOCATED(U_matrix_opt))) THEN
         WRITE(*,*) 'MLWF_ROTATE_ORBITALS_FULLK: ERROR: rotation matrices not available'
         STOP
      ENDIF

      NW=SIZE(U_matrix,1)
      NB=SIZE(U_matrix_opt,1)
      IF (SIZE(lexclude_band)/=WDES%NB_TOT) THEN
         WRITE(*,*) 'MLWF_ROTATE_ORBITALS_FULLK: ERROR: lexclude_band array is corrupt:', &
        &   SIZE(lexclude_band),WDES%NB_TOT
         STOP
      ENDIF
      NBEXCL=0
      DO I=1,SIZE(lexclude_band)
         IF (lexclude_band(I)) NBEXCL=NBEXCL+1
      ENDDO
      IF ((NB+NBEXCL)/=WDES%NB_TOT) THEN
         WRITE(*,*) 'MLWF_ROTATE_ORBITALS_FULLK: ERROR: inconsistent number of bands:', &
        &   NB+NBEXCL,WDES%NB_TOT
      ENDIF

      ALLOCATE(U(WDES%NB_TOT,WDES%NB_TOT))

      spin: DO ISP=1,WDES%ISPIN
         kpoints: DO NK=1,KPOINTS_FULL%NKPTS
            ! find the corresponding entry in KPOINTS_FULL_ORIG, that contains
            ! the set of k-points used in the computation of the rotation matrices
            DO NKP=1,KPOINTS_FULL_ORIG%NKPTS
               IF (LIDENTICAL_KPOINT(KPOINTS_FULL%VKPT(:,NK),KPOINTS_FULL_ORIG%VKPT(:,NKP))) exit
            ENDDO
            IF (NKP>KPOINTS_FULL_ORIG%NKPTS) THEN
               WRITE(*,*) 'MLWF_ROTATE_ORBITALS_FULLK: ERROR: no matching k-point found in KPOINTS_FULL_ORIG',NK
               STOP
            ENDIF
            ! setup the effective rotation matrix U at the present k-point
            U=0._q; IP=0; IPP=0
            bands: DO I=1,WDES%NB_TOT
               IF (lexclude_band(I)) CYCLE bands
               IP=IP+1
               IF (.NOT.lwindow(IP,NKP,ISP)) CYCLE bands
               IPP=IPP+1
               DO J=1,NW
                  CTMP=0._q
                  DO K=1,NW
                     CTMP=CTMP+U_matrix_opt(IPP,K,NKP,ISP)*U_matrix(K,J,NKP,ISP)
                  ENDDO
                  U(I,J)=CTMP
               ENDDO
            ENDDO bands

            CALL SETWDES(WDES,WDES1,NK)
            WA=ELEMENTS(W,WDES1,ISP)

            ! redistribute over plane wave coefficients
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF

            ! build linear combinations of the wave functions in
            ! accordance with the effective rotation matrix U
            CALL LINCOM('F',WA%CW_RED,WA%CPROJ_RED,U(1,1), &
           &     WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), &
           &     WDES1%NPL_RED,WDES1%NPRO_RED,WDES%NRPLWV_RED,WDES%NPROD_RED,WDES%NB_TOT, &
           &     WA%CW_RED,WA%CPROJ_RED)

            ! redistribute back over bands
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF
         ENDDO kpoints
      ENDDO spin

      DEALLOCATE(U)

      RETURN
      END SUBROUTINE MLWF_ROTATE_ORBITALS_FULLK


!******************** SUBROUTINE MLWF_ROTATE_ORBITALS_NOSYMM ***********
!
!***********************************************************************
      SUBROUTINE MLWF_ROTATE_ORBITALS_NOSYMM(WDES,W)
      USE dfast
      USE wave_high
      USE full_kpoints
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      ! local variables
      TYPE(wavedes1) WDES1
      TYPE(wavefuna) WA
      INTEGER I,J,K,IP,IPP
      INTEGER ISP,NK,NKP,NW,NB,NBEXCL,NBwin
      COMPLEX(q) CTMP
      GDEF, ALLOCATABLE :: U(:,:)

      CALL CHECK_FULL_KPOINTS

      IF ((.NOT.ALLOCATED(U_matrix)).OR.(.NOT.ALLOCATED(U_matrix_opt))) THEN
         WRITE(*,*) 'MLWF_ROTATE_ORBITALS_NOSYMM: ERROR: rotation matrices not available'
         STOP
      ENDIF

      NW=SIZE(U_matrix,1)
      NB=SIZE(U_matrix_opt,1)
      IF (SIZE(lexclude_band)/=WDES%NB_TOT) THEN
         WRITE(*,*) 'MLWF_ROTATE_ORBITALS_NOSYMM: ERROR: lexclude_band array is corrupt:', &
        &   SIZE(lexclude_band),WDES%NB_TOT
         STOP
      ENDIF
      NBEXCL=0
      DO I=1,SIZE(lexclude_band)
         IF (lexclude_band(I)) NBEXCL=NBEXCL+1
      ENDDO
      IF ((NB+NBEXCL)/=WDES%NB_TOT) THEN
         WRITE(*,*) 'MLWF_ROTATE_ORBITALS_NOSYMM: ERROR: inconsistent number of bands:', &
        &   NB+NBEXCL,WDES%NB_TOT
      ENDIF

      ALLOCATE(U(WDES%NB_TOT,WDES%NB_TOT))

!      WRITE(*,*) KPOINTS_FULL%NKPTS
!      WRITE(*,'(3F14.7)') KPOINTS_FULL%VKPT

      spin: DO ISP=1,WDES%ISPIN
         kpoints: DO NK=1,WDES%NKPTS
            ! find the corresponding entry in KPOINTS_FULL, that contains
            ! the set of k-points used in the computation of the rotation matrices
            DO NKP=1,KPOINTS_FULL%NKPTS
               IF (LIDENTICAL_KPOINT(WDES%VKPT(:,NK),KPOINTS_FULL%VKPT(:,NKP))) exit
            ENDDO
            IF (NKP>KPOINTS_FULL%NKPTS) THEN
               WRITE(*,*) 'MLWF_ROTATE_ORBITALS: ERROR: no matching k-point found in KPOINTS_FULL',NK
               STOP
            ENDIF
            ! setup the effective rotation matrix U at the present k-point
            U=0._q; IP=0; IPP=0
            bands: DO I=1,WDES%NB_TOT
               IF (lexclude_band(I)) CYCLE bands
               IP=IP+1
               IF (.NOT.lwindow(IP,NKP,ISP)) CYCLE bands
               IPP=IPP+1
               DO J=1,NW
                  CTMP=0._q
                  DO K=1,NW
                     CTMP=CTMP+U_matrix_opt(IPP,K,NKP,ISP)*U_matrix(K,J,NKP,ISP)
                  ENDDO
                  U(I,J)=CTMP
               ENDDO
            ENDDO bands

            CALL SETWDES(WDES,WDES1,NK)
            WA=ELEMENTS(W,WDES1,ISP)

            ! redistribute over plane wave coefficients
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF

!           CALL ROTATE_CDER_NOINTER(WDES, U, SIZE(U,1), NK, ISP)

            ! build linear combinations of the wave functions in
            ! accordance with the effective rotation matrix U
            CALL LINCOM('F',WA%CW_RED,WA%CPROJ_RED,U(1,1), &
           &     WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), &
           &     WDES1%NPL_RED,WDES1%NPRO_RED,WDES%NRPLWV_RED,WDES%NPROD_RED,WDES%NB_TOT, &
           &     WA%CW_RED,WA%CPROJ_RED)

            ! redistribute back over bands
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF
         ENDDO kpoints
      ENDDO spin

      DEALLOCATE(U)

      RETURN
      END SUBROUTINE MLWF_ROTATE_ORBITALS_NOSYMM


      SUBROUTINE MLWF_ROTATE_ORBITALS_NOSYMM_U(WDES,W,U)
      USE wave_high
      USE full_kpoints
      USE dfast
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      GDEF :: U(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)
      ! local variables
      TYPE(wavedes1) WDES1
      TYPE(wavefuna) WA
      INTEGER I,J,K,IP,IPP
      INTEGER ISP,NK,NKP,NW,NB,NBEXCL,NBwin
      GDEF :: UU(WDES%NB_TOT,WDES%NB_TOT)

      spin: DO ISP=1,WDES%ISPIN
         kpoints: DO NK=1,WDES%NKPTS
            UU=U(:,:,NK,ISP)

            CALL SETWDES(WDES,WDES1,NK)
            WA=ELEMENTS(W,WDES1,ISP)

            ! redistribute over plane wave coefficients
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF

            CALL ROTATE_CDER_NOINTER(WDES, UU, SIZE(UU,1), NK, ISP)

            ! build linear combinations of the wave functions in
            ! accordance with the effective rotation matrix U
            CALL LINCOM('F',WA%CW_RED,WA%CPROJ_RED,UU, &
           &     WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), &
           &     WDES1%NPL_RED,WDES1%NPRO_RED,WDES%NRPLWV_RED,WDES%NPROD_RED,WDES%NB_TOT, &
           &     WA%CW_RED,WA%CPROJ_RED)

            ! redistribute back over bands
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF
         ENDDO kpoints
      ENDDO spin

      RETURN
    END SUBROUTINE MLWF_ROTATE_ORBITALS_NOSYMM_U


!******************** SUBROUTINE MLWF_DPSI_DK **************************
!
!***********************************************************************
      SUBROUTINE MLWF_DPSI_DK(W,KPOINTS,LATT_CUR,RPHI)
      USE constant
      USE pseudo
      USE lattice
      USE mkpoints
      USE full_kpoints
      USE wave_high
      USE twoelectron4o
      IMPLICIT NONE
      TYPE(wavespin) W
      TYPE(latt) LATT_CUR
      TYPE(kpoints_struct) KPOINTS
      COMPLEX(qs) RPHI(W%WDES%NRPLWV,W%WDES%NBANDS,KPOINTS%NKPTS,W%WDES%ISPIN,3)
      ! local variables
      TYPE (wavedes1) WDESK,WDESIK
      TYPE (wavefun1) WK
      INTEGER ISP,IK,NK,IB,ISPINOR,IBG
      INTEGER IDIR,J,IDELTA,IORDER,ISGN
      REAL(q) K0(3),K(3),DK(3),KP(3)
      COMPLEX(q), ALLOCATABLE :: CPHASE(:)
      LOGICAL LPHASE

      IORDER=4

      DK(1)=1._q/REAL(KPOINTS_FULL%NKPX,KIND=q)
      DK(2)=1._q/REAL(KPOINTS_FULL%NKPY,KIND=q)
      DK(3)=1._q/REAL(KPOINTS_FULL%NKPZ,KIND=q)

      RPHI=0
      spin: DO ISP=1,W%WDES%ISPIN
      kpoint: DO IK=1,KPOINTS%NKPTS
         CALL SETWDES(W%WDES,WDESIK,IK)
         K0(:)=KPOINTS%VKPT(:,IK)

         dir: DO IDIR=1,3
         delta: DO IDELTA=1,IORDER
         sgn: DO ISGN=-1,1,2

            K(:)=K0(:)
            K(IDIR)=K0(IDIR)+ISGN*IDELTA*DK(IDIR)
            NK=KPOINT_IN_FULL_GRID(K,KPOINTS_FULL)
            KP(:)=KPOINTS_FULL%VKPT(:,NK)

            CALL SETWDES(W%WDES,WDESK,NK)
            CALL NEWWAV(WK,WDESK,.TRUE.)

            ALLOCATE(CPHASE(WDESK%GRID%MPLWV))
            CALL SETPHASE( KP(:)-K(:),WDESK%GRID,CPHASE,LPHASE)

!           bands: DO IB=1,MLWF_num_wann
            bands: DO IB=1,W%WDES%NBANDS
               ! global band index
#ifdef MPI
               IBG=(IB-1)*W%WDES%NB_PAR+W%WDES%COMM_INTER%NODE_ME
               IF(IBG>MLWF_num_wann) EXIT
#endif

               CALL W1_COPY(ELEMENT(W,WDESK,IB,ISP),WK)
               CALL FFTWAV_W1(WK)

               IF (LPHASE) CALL APPLY_PHASE(WDESK%GRID,CPHASE,WK%CR,WK%CR)

               WK%CR=WK%CR/WDESK%GRID%NPLWV
               DO ISPINOR=0,WDESK%NRSPINORS-1
                  CALL FFTEXT(WDESIK%NGVECTOR,WDESIK%NINDPW(1),WK%CR(1+ISPINOR*WDESK%GRID%MPLWV),WK%CW(1+ISPINOR*WDESIK%NGVECTOR),WDESK%GRID,.FALSE.)
               ENDDO

               RPHI(1:SIZE(WK%CW),IB,IK,ISP,IDIR)=RPHI(1:SIZE(WK%CW),IB,IK,ISP,IDIR)+CMPLX(ISGN*WK%CW(1:SIZE(WK%CW))*FAC(IDELTA,IORDER)/DK(IDIR)/TPI/2._q,KIND=qs)
!               DO J=1,3
!                  RPHI(1:SIZE(WK%CW),IB,IK,ISP,J)=RPHI(1:SIZE(WK%CW),IB,IK,ISP,J)+CMPLX(ISGN*LATT_CUR%A(J,IDIR)*WK%CW(1:SIZE(WK%CW))*FAC(IDELTA,IORDER)/DK(IDIR)/TPI/2._q,KIND=qs)
!               ENDDO
            ENDDO bands

            DEALLOCATE(CPHASE)
            CALL DELWAV(WK,.TRUE.)

         ENDDO sgn
         ENDDO delta
         ENDDO dir

      ENDDO kpoint
      ENDDO spin

      RETURN
      END SUBROUTINE MLWF_DPSI_DK

!******************** SUBROUTINE WANNIER90_ORBITAL_DEFINITIONS *********
!
! VASP sets up the YLM functions in the following order
! (see for instance SETYLM in asa.F)
! YLM(:,1)     -> s
! YLM(:,2:4)   -> p:= y, z, x
! YLM(:,5:9)   -> d:= xy, yz, z2, xz, x2
! YLM(:,10:16) -> f:= y(3x2-y2), xyz, yz2, z3, xz2, z(x2-y2), x(x2-3y2)
!
! This routine provides a translation between the aforementioned and
! the orbital definitions used in wannier90
!
!***********************************************************************
      SUBROUTINE WANNIER90_ORBITAL_DEFINITIONS( &
     &   L,M,HYBRID_ORBITAL &
     &)
      IMPLICIT NONE
      INTEGER L,M
      REAL(q) HYBRID_ORBITAL(:)
      ! local variables

      ! we should be able to deal with anything up to and including L=3
      IF (SIZE(HYBRID_ORBITAL)<16) THEN
         WRITE(*,*) 'WANNIER90_ORBITAL_DEFINITIONS: ERROR: HYBRID_ORBITAL array too small', &
        &   SIZE(HYBRID_ORBITAL)
         STOP
      ENDIF

      HYBRID_ORBITAL=0

      SELECT CASE(L)
         CASE(0)
            ! s-function
            IF (M==1) THEN
               HYBRID_ORBITAL(1)=1._q
            ELSE
               WRITE(*,'(A,I2,A,I2,A)')  &
              &   'WANNIER90_ORBITAL_DEFINITIONS: ERROR: L=',L,' M=',M,' not implemented'
               STOP
            ENDIF
         CASE(1)
            ! p-functions
            IF (M==1) THEN
               HYBRID_ORBITAL(3)=1._q
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(4)=1._q
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(2)=1._q
            ELSE
               WRITE(*,'(A,I2,A,I2,A)')  &
              &   'WANNIER90_ORBITAL_DEFINITIONS: ERROR: L=',L,' M=',M,' not implemented'
               STOP
            ENDIF
         CASE(2)
            ! d-functions
            IF (M==1) THEN
               HYBRID_ORBITAL(7)=1._q
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(8)=1._q
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(6)=1._q
            ELSEIF (M==4) THEN
               HYBRID_ORBITAL(9)=1._q
            ELSEIF (M==5) THEN
               HYBRID_ORBITAL(5)=1._q
            ELSE
               WRITE(*,'(A,I2,A,I2,A)')  &
              &   'WANNIER90_ORBITAL_DEFINITIONS: ERROR: L=',L,' M=',M,' not implemented'
               STOP
            ENDIF
         CASE(3)
            ! f-functions
            IF (M==1) THEN
               HYBRID_ORBITAL(13)=1._q
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(14)=1._q
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(12)=1._q
            ELSEIF (M==4) THEN
               HYBRID_ORBITAL(15)=1._q
            ELSEIF (M==5) THEN
               HYBRID_ORBITAL(11)=1._q
            ELSEIF (M==6) THEN
               HYBRID_ORBITAL(16)=1._q
            ELSEIF (M==7) THEN
               HYBRID_ORBITAL(10)=1._q
            ELSE
               WRITE(*,'(A,I2,A,I2,A)')  &
              &   'WANNIER90_ORBITAL_DEFINITIONS: ERROR: L=',L,' M=',M,' not implemented'
               STOP
            ENDIF
         CASE(4:)
            ! nothing beyond L=3 yet
            WRITE(*,'(A,I2,A,I2,A)')  &
           &   'WANNIER90_ORBITAL_DEFINITIONS: ERROR: L=',L,' M=',M,' not implemented'
            STOP
         CASE(-1)
            ! sp-hybrids
            IF (M==1) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(2._q)
               HYBRID_ORBITAL(4)= 1._q/SQRT(2._q)
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(2._q)
               HYBRID_ORBITAL(4)=-1._q/SQRT(2._q)
            ELSE
               WRITE(*,'(A,I2,A,I2,A)')  &
              &   'WANNIER90_ORBITAL_DEFINITIONS: ERROR: L=',L,' M=',M,' not implemented'
               STOP
            ENDIF
         CASE(-2)
            ! sp2-hybrids
            IF (M==1) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(3._q)
               HYBRID_ORBITAL(4)=-1._q/SQRT(6._q)
               HYBRID_ORBITAL(2)= 1._q/SQRT(2._q)
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(3._q)
               HYBRID_ORBITAL(4)=-1._q/SQRT(6._q)
               HYBRID_ORBITAL(2)=-1._q/SQRT(2._q)
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(3._q)
               HYBRID_ORBITAL(4)= 2._q/SQRT(6._q)
            ELSE
               WRITE(*,'(A,I2,A,I2,A)')  &
              &   'WANNIER90_ORBITAL_DEFINITIONS: ERROR: L=',L,' M=',M,' not implemented'
               STOP
            ENDIF
         CASE(-3)
            ! sp3-hybrids
            IF (M==1) THEN
               HYBRID_ORBITAL(1)= 1._q/2._q
               HYBRID_ORBITAL(4)= 1._q/2._q
               HYBRID_ORBITAL(2)= 1._q/2._q
               HYBRID_ORBITAL(3)= 1._q/2._q
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(1)= 1._q/2._q
               HYBRID_ORBITAL(4)= 1._q/2._q
               HYBRID_ORBITAL(2)=-1._q/2._q
               HYBRID_ORBITAL(3)=-1._q/2._q
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(1)= 1._q/2._q
               HYBRID_ORBITAL(4)=-1._q/2._q
               HYBRID_ORBITAL(2)= 1._q/2._q
               HYBRID_ORBITAL(3)=-1._q/2._q
            ELSEIF (M==4) THEN
               HYBRID_ORBITAL(1)= 1._q/2._q
               HYBRID_ORBITAL(4)=-1._q/2._q
               HYBRID_ORBITAL(2)=-1._q/2._q
               HYBRID_ORBITAL(3)= 1._q/2._q
            ELSE
               WRITE(*,'(A,I2,A,I2,A)')  &
              &   'WANNIER90_ORBITAL_DEFINITIONS: ERROR: L=',L,' M=',M,' not implemented'
               STOP
            ENDIF
         CASE(-4)
            ! spd3d-hybrids
            IF (M==1) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(3._q)
               HYBRID_ORBITAL(4)=-1._q/SQRT(6._q)
               HYBRID_ORBITAL(2)= 1._q/SQRT(2._q)
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(3._q)
               HYBRID_ORBITAL(4)=-1._q/SQRT(6._q)
               HYBRID_ORBITAL(2)=-1._q/SQRT(2._q)
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(3._q)
               HYBRID_ORBITAL(4)= 2._q/SQRT(6._q)
            ELSEIF (M==4) THEN
               HYBRID_ORBITAL(3)= 1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)= 1._q/SQRT(2._q)
            ELSEIF (M==5) THEN
               HYBRID_ORBITAL(3)=-1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)= 1._q/SQRT(2._q)
            ELSE
               WRITE(*,'(A,I2,A,I2,A)')  &
              &   'WANNIER90_ORBITAL_DEFINITIONS: ERROR: L=',L,' M=',M,' not implemented'
               STOP
            ENDIF
         CASE(-5)
            ! sp3d2-hybrids
            IF (M==1) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(6._q)
               HYBRID_ORBITAL(4)=-1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)=-1._q/SQRT(12._q)
               HYBRID_ORBITAL(9)= 1._q/2._q
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(6._q)
               HYBRID_ORBITAL(4)= 1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)=-1._q/SQRT(12._q)
               HYBRID_ORBITAL(9)= 1._q/2._q
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(6._q)
               HYBRID_ORBITAL(2)=-1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)=-1._q/SQRT(12._q)
               HYBRID_ORBITAL(9)=-1._q/2._q
            ELSEIF (M==4) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(6._q)
               HYBRID_ORBITAL(2)= 1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)=-1._q/SQRT(12._q)
               HYBRID_ORBITAL(9)=-1._q/2._q
            ELSEIF (M==5) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(6._q)
               HYBRID_ORBITAL(3)=-1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)= 1._q/SQRT(3._q)
            ELSEIF (M==6) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(6._q)
               HYBRID_ORBITAL(3)= 1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)= 1._q/SQRT(3._q)
            ELSE
               WRITE(*,'(A,I2,A,I2,A)')  &
              &   'WANNIER90_ORBITAL_DEFINITIONS: ERROR: L=',L,' M=',M,' not implemented'
               STOP
            ENDIF
         CASE(:-6)
            WRITE(*,'(A,I2,A,I2,A)')  &
           &   'WANNIER90_ORBITAL_DEFINITIONS: ERROR: L=',L,' M=',M,' not implemented'
            STOP

      END SELECT

      RETURN
      END SUBROUTINE WANNIER90_ORBITAL_DEFINITIONS


!******************** SUBROUTINE SETROTYLM *****************************
!
!***********************************************************************
      SUBROUTINE SETROTYLM( &
     & XIN,ZIN,LMAX,ROTYLM &
     &)
      USE asa
      USE constant
      IMPLICIT NONE
      INTEGER LMAX
      REAL(q) XIN(3),ZIN(3)
      REAL(q) ROTYLM(:,:)
      ! local variables
      INTEGER LMMAX
      REAL(q) XN,ZN,XDOTZ
      REAL(q) X(3),Y(3),Z(3),U(3,3)
      REAL(q), PARAMETER :: TINY=1.0E-6_q

      INTEGER I,J,IFAIL,LM,LMP
      INTEGER PHPTS,THPTS,NPTS,NP
      REAL(q) SCALE,DELTAPHI,SIM_FAKT
      REAL(q), ALLOCATABLE :: RADPTS(:,:),XYZPTS(:,:),UXYZPTS(:,:)
      REAL(q), ALLOCATABLE :: YLM(:,:),YLMP(:,:)
      REAL(q), ALLOCATABLE :: WEIGHT(:),ABSCIS(:)
      EXTERNAL GAUSSI2

      X=XIN; Z=ZIN

      LMMAX=(LMAX+1)**2

      ! check size of ROTYLM
      IF (SIZE(ROTYLM,1)<LMMAX.OR.SIZE(ROTYLM,2)<LMMAX) THEN
         WRITE(*,*) 'SETROTYLM: ERROR: ROTYLM too small:',LMMAX,SIZE(ROTYLM,1),SIZE(ROTYLM,2)
         STOP
      ENDIF
      ! check size of X and Z
      XN=SQRT(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
      ZN=SQRT(Z(1)*Z(1)+Z(2)*Z(2)+Z(3)*Z(3))
      IF (XN<TINY.OR.ZN<TINY) THEN
         WRITE(*,*) 'SETROTYLM: ERROR: |X| or |Z| very small:',XN,ZN
         STOP
      ENDIF
      ! and normalize to 1
      X=X/XN; Z=Z/ZN
      ! check orthogonality of X and Z
      XDOTZ=X(1)*Z(1)+X(2)*Z(2)+X(3)*Z(3)
      IF (ABS(XDOTZ)>TINY) THEN
         WRITE(*,*) 'SETROTYLM: ERROR: X and Z are not orthogonal (enough):',ABS(XDOTZ)
         STOP
      ENDIF
      ! y=Z \times X
      Y(1)=(Z(2)*X(3)-X(2)*Z(3))
      Y(2)=(Z(3)*X(1)-X(3)*Z(1))
      Y(3)=(Z(1)*X(2)-X(1)*Z(2))
      ! transformation matrix
      U(1,:)=X(:); U(2,:)=Y(:); U(3,:)=Z(:)
! test
!     WRITE(*,*) 'transformation matrix'
!     WRITE(*,'(3F10.5)') U(1,:)
!     WRITE(*,'(3F10.5)') U(2,:)
!     WRITE(*,'(3F10.5)') U(3,:)
! test
      SCALE=2*SQRT(PI) ! 1/Y00

!========================================================================
! number of theta and phi pivot points to perform angular integration
! since Exc=f(a*Yllmax,m) we need more pivot points than theoretically
! needed to integrate Yllmax,m.
! the factor 2 is the minium, 3 is more accurate
!========================================================================
      PHPTS=3*(LMAX+1)
      THPTS=3*FLOOR(REAL(LMAX/2+1,KIND=q))
      NPTS=PHPTS*THPTS
      DELTAPHI=REAL(2_q*PI/PHPTS,KIND=q)
      ! allocate arrays
      ALLOCATE(XYZPTS(NPTS,3),UXYZPTS(NPTS,3),RADPTS(NPTS,2),WEIGHT(THPTS),ABSCIS(THPTS), &
     &     YLM(NPTS,LMMAX),YLMP(NPTS,LMMAX))

      ! set phi positions, equally spaced
      RADPTS=0; WEIGHT=0; ABSCIS=0
      DO I=1,PHPTS
         DO J=1,THPTS
            RADPTS((J-1)*PHPTS+I,2)=(I-1)*DELTAPHI
         ENDDO
      ENDDO
     ! get theta positions (actually get cos(theta)) (Gauss integration)
      CALL GAUSSI(GAUSSI2,-1._q,1._q,0,THPTS,WEIGHT,ABSCIS,IFAIL)
      DO I=1,THPTS
         RADPTS((I-1)*PHPTS+1:I*PHPTS,1)=ABSCIS(I)
      ENDDO
      ! convert radial to cartesian coordinates
      DO I=1,NPTS
         XYZPTS(I,1)=COS(RADPTS(I,2))*SQRT(1_q-RADPTS(I,1)**2_q) ! x
         XYZPTS(I,2)=SIN(RADPTS(I,2))*SQRT(1_q-RADPTS(I,1)**2_q) ! y
         XYZPTS(I,3)=RADPTS(I,1)                                 ! z
      ENDDO

      YLM=0
      CALL SETYLM(LMAX,NPTS,YLM,XYZPTS(:,1),XYZPTS(:,2),XYZPTS(:,3))

      ! rotate XYZPTS, using accordance with U
      DO I=1,NPTS
         UXYZPTS(I,1)=U(1,1)*XYZPTS(I,1)+U(1,2)*XYZPTS(I,2)+U(1,3)*XYZPTS(I,3)
         UXYZPTS(I,2)=U(2,1)*XYZPTS(I,1)+U(2,2)*XYZPTS(I,2)+U(2,3)*XYZPTS(I,3)
         UXYZPTS(I,3)=U(3,1)*XYZPTS(I,1)+U(3,2)*XYZPTS(I,2)+U(3,3)*XYZPTS(I,3)
      ENDDO

      YLMP=0
      CALL SETYLM(LMAX,NPTS,YLMP,UXYZPTS(:,1),UXYZPTS(:,2),UXYZPTS(:,3))

      ROTYLM=0
      ! loop over all points in the angular grid
      points: DO NP=1,NPTS

         ! weight of this points
         SIM_FAKT=DELTAPHI*WEIGHT((INT((NP-1)/PHPTS)+1))

         DO LM=1,LMMAX
         DO LMP=1,LMMAX
            ROTYLM(LMP,LM)=ROTYLM(LMP,LM)+YLMP(NP,LMP)*YLM(NP,LM)*SIM_FAKT
         ENDDO
         ENDDO

      ENDDO points

      DEALLOCATE(XYZPTS,UXYZPTS,RADPTS,WEIGHT,ABSCIS,YLM,YLMP)
      RETURN
      END SUBROUTINE SETROTYLM


!******************** SUBROUTINE SETRGRID ******************************
!
!***********************************************************************
      SUBROUTINE SETRGRID( &
     & RSTART,REND,H,R &
     &)
      USE radial
      IMPLICIT NONE
      TYPE(rgrid) R
      REAL(q) RSTART,REND,H
      ! local variables
      INTEGER I

      I=0
      DO
        I=I+1
        IF (RSTART*EXP(H*(I-1))>=REND) EXIT
      ENDDO
      ! make sure I is uneven
      I=I+(MODULO(I,2)+1)
      R%NMAX=I
      R%RSTART=RSTART; R%H=H; R%D=EXP(H)
      IF (ASSOCIATED(R%R)) THEN
         DEALLOCATE(R%R); NULLIFY(R%R)
      ENDIF
      ALLOCATE(R%R(I))
      DO I=1,R%NMAX
         R%R(I)=R%RSTART*EXP(H*(I-1))
      ENDDO
      R%REND=R%R(R%NMAX)

      IF (ASSOCIATED(R%SI)) THEN
         DEALLOCATE(R%SI); NULLIFY(R%SI)
      ENDIF
      CALL SET_SIMP(R)

      RETURN
      END SUBROUTINE SETRGRID


!******************** SUBROUTINE RADIAL_FUNCTION ***********************
!
!***********************************************************************
      SUBROUTINE RADIAL_FUNCTION( &
     & ITYP,R,ZA,FR &
     &)
      USE radial
      IMPLICIT NONE
      TYPE(rgrid) R
      INTEGER ITYP
      REAL(q) ZA
      REAL(q) FR(:)

      IF (SIZE(FR)<R%NMAX) THEN
         WRITE(*,*) 'RADIAL_FUNCTION: ERROR: FR too small:',SIZE(FR),R%NMAX
      ENDIF
      IF (ITYP==1) THEN
         FR(:)=2._q*ZA**(3._q/2._q)*EXP(-ZA*R%R(:))
      ELSEIF (ITYP==2) THEN
         FR(:)=1._q/(2._q*SQRT(2._q))*ZA**(3._q/2._q)*(2._q-ZA*R%R(:))*EXP(-ZA*R%R(:)/2._q)
      ELSEIF (ITYP==3) THEN
         FR(:)=SQRT(4._q/27._q)*ZA**(3._q/2._q)* &
        &   (1._q-2._q/3._q*ZA*R%R(:)+2._q/27._q*ZA*ZA*R%R(:)*R%R(:))*EXP(-ZA*R%R(:)/3._q)
      ELSE
         WRITE(*,*) 'RADIAL_FUNCTION: TYPE does not exist',ITYP
      ENDIF

      RETURN
      END SUBROUTINE RADIAL_FUNCTION


!******************** SUBROUTINE CALC_OVERLAP_GN ***********************
!
!***********************************************************************
      SUBROUTINE CALC_OVERLAP_GN_band( &
     &   nbmin,nbmax,L,FG,POS,W,K,ISP,ISPINOR,P,CQIJ,LATT_CUR,T_INFO,S &
     &)
      USE pead
      USE poscar
      USE pseudo
      USE lattice
      USE full_kpoints
      USE wave_high
      USE nonl_high
      IMPLICIT NONE
      TYPE(wavespin) W
      TYPE(potcar) P(:)
      TYPE(latt) LATT_CUR
      TYPE(type_info) T_INFO
      INTEGER L
      INTEGER ISP
      INTEGER ISPINOR
      INTEGER nbmin
      INTEGER nbmax
      REAL(q) K(3)
      REAL(q) FG(:,:,:)
      REAL(q) POS(3)
      OVERLAP CQIJ(:,:,:,:)
      GDEF S(W%WDES%NB_TOT,(L+1)**2)
      ! local variables
      TYPE(wavespin) WP
      TYPE(wavefuna) WK,WRYlm
      TYPE(wavedes1), TARGET :: WDESK
      TYPE(nonl_struct) NONL_S

      TYPE(rotation_handle), POINTER :: ROT_HANDLE

      GDEF C
      REAL(q) WSCAL
      INTEGER NK,NB,N,NYLM

      WP=W
      WP%WDES=>WDES_FULL_PEAD

      CALL CHECK_FULL_KPOINTS

      NULLIFY(ROT_HANDLE)

      NYLM=(L+1)**2

      ! search for kpoint k in BZ
      NK=KPOINT_IN_FULL_GRID(K,KPOINTS_FULL)
      CALL SETWDES(WP%WDES,WDESK,NK)
      IF (NK==KPOINTS_FULL%NEQUIV(NK)) THEN
         ! k is a kpoint in the IBZ
         WK=ELEMENTS(WP,WDESK,ISP)
      ELSE
         ! k is not a kpoint in the IBZ
         CALL NEWWAVA(WK,WDESK,WDESK%NBANDS)
         CALL PEAD_WA_ROTATE(WP,P,LATT_CUR,ISP,WK)
      ENDIF

      CALL NONL_ALLOC(NONL_S,T_INFO,P,WP%WDES,.FALSE.)
      CALL SPHER(WP%WDES%GRID,NONL_S,P,WP%WDES,LATT_CUR,1,NK)
      CALL PHASE(WP%WDES,NONL_S,NK)

      CALL NEWWAVA(WRYlm,WDESK,NYLM)
      WRYlm%CW=0
      CALL CONSTRUCT_FUNCTION_RYlm(L,FG,LATT_CUR,POS,NONL_S,ISPINOR,WRYlm)
      ! and normalize the functions WRYlm
      DO N=1,NYLM
         CALL CNORMN(ELEMENT(WRYlm,N),CQIJ,1,WSCAL)
      ENDDO

      ! calculate overlap between Wk and WRYlm: < w_{m,k1} | S | RYlm >
      !WRITE(*,*) WDESK%NB_LOW,WDESK%NB_PAR,WP%WDES%NBANDS
      !STOP
      S=0
      DO NB=1,WP%WDES%NBANDS
         IF (WDESK%NB_LOW+WDESK%NB_PAR*(NB-1).LT.nbmin.OR.WDESK%NB_LOW+WDESK%NB_PAR*(NB-1).GT.nbmax) CYCLE
         DO N=1,NYLM
            C=W1_DOT(ELEMENT(WK,NB),ELEMENT(WRYlm,N),CQIJ)
            S(WDESK%NB_LOW+WDESK%NB_PAR*(NB-1),N)=C
!           WRITE(*,*) WDESK%NB_LOW+WDESK%NB_PAR*(NB-1),N,C
         ENDDO
      ENDDO

      CALLMPI( M_sum_g(WDESK%COMM_INTER,S(1,1),WDESK%NB_TOT*NYLM) )

      ! some deallocation to be done
      CALL DEALLOCATE_ROT_HANDLE(ROT_HANDLE)

      CALL DELWAVA(WRYlm)
      IF (NK/=KPOINTS_FULL%NEQUIV(NK)) CALL DELWAVA(WK)

      CALL NONL_DEALLOC(NONL_S)

      RETURN
      END SUBROUTINE CALC_OVERLAP_GN_band
      SUBROUTINE CALC_OVERLAP_GN( &
     &   L,FG,POS,W,K,ISP,ISPINOR,P,CQIJ,LATT_CUR,T_INFO,S &
     &)
      USE pead
      USE poscar
      USE pseudo
      USE lattice
      USE full_kpoints
      USE wave_high
      USE nonl_high
      IMPLICIT NONE
      TYPE(wavespin) W
      TYPE(potcar) P(:)
      TYPE(latt) LATT_CUR
      TYPE(type_info) T_INFO
      INTEGER L
      INTEGER ISP
      INTEGER ISPINOR
      REAL(q) K(3)
      REAL(q) FG(:,:,:)
      REAL(q) POS(3)
      OVERLAP CQIJ(:,:,:,:)
      GDEF S(W%WDES%NB_TOT,(L+1)**2)
      ! local variables
      TYPE(wavespin) WP
      TYPE(wavefuna) WK,WRYlm
      TYPE(wavedes1), TARGET :: WDESK
      TYPE(nonl_struct) NONL_S

      TYPE(rotation_handle), POINTER :: ROT_HANDLE

      GDEF C
      REAL(q) WSCAL
      INTEGER NK,NB,N,NYLM

      WP=W
      WP%WDES=>WDES_FULL_PEAD

      CALL CHECK_FULL_KPOINTS

      NULLIFY(ROT_HANDLE)

      NYLM=(L+1)**2

      ! search for kpoint k in BZ
      NK=KPOINT_IN_FULL_GRID(K,KPOINTS_FULL)
      CALL SETWDES(WP%WDES,WDESK,NK)
      IF (NK==KPOINTS_FULL%NEQUIV(NK)) THEN
         ! k is a kpoint in the IBZ
         WK=ELEMENTS(WP,WDESK,ISP)
      ELSE
         ! k is not a kpoint in the IBZ
         CALL NEWWAVA(WK,WDESK,WDESK%NBANDS)
         CALL PEAD_WA_ROTATE(WP,P,LATT_CUR,ISP,WK)
      ENDIF

      CALL NONL_ALLOC(NONL_S,T_INFO,P,WP%WDES,.FALSE.)
      CALL SPHER(WP%WDES%GRID,NONL_S,P,WP%WDES,LATT_CUR,1,NK)
      CALL PHASE(WP%WDES,NONL_S,NK)

      CALL NEWWAVA(WRYlm,WDESK,NYLM)
      WRYlm%CW=0
      CALL CONSTRUCT_FUNCTION_RYlm(L,FG,LATT_CUR,POS,NONL_S,ISPINOR,WRYlm)
      ! and normalize the functions WRYlm
      DO N=1,NYLM
         CALL CNORMN(ELEMENT(WRYlm,N),CQIJ,1,WSCAL)
      ENDDO

      ! calculate overlap between Wk and WRYlm: < w_{m,k1} | S | RYlm >
      S=0
      DO NB=1,WP%WDES%NBANDS
         DO N=1,NYLM
            C=W1_DOT(ELEMENT(WK,NB),ELEMENT(WRYlm,N),CQIJ)
            S(WDESK%NB_LOW+WDESK%NB_PAR*(NB-1),N)=C
!           WRITE(*,*) WDESK%NB_LOW+WDESK%NB_PAR*(NB-1),N,C
         ENDDO
      ENDDO

      CALLMPI( M_sum_g(WDESK%COMM_INTER,S(1,1),WDESK%NB_TOT*NYLM) )

      ! some deallocation to be done
      CALL DEALLOCATE_ROT_HANDLE(ROT_HANDLE)

      CALL DELWAVA(WRYlm)
      IF (NK/=KPOINTS_FULL%NEQUIV(NK)) CALL DELWAVA(WK)

      CALL NONL_DEALLOC(NONL_S)

      RETURN
      END SUBROUTINE CALC_OVERLAP_GN


!******************** SUBROUTINE CONSTRUCT_FUNCTION_RYlm ***************
!
!***********************************************************************
      SUBROUTINE CONSTRUCT_FUNCTION_RYlm( &
     &   LMAX,F,LATT_CUR,POS,NONL_S,ISPINOR,WRYLM &
     &)
      USE ini
      USE asa
      USE pead
      USE lattice
      USE constant
      USE wave_high
      USE nonl_high
      IMPLICIT NONE
      TYPE(latt) LATT_CUR
      TYPE(wavefuna) WRYLM
      TYPE (nonl_struct) NONL_S
      REAL(q) F(:,:,:)
      REAL(q) POS(3)
      INTEGER LMAX
      INTEGER ISPINOR
      ! local variables
      TYPE(wavefun1) W1
      INTEGER N1,N2,N3,IND,NPL,LMMAX,L,M,LM,IG
      REAL(q) G1,G2,G3,GKX,GKY,GKZ,FACTM,FAKT,FDER
      REAL(q), ALLOCATABLE :: XS(:),YS(:),ZS(:),YLM(:,:)
      REAL(q), ALLOCATABLE :: G(:),FG(:)
      COMPLEX(q) CSET,CGDR
      COMPLEX(q), ALLOCATABLE :: CFAKTX(:)

      LMMAX=(LMAX+1)**2

      NPL=WRYLM%WDES1%NGVECTOR

      IF (ISPINOR/=1.AND.(.NOT.WRYLM%WDES1%LNONCOLLINEAR)) THEN
         WRITE(*,*) 'CONSTRUCT_FUNCTION_RYlm: ERROR: ISPINOR=',ISPINOR,' but LNONCOLLINEAR=.FALSE.'
         STOP
      ENDIF

      ALLOCATE(G(NPL),FG(NPL))
      ALLOCATE(XS(NPL),YS(NPL),ZS(NPL),CFAKTX(NPL))

      ! loop over all G-vectors in the basis at this k-point
      DO IND=1,WRYLM%WDES1%NGVECTOR
         N1=MOD(WRYLM%WDES1%IGX(IND)+WRYLM%WDES1%GRID%NGX,WRYLM%WDES1%GRID%NGX)+1
         N2=MOD(WRYLM%WDES1%IGY(IND)+WRYLM%WDES1%GRID%NGY,WRYLM%WDES1%GRID%NGY)+1
         N3=MOD(WRYLM%WDES1%IGZ(IND)+WRYLM%WDES1%GRID%NGZ,WRYLM%WDES1%GRID%NGZ)+1

         G1=(WRYLM%WDES1%GRID%LPCTX(N1)+WRYLM%WDES1%VKPT(1))
         G2=(WRYLM%WDES1%GRID%LPCTY(N2)+WRYLM%WDES1%VKPT(2))
         G3=(WRYLM%WDES1%GRID%LPCTZ(N3)+WRYLM%WDES1%VKPT(3))

         FACTM=1._q
         IF (WRYLM%WDES1%LGAMMA .AND. (N1/=1 .OR. N2/=1 .OR. N3/=1)) FACTM=SQRT(2._q)

         GKX=(G1*LATT_CUR%B(1,1)+G2*LATT_CUR%B(1,2)+G3*LATT_CUR%B(1,3))*TPI
         GKY=(G1*LATT_CUR%B(2,1)+G2*LATT_CUR%B(2,2)+G3*LATT_CUR%B(2,3))*TPI
         GKZ=(G1*LATT_CUR%B(3,1)+G2*LATT_CUR%B(3,2)+G3*LATT_CUR%B(3,3))*TPI

         G(IND)=MAX(SQRT(GKX*GKX+GKY*GKY+GKZ*GKZ),1E-10_q)

         ! phase factor e^{-i(k+G)R} where R is the origin
         ! of the localized function
         CGDR=CITPI*(G1*POS(1)+G2*POS(2)+G3*POS(3))
!        CGDR=0
         CFAKTX(IND)=FACTM*EXP(-CGDR)

         XS(IND)  =GKX/G(IND)
         YS(IND)  =GKY/G(IND)
         ZS(IND)  =GKZ/G(IND)
      ENDDO

      ALLOCATE(YLM(NPL,LMMAX))
      ! get me all the Y_lm up to and including l=LMAX
      CALL SETYLM(LMAX,NPL,YLM,XS,YS,ZS)

      ! Setup the plane wave part of the desired function
      FAKT= 1/SQRT(LATT_CUR%OMEGA)
      CSET=CMPLX(0._q,-1._q,q)

!     WRYLM%CW=0
      LM=1
      DO L=0,LMAX
         ! get me the Bessel transform of the radial function
         DO IG=1,SIZE(G)
            CALL SPLVAL(G(IG),FG(IG),FDER,F(:,:,L+1),SIZE(F,1),SIZE(F,1))
         ENDDO
         DO IND=1,NPL
            DO M=1,2*L+1
               WRYLM%CW(IND+(ISPINOR-1)*NPL,LM+M-1)= &
              &   FAKT*(CSET**L)*CFAKTX(IND)*FG(IND)*YLM(IND,LM+M-1)
            ENDDO
         ENDDO
         LM=LM+2*L+1
      ENDDO

      ! and get the projections of RYlm onto the PAW projectors
      WRYLM%CPROJ=0
      DO LM=1,LMMAX
         W1=ELEMENT(WRYLM,LM)
         CALL PROJ1(NONL_S,WRYLM%WDES1,W1)
      ENDDO

      DEALLOCATE(G,FG,XS,YS,ZS,CFAKTX,YLM)

      RETURN
      END SUBROUTINE CONSTRUCT_FUNCTION_RYlm


!******************** SUBROUTINE WRITE_WAVE_FUNCTIONS ******************
!
!***********************************************************************
      SUBROUTINE WRITE_WAVE_FUNCTIONS(W,K,ISP,EXCLUDE_BAND,P,LATT_CUR,IU)
      USE pead
      USE pseudo
      USE lattice
      USE wave_high
      USE full_kpoints
      IMPLICIT NONE
      TYPE(wavespin) W
      TYPE(potcar) P(:)
      TYPE(latt) LATT_CUR
      INTEGER ISP,IU
      REAL(q) K(3)
      LOGICAL EXCLUDE_BAND(W%WDES%NB_TOT)
     ! local variables
      TYPE(wavespin) WP
      TYPE(wavefuna) WK
      TYPE(wavedes1), TARGET :: WDESK
      TYPE (wavefun1), ALLOCATABLE :: WSTRIP(:),WCOLLECT(:)

      TYPE(rotation_handle), POINTER :: ROT_HANDLE

      INTEGER NK,NB,ISTRIP
      INTEGER NX,NY,NZ,NC,NGX,NGY,NGZ,IND,I,NWRITTEN
      INTEGER, PARAMETER :: NSTRIP=1

      GDEF, ALLOCATABLE :: WVFN(:)

      WP=W
      WP%WDES=>WDES_FULL_PEAD
      NGX=WP%WDES%GRID%NGX; NGY=WP%WDES%GRID%NGY; NGZ=WP%WDES%GRID%NGZ

      CALL CHECK_FULL_KPOINTS

      NULLIFY(ROT_HANDLE)

      ! search for kpoint k in BZ
      NK=KPOINT_IN_FULL_GRID(K,KPOINTS_FULL)
      CALL SETWDES(WP%WDES,WDESK,NK)
      IF (NK==KPOINTS_FULL%NEQUIV(NK)) THEN
         ! k is a kpoint in the IBZ
         WK=ELEMENTS(WP,WDESK,ISP)
      ELSE
         ! k is not a kpoint in the IBZ
         CALL NEWWAVA(WK,WDESK,WDESK%NBANDS)
         CALL PEAD_WA_ROTATE(WP,P,LATT_CUR,ISP,WK)
      ENDIF

      ALLOCATE(WSTRIP(NSTRIP))

      ALLOCATE(WCOLLECT(NSTRIP*WP%WDES%NB_PAR))
      DO NB=1,NSTRIP*WP%WDES%NB_PAR
         CALL NEWWAV(WCOLLECT(NB),WDESK,.TRUE.)
      ENDDO

      ALLOCATE(WVFN(NGX*NGY*NGZ))

      NWRITTEN=0

      DO NB=1,WDESK%NBANDS,NSTRIP
         DO ISTRIP=NB,MIN(WDESK%NBANDS,NB+NSTRIP-1)
            WSTRIP(ISTRIP-NB+1)=ELEMENT(WK,ISTRIP)
         ENDDO
         CALL W1_GATHER_W1(WP,MIN(WDESK%NBANDS-NB+1,NSTRIP),WSTRIP,WCOLLECT)
         DO ISTRIP=1,MIN(WDESK%NBANDS-NB+1,NSTRIP)*WP%WDES%NB_PAR
            ! skip this band?
            IF (EXCLUDE_BAND((NB-1)*WP%WDES%NB_PAR+ISTRIP)) CYCLE
            ! if not, proceed ...
            I=0; WVFN=0
            ! might seem strange to duplicate these loops,
            ! but I don't want this conditional statement nested inside

            IF (WP%WDES%GRID%RL%NFAST==1) THEN
               DO NC=1,WP%WDES%GRID%RL%NCOL
                  NY=WP%WDES%GRID%RL%I2(NC)
                  NZ=WP%WDES%GRID%RL%I3(NC)
                  DO NX=1,WP%WDES%GRID%RL%NROW
                     IND=NX+(NY-1)*NGX+(NZ-1)*NGX*NGY
#ifdef gammareal
                     WVFN(IND)=REAL(WCOLLECT(ISTRIP)%CR(NX+(NC-1)*WP%WDES%GRID%RL%NROW),KIND=q)
#else
                     WVFN(IND)=WCOLLECT(ISTRIP)%CR(NX+(NC-1)*WP%WDES%GRID%RL%NROW)
#endif
                  ENDDO
               ENDDO
            ELSEIF (WP%WDES%GRID%RL%NFAST==3) THEN
               DO NC=1,WP%WDES%GRID%RL%NCOL
                  NX=WP%WDES%GRID%RL%I2(NC)
                  NY=WP%WDES%GRID%RL%I3(NC)
                  DO NZ=1,WP%WDES%GRID%RL%NROW
                     IND=NX+(NY-1)*NGX+(NZ-1)*NGX*NGY
#ifdef gammareal
                     WVFN(IND)=REAL(WCOLLECT(ISTRIP)%CR(NZ+(NC-1)*WP%WDES%GRID%RL%NROW),KIND=q)
#else
                     WVFN(IND)=WCOLLECT(ISTRIP)%CR(NZ+(NC-1)*WP%WDES%GRID%RL%NROW)
#endif
                  ENDDO
               ENDDO
            ELSE
               WRITE(*,'(A)') 'WRITE_WAVE_FUNCTIONS: ERROR: W1 grid not set'
               STOP
            ENDIF
            CALLMPI(M_sum_g(WP%WDES%COMM_INB,WVFN(1),NGX*NGY*NGZ))

            ! write WVFN to file
#ifdef MPI
            IF (WP%WDES%COMM%NODE_ME==WP%WDES%COMM%IONODE) THEN
#endif
               WRITE(IU) (WVFN(I),I=1,NGX*NGY*NGZ)
#ifdef MPI
            ENDIF
#endif
            NWRITTEN=NWRITTEN+1
         ENDDO
      ENDDO

      ! some deallocation to be done
      CALL DEALLOCATE_ROT_HANDLE(ROT_HANDLE)
      IF (NK/=KPOINTS_FULL%NEQUIV(NK)) CALL DELWAVA(WK)
      DO NB=1,NSTRIP*WP%WDES%NB_PAR
         CALL DELWAV(WCOLLECT(NB),.TRUE.)
      ENDDO
      DEALLOCATE(WSTRIP,WCOLLECT)
      DEALLOCATE(WVFN)

      RETURN
      END SUBROUTINE WRITE_WAVE_FUNCTIONS


!******************** SUBROUTINE OCCURS_IN_FILE ************************
!
!***********************************************************************
      SUBROUTINE OCCURS_IN_FILE(IUNIT,PATTERN,N)
      IMPLICIT NONE
      INTEGER IUNIT,N
      CHARACTER(*) PATTERN
      ! local variables
      INTEGER NOCCUR,NN
      CHARACTER(LEN=20) PATT
      CHARACTER(LEN=255) BUFLIN

      EXTERNAL NOCCUR
      REWIND(IUNIT)

      WRITE(PATT,'(A20)') PATTERN
      CALL STRIP(PATT,NN,'L')
      CALL UPPER(PATT)

      N=0
      DO
         READ(IUNIT,'(A)',ERR=100,END=100) BUFLIN
         CALL UPPER(BUFLIN)
         N=N+NOCCUR(BUFLIN,PATT,-1)
      ENDDO

 100  CONTINUE
      BACKSPACE(IUNIT)
      RETURN
      END SUBROUTINE OCCURS_IN_FILE


!******************** SUBROUTINE COPYFILE ******************************
!
!***********************************************************************
      SUBROUTINE COPYFILE(SRC,TRG)
      IMPLICIT NONE
      CHARACTER(*) :: SRC,TRG
      CHARACTER(LEN=1) :: CHARAC
      OPEN(UNIT=99,FILE=SRC,STATUS='OLD')
      OPEN(UNIT=98,FILE=TRG,STATUS='REPLACE')
  10  READ(99,'(A)',ADVANCE='No',EOR=20,END=30) CHARAC
      WRITE(98,'(A)',ADVANCE='No') CHARAC
      GOTO 10
  20  WRITE(98,*)
      GOTO 10
  30  CLOSE(98); CLOSE(99)
      RETURN
      END SUBROUTINE COPYFILE


!************************* SETYLM_TST ********************************
!
! calculate the real spherical harmonics
! for a set of grid points up to LMAX
!
! YLM(:,1:LMAX) is set to Y_lm(hat x)
!
! hat x must be lying on the unit sphere
!
! written by Georg Kresse
!
!*********************************************************************

    SUBROUTINE SETYLM_TST(LYDIM,INDMAX,YLM,X,Y,Z)
      USE prec
      USE asa
      USE constant
      IMPLICIT NONE
      INTEGER LYDIM           ! maximum L
      INTEGER INDMAX          ! number of points (X,Y,Z)
      REAL(q) YLM(:,:)        ! spherical harmonics
      REAL(q) X(:),Y(:),Z(:)  ! x,y and z coordinates

! local variables
      REAL(q) FAK
      INTEGER IND,LSET,LM,LP,LMINDX,ISTART,IEND,LNEW,L,M,MP,IC
!-----------------------------------------------------------------------
! runtime check of workspace
!-----------------------------------------------------------------------
      IF ( UBOUND(YLM,2) < (LYDIM+1)**2) THEN
         WRITE(0,*)'internal ERROR: SETYLM, insufficient L workspace'
         STOP
      ENDIF

      IF ( UBOUND(YLM,1) < INDMAX) THEN
         WRITE(0,*)'internal ERROR: SETYLM, insufficient INDMAX workspace'
         STOP
      ENDIF

      FAK=1/(2._q * SQRT(PI))
!-----------------------------------------------------------------------
! here is the code for L=0, hard coded
!-----------------------------------------------------------------------
      IF (LYDIM <0) GOTO 100
!DIR$ IVDEP
!OCL NOVREC
      DO IND=1,INDMAX
        YLM(IND,1)=FAK
      ENDDO
!-----------------------------------------------------------------------
! here is the code for L=1, once again hard coded
!-----------------------------------------------------------------------
      IF (LYDIM <1) GOTO 100
!DIR$ IVDEP
!OCL NOVREC
      DO IND=1,INDMAX
        YLM(IND,2)  = (FAK*SQRT(3._q))*Y(IND)
        YLM(IND,3)  = (FAK*SQRT(3._q))*Z(IND)
        YLM(IND,4)  = (FAK*SQRT(3._q))*X(IND)
      ENDDO
!-----------------------------------------------------------------------
! code for L=2,
!-----------------------------------------------------------------------
      IF (LYDIM <2) GOTO 100
!DIR$ IVDEP
!OCL NOVREC
      DO IND=1,INDMAX
        YLM(IND,5)= (FAK*SQRT(15._q))  *X(IND)*Y(IND)
        YLM(IND,6)= (FAK*SQRT(15._q))  *Y(IND)*Z(IND)
        YLM(IND,7)= (FAK*SQRT(5._q)/2._q)*(3*Z(IND)*Z(IND)-1)
        YLM(IND,8)= (FAK*SQRT(15._q))  *X(IND)*Z(IND)
        YLM(IND,9)= (FAK*SQRT(15._q)/2._q)*(X(IND)*X(IND)-Y(IND)*Y(IND))
      ENDDO
!-----------------------------------------------------------------------
! code for L=3,
!-----------------------------------------------------------------------
      IF (LYDIM <3) GOTO 100

      DO IND=1,INDMAX
        YLM(IND,10)=(FAK*SQRT(7._q)/2._q) *Z(IND)*(2*Z(IND)*Z(IND)-3*X(IND)*X(IND)-3*Y(IND)*Y(IND))
        YLM(IND,11)=(FAK*SQRT(35._q/2._q)/2._q) *Y(IND)*(3*X(IND)*X(IND)-Y(IND)*Y(IND))
        YLM(IND,12)=(FAK*SQRT(35._q/2._q)/2._q) *X(IND)*(X(IND)*X(IND)-3*Y(IND)*Y(IND))
        YLM(IND,13)=(FAK*SQRT(105._q)/2._q) *Z(IND)*(X(IND)*X(IND)-Y(IND)*Y(IND))
        YLM(IND,14)=(FAK*SQRT(105._q)) *X(IND)*Y(IND)*Z(IND)
        YLM(IND,15)=(FAK*SQRT(21._q/2._q)/2._q) *Y(IND)*(4*Z(IND)*Z(IND)-X(IND)*X(IND)-Y(IND)*Y(IND))
        YLM(IND,16)=(FAK*SQRT(21._q/2._q)/2._q) *X(IND)*(4*Z(IND)*Z(IND)-X(IND)*X(IND)-Y(IND)*Y(IND))
      ENDDO

 100  CONTINUE

    END SUBROUTINE SETYLM_TST


!******************** SUBROUTINE MLWF_DHAM_DK ******************
!
!***********************************************************************

      SUBROUTINE MLWF_SETUP_HAM_WANNIER(WDES,W,HAM)
        USE wave_high
        USE full_kpoints
        IMPLICIT NONE
        TYPE(wavedes) WDES
        TYPE(wavespin) W
        COMPLEX(q) HAM(W%WDES%NB_TOT,W%WDES%NB_TOT,WDES%NKPTS,W%WDES%ISPIN)
        ! local variables
        GDEF, ALLOCATABLE :: U(:,:,:,:)
        TYPE(wavedes1) WDES1
        TYPE(wavefuna) WA
        INTEGER I,J,K,IP,IPP
        INTEGER ISP,NK,NKP,NW,NB,NBEXCL,NBwin,IK,B1,B2,IB
        COMPLEX(q) CTMP

        CALL CHECK_FULL_KPOINTS

        IF ((.NOT.ALLOCATED(U_matrix)).OR.(.NOT.ALLOCATED(U_matrix_opt))) THEN
           WRITE(*,*) 'MLWF_GET_ROTATION_MATRICES: ERROR: rotation matrices not available'
           STOP
        ENDIF

        NW=SIZE(U_matrix,1)
        NB=SIZE(U_matrix_opt,1)
        IF (SIZE(lexclude_band)/=WDES%NB_TOT) THEN
           WRITE(*,*) 'MLWF_GET_ROTATION_MATRICES: ERROR: lexclude_band array is corrupt:', &
                &   SIZE(lexclude_band),WDES%NB_TOT
           STOP
        ENDIF

        NBEXCL=0
        DO I=1,SIZE(lexclude_band)
           IF (lexclude_band(I)) NBEXCL=NBEXCL+1
        ENDDO
        IF ((NB+NBEXCL)/=WDES%NB_TOT) THEN
           WRITE(*,*) 'MLWF_GET_ROTATION_MATRICES: ERROR: inconsistent number of bands:', &
                &   NB+NBEXCL,WDES%NB_TOT
        ENDIF

        ALLOCATE(U(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN))

        spin: DO ISP=1,WDES%ISPIN
           kpoint: DO NK=1,WDES%NKPTS
              ! find the corresponding entry in KPOINTS_FULL, that contains
              ! the set of k-points used in the computation of the rotation matrices

              DO NKP=1,KPOINTS_FULL%NKPTS
                 IF (LIDENTICAL_KPOINT(WDES%VKPT(:,NK),KPOINTS_FULL%VKPT(:,NKP))) exit
              ENDDO
              IF (NKP>KPOINTS_FULL%NKPTS) THEN
                 WRITE(*,*) 'MLWF_GET_ROTATION_MATRICES: ERROR: no matching k-point found in KPOINTS_FULL',NK,WDES%VKPT(:,NK)
                 STOP
              ENDIF

              ! setup the effective rotation matrix U at the present k-point
              U(:,:,NK,ISP)=0._q; IP=0; IPP=0
              band: DO I=1,WDES%NB_TOT
                 IF (lexclude_band(I)) CYCLE band
                 IP=IP+1
                 IF (.NOT.lwindow(IP,NKP,ISP)) CYCLE band
                 IPP=IPP+1
                 DO J=1,NW
                    CTMP=0._q
                    DO K=1,NW
                       CTMP=CTMP+U_matrix_opt(IPP,K,NKP,ISP)*U_matrix(K,J,NKP,ISP)
                    ENDDO
                    U(I,J,NK,ISP)=CTMP
                 ENDDO
              ENDDO band

           ENDDO kpoint
        ENDDO spin

        ! calculate Hamiltonian matrices in Wannier basis at each k-point
        ! setup matrix here
        ! The Hamiltonian in the Wannier band basis
        !
        ! w .. Wannier bands
        ! b .. Bloch bands (cell-periodic part only)
        ! < w_m |H| w_n > =: H
        ! w_m = U_jm b_j
        ! eps = eps_j delta_ij
        !
        ! therefore:
        ! H_mn = U*_jm U_in < b_j|H|b_i > = U*_jm U_in eps_j delta_ji = U*_jm U_jn eps_j
        HAM=0
        spin1: DO ISP=1,W%WDES%ISPIN
           kpoint1: DO IK=1,WDES%NKPTS

              DO B1=1,W%WDES%NB_TOTK(IK,ISP)
                 DO B2=1,W%WDES%NB_TOTK(IK,ISP)
                    DO IB=1,W%WDES%NB_TOTK(IK,ISP)
                       HAM(B2,B1,IK,ISP)=HAM(B2,B1,IK,ISP)+GCONJG(U(IB,B2,IK,ISP))*W%CELTOT(IB,IK,ISP)*U(IB,B1,IK,ISP)
                    END DO
                 END DO
              END DO
           END DO kpoint1
        END DO spin1

        DEALLOCATE(U)
        RETURN
      END SUBROUTINE MLWF_SETUP_HAM_WANNIER


!******************** SUBROUTINE MLWF_DHAM_DK ******************
!
!***********************************************************************
      SUBROUTINE MLWF_DHAM_DK(W,KPOINTS,LATT_CUR,HAM,DHAM)
        USE constant
        USE pseudo
        USE lattice
        USE mkpoints
        USE full_kpoints
        USE wave_high
        IMPLICIT NONE
        TYPE(wavespin) W
        TYPE(latt) LATT_CUR
        TYPE(kpoints_struct) KPOINTS
        COMPLEX(q) DHAM(W%WDES%NB_TOT,W%WDES%NB_TOT,KPOINTS%NKPTS,W%WDES%ISPIN,3), &
             HAM(W%WDES%NB_TOT,W%WDES%NB_TOT,KPOINTS_FULL%NKPTS,W%WDES%ISPIN)
        ! local variables
        TYPE (wavedes1) WDESK
        TYPE (wavefun1) WK
        INTEGER ISP,IK,NK,IB,ISPINOR,B1,B2
        INTEGER IDIR,J,IDELTA,IORDER,ISGN
        REAL(q) K0(3),K(3),DK(3),KP(3)

        ! calculate derivative of Hamiltonian at each k-point
        IORDER=4

        DK(1)=1._q/REAL(KPOINTS_FULL%NKPX,KIND=q)
        DK(2)=1._q/REAL(KPOINTS_FULL%NKPY,KIND=q)
        DK(3)=1._q/REAL(KPOINTS_FULL%NKPZ,KIND=q)

        DHAM=0

        spin2: DO ISP=1,W%WDES%ISPIN
           kpoint2: DO IK=1,KPOINTS%NKPTS

              K0(:)=KPOINTS%VKPT(:,IK)

              dir: DO IDIR=1,3
                 delta: DO IDELTA=1,IORDER
                    sgn: DO ISGN=-1,1,2
                       K(:)=K0(:)
                       K(IDIR)=K0(IDIR)+ISGN*IDELTA*DK(IDIR)
                       NK=KPOINT_IN_FULL_GRID(K,KPOINTS_FULL)
                       KP(:)=KPOINTS_FULL%VKPT(:,NK)

                       ! non-cartesian coordinates !!!
                       DHAM(:,:,IK,ISP,IDIR)=DHAM(:,:,IK,ISP,IDIR)+CMPLX(ISGN*HAM(:,:,NK,ISP)*FAC(IDELTA,IORDER)/DK(IDIR)/TPI/2._q,KIND=q)
                       DO J=1,3
!                          DHAM(:,:,IK,ISP,J)=DHAM(:,:,IK,ISP,J)+CMPLX(ISGN*LATT_CUR%A(J,IDIR)*HAM(:,:,NK,ISP)*FAC(IDELTA,IORDER)/DK(IDIR)/TPI/2._q,KIND=q)
                          ! derivative w.r.t. to dimensionless coordinate
!                          DHAM(:,:,IK,ISP,J)=DHAM(:,:,IK,ISP,J)+CMPLX(ISGN*HAM(:,:,NK,ISP)*FAC(IDELTA,IORDER)/DK(IDIR)/TPI/2._q,KIND=q)


                       ENDDO

                    ENDDO sgn
                 ENDDO delta
              ENDDO dir

           ENDDO kpoint2
        ENDDO spin2

        RETURN
      END SUBROUTINE MLWF_DHAM_DK


      SUBROUTINE MLWF_GET_U(WDES,U)
        USE wave_high
        USE full_kpoints
        IMPLICIT NONE
        TYPE(wavedes) WDES
        GDEF :: U(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)
        ! local variables
        INTEGER I,J,K,IP,IPP
        INTEGER ISP,NK,NKP,NW,NB,NBEXCL
        COMPLEX(q) CTMP

        CALL CHECK_FULL_KPOINTS

        IF ((.NOT.ALLOCATED(U_matrix)).OR.(.NOT.ALLOCATED(U_matrix_opt))) THEN
           WRITE(*,*) 'MLWF_GET_ROTATION_MATRICES: ERROR: rotation matrices not available'
           STOP
        ENDIF

        NW=SIZE(U_matrix,1)
        NB=SIZE(U_matrix_opt,1)
        IF (SIZE(lexclude_band)/=WDES%NB_TOT) THEN
           WRITE(*,*) 'MLWF_GET_ROTATION_MATRICES: ERROR: lexclude_band array is corrupt:', &
                &   SIZE(lexclude_band),WDES%NB_TOT
           STOP
        ENDIF

        NBEXCL=0
        DO I=1,SIZE(lexclude_band)
           IF (lexclude_band(I)) NBEXCL=NBEXCL+1
        ENDDO
        IF ((NB+NBEXCL)/=WDES%NB_TOT) THEN
           WRITE(*,*) 'MLWF_GET_ROTATION_MATRICES: ERROR: inconsistent number of bands:', &
                &   NB+NBEXCL,WDES%NB_TOT
        ENDIF

        spin: DO ISP=1,WDES%ISPIN
           kpoint: DO NK=1,WDES%NKPTS
              ! find the corresponding entry in KPOINTS_FULL, that contains
              ! the set of k-points used in the computation of the rotation matrices

              DO NKP=1,KPOINTS_FULL%NKPTS
                 IF (LIDENTICAL_KPOINT(WDES%VKPT(:,NK),KPOINTS_FULL%VKPT(:,NKP))) exit
              ENDDO
              IF (NKP>KPOINTS_FULL%NKPTS) THEN
                 WRITE(*,*) 'MLWF_GET_ROTATION_MATRICES: ERROR: no matching k-point found in KPOINTS_FULL',NK,WDES%VKPT(:,NK)
                 STOP
              ENDIF

              ! setup the effective rotation matrix U at the present k-point
              U(:,:,NK,ISP)=0._q; IP=0; IPP=0
              band: DO I=1,WDES%NB_TOT
                 IF (lexclude_band(I)) CYCLE band
                 IP=IP+1
                 IF (.NOT.lwindow(IP,NKP,ISP)) CYCLE band
                 IPP=IPP+1
                 DO J=1,NW
                    CTMP=0._q
                    DO K=1,NW
                       CTMP=CTMP+U_matrix_opt(IPP,K,NKP,ISP)*U_matrix(K,J,NKP,ISP)
                    ENDDO
                    U(I,J,NK,ISP)=CTMP
                 ENDDO
              ENDDO band
           ENDDO kpoint
        ENDDO spin
      END SUBROUTINE MLWF_GET_U


      SUBROUTINE MLWF_ROTATE_HAM(WDES,W,HAM,U)
        USE wave_high
        USE full_kpoints
        IMPLICIT NONE
        TYPE(wavedes) WDES
        TYPE(wavespin) W
        COMPLEX(q) HAM(W%WDES%NB_TOT,W%WDES%NB_TOT,WDES%NKPTS,W%WDES%ISPIN)
        GDEF :: U(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)
        ! local variables
        INTEGER ISP,IK,B1,B2,IB
        ! calculate Hamiltonian matrices in Wannier basis at each k-point
        ! setup matrix here
        ! The Hamiltonian in the Wannier band basis
        !
        ! w .. Wannier bands
        ! b .. Bloch bands (cell-periodic part only)
        ! < w_m |H| w_n > =: H
        ! w_m = U_jm b_j
        ! eps = eps_j delta_ij
        !
        ! therefore:
        ! H_mn = U*_jm U_in < b_j|H|b_i > = U*_jm U_in eps_j delta_ji = U*_jm U_jn eps_j
        HAM=0
        spin1: DO ISP=1,W%WDES%ISPIN
           kpoint1: DO IK=1,WDES%NKPTS

              DO B1=1,W%WDES%NB_TOTK(IK,ISP)
                 DO B2=1,W%WDES%NB_TOTK(IK,ISP)
                    DO IB=1,W%WDES%NB_TOTK(IK,ISP)
                       HAM(B2,B1,IK,ISP)=HAM(B2,B1,IK,ISP)+GCONJG(U(IB,B2,IK,ISP))*W%CELTOT(IB,IK,ISP)*U(IB,B1,IK,ISP)
                    END DO
                 END DO
              END DO
           END DO kpoint1
        END DO spin1
        RETURN
      END SUBROUTINE MLWF_ROTATE_HAM

!***********************************************************************
! writes U_(nm)^K files
!***********************************************************************

      SUBROUTINE WRITE_WANPROJ_FILE(WDES)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         COMPLEX(q), ALLOCATABLE :: U(:,:)
         INTEGER I,J,K,IP,IPP
         INTEGER ISP,NK,NKP,NW,NB
         COMPLEX(q) CTMP

         CALL CHECK_FULL_KPOINTS
         IF ((.NOT.ALLOCATED(U_matrix)).OR.&
            (.NOT.ALLOCATED(U_matrix_opt))) THEN
            WRITE(*,*) &
               'MLWF_ROTATE_ORBITALS_FULLK: ERROR: rotation '//&
               'matrices not available'
            STOP
         ENDIF
         NW=SIZE(U_matrix,1)
         NB=SIZE(U_matrix_opt,1)
         ALLOCATE(U(WDES%NB_TOT,WDES%NB_TOT))
         U=0

         OPEN(UNIT=99, FILE=WANPROJFILE, STATUS='REPLACE')

         WRITE(99,*)WDES%ISPIN,WDES%NKPTS,WDES%NB_TOT,NW
         DO NK=1,WDES%NKPTS
            WRITE(99,999)NK,WDES%VKPT(:,NK)
         ENDDO
         spin: DO ISP=1,WDES%ISPIN
            kpts: DO NK=1,WDES%NKPTS
               ! find the corresponding entry in kpt_latt, that contains
               ! the set of k-points used in the computation of the rotation matrices
               DO NKP=1,SIZE(kpt_latt,2)
                  IF (LIDENTICAL_KPOINT(WDES%VKPT(:,NK),kpt_latt(:,NKP))) exit
               ENDDO
               IF (NKP>KPOINTS_FULL%NKPTS) THEN
                  WRITE(*,*) &
                      'ERROR: no matching k-point found in kpt_latt',NK
                  STOP
               ENDIF

               ! setup the effective rotation matrix U at the present k-point
               IP=0; IPP=0
               bands1: DO I=1,WDES%NB_TOT
                  IF (lexclude_band(I)) CYCLE bands1
                  IP=IP+1
                  IF (.NOT.lwindow(IP,NKP,ISP)) CYCLE bands1
                  IPP=IPP+1
                  DO J=1,NW
                     CTMP=0._q
                     DO K=1,NW
                        CTMP=CTMP+U_matrix_opt(IPP,K,NKP,ISP)*U_matrix(K,J,NKP,ISP)
                     ENDDO
                  ENDDO
               ENDDO bands1

               WRITE(99,1000)ISP,IPP,WDES%VKPT(:,NK)

               ! setup the effective rotation matrix U at the present k-point
               IP=0; IPP=0
               bands: DO I=1,WDES%NB_TOT
                  IF (lexclude_band(I)) CYCLE bands
                  IP=IP+1
                  IF (.NOT.lwindow(IP,NKP,ISP)) CYCLE bands
                  IPP=IPP+1
                  DO J=1,NW
                     CTMP=0._q
                     DO K=1,NW
                        CTMP=CTMP+U_matrix_opt(IPP,K,NKP,ISP)*U_matrix(K,J,NKP,ISP)
                     ENDDO
                     WRITE(99,1001)I,J,CTMP
                  ENDDO
               ENDDO bands


            ENDDO kpts
         ENDDO spin
         CLOSE(99)

999      FORMAT(I6, 3F26.16)
1000     FORMAT(2I6,3F26.16)
1001     FORMAT(2I6,2F26.16)
      ENDSUBROUTINE WRITE_WANPROJ_FILE
      END MODULE mlwf
