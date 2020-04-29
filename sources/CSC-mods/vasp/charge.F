#include "symbol.inc"
!************************* SUBROUTINE FFT_RC_SCALE *********************
! RCS:  $Id: charge.F,v 1.5 2002/08/14 13:59:37 kresse Exp $
!
! subroutine transforms a real space chargedensity to
! reciprocal space applying the real to complex FFT transformation
!
!***********************************************************************

  SUBROUTINE FFT_RC_SCALE(CHDENR,CHDEN,GRID)
    USE prec
    USE mpimy
    USE mgrid
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)
    
    TYPE (grid_3d)     GRID
    
    COMPLEX(q) CHDEN(GRID%MPLWV)
    RGRID      CHDENR(GRID%RL%NP)
    
    RINPLW=1.0_q/GRID%NPLWV
    CALL  RL_ADD(CHDENR,RINPLW,CHDENR,0.0_q,CHDEN,GRID)
    CALL FFT3D(CHDEN,GRID,-1)
    RETURN
  END SUBROUTINE FFT_RC_SCALE
  

MODULE charge
  USE prec
  
CONTAINS

!********************* SOFT_CHARGE *************************************
!
! subroutine constructs the electronic charge density according
! to the current wavefunctions and fermi-weights
!  rho(r) = \sum_n,k phi_nk(r) phi_nk^*(r) * f_nk
! where f_nk are the single electron occupancies (FERWE)
! the phi_nk are the pseudowavefunctions 
!
!***********************************************************************

  SUBROUTINE SOFT_CHARGE_DMFT_IBZ(GRID,GRID_SOFT,W,WDES, CHDEN)
    USE wave_high
    USE gridq
    USE hamil
    USE mlwf
    USE pseudo
    USE lattice
    USE sym_prec

    IMPLICIT COMPLEX(q) (C)

    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
    TYPE (grid_3d)     GRID,GRID_SOFT
    COMPLEX(q)   CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
  ! local
    TYPE (wavedes1)     WDES1
    TYPE (wavefun1)     W1
    INTEGER ISPINOR
    INTEGER, PARAMETER :: NSTRIPD=2
    TYPE (REDIS_PW_CTR),POINTER :: H_PW
    TYPE (GRIDQUANT) :: CHARGE_REAL_SPACE, CHARGE
    INTEGER I
    INTEGER :: ierr,N2,NBG,ISP, NK, N
    REAL(q) :: WEIGHT
    COMPLEX(q),ALLOCATABLE  :: CR_TOT(:,:),CR_TOT2(:,:)

    PROFILING_START('soft_charge')

    IF (W%OVER_BAND) THEN
#ifdef MPI
       NCPU=WDES%COMM_INTER%NCPU ! number of procs involved in band dis.
#else
       NCPU=1
#endif
       NSTRIP=MIN(NSTRIPD,WDES%NBANDS)
       CALL REDIS_PW_ALLOC(WDES, NSTRIP, H_PW)
    ENDIF

    CALL SETWDES(WDES, WDES1, 0)
    CALL NEWWAV_R(W1, WDES1)

    CALL GENERATE_GRID_QUANTITY(CHARGE, GRID_SOFT, CHDEN) 
    CALL ALLOCATE_GRID_QUANTITY_FORCE_RL(CHARGE_REAL_SPACE, GRID, GRID_SOFT, WDES%NCDIJ)

    CHARGE_REAL_SPACE=0.0_q
    CHARGE_REAL_SPACE%REALSPACE =.TRUE.

    spin: DO ISP=1,WDES%ISPIN

       kpoints: DO NK=1,WDES%NKPTS
#ifdef MPI
          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif

          CALL SETWDES(WDES, WDES1, NK)

          IF (W%OVER_BAND) THEN
             DO N=1,NSTRIP
                CALL REDIS_PW_START(WDES, W%CW(1,N,NK,ISP), N, H_PW)
             ENDDO
          ENDIF

          !write(*,*) NK
          !write(*,*) W%WDES%NB_LOW, WDES1%NB_LOW
          num_band_max=band_win_(2,NK)-band_win_(1,NK)+1
          ALLOCATE(CR_TOT(WDES1%GRID%MPLWV,num_band_max))
          !ALLOCATE(CR_TOT2(WDES1%GRID%MPLWV,num_band_max))
          CR_TOT=0._q;!CR_TOT2=0._q;
          DO N=1,WDES1%NBANDS
             CALL SETWAV(W,W1,WDES1,N,ISP)
             IF (W%OVER_BAND) THEN
                CALL REDIS_PW_STOP (WDES, W%CW(1,N,NK,ISP), N, H_PW)
                IF (N+NSTRIP<=WDES%NBANDS) &
                     CALL REDIS_PW_START(WDES, W%CW(1,N+NSTRIP,NK,ISP), N+NSTRIP, H_PW)
             ENDIF
             CALL W1_COPY(ELEMENT(W, WDES1, N, ISP), W1)
             !write(*,*) 'a',NK
             CALL FFTWAV_W1(W1)
             !write(*,*) N, band_win_(1,NK), band_win_(2,NK)
             NBG=W%WDES%NB_LOW+(N-1)*W%WDES%NB_PAR
             if (NBG .GE. band_win_(1,NK) .AND. NBG .LE. band_win_(2,NK)) then
               CR_TOT(:,NBG-band_win_(1,NK)+1)=W1%CR(:)
             endif
          ENDDO
          !write(*,*) 'b',NK
          !write(*,*) WDES1%GRID%MPLWV,num_band_max,CR_TOT(1,1),CR_TOT2(1,1)
          !CALL MPI_ALLREDUCE(CR_TOT,CR_TOT2,WDES1%GRID%MPLWV*num_band_max,MPI_DOUBLE_COMPLEX,MPI_SUM,WDES%COMM%MPI_COMM,ierr)
          !DEALLOCATE(CR_TOT)
          !write(*,*) 'c',NK
!
!#ifdef MPI
!          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) DEALLOCATE(CR_TOT2); CYCLE
!#endif
          !write(*,*) 'd',NK

          band: DO N=1,WDES%NBANDS
             CALL SETWAV(W,W1,WDES1,N,ISP)

             IF (W%OVER_BAND) THEN
                CALL REDIS_PW_STOP (WDES, W%CW(1,N,NK,ISP), N, H_PW)
                IF (N+NSTRIP<=WDES%NBANDS) &
                     CALL REDIS_PW_START(WDES, W%CW(1,N+NSTRIP,NK,ISP), N+NSTRIP, H_PW)
             ENDIF

             CALL W1_COPY(ELEMENT(W, WDES1, N, ISP), W1)
             CALL FFTWAV_W1(W1)
             NBG=W%WDES%NB_LOW+(N-1)*W%WDES%NB_PAR
             if (NBG .GE. band_win_(1,NK) .AND. NBG .LE. band_win_(2,NK)) then
               W1%CR=0._q
               DO N2=1,num_band_max
                 W1%CR(:)=W1%CR(:)+CR_TOT(:,N2)*DMFT_EVEC_(N2,NBG-band_win_(1,NK)+1,NK)
               ENDDO
               WEIGHT=W%WDES%RSPIN*W%WDES%WTKPT(NK)*DMFT_EVAL_(NBG-band_win_(1,NK)+1,NK)
             else
               WEIGHT=W%WDES%RSPIN*W%WDES%WTKPT(NK)*W%FERWE(N,NK,ISP)
             endif

             !WEIGHT=WDES%RSPIN*WDES%WTKPT(NK)*W%FERWE(N,NK,ISP)
             IF (WEIGHT==0) CYCLE
             CALL PW_CHARGE(WDES1, CHARGE_REAL_SPACE%RG(1,ISP), SIZE(CHARGE_REAL_SPACE%RG,1), W1%CR(1), W1%CR(1), WEIGHT )
          ENDDO band
          DEALLOCATE(CR_TOT)
       ENDDO kpoints
    ENDDO spin

    IF (W%OVER_BAND) THEN
       W%OVER_BAND=.FALSE.
       CALL REDIS_PW_DEALLOC(H_PW)
    ENDIF

    ! fourier-transformation of charge-density using GRID_SOFT
    ! only input data from first in-band-group is used, and all nodes
    ! are involved in the FFT, final result is distributed among nodes
    ! (see SET_RL_GRID() in mgrid.F, and M_divide() in mpi.F)

    ! merge charge from all bands
#ifdef MPI
    IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
       DO I=1,CHARGE_REAL_SPACE%NCDIJ
#ifdef realmode
          CALLMPI( M_sum_d(WDES%COMM_KINTER, CHARGE_REAL_SPACE%RG(1,I),CHARGE_REAL_SPACE%GRID%RL%NP) )
#else
          CALLMPI( M_sum_z(WDES%COMM_KINTER, CHARGE_REAL_SPACE%RG(1,I),CHARGE_REAL_SPACE%GRID%RL%NP) )
#endif
       END DO
    END IF
#endif

    CALL SUMRL_GQ( CHARGE_REAL_SPACE, CHARGE, WDES%COMM_INTER)
    
    CALL FFT_GQ(CHARGE)
      
    ! set the charge-density of unbalanced lattice vectors to 0
    CALL SETUNB_GQ(CHARGE)

    CALL DELWAV_R(W1)
    CALL DEALLOCATE_GRID_QUANTITY(CHARGE_REAL_SPACE)

    PROFILING_STOP('soft_charge')

  END SUBROUTINE SOFT_CHARGE_DMFT_IBZ

  SUBROUTINE SOFT_CHARGE(GRID,GRID_SOFT,W,WDES, CHDEN)
    USE wave_high
    USE gridq
    USE hamil
    IMPLICIT COMPLEX(q) (C)

    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
    TYPE (grid_3d)     GRID,GRID_SOFT
    COMPLEX(q)   CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
  ! local
    TYPE (wavedes1)     WDES1
    TYPE (wavefun1)     W1
    INTEGER ISPINOR
    INTEGER, PARAMETER :: NSTRIPD=2
    TYPE (REDIS_PW_CTR),POINTER :: H_PW
    TYPE (GRIDQUANT) :: CHARGE_REAL_SPACE, CHARGE
    INTEGER I

    PROFILING_START('soft_charge')

    IF (W%OVER_BAND) THEN
#ifdef MPI
       NCPU=WDES%COMM_INTER%NCPU ! number of procs involved in band dis.
#else
       NCPU=1
#endif
       NSTRIP=MIN(NSTRIPD,WDES%NBANDS)
       CALL REDIS_PW_ALLOC(WDES, NSTRIP, H_PW)
    ENDIF

    CALL SETWDES(WDES, WDES1, 0)
    CALL NEWWAV_R(W1, WDES1)

    CALL GENERATE_GRID_QUANTITY(CHARGE, GRID_SOFT, CHDEN) 
    CALL ALLOCATE_GRID_QUANTITY_FORCE_RL(CHARGE_REAL_SPACE, GRID, GRID_SOFT, WDES%NCDIJ)

    CHARGE_REAL_SPACE=0.0_q
    CHARGE_REAL_SPACE%REALSPACE =.TRUE.

    spin: DO ISP=1,WDES%ISPIN

       kpoints: DO NK=1,WDES%NKPTS
#ifdef MPI
          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
          CALL SETWDES(WDES, WDES1, NK)
          

          IF (W%OVER_BAND) THEN
             DO N=1,NSTRIP
                CALL REDIS_PW_START(WDES, W%CW(1,N,NK,ISP), N, H_PW)
             ENDDO
          ENDIF

          band: DO N=1,WDES%NBANDS
             CALL SETWAV(W,W1,WDES1,N,ISP)

             IF (W%OVER_BAND) THEN
                CALL REDIS_PW_STOP (WDES, W%CW(1,N,NK,ISP), N, H_PW)
                IF (N+NSTRIP<=WDES%NBANDS) &
                     CALL REDIS_PW_START(WDES, W%CW(1,N+NSTRIP,NK,ISP), N+NSTRIP, H_PW)
             ENDIF

             WEIGHT=WDES%RSPIN*WDES%WTKPT(NK)*W%FERWE(N,NK,ISP)
             IF (WEIGHT==0) CYCLE
             CALL FFTWAV_W1(W1)
             CALL PW_CHARGE(WDES1, CHARGE_REAL_SPACE%RG(1,ISP), SIZE(CHARGE_REAL_SPACE%RG,1), W1%CR(1), W1%CR(1), WEIGHT )
          ENDDO band
       ENDDO kpoints
    ENDDO spin

    IF (W%OVER_BAND) THEN
       W%OVER_BAND=.FALSE.
       CALL REDIS_PW_DEALLOC(H_PW)
    ENDIF

    ! fourier-transformation of charge-density using GRID_SOFT
    ! only input data from first in-band-group is used, and all nodes
    ! are involved in the FFT, final result is distributed among nodes
    ! (see SET_RL_GRID() in mgrid.F, and M_divide() in mpi.F)

    ! merge charge from all bands
#ifdef MPI
    IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
       DO I=1,CHARGE_REAL_SPACE%NCDIJ
#ifdef realmode
          CALLMPI( M_sum_d(WDES%COMM_KINTER, CHARGE_REAL_SPACE%RG(1,I),CHARGE_REAL_SPACE%GRID%RL%NP) )
#else
          CALLMPI( M_sum_z(WDES%COMM_KINTER, CHARGE_REAL_SPACE%RG(1,I),CHARGE_REAL_SPACE%GRID%RL%NP) )
#endif
       END DO
    END IF
#endif

    CALL SUMRL_GQ( CHARGE_REAL_SPACE, CHARGE, WDES%COMM_INTER)
    
    CALL FFT_GQ(CHARGE)
      
    ! set the charge-density of unbalanced lattice vectors to 0
    CALL SETUNB_GQ(CHARGE)

    CALL DELWAV_R(W1)
    CALL DEALLOCATE_GRID_QUANTITY(CHARGE_REAL_SPACE)

    PROFILING_STOP('soft_charge')

  END SUBROUTINE SOFT_CHARGE

  SUBROUTINE SOFT_CHARGE_DMFT(GRID,GRID_SOFT,W,WDES, CHDEN ,P, LATT_CUR)
    USE wave_high
    USE gridq
    USE hamil
    USE full_kpoints
    USE pseudo
    USE lattice
    USE sym_prec
    USE mlwf
    USE pead

    IMPLICIT NONE

    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
    TYPE (grid_3d)     GRID,GRID_SOFT
    COMPLEX(q)   CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
  ! local
    INTEGER :: ISP, NK, N, ISP_IRZ
    REAL(q) :: WEIGHT
    LOGICAL :: LSHIFT
    TYPE (wavedes1)     WDES1, WDES1_IRZ
    TYPE (wavefun1)     W1
    INTEGER ISPINOR
    TYPE (REDIS_PW_CTR),POINTER :: H_PW
    TYPE (GRIDQUANT) :: CHARGE_REAL_SPACE, CHARGE
    INTEGER I
    INTEGER :: ierr,N2,NBG
    COMPLEX(q),ALLOCATABLE  :: CR_TOT(:,:),CR_TOT2(:,:)
    TYPE(rotation_handle), POINTER :: ROT_HANDLE
    TYPE(potcar) P(:)
    TYPE (latt)     LATT_CUR
    TYPE (wavespin)    WORIG
    TYPE(wavefuna) WK


    PROFILING_START('soft_charge_sym')

    IF (W%OVER_BAND) THEN
       W%OVER_BAND=.FALSE.
       WRITE(0,*) "internal error in SOFT_CHARGE_SYM: W%OVER_BAND is not supported"
       STOP
    ENDIF
    WORIG=W
    WORIG%WDES=>WDES_FULL_PEAD
    !write(*,*) WDES_FULL_PEAD%NKPTS

    CALL CHECK_FULL_KPOINTS

    CALL SETWDES(WORIG%WDES, WDES1, 0)
    CALL NEWWAV(W1 , WDES1,.TRUE.)

    CALL GENERATE_GRID_QUANTITY(CHARGE, GRID_SOFT, CHDEN) 
    CALL ALLOCATE_GRID_QUANTITY_FORCE_RL(CHARGE_REAL_SPACE, GRID, GRID_SOFT, WDES%NCDIJ)

    CHARGE_REAL_SPACE=0.0_q
    CHARGE_REAL_SPACE%REALSPACE =.TRUE.
    !write(*,*) 'hihi'

    spin: DO ISP=1,WORIG%WDES%ISPIN

       kpoints: DO NK=1,KPOINTS_FULL%NKPTS
#ifdef MPI
          IF (MOD(NK-1,WORIG%WDES%COMM_KINTER%NCPU).NE.WORIG%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
          !write(*,*) 'hihi'
          NULLIFY(ROT_HANDLE)
          CALL SETWDES(WORIG%WDES,WDES1,NK)
          !write(*,*) 'hihi2'
          IF (NK==KPOINTS_FULL%NEQUIV(NK)) THEN
             ! k is a kpoint in the IBZ
             WK=ELEMENTS(WORIG,WDES1,ISP)
          ELSE
             ! k is not a kpoint in the IBZ
             CALL NEWWAVA(WK,WDES1,WDES1%NBANDS)
             CALL PEAD_WA_ROTATE(WORIG,P,LATT_CUR,ISP,WK)
          ENDIF

          !CALL SETWDES(WDES, WDES1, NK)
          !!CALL SETWDES(WDES, WDES1_IRZ, KPOINTS_FULL%NEQUIV(NK))
          !
          !ISP_IRZ=ISP
          !IF (KPOINTS_FULL%SPINFLIP(NK)==1) THEN
          !   ISP_IRZ=3-ISP
          !ENDIF

          num_band_max=band_win_(2,NK)-band_win_(1,NK)+1
          ALLOCATE(CR_TOT(WDES1%GRID%MPLWV,num_band_max))
          !ALLOCATE(CR_TOT2(WDES1%GRID%MPLWV,num_band_max2))
          CR_TOT=0._q!;CR_TOT2=0._q;
          DO N=1,WDES1%NBANDS
             CALL W1_COPY(ELEMENT(WK, N), W1)
             CALL FFTWAV_W1(W1)
             NBG=WORIG%WDES%NB_LOW+(N-1)*WORIG%WDES%NB_PAR
             if (NBG .GE. band_win_(1,NK) .AND. NBG .LE. band_win_(2,NK)) then
               CR_TOT(:,NBG-band_win_(1,NK)+1)=W1%CR(:)
             endif
          ENDDO
          !write(*,*) 'hihi3'
         
          band: DO N=1,WORIG%WDES%NBANDS
             CALL W1_COPY(ELEMENT(WK, N), W1)
             CALL FFTWAV_W1(W1)
             NBG=WORIG%WDES%NB_LOW+(N-1)*WORIG%WDES%NB_PAR
             if (NBG .GE. band_win_(1,NK) .AND. NBG .LE. band_win_(2,NK)) then
               W1%CR=0._q
               DO N2=1,num_band_max
                 W1%CR(:)=W1%CR(:)+CR_TOT(:,N2)*DMFT_EVEC_(N2,NBG-band_win_(1,NK)+1,NK)
               ENDDO
               WEIGHT=WORIG%WDES%RSPIN*KPOINTS_FULL%WTKPT(NK)*DMFT_EVAL_(NBG-band_win_(1,NK)+1,NK)
             else
               WEIGHT=WORIG%WDES%RSPIN*KPOINTS_FULL%WTKPT(NK)*W%FERWE(N,KPOINTS_FULL%NEQUIV(NK),ISP)
             endif

             !WEIGHT=WDES%RSPIN*KPOINTS_FULL%WTKPT(NK)*W%FERWE(N,KPOINTS_FULL%NEQUIV(NK),ISP)

             IF (WEIGHT==0) CYCLE

             !IF (NK <= WDES%NKPTS) THEN
             !   CALL W1_COPY(ELEMENT(W, WDES1, N, ISP), W1)
             !   CALL FFTWAV_W1(W1)
             !ELSE
#ifndef gammareal
           ! 
           ! symmetry must be considered if the wavefunctions for this
           ! k-point NK (containing all k-points in the entire BZ) 
           ! are not stored in W
           !
                !LSHIFT=.FALSE.
                !IF ((ABS(KPOINTS_FULL%TRANS(1,NK)) >TINY) .OR. &
                !     (ABS(KPOINTS_FULL%TRANS(2,NK))>TINY) .OR. &
                !     (ABS(KPOINTS_FULL%TRANS(3,NK))>TINY)) LSHIFT=.TRUE.
                !CALL W1_ROTATE_AND_FFT_NO_PROJ(W1, ELEMENT(W, WDES1_IRZ, N, ISP_IRZ), LSHIFT)
#endif
             !ENDIF


             CALL PW_CHARGE(WDES1, CHARGE_REAL_SPACE%RG(1,ISP), SIZE(CHARGE_REAL_SPACE%RG,1), W1%CR(1), W1%CR(1), WEIGHT )
          ENDDO band
          CALL DEALLOCATE_ROT_HANDLE(ROT_HANDLE)
          IF (NK/=KPOINTS_FULL%NEQUIV(NK)) CALL DELWAVA(WK)
          !CALL DELWAVA(WK)
          DEALLOCATE(CR_TOT)
       ENDDO kpoints
       

    ENDDO spin


    ! fourier-transformation of charge-density using GRID_SOFT
    ! only input data from first in-band-group is used, and all nodes
    ! are involved in the FFT, final result is distributed among nodes
    ! (see SET_RL_GRID() in mgrid.F, and M_divide() in mpi.F)

    ! merge charge from all bands
#ifdef MPI
    IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
       DO I=1,CHARGE_REAL_SPACE%NCDIJ
#ifdef realmode
          CALLMPI( M_sum_d(WDES%COMM_KINTER, CHARGE_REAL_SPACE%RG(1,I),CHARGE_REAL_SPACE%GRID%RL%NP) )
#else
          CALLMPI( M_sum_z(WDES%COMM_KINTER, CHARGE_REAL_SPACE%RG(1,I),CHARGE_REAL_SPACE%GRID%RL%NP) )
#endif
       END DO
    END IF
#endif

    CALL SUMRL_GQ( CHARGE_REAL_SPACE, CHARGE, WDES%COMM_INTER)
    
    CALL FFT_GQ(CHARGE)
      
    ! set the charge-density of unbalanced lattice vectors to 0
    CALL SETUNB_GQ(CHARGE)

    CALL DELWAV(W1, .TRUE.)
    CALL DEALLOCATE_GRID_QUANTITY(CHARGE_REAL_SPACE)

    PROFILING_STOP('soft_charge_sym')

  END SUBROUTINE SOFT_CHARGE_DMFT

!********************* SOFT_CHARGE *************************************
!
! this subroutine constructs the electronic charge density applying
! symmetry to determine the wavefunctions at each k-point
! it should be fully compatible with the Fock like routines but
! of course is slightly slower than SOFT_CHARGE
! if this routine is used to construct the charge density no
! a posteriori symmetrization of the charge density should be required
!
!***********************************************************************

  SUBROUTINE SOFT_CHARGE_SYM(GRID,GRID_SOFT,W,WDES, CHDEN )
    USE wave_high
    USE gridq
    USE hamil
    USE full_kpoints
    USE pseudo
    USE lattice
    USE sym_prec
    IMPLICIT NONE

    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
    TYPE (grid_3d)     GRID,GRID_SOFT
    COMPLEX(q)   CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
  ! local
    INTEGER :: ISP, NK, N, ISP_IRZ
    REAL(q) :: WEIGHT
    LOGICAL :: LSHIFT
    TYPE (wavedes1)     WDES1, WDES1_IRZ
    TYPE (wavefun1)     W1
    INTEGER ISPINOR
    TYPE (REDIS_PW_CTR),POINTER :: H_PW
    TYPE (GRIDQUANT) :: CHARGE_REAL_SPACE, CHARGE
    INTEGER I

    PROFILING_START('soft_charge_sym')

    IF (W%OVER_BAND) THEN
       W%OVER_BAND=.FALSE.
       WRITE(0,*) "internal error in SOFT_CHARGE_SYM: W%OVER_BAND is not supported"
       STOP
    ENDIF
    CALL CHECK_FULL_KPOINTS

    CALL SETWDES(WDES, WDES1, 0)
    CALL NEWWAV(W1 , WDES1,.TRUE.)

    CALL GENERATE_GRID_QUANTITY(CHARGE, GRID_SOFT, CHDEN) 
    CALL ALLOCATE_GRID_QUANTITY_FORCE_RL(CHARGE_REAL_SPACE, GRID, GRID_SOFT, WDES%NCDIJ)

    CHARGE_REAL_SPACE=0.0_q
    CHARGE_REAL_SPACE%REALSPACE =.TRUE.

    spin: DO ISP=1,WDES%ISPIN

       kpoints: DO NK=1,KPOINTS_FULL%NKPTS
#ifdef MPI
          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
          CALL SETWDES(WDES, WDES1, NK)
          CALL SETWDES(WDES, WDES1_IRZ, KPOINTS_FULL%NEQUIV(NK))
          
          ISP_IRZ=ISP
          IF (KPOINTS_FULL%SPINFLIP(NK)==1) THEN
             ISP_IRZ=3-ISP
          ENDIF
         
          band: DO N=1,WDES%NBANDS
             WEIGHT=WDES%RSPIN*KPOINTS_FULL%WTKPT(NK)*W%FERWE(N,KPOINTS_FULL%NEQUIV(NK),ISP)

             IF (WEIGHT==0) CYCLE

             IF (NK <= WDES%NKPTS) THEN
                CALL W1_COPY(ELEMENT(W, WDES1, N, ISP), W1)
                CALL FFTWAV_W1(W1)
             ELSE
#ifndef gammareal
           ! 
           ! symmetry must be considered if the wavefunctions for this
           ! k-point NK (containing all k-points in the entire BZ) 
           ! are not stored in W
           !
                LSHIFT=.FALSE.
                IF ((ABS(KPOINTS_FULL%TRANS(1,NK)) >TINY) .OR. &
                     (ABS(KPOINTS_FULL%TRANS(2,NK))>TINY) .OR. &
                     (ABS(KPOINTS_FULL%TRANS(3,NK))>TINY)) LSHIFT=.TRUE.
                CALL W1_ROTATE_AND_FFT_NO_PROJ(W1, ELEMENT(W, WDES1_IRZ, N, ISP_IRZ), LSHIFT)
#endif
             ENDIF


             CALL PW_CHARGE(WDES1, CHARGE_REAL_SPACE%RG(1,ISP), SIZE(CHARGE_REAL_SPACE%RG,1), W1%CR(1), W1%CR(1), WEIGHT )
          ENDDO band
       ENDDO kpoints
       

    ENDDO spin


    ! fourier-transformation of charge-density using GRID_SOFT
    ! only input data from first in-band-group is used, and all nodes
    ! are involved in the FFT, final result is distributed among nodes
    ! (see SET_RL_GRID() in mgrid.F, and M_divide() in mpi.F)

    ! merge charge from all bands
#ifdef MPI
    IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
       DO I=1,CHARGE_REAL_SPACE%NCDIJ
#ifdef realmode
          CALLMPI( M_sum_d(WDES%COMM_KINTER, CHARGE_REAL_SPACE%RG(1,I),CHARGE_REAL_SPACE%GRID%RL%NP) )
#else
          CALLMPI( M_sum_z(WDES%COMM_KINTER, CHARGE_REAL_SPACE%RG(1,I),CHARGE_REAL_SPACE%GRID%RL%NP) )
#endif
       END DO
    END IF
#endif

    CALL SUMRL_GQ( CHARGE_REAL_SPACE, CHARGE, WDES%COMM_INTER)
    
    CALL FFT_GQ(CHARGE)
      
    ! set the charge-density of unbalanced lattice vectors to 0
    CALL SETUNB_GQ(CHARGE)

    CALL DELWAV(W1, .TRUE.)
    CALL DEALLOCATE_GRID_QUANTITY(CHARGE_REAL_SPACE)

    PROFILING_STOP('soft_charge_sym')

  END SUBROUTINE SOFT_CHARGE_SYM


!*********************** SOFT_CHARGE1 *********************************
!
! subroutine SOFT_CHARGE1 constructs the electronic charge density from
! two set of wavefunctions i.e.
!  rho(r) = sum_n,k (phi(0)_nk(r) phi(1)_nk*(r) + c.c. )* f(0)_nk
!                  + phi(0)_nk(r) phi(0)_nk^*(r) * f(1)_nk
! this is for linear response theory
! f(0) (W0%FERWE) are the initial occupancies
! f(1) (W1%FERWE) are the first order change of the single electron
!                 occupancies
!
!***********************************************************************

      SUBROUTINE SOFT_CHARGE1(GRID,GRID_SOFT,W0,W1,WDES, CHDEN)
      USE prec
      USE mpimy
      USE mgrid
      USE wave
      USE wave_mpi
      IMPLICIT NONE 

      TYPE (wavedes)     WDES
      TYPE (wavespin)    W0
      TYPE (wavespin)    W1
      TYPE (grid_3d)     GRID,GRID_SOFT
      COMPLEX(q)   CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)

    ! work arrays
      RGRID,   ALLOCATABLE :: CDWORK(:,:)
      COMPLEX(q), ALLOCATABLE :: CPTDUM0(:), CPTDUM1(:)
      INTEGER ISPINOR, ISPINOR_, ISP, NK, N, NPL, M, MM, MM_, I
      REAL(q) :: WEIGHT0, WEIGHT1

      ! MPLWV is the allocation in complex words
      ! hence if CDWORK is REAL one needs to double the allocation 
      ALLOCATE(CDWORK(GRID%MPLWV*2,WDES%NCDIJ), &
           CPTDUM0(GRID%MPLWV*WDES%NRSPINORS),CPTDUM1(GRID%MPLWV*WDES%NRSPINORS))

      CDWORK=0
!=======================================================================
! loop over k-points and bands
!=======================================================================
      spin:    DO ISP=1,WDES%ISPIN
      kpoints: DO NK=1,WDES%NKPTS
#ifdef MPI
               IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
      band:    DO N=1,WDES%NBANDS

         WEIGHT0=WDES%RSPIN*WDES%WTKPT(NK)*W0%FERWE(N,NK,ISP)
         WEIGHT1=WDES%RSPIN*WDES%WTKPT(NK)*W1%FERWE(N,NK,ISP)
         IF (WEIGHT0==0 .AND. WEIGHT1==0) CYCLE

         NPL=WDES%NGVECTOR(NK)
!=======================================================================
! fourier-transformation of wave-function
! sum up  real space charge density
!=======================================================================
         DO ISPINOR=0,WDES%NRSPINORS-1
            CALL FFTWAV(NPL,WDES%NINDPW(1,NK),CPTDUM0(1+ISPINOR*GRID%MPLWV),W0%CPTWFP(1+ISPINOR*NPL,N,NK,ISP),GRID)
            CALL FFTWAV(NPL,WDES%NINDPW(1,NK),CPTDUM1(1+ISPINOR*GRID%MPLWV),W1%CPTWFP(1+ISPINOR*NPL,N,NK,ISP),GRID)
         ENDDO
         spinor: DO ISPINOR=0,WDES%NRSPINORS-1 
         DO ISPINOR_=0,WDES%NRSPINORS-1 
            DO M=1,GRID%RL%NP
              MM =M+ISPINOR *GRID%MPLWV
              MM_=M+ISPINOR_*GRID%MPLWV
              CDWORK(M,ISP+ISPINOR_+2*ISPINOR)=CDWORK(M,ISP+ISPINOR_+2*ISPINOR)+ &
                   CPTDUM0(MM)*CONJG(CPTDUM0(MM_))*WEIGHT1+ &
                   CPTDUM0(MM)*CONJG(CPTDUM1(MM_))*WEIGHT0+ &
                   CPTDUM1(MM)*CONJG(CPTDUM0(MM_))*WEIGHT0
           ENDDO
        ENDDO
        ENDDO spinor
      ENDDO band
      ENDDO kpoints
      ENDDO spin

!=======================================================================
! fourier-transformation of charge-density using GRID_SOFT (see above)
!=======================================================================
      DO I=1,WDES%NCDIJ
! now merge the chargedensity from all nodes

!PK Reduce onto first KINTER nodes then broadcast
#ifdef realmode
         CALLMPI( M_sum_d(WDES%COMM_INTER, CDWORK(1,I), GRID%RL%NP))
         CALLMPI( M_sum_d(WDES%COMM_KINTER, CDWORK(1,I), GRID%RL%NP))
#else
         CALLMPI( M_sum_z(WDES%COMM_INTER, CDWORK(1,I), GRID%RL%NP))
         CALLMPI( M_sum_z(WDES%COMM_KINTER, CDWORK(1,I), GRID%RL%NP))
#endif
         CALL FFT_RC_SCALE(CDWORK(1,I),CHDEN(1,I),GRID_SOFT)
! set the charge-density of unbalanced lattic-vectors to 0
         CALL SETUNB(CHDEN(1,I),GRID_SOFT)
      ENDDO

      DEALLOCATE(CDWORK, CPTDUM0, CPTDUM1)

      RETURN
      END SUBROUTINE



!*************************SUBROUTINE RHOAT0 ****************************
!
!  This routine calculates the term G^2 in the Taylor expansion
!  of the atomic chargedensities multiplies this term with
!  the number of atoms per type and the electronic field constant
!
!***********************************************************************

      SUBROUTINE RHOAT0(P,T_INFO, BETATO,OMEGA)
      USE prec
      USE pseudo
      USE poscar
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P (T_INFO%NTYP)

      BETATO=0._q
      DO 200 NT=1,T_INFO%NTYP
        IF (P(NT)%PSPRHO(1)==0) GOTO 200
        IMAX=4
        DQ=P(NT)%PSGMAX/NPSPTS
        B1=0
        B2=0
        DO 100 I=1,IMAX
          A=(P(NT)%PSPRHO(I+1)-P(NT)%PSPRHO(1))/(DQ*I)**2
          B1=B1+A
          B2=B2+A*A
  100   CONTINUE
        B1=B1/IMAX
        B2=B2/IMAX
        BETATO=BETATO+B1*T_INFO%NITYP(NT)*T_INFO%VCA(NT)*EDEPS/OMEGA

  200 CONTINUE
      RETURN
      END SUBROUTINE


!*************************SUBROUTINE RHOATO ****************************
!  This routine calculates the chargedensity and its derivatives
!  corresponding to overlapping  atoms
!  LFOUR
!   .FALSE. set up charge density in reciprocal space
!   .TRUE.  set up charge density in real space
!  LPAR
!   .FALSE. use PSPRHO
!   .TRUE.  use PSPCOR
!***********************************************************************

      SUBROUTINE RHOATO(LFOUR,LPAR,GRIDC,T_INFO,B,P,CSTRF,CHTOT,CHDER)
      USE prec
      USE mpimy
      USE mgrid
      USE poscar
      USE pseudo
      USE constant
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC
      COMPLEX(q)   CHTOT(GRIDC%RC%NP),CHDER(GRIDC%RC%NP)
      COMPLEX(q)   CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
      REAL(q)      B(3,3)
      LOGICAL LFOUR,LPAR
! local variables
      REAL(q), POINTER :: PRHO(:)

      CHTOT=0
      CHDER=0
!=======================================================================
! loop over all types of atoms
! if no pseudocharge for this ion next ion
!=======================================================================
      typ: DO NT=1,T_INFO%NTYP
       IF (LPAR) THEN
         PRHO=>P(NT)%PSPCOR
       ELSE
         PRHO=>P(NT)%PSPRHO
       ENDIF
       IF (.NOT.ASSOCIATED(PRHO)) GOTO 200
!=======================================================================
! calculate the scaling factor ARGSC that converts the magnitude of a
! reciprocal lattice vector to the correponding position in the
! pseudopotential arrays
!=======================================================================
      ARGSC=NPSPTS/P(NT)%PSGMAX
      PSGMA2=P(NT)%PSGMAX-P(NT)%PSGMAX/NPSPTS*3

      DO N=1,GRIDC%RC%NP
        N1= MOD((N-1),GRIDC%RC%NROW) +1
        NC= (N-1)/GRIDC%RC%NROW+1
        N2= GRIDC%RC%I2(NC)
        N3= GRIDC%RC%I3(NC)
!=======================================================================
! calculate the magnitude of the reciprocal lattice vector
!=======================================================================
        GX= GRIDC%LPCTX(N1)*B(1,1)+GRIDC%LPCTY(N2)*B(1,2)+GRIDC%LPCTZ(N3)*B(1,3)
        GY= GRIDC%LPCTX(N1)*B(2,1)+GRIDC%LPCTY(N2)*B(2,2)+GRIDC%LPCTZ(N3)*B(2,3)
        GZ= GRIDC%LPCTX(N1)*B(3,1)+GRIDC%LPCTY(N2)*B(3,2)+GRIDC%LPCTZ(N3)*B(3,3)

        G=SQRT(GX**2+GY**2+GZ**2)*TPI

        IF (G/=0 .AND.G<PSGMA2) THEN
!=======================================================================
! convert the magnitude of the reciprocal latice vector to a position
! in the charge-dens. array  and interpolate the atomic-chargedensity
!=======================================================================
          ARG=(G*ARGSC)+1
          NADDR=MAX(INT(ARG),2)
          REM=ARG-NADDR
          V1=PRHO(NADDR-1)
          V2=PRHO(NADDR)
          V3=PRHO(NADDR+1)
          V4=PRHO(NADDR+2)
          T0=V2
          T1=((6*V3)-(2*V1)-(3*V2)-V4)/6._q
          T2=(V1+V3-(2*V2))/2._q
          T3=(V4-V1+(3*(V2-V3)))/6._q
          CHTOT(N)=CHTOT(N)+(T0+REM*(T1+REM*(T2+REM*T3))) *CSTRF(N,NT)
          CHDER(N)=CHDER(N)+(T1+REM*(2*T2+REM*3*T3))*ARGSC*CSTRF(N,NT)
        ELSE IF (G==0) THEN
          CHTOT(N)=CHTOT(N)+PRHO(1)*CSTRF(N,NT)
          CHDER(N)=0
        ENDIF
      ENDDO
  200 CONTINUE
      ENDDO typ
!=======================================================================
! set the charge-density of unbalanced lattice-vectors to 0
! and transform the charge-density to real space
!=======================================================================
      CALL SETUNB(CHTOT(1),GRIDC)
      CALL SETUNB(CHDER(1),GRIDC)

      IF (LFOUR) THEN
        CALL FFT3D(CHTOT,GRIDC,1)
      ENDIF

      RETURN
      END SUBROUTINE
!
!  small subroutine which allocates the charge derivative array 
!  CWORK explicitly and calls the subroutine above
! 
      SUBROUTINE RHOATO_WORK(LFOUR,LPAR,GRIDC,T_INFO,B,P,CSTRF,CHTOT)
      USE prec
      USE mpimy
      USE mgrid
      USE poscar
      USE pseudo
      IMPLICIT NONE


      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC
      COMPLEX(q)   CHTOT(GRIDC%MPLWV)
      COMPLEX(q)   CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
      REAL(q)      B(3,3)
      LOGICAL LFOUR,LPAR
! work arrays
      COMPLEX(q),ALLOCATABLE :: CWORK(:)

      ALLOCATE(CWORK(1:GRIDC%MPLWV))
      CALL RHOATO(LFOUR,LPAR,GRIDC,T_INFO,B,P,CSTRF,CHTOT,CWORK)
      DEALLOCATE(CWORK)
      
      RETURN
      END SUBROUTINE


!*************************SUBROUTINE RHOATO_PARTICLE_MESH **************
!  This routine calculates the charge density
!  corresponding to overlapping  atoms using the particle mesh method
!  LFOUR
!   .FALSE. set up charge density in reciprocal space
!   .TRUE.  set up charge density in real space
!  LPAR
!   .FALSE. use PSPRHO
!   .TRUE.  use PSPCOR
!
!  NOTE: Does NOT calculate the derivatives!
!***********************************************************************
      SUBROUTINE RHOATO_PARTICLE_MESH(LFOUR,LPAR,GRIDC,LATT_CUR,T_INFO,INFO,P,CHTOT,IU6)
      USE base
      USE prec
      USE mpimy
      USE mgrid
      USE poscar
      USE pseudo
      USE constant
      USE lattice
      IMPLICIT NONE
      LOGICAL LFOUR,LPAR
      TYPE (grid_3d) GRIDC
      TYPE (latt) LATT_CUR
      TYPE (type_info) T_INFO
      TYPE (info_struct) INFO
      TYPE (potcar) P(T_INFO%NTYP)
      COMPLEX(q) CHTOT(GRIDC%RC%NP)
      INTEGER IU6
      ! local variables
      INTEGER NT
      REAL(q) G1,G2,GMAX
      REAL(q), POINTER :: PRHO(:)

      INTEGER, PARAMETER :: NFFT=32768

#ifdef debug
      REAL(q) RHOT
      REAL(q), EXTERNAL :: RHO0
      RHOT=RHO0(GRIDC,CHTOT(1))
      WRITE(*,'(A,F20.10)') 'RHOATO_PARTICLE_MESH: Total number of electrons =',RHOT
#endif

      DO NT=1,T_INFO%NTYP
         NULLIFY(P(NT)%USESPL)
         IF (LPAR) THEN
            PRHO=>P(NT)%PSPCOR
         ELSE
            PRHO=>P(NT)%PSPRHO
         ENDIF
         IF (.NOT.ASSOCIATED(PRHO)) CYCLE
         ! momentum cutoffs for Fourier filtering
         GMAX=SQRT(INFO%ENMAX/HSQDTM)
         IF (INFO%SZPREC(1:1)=='s') THEN
            G1=GMAX
            G2=MIN(3.0_q*GMAX,P(NT)%PSGMAX)/2
         ELSE
            G1=2.0_q*GMAX
            G2=MIN(3.0_q*GMAX,P(NT)%PSGMAX)
         ENDIF
         ! Fourier-filter atomic charge density, transform to real space
         ! and calculate spline interpolation
         ALLOCATE(P(NT)%ATOSPL(NFFT/2,5))
         CALL ATORHO_FILTER(G1,G2,NPSPTS,PRHO,P(NT)%PSGMAX/NPSPTS,NFFT,P(NT)%ATOSPL,P(NT)%RCUTATO,-1)
         ! set the arrays and cutoff that are used by RHOADD
         P(NT)%USESPL=>P(NT)%ATOSPL
         P(NT)%USEZ=PRHO(1)
         P(NT)%USECUT=P(NT)%RCUTATO
      ENDDO

      CHTOT=0
      ! put the atomic charge density of all ions on the grid 
      CALL RHOADD(T_INFO,LATT_CUR,P,GRIDC,CHTOT)

#ifdef debug
      RHOT=RHO0(GRIDC,CHTOT(1))
      WRITE(*,'(A,F20.10)') 'RHOATO_PARTICLE_MESH: Total number of electrons =',RHOT
#endif
      ! if requested, take the charge density to real space
      IF (LFOUR) THEN
         CALL FFT3D(CHTOT,GRIDC,1)
      ENDIF

      DO NT=1,T_INFO%NTYP
         IF (ASSOCIATED(P(NT)%ATOSPL)) THEN
            DEALLOCATE(P(NT)%ATOSPL)
            NULLIFY(P(NT)%ATOSPL)
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE RHOATO_PARTICLE_MESH


!******************** SUBROUTINE RHOADD ********************************
!
!***********************************************************************
      SUBROUTINE RHOADD(T_INFO,LATT_CUR,P,GRIDC,CHTOT,QTOT)
      USE ini
      USE prec
      USE poscar
      USE lattice
      USE mgrid
      USE pseudo
      USE constant
      IMPLICIT NONE
      TYPE (type_info) T_INFO
      TYPE (latt) LATT_CUR
      TYPE (potcar) P(T_INFO%NTYP)
      TYPE (grid_3d) GRIDC
      COMPLEX(q) CHTOT(GRIDC%MPLWV)
      REAL(q), OPTIONAL :: QTOT(T_INFO%NIONS)
      ! local variables
      REAL(q) QTOT_(T_INFO%NIONS)

#ifdef debug
      REAL(q) RHOT
      REAL(q), EXTERNAL :: RHO0
      RHOT=RHO0(GRIDC,CHTOT(1))
      WRITE(*,'(A,F20.10)') 'Total number of electrons =',RHOT
#endif

      ! take total charge to real space
      CALL FFT3D(CHTOT(1),GRIDC,1)
      
      CALL RHOADD_RL(T_INFO,LATT_CUR,P,GRIDC,QTOT_,CHTOT(1))

      ! and back to reciprocal space
      CALL FFT_RC_SCALE(CHTOT(1),CHTOT(1),GRIDC)
      CALL SETUNB_COMPAT(CHTOT(1),GRIDC)

#ifdef debug
      RHOT=RHO0(GRIDC,CHTOT(1))
      WRITE(*,'(A,F20.10)') 'Total number of electrons =',RHOT
#endif
      IF (PRESENT(QTOT)) QTOT=QTOT_

      RETURN
      END SUBROUTINE RHOADD


!*************************SUBROUTINE RHOATR ****************************
!  subroutine to calculate the charge-density (or partial-core
!  charge-densities) from overlapping atoms in real space
!  and store it in an real array
!  NOTE: the DENCOR array is not necessarily suffiently large
!  to allow for 3d-FFT, hence a work array is temporarily allocated
!***********************************************************************

      SUBROUTINE RHOPAR(GRIDC,T_INFO,INFO,LATT_CUR,P,CSTRF,DENS,IU6)
      USE prec
      USE base
      USE mpimy
      USE mgrid
      USE poscar
      USE pseudo
      USE lattice
      IMPLICIT NONE
      TYPE (type_info) T_INFO
      TYPE (info_struct) INFO
      TYPE (potcar) P(T_INFO%NTYP)
      TYPE (grid_3d) GRIDC
      TYPE (latt) LATT_CUR
      RGRID DENS(GRIDC%RL%NP)
      COMPLEX(q) CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
      INTEGER IU6
! work arrays
      COMPLEX(q),ALLOCATABLE:: CWORK1(:),CWORK2(:)

      ALLOCATE(CWORK1(GRIDC%MPLWV),CWORK2(GRIDC%MPLWV))
      IF(INFO%TURBO==0)THEN
         CALL RHOATO(.TRUE.,.TRUE.,GRIDC,T_INFO,LATT_CUR%B,P,CSTRF,CWORK1,CWORK2)
      ELSE
         CALL RHOATO_PARTICLE_MESH(.TRUE.,.TRUE.,GRIDC,LATT_CUR,T_INFO,INFO,P,CWORK1,IU6)
      END IF
      CALL RL_ADD(CWORK1,1.0_q,CWORK1,0.0_q,DENS,GRIDC)

      DEALLOCATE(CWORK1,CWORK2)

      RETURN
      END SUBROUTINE
      END MODULE


!************************ SUBROUTINE MRHOATO ***************************
!
!  This routine calculates the magnetization density and its derivatives
!  corresponding to overlapping  atoms
!  LFOUR  set up charge density in reciprocal space and store in CHTOT
!
!***********************************************************************
      SUBROUTINE MRHOATO(LFOUR,GRIDC,T_INFO,B,P, CHTOT, ISPIN)
      USE prec
      USE mpimy
      USE mgrid
      USE poscar
      USE pseudo
      USE constant
      IMPLICIT NONE

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC
      INTEGER ISPIN
      COMPLEX(q)   CHTOT(GRIDC%MPLWV,ISPIN)
      REAL(q)      B(3,3)
      LOGICAL   LFOUR
      INTEGER ISP
! local variables
      INTEGER NT,NADDR,NI,N,N3,N2,N1,NA(T_INFO%NIONS),NIS,NC
      REAL(q) V1,V2,V3,V4,T0,T1,T2,T3,REM,ARGSC,PSGMA2, &
              GX,GY,GZ,G,ARG,AMAG,FNORM
! work arrays
      COMPLEX(q),ALLOCATABLE:: CSTRF(:)

      IF (ISPIN==0) RETURN

      ALLOCATE(CSTRF(GRIDC%RC%NP))

      CHTOT=0
!=======================================================================
! loop over all types of atoms
! if no pseudocharge for this ion next ion
!=======================================================================
      NIS=1
      typ: DO NT=1,T_INFO%NTYP
       IF (P(NT)%PSPRHO(1)==0) GOTO 200

      ARGSC=NPSPTS/P(NT)%PSGMAX
      PSGMA2=P(NT)%PSGMAX-P(NT)%PSGMAX/NPSPTS*3
      FNORM=1._q/P(NT)%ZVALF

 ion: DO NI=NIS,T_INFO%NITYP(NT)+NIS-1

    ! set up phase factor for this ion
      CALL STUFAK_ONE(GRIDC,1,T_INFO%POSION(1,NI),T_INFO%VCA(NT),CSTRF)

    ! AMAG is the magnitude of the  magnetization for this ion
    ! and for the present direction
spin: DO ISP=1,ISPIN

      AMAG=T_INFO%ATOMOM(ISP+(NI-1)*ISPIN)*FNORM

      DO N=1,GRIDC%RC%NP
        N1= MOD((N-1),GRIDC%RC%NROW) +1
        NC= (N-1)/GRIDC%RC%NROW+1
        N2= GRIDC%RC%I2(NC)
        N3= GRIDC%RC%I3(NC)
        ! calculate the magnitude of the reciprocal lattice vector
        GX= GRIDC%LPCTX(N1)*B(1,1)+GRIDC%LPCTY(N2)*B(1,2)+GRIDC%LPCTZ(N3)*B(1,3)
        GY= GRIDC%LPCTX(N1)*B(2,1)+GRIDC%LPCTY(N2)*B(2,2)+GRIDC%LPCTZ(N3)*B(2,3)
        GZ= GRIDC%LPCTX(N1)*B(3,1)+GRIDC%LPCTY(N2)*B(3,2)+GRIDC%LPCTZ(N3)*B(3,3)

        G=SQRT(GX**2+GY**2+GZ**2)*TPI

        IF (G/=0 .AND.G<PSGMA2) THEN
        ! convert the magnitude of the reciprocal latice vector to a position
        ! in the charge-dens. array  and interpolate the atomic-chargedensity
          ARG=(G*ARGSC)+1
          NADDR=MAX(INT(ARG),2)
          REM=ARG-NADDR
          V1=P(NT)%PSPRHO(NADDR-1)
          V2=P(NT)%PSPRHO(NADDR)
          V3=P(NT)%PSPRHO(NADDR+1)
          V4=P(NT)%PSPRHO(NADDR+2)
          T0=V2
          T1=((6*V3)-(2*V1)-(3*V2)-V4)/6._q
          T2=(V1+V3-(2*V2))/2._q
          T3=(V4-V1+(3*(V2-V3)))/6._q
          ! voila: phase factor  * charge * AMAG
          CHTOT(N,ISP)=CHTOT(N,ISP)+(T0+REM*(T1+REM*(T2+REM*T3))) *CSTRF(N)*AMAG
        ELSE IF (G==0) THEN
          CHTOT(N,ISP)=CHTOT(N,ISP)+P(NT)%PSPRHO(1)*CSTRF(N)*AMAG
        ENDIF
      ENDDO
      ENDDO spin
      ENDDO ion
  200 NIS=NIS+T_INFO%NITYP(NT)
!=======================================================================
      ENDDO typ
!=======================================================================

! set the charge-density of unbalanced lattic-vectors to 0
! and transform the charge-density to real space

      DO ISP=1,ISPIN
         CALL SETUNB(CHTOT(1,ISP),GRIDC)
         IF (LFOUR) THEN
            CALL FFT3D(CHTOT(1,ISP),GRIDC,1)
         ENDIF
      ENDDO

      DEALLOCATE( CSTRF )

      RETURN
      END SUBROUTINE


!******************** SUBROUTINE RHOADD_RL *****************************
! Evaluate the atomic or pseudo or ... charge density IN P(:)%USESPL
! on the grid GRIDC,
! cutting of at a cube of length 2*P(:)%USECUT, and
! normalizing to P(:)%USEZ for each ion
!***********************************************************************
      SUBROUTINE RHOADD_RL(T_INFO,LATT_CUR,P,GRIDC,QITOT,DHTOT)
      USE prec
      USE constant
      USE poscar
      USE lattice
      USE mgrid
      USE mpimy
      USE pseudo
      USE radial
      IMPLICIT NONE
      TYPE (type_info) T_INFO
      TYPE (latt) LATT_CUR
      TYPE (potcar) P(T_INFO%NTYP)
      TYPE (grid_3d) GRIDC
      REAL(q) QITOT(T_INFO%NIONS)
      RGRID DHTOT(DIMREAL(GRIDC%MPLWV))
      ! local variables
      INTEGER I
      INTEGER NIS,NT,NI,K
      INTEGER N1LOW,N2LOW,N3LOW,N1HI,N2HI,N3HI
      INTEGER N1,N2,N3,N1P,N2P,N3P,NCOL,IND,IND1
      REAL(q) F1,F2,F3,D1,D2,D3,X1,X2,X3
      REAL(q) XC,YC,ZC,D,DX
      REAL(q) QR,BJ,QSUM,DQ,QN,RHOT,QSM1
      INTEGER NTOT(T_INFO%NIONS)
      REAL(q) QTOT(T_INFO%NIONS)
      REAL(q) QSPH(T_INFO%NIONS)
      INTEGER, PARAMETER :: NQ=2
      REAL(q) QQ(NQ),A(NQ),VAL,VALDER,RCUT
      REAL(q), PARAMETER :: TINY = 1.E-5_q

      NTOT=0; QTOT=0; QSPH=0

      NIS=1

      type1: DO NT=1,T_INFO%NTYP
         IF (.NOT.ASSOCIATED(P(NT)%USESPL)) THEN
            NIS = NIS+T_INFO%NITYP(NT); CYCLE type1
         ENDIF

         ions1: DO NI=NIS,T_INFO%NITYP(NT)+NIS-1

            !=======================================================================
            ! find lattice points contained within the cutoff-sphere
            !=======================================================================
            F1=1._q/GRIDC%NGX
            F2=1._q/GRIDC%NGY
            F3=1._q/GRIDC%NGZ

            !-----------------------------------------------------------------------
            ! restrict loop to points contained within a cubus around the ion
            !-----------------------------------------------------------------------
            D1= P(NT)%USECUT*LATT_CUR%BNORM(1)*GRIDC%NGX
            D2= P(NT)%USECUT*LATT_CUR%BNORM(2)*GRIDC%NGY
            D3= P(NT)%USECUT*LATT_CUR%BNORM(3)*GRIDC%NGZ

            N3LOW= INT(T_INFO%POSION(3,NI)*GRIDC%NGZ-D3+10*GRIDC%NGZ+.99_q)-10*GRIDC%NGZ
            N2LOW= INT(T_INFO%POSION(2,NI)*GRIDC%NGY-D2+10*GRIDC%NGY+.99_q)-10*GRIDC%NGY
            N1LOW= INT(T_INFO%POSION(1,NI)*GRIDC%NGX-D1+10*GRIDC%NGX+.99_q)-10*GRIDC%NGX
            
            N3HI = INT(T_INFO%POSION(3,NI)*GRIDC%NGZ+D3+10*GRIDC%NGZ)-10*GRIDC%NGZ
            N2HI = INT(T_INFO%POSION(2,NI)*GRIDC%NGY+D2+10*GRIDC%NGY)-10*GRIDC%NGY
            N1HI = INT(T_INFO%POSION(1,NI)*GRIDC%NGX+D1+10*GRIDC%NGX)-10*GRIDC%NGX
            
            IND=0; QSUM=0; QSM1=0

#ifdef MPI
            !-----------------------------------------------------------------------
            ! loop over cubus
            ! MPI version z ist the fast index
            !-----------------------------------------------------------------------
            DO N2=N2LOW,N2HI
               X2=(N2*F2-T_INFO%POSION(2,NI))
               N2P=MOD(N2+10*GRIDC%NGY,GRIDC%NGY)
               
               DO N1=N1LOW,N1HI
                  X1=(N1*F1-T_INFO%POSION(1,NI))
                  N1P=MOD(N1+10*GRIDC%NGX,GRIDC%NGX)
                  
                  NCOL=GRIDC%RL%INDEX(N1P,N2P)
                  IF (NCOL==0) CYCLE ! not on local node, move on
                  IF (GRIDC%RL%I2(NCOL) /= N1P+1 .OR. GRIDC%RL%I3(NCOL) /= N2P+1) THEN
                     WRITE(*,*)'RHOADD: internal ERROR:', &
                          GRIDC%RL%I2(NCOL),N1P+1, GRIDC%RL%I3(NCOL),N2P+1
                     STOP
                  ENDIF
                  !OCL SCALAR
                  DO N3=N3LOW,N3HI
                     X3=(N3*F3-T_INFO%POSION(3,NI))
                     N3P=MOD(N3+10*GRIDC%NGZ,GRIDC%NGZ)
                     IND1=N3P+(NCOL-1)*GRIDC%NGZ+1          
#else
            !-----------------------------------------------------------------------
            ! loop over cubus around one ion
            ! conventional version x is fast index
            !-----------------------------------------------------------------------
            DO N3=N3LOW,N3HI
               X3=(N3*F3-T_INFO%POSION(3,NI))
               N3P=MOD(N3+10*GRIDC%NGZ,GRIDC%NGZ)
               
               DO N2=N2LOW,N2HI
                  X2=(N2*F2-T_INFO%POSION(2,NI))
                  N2P=MOD(N2+10*GRIDC%NGY,GRIDC%NGY)
                  
                  NCOL=GRIDC%RL%INDEX(N2P,N3P)
                  
                  DO N1=N1LOW,N1HI
                     X1=(N1*F1-T_INFO%POSION(1,NI))
                     N1P=MOD(N1+10*GRIDC%NGX,GRIDC%NGX)
                     IND1=N1P+(NCOL-1)*GRIDC%NGX+1
#endif
            
                     XC= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
                     YC= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
                     ZC= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)
                     
                     D=SQRT(XC*XC+YC*YC+ZC*ZC)
                     
                     IF (D<=P(NT)%USECUT) THEN
                        IND=IND+1
                                                 
                        K=D/P(NT)%USESPL(2,1)+1
                        DX=D-P(NT)%USESPL(K,1)
                        VAL   =((P(NT)%USESPL(K,5)*DX+P(NT)%USESPL(K,4))*DX+P(NT)%USESPL(K,3))*DX+P(NT)%USESPL(K,2)
                        VALDER=(3.0_q*P(NT)%USESPL(K,5)*DX+2.0_q*P(NT)%USESPL(K,4))*DX+P(NT)%USESPL(K,3)

                        QSUM=QSUM+VAL
            
                        ! sum total charge in sphere
                        QSM1=QSM1+DHTOT(IND1)
                     ENDIF
                     
                  ENDDO
               ENDDO
            ENDDO
                     
            QTOT(NI)=QSUM
            NTOT(NI)=IND
            QSPH(NI)=QSM1

         ENDDO ions1
         NIS = NIS+T_INFO%NITYP(NT)
      ENDDO type1

      CALLMPI( M_sum_i(GRIDC%COMM, NTOT, T_INFO%NIONS))
      CALLMPI( M_sum_d(GRIDC%COMM, QTOT, T_INFO%NIONS))
      CALLMPI( M_sum_d(GRIDC%COMM, QSPH, T_INFO%NIONS))
      QITOT=QTOT

#ifdef debug
      DO NI=1,T_INFO%NIONS
         ! WRITE(*,'(I6,I8,F14.7)') NI,NTOT(NI),QTOT(NI)/NTOT(NI)!*LATT_CUR%OMEGA*F1*F2*F3
         WRITE(6,'(a,I6,10F14.7)') 'tot',NI,QTOT(NI)*LATT_CUR%OMEGA*F1*F2*F3,QSPH(NI)/GRIDC%NPLWV
      ENDDO
#endif

      NIS=1
      type2: DO NT=1,T_INFO%NTYP

         IF (.NOT.ASSOCIATED(P(NT)%USESPL)) THEN
            NIS = NIS+T_INFO%NITYP(NT); CYCLE type2
         ENDIF

         ions2: DO NI=NIS,T_INFO%NITYP(NT)+NIS-1

            !-----------------------------------------------------------------------
            ! restrict loop to points contained within a cubus around the ion
            !-----------------------------------------------------------------------
            D1= P(NT)%USECUT*LATT_CUR%BNORM(1)*GRIDC%NGX
            D2= P(NT)%USECUT*LATT_CUR%BNORM(2)*GRIDC%NGY
            D3= P(NT)%USECUT*LATT_CUR%BNORM(3)*GRIDC%NGZ
            
            N3LOW= INT(T_INFO%POSION(3,NI)*GRIDC%NGZ-D3+10*GRIDC%NGZ+.99_q)-10*GRIDC%NGZ
            N2LOW= INT(T_INFO%POSION(2,NI)*GRIDC%NGY-D2+10*GRIDC%NGY+.99_q)-10*GRIDC%NGY
            N1LOW= INT(T_INFO%POSION(1,NI)*GRIDC%NGX-D1+10*GRIDC%NGX+.99_q)-10*GRIDC%NGX
            
            N3HI = INT(T_INFO%POSION(3,NI)*GRIDC%NGZ+D3+10*GRIDC%NGZ)-10*GRIDC%NGZ
            N2HI = INT(T_INFO%POSION(2,NI)*GRIDC%NGY+D2+10*GRIDC%NGY)-10*GRIDC%NGY
            N1HI = INT(T_INFO%POSION(1,NI)*GRIDC%NGX+D1+10*GRIDC%NGX)-10*GRIDC%NGX
            
            IND=0; QSUM=0
            !-----------------------------------------------------------------------
            ! loop over cubus
            ! MPI version z ist the fast index
            !-----------------------------------------------------------------------
#ifdef MPI
            DO N2=N2LOW,N2HI
               X2=(N2*F2-T_INFO%POSION(2,NI))
               N2P=MOD(N2+10*GRIDC%NGY,GRIDC%NGY)
               
               DO N1=N1LOW,N1HI
                  X1=(N1*F1-T_INFO%POSION(1,NI))
                  N1P=MOD(N1+10*GRIDC%NGX,GRIDC%NGX)
                  
                  NCOL=GRIDC%RL%INDEX(N1P,N2P)
                  IF (NCOL==0) CYCLE ! not on local node, move on
                  IF (GRIDC%RL%I2(NCOL) /= N1P+1 .OR. GRIDC%RL%I3(NCOL) /= N2P+1) THEN
                     WRITE(*,*)'RHOADD: internal ERROR:', &
                          GRIDC%RL%I2(NCOL),N1P+1, GRIDC%RL%I3(NCOL),N2P+1
                     STOP
                  ENDIF
                  !OCL SCALAR
                  DO N3=N3LOW,N3HI
                     X3=(N3*F3-T_INFO%POSION(3,NI))
                     N3P=MOD(N3+10*GRIDC%NGZ,GRIDC%NGZ)
                     IND=N3P+(NCOL-1)*GRIDC%NGZ+1
#else
            !-----------------------------------------------------------------------
            ! loop over cubus around one ion
            ! conventional version x is fast index
            !-----------------------------------------------------------------------
            DO N3=N3LOW,N3HI
               X3=(N3*F3-T_INFO%POSION(3,NI))
               N3P=MOD(N3+10*GRIDC%NGZ,GRIDC%NGZ)
               
               DO N2=N2LOW,N2HI
                  X2=(N2*F2-T_INFO%POSION(2,NI))
                  N2P=MOD(N2+10*GRIDC%NGY,GRIDC%NGY)
                  
                  NCOL=GRIDC%RL%INDEX(N2P,N3P)
                  
                  DO N1=N1LOW,N1HI
                     X1=(N1*F1-T_INFO%POSION(1,NI))
                     N1P=MOD(N1+10*GRIDC%NGX,GRIDC%NGX)
                     IND=N1P+(NCOL-1)*GRIDC%NGX+1
#endif        
             
                     XC= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
                     YC= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
                     ZC= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)
                     
                     D=SQRT(XC*XC+YC*YC+ZC*ZC)
                     
                     IF (D<=P(NT)%USECUT) THEN

                        K=D/P(NT)%USESPL(2,1)+1
                        DX=D-P(NT)%USESPL(K,1)
                        VAL   =((P(NT)%USESPL(K,5)*DX+P(NT)%USESPL(K,4))*DX+P(NT)%USESPL(K,3))*DX+P(NT)%USESPL(K,2)
                        VALDER=(3.0_q*P(NT)%USESPL(K,5)*DX+2.0_q*P(NT)%USESPL(K,4))*DX+P(NT)%USESPL(K,3)

                        QSUM=QSUM+VAL/QTOT(NI)
                        DHTOT(IND)=DHTOT(IND)+VAL/QTOT(NI)*P(NT)%USEZ*GRIDC%NPLWV
                     ENDIF
                     
                  ENDDO
               ENDDO
            ENDDO
                                 
            QTOT(NI)=QSUM
            
         ENDDO ions2
         NIS = NIS+T_INFO%NITYP(NT)
      ENDDO type2
      
      CALLMPI( M_sum_i(GRIDC%COMM, NTOT, T_INFO%NIONS))
      CALLMPI( M_sum_d(GRIDC%COMM, QTOT, T_INFO%NIONS))

#ifdef debug
      RHOT=0
      DO IND=1,GRIDC%RL%NP
         RHOT=RHOT+DHTOT(IND)
      ENDDO
      CALLMPI( M_sum_d(GRIDC%COMM, RHOT, 1))
      
      DO NI=1,T_INFO%NIONS
         WRITE(*,'(I6,I8,F14.7)') NI,NTOT(NI),QTOT(NI)!*LATT_CUR%OMEGA*F1*F2*F3
      ENDDO
      WRITE(*,'(A,F20.10)') 'Total number of electrons =',RHOT/GRIDC%NPLWV
#endif

      RETURN
      END SUBROUTINE RHOADD_RL


!******************** SUBROUTINE RHODER *****************************
! Calculates the derivative of the charge density with respect to
! one ionic coordinate r_NI:
!        f_i(r)=d rho/dr*(r-R_i)/|r-R_i|
! (rho .. P%USESPL)
! and integrate
!        F_i = \int dr f_i(r) V(r)
! (V .. POTTOT, F_i .. CIFOR)
! to get the force acting on ion NI
!********************************************************************
      SUBROUTINE RHODER(T_INFO,LATT_CUR,P,NI,GRIDC,QTOT,POTTOT,IFOR)
      USE prec
      USE constant
      USE poscar
      USE lattice
      USE mgrid
      USE mpimy
      USE pseudo
      USE radial
      IMPLICIT NONE
      TYPE (type_info) T_INFO
      TYPE (latt) LATT_CUR
      TYPE (potcar) P
      INTEGER NI
      TYPE (grid_3d) GRIDC
      REAL(q) QTOT(T_INFO%NIONS)
      RGRID POTTOT(DIMREAL(GRIDC%MPLWV))
      REAL(q) IFOR(3)
      ! local variables
      INTEGER I,K
      INTEGER N1LOW,N2LOW,N3LOW,N1HI,N2HI,N3HI
      INTEGER N1,N2,N3,N1P,N2P,N3P,NCOL,IND
      REAL(q) F1,F2,F3,D1,D2,D3,X1,X2,X3
      REAL(q) XC,YC,ZC,D,DX
      REAL(q) VAL,VALDER,RCUT
      REAL(q) :: DHDER(3)
      COMPLEX(q) CIFOR(3)

      F1=1._q/GRIDC%NGX
      F2=1._q/GRIDC%NGY
      F3=1._q/GRIDC%NGZ

      !-----------------------------------------------------------------------
      ! restrict loop to points contained within a cubus around the ion
      !-----------------------------------------------------------------------
      D1= P%USECUT*LATT_CUR%BNORM(1)*GRIDC%NGX
      D2= P%USECUT*LATT_CUR%BNORM(2)*GRIDC%NGY
      D3= P%USECUT*LATT_CUR%BNORM(3)*GRIDC%NGZ

      N3LOW= INT(T_INFO%POSION(3,NI)*GRIDC%NGZ-D3+10*GRIDC%NGZ+.99_q)-10*GRIDC%NGZ
      N2LOW= INT(T_INFO%POSION(2,NI)*GRIDC%NGY-D2+10*GRIDC%NGY+.99_q)-10*GRIDC%NGY
      N1LOW= INT(T_INFO%POSION(1,NI)*GRIDC%NGX-D1+10*GRIDC%NGX+.99_q)-10*GRIDC%NGX

      N3HI = INT(T_INFO%POSION(3,NI)*GRIDC%NGZ+D3+10*GRIDC%NGZ)-10*GRIDC%NGZ
      N2HI = INT(T_INFO%POSION(2,NI)*GRIDC%NGY+D2+10*GRIDC%NGY)-10*GRIDC%NGY
      N1HI = INT(T_INFO%POSION(1,NI)*GRIDC%NGX+D1+10*GRIDC%NGX)-10*GRIDC%NGX

      IND=0
      DHDER=0
      CIFOR=0
      !-----------------------------------------------------------------------
      ! loop over cubus
      ! MPI version z ist the fast index
      !-----------------------------------------------------------------------
#ifdef MPI
      DO N2=N2LOW,N2HI
         X2=(N2*F2-T_INFO%POSION(2,NI))
         N2P=MOD(N2+10*GRIDC%NGY,GRIDC%NGY)

         DO N1=N1LOW,N1HI
            X1=(N1*F1-T_INFO%POSION(1,NI))
            N1P=MOD(N1+10*GRIDC%NGX,GRIDC%NGX)

            NCOL=GRIDC%RL%INDEX(N1P,N2P)
            IF (NCOL==0) CYCLE ! not on local node, move on
            IF (GRIDC%RL%I2(NCOL) /= N1P+1 .OR. GRIDC%RL%I3(NCOL) /= N2P+1) THEN
               WRITE(*,*)'RHOADD: internal ERROR:', &
                    GRIDC%RL%I2(NCOL),N1P+1, GRIDC%RL%I3(NCOL),N2P+1
               STOP
            ENDIF
            !OCL SCALAR
            DO N3=N3LOW,N3HI
               X3=(N3*F3-T_INFO%POSION(3,NI))
               N3P=MOD(N3+10*GRIDC%NGZ,GRIDC%NGZ)
               IND=N3P+(NCOL-1)*GRIDC%NGZ+1
#else
      !-----------------------------------------------------------------------
      ! loop over cubus around one ion
      ! conventional version x is fast index
      !-----------------------------------------------------------------------
      DO N3=N3LOW,N3HI
         X3=(N3*F3-T_INFO%POSION(3,NI))
         N3P=MOD(N3+10*GRIDC%NGZ,GRIDC%NGZ)
         
         DO N2=N2LOW,N2HI
            X2=(N2*F2-T_INFO%POSION(2,NI))
            N2P=MOD(N2+10*GRIDC%NGY,GRIDC%NGY)
            
            NCOL=GRIDC%RL%INDEX(N2P,N3P)
            
            DO N1=N1LOW,N1HI
               X1=(N1*F1-T_INFO%POSION(1,NI))
               N1P=MOD(N1+10*GRIDC%NGX,GRIDC%NGX)
               IND=N1P+(NCOL-1)*GRIDC%NGX+1
#endif        
               
               XC= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
               YC= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
               ZC= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)
               
               D=SQRT(XC*XC+YC*YC+ZC*ZC)
               
               IF (D<=P%USECUT.AND.D>0) THEN
                  ! derivative DHDER(a) = d rho/d R_NI,a

                  K=D/P%USESPL(2,1)+1
                  DX=D-P%USESPL(K,1)
                  VALDER=(3.0_q*P%USESPL(K,5)*DX+2.0_q*P%USESPL(K,4))*DX+P%USESPL(K,3)

                  VAL=VALDER/QTOT(NI)*P%USEZ*GRIDC%NPLWV/D
                  DHDER(1)=VAL*XC
                  DHDER(2)=VAL*YC
                  DHDER(3)=VAL*ZC
                  ! force CIFOR(a) = \int d rho/d R_NI,a (r) V(r)
                  CIFOR=CIFOR+DHDER*POTTOT(IND)
               ENDIF
               
            ENDDO
         ENDDO
      ENDDO

      CALLMPI( M_sum_z(GRIDC%COMM, CIFOR, 3) )
      CIFOR=CIFOR/LATT_CUR%OMEGA*F1*F2*F3

      IFOR=REAL(CIFOR,kind=q)

      RETURN
      END SUBROUTINE RHODER


!*************************SUBROUTINE RHO0  *****************************
! this subroutine calculates the total number of electrons in recip space
! (not quite trival in parallel mode)
!***********************************************************************

      FUNCTION RHO0(GRID, CHTOT)
      USE prec
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d) GRID
      COMPLEX(q) CHTOT(GRID%RC%NROW,GRID%RC%NCOL)

      N2= 1
      N3= 1
      N1= 1

      RHO_SUM=0
      DO NC=1,GRID%RC%NCOL
        IF (GRID%RC%I2(NC)==N2 .AND. GRID%RC%I3(NC)==N3) THEN
          RHO_SUM=CHTOT(NC,N1)
        ENDIF
      ENDDO

      CALLMPI( M_sum_d( GRID%COMM, RHO_SUM, 1))
      RHO0= RHO_SUM
      RETURN
      END FUNCTION


      FUNCTION CRHO0(GRID, CHTOT)
      USE prec
      USE mgrid
      IMPLICIT NONE

      TYPE (grid_3d) GRID
      COMPLEX(q) CHTOT(GRID%RC%NROW,GRID%RC%NCOL)
      COMPLEX(q) CRHO0
    ! local
      INTEGER N1, N2, N3, NC
      COMPLEX(q) CRHO_SUM

      N2= 1
      N3= 1
      N1= 1

      CRHO_SUM=0
      DO NC=1,GRID%RC%NCOL
        IF (GRID%RC%I2(NC)==N2 .AND. GRID%RC%I3(NC)==N3) THEN
          CRHO_SUM=CHTOT(NC,N1)
        ENDIF
      ENDDO

      CALLMPI( M_sum_z( GRID%COMM, CRHO_SUM, 1))
      CRHO0= CRHO_SUM
      RETURN
      END FUNCTION


!*************************SUBROUTINE RHO0  *****************************
! this subroutine calculates the total number of electrons in recip space
! (not quite trival in parallel mode)
!***********************************************************************

      SUBROUTINE GET_RHO0(GRID, CHTOT, RHO0)
      USE prec
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d) GRID
      COMPLEX(q) CHTOT(GRID%RC%NROW,GRID%RC%NCOL)
      N2= 1
      N3= 1
      N1= 1
      RHO=0
      DO NC=1,GRID%RC%NCOL
        IF (GRID%RC%I2(NC)==N2 .AND. GRID%RC%I3(NC)==N3) THEN
          RHO=CHTOT(NC,N1)
        ENDIF
      ENDDO

      CALLMPI( M_sum_d( GRID%COMM, RHO, 1))
      RHO0= RHO
      RETURN
      END SUBROUTINE

!*************************SUBROUTINE SET_RHO0  *************************
! this subroutine sets the  total number of electrons in recip space
! (not quite trival in parallel mode)
!***********************************************************************

      SUBROUTINE SET_RHO0(GRID, CHTOT, RHO_SOLL)
      USE prec
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d) GRID
      COMPLEX(q) CHTOT(GRID%RC%NROW,GRID%RC%NCOL)

      N2= 1
      N3= 1
      N1= 1

      RHO_SUM=0
      DO NC=1,GRID%RC%NCOL
        IF (GRID%RC%I2(NC)==N2 .AND. GRID%RC%I3(NC)==N3) THEN
          CHTOT(NC,N1)=RHO_SOLL
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE


!*************************SUBROUTINE TRUNC_HIGH_FREQU  *****************
!
! this subroutine truncates the high frequency components
! of an array by imposing a spherical cutoff
! the spherical cutoff is chosen in such a way that the entire sphere
! fits into the parallelepided spanned by the plane wave basis set
! this subroutine is required for GGA calculations to maintain the
! full symmetry of the exchange correlation potential
!
!***********************************************************************

      MODULE compat_gga
        LOGICAL GGA_COMPAT
      END MODULE compat_gga

      SUBROUTINE TRUNC_HIGH_FREQU(LATT_CUR, GRID, C)
      USE prec
      USE mgrid
      USE lattice
      USE compat_gga
      IMPLICIT NONE

      TYPE (grid_3d) GRID
      TYPE (latt)    LATT_CUR
      COMPLEX(q)     C(GRID%MPLWV)
      REAL(q) GX, GY, GZ, G2, GMIN2
      INTEGER I, N1, N2, N3, NC, NZERO
!
! if you want to revert to the behaviour of vasp before vasp.4.6.16
! comment in the RETURN statment below
!
      IF (GGA_COMPAT) RETURN

      GX=GRID%NGX/2/LATT_CUR%ANORM(1)
      GY=GRID%NGY/2/LATT_CUR%ANORM(2)
      GZ=GRID%NGZ/2/LATT_CUR%ANORM(3)

      GMIN2=MIN(GX,GY,GZ)**2
      NZERO=0

      DO I=1,GRID%RC%NP
         N1= MOD((I-1),GRID%RC%NROW) +1
         NC= (I-1)/GRID%RC%NROW+1
         N2= GRID%RC%I2(NC)
         N3= GRID%RC%I3(NC)
         GX=(GRID%LPCTX(N1)*LATT_CUR%B(1,1)+GRID%LPCTY(N2)*LATT_CUR%B(1,2)+GRID%LPCTZ(N3)*LATT_CUR%B(1,3))
         GY=(GRID%LPCTX(N1)*LATT_CUR%B(2,1)+GRID%LPCTY(N2)*LATT_CUR%B(2,2)+GRID%LPCTZ(N3)*LATT_CUR%B(2,3))
         GZ=(GRID%LPCTX(N1)*LATT_CUR%B(3,1)+GRID%LPCTY(N2)*LATT_CUR%B(3,2)+GRID%LPCTZ(N3)*LATT_CUR%B(3,3))
         G2=GX*GX+GY*GY+GZ*GZ

         IF (G2>GMIN2) THEN
            NZERO=NZERO+1
            C(I)=0
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE





!****************** subroutine GGA_COMPAT_MODE ***********************
!
! If GGA_COMPAT is .TRUE. the vasp.4.4-4.6 behavior is used whereas
! for GGA_COMPAT .FALSE. the new corrected version is used
!
! GGA_COMPAT = .TRUE.   is the default for vasp.4.6
! GGA_COMPAT = .TRUE.   is the default for vasp.5.0
!
!***********************************************************************

      SUBROUTINE GGA_COMPAT_MODE(IU5, IU0, LCOMPAT)
        USE prec
        USE base
        USE compat_gga
        IMPLICIT NONE 
        
        INTEGER               :: IU5,IU0    ! input unit
        LOGICAL :: LCOMPAT

        LOGICAL :: LOPEN,LDUM
        INTEGER :: IDUM, N, IERR
        REAL(q) :: RDUM
        COMPLEX(q)  :: CDUM
        CHARACTER (LEN=1) :: CHARAC
        
        LOPEN=.FALSE.
        OPEN(UNIT=IU5,FILE=INCAR,STATUS='OLD')

        GGA_COMPAT=.TRUE.
        CALL RDATAB(LOPEN,INCAR,IU5,'GGA_COMPAT','=','#',';','L', &
             &  IDUM,RDUM,CDUM,GGA_COMPAT,CHARAC,N,1,IERR)

        IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
           IF (IU0>=0) THEN
              WRITE(IU0,*)'Error reading item ''GGA_COMPAT'' from file INCAR.'
           ENDIF
           GGA_COMPAT=.FALSE.
        ENDIF

        CALL XML_INCAR('GGA_COMPAT','L',IDUM,RDUM,CDUM,GGA_COMPAT,CHARAC,N)

      END SUBROUTINE GGA_COMPAT_MODE

      SUBROUTINE XML_WRITE_GGA_COMPAT_MODE
        USE prec
        USE compat_gga
        IMPLICIT NONE
        INTEGER :: IDUM
        REAL(q) :: RDUM
        COMPLEX(q)  :: CDUM
        CHARACTER (LEN=1) :: CHARAC

        CALL XML_INCAR('GGA_COMPAT','L',IDUM,RDUM,CDUM,GGA_COMPAT,CHARAC,1)

      END SUBROUTINE XML_WRITE_GGA_COMPAT_MODE
