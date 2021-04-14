
#include "symbol.inc"
      MODULE us
      USE prec
!***********************************************************************
! RCS:  $Id: us.F,v 1.6 2003/06/27 13:22:23 kresse Exp kresse $
!
! a few comments concerning implementation:
! It would be simplest to have assumed shape arrays for CDIJ...
! But I have not done this up do now.
! One problem (which I really hate) is that VASP uses complex
! to real FFT. This requires that the same array is in some
! places addressed as a real, and in other places as a complex
! array. Because no cast exists in f90, I have to use the
! strange INTERFACE construct below to get at least some error
! checking.
! This fact also keeps me from using assumed shape arrays (gK)
!***********************************************************************
!***********************************************************************
!
! specify the interface for SETDIJ and DEPLE
! the dummy argument CVTOT differs here from the actual definition
! in the subroutine used below
!
!***********************************************************************
  INTERFACE
    SUBROUTINE SETDIJ(WDES, GRIDC_,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO, LOVERL, &
          LMDIM,CDIJ,CQIJ, CVTOT_, IRDMAA,IRDMAX)
      USE prec
      USE pseudo
      USE poscar
      USE mpimy
      USE mgrid
      USE lattice
      USE wave
      USE asa
      USE paw
      USE constant
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC_,GRIDUS
      TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES

      INTEGER IRDMAX      ! allocation required for augmentation
      INTEGER IRDMAA      ! actual maximum augmentation index
      OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
              CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q) CVTOT_(GRIDC_%MPLWV,WDES%NCDIJ)
      LOGICAL LOVERL
    END SUBROUTINE SETDIJ
  END INTERFACE
  INTERFACE
    SUBROUTINE SETDIJ_(WDES, GRIDC_,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO, LOVERL, &
        LMDIM,CDIJ,CQIJ, CVTOT_, LDIAGONAL_TERMS, IRDMAA,IRDMAX, DISPL)
      USE prec
      USE pseudo
      USE poscar
      USE mpimy
      USE mgrid
      USE lattice
      USE wave
      USE asa
      USE paw
      USE constant
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d),TARGET  :: GRIDC_,GRIDUS
      TYPE (grid_3d),POINTER :: GRIDC
      TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES

      INTEGER  IRDMAX      ! allocation required for augmentation
      INTEGER  IRDMAA      ! actual maximum augmentation index
      OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
               CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q),TARGET :: CVTOT_(GRIDC_%MPLWV,WDES%NCDIJ)
      LOGICAL  LOVERL
      REAL(q)  DISPL(3,T_INFO%NIONS)
      LOGICAL  LDIAGONAL_TERMS
    END SUBROUTINE SETDIJ_
  END INTERFACE
  INTERFACE
    SUBROUTINE DEPLE(WDES, GRID_SOFT,GRIDC_,GRIDUS,C_TO_US, &
          LATT_CUR,P,T_INFO,SYMM, LOVERL, SOFT_TO_C,&
          LMDIM,CRHODE, CHTOT_,CHDEN, IRDMAX )
      USE prec
      USE base
      USE charge
      USE pseudo
      USE poscar
      USE mpimy
      USE mgrid
      USE lattice
      USE wave
      USE asa
      USE paw

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRID_SOFT,GRIDC_,GRIDUS
      TYPE (transit)     C_TO_US
      TYPE (transit)     SOFT_TO_C
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES
      TYPE (symmetry)    SYMM

      INTEGER IRDMAX,ISP      ! allocation required for augmentation
      OVERLAP CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q) CHTOT_(GRIDC_%MPLWV,WDES%NCDIJ)
      COMPLEX(q) CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
      LOGICAL LOVERL
    END SUBROUTINE DEPLE
  END INTERFACE
  INTERFACE
    SUBROUTINE SETDIJ_AVEC_(WDES, GRIDC_,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO, LOVERL, &
          LMDIM,CDIJ,AVTOT_, NONLR_S, NONL_S, IRDMAX, DISPL )
      USE prec
      USE pseudo
      USE poscar
      USE mgrid
      USE lattice
      USE wave
      USE paw
      USE nonl_high
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)!

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d),TARGET  :: GRIDC_,GRIDUS
      TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES
      INTEGER  IRDMAX      ! allocation required for augmentation
      OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q) :: AVTOT_(GRIDC_%MPLWV,3)
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct)  NONL_S
      LOGICAL  LOVERL
      REAL(q)  DISPL(3,T_INFO%NIONS)
    END SUBROUTINE SETDIJ_AVEC_
  END INTERFACE
  INTERFACE
    SUBROUTINE SETDIJ_AVEC_ONE_CENTER(WDES, GRIDC_,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO, LOVERL, &
        LMDIM,CDIJ,AVTOT_, NONLR_S, NONL_S, IRDMAX, DISPL,  A_ONE_CENTER, LMMAX)
      USE prec
      USE pseudo
      USE poscar
      USE mgrid
      USE lattice
      USE wave
      USE paw
      USE nonl_high
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)!

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d),TARGET  :: GRIDC_,GRIDUS
      TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES
      INTEGER  IRDMAX      ! allocation required for augmentation
      OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q) :: AVTOT_(GRIDC_%MPLWV,3)
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct)  NONL_S
      LOGICAL  LOVERL
      INTEGER  LMMAX
      REAL(q)  DISPL(3,T_INFO%NIONS)
      REAL(q)  A_ONE_CENTER(3,LMMAX,T_INFO%NIONS)
    END SUBROUTINE SETDIJ_AVEC_ONE_CENTER
  END INTERFACE
  INTERFACE
    SUBROUTINE AUGMENTATION_CHARGE( &
           WDES, GRIDC_, GRIDUS, C_TO_US, &
           LATT_CUR, P, T_INFO, SYMM, LOVERL, &
           LMDIM, CRHODE, CHTOT_, IRDMAX, DISPL)
      USE prec
      USE base
      USE charge
      USE pseudo
      USE poscar
      USE mpimy
      USE mgrid
      USE lattice
      USE wave
      USE asa
      USE paw
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d),TARGET ::  GRIDC_,GRIDUS
      TYPE (grid_3d),POINTER :: GRIDC
      TYPE (transit)     C_TO_US
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES
      TYPE (symmetry)    SYMM

      INTEGER   IRDMAX         ! allocation required for augmentation
      INTEGER   LMDIM
      OVERLAP   CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      RGRID,TARGET  :: CHTOT_(DIMREAL(GRIDC_%MPLWV),WDES%NCDIJ)
      LOGICAL   LOVERL,LADDITIONAL
      REAL(q)   DISPL(3,T_INFO%NIONS)
    END SUBROUTINE AUGMENTATION_CHARGE
  END INTERFACE

   INTERFACE
      SUBROUTINE CURRENT_AUGMENTATION( &
           WDES, GRIDC_, GRIDUS, C_TO_US, &
           LATT_CUR, P, T_INFO, &
           LMDIM, CRHODE, JTOT_, IRDMAX, DISPL, LMONOPOLE )
      USE prec
      USE base
      USE pseudo
      USE poscar
      USE lattice
      USE mgrid
      USE wave

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d),TARGET ::  GRIDC_,GRIDUS
      TYPE (grid_3d),POINTER :: GRIDC
      TYPE (transit)     C_TO_US
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES

      INTEGER   IRDMAX         ! allocation required for augmentation
      INTEGER   LMDIM
      OVERLAP   CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q)   JTOT_(GRIDC_%MPLWV,3)
      LOGICAL   LADDITIONAL
      REAL(q)   DISPL(3,T_INFO%NIONS)
      LOGICAL :: LMONOPOLE
      END SUBROUTINE CURRENT_AUGMENTATION
   END INTERFACE

   INTERFACE
      SUBROUTINE CURRENT_AUGMENTATION_CMPLX( &
           WDES, GRIDC_, GRIDUS, C_TO_US, &
           LATT_CUR, P, T_INFO, &
           LMDIM, CRHODE, JTOT_, IRDMAX, DISPL, LMONOPOLE )
      USE prec
      USE base
      USE pseudo
      USE poscar
      USE lattice
      USE mgrid
      USE wave

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d),TARGET ::  GRIDC_,GRIDUS
      TYPE (grid_3d),POINTER :: GRIDC
      TYPE (transit)     C_TO_US
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES

      INTEGER   IRDMAX         ! allocation required for augmentation
      INTEGER   LMDIM
      COMPLEX(q)   CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q)   JTOT_(GRIDC_%MPLWV,3)
      LOGICAL   LADDITIONAL
      REAL(q)   DISPL(3,T_INFO%NIONS)
      LOGICAL :: LMONOPOLE
      END SUBROUTINE CURRENT_AUGMENTATION_CMPLX
   END INTERFACE

      CONTAINS

      !DMFTwDFT
      SUBROUTINE DEPSUM_DMFT_IBZ(W,WDES, LMDIM, CRHODE, LOVERL)
      USE prec
      USE wave
      USE mlwf
      USE wave_high

      IMPLICIT NONE

      TYPE (wavespin)    W
      TYPE (wavedes)     WDES
      LOGICAL LOVERL
      INTEGER LMDIM
      OVERLAP CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)

      INTEGER ISP, NT, NK, N, ISPINOR, ISPINOR_, LMBASE, LMBASE_, NIS, &
           LMMAXC, NI, L, LP
      REAL(q) WEIGHT, norm, tot_n
      TYPE (wavedes1)     WDES1
      TYPE (wavefun1)     W1
      INTEGER :: ierr,N2,NBG
      COMPLEX(q),ALLOCATABLE  :: CPROJ_TOT(:,:)!,CPROJ_TOT2(:,:)

      PROFILING_START('depsum')

      IF (.NOT.LOVERL) THEN
         PROFILING_STOP('depsum')
         RETURN
      ENDIF

      !!!!!!ADDED by Hyowon Park!!!!!!!!!!
      CALL SETWDES(WDES, WDES1, 0)
      CALL NEWWAV(W1 , WDES1, .FALSE.)

!=======================================================================
! initialise to zero
!=======================================================================
      CRHODE=0

!=======================================================================
! loop over all bands and k-points
!=======================================================================
      spin:   DO ISP=1,WDES%ISPIN
      kpoint: DO NK=1,WDES%NKPTS
#ifdef MPI
      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
      CALL SETWDES(W%WDES,WDES1,NK)

      num_band_max=band_win_(2,NK)-band_win_(1,NK)+1
      ALLOCATE(CPROJ_TOT(WDES%NPROD,num_band_max))
      !ALLOCATE(CPROJ_TOT2(WDES%NPROD,num_band_max))
      CPROJ_TOT=0._q!;CPROJ_TOT2=0._q;
      DO N=1,WDES%NBANDS
         CALL W1_COPY(ELEMENT(W, WDES1, N, ISP), W1)
         NBG=WDES%NB_LOW+(N-1)*WDES%NB_PAR
         if (NBG .GE. band_win_(1,NK) .AND. NBG .LE. band_win_(2,NK)) then
           CPROJ_TOT(:,NBG-band_win_(1,NK)+1)=W1%CPROJ(:)
         endif
      ENDDO
      !CALL
      !MPI_ALLREDUCE(CPROJ_TOT,CPROJ_TOT2,WDES%NPROD*num_band_max,MPI_DOUBLE_COMPLEX,MPI_SUM,WDES%COMM%MPI_COMM,ierr)

      band:   DO N=1,WDES%NBANDS
      CALL W1_COPY(ELEMENT(W, WDES1, N, ISP), W1)
      NBG=WDES%NB_LOW+(N-1)*WDES%NB_PAR
      !if (NBG .GE. band_win_(1,NK) .AND. NBG .LE. band_win_(2,NK)) then
      if (NBG .GE. band_win_(1,NK) .AND. NBG .LE. band_win_(1,NK)+num_wann_-1) then
         !print *, NBG
         norm=0._q
         W1%CPROJ=0._q
         !DO N2=1,num_band_max
         DO N2=1,num_wann_
         !   norm=norm+real(DMFT_evec_(N2,NBG-band_win_(1,NK)+1,NK)*conjg(DMFT_evec_(N2,NBG-band_win_(1,NK)+1,NK)))
            W1%CPROJ(:)=W1%CPROJ(:)+CPROJ_TOT(:,N2)*DMFT_evec_(N2,NBG-band_win_(1,NK)+1,NK)
         ENDDO
         WEIGHT=WDES%RSPIN*WDES%WTKPT(NK)*DMFT_eval_(NBG-band_win_(1,NK)+1,NK)
         !if (NBG .GE. band_win_(1,NK) .AND. NBG .LE. band_win_(1,NK)+31) then
         !else
         !   WEIGHT=0
         !endif
         !WEIGHT=WDES%RSPIN*WDES%WTKPT(NK)
         !print *, NK,NBG,DMFT_eval_(NBG-band_win_(1,NK)+1,NK)
         !if (abs(norm-1.0).GE.1e-6) then 
         !   print *, NK, NBG, norm !DMFT_eval_(NBG-band_win_(1,NK)+1,NK),W%FERWE(N,NK,ISP)
         !endif
      else if (NBG .LT. band_win_(1,NK)) then
         !WEIGHT=0._q
         WEIGHT=WDES%RSPIN*WDES%WTKPT(NK)*W%FERWE(N,NK,ISP)
      else
         WEIGHT=0._q
      endif
      tot_n=tot_n+WEIGHT
      !WEIGHT=WDES%RSPIN*W%FERWE(N,NK,ISP)*WDES%WTKPT(NK)

      spinor: DO ISPINOR =0,WDES%NRSPINORS-1
      DO ISPINOR_=0,WDES%NRSPINORS-1

      LMBASE =ISPINOR *(WDES%NPRO/2)
      LMBASE_=ISPINOR_*(WDES%NPRO/2)
      NIS   =1
      typ:  DO NT=1,WDES%NTYP
      LMMAXC=WDES%LMMAX(NT)
      IF (LMMAXC==0) GOTO 210

      ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1
      
!DIR$ IVDEP
!OCL NOVREC
        DO L=1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
        DO LP=1,LMMAXC
           CRHODE(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)=CRHODE(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)+&
              !WEIGHT*W%CPROJ(L+LMBASE,N,NK,ISP)*GCONJG(W%CPROJ(LP+LMBASE_,N,NK,ISP))
              WEIGHT*W1%CPROJ(L+LMBASE)*GCONJG(W1%CPROJ(LP+LMBASE_))
        ENDDO
        ENDDO

      LMBASE = LMMAXC+LMBASE
      LMBASE_= LMMAXC+LMBASE_

      ENDDO ion

  210 NIS = NIS+WDES%NITYP(NT)
      ENDDO typ
      ENDDO
      ENDDO spinor

      ENDDO band
      DEALLOCATE(CPROJ_TOT)
      ENDDO kpoint
      ENDDO spin
      !STOP
      !print *, tot_n
      ! sum over all bands
#ifdef realmode
      CALLMPI(M_sum_d(WDES%COMM_INTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ))
      CALLMPI(M_sum_d(WDES%COMM_KINTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ))
#else
      CALLMPI(M_sum_z(WDES%COMM_INTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ))
      CALLMPI(M_sum_z(WDES%COMM_KINTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ))
#endif

      PROFILING_STOP('depsum')

      RETURN
      END SUBROUTINE


!************************ SUBROUTINE DEPSUM ****************************
!
! this subroutine calculates  the total one center (on site) "occupancy" 
! matrix of each ll'mm'augmentation channel from the FERMI weights, 
! the weight of each k-point WTKPT and  the wavefunction character 
! of all bands
! result is stored in CRHODE (thesis gK 10.32)
!
!
!***********************************************************************


      SUBROUTINE DEPSUM(W,WDES, LMDIM, CRHODE, LOVERL)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavespin)    W
      TYPE (wavedes)     WDES
      LOGICAL LOVERL
      INTEGER LMDIM
      OVERLAP CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)

      INTEGER ISP, NT, NK, N, ISPINOR, ISPINOR_, LMBASE, LMBASE_, NIS, &
           LMMAXC, NI, L, LP
      REAL(q) WEIGHT

      PROFILING_START('depsum')

      IF (.NOT.LOVERL) THEN
         PROFILING_STOP('depsum')
         RETURN
      ENDIF

!=======================================================================
! initialise to zero
!=======================================================================
      CRHODE=0
!=======================================================================
! loop over all bands and k-points
!=======================================================================
      spin:   DO ISP=1,WDES%ISPIN
      kpoint: DO NK=1,WDES%NKPTS
#ifdef MPI
      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
      band:   DO N=1,WDES%NBANDS

      WEIGHT=WDES%RSPIN*W%FERWE(N,NK,ISP)*WDES%WTKPT(NK)

      spinor: DO ISPINOR =0,WDES%NRSPINORS-1
      DO ISPINOR_=0,WDES%NRSPINORS-1

      LMBASE =ISPINOR *(WDES%NPRO/2)
      LMBASE_=ISPINOR_*(WDES%NPRO/2)
      NIS   =1
      typ:  DO NT=1,WDES%NTYP
      LMMAXC=WDES%LMMAX(NT)
      IF (LMMAXC==0) GOTO 210

      ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1

!DIR$ IVDEP
!OCL NOVREC
        DO L=1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
        DO LP=1,LMMAXC
           CRHODE(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)=CRHODE(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)+ &
              WEIGHT*W%CPROJ(L+LMBASE,N,NK,ISP)*GCONJG(W%CPROJ(LP+LMBASE_,N,NK,ISP))
        ENDDO
        ENDDO
      
      LMBASE = LMMAXC+LMBASE
      LMBASE_= LMMAXC+LMBASE_

      ENDDO ion

  210 NIS = NIS+WDES%NITYP(NT)
      ENDDO typ
      ENDDO
      ENDDO spinor

      ENDDO band
      ENDDO kpoint
      ENDDO spin
      ! sum over all bands
#ifdef realmode
      CALLMPI( M_sum_d(WDES%COMM_INTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ))
      CALLMPI( M_sum_d(WDES%COMM_KINTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ))
#else
      CALLMPI( M_sum_z(WDES%COMM_INTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ))
      CALLMPI( M_sum_z(WDES%COMM_KINTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ))
#endif

      PROFILING_STOP('depsum')

      RETURN
      END SUBROUTINE


!************************ SUBROUTINE DEPSUM_SYM ************************
!
! this subroutine calculates  the total one center (on site) "occupancy" 
! applying 
! symmetry to determine the wavefunctions at each k-point
! it should be fully compatible with the Fock like routines but
! of course is slightly slower than DEPSUM
! if this routine is used to construct CRHODE no
! a posteriori symmetrization of the charge density should be required
!
!***********************************************************************

    SUBROUTINE DEPSUM_SYM(W,WDES, LMDIM, CRHODE, LOVERL, P, LATT_CUR)
      USE prec
      USE sym_prec
      USE wave_high
      USE full_kpoints
      USE pseudo
      USE lattice
      IMPLICIT NONE

      TYPE (wavespin)    W
      TYPE (wavedes)     WDES
      LOGICAL LOVERL
      INTEGER LMDIM
      OVERLAP CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      TYPE (potcar)       P(:)
      TYPE (latt) LATT_CUR
      
   ! local
      INTEGER ISP, ISP_IRZ, NT, NK, N, ISPINOR, ISPINOR_, LMBASE, LMBASE_, NIS, &
           LMMAXC, NI, L, LP
      LOGICAL :: LSHIFT
      TYPE( rotation_handle), POINTER :: ROT_HANDLE
      TYPE (wavedes1)     WDES1, WDES1_IRZ
      TYPE (wavefun1)     W1
      REAL(q) WEIGHT

      PROFILING_START('depsum_sym')

      IF (.NOT.LOVERL) THEN
         PROFILING_STOP('depsum_sym')
         RETURN
      ENDIF

      CALL CHECK_FULL_KPOINTS
      NULLIFY(ROT_HANDLE)

      CALL SETWDES(WDES, WDES1, 0)
      CALL NEWWAV(W1 , WDES1, .FALSE.)
!=======================================================================
! initialise to zero
!=======================================================================
      CRHODE=0
!=======================================================================
! loop over all bands and k-points
!=======================================================================
      spin:   DO ISP=1,WDES%ISPIN
      kpoint: DO NK=1,KPOINTS_FULL%NKPTS
#ifdef MPI
         IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
         CALL SETWDES(WDES, WDES1, NK)
         CALL SETWDES(WDES, WDES1_IRZ, KPOINTS_FULL%NEQUIV(NK))
         
         ISP_IRZ=ISP
         IF (KPOINTS_FULL%SPINFLIP(NK)==1) THEN
            ISP_IRZ=3-ISP
         ENDIF
         
         band:   DO N=1,WDES%NBANDS

            WEIGHT=WDES%RSPIN*KPOINTS_FULL%WTKPT(NK)*W%FERWE(N,KPOINTS_FULL%NEQUIV(NK),ISP)

            IF (NK <= WDES%NKPTS) THEN
               CALL W1_COPY(ELEMENT(W, WDES1, N, ISP), W1)
            ELSE
#ifndef gammareal
           ! 
           ! use symmetry to construct wave function charakter at k
               CALL ROTATE_WAVE_CHARACTER(ROT_HANDLE, P, LATT_CUR, WDES1, ELEMENT(W, WDES1_IRZ, N, ISP_IRZ), W1)
#endif
            ENDIF

            spinor: DO ISPINOR =0,WDES%NRSPINORS-1
               DO ISPINOR_=0,WDES%NRSPINORS-1

                  LMBASE =ISPINOR *(WDES%NPRO/2)
                  LMBASE_=ISPINOR_*(WDES%NPRO/2)
                  
                  NIS   =1
                  typ:  DO NT=1,WDES%NTYP
                     LMMAXC=WDES%LMMAX(NT)
                     IF (LMMAXC==0) GOTO 210

                     ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1
!DIR$ IVDEP
!OCL NOVREC
                        DO L=1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
                           DO LP=1,LMMAXC
                              CRHODE(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)=CRHODE(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)+ &
                                   WEIGHT*W1%CPROJ(L+LMBASE)*GCONJG(W1%CPROJ(LP+LMBASE_))
                           ENDDO
                        ENDDO
                        
                        LMBASE = LMMAXC+LMBASE
                        LMBASE_= LMMAXC+LMBASE_

                     ENDDO ion

210                  NIS = NIS+WDES%NITYP(NT)
                  ENDDO typ
               ENDDO
            ENDDO spinor

         ENDDO band
      ENDDO kpoint
      ENDDO spin
      ! sum over all bands
#ifdef realmode
      CALLMPI( M_sum_d(WDES%COMM_INTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ))
      CALLMPI( M_sum_d(WDES%COMM_KINTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ))
#else
      CALLMPI( M_sum_z(WDES%COMM_INTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ))
      CALLMPI( M_sum_z(WDES%COMM_KINTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ))
#endif

      CALL DEALLOCATE_ROT_HANDLE(ROT_HANDLE)
      CALL DELWAV(W1, .FALSE.)

      PROFILING_STOP('depsum_sym')

      RETURN
    END SUBROUTINE DEPSUM_SYM



!************************ SUBROUTINE DEPSUM1 ****************************
!
! this subroutine calculates  the first order change of the 
! one center (on site) "occupancy" 
! matrix of each ll'mm'augmentation channel from the FERMI weights, 
! the weight of each k-point WTKPT and  the wavefunction character 
! of all bands
! result is stored in CRHODE1
!
!***********************************************************************

      SUBROUTINE DEPSUM1(W0, W1, WDES, LMDIM, CRHODE, LOVERL)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavespin)    W0
      TYPE (wavespin)    W1
      TYPE (wavedes)     WDES
      LOGICAL LOVERL
      INTEGER LMDIM
      OVERLAP CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)

      INTEGER ISP, NT, NK, N, ISPINOR, ISPINOR_, LMBASE, LMBASE_, NIS, &
           LMMAXC, NI, L, LP
      REAL(q) WEIGHT0, WEIGHT1


      IF (.NOT.LOVERL) RETURN
!=======================================================================
! initialise to zero
!=======================================================================
      CRHODE=0
!=======================================================================
! loop over all bands and k-points
!=======================================================================
      spin:   DO ISP=1,WDES%ISPIN
      kpoint: DO NK=1,WDES%NKPTS
#ifdef MPI
      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
      band:   DO N=1,WDES%NBANDS

      WEIGHT0=WDES%RSPIN*W0%FERWE(N,NK,ISP)*WDES%WTKPT(NK)
      WEIGHT1=WDES%RSPIN*W1%FERWE(N,NK,ISP)*WDES%WTKPT(NK)

      spinor: DO ISPINOR =0,WDES%NRSPINORS-1
      DO ISPINOR_=0,WDES%NRSPINORS-1

      LMBASE =ISPINOR *(WDES%NPRO/2)
      LMBASE_=ISPINOR_*(WDES%NPRO/2)

      NIS   =1
      typ:  DO NT=1,WDES%NTYP
      LMMAXC=WDES%LMMAX(NT)
      IF (LMMAXC==0) GOTO 210

      ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1

!DIR$ IVDEP
!OCL NOVREC
        DO L=1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
        DO LP=1,LMMAXC
           CRHODE(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)=CRHODE(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)+ &
                WEIGHT1*W0%CPROJ(L+LMBASE,N,NK,ISP)*GCONJG(W0%CPROJ(LP+LMBASE_,N,NK,ISP))+ &
                WEIGHT0*W1%CPROJ(L+LMBASE,N,NK,ISP)*GCONJG(W0%CPROJ(LP+LMBASE_,N,NK,ISP))+ &
                WEIGHT0*W0%CPROJ(L+LMBASE,N,NK,ISP)*GCONJG(W1%CPROJ(LP+LMBASE_,N,NK,ISP))
        ENDDO
        ENDDO
      
      LMBASE = LMMAXC+LMBASE
      LMBASE_= LMMAXC+LMBASE_

      ENDDO ion

  210 NIS = NIS+WDES%NITYP(NT)
      ENDDO typ
      ENDDO
      ENDDO spinor

      ENDDO band
      ENDDO kpoint
      ENDDO spin
      ! sum over all bands
#ifdef realmode
      CALLMPI( M_sum_d(WDES%COMM_INTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ))
      CALLMPI( M_sum_d(WDES%COMM_KINTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ))
#else
      CALLMPI( M_sum_d(WDES%COMM_INTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ*2))
      CALLMPI( M_sum_d(WDES%COMM_KINTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ*2))
#endif


      RETURN
      END SUBROUTINE



!************************ SUBROUTINE DEPATO ****************************
!
! set the "occupancy" matrix CRHODE to the values during
! the pseudopotential generation (QATO)
!
! in the spinpolarized collinear case the CRHODE(:,2) is set
! to CRHODE(:,1) / ZVAL * ATOMOM
! in the non collinear case the magnetization in each direction 
! x,y and z is set such that the moment is paralle to ATOMOM
!
!***********************************************************************

      SUBROUTINE DEPATO(WDES, LMDIM, CRHODE, LOVERL, P, T_INFO)
      USE prec
      USE pseudo
      USE wave
      USE poscar
      IMPLICIT NONE

      TYPE (wavedes)     WDES
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      LOGICAL LOVERL
      INTEGER LMDIM
      OVERLAP CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
! local variables
      INTEGER NI,NIP,LOW,LM,NT,LL,LHI,MMAX,L,LP,M,ISP
      REAL(q) :: VALUE(WDES%NCDIJ)

      IF (.NOT.LOVERL) RETURN
      CRHODE=0

  ion: DO NI=1,T_INFO%NIONS
      NIP=NI_LOCAL(NI, WDES%COMM_INB)
      IF (NIP==0) CYCLE ion ! not on local node

      LOW=1
      LM =1
      NT=T_INFO%ITYP(NI)

      IF (.NOT.ASSOCIATED(P(NT)%QATO)) CYCLE ion
      
      !
      ! this little trick does pushes the PAW occupancies in the right
      ! direction for core level shift calculations
      !
      VALUE(1) = P(NT)%ZVALF/P(NT)%ZVALF_ORIG
      IF (WDES%NCDIJ==2) THEN
         VALUE(2)= T_INFO%ATOMOM(NI)/P(NT)%ZVALF
      ELSE IF (WDES%NCDIJ == 4) THEN
         VALUE(2)= T_INFO%ATOMOM(1+(NI-1)*3)/P(NT)%ZVALF
         VALUE(3)= T_INFO%ATOMOM(2+(NI-1)*3)/P(NT)%ZVALF
         VALUE(4)= T_INFO%ATOMOM(3+(NI-1)*3)/P(NT)%ZVALF
      ENDIF

      block: DO
         LL=P(NT)%LPS(LOW)
         ! search block with same L
         DO LHI=LOW,P(NT)%LMAX
            IF (LL/=P(NT)%LPS(LHI)) EXIT
         ENDDO
         LHI=LHI-1
         MMAX=2*LL+1

         DO ISP=1,WDES%NCDIJ

         DO L =LOW,LHI
         DO LP=LOW,LHI
         DO M =0,MMAX-1
            CRHODE(LM+(L-LOW)*MMAX+M,LM+(LP-LOW)*MMAX+M,NIP,ISP)= &
                    P(NT)%QATO(L,LP)*VALUE(ISP)
         ENDDO
         ENDDO
         ENDDO

         ENDDO

         ! set new LOW value and LM value and go on
         LM=LM+(LHI-LOW+1)*MMAX
         LOW=LHI+1
         IF (LOW > P(NT)%LMAX) EXIT block
      ENDDO block

      ENDDO ion
      RETURN
      END SUBROUTINE


!************************* SUBROUTINE FORDEP ***************************
!
! this subroutine calculates the forces, which are related to the
! the change in the position of the augmentation charges in a fixed
! local potential
! finite central differences are usually used
! the accuracy of the routine is at least 7 significant digits
! Equ. (54) of PRB 59, 1758
!
! NOTE: the routine resets the CDIJ array 
! CDIJ must be restored after calling the routine
!
!
!***********************************************************************

 
      SUBROUTINE FORDEP(WDES, GRIDC,GRIDUS,C_TO_US, &
        LATT_CUR,P,T_INFO, LOVERL, &
        LMDIM, CDIJ, CQIJ, CRHODE, CVTOT, IRDMAX, DISPL0, FORNL)
      USE prec

      USE pseudo
      USE poscar
      USE mpimy
      USE mgrid
      USE lattice
      USE wave
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavedes)     WDES
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC,GRIDUS
      TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
      TYPE (latt)        LATT_CUR

      INTEGER  IRDMAX      ! allocation required for augmentation
      INTEGER  IRDMAA      ! actual maximum augmentation index
      OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
               CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
               CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q) CVTOT(GRIDC%MPLWV,WDES%NCDIJ)
      REAL(q)   DISPL0(3,T_INFO%NIONS)
      REAL(q)   FORNL(3,T_INFO%NIONS)
      LOGICAL   LOVERL
! work array
      REAL(q)   DISPL(3,T_INFO%NIONS)
      REAL(q)   EDEP(T_INFO%NIONS)
      REAL(q)   FORAUG(3,T_INFO%NIONS)
      REAL(q):: DIS=1E-4_q

      PROFILING_START('fordep')

      IF (.NOT. LOVERL) THEN
         PROFILING_STOP('fordep')
         RETURN
      ENDIF

      ! we need the spinor representation of CRHODE at this point
      CALL US_FLIP(WDES, LMDIM, CRHODE, LOVERL, .TRUE.)
!TEST
!     WRITE(*,'(4E14.7)') DIS
!     DIS=1E-2
!1000 DIS=DIS/2
!TEST
      FORAUG=0
!=======================================================================
! calculate the contribution to the energy from the augmentation charges
! for displacement X
! (and displacement -X if central differences are used)
!=======================================================================

      dir: DO IDIR=1,3
      DISPL=DISPL0
      DISPL(IDIR,:)=DISPL0(IDIR,:)-DIS/2
      CALL SETDIJ_(WDES, GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO, LOVERL, &
                  LMDIM,CDIJ,CQIJ,CVTOT,.FALSE.,IRDMAA,IRDMAX,DISPL)
      NIS=1
      DO NT=1,WDES%NTYP
        LMMAXC=WDES%LMMAX(NT)
        DO NI=NIS,WDES%NITYP(NT)+NIS-1
          CALL FORDE1(LMDIM,WDES%NIONS,WDES%NCDIJ,LMMAXC,NI,CDIJ,CRHODE,ADD)
          EDEP(NI)=ADD
        ENDDO
      NIS = NIS+WDES%NITYP(NT)
      ENDDO
      DISPL=DISPL0
      DISPL(IDIR,:)=DISPL0(IDIR,:)+DIS/2
      CALL SETDIJ_(WDES, GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,LOVERL, &
                  LMDIM,CDIJ,CQIJ,CVTOT,.FALSE.,IRDMAA,IRDMAX,DISPL)
      NIS=1
      DO NT=1,WDES%NTYP
        LMMAXC=WDES%LMMAX(NT)
        DO NI=NIS,WDES%NITYP(NT)+NIS-1
          CALL FORDE1(LMDIM,WDES%NIONS,WDES%NCDIJ,LMMAXC,NI,CDIJ,CRHODE,ADD)
          NIP=NI_GLOBAL(NI, WDES%COMM_INB)     !  local storage index
          FORAUG(IDIR,NIP)=FORAUG(IDIR,NIP)-(ADD-EDEP(NI))/DIS
        ENDDO
        NIS = NIS+WDES%NITYP(NT)
      ENDDO
      ENDDO dir

      CALLMPI( M_sum_d(WDES%COMM_INB, FORAUG(1,1),T_INFO%NIONS*3))

!TEST
!      WRITE(*,'(4E20.12)') DIS
!      WRITE(*,'("a",3F15.9 )') FORAUG

!      GOTO 1000
!TEST
      FORNL=FORNL+FORAUG

      ! back to original representation
      CALL US_FLIP(WDES, LMDIM, CRHODE, LOVERL, .FALSE.)

      PROFILING_STOP('fordep')

      RETURN
      END SUBROUTINE

!
! small helper routine
!
      SUBROUTINE FORDE1(LMDIM,NIONS,ISPIN,LMMAXC,NI,CDIJ,CRHODE,ADD)
      USE prec
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      OVERLAP CRHODE(LMDIM,LMDIM,NIONS,ISPIN)
      OVERLAP CDIJ(LMDIM,LMDIM,NIONS,ISPIN)

      ADD=0

      DO ISP=1,ISPIN
      DO L=1,LMMAXC
      DO LP=1,LMMAXC
#ifdef realmode
        ADD=ADD+ CDIJ(LP,L,NI,ISP)*CRHODE(LP,L,NI,ISP)
#else
        ADD=ADD+ CDIJ(LP,L,NI,ISP)*CONJG(CRHODE(LP,L,NI,ISP))
#endif
      ENDDO; ENDDO; ENDDO

      RETURN

      END SUBROUTINE


!************************* SUBROUTINE STRDEP ***************************
!
!  subroutine for calculating the contributions to the stress,
!  related to the augmentation charges
!  uncomment CTEST-lines if you want to test finite differences
!
!***********************************************************************

      SUBROUTINE STRDEP(WDES, GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO, LOVERL, &
       LMDIM,CDIJ,CQIJ,CRHODE, CVTOT, IRDMAX, ISIF,AUGSIF)
      USE prec

      USE pseudo
      USE poscar
      USE mpimy
      USE mgrid
      USE lattice
      USE wave
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavedes)     WDES
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRIDC,GRIDUS
      TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
      TYPE (latt)        LATT_CUR,LATT_FIN

      INTEGER  IRDMAX      ! allocation required for augmentation
      INTEGER  IRDMAA      ! actual maximum augmentation index
      OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
               CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
               CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q)  CVTOT(GRIDC%MPLWV,WDES%NCDIJ)
      LOGICAL  LOVERL
      REAL(q)     AUGSIF(3,3)

      IF (.NOT. LOVERL) RETURN

      ! we need the spinor representation of CRHODE at this point
      CALL US_FLIP(WDES, LMDIM, CRHODE, LOVERL, .TRUE.)

      DIS=fd_displacement
!TEST which step should be used in the finite differences
!      DIS=1E-3
! 1000 DIS=DIS/2
!TEST
!=======================================================================
! calculate the contribution to the energy from the augmentation-hole
! for undistorted lattice
!=======================================================================
      AUGSIF=0
      CALL LATTIC(LATT_CUR)

      CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,LOVERL, &
                  LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

      EAUG=0
      NIS=1

      DO NT=1,WDES%NTYP
        LMMAXC=WDES%LMMAX(NT)
        DO NI=NIS,WDES%NITYP(NT)+NIS-1
        DO ISP=1,WDES%NCDIJ
        DO L=1,LMMAXC
        DO LP=1,LMMAXC
#ifdef realmode
          EAUG=EAUG+ REAL( CDIJ(LP,L,NI,ISP)*CRHODE(LP,L,NI,ISP),KIND=q)
#else
          EAUG=EAUG+ REAL( CDIJ(LP,L,NI,ISP)*CONJG(CRHODE(LP,L,NI,ISP)),KIND=q)
#endif
        ENDDO; ENDDO; ENDDO; ENDDO

        NIS = NIS+WDES%NITYP(NT)
      ENDDO
!=======================================================================
! calculate the contribution to the energy from the augmentation charges
! for distortion X
! 1. order finite differences are used
! (for stress high precision is not really required)
!=======================================================================
      DO IDIR=1,3
      DO JDIR=1,3

      LATT_FIN=LATT_CUR
      IF (ISIF==1) THEN
!  only isotrop pressure
        DO I=1,3
           DO J=1,3
              LATT_FIN%A(I,J)=LATT_CUR%A(I,J)*(1+DIS/3)
           ENDDO
        ENDDO
      ELSE
!  all directions
         DO I=1,3
            LATT_FIN%A(IDIR,I)=LATT_CUR%A(IDIR,I)+DIS*LATT_CUR%A(JDIR,I)
         ENDDO
      ENDIF
      CALL LATTIC(LATT_FIN)

      CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_FIN,P,T_INFO,LOVERL, &
           LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

      EAUGD=0
      NIS=1
      DO  NT=1,WDES%NTYP
        LMMAXC=WDES%LMMAX(NT)
        DO NI=NIS,WDES%NITYP(NT)+NIS-1
        DO ISP=1,WDES%NCDIJ
        DO L=1,LMMAXC
        DO LP=1,LMMAXC
#ifdef realmode
          EAUGD=EAUGD+ REAL( CDIJ(LP,L,NI,ISP)*CRHODE(LP,L,NI,ISP),KIND=q)
#else
          EAUGD=EAUGD+ REAL( CDIJ(LP,L,NI,ISP)*CONJG(CRHODE(LP,L,NI,ISP)),KIND=q)
#endif
        ENDDO; ENDDO; ENDDO; ENDDO

        NIS = NIS+WDES%NITYP(NT)
      ENDDO

      LATT_FIN=LATT_CUR
      IF (ISIF==1) THEN
!  only isotrop pressure
        DO I=1,3
           DO J=1,3
              LATT_FIN%A(I,J)=LATT_CUR%A(I,J)*(1-DIS/3)
           ENDDO
        ENDDO
      ELSE
!  all directions
         DO I=1,3
            LATT_FIN%A(IDIR,I)=LATT_CUR%A(IDIR,I)-DIS*LATT_CUR%A(JDIR,I)
         ENDDO
      ENDIF
      CALL LATTIC(LATT_FIN)

      CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_FIN,P,T_INFO,LOVERL, &
           LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

      EAUG=0
      NIS=1
      DO  NT=1,WDES%NTYP
        LMMAXC=WDES%LMMAX(NT)
        DO NI=NIS,WDES%NITYP(NT)+NIS-1
        DO ISP=1,WDES%NCDIJ
        DO L=1,LMMAXC
        DO LP=1,LMMAXC
#ifdef realmode
          EAUG=EAUG+ REAL( CDIJ(LP,L,NI,ISP)*CRHODE(LP,L,NI,ISP),KIND=q)
#else
          EAUG=EAUG+ REAL( CDIJ(LP,L,NI,ISP)*CONJG(CRHODE(LP,L,NI,ISP)),KIND=q)
#endif
        ENDDO; ENDDO; ENDDO; ENDDO

        NIS = NIS+WDES%NITYP(NT)
      ENDDO

      AUGSIF(IDIR,JDIR)=(EAUG-EAUGD)/2

!
!  only isotrop pressure terminate loop (IDIR=1 and JDIR=1)
!
      IF (ISIF==1) THEN
        AUGSIF(2,2)= AUGSIF(1,1)
        AUGSIF(3,3)= AUGSIF(1,1)
        GOTO 400 ! terminate (not very clean but who cares)
      ENDIF

      ENDDO
      ENDDO
!=======================================================================
! calculation finished  scale stress
!=======================================================================
  400 CONTINUE
      CALLMPI( M_sum_d(WDES%COMM_INB, AUGSIF, 9))

      AUGSIF=AUGSIF/DIS

#ifdef debug
!TEST
!      WRITE(*,'(E10.3,"  ",7E14.7)')DIS,AUGSIF
!      IF (DIS>1E-10) GOTO 1000
!      STOP
!TEST
#endif
      CALL US_FLIP(WDES, LMDIM, CRHODE, LOVERL, .FALSE.)

      RETURN
      END SUBROUTINE

!************************ SUBROUTINE CHARGE     ************************
!
! calculate the charge density on the soft and fine grid
! results are returned in the convention
! (total, magnetization (x,y,z))
! for all considered arrays
!
! CHDEN          charge density on soft grid
! CHTOT          total charge density on fine grid
! CRHODE         one center occupancies
!
!***********************************************************************

      SUBROUTINE SET_CHARGE(W, WDES, LOVERL, &
                  GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
                  LATT_CUR, P, SYMM, T_INFO, &
                  CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

      USE paw
      USE prec
      USE base
      USE charge
      USE pseudo
      USE poscar
      USE mpimy
      USE mgrid
      USE lattice
      USE wave
      USE mlwf

      IMPLICIT NONE

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRID,GRIDC,GRID_SOFT,GRIDUS
      TYPE (transit)     C_TO_US,SOFT_TO_C
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES
      TYPE (symmetry)    SYMM
      TYPE (wavespin)    W

      INTEGER LMDIM
      INTEGER IRDMAX
      INTEGER N_MIX_PAW
      REAL(q)  RHOLM(N_MIX_PAW,WDES%NCDIJ)
! on return the following arrays are set
      OVERLAP     CRHODE_cor(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ) ! occupancy of augmentation channels
      OVERLAP     CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ) ! occupancy of augmentation channels
      COMPLEX(q)  CHDEN_cor(GRID_SOFT%MPLWV,WDES%NCDIJ)         ! soft pseudo charge
      COMPLEX(q)  CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)         ! soft pseudo charge
      COMPLEX(q)  CHTOT(GRIDC%MPLWV,WDES%NCDIJ)             ! total charge
      LOGICAL LOVERL
      INTEGER ISP

      PROFILING_START('set_charge')

      IF (SYMM%ISYM==3) THEN
         CALL SOFT_CHARGE_SYM(GRID,GRID_SOFT,W, WDES, CHDEN)
      ELSE
        !DMFTwDFT
         if (LDMFT.eq..TRUE.) then
            !CALL SOFT_CHARGE_DMFT(GRID,GRID_SOFT,W, WDES, CHDEN, P, LATT_CUR)
            !CALL SOFT_CHARGE(GRID,GRID_SOFT,W, WDES, CHDEN)
            CALL SOFT_CHARGE_DMFT_IBZ(GRID,GRID_SOFT,W, WDES, CHDEN)
            !CHDEN(:,:)=CHDEN(:,:)+CHDEN_cor(:,:)
         else
            CALL SOFT_CHARGE(GRID,GRID_SOFT,W, WDES, CHDEN)
         endif
         !CALL SOFT_CHARGE(GRID,GRID_SOFT,W, WDES, CHDEN)
      ENDIF

      ! change storage convention to (total, magnetization)
      CALL RC_FLIP(CHDEN,GRID_SOFT,WDES%NCDIJ,.FALSE.)

      IF (SYMM%ISYM==2) THEN
         IF (WDES%LNONCOLLINEAR) THEN
            CALL RHOSYM(CHDEN(1,1),GRID_SOFT,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,1)
            IF (.NOT.WDES%LSPIRAL) &
           &   CALL SYMFIELD(CHDEN(1,2),GRID_SOFT,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,WDES%SAXIS,LATT_CUR)
         ELSE
            DO ISP=1,WDES%ISPIN
               CALL RHOSYM(CHDEN(1,ISP),GRID_SOFT,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,ISP)
            ENDDO
         ENDIF
      ENDIF

      IF (SYMM%ISYM==3) THEN
         CALL DEPSUM_SYM(W, WDES, LMDIM, CRHODE, LOVERL, P, LATT_CUR)
      ELSE
         !DMFTwDFT
         if (LDMFT.eq..TRUE.) then
            !CALL DEPSUM_DMFT(W, WDES, LMDIM, CRHODE, LOVERL,P,LATT_CUR)
            !CALL DEPSUM(W, WDES, LMDIM, CRHODE, LOVERL)
            CALL DEPSUM_DMFT_IBZ(W, WDES, LMDIM, CRHODE, LOVERL)
            !CRHODE(:,:,:,:)=CRHODE(:,:,:,:)+CRHODE_cor(:,:,:,:)
         else
            CALL DEPSUM(W, WDES, LMDIM, CRHODE, LOVERL)
         endif
         !CALL DEPSUM(W, WDES, LMDIM, CRHODE, LOVERL)
      ENDIF

      ! change storage convention to (total, magnetization)
      CALL US_FLIP(WDES, LMDIM, CRHODE, LOVERL, .FALSE.)

      CALL DEPLE(WDES,GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
                 LATT_CUR,P,T_INFO,SYMM, LOVERL, SOFT_TO_C,&
                 LMDIM, CRHODE, CHTOT,CHDEN, IRDMAX)

      CALL SET_RHO_PAW(WDES, P, T_INFO, LOVERL, WDES%NCDIJ, LMDIM, &
           CRHODE, RHOLM)

      IF (SYMM%ISYM==1) THEN
         IF (WDES%LNONCOLLINEAR) THEN
            CALL RHOSYM(CHTOT(1,1),GRIDC,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,1)
            IF (.NOT.WDES%LSPIRAL) &
           &   CALL SYMFIELD(CHTOT(1,2),GRIDC,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,WDES%SAXIS,LATT_CUR)
         ELSE
            DO ISP=1,WDES%ISPIN
               CALL RHOSYM(CHTOT(1,ISP),GRIDC,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,ISP)
            ENDDO
         ENDIF
      ENDIF

      PROFILING_STOP('set_charge')

    END SUBROUTINE SET_CHARGE

END MODULE

!************************ SUBROUTINE SETYLM_AUG ************************
!
! this subroutine performes the following tasks
! ) finds the points, which are within a certain cutoff around one ion
! ) calculates the distance of each of this points from the ion
! ) calculates the spherical harmonics Y_lm(Omega(r-R(ion))
! DISX,Y,Z are additional displacements of the ions
!
! mind that the cutoff-sphere extends up to PSDMAX*(NPSRNL-1)/NPSRNL
!
!***********************************************************************

    SUBROUTINE SETYLM_AUG(GRID,LATT_CUR,POSION,PSDMAX,NPSRNL, &
             LMYDIM,LYDIM,YLM,IRMAX,INDMAX,DISX,DISY,DISZ,DIST,NLI, &
             XS,YS,ZS)
      USE prec

      USE mpimy
      USE mgrid
      USE lattice
      USE asa
      USE constant
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)  GRID
      TYPE (latt)     LATT_CUR

      DIMENSION POSION(3)

      DIMENSION DIST(IRMAX)
      DIMENSION YLM(IRMAX,LMYDIM)
! work-arrays
      DIMENSION NLI(IRMAX)
      REAL(q) XS(IRMAX),YS(IRMAX),ZS(IRMAX)
      XS=0;YS=0;ZS=0
      IF ((LYDIM+1)**2 > LMYDIM) THEN
         WRITE(0,*)'internal error: LMYDIM is too small',LYDIM,LMYDIM
         STOP
      ENDIF
!=======================================================================
! find lattice points contained within the cutoff-sphere
! mind that the cutoff-sphere extends up to PSDMAX*(NPSRNL-1)/NPSRNL
! which is a somewhat strange convention
!=======================================================================
      F1=1._q/GRID%NGX
      F2=1._q/GRID%NGY
      F3=1._q/GRID%NGZ

      ARGSC=NPSRNL/PSDMAX
!-----------------------------------------------------------------------
! restrict loop to points contained within a cubus around the ion
!-----------------------------------------------------------------------
      D1= PSDMAX*LATT_CUR%BNORM(1)*GRID%NGX
      D2= PSDMAX*LATT_CUR%BNORM(2)*GRID%NGY
      D3= PSDMAX*LATT_CUR%BNORM(3)*GRID%NGZ

      N3LOW= INT(POSION(3)*GRID%NGZ-D3+10*GRID%NGZ+.99_q)-10*GRID%NGZ
      N2LOW= INT(POSION(2)*GRID%NGY-D2+10*GRID%NGY+.99_q)-10*GRID%NGY
      N1LOW= INT(POSION(1)*GRID%NGX-D1+10*GRID%NGX+.99_q)-10*GRID%NGX

      N3HI = INT(POSION(3)*GRID%NGZ+D3+10*GRID%NGZ)-10*GRID%NGZ
      N2HI = INT(POSION(2)*GRID%NGY+D2+10*GRID%NGY)-10*GRID%NGY
      N1HI = INT(POSION(1)*GRID%NGX+D1+10*GRID%NGX)-10*GRID%NGX

!-----------------------------------------------------------------------
! MPI version z ist the fast index
!-----------------------------------------------------------------------
#ifdef MPI
      IND=1

      DO N2=N2LOW,N2HI
      X2=(N2*F2-POSION(2))
      N2P=MOD(N2+10*GRID%NGY,GRID%NGY)

      DO N1=N1LOW,N1HI
      X1=(N1*F1-POSION(1))
      N1P=MOD(N1+10*GRID%NGX,GRID%NGX)

      NCOL=GRID%RL%INDEX(N1P,N2P)
      IF (NCOL==0) CYCLE
      IF (GRID%RL%I2(NCOL) /= N1P+1 .OR. GRID%RL%I3(NCOL) /= N2P+1) THEN
        WRITE(*,*)'internal ERROR SETYLM_AUG:',N1P+1,N2P+1,NCOL, &
           GRID%RL%I2(NCOL),N1P , GRID%RL%I3(NCOL),N2P
        STOP
      ENDIF

!OCL SCALAR
      DO N3=N3LOW,N3HI
        X3=(N3*F3-POSION(3))

        X= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
        Y= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
        Z= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)

        D=SQRT(X*X+Y*Y+Z*Z)
        ARG=(D*ARGSC)+1
        NADDR=INT(ARG)

        IF (NADDR<NPSRNL) THEN
          N3P=MOD(N3+10*GRID%NGZ,GRID%NGZ)
          NLI (IND)=1+N3P+ GRID%NGZ*(NCOL-1)

          ZZ=Z-DISZ
          YY=Y-DISY
          XX=X-DISX
! the calculation of the | R(ion)-R(mesh)+d | for displaced ions
! was done using  | R+d | = | R | + d . R/|R|
!          IF (D<1E-4_q) THEN
!            DIST(IND)=1E-4_q
!          ELSE
!            DIST(IND)=MAX(D-(DISX*X+DISY*Y+DISZ*Z)/D,1E-10_q)
!          ENDIF
! new version correct to second order
          DIST(IND)=MAX(SQRT(XX*XX+YY*YY+ZZ*ZZ),1E-15_q)

          XS(IND)  =XX/DIST(IND)
          YS(IND)  =YY/DIST(IND)
          ZS(IND)  =ZZ/DIST(IND)

          IND=IND+1
        ENDIF
      ENDDO
      ENDDO; ENDDO
#else
!-----------------------------------------------------------------------
! conventional version x is fast index
!-----------------------------------------------------------------------
      IND=1

      DO N3=N3LOW,N3HI
      X3=(N3*F3-POSION(3))
      N3P=MOD(N3+10*GRID%NGZ,GRID%NGZ)

      DO N2=N2LOW,N2HI
      X2=(N2*F2-POSION(2))
      N2P=MOD(N2+10*GRID%NGY,GRID%NGY)

      NCOL=GRID%RL%INDEX(N2P,N3P)

      DO N1=N1LOW,N1HI
      X1=(N1*F1-POSION(1))

      X= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
      Y= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
      Z= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)

      D=SQRT(X*X+Y*Y+Z*Z)
      ARG=(D*ARGSC)+1
      NADDR=INT(ARG)

      IF (NADDR<NPSRNL) THEN
        N1P=MOD(N1+10*GRID%NGX,GRID%NGX)
        NLI (IND) = N1P+(NCOL-1)*GRID%NGX+1
        IF (NLI(IND) /= 1+N1P+ GRID%NGX*( N2P + GRID%NGY*N3P )) THEN
          WRITE(*,*)'SETYLM_AUG ERROR:',N1P,N2P,N3P, NCOL
          STOP
        ENDIF
        ZZ=Z-DISZ
        YY=Y-DISY
        XX=X-DISX
! the calculation of the | R(ion)-R(mesh)+d | for displaced ions
! was done using  | R+d | = | R | + d . R/|R|
!        IF (D<1E-4_q) THEN
!          DIST(IND)=1E-4_q
!        ELSE
!          DIST(IND)=MAX(D-(DISX*X+DISY*Y+DISZ*Z)/D,1E-10_q)
!        ENDIF
! new version correct to second order
        DIST(IND)=MAX(SQRT(XX*XX+YY*YY+ZZ*ZZ),1E-10_q)

        XS(IND)  =XX/DIST(IND)
        YS(IND)  =YY/DIST(IND)
        ZS(IND)  =ZZ/DIST(IND)

        IND=IND+1
      ENDIF
      ENDDO; ENDDO; ENDDO
#endif
!-----------------------------------------------------------------------
!  compare maximum index with INDMAX
!-----------------------------------------------------------------------
      INDMAX=IND-1
      IF (INDMAX>IRMAX) THEN
        WRITE(*,*) &
     &  'internal ERROR: DEPLE:  IRDMAX must be increased to',INT(INDMAX*1.1)
        STOP
      ENDIF
      CALL SETYLM(LYDIM,INDMAX,YLM,XS,YS,ZS)
! spin spirals
! phase shifts in SETDIJ
      DO IND=1,INDMAX
         XS(IND)=XS(IND)*DIST(IND)
         YS(IND)=YS(IND)*DIST(IND)
         ZS(IND)=ZS(IND)*DIST(IND)
      ENDDO

      RETURN
    END SUBROUTINE


!************************ SUBROUTINE MAXL ******************************
!
! calculate the maximum L quantum number for the augmentation charges
! (required to allocate the tables for e.g. the spherical harmonics)
! this quantity differs for US and PAW potentials
! for the former it is MAXL, whereas for the later it is 2*MAXL
!
!***********************************************************************

   FUNCTION MAXL_AUG(NTYP,P)
      USE prec
      USE pseudo
      IMPLICIT NONE
      INTEGER MAXL_AUG,NTYP
      TYPE (potcar) P(NTYP)
   ! local varibale
      INTEGER I,NT,LTMP,CHANNELS

      MAXL_AUG=0

      DO NT=1,NTYP
         CHANNELS=P(NT)%LMAX
         LTMP=0
         DO I=1,CHANNELS
            LTMP=MAX( P(NT)%LPS(I),LTMP )
         ENDDO
         ! paw requires 2*L for the augmentation charges
         IF ( ASSOCIATED( P(NT)%QPAW) ) THEN
            LTMP=LTMP*2
         ENDIF
         MAXL_AUG=MAX( MAXL_AUG, LTMP)
      END DO

    END FUNCTION MAXL_AUG

!************************ SUBROUTINE MAXL1 *****************************
!
! calculate the maximum L quantum number for one particular type
!
!***********************************************************************

   FUNCTION MAXL1(P)
      USE prec
      USE pseudo
      IMPLICIT NONE
      INTEGER MAXL1
      TYPE (potcar) P
   ! local varibale
      INTEGER I,LTMP,CHANNELS

      MAXL1=0

      CHANNELS=P%LMAX
      LTMP=0
      DO I=1,CHANNELS
         LTMP=MAX( P%LPS(I),LTMP )
      ENDDO
      MAXL1=LTMP

    END FUNCTION MAXL1

!************************ SUBROUTINE MAXL ******************************
!
! calculate the maximum L quantum number found among all pseudopotentials
!
!***********************************************************************

   FUNCTION MAXL(NTYP,P)
      USE prec
      USE pseudo
      IMPLICIT NONE
      INTEGER MAXL,NTYP
      TYPE (potcar) P(NTYP)
   ! local varibale
      INTEGER I,NT,LTMP,CHANNELS

      MAXL=0

      DO NT=1,NTYP
         CHANNELS=P(NT)%LMAX
         LTMP=0
         DO I=1,CHANNELS
            LTMP=MAX( P(NT)%LPS(I),LTMP )
         ENDDO
         MAXL=MAX( MAXL, LTMP)
      END DO

    END FUNCTION MAXL

!************************ SUBROUTINE SETDEP ****************************
!
! this subroutine interpolates the augmentation charge on the grid
! around one ion using a cubic spline interpolation
! result is returned in DEP
!***********************************************************************

    SUBROUTINE SETDEP(QDEP,PSDMAX,NPSRNL,OMEGA,INDMAX,DIST,DEP)
      USE prec
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)


      DIMENSION QDEP(NPSRNL,5)
      DIMENSION DIST(INDMAX)
      DIMENSION DEP(INDMAX)

      FAKT= OMEGA

      ARGSC=NPSRNL/PSDMAX
!DIR$ IVDEP
!OCL NOVREC
      DO IND=1,INDMAX
        I  =MIN(INT(DIST(IND)*ARGSC)+1,NPSRNL-1)
        REM=DIST(IND)-QDEP(I,1)
        DEP(IND)=(QDEP(I,2)+REM*(QDEP(I,3)+ &
     &               REM*(QDEP(I,4)+REM*QDEP(I,5))))*FAKT
      ENDDO

      RETURN
    END SUBROUTINE SETDEP


!************************ SUBROUTINE SETDIJ  ****************************
!
! this subroutine calculates the corrections DION(I,J) corresponding
! to the integral over the augmentation-holes * total local potential
! as input it requires the total potential CVTOT (6.16) (6.12)
! the subroutine also sets up CQIJ(I,J) (i.e. Q(i,j) in thesis)
! DISX,DISY,DISZ are additional displacements of the ions
! small changes required
!
!***********************************************************************


    SUBROUTINE SETDIJ(WDES, GRIDC_,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO, LOVERL, &
        LMDIM,CDIJ,CQIJ, CVTOT_, IRDMAA,IRDMAX)
      USE prec
      USE pseudo
      USE poscar
      USE mpimy
      USE mgrid
      USE lattice
      USE wave
      USE asa
      USE paw
      USE constant
      USE us, ONLY: SETDIJ_
      IMPLICIT NONE

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d),TARGET  :: GRIDC_,GRIDUS
      TYPE (grid_3d),POINTER :: GRIDC
      TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES

      INTEGER  IRDMAX      ! allocation required for augmentation
      INTEGER  IRDMAA      ! actual maximum augmentation index
      INTEGER  LMDIM
      OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
               CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q),TARGET :: CVTOT_(GRIDC_%MPLWV,WDES%NCDIJ)
      LOGICAL  LOVERL
      REAL(q)  DISPL(3,T_INFO%NIONS)

      DISPL=0
      
      CALL SETDIJ_(WDES, GRIDC_,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO, LOVERL, &
        LMDIM,CDIJ,CQIJ, CVTOT_, .TRUE., IRDMAA,IRDMAX, DISPL)
    END SUBROUTINE SETDIJ


!
    SUBROUTINE SETDIJ_(WDES, GRIDC_,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO, LOVERL, &
        LMDIM,CDIJ,CQIJ, CVTOT_, LDIAGONAL_TERMS, IRDMAA,IRDMAX, DISPL)
      USE prec
      USE pseudo
      USE poscar
      USE mpimy
      USE mgrid
      USE lattice
      USE wave
      USE asa
      USE paw
      USE constant
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d),TARGET  :: GRIDC_,GRIDUS
      TYPE (grid_3d),POINTER :: GRIDC
      TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES

      INTEGER  IRDMAX      ! allocation required for augmentation
      INTEGER  IRDMAA      ! actual maximum augmentation index
      OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
               CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      RGRID,TARGET :: CVTOT_(DIMREAL(GRIDC_%MPLWV),WDES%NCDIJ)
      LOGICAL  LOVERL
      REAL(q)  DISPL(3,T_INFO%NIONS)
      LOGICAL  LDIAGONAL_TERMS
                           ! add atomic reference diagonal terms
!  work arrays
      LOGICAL LADDITIONAL
      REAL(q)   DLM(256)
      REAL(q)   ,ALLOCATABLE ::   DIST(:),DEP(:),POT(:),YLM(:,:)
      INTEGER,ALLOCATABLE ::   NLI(:)
      RGRID,POINTER :: CVTOT(:),CWORK(:)
      REAL(q) QVEC(3),QR
      REAL(q),ALLOCATABLE :: XS(:),YS(:),ZS(:)
      INTEGER :: LYDIM, LMYDIM, ISP, ISP_INC, NI, NIP, LOW, NT, LYMAX, MPLOW, & 
                 LL, LLP, LMP, LM, N, INDMAX, LDEP_INDEX, L, LP, M, MP, & 
                 INDYLM, INDPYL, IND, QDEP_FOCK_INDEX
      TYPE (potcar), POINTER :: PP
      INTEGER, EXTERNAL :: ONE_CENTER_NMAX_FOCKAE
#ifdef MPI
      ! mind: in the MPI version CTMP holds all elements
      ! to achieve good load balancing CDIJ is calculated locally on all nodes,
      ! merged, and finally distributed among nodes
      OVERLAP,ALLOCATABLE :: CTMP(:,:,:,:)

      PROFILING_START('setdij_')

      ALLOCATE(CTMP(LMDIM,LMDIM,T_INFO%NIONS,WDES%NCDIJ))
#else
#define CTMP CDIJ
      PROFILING_START('setdij_')
#endif
      ! LADDITIONAL uses an even finer grid for
      ! calculating the augmentation charges
      LADDITIONAL=(GRIDUS%NGX/=GRIDC_%NGX) .OR. &
                  (GRIDUS%NGY/=GRIDC_%NGY) .OR. &
                  (GRIDUS%NGZ/=GRIDC_%NGZ)
      CTMP  =0
      IRDMAA=0
      
 overl: IF (LOVERL) THEN
      ! storage convention to (total,magnetization)
      CALL RL_FLIP(CVTOT_(1,1),GRIDC_,WDES%NCDIJ,.FALSE.)

      LYDIM=MAXL_AUG(T_INFO%NTYP,P)
      LMYDIM=(LYDIM+1)**2          ! number of lm pairs

      ALLOCATE( DIST(IRDMAX),DEP(IRDMAX),POT(IRDMAX), &
     &          YLM(IRDMAX,LMYDIM),NLI(IRDMAX))
      IF (LADDITIONAL) THEN
         ALLOCATE(CVTOT(DIMREAL(GRIDUS%MPLWV)),CWORK(DIMREAL(GRIDC_%MPLWV)))
      ENDIF
      
      ALLOCATE(XS(IRDMAX),YS(IRDMAX),ZS(IRDMAX))
      IF (WDES%LSPIRAL) THEN
      ! Take QSPIRAL from direct to cartesian coordinates
         QVEC(1)=WDES%QSPIRAL(1)*LATT_CUR%B(1,1)+WDES%QSPIRAL(2)*LATT_CUR%B(1,2)+WDES%QSPIRAL(3)*LATT_CUR%B(1,3)
         QVEC(2)=WDES%QSPIRAL(1)*LATT_CUR%B(2,1)+WDES%QSPIRAL(2)*LATT_CUR%B(2,2)+WDES%QSPIRAL(3)*LATT_CUR%B(2,3)
         QVEC(3)=WDES%QSPIRAL(1)*LATT_CUR%B(3,1)+WDES%QSPIRAL(2)*LATT_CUR%B(3,2)+WDES%QSPIRAL(3)*LATT_CUR%B(3,3)
      ENDIF
 
 spin:DO ISP=1,WDES%NCDIJ

      IF (LADDITIONAL) THEN
         RINPL=1._q/GRIDC_%NPLWV
         CALL RL_ADD(CVTOT_(1,ISP),RINPL,CWORK,0.0_q,CWORK,GRIDC_)
         CALL FFT3D(CWORK(1),GRIDC_,-1)

         CVTOT=0
         CALL CPB_GRID(GRIDUS,GRIDC_,C_TO_US,CWORK(1),CVTOT(1))
         CALL FFT3D(CVTOT(1),GRIDUS,1)
         GRIDC => GRIDUS
      ELSE
         CVTOT => CVTOT_(:,ISP)
         GRIDC => GRIDC_
      ENDIF
      RINPL=1._q/GRIDC%NPLWV

!=======================================================================
! loop over all ions
!=======================================================================
      ion: DO NI=1,T_INFO%NIONS
      NT=T_INFO%ITYP(NI)
      PP=>PP_POINTER(P, NI, NT)

      ! for this ion (this type of ion) no depletion charge
      IF (PP%PSDMAX==0) CYCLE
!-----------------------------------------------------------------------
! calculate the spherical harmonics YLM and the distance DIST  between
! grid-points and central atom (DEP and POT are work-arrays)
!-----------------------------------------------------------------------
      LYMAX=MAXL1(PP)
      IF ( ASSOCIATED(PP%QPAW) ) THEN
! old VASP version: in paw method we truncate the augmentation charge at L=4
!         LYMAX=MIN(4,LYMAX*2)
         LYMAX=LYMAX*2
      ENDIF

      CALL SETYLM_AUG(GRIDC,LATT_CUR,T_INFO%POSION(1,NI),PP%PSDMAX,NPSRNL, &
     &        LMYDIM,LYMAX,YLM(1,1),IRDMAX,INDMAX, &
     &        DISPL(1,NI),DISPL(2,NI), DISPL(3,NI),DIST(1),NLI(1),XS(1),YS(1),ZS(1))

      IRDMAA=MAX(IRDMAA,INDMAX)

      DO N=1,INDMAX
         POT(N)=CVTOT(NLI(N))
      ENDDO
!=======================================================================
! US-PP
!=======================================================================
  lpaw: IF ( .NOT. ASSOCIATED(PP%QPAW) ) THEN
      LDEP_INDEX=1

    ! loop over all channels (l,epsilon)
      LM =1
      l_loop:  DO L =1,PP%LMAX
      LMP=LM
      lp_loop: DO LP=L,PP%LMAX
      IF (PP%NDEP(L,LP)==0) GOTO 510

      CALL SETDEP(PP%QDEP(1,1,LDEP_INDEX),PP%PSDMAX,NPSRNL, &
     &            LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))
      LDEP_INDEX=LDEP_INDEX+ABS(PP%NDEP(L,LP))

    ! quantum numbers l and lp of these two channels
      LL =PP%LPS(L )
      LLP=PP%LPS(LP)

    ! loop over all m mp
      m_loop:  DO M=1,2*LL+1
      MPLOW=1
      IF (L==LP) MPLOW=M
      mp_loop: DO MP=MPLOW,2*LLP+1

      !   calculate the indices into the array containing the spherical
      !   harmonics
         INDYLM =LL**2   +M
         INDPYL =LLP**2  +MP

         SUM=0

      !   sum over all r
!DIR$ IVDEP
!OCL NOVREC
         DO IND=1,INDMAX
            SUM=SUM+POT(IND)*DEP(IND)*YLM(IND,INDYLM)*YLM(IND,INDPYL)
         ENDDO

         !   add to array CTMP and make symmetric
         CTMP(LM+M-1,LMP+MP-1,NI,ISP)=SUM*RINPL
         CTMP(LMP+MP-1,LM+M-1,NI,ISP)=SUM*RINPL

      ENDDO mp_loop
      ENDDO m_loop

 510  LMP=LMP+2*LLP+1
      ENDDO lp_loop
      LM =LM +2*LL +1
      ENDDO l_loop
   ELSE lpaw
!=======================================================================
! PAW approach
! calculate first the integral int V Y(L,M) Q(L)
! (Q(L) are the L dependent compensation charges in the PAW method)
! around one atom
! then transform to  L,L',M,M' use Clebsch-Gordan coefficients
!=======================================================================
      IF (WDES%LSPIRAL .AND. ISP==2) THEN
! V_x -> cos(qr)V_x - sin(qr)V_y
! corresponds to the multiplication of V_12 and V_21 in the spinor 
! representation of the potential by exp(-iqr) and exp(+iqr), respectively
      DLM=0
      DO L =0,LYMAX
         CALL SETDEP(PP%QDEP(1,1,L),PP%PSDMAX,NPSRNL, &
              LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))
         DO M=1,(L*2)+1
            INDYLM =L*L  +M
            SUM=0
!DIR$ IVDEP
!OCL NOVREC
            DO IND=1,INDMAX
            ! calculate q dot r, both given in cartesian coordinates
               QR=TPI*(QVEC(1)*XS(IND)+QVEC(2)*YS(IND)+QVEC(3)*ZS(IND))
               SUM=SUM+DEP(IND)*YLM(IND,INDYLM)* &
              &     (CVTOT_(NLI(IND),2)*COS(QR)-CVTOT_(NLI(IND),3)*SIN(QR))
            ENDDO
            DLM(INDYLM)=SUM*RINPL
         ENDDO
         ! WRITE(0,'("DLM",I2,20F7.4)') L,(DLM(L**2+M),M=1,(L*2)+1)
      ENDDO
      CALL CALC_DLLMM( CTMP(:,:,NI,ISP), DLM, PP)
      ENDIF
      
      IF (WDES%LSPIRAL .AND. ISP==3) THEN
! V_y -> cos(qr)V_y + sin(qr)V_x
! corresponds to the multiplication of V_12 and V_21 in the spinor 
! representation of the potential by exp(-iqr) and exp(+iqr), respectively
      DLM=0
      DO L =0,LYMAX
         CALL SETDEP(PP%QDEP(1,1,L),PP%PSDMAX,NPSRNL, &
              LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))
         DO M=1,(L*2)+1
            INDYLM =L*L  +M
            SUM=0
!DIR$ IVDEP
!OCL NOVREC
            DO IND=1,INDMAX
            ! calculate q dot r, both given in cartesian coordinates
               QR=TPI*(QVEC(1)*XS(IND)+QVEC(2)*YS(IND)+QVEC(3)*ZS(IND))
               SUM=SUM+DEP(IND)*YLM(IND,INDYLM)* &
              &     (CVTOT_(NLI(IND),3)*COS(QR)+CVTOT_(NLI(IND),2)*SIN(QR))
            ENDDO
            DLM(INDYLM)=SUM*RINPL
         ENDDO
      !   WRITE(0,'("DLM",I2,20F7.4)') L,(DLM(L**2+M),M=1,(L*2)+1)
      ENDDO
      CALL CALC_DLLMM( CTMP(:,:,NI,ISP), DLM, PP)
      ENDIF
      
      IF (.NOT.WDES%LSPIRAL .OR. ISP==1 .OR. ISP==4) THEN
! no phase factor for total charge and m_z or if LSPIRAL=.FALSE.
      DLM=0
      DO L =0,LYMAX
         CALL SETDEP(PP%QDEP(1,1,L),PP%PSDMAX,NPSRNL, &
              LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))
         DO M=1,(L*2)+1
            INDYLM =L*L  +M
            SUM=0
            SUMN=0
!DIR$ IVDEP
!OCL NOVREC
            DO IND=1,INDMAX
               SUM=SUM+POT(IND)*DEP(IND)*YLM(IND,INDYLM)
               SUMN=SUMN+DEP(IND)*YLM(IND,INDYLM)
            ENDDO
            DLM(INDYLM)=SUM*RINPL
         ENDDO
       !  IF (L<=1) WRITE(0,'("DLM",I2,10F7.4)') L,(1E6*DLM(L**2+M),M=1,(L*2)+1)
      ENDDO
      CALL CALC_DLLMM( CTMP(:,:,NI,ISP), DLM, PP)
!=======================================================================
! PAW contributions from accurate augmentation 
! currently spin spirals are no supported
!=======================================================================
      onecAE: IF (ONE_CENTER_NMAX_FOCKAE()>0) THEN

      DLM=0
      NMAX_FOCKAE=SIZE(PP%QPAW_FOCK,4)
      LMAX_FOCKAE=SIZE(PP%QPAW_FOCK,3)-1
      ! proper indexing of QDEP_FOCK is nasty see fast_aug.F
      QDEP_FOCK_INDEX=SIZE(PP%QDEP_FOCK,3)-NMAX_FOCKAE*(LMAX_FOCKAE+1)-1
      DO NAE=1,NMAX_FOCKAE
      DO L  =0,LMAX_FOCKAE
         QDEP_FOCK_INDEX=QDEP_FOCK_INDEX+1

         CALL SETDEP(PP%QDEP_FOCK(1,1,QDEP_FOCK_INDEX),PP%PSDMAX,NPSRNL, &
              LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))
         DO M=1,(L*2)+1
            INDYLM =L*L  +M
            SUM=0
!DIR$ IVDEP
!OCL NOVREC
            DO IND=1,INDMAX
               SUM=SUM+POT(IND)*DEP(IND)*YLM(IND,INDYLM)
            ENDDO
            DLM(INDYLM)=SUM*RINPL
         ENDDO
       !  IF (L<=1) WRITE(0,'("DLM",I2,10F7.4)') L,(1E6*DLM(L**2+M),M=1,(L*2)+1)
      ENDDO
!     WRITE(*,*) '1',DLM(1)
      CALL CALC_DLLMM_FOCK( CTMP(:,:,NI,ISP), DLM, PP, NAE)
      ENDDO

      IF (QDEP_FOCK_INDEX /= SIZE(PP%QDEP_FOCK,3)-1) THEN
         WRITE(0,*) 'internal error in SETDIJ_FOCKAE: size of QDEP_FOCK seems to be incorrect',SIZE(PP%AUG_SOFT,2),QDEP_FOCK_INDEX,SIZE(PP%QDEP_FOCK,3)-1
         STOP
      ENDIF

      ENDIF onecAE

      ENDIF
      ENDIF lpaw
!=======================================================================
      ENDDO ion
!-----------------------------------------------------------------------
      ENDDO spin

      DEALLOCATE(DIST,DEP,POT,YLM,NLI,XS,YS,ZS)
      IF (LADDITIONAL) DEALLOCATE(CVTOT,CWORK)
! reduce CTMP
#ifdef realmode
      CALLMPI( M_sum_d(GRIDC%COMM, CTMP, LMDIM*LMDIM*T_INFO%NIONS*WDES%NCDIJ))
#else
      CALLMPI( M_sum_d(GRIDC%COMM, CTMP, LMDIM*LMDIM*T_INFO%NIONS*WDES%NCDIJ*2))
#endif
      ! back to spinor representation
      CALL RL_FLIP(CVTOT_(1,1),GRIDC_,WDES%NCDIJ,.TRUE.)

      ENDIF overl
!-----------------------------------------------------------------------
! now set up CQIJ and add diagonal part to CDIJ
! find blocks with same quantum number L
! the routine is somewhat complicated
! only terms with same quantum numbers L L' and M M' are non-zero
!-----------------------------------------------------------------------
      ion2: DO NI=1,T_INFO%NIONS
      NIP=NI_LOCAL(NI, WDES%COMM_INB)
      IF (NIP==0) CYCLE ion2
        DO ISP=1,WDES%NCDIJ
           CQIJ(:,:,NIP,ISP)=0
#ifdef MPI
           CDIJ(:,:,NIP,ISP)=CTMP(:,:,NI,ISP)
#endif
        ENDDO

        LOW=1
        LM =1 
        NT=T_INFO%ITYP(NI)
        PP=>PP_POINTER(P, NI, NT)

        IF (WDES%LNONCOLLINEAR) THEN 
           ISP_INC=3
           FAKT   =2
        ELSE
           ISP_INC=1
           FAKT=WDES%ISPIN
        ENDIF

        block: DO
        IF (PP%LMAX==0) EXIT
           LL=PP%LPS(LOW)
           ! search block with same L
           DO LHI=LOW,PP%LMAX
              IF (LL/=PP%LPS(LHI)) EXIT
           ENDDO
           LHI=LHI-1
           MMAX=2*LL+1

           DO ISP=1,WDES%NCDIJ,ISP_INC
           DO L =LOW,LHI
           DO LP=LOW,LHI
           DO M =0,MMAX-1
              CQIJ(LM+(L-LOW)*MMAX+M,LM+(LP-LOW)*MMAX+M,NIP,ISP)=PP%QION(L,LP)
           ENDDO
           ENDDO
           ENDDO
           ENDDO

           IF (LDIAGONAL_TERMS) THEN
           ISP=1
           DO L =LOW,LHI
           DO LP=LOW,LHI
           DO M =0,MMAX-1
              CDIJ(LM+(L-LOW)*MMAX+M,LM+(LP-LOW)*MMAX+M,NIP,ISP)=PP%DION(L,LP)*FAKT &
                   &   +CDIJ(LM+(L-LOW)*MMAX+M,LM+(LP-LOW)*MMAX+M,NIP,ISP)
           ENDDO
           ENDDO
           ENDDO
           ENDIF
        ! set new LOW value and LM value and goon
           LM=LM+(LHI-LOW+1)*MMAX
           LOW=LHI+1
           IF (LOW > PP%LMAX) EXIT block
        ENDDO block

        IF (T_INFO%VCA(NT)/=1.0) THEN
           DO ISP=1,WDES%NCDIJ
              CDIJ(:,:,NIP,ISP)=CDIJ(:,:,NIP,ISP)*T_INFO%VCA(NT)
              CQIJ(:,:,NIP,ISP)=CQIJ(:,:,NIP,ISP)*T_INFO%VCA(NT)
           ENDDO
        ENDIF

!       IF (GRIDC%COMM%NODE_ME == GRIDC%COMM%IONODE) THEN
!       DO ISP=1,WDES%NCDIJ
!       WRITE(*,*)'ion',NI,NIP,ISP
!       DO LP=1,P(1)%LMMAX
!         WRITE(*,'(16(F10.3,1X))') (CDIJ(L,LP,NIP,ISP),L=1,MIN(8,P(1)%LMMAX))
!       ENDDO
!       WRITE(*,*)
!       ENDDO
!       ENDIF

      ENDDO ion2
      ! go to spin up and down presentation for CDIJ
      CALL US_FLIP(WDES, LMDIM, CDIJ, .TRUE., .TRUE.)

#ifdef MPI
      DEALLOCATE(CTMP)
#else
#undef CTMP
#endif
      PROFILING_STOP('setdij_')

      RETURN
    END SUBROUTINE SETDIJ_


!************************ SUBROUTINE SETDIJ_AVEC_ **********************
!
! this subroutine calculates the vector potential related 
! contributions to the non-local pseudopotential strength parameters
!
! here we have to calculate an augmentation contribution to the non-
! local strength parameter
!  D_ij = \int(r) j_ij(r) A(r) d^3 r
! where j_ij(r) is the augmentation current in each sphere for pairs
! if states ij, and A(r)
! the vector potential
! if the projectors are gauge twisted the equation becomes
!  D_ij = \int(r) j_ij(r) [A(r)-A(R)] d^3 r
!
!***********************************************************************


    SUBROUTINE SETDIJ_AVEC_(WDES, GRIDC_,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO, LOVERL, &
        LMDIM,CDIJ,AVTOT_, NONLR_S, NONL_S, IRDMAX, DISPL )
      USE prec
      USE pseudo
      USE poscar
      USE mpimy
      USE mgrid
      USE lattice
      USE wave
      USE asa
      USE paw
      USE constant
      USE nonl_high
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d),TARGET  :: GRIDC_,GRIDUS
      TYPE (grid_3d),POINTER :: GRIDC
      TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES
      INTEGER  IRDMAX      ! allocation required for augmentation
      OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      RGRID,TARGET :: AVTOT_(DIMREAL(GRIDC_%MPLWV),3)
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct)  NONL_S
      LOGICAL  LOVERL
      REAL(q)  DISPL(3,T_INFO%NIONS)
!  work arrays
      LOGICAL LADDITIONAL
      OVERLAP  ::   DLM(256),SUM
      REAL(q), ALLOCATABLE ::   DIST(:),DEP(:),YLM(:,:),XS(:),YS(:),ZS(:)
      RGRID,   ALLOCATABLE ::   POT(:)
      INTEGER, ALLOCATABLE ::   NLI(:)
      RGRID,POINTER :: AVTOT(:),CWORK(:)
      TYPE (potcar), POINTER :: PP
! CTMP holds the contribution to the strength parameter from 
! the vector field for all ions
! for the serial version, storing one ion would be sufficient
      OVERLAP,ALLOCATABLE :: CTMP(:,:,:)
!-----------------------------------------------------------------------
! allocation
!-----------------------------------------------------------------------
      ALLOCATE(CTMP(LMDIM,LMDIM,T_INFO%NIONS))

      ! LADDITIONAL uses an even finer grid for
      ! calculating the augmentation charges
      LADDITIONAL=(GRIDUS%NGX/=GRIDC_%NGX) .OR. &
                  (GRIDUS%NGY/=GRIDC_%NGY) .OR. &
                  (GRIDUS%NGZ/=GRIDC_%NGZ)
      CTMP  =0
      
 overl: IF (LOVERL) THEN
      ! go to total and magnetization presentation for CDIJ
      CALL CDIJ_FLIP(WDES, LMDIM, CDIJ, LOVERL, .FALSE.)


      LYDIM=MAXL_AUG(T_INFO%NTYP,P)+1
      LMYDIM=(LYDIM+1)**2          ! number of lm pairs

      ALLOCATE( DIST(IRDMAX),DEP(IRDMAX),POT(IRDMAX), &
                YLM(IRDMAX,LMYDIM),NLI(IRDMAX),XS(IRDMAX),YS(IRDMAX),ZS(IRDMAX))
      IF (LADDITIONAL) THEN
         ALLOCATE(AVTOT(DIMREAL(GRIDUS%MPLWV)),CWORK(DIMREAL(GRIDC_%MPLWV)))
      ENDIF

!-----------------------------------------------------------------------
! loop over all ions
!-----------------------------------------------------------------------
      IF (NONLR_S%LREAL .AND. ASSOCIATED(NONLR_S%VKPT_SHIFT)) THEN
         NONLR_S%VKPT_SHIFT=0
      ELSE IF (ASSOCIATED(NONL_S%VKPT_SHIFT)) THEN
         NONL_S%VKPT_SHIFT=0
      ENDIF

 cart:DO ICART=1,3

      IF (LADDITIONAL) THEN
         RINPL=1._q/GRIDC_%NPLWV
         CALL RL_ADD(AVTOT_(1,ICART),RINPL,CWORK,0.0_q,CWORK,GRIDC_)
         CALL FFT3D(CWORK(1),GRIDC_,-1)

         AVTOT=0
         CALL CPB_GRID(GRIDUS,GRIDC_,C_TO_US,CWORK(1),AVTOT(1))
         CALL FFT3D(AVTOT(1),GRIDUS,1)
         GRIDC => GRIDUS
      ELSE
         AVTOT => AVTOT_(:,ICART)
         GRIDC => GRIDC_
      ENDIF
      RINPL=1._q/GRIDC%NPLWV
!-----------------------------------------------------------------------
! loop over all ions
!-----------------------------------------------------------------------
      ion: DO NI=1,T_INFO%NIONS
      NT=T_INFO%ITYP(NI)
      PP=>PP_POINTER(P, NI, NT)
      ! for this ion (this type of ion) no depletion charge
      IF (PP%PSDMAX==0 .OR.  .NOT. ASSOCIATED(PP%JPAW)) CYCLE
!-----------------------------------------------------------------------
! calculate the spherical harmonics YLM and the distance DIST  between
! grid-points and central atom (DEP and POT are work-arrays)
!-----------------------------------------------------------------------
      ! nabla operator increases the maximum L one-center number by one
      LYMAX=MAXL1(PP)*2+1

      CALL SETYLM_AUG(GRIDC,LATT_CUR,T_INFO%POSION(1,NI),PP%PSDMAX,NPSRNL, &
     &        LMYDIM,LYMAX,YLM(1,1),IRDMAX,INDMAX, &
     &        DISPL(1,NI),DISPL(2,NI), DISPL(3,NI),DIST(1),NLI(1),XS(1),YS(1),ZS(1))

      DO N=1,INDMAX
         POT(N)=AVTOT(NLI(N))
      ENDDO

      DLM=0
      DO L =0,LYMAX
         CALL SETDEP(PP%QDEP(1,1,L),PP%PSDMAX,NPSRNL, &
              LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))
         DO M=1,(L*2)+1
            INDYLM =L*L  +M
            SUM=0
!DIR$ IVDEP
!OCL NOVREC
            DO IND=1,INDMAX
               SUM=SUM+POT(IND)*DEP(IND)*YLM(IND,INDYLM)
            ENDDO
            DLM(INDYLM)=SUM*RINPL
         ENDDO
!         IF (L<=1) WRITE(*,'("AVEC",3I4,30F10.4)') ICART,NI,L,(DLM(L**2+M)*1000,M=1,(L*2)+1)
      ENDDO
      WRITE(*,'("AVEC",2I4,30F10.4)') ICART,NI,(DLM(1:4)*1000)

! important, if the projectors are phase twisted the monopole terms must not be
! included. Essentially they are handled by the phase twist or in other words
! A(r) in the sphere needs to be replaced by A(r)-A(R)
! the phase twist is passed on the NONLR_S%VKPT_SHIFT below
      IF (NONLR_S%LREAL .AND. ASSOCIATED(NONLR_S%VKPT_SHIFT)) THEN
         ! store phase shift in the non-local projectors
         NONLR_S%VKPT_SHIFT(ICART,NI) =DLM(1)
         DLM(1)=0
      ELSE IF (ASSOCIATED(NONL_S%VKPT_SHIFT)) THEN
         ! store phase shift in the non-local projectors
         NONL_S%VKPT_SHIFT(ICART,NI)  =DLM(1)
         DLM(1)=0
      ENDIF

      CALL CALC_DLLMM_AVEC( CTMP(:,:,NI), DLM, PP, ICART)
!-----------------------------------------------------------------------
      ENDDO ion
      ENDDO cart

      DEALLOCATE(DIST,DEP,POT,YLM,NLI,XS,YS,ZS)
      IF (LADDITIONAL) DEALLOCATE(AVTOT,CWORK)
!-----------------------------------------------------------------------
! reduce CTMP and bring results to all nodes
!-----------------------------------------------------------------------
#ifdef realmode
      CALLMPI( M_sum_d(GRIDC%COMM, CTMP, LMDIM*LMDIM*T_INFO%NIONS))
#else
      CALLMPI( M_sum_z(GRIDC%COMM, CTMP, LMDIM*LMDIM*T_INFO%NIONS))
#endif

      IF (NONLR_S%LREAL .AND. ASSOCIATED(NONLR_S%VKPT_SHIFT)) THEN
         CALLMPI( M_sum_d(GRIDC%COMM, NONLR_S%VKPT_SHIFT, SIZE(NONLR_S%VKPT_SHIFT)))
         NONLR_S%VKPT_SHIFT=NONLR_S%VKPT_SHIFT  *MOMTOMOM/MAGMOMTOENERGY/SQRT(4 *PI)/TPI
         CALL KARDIR( SIZE(NONLR_S%VKPT_SHIFT,2), NONLR_S%VKPT_SHIFT , LATT_CUR%A)

      ELSE IF (ASSOCIATED(NONL_S%VKPT_SHIFT)) THEN
         CALLMPI( M_sum_d(GRIDC%COMM, NONL_S%VKPT_SHIFT, SIZE(NONL_S%VKPT_SHIFT)))
         NONL_S%VKPT_SHIFT=NONL_S%VKPT_SHIFT    *MOMTOMOM/MAGMOMTOENERGY/SQRT(4 *PI)/TPI
         CALL KARDIR( SIZE(NONL_S%VKPT_SHIFT,2), NONL_S%VKPT_SHIFT , LATT_CUR%A)
      ENDIF

      ion2: DO NI=1,T_INFO%NIONS
         NIP=NI_LOCAL(NI, WDES%COMM_INB)
         IF (NIP==0) CYCLE ion2
         ! action of nabla is stored in CTMP
         ! what is lacking is the factor -i
         CDIJ(:,:,NIP,1)=CDIJ(:,:,NIP,1)-CTMP(:,:,NI)*(0.0_q,1.0_q)
      ENDDO ion2

      ! back to spinor representation for CDIJ
      CALL CDIJ_FLIP(WDES, LMDIM, CDIJ, LOVERL, .TRUE.)

     ENDIF overl

!    VCA is not implemented here

      DEALLOCATE(CTMP)

      RETURN
    END SUBROUTINE SETDIJ_AVEC_


! identical routine but with one more argument A_ONE_CENTER
! we can not have an explict interface for this routines, because
! of issues with POINTERS (AVEC) passing to f77 line routines

    SUBROUTINE SETDIJ_AVEC_ONE_CENTER(WDES, GRIDC_,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO, LOVERL, &
        LMDIM,CDIJ,AVTOT_, NONLR_S, NONL_S, IRDMAX, DISPL,  A_ONE_CENTER, LMMAX)
      USE prec
      USE pseudo
      USE poscar
      USE mpimy
      USE mgrid
      USE lattice
      USE wave
      USE asa
      USE paw
      USE constant
      USE nonl_high
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d),TARGET  :: GRIDC_,GRIDUS
      TYPE (grid_3d),POINTER :: GRIDC
      TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES
      INTEGER  IRDMAX      ! allocation required for augmentation
      OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      RGRID,TARGET :: AVTOT_(DIMREAL(GRIDC_%MPLWV),3)
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct)  NONL_S
      LOGICAL  LOVERL
      INTEGER  LMMAX
      REAL(q)  DISPL(3,T_INFO%NIONS)
      REAL(q)  A_ONE_CENTER(3,LMMAX,T_INFO%NIONS)
!  work arrays
      LOGICAL LADDITIONAL
      OVERLAP  ::   DLM(256),SUM
      REAL(q), ALLOCATABLE ::   DIST(:),DEP(:),YLM(:,:),XS(:),YS(:),ZS(:)
      RGRID,   ALLOCATABLE ::   POT(:)
      INTEGER, ALLOCATABLE ::   NLI(:)
      RGRID,POINTER :: AVTOT(:),CWORK(:)
      TYPE (potcar), POINTER :: PP
! CTMP holds the contribution to the strength parameter from 
! the vector field for all ions
! for the serial version, storing one ion would be sufficient
      OVERLAP,ALLOCATABLE :: CTMP(:,:,:)
!-----------------------------------------------------------------------
! allocation
!-----------------------------------------------------------------------
      ALLOCATE(CTMP(LMDIM,LMDIM,T_INFO%NIONS))

      ! LADDITIONAL uses an even finer grid for
      ! calculating the augmentation charges
      LADDITIONAL=(GRIDUS%NGX/=GRIDC_%NGX) .OR. &
                  (GRIDUS%NGY/=GRIDC_%NGY) .OR. &
                  (GRIDUS%NGZ/=GRIDC_%NGZ)
      CTMP  =0
      
 overl: IF (LOVERL) THEN
      ! go to total and magnetization presentation for CDIJ
      CALL CDIJ_FLIP(WDES, LMDIM, CDIJ, LOVERL, .FALSE.)


      LYDIM=MAXL_AUG(T_INFO%NTYP,P)+1
      LMYDIM=(LYDIM+1)**2          ! number of lm pairs

      ALLOCATE( DIST(IRDMAX),DEP(IRDMAX),POT(IRDMAX), &
                YLM(IRDMAX,LMYDIM),NLI(IRDMAX),XS(IRDMAX),YS(IRDMAX),ZS(IRDMAX))
      IF (LADDITIONAL) THEN
         ALLOCATE(AVTOT(DIMREAL(GRIDUS%MPLWV)),CWORK(DIMREAL(GRIDC_%MPLWV)))
      ENDIF

!-----------------------------------------------------------------------
! loop over all ions
!-----------------------------------------------------------------------
      IF (NONLR_S%LREAL .AND. ASSOCIATED(NONLR_S%VKPT_SHIFT)) THEN
         NONLR_S%VKPT_SHIFT=0
      ELSE IF (ASSOCIATED(NONL_S%VKPT_SHIFT)) THEN
         NONL_S%VKPT_SHIFT=0
      ENDIF

 cart:DO ICART=1,3

      IF (LADDITIONAL) THEN
         RINPL=1._q/GRIDC_%NPLWV
         CALL RL_ADD(AVTOT_(1,ICART),RINPL,CWORK,0.0_q,CWORK,GRIDC_)
         CALL FFT3D(CWORK(1),GRIDC_,-1)

         AVTOT=0
         CALL CPB_GRID(GRIDUS,GRIDC_,C_TO_US,CWORK(1),AVTOT(1))
         CALL FFT3D(AVTOT(1),GRIDUS,1)
         GRIDC => GRIDUS
      ELSE
         AVTOT => AVTOT_(:,ICART)
         GRIDC => GRIDC_
      ENDIF
      RINPL=1._q/GRIDC%NPLWV
!-----------------------------------------------------------------------
! loop over all ions
!-----------------------------------------------------------------------
      ion: DO NI=1,T_INFO%NIONS
      NT=T_INFO%ITYP(NI)
      PP=>PP_POINTER(P, NI, NT)
      ! for this ion (this type of ion) no depletion charge
      IF (PP%PSDMAX==0 .OR.  .NOT. ASSOCIATED(PP%JPAW)) CYCLE
!-----------------------------------------------------------------------
! calculate the spherical harmonics YLM and the distance DIST  between
! grid-points and central atom (DEP and POT are work-arrays)
!-----------------------------------------------------------------------
      ! nabla operator increases the maximum L one-center number by one
      LYMAX=MAXL1(PP)*2+1

      CALL SETYLM_AUG(GRIDC,LATT_CUR,T_INFO%POSION(1,NI),PP%PSDMAX,NPSRNL, &
     &        LMYDIM,LYMAX,YLM(1,1),IRDMAX,INDMAX, &
     &        DISPL(1,NI),DISPL(2,NI), DISPL(3,NI),DIST(1),NLI(1),XS(1),YS(1),ZS(1))

      DO N=1,INDMAX
         POT(N)=AVTOT(NLI(N))
      ENDDO

      DLM=0
      DO L =0,LYMAX
         CALL SETDEP(PP%QDEP(1,1,L),PP%PSDMAX,NPSRNL, &
              LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))
         DO M=1,(L*2)+1
            INDYLM =L*L  +M
            SUM=0
!DIR$ IVDEP
!OCL NOVREC
            DO IND=1,INDMAX
               SUM=SUM+POT(IND)*DEP(IND)*YLM(IND,INDYLM)
            ENDDO
            DLM(INDYLM)=SUM*RINPL
         ENDDO
!         IF (L<=1) WRITE(*,'("AVEC",3I4,30F10.4)') ICART,NI,L,(DLM(L**2+M)*1000,M=1,(L*2)+1)
      ENDDO

! added line
      A_ONE_CENTER(ICART,:,NI)=DLM(1:SIZE(A_ONE_CENTER,2))

      IF (NONLR_S%LREAL .AND. ASSOCIATED(NONLR_S%VKPT_SHIFT)) THEN
         NONLR_S%VKPT_SHIFT(ICART,NI) =DLM(1)
         DLM(1)=0
      ELSE IF (ASSOCIATED(NONL_S%VKPT_SHIFT)) THEN
         NONL_S%VKPT_SHIFT(ICART,NI)  =DLM(1)
         DLM(1)=0
      ENDIF
!      WRITE(*,'("AVEC",2I4,30F10.4)') ICART,NI,(DLM(1:4)*1000)

      CALL CALC_DLLMM_AVEC( CTMP(:,:,NI), DLM, PP, ICART)
!-----------------------------------------------------------------------
      ENDDO ion
      ENDDO cart

      DEALLOCATE(DIST,DEP,POT,YLM,NLI,XS,YS,ZS)
      IF (LADDITIONAL) DEALLOCATE(AVTOT,CWORK)
!-----------------------------------------------------------------------
! reduce CTMP and bring results to all nodes
!-----------------------------------------------------------------------
#ifdef realmode
      CALLMPI( M_sum_d(GRIDC%COMM, CTMP, LMDIM*LMDIM*T_INFO%NIONS))
#else
      CALLMPI( M_sum_z(GRIDC%COMM, CTMP, LMDIM*LMDIM*T_INFO%NIONS))
#endif

!added line
      CALLMPI( M_sum_d(GRIDC%COMM, A_ONE_CENTER, SIZE(A_ONE_CENTER)))

      IF (NONLR_S%LREAL .AND. ASSOCIATED(NONLR_S%VKPT_SHIFT)) THEN
         CALLMPI( M_sum_d(GRIDC%COMM, NONLR_S%VKPT_SHIFT, SIZE(NONLR_S%VKPT_SHIFT)))
         NONLR_S%VKPT_SHIFT=NONLR_S%VKPT_SHIFT  *MOMTOMOM/MAGMOMTOENERGY/SQRT(4 *PI)/TPI
         CALL KARDIR( SIZE(NONLR_S%VKPT_SHIFT,2), NONLR_S%VKPT_SHIFT , LATT_CUR%A)
         

      ELSE IF (ASSOCIATED(NONL_S%VKPT_SHIFT)) THEN
         CALLMPI( M_sum_d(GRIDC%COMM, NONL_S%VKPT_SHIFT, SIZE(NONL_S%VKPT_SHIFT)))
         NONL_S%VKPT_SHIFT=NONL_S%VKPT_SHIFT    *MOMTOMOM/MAGMOMTOENERGY/SQRT(4 *PI)/TPI
         CALL KARDIR( SIZE(NONL_S%VKPT_SHIFT,2), NONL_S%VKPT_SHIFT , LATT_CUR%A)
      ENDIF

      ion2: DO NI=1,T_INFO%NIONS
         NIP=NI_LOCAL(NI, WDES%COMM_INB)
         IF (NIP==0) CYCLE ion2
         ! action of nabla is stored in JPAW
         ! what is lacking is the factor -i
         CDIJ(:,:,NIP,1)=CDIJ(:,:,NIP,1)-CTMP(:,:,NI)*(0.0_q,1.0_q)
      ENDDO ion2

      ! back to spinor representation for CDIJ
      CALL CDIJ_FLIP(WDES, LMDIM, CDIJ, LOVERL, .TRUE.)

     ENDIF overl

!    VCA is not implemented here

      DEALLOCATE(CTMP)

      RETURN
    END SUBROUTINE SETDIJ_AVEC_ONE_CENTER

!************************ SUBROUTINE US_FLIP ***************************
!
! rearranges the storage mode for spin components of array CRHODE:
! given crhode_up and crhode_down on input the quantities
! (crhode_up+crhode_down) and (crhode_up-crhode_down) = total charge
! and magnetization are returned
! also the reverse operation is possible if setting LBACK=.TRUE.
!
!***********************************************************************

    SUBROUTINE US_FLIP(WDES, LMDIM, CRHODE, LOVERL, LBACK)
      USE prec
      USE wave

      IMPLICIT NONE

      TYPE (wavedes)     WDES
      LOGICAL LOVERL, LBACK
      INTEGER LMDIM,LMMAXC,NT,NI,NIS,L,LP,LMBASE
      REAL(q) FAC
      OVERLAP :: CQU,CQD,C01,C10
      OVERLAP :: C11,C00,CX,CY,CZ
      OVERLAP :: CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)


      IF (.NOT.LOVERL) RETURN

      IF (WDES%NCDIJ==2 ) THEN
!=======================================================================
         FAC=1._q
         IF (LBACK) FAC=0.5_q
      
         NIS=1
         LMBASE=0
         DO NT=1,WDES%NTYP
            LMMAXC=WDES%LMMAX(NT)
            IF (LMMAXC==0) GOTO 100

            DO NI=NIS,WDES%NITYP(NT)+NIS-1

            DO L=1,LMMAXC
               DO LP=1,LMMAXC
                  CQU=CRHODE(L,LP,NI,1)
                  CQD=CRHODE(L,LP,NI,2)
                  CRHODE(L,LP,NI,1)=FAC*(CQU+CQD)
                  CRHODE(L,LP,NI,2)=FAC*(CQU-CQD)
               ENDDO
            ENDDO
         
            LMBASE= LMMAXC+LMBASE
            ENDDO

 100     NIS = NIS+WDES%NITYP(NT)
         ENDDO
      ELSE IF (  WDES%NCDIJ==4 .AND. .NOT. LBACK) THEN
!=======================================================================
         NIS=1
         LMBASE=0
         DO NT=1,WDES%NTYP
            LMMAXC=WDES%LMMAX(NT)
            IF (LMMAXC==0) GOTO 200

            DO NI=NIS,WDES%NITYP(NT)+NIS-1

            DO L=1,LMMAXC
               DO LP=1,LMMAXC
                  C00=CRHODE(L,LP,NI,1)
                  C01=CRHODE(L,LP,NI,2)
                  C10=CRHODE(L,LP,NI,3)
                  C11=CRHODE(L,LP,NI,4)

                  CRHODE(L,LP,NI,1)= C00+C11
                  CRHODE(L,LP,NI,2)= C01+C10
                  CRHODE(L,LP,NI,3)=(C01-C10)*(0._q,1._q)
                  CRHODE(L,LP,NI,4)= C00-C11             
               ENDDO
            ENDDO
         
            LMBASE= LMMAXC+LMBASE
            ENDDO

 200     NIS = NIS+WDES%NITYP(NT)
         ENDDO
      ELSE IF (  WDES%NCDIJ==4 .AND. LBACK) THEN
!=======================================================================
         FAC=0.5_q
         NIS=1
         LMBASE=0
         typ:  DO NT=1,WDES%NTYP
            LMMAXC=WDES%LMMAX(NT)
            IF (LMMAXC==0) GOTO 300

            ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1

            DO L=1,LMMAXC
               DO LP=1,LMMAXC
                  C00=CRHODE(L,LP,NI,1)
                  CX =CRHODE(L,LP,NI,2)
                  CY =CRHODE(L,LP,NI,3)
                  CZ =CRHODE(L,LP,NI,4)

                  CRHODE(L,LP,NI,1)= (C00+CZ)*FAC
                  CRHODE(L,LP,NI,2)= (CX-CY*(0._q,1._q))*FAC
                  CRHODE(L,LP,NI,3)= (CX+CY*(0._q,1._q))*FAC
                  CRHODE(L,LP,NI,4)= (C00-CZ)*FAC
               ENDDO
            ENDDO
         
            LMBASE= LMMAXC+LMBASE
            ENDDO ion

 300     NIS = NIS+WDES%NITYP(NT)
         ENDDO typ
      ENDIF

    END SUBROUTINE US_FLIP

!************************ SUBROUTINE US_FLIP_CMPLX  ********************
!
! rearranges the storage mode for spin components of array CRHODE:
! given crhode_up and crhode_down on input the quantities
! (crhode_up+crhode_down) and (crhode_up-crhode_down) = total charge
! and magnetization are returned
! also the reverse operation is possible if setting LBACK=.TRUE.
!
!***********************************************************************

    SUBROUTINE US_FLIP_CMPLX(WDES, LMDIM, CRHODE, LOVERL, LBACK)
      USE prec
      USE wave

      IMPLICIT NONE

      TYPE (wavedes)     WDES
      LOGICAL LOVERL, LBACK
      INTEGER LMDIM,LMMAXC,NT,NI,NIS,L,LP,LMBASE
      REAL(q) FAC
      COMPLEX(q) :: CQU,CQD,C01,C10
      COMPLEX(q) :: C11,C00,CX,CY,CZ
      COMPLEX(q) :: CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)


      IF (.NOT.LOVERL) RETURN

      IF (WDES%NCDIJ==2 ) THEN
!=======================================================================
         FAC=1._q
         IF (LBACK) FAC=0.5_q
      
         NIS=1
         LMBASE=0
         DO NT=1,WDES%NTYP
            LMMAXC=WDES%LMMAX(NT)
            IF (LMMAXC==0) GOTO 100

            DO NI=NIS,WDES%NITYP(NT)+NIS-1

            DO L=1,LMMAXC
               DO LP=1,LMMAXC
                  CQU=CRHODE(L,LP,NI,1)
                  CQD=CRHODE(L,LP,NI,2)
                  CRHODE(L,LP,NI,1)=FAC*(CQU+CQD)
                  CRHODE(L,LP,NI,2)=FAC*(CQU-CQD)
               ENDDO
            ENDDO
         
            LMBASE= LMMAXC+LMBASE
            ENDDO

 100     NIS = NIS+WDES%NITYP(NT)
         ENDDO
      ELSE IF (  WDES%NCDIJ==4 .AND. .NOT. LBACK) THEN
!=======================================================================
         NIS=1
         LMBASE=0
         DO NT=1,WDES%NTYP
            LMMAXC=WDES%LMMAX(NT)
            IF (LMMAXC==0) GOTO 200

            DO NI=NIS,WDES%NITYP(NT)+NIS-1

            DO L=1,LMMAXC
               DO LP=1,LMMAXC
                  C00=CRHODE(L,LP,NI,1)
                  C01=CRHODE(L,LP,NI,2)
                  C10=CRHODE(L,LP,NI,3)
                  C11=CRHODE(L,LP,NI,4)

                  CRHODE(L,LP,NI,1)= C00+C11
                  CRHODE(L,LP,NI,2)= C01+C10
                  CRHODE(L,LP,NI,3)=(C01-C10)*(0._q,1._q)
                  CRHODE(L,LP,NI,4)= C00-C11             
               ENDDO
            ENDDO
         
            LMBASE= LMMAXC+LMBASE
            ENDDO

 200     NIS = NIS+WDES%NITYP(NT)
         ENDDO
      ELSE IF (  WDES%NCDIJ==4 .AND. LBACK) THEN
!=======================================================================
         FAC=0.5_q
         NIS=1
         LMBASE=0
         typ:  DO NT=1,WDES%NTYP
            LMMAXC=WDES%LMMAX(NT)
            IF (LMMAXC==0) GOTO 300

            ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1

            DO L=1,LMMAXC
               DO LP=1,LMMAXC
                  C00=CRHODE(L,LP,NI,1)
                  CX =CRHODE(L,LP,NI,2)
                  CY =CRHODE(L,LP,NI,3)
                  CZ =CRHODE(L,LP,NI,4)

                  CRHODE(L,LP,NI,1)= (C00+CZ)*FAC
                  CRHODE(L,LP,NI,2)= (CX-CY*(0._q,1._q))*FAC
                  CRHODE(L,LP,NI,3)= (CX+CY*(0._q,1._q))*FAC
                  CRHODE(L,LP,NI,4)= (C00-CZ)*FAC
               ENDDO
            ENDDO
         
            LMBASE= LMMAXC+LMBASE
            ENDDO ion

 300     NIS = NIS+WDES%NITYP(NT)
         ENDDO typ
      ENDIF

    END SUBROUTINE US_FLIP_CMPLX

!************************ SUBROUTINE US_FLIP ***************************
!
! rearranges the storage mode for spin components of array CDIJ:
! given CDIJ_up and CDIJ_down on input the quantities
! (CDIJ_up+CDIJ_down)/2 and (CDIJ_up-CDIJ_down)/2 = total charge
! potential and difference potential are returned
! also the reverse operation is possible if setting LBACK=.TRUE.
!
!***********************************************************************

    SUBROUTINE CDIJ_FLIP(WDES, LMDIM, CDIJ, LOVERL, LBACK)
      USE prec
      USE wave

      IMPLICIT NONE

      TYPE (wavedes)     WDES
      LOGICAL LOVERL, LBACK
      INTEGER LMDIM,LMMAXC,NT,NI,NIS,L,LP,LMBASE
      REAL(q) FAC
      OVERLAP :: CQU,CQD,C01,C10
      OVERLAP :: C11,C00,CX,CY,CZ
      OVERLAP :: CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)


      IF (.NOT.LOVERL) RETURN

      IF (WDES%NCDIJ==2 ) THEN
!=======================================================================
         FAC=0.5_q
         IF (LBACK) FAC=1.0_q
      
         NIS=1
         LMBASE=0
         DO NT=1,WDES%NTYP
            LMMAXC=WDES%LMMAX(NT)
            IF (LMMAXC==0) GOTO 100

            DO NI=NIS,WDES%NITYP(NT)+NIS-1

            DO L=1,LMMAXC
               DO LP=1,LMMAXC
                  CQU=CDIJ(L,LP,NI,1)
                  CQD=CDIJ(L,LP,NI,2)
                  CDIJ(L,LP,NI,1)=FAC*(CQU+CQD)
                  CDIJ(L,LP,NI,2)=FAC*(CQU-CQD)
               ENDDO
            ENDDO
         
            LMBASE= LMMAXC+LMBASE
            ENDDO

 100     NIS = NIS+WDES%NITYP(NT)
         ENDDO
      ELSE IF (  WDES%NCDIJ==4 .AND. .NOT. LBACK) THEN
!=======================================================================
         FAC=0.5_q
         NIS=1
         LMBASE=0
         DO NT=1,WDES%NTYP
            LMMAXC=WDES%LMMAX(NT)
            IF (LMMAXC==0) GOTO 200

            DO NI=NIS,WDES%NITYP(NT)+NIS-1

            DO L=1,LMMAXC
               DO LP=1,LMMAXC
                  C00=CDIJ(L,LP,NI,1)
                  C01=CDIJ(L,LP,NI,2)
                  C10=CDIJ(L,LP,NI,3)
                  C11=CDIJ(L,LP,NI,4)

                  CDIJ(L,LP,NI,1)= (C00+C11)*FAC
                  CDIJ(L,LP,NI,2)= (C01+C10)*FAC
                  CDIJ(L,LP,NI,3)=((C01-C10)*(0._q,1._q))*FAC
                  CDIJ(L,LP,NI,4)= (C00-C11)*FAC
               ENDDO
            ENDDO
         
            LMBASE= LMMAXC+LMBASE
            ENDDO

 200     NIS = NIS+WDES%NITYP(NT)
         ENDDO
      ELSE IF (  WDES%NCDIJ==4 .AND. LBACK) THEN
!=======================================================================
         NIS=1
         LMBASE=0
         typ:  DO NT=1,WDES%NTYP
            LMMAXC=WDES%LMMAX(NT)
            IF (LMMAXC==0) GOTO 300

            ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1

            DO L=1,LMMAXC
               DO LP=1,LMMAXC
                  C00=CDIJ(L,LP,NI,1)
                  CX =CDIJ(L,LP,NI,2)
                  CY =CDIJ(L,LP,NI,3)
                  CZ =CDIJ(L,LP,NI,4)

                  CDIJ(L,LP,NI,1)= (C00+CZ)
                  CDIJ(L,LP,NI,2)= (CX-CY*(0._q,1._q))
                  CDIJ(L,LP,NI,3)= (CX+CY*(0._q,1._q))
                  CDIJ(L,LP,NI,4)= (C00-CZ)
               ENDDO
            ENDDO
         
            LMBASE= LMMAXC+LMBASE
            ENDDO ion

 300     NIS = NIS+WDES%NITYP(NT)
         ENDDO typ
      ENDIF

    END SUBROUTINE CDIJ_FLIP

!************************ SUBROUTINE DEPLE  ****************************
!
! this subroutine calculates  the augmentation 
! charge-density-distribution in real space and adds the 
! pseudo charge density CHDEN
! the result is returned in CHTOT
!
! the second version DEPLE_ADD adds the charge density to CHTOT
!
!***********************************************************************

    SUBROUTINE DEPLE(WDES, GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
        LATT_CUR,P,T_INFO,SYMM, LOVERL, SOFT_TO_C,&
        LMDIM,CRHODE, CHTOT,CHDEN, IRDMAX )
      USE base
      USE pseudo
      USE poscar
      USE mgrid
      USE lattice
      USE wave
      USE us, ONLY: AUGMENTATION_CHARGE
      IMPLICIT NONE

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRID_SOFT,GRIDC,GRIDUS
      TYPE (transit)     C_TO_US
      TYPE (transit)     SOFT_TO_C
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES
      TYPE (symmetry)    SYMM

      INTEGER    LMDIM          
      INTEGER    IRDMAX
      OVERLAP    CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      RGRID      CHTOT(DIMREAL(GRIDC%MPLWV),WDES%NCDIJ)
      COMPLEX(q) CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
      LOGICAL    LOVERL,LADDITIONAL
      REAL(q)    DISPL(3,T_INFO%NIONS)
  ! local
      INTEGER    ISP
!=======================================================================
! if no overlap copy CHDEN to CHTOT and thats it
!=======================================================================
      overl: IF (.NOT.LOVERL) THEN

         DO ISP=1,WDES%NCDIJ
            CALL RC_ADD(CHDEN(1,ISP),1.0_q,CHDEN(1,ISP),0.0_q,CHTOT(1,ISP),GRID_SOFT)
         ENDDO

      ELSE overl
!=======================================================================
! calculate augmentation charge 
!=======================================================================
      CHTOT=0
      DISPL=0
      CALL AUGMENTATION_CHARGE( &
           WDES, GRIDC, GRIDUS, C_TO_US, &
           LATT_CUR, P, T_INFO, SYMM, LOVERL, &
           LMDIM, CRHODE, CHTOT, IRDMAX, DISPL)

      DO ISP=1,WDES%NCDIJ

#ifdef debug
         RHO_AUG =RHO0(GRIDC, CHTOT(1,ISP))
         RHO_SOFT=RHO0(GRID_SOFT , CHDEN)
         DWRITE "augmentation electrons", RHO_AUG
         DWRITE "soft         electrons", RHO_SOFT
#endif
         CALL ADD_GRID(GRIDC, GRID_SOFT, SOFT_TO_C,CHDEN(1,ISP),CHTOT(1,ISP))

         CALL SETUNB_COMPAT(CHTOT(1,ISP),GRIDC)
#ifdef debug
         RHO_AUG =RHO0(GRIDC, CHTOT(1,ISP))
         DWRITE "total        electrons", RHO_AUG
#endif
      ENDDO
! hard set of m_z to zero; this should not be necessary
      IF (WDES%LSPIRAL.AND.WDES%LZEROZ) CHTOT(:,4)=0
      
      ENDIF overl
      RETURN
      END SUBROUTINE


      SUBROUTINE DEPLE_ADD(WDES, GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
        LATT_CUR,P,T_INFO,SYMM, LOVERL, SOFT_TO_C,&
        LMDIM,CRHODE, CHTOT,CHDEN, IRDMAX )
      USE base
      USE pseudo
      USE poscar
      USE mgrid
      USE lattice
      USE wave
      USE us, ONLY: AUGMENTATION_CHARGE
      IMPLICIT NONE

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d)     GRID_SOFT,GRIDC,GRIDUS
      TYPE (transit)     C_TO_US
      TYPE (transit)     SOFT_TO_C
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES
      TYPE (symmetry)    SYMM

      INTEGER    LMDIM          
      INTEGER    IRDMAX
      OVERLAP    CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      RGRID      CHTOT(DIMREAL(GRIDC%MPLWV),WDES%NCDIJ)
      COMPLEX(q) CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
      LOGICAL    LOVERL,LADDITIONAL
      REAL(q)    DISPL(3,T_INFO%NIONS)

  ! local
      INTEGER    ISP
!=======================================================================
! if no overlap copy CHDEN to CHTOT and thats it
!=======================================================================
      overl: IF (.NOT.LOVERL) THEN

         DO ISP=1,WDES%NCDIJ
            CALL RC_ADD(CHTOT(1,ISP),1.0_q,CHDEN(1,ISP),1.0_q,CHTOT(1,ISP),GRID_SOFT)
         ENDDO

      ELSE overl
!=======================================================================
! calculate augmentation charge 
!=======================================================================
      DISPL=0
      CALL AUGMENTATION_CHARGE( &
           WDES, GRIDC, GRIDUS, C_TO_US, &
           LATT_CUR, P, T_INFO, SYMM, LOVERL, &
           LMDIM, CRHODE, CHTOT, IRDMAX, DISPL)

      DO ISP=1,WDES%NCDIJ
         CALL ADD_GRID(GRIDC, GRID_SOFT, SOFT_TO_C,CHDEN(1,ISP),CHTOT(1,ISP))

         CALL SETUNB_COMPAT(CHTOT(1,ISP),GRIDC)
      ENDDO
! hard set of m_z to zero; this should not be necessary
      IF (WDES%LSPIRAL.AND.WDES%LZEROZ) CHTOT(:,4)=0
      
      ENDIF overl
      RETURN
    END SUBROUTINE DEPLE_ADD


!************************ SUBROUTINE AUGMENTATION_CHARGE  **************
!
! this subroutine calculates  the augmentation charge-density-
! distribution in reciprocal space and adds it to the array
! CHTOT
! as input it requires  CRHODE(LM,LMP,ION,ISP)
! and a possible displacement vector DISPL
!
!***********************************************************************

    SUBROUTINE AUGMENTATION_CHARGE( &
           WDES, GRIDC_, GRIDUS, C_TO_US, &
           LATT_CUR, P, T_INFO, SYMM, LOVERL, &
           LMDIM, CRHODE, CHTOT_, IRDMAX, DISPL)
      USE prec
      USE base
      USE charge
      USE pseudo
      USE poscar
      USE mpimy
      USE mgrid
      USE lattice
      USE wave
      USE asa
      USE paw
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d),TARGET ::  GRIDC_,GRIDUS
      TYPE (grid_3d),POINTER :: GRIDC
      TYPE (transit)     C_TO_US
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES
      TYPE (symmetry)    SYMM

      INTEGER   IRDMAX         ! allocation required for augmentation
      INTEGER   LMDIM
      OVERLAP   CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      RGRID,TARGET  :: CHTOT_(DIMREAL(GRIDC_%MPLWV),WDES%NCDIJ)
      LOGICAL   LOVERL,LADDITIONAL
      REAL(q)   DISPL(3,T_INFO%NIONS)
!  work arrays
      TYPE (potcar), POINTER :: PP
      RGRID,POINTER :: CHTOT(:)
      REAL(q)   RHOLM(256)
      REAL(q),ALLOCATABLE ::   DIST(:),DEP(:),SUM(:),YLM(:,:)
      INTEGER,ALLOCATABLE ::   NLI(:)
      INTEGER :: QDEP_FOCK_INDEX
      LOGICAL L_SYM
!  spin spiral stuff
      REAL(q)   QVEC(3),QR
      REAL(q),ALLOCATABLE :: XS(:),YS(:),ZS(:)
      REAL(q)   RHOLMX(256),RHOLMY(256)
      INTEGER, EXTERNAL :: ONE_CENTER_NMAX_FOCKAE
#ifdef MPI
      ! in the MPI version CRHODE holds the contribution to the augmentation
      ! occupation only for ions and bands which are local
      ! CTMP holds all elements (merged)
      ! to achive good load balancing, CTMP is first merged, and then each node
      ! calculates augmentation charges for his columns
      OVERLAP,ALLOCATABLE :: CTMP(:,:,:,:)
      ALLOCATE(CTMP(LMDIM,LMDIM,T_INFO%NIONS,WDES%NCDIJ))
#else
#define CTMP CRHODE
#endif
      ! LADDITIONAL uses an even finer grid for
      ! calculating the augmentation charges
      LADDITIONAL=(GRIDUS%NGX/=GRIDC_%NGX) .OR. &
                  (GRIDUS%NGY/=GRIDC_%NGY) .OR. &
                  (GRIDUS%NGZ/=GRIDC_%NGZ)

      ALLOCATE(XS(IRDMAX),YS(IRDMAX),ZS(IRDMAX))
      IF (WDES%LSPIRAL) THEN
      ! Take QSPIRAL from direct to cartesian coordinates
         QVEC(1)=WDES%QSPIRAL(1)*LATT_CUR%B(1,1)+WDES%QSPIRAL(2)*LATT_CUR%B(1,2)+WDES%QSPIRAL(3)*LATT_CUR%B(1,3)
         QVEC(2)=WDES%QSPIRAL(1)*LATT_CUR%B(2,1)+WDES%QSPIRAL(2)*LATT_CUR%B(2,2)+WDES%QSPIRAL(3)*LATT_CUR%B(2,3)
         QVEC(3)=WDES%QSPIRAL(1)*LATT_CUR%B(3,1)+WDES%QSPIRAL(2)*LATT_CUR%B(3,2)+WDES%QSPIRAL(3)*LATT_CUR%B(3,3)
      ENDIF

! find the maximum L for augmentation charge (usually just 2 l)
      LYDIM=MAXL_AUG(T_INFO%NTYP,P)
      LMYDIM=(LYDIM+1)**2          ! number of lm pairs

      ALLOCATE( &
     &          DIST(IRDMAX),DEP(IRDMAX),SUM(IRDMAX),YLM(IRDMAX,LMYDIM), &
     &          NLI(IRDMAX))

      IF (LADDITIONAL) THEN
         GRIDC => GRIDUS
      ELSE
         GRIDC => GRIDC_
      ENDIF

      ALLOCATE(CHTOT(DIMREAL(GRIDC%MPLWV)))

#ifdef MPI
!=======================================================================
! merge CRHODE from all nodes
! for simplicity I do this with M_sum_d but there are of course better
! ways to do this
!=======================================================================
      CTMP=0
      DO ISP=1,WDES%NCDIJ
         DO NI=1,T_INFO%NIONS
            NIP=NI_LOCAL(NI, WDES%COMM_INB)
            IF (NIP/=0) THEN
               CTMP(:,:,NI,ISP)=CRHODE(:,:,NIP,ISP)
            ENDIF
         ENDDO
      ENDDO
#ifdef realmode
      CALLMPI( M_sum_d(WDES%COMM_INB,CTMP,LMDIM*LMDIM*T_INFO%NIONS*WDES%NCDIJ))
#else
      CALLMPI( M_sum_d(WDES%COMM_INB,CTMP,LMDIM*LMDIM*T_INFO%NIONS*WDES%NCDIJ*2))
#endif
#endif
!-----------------------------------------------------------------------
! do  symmetrization of the CRHODE
! (in MPI version this is the only position where I can do that
!  without additional communication)
!-----------------------------------------------------------------------
! if PAW is selected, do symmetrization in any case
      L_SYM=.FALSE.
      DO NT=1,T_INFO%NTYP
        IF ( ASSOCIATED(P(NT)%QPAW) ) L_SYM=.TRUE.
      ENDDO
! no symmetry used, well switch it off
      IF (SYMM%ISYM<=0) L_SYM=.FALSE.
! switch it off for ISYM=3 as well
      IF (SYMM%ISYM==3 .OR. SYMM%ISYM==4) L_SYM=.FALSE.
! CHDEN is symmetrized and not CHTOT, in that case do symmetrization in any case
      IF (SYMM%ISYM==2) L_SYM=.TRUE.
! now do the symmetrization
      IF (L_SYM) THEN
         IF (WDES%LNONCOLLINEAR) THEN
#ifdef MPI
           CALL AUGSYM_(P,LMDIM,T_INFO%NIONS,T_INFO%NIOND,T_INFO%NTYP,T_INFO%NITYP, &
                CTMP(1,1,1,1), SYMM%ROTMAP(1,1,1), SYMM%MAGROT(1,1), LATT_CUR%A, LATT_CUR%B, 1)
! symmetrize the vectors (DX,DY,DZ)
           IF (.NOT.WDES%LSPIRAL) &
          &   CALL AUGSYM_NONCOL_(P,LMDIM,T_INFO%NIONS,T_INFO%NIOND,T_INFO%NTYP,T_INFO%NITYP, &
                CTMP(1,1,1,2), SYMM%ROTMAP(1,1,1), SYMM%MAGROT(1,1), WDES%SAXIS, LATT_CUR%A, LATT_CUR%B)
! store result back to original storage position
           DO NI=1,T_INFO%NIONS
              NIP=NI_LOCAL(NI, WDES%COMM_INB)
              IF (NIP/=0) THEN
                 CRHODE(:,:,NIP,:)=CTMP(:,:,NI,:)
              ENDIF
           ENDDO
#else
           CALL AUGSYM_(P,LMDIM,T_INFO%NIONS,T_INFO%NIOND,T_INFO%NTYP,T_INFO%NITYP(1), &
              CRHODE(1,1,1,1), SYMM%ROTMAP(1,1,1), SYMM%MAGROT(1,1), LATT_CUR%A, LATT_CUR%B, 1)
! symmetrize the vectors (DX,DY,DZ)
           IF (.NOT.WDES%LSPIRAL) &
          &   CALL AUGSYM_NONCOL_(P,LMDIM,T_INFO%NIONS,T_INFO%NIOND,T_INFO%NTYP,T_INFO%NITYP, &
                CTMP(1,1,1,2), SYMM%ROTMAP(1,1,1), SYMM%MAGROT(1,1), WDES%SAXIS, LATT_CUR%A, LATT_CUR%B)
#endif
         ELSE
         DO ISP=1,WDES%NCDIJ
#ifdef MPI
           CALL AUGSYM_(P,LMDIM,T_INFO%NIONS,T_INFO%NIOND,T_INFO%NTYP,T_INFO%NITYP, &
                CTMP(1,1,1,ISP), SYMM%ROTMAP(1,1,1), SYMM%MAGROT(1,1), LATT_CUR%A, LATT_CUR%B, ISP)
        ! store result back to original storage position
           DO NI=1,T_INFO%NIONS
              NIP=NI_LOCAL(NI, WDES%COMM_INB)
              IF (NIP/=0) THEN
                 CRHODE(:,:,NIP,ISP)=CTMP(:,:,NI,ISP)
              ENDIF
           ENDDO
#else
           CALL AUGSYM_(P,LMDIM,T_INFO%NIONS,T_INFO%NIOND,T_INFO%NTYP,T_INFO%NITYP(1), &
              CRHODE(1,1,1,ISP), SYMM%ROTMAP(1,1,1), SYMM%MAGROT(1,1), LATT_CUR%A, LATT_CUR%B, ISP)
#endif
         ENDDO
        ENDIF
     ENDIF

!=======================================================================
! now the actual work starts
!=======================================================================
     spin:DO ISP=1,WDES%NCDIJ
      CHTOT = 0
!=======================================================================
! loop over all ions
!=======================================================================
      ion: DO NI=1,T_INFO%NIONS
      NT=T_INFO%ITYP(NI)
      PP=>PP_POINTER(P, NI, NT)
!-----------------------------------------------------------------------
! for this ion (this type of ion) no depletion charge
!-----------------------------------------------------------------------
      IF (PP%PSDMAX==0) CYCLE
!-----------------------------------------------------------------------
! calulate the spherical harmonics (DEP is Work-arrays)
!-----------------------------------------------------------------------
      LYMAX=MAXL1(PP)
      IF ( ASSOCIATED(PP%QPAW) ) THEN
! old VASP version: in paw method we truncate the augmentation charge at L=4
!         LYMAX=MIN(4,LYMAX*2)
         LYMAX=LYMAX*2
      ENDIF

      CALL SETYLM_AUG(GRIDC,LATT_CUR,T_INFO%POSION(1,NI),PP%PSDMAX,NPSRNL, &
     &        LMYDIM,LYMAX,YLM(1,1),IRDMAX,INDMAX, &
     &        DISPL(1,NI),DISPL(2,NI), DISPL(3,NI),DIST(1),NLI(1),XS,YS,ZS)

      SUM=0
!=======================================================================
! US-PP
! now loop over pseudopotential indexes L and LP
!=======================================================================
  lpaw: IF ( .NOT. ASSOCIATED(PP%QPAW) ) THEN
      LDEP_INDEX=1

    ! loop over all channels (l,epsilon)
      LM=1
      l_loop:  DO L =1,PP%LMAX
      LMP=LM
      lp_loop: DO LP=L,PP%LMAX
      IF (PP%NDEP(L,LP)==0) GOTO 510

      CALL SETDEP(PP%QDEP(1,1,LDEP_INDEX),PP%PSDMAX,NPSRNL, &
           LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))
      LDEP_INDEX=LDEP_INDEX+ABS(PP%NDEP(L,LP))

    ! quantum numbers l and lp of these two channels
      LL =PP%LPS(L )
      LLP=PP%LPS(LP)

    ! loop over all m mp
      m_loop:  DO M=1,2*LL+1
      MPLOW=1
      IF (L==LP) MPLOW=M
      mp_loop: DO MP=MPLOW,2*LLP+1
         FAKT=1
         IF (LMP+MP/=LM+M) FAKT=2

      !   calculate the indices into the array containing the spherical
      !   harmonics
         INDYLM =LL **2  +M
         INDPYL =LLP**2  +MP

         TFAKT=CTMP(LM+M-1,LMP+MP-1,NI,ISP)*FAKT
         IF (T_INFO%VCA(NT)/=1.0) THEN
            TFAKT=TFAKT*T_INFO%VCA(NT)
         ENDIF
      !   add augmentation charge (augmentation charge is real)
!DIR$ IVDEP
!OCL NOVREC
         DO IND=1,INDMAX
            SUM(IND)=SUM(IND)+DEP(IND)*YLM(IND,INDYLM)*YLM(IND,INDPYL)*TFAKT
         ENDDO

      ENDDO mp_loop
      ENDDO m_loop
  510 LMP=LMP+2*LLP+1
      ENDDO lp_loop
      LM =LM +2*LL +1
      ENDDO l_loop
   ELSE lpaw
!=======================================================================
! PAW
! transform CRHODE from the basis L,L',M,M'  to
! the basis L,LP,Lmain,Mmain using Clebsch-Gordan coefficients
! then add the compensation charge to grid
!=======================================================================
      IF (WDES%LSPIRAL .AND. ISP==2) THEN
! calculate cell periodic part of the augmentation magnetization density
! by rotating rho_x and rho_y against the spiral
! rho_x -> cos(qr)rho_x + sin(qr)rho_y
      CALL CALC_RHOLM( LYMAX, CTMP(:,:,NI,2) , RHOLMX, PP)
      CALL CALC_RHOLM( LYMAX, CTMP(:,:,NI,3) , RHOLMY, PP)

      DO L =0,LYMAX
         ! WRITE(0,'("RHOLM",I2,10F10.6)') L,(RHOLM(L**2+M),M=1,(L*2)+1)
         CALL SETDEP(PP%QDEP(1,1,L),PP%PSDMAX,NPSRNL, &
              LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))
         DO M=1,(L*2)+1
            INDYLM =L **2  +M
            TFAKT=RHOLM(INDYLM)
            IF (T_INFO%VCA(NT)/=1.0) THEN
               TFAKT=TFAKT*T_INFO%VCA(NT)
            ENDIF
!DIR$ IVDEP
!OCL NOVREC
            DO IND=1,INDMAX
               QR=TPI*(QVEC(1)*XS(IND)+QVEC(2)*YS(IND)+QVEC(3)*ZS(IND))
               SUM(IND)=SUM(IND)+DEP(IND)*YLM(IND,INDYLM)*(RHOLMX(INDYLM)*COS(QR)+RHOLMY(INDYLM)*SIN(QR))
            ENDDO

         ENDDO
      ENDDO
      ENDIF
      
      IF (WDES%LSPIRAL .AND. ISP==3) THEN
! calculate cell periodic part of the augmentation magnetization density
! by rotating rho_x and rho_y against the spiral
! rho_y -> cos(qr)rho_y - sin(qr)rho_x
      CALL CALC_RHOLM( LYMAX, CTMP(:,:,NI,2) , RHOLMX, PP)
      CALL CALC_RHOLM( LYMAX, CTMP(:,:,NI,3) , RHOLMY, PP)

      DO L =0,LYMAX
         ! WRITE(0,'("RHOLM",I2,10F10.6)') L,(RHOLM(L**2+M),M=1,(L*2)+1)
         CALL SETDEP(PP%QDEP(1,1,L),PP%PSDMAX,NPSRNL, &
              LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))
         DO M=1,(L*2)+1
            INDYLM =L **2  +M
            TFAKT=RHOLM(INDYLM)
            IF (T_INFO%VCA(NT)/=1.0) THEN
               TFAKT=TFAKT*T_INFO%VCA(NT)
            ENDIF
!DIR$ IVDEP
!OCL NOVREC
            DO IND=1,INDMAX
               QR=TPI*(QVEC(1)*XS(IND)+QVEC(2)*YS(IND)+QVEC(3)*ZS(IND))
               SUM(IND)=SUM(IND)+DEP(IND)*YLM(IND,INDYLM)*(-RHOLMX(INDYLM)*SIN(QR)+RHOLMY(INDYLM)*COS(QR))
            ENDDO

         ENDDO
      ENDDO
      ENDIF
      
      IF (.NOT.WDES%LSPIRAL .OR. ISP==1 .OR. ISP==4) THEN
! no phase factor for total charge and m_z or if LSPIRAL=.FALSE.
      CALL CALC_RHOLM( LYMAX, CTMP(:,:,NI,ISP) , RHOLM, PP)

      DO L =0,LYMAX
         ! WRITE(0,'("RHOLM",I2,10F10.6)') L,(RHOLM(L**2+M),M=1,(L*2)+1)
         CALL SETDEP(PP%QDEP(1,1,L),PP%PSDMAX,NPSRNL, &
              LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))
         DO M=1,(L*2)+1
            INDYLM =L **2  +M
            TFAKT=RHOLM(INDYLM)
            IF (T_INFO%VCA(NT)/=1.0) THEN
               TFAKT=TFAKT*T_INFO%VCA(NT)
            ENDIF
!DIR$ IVDEP
!OCL NOVREC
            DO IND=1,INDMAX
               SUM(IND)=SUM(IND)+DEP(IND)*YLM(IND,INDYLM)*TFAKT
            ENDDO

         ENDDO
      ENDDO
!=======================================================================
! PAW contributions from accurate augmentation 
! currently spin spirals are no supported
!=======================================================================
      onecAE: IF (ONE_CENTER_NMAX_FOCKAE()>0) THEN
      NMAX_FOCKAE=SIZE(PP%QPAW_FOCK,4)
      LMAX_FOCKAE=SIZE(PP%QPAW_FOCK,3)-1
      ! proper indexing of QDEP_FOCK is nasty see fast_aug.F
      QDEP_FOCK_INDEX=SIZE(PP%QDEP_FOCK,3)-NMAX_FOCKAE*(LMAX_FOCKAE+1)-1

      DO NAE=1,NMAX_FOCKAE
      CALL CALC_RHOLM_FOCK( LYMAX, CTMP(:,:,NI,ISP) , RHOLM, PP, NAE)

      DO L  =0,LMAX_FOCKAE
         QDEP_FOCK_INDEX=QDEP_FOCK_INDEX+1

         CALL SETDEP(PP%QDEP_FOCK(1,1,QDEP_FOCK_INDEX),PP%PSDMAX,NPSRNL, &
              LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))
         DO M=1,(L*2)+1
            INDYLM =L **2  +M
            TFAKT=RHOLM(INDYLM)
            IF (T_INFO%VCA(NT)/=1.0) THEN
               TFAKT=TFAKT*T_INFO%VCA(NT)
            ENDIF
!DIR$ IVDEP
!OCL NOVREC
            DO IND=1,INDMAX
               SUM(IND)=SUM(IND)+DEP(IND)*YLM(IND,INDYLM)*TFAKT
            ENDDO
         ENDDO
      ENDDO
      ENDDO

    ENDIF onecAE
    ENDIF
    ENDIF lpaw
!=======================================================================
! add the calculated augmentation charge to the total charge
!=======================================================================
      SUMN=0
      DO IND=1,INDMAX
        CHTOT(NLI(IND))=CHTOT(NLI(IND))+SUM(IND)
        SUMN=SUMN+SUM(IND)
      ENDDO
!-----------------------------------------------------------------------
      ENDDO ion
!-----------------------------------------------------------------------
      CALL FFT_RC_SCALE(CHTOT(1),CHTOT(1),GRIDC)

      IF (LADDITIONAL) THEN
         CALL CP_ADD_GRID(GRIDUS,GRIDC_,C_TO_US,CHTOT(1),CHTOT_(1,ISP))
      ELSE
         CALL RC_ADD(CHTOT_(1,ISP), 1.0_q, CHTOT(1), 1.0_q, CHTOT_(1,ISP), GRIDC_)
      ENDIF
      ENDDO spin

      DEALLOCATE(CHTOT)
      DEALLOCATE(DIST,DEP,SUM,YLM,NLI,XS,YS,ZS)


#ifdef MPI
      DEALLOCATE(CTMP)
#else
#undef CTMP
#endif


      RETURN
    END SUBROUTINE AUGMENTATION_CHARGE


!************************ SUBROUTINE CURRENT_AUGMENTATION **************
!
! this subroutine calculates  the augmentation contribution to 
! the current density JTOT
! as input it requires  CRHODE(LM,LMP,ION,ISP)
! and a possible displacement vector DISPL
!
!***********************************************************************


    SUBROUTINE CURRENT_AUGMENTATION( &
           WDES, GRIDC_, GRIDUS, C_TO_US, &
           LATT_CUR, P, T_INFO, &
           LMDIM, CRHODE, JTOT_, IRDMAX, DISPL, LMONOPOLE)
      USE prec
      USE base
      USE charge
      USE pseudo
      USE poscar
      USE mpimy
      USE mgrid
      USE lattice
      USE wave
      USE asa
      USE paw
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d),TARGET ::  GRIDC_,GRIDUS
      TYPE (grid_3d),POINTER :: GRIDC
      TYPE (transit)     C_TO_US
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES

      INTEGER   IRDMAX         ! allocation required for augmentation
      INTEGER   LMDIM
      OVERLAP   CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      RGRID,TARGET  :: JTOT_(DIMREAL(GRIDC_%MPLWV),3)
      REAL(q), POINTER :: A_ONE_CENTER(:,:)
      REAL(q)   DISPL(3,T_INFO%NIONS)
!  work arrays
      INTEGER       :: LYMAX, LYMAXP1
      TYPE (potcar), POINTER :: PP
      RGRID,POINTER :: JTOT(:)
      OVERLAP       :: JLM(256)
      REAL(q)       :: RHOLM(256)
      REAL(q),ALLOCATABLE ::   DIST(:),DEP(:),YLM(:,:),XS(:),YS(:),ZS(:)
      OVERLAP,ALLOCATABLE ::   SUM(:)
      INTEGER,ALLOCATABLE ::   NLI(:)
      OVERLAP,ALLOCATABLE ::   CTMP(:,:,:)
      LOGICAL :: LMONOPOLE
!=======================================================================
! allocation of all required quantities
!=======================================================================

! find the maximum L for augmentation charge (usually just 2 maximum l)
      LYDIM=MAXL_AUG(T_INFO%NTYP,P)+1
      LMYDIM=(LYDIM+1)**2          ! number of lm pairs

      ALLOCATE( DIST(IRDMAX),DEP(IRDMAX),SUM(IRDMAX),YLM(IRDMAX,LMYDIM), &
                NLI(IRDMAX),XS(IRDMAX),YS(IRDMAX),ZS(IRDMAX))

      GRIDC => GRIDC_

      ALLOCATE(JTOT(DIMREAL(GRIDC%MPLWV)),CTMP(LMDIM,LMDIM,T_INFO%NIONS))
!-----------------------------------------------------------------------
! merge augmentation occupancies CRHODE from all nodes
! for simplicity I do this using M_sum_d but there are of course better
! ways to do this
!-----------------------------------------------------------------------

! hopefully CRHODE is in the representation such that the *total*
! charge is stored in channel one
      CTMP=0
      DO NI=1,T_INFO%NIONS
         NIP=NI_LOCAL(NI, WDES%COMM_INB)
         IF (NIP/=0) THEN
            CTMP(:,:,NI)=CRHODE(:,:,NIP,1)
         ENDIF
      ENDDO
#ifdef realmode
      CALLMPI( M_sum_d(WDES%COMM_INB,CTMP,LMDIM*LMDIM*T_INFO%NIONS))
#else
      CALLMPI( M_sum_z(WDES%COMM_INB,CTMP,LMDIM*LMDIM*T_INFO%NIONS))
#endif
!-----------------------------------------------------------------------
! loop over cartesian index and ion
!-----------------------------------------------------------------------
   cart: DO ICART=1,3
      JTOT = 0
      ion: DO NI=1,T_INFO%NIONS
      NT=T_INFO%ITYP(NI)
      PP=>PP_POINTER(P, NI, NT)
!-----------------------------------------------------------------------
! for this ion (this type of ion) no depletion charge
!-----------------------------------------------------------------------
      IF (PP%PSDMAX==0 .OR.  .NOT. ASSOCIATED(PP%JPAW) ) CYCLE
!-----------------------------------------------------------------------
! calulate the spherical harmonics (DEP is Work-arrays)
!-----------------------------------------------------------------------
      ! nabla operator increases the maximum L one-center number by one
      LYMAX  =MAXL1(PP)*2
      LYMAXP1=LYMAX+1

      CALL SETYLM_AUG(GRIDC,LATT_CUR,T_INFO%POSION(1,NI),PP%PSDMAX,NPSRNL, &
     &        LMYDIM,LYMAXP1,YLM(1,1),IRDMAX,INDMAX, &
     &        DISPL(1,NI),DISPL(2,NI), DISPL(3,NI),DIST(1),NLI(1),XS,YS,ZS)

      JLM=0
      CALL CALC_RHOLM_JVEC( CTMP(:,:,NI) , JLM, PP, ICART)

      IF (.NOT. LMONOPOLE) THEN
        JLM(1)=0
      ENDIF

      SUM=0

      DO L =0,LYMAXP1
!         IF (L<=1.AND.WDES%COMM%NODE_ME==WDES%COMM%IONODE) WRITE(*,'("JAUG",3I4,10F10.6)') ICART,NI,L,(AIMAG(JLM(L**2+M))*1E6,M=1,(L*2)+1)
         CALL SETDEP(PP%QDEP(1,1,L),PP%PSDMAX,NPSRNL, &
              LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))

         DO M=1,(L*2)+1
            INDYLM =L **2  +M
!DIR$ IVDEP
!OCL NOVREC
            DO IND=1,INDMAX
               SUM(IND)=SUM(IND)+DEP(IND)*YLM(IND,INDYLM)*JLM(INDYLM)
            ENDDO
         ENDDO
      ENDDO
!-----------------------------------------------------------------------
! add the calculated augmentation currents to the current density
!-----------------------------------------------------------------------
      DO IND=1,INDMAX
        JTOT(NLI(IND))=JTOT(NLI(IND))-AIMAG(SUM(IND)*(1.0_q,0.0_q))
      ENDDO
!-----------------------------------------------------------------------
      ENDDO ion
!-----------------------------------------------------------------------
      CALL FFT_RC_SCALE(JTOT(1),JTOT(1),GRIDC)

      CALL RC_ADD(JTOT_(1,ICART), 1.0_q, JTOT(1), 1.0_q, JTOT_(1,ICART), GRIDC_)
      ENDDO cart

      DEALLOCATE(DIST,DEP,SUM,YLM,NLI,XS,YS,ZS)
      DEALLOCATE(JTOT,CTMP)

      RETURN
   END SUBROUTINE CURRENT_AUGMENTATION

   SUBROUTINE CURRENT_AUGMENTATION_CMPLX( &
           WDES, GRIDC_, GRIDUS, C_TO_US, &
           LATT_CUR, P, T_INFO, &
           LMDIM, CRHODE, JTOT_, IRDMAX, DISPL, LMONOPOLE)
      USE prec
      USE base
      USE charge
      USE pseudo
      USE poscar
      USE mpimy
      USE mgrid
      USE lattice
      USE wave
      USE asa
      USE paw
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (grid_3d),TARGET ::  GRIDC_,GRIDUS
      TYPE (grid_3d),POINTER :: GRIDC
      TYPE (transit)     C_TO_US
      TYPE (latt)        LATT_CUR
      TYPE (wavedes)     WDES

      INTEGER   IRDMAX         ! allocation required for augmentation
      INTEGER   LMDIM
      COMPLEX(q) CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)                   !GIL
      RGRID,TARGET  :: JTOT_(DIMREAL(GRIDC_%MPLWV),3)
      REAL(q), POINTER :: A_ONE_CENTER(:,:)
      REAL(q)   DISPL(3,T_INFO%NIONS)
!  work arrays
      INTEGER       :: LYMAX, LYMAXP1
      TYPE (potcar), POINTER :: PP
      RGRID,POINTER :: JTOT(:)
      OVERLAP       :: JLM(256)
      REAL(q)       :: RHOLM(256)
      REAL(q),ALLOCATABLE ::   DIST(:),DEP(:),YLM(:,:),XS(:),YS(:),ZS(:)
      OVERLAP,ALLOCATABLE ::   SUM(:)
      INTEGER,ALLOCATABLE ::   NLI(:)
      OVERLAP,ALLOCATABLE ::   CTMP(:,:,:)
      LOGICAL :: LMONOPOLE
!=======================================================================
! allocation of all required quantities
!=======================================================================

! find the maximum L for augmentation charge (usually just 2 maximum l)
      LYDIM=MAXL_AUG(T_INFO%NTYP,P)+1
      LMYDIM=(LYDIM+1)**2          ! number of lm pairs

      ALLOCATE( DIST(IRDMAX),DEP(IRDMAX),SUM(IRDMAX),YLM(IRDMAX,LMYDIM), &
                NLI(IRDMAX),XS(IRDMAX),YS(IRDMAX),ZS(IRDMAX))

      GRIDC => GRIDC_

      ALLOCATE(JTOT(DIMREAL(GRIDC%MPLWV)),CTMP(LMDIM,LMDIM,T_INFO%NIONS))
!-----------------------------------------------------------------------
! merge augmentation occupancies CRHODE from all nodes
! for simplicity I do this using M_sum_d but there are of course better
! ways to do this
!-----------------------------------------------------------------------

! hopefully CRHODE is in the representation such that the *total*
! charge is stored in channel one
      CTMP=0
      DO NI=1,T_INFO%NIONS
         NIP=NI_LOCAL(NI, WDES%COMM_INB)
         IF (NIP/=0) THEN
            CTMP(:,:,NI)=-AIMAG(CRHODE(:,:,NIP,1))                            !GIL
         ENDIF
      ENDDO
#ifdef realmode
      CALLMPI( M_sum_d(WDES%COMM_INB,CTMP,LMDIM*LMDIM*T_INFO%NIONS))
#else
      CALLMPI( M_sum_z(WDES%COMM_INB,CTMP,LMDIM*LMDIM*T_INFO%NIONS))
#endif
!-----------------------------------------------------------------------
! loop over cartesian index and ion
!-----------------------------------------------------------------------
   cart: DO ICART=1,3
      JTOT = 0
      ion: DO NI=1,T_INFO%NIONS
      NT=T_INFO%ITYP(NI)
      PP=>PP_POINTER(P, NI, NT)
!-----------------------------------------------------------------------
! for this ion (this type of ion) no depletion charge
!-----------------------------------------------------------------------
      IF (PP%PSDMAX==0 .OR.  .NOT. ASSOCIATED(PP%JPAW) ) CYCLE
!-----------------------------------------------------------------------
! calulate the spherical harmonics (DEP is Work-arrays)
!-----------------------------------------------------------------------
      ! nabla operator increases the maximum L one-center number by one
      LYMAX  =MAXL1(PP)*2
      LYMAXP1=LYMAX+1

      CALL SETYLM_AUG(GRIDC,LATT_CUR,T_INFO%POSION(1,NI),PP%PSDMAX,NPSRNL, &
     &        LMYDIM,LYMAXP1,YLM(1,1),IRDMAX,INDMAX, &
     &        DISPL(1,NI),DISPL(2,NI), DISPL(3,NI),DIST(1),NLI(1),XS,YS,ZS)

      JLM=0
      CALL CALC_RHOLM_JVEC( CTMP(:,:,NI) , JLM, PP, ICART)

      IF (.NOT. LMONOPOLE) THEN
        JLM(1)=0
      ENDIF

      SUM=0

      DO L =0,LYMAXP1
!         IF (L<=1.AND.WDES%COMM%NODE_ME==WDES%COMM%IONODE) WRITE(*,'("JAUG",3I4,10F10.6)') ICART,NI,L,(AIMAG(JLM(L**2+M))*1E6,M=1,(L*2)+1)
         CALL SETDEP(PP%QDEP(1,1,L),PP%PSDMAX,NPSRNL, &
              LATT_CUR%OMEGA,INDMAX,DIST(1),DEP(1))

         DO M=1,(L*2)+1
            INDYLM =L **2  +M
!DIR$ IVDEP
!OCL NOVREC
            DO IND=1,INDMAX
               SUM(IND)=SUM(IND)+DEP(IND)*YLM(IND,INDYLM)*JLM(INDYLM)
            ENDDO
         ENDDO
      ENDDO
!-----------------------------------------------------------------------
! add the calculated augmentation currents to the current density
!-----------------------------------------------------------------------
      DO IND=1,INDMAX
        JTOT(NLI(IND))=JTOT(NLI(IND))-(SUM(IND)*(1.0_q,0.0_q))                      !GIL
      ENDDO
!-----------------------------------------------------------------------
      ENDDO ion
!-----------------------------------------------------------------------
      CALL FFT_RC_SCALE(JTOT(1),JTOT(1),GRIDC)                                      !OEI

      CALL RC_ADD(JTOT_(1,ICART), 1.0_q, JTOT(1), 1.0_q, JTOT_(1,ICART), GRIDC_)
      ENDDO cart

      DEALLOCATE(DIST,DEP,SUM,YLM,NLI,XS,YS,ZS)
      DEALLOCATE(JTOT,CTMP)

      RETURN
   END SUBROUTINE CURRENT_AUGMENTATION_CMPLX

