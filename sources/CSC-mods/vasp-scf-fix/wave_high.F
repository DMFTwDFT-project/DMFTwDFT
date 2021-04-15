#include "symbol.inc"
#define boundary_check
#define MPI_barrier_after_bcast
!***********************************************************************
!
! this module implements high level routines to operate on wavefunctions
! care is taken to avoid indexing pointer arrays since
! this incures performance penalties on all Intel compilers 
! (but particularly on the efc compiler)
! (constructs such as W1%CW(I,..) are avoided)
! instead low level F77 routines are used, which are implemented
! at the end of the routine
!
!***********************************************************************

MODULE wave_high
  USE prec
  USE wave
  USE wave_mpi
!***********************************************************************
!
! interfaces for functions that accept only F77 style arrays
! these F77 routines would cause a problem if the first element
! of a pointer array is passed to the F77 routines
!
!***********************************************************************

  INTERFACE
     SUBROUTINE ARRAY_TO_W1( W1, C, CPROJ)
       USE wave
       TYPE (wavefun1)    W1
       COMPLEX(q):: C
       GDEF, OPTIONAL :: CPROJ
     END SUBROUTINE ARRAY_TO_W1
  END INTERFACE
  
  INTERFACE
     SUBROUTINE W1_TO_ARRAY( W1, C, CPROJ)
       USE wave
       TYPE (wavefun1)    W1
       COMPLEX(q):: C
       GDEF, OPTIONAL :: CPROJ
     END SUBROUTINE W1_TO_ARRAY
  END INTERFACE

  INTERFACE
     SUBROUTINE ECCP_NL(LMDIM,LMMAXC,CDIJ,CPROJ1,CPROJ2,CNL)
       USE prec
       GDEF      CNL
       INTEGER LMDIM, LMMAXC
       OVERLAP CDIJ
       GDEF CPROJ1,CPROJ2
     END SUBROUTINE ECCP_NL
  END INTERFACE

  INTERFACE
     SUBROUTINE OVERL(WDES1, LOVERL, LMDIM, CQIJ, CPROF, CRESUL)
       USE wave
       TYPE (wavedes1) WDES1
       LOGICAL LOVERL
       INTEGER LMDIM
       OVERLAP CQIJ,CDIJ
       GDEF CRESUL,CPROF
     END SUBROUTINE OVERL
  END INTERFACE

  INTERFACE
     SUBROUTINE OVERL1(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF,CRESUL)
       USE wave
       TYPE (wavedes1) WDES1
       INTEGER LMDIM
       OVERLAP CQIJ,CDIJ
       REAL(q) :: EVALUE
       GDEF CRESUL,CPROF
     END SUBROUTINE OVERL1
  END INTERFACE

  INTERFACE
     SUBROUTINE OVERL1_C(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF,CRESUL)
       USE wave
       TYPE (wavedes1) WDES1
       INTEGER LMDIM
       OVERLAP CQIJ,CDIJ
       COMPLEX(q) EVALUE
       GDEF CRESUL,CPROF
     END SUBROUTINE OVERL1_C
  END INTERFACE

  INTERFACE
     SUBROUTINE OVERL1_CCDIJ(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF,CRESUL)
       USE wave
       TYPE (wavedes1) WDES1
       INTEGER LMDIM
       COMPLEX(q) CQIJ,CDIJ
       REAL(q) :: EVALUE
       GDEF CRESUL,CPROF
     END SUBROUTINE OVERL1_CCDIJ
  END INTERFACE


  INTERFACE ELEMENT
     MODULE PROCEDURE W1_FROM_W
     MODULE PROCEDURE W1_FROM_WA
     MODULE PROCEDURE W1_FROM_WA2
  END INTERFACE

  INTERFACE ELEMENTS
     MODULE PROCEDURE WA_FROM_W
     MODULE PROCEDURE WA_FROM_WA
     MODULE PROCEDURE WA_FROM_WA2
  END INTERFACE

  INTERFACE REDISTRIBUTE_PROJ
     MODULE PROCEDURE W_REDIS_PROJ
     MODULE PROCEDURE WA_REDIS_PROJ
     MODULE PROCEDURE W1_REDIS_PROJ
  END INTERFACE

  INTERFACE REDISTRIBUTE_PW
     MODULE PROCEDURE W_REDIS_PW
     MODULE PROCEDURE WA_REDIS_PW
     MODULE PROCEDURE W1_REDIS_PW
  END INTERFACE

  INTERFACE ASSIGNMENT
     MODULE PROCEDURE W1_COPY_REVERSE_ARG
  END INTERFACE
  
  CONTAINS

!***********************************************************************
!
! FFT of a wavefunction to real space
!      
!***********************************************************************

    SUBROUTINE FFTWAV_W1( W1)
      IMPLICIT NONE
      INTEGER ISPINOR
      TYPE (wavefun1)    W1

      IF (W1%WDES1%NK==0) THEN
         WRITE(0,*) 'internal error in FFTWAV_W1: NK is set to zero: WDES not set up properly'
         STOP
      ENDIF
      DO ISPINOR=0,W1%WDES1%NRSPINORS-1
         CALL FFTWAV(W1%WDES1%NGVECTOR, W1%WDES1%NINDPW(1), & 
              W1%CR(1+ISPINOR*W1%WDES1%GRID%MPLWV), & 
              W1%CW(1+ISPINOR*W1%WDES1%NGVECTOR),W1%WDES1%GRID)
      ENDDO
    END SUBROUTINE FFTWAV_W1


!***********************************************************************
!
! This routine accomplishes the same thing as the assignment:
!
!   W1_LHS = W1_RHS
!
! would, but allows for the use of "fast"-memory (in the sense that
! the fm_* arrays are excluded from the assignment).
!
!***********************************************************************

  SUBROUTINE W1_ASSIGN(W1_LHS,W1_RHS)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1_RHS
    TYPE (wavefun1) :: W1_LHS

    W1_LHS%WDES1  =>W1_RHS%WDES1
    W1_LHS%FERWE  = W1_RHS%FERWE
    W1_LHS%AUX    = W1_RHS%AUX
    W1_LHS%CELEN  = W1_RHS%CELEN
    W1_LHS%NB     = W1_RHS%NB
    W1_LHS%ISP    = W1_RHS%ISP  
    W1_LHS%LDO    = W1_RHS%LDO

    W1_LHS%CPTWFP =>W1_RHS%CPTWFP
    W1_LHS%CPROJ  =>W1_RHS%CPROJ    
    W1_LHS%CR     =>W1_RHS%CR

  END SUBROUTINE W1_ASSIGN


!***********************************************************************
!
! copy a W1 structure
!  W2 = W1  (W1 -> W2)
! the argument arrangement is similar to  DCOPY, ZCOPY in BLAS level 1
! for syntactic shugar the reverse operation is also present
!
!
!***********************************************************************

  SUBROUTINE W1_COPY( W1, W2)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1    
    TYPE (wavefun1) :: W2

    PROFILING_START('w1_copy')

#if defined (use_fastbcopy)
    CALL FASTBCOPY8(W1%WDES1%NRPLWV*2 ,W1%CW(1), W2%CW(1))
#else
    CALL ZCOPY( W1%WDES1%NRPLWV, W1%CW(1), 1, W2%CW(1), 1)
#endif

    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( W1%WDES1%NPROD, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
       CALL ZCOPY( W1%WDES1%NPROD, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF

    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          WRITE(0,*) 'internal error in W1_COPY: real space allocation is not correct'
          STOP
       ENDIF
       IF (SIZE(W1%CR) /=SIZE(W2%CR)) THEN
          WRITE(0,*) 'internal error in W1_COPY: real space allocation is different',SIZE(W1%CR), SIZE(W2%CR)
          STOP
       ENDIF
       CALL ZCOPY( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS, W1%CR(1), 1, W2%CR(1), 1)
    END IF

    PROFILING_STOP('w1_copy')

  END SUBROUTINE W1_COPY

  SUBROUTINE W1_COPY_NOCR( W1, W2)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1    
    TYPE (wavefun1) :: W2

#if defined (use_fastbcopy)
    CALL FASTBCOPY8(W1%WDES1%NRPLWV*2 ,W1%CW(1), W2%CW(1))
#else
    CALL ZCOPY( W1%WDES1%NRPLWV, W1%CW(1), 1, W2%CW(1), 1)
#endif

    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( W1%WDES1%NPROD, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
       CALL ZCOPY( W1%WDES1%NPROD, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF
  END SUBROUTINE W1_COPY_NOCR

  SUBROUTINE W1_COPY_CPROJ( W1, W2)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1
    TYPE (wavefun1) :: W2
    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( W1%WDES1%NPROD, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
       CALL ZCOPY( W1%WDES1%NPROD, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF
  END SUBROUTINE W1_COPY_CPROJ

  SUBROUTINE W1_COPY_REVERSE_ARG( W2, W1 )
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1
    TYPE (wavefun1) :: W2

#if defined (use_fastbcopy)
    CALL FASTBCOPY8(W1%WDES1%NRPLWV*2 ,W1%CW(1), W2%CW(1))
#else
    CALL ZCOPY( W1%WDES1%NRPLWV, W1%CW(1), 1, W2%CW(1), 1)
#endif
    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( W1%WDES1%NPROD, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
       CALL ZCOPY( W1%WDES1%NPROD, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF

    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          WRITE(0,*) 'internal error in W1_COPY: real space allocation is not correct'
          STOP
       ENDIF
       IF (SIZE(W1%CR) /=SIZE(W2%CR)) THEN
          WRITE(0,*) 'internal error in W1_COPY: real space allocation is different',SIZE(W1%CR), SIZE(W2%CR)
          STOP
       ENDIF

       CALL ZCOPY( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS, W1%CR(1), 1, W2%CR(1), 1)
    END IF

    
  END SUBROUTINE W1_COPY_REVERSE_ARG


!***********************************************************************
!
! copy a WA structure
!  W2 = W1  (W1 -> W2)
! the argument arrangement is similar to  DCOPY, ZCOPY in BLAS level 1
! the redistributed wavefunctions pointers in the destination
! are also properly set
!
!***********************************************************************

  SUBROUTINE WA_COPY( W1, W2 )
    IMPLICIT NONE
    TYPE (wavefuna), INTENT(IN) :: W1    
    TYPE (wavefuna) :: W2

    IF (SIZE(W1%CW) /= SIZE(W2%CW)) THEN
       WRITE(0,*)'internal error in WA_COPY: size mismatch in CW',SIZE(W1%CW), SIZE(W2%CW)
       STOP
    ENDIF

    IF (SIZE(W1%CPROJ) /= SIZE(W2%CPROJ)) THEN
       WRITE(0,*)'internal error in WA_COPY: size mismatch in CPROJ',SIZE(W1%CPROJ), SIZE(W2%CPROJ)
       STOP
    ENDIF

    CALL ZCOPY( SIZE(W1%CW), W1%CW(1,1), 1, W2%CW(1,1), 1)

    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( SIZE(W1%CPROJ), W1%CPROJ(1,1), 1,  W2%CPROJ(1,1), 1)
    ELSE
       CALL ZCOPY( SIZE(W1%CPROJ), W1%CPROJ(1,1), 1,  W2%CPROJ(1,1), 1)
    ENDIF
    ! remember WDES1
    W2%WDES1 =>W1%WDES1
    ! one dimensional indexing assumed
    W2%FIRST_DIM=0
    ! remember spin index
    W2%ISP =W1%ISP
    ! set redistributed wavefunction indices
    IF (W2%WDES1%DO_REDIS) THEN
       CALL SET_WPOINTER(W2%CW_RED,    W2%WDES1%NRPLWV_RED, W2%WDES1%NB_TOT, W2%CW(1,1))
       CALL SET_GPOINTER(W2%CPROJ_RED, W2%WDES1%NPROD_RED,  W2%WDES1%NB_TOT, W2%CPROJ(1,1))
    ELSE
       W2%CW_RED=>W2%CPTWFP
       W2%CPROJ_RED=>W2%CPROJ
    ENDIF


  END SUBROUTINE WA_COPY


  SUBROUTINE WA_COPY_CPROJ( W1, W2 )
    IMPLICIT NONE
    TYPE (wavefuna), INTENT(IN) :: W1    
    TYPE (wavefuna) :: W2

    IF (SIZE(W1%CPROJ) /= SIZE(W2%CPROJ)) THEN
       WRITE(0,*)'internal error in WA_COPY: size mismatch in CPROJ',SIZE(W1%CPROJ), SIZE(W2%CPROJ)
       STOP
    ENDIF

    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( SIZE(W1%CPROJ), W1%CPROJ(1,1), 1,  W2%CPROJ(1,1), 1)
    ELSE
       CALL ZCOPY( SIZE(W1%CPROJ), W1%CPROJ(1,1), 1,  W2%CPROJ(1,1), 1)
      ENDIF
    ! remember WDES1
    W2%WDES1 =>W1%WDES1
    ! one dimensional indexing assumed
    W2%FIRST_DIM=0
    ! remember spin index
    W2%ISP =W1%ISP
    ! set redistributed wavefunction indices
    IF (W2%WDES1%DO_REDIS) THEN
       CALL SET_GPOINTER(W2%CPROJ_RED, W2%WDES1%NPROD_RED,  W2%WDES1%NB_TOT, W2%CPROJ(1,1))
    ELSE
       W2%CPROJ_RED=>W2%CPROJ
    ENDIF


  END SUBROUTINE WA_COPY_CPROJ


!***********************************************************************
!
! redistribute the wavefunction character for a W1 or WA array
!
!***********************************************************************

  SUBROUTINE W1_REDIS_PROJ( W1)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1

    IF (W1%WDES1%DO_REDIS) CALL REDIS_PROJ(W1%WDES1, 1, W1%CPROJ(1))
  END SUBROUTINE W1_REDIS_PROJ

  SUBROUTINE WA_REDIS_PROJ( WA)
    IMPLICIT NONE
    TYPE (wavefuna), INTENT(IN) :: WA

    IF (WA%WDES1%DO_REDIS) CALL REDIS_PROJ(WA%WDES1, SIZE(WA%CPROJ,2), WA%CPROJ(1,1))
  END SUBROUTINE WA_REDIS_PROJ

  SUBROUTINE W_REDIS_PROJ( W)
    IMPLICIT NONE
    TYPE (wavespin), INTENT(IN) :: W
    TYPE (wavedes1)    WDES1          ! descriptor for one k-point
    INTEGER :: K1, ISP

    IF (W%WDES%DO_REDIS) THEN
       DO K1=1,W%WDES%NKPTS
#ifdef MPI
          IF (MOD(K1-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
          CALL SETWDES(W%WDES,WDES1,K1)
          DO ISP=1,W%WDES%ISPIN
             CALL WA_REDIS_PROJ( ELEMENTS( W, WDES1, ISP))
          ENDDO
       ENDDO
    ENDIF
  END SUBROUTINE W_REDIS_PROJ



!***********************************************************************
!
! redistribute the plane wave coefficients for a W1 or WA array
!
!***********************************************************************

  SUBROUTINE W1_REDIS_PW( W1)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1

    IF (W1%WDES1%DO_REDIS) CALL REDIS_PW(W1%WDES1, 1, W1%CW(1))
  END SUBROUTINE W1_REDIS_PW

  SUBROUTINE WA_REDIS_PW( WA)
    IMPLICIT NONE
    TYPE (wavefuna), INTENT(IN) :: WA

    IF (WA%WDES1%DO_REDIS) CALL REDIS_PW(WA%WDES1, SIZE(WA%CW,2), WA%CW(1,1))
  END SUBROUTINE WA_REDIS_PW


  SUBROUTINE W_REDIS_PW( W)
    IMPLICIT NONE
    TYPE (wavespin), INTENT(IN) :: W
    TYPE (wavedes1)    WDES1          ! descriptor for one k-point
    INTEGER :: K1, ISP

    IF (W%WDES%DO_REDIS) THEN
       DO K1=1,W%WDES%NKPTS
#ifdef MPI
          IF (MOD(K1-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
          CALL SETWDES(W%WDES,WDES1,K1)
          DO ISP=1,W%WDES%ISPIN
             CALL WA_REDIS_PW( ELEMENTS( W, WDES1, ISP))
          ENDDO
       ENDDO
    ENDIF
  END SUBROUTINE W_REDIS_PW


!***********************************************************************
!
! update of vector
!  W2 = W1*a + W2
!
!
!***********************************************************************


  SUBROUTINE W1_DAXPY( W1, SCALE, W2)
    IMPLICIT NONE
    REAL(q) SCALE
    TYPE (wavefun1)    W1, W2

    PROFILING_START('w1_daxpy')

    CALL DAXPY( W1%WDES1%NPL*2, SCALE, W1%CW(1), 1, W2%CW(1), 1)

    IF (W1%WDES1%LGAMMA) THEN
       CALL DAXPY( W1%WDES1%NPRO, SCALE, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
       CALL DAXPY( W1%WDES1%NPRO*2, SCALE, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF

    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          WRITE(0,*) 'internal error in W1_COPY: real space allocation is not correct'
          STOP
       ENDIF

       ! real space wavefunction complex, SCALE real
       CALL DAXPY( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS*2, SCALE, W1%CR(1), 1, W2%CR(1), 1)
    END IF

    PROFILING_STOP('w1_daxpy')

  END SUBROUTINE W1_DAXPY

!
!
  SUBROUTINE W1_GAXPY( W1, SCALE, W2)
    IMPLICIT NONE
    GDEF SCALE
    TYPE (wavefun1)    W1, W2

    PROFILING_START('w1_gaxpy')

    IF (W1%WDES1%LGAMMA) THEN
       ! wavefunction complex, SCALE real
       CALL DAXPY( W1%WDES1%NPL*2, SCALE, W1%CW(1), 1, W2%CW(1), 1)
       CALL DAXPY( W1%WDES1%NPRO, SCALE, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
       ! wavefunction complex, SCALE complex
       CALL ZAXPY( W1%WDES1%NPL, SCALE, W1%CW(1), 1, W2%CW(1), 1)
       CALL ZAXPY( W1%WDES1%NPRO, SCALE, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF

    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          WRITE(0,*) 'internal error in W1_COPY: real space allocation is not correct'
          STOP
       ENDIF

       IF (W1%WDES1%LGAMMA) THEN
          ! real space wavefunction complex, SCALE real
          CALL DAXPY( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS*2, SCALE, W1%CR(1), 1, W2%CR(1), 1)
       ELSE
          ! real space wavefunction complex, SCALE complex
          CALL ZAXPY( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS, SCALE, W1%CR(1), 1, W2%CR(1), 1)
       ENDIF
             
    END IF

    PROFILING_STOP('w1_gaxpy')

  END SUBROUTINE W1_GAXPY

!***********************************************************************
!
! update of vector
!  W1 = W1*a
!
!
!***********************************************************************


  SUBROUTINE W1_DSCAL( W1, SCALE)
    IMPLICIT NONE
    REAL(q) SCALE
    TYPE (wavefun1)    W1

    PROFILING_START('w1_dscal')

    ! since this function can be used to zero out an array we operate
    ! on all elements (dimension) and not only on those that are
    ! actually used
    CALL DSCAL( W1%WDES1%NRPLWV*2, SCALE, W1%CW(1), 1)

    IF (W1%WDES1%LGAMMA) THEN
       CALL DSCAL( W1%WDES1%NPROD, SCALE, W1%CPROJ(1), 1)
    ELSE
       CALL DSCAL( W1%WDES1%NPROD*2, SCALE, W1%CPROJ(1),  1)
    ENDIF

    IF (ASSOCIATED(W1%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          WRITE(0,*) 'internal error in W1_COPY: real space allocation is not correct'
          STOP
       ENDIF
       ! real space wavefunction complex, SCALE real
       CALL DSCAL( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS*2, SCALE, W1%CR(1), 1)
    END IF

    PROFILING_STOP('w1_dscal')

  END SUBROUTINE W1_DSCAL

!***********************************************************************
!
! calculate the dot product between two wavefunctions
!  C=   W1^* x dW2
!
!
!***********************************************************************

  FUNCTION W1_DOT2( W1, W2, CQIJ, CPROJ1, WDES, LATT_CUR, NIONS, DIREC) RESULT (C)
    USE constant
    USE lattice
    IMPLICIT NONE
    TYPE (wavefun1)    W1, W2
    TYPE (wavedes1)    WDES
    GDEF :: C(NIONS+1)
    INTEGER NC                      ! stride of C
    OVERLAP, OPTIONAL :: CQIJ(:,:,:,:) ! optional overlap operator
    COMPLEX(q), EXTERNAL :: ZDOTC
    REAL(q), EXTERNAL ::  DDOT
  ! local    
    GDEF :: CNL
    INTEGER :: LMDIM, NPRO, NPRO_, ISPINOR, ISPINOR_, LMMAXC, NT, NI, NIS, NIONS
    GDEF :: CPROJ1(W1%WDES1%NPROD)
    GDEF  :: CWI(WDES%NPL)
    TYPE (latt)        LATT_CUR
    INTEGER :: DIREC, K
    GDEF  :: CX_, CY_, CZ_

    PROFILING_START('w1_dot2')

    C=0
    DO K=1,WDES%NPL
       CX_=-1.0*(WDES%IGX(K))*W2%CW(K)*CITPI
       CY_=-1.0*(WDES%IGY(K))*W2%CW(K)*CITPI
       CZ_=-1.0*(WDES%IGZ(K))*W2%CW(K)*CITPI
       CWI(K)=CX_*LATT_CUR%B(DIREC,1)+CY_*LATT_CUR%B(DIREC,2)+CZ_*LATT_CUR%B(DIREC,3)
    ENDDO

    IF (W1%WDES1%LGAMMA) THEN
       C(1)=DDOT( 2*WDES%NPL,W1%CW(1),1, CWI(1),1)
    ELSE
       C(1)=ZDOTC( WDES%NPL,W1%CW(1),1, CWI(1),1)
    ENDIF

    IF (PRESENT(CQIJ) .AND. W1%WDES1%LOVERL .AND.  W1%WDES1%NPROD>0 ) THEN
       !CNL =0
       LMDIM = SIZE(CQIJ,1)

       spinor: DO ISPINOR=0,W1%WDES1%NRSPINORS-1
       DO ISPINOR_=0,W1%WDES1%NRSPINORS-1
             
          NPRO =ISPINOR *(W1%WDES1%NPRO/2)
          NPRO_=ISPINOR_*(W1%WDES1%NPRO/2)
          NIS =1
          DO NT=1,W1%WDES1%NTYP
             LMMAXC=W1%WDES1%LMMAX(NT)
             IF (LMMAXC/=0) THEN
                DO NI=NIS,W1%WDES1%NITYP(NT)+NIS-1
                   CNL =0
                   !CALL ECCP_NL(LMDIM,LMMAXC,CQIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W2%CPROJ(NPRO_+1),W1%CPROJ(NPRO+1),CNL)
                   CALL ECCP_NL(LMDIM,LMMAXC,CQIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W2%CPROJ(NPRO_+1),CPROJ1(NPRO+1),CNL)
                   NPRO = LMMAXC+NPRO
                   NPRO_= LMMAXC+NPRO_
                   C(NI+1)=C(NI+1)+CNL
                ENDDO
             ENDIF
             NIS = NIS+W1%WDES1%NITYP(NT)
          ENDDO
       ENDDO
       ENDDO spinor
       !C=C+CNL
    ENDIF

    CALLMPI( M_sum_g(W1%WDES1%COMM_INB, C, NIONS+1))

    PROFILING_STOP('w1_dot2')
      
  END FUNCTION W1_DOT2

!***********************************************************************
!
! calculate the dot product between two wavefunctions
!  C=   W1^* x W2
! this is a substitue for the routine CINDPROD but mind
! the W1 and W2 are interchanged
!
!
!***********************************************************************

  FUNCTION W1_DOT( W1, W2, CQIJ) RESULT (C)
    IMPLICIT NONE
    TYPE (wavefun1)    W1, W2
    GDEF :: C
    INTEGER NC                      ! stride of C
    OVERLAP, OPTIONAL :: CQIJ(:,:,:,:) ! optional overlap operator
    COMPLEX(q), EXTERNAL :: ZDOTC
    REAL(q), EXTERNAL ::  DDOT
  ! local    
    GDEF :: CNL
    INTEGER :: LMDIM, NPRO, NPRO_, ISPINOR, ISPINOR_, LMMAXC, NT, NI, NIS

    PROFILING_START('w1_dot')

    IF (W1%WDES1%LGAMMA) THEN
       C=DDOT( 2* W1%WDES1%NPL,W1%CW(1),1,W2%CW(1),1)
    ELSE
       C=ZDOTC( W1%WDES1%NPL,W1%CW(1),1,W2%CW(1),1)
    ENDIF

    IF (PRESENT(CQIJ) .AND. W1%WDES1%LOVERL .AND.  W1%WDES1%NPROD>0 ) THEN
       CNL =0
       LMDIM = SIZE(CQIJ,1)

       spinor: DO ISPINOR=0,W1%WDES1%NRSPINORS-1
       DO ISPINOR_=0,W1%WDES1%NRSPINORS-1
             
          NPRO =ISPINOR *(W1%WDES1%NPRO/2)
          NPRO_=ISPINOR_*(W1%WDES1%NPRO/2)
          NIS =1
          DO NT=1,W1%WDES1%NTYP
             LMMAXC=W1%WDES1%LMMAX(NT)
             IF (LMMAXC/=0) THEN
                DO NI=NIS,W1%WDES1%NITYP(NT)+NIS-1
                   CALL ECCP_NL(LMDIM,LMMAXC,CQIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W2%CPROJ(NPRO_+1),W1%CPROJ(NPRO+1),CNL)
                   NPRO = LMMAXC+NPRO
                   NPRO_= LMMAXC+NPRO_
                ENDDO
             ENDIF
             NIS = NIS+W1%WDES1%NITYP(NT)
          ENDDO
       ENDDO
       ENDDO spinor
       C=C+CNL
    ENDIF

    CALLMPI( M_sum_g(W1%WDES1%COMM_INB, C, 1))

    PROFILING_STOP('w1_dot')
      
  END FUNCTION W1_DOT


!***********************************************************************
!
! calculate the inproduct between one wavefunction and 
! a set of wavefunctions with possible scales and add the result
! to a third vector
!  C=   WA^* x W1 * SCALEA + C * SCALEC
!     
!***********************************************************************

  SUBROUTINE W1_GEMV( SCALEA, WA, W1, SCALEC, C, NC , CQIJ)
    IMPLICIT NONE
    GDEF :: SCALEA                  ! scaleing constant for WA x W1
    TYPE (wavefuna)    WA
    TYPE (wavefun1)    W1
    GDEF :: SCALEC                  ! scaling constant for C
    GDEF :: C(*)                    ! result
    INTEGER NC                      ! stride of C
    OVERLAP, OPTIONAL :: CQIJ(:,:,:,:) ! optional overlap operator
  ! local
    INTEGER :: LMDIM, I
    GDEF CRESUL(WA%WDES1%NPRO)
    GDEF CTMP(SIZE(WA%CW,2))

    PROFILING_START('w1_gemv')

    IF (W1%WDES1%LGAMMA) THEN
       CALL DGEMV( dtrans, 2* WA%WDES1%NPL, SIZE(WA%CW,2) , done , WA%CW(1,1) , &
            &           2* WA%WDES1%NRPLWV, W1%CW(1) , 1 , dzero,  CTMP(1), 1)
    ELSE
       CALL ZGEMV( ztrans,  WA%WDES1%NPL, SIZE(WA%CW,2) , zone , WA%CW(1,1) , &
            &            WA%WDES1%NRPLWV, W1%CW(1) , 1 , zzero,  CTMP(1), 1)
    ENDIF

    IF (PRESENT(CQIJ) .AND. WA%WDES1%LOVERL .AND.  WA%WDES1%NPROD>0 ) THEN
       LMDIM = SIZE(CQIJ,1)
       CALL OVERL1(WA%WDES1, LMDIM, CQIJ(1,1,1,1), CQIJ(1,1,1,1), 0.0_q, WA%CPROJ(1,1), CRESUL(1))

       IF (W1%WDES1%LGAMMA) THEN
          CALL DGEMV( dtrans, WA%WDES1%NPRO,  SIZE(WA%CW,2), done ,WA%CPROJ(1,1), &
               &           WA%WDES1%NPROD, CRESUL , 1 , done,  CTMP(1), 1)
       ELSE
          CALL ZGEMV( ztrans, WA%WDES1%NPRO,  SIZE(WA%CW,2), zone ,WA%CPROJ(1,1), &
               &           WA%WDES1%NPROD, CRESUL , 1 , zone,  CTMP(1), 1)
       ENDIF

    ENDIF

    CALLMPI( M_sum_g(WA%WDES1%COMM_INB, CTMP, SIZE(WA%CW,2)))

    IF (SCALEC==0) THEN
       DO I=0,SIZE(WA%CW,2)-1
          C(I*NC+1)=CTMP(I+1)*SCALEA
       ENDDO
    ELSE
       DO I=0,SIZE(WA%CW,2)-1
          C(I*NC+1)=C(I*NC+1)*SCALEC+CTMP(I+1)*SCALEA
       ENDDO
    ENDIF

    PROFILING_STOP('w1_gemv')

  END SUBROUTINE W1_GEMV


!***********************************************************************
!
! W1 descriptor from W array
! alternative to SETWAV (returns a W1 descriptor)
! it is somewhat slimmed down to optimize performance
! and returns only the wavefunction and wavefunction character pointers
! the full version remains SETWAV
! 
!***********************************************************************

  FUNCTION W1_FROM_W( W, WDES1, NB, ISP) RESULT (W1)
    IMPLICIT NONE
    INTEGER NB, ISP
    TYPE (wavespin) W
    TYPE (wavefun1) W1
    TYPE (wavedes1), TARGET :: WDES1
    INTEGER NK
    
    NK=WDES1%NK
#ifdef boundary_check
    IF (NB<=0 .OR. NB> SIZE(W%CPTWFP,2)) THEN
       WRITE(0,*) 'internal error in W1_FROM_W: bounds exceed ',NB, SIZE(W%CPTWFP,2)
       STOP
    ENDIF
#endif    
    W1%CPTWFP=>W%CPTWFP(:,NB,NK,ISP)
    W1%CPROJ =>W%CPROJ(:,NB,NK,ISP)
    W1%WDES1 => WDES1
    NULLIFY(W1%CR)
    W1%LDO=.TRUE.
  END FUNCTION W1_FROM_W

!***********************************************************************
!
! W1 descriptor from WA descriptor
!     
!***********************************************************************

  FUNCTION W1_FROM_WA( WA, N1) RESULT (W1)
    IMPLICIT NONE
    INTEGER N1
    TYPE (wavefun1) W1
    TYPE (wavefuna) WA
#ifdef boundary_check
    IF (N1<=0 .OR. N1> SIZE(WA%CPTWFP,2)) THEN
       WRITE(0,*) 'internal error in W1_FROM_WA: bounds exceed ',N1, SIZE(WA%CPTWFP,2)
       STOP
    ENDIF
#endif    
    W1%CPTWFP=>WA%CPTWFP(:,N1)
    W1%CPROJ =>WA%CPROJ(:,N1)
    W1%WDES1 =>WA%WDES1
    NULLIFY(W1%CR)
    W1%LDO=.TRUE.
  END FUNCTION W1_FROM_WA

!***********************************************************************
!
! subindex a WA array
! return WA(N1:N2)
!     
!***********************************************************************

  FUNCTION WA_FROM_WA( WA, N1, N2) RESULT (W1)
    IMPLICIT NONE
    INTEGER N1, N2
    TYPE (wavefuna) W1
    TYPE (wavefuna) WA
#ifdef boundary_check
    IF (N1<=0 .OR. N1> SIZE(WA%CPTWFP,2)) THEN
       WRITE(0,*) 'internal error in WA_FROM_WA: bounds exceed ',N1, SIZE(WA%CPTWFP,2)
       STOP
    ENDIF
    IF (N2<N1 .OR. N2> SIZE(WA%CPTWFP,2)) THEN
       WRITE(0,*) 'internal error in WA_FROM_WA: bounds exceed ',N1,N2, SIZE(WA%CPTWFP,2)
       STOP
    ENDIF
#endif    
    W1%CPTWFP=>WA%CPTWFP(:,N1:N2)
    W1%CPROJ =>WA%CPROJ(:,N1:N2)
    W1%WDES1 =>WA%WDES1
  END FUNCTION WA_FROM_WA

!***********************************************************************
!
! one element from a WA array using two indices
! return WA(N1,N2)
!     
!***********************************************************************

  FUNCTION W1_FROM_WA2( WA, N1, N2) RESULT (W1)
    IMPLICIT NONE
    INTEGER N1, N2
    TYPE (wavefun1) W1
    TYPE (wavefuna) WA
#ifdef boundary_check
    IF (N1<=0 .OR. N1> WA%FIRST_DIM) THEN
       WRITE(0,*) 'internal error in W1_FROM_WA: bounds exceed ',N1, WA%FIRST_DIM
       STOP
    ENDIF
    IF (N1+ (N2-1)*WA%FIRST_DIM<=0 .OR. N1+ (N2-1)*WA%FIRST_DIM> SIZE(WA%CPTWFP,2)) THEN
       WRITE(0,*) 'internal error in W1_FROM_WA: bounds exceed ',N1+ (N2-1)*WA%FIRST_DIM,SIZE(WA%CPTWFP,2)
       STOP
    ENDIF
#endif
    W1%CPTWFP=>WA%CPTWFP(:,N1+ (N2-1)*WA%FIRST_DIM)
    W1%CPROJ =>WA%CPROJ (:,N1+ (N2-1)*WA%FIRST_DIM)
    W1%WDES1 =>WA%WDES1
    NULLIFY(W1%CR)
    W1%LDO=.TRUE.
  END FUNCTION W1_FROM_WA2

!***********************************************************************
!
! subindex a WA array using two indices
! return WA(N1:N12,N2)
!     
!***********************************************************************

  FUNCTION WA_FROM_WA2( WA, N1, N12, N2) RESULT (W1)
    IMPLICIT NONE
    INTEGER N1, N12, N2
    TYPE (wavefuna) W1
    TYPE (wavefuna) WA
#ifdef boundary_check
    IF (N1<=0 .OR. N1> WA%FIRST_DIM) THEN
       WRITE(0,*) 'internal error in W1_FROM_WA: bounds exceed ',N1, WA%FIRST_DIM
       STOP
    ENDIF
    IF (N1+ (N2-1)*WA%FIRST_DIM<=0 .OR. N1+ (N2-1)*WA%FIRST_DIM> SIZE(WA%CPTWFP,2)) THEN
       WRITE(0,*) 'internal error in W1_FROM_WA: bounds exceed ',N1+ (N2-1)*WA%FIRST_DIM,SIZE(WA%CPTWFP,2)
       STOP
    ENDIF
    IF (N12<N1 .OR. N12> WA%FIRST_DIM) THEN
       WRITE(0,*) 'internal error in W1_FROM_WA: bounds exceed ',N1,N12, WA%FIRST_DIM
       STOP
    ENDIF
    IF (N12+ (N2-1)*WA%FIRST_DIM<N1+ (N2-1)*WA%FIRST_DIM .OR. N12+ (N2-1)*WA%FIRST_DIM> SIZE(WA%CPTWFP,2)) THEN
       WRITE(0,*) 'internal error in W1_FROM_WA: bounds exceed ',N12+ (N2-1)*WA%FIRST_DIM,SIZE(WA%CPTWFP,2)
       STOP
    ENDIF

#endif
    W1%CPTWFP=>WA%CPTWFP(:,N1+ (N2-1)*WA%FIRST_DIM:N12+ (N2-1)*WA%FIRST_DIM)
    W1%CPROJ =>WA%CPROJ (:,N1+ (N2-1)*WA%FIRST_DIM:N12+ (N2-1)*WA%FIRST_DIM)
    W1%WDES1 =>WA%WDES1
  END FUNCTION WA_FROM_WA2


!***********************************************************************
!
! get a WA structure from a wavefunction array
! this is equivalent to the low level routine SETWAVA
!
!***********************************************************************

  FUNCTION WA_FROM_W( W, WDES1, ISP ) RESULT (WA)
    IMPLICIT NONE
    TYPE (wavespin), INTENT(IN) :: W
    TYPE (wavedes1), TARGET ::  WDES1
    INTEGER ISP
    TYPE (wavefuna) ::  WA

    WA%CPTWFP=>W%CPTWFP(:,:,WDES1%NK,ISP)
    WA%CPROJ =>W%CPROJ(:,:,WDES1%NK,ISP)
    WA%FERWE =>W%FERWE(:,WDES1%NK,ISP)
    WA%AUX   =>W%AUX  (:,WDES1%NK,ISP)
    WA%CELEN =>W%CELEN(:,WDES1%NK,ISP)
    WA%WDES1 =>WDES1
    ! remember spin index
    WA%ISP =ISP
    ! set redistributed wavefunction indices
    IF (WDES1%DO_REDIS) THEN
       CALL SET_WPOINTER(WA%CW_RED,    WDES1%NRPLWV_RED, W%WDES%NB_TOT, WA%CW(1,1))
       CALL SET_GPOINTER(WA%CPROJ_RED, WDES1%NPROD_RED,  W%WDES%NB_TOT, WA%CPROJ(1,1))
    ELSE
       WA%CW_RED=>WA%CPTWFP
       WA%CPROJ_RED=>WA%CPROJ
    ENDIF
  END FUNCTION WA_FROM_W


!=======================================================================
!  create storage for a wavefunction array WA
!=======================================================================

  SUBROUTINE NEWWAVA(WA, WDES1, NDIM, NDIM2)
    USE prec
    IMPLICIT NONE
    TYPE (wavefuna), TARGET :: WA
    TYPE (wavedes1), TARGET :: WDES1
    INTEGER :: NDIM
    INTEGER, OPTIONAL ::  NDIM2
    IF (PRESENT(NDIM2)) THEN
       ALLOCATE(WA%CPTWFP(WDES1%NRPLWV,NDIM*NDIM2),WA%CPROJ(WDES1%NPROD,NDIM*NDIM2), &
            WA%CELEN(NDIM*NDIM2), WA%FERWE(NDIM*NDIM2), WA%AUX(NDIM*NDIM2))
       WA%FIRST_DIM=NDIM
       WA%CPTWFP=0
       WA%CPROJ=0
    ELSE
       ALLOCATE(WA%CPTWFP(WDES1%NRPLWV,NDIM),WA%CPROJ(WDES1%NPROD,NDIM), &
            WA%CELEN(NDIM), WA%FERWE(NDIM), WA%AUX(NDIM))
       WA%FIRST_DIM=0
    ENDIF
    WA%ISP=-1
    ! set redistributed wavefunction indices
    IF (WDES1%DO_REDIS) THEN
       CALL SET_WPOINTER(WA%CW_RED,    WDES1%NRPLWV_RED, NDIM*WDES1%NB_PAR, WA%CW(1,1))
       CALL SET_GPOINTER(WA%CPROJ_RED, WDES1%NPROD_RED,  NDIM*WDES1%NB_PAR, WA%CPROJ(1,1))
    ELSE
       WA%CW_RED=>WA%CPTWFP
       WA%CPROJ_RED=>WA%CPROJ
    ENDIF
    WA%WDES1=>WDES1

  END SUBROUTINE NEWWAVA


!=======================================================================
!  create storage for one wavefunction W array
!  to store non local part only
!=======================================================================

  SUBROUTINE NEWWAVA_PROJ(WA, WDES1, NDIM)
    USE prec
    USE wave_mpi
    IMPLICIT NONE
    TYPE (wavefuna), TARGET :: WA
    TYPE (wavedes1), TARGET :: WDES1
    INTEGER, OPTIONAL :: NDIM
    IF (PRESENT(NDIM)) THEN
       ALLOCATE(WA%CPROJ(WDES1%NPROD,NDIM))
    ELSE
       ALLOCATE(WA%CPROJ(WDES1%NPROD,WDES1%NBANDS))
    ENDIF

    WA%FIRST_DIM=0
    WA%ISP=-1
    NULLIFY(WA%CPTWFP)
    NULLIFY(WA%CW_RED)
    NULLIFY(WA%FERWE)
    NULLIFY(WA%AUX  )
    NULLIFY(WA%CELEN)

    WA%WDES1=>WDES1

    IF (WDES1%DO_REDIS) THEN
       IF (PRESENT(NDIM)) THEN
          CALL SET_GPOINTER(WA%CPROJ_RED, WDES1%NPROD_RED,  NDIM*WDES1%NB_PAR, WA%CPROJ(1,1))
       ELSE
          CALL SET_GPOINTER(WA%CPROJ_RED, WDES1%NPROD_RED,  WDES1%NB_TOT, WA%CPROJ(1,1))
       ENDIF
    ELSE
       WA%CPROJ_RED=>WA%CPROJ
    ENDIF

  END SUBROUTINE NEWWAVA_PROJ


!=======================================================================
!  set one single wavefunction array (WA) from an array of wavefunctions
!=======================================================================

  SUBROUTINE SETWAVA(W, WA, WDES1, ISP)
    USE prec
    IMPLICIT NONE
    INTEGER NB,ISP
    TYPE (wavespin) W
    TYPE (wavefuna) WA
    TYPE (wavedes1), TARGET :: WDES1
    INTEGER NK

    NK=WDES1%NK

    WA%CPTWFP=>W%CPTWFP(:,:,NK,ISP)
    WA%CPROJ =>W%CPROJ(:,:,NK,ISP)
    WA%FERWE =>W%FERWE(:,NK,ISP)
    WA%AUX   =>W%AUX  (:,NK,ISP)
    WA%CELEN =>W%CELEN(:,NK,ISP)
    WA%WDES1 =>WDES1
    ! one dimensional indexing assumed
    WA%FIRST_DIM=0
    ! remember spin index
    WA%ISP =ISP
    ! set redistributed wavefunction indices
    IF (WDES1%DO_REDIS) THEN
       CALL SET_WPOINTER(WA%CW_RED,    WA%WDES1%NRPLWV_RED, W%WDES%NB_TOT, WA%CW(1,1))
       CALL SET_GPOINTER(WA%CPROJ_RED, WA%WDES1%NPROD_RED,  W%WDES%NB_TOT, WA%CPROJ(1,1))
    ELSE
       WA%CW_RED=>WA%CPTWFP
       WA%CPROJ_RED=>WA%CPROJ
    ENDIF

  END SUBROUTINE SETWAVA

!=======================================================================
!  destroy storage for a wavefunctionarray WA
!=======================================================================

      SUBROUTINE DELWAVA(WA)
      USE prec
      IMPLICIT NONE
      TYPE (wavefuna) WA

      IF (.NOT. ASSOCIATED(WA%CPTWFP) .OR. .NOT. ASSOCIATED(WA%FERWE) .OR. &
          .NOT. ASSOCIATED(WA%CELEN) .OR. .NOT. ASSOCIATED(WA%AUX)) THEN
         WRITE(*,*)'internal error in DELWAVA: not all enities are associated, try DELWAVA_PROJ'
         STOP
      ENDIF
      DEALLOCATE(WA%CPTWFP,WA%CPROJ,WA%CELEN, WA%FERWE, WA%AUX)

      END SUBROUTINE DELWAVA


!=======================================================================
!  destroy storage for a wavefunctionarray WA
!  wavefunction character  only
!=======================================================================

      SUBROUTINE DELWAVA_PROJ(WA)
      USE prec
      IMPLICIT NONE
      TYPE (wavefuna) WA

      IF (ASSOCIATED(WA%CPTWFP) .OR. ASSOCIATED(WA%FERWE) .OR. &
          ASSOCIATED(WA%CELEN) .OR. ASSOCIATED(WA%AUX)) THEN
         WRITE(*,*)'internal error in DELWAVA_PROJ: enities are associated, try DELWAVA'
         STOP
      ENDIF
      DEALLOCATE(WA%CPROJ)

    END SUBROUTINE DELWAVA_PROJ


!***********************************************************************
!
! index the wavefunction array or character array in WA
!     
!***********************************************************************

  FUNCTION PCW( WA, N1, N2)
    INTEGER N1, N2
    TYPE (wavefuna) WA
    COMPLEX(q), POINTER :: PCW(:)
#ifdef boundary_check
    IF (N1<=0 .OR. N1> WA%FIRST_DIM) THEN
       WRITE(0,*) 'internal error in W1_FROM_WA: bounds exceed ',N1, WA%FIRST_DIM
       STOP
    ENDIF
    IF (N1+ (N2-1)*WA%FIRST_DIM<=0 .OR. N1+ (N2-1)*WA%FIRST_DIM> SIZE(WA%CPTWFP,2)) THEN
       WRITE(0,*) 'internal error in W1_FROM_WA: bounds exceed ',N1+ (N2-1)*WA%FIRST_DIM,SIZE(WA%CPTWFP,2)
       STOP
    ENDIF
#endif
    PCW=>WA%CPTWFP(:,N1+ (N2-1)*WA%FIRST_DIM)
  END FUNCTION PCW

  FUNCTION PCPROJ( WA, N1, N2)
    INTEGER N1, N2
    TYPE (wavefuna) WA
    GDEF, POINTER :: PCPROJ(:)
#ifdef boundary_check
    IF (N1<=0 .OR. N1> WA%FIRST_DIM) THEN
       WRITE(0,*) 'internal error in W1_FROM_WA: bounds exceed ',N1, WA%FIRST_DIM
       STOP
    ENDIF
    IF (N1+ (N2-1)*WA%FIRST_DIM<=0 .OR. N1+ (N2-1)*WA%FIRST_DIM> SIZE(WA%CPTWFP,2)) THEN
       WRITE(0,*) 'internal error in W1_FROM_WA: bounds exceed ',N1+ (N2-1)*WA%FIRST_DIM,SIZE(WA%CPTWFP,2)
       STOP
    ENDIF
#endif
    PCPROJ=>WA%CPROJ(:,N1+ (N2-1)*WA%FIRST_DIM)
  END FUNCTION PCPROJ



!************************* SUBROUTINE ORTHON ***************************
!
! orthogonalize a wavefunction W1 to all other bands
! including the current band
! the subroutine uses BLAS 3 calls,
!
!***********************************************************************
    
  SUBROUTINE ORTHON(NK, W, W1, CQIJ, ISP)
    IMPLICIT NONE

    INTEGER NK
    TYPE (wavespin)   W
    TYPE (wavefun1)   W1
    LOGICAL LOVERL
    OVERLAP CQIJ(:,:,:,:)
    INTEGER ISP
  ! local
    GDEF :: CPRO(W%WDES%NBANDS),CWORK(W%WDES%NPRO)
    REAL(q) :: WFMAG
    INTEGER :: I

#ifdef MPI
    IF (W%WDES%COMM_KIN%NCPU /= W%WDES%COMM_INB%NCPU) THEN
       WRITE(*,*)'internal error: ORTHON does not support band-par.'
       STOP
    ENDIF
#endif
    IF (W1%WDES1%LOVERL) THEN
       CALL OVERL1(W1%WDES1, SIZE(CQIJ,1),CQIJ(1,1,1,ISP),CQIJ(1,1,1,ISP), 0.0_q, W1%CPROJ(1),CWORK(1))
    ENDIF

    CALL GGEMV( trans , m_ W%WDES%NPLWKP(NK) , W%WDES%NBANDS ,one , W%CW(1,1,NK,ISP), &
         &             m_ W%WDES%NRPLWV, W1%CW(1) , 1 , zero ,  CPRO(1), 1)

    IF (W1%WDES1%LOVERL) THEN
       IF (W%WDES%NPRO /= 0) &
            CALL GGEMV( trans ,  W%WDES%NPRO , W%WDES%NBANDS ,one , W%CPROJ(1,1,NK,ISP) , &
            W%WDES%NPROD, CWORK(1), 1 , one ,  CPRO(1), 1)
    ENDIF

    CALLMPI( M_sum_g(W%WDES%COMM_KIN, CPRO(1), W%WDES%NBANDS))

    CALL GGEMM( 'N', 'N' , m_ W%WDES%NPLWKP(NK) , 1 , W%WDES%NBANDS , -one , &
         W%CW(1,1,NK,ISP), m_ W%WDES%NRPLWV , CPRO(1) , W%WDES%NBANDS , &
         one , W1%CW(1) , m_ W%WDES%NRPLWV )

    IF (W%WDES%NPRO /= 0) &
         CALL GGEMM( 'N', 'N' ,  W%WDES%NPRO , 1 , W%WDES%NBANDS  , -one , &
         W%CPROJ(1,1,NK,ISP) ,  W%WDES%NPROD , CPRO(1) , W%WDES%NBANDS , &
         one , W1%CPROJ(1) ,  W%WDES%NPROD  )

  END SUBROUTINE ORTHON

  SUBROUTINE ORTHON1P(NK, W, W1, CQIJ, ISP, NB)
    IMPLICIT NONE

    INTEGER NK
    TYPE (wavespin)   W
    TYPE (wavefun1)   W1(:)
    LOGICAL LOVERL
    OVERLAP CQIJ(:,:,:,:)
    INTEGER ISP,NB(:)
  ! local
    GDEF :: CPRO(W%WDES%NBANDS,SIZE(W1)),CWORK(W%WDES%NPRO,SIZE(W1))
    GDEF,ALLOCATABLE :: CW1(:,:),CPROJ1(:,:)
    REAL(q) :: WFMAG
    INTEGER :: I,BW

    BW=128

#ifdef MPI
    IF (W%WDES%COMM_KIN%NCPU /= W%WDES%COMM_INB%NCPU) THEN
       WRITE(*,*)'internal error: ORTHON does not support band-par.'
       STOP
    ENDIF
#endif

    ALLOCATE(CW1(SIZE(W1(1)%CW,1),SIZE(W1)))
    DO I=1,SIZE(W1)
       CW1(:,I)=W1(I)%CW
    END DO

    DO I=1,SIZE(W1)
       IF (W1(I)%WDES1%LOVERL)THEN
          CALL OVERL1(W1(I)%WDES1, SIZE(CQIJ,1),CQIJ(1,1,1,ISP),CQIJ(1,1,1,ISP), 0.0_q, W1(I)%CPROJ(1),CWORK(1,I))
       ELSE
          CWORK(:,I)=0.0_q
       END IF
    END DO

    CALL GGEMM( trans , 'N', W%WDES%NBANDS, size(w1), m_ W%WDES%NPLWKP(NK), &
         & one , W%CW(1,1,NK,ISP), m_ W%WDES%NRPLWV, CW1(1,1) , m_ W%WDES%NRPLWV, zero ,  CPRO(1,1), W%WDES%NBANDS)


    IF (W%WDES%NPRO /= 0) &
         CALL GGEMM( trans , 'N', W%WDES%NBANDS , size(w1), W%WDES%NPRO , one ,&
         & W%CPROJ(1,1,NK,ISP) , W%WDES%NPROD, CWORK(1,1), W%WDES%NPRO, one,  CPRO(1,1), W%WDES%NBANDS)

    CALLMPI( M_sum_g(W%WDES%COMM_KIN, CPRO(1,1), W%WDES%NBANDS*SIZE(W1)))

    DO I=1,SIZE(W1)
       CPRO(NB(I),I)=0.d0
       IF(NB(I)-BW>=1)CPRO(1:NB(I)-BW,I)=0.d0
       IF(NB(I)+BW<=SIZE(CPRO,1))CPRO(NB(I)+BW:SIZE(CPRO,1),I)=0.d0
    END DO

    CALL GGEMM( 'N', 'N' , m_ W%WDES%NPLWKP(NK) , SIZE(W1) , W%WDES%NBANDS , -one , &
         W%CW(1,1,NK,ISP), m_ W%WDES%NRPLWV , CPRO(1,1) , W%WDES%NBANDS , &
         one , CW1(1,1) , m_ W%WDES%NRPLWV )

    IF (W%WDES%NPRO /= 0)THEN
       DO I=1,SIZE(W1)
          CWORK(:,i)=W1(i)%CPROJ(:)
        END DO
       CALL GGEMM( 'N', 'N' ,  W%WDES%NPRO , SIZE(W1) , W%WDES%NBANDS  , -one , &
            W%CPROJ(1,1,NK,ISP) ,  W%WDES%NPROD , CPRO(1,1) , W%WDES%NBANDS , &
            one , CWORK(1,1) ,  W%WDES%NPROD  )
    end IF

    DO I=1,SIZE(W1)
       W1(I)%CW=CW1(:,I)
       W1(I)%CPROJ=CWORK(:,I)
    END DO

    DEALLOCATE(CW1)

  END SUBROUTINE ORTHON1P

!************************* SUBROUTINE CNORMN  **************************
!
! this subroutine normalises a wavefunction
! subroutine is not important for performance
!
!***********************************************************************

  SUBROUTINE CNORMN(W, CQIJ, ISP, WSCAL)
    IMPLICIT NONE
    TYPE (wavefun1)    W

    OVERLAP CQIJ(:,:,:,:)
    INTEGER   ISP
    REAL(q) WSCAL
  ! local
    GDEF      CP
    COMPLEX(q), EXTERNAL ::  ZDOTC
    REAL(q), EXTERNAL ::  DDOT
    REAL(q) WFMAG
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NIS, NI, LMMAXC

    WFMAG=ZDOTC(W%WDES1%NPL,W%CW(1),1,W%CW(1),1)
!=======================================================================
! if necessary caclulate <w| P |w>
!=======================================================================
    IF (W%WDES1%LOVERL) THEN
       CP  =0

       spinor: DO ISPINOR=0,W%WDES1%NRSPINORS-1
          DO ISPINOR_=0,W%WDES1%NRSPINORS-1

             NPRO =ISPINOR *(W%WDES1%NPRO/2)
             NPRO_=ISPINOR_*(W%WDES1%NPRO/2)

             NIS =1
             DO NT=1,W%WDES1%NTYP
                LMMAXC=W%WDES1%LMMAX(NT)
                IF (LMMAXC==0) GOTO 230

                DO NI=NIS,W%WDES1%NITYP(NT)+NIS-1
                   CALL ECCP_NL(SIZE(CQIJ,1),LMMAXC,CQIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W%CPROJ(NPRO_+1),W%CPROJ(NPRO+1),CP)
                   NPRO = LMMAXC+NPRO
                   NPRO_= LMMAXC+NPRO_
                ENDDO
230             NIS = NIS+W%WDES1%NITYP(NT)
             ENDDO
          ENDDO
       ENDDO spinor

       WFMAG=WFMAG+CP
    ENDIF

    CALLMPI( M_sum_d(W%WDES1%COMM_INB, WFMAG, 1))

    !-----check that it is non-zero
    IF(WFMAG<=0) THEN
!=======================================================================
! if it is smaller zero write a warning
!=======================================================================
#ifdef MPI
       IF (W%WDES1%COMM_INB%NODE_ME == W%WDES1%COMM_INB%IONODE) THEN
#endif
          WRITE(*,*)'WARNING: CNORMN: search vector ill defined'
#ifdef MPI
       ENDIF
#endif
       WSCAL= -1._q/SQRT(-WFMAG)
    ELSE
       WSCAL= 1._q/SQRT(WFMAG)
    ENDIF
    CALL ZDSCAL( W%WDES1%NPL ,WSCAL,W%CW(1),1)
    CALL GDSCAL( W%WDES1%NPRO,WSCAL,W%CPROJ(1),1)

  END SUBROUTINE CNORMN


!************************* SUBROUTINE CNORMN_REAL **********************
!
! performs operations on real space part of wavefunction
! after a call to CPROJCN and CNORMN
!
!***********************************************************************

  SUBROUTINE CNORMN_REAL(W, W1, ISP, WSCAL, CSCPD )
    IMPLICIT NONE
    TYPE (wavefun1)    W, W1

    INTEGER   ISP
    REAL(q) WSCAL
    COMPLEX(q) :: CSCPD
  ! local
    INTEGER ISPINOR, K, KK

    DO ISPINOR=0,W%WDES1%NRSPINORS-1
       DO K=1,W%WDES1%GRID%RL%NP
          KK=K+ISPINOR*W%WDES1%GRID%MPLWV
          W%CR(KK)=(W%CR(KK)-CSCPD*W1%CR(KK))*WSCAL
       ENDDO
    ENDDO

  END SUBROUTINE CNORMN_REAL


!************************* SUBROUTINE CNORMA  **************************
!
! this subroutine calculates the norm of a wavefunction
!
!
!***********************************************************************

  SUBROUTINE CNORMA(W, CQIJ, ISP, WSCAL)
    IMPLICIT NONE

    TYPE (wavefun1)  W
    OVERLAP   CQIJ(:,:,:,:)
    INTEGER   ISP
    REAL(q) :: WSCAL
  ! local
    GDEF      CP
    REAL(q) :: WFMAG
    COMPLEX(q), EXTERNAL ::  ZDOTC
    REAL(q), EXTERNAL ::  DDOT
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NIS, NI, LMMAXC

    PROFILING_START('cnorma')

    WFMAG=ZDOTC(W%WDES1%NPL,W%CW(1),1,W%CW(1),1)
    CP  =0
!=======================================================================
! if necessary caclulate <w| P |w>
!=======================================================================
    IF (W%WDES1%LOVERL) THEN

       NPRO=0

       spinor: DO ISPINOR=0,W%WDES1%NRSPINORS-1
          DO ISPINOR_=0,W%WDES1%NRSPINORS-1

             NPRO =ISPINOR *(W%WDES1%NPRO/2)
             NPRO_=ISPINOR_*(W%WDES1%NPRO/2)
             NIS =1
             DO NT=1,W%WDES1%NTYP
                LMMAXC=W%WDES1%LMMAX(NT)
                IF (LMMAXC==0) GOTO 230

                DO NI=NIS,W%WDES1%NITYP(NT)+NIS-1
                   CALL ECCP_NL(SIZE(CQIJ,1),LMMAXC,CQIJ(1,1,NI,ISP+ISPINOR_+2*ISPINOR),W%CPROJ(NPRO_+1),W%CPROJ(NPRO+1),CP)
                   NPRO = LMMAXC+NPRO
                   NPRO_= LMMAXC+NPRO_
                ENDDO
230             NIS = NIS+W%WDES1%NITYP(NT)
             ENDDO
          ENDDO
       ENDDO spinor

       WFMAG=WFMAG+CP
    ENDIF

    CALLMPI( M_sum_d(W%WDES1%COMM_INB, WFMAG, 1))
    !-----check that it is non-zero
    IF(WFMAG<=0) THEN
!=======================================================================
! if it is smaller zero write a warning
!=======================================================================
#ifdef MPI
       IF (W%WDES1%COMM_INB%NODE_ME == W%WDES1%COMM_INB%IONODE) THEN
#endif
          WRITE(*,*)'WARNING: CNORMN: search vector ill defined'
#ifdef MPI
       ENDIF
#endif
       WSCAL= -1._q/SQRT(-WFMAG)
    ELSE
       WSCAL= 1._q/SQRT(WFMAG)
    ENDIF

    PROFILING_STOP('cnorma')

  END SUBROUTINE CNORMA

!************************* SUBROUTINE CINPROD  *************************
!
! this subroutine calculates the inproduct between two wavefunctions
!  <W2 | S | W1>      = W2(G)*  W1(G) +  \sum_i W2_i* Q_ij W1_j
!***********************************************************************

  SUBROUTINE CINPROD(W1,W2,CQIJ,CWFMAG)
    IMPLICIT NONE

    TYPE (wavefun1)    W1
    TYPE (wavefun1)    W2

    OVERLAP   CQIJ(:,:,:,:)
    COMPLEX(q) CWFMAG
  ! local
    GDEF      CP
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NIS, NI, LMMAXC
    COMPLEX(q), EXTERNAL :: ZDOTC
    REAL(q), EXTERNAL ::  DDOT

    CWFMAG=ZDOTC(W1%WDES1%NPL,W2%CW(1),1,W1%CW(1),1)
    CP  =0
!=======================================================================
! if necessary caclulate <w| P |w>
!=======================================================================
    IF (W1%WDES1%LOVERL) THEN
       NPRO=0

       spinor: DO ISPINOR=0,W1%WDES1%NRSPINORS-1
          DO ISPINOR_=0,W1%WDES1%NRSPINORS-1

             NPRO =ISPINOR *(W1%WDES1%NPRO/2)
             NPRO_=ISPINOR_*(W1%WDES1%NPRO/2)

             NIS =1
             DO NT=1,W1%WDES1%NTYP
                LMMAXC=W1%WDES1%LMMAX(NT)
                IF (LMMAXC==0) GOTO 230

                DO NI=NIS,W1%WDES1%NITYP(NT)+NIS-1
                   CALL ECCP_NL(SIZE(CQIJ,1),LMMAXC,CQIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CP)
                   NPRO = LMMAXC+NPRO
                   NPRO_= LMMAXC+NPRO_
                ENDDO
230             NIS = NIS+W1%WDES1%NITYP(NT)
             ENDDO
          ENDDO
       ENDDO spinor

       CWFMAG=CWFMAG+CP
    ENDIF

    CALLMPI( M_sum_d(W1%WDES1%COMM_INB, CWFMAG, 2))
  END SUBROUTINE CINPROD


!************************* SUBROUTINE PROJCN ***************************
!
! this subroutine projects out from one wavefunction
! CF another wavefunction  CPRO
! subroutine is not important for performance
!
!***********************************************************************

  SUBROUTINE PROJCN(W1, W2, CQIJ, ISP, CSCPD)
    IMPLICIT NONE

    TYPE (wavefun1)    W1,W2
    GDEF    CADD
    OVERLAP CQIJ(:,:,:,:)
    INTEGER :: ISP
    COMPLEX(q) :: CSCPD
    ! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NIS, NI, LMMAXC
    COMPLEX(q), EXTERNAL ::  ZDOTC
    REAL(q), EXTERNAL ::  DDOT

    CSCPD= GREAL(ZDOTC(W1%WDES1%NPL,W2%CW(1),1,W1%CW(1),1))
!=======================================================================
! if necessary caclulate <p| P |w>
!=======================================================================
    IF (W1%WDES1%LOVERL) THEN
       CADD=0
       NPRO=0
       spinor: DO ISPINOR=0,W1%WDES1%NRSPINORS-1
          DO ISPINOR_=0,W1%WDES1%NRSPINORS-1

             NPRO =ISPINOR *(W1%WDES1%NPRO/2)
             NPRO_=ISPINOR_*(W1%WDES1%NPRO/2)

             NIS =1
             DO  NT=1,W1%WDES1%NTYP
                LMMAXC=W1%WDES1%LMMAX(NT)
                IF (LMMAXC==0) GOTO 230

                DO NI=NIS,W1%WDES1%NITYP(NT)+NIS-1
                   CALL ECCP_NL(SIZE(CQIJ,1),LMMAXC,CQIJ(1,1,NI,ISP+ISPINOR_+2*ISPINOR),W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CADD)
                   NPRO = LMMAXC+NPRO
                   NPRO_= LMMAXC+NPRO_
                ENDDO
230             NIS = NIS+W1%WDES1%NITYP(NT)
             ENDDO
          ENDDO
       ENDDO spinor

       CSCPD=GREAL(CSCPD+CADD)
    ENDIF
!=======================================================================
! performe orthogonalisations
!=======================================================================
    CALLMPI( M_sum_z(W1%WDES1%COMM_INB, CSCPD, 1))

    CALL ZAXPY(W1%WDES1%NPL ,-CSCPD,W2%CW(1)   ,1,W1%CW(1)   ,1)
    IF (W1%WDES1%LGAMMA) THEN
       CALL DAXPY(W1%WDES1%NPRO,-CSCPD,W2%CPROJ(1),1,W1%CPROJ(1),1)
    ELSE
       CALL ZAXPY(W1%WDES1%NPRO,-CSCPD,W2%CPROJ(1),1,W1%CPROJ(1),1)
    ENDIF

    RETURN
  END SUBROUTINE PROJCN


!************************ SUBROUTINE W1_GATHER ************************
!
! This subroutine gathers a set of wavefunctions starting 
! from band NB1 until NB2 to all nodes
! 
!**********************************************************************

  SUBROUTINE W1_GATHER( W, NB1, NB2, ISP, W1)
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed

  ! local
    INTEGER :: NN, N, NLOC, NCPU

    NCPU=W%WDES%NB_PAR

    DO N=NB1,NB2
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       CALL W1_COPY( ELEMENT( W, W1(NN)%WDES1, N, ISP), W1(NN) )
       CALL FFTWAV_W1( W1(NN))
    ENDDO

    NLOC=(NB2-NB1+1)*W%WDES%NB_PAR

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO NN=1,NLOC
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CW(1), &
               SIZE(W1(NN)%CW),MOD(NN-1,NCPU)+1)
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(NN-1,NCPU)+1)
#ifndef gammareal
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
               W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
               W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#endif
       ENDDO

#ifdef MPI_barrier_after_bcast
      CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif
    ENDIF
#endif
  END SUBROUTINE W1_GATHER

!
!  new version that uses MPI_Allgather
!  this seems to be even slower then the version above
!
  SUBROUTINE W1_GATHER_NEW( W, NB1, NB2, ISP, W1)
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed

    COMPLEX(q) :: CWBUFF(W%WDES%NRPLWV,W%WDES%NB_PAR)        ! receive buffer
    GDEF       :: CPROJBUFF(W%WDES%NPROD,W%WDES%NB_PAR)      ! receive buffer for projectors
    COMPLEX(q) :: CRBUFF(W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,W%WDES%NB_PAR) ! receive buffer for real space part
  ! local
    INTEGER :: NN, N, NLOC, NCPU
    INTEGER :: ierr          ! error for mpi

    NCPU=W%WDES%NB_PAR
    
    IF (SIZE(W1(1)%CW) /= W%WDES%NRPLWV) THEN
       WRITE (*,*) 'internal error in W1_GATHER: sizes inconsistent ', SIZE(W1(1)%CW) /= W%WDES%NRPLWV
       STOP
    ENDIF

    DO N=NB1,NB2
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       CALL W1_COPY( ELEMENT( W, W1(NN)%WDES1, N, ISP), W1(NN) )
       CALL FFTWAV_W1( W1(NN))

    ! distribute W1 to all nodes
#ifdef MPI
       IF (W%WDES%DO_REDIS) THEN
          CALL MPI_Allgather(W1(NN)%CR(1), W%WDES%GRID%MPLWV*W%WDES%NRSPINORS, MPI_double_complex, CRBUFF(1,1), W%WDES%GRID%MPLWV*W%WDES%NRSPINORS, & 
               MPI_double_complex, W%WDES%COMM_INTER%MPI_COMM, ierr)
          IF ( ierr /= MPI_success ) CALL M_stop_ierr('ERROR: M_Allgather returns',ierror)

          CALL MPI_Allgather(W1(NN)%CW(1), W%WDES%NRPLWV , MPI_double_complex, CWBUFF(1,1), W%WDES%NRPLWV , & 
               MPI_double_complex, W%WDES%COMM_INTER%MPI_COMM, ierr)
          IF ( ierr /= MPI_success ) CALL M_stop_ierr('ERROR: M_Allgather returns',ierror)

#ifdef gammareal
          IF (W%WDES%LOVERL) CALL MPI_Allgather(W1(NN)%CPROJ(1), W%WDES%NPROD, MPI_double_precision, &
               CPROJBUFF(1,1), W%WDES%NPROD, MPI_double_precision, W%WDES%COMM_INTER%MPI_COMM, ierr)
#else 
          IF (W%WDES%LOVERL) CALL MPI_Allgather(W1(NN)%CPROJ(1), W%WDES%NPROD, MPI_double_complex, &
               CPROJBUFF(1,1), W%WDES%NPROD, MPI_double_complex, W%WDES%COMM_INTER%MPI_COMM, ierr)
#endif
          IF ( ierr /= MPI_success ) CALL M_stop_ierr('ERROR: M_Allgather returns',ierror)

          ! store results back in correct storage position
          DO NN=1,W%WDES%NB_PAR
             W1((N-NB1)*NCPU+NN)%CW(:)   =CWBUFF(:,NN)
             W1((N-NB1)*NCPU+NN)%CR(:)   =CRBUFF(:,NN)
             W1((N-NB1)*NCPU+NN)%CPROJ(:)=CPROJBUFF(:,NN)
          ENDDO
       END IF
#endif
    ENDDO
  END SUBROUTINE W1_GATHER_NEW

!************************ SUBROUTINE W1_GATHER ************************
!
! This subroutine gathers a set of wavefunctions starting 
! from band NB1 until NB2 to all nodes
! compared to the previous routine only up to NLOC bands are collected
! 
!**********************************************************************

  SUBROUTINE W1_GATHER_N( W, NB1, NB2, ISP, W1, NLOC)
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed
    INTEGER :: NLOC          ! total number of bands to be collected

  ! local
    INTEGER :: NN, N, NCPU

    NCPU=W%WDES%NB_PAR

    DO N=NB1,NB2
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       IF (NN>NLOC) EXIT
       CALL W1_COPY( ELEMENT( W, W1(NN)%WDES1, N, ISP), W1(NN) )
       CALL FFTWAV_W1( W1(NN))
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO NN=1,NLOC
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CW(1), &
               SIZE(W1(NN)%CW),MOD(NN-1,NCPU)+1)
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(NN-1,NCPU)+1)
#ifndef gammareal
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
               W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
               W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#endif
       ENDDO

#ifdef MPI_barrier_after_bcast
      CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif

    ENDIF
#endif

  END SUBROUTINE W1_GATHER_N


!************************ SUBROUTINE W1_GATHER_GLB ********************
!
! This subroutine gathers a set of wavefunctions starting 
! from band NB1 until NB2 to all nodes
! compared to the previous version the global instead of local band 
! indices are supplied
! 
!**********************************************************************

  SUBROUTINE W1_GATHER_GLB( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed

  ! local
    INTEGER :: N_INTO_TOT, N, NCPU
    INTEGER :: ierror

    PROFILING_START('w1_gather_glb')

    NCPU=W%WDES%NB_PAR
    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(N_INTO_TOT-NB1+1)%WDES1, N, ISP), W1(N_INTO_TOT-NB1+1) )
          CALL FFTWAV_W1( W1(N_INTO_TOT-NB1+1))
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CW(1), &
               SIZE(W1(N)%CW),MOD(N_INTO_TOT-1,NCPU)+1)
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(N_INTO_TOT-1,NCPU)+1)
#ifndef gammareal
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
               W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
               W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#endif
       ENDDO

#ifdef MPI_barrier_after_bcast
      CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif

    ENDIF
#endif
    PROFILING_STOP('w1_gather_glb')

  END SUBROUTINE W1_GATHER_GLB

  SUBROUTINE W1_GATHER_GLB_( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed

  ! local
    INTEGER :: N_INTO_TOT, N, NB_LOCAL

    PROFILING_START('w1_gather_glb_')

    DO N=NB1,NB2
       N_INTO_TOT=N-NB1+1
       IF (MOD(N-1,W%WDES%NB_PAR)+1==W%WDES%NB_LOW) THEN
          NB_LOCAL=1+(N-1)/W%WDES%NB_PAR
          CALL W1_COPY(ELEMENT(W,W1(N_INTO_TOT)%WDES1,NB_LOCAL,ISP),W1(N_INTO_TOT))
          CALL FFTWAV_W1(W1(N_INTO_TOT))
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO N=NB1,NB2
          N_INTO_TOT=N-NB1+1
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N_INTO_TOT)%CW(1), &
               SIZE(W1(N_INTO_TOT)%CW),MOD(N-1,W%WDES%NB_PAR)+1)         
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N_INTO_TOT)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(N-1,W%WDES%NB_PAR)+1)
#ifndef gammareal
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N_INTO_TOT)%CPROJ(1), &
               W%WDES%NPROD,MOD(N-1,W%WDES%NB_PAR)+1)
#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(N_INTO_TOT)%CPROJ(1), &
               W%WDES%NPROD,MOD(N-1,W%WDES%NB_PAR)+1)
#endif
       ENDDO
#ifdef MPI_barrier_after_bcast
       CALLMPI( M_barrier(W%WDES%COMM_INTER))
#endif
    ENDIF
#endif
    PROFILING_STOP('w1_gather_glb_')

  END SUBROUTINE W1_GATHER_GLB_


!************************ SUBROUTINE W1_IGATHER_GLB *******************
!
! Gathers a set of wavefunctions starting from band NB1 until NB2
! (global band indices) to all nodes using non-blocking bcast_from.
!
! In case the code is compiled with -Dshmem_bcast_buffer communication
! will be between ranks within COMM_inter_node, i.e., only those ranks
! within COMM_INTER that are NOT on the same physical node will talk to
! eachother. MPI communication between ranks within COMM_INTER that
! reside on the same node is not necessary since they access a common
! shared memory segment.
!
! In case the code is NOT compiled with -Dshmem_bcast_buffer, the
! communication will be between all ranks in COMM_INTER.
! 
!**********************************************************************
  SUBROUTINE W1_IGATHER_GLB( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed

  ! local
    INTEGER :: NI, N, NB_LOCAL
    INTEGER :: ierror

    INTEGER :: nrequests,requests(3*(NB2-NB1+1))

    PROFILING_START('w1_igather_glb')

    DO N=NB1,NB2
       NI=N-NB1+1
       IF (MOD(N-1,W%WDES%NB_PAR)+1==W%WDES%NB_LOW) THEN
          NB_LOCAL=1+(N-1)/W%WDES%NB_PAR
          CALL W1_COPY(ELEMENT(W,W1(NI)%WDES1,NB_LOCAL,ISP),W1(NI))
          CALL FFTWAV_W1(W1(NI))
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%COMM_INTER%NCPU>1) THEN

       nrequests=0
       DO N=NB1,NB2
          NI=N-NB1+1

          nrequests=nrequests+1
          CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CW(1), &
         &     SIZE(W1(NI)%CW),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
          nrequests=nrequests+1
          CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CR(1), &
         &     SIZE(W1(NI)%CR),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
          IF (W%WDES%LOVERL) THEN
             nrequests=nrequests+1
#ifndef gammareal
             CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CPROJ(1), &
            &     SIZE(W1(NI)%CPROJ),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
#else
             CALL M_ibcast_d_from(W%WDES%COMM_INTER,W1(NI)%CPROJ(1), &
            &     SIZE(W1(NI)%CPROJ),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
#endif
          ENDIF
       ENDDO

       CALL M_waitall(nrequests,requests(1))

    ENDIF
#endif
    PROFILING_STOP('w1_igather_glb')

  END SUBROUTINE W1_IGATHER_GLB

!************************ SUBROUTINE W1_REDUCE_GLB ********************
!
! Take the sum of W1(i)%CW, W1(i)%CPROJ, for i=1,NB2-NB1+1 over all
! ranks in COMM_inter_node, using reduce_to.
! The result is reduced onto the MPI-rank in COMM_inter_node that would
! normally own (part of) the band with the global index N=NB1,...,NB2.
!
!**********************************************************************

  SUBROUTINE W1_REDUCE_GLB(WDES1,W1,NB1,NB2)
    IMPLICIT NONE
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1(:) ! array into which the merge is performed
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band

  ! local
    INTEGER :: N_INTO_TOT, N
    INTEGER :: ierror

    INTEGER :: IDO(NB2-NB1+1)

    PROFILING_START('w1_reduce_glb')
#ifdef MPI
    DO N=NB1,NB2
       N_INTO_TOT=N-NB1+1
       IF (MOD(N-1,WDES1%NB_PAR)+1==WDES1%NB_LOW) THEN
          IDO(N_INTO_TOT)=WDES1%COMM_inter_node%NODE_ME
       ELSE
          IDO(N_INTO_TOT)=0
       ENDIF
    ENDDO

    CALLMPI( M_sum_i(WDES1%COMM_inter_node,IDO,NB2-NB1+1))

    ! reduce W1
    IF (WDES1%DO_REDIS) THEN

       DO N=NB1,NB2
          N_INTO_TOT=N-NB1+1
          IF (IDO(N_INTO_TOT)>0) THEN
             CALL M_reduce_z_to(WDES1%COMM_inter_node,W1(N_INTO_TOT)%CW(1), &
            &     SIZE(W1(N_INTO_TOT)%CW),IDO(N_INTO_TOT))
#ifndef gammareal
             IF (WDES1%LOVERL) &
                CALL M_reduce_z_to(WDES1%COMM_inter_node,W1(N_INTO_TOT)%CPROJ(1), &
            &     SIZE(W1(N_INTO_TOT)%CPROJ),IDO(N_INTO_TOT))
#else
             IF (WDES1%LOVERL) &
                CALL M_reduce_d_to(WDES1%COMM_inter_node,W1(N_INTO_TOT)%CPROJ(1), &
            &     SIZE(W1(N_INTO_TOT)%CPROJ),IDO(N_INTO_TOT))
#endif
          ENDIF
       ENDDO

    ENDIF
#endif
    PROFILING_STOP('w1_reduce_glb')

  END SUBROUTINE W1_REDUCE_GLB


!************************ SUBROUTINE W1_IREDUCE_GLB *******************
!
! Take the sum of W1(i)%CW, W1(i)%CPROJ, for i=1,NB2-NB1+1 over all
! ranks in COMM_INTER, using non-blocking reduce_to.
! The result is reduced onto the MPI-rank in COMM_INTER that would
! normally own (part of) the band with the global index N=NB1,...,NB2.
!
!**********************************************************************

  SUBROUTINE W1_IREDUCE_GLB(WDES1,W1,NB1,NB2)
    IMPLICIT NONE
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1(:) ! array into which the merge is performed
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band

  ! local
    INTEGER :: NI, N
    INTEGER :: requests(2*(NB2-NB1+1))
    INTEGER :: nrequests

    PROFILING_START('w1_ireduce_glb')

    ! reduce W1
#ifdef MPI
    IF (WDES1%DO_REDIS) THEN

       nrequests=0
       DO N=NB1,NB2
          NI=N-NB1+1

          nrequests=nrequests+1
          CALL M_ireduce_z_to(WDES1%COMM_INTER,W1(NI)%CW(1), &
         &     SIZE(W1(NI)%CW),MOD(N-1,WDES1%NB_PAR)+1,requests(nrequests))

          IF (WDES1%LOVERL) THEN
             nrequests=nrequests+1
#ifndef gammareal
             CALL M_ireduce_z_to(WDES1%COMM_INTER,W1(NI)%CPROJ(1), &
            &     SIZE(W1(NI)%CPROJ),MOD(N-1,WDES1%NB_PAR)+1,requests(nrequests))
#else
             CALL M_ireduce_d_to(WDES1%COMM_INTER,W1(NI)%CPROJ(1), &
            &     SIZE(W1(NI)%CPROJ),MOD(N-1,WDES1%NB_PAR)+1,requests(nrequests))
#endif
          ENDIF
       ENDDO

       CALL M_waitall(nrequests,requests(1))

    ENDIF
#endif
    PROFILING_STOP('w1_ireduce_glb')

  END SUBROUTINE W1_IREDUCE_GLB

!************************ SUBROUTINE W1_GATHER_GLB_NOCR ***************
!
! This subroutine gathers a set of wavefunctions starting 
! from band NB1 until NB2 to all nodes
! compared to the previous version the global instead of local band 
! indices are supplied and the real space part is not
! communicated
! 
!**********************************************************************

  SUBROUTINE W1_GATHER_GLB_NOCR( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed

  ! local
    INTEGER :: N_INTO_TOT, N, NCPU

    PROFILING_START('w1_gather_glb_nocr')

    NCPU=W%WDES%NB_PAR
    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(N_INTO_TOT-NB1+1)%WDES1, N, ISP), W1(N_INTO_TOT-NB1+1) )
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CW(1), &
               SIZE(W1(N)%CW),MOD(N_INTO_TOT-1,NCPU)+1)
#ifndef gammareal
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
               W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
               W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#endif
       ENDDO
    ENDIF
#endif
    PROFILING_STOP('w1_gather_glb_nocr')

  END SUBROUTINE W1_GATHER_GLB_NOCR

#ifdef use_shmem
!************************ SUBROUTINE W1_GATHER_GLB_NOCR_SHMEM *********
!
! This subroutine gathers a set of wavefunctions starting 
! from band NB1 until NB2 to all nodes
! compared to the previous version the global instead of local band 
! indices are supplied
! 
!**********************************************************************

  SUBROUTINE W1_GATHER_GLB_NOCR_SHMEM(W,NB1,NB2,ISP,W1)
    IMPLICIT NONE
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1) :: W1(:) ! array into which the merge is performed
    ! local variables
    COMPLEX(q) :: CWORK(W1(1)%WDES1%NRPLWV) 
    INTEGER :: N_INTO_TOT,N,NCPU

    NCPU=W%WDES%NB_PAR
    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(N_INTO_TOT-NB1+1)%WDES1, N, ISP), W1(N_INTO_TOT-NB1+1) )
!         CALL FFTWAV_W1( W1(N_INTO_TOT-NB1+1))
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          ! the node that owns this band stores it onto CWORK
          IF (MOD(N_INTO_TOT-1,NCPU)+1==W%WDES%COMM_INTER%NODE_ME) CWORK=W1(N)%CW
          ! and broadcasts CWORK to all other nodes
          CALL M_bcast_z_from(W%WDES%COMM_INTER,CWORK(1), &
               SIZE(CWORK),MOD(N_INTO_TOT-1,NCPU)+1)
          ! but only the root nodes in COMM_SHMEM need to store CWORK back into W1
          IF (MOD(N_INTO_TOT-1,NCPU)+1/=W%WDES%COMM_INTER%NODE_ME &
             .AND.W%WDES%COMM_SHMEM%NODE_ME==1) W1(N)%CW=CWORK
#ifndef gammareal
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
               W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
               W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#endif
       ENDDO
    ENDIF
#endif

  END SUBROUTINE W1_GATHER_GLB_NOCR_SHMEM


!************************ SUBROUTINE W1_GATHER_GLB_ALLK_SHMEM *********
!
! Gathers a set of wavefunctions including band NB1 to NB2 over
! all k-points to all nodes using shared memory
!
!**********************************************************************

  SUBROUTINE W1_GATHER_GLB_ALLK_SHMEM( W, NB1, NB2, ISP, WF)
    IMPLICIT NONE
    TYPE (wavespin) W          ! wavefunction
    INTEGER :: NB1             ! starting band
    INTEGER :: NB2             ! final band
    INTEGER :: ISP             ! spin
    TYPE (wavefun1) :: WF(:,:) ! array into which the merge is performed
    ! local variables
    TYPE(wavefun1):: WAUX  
    TYPE(wavedes1), TARGET :: WDESAUX
    COMPLEX(q) :: CWORK(W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)
    GDEF CPROJ(W%WDES%NPROD)
    INTEGER :: N_INTO_TOT,N,NCPU,IK,NKPTS

    CALL SETWDES(W%WDES,WDESAUX,0)
    CALL NEWWAV(WAUX,WDESAUX,.TRUE.)
    NCPU=W%WDES%NB_PAR
    NKPTS=W%WDES%NKPTS

    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          DO IK=1,NKPTS
             ! set wave descriptor to current k-point for FFT
             CALL SETWDES(W%WDES,WDESAUX,IK)
             CALL W1_COPY_NOCR( ELEMENT( W, WDESAUX, N, ISP), WAUX)
             CALL FFTWAV_W1(WAUX)
             CALL ZCOPY( WDESAUX%GRID%MPLWV*WDESAUX%NRSPINORS, WAUX%CR(1), 1, WF(N_INTO_TOT-NB1+1,IK)%CR(1), 1)
             IF (W%WDES%LGAMMA) THEN
                CALL DCOPY( W%WDES%NPROD, WAUX%CPROJ(1), 1,  WF(N_INTO_TOT-NB1+1,IK)%CPROJ(1), 1)
             ELSE
                CALL ZCOPY( W%WDES%NPROD, WAUX%CPROJ(1), 1,  WF(N_INTO_TOT-NB1+1,IK)%CPROJ(1), 1)
             ENDIF
          END DO
       ENDIF
    ENDDO
    CALL DELWAV(WAUX,.TRUE.)
    ! distribute WF to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          DO IK=1,NKPTS
             ! the node that owns this band stores it onto CWORK
             IF (MOD(N_INTO_TOT-1,NCPU)+1==W%WDES%COMM_INTER%NODE_ME) THEN
                CWORK(1:SIZE(WF(N,IK)%CR))=WF(N,IK)%CR(:); CPROJ(1:SIZE(WF(N,IK)%CPROJ))=WF(N,IK)%CPROJ(:)
             ENDIF
             ! and broadcasts CWORK to all other nodes
             CALL M_bcast_z_from(W%WDES%COMM_INTER,CWORK(1),SIZE(CWORK),MOD(N_INTO_TOT-1,NCPU)+1)
             IF (W%WDES%LOVERL) THEN
#ifdef gammareal
                CALL M_bcast_d_from(W%WDES%COMM_INTER,CPROJ(1),SIZE(CPROJ),MOD(N_INTO_TOT-1,NCPU)+1)
#else
                CALL M_bcast_z_from(W%WDES%COMM_INTER,CPROJ(1),SIZE(CPROJ),MOD(N_INTO_TOT-1,NCPU)+1)         
#endif
             ENDIF
             ! but only the root nodes in COMM_SHMEM need to store CWORK back into WF
             IF (MOD(N_INTO_TOT-1,NCPU)+1/=W%WDES%COMM_INTER%NODE_ME.AND.W%WDES%COMM_SHMEM%NODE_ME==1) THEN
                WF(N,IK)%CR=CWORK(1:SIZE(WF(N,IK)%CR)); WF(N,IK)%CPROJ(:)=CPROJ(1:SIZE(WF(N,IK)%CPROJ))
             ENDIF
          ENDDO
       ENDDO
#ifdef MPI_barrier_after_bcast
       CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif
    ENDIF
#endif
  END SUBROUTINE W1_GATHER_GLB_ALLK_SHMEM
#endif

!************************ SUBROUTINE W1_GATHER_GLB_ALLK ***************
!
! Gathers a set of wavefunctions including band NB1 to NB2 over
! all k-points to all nodes
!
!**********************************************************************

  SUBROUTINE W1_GATHER_GLB_ALLK( W, NB1, NB2, ISP, WF)
    IMPLICIT NONE
    TYPE (wavespin) W          ! wavefunction
    INTEGER :: NB1             ! starting band
    INTEGER :: NB2             ! final band
    INTEGER :: ISP             ! spin
    TYPE (wavefun1) :: WF(:,:) ! array into which the merge is performed
    ! local
    TYPE (wavefun1):: WAUX  
    TYPE(wavedes1), TARGET :: WDESAUX
    INTEGER :: N_INTO_TOT,N,NCPU,IK,NKPTS


    CALL SETWDES(W%WDES,WDESAUX,0)
    CALL NEWWAV(WAUX,WDESAUX,.TRUE.)
    NCPU=W%WDES%NB_PAR
    NKPTS=W%WDES%NKPTS

    ! W1 contains bands NB1 to NB2 for all k-points
    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          DO IK=1,NKPTS
             ! set wave descriptor to current k-point for FFT
             CALL SETWDES(W%WDES,WDESAUX,IK)
             CALL W1_COPY_NOCR( ELEMENT( W, WDESAUX, N, ISP), WAUX)
             CALL FFTWAV_W1(WAUX)
             CALL ZCOPY( WDESAUX%GRID%MPLWV*WDESAUX%NRSPINORS, WAUX%CR(1), 1, WF(N_INTO_TOT-NB1+1,IK)%CR(1), 1)
             IF (W%WDES%LGAMMA) THEN
                CALL DCOPY( W%WDES%NPROD, WAUX%CPROJ(1), 1,  WF(N_INTO_TOT-NB1+1,IK)%CPROJ(1), 1)
             ELSE
                CALL ZCOPY( W%WDES%NPROD, WAUX%CPROJ(1), 1,  WF(N_INTO_TOT-NB1+1,IK)%CPROJ(1), 1)
             ENDIF
          END DO
       ENDIF
    ENDDO
    CALL DELWAV(WAUX,.TRUE.)
    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          DO IK=1,NKPTS
             CALL M_bcast_z_from(W%WDES%COMM_INTER,WF(N,IK)%CR(1), &
                  W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(N_INTO_TOT-1,NCPU)+1)
#ifndef gammareal
             IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,WF(N,IK)%CPROJ(1), &
                  W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#else
             IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,WF(N,IK)%CPROJ(1), &
                  W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#endif
          ENDDO
       END DO

#ifdef MPI_barrier_after_bcast
       CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif

    ENDIF
#endif
  END SUBROUTINE W1_GATHER_GLB_ALLK


!************************ SUBROUTINE W1_GATHER ************************
!
! This subroutine gathers a set of orbitals
! from band NB1 until NB2 to all nodes
! compared to previous versions the set in collected
! into work arrays 
! 
!**********************************************************************

  SUBROUTINE W1_GATHER_ARRAY( W, NB1, NB2, ISP, W1, CR, CPROJ)
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1     ! array for FFT 
    GDEF    :: CR(:,:)       ! collected real space orbitals
    GDEF    :: CPROJ(:,:)    ! collected projected orbitals

  ! local
    INTEGER :: NN, N, NLOC, NCPU

    IF (W%WDES%NRSPINORS/=1) THEN
       WRITE(*,*) 'internal error in W1_GATHER_ARRAY: at the moment spinors are not supported'
       STOP
    ENDIF

    NCPU=W%WDES%NB_PAR

    DO N=NB1,NB2
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       CALL W1_COPY( ELEMENT( W, W1%WDES1, N, ISP), W1 )
       CALL FFTWAV_W1( W1)

       IF (NN> SIZE(CR,2) .OR. NN >SIZE(CPROJ,2)) THEN
          WRITE(*,*) 'internal error in W1_GATHER_ARRAY: bound exceed ',NN, SIZE(CR,2),SIZE(CPROJ,2)
          STOP
       ENDIF
       IF (SIZE(CR,1)/=W1%WDES1%GRID%NPLWV) THEN
          WRITE(*,*) 'internal error in W1_GATHER_ARRAY: size mismatch ',SIZE(CR,1),W1%WDES1%GRID%NPLWV
       ENDIF

       CR(1:W1%WDES1%GRID%RL%NP, NN)=W1%CR(1:W1%WDES1%GRID%RL%NP)
       ! pad CR with zeros (just in case)
       CR(W1%WDES1%GRID%RL%NP+1:SIZE(CR,1), NN)=0
       CPROJ(:,NN)=W1%CPROJ(:)

    ! distribute W1 to all nodes
#ifdef MPI
       IF (W%WDES%DO_REDIS) THEN
        DO NN=1,NCPU
#ifndef gammareal
          CALL M_bcast_z_from(W%WDES%COMM_INTER,CR(1,NN+(N-NB1)*NCPU), &
               SIZE(CR,1),MOD(NN-1,NCPU)+1)
#else
          CALL M_bcast_d_from(W%WDES%COMM_INTER,CR(1,NN+(N-NB1)*NCPU), &
               SIZE(CR,1),MOD(NN-1,NCPU)+1)
#endif

#ifndef gammareal
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,CPROJ(1,NN+(N-NB1)*NCPU), &
               W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,CPROJ(1,NN+(N-NB1)*NCPU), &
               W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#endif
         ENDDO
       END IF
#endif
      ENDDO

  END SUBROUTINE W1_GATHER_ARRAY


!************************ SUBROUTINE W1_GATHER ************************
!
! This subroutine gathers a set of orbitals
! from band NB1 until NB2 to all nodes
! this version collects the place wave coefficients and the
! CPROJ coefficients
! 
!**********************************************************************

  SUBROUTINE W1_GATHER_ARRAY_RECIPROCAL( W, NB1, NB2, ISP, W1, CG, CPROJ)
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1     ! array for FFT 
    COMPLEX(q) :: CG(:,:)    ! collected real space orbitals
    GDEF       :: CPROJ(:,:) ! collected projected orbitals

    COMPLEX(q) :: CWBUFF(SIZE(CG,1))        ! send buffer
    GDEF       :: CPROJBUFF(SIZE(CPROJ,1))  ! send buffer for projectors

  ! local
    INTEGER :: NN, N, NLOC, NCPU
    INTEGER :: ierr          ! error for mpi

    IF (W%WDES%NRSPINORS/=1) THEN
       WRITE(*,*) 'internal error in W1_GATHER_ARRAY: at the moment spinors are not supported'
       STOP
    ENDIF
    IF (SIZE(CPROJ,1) > SIZE(W1%CPROJ,1)) THEN
       WRITE (*,*) 'internal error in W1_GATHER_ARRAY_RECIPROCAL: CPROJ size inconsistent ', SIZE(CPROJ,1), SIZE(W1%CPROJ,1)
       STOP
    ENDIF

    NCPU=W%WDES%NB_PAR

    ! NB1 and NB2 are local indices for collecting the bands
    DO N=NB1,NB2
       !NB_LOW is the off-set of the local node
       !NN is index into CG array
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       !this copies orbital N (reciprocal and real (if allocated) part to W1)
       !these are local indices
       CALL W1_COPY( ELEMENT( W, W1%WDES1, N, ISP), W1 )

       IF (NN> SIZE(CG,2) .OR. NN >SIZE(CPROJ,2)) THEN
          WRITE(*,*) 'internal error in W1_GATHER_ARRAY_RECIPROCAL: bound exceed ',NN, SIZE(CG,2),SIZE(CPROJ,2)
          STOP
       ENDIF
  
       ! copy data over to return arrays: CG and CPROJ
       CG(1:W1%WDES1%NPL ,NN)=W1%CW(1:W1%WDES1%NPL)
       CG(W1%WDES1%NPL+1:,NN)=0.0_q                 ! pad with zero

       CPROJ(:,NN)=W1%CPROJ(1:SIZE(CPROJ,1))
      ! distribute W1 to all nodes, CG(1,(N-NB)*NCPU)
#ifdef MPI
       ! MPI_IN_PLACE was tested but found to be much slower (code has been removed)
       IF (W%WDES%DO_REDIS) THEN
          CPROJBUFF(:)=CPROJ(:,NN)
          CWBUFF(:)   =CG(:,NN)
          CALL MPI_Allgather(CWBUFF(1), SIZE(CWBUFF), MPI_double_complex, CG(1,(N-NB1)*NCPU+1), SIZE(CG,1), & 
               MPI_double_complex, W%WDES%COMM_INTER%MPI_COMM, ierr)
          IF ( ierr /= MPI_success ) CALL M_stop_ierr('ERROR: M_Allgather returns',ierror)
#ifdef gammareal
          IF (W%WDES%LOVERL) CALL MPI_Allgather(CPROJBUFF(1), SIZE(CPROJBUFF), MPI_double_precision, &
               CPROJ(1,(N-NB1)*NCPU+1), SIZE(CPROJ,1), MPI_double_precision, W%WDES%COMM_INTER%MPI_COMM, ierr)
#else 
          IF (W%WDES%LOVERL) CALL MPI_Allgather(CPROJBUFF(1), SIZE(CPROJBUFF), MPI_double_complex,&
               CPROJ(1,(N-NB1)*NCPU+1), SIZE(CPROJ,1), MPI_double_complex, W%WDES%COMM_INTER%MPI_COMM, ierr)
#endif
          IF ( ierr /= MPI_success ) CALL M_stop_ierr('ERROR: M_Allgather returns',ierror)
       END IF
#endif
    ENDDO

  END SUBROUTINE W1_GATHER_ARRAY_RECIPROCAL

!
! old version, uses bcast and is much slower
!
  SUBROUTINE W1_GATHER_ARRAY_RECIPROCAL_OLD( W, NB1, NB2, ISP, W1, CG, CPROJ)
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1     ! array for FFT 
    COMPLEX(q) :: CG(:,:)    ! collected real space orbitals
    GDEF       :: CPROJ(:,:) ! collected projected orbitals

  ! local
    INTEGER :: NN, N, NLOC, NCPU

    IF (W%WDES%NRSPINORS/=1) THEN
       WRITE(*,*) 'internal error in W1_GATHER_ARRAY: at the moment spinors are not supported'
       STOP
    ENDIF

    NCPU=W%WDES%NB_PAR

    DO N=NB1,NB2
       !NB_LOW will be the off-set of the local node
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       !this copies orbital N (reciprocal and real (if allocated) part to W1)
       !these are local indices
       CALL W1_COPY( ELEMENT( W, W1%WDES1, N, ISP), W1 )

       IF (NN> SIZE(CG,2) .OR. NN >SIZE(CPROJ,2)) THEN
          WRITE(*,*) 'internal error in W1_GATHER_ARRAY: bound exceed ',NN, SIZE(CG,2),SIZE(CPROJ,2)
          STOP
       ENDIF
  
       !copy the local data to the position with global index in the strip
       CG(1:W1%WDES1%NPL, NN)=W1%CW(1:W1%WDES1%NPL)
       CG(W1%WDES1%NPL+1: SIZE(CG,1), NN)=0
       !the first index CPROJ size can be smaller than the second one (less cpus/tau than total number)
       !but this copies only the data for valid indices in the first array (or?)
       !!!test
       !IF (SIZE(CPROJ,1).LT.SIZE(W1%CPROJ,1)) WRITE (*,*) NN, W1%CPROJ(SIZE(CPROJ,1)+1:)
       IF (SIZE(CPROJ,1).GT.SIZE(W1%CPROJ,1)) THEN
         WRITE (*,*) 'internal error in W1_GATHER_ARRAY: Sizes inconsistent ', SIZE(CPROJ,1), SIZE(W1%CPROJ,1)
       ENDIF
       CPROJ(:,NN)=W1%CPROJ(1:SIZE(CPROJ,1))
       !IF (N.eq.NB1) WRITE (*,*) 'in gather array', SIZE(CPROJ,1), SIZE(CPROJ,2), SIZE(W1%CPROJ,1)
       !!!test

    ! distribute W1 to all nodes
#ifdef MPI
       IF (W%WDES%DO_REDIS) THEN
       DO NN=1,NCPU
          !so broadcast the data from node holding index NN to all other nodes 
          !I wonder how efficient this is, might be all right 10000 coefficients is 160kB of data, 
          !a bit on the low side but fine
          CALL M_bcast_z_from(W%WDES%COMM_INTER,CG(1,NN+(N-NB1)*NCPU), &
               SIZE(CG,1),MOD(NN-1,NCPU)+1)
#ifndef gammareal
        ! OK, the point is that the W%WDES%NPROD array can be longer than the size of CPROJ
        ! since both are a multiple of NCPU but once the total number and in the second case per tau
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,CPROJ(1,NN+(N-NB1)*NCPU), &
               SIZE(CPROJ,1),MOD(NN-1,NCPU)+1)

#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,CPROJ(1,NN+(N-NB1)*NCPU), &
               SIZE(CPROJ,1),MOD(NN-1,NCPU)+1)
#endif
       ENDDO
       END IF
#endif
    ENDDO

  END SUBROUTINE W1_GATHER_ARRAY_RECIPROCAL_OLD


!************************ SUBROUTINE W1_GATHER_DISTR ******************
!
! This subroutine gathers a set of wavefunctions starting 
! from band NB1 until NB2 to all nodes
! compared to the previous version the global instead of local band 
! indices are supplied
! 
!**********************************************************************

  SUBROUTINE W1_GATHER_DISTR( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed
    ! local variables
    INTEGER :: N_INTO_TOT, N, NCPU
    INTEGER :: NBMIN, NBMAX
    TYPE (wavefun1):: WTMP
    
    NCPU=W%WDES%NB_PAR
    
    CALL NEWWAV(WTMP, W1(1)%WDES1, .FALSE.)
        
    ! establish global band index interval
    NBMIN=NB1
    NBMAX=NB2
#ifdef MPI
    DO N=1,NCPU
       CALL M_bcast_i_from(W%WDES%COMM_INTER, NBMIN, 1, n)
       CALL M_bcast_i_from(W%WDES%COMM_INTER, NBMAX, 1, n)
       IF (NB1<NBMIN) NBMIN=NB1
       IF (NB2>NBMAX) NBMAX=NB2
    ENDDO
#endif
    
    DO N_INTO_TOT=NBMIN,NBMAX
    
       ! if band resides on this node copy it to WTMP
       IF (MOD(N_INTO_TOT-W%WDES%NB_LOW,NCPU)==0) THEN
          N=(N_INTO_TOT-W%WDES%NB_LOW)/NCPU+1
          CALL W1_COPY( ELEMENT( W, WTMP%WDES1, N, ISP), WTMP )
       ENDIF

#ifdef MPI
       ! broadcast WTMP from the node where it resides
       CALL M_bcast_z_from(W%WDES%COMM_INTER,WTMP%CW(1), &
            SIZE(WTMP%CW),MOD(N_INTO_TOT-1,NCPU)+1)
#ifndef gammareal
       IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,WTMP%CPROJ(1), &
            W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#else
       IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,WTMP%CPROJ(1), &
            W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#endif
#endif
       ! if this band is targeted to reside on this node
       ! copy WTMP to W1(N_INTO_TOT-NB1+1)
       IF (N_INTO_TOT>=NB1 .AND. N_INTO_TOT<=NB2) THEN
          N=N_INTO_TOT-NB1+1
          CALL W1_COPY(WTMP,W1(N))
       ENDIF
    ENDDO
    
    ! FFT to real space
    DO N=1,(NB2-NB1)+1
       CALL FFTWAV_W1(W1(N))
    ENDDO

    CALL DELWAV(WTMP, .FALSE.)

  END SUBROUTINE W1_GATHER_DISTR


!************************ SUBROUTINE W1_GATHER_KSEL *******************
!
! This subroutine gathers a set of wavefunctions starting 
! from band NB1 until NB2 and distributes the data over k in a round
! robin fashion
! the data distribution is based on an index (k-point index) supplied as
! the last argument
!
!**********************************************************************

  SUBROUTINE W1_GATHER_KSEL( W, NB1, NB2, ISP, W1, NK)
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed
    INTEGER :: NK            ! k-points intex

   ! local
    TYPE (wavefun1):: WTMP   ! temporary for FFT
    INTEGER :: NN, N, NLOC, NCPU

    NCPU=W%WDES%NB_PAR

    CALL NEWWAV(WTMP, W1(1)%WDES1, .TRUE.)

    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(1)%WDES1, N, ISP), WTMP )
          CALL FFTWAV_W1( WTMP)
       ENDIF
#ifdef MPI
       IF (MOD(NK-1,NCPU)+1 ==  W%WDES%NB_LOW) THEN
          ! receive from all other nodes and local copy
          DO NN=1,NCPU
             N_INTO_TOT=(N-1)*NCPU+NN
             IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
                IF (NN==W%WDES%NB_LOW) THEN
                   CALL W1_COPY( WTMP, W1(N_INTO_TOT-NB1+1) )
                ELSE
!                   WRITE(*,*) W%WDES%NB_LOW, 'receive', N_INTO_TOT-NB1+1, 'from' , NN
                   CALL M_recv_z(W%WDES%COMM_INTER, NN , & 
                        W1(N_INTO_TOT-NB1+1)%CR(1), W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)
#ifndef gammareal
                   IF (W%WDES%LOVERL) CALL M_recv_z(W%WDES%COMM_INTER, NN, & 
                        W1(N_INTO_TOT-NB1+1)%CPROJ(1), W%WDES%NPROD)
#else
                   IF (W%WDES%LOVERL) CALL M_recv_d(W%WDES%COMM_INTER, NN, & 
                        W1(N_INTO_TOT-NB1+1)%CPROJ(1), W%WDES%NPROD)
#endif
                ENDIF
             ENDIF
          END DO
       ELSE
          IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
!             WRITE(*,*) W%WDES%NB_LOW,'send', N_INTO_TOT, 'to' , MOD(NK-1,NCPU)+1
             CALL M_send_z(W%WDES%COMM_INTER, MOD(NK-1,NCPU)+1, & 
                  WTMP%CR(1), W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)
#ifndef gammareal
             IF (W%WDES%LOVERL) CALL M_send_z(W%WDES%COMM_INTER, MOD(NK-1,NCPU)+1, & 
                  WTMP%CPROJ(1), W%WDES%NPROD)
#else
             IF (W%WDES%LOVERL) CALL M_send_d(W%WDES%COMM_INTER, MOD(NK-1,NCPU)+1, & 
                  WTMP%CPROJ(1), W%WDES%NPROD)
#endif
          ENDIF
       ENDIF
#endif
    ENDDO

    CALL DELWAV(WTMP, .TRUE.)

  END SUBROUTINE W1_GATHER_KSEL


!************************ SUBROUTINE W1_GATHER_KNODESEL ***************
!
! This subroutine is similar to W1_GATHER_KSEL, but at variance with
! it, one explicitly select the node on which you gather the wave
! fucntion (instead of distributing k-points in round robin fashion).
!
!**********************************************************************

  SUBROUTINE W1_GATHER_KNODESEL( W, NB1, NB2, ISP, W1, NODE)
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed
    INTEGER :: NODE          ! which node receives

   ! local
    TYPE (wavefun1):: WTMP   ! temporary for FFT
    INTEGER :: NN, N, NLOC, NCPU

    NCPU=W%WDES%NB_PAR

    CALL NEWWAV(WTMP, W1(1)%WDES1, .TRUE.)

    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(1)%WDES1, N, ISP), WTMP )
          CALL FFTWAV_W1( WTMP)
       ENDIF
#ifdef MPI
       IF (NODE==W%WDES%NB_LOW) THEN
          ! receive from all other nodes and local copy
          DO NN=1,NCPU
             N_INTO_TOT=(N-1)*NCPU+NN
             IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
                IF (NN==W%WDES%NB_LOW) THEN
                   CALL W1_COPY( WTMP, W1(N_INTO_TOT-NB1+1) )
                ELSE
                   CALL M_recv_z(W%WDES%COMM_INTER, NN , & 
                        W1(N_INTO_TOT-NB1+1)%CR(1), W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)
#ifndef gammareal
                   IF (W%WDES%LOVERL) CALL M_recv_z(W%WDES%COMM_INTER, NN, & 
                        W1(N_INTO_TOT-NB1+1)%CPROJ(1), W%WDES%NPROD)
#else
                   IF (W%WDES%LOVERL) CALL M_recv_d(W%WDES%COMM_INTER, NN, & 
                        W1(N_INTO_TOT-NB1+1)%CPROJ(1), W%WDES%NPROD)
#endif
                ENDIF
             ENDIF
          END DO
       ELSE
          IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
             CALL M_send_z(W%WDES%COMM_INTER,NODE, & 
                  WTMP%CR(1), W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)
#ifndef gammareal
             IF (W%WDES%LOVERL) CALL M_send_z(W%WDES%COMM_INTER,NODE, & 
                  WTMP%CPROJ(1), W%WDES%NPROD)
#else
             IF (W%WDES%LOVERL) CALL M_send_d(W%WDES%COMM_INTER,NODE, & 
                  WTMP%CPROJ(1), W%WDES%NPROD)
#endif
          ENDIF
       ENDIF
#endif
    ENDDO
    CALL DELWAV(WTMP, .TRUE.)
  END SUBROUTINE W1_GATHER_KNODESEL


!************************ SUBROUTINE W1_GATHER_W1 *********************
!
! This subroutine gathers from a wavefunction array
! instead of W
! 
!**********************************************************************

  SUBROUTINE W1_GATHER_W1( W, NB2, W1_ORIG, W1)
    TYPE (wavespin):: W       ! wavefunction
    TYPE (wavefun1):: W1_ORIG(:) ! wavefunction array to be merged
    TYPE (wavefun1):: W1(:)   ! array into which the merge is performed
    INTEGER :: NB2            ! final band

  ! local
    INTEGER :: NN, N, NLOC, NCPU

    NCPU=W%WDES%NB_PAR
    IF (SIZE(W1_ORIG)< NB2) THEN
       WRITE(*,*) 'internal error in W1_GATHER_W1: W1_ORIG is not sufficiently large'
       STOP
    ENDIF
       
    DO N=1,NB2
       NN=(N-1)*NCPU+W%WDES%NB_LOW
       CALL W1_COPY_NOCR( W1_ORIG(N), W1(NN) )
       CALL FFTWAV_W1( W1(NN))
    ENDDO

    NLOC=NB2*W%WDES%NB_PAR

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO NN=1,NLOC
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(NN-1,NCPU)+1)
#ifndef gammareal
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
               W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
               W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#endif
       ENDDO

#ifdef MPI_barrier_after_bcast
      CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif

    ENDIF
#endif

  END SUBROUTINE W1_GATHER_W1


!************************ SUBROUTINE W1_GATHER_STRIP ******************
!
! This subroutine gathers from a wavefunction array instead of W.
! Collected are the bands W1_ORIG(N): NB1<=NB_LOW+(N-1)*NB_PAR<=NB2,
! into W1(1:NB2-NB1+1)
! 
!**********************************************************************

  SUBROUTINE W1_GATHER_STRIP( W, NB1, NB2, W1_ORIG, W1)
    TYPE (wavespin):: W       ! wavefunction
    TYPE (wavefun1):: W1_ORIG(:) ! wavefunction array to be merged
    TYPE (wavefun1):: W1(:)   ! array into which the merge is performed
    INTEGER :: NB1, NB2       ! final band

  ! local
    INTEGER :: NN, N, NCPU

    PROFILING_START('w1_gather_strip')

    NCPU=W%WDES%NB_PAR
    IF (SIZE(W1_ORIG)*NCPU<NB2) THEN
       WRITE(*,*) 'internal error in W1_GATHER_STRIP: W1_ORIG is not sufficiently large'
       STOP
    ENDIF

    IF (SIZE(W1)<NB2-NB1+1) THEN
       WRITE(*,*) 'internal error in W1_GATHER_STRIP: W1 is not sufficiently large'
       STOP
    ENDIF
       
    DO N=1,SIZE(W1_ORIG)
       NN=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NN<NB1.OR.NN>NB2) CYCLE
       CALL W1_COPY_NOCR( W1_ORIG(N), W1(NN-NB1+1) )
       CALL FFTWAV_W1( W1(NN-NB1+1))
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO NN=NB1,NB2
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN-NB1+1)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(NN-1,NCPU)+1)
#ifndef gammareal
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN-NB1+1)%CPROJ(1), &
               W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(NN-NB1+1)%CPROJ(1), &
               W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#endif
       ENDDO

#ifdef MPI_barrier_after_bcast
      CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif

    ENDIF
#endif
    PROFILING_STOP('w1_gather_strip')

  END SUBROUTINE W1_GATHER_STRIP


!************************ SUBROUTINE W1_IGATHER_STRIP *****************
!
! Gathers a set of wavefunctions starting from band NB1 until NB2
! (global band indices) to all nodes using non-blocking bcast_from.
! This subroutine gathers from a wavefunction array W1_ORIG instead
! of W.
!
! In case the code is compiled with -Dshmem_bcast_buffer communication
! will be between ranks within COMM_inter_node, i.e., only those ranks
! within COMM_INTER that are NOT on the same physical node will talk to
! eachother. MPI communication between ranks within COMM_INTER that
! reside on the same node is not necessary since they access a common
! shared memory segment.
!
! In case the code is NOT compiled with -Dshmem_bcast_buffer, the
! communication will be between all ranks in COMM_INTER.
! 
!**********************************************************************
  SUBROUTINE W1_IGATHER_STRIP( W, NB1, NB2, W1_ORIG, W1)
    TYPE (wavespin):: W          ! wavefunction
    TYPE (wavefun1):: W1_ORIG(:) ! wavefunction array to be merged
    TYPE (wavefun1):: W1(:)      ! array into which the merge is performed
    INTEGER :: NB1,NB2           ! first and final band

  ! local
    INTEGER :: NI,N,NB_LOCAL
    INTEGER :: ierror

    INTEGER :: nrequests,requests(3*(NB2-NB1+1))


    PROFILING_START('w1_igather_strip')

    IF (SIZE(W1_ORIG)*W%WDES%NB_PAR<NB2) THEN
       WRITE(*,*) 'internal error in W1_IGATHER_STRIP: W1_ORIG is not sufficiently large'
       STOP
    ENDIF

    IF (SIZE(W1)<NB2-NB1+1) THEN
       WRITE(*,*) 'internal error in W1_IGATHER_STRIP: W1 is not sufficiently large'
       STOP
    ENDIF
       
    DO N=NB1,NB2
       NI=N-NB1+1
       IF (MOD(N-1,W%WDES%NB_PAR)+1==W%WDES%NB_LOW) THEN
          NB_LOCAL=1+(N-1)/W%WDES%NB_PAR
          CALL W1_COPY_NOCR(W1_ORIG(NB_LOCAL),W1(NI))
          CALL FFTWAV_W1(W1(NI))
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%COMM_INTER%NCPU>1) THEN

       nrequests=0
       DO N=NB1,NB2
          NI=N-NB1+1

          nrequests=nrequests+1
          CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CW(1), &
         &     SIZE(W1(NI)%CW),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
          nrequests=nrequests+1
          CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CR(1), &
         &     SIZE(W1(NI)%CR),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
          IF (W%WDES%LOVERL) THEN
             nrequests=nrequests+1
#ifndef gammareal
             CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CPROJ(1), &
            &     SIZE(W1(NI)%CPROJ),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
#else
             CALL M_ibcast_d_from(W%WDES%COMM_INTER,W1(NI)%CPROJ(1), &
            &     SIZE(W1(NI)%CPROJ),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
#endif
          ENDIF
       ENDDO

       CALL M_waitall(nrequests,requests(1))

    ENDIF
#endif
    PROFILING_STOP('w1_igather_strip')

  END SUBROUTINE W1_IGATHER_STRIP

!************************* SUBROUTINE WVREAL_PRECISE *******************
!
! this subroutine forces the wavefunction to be real at the Gamma-point
! it is required for the gamma point only mode
! to avoid that small non real components develop
! this version is exact and works through an FFT to real space
! and the forcing the wavefunction to be real
! the routine is required if only subspace rotations are performed
! since the subspace rotation routine can not force wavefunctions
! to become real
!
! if LORBITALREAL is .TRUE., the routine also tries to make the
! orbitals real (in real space) by calculating the orbital
!   phi(r) = u_k(r) e^ikr
! taking the real part, storing the coefficients back
! this is possible at special k-points such as points at the BZ boundary
! e.g. k=(+-0.5,+-0.5,+-0.5)
!      k=(+-0.5,+-0.5,  0)
!      k=(+-0.5,  0  ,  0)    and permutations thereof
! the routine will not work properly at any other k-points
! (and is in fact by-passed in this case)
! ideally no subspace rotation should be called after calling this
! routine (though orthogonalization is fine)
!
!***********************************************************************

  SUBROUTINE WVREAL_PRECISE(W)
    USE constant
    IMPLICIT NONE
    TYPE (wavespin) W
  ! local
    INTEGER :: NK, ISP, NB, ISPINOR, N
    TYPE (wavedes1)    WDES1          ! descriptor for one k-point
    TYPE (wavefun1)    W1             ! current wavefunction
    LOGICAL :: LPHASE
    COMPLEX(q) :: CPHASE(W%WDES%GRID%MPLWV)
    COMPLEX(q) :: AVERAGE_PHASE
    REAL(q), EXTERNAL :: RANE

#ifdef gammareal
    CALL SETWDES(W%WDES,WDES1,0)
    CALL NEWWAV(W1, WDES1, .TRUE.)
    !
    ! search for position of G=0 zero element
    !
    DO NK  =1,W%WDES%NKPTS
#ifdef MPI
       IF (MOD(NK-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
       IF (ABS(SUM(W%WDES%VKPT(:,NK)*W%WDES%VKPT(:,NK)))< 1E-12) THEN
          CALL SETWDES(W%WDES,WDES1,NK)
          DO ISP=1,W%WDES%ISPIN
             DO NB=1,W%WDES%NBANDS
                CALL W1_COPY(ELEMENT(W, WDES1, NB, ISP), W1)
                CALL FFTWAV_W1(W1)
                DO ISPINOR =0,WDES1%NRSPINORS-1
                   DO N=1,WDES1%GRID%RL%NP
                      W1%CR(N+ISPINOR*WDES1%GRID%MPLWV)=REAL(W1%CR(N+ISPINOR*WDES1%GRID%MPLWV),q) & 
                          *(1.0_q/WDES1%GRID%NPLWV)
                   ENDDO
                   CALL FFTEXT(WDES1%NGVECTOR, WDES1%NINDPW(1), W1%CR(1+ISPINOR*WDES1%GRID%MPLWV),W%CW(1+ISPINOR*WDES1%NGVECTOR,NB,NK,ISP),WDES1%GRID,.FALSE.)
                ENDDO
             ENDDO
          ENDDO
       END IF
    ENDDO

    CALL DELWAV(W1, .TRUE.)
#else
    IF (W%WDES%LORBITALREAL) THEN
    CALL SETWDES(W%WDES,WDES1,0)
    CALL NEWWAV(W1, WDES1, .TRUE.)
!
! force orbitals to be real for non-gamma point only method
!
    DO NK  =1,W%WDES%NKPTS
#ifdef MPI
       IF (MOD(NK-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
       IF (ABS(MOD(W%WDES%VKPT(1,NK)*2+100,1._q))< 1E-6 .AND. &
           ABS(MOD(W%WDES%VKPT(2,NK)*2+100,1._q))< 1E-6 .AND. &
           ABS(MOD(W%WDES%VKPT(3,NK)*2+100,1._q))< 1E-6) THEN
       CALL SETPHASE_WVREAL(W%WDES%VKPT(:,NK), W%WDES%GRID, CPHASE, LPHASE)
       
       CALL SETWDES(W%WDES,WDES1,NK)
       DO ISP=1,W%WDES%ISPIN
          DO NB=1,W%WDES%NBANDS
             CALL W1_COPY(ELEMENT(W, WDES1, NB, ISP), W1)
             CALL FFTWAV_W1(W1)
             DO ISPINOR =0,WDES1%NRSPINORS-1
                ! chose a random complex direction and project onto this
                AVERAGE_PHASE=EXP(CITPI*RANE())
                CALLMPI (M_bcast_z( WDES1%COMM_INB, AVERAGE_PHASE, 1))
                DO N=1,WDES1%GRID%RL%NP
                   IF (LPHASE) THEN
                      W1%CR(N+ISPINOR*WDES1%GRID%MPLWV)=CONJG(CPHASE(N))*REAL(W1%CR(N+ISPINOR*WDES1%GRID%MPLWV)*CPHASE(N)*AVERAGE_PHASE,q) & 
                      *(1.0_q/WDES1%GRID%NPLWV)
                   ELSE
                      W1%CR(N+ISPINOR*WDES1%GRID%MPLWV)=REAL(W1%CR(N+ISPINOR*WDES1%GRID%MPLWV)*AVERAGE_PHASE,q) & 
                      *(1.0_q/WDES1%GRID%NPLWV)
                   ENDIF
                ENDDO
                CALL FFTEXT(WDES1%NGVECTOR, WDES1%NINDPW(1), W1%CR(1+ISPINOR*WDES1%GRID%MPLWV),W%CW(1+ISPINOR*WDES1%NGVECTOR,NB,NK,ISP),WDES1%GRID,.FALSE.)
             ENDDO
          ENDDO
       ENDDO
       ENDIF
    ENDDO

    CALL DELWAV(W1, .TRUE.)
    ENDIF
#endif

    RETURN
  END SUBROUTINE WVREAL_PRECISE


  SUBROUTINE SETPHASE_WVREAL(VKPT, GRID, CPHASE, LPHASE)
    USE constant

    IMPLICIT NONE

    REAL(q) :: VKPT(3)
    TYPE (grid_3d) GRID
    COMPLEX(q) :: CPHASE(GRID%MPLWV)
    LOGICAL LPHASE
  ! local
    REAL(q),PARAMETER :: TINY=1E-6_q
    REAL(q) F1, F2, F3
    INTEGER NC, N, IND
    COMPLEX(q) C, CD, CSUM
    CSUM=0
#ifdef gammareal
    LPHASE=.FALSE.
    RETURN
#endif
    
    IF (ABS(VKPT(1))>TINY .OR. ABS(VKPT(2))>TINY .OR. ABS(VKPT(3))>TINY) THEN
       LPHASE=.TRUE.
       F1=TPI/GRID%NGX*VKPT(1)
       F2=TPI/GRID%NGY*VKPT(2)
       F3=TPI/GRID%NGZ*VKPT(3)

       IF (GRID%RL%NFAST==3) THEN
          CD=EXP(CMPLX(0,F3,q))
          IND=0
          DO NC=1,GRID%RL%NCOL
             C=EXP(CMPLX(0,F1*(GRID%RL%I2(NC)-1)+F2*(GRID%RL%I3(NC)-1),q))
             DO N=1,GRID%RL%NROW
                IND=IND+1
                CPHASE(IND)=C
                C=C*CD
             ENDDO
          ENDDO
       ELSE
          CD=EXP(CMPLX(0,F1,q))
          IND=0
          DO NC=1,GRID%RL%NCOL
             C=EXP(CMPLX(0,F2*(GRID%RL%I2(NC)-1)+F3*(GRID%RL%I3(NC)-1),q))
             DO N=1,GRID%RL%NROW
                IND=IND+1
                CPHASE(IND)=C
                CSUM=CSUM+C
                C=C*CD
             ENDDO
          ENDDO
       ENDIF
       LPHASE=.TRUE.
    ELSE
       LPHASE=.FALSE.
    ENDIF
  END SUBROUTINE SETPHASE_WVREAL
      
END MODULE wave_high


!***********************************************************************
!
! assign a a W1 structure an initial value from a wavefunction 
! array
! the plane wave coefficients and the wave function character (optional)
! must be supplied
! for performance reason no runtime checking on anything is performed
!      
!***********************************************************************

  SUBROUTINE ARRAY_TO_W1( W1, C, CPROJ)
    USE prec
    USE wave
    IMPLICIT NONE
    TYPE (wavefun1)    W1
    COMPLEX(q):: C(*)
    GDEF, OPTIONAL :: CPROJ(*)

#if defined (use_fastbcopy)
    CALL FASTBCOPY8(W1%WDES1%NRPLWV*2 ,C, W1%CW(1))
#else
    CALL ZCOPY( W1%WDES1%NRPLWV, C, 1, W1%CW(1), 1)
#endif
    IF (PRESENT(CPROJ)) THEN
       IF (W1%WDES1%LGAMMA) THEN
          CALL DCOPY( W1%WDES1%NPROD, CPROJ(1), 1, W1%CPROJ(1), 1)
       ELSE
          CALL ZCOPY( W1%WDES1%NPROD, CPROJ(1), 1, W1%CPROJ(1), 1)
       ENDIF
    END IF


  END SUBROUTINE ARRAY_TO_W1


!***********************************************************************
!
! assign a W1 structure to an array (reverse of the previous operation)
!      
!***********************************************************************

  SUBROUTINE W1_TO_ARRAY( W1, C,  CPROJ)
    USE prec
    USE wave
    IMPLICIT NONE
    TYPE (wavefun1)    W1
    COMPLEX(q) :: C(*)
    GDEF, OPTIONAL :: CPROJ(*)

#if defined (use_fastbcopy)
    CALL FASTBCOPY8(W1%WDES1%NRPLWV*2 , W1%CW(1), C)
#else
    CALL ZCOPY( W1%WDES1%NRPLWV, W1%CW(1), 1, C, 1)
#endif
    IF (PRESENT(CPROJ)) THEN
       IF (W1%WDES1%LGAMMA) THEN
          CALL DCOPY( W1%WDES1%NPROD,  W1%CPROJ(1), 1, CPROJ(1), 1)
       ELSE
          CALL ZCOPY( W1%WDES1%NPROD,  W1%CPROJ(1), 1, CPROJ(1), 1)
       ENDIF
    END IF

  END SUBROUTINE W1_TO_ARRAY


!************************* SUBROUTINE ECCP_NL   ************************
!
! this subroutine calculates the expectation value of <c|H|cp>
! where c and cp are two wavefunctions; non local part only
! for one ion only
! I have put this in a separate routine because optimization
! is than easier
!
!***********************************************************************

  SUBROUTINE ECCP_NL(LMDIM,LMMAXC,CDIJ,CPROJ1,CPROJ2,CNL)
    USE prec
    IMPLICIT NONE
    GDEF      CNL
    INTEGER LMDIM, LMMAXC
    OVERLAP CDIJ(LMDIM,LMDIM)
    GDEF CPROJ1(LMMAXC),CPROJ2(LMMAXC)
 ! local
    INTEGER L, LP

!   PROFILING_START('eccp_nl')

!DIR$ IVDEP
!OCL NOVREL
    DO L=1,LMMAXC
       DO LP=1,LMMAXC
          CNL=CNL+CDIJ(LP,L)*CPROJ1(LP)*GCONJG(CPROJ2(L))
       ENDDO
    ENDDO

!   PROFILING_STOP('eccp_nl')

  END SUBROUTINE ECCP_NL


!************************* SUBROUTINE OVERL ***************************
!
! F77 low level routine
! calculate the result of the overlap-operator acting onto a set of
! wave function characters
!
!' CRESUL^n_N,nlm = sum_n'l'm' D_N,n'l'm',nlm CPROF^n_N,n'l'm'  
!  CRESUL^n_N,i   = sum_j D_N,j,i CPROF^n_N,j  
!
! n is the band index, N the ion index, nlm=i the index for the 
! one centre partial waves
!
!
!**********************************************************************

  SUBROUTINE OVERL(WDES1, LOVERL, LMDIM, CQIJ, CPROF, CRESUL)
    USE wave
    IMPLICIT NONE
    TYPE (wavedes1) WDES1
    LOGICAL LOVERL
    INTEGER LMDIM
    OVERLAP CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    GDEF CRESUL(WDES1%NPROD,WDES1%NBANDS),CPROF(WDES1%NPROD,WDES1%NBANDS)
  ! local
    INTEGER NB, NP, ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NIS, NI, LMMAXC, L, LP

    PROFILING_START('overl')

    IF (LOVERL) THEN

       bands: DO NB=1,WDES1%NBANDS

          DO NP=1,WDES1%NPRO
             CRESUL(NP,NB)=0
          ENDDO

          spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
          DO ISPINOR_=0,WDES1%NRSPINORS-1

             NPRO =ISPINOR *(WDES1%NPRO/2)
             NPRO_=ISPINOR_*(WDES1%NPRO/2)
!#ifndef _OPENMP
             NIS =1
             DO NT=1,WDES1%NTYP
                LMMAXC=WDES1%LMMAX(NT)
                IF (LMMAXC/=0) THEN
                   DO NI=NIS,WDES1%NITYP(NT)+NIS-1
                      DO L =1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
                      DO LP=1,LMMAXC
                         CRESUL(L+NPRO,NB)=CRESUL(L+NPRO,NB)+ &
                              CQIJ(LP,L,NI,1+ISPINOR_+2*ISPINOR)*CPROF(LP+NPRO_,NB)
                      ENDDO
                      ENDDO
                      NPRO = LMMAXC+NPRO
                      NPRO_= LMMAXC+NPRO_
                   ENDDO
                ENDIF
                NIS = NIS+WDES1%NITYP(NT)
             ENDDO
!#else
!             DO NI=1,WDES1%NIONS
!                NT=WDES1%ITYP(NI)
!                LMMAXC=WDES1%LMMAX(NT)
!                NPRO2 =WDES1%LMBASE(NI)+NPRO
!                NPRO2_=WDES1%LMBASE(NI)+NPRO_
!                IF (LMMAXC==0) CYCLE
!                DO L =1,LMMAXC
!!DIR$ IVDEP
!!OCL NOVREC
!                DO LP=1,LMMAXC
!                   CRESUL(L+NPRO2,NB)=CRESUL(L+NPRO2,NB)+ &
!                        CQIJ(LP,L,NI,1+ISPINOR_+2*ISPINOR)*CPROF(LP+NPRO2_,NB)
!                ENDDO
!                ENDDO
!             ENDDO
!#endif
          ENDDO
          ENDDO spinor

       ENDDO bands
    ENDIF

    PROFILING_STOP('overl')

  END SUBROUTINE OVERL


!************************* SUBROUTINE OVERL1 **************************
!
! F77 low level routine
! calculate the result of the overlap-operator acting onto one
! wave function character
!
!  CRESUL_N,nlm = sum_n'l'm' (D_N,n'l'm',nlm-e Q_N,n'l'm',nlm) CPROF_N,n'l'm'  
!  CRESUL_N,i   = sum_j (D_N,j,i -e Q_N,j,i) CPROF_N,j  
!
! N the ion index, nlm=i the index for the one centre partial waves
!
!
!**********************************************************************

!
! scheduled for removal to be replaced by OVERL1_
!
  SUBROUTINE OVERL1(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF,CRESUL)
    USE wave
    IMPLICIT NONE

    TYPE (wavedes1) WDES1
    INTEGER LMDIM
    OVERLAP CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
         CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    REAL(q) :: EVALUE
    GDEF CRESUL(WDES1%NPRO),CPROF(WDES1%NPRO)
  ! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NIS, NI, LMMAXC, L, LP

    PROFILING_START('overl1')

    CRESUL=0
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
    DO ISPINOR_=0,WDES1%NRSPINORS-1

       NPRO =ISPINOR *(WDES1%NPRO/2)
       NPRO_=ISPINOR_*(WDES1%NPRO/2)
       NIS =1
       DO NT=1,WDES1%NTYP
          LMMAXC=WDES1%LMMAX(NT)
          IF (LMMAXC/=0) THEN
             DO NI=NIS,WDES1%NITYP(NT)+NIS-1
                IF (EVALUE==0) THEN
                   DO L =1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
                   DO LP=1,LMMAXC
                      CRESUL(L+NPRO)=CRESUL(L+NPRO)+ &
                           CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)*CPROF(LP+NPRO_)
                   ENDDO
                   ENDDO
                ELSE
                   DO L =1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
                   DO LP=1,LMMAXC
                      CRESUL(L+NPRO)=CRESUL(L+NPRO)+ &
                           (CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)- &
                           EVALUE*CQIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)) * CPROF(LP+NPRO_)
                   ENDDO
                   ENDDO
                ENDIF
                NPRO = LMMAXC+NPRO
                NPRO_= LMMAXC+NPRO_
             ENDDO
          ENDIF
          NIS = NIS+WDES1%NITYP(NT)
       ENDDO
    ENDDO
    ENDDO spinor

    PROFILING_STOP('overl1')

  END SUBROUTINE OVERL1
!
! identical to previous version but with complex EVALUE
!
!
  SUBROUTINE OVERL1_C(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF,CRESUL)
    USE wave
    IMPLICIT NONE

    TYPE (wavedes1) WDES1
    INTEGER LMDIM
    OVERLAP CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
         CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    COMPLEX(q) EVALUE
    GDEF CRESUL(WDES1%NPRO),CPROF(WDES1%NPRO)
  ! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NIS, NI, LMMAXC, L, LP

    PROFILING_START('overl1_c')

    CRESUL=0

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
    DO ISPINOR_=0,WDES1%NRSPINORS-1

       NPRO =ISPINOR *(WDES1%NPRO/2)
       NPRO_=ISPINOR_*(WDES1%NPRO/2)
       NIS =1
       DO NT=1,WDES1%NTYP
          LMMAXC=WDES1%LMMAX(NT)
          IF (LMMAXC/=0) THEN
             DO NI=NIS,WDES1%NITYP(NT)+NIS-1

                DO L =1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
                DO LP=1,LMMAXC
                   CRESUL(L+NPRO)=CRESUL(L+NPRO)+(CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)- &
                        EVALUE*CQIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)) * CPROF(LP+NPRO_)
                ENDDO
                ENDDO
                NPRO = LMMAXC+NPRO
                NPRO_= LMMAXC+NPRO_
             ENDDO
          ENDIF
          NIS = NIS+WDES1%NITYP(NT)
       ENDDO
    ENDDO
    ENDDO spinor

    PROFILING_STOP('overl1_c')

  END SUBROUTINE OVERL1_C

!
! indentical to OVERL1 but with CDIJ and CQIJ always complex
!
  SUBROUTINE OVERL1_CCDIJ(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF,CRESUL)
    USE wave
    IMPLICIT NONE

    TYPE (wavedes1) WDES1
    INTEGER LMDIM
    COMPLEX(q) CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
         CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    REAL(q) :: EVALUE
    GDEF CRESUL(WDES1%NPRO),CPROF(WDES1%NPRO)
  ! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NIS, NI, LMMAXC, L, LP

    CRESUL=0
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
    DO ISPINOR_=0,WDES1%NRSPINORS-1

       NPRO =ISPINOR *(WDES1%NPRO/2)
       NPRO_=ISPINOR_*(WDES1%NPRO/2)

       NIS =1
       DO NT=1,WDES1%NTYP
          LMMAXC=WDES1%LMMAX(NT)
          IF (LMMAXC/=0) THEN
             DO NI=NIS,WDES1%NITYP(NT)+NIS-1
                IF (EVALUE==0) THEN
                   DO L =1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
                   DO LP=1,LMMAXC
                      CRESUL(L+NPRO)=CRESUL(L+NPRO)+ &
                           CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)*CPROF(LP+NPRO_)
                   ENDDO
                   ENDDO
                ELSE
                   DO L =1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
                   DO LP=1,LMMAXC
                      CRESUL(L+NPRO)=CRESUL(L+NPRO)+ &
                           (CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)- &
                           EVALUE*CQIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)) * CPROF(LP+NPRO_)
                   ENDDO
                   ENDDO
                ENDIF
                NPRO = LMMAXC+NPRO
                NPRO_= LMMAXC+NPRO_
             ENDDO
          ENDIF
          NIS = NIS+WDES1%NITYP(NT)
       ENDDO
    ENDDO
    ENDDO spinor

  END SUBROUTINE OVERL1_CCDIJ
