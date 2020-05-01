program test 
  implicit none

  include 'mpif.h'

  integer, parameter :: dp = kind(1.0d0)
  integer rank, size, ierror, tag, status(MPI_STATUS_SIZE)
  integer :: num_wann,nkpts,nkp,c1,c2,c3,qx,qy,qz
  real(kind=dp),allocatable :: kpts(:,:) ! Input for dmft ksum library mode
  real(kind=dp),allocatable :: wght(:)   ! Input for dmft ksum library mode
  integer,allocatable :: band_win(:,:) ! Output for charge update
  real(kind=dp),allocatable :: DMFT_eval(:,:) ! Output for charge update
  complex(kind=dp),allocatable :: DMFT_evec(:,:,:) !Output for charge update
   
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
  qx=8;qy=8;qz=8;num_wann=28
  nkpts=qx*qy*qz
  allocate(band_win(2,nkpts))
  allocate(kpts(3,nkpts))
  allocate(wght(nkpts))
  wght=0.0_dp; kpts=0.0_dp; band_win=0
  
  !!!!!!! uniform mesh !!!!!!!!!!!!
  nkp=0
  do c1=1,qx
    do c2=1,qy
      do c3=1,qz
        nkp=nkp+1
        kpts(1,nkp)=dfloat(c1-qx/2-1)/dfloat(qx)
        kpts(2,nkp)=dfloat(c2-qy/2-1)/dfloat(qy)
        kpts(3,nkp)=dfloat(c3-qz/2-1)/dfloat(qz)
        wght(nkp)=1.0/dfloat(qx*qy*qz)
      enddo
    enddo
  enddo

  allocate(DMFT_eval(num_wann,nkpts))
  allocate(DMFT_evec(num_wann,num_wann,nkpts))
  DMFT_eval=0.0_dp;DMFT_evec=(0.0_dp,0.0_dp)
  call Compute_DMFT(nkpts,num_wann,kpts,wght,band_win,DMFT_eval,DMFT_evec)
  if (rank.eq.0) then 
     write(*,*) "k-point,   band,   DMFT_eval,  DMFT_evec"
     do c1=1,nkpts
        do c2=1,num_wann
           write(*,*) c1, band_win(1,nkpts)+c2-1, DMFT_eval(c2,c1),DMFT_evec(1,c2,c1)
        enddo
     enddo
  endif
  call MPI_FINALIZE(ierror)
end program test   
