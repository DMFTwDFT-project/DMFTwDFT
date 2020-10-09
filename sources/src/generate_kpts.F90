
module generate_kpts
  !! This module contains parameters to control the actions of wannier90.
  !! Also routines to read the parameters and write them out again.

  use constants, only: dp
  use io, only: stdout, maxlen
  use read_inputs 

  implicit none

  integer, save :: n_kpts
  real(kind=dp),allocatable, save :: kpts(:,:)
  real(kind=dp),allocatable, save :: weight(:)
  integer, save :: num_new_kpts, nfine(3)
  

contains
 
  subroutine generate_uniform_kmesh()
    use constants
    use io
    use read_inputs, only: qx,qy,qz,nspin
   
    implicit none

    integer :: nkp,c1,c2,c3 
!    integer :: qx,qy,qz 


!    qx=24;qy=24;qz=24
    num_new_kpts=qx*qx*(qz/2+1)
    allocate(kpts(3,num_new_kpts))
    allocate(weight(num_new_kpts))
    kpts=0.0_dp
    weight=0.0_dp
    nkp=0
    if (MOD(qz,2) .eq. 0) then
      do c1=1,qx
        do c2=1,qy
          do c3=1,qz/2+1
            nkp=nkp+1
            kpts(1,nkp)=dfloat(c1-qx/2-1)/dfloat(qx)
            kpts(2,nkp)=dfloat(c2-qy/2-1)/dfloat(qy)
            kpts(3,nkp)=dfloat(c3-qz/2-1)/dfloat(qz)
            if (c3.eq.1 .OR. c3.eq.(qz/2+1)) then
              weight(nkp)=1.0/dfloat(qx*qy*qz)/nspin
            else
              weight(nkp)=2.0/dfloat(qx*qy*qz)/nspin
            endif
          enddo
        enddo
      enddo
    else
      do c1=1,qx
        do c2=1,qy
          do c3=1,qz/2+1
            nkp=nkp+1
            kpts(1,nkp)=dfloat(c1-qx/2-1)/dfloat(qx)
            kpts(2,nkp)=dfloat(c2-qy/2-1)/dfloat(qy)
            kpts(3,nkp)=dfloat(c3-qz/2-1)/dfloat(qz)
            if (c3.eq.(qz/2+1)) then
              weight(nkp)=1.0/dfloat(qx*qy*qz)/nspin
            else
              weight(nkp)=2.0/dfloat(qx*qy*qz)/nspin
            endif
          enddo
        enddo
      enddo
    endif
  end subroutine generate_uniform_kmesh

  subroutine generate_dense_kmesh(n_kpts_loc,kpt_dft,wght_dft)
    use constants
    use io
    use read_inputs, only: mp_grid,qx,qy,qz,nspin
   
    implicit none

    integer, intent(in) :: n_kpts_loc
    real(kind=dp), intent(in) :: kpt_dft(3,n_kpts_loc)
    real(kind=dp), intent(in) :: wght_dft(n_kpts_loc)
    integer :: nc,nk,nkp,c1,c2,c3 
!    integer :: qx,qy,qz 


    nfine(1)=int(qx/mp_grid(1))
    nfine(2)=int(qy/mp_grid(2))
    nfine(3)=int(qz/mp_grid(3))
    !write(*,*) kpts
    num_new_kpts=nfine(1)*nfine(2)*nfine(3)*n_kpts_loc
    allocate(kpts(3,num_new_kpts))
    allocate(weight(num_new_kpts))
    kpts=0.0_dp
    weight=0.0_dp
    nkp=0
    do nk=1,n_kpts_loc
      do c1=1,nfine(3)
        do c2=1,nfine(2)
          do c3=1,nfine(1)
            nkp=nkp+1
            if (mod(nfine(1),2) .eq. 0) then
              nc=nfine(1)/2
              kpts(1,nkp)=kpt_dft(1,nk)+(dfloat(c3-nc)-0.5)/dfloat(nfine(1)*mp_grid(1))
            else
              nc=(nfine(1)+1)/2
              kpts(1,nkp)=kpt_dft(1,nk)+dfloat(c3-nc)/dfloat(nfine(1)*mp_grid(1))
            !  write(*,*) kpts(1,nkp)
            endif
            if (mod(nfine(2),2) .eq. 0) then
              nc=nfine(2)/2
              kpts(2,nkp)=kpt_dft(2,nk)+(dfloat(c2-nc)-0.5)/dfloat(nfine(2)*mp_grid(2))
            else
              nc=(nfine(2)+1)/2
              kpts(2,nkp)=kpt_dft(2,nk)+dfloat(c2-nc)/dfloat(nfine(2)*mp_grid(2))
            !  write(*,*) kpts(2,nkp)
            endif
            if (mod(nfine(3),2) .eq. 0) then
              nc=nfine(3)/2
              kpts(3,nkp)=kpt_dft(3,nk)+(dfloat(c1-nc)-0.5)/dfloat(nfine(3)*mp_grid(3))
            else
              nc=(nfine(3)+1)/2
              kpts(3,nkp)=kpt_dft(3,nk)+dfloat(c1-nc)/dfloat(nfine(3)*mp_grid(3))
            !  write(*,*) kpts(3,nkp)
            endif
            weight(nkp)=wght_dft(nk)/dfloat(nfine(1)*nfine(2)*nfine(3))/nspin
            !write(*,*) weight(nkp)
          enddo
        enddo
      enddo
    enddo
    !write(*,*) kpts
    !write(*,*) weight
    n_kpts=n_kpts_loc
  end subroutine generate_dense_kmesh

end module generate_kpts
 
