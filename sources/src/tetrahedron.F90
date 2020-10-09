module tetrahedron 
  !! This module contains parameters to control the actions of wannier90.
  !! Also routines to read the parameters and write them out again.

  use constants, only: dp
  use io, only: stdout, maxlen
  use read_inputs 
  use utility 

  implicit none

!  integer, save :: n_kpts
!  real(kind=dp),allocatable, save :: kpts(:,:)
!  real(kind=dp),allocatable, save :: weight(:)
  integer, save :: num_new_kpts!, nfine(3)
  integer, allocatable, save :: tetptr(:,:)
  integer, allocatable, save :: tet_idx(:) 
  integer, allocatable, save :: kibz(:,:) !diagonal directions
  

contains
  function dist_sq_int(a,b,dim)
    use constants, only: dp
    integer :: dim, c0 
    integer :: a(dim) 
    integer :: b(dim) 
    integer :: dist_sq_int

    dist_sq_int=0
    do c0=1,dim
      dist_sq_int=dist_sq_int+(a(c0)-b(c0))**2
    enddo
  end function dist_sq_int

  function dist_sq(a,b,dim)
    use constants, only: dp
    integer :: dim, c0 
    real(kind=dp) :: a(dim) 
    real(kind=dp) :: b(dim) 
    real(kind=dp) :: dist_sq

    dist_sq=0.0_dp
    do c0=1,dim
      dist_sq=dist_sq+(a(c0)-b(c0))**2
    enddo
  end function dist_sq

  function comp(a,b)
    use constants, only: dp
    integer :: comp
    !real(kind=dp) :: a
    !real(kind=dp) :: b
    integer :: a
    integer :: b

    if ( a .lt. b ) comp = -1
    if ( a .eq. b ) comp = 0
    if ( a .gt. b ) comp = 1
  end function comp

  SUBROUTINE sort(X,N)
    ! this sorts in descending order...
    IMPLICIT NONE
    INTEGER N, I, J
    REAL(kind=dp) X(N)

    DO 100 I=2,N
      IF ( X(I).GT.X(I-1) ) THEN
        DO 50 J=I-2,1,-1
          IF(X(I).LT.X(J)) go to 70
50        CONTINUE
        J=0
70      x(j+1:i) = cshift(x(j+1:i),-1)
        !iy(j+1:i) = cshift(iy(j+1:i),-1)
      ENDIF
100 CONTINUE
    RETURN
  END SUBROUTINE

  SUBROUTINE sortr(X,iy,N)
      ! sorts a real list in ascending order...
    IMPLICIT NONE
    INTEGER N, I, J, iy(n)
    REAL(kind=dp) X(N)

    do i=1,n;iy(i)=i;enddo
    DO 100 I=2,N
      IF ( X(I).LT.X(I-1) ) THEN
        DO 50 J=I-2,1,-1
          IF(X(I).GT.X(J)) go to 70
50        CONTINUE
        J=0
70      x(j+1:i) = cshift(x(j+1:i),-1)
        iy(j+1:i) = cshift(iy(j+1:i),-1)
      ENDIF
100 CONTINUE
    RETURN
  END SUBROUTINE

  SUBROUTINE sorti(X,iy,N)
    ! sort in increasing order...
    IMPLICIT NONE
    INTEGER N, I, J
    integer X(N),iy(n)

    do i=1,n;iy(i)=i;enddo
    DO 100 I=2,N
       IF ( X(I).LT.X(I-1) ) THEN
          DO 50 J=I-2,1,-1
            IF(X(I).GT.X(J)) go to 70
50          CONTINUE
          J=0
70        x(j+1:i)  = cshift(x(j+1:i),-1)
          iy(j+1:i) = cshift(iy(j+1:i),-1)
       ENDIF
100 CONTINUE
    RETURN
  END SUBROUTINE
 
  subroutine get_tetra()
   
    implicit none

    integer :: verts(8,3) !vertices
    integer :: diags(4,2) !diagonal directions
    integer :: temp_verts(6,3) !diagonal directions
    integer :: ptet(6,4,3) !diagonal directions
    integer :: pptet(4,3) !diagonal directions
    integer :: sym(1,3,3) !diagonal directions
    integer :: kt((mp_grid(1)+1)*(mp_grid(2)+1)*(mp_grid(3)+1),3) !diagonal directions
    integer :: kptr((mp_grid(1)+1)*(mp_grid(2)+1)*(mp_grid(3)+1),3) !diagonal directions
    integer :: gptr((mp_grid(1)+1)*(mp_grid(2)+1)*(mp_grid(3)+1)) !diagonal directions
    integer :: tet((mp_grid(1))*(mp_grid(2))*(mp_grid(3))*6,4) 
    real(kind=dp) :: mdiag(4) !diagonal directions
    real(kind=dp) :: mdiag_orig(4) !diagonal directions
!    real(kind=dp), intent(in) :: kpt_dft(3,n_kpts_loc)
!    real(kind=dp), intent(in) :: wght_dft(n_kpts_loc)
    integer :: ntet,nk,nkp,c1,c2,c3 
    integer :: idx,cidx,qx,qy,qz 

    qx=mp_grid(1); qy=mp_grid(2); qz=mp_grid(3)
    if (.not. allocated(tetptr)) allocate(tetptr(qx*qy*qz*6,4))
    if (.not. allocated(tet_idx)) allocate(tet_idx(qx*qy*qz))
    if (.not. allocated(kibz)) allocate(kibz(qx*qy*qz,3))

    num_new_kpts=qx*qy*qz
    !recip_latt(1,:)=(/ 2.0,2.0,2.0 /)
    !recip_latt(2,:)=(/ 3.0,1.0,3.0 /)
    !recip_latt(3,:)=(/ 4.0,3.0,6.0 /)
    verts(1,:)=(/ 0,0,0 /)
    verts(2,:)=(/ 1,0,0 /)
    verts(3,:)=(/ 0,1,0 /)
    verts(4,:)=(/ 1,1,0 /)
    verts(5,:)=(/ 0,0,1 /)
    verts(6,:)=(/ 1,0,1 /)
    verts(7,:)=(/ 0,1,1 /)
    verts(8,:)=(/ 1,1,1 /)
    diags(1,:)=(/ 1,8 /)
    diags(2,:)=(/ 2,7 /)
    diags(3,:)=(/ 3,6 /)
    diags(4,:)=(/ 4,5 /)
    sym(1,1,:)=(/ 1,0,0 /)
    sym(1,2,:)=(/ 0,1,0 /)
    sym(1,3,:)=(/ 0,0,1 /)

    do c1=1,4
      mdiag(c1)=dist_sq(matmul(recip_latt,verts(diags(c1,1),:)),matmul(recip_latt,verts(diags(c1,2),:)),3)
      mdiag_orig(c1)=mdiag(c1)
    enddo
!    print *, mdiag
    call quicksort(mdiag)
    do c1=1,4
      if (abs(mdiag_orig(c1)-mdiag(1)) .lt. 1e-4) then
        cidx=c1
        exit
      endif
    enddo
    idx=1
    do c1=1,8
      if (c1 .ne. diags(cidx,1) .AND. c1 .ne. diags(cidx,2)) then
        temp_verts(idx,:)=verts(c1,:)
        idx=idx+1
      endif
    enddo
    idx=1
    do c1=1,6
      do c2=c1+1,6
        if (dist_sq_int(temp_verts(c1,:),temp_verts(c2,:),3) .eq. 1) then
          pptet(1,:)=verts(diags(cidx,1),:) 
          !print *, verts(diags(cidx,1),:)
          pptet(2,:)=verts(diags(cidx,2),:) 
          !print *, verts(diags(cidx,2),:)
          pptet(3,:)=temp_verts(c1,:)
          !print *, temp_verts(c1,:)        
          pptet(4,:)=temp_verts(c2,:)
          !print *, temp_verts(c2,:)          
          !call quicksort(pptet)
          ! I am not sure if I have to sort here since the final result is sorted anyway
          ptet(idx,:,:)=pptet(:,:)
          idx=idx+1
        endif
      enddo
    enddo

    idx=1
    ntet=qx*qy*qz*6
    do c3=1,qz+1
      do c2=1,qy+1
        do c1=1,qx+1
          kt(idx,:)=(/ c1-1,c2-1,c3-1 /)
          idx=idx+1
        enddo
      enddo
    enddo
    
    call get_kptr(mp_grid,1,(qx+1)*(qy+1)*(qz+1),sym,kt,kptr,gptr)
    
    ! If sym=0, kibz is the same as full kpoints
    idx=1
    do c1=1,qx
      do c2=1,qy
        do c3=1,qz
          kibz(idx,:)=(/ c1-1,c2-1,c3-1 /)
          idx=idx+1
        enddo
      enddo
    enddo
   
    do c1=1,qx*qy*qz
      do c2=1,qx*qy*qz
        if (kibz(c1,1).eq.ikpt(1,c2) .and. kibz(c1,2).eq.ikpt(2,c2) .and. kibz(c1,3).eq.ikpt(3,c2)) then
          tet_idx(c1)=c2
          exit
        endif
      enddo
    enddo

    call get_itet(mp_grid,ptet,qx*qy*qz,(qx+1)*(qy+1)*(qz+1),kibz,kptr,ntet,tetptr,tet)
!    print *, tetptr(1,:)
!    print *, tetptr(2,:)
!    print *, tetptr(ntet-1,:)
!    print *, tetptr(ntet,:)
!    print *, ikpt
!    print *, tet_idx

  end subroutine get_tetra

  subroutine get_kptr(nk,nop,nkt,symm,kt,kptr,gptr)
    ! this routine returns a pointer from each kpoint to the ikp...
    implicit none
    integer, intent(in) :: nk(3),symm(nop,3,3),kt(nkt,3)
    integer, intent(out) :: kptr(nkt,3),gptr(nkt)
    integer, intent(in) :: nop,nkt
    integer :: g,k,temp(3),tt,i

    do k=1,nkt
      tt=kt(k,1)+(nk(1)+1)*(kt(k,2)+(nk(2)+1)*kt(k,3))
      !mo=kt(k,:)
      kptr(k,:)=kt(k,:)
      gptr(k)=0
      do g=1,nop
        temp=matmul(symm(g,:,:),kt(k,:))
        do i=1,3
          if (temp(i)<0)then; temp(i)=temp(i)+nk(i)
          elseif (temp(i) .ge. nk(i))then; temp(i)=temp(i)-nk(i)
          endif
        enddo
        !print *, 'temp=',temp
        if (temp(1)+(nk(1)+1)*(temp(2)+(nk(2)+1)*temp(3)) .le.tt)then
          tt=temp(1)+(nk(1)+1)*(temp(2)+(nk(2)+1)*temp(3))
          !gptr(k)=g-1   ! this is the python convention... not using.
          gptr(k)=g
          kptr(k,:)=temp(:)
        endif
      enddo
      if (gptr(k).eq.0) then
        print*,'Problem in get_kptr...'
        STOP
      endif
    enddo

  end subroutine get_kptr


  subroutine get_itet(nk,ptet,nibz,nkptr,kibz,kptr,ntet,tetptr,tet)
    implicit none
    integer, intent(in) :: nk(3),ptet(6,4,3),kibz(nibz,3), nibz,ntet
    integer, intent(out) :: tetptr(ntet,4),tet(ntet,4) 
    integer, intent(in) ::  kptr(nkptr,3),nkptr
    integer :: order(4)
    integer i,j,k,l,p,t,temp(3),cc,temp2(4),temp3(4)

    cc=1
    ! i,j,k loop over kpoints...
    do i=0,nk(1)-1
      do j=0,nk(2)-1
        do k=0,nk(3)-1
          ! p loops over the primitive tetra...
          do p=1,6
            ! t loops over the vertices of the tetra...
            do t=1,4
              temp=0
              ! temp is a given vertex of a tetra...
              temp(1)=ptet(p,t,1)+i
              temp(2)=ptet(p,t,2)+j
              temp(3)=ptet(p,t,3)+k
              ! store vertex for use below...
              temp3(t)=1+temp(1)+(nk(1)+1)*(temp(2)+(nk(2)+1)*temp(3))
              ! now get the corresponding irreducible kpoint...
              temp=kptr(temp3(t),:)
              ! now find the INDEX of that irreducible kpoint...
              do l=1,nibz
                if (temp(1) .eq. kibz(l,1) .and. temp(2) .eq. kibz(l,2) .and. temp(3) .eq. kibz(l,3))then
                  !tetptr(cc,t)=l-1
                  tetptr(cc,t)=l
                  exit
                elseif (l==nibz) then; stop 'Didnt find kpoint in IBZ'
                endif
              enddo
            enddo
            temp2=tetptr(cc,:)
            ! now sort in order of ikp...
            call sorti(temp2,order,4)
            ! now we have the sorted pointer to the INDEX of the ikp...
            tetptr(cc,:)=temp2
            ! now reorder to original tet so it corresponds with ikp and store
            ! for output...
            tet(cc,1)=temp3(order(1))
            tet(cc,2)=temp3(order(2))
            tet(cc,3)=temp3(order(3))
            tet(cc,4)=temp3(order(4))
            cc=cc+1
          enddo
        enddo
      enddo
    enddo

  end subroutine get_itet

  subroutine tetra_totdos(nbin,norb,nki,ntet,tet_min,tet_max,kibz,itet,tet_idx,tval,emin,emax,tdos,cdos)
    implicit none
    integer, intent(in) :: nbin,norb,nki,ntet,kibz(nki,3),itet(ntet,4), tet_idx(ntet), tet_min,tet_max
    real(kind=dp), intent(in) :: tval(nki,norb), emin, emax
    real(kind=dp), intent(out) :: tdos(nbin+1), cdos(nbin+1)
    real(kind=dp) :: temp(4),e1,e2,e3,e4,en,ttemp,ctemp,deps
    integer te,eps,b

    ! define the deps...
    if (nbin==0) then
      deps=0
    else
      deps=(emax-emin)/nbin
    endif

    tdos=0.0_dp
    cdos=0.0_dp
    ! loop over all irred tetra and compute DOS...
    do eps=1,nbin+1
      en=emin+(eps-1)*deps
      !do te=1,ntet
      do te=tet_min,tet_max
        do b=1,norb
          !print *, itet(te,1)
          !print *, itet(te,2)
          !print *, itet(te,3)
          !print *, itet(te,4)
          temp(1)=tval(tet_idx(itet(te,1)),b)
          temp(2)=tval(tet_idx(itet(te,2)),b)
          temp(3)=tval(tet_idx(itet(te,3)),b)
          temp(4)=tval(tet_idx(itet(te,4)),b)
          call sort(temp,4)
          !print *, temp
          e1=temp(4);e2=temp(3);e3=temp(2);e4=temp(1)
          ctemp=0;ttemp=0
          if (en>e1 .and. en<e2) then
            ttemp=3*(en-e1)**2/((e2-e1)*(e3-e1)*(e4-e1))
            ctemp=(en-e1)**3/((e2-e1)*(e3-e1)*(e4-e1))
          elseif (en>e2 .and. en<e3) then
            ttemp=1/((e3-e1)*(e4-e1))*(3*(e2-e1)+6*(en-e2)-3*((e3-e1+e4-e2)*(en-e2)**2)/((e3-e2)*(e4-e2)))
            ctemp=1/((e3-e1)*(e4-e1))*((e2-e1)**2+3*(e2-e1)*(en-e2)+3*(en-e2)**2-(e3-e1+e4-e2)/((e3-e2)*(e4-e2))*(en-e2)**3)
          elseif (en>e3 .and. en<e4) then
            ttemp=3*(e4-en)**2/((e4-e1)*(e4-e2)*(e4-e3))
            ctemp=1-(e4-en)**3/((e4-e1)*(e4-e2)*(e4-e3))
          elseif (en>e4) then
            ctemp=1
          endif
          !tdos(eps)=tdos(eps)+ttemp*tetmult(te)
          !cdos(eps)=cdos(eps)+ctemp*tetmult(te)
          tdos(eps)=tdos(eps)+ttemp
          cdos(eps)=cdos(eps)+ctemp
        enddo
      enddo
      !tdos(eps)=tdos(eps)/tottet ! put in vt/vg
      !cdos(eps)=cdos(eps)/tottet ! put in vt/vg
      tdos(eps)=tdos(eps)/ntet ! put in vt/vg
      cdos(eps)=cdos(eps)/ntet ! put in vt/vg
    enddo
  end subroutine

  subroutine tetra_poccup(norb,nki,ntet,tet_min,tet_max,kibz,tetkptr,tet_idx,tval,en,cpdos)
    implicit none
    integer, intent(in) :: norb,nki,ntet,kibz(nki,3),tetkptr(ntet,4), tet_idx(ntet), tet_min,tet_max
    real(kind=dp), intent(in) :: tval(nki,norb), en
    real(kind=dp), intent(out) :: cpdos(nki,norb)
    real(kind=dp) :: temp(4),e1,e2,e3,e4,ttemp,ctemp,deps
    real(kind=dp) :: c,c1,c2,c3,cw1,cw2,cw3,cw4
    integer te,eps,b,i1,i2,i3,i4, order(4)

    cpdos=0.0_dp
        ! loop over the tetra...
    do t=tet_min,tet_max
      do b=1,norb
        temp(1)=tval(tet_idx(tetkptr(t,1)),b)
        temp(2)=tval(tet_idx(tetkptr(t,2)),b)
        temp(3)=tval(tet_idx(tetkptr(t,3)),b)
        temp(4)=tval(tet_idx(tetkptr(t,4)),b)
        call sortr(temp,order,4)
        e1=temp(1);e2=temp(2);e3=temp(3);e4=temp(4)
          ! store the projections with proper order for given band...
        i1=tet_idx(tetkptr(t,order(1)))
        i2=tet_idx(tetkptr(t,order(2)))
        i3=tet_idx(tetkptr(t,order(3)))
        i4=tet_idx(tetkptr(t,order(4)))
        cw1=0.0_dp;cw2=0.0_dp;cw3=0.0_dp;cw4=0.0_dp
        if (en>e1 .and. en<e2) then
          c=0.25*(en-e1)**3/((e2-e1)*(e3-e1)*(e4-e1))
          cw1=c*(4-(en-e1)*(1/(e2-e1)+1/(e3-e1)+1/(e4-e1)))
          cw2=c*(en-e1)/(e2-e1)
          cw3=c*(en-e1)/(e3-e1)
          cw4=c*(en-e1)/(e4-e1)
        elseif (en>e2 .and. en<e3) then
          c1=0.25*(en-e1)**2/((e4-e1)*(e3-e1))
          c2=0.25*((en-e1)*(en-e2)*(e3-en))/((e4-e1)*(e3-e2)* &
                   (e3-e1))
          c3=0.25*(en-e2)**2*(e4-en)/((e4-e2)*(e3-e2)*(e4-e1))
          cw1=c1+(c1+c2)*(e3-en)/(e3-e1)+(c1+c2+c3)*(e4-en)/(e4-e1)
          cw2=c1+c2+c3+(c2+c3)*(e3-en)/(e3-e2)+c3*(e4-en)/(e4-e2)
          cw3=(c1+c2)*(en-e1)/(e3-e1)+(c2+c3)*(en-e2)/(e3-e2)
          cw4=(c1+c2+c3)*(en-e1)/(e4-e1)+c3*(en-e2)/(e4-e2)
        elseif (en>e3 .and. en<e4) then
          c=0.25*(e4-en)**3/((e4-e1)*(e4-e2)*(e4-e3))
          cw1=0.25-c*(e4-en)/(e4-e1)
          cw2=0.25-c*(e4-en)/(e4-e2)
          cw3=0.25-c*(e4-en)/(e4-e3)
          cw4=0.25-c*(4-(1/(e4-e1)+1/(e4-e2)+1/(e4-e3))*(e4-en))
        elseif (en>e4) then
          cw1=0.25;cw2=0.25;cw3=0.25;cw4=0.25
        endif

        cpdos(i1,b)=cpdos(i1,b)+cw1/6.0
        cpdos(i2,b)=cpdos(i2,b)+cw2/6.0
        cpdos(i3,b)=cpdos(i3,b)+cw3/6.0
        cpdos(i4,b)=cpdos(i4,b)+cw4/6.0
      enddo
    enddo
  end subroutine

  subroutine tetra_pdos(nbin,norb,nki,ntet,tet_min,tet_max,kibz,tetkptr,tet_idx,tval,emin,emax,pdos,cpdos)
    implicit none
    integer, intent(in) :: nbin,norb,nki,ntet,kibz(nki,3),tetkptr(ntet,4), tet_idx(ntet), tet_min,tet_max
    real(kind=dp), intent(in) :: tval(nki,norb), emin, emax
    real(kind=dp), intent(out) :: pdos(nbin+1,nki,norb), cpdos(nbin+1,nki,norb)
    real(kind=dp) :: temp(4),e1,e2,e3,e4,en,ttemp,ctemp,deps
    real(kind=dp) :: c,c1,c2,c3,cw1,cw2,cw3,cw4,w1,w2,w3,w4
    integer te,eps,b,i1,i2,i3,i4, order(4)

      ! define the deps...
    if (nbin==0) then
      deps=0
    else
      deps=(emax-emin)/nbin
    endif

    pdos=0.0_dp
    cpdos=0.0_dp
      ! start loop over energy...
    do eps=1,nbin+1
      en=emin+(eps-1)*deps
        ! loop over the tetra...
      do t=tet_min,tet_max
!          do i=1,4
!            ! compute projections for each kpoint in tetra...
!            ctemp=matmul(transpose(gred(gptr(tet(t,i)),:,:)),
!     1                 tvec(tetkptr(t,i),:,:))
!            !ctemp=tvec(tetkptr(t,i),:,:)
!            ptemp0(i,:,:)=real(conjg(ctemp)*ctemp)
!          enddo
        do b=1,norb
          temp(1)=tval(tet_idx(tetkptr(t,1)),b)
          temp(2)=tval(tet_idx(tetkptr(t,2)),b)
          temp(3)=tval(tet_idx(tetkptr(t,3)),b)
          temp(4)=tval(tet_idx(tetkptr(t,4)),b)
          call sortr(temp,order,4)
          e1=temp(1);e2=temp(2);e3=temp(3);e4=temp(4)
            ! store the projections with proper order for given band...
          i1=tet_idx(tetkptr(t,order(1)))
          i2=tet_idx(tetkptr(t,order(2)))
          i3=tet_idx(tetkptr(t,order(3)))
          i4=tet_idx(tetkptr(t,order(4)))
  !          ptemp(1,:)=ptemp0(order(1),:,b)
  !          ptemp(2,:)=ptemp0(order(2),:,b)
  !          ptemp(3,:)=ptemp0(order(3),:,b)
  !          ptemp(4,:)=ptemp0(order(4),:,b)
            ! now compute the weights...
          w1=0;w2=0;w3=0;w4=0;cw1=0;cw2=0;cw3=0;cw4=0
          if (en>e1 .and. en<e2) then
            c=0.25*(en-e1)**3/((e2-e1)*(e3-e1)*(e4-e1))
            cw1=c*(4-(en-e1)*(1/(e2-e1)+1/(e3-e1)+1/(e4-e1)))
            cw2=c*(en-e1)/(e2-e1)
            cw3=c*(en-e1)/(e3-e1)
            cw4=c*(en-e1)/(e4-e1)
        w1=0.25*(en-e1)**3*(-1/(e2-e1)-1/(e3-e1)-1/(e4-e1))/((e2-e1)* & 
          (e3-e1)*(e4-e1))+0.75*(en-e1)**2*(4-(en-e1)*(1/(e2-e1)+1/   &
          (e3-e1)+1/(e4-e1)))/((e2-e1)*(e3-e1)*(e4-e1))
        w2=(en-e1)**3/((e2-e1)**2*(e3-e1)*(e4-e1))
        w3=(en-e1)**3/((e2-e1)*(e3-e1)**2*(e4-e1))
        w4=(en-e1)**3/((e2-e1)*(e3-e1)*(e4-e1)**2)
          elseif (en>e2 .and. en<e3) then
            c1=0.25*(en-e1)**2/((e4-e1)*(e3-e1))
            c2=0.25*((en-e1)*(en-e2)*(e3-en))/((e4-e1)*(e3-e2)* &
                     (e3-e1))
            c3=0.25*(en-e2)**2*(e4-en)/((e4-e2)*(e3-e2)*(e4-e1))
            cw1=c1+(c1+c2)*(e3-en)/(e3-e1)+(c1+c2+c3)*(e4-en)/(e4-e1)
            cw2=c1+c2+c3+(c2+c3)*(e3-en)/(e3-e2)+c3*(e4-en)/(e4-e2)
            cw3=(c1+c2)*(en-e1)/(e3-e1)+(c2+c3)*(en-e2)/(e3-e2)
            cw4=(c1+c2+c3)*(en-e1)/(e4-e1)+c3*(en-e2)/(e4-e2)
         w1=-(0.25*(e3-en)*(en-e1)*(en-e2)/((e3-e1)*(e3-e2)*(e4-e1))+   &
           0.25*(en-e1)**2/((e3-e1)*(e4-e1)))/(e3-e1)-(0.25*(en-e2)**2  & 
           *(e4-en)/((e3-e2)*(e4-e1)*(e4-e2))+0.25*(e3-en)*(en-e1)*     &
           (en-e2)/((e3-e1)*(e3-e2)*(e4-e1))+0.25*(en-e1)**2/((e3-e1)*  &
           (e4-e1)))/(e4-e1)+(e3-en)*(0.25*(-2*e1+2*en)/((e3-e1)*       &
           (e4-e1))+0.25*(e3-en)*(en-e1)/((e3-e1)*(e3-e2)*(e4-e1))+     &
           0.25*(e3-en)*(en-e2)/((e3-e1)*(e3-e2)*(e4-e1))-0.25*(en-e1)* &
           (en-e2)/((e3-e1)*(e3-e2)*(e4-e1)))/(e3-e1)+(e4-en)*(0.25*    &
           (-2*e1+2*en)/((e3-e1)*(e4-e1))+0.25*(e3-en)*(en-e1)/((e3-e1) &
           *(e3-e2)*(e4-e1))+0.25*(e3-en)*(en-e2)/((e3-e1)*(e3-e2)*     &
           (e4-e1))+0.25*(e4-en)*(-2*e2+2*en)/((e3-e2)*(e4-e1)*(e4-e2)) &
           -0.25*(en-e1)*(en-e2)/((e3-e1)*(e3-e2)*(e4-e1))-0.25*        &
           (en-e2)**2/((e3-e2)*(e4-e1)*(e4-e2)))/(e4-e1)+0.25*(-2*e1+   &
           2*en)/((e3-e1)*(e4-e1))
        w2=-(0.25*(en-e2)**2*(e4-en)/((e3-e2)*(e4-e1)*(e4-e2))+0.25*    &
           (e3-en)*(en-e1)*(en-e2)/((e3-e1)*(e3-e2)*(e4-e1)))/(e3-e2)+  &
           (e3-en)*(0.25*(e3-en)*(en-e1)/((e3-e1)*(e3-e2)*(e4-e1))+     &
           0.25*(e3-en)*(en-e2)/((e3-e1)*(e3-e2)*(e4-e1))+0.25*(e4-en)* &
           (-2*e2+2*en)/((e3-e2)*(e4-e1)*(e4-e2))-0.25*(en-e1)*(en-e2)  &
           /((e3-e1)*(e3-e2)*(e4-e1))-0.25*(en-e2)**2/((e3-e2)*(e4-e1)  &
           *(e4-e2)))/(e3-e2)+0.25*(-2*e1+2*en)/((e3-e1)*(e4-e1))+0.25  &
           *(e3-en)*(en-e1)/((e3-e1)*(e3-e2)*(e4-e1))+0.25*(e3-en)*     &
           (en-e2)/((e3-e1)*(e3-e2)*(e4-e1))+0.25*(e4-en)*(-2*e2+2*en)  &
           /((e3-e2)*(e4-e1)*(e4-e2))+0.25*(e4-en)**2*(-2*e2+2*en)/     &
           ((e3-e2)*(e4-e1)*(e4-e2)**2)+0.25*(en-e2)**2*(-2*e4+2*en)/   &
           ((e3-e2)*(e4-e1)*(e4-e2)**2)-0.25*(en-e1)*(en-e2)/((e3-e1)*  &
           (e3-e2)*(e4-e1))-0.25*(en-e2)**2/((e3-e2)*(e4-e1)*(e4-e2))   
        w3=(0.25*(e3-en)*(en-e1)*(en-e2)/((e3-e1)*(e3-e2)*(e4-e1))+     &
           0.25*(en-e1)**2/((e3-e1)*(e4-e1)))/(e3-e1)+(0.25*(en-e2)**2  &
           *(e4-en)/((e3-e2)*(e4-e1)*(e4-e2))+0.25*(e3-en)*(en-e1)*     &
           (en-e2)/((e3-e1)*(e3-e2)*(e4-e1)))/(e3-e2)+(en-e1)*(0.25*    &
           (-2*e1+2*en)/((e3-e1)*(e4-e1))+0.25*(e3-en)*(en-e1)/         &
           ((e3-e1)*(e3-e2)*(e4-e1))+0.25*(e3-en)*(en-e2)/((e3-e1)*     &
           (e3-e2)*(e4-e1))-0.25*(en-e1)*(en-e2)/((e3-e1)*(e3-e2)*      &
           (e4-e1)))/(e3-e1)+(en-e2)*(0.25*(e3-en)*(en-e1)/((e3-e1)*    &
           (e3-e2)*(e4-e1))+0.25*(e3-en)*(en-e2)/((e3-e1)*(e3-e2)*      &
           (e4-e1))+0.25*(e4-en)*(-2*e2+2*en)/((e3-e2)*(e4-e1)*         &
           (e4-e2))-0.25*(en-e1)*(en-e2)/((e3-e1)*(e3-e2)*(e4-e1))-     &
           0.25*(en-e2)**2/((e3-e2)*(e4-e1)*(e4-e2)))/(e3-e2)         
        w4=(0.25*(en-e2)**2*(e4-en)/((e3-e2)*(e4-e1)*(e4-e2))+0.25*     &
           (e3-en)*(en-e1)*(en-e2)/((e3-e1)*(e3-e2)*(e4-e1))+0.25*      &
           (en-e1)**2/((e3-e1)*(e4-e1)))/(e4-e1)+(en-e1)*(0.25*         &
           (-2*e1+2*en)/((e3-e1)*(e4-e1))+0.25*(e3-en)*(en-e1)/((e3-e1) &
           *(e3-e2)*(e4-e1))+0.25*(e3-en)*(en-e2)/((e3-e1)*(e3-e2)*     &
           (e4-e1))+0.25*(e4-en)*(-2*e2+2*en)/((e3-e2)*(e4-e1)*         &
           (e4-e2))-0.25*(en-e1)*(en-e2)/((e3-e1)*(e3-e2)*(e4-e1))-     &
           0.25*(en-e2)**2/((e3-e2)*(e4-e1)*(e4-e2)))/(e4-e1)+0.75*     &
           (en-e2)**2*(e4-en)/((e3-e2)*(e4-e1)*(e4-e2)**2)-0.25*        &
           (en-e2)**3/((e3-e2)*(e4-e1)*(e4-e2)**2)
          elseif (en>e3 .and. en<e4) then
            c=0.25*(e4-en)**3/((e4-e1)*(e4-e2)*(e4-e3))
            cw1=0.25-c*(e4-en)/(e4-e1)
            cw2=0.25-c*(e4-en)/(e4-e2)
            cw3=0.25-c*(e4-en)/(e4-e3)
            cw4=0.25-c*(4-(1/(e4-e1)+1/(e4-e2)+1/(e4-e3))*(e4-en))
        w1=(e4-en)**3/((e4-e1)**2*(e4-e2)*(e4-e3))
        w2=(e4-en)**3/((e4-e1)*(e4-e2)**2*(e4-e3))
        w3=(e4-en)**3/((e4-e1)*(e4-e2)*(e4-e3)**2)
        w4=0.75*(e4-en)**2*(4-(e4-en)*(1/(e4-e1)+1/(e4-e2)+1/(e4-e3)))/  &
        ((e4-e1)*(e4-e2)*(e4-e3))-0.25*(e4-en)**3*(1/(e4-e1)+1/          &
        (e4-e2)+1/(e4-e3))/((e4-e1)*(e4-e2)*(e4-e3))
          elseif (en>e4) then
            cw1=0.25;cw2=0.25;cw3=0.25;cw4=0.25
          endif

          pdos(eps,i1,b)=pdos(eps,i1,b)+w1/6.0
          pdos(eps,i2,b)=pdos(eps,i2,b)+w2/6.0
          pdos(eps,i3,b)=pdos(eps,i3,b)+w3/6.0
          pdos(eps,i4,b)=pdos(eps,i4,b)+w4/6.0
          cpdos(eps,i1,b)=cpdos(eps,i1,b)+cw1/6.0
          cpdos(eps,i2,b)=cpdos(eps,i2,b)+cw2/6.0
          cpdos(eps,i3,b)=cpdos(eps,i3,b)+cw3/6.0
          cpdos(eps,i4,b)=cpdos(eps,i4,b)+cw4/6.0
!          pdos(eps,:)=pdos(eps,:)+w1*ptemp(1,:)+w2*ptemp(2,:)+
!     1                              w3*ptemp(3,:)+w4*ptemp(4,:)
!          cpdos(eps,:)=cpdos(eps,:)+cw1*ptemp(1,:)+cw2*ptemp(2,:)+
!     1                              cw3*ptemp(3,:)+cw4*ptemp(4,:)
        enddo
      enddo
!      pdos(eps,:)=pdos(eps,:)/ntet
!      cpdos(eps,:)=cpdos(eps,:)/ntet
    enddo

  end subroutine


end module tetrahedron  
 
