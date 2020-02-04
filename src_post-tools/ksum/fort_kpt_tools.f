      subroutine get_kptr(nk,nop,nkt,symm,kt,kptr,gptr)
      ! this routine returns a pointer from each kpoint to the ikp...
      ! f2py --link-lapack_opt -c fort_kpt_tools.f -m  fort_kpt_tools
      implicit none
      integer nk(3),symm(nop,3,3),kt(nkt,3),kptr(nkt,3),gptr(nkt)
      integer nop,nkt
      integer g,k,temp(3),tt,i
Cf2py intent(in) nkt,nibz,nop,nk,symm,kt
Cf2py intent(out) kptr,gptr


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
          if (temp(1)+(nk(1)+1)*(temp(2)+(nk(2)+1)*temp(3)) .le.tt)then 
            tt=temp(1)+(nk(1)+1)*(temp(2)+(nk(2)+1)*temp(3))
            !gptr(k)=g-1   ! this is the python convention... not using.
            gptr(k)=g
            kptr(k,:)=temp
          endif
        enddo
        if (gptr(k).eq.0) then
          print*,'Problem in get_kptr...'
          STOP
        endif
      enddo
      
      end             

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
      END 

      subroutine get_itet(nk,ptet,nibz,nkptr,kibz,kptr,ntet,tetptr,tet)
      implicit none
      integer nk(3),ptet(6,4,3),kibz(nibz,3), nibz,tetptr(ntet,4),ntet
      integer kptr(nkptr,3),nkptr,order(4),tet(ntet,4)
      integer i,j,k,l,p,t,temp(3),cc,temp2(4),temp3(4)
Cf2py intent(in) nk,ptet,kibz,nibz,kptr,nkptr,ntet
Cf2py intent(out) tetptr,tet

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
                temp3(t)=1+temp(1)+(nk(1)+1)*(temp(2)+
     1                      (nk(2)+1)*temp(3))
                ! now get the corresponding irreducible kpoint...
                temp=kptr(temp3(t),:)
                ! now find the INDEX of that irreducible kpoint...
                do l=1,nibz
                  if (temp(1) .eq. kibz(l,1) .and. temp(2) .eq. 
     1                kibz(l,2) .and. temp(3) .eq. kibz(l,3))then
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
              ! now reorder to original tet so it corresponds with ikp and store for output...
              tet(cc,1)=temp3(order(1))
              tet(cc,2)=temp3(order(2))
              tet(cc,3)=temp3(order(3))
              tet(cc,4)=temp3(order(4))
              cc=cc+1
            enddo
          enddo
        enddo
      enddo

      end

      subroutine get_tetmult(ntet,nitet,tet,itet,mult)
      implicit none
      integer  ntet,nitet,tet(ntet,4),itet(nitet,4),mult(nitet)
      integer  i,j
Cf2py intent(in) ntet,nitet,tet,itet
Cf2py intent(out) mult
      
      mult=0.d0
      do i=1,nitet
        do j=1,ntet
          if (itet(i,1)==tet(j,1) .and. itet(i,2)==tet(j,2) .and. 
     1        itet(i,3)==tet(j,3) .and. itet(i,4)==tet(j,4) )
     1        mult(i)=mult(i)+1
        enddo
      enddo

      end

      subroutine tetra_pnf(norb,nki,ntet,tet_idx,tetkptr,eigval,cpdos)
      ! f2py --link-lapack_opt -c tempfort.f -m  tempfort
      implicit none
      integer norb,i,nki,t
      real*8 eigval(nki,norb)
      real*8 cpdos(nki,norb)
      real*8 temp(4)
      integer b,ntet,order(4),i1,i2,i3,i4
      integer tetkptr(ntet,4),tet_idx(nki)
      real*8 e1,e2,e3,e4,ef
      real*8 c,c1,c2,c3,cw1,cw2,cw3,cw4
Cf2py intent(in) norb,nki,ntet,tet_idx,tetkptr,eigval
Cf2py intent(out) cpdos

      ! loop over the tetra...
      cpdos=0
      do t=1,ntet
        do b=1,norb
          temp(1)=eigval(tet_idx(tetkptr(t,1)),b)
          temp(2)=eigval(tet_idx(tetkptr(t,2)),b)
          temp(3)=eigval(tet_idx(tetkptr(t,3)),b)
          temp(4)=eigval(tet_idx(tetkptr(t,4)),b)
          call sortr(temp,order,4)
          i1=tet_idx(tetkptr(t,order(1)))
          i2=tet_idx(tetkptr(t,order(2)))
          i3=tet_idx(tetkptr(t,order(3)))
          i4=tet_idx(tetkptr(t,order(4)))
          e1=temp(1);e2=temp(2);e3=temp(3);e4=temp(4)
          ! now compute the weights...
          cw1=0;cw2=0;cw3=0;cw4=0
          if (ef>e1 .and. ef<e2) then
            c=0.25*(ef-e1)**3/((e2-e1)*(e3-e1)*(e4-e1))
            cw1=c*(4-(ef-e1)*(1/(e2-e1)+1/(e3-e1)+1/(e4-e1)))
            cw2=c*(ef-e1)/(e2-e1)
            cw3=c*(ef-e1)/(e3-e1)
            cw4=c*(ef-e1)/(e4-e1)
          elseif (ef>e2 .and. ef<e3) then
            c1=0.25*(ef-e1)**2/((e4-e1)*(e3-e1))
            c2=0.25*((ef-e1)*(ef-e2)*(e3-ef))/((e4-e1)*(e3-e2)*
     1              (e3-e1))
            c3=0.25*(ef-e2)**2*(e4-ef)/((e4-e2)*(e3-e2)*(e4-e1))
            cw1=c1+(c1+c2)*(e3-ef)/(e3-e1)+(c1+c2+c3)*(e4-ef)/(e4-e1)
            cw2=c1+c2+c3+(c2+c3)*(e3-ef)/(e3-e2)+c3*(e4-ef)/(e4-e2)
            cw3=(c1+c2)*(ef-e1)/(e3-e1)+(c2+c3)*(ef-e2)/(e3-e2)
            cw4=(c1+c2+c3)*(ef-e1)/(e4-e1)+c3*(ef-e2)/(e4-e2)
          elseif (ef>e3 .and. ef<e4) then
            c=0.25*(e4-ef)**3/((e4-e1)*(e4-e2)*(e4-e3))
            cw1=0.25-c*(e4-ef)/(e4-e1)
            cw2=0.25-c*(e4-ef)/(e4-e2)
            cw3=0.25-c*(e4-ef)/(e4-e3)
            cw4=0.25-c*(4-(1/(e4-e1)+1/(e4-e2)+1/(e4-e3))*(e4-ef))
          elseif (ef>e4) then
            cw1=0.25;cw2=0.25;cw3=0.25;cw4=0.25
          endif

          cpdos(i1,b)=cpdos(i1,b)+cw1/6
          cpdos(i2,b)=cpdos(i2,b)+cw2/6
          cpdos(i3,b)=cpdos(i3,b)+cw3/6
          cpdos(i4,b)=cpdos(i4,b)+cw4/6
        enddo
      enddo

      end

      SUBROUTINE sortr(X,iy,N)
      ! sorts a real list in ascending order...
      IMPLICIT NONE
      INTEGER N, I, J, iy(n)
      REAL*8 X(N)

      do i=1,n;iy(i)=i;enddo
      DO 100 I=2,N
         IF ( X(I).LT.X(I-1) ) THEN
            DO 50 J=I-2,1,-1
              IF(X(I).GT.X(J)) go to 70
  50          CONTINUE
            J=0
  70        x(j+1:i) = cshift(x(j+1:i),-1)
            iy(j+1:i) = cshift(iy(j+1:i),-1)
         ENDIF
  100 CONTINUE
      RETURN
      END
