module OrderParameter

use omp_lib
use Setup
use Geometry
use TightBinding
use lapack_routines
implicit none

contains

!!!!!!!!!!!!!!!

subroutine NearestNeighbours(NearestNeighborInd,NearestNeighborCell,Coords,ndim,tn)

    integer(dp), intent(in)    :: ndim
    integer(dp), intent(inout) :: NearestNeighborInd(ndim,3),NearestNeighborCell(ndim,3)
    real(dp), intent(in)    :: tn(6,2), Coords(ndim,3)
    
    
    integer(dp) :: i,j,n,ncount,ntemp,ntemp1,ntemp2,ntempt,ntempt1,ntempt2
    real(dp) :: rij,x,y,phi,temp,temp1,temp2,tempt,tempt1,tempt2,rot(ndim,3) 
    
    
    do i=1,ndim
        ncount=0
        do j=1,ndim
    
            if((Coords(i,3)-Coords(j,3))**2 .LT. deltaR)then
                
    
                rij=(Coords(i,1)-Coords(j,1))**2+(Coords(i,2)-Coords(j,2))**2
    
                rij=sqrt(rij)
    
                if(rij.GT.deltaR.AND.rij.LT.a0*1.2)then
                    ncount=ncount+1
                    NearestNeighborInd(i,ncount)=j
                    NearestNeighborCell(i,ncount)=0
    
                    x=Coords(i,1)-Coords(j,1)
                    y=Coords(i,2)-Coords(j,2)
    
                    if(x.GT.0)then
                        phi=atan(y/x)
                    elseif(x.LT.0.AND.y.GT.0)then
                        phi=atan(y/x)+pi
                    else
                        phi=-pi+atan(y/x)   
                    endif
                    if(x.EQ.0)then
                        if(y.GT.0)then
                            phi=pi/2.
                        else
                            phi=-pi/2.
                        endif
                    endif
                    ! write(*,*) x,y,phi*180/pi
    
                    rot(i,ncount)=phi
                    ! write(*,*) i,ncount
                    ! write(*,*) Coords(i,3), Coords(j,3)
                    ! write(*,*) rij
                endif
    
                !!!!!!! Coupling to adjacent Wigner-Seitz cell
    
                do n=1,6
    
                    rij=(Coords(i,1)-Coords(j,1)+tn(n,1))**2+(Coords(i,2)-Coords(j,2)+tn(n,2))**2
    
                    rij=sqrt(rij)
    
                    if(rij.LT.a0*1.2)then
                        ncount=ncount+1
                        NearestNeighborInd(i,ncount)=j
                        NearestNeighborCell(i,ncount)=n
    
                        x=Coords(i,1)-Coords(j,1)+tn(n,1)
                        y=Coords(i,2)-Coords(j,2)+tn(n,2)

                        if(x.GT.0)then
                            phi=atan(y/x)
                        elseif(x.LT.0.AND.y.GT.0)then
                            phi=atan(y/x)+pi
                        else
                            phi=-pi+atan(y/x)   
                        endif
                        if(x.EQ.0)then
                            if(y.GT.0)then
                                phi=pi/2.
                            else
                                phi=-pi/2.
                            endif
                        endif
    
    
                        rot(i,ncount)=phi
                        ! write(*,*) i,ncount
                        ! write(*,*) Coords(i,3), Coords(j,3)
                        ! write(*,*) rij
                    endif
    
                enddo
    
            endif
    
        enddo
    
        if(ncount.NE.3)then
            write(*,*) 'ERROR NEARESTNEIGHBOUR', i,ncount
        endif
    
        !!!!!!!!!!!!!!!!!!!
        !!!!!! Order cycle
        !!!!!!!!!!!!!!!!!!!
    
    
        if(rot(i,1).LT.rot(i,2))then
            temp=rot(i,1)
            ntemp=NearestNeighborInd(i,1)
            ntempt=NearestNeighborCell(i,1)
            rot(i,1)=rot(i,2)
            NearestNeighborInd(i,1)=NearestNeighborInd(i,2)
            NearestNeighborCell(i,1)=NearestNeighborCell(i,2)
            rot(i,2)=temp
            NearestNeighborInd(i,2)=ntemp
            NearestNeighborCell(i,2)=ntempt
        endif
    
        if(rot(i,1).LT.rot(i,3))then
            temp1=rot(i,1)
            ntemp1=NearestNeighborInd(i,1)
            ntempt1=NearestNeighborCell(i,1)
            temp2=rot(i,2)
            ntemp2=NearestNeighborInd(i,2)
            ntempt2=NearestNeighborCell(i,2)
    
            rot(i,1)=rot(i,3)
            NearestNeighborInd(i,1)=NearestNeighborInd(i,3)
            NearestNeighborCell(i,1)=NearestNeighborCell(i,3)
            rot(i,2)=temp1
            NearestNeighborInd(i,2)=ntemp1
            NearestNeighborCell(i,2)=ntempt1
            rot(i,3)=temp2
            NearestNeighborInd(i,3)=ntemp2
            NearestNeighborCell(i,3)=ntempt2
            else
            if(rot(i,2).LT.rot(i,3))then
                temp=rot(i,2)
                ntemp=NearestNeighborInd(i,2)
                ntempt=NearestNeighborCell(i,2)
                rot(i,2)=rot(i,3)
                NearestNeighborInd(i,2)=NearestNeighborInd(i,3)
                NearestNeighborCell(i,2)=NearestNeighborCell(i,3)
                rot(i,3)=temp
                NearestNeighborInd(i,3)=ntemp
                NearestNeighborCell(i,3)=ntempt
            endif
        endif
    
    
    enddo
     
    
    return
end subroutine NearestNeighbours

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine HexagonalLoops(KekuleLoops,Coords,ndim,a1,a2,RotMatrix,tn)

    integer(dp),  intent(in) :: ndim
    integer(dp), intent(out) :: KekuleLoops(ndim/2,2,3)
    real(dp), intent(in) :: Coords(ndim,3), RotMatrix(2,2), a1(2), a2(2), tn(6,2)

    integer(dp) :: i,n,m, counter1, counter2
    real(dp) :: a11(2), a21(2), a12(2), a22(2), r0(2),r1(2),r2(2),r3(2),r4(2),r5(2),r6(2),tm(0:6,2)

    tm(1:6,:) = tn(:,:)
    tm(0,:) = [0.0_dp, 0.0_dp]
    a11 = matmul(RotMatrix,a1)
    a21 = matmul(RotMatrix,a2)
    a12 = matmul(transpose(RotMatrix),a1)
    a22 = matmul(transpose(RotMatrix),a2)

    KekuleLoops(:,:,:) = 0_dp

    do n=1,ndim/4
        r0 = Coords(n,1:2)

        r1 = r0 + a11 - a21
        r2 = r0 + a21
        r3 = r0 - a11
        r4 = r0 + a11
        r5 = r0 + a21 - a11
        r6 = r0 - a21

        counter1 = 0_dp
        counter2 = 0_dp

        do i=1,ndim/4

            do m=0,6
                if(norm2(Coords(i,1:2) + tm(m,:) - r1).lt.a0*.1_dp) then
                    counter1 = counter1 + 1_dp
                    KekuleLoops(n,1,counter1) = i
                endif
            enddo

            do m=0,6
                if(norm2(Coords(i,1:2) + tm(m,:) - r2).lt.a0*.1_dp) then
                    counter1 = counter1 + 1_dp
                    KekuleLoops(n,1,counter1) = i
                endif
            enddo

            do m=0,6
                if(norm2(Coords(i,1:2)+ tm(m,:) - r3).lt.a0*.1_dp) then
                    counter1 = counter1 + 1_dp
                    KekuleLoops(n,1,counter1) = i
                endif
            enddo

            do m=0,6
                if(norm2(Coords(i,1:2)+ tm(m,:) - r4).lt.a0*.1_dp) then
                    counter2 = counter2 + 1_dp
                    KekuleLoops(n,2,counter2) = i
                endif
            enddo

            do m=0,6
                if(norm2(Coords(i,1:2)+ tm(m,:) - r5).lt.a0*.1_dp) then
                    counter2 = counter2 + 1_dp
                    KekuleLoops(n,2,counter2) = i
                endif
            enddo

            do m=0,6
                if(norm2(Coords(i,1:2)+ tm(m,:) - r6).lt.a0*.1_dp) then
                    counter2 = counter2 + 1_dp
                    KekuleLoops(n,2,counter2) = i
                endif
            enddo

        enddo
    enddo

    do n=ndim/2+1,3*ndim/4
        r0 = Coords(n,1:2)

        r1 = r0 + a12 - a22
        r2 = r0 + a22
        r3 = r0 - a12
        r4 = r0 + a12
        r5 = r0 + a22 - a12
        r6 = r0 - a22

        counter1 = 0_dp
        counter2 = 0_dp

        do i=ndim/2+1,3*ndim/4

            do m=0,6
            if(norm2(Coords(i,1:2)+tm(m,:) - r1).lt.a0*.1_dp) then
                counter1 = counter1 + 1_dp
                KekuleLoops(n-ndim/4,1,counter1) = i
            endif
            enddo
            
            do m=0,6
            if(norm2(Coords(i,1:2)+tm(m,:) - r2).lt.a0*.1_dp) then
                counter1 = counter1 + 1_dp
                KekuleLoops(n-ndim/4,1,counter1) = i
            endif
            enddo

            do m=0,6
            if(norm2(Coords(i,1:2)+tm(m,:) - r3).lt.a0*.1_dp) then
                counter1 = counter1 + 1_dp
                KekuleLoops(n-ndim/4,1,counter1) = i
            endif
            enddo

            do m=0,6
            if(norm2(Coords(i,1:2)+tm(m,:) - r4).lt.a0*.1_dp) then
                counter2 = counter2 + 1_dp
                KekuleLoops(n-ndim/4,2,counter2) = i
            endif
            enddo

            do m=0,6
            if(norm2(Coords(i,1:2)+tm(m,:) - r5).lt.a0*.1_dp) then
                counter2 = counter2 + 1_dp
                KekuleLoops(n-ndim/4,2,counter2) = i
            endif
            enddo

            do m=0,6
            if(norm2(Coords(i,1:2)+tm(m,:) - r6).lt.a0*.1_dp) then
                counter2 = counter2 + 1_dp
                KekuleLoops(n-ndim/4,2,counter2) = i
            endif
            enddo

        enddo
    enddo

end subroutine HexagonalLoops

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine KekuleNeighbors(Coords,ndim,RotMatrix,a1,a2,tnn,NeighborsKekule)

    integer(dp), intent(in) :: ndim
    real(dp), intent(in) :: Coords(ndim,3), RotMatrix(2,2), a1(2), a2(2), tnn(6,2)
    integer(dp), intent(inout) :: NeighborsKekule(ndim/2,6,2)

    integer(dp) :: n,h, i
    real(dp) :: a11(2), a21(2), a12(2), a22(2), nCoords(2), tn(0:6,2)

    ndim = ndim
    tn(1:6,:) = tnn(:,:)
    tn(0,1) = 0.0_dp
    tn(0,2) = 0.0_dp

    NeighborsKekule(:,:,:) = 0_dp

    a11 = matmul(RotMatrix,a1)
    a21 = matmul(RotMatrix,a2)
    a12 = matmul(transpose(RotMatrix),a1)
    a22 = matmul(transpose(RotMatrix),a2)

    do n=1,ndim/4

        NeighborsKekule(n,1,1) = n
        NeighborsKekule(n,1,2) = 0

        nCoords = Coords(n,1:2) + a11/3.0 - 2.0*a21/3.0
        do h = ndim/4+1,ndim/2
            do i=0,6
                if (norm2(Coords(h,1:2) + tn(i,:) - nCoords).lt.a0*.1) then
                    
                    NeighborsKekule(n,2,1) = h
                    NeighborsKekule(n,2,2) = i
                endif
            enddo
        enddo
        
        nCoords = Coords(n,1:2) + a11 - a21
        do h = 1,ndim/4
            do i=0,6
                if (norm2(Coords(h,1:2) + tn(i,:) - nCoords).lt.a0*.1) then

                    NeighborsKekule(n,3,1) = h
                    NeighborsKekule(n,3,2) = i
                endif
            enddo
        enddo
        
        nCoords = Coords(n,1:2) + 4*a11/3.0 - 2*a21/3.0
        do h=ndim/4+1,ndim/2
            do i=0,6
                if (norm2(Coords(h,1:2) + tn(i,:) - nCoords).lt.a0*.1) then

                    NeighborsKekule(n,4,1) = h
                    NeighborsKekule(n,4,2) = i
                endif
            enddo
        enddo
        
        nCoords = Coords(n,1:2) + a11
        do h=1,ndim/4
            do i=0,6
                if (norm2(Coords(h,1:2) + tn(i,:) - nCoords).lt.a0*.1) then

                    NeighborsKekule(n,5,1) = h
                    NeighborsKekule(n,5,2) = i
                endif
            enddo
        enddo
        
        nCoords = Coords(n,1:2) + a11/3.0 + a21/3.0
        do h=ndim/4+1,ndim/2
            do i=0,6
                if (norm2(Coords(h,1:2) + tn(i,:) - nCoords).lt.a0*.1) then

                    NeighborsKekule(n,6,1) = h
                    NeighborsKekule(n,6,2) = i
                endif
            enddo
        enddo

        if(NeighborsKekule(n,2,1)*NeighborsKekule(n,3,1)*NeighborsKekule(n,4,1)*NeighborsKekule(n,5,1)*NeighborsKekule(n,6,1).eq.0)then
            NeighborsKekule(n,2:6,:) = 0_dp
            write(*,*) 'ERROR KEKULELATTICE'
        endif


    enddo

    do n=ndim/2+1,3*ndim/4

        NeighborsKekule(n-ndim/4,1,1) = n
        NeighborsKekule(n-ndim/4,1,2) = 0


        nCoords = Coords(n,1:2) + a12/3.0 - 2*a22/3.0
        do h = 3*ndim/4+1,ndim
            do i=0,6
                if (norm2(Coords(h,1:2) + tn(i,:) - nCoords).lt.a0*.1) then
                    NeighborsKekule(n-ndim/4,2,1) = h
                    NeighborsKekule(n-ndim/4,2,2) = i
                endif
            enddo
        enddo
        
        nCoords = Coords(n,1:2) + a12 - a22
        do h=ndim/2+1,3*ndim/4
            do i=0,6
                if (norm2(Coords(h,1:2) + tn(i,:) - nCoords).lt.a0*.1) then
                    NeighborsKekule(n-ndim/4,3,1) = h
                    NeighborsKekule(n-ndim/4,3,2) = i
                endif
            enddo
        enddo
        
        nCoords = Coords(n,1:2) + 4.0*a12/3.0 - 2.0*a22/3.0
        do h=3*ndim/4+1,ndim
            do i=0,6
                if (norm2(Coords(h,1:2) + tn(i,:) - nCoords).lt.a0*.1) then
                    NeighborsKekule(n-ndim/4,4,1) = h
                    NeighborsKekule(n-ndim/4,4,2) = i
                endif
            enddo
        enddo
        
        nCoords = Coords(n,1:2) + a12
        do h=ndim/2+1,3*ndim/4
            do i=0,6
                if (norm2(Coords(h,1:2) + tn(i,:) - nCoords).lt.a0*.1) then
                    NeighborsKekule(n-ndim/4,5,1) = h
                    NeighborsKekule(n-ndim/4,5,2) = i
                endif
            enddo
        enddo
        
        nCoords = Coords(n,1:2) + a12/3.0 + a22/3.0
        do h=3*ndim/4+1,ndim
            do i=0,6
                if (norm2(Coords(h,1:2) + tn(i,:) - nCoords).lt.a0*.1) then
                    NeighborsKekule(n-ndim/4,6,1) = h
                    NeighborsKekule(n-ndim/4,6,2) = i
                endif
            enddo
        enddo

        if(NeighborsKekule(n-ndim/4,2,1)*NeighborsKekule(n-ndim/4,3,1)* &
            NeighborsKekule(n-ndim/4,4,1)*NeighborsKekule(n-ndim/4,5,1)*NeighborsKekule(n-ndim/4,6,1).eq.0)then
            NeighborsKekule(n-ndim/4,2:6,:) = 0
            write(*,*) 'ERROR KEKULELATTICE'
        endif

    enddo

end subroutine KekuleNeighbors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine IntraSubIntraVal(densvp,densnorm,ndim,NearestNeighborInd,NearestNeighborCell,numNeighborCells,nUnitCell_1,nUnitCell_2,tnton1n2,zFock)

    integer(dp), intent(in) :: ndim, numNeighborCells
    integer(dp), intent(in) :: NearestNeighborInd(ndim,3),NearestNeighborCell(ndim,3), nUnitCell_1(numNeighborCells), nUnitCell_2(numNeighborCells), tnton1n2(0:6,2)
    real(dp), intent(out) :: densvp(ndim), densnorm(ndim)
    complex(dp), intent(in) :: zFock(ndim,ndim,numNeighborCells)
    
    integer(dp) :: n,m, index
    real(dp) :: phi0
    real(dp) :: hA1,hB1,hA2,hB2
    complex(dp) :: ztemp,z1,z2,z3, densChern(ndim)
    
    densvp(:) = 0.0_dp
    densnorm(:) = 0.0_dp
    densChern(:) = cmplx(0.0_dp,0.0_dp,dp)

    do n=1,ndim

        if(numI.eq.0)then
            z1 =  zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),1)
            z2 = zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),1)
            z3 = zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),1)
        
        else

            call n1n2toind(index, -tnton1n2(NearestNeighborCell(n,1),1) + tnton1n2(NearestNeighborCell(n,3),1),&
            -tnton1n2(NearestNeighborCell(n,1),2) + tnton1n2(NearestNeighborCell(n,3),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (index.eq.-1_dp)then
                call n1n2toind(index, +tnton1n2(NearestNeighborCell(n,1),1) - tnton1n2(NearestNeighborCell(n,3),1),&
                +tnton1n2(NearestNeighborCell(n,1),2) - tnton1n2(NearestNeighborCell(n,3),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z1 = conjg(zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,1),index))
            else
                
                z1 = zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),index)
            endif
            
            ! if((n.eq.3*ndim/4).or.(n.eq.ndim))then
            ! write(*,*) n,index,z1,zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),1)
            ! endif

            call n1n2toind(index, -tnton1n2(NearestNeighborCell(n,3),1) + tnton1n2(NearestNeighborCell(n,2),1),&
            -tnton1n2(NearestNeighborCell(n,3),2) + tnton1n2(NearestNeighborCell(n,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (index.eq.-1_dp)then
                call n1n2toind(index, +tnton1n2(NearestNeighborCell(n,3),1) - tnton1n2(NearestNeighborCell(n,2),1),&
                +tnton1n2(NearestNeighborCell(n,3),2) - tnton1n2(NearestNeighborCell(n,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z2 = conjg(zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,3),index))
            else
                
                z2 = zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),index)
            endif
            !write(*,*) n,index,z2,zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),1)
            ! if((n.eq.3*ndim/4).or.(n.eq.ndim))then
            ! write(*,*) n,index,z2,zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),1)
            ! endif
            
            call n1n2toind(index, -tnton1n2(NearestNeighborCell(n,2),1) + tnton1n2(NearestNeighborCell(n,1),1),&
            -tnton1n2(NearestNeighborCell(n,2),2) + tnton1n2(NearestNeighborCell(n,1),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (index.eq.-1_dp)then
                call n1n2toind(index, +tnton1n2(NearestNeighborCell(n,2),1) - tnton1n2(NearestNeighborCell(n,1),1),&
                +tnton1n2(NearestNeighborCell(n,2),2) - tnton1n2(NearestNeighborCell(n,1),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z3 = conjg(zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,2),index))
            else
                
                z3 = zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),index)
            endif
            !write(*,*) n,index,z3,zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),1)
            ! if((n.eq.3*ndim/4).or.(n.eq.ndim))then
            ! write(*,*) n,index,z3,zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),1),NearestNeighborInd(n,2),NearestNeighborInd(n,1),NearestNeighborCell(n,2),NearestNeighborCell(n,1)
            ! endif

        endif

        ztemp = z1 + z2 + z3
        densChern(n) = ztemp

    enddo
    
    do n=1,ndim/4
        m=n
        densvp(m) = 2.0_dp/3.0_dp/sqrt(3.0_dp)*aimag(densChern(m))
        densnorm(m) = -2.0_dp/3.0_dp*real(densChern(m))

        m=n+ndim/4
        densvp(m) = -2.0_dp/3.0_dp/sqrt(3.0_dp)*aimag(densChern(m))
        densnorm(m) = -2.0_dp/3.0_dp*real(densChern(m))

        m=n+ndim/2
        densvp(m) = 2.0_dp/3.0_dp/sqrt(3.0_dp)*aimag(densChern(m))
        densnorm(m) = -2.0_dp/3.0_dp*real(densChern(m))

        m=n+3*ndim/4
        densvp(m) = -2.0_dp/3.0_dp/sqrt(3.0_dp)*aimag(densChern(m))
        densnorm(m) = -2.0_dp/3.0_dp*real(densChern(m))

    enddo


    return
end subroutine IntraSubIntraVal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine IntraSubInterVal(fkpk,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,NearestNeighborInd,NearestNeighborCell,tnton1n2,zFock)

    integer(dp), intent(in) :: ndim, numNeighborCells
    integer(dp), intent(in) :: NearestNeighborInd(ndim,3),NearestNeighborCell(ndim,3), tnton1n2(0:6,2), nUnitCell_1(numNeighborCells), nUnitCell_2(numNeighborCells)
    complex(dp), intent(out) :: fkpk(ndim)
    complex(dp), intent(in) :: zFock(ndim,ndim,numNeighborCells)
    
    integer(dp) :: n,m, index
    complex(dp) :: ztemp,z1,z2,z3
    complex(dp) :: zi,zPi3,z2Pi3
    
    zi=cmplx(0._dp,1._dp,dp)
    zPi3=exp(cmplx(0._dp,1._dp,dp)*pi/3._dp)
    z2Pi3=exp(cmplx(0.0_dp,2._dp*pi/3._dp,dp))

    fkpk(:) = cmplx(0.0_dp,0.0_dp,dp)

    do n=1,ndim/4

        if(numI.EQ.0)then
            z1 = zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),1)*z2Pi3
            z2 = zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),1)/z2Pi3
            z3 = zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),1)
            ztemp= z1 + z2 + z3
        else

            call n1n2toind(index, -tnton1n2(NearestNeighborCell(n,1),1) + tnton1n2(NearestNeighborCell(n,3),1),&
            -tnton1n2(NearestNeighborCell(n,1),2) + tnton1n2(NearestNeighborCell(n,3),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (index.eq.-1_dp)then
                call n1n2toind(index, +tnton1n2(NearestNeighborCell(n,1),1) - tnton1n2(NearestNeighborCell(n,3),1),&
                +tnton1n2(NearestNeighborCell(n,1),2) - tnton1n2(NearestNeighborCell(n,3),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z1 = conjg(zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,1),index))
            else
                
                z1 = zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),index)
            endif
            !write(*,*) n,index,z1,zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),1)

            call n1n2toind(index, -tnton1n2(NearestNeighborCell(n,3),1) + tnton1n2(NearestNeighborCell(n,2),1),&
            -tnton1n2(NearestNeighborCell(n,3),2) + tnton1n2(NearestNeighborCell(n,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (index.eq.-1_dp)then
                call n1n2toind(index, +tnton1n2(NearestNeighborCell(n,3),1) - tnton1n2(NearestNeighborCell(n,2),1),&
                +tnton1n2(NearestNeighborCell(n,3),2) - tnton1n2(NearestNeighborCell(n,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z2 = conjg(zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,3),index))
            else
                
                z2 = zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),index)
            endif
            !write(*,*) n,index,z2,zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),1)

            call n1n2toind(index, -tnton1n2(NearestNeighborCell(n,2),1) + tnton1n2(NearestNeighborCell(n,1),1),&
            -tnton1n2(NearestNeighborCell(n,2),2) + tnton1n2(NearestNeighborCell(n,1),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (index.eq.-1_dp)then
                call n1n2toind(index, +tnton1n2(NearestNeighborCell(n,2),1) - tnton1n2(NearestNeighborCell(n,1),1),&
                +tnton1n2(NearestNeighborCell(n,2),2) - tnton1n2(NearestNeighborCell(n,1),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z3 = conjg(zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,2),index))
            else
                
                z3 = zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),index)
            endif
            !write(*,*) n,index,z3,zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),1)

            ztemp = z1*z2Pi3 + z2/z2Pi3 + z3
        endif

        fkpk(n) = ztemp/3.0_dp

    enddo

    do n=ndim/4+1,ndim/2

        if(numI.eq.0)then
            z1 = zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),1)*z2Pi3
            z2 = zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),1)
            z3 = zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),1)/z2Pi3
            ztemp= z1 + z2 + z3
        else

                call n1n2toind(index, -tnton1n2(NearestNeighborCell(n,1),1) + tnton1n2(NearestNeighborCell(n,3),1),&
            -tnton1n2(NearestNeighborCell(n,1),2) + tnton1n2(NearestNeighborCell(n,3),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (index.eq.-1_dp)then
                call n1n2toind(index, +tnton1n2(NearestNeighborCell(n,1),1) - tnton1n2(NearestNeighborCell(n,3),1),&
                +tnton1n2(NearestNeighborCell(n,1),2) - tnton1n2(NearestNeighborCell(n,3),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z1 = conjg(zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,1),index))
            else
                
                z1 = zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),index)
            endif
            !write(*,*) n,index,z1,zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),1)

            call n1n2toind(index, -tnton1n2(NearestNeighborCell(n,3),1) + tnton1n2(NearestNeighborCell(n,2),1),&
            -tnton1n2(NearestNeighborCell(n,3),2) + tnton1n2(NearestNeighborCell(n,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (index.eq.-1_dp)then
                call n1n2toind(index, +tnton1n2(NearestNeighborCell(n,3),1) - tnton1n2(NearestNeighborCell(n,2),1),&
                +tnton1n2(NearestNeighborCell(n,3),2) - tnton1n2(NearestNeighborCell(n,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z2 = conjg(zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,3),index))
            else
                
                z2 = zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),index)
            endif
            !write(*,*) n,index,z2,zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),1)

            call n1n2toind(index, -tnton1n2(NearestNeighborCell(n,2),1) + tnton1n2(NearestNeighborCell(n,1),1),&
            -tnton1n2(NearestNeighborCell(n,2),2) + tnton1n2(NearestNeighborCell(n,1),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (index.eq.-1_dp)then
                call n1n2toind(index, +tnton1n2(NearestNeighborCell(n,2),1) - tnton1n2(NearestNeighborCell(n,1),1),&
                +tnton1n2(NearestNeighborCell(n,2),2) - tnton1n2(NearestNeighborCell(n,1),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z3 = conjg(zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,2),index))
            else
                
                z3 = zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),index)
            endif
            !write(*,*) n,index,z3,zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),1)

            ztemp = z1*z2Pi3 + z2 + z3/z2Pi3

        endif

        fkpk(n) = ztemp/3.0_dp

    enddo

    do n=ndim/2+1,3*ndim/4

        if(numI.eq.0)then
            z1 = zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),1)*z2Pi3
            z2 = zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),1)/z2Pi3
            z3 = zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),1)
            ztemp = z1 + z2 + z3

        else
            call n1n2toind(index, -tnton1n2(NearestNeighborCell(n,1),1) + tnton1n2(NearestNeighborCell(n,3),1),&
            -tnton1n2(NearestNeighborCell(n,1),2) + tnton1n2(NearestNeighborCell(n,3),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (index.eq.-1_dp)then
                call n1n2toind(index, +tnton1n2(NearestNeighborCell(n,1),1) - tnton1n2(NearestNeighborCell(n,3),1),&
                +tnton1n2(NearestNeighborCell(n,1),2) - tnton1n2(NearestNeighborCell(n,3),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z1 = conjg(zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,1),index))
            else
                
                z1 = zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),index)
            endif
            !write(*,*) n,index,z1,zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),1)

            call n1n2toind(index, -tnton1n2(NearestNeighborCell(n,3),1) + tnton1n2(NearestNeighborCell(n,2),1),&
            -tnton1n2(NearestNeighborCell(n,3),2) + tnton1n2(NearestNeighborCell(n,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (index.eq.-1_dp)then
                call n1n2toind(index, +tnton1n2(NearestNeighborCell(n,3),1) - tnton1n2(NearestNeighborCell(n,2),1),&
                +tnton1n2(NearestNeighborCell(n,3),2) - tnton1n2(NearestNeighborCell(n,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z2 = conjg(zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,3),index))
            else
                
                z2 = zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),index)
            endif
            !write(*,*) n,index,z2,zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),1)

            call n1n2toind(index, -tnton1n2(NearestNeighborCell(n,2),1) + tnton1n2(NearestNeighborCell(n,1),1),&
            -tnton1n2(NearestNeighborCell(n,2),2) + tnton1n2(NearestNeighborCell(n,1),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (index.eq.-1_dp)then
                call n1n2toind(index, +tnton1n2(NearestNeighborCell(n,2),1) - tnton1n2(NearestNeighborCell(n,1),1),&
                +tnton1n2(NearestNeighborCell(n,2),2) - tnton1n2(NearestNeighborCell(n,1),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z3 = conjg(zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,2),index))
            else
                
                z3 = zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),index)
            endif
            !write(*,*) n,index,z3,zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),1)

            ztemp = z1*z2Pi3 + z2/z2Pi3 + z3
        endif
        
        fkpk(n) = ztemp/3.0_dp
    enddo

    do n=3*ndim/4+1,ndim

        if (numI.eq.0) then
            z1 = zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),1)*z2Pi3
            z2 = zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),1)
            z3 = zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),1)/z2Pi3
            ztemp = z1 + z2 + z3
        else

            call n1n2toind(index, -tnton1n2(NearestNeighborCell(n,1),1) + tnton1n2(NearestNeighborCell(n,3),1),&
            -tnton1n2(NearestNeighborCell(n,1),2) + tnton1n2(NearestNeighborCell(n,3),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (index.eq.-1_dp)then
                call n1n2toind(index, +tnton1n2(NearestNeighborCell(n,1),1) - tnton1n2(NearestNeighborCell(n,3),1),&
                +tnton1n2(NearestNeighborCell(n,1),2) - tnton1n2(NearestNeighborCell(n,3),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z1 = conjg(zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,1),index))
            else
                
                z1 = zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),index)
            endif
            !write(*,*) n,index,z1,zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),1)

            call n1n2toind(index, -tnton1n2(NearestNeighborCell(n,3),1) + tnton1n2(NearestNeighborCell(n,2),1),&
            -tnton1n2(NearestNeighborCell(n,3),2) + tnton1n2(NearestNeighborCell(n,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (index.eq.-1_dp)then
                call n1n2toind(index, +tnton1n2(NearestNeighborCell(n,3),1) - tnton1n2(NearestNeighborCell(n,2),1),&
                +tnton1n2(NearestNeighborCell(n,3),2) - tnton1n2(NearestNeighborCell(n,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z2 = conjg(zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,3),index))
            else
                
                z2 = zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),index)
            endif
            !write(*,*) n,index,z2,zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),1)

            call n1n2toind(index, -tnton1n2(NearestNeighborCell(n,2),1) + tnton1n2(NearestNeighborCell(n,1),1),&
            -tnton1n2(NearestNeighborCell(n,2),2) + tnton1n2(NearestNeighborCell(n,1),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (index.eq.-1_dp)then
                call n1n2toind(index, +tnton1n2(NearestNeighborCell(n,2),1) - tnton1n2(NearestNeighborCell(n,1),1),&
                +tnton1n2(NearestNeighborCell(n,2),2) - tnton1n2(NearestNeighborCell(n,1),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z3 = conjg(zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,2),index))
            else
                
                z3 = zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),index)
            endif
            !write(*,*) n,index,z3,zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),1)
            
            ztemp = z1*z2Pi3 + z2 + z3/z2Pi3
        endif

        fkpk(n) = ztemp/3.0_dp
    enddo
        
    return
end subroutine IntraSubInterVal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InterSubInterVal(fkakpb,fkbkpa,NeighborsKekule,KekuleLoops,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,tnton1n2,zFock)

    integer(dp), intent(in) :: ndim, numNeighborCells
    integer(dp), intent(in) :: nUnitCell_1(numNeighborCells), nUnitCell_2(numNeighborCells), tnton1n2(0:6,2)
    integer(dp), intent(in) :: NeighborsKekule(ndim/2,6,2), KekuleLoops(ndim/2,2,3)
    complex(dp), intent(out) :: fkakpb(ndim/2), fkbkpa(ndim/2)
    complex(dp), intent(in) :: zFock(ndim,ndim,numNeighborCells)
    
    integer(dp) :: n, m, index
    complex(dp) :: zi,zPi3,z2Pi3,z1,z2,z3,z4,z5,z6,delta0,delta1,delta2
    complex(dp) :: current(ndim/2)

    zi=cmplx(0.0_dp,1._dp,dp)
    zPi3=exp(cmplx(0.0_dp,1._dp,dp)*pi/3._dp)
    z2Pi3=exp(cmplx(0.0_dp,2._dp*pi/3._dp,dp))

    current(:) = cmplx(0.0_dp,0.0_dp,dp)
    fkakpb(:) = cmplx(0.0_dp,0.0_dp,dp)
    fkbkpa(:) = cmplx(0.0_dp,0.0_dp,dp)

    do n=1,ndim/2

        if(numI.eq.0)then

            z1 = zFock(NeighborsKekule(n,1,1),NeighborsKekule(n,2,1),1)
            z2 = zFock(NeighborsKekule(n,2,1),NeighborsKekule(n,3,1),1)
            z3 = zFock(NeighborsKekule(n,3,1),NeighborsKekule(n,4,1),1)
            z4 = zFock(NeighborsKekule(n,4,1),NeighborsKekule(n,5,1),1)
            z5 = zFock(NeighborsKekule(n,5,1),NeighborsKekule(n,6,1),1)
            z6 = zFock(NeighborsKekule(n,6,1),NeighborsKekule(n,1,1),1)

        else

            call n1n2toind(index, tnton1n2(NeighborsKekule(n,1,2),1) - tnton1n2(NeighborsKekule(n,2,2),1),&
                tnton1n2(NeighborsKekule(n,1,2),2) - tnton1n2(NeighborsKekule(n,2,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(index.eq.-1_dp)then
                call n1n2toind(index, -tnton1n2(NeighborsKekule(n,1,2),1) + tnton1n2(NeighborsKekule(n,2,2),1),&
                -tnton1n2(NeighborsKekule(n,1,2),2) + tnton1n2(NeighborsKekule(n,2,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z1 = conjg(zFock(NeighborsKekule(n,2,1),NeighborsKekule(n,1,1),index))
            else

                z1 = zFock(NeighborsKekule(n,1,1),NeighborsKekule(n,2,1),index)
            endif

            call n1n2toind(index, tnton1n2(NeighborsKekule(n,2,2),1) - tnton1n2(NeighborsKekule(n,3,2),1),&
                tnton1n2(NeighborsKekule(n,2,2),2) - tnton1n2(NeighborsKekule(n,3,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(index.eq.-1_dp)then
                call n1n2toind(index, -tnton1n2(NeighborsKekule(n,2,2),1) + tnton1n2(NeighborsKekule(n,3,2),1),&
                -tnton1n2(NeighborsKekule(n,2,2),2) + tnton1n2(NeighborsKekule(n,3,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z2 = conjg(zFock(NeighborsKekule(n,3,1),NeighborsKekule(n,2,1),index))
            else

                z2 = zFock(NeighborsKekule(n,2,1),NeighborsKekule(n,3,1),index)
            endif
            
            call n1n2toind(index, tnton1n2(NeighborsKekule(n,3,2),1) - tnton1n2(NeighborsKekule(n,4,2),1),&
                tnton1n2(NeighborsKekule(n,3,2),2) - tnton1n2(NeighborsKekule(n,4,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(index.eq.-1_dp)then
                call n1n2toind(index, -tnton1n2(NeighborsKekule(n,3,2),1) + tnton1n2(NeighborsKekule(n,4,2),1),&
                -tnton1n2(NeighborsKekule(n,3,2),2) + tnton1n2(NeighborsKekule(n,4,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z3 = conjg(zFock(NeighborsKekule(n,4,1),NeighborsKekule(n,3,1),index))
            else

                z3 = zFock(NeighborsKekule(n,3,1),NeighborsKekule(n,4,1),index)
            endif
            
            call n1n2toind(index, tnton1n2(NeighborsKekule(n,4,2),1) - tnton1n2(NeighborsKekule(n,5,2),1),&
            tnton1n2(NeighborsKekule(n,4,2),2) - tnton1n2(NeighborsKekule(n,5,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(index.eq.-1_dp)then
                call n1n2toind(index, -tnton1n2(NeighborsKekule(n,4,2),1) + tnton1n2(NeighborsKekule(n,5,2),1),&
                -tnton1n2(NeighborsKekule(n,4,2),2) + tnton1n2(NeighborsKekule(n,5,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z4 = conjg(zFock(NeighborsKekule(n,5,1),NeighborsKekule(n,4,1),index))
            else

                z4 = zFock(NeighborsKekule(n,4,1),NeighborsKekule(n,5,1),index)
            endif
            
            call n1n2toind(index, tnton1n2(NeighborsKekule(n,5,2),1) - tnton1n2(NeighborsKekule(n,6,2),1),&
            tnton1n2(NeighborsKekule(n,5,2),2) - tnton1n2(NeighborsKekule(n,6,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(index.eq.-1_dp)then
                call n1n2toind(index, -tnton1n2(NeighborsKekule(n,5,2),1) + tnton1n2(NeighborsKekule(n,6,2),1),&
                -tnton1n2(NeighborsKekule(n,5,2),2) + tnton1n2(NeighborsKekule(n,6,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z5 = conjg(zFock(NeighborsKekule(n,6,1),NeighborsKekule(n,5,1),index))
            else

                z5 = zFock(NeighborsKekule(n,5,1),NeighborsKekule(n,6,1),index)
            endif
            
            call n1n2toind(index, tnton1n2(NeighborsKekule(n,6,2),1) - tnton1n2(NeighborsKekule(n,1,2),1),&
            tnton1n2(NeighborsKekule(n,6,2),2) - tnton1n2(NeighborsKekule(n,1,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(index.eq.-1_dp)then
                call n1n2toind(index, -tnton1n2(NeighborsKekule(n,6,2),1) + tnton1n2(NeighborsKekule(n,1,2),1),&
                -tnton1n2(NeighborsKekule(n,6,2),2) + tnton1n2(NeighborsKekule(n,1,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z6 = conjg(zFock(NeighborsKekule(n,1,1),NeighborsKekule(n,6,1),index))
            else

                z6 = zFock(NeighborsKekule(n,6,1),NeighborsKekule(n,1,1),index)
            endif

        endif

        current(n) = z1+z2+z3+z4+z5+z6


    enddo

    do n=1,ndim/4

        !if (product(KekuleLoops(n,:,:)).ne.0_dp)then
        
        delta0 = current(n)
        delta1 = (current(KekuleLoops(n,1,1)) + current(KekuleLoops(n,1,2)) + current(KekuleLoops(n,1,3)))/3.0_dp 
        delta2 = (current(KekuleLoops(n,2,1)) + current(KekuleLoops(n,2,2)) + current(KekuleLoops(n,2,3)))/3.0_dp 

        fkakpb(n) = (zPi3*(conjg(delta0) + conjg(delta1)/z2Pi3 + conjg(delta2)*z2Pi3) +&
                    (delta0 + delta1/z2Pi3 + delta2*z2Pi3))*cmplx(0.0_dp,-1.0_dp,dp)/9.0_dp/sqrt(3.0_dp)

        fkbkpa(n) = ((conjg(delta0) + conjg(delta1)/z2Pi3 + conjg(delta2)*z2Pi3) +&
                    zPi3*(delta0 + delta1/z2Pi3 + delta2*z2Pi3))*cmplx(0.0_dp,-1.0_dp,dp)/9.0_dp/sqrt(3.0_dp)
            
        !endif
                
    enddo

    do n=ndim/4+1,ndim/2

        !if (sum(KekuleLoops(n,:,:)).ne.0_dp)then
        
        delta0 = current(n)
        delta1 = (current(KekuleLoops(n,1,1)-ndim/4) + current(KekuleLoops(n,1,2)-ndim/4)&
                    + current(KekuleLoops(n,1,3)-ndim/4))/3.0_dp 
        delta2 = (current(KekuleLoops(n,2,1)-ndim/4) + current(KekuleLoops(n,2,2)-ndim/4)&
                    + current(KekuleLoops(n,2,3)-ndim/4))/3.0_dp 

        fkakpb(n) = (zPi3*(conjg(delta0) + conjg(delta1)/z2Pi3 + conjg(delta2)*z2Pi3) +&
                    (delta0 + delta1/z2Pi3 + delta2*z2Pi3))*cmplx(0.0_dp,-1.0_dp,dp)/9.0_dp/sqrt(3.0_dp)

        fkbkpa(n) = ((conjg(delta0) + conjg(delta1)/z2Pi3 + conjg(delta2)*z2Pi3) +&
                    zPi3*(delta0 + delta1/z2Pi3 + delta2*z2Pi3))*cmplx(0.0_dp,-1.0_dp,dp)/9.0_dp/sqrt(3.0_dp)

            
        !endif
                
    enddo


    return
end subroutine InterSubInterVal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InterSubIntraVal_old(fkakb,fkpakpb,NeighborsKekule,KekuleLoops,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,tnton1n2,zFock)

    integer(dp), intent(in) :: ndim,numNeighborCells
    integer(dp), intent(in) :: nUnitCell_1(numNeighborCells), nUnitCell_2(numNeighborCells)
    integer(dp), intent(in) :: NeighborsKekule(ndim/2,6,2), KekuleLoops(ndim/2,2,3), tnton1n2(0:6,2)
    complex(dp), intent(out) :: fkakb(ndim/2), fkpakpb(ndim/2)
    complex(dp), intent(in) :: zFock(ndim,ndim,numNeighborCells)
    
    integer(dp) :: n,m,index
    complex(dp) :: zi,zPi3,z2Pi3,z1,z2,z3,z4,z5,z6
    complex(dp) :: current(ndim/2,4)
    
    zi=cmplx(0.0_dp,1._dp,dp)
    zPi3=exp(cmplx(0.0_dp,1._dp,dp)*pi/3._dp)
    z2Pi3=exp(cmplx(0.0_dp,2._dp*pi/3._dp,dp))

    current(:,:) = cmplx(0.0_dp,0.0_dp,dp)
    fkakb(:) = cmplx(0.0_dp,0.0_dp,dp)
    fkpakpb(:) = cmplx(0.0_dp,0.0_dp,dp)

    do n=1,ndim/2

        if(numI.eq.0)then

            z1 = zFock(NeighborsKekule(n,1,1),NeighborsKekule(n,2,1),1)
            z2 = zFock(NeighborsKekule(n,2,1),NeighborsKekule(n,3,1),1)
            z3 = zFock(NeighborsKekule(n,3,1),NeighborsKekule(n,4,1),1)
            z4 = zFock(NeighborsKekule(n,4,1),NeighborsKekule(n,5,1),1)
            z5 = zFock(NeighborsKekule(n,5,1),NeighborsKekule(n,6,1),1)
            z6 = zFock(NeighborsKekule(n,6,1),NeighborsKekule(n,1,1),1)
        else

            call n1n2toind(index, tnton1n2(NeighborsKekule(n,1,2),1) - tnton1n2(NeighborsKekule(n,2,2),1),&
                tnton1n2(NeighborsKekule(n,1,2),2) - tnton1n2(NeighborsKekule(n,2,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(index.eq.-1_dp)then
                call n1n2toind(index, -tnton1n2(NeighborsKekule(n,1,2),1) + tnton1n2(NeighborsKekule(n,2,2),1),&
                -tnton1n2(NeighborsKekule(n,1,2),2) + tnton1n2(NeighborsKekule(n,2,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z1 = conjg(zFock(NeighborsKekule(n,2,1),NeighborsKekule(n,1,1),index))
            else

                z1 = zFock(NeighborsKekule(n,1,1),NeighborsKekule(n,2,1),index)
            endif
            ! write(*,*) index

            call n1n2toind(index, tnton1n2(NeighborsKekule(n,2,2),1) - tnton1n2(NeighborsKekule(n,3,2),1),&
                tnton1n2(NeighborsKekule(n,2,2),2) - tnton1n2(NeighborsKekule(n,3,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(index.eq.-1_dp)then
                call n1n2toind(index, -tnton1n2(NeighborsKekule(n,2,2),1) + tnton1n2(NeighborsKekule(n,3,2),1),&
                -tnton1n2(NeighborsKekule(n,2,2),2) + tnton1n2(NeighborsKekule(n,3,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z2 = conjg(zFock(NeighborsKekule(n,3,1),NeighborsKekule(n,2,1),index))
            else

                z2 = zFock(NeighborsKekule(n,2,1),NeighborsKekule(n,3,1),index)
            endif
            ! write(*,*) index

            
            call n1n2toind(index, tnton1n2(NeighborsKekule(n,3,2),1) - tnton1n2(NeighborsKekule(n,4,2),1),&
                tnton1n2(NeighborsKekule(n,3,2),2) - tnton1n2(NeighborsKekule(n,4,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(index.eq.-1_dp)then
                call n1n2toind(index, -tnton1n2(NeighborsKekule(n,3,2),1) + tnton1n2(NeighborsKekule(n,4,2),1),&
                -tnton1n2(NeighborsKekule(n,3,2),2) + tnton1n2(NeighborsKekule(n,4,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z3 = conjg(zFock(NeighborsKekule(n,4,1),NeighborsKekule(n,3,1),index))
            else

                z3 = zFock(NeighborsKekule(n,3,1),NeighborsKekule(n,4,1),index)
            endif
            ! write(*,*) index

            
            call n1n2toind(index, tnton1n2(NeighborsKekule(n,4,2),1) - tnton1n2(NeighborsKekule(n,5,2),1),&
            tnton1n2(NeighborsKekule(n,4,2),2) - tnton1n2(NeighborsKekule(n,5,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(index.eq.-1_dp)then
                call n1n2toind(index, -tnton1n2(NeighborsKekule(n,4,2),1) + tnton1n2(NeighborsKekule(n,5,2),1),&
                -tnton1n2(NeighborsKekule(n,4,2),2) + tnton1n2(NeighborsKekule(n,5,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z4 = conjg(zFock(NeighborsKekule(n,5,1),NeighborsKekule(n,4,1),index))
            else

                z4 = zFock(NeighborsKekule(n,4,1),NeighborsKekule(n,5,1),index)
            endif
            ! write(*,*) index

            
            call n1n2toind(index, tnton1n2(NeighborsKekule(n,5,2),1) - tnton1n2(NeighborsKekule(n,6,2),1),&
            tnton1n2(NeighborsKekule(n,5,2),2) - tnton1n2(NeighborsKekule(n,6,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(index.eq.-1_dp)then
                call n1n2toind(index, -tnton1n2(NeighborsKekule(n,5,2),1) + tnton1n2(NeighborsKekule(n,6,2),1),&
                -tnton1n2(NeighborsKekule(n,5,2),2) + tnton1n2(NeighborsKekule(n,6,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z5 = conjg(zFock(NeighborsKekule(n,6,1),NeighborsKekule(n,5,1),index))
            else

                z5 = zFock(NeighborsKekule(n,5,1),NeighborsKekule(n,6,1),index)
            endif
            ! write(*,*) index

            
            call n1n2toind(index, tnton1n2(NeighborsKekule(n,6,2),1) - tnton1n2(NeighborsKekule(n,1,2),1),&
            tnton1n2(NeighborsKekule(n,6,2),2) - tnton1n2(NeighborsKekule(n,1,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(index.eq.-1_dp)then
                call n1n2toind(index, -tnton1n2(NeighborsKekule(n,6,2),1) + tnton1n2(NeighborsKekule(n,1,2),1),&
                -tnton1n2(NeighborsKekule(n,6,2),2) + tnton1n2(NeighborsKekule(n,1,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z6 = conjg(zFock(NeighborsKekule(n,1,1),NeighborsKekule(n,6,1),index))
            else

                z6 = zFock(NeighborsKekule(n,6,1),NeighborsKekule(n,1,1),index)
            endif
            ! write(*,*) index

        endif

        current(n,1) = z1*z2Pi3 + z2*z2Pi3 + z3 + z4/z2Pi3 + z5/z2Pi3 + z6
        current(n,2) = z1 + z2/z2Pi3 + z3/z2Pi3 + z4 + z5*z2Pi3 + z6*z2Pi3
        current(n,3) = z1/z2Pi3 + z2/z2Pi3 + z3 + z4*z2Pi3 + z5*z2Pi3 + z6
        current(n,4) = z1*z2Pi3 + z2 + z3/z2Pi3 + z4/z2Pi3 + z5 + z6*z2Pi3


    enddo

    do n=1,ndim/2

        ! if (product(KekuleLoops(n,:,:)).ne.0_dp)then

            fkakb(n) = (current(n,1)*z2Pi3 - current(n,2))/3.0_dp/(z2Pi3 - 1/z2Pi3)
            ! fkakb(n) = current(n,1)
            fkpakpb(n) = (current(n,3)*z2Pi3 - current(n,4))/3.0_dp/(z2Pi3 - 1/z2Pi3)
            ! fkpakpb(n) = current(n,2)


        ! endif
                
    enddo

    return
end subroutine InterSubIntraVal_old

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InterSubIntraVal(fkakb,fkpakpb,NeighborsKekule,KekuleLoops,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,tnton1n2,zFock)

    integer(dp), intent(in) :: ndim,numNeighborCells
    integer(dp), intent(in) :: nUnitCell_1(numNeighborCells), nUnitCell_2(numNeighborCells)
    integer(dp), intent(in) :: NeighborsKekule(ndim/2,6,2), KekuleLoops(ndim/2,2,3), tnton1n2(0:6,2)
    complex(dp), intent(out) :: fkakb(ndim/2), fkpakpb(ndim/2)
    complex(dp), intent(in) :: zFock(ndim,ndim,numNeighborCells)
    
    integer(dp) :: n,m,index
    complex(dp) :: zi,zPi3,z2Pi3,z1,z2,z3,z4,z5,z6
    complex(dp) :: current(ndim/2,4)
    
    zi=cmplx(0.0_dp,1._dp,dp)
    zPi3=exp(cmplx(0.0_dp,1._dp,dp)*pi/3._dp)
    z2Pi3=exp(cmplx(0.0_dp,2._dp*pi/3._dp,dp))

    current(:,:) = cmplx(0.0_dp,0.0_dp,dp)
    fkakb(:) = cmplx(0.0_dp,0.0_dp,dp)
    fkpakpb(:) = cmplx(0.0_dp,0.0_dp,dp)

    do n=1,ndim/2

        if(numI.eq.0)then

            z1 = zFock(NeighborsKekule(n,1,1),NeighborsKekule(n,2,1),1)
            z2 = zFock(NeighborsKekule(n,2,1),NeighborsKekule(n,3,1),1)
            z3 = zFock(NeighborsKekule(n,3,1),NeighborsKekule(n,4,1),1)
            z4 = zFock(NeighborsKekule(n,4,1),NeighborsKekule(n,5,1),1)
            z5 = zFock(NeighborsKekule(n,5,1),NeighborsKekule(n,6,1),1)
            z6 = zFock(NeighborsKekule(n,6,1),NeighborsKekule(n,1,1),1)
        else

            call n1n2toind(index, tnton1n2(NeighborsKekule(n,1,2),1) - tnton1n2(NeighborsKekule(n,2,2),1),&
                tnton1n2(NeighborsKekule(n,1,2),2) - tnton1n2(NeighborsKekule(n,2,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(index.eq.-1_dp)then
                call n1n2toind(index, -tnton1n2(NeighborsKekule(n,1,2),1) + tnton1n2(NeighborsKekule(n,2,2),1),&
                -tnton1n2(NeighborsKekule(n,1,2),2) + tnton1n2(NeighborsKekule(n,2,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z1 = conjg(zFock(NeighborsKekule(n,2,1),NeighborsKekule(n,1,1),index))
            else

                z1 = zFock(NeighborsKekule(n,1,1),NeighborsKekule(n,2,1),index)
            endif
            ! write(*,*) index

            call n1n2toind(index, tnton1n2(NeighborsKekule(n,2,2),1) - tnton1n2(NeighborsKekule(n,3,2),1),&
                tnton1n2(NeighborsKekule(n,2,2),2) - tnton1n2(NeighborsKekule(n,3,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(index.eq.-1_dp)then
                call n1n2toind(index, -tnton1n2(NeighborsKekule(n,2,2),1) + tnton1n2(NeighborsKekule(n,3,2),1),&
                -tnton1n2(NeighborsKekule(n,2,2),2) + tnton1n2(NeighborsKekule(n,3,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z2 = conjg(zFock(NeighborsKekule(n,3,1),NeighborsKekule(n,2,1),index))
            else

                z2 = zFock(NeighborsKekule(n,2,1),NeighborsKekule(n,3,1),index)
            endif
            ! write(*,*) index

            
            call n1n2toind(index, tnton1n2(NeighborsKekule(n,3,2),1) - tnton1n2(NeighborsKekule(n,4,2),1),&
                tnton1n2(NeighborsKekule(n,3,2),2) - tnton1n2(NeighborsKekule(n,4,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(index.eq.-1_dp)then
                call n1n2toind(index, -tnton1n2(NeighborsKekule(n,3,2),1) + tnton1n2(NeighborsKekule(n,4,2),1),&
                -tnton1n2(NeighborsKekule(n,3,2),2) + tnton1n2(NeighborsKekule(n,4,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z3 = conjg(zFock(NeighborsKekule(n,4,1),NeighborsKekule(n,3,1),index))
            else

                z3 = zFock(NeighborsKekule(n,3,1),NeighborsKekule(n,4,1),index)
            endif
            ! write(*,*) index

            
            call n1n2toind(index, tnton1n2(NeighborsKekule(n,4,2),1) - tnton1n2(NeighborsKekule(n,5,2),1),&
            tnton1n2(NeighborsKekule(n,4,2),2) - tnton1n2(NeighborsKekule(n,5,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(index.eq.-1_dp)then
                call n1n2toind(index, -tnton1n2(NeighborsKekule(n,4,2),1) + tnton1n2(NeighborsKekule(n,5,2),1),&
                -tnton1n2(NeighborsKekule(n,4,2),2) + tnton1n2(NeighborsKekule(n,5,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z4 = conjg(zFock(NeighborsKekule(n,5,1),NeighborsKekule(n,4,1),index))
            else

                z4 = zFock(NeighborsKekule(n,4,1),NeighborsKekule(n,5,1),index)
            endif
            ! write(*,*) index

            
            call n1n2toind(index, tnton1n2(NeighborsKekule(n,5,2),1) - tnton1n2(NeighborsKekule(n,6,2),1),&
            tnton1n2(NeighborsKekule(n,5,2),2) - tnton1n2(NeighborsKekule(n,6,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(index.eq.-1_dp)then
                call n1n2toind(index, -tnton1n2(NeighborsKekule(n,5,2),1) + tnton1n2(NeighborsKekule(n,6,2),1),&
                -tnton1n2(NeighborsKekule(n,5,2),2) + tnton1n2(NeighborsKekule(n,6,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z5 = conjg(zFock(NeighborsKekule(n,6,1),NeighborsKekule(n,5,1),index))
            else

                z5 = zFock(NeighborsKekule(n,5,1),NeighborsKekule(n,6,1),index)
            endif
            ! write(*,*) index

            
            call n1n2toind(index, tnton1n2(NeighborsKekule(n,6,2),1) - tnton1n2(NeighborsKekule(n,1,2),1),&
            tnton1n2(NeighborsKekule(n,6,2),2) - tnton1n2(NeighborsKekule(n,1,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(index.eq.-1_dp)then
                call n1n2toind(index, -tnton1n2(NeighborsKekule(n,6,2),1) + tnton1n2(NeighborsKekule(n,1,2),1),&
                -tnton1n2(NeighborsKekule(n,6,2),2) + tnton1n2(NeighborsKekule(n,1,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z6 = conjg(zFock(NeighborsKekule(n,1,1),NeighborsKekule(n,6,1),index))
            else

                z6 = zFock(NeighborsKekule(n,6,1),NeighborsKekule(n,1,1),index)
            endif
            ! write(*,*) index

        endif

        current(n,1) = z1*z2Pi3 + z2/z2Pi3 + z3 + z4*z2Pi3 + z5/z2Pi3 + z6
        current(n,2) = z1/z2Pi3 + z2*z2Pi3 + z3*z2Pi3 + z4 + z5 + z6/z2Pi3

    enddo

    do n=1,ndim/2

        ! if (product(KekuleLoops(n,:,:)).ne.0_dp)then

            fkakb(n) = (current(n,2) - current(n,1)/z2Pi3)/3.0_dp*cmplx(0.0_dp,-1.0_dp,dp)/sqrt(3.0_dp)
            ! fkakb(n) = current(n,1)
            fkpakpb(n) = conjg((current(n,2) - current(n,1)*z2Pi3)/3.0_dp*cmplx(0.0_dp,1.0_dp,dp)/sqrt(3.0_dp))
            ! fkpakpb(n) = current(n,2)


        ! endif
                
    enddo

    return
end subroutine InterSubIntraVal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real function fphi(x,y)

    real(dp), intent(in) :: x,y
    
    if(x.GT.0)then
    fphi=atan(y/x)
    endif
    if(x.LT.0)then
    fphi=atan(y/x)+pi    
    endif
    if(x.EQ.0)then
    if(y.GT.0)then
    fphi=pi/2.
    else
    fphi=-pi/2.
    endif
    endif

end function fphi

end module OrderParameter
