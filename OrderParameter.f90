module OrderPlot

use omp_lib
use Setup
use lapack_routines

implicit none

contains

!!!!!!!!!!!!!!!


  subroutine nUnitCell12ToInd(index,n1,n2,nUnitCell_1,nUnitCell_2,ind)

    integer(dp), intent(in) :: ind
    integer(dp), intent(in) :: n1, n2, nUnitCell_1(ind), nUnitCell_2(ind)
    integer(dp), intent(inout) :: index

    integer(dp) :: n
    
    index = -1_dp 
    do n=1,ind
        if((abs(real(nUnitCell_1(n))-real(n1)).lt.1e-3).and.(abs(real(nUnitCell_2(n))-real(n2)).lt.1e-3))then
            index = n
        endif
    enddo

end subroutine nUnitCell12ToInd



subroutine getKekuleNeighbors(KekuleNeighbors,Coords,ndim,a1,a2,RotMatrix,tn)

    integer(dp),  intent(in) :: ndim
    integer(dp), intent(out) :: KekuleNeighbors(ndim/2,2,3)
    real(dp), intent(in) :: Coords(ndim,3), RotMatrix(2,2), a1(2), a2(2), tn(6,2)

    integer(dp) :: i,n,m, counter1, counter2
    real(dp) :: a11(2), a21(2), a12(2), a22(2), r0(2),r1(2),r2(2),r3(2),r4(2),r5(2),r6(2),tm(0:6,2)

    tm(1:6,:) = tn(:,:)
    tm(0,:) = [0.0_dp, 0.0_dp]
    a11 = matmul(RotMatrix,a1)
    a21 = matmul(RotMatrix,a2)
    a12 = matmul(transpose(RotMatrix),a1)
    a22 = matmul(transpose(RotMatrix),a2)

    KekuleNeighbors(:,:,:) = 0_dp

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
                    KekuleNeighbors(n,1,counter1) = i
                endif
            enddo

            do m=0,6
                if(norm2(Coords(i,1:2) + tm(m,:) - r2).lt.a0*.1_dp) then
                    counter1 = counter1 + 1_dp
                    KekuleNeighbors(n,1,counter1) = i
                endif
            enddo

            do m=0,6
                if(norm2(Coords(i,1:2)+ tm(m,:) - r3).lt.a0*.1_dp) then
                    counter1 = counter1 + 1_dp
                    KekuleNeighbors(n,1,counter1) = i
                endif
            enddo

            do m=0,6
                if(norm2(Coords(i,1:2)+ tm(m,:) - r4).lt.a0*.1_dp) then
                    counter2 = counter2 + 1_dp
                    KekuleNeighbors(n,2,counter2) = i
                endif
            enddo

            do m=0,6
                if(norm2(Coords(i,1:2)+ tm(m,:) - r5).lt.a0*.1_dp) then
                    counter2 = counter2 + 1_dp
                    KekuleNeighbors(n,2,counter2) = i
                endif
            enddo

            do m=0,6
                if(norm2(Coords(i,1:2)+ tm(m,:) - r6).lt.a0*.1_dp) then
                    counter2 = counter2 + 1_dp
                    KekuleNeighbors(n,2,counter2) = i
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
                KekuleNeighbors(n-ndim/4,1,counter1) = i
            endif
            enddo
            
            do m=0,6
            if(norm2(Coords(i,1:2)+tm(m,:) - r2).lt.a0*.1_dp) then
                counter1 = counter1 + 1_dp
                KekuleNeighbors(n-ndim/4,1,counter1) = i
            endif
            enddo

            do m=0,6
            if(norm2(Coords(i,1:2)+tm(m,:) - r3).lt.a0*.1_dp) then
                counter1 = counter1 + 1_dp
                KekuleNeighbors(n-ndim/4,1,counter1) = i
            endif
            enddo

            do m=0,6
            if(norm2(Coords(i,1:2)+tm(m,:) - r4).lt.a0*.1_dp) then
                counter2 = counter2 + 1_dp
                KekuleNeighbors(n-ndim/4,2,counter2) = i
            endif
            enddo

            do m=0,6
            if(norm2(Coords(i,1:2)+tm(m,:) - r5).lt.a0*.1_dp) then
                counter2 = counter2 + 1_dp
                KekuleNeighbors(n-ndim/4,2,counter2) = i
            endif
            enddo

            do m=0,6
            if(norm2(Coords(i,1:2)+tm(m,:) - r6).lt.a0*.1_dp) then
                counter2 = counter2 + 1_dp
                KekuleNeighbors(n-ndim/4,2,counter2) = i
            endif
            enddo

        enddo
    enddo

end subroutine getKekuleNeighbors

subroutine getKekuleLattice(Coords,ndim,RotMatrix,a1,a2,tnn,KekuleLattice)

    integer(dp), intent(in) :: ndim
    real(dp), intent(in) :: Coords(ndim,3), RotMatrix(2,2), a1(2), a2(2), tnn(6,2)
    integer(dp), intent(inout) :: KekuleLattice(ndim/2,6,2)

    integer(dp) :: n,h, i, ndimq
    real(dp) :: a11(2), a21(2), a12(2), a22(2), nCoords(2), tn(0:6,2)

    ndimq = ndim
    tn(1:6,:) = tnn(:,:)
    tn(0,1) = 0.0_dp
    tn(0,2) = 0.0_dp

    KekuleLattice(:,:,:) = 0_dp

    a11 = matmul(RotMatrix,a1)
    a21 = matmul(RotMatrix,a2)
    a12 = matmul(transpose(RotMatrix),a1)
    a22 = matmul(transpose(RotMatrix),a2)

    do n=1,ndim/4

        KekuleLattice(n,1,1) = n
        KekuleLattice(n,1,2) = 0

        nCoords = Coords(n,1:2) + a11/3.0 - 2.0*a21/3.0
        do h = ndimq/4+1,ndimq/2
            do i=0,6
                if (norm2(Coords(h,1:2) + tn(i,:) - nCoords).lt.a0*.1) then
                    
                    KekuleLattice(n,2,1) = h
                    KekuleLattice(n,2,2) = i
                endif
            enddo
        enddo
        
        nCoords = Coords(n,1:2) + a11 - a21
        do h = 1,ndimq/4
            do i=0,6
                if (norm2(Coords(h,1:2) + tn(i,:) - nCoords).lt.a0*.1) then

                    KekuleLattice(n,3,1) = h
                    KekuleLattice(n,3,2) = i
                endif
            enddo
        enddo
        
        nCoords = Coords(n,1:2) + 4*a11/3.0 - 2*a21/3.0
        do h=ndimq/4+1,ndimq/2
            do i=0,6
                if (norm2(Coords(h,1:2) + tn(i,:) - nCoords).lt.a0*.1) then

                    KekuleLattice(n,4,1) = h
                    KekuleLattice(n,4,2) = i
                endif
            enddo
        enddo
        
        nCoords = Coords(n,1:2) + a11
        do h=1,ndimq/4
            do i=0,6
                if (norm2(Coords(h,1:2) + tn(i,:) - nCoords).lt.a0*.1) then

                    KekuleLattice(n,5,1) = h
                    KekuleLattice(n,5,2) = i
                endif
            enddo
        enddo
        
        nCoords = Coords(n,1:2) + a11/3.0 + a21/3.0
        do h=ndimq/4+1,ndimq/2
            do i=0,6
                if (norm2(Coords(h,1:2) + tn(i,:) - nCoords).lt.a0*.1) then

                    KekuleLattice(n,6,1) = h
                    KekuleLattice(n,6,2) = i
                endif
            enddo
        enddo

        if(KekuleLattice(n,2,1)*KekuleLattice(n,3,1)*KekuleLattice(n,4,1)*KekuleLattice(n,5,1)*KekuleLattice(n,6,1).eq.0)then
            KekuleLattice(n,2:6,:) = 0_dp
            write(*,*) 'ERROR getKekuleLattice'
        endif


    enddo

    do n=ndim/2+1,3*ndim/4

        KekuleLattice(n-ndim/4,1,1) = n
        KekuleLattice(n-ndim/4,1,2) = 0


        nCoords = Coords(n,1:2) + a12/3.0 - 2*a22/3.0
        do h = 3*ndimq/4+1,ndimq
            do i=0,6
                if (norm2(Coords(h,1:2) + tn(i,:) - nCoords).lt.a0*.1) then
                    KekuleLattice(n-ndim/4,2,1) = h
                    KekuleLattice(n-ndim/4,2,2) = i
                endif
            enddo
        enddo
        
        nCoords = Coords(n,1:2) + a12 - a22
        do h=ndimq/2+1,3*ndimq/4
            do i=0,6
                if (norm2(Coords(h,1:2) + tn(i,:) - nCoords).lt.a0*.1) then
                    KekuleLattice(n-ndim/4,3,1) = h
                    KekuleLattice(n-ndim/4,3,2) = i
                endif
            enddo
        enddo
        
        nCoords = Coords(n,1:2) + 4.0*a12/3.0 - 2.0*a22/3.0
        do h=3*ndimq/4+1,ndimq
            do i=0,6
                if (norm2(Coords(h,1:2) + tn(i,:) - nCoords).lt.a0*.1) then
                    KekuleLattice(n-ndim/4,4,1) = h
                    KekuleLattice(n-ndim/4,4,2) = i
                endif
            enddo
        enddo
        
        nCoords = Coords(n,1:2) + a12
        do h=ndimq/2+1,3*ndimq/4
            do i=0,6
                if (norm2(Coords(h,1:2) + tn(i,:) - nCoords).lt.a0*.1) then
                    KekuleLattice(n-ndim/4,5,1) = h
                    KekuleLattice(n-ndim/4,5,2) = i
                endif
            enddo
        enddo
        
        nCoords = Coords(n,1:2) + a12/3.0 + a22/3.0
        do h=3*ndimq/4+1,ndimq
            do i=0,6
                if (norm2(Coords(h,1:2) + tn(i,:) - nCoords).lt.a0*.1) then
                    KekuleLattice(n-ndim/4,6,1) = h
                    KekuleLattice(n-ndim/4,6,2) = i
                endif
            enddo
        enddo

        if(KekuleLattice(n-ndim/4,2,1)*KekuleLattice(n-ndim/4,3,1)* &
            KekuleLattice(n-ndim/4,4,1)*KekuleLattice(n-ndim/4,5,1)*KekuleLattice(n-ndim/4,6,1).eq.0)then
            KekuleLattice(n-ndim/4,2:6,:) = 0
            write(*,*) 'ERROR getKekuleLattice'
        endif

    enddo

end subroutine getKekuleLattice


subroutine IntraSubIntraValAlt(ValleyPol,NormSquared,ndim,NearestNeighborsUC,NearestNeighborsT,ind,nUnitCell_1,nUnitCell_2,TnTonUnitCell12,zFock)

    integer(dp), intent(in) :: ndim, ind
    integer(dp), intent(in) :: NearestNeighborsUC(ndim,3),NearestNeighborsT(ndim,3), nUnitCell_1(ind), nUnitCell_2(ind), TnTonUnitCell12(0:6,2)
    real(dp), intent(out) :: ValleyPol(ndim), NormSquared(ndim)
    complex(dp), intent(in) :: zFock(ndim,ndim,ind)
    
    integer(dp) :: n,m, index
    real(dp) :: phi0
    real(dp) :: hA1,hB1,hA2,hB2
    complex(dp) :: ztemp,z1,z2,z3, TempLoop(ndim)
    
    ValleyPol(:) = 0.0_dp
    NormSquared(:) = 0.0_dp
    TempLoop(:) = cmplx(0.0_dp,0.0_dp,dp)

    do n=1,ndim

        if(numI.eq.0)then
            z1 =  zFock(NearestNeighborsUC(n,1),NearestNeighborsUC(n,3),1)
            z2 = zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,2),1)
            z3 = zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,1),1)
        
        else

            call nUnitCell12ToInd(index, -TnTonUnitCell12(NearestNeighborsT(n,1),1) + TnTonUnitCell12(NearestNeighborsT(n,3),1),&
            -TnTonUnitCell12(NearestNeighborsT(n,1),2) + TnTonUnitCell12(NearestNeighborsT(n,3),2), nUnitCell_1, nUnitCell_2, ind)  
            if (index.eq.-1_dp)then
                call nUnitCell12ToInd(index, +TnTonUnitCell12(NearestNeighborsT(n,1),1) - TnTonUnitCell12(NearestNeighborsT(n,3),1),&
                +TnTonUnitCell12(NearestNeighborsT(n,1),2) - TnTonUnitCell12(NearestNeighborsT(n,3),2), nUnitCell_1, nUnitCell_2, ind)
                
                z1 = conjg(zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,1),index))
            else
                
                z1 = zFock(NearestNeighborsUC(n,1),NearestNeighborsUC(n,3),index)
            endif
            
            ! if((n.eq.3*ndim/4).or.(n.eq.ndim))then
            ! write(*,*) n,index,z1,zFock(NearestNeighborsUC(n,1),NearestNeighborsUC(n,3),1)
            ! endif

            call nUnitCell12ToInd(index, -TnTonUnitCell12(NearestNeighborsT(n,3),1) + TnTonUnitCell12(NearestNeighborsT(n,2),1),&
            -TnTonUnitCell12(NearestNeighborsT(n,3),2) + TnTonUnitCell12(NearestNeighborsT(n,2),2), nUnitCell_1, nUnitCell_2, ind)  
            if (index.eq.-1_dp)then
                call nUnitCell12ToInd(index, +TnTonUnitCell12(NearestNeighborsT(n,3),1) - TnTonUnitCell12(NearestNeighborsT(n,2),1),&
                +TnTonUnitCell12(NearestNeighborsT(n,3),2) - TnTonUnitCell12(NearestNeighborsT(n,2),2), nUnitCell_1, nUnitCell_2, ind)
                
                z2 = conjg(zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,3),index))
            else
                
                z2 = zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,2),index)
            endif
            !write(*,*) n,index,z2,zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,2),1)
            ! if((n.eq.3*ndim/4).or.(n.eq.ndim))then
            ! write(*,*) n,index,z2,zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,2),1)
            ! endif
            
            call nUnitCell12ToInd(index, -TnTonUnitCell12(NearestNeighborsT(n,2),1) + TnTonUnitCell12(NearestNeighborsT(n,1),1),&
            -TnTonUnitCell12(NearestNeighborsT(n,2),2) + TnTonUnitCell12(NearestNeighborsT(n,1),2), nUnitCell_1, nUnitCell_2, ind)  
            if (index.eq.-1_dp)then
                call nUnitCell12ToInd(index, +TnTonUnitCell12(NearestNeighborsT(n,2),1) - TnTonUnitCell12(NearestNeighborsT(n,1),1),&
                +TnTonUnitCell12(NearestNeighborsT(n,2),2) - TnTonUnitCell12(NearestNeighborsT(n,1),2), nUnitCell_1, nUnitCell_2, ind)
                
                z3 = conjg(zFock(NearestNeighborsUC(n,1),NearestNeighborsUC(n,2),index))
            else
                
                z3 = zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,1),index)
            endif
            !write(*,*) n,index,z3,zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,1),1)
            ! if((n.eq.3*ndim/4).or.(n.eq.ndim))then
            ! write(*,*) n,index,z3,zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,1),1),NearestNeighborsUC(n,2),NearestNeighborsUC(n,1),NearestNeighborsT(n,2),NearestNeighborsT(n,1)
            ! endif

        endif

        ztemp = z1 + z2 + z3
        TempLoop(n) = ztemp

    enddo
    
    do n=1,ndim/4
        m=n
        ValleyPol(m) = 2.0_dp/3.0_dp/sqrt(3.0_dp)*aimag(TempLoop(m))
        NormSquared(m) = -2.0_dp/3.0_dp*real(TempLoop(m))

        m=n+ndim/4
        ValleyPol(m) = -2.0_dp/3.0_dp/sqrt(3.0_dp)*aimag(TempLoop(m))
        NormSquared(m) = -2.0_dp/3.0_dp*real(TempLoop(m))

        m=n+ndim/2
        ValleyPol(m) = 2.0_dp/3.0_dp/sqrt(3.0_dp)*aimag(TempLoop(m))
        NormSquared(m) = -2.0_dp/3.0_dp*real(TempLoop(m))

        m=n+3*ndim/4
        ValleyPol(m) = -2.0_dp/3.0_dp/sqrt(3.0_dp)*aimag(TempLoop(m))
        NormSquared(m) = -2.0_dp/3.0_dp*real(TempLoop(m))

    enddo


    return
end subroutine IntraSubIntraValAlt

subroutine IntraSubInterValAlt(fKp_conjxfK,ndim,ind,nUnitCell_1,nUnitCell_2,NearestNeighborsUC,NearestNeighborsT,TnTonUnitCell12,zFock)

    integer(dp), intent(in) :: ndim, ind
    integer(dp), intent(in) :: NearestNeighborsUC(ndim,3),NearestNeighborsT(ndim,3), TnTonUnitCell12(0:6,2), nUnitCell_1(ind), nUnitCell_2(ind)
    complex(dp), intent(out) :: fKp_conjxfK(ndim)
    complex(dp), intent(in) :: zFock(ndim,ndim,ind)
    
    integer(dp) :: n,m, index
    complex(dp) :: ztemp,z1,z2,z3
    complex(dp) :: zi,zPi3,z2Pi3
    
    zi=cmplx(0._dp,1._dp,dp)
    zPi3=exp(cmplx(0._dp,1._dp,dp)*pi/3._dp)
    z2Pi3=exp(cmplx(0.0_dp,2._dp*pi/3._dp,dp))

    fKp_conjxfK(:) = cmplx(0.0_dp,0.0_dp,dp)

    do n=1,ndim/4

        if(numI.EQ.0)then
            z1 = zFock(NearestNeighborsUC(n,1),NearestNeighborsUC(n,3),1)*z2Pi3
            z2 = zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,2),1)/z2Pi3
            z3 = zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,1),1)
            ztemp= z1 + z2 + z3
        else

            call nUnitCell12ToInd(index, -TnTonUnitCell12(NearestNeighborsT(n,1),1) + TnTonUnitCell12(NearestNeighborsT(n,3),1),&
            -TnTonUnitCell12(NearestNeighborsT(n,1),2) + TnTonUnitCell12(NearestNeighborsT(n,3),2), nUnitCell_1, nUnitCell_2, ind)  
            if (index.eq.-1_dp)then
                call nUnitCell12ToInd(index, +TnTonUnitCell12(NearestNeighborsT(n,1),1) - TnTonUnitCell12(NearestNeighborsT(n,3),1),&
                +TnTonUnitCell12(NearestNeighborsT(n,1),2) - TnTonUnitCell12(NearestNeighborsT(n,3),2), nUnitCell_1, nUnitCell_2, ind)
                
                z1 = conjg(zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,1),index))
            else
                
                z1 = zFock(NearestNeighborsUC(n,1),NearestNeighborsUC(n,3),index)
            endif
            !write(*,*) n,index,z1,zFock(NearestNeighborsUC(n,1),NearestNeighborsUC(n,3),1)

            call nUnitCell12ToInd(index, -TnTonUnitCell12(NearestNeighborsT(n,3),1) + TnTonUnitCell12(NearestNeighborsT(n,2),1),&
            -TnTonUnitCell12(NearestNeighborsT(n,3),2) + TnTonUnitCell12(NearestNeighborsT(n,2),2), nUnitCell_1, nUnitCell_2, ind)  
            if (index.eq.-1_dp)then
                call nUnitCell12ToInd(index, +TnTonUnitCell12(NearestNeighborsT(n,3),1) - TnTonUnitCell12(NearestNeighborsT(n,2),1),&
                +TnTonUnitCell12(NearestNeighborsT(n,3),2) - TnTonUnitCell12(NearestNeighborsT(n,2),2), nUnitCell_1, nUnitCell_2, ind)
                
                z2 = conjg(zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,3),index))
            else
                
                z2 = zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,2),index)
            endif
            !write(*,*) n,index,z2,zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,2),1)

            call nUnitCell12ToInd(index, -TnTonUnitCell12(NearestNeighborsT(n,2),1) + TnTonUnitCell12(NearestNeighborsT(n,1),1),&
            -TnTonUnitCell12(NearestNeighborsT(n,2),2) + TnTonUnitCell12(NearestNeighborsT(n,1),2), nUnitCell_1, nUnitCell_2, ind)  
            if (index.eq.-1_dp)then
                call nUnitCell12ToInd(index, +TnTonUnitCell12(NearestNeighborsT(n,2),1) - TnTonUnitCell12(NearestNeighborsT(n,1),1),&
                +TnTonUnitCell12(NearestNeighborsT(n,2),2) - TnTonUnitCell12(NearestNeighborsT(n,1),2), nUnitCell_1, nUnitCell_2, ind)
                
                z3 = conjg(zFock(NearestNeighborsUC(n,1),NearestNeighborsUC(n,2),index))
            else
                
                z3 = zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,1),index)
            endif
            !write(*,*) n,index,z3,zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,1),1)

            ztemp = z1*z2Pi3 + z2/z2Pi3 + z3
        endif

        fKp_conjxfK(n) = ztemp/3.0_dp*1.0 

    enddo

    do n=ndim/4+1,ndim/2

        if(numI.eq.0)then
            z1 = zFock(NearestNeighborsUC(n,1),NearestNeighborsUC(n,3),1)*z2Pi3
            z2 = zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,2),1)
            z3 = zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,1),1)/z2Pi3
            ztemp= z1 + z2 + z3
        else

                call nUnitCell12ToInd(index, -TnTonUnitCell12(NearestNeighborsT(n,1),1) + TnTonUnitCell12(NearestNeighborsT(n,3),1),&
            -TnTonUnitCell12(NearestNeighborsT(n,1),2) + TnTonUnitCell12(NearestNeighborsT(n,3),2), nUnitCell_1, nUnitCell_2, ind)  
            if (index.eq.-1_dp)then
                call nUnitCell12ToInd(index, +TnTonUnitCell12(NearestNeighborsT(n,1),1) - TnTonUnitCell12(NearestNeighborsT(n,3),1),&
                +TnTonUnitCell12(NearestNeighborsT(n,1),2) - TnTonUnitCell12(NearestNeighborsT(n,3),2), nUnitCell_1, nUnitCell_2, ind)
                
                z1 = conjg(zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,1),index))
            else
                
                z1 = zFock(NearestNeighborsUC(n,1),NearestNeighborsUC(n,3),index)
            endif
            !write(*,*) n,index,z1,zFock(NearestNeighborsUC(n,1),NearestNeighborsUC(n,3),1)

            call nUnitCell12ToInd(index, -TnTonUnitCell12(NearestNeighborsT(n,3),1) + TnTonUnitCell12(NearestNeighborsT(n,2),1),&
            -TnTonUnitCell12(NearestNeighborsT(n,3),2) + TnTonUnitCell12(NearestNeighborsT(n,2),2), nUnitCell_1, nUnitCell_2, ind)  
            if (index.eq.-1_dp)then
                call nUnitCell12ToInd(index, +TnTonUnitCell12(NearestNeighborsT(n,3),1) - TnTonUnitCell12(NearestNeighborsT(n,2),1),&
                +TnTonUnitCell12(NearestNeighborsT(n,3),2) - TnTonUnitCell12(NearestNeighborsT(n,2),2), nUnitCell_1, nUnitCell_2, ind)
                
                z2 = conjg(zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,3),index))
            else
                
                z2 = zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,2),index)
            endif
            !write(*,*) n,index,z2,zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,2),1)

            call nUnitCell12ToInd(index, -TnTonUnitCell12(NearestNeighborsT(n,2),1) + TnTonUnitCell12(NearestNeighborsT(n,1),1),&
            -TnTonUnitCell12(NearestNeighborsT(n,2),2) + TnTonUnitCell12(NearestNeighborsT(n,1),2), nUnitCell_1, nUnitCell_2, ind)  
            if (index.eq.-1_dp)then
                call nUnitCell12ToInd(index, +TnTonUnitCell12(NearestNeighborsT(n,2),1) - TnTonUnitCell12(NearestNeighborsT(n,1),1),&
                +TnTonUnitCell12(NearestNeighborsT(n,2),2) - TnTonUnitCell12(NearestNeighborsT(n,1),2), nUnitCell_1, nUnitCell_2, ind)
                
                z3 = conjg(zFock(NearestNeighborsUC(n,1),NearestNeighborsUC(n,2),index))
            else
                
                z3 = zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,1),index)
            endif
            !write(*,*) n,index,z3,zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,1),1)

            ztemp = z1*z2Pi3 + z2 + z3/z2Pi3

        endif

        fKp_conjxfK(n) = ztemp/3.0_dp*1.0 

    enddo

    do n=ndim/2+1,3*ndim/4

        if(numI.eq.0)then
            z1 = zFock(NearestNeighborsUC(n,1),NearestNeighborsUC(n,3),1)*z2Pi3
            z2 = zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,2),1)/z2Pi3
            z3 = zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,1),1)
            ztemp = z1 + z2 + z3

        else
            call nUnitCell12ToInd(index, -TnTonUnitCell12(NearestNeighborsT(n,1),1) + TnTonUnitCell12(NearestNeighborsT(n,3),1),&
            -TnTonUnitCell12(NearestNeighborsT(n,1),2) + TnTonUnitCell12(NearestNeighborsT(n,3),2), nUnitCell_1, nUnitCell_2, ind)  
            if (index.eq.-1_dp)then
                call nUnitCell12ToInd(index, +TnTonUnitCell12(NearestNeighborsT(n,1),1) - TnTonUnitCell12(NearestNeighborsT(n,3),1),&
                +TnTonUnitCell12(NearestNeighborsT(n,1),2) - TnTonUnitCell12(NearestNeighborsT(n,3),2), nUnitCell_1, nUnitCell_2, ind)
                
                z1 = conjg(zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,1),index))
            else
                
                z1 = zFock(NearestNeighborsUC(n,1),NearestNeighborsUC(n,3),index)
            endif
            !write(*,*) n,index,z1,zFock(NearestNeighborsUC(n,1),NearestNeighborsUC(n,3),1)

            call nUnitCell12ToInd(index, -TnTonUnitCell12(NearestNeighborsT(n,3),1) + TnTonUnitCell12(NearestNeighborsT(n,2),1),&
            -TnTonUnitCell12(NearestNeighborsT(n,3),2) + TnTonUnitCell12(NearestNeighborsT(n,2),2), nUnitCell_1, nUnitCell_2, ind)  
            if (index.eq.-1_dp)then
                call nUnitCell12ToInd(index, +TnTonUnitCell12(NearestNeighborsT(n,3),1) - TnTonUnitCell12(NearestNeighborsT(n,2),1),&
                +TnTonUnitCell12(NearestNeighborsT(n,3),2) - TnTonUnitCell12(NearestNeighborsT(n,2),2), nUnitCell_1, nUnitCell_2, ind)
                
                z2 = conjg(zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,3),index))
            else
                
                z2 = zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,2),index)
            endif
            !write(*,*) n,index,z2,zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,2),1)

            call nUnitCell12ToInd(index, -TnTonUnitCell12(NearestNeighborsT(n,2),1) + TnTonUnitCell12(NearestNeighborsT(n,1),1),&
            -TnTonUnitCell12(NearestNeighborsT(n,2),2) + TnTonUnitCell12(NearestNeighborsT(n,1),2), nUnitCell_1, nUnitCell_2, ind)  
            if (index.eq.-1_dp)then
                call nUnitCell12ToInd(index, +TnTonUnitCell12(NearestNeighborsT(n,2),1) - TnTonUnitCell12(NearestNeighborsT(n,1),1),&
                +TnTonUnitCell12(NearestNeighborsT(n,2),2) - TnTonUnitCell12(NearestNeighborsT(n,1),2), nUnitCell_1, nUnitCell_2, ind)
                
                z3 = conjg(zFock(NearestNeighborsUC(n,1),NearestNeighborsUC(n,2),index))
            else
                
                z3 = zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,1),index)
            endif
            !write(*,*) n,index,z3,zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,1),1)

            ztemp = z1*z2Pi3 + z2/z2Pi3 + z3
        endif
        
        fKp_conjxfK(n) = ztemp/3.0_dp*1.0 
    enddo

    do n=3*ndim/4+1,ndim

        if (numI.eq.0) then
            z1 = zFock(NearestNeighborsUC(n,1),NearestNeighborsUC(n,3),1)*z2Pi3
            z2 = zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,2),1)
            z3 = zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,1),1)/z2Pi3
            ztemp = z1 + z2 + z3
        else

            call nUnitCell12ToInd(index, -TnTonUnitCell12(NearestNeighborsT(n,1),1) + TnTonUnitCell12(NearestNeighborsT(n,3),1),&
            -TnTonUnitCell12(NearestNeighborsT(n,1),2) + TnTonUnitCell12(NearestNeighborsT(n,3),2), nUnitCell_1, nUnitCell_2, ind)  
            if (index.eq.-1_dp)then
                call nUnitCell12ToInd(index, +TnTonUnitCell12(NearestNeighborsT(n,1),1) - TnTonUnitCell12(NearestNeighborsT(n,3),1),&
                +TnTonUnitCell12(NearestNeighborsT(n,1),2) - TnTonUnitCell12(NearestNeighborsT(n,3),2), nUnitCell_1, nUnitCell_2, ind)
                
                z1 = conjg(zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,1),index))
            else
                
                z1 = zFock(NearestNeighborsUC(n,1),NearestNeighborsUC(n,3),index)
            endif
            !write(*,*) n,index,z1,zFock(NearestNeighborsUC(n,1),NearestNeighborsUC(n,3),1)

            call nUnitCell12ToInd(index, -TnTonUnitCell12(NearestNeighborsT(n,3),1) + TnTonUnitCell12(NearestNeighborsT(n,2),1),&
            -TnTonUnitCell12(NearestNeighborsT(n,3),2) + TnTonUnitCell12(NearestNeighborsT(n,2),2), nUnitCell_1, nUnitCell_2, ind)  
            if (index.eq.-1_dp)then
                call nUnitCell12ToInd(index, +TnTonUnitCell12(NearestNeighborsT(n,3),1) - TnTonUnitCell12(NearestNeighborsT(n,2),1),&
                +TnTonUnitCell12(NearestNeighborsT(n,3),2) - TnTonUnitCell12(NearestNeighborsT(n,2),2), nUnitCell_1, nUnitCell_2, ind)
                
                z2 = conjg(zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,3),index))
            else
                
                z2 = zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,2),index)
            endif
            !write(*,*) n,index,z2,zFock(NearestNeighborsUC(n,3),NearestNeighborsUC(n,2),1)

            call nUnitCell12ToInd(index, -TnTonUnitCell12(NearestNeighborsT(n,2),1) + TnTonUnitCell12(NearestNeighborsT(n,1),1),&
            -TnTonUnitCell12(NearestNeighborsT(n,2),2) + TnTonUnitCell12(NearestNeighborsT(n,1),2), nUnitCell_1, nUnitCell_2, ind)  
            if (index.eq.-1_dp)then
                call nUnitCell12ToInd(index, +TnTonUnitCell12(NearestNeighborsT(n,2),1) - TnTonUnitCell12(NearestNeighborsT(n,1),1),&
                +TnTonUnitCell12(NearestNeighborsT(n,2),2) - TnTonUnitCell12(NearestNeighborsT(n,1),2), nUnitCell_1, nUnitCell_2, ind)
                
                z3 = conjg(zFock(NearestNeighborsUC(n,1),NearestNeighborsUC(n,2),index))
            else
                
                z3 = zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,1),index)
            endif
            !write(*,*) n,index,z3,zFock(NearestNeighborsUC(n,2),NearestNeighborsUC(n,1),1)
            
            ztemp = z1*z2Pi3 + z2 + z3/z2Pi3
        endif

        fKp_conjxfK(n) = ztemp/3.0_dp*1.0 
    enddo
        
    return
end subroutine IntraSubInterValAlt

subroutine InterSubInterValAlt(fKA_conjxfKpB,fKB_conjxfKpA,KekuleLattice,KekuleNeighbors,ndim,ind,nUnitCell_1,nUnitCell_2,TnTonUnitCell12,zFock)

    integer(dp), intent(in) :: ndim, ind
    integer(dp), intent(in) :: nUnitCell_1(ind), nUnitCell_2(ind), TnTonUnitCell12(0:6,2)
    integer(dp), intent(in) :: KekuleLattice(ndim/2,6,2), KekuleNeighbors(ndim/2,2,3)
    complex(dp), intent(out) :: fKA_conjxfKpB(ndim/2), fKB_conjxfKpA(ndim/2)
    complex(dp), intent(in) :: zFock(ndim,ndim,ind)
    
    integer(dp) :: n, m, index
    complex(dp) :: zi,zPi3,z2Pi3,z1,z2,z3,z4,z5,z6,delta0,delta1,delta2
    complex(dp) :: CurrentTemp(ndim/2)

    zi=cmplx(0.0_dp,1._dp,dp)
    zPi3=exp(cmplx(0.0_dp,1._dp,dp)*pi/3._dp)
    z2Pi3=exp(cmplx(0.0_dp,2._dp*pi/3._dp,dp))

    CurrentTemp(:) = cmplx(0.0_dp,0.0_dp,dp)
    fKA_conjxfKpB(:) = cmplx(0.0_dp,0.0_dp,dp)
    fKB_conjxfKpA(:) = cmplx(0.0_dp,0.0_dp,dp)

    do n=1,ndim/2

        if(numI.eq.0)then

            z1 = zFock(KekuleLattice(n,1,1),KekuleLattice(n,2,1),1)
            z2 = zFock(KekuleLattice(n,2,1),KekuleLattice(n,3,1),1)
            z3 = zFock(KekuleLattice(n,3,1),KekuleLattice(n,4,1),1)
            z4 = zFock(KekuleLattice(n,4,1),KekuleLattice(n,5,1),1)
            z5 = zFock(KekuleLattice(n,5,1),KekuleLattice(n,6,1),1)
            z6 = zFock(KekuleLattice(n,6,1),KekuleLattice(n,1,1),1)
           
        else

            call nUnitCell12ToInd(index, TnTonUnitCell12(KekuleLattice(n,1,2),1) - TnTonUnitCell12(KekuleLattice(n,2,2),1),&
                TnTonUnitCell12(KekuleLattice(n,1,2),2) - TnTonUnitCell12(KekuleLattice(n,2,2),2), nUnitCell_1, nUnitCell_2, ind)  
            if(index.eq.-1_dp)then
                call nUnitCell12ToInd(index, -TnTonUnitCell12(KekuleLattice(n,1,2),1) + TnTonUnitCell12(KekuleLattice(n,2,2),1),&
                -TnTonUnitCell12(KekuleLattice(n,1,2),2) + TnTonUnitCell12(KekuleLattice(n,2,2),2), nUnitCell_1, nUnitCell_2, ind) 

                z1 = conjg(zFock(KekuleLattice(n,2,1),KekuleLattice(n,1,1),index))
            else

                z1 = zFock(KekuleLattice(n,1,1),KekuleLattice(n,2,1),index)
            endif

            call nUnitCell12ToInd(index, TnTonUnitCell12(KekuleLattice(n,2,2),1) - TnTonUnitCell12(KekuleLattice(n,3,2),1),&
                TnTonUnitCell12(KekuleLattice(n,2,2),2) - TnTonUnitCell12(KekuleLattice(n,3,2),2), nUnitCell_1, nUnitCell_2, ind)  
            if(index.eq.-1_dp)then
                call nUnitCell12ToInd(index, -TnTonUnitCell12(KekuleLattice(n,2,2),1) + TnTonUnitCell12(KekuleLattice(n,3,2),1),&
                -TnTonUnitCell12(KekuleLattice(n,2,2),2) + TnTonUnitCell12(KekuleLattice(n,3,2),2), nUnitCell_1, nUnitCell_2, ind) 

                z2 = conjg(zFock(KekuleLattice(n,3,1),KekuleLattice(n,2,1),index))
            else

                z2 = zFock(KekuleLattice(n,2,1),KekuleLattice(n,3,1),index)
            endif
            
            call nUnitCell12ToInd(index, TnTonUnitCell12(KekuleLattice(n,3,2),1) - TnTonUnitCell12(KekuleLattice(n,4,2),1),&
                TnTonUnitCell12(KekuleLattice(n,3,2),2) - TnTonUnitCell12(KekuleLattice(n,4,2),2), nUnitCell_1, nUnitCell_2, ind)  
            if(index.eq.-1_dp)then
                call nUnitCell12ToInd(index, -TnTonUnitCell12(KekuleLattice(n,3,2),1) + TnTonUnitCell12(KekuleLattice(n,4,2),1),&
                -TnTonUnitCell12(KekuleLattice(n,3,2),2) + TnTonUnitCell12(KekuleLattice(n,4,2),2), nUnitCell_1, nUnitCell_2, ind) 

                z3 = conjg(zFock(KekuleLattice(n,4,1),KekuleLattice(n,3,1),index))
            else

                z3 = zFock(KekuleLattice(n,3,1),KekuleLattice(n,4,1),index)
            endif
            
            call nUnitCell12ToInd(index, TnTonUnitCell12(KekuleLattice(n,4,2),1) - TnTonUnitCell12(KekuleLattice(n,5,2),1),&
            TnTonUnitCell12(KekuleLattice(n,4,2),2) - TnTonUnitCell12(KekuleLattice(n,5,2),2), nUnitCell_1, nUnitCell_2, ind)  
            if(index.eq.-1_dp)then
                call nUnitCell12ToInd(index, -TnTonUnitCell12(KekuleLattice(n,4,2),1) + TnTonUnitCell12(KekuleLattice(n,5,2),1),&
                -TnTonUnitCell12(KekuleLattice(n,4,2),2) + TnTonUnitCell12(KekuleLattice(n,5,2),2), nUnitCell_1, nUnitCell_2, ind) 

                z4 = conjg(zFock(KekuleLattice(n,5,1),KekuleLattice(n,4,1),index))
            else

                z4 = zFock(KekuleLattice(n,4,1),KekuleLattice(n,5,1),index)
            endif
            
            call nUnitCell12ToInd(index, TnTonUnitCell12(KekuleLattice(n,5,2),1) - TnTonUnitCell12(KekuleLattice(n,6,2),1),&
            TnTonUnitCell12(KekuleLattice(n,5,2),2) - TnTonUnitCell12(KekuleLattice(n,6,2),2), nUnitCell_1, nUnitCell_2, ind)  
            if(index.eq.-1_dp)then
                call nUnitCell12ToInd(index, -TnTonUnitCell12(KekuleLattice(n,5,2),1) + TnTonUnitCell12(KekuleLattice(n,6,2),1),&
                -TnTonUnitCell12(KekuleLattice(n,5,2),2) + TnTonUnitCell12(KekuleLattice(n,6,2),2), nUnitCell_1, nUnitCell_2, ind) 

                z5 = conjg(zFock(KekuleLattice(n,6,1),KekuleLattice(n,5,1),index))
            else

                z5 = zFock(KekuleLattice(n,5,1),KekuleLattice(n,6,1),index)
            endif
            
            call nUnitCell12ToInd(index, TnTonUnitCell12(KekuleLattice(n,6,2),1) - TnTonUnitCell12(KekuleLattice(n,1,2),1),&
            TnTonUnitCell12(KekuleLattice(n,6,2),2) - TnTonUnitCell12(KekuleLattice(n,1,2),2), nUnitCell_1, nUnitCell_2, ind)  
            if(index.eq.-1_dp)then
                call nUnitCell12ToInd(index, -TnTonUnitCell12(KekuleLattice(n,6,2),1) + TnTonUnitCell12(KekuleLattice(n,1,2),1),&
                -TnTonUnitCell12(KekuleLattice(n,6,2),2) + TnTonUnitCell12(KekuleLattice(n,1,2),2), nUnitCell_1, nUnitCell_2, ind) 

                z6 = conjg(zFock(KekuleLattice(n,1,1),KekuleLattice(n,6,1),index))
            else

                z6 = zFock(KekuleLattice(n,6,1),KekuleLattice(n,1,1),index)
            endif

        endif

        CurrentTemp(n) = z1+z2+z3+z4+z5+z6


    enddo

    do n=1,ndim/4

        !if (product(KekuleNeighbors(n,:,:)).ne.0_dp)then
        
        delta0 = CurrentTemp(n)
        delta1 = (CurrentTemp(KekuleNeighbors(n,1,1)) + CurrentTemp(KekuleNeighbors(n,1,2)) + CurrentTemp(KekuleNeighbors(n,1,3)))/3.0_dp 
        delta2 = (CurrentTemp(KekuleNeighbors(n,2,1)) + CurrentTemp(KekuleNeighbors(n,2,2)) + CurrentTemp(KekuleNeighbors(n,2,3)))/3.0_dp 

        fKA_conjxfKpB(n) = (zPi3*(conjg(delta0) + conjg(delta1)/z2Pi3 + conjg(delta2)*z2Pi3) +&
                    (delta0 + delta1/z2Pi3 + delta2*z2Pi3))*cmplx(0.0_dp,-1.0_dp,dp)/9.0_dp/sqrt(3.0_dp)

        fKB_conjxfKpA(n) = ((conjg(delta0) + conjg(delta1)/z2Pi3 + conjg(delta2)*z2Pi3) +&
                    zPi3*(delta0 + delta1/z2Pi3 + delta2*z2Pi3))*cmplx(0.0_dp,-1.0_dp,dp)/9.0_dp/sqrt(3.0_dp)
            
        !endif
                
    enddo

    do n=ndim/4+1,ndim/2

        !if (sum(KekuleNeighbors(n,:,:)).ne.0_dp)then
        
        delta0 = CurrentTemp(n)
        delta1 = (CurrentTemp(KekuleNeighbors(n,1,1)-ndim/4) + CurrentTemp(KekuleNeighbors(n,1,2)-ndim/4)&
                    + CurrentTemp(KekuleNeighbors(n,1,3)-ndim/4))/3.0_dp 
        delta2 = (CurrentTemp(KekuleNeighbors(n,2,1)-ndim/4) + CurrentTemp(KekuleNeighbors(n,2,2)-ndim/4)&
                    + CurrentTemp(KekuleNeighbors(n,2,3)-ndim/4))/3.0_dp 

        fKA_conjxfKpB(n) = (zPi3*(conjg(delta0) + conjg(delta1)/z2Pi3 + conjg(delta2)*z2Pi3) +&
                    (delta0 + delta1/z2Pi3 + delta2*z2Pi3))*cmplx(0.0_dp,-1.0_dp,dp)/9.0_dp/sqrt(3.0_dp)

        fKB_conjxfKpA(n) = ((conjg(delta0) + conjg(delta1)/z2Pi3 + conjg(delta2)*z2Pi3) +&
                    zPi3*(delta0 + delta1/z2Pi3 + delta2*z2Pi3))*cmplx(0.0_dp,-1.0_dp,dp)/9.0_dp/sqrt(3.0_dp)

            
        !endif
                
    enddo


    return
end subroutine InterSubInterValAlt


subroutine InterSubIntraValAlt(fKA_conjxfKB,fKpA_conjxfKpB,KekuleLattice,KekuleNeighbors,ndim,ind,nUnitCell_1,nUnitCell_2,TnTonUnitCell12,zFock)

    integer(dp), intent(in) :: ndim,ind
    integer(dp), intent(in) :: nUnitCell_1(ind), nUnitCell_2(ind)
    integer(dp), intent(in) :: KekuleLattice(ndim/2,6,2), KekuleNeighbors(ndim/2,2,3), TnTonUnitCell12(0:6,2)
    complex(dp), intent(out) :: fKA_conjxfKB(ndim/2), fKpA_conjxfKpB(ndim/2)
    complex(dp), intent(in) :: zFock(ndim,ndim,ind)
    
    integer(dp) :: n,m,index
    complex(dp) :: zi,zPi3,z2Pi3,z1,z2,z3,z4,z5,z6
    complex(dp) :: CurrentTemp(ndim/2,4)
    
    zi=cmplx(0.0_dp,1._dp,dp)
    zPi3=exp(cmplx(0.0_dp,1._dp,dp)*pi/3._dp)
    z2Pi3=exp(cmplx(0.0_dp,2._dp*pi/3._dp,dp))

    CurrentTemp(:,:) = cmplx(0.0_dp,0.0_dp,dp)
    fKA_conjxfKB(:) = cmplx(0.0_dp,0.0_dp,dp)
    fKpA_conjxfKpB(:) = cmplx(0.0_dp,0.0_dp,dp)

    do n=1,ndim/2

        if(numI.eq.0)then

            z1 = zFock(KekuleLattice(n,1,1),KekuleLattice(n,2,1),1)
            z2 = zFock(KekuleLattice(n,2,1),KekuleLattice(n,3,1),1)
            z3 = zFock(KekuleLattice(n,3,1),KekuleLattice(n,4,1),1)
            z4 = zFock(KekuleLattice(n,4,1),KekuleLattice(n,5,1),1)
            z5 = zFock(KekuleLattice(n,5,1),KekuleLattice(n,6,1),1)
            z6 = zFock(KekuleLattice(n,6,1),KekuleLattice(n,1,1),1)
        else

            call nUnitCell12ToInd(index, TnTonUnitCell12(KekuleLattice(n,1,2),1) - TnTonUnitCell12(KekuleLattice(n,2,2),1),&
                TnTonUnitCell12(KekuleLattice(n,1,2),2) - TnTonUnitCell12(KekuleLattice(n,2,2),2), nUnitCell_1, nUnitCell_2, ind)  
            if(index.eq.-1_dp)then
                call nUnitCell12ToInd(index, -TnTonUnitCell12(KekuleLattice(n,1,2),1) + TnTonUnitCell12(KekuleLattice(n,2,2),1),&
                -TnTonUnitCell12(KekuleLattice(n,1,2),2) + TnTonUnitCell12(KekuleLattice(n,2,2),2), nUnitCell_1, nUnitCell_2, ind) 

                z1 = conjg(zFock(KekuleLattice(n,2,1),KekuleLattice(n,1,1),index))
            else

                z1 = zFock(KekuleLattice(n,1,1),KekuleLattice(n,2,1),index)
            endif
            ! write(*,*) index

            call nUnitCell12ToInd(index, TnTonUnitCell12(KekuleLattice(n,2,2),1) - TnTonUnitCell12(KekuleLattice(n,3,2),1),&
                TnTonUnitCell12(KekuleLattice(n,2,2),2) - TnTonUnitCell12(KekuleLattice(n,3,2),2), nUnitCell_1, nUnitCell_2, ind)  
            if(index.eq.-1_dp)then
                call nUnitCell12ToInd(index, -TnTonUnitCell12(KekuleLattice(n,2,2),1) + TnTonUnitCell12(KekuleLattice(n,3,2),1),&
                -TnTonUnitCell12(KekuleLattice(n,2,2),2) + TnTonUnitCell12(KekuleLattice(n,3,2),2), nUnitCell_1, nUnitCell_2, ind) 

                z2 = conjg(zFock(KekuleLattice(n,3,1),KekuleLattice(n,2,1),index))
            else

                z2 = zFock(KekuleLattice(n,2,1),KekuleLattice(n,3,1),index)
            endif
            ! write(*,*) index

            
            call nUnitCell12ToInd(index, TnTonUnitCell12(KekuleLattice(n,3,2),1) - TnTonUnitCell12(KekuleLattice(n,4,2),1),&
                TnTonUnitCell12(KekuleLattice(n,3,2),2) - TnTonUnitCell12(KekuleLattice(n,4,2),2), nUnitCell_1, nUnitCell_2, ind)  
            if(index.eq.-1_dp)then
                call nUnitCell12ToInd(index, -TnTonUnitCell12(KekuleLattice(n,3,2),1) + TnTonUnitCell12(KekuleLattice(n,4,2),1),&
                -TnTonUnitCell12(KekuleLattice(n,3,2),2) + TnTonUnitCell12(KekuleLattice(n,4,2),2), nUnitCell_1, nUnitCell_2, ind) 

                z3 = conjg(zFock(KekuleLattice(n,4,1),KekuleLattice(n,3,1),index))
            else

                z3 = zFock(KekuleLattice(n,3,1),KekuleLattice(n,4,1),index)
            endif
            ! write(*,*) index

            
            call nUnitCell12ToInd(index, TnTonUnitCell12(KekuleLattice(n,4,2),1) - TnTonUnitCell12(KekuleLattice(n,5,2),1),&
            TnTonUnitCell12(KekuleLattice(n,4,2),2) - TnTonUnitCell12(KekuleLattice(n,5,2),2), nUnitCell_1, nUnitCell_2, ind)  
            if(index.eq.-1_dp)then
                call nUnitCell12ToInd(index, -TnTonUnitCell12(KekuleLattice(n,4,2),1) + TnTonUnitCell12(KekuleLattice(n,5,2),1),&
                -TnTonUnitCell12(KekuleLattice(n,4,2),2) + TnTonUnitCell12(KekuleLattice(n,5,2),2), nUnitCell_1, nUnitCell_2, ind) 

                z4 = conjg(zFock(KekuleLattice(n,5,1),KekuleLattice(n,4,1),index))
            else

                z4 = zFock(KekuleLattice(n,4,1),KekuleLattice(n,5,1),index)
            endif
            ! write(*,*) index

            
            call nUnitCell12ToInd(index, TnTonUnitCell12(KekuleLattice(n,5,2),1) - TnTonUnitCell12(KekuleLattice(n,6,2),1),&
            TnTonUnitCell12(KekuleLattice(n,5,2),2) - TnTonUnitCell12(KekuleLattice(n,6,2),2), nUnitCell_1, nUnitCell_2, ind)  
            if(index.eq.-1_dp)then
                call nUnitCell12ToInd(index, -TnTonUnitCell12(KekuleLattice(n,5,2),1) + TnTonUnitCell12(KekuleLattice(n,6,2),1),&
                -TnTonUnitCell12(KekuleLattice(n,5,2),2) + TnTonUnitCell12(KekuleLattice(n,6,2),2), nUnitCell_1, nUnitCell_2, ind) 

                z5 = conjg(zFock(KekuleLattice(n,6,1),KekuleLattice(n,5,1),index))
            else

                z5 = zFock(KekuleLattice(n,5,1),KekuleLattice(n,6,1),index)
            endif
            ! write(*,*) index

            
            call nUnitCell12ToInd(index, TnTonUnitCell12(KekuleLattice(n,6,2),1) - TnTonUnitCell12(KekuleLattice(n,1,2),1),&
            TnTonUnitCell12(KekuleLattice(n,6,2),2) - TnTonUnitCell12(KekuleLattice(n,1,2),2), nUnitCell_1, nUnitCell_2, ind)  
            if(index.eq.-1_dp)then
                call nUnitCell12ToInd(index, -TnTonUnitCell12(KekuleLattice(n,6,2),1) + TnTonUnitCell12(KekuleLattice(n,1,2),1),&
                -TnTonUnitCell12(KekuleLattice(n,6,2),2) + TnTonUnitCell12(KekuleLattice(n,1,2),2), nUnitCell_1, nUnitCell_2, ind) 

                z6 = conjg(zFock(KekuleLattice(n,1,1),KekuleLattice(n,6,1),index))
            else

                z6 = zFock(KekuleLattice(n,6,1),KekuleLattice(n,1,1),index)
            endif
            ! write(*,*) index

        endif

        CurrentTemp(n,1) = z1*z2Pi3 + z2/z2Pi3 + z3 + z4*z2Pi3 + z5/z2Pi3 + z6
        CurrentTemp(n,2) = z1/z2Pi3 + z2*z2Pi3 + z3*z2Pi3 + z4 + z5 + z6/z2Pi3

    enddo

    do n=1,ndim/2

        ! if (product(KekuleNeighbors(n,:,:)).ne.0_dp)then

            fKA_conjxfKB(n) = (CurrentTemp(n,2) - CurrentTemp(n,1)/z2Pi3)/3.0_dp*cmplx(0.0_dp,-1.0_dp,dp)/sqrt(3.0_dp)
            ! fKA_conjxfKB(n) = CurrentTemp(n,1)
            fKpA_conjxfKpB(n) = conjg((CurrentTemp(n,2) - CurrentTemp(n,1)*z2Pi3)/3.0_dp*cmplx(0.0_dp,1.0_dp,dp)/sqrt(3.0_dp))
            ! fKpA_conjxfKpB(n) = CurrentTemp(n,2)


        ! endif
                
    enddo

    return
end subroutine InterSubIntraValAlt

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

end module OrderPlot
