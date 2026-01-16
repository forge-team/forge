module OrderParameter

use omp_lib
use Setup
use Geometry
use TightBinding
use lapack_routines
implicit none

contains

!!!!!!!!!!!!!!!

subroutine NearestNeighbours(NearestNeighborInd,NearestNeighborCell,Coords,ndim,LatticeVectors)

    integer(dp), intent(in)    :: ndim
    integer(dp), intent(inout) :: NearestNeighborInd(ndim,3),NearestNeighborCell(ndim,3)
    real(dp), intent(in)    :: LatticeVectors(0:6,2), Coords(ndim,3)
    
    
    integer(dp) :: nlayer,i,j,n,ncount,ntemp,ntemp1,ntemp2,ntempt,ntempt1,ntempt2
    real(dp) :: rij,x,y,phi,temp,temp1,temp2,tempt,tempt1,tempt2,angleNN(ndim,3) 
    
    
    do nlayer=1,nlayers

        do i=ndim/nlayers*(nlayer-1)+1, ndim/nlayers*nlayer
            ncount=0
            do j=ndim/nlayers*(nlayer-1)+1, ndim/nlayers*nlayer
                            
                rij=(Coords(i,1)-Coords(j,1))**2+(Coords(i,2)-Coords(j,2))**2
                rij=sqrt(rij)
    
                if(rij.GT.(a0*1.2).AND.rij.LT.(a0*1.2))then
                    ncount=ncount+1
                    NearestNeighborInd(i,ncount)=j
                    NearestNeighborCell(i,ncount)=0
    
                    x=Coords(i,1)-Coords(j,1)
                    y=Coords(i,2)-Coords(j,2)
    
                    angleNN(i,ncount)=Angle(x,y)

                endif
    
                !!!!!!! Coupling to adjacent Wigner-Seitz cell
                do n=1,6
    
                    rij=(Coords(i,1)-Coords(j,1)+LatticeVectors(n,1))**2+(Coords(i,2)-Coords(j,2)+LatticeVectors(n,2))**2
                    rij=sqrt(rij)
    
                    if(rij.LT.(a0*1.2))then
                        ncount=ncount+1
                        NearestNeighborInd(i,ncount)=j
                        NearestNeighborCell(i,ncount)=n
    
                        x=Coords(i,1)-Coords(j,1)+LatticeVectors(n,1)
                        y=Coords(i,2)-Coords(j,2)+LatticeVectors(n,2)

                        angleNN(i,ncount)=Angle(x,y)

                    endif
    
                enddo
                
            enddo
        
            if(ncount.NE.3)then
                write(*,*) 'ERROR NEARESTNEIGHBOUR', i,ncount
            endif
        
            !!!!!!!!!!!!!!!!!!!
            !!!!!! Order cycle
            !!!!!!!!!!!!!!!!!!!
            if(angleNN(i,1).LT.angleNN(i,2))then
                temp=angleNN(i,1)
                ntemp=NearestNeighborInd(i,1)
                ntempt=NearestNeighborCell(i,1)
                angleNN(i,1)=angleNN(i,2)
                NearestNeighborInd(i,1)=NearestNeighborInd(i,2)
                NearestNeighborCell(i,1)=NearestNeighborCell(i,2)
                angleNN(i,2)=temp
                NearestNeighborInd(i,2)=ntemp
                NearestNeighborCell(i,2)=ntempt
            endif
        
            if(angleNN(i,1).LT.angleNN(i,3))then
                temp1=angleNN(i,1)
                ntemp1=NearestNeighborInd(i,1)
                ntempt1=NearestNeighborCell(i,1)
                temp2=angleNN(i,2)
                ntemp2=NearestNeighborInd(i,2)
                ntempt2=NearestNeighborCell(i,2)
        
                angleNN(i,1)=angleNN(i,3)
                NearestNeighborInd(i,1)=NearestNeighborInd(i,3)
                NearestNeighborCell(i,1)=NearestNeighborCell(i,3)
                angleNN(i,2)=temp1
                NearestNeighborInd(i,2)=ntemp1
                NearestNeighborCell(i,2)=ntempt1
                angleNN(i,3)=temp2
                NearestNeighborInd(i,3)=ntemp2
                NearestNeighborCell(i,3)=ntempt2
                else
                if(angleNN(i,2).LT.angleNN(i,3))then
                    temp=angleNN(i,2)
                    ntemp=NearestNeighborInd(i,2)
                    ntempt=NearestNeighborCell(i,2)
                    angleNN(i,2)=angleNN(i,3)
                    NearestNeighborInd(i,2)=NearestNeighborInd(i,3)
                    NearestNeighborCell(i,2)=NearestNeighborCell(i,3)
                    angleNN(i,3)=temp
                    NearestNeighborInd(i,3)=ntemp
                    NearestNeighborCell(i,3)=ntempt
                endif
            endif
        
        
        enddo
    enddo
     
    
end subroutine NearestNeighbours

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine HexagonalLoops(KekuleLoops,Coords,ndim,a1,a2,RotMatrix,LatticeVectors)

    integer(dp),  intent(in) :: ndim
    integer(dp), intent(out) :: KekuleLoops(ndim/2,2,3)
    real(dp), intent(in) :: Coords(ndim,3), RotMatrix(2,2), a1(2), a2(2), LatticeVectors(0:6,2)

    integer(dp) :: i,n,m,counter1,counter2,nlayer
    real(dp) :: a1r(2),a2r(2),r0(2),r1(2),r2(2),r3(2),r4(2),r5(2),r6(2)
    
    KekuleLoops(:,:,:) = 0_dp

    do nlayer=1,nlayers

        do n = ndim/nlayers*(nlayer-1)+1, ndim/nlayers*(nlayer-1) + ndim/nlayers*(nlayer-1) +  ndim/nlayers/2
            
            if (RotateLayers(nlayer).eq.-1)then
                a1 = matmul(RotMatrix,a1)
                a2 = matmul(RotMatrix,a2)
            if (RotateLayers(nlayer).eq.1)then
                a1 = matmul(transpose(RotMatrix),a1)
                a2 = matmul(transpose(RotMatrix),a2)
            endif

            r0 = Coords(n,1:2)

            r1 = r0 + a1r - a2r
            r2 = r0 + a2r
            r3 = r0 - a1r
            r4 = r0 + a1r
            r5 = r0 + a2r - a1r
            r6 = r0 - a2r

            counter1 = 0_dp
            counter2 = 0_dp

            do i = ndim/nlayers*(nlayer-1)+1, ndim/nlayers*(nlayer-1) + ndim/nlayers*(nlayer-1) +  ndim/nlayers/2

                do m=0,6
                if(norm2(Coords(i,1:2)+LatticeVectors(m,:) - r1).lt.(a0*.1_dp)) then
                    counter1 = counter1 + 1_dp
                    KekuleLoops(n-ndim/nlayers/2*(nlayer-1),1,counter1) = i
                endif
                enddo
                
                do m=0,6
                if(norm2(Coords(i,1:2)+LatticeVectors(m,:) - r2).lt.(a0*.1_dp)) then
                    counter1 = counter1 + 1_dp
                    KekuleLoops(n-ndim/nlayers/2*(nlayer-1),1,counter1) = i
                endif
                enddo

                do m=0,6
                if(norm2(Coords(i,1:2)+LatticeVectors(m,:) - r3).lt.(a0*.1_dp)) then
                    counter1 = counter1 + 1_dp
                    KekuleLoops(n-ndim/nlayers/2*(nlayer-1),1,counter1) = i
                endif
                enddo

                do m=0,6
                if(norm2(Coords(i,1:2)+LatticeVectors(m,:) - r4).lt.(a0*.1_dp)) then
                    counter2 = counter2 + 1_dp
                    KekuleLoops(n-ndim/nlayers/2*(nlayer-1),2,counter2) = i
                endif
                enddo

                do m=0,6
                if(norm2(Coords(i,1:2)+LatticeVectors(m,:) - r5).lt.(a0*.1_dp)) then
                    counter2 = counter2 + 1_dp
                    KekuleLoops(n-ndim/nlayers/2*(nlayer-1),2,counter2) = i
                endif
                enddo

                do m=0,6
                if(norm2(Coords(i,1:2)+LatticeVectors(m,:) - r6).lt.(a0*.1_dp)) then
                    counter2 = counter2 + 1_dp
                    KekuleLoops(n-ndim/nlayers/2*(nlayer-1),2,counter2) = i
                endif
                enddo

            enddo
        enddo
    enddo

end subroutine HexagonalLoops

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine KekuleNeighbors(Coords,ndim,RotMatrix,a1,a2,LatticeVectors,NeighborsKekule)

    integer(dp), intent(in) :: ndim
    real(dp), intent(in) :: Coords(ndim,3), RotMatrix(2,2), a1(2), a2(2), LatticeVectors(0:6,2)
    integer(dp), intent(inout) :: NeighborsKekule(ndim/2,6,2)

    integer(dp) :: n,h,i
    real(dp) :: a1r(2), a2r(2), newCoords(2), LatticeVectors(0:6,2)

    NeighborsKekule(:,:,:) = 0_dp

    do nlayer=1,nlayers

        if (RotateLayers(nlayer).eq.-1)then
            a1r = matmul(RotMatrix,a1)
            a2r = matmul(RotMatrix,a2)
        if (RotateLayers(nlayer).eq.1)then
            a1r = matmul(transpose(RotMatrix),a1)
            a2r = matmul(transpose(RotMatrix),a2)
        endif

        do n= ndim/nlayers*(nlayer-1)+1,ndim/nlayers*(nlayer-1)+ndim/nlayer/2

            NeighborsKekule(n - ndim/nlayers*(nlayer-1) + ndim/nlayers/2*(nlayer-1),1,1) = n
            NeighborsKekule(n - ndim/nlayers*(nlayer-1) + ndim/nlayers/2*(nlayer-1),1,2) = 0


            newCoords = Coords(n,1:2) + a1r/3.0 - 2*a2r/3.0
            do h = = ndim/nlayers*(nlayer-1)+ndim/nlayer/2+1, ndim/nlayers*(nlayer) 
                do i=0,6
                    if (norm2(Coords(h,1:2) + LatticeVectors(i,:) - newCoords).lt.a0*.1) then
                        NeighborsKekule(n - ndim/nlayers*(nlayer-1) + ndim/nlayers/2*(nlayer-1),2,1) = h
                        NeighborsKekule(n - ndim/nlayers*(nlayer-1) + ndim/nlayers/2*(nlayer-1),2,2) = i
                    endif
                enddo
            enddo
            
            newCoords = Coords(n,1:2) + a1r - a2r
            do h= ndim/nlayers*(nlayer-1)+1,ndim/nlayers*(nlayer-1)+ndim/nlayer/2
                do i=0,6
                    if (norm2(Coords(h,1:2) + LatticeVectors(i,:) - newCoords).lt.a0*.1) then
                        NeighborsKekule(n - ndim/nlayers*(nlayer-1) + ndim/nlayers/2*(nlayer-1),3,1) = h
                        NeighborsKekule(n - ndim/nlayers*(nlayer-1) + ndim/nlayers/2*(nlayer-1),3,2) = i
                    endif
                enddo
            enddo
            
            newCoords = Coords(n,1:2) + 4.0*a1r/3.0 - 2.0*a2r/3.0
            do h = = ndim/nlayers*(nlayer-1)+ndim/nlayer/2+1, ndim/nlayers*(nlayer) 
                do i=0,6
                    if (norm2(Coords(h,1:2) + LatticeVectors(i,:) - newCoords).lt.a0*.1) then
                        NeighborsKekule(n - ndim/nlayers*(nlayer-1) + ndim/nlayers/2*(nlayer-1),4,1) = h
                        NeighborsKekule(n - ndim/nlayers*(nlayer-1) + ndim/nlayers/2*(nlayer-1),4,2) = i
                    endif
                enddo
            enddo
            
            newCoords = Coords(n,1:2) + a1r
            do h= ndim/nlayers*(nlayer-1)+1,ndim/nlayers*(nlayer-1)+ndim/nlayer/2
                do i=0,6
                    if (norm2(Coords(h,1:2) + LatticeVectors(i,:) - newCoords).lt.a0*.1) then
                        NeighborsKekule(n - ndim/nlayers*(nlayer-1) + ndim/nlayers/2*(nlayer-1),5,1) = h
                        NeighborsKekule(n - ndim/nlayers*(nlayer-1) + ndim/nlayers/2*(nlayer-1),5,2) = i
                    endif
                enddo
            enddo
            
            newCoords = Coords(n,1:2) + a1r/3.0 + a2r/3.0
            do h = = ndim/nlayers*(nlayer-1)+ndim/nlayer/2+1, ndim/nlayers*(nlayer) 
                do i=0,6
                    if (norm2(Coords(h,1:2) + LatticeVectors(i,:) - newCoords).lt.a0*.1) then
                        NeighborsKekule(n - ndim/nlayers*(nlayer-1) + ndim/nlayers/2*(nlayer-1),6,1) = h
                        NeighborsKekule(n - ndim/nlayers*(nlayer-1) + ndim/nlayers/2*(nlayer-1),6,2) = i
                    endif
                enddo
            enddo

            if(product(NeighborsKekule(n - ndim/nlayers*(nlayer-1) + ndim/nlayers/2*(nlayer-1),:,1)).eq.0)then
                NeighborsKekule(n-ndim/4,2:6,:) = 0
                write(*,*) 'ERROR KEKULENEIGHBORS'
            endif

        enddo
    enddo

end subroutine KekuleNeighbors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine IntraSubIntraVal(ValleyPol,TotalDens,ndim,NearestNeighborInd,NearestNeighborCell,numNeighborCells,nUnitCell_1,nUnitCell_2,FockIndexToVector,zFock)

    integer(dp), intent(in) :: ndim, numNeighborCells
    integer(dp), intent(in) :: NearestNeighborInd(ndim,3),NearestNeighborCell(ndim,3), nUnitCell_1(numNeighborCells), nUnitCell_2(numNeighborCells), FockIndexToVector(0:6,2)
    real(dp), intent(out) :: ValleyPol(ndim), TotalDens(ndim)
    complex(dp), intent(in) :: zFock(ndim,ndim,numNeighborCells)
    
    integer(dp) :: n,m, FockIndex, nlayer
    complex(dp) :: z1,z2,z3, temp(ndim)
    
    ValleyPol(:) = 0.0_dp
    TotalDens(:) = 0.0_dp
    temp(:) = cmplx(0.0_dp,0.0_dp,dp)

    do n=1,ndim

        if(numI.eq.0)then
            z1 =  zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),1)
            z2 = zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),1)
            z3 = zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),1)
        
        else

            call VectorToFockIndex(FockIndex, -FockIndexToVector(NearestNeighborCell(n,1),1) + FockIndexToVector(NearestNeighborCell(n,3),1),&
            -FockIndexToVector(NearestNeighborCell(n,1),2) + FockIndexToVector(NearestNeighborCell(n,3),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, +FockIndexToVector(NearestNeighborCell(n,1),1) - FockIndexToVector(NearestNeighborCell(n,3),1),&
                +FockIndexToVector(NearestNeighborCell(n,1),2) - FockIndexToVector(NearestNeighborCell(n,3),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z1 = conjg(zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,1),FockIndex))
            else
                
                z1 = zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),FockIndex)
            endif


            call VectorToFockIndex(FockIndex, -FockIndexToVector(NearestNeighborCell(n,3),1) + FockIndexToVector(NearestNeighborCell(n,2),1),&
            -FockIndexToVector(NearestNeighborCell(n,3),2) + FockIndexToVector(NearestNeighborCell(n,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, +FockIndexToVector(NearestNeighborCell(n,3),1) - FockIndexToVector(NearestNeighborCell(n,2),1),&
                +FockIndexToVector(NearestNeighborCell(n,3),2) - FockIndexToVector(NearestNeighborCell(n,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z2 = conjg(zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,3),FockIndex))
            else
                
                z2 = zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),FockIndex)
            endif
      
            call VectorToFockIndex(FockIndex, -FockIndexToVector(NearestNeighborCell(n,2),1) + FockIndexToVector(NearestNeighborCell(n,1),1),&
            -FockIndexToVector(NearestNeighborCell(n,2),2) + FockIndexToVector(NearestNeighborCell(n,1),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, +FockIndexToVector(NearestNeighborCell(n,2),1) - FockIndexToVector(NearestNeighborCell(n,1),1),&
                +FockIndexToVector(NearestNeighborCell(n,2),2) - FockIndexToVector(NearestNeighborCell(n,1),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z3 = conjg(zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,2),FockIndex))
            else
                
                z3 = zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),FockIndex)
            endif

        endif

        temp(n) = z1 + z2 + z3

    enddo
    
    do nlayer=1,nlayers
        do n=1,ndim/nlayers/2

            m = ndim/nlayers*(nlayer-1) + n
            ValleyPol(m) = 2.0_dp/3.0_dp/sqrt(3.0_dp)*aimag(temp(m))
            TotalDens(m) = -2.0_dp/3.0_dp*real(temp(m))

            m = ndim/nlayers*(nlayer-1) + ndim/nlayers/2 + n
            ValleyPol(m) = -2.0_dp/3.0_dp/sqrt(3.0_dp)*aimag(temp(m))
            TotalDens(m) = -2.0_dp/3.0_dp*real(temp(m))

        enddo
    enddo


    return
end subroutine IntraSubIntraVal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine IntraSubInterVal(fkpk,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,NearestNeighborInd,NearestNeighborCell,FockIndexToVector,zFock)

    integer(dp), intent(in) :: ndim, numNeighborCells
    integer(dp), intent(in) :: NearestNeighborInd(ndim,3),NearestNeighborCell(ndim,3), FockIndexToVector(0:6,2), nUnitCell_1(numNeighborCells), nUnitCell_2(numNeighborCells)
    complex(dp), intent(out) :: fkpk(ndim)
    complex(dp), intent(in) :: zFock(ndim,ndim,numNeighborCells)
    
    integer(dp) :: n,m, FockIndex
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

            call VectorToFockIndex(FockIndex, -FockIndexToVector(NearestNeighborCell(n,1),1) + FockIndexToVector(NearestNeighborCell(n,3),1),&
            -FockIndexToVector(NearestNeighborCell(n,1),2) + FockIndexToVector(NearestNeighborCell(n,3),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, +FockIndexToVector(NearestNeighborCell(n,1),1) - FockIndexToVector(NearestNeighborCell(n,3),1),&
                +FockIndexToVector(NearestNeighborCell(n,1),2) - FockIndexToVector(NearestNeighborCell(n,3),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z1 = conjg(zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,1),FockIndex))
            else
                
                z1 = zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),FockIndex)
            endif
            !write(*,*) n,FockIndex,z1,zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),1)

            call VectorToFockIndex(FockIndex, -FockIndexToVector(NearestNeighborCell(n,3),1) + FockIndexToVector(NearestNeighborCell(n,2),1),&
            -FockIndexToVector(NearestNeighborCell(n,3),2) + FockIndexToVector(NearestNeighborCell(n,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, +FockIndexToVector(NearestNeighborCell(n,3),1) - FockIndexToVector(NearestNeighborCell(n,2),1),&
                +FockIndexToVector(NearestNeighborCell(n,3),2) - FockIndexToVector(NearestNeighborCell(n,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z2 = conjg(zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,3),FockIndex))
            else
                
                z2 = zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),FockIndex)
            endif
            !write(*,*) n,FockIndex,z2,zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),1)

            call VectorToFockIndex(FockIndex, -FockIndexToVector(NearestNeighborCell(n,2),1) + FockIndexToVector(NearestNeighborCell(n,1),1),&
            -FockIndexToVector(NearestNeighborCell(n,2),2) + FockIndexToVector(NearestNeighborCell(n,1),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, +FockIndexToVector(NearestNeighborCell(n,2),1) - FockIndexToVector(NearestNeighborCell(n,1),1),&
                +FockIndexToVector(NearestNeighborCell(n,2),2) - FockIndexToVector(NearestNeighborCell(n,1),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z3 = conjg(zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,2),FockIndex))
            else
                
                z3 = zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),FockIndex)
            endif
            !write(*,*) n,FockIndex,z3,zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),1)

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

                call VectorToFockIndex(FockIndex, -FockIndexToVector(NearestNeighborCell(n,1),1) + FockIndexToVector(NearestNeighborCell(n,3),1),&
            -FockIndexToVector(NearestNeighborCell(n,1),2) + FockIndexToVector(NearestNeighborCell(n,3),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, +FockIndexToVector(NearestNeighborCell(n,1),1) - FockIndexToVector(NearestNeighborCell(n,3),1),&
                +FockIndexToVector(NearestNeighborCell(n,1),2) - FockIndexToVector(NearestNeighborCell(n,3),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z1 = conjg(zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,1),FockIndex))
            else
                
                z1 = zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),FockIndex)
            endif
            !write(*,*) n,FockIndex,z1,zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),1)

            call VectorToFockIndex(FockIndex, -FockIndexToVector(NearestNeighborCell(n,3),1) + FockIndexToVector(NearestNeighborCell(n,2),1),&
            -FockIndexToVector(NearestNeighborCell(n,3),2) + FockIndexToVector(NearestNeighborCell(n,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, +FockIndexToVector(NearestNeighborCell(n,3),1) - FockIndexToVector(NearestNeighborCell(n,2),1),&
                +FockIndexToVector(NearestNeighborCell(n,3),2) - FockIndexToVector(NearestNeighborCell(n,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z2 = conjg(zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,3),FockIndex))
            else
                
                z2 = zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),FockIndex)
            endif
            !write(*,*) n,FockIndex,z2,zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),1)

            call VectorToFockIndex(FockIndex, -FockIndexToVector(NearestNeighborCell(n,2),1) + FockIndexToVector(NearestNeighborCell(n,1),1),&
            -FockIndexToVector(NearestNeighborCell(n,2),2) + FockIndexToVector(NearestNeighborCell(n,1),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, +FockIndexToVector(NearestNeighborCell(n,2),1) - FockIndexToVector(NearestNeighborCell(n,1),1),&
                +FockIndexToVector(NearestNeighborCell(n,2),2) - FockIndexToVector(NearestNeighborCell(n,1),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z3 = conjg(zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,2),FockIndex))
            else
                
                z3 = zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),FockIndex)
            endif
            !write(*,*) n,FockIndex,z3,zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),1)

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
            call VectorToFockIndex(FockIndex, -FockIndexToVector(NearestNeighborCell(n,1),1) + FockIndexToVector(NearestNeighborCell(n,3),1),&
            -FockIndexToVector(NearestNeighborCell(n,1),2) + FockIndexToVector(NearestNeighborCell(n,3),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, +FockIndexToVector(NearestNeighborCell(n,1),1) - FockIndexToVector(NearestNeighborCell(n,3),1),&
                +FockIndexToVector(NearestNeighborCell(n,1),2) - FockIndexToVector(NearestNeighborCell(n,3),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z1 = conjg(zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,1),FockIndex))
            else
                
                z1 = zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),FockIndex)
            endif
            !write(*,*) n,FockIndex,z1,zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),1)

            call VectorToFockIndex(FockIndex, -FockIndexToVector(NearestNeighborCell(n,3),1) + FockIndexToVector(NearestNeighborCell(n,2),1),&
            -FockIndexToVector(NearestNeighborCell(n,3),2) + FockIndexToVector(NearestNeighborCell(n,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, +FockIndexToVector(NearestNeighborCell(n,3),1) - FockIndexToVector(NearestNeighborCell(n,2),1),&
                +FockIndexToVector(NearestNeighborCell(n,3),2) - FockIndexToVector(NearestNeighborCell(n,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z2 = conjg(zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,3),FockIndex))
            else
                
                z2 = zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),FockIndex)
            endif
            !write(*,*) n,FockIndex,z2,zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),1)

            call VectorToFockIndex(FockIndex, -FockIndexToVector(NearestNeighborCell(n,2),1) + FockIndexToVector(NearestNeighborCell(n,1),1),&
            -FockIndexToVector(NearestNeighborCell(n,2),2) + FockIndexToVector(NearestNeighborCell(n,1),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, +FockIndexToVector(NearestNeighborCell(n,2),1) - FockIndexToVector(NearestNeighborCell(n,1),1),&
                +FockIndexToVector(NearestNeighborCell(n,2),2) - FockIndexToVector(NearestNeighborCell(n,1),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z3 = conjg(zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,2),FockIndex))
            else
                
                z3 = zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),FockIndex)
            endif
            !write(*,*) n,FockIndex,z3,zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),1)

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

            call VectorToFockIndex(FockIndex, -FockIndexToVector(NearestNeighborCell(n,1),1) + FockIndexToVector(NearestNeighborCell(n,3),1),&
            -FockIndexToVector(NearestNeighborCell(n,1),2) + FockIndexToVector(NearestNeighborCell(n,3),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, +FockIndexToVector(NearestNeighborCell(n,1),1) - FockIndexToVector(NearestNeighborCell(n,3),1),&
                +FockIndexToVector(NearestNeighborCell(n,1),2) - FockIndexToVector(NearestNeighborCell(n,3),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z1 = conjg(zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,1),FockIndex))
            else
                
                z1 = zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),FockIndex)
            endif
            !write(*,*) n,FockIndex,z1,zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,3),1)

            call VectorToFockIndex(FockIndex, -FockIndexToVector(NearestNeighborCell(n,3),1) + FockIndexToVector(NearestNeighborCell(n,2),1),&
            -FockIndexToVector(NearestNeighborCell(n,3),2) + FockIndexToVector(NearestNeighborCell(n,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, +FockIndexToVector(NearestNeighborCell(n,3),1) - FockIndexToVector(NearestNeighborCell(n,2),1),&
                +FockIndexToVector(NearestNeighborCell(n,3),2) - FockIndexToVector(NearestNeighborCell(n,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z2 = conjg(zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,3),FockIndex))
            else
                
                z2 = zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),FockIndex)
            endif
            !write(*,*) n,FockIndex,z2,zFock(NearestNeighborInd(n,3),NearestNeighborInd(n,2),1)

            call VectorToFockIndex(FockIndex, -FockIndexToVector(NearestNeighborCell(n,2),1) + FockIndexToVector(NearestNeighborCell(n,1),1),&
            -FockIndexToVector(NearestNeighborCell(n,2),2) + FockIndexToVector(NearestNeighborCell(n,1),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if (FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, +FockIndexToVector(NearestNeighborCell(n,2),1) - FockIndexToVector(NearestNeighborCell(n,1),1),&
                +FockIndexToVector(NearestNeighborCell(n,2),2) - FockIndexToVector(NearestNeighborCell(n,1),2), nUnitCell_1, nUnitCell_2, numNeighborCells)
                
                z3 = conjg(zFock(NearestNeighborInd(n,1),NearestNeighborInd(n,2),FockIndex))
            else
                
                z3 = zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),FockIndex)
            endif
            !write(*,*) n,FockIndex,z3,zFock(NearestNeighborInd(n,2),NearestNeighborInd(n,1),1)
            
            ztemp = z1*z2Pi3 + z2 + z3/z2Pi3
        endif

        fkpk(n) = ztemp/3.0_dp
    enddo
        
    return
end subroutine IntraSubInterVal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InterSubInterVal(fkakpb,fkbkpa,NeighborsKekule,KekuleLoops,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,FockIndexToVector,zFock)

    integer(dp), intent(in) :: ndim, numNeighborCells
    integer(dp), intent(in) :: nUnitCell_1(numNeighborCells), nUnitCell_2(numNeighborCells), FockIndexToVector(0:6,2)
    integer(dp), intent(in) :: NeighborsKekule(ndim/2,6,2), KekuleLoops(ndim/2,2,3)
    complex(dp), intent(out) :: fkakpb(ndim/2), fkbkpa(ndim/2)
    complex(dp), intent(in) :: zFock(ndim,ndim,numNeighborCells)
    
    integer(dp) :: n, m, FockIndex
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

            call VectorToFockIndex(FockIndex, FockIndexToVector(NeighborsKekule(n,1,2),1) - FockIndexToVector(NeighborsKekule(n,2,2),1),&
                FockIndexToVector(NeighborsKekule(n,1,2),2) - FockIndexToVector(NeighborsKekule(n,2,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, -FockIndexToVector(NeighborsKekule(n,1,2),1) + FockIndexToVector(NeighborsKekule(n,2,2),1),&
                -FockIndexToVector(NeighborsKekule(n,1,2),2) + FockIndexToVector(NeighborsKekule(n,2,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z1 = conjg(zFock(NeighborsKekule(n,2,1),NeighborsKekule(n,1,1),FockIndex))
            else

                z1 = zFock(NeighborsKekule(n,1,1),NeighborsKekule(n,2,1),FockIndex)
            endif

            call VectorToFockIndex(FockIndex, FockIndexToVector(NeighborsKekule(n,2,2),1) - FockIndexToVector(NeighborsKekule(n,3,2),1),&
                FockIndexToVector(NeighborsKekule(n,2,2),2) - FockIndexToVector(NeighborsKekule(n,3,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, -FockIndexToVector(NeighborsKekule(n,2,2),1) + FockIndexToVector(NeighborsKekule(n,3,2),1),&
                -FockIndexToVector(NeighborsKekule(n,2,2),2) + FockIndexToVector(NeighborsKekule(n,3,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z2 = conjg(zFock(NeighborsKekule(n,3,1),NeighborsKekule(n,2,1),FockIndex))
            else

                z2 = zFock(NeighborsKekule(n,2,1),NeighborsKekule(n,3,1),FockIndex)
            endif
            
            call VectorToFockIndex(FockIndex, FockIndexToVector(NeighborsKekule(n,3,2),1) - FockIndexToVector(NeighborsKekule(n,4,2),1),&
                FockIndexToVector(NeighborsKekule(n,3,2),2) - FockIndexToVector(NeighborsKekule(n,4,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, -FockIndexToVector(NeighborsKekule(n,3,2),1) + FockIndexToVector(NeighborsKekule(n,4,2),1),&
                -FockIndexToVector(NeighborsKekule(n,3,2),2) + FockIndexToVector(NeighborsKekule(n,4,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z3 = conjg(zFock(NeighborsKekule(n,4,1),NeighborsKekule(n,3,1),FockIndex))
            else

                z3 = zFock(NeighborsKekule(n,3,1),NeighborsKekule(n,4,1),FockIndex)
            endif
            
            call VectorToFockIndex(FockIndex, FockIndexToVector(NeighborsKekule(n,4,2),1) - FockIndexToVector(NeighborsKekule(n,5,2),1),&
            FockIndexToVector(NeighborsKekule(n,4,2),2) - FockIndexToVector(NeighborsKekule(n,5,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, -FockIndexToVector(NeighborsKekule(n,4,2),1) + FockIndexToVector(NeighborsKekule(n,5,2),1),&
                -FockIndexToVector(NeighborsKekule(n,4,2),2) + FockIndexToVector(NeighborsKekule(n,5,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z4 = conjg(zFock(NeighborsKekule(n,5,1),NeighborsKekule(n,4,1),FockIndex))
            else

                z4 = zFock(NeighborsKekule(n,4,1),NeighborsKekule(n,5,1),FockIndex)
            endif
            
            call VectorToFockIndex(FockIndex, FockIndexToVector(NeighborsKekule(n,5,2),1) - FockIndexToVector(NeighborsKekule(n,6,2),1),&
            FockIndexToVector(NeighborsKekule(n,5,2),2) - FockIndexToVector(NeighborsKekule(n,6,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, -FockIndexToVector(NeighborsKekule(n,5,2),1) + FockIndexToVector(NeighborsKekule(n,6,2),1),&
                -FockIndexToVector(NeighborsKekule(n,5,2),2) + FockIndexToVector(NeighborsKekule(n,6,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z5 = conjg(zFock(NeighborsKekule(n,6,1),NeighborsKekule(n,5,1),FockIndex))
            else

                z5 = zFock(NeighborsKekule(n,5,1),NeighborsKekule(n,6,1),FockIndex)
            endif
            
            call VectorToFockIndex(FockIndex, FockIndexToVector(NeighborsKekule(n,6,2),1) - FockIndexToVector(NeighborsKekule(n,1,2),1),&
            FockIndexToVector(NeighborsKekule(n,6,2),2) - FockIndexToVector(NeighborsKekule(n,1,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, -FockIndexToVector(NeighborsKekule(n,6,2),1) + FockIndexToVector(NeighborsKekule(n,1,2),1),&
                -FockIndexToVector(NeighborsKekule(n,6,2),2) + FockIndexToVector(NeighborsKekule(n,1,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z6 = conjg(zFock(NeighborsKekule(n,1,1),NeighborsKekule(n,6,1),FockIndex))
            else

                z6 = zFock(NeighborsKekule(n,6,1),NeighborsKekule(n,1,1),FockIndex)
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

subroutine InterSubIntraVal_old(fkakb,fkpakpb,NeighborsKekule,KekuleLoops,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,FockIndexToVector,zFock)

    integer(dp), intent(in) :: ndim,numNeighborCells
    integer(dp), intent(in) :: nUnitCell_1(numNeighborCells), nUnitCell_2(numNeighborCells)
    integer(dp), intent(in) :: NeighborsKekule(ndim/2,6,2), KekuleLoops(ndim/2,2,3), FockIndexToVector(0:6,2)
    complex(dp), intent(out) :: fkakb(ndim/2), fkpakpb(ndim/2)
    complex(dp), intent(in) :: zFock(ndim,ndim,numNeighborCells)
    
    integer(dp) :: n,m,FockIndex
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

            call VectorToFockIndex(FockIndex, FockIndexToVector(NeighborsKekule(n,1,2),1) - FockIndexToVector(NeighborsKekule(n,2,2),1),&
                FockIndexToVector(NeighborsKekule(n,1,2),2) - FockIndexToVector(NeighborsKekule(n,2,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, -FockIndexToVector(NeighborsKekule(n,1,2),1) + FockIndexToVector(NeighborsKekule(n,2,2),1),&
                -FockIndexToVector(NeighborsKekule(n,1,2),2) + FockIndexToVector(NeighborsKekule(n,2,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z1 = conjg(zFock(NeighborsKekule(n,2,1),NeighborsKekule(n,1,1),FockIndex))
            else

                z1 = zFock(NeighborsKekule(n,1,1),NeighborsKekule(n,2,1),FockIndex)
            endif

            call VectorToFockIndex(FockIndex, FockIndexToVector(NeighborsKekule(n,2,2),1) - FockIndexToVector(NeighborsKekule(n,3,2),1),&
                FockIndexToVector(NeighborsKekule(n,2,2),2) - FockIndexToVector(NeighborsKekule(n,3,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, -FockIndexToVector(NeighborsKekule(n,2,2),1) + FockIndexToVector(NeighborsKekule(n,3,2),1),&
                -FockIndexToVector(NeighborsKekule(n,2,2),2) + FockIndexToVector(NeighborsKekule(n,3,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z2 = conjg(zFock(NeighborsKekule(n,3,1),NeighborsKekule(n,2,1),FockIndex))
            else

                z2 = zFock(NeighborsKekule(n,2,1),NeighborsKekule(n,3,1),FockIndex)
            endif

            
            call VectorToFockIndex(FockIndex, FockIndexToVector(NeighborsKekule(n,3,2),1) - FockIndexToVector(NeighborsKekule(n,4,2),1),&
                FockIndexToVector(NeighborsKekule(n,3,2),2) - FockIndexToVector(NeighborsKekule(n,4,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, -FockIndexToVector(NeighborsKekule(n,3,2),1) + FockIndexToVector(NeighborsKekule(n,4,2),1),&
                -FockIndexToVector(NeighborsKekule(n,3,2),2) + FockIndexToVector(NeighborsKekule(n,4,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z3 = conjg(zFock(NeighborsKekule(n,4,1),NeighborsKekule(n,3,1),FockIndex))
            else

                z3 = zFock(NeighborsKekule(n,3,1),NeighborsKekule(n,4,1),FockIndex)
            endif

            
            call VectorToFockIndex(FockIndex, FockIndexToVector(NeighborsKekule(n,4,2),1) - FockIndexToVector(NeighborsKekule(n,5,2),1),&
            FockIndexToVector(NeighborsKekule(n,4,2),2) - FockIndexToVector(NeighborsKekule(n,5,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, -FockIndexToVector(NeighborsKekule(n,4,2),1) + FockIndexToVector(NeighborsKekule(n,5,2),1),&
                -FockIndexToVector(NeighborsKekule(n,4,2),2) + FockIndexToVector(NeighborsKekule(n,5,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z4 = conjg(zFock(NeighborsKekule(n,5,1),NeighborsKekule(n,4,1),FockIndex))
            else

                z4 = zFock(NeighborsKekule(n,4,1),NeighborsKekule(n,5,1),FockIndex)
            endif

            
            call VectorToFockIndex(FockIndex, FockIndexToVector(NeighborsKekule(n,5,2),1) - FockIndexToVector(NeighborsKekule(n,6,2),1),&
            FockIndexToVector(NeighborsKekule(n,5,2),2) - FockIndexToVector(NeighborsKekule(n,6,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, -FockIndexToVector(NeighborsKekule(n,5,2),1) + FockIndexToVector(NeighborsKekule(n,6,2),1),&
                -FockIndexToVector(NeighborsKekule(n,5,2),2) + FockIndexToVector(NeighborsKekule(n,6,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z5 = conjg(zFock(NeighborsKekule(n,6,1),NeighborsKekule(n,5,1),FockIndex))
            else

                z5 = zFock(NeighborsKekule(n,5,1),NeighborsKekule(n,6,1),FockIndex)
            endif

            
            call VectorToFockIndex(FockIndex, FockIndexToVector(NeighborsKekule(n,6,2),1) - FockIndexToVector(NeighborsKekule(n,1,2),1),&
            FockIndexToVector(NeighborsKekule(n,6,2),2) - FockIndexToVector(NeighborsKekule(n,1,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, -FockIndexToVector(NeighborsKekule(n,6,2),1) + FockIndexToVector(NeighborsKekule(n,1,2),1),&
                -FockIndexToVector(NeighborsKekule(n,6,2),2) + FockIndexToVector(NeighborsKekule(n,1,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z6 = conjg(zFock(NeighborsKekule(n,1,1),NeighborsKekule(n,6,1),FockIndex))
            else

                z6 = zFock(NeighborsKekule(n,6,1),NeighborsKekule(n,1,1),FockIndex)
            endif

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

subroutine InterSubIntraVal(fkakb,fkpakpb,NeighborsKekule,KekuleLoops,ndim,numNeighborCells,nUnitCell_1,nUnitCell_2,FockIndexToVector,zFock)

    integer(dp), intent(in) :: ndim,numNeighborCells
    integer(dp), intent(in) :: nUnitCell_1(numNeighborCells), nUnitCell_2(numNeighborCells)
    integer(dp), intent(in) :: NeighborsKekule(ndim/2,6,2), KekuleLoops(ndim/2,2,3), FockIndexToVector(0:6,2)
    complex(dp), intent(out) :: fkakb(ndim/2), fkpakpb(ndim/2)
    complex(dp), intent(in) :: zFock(ndim,ndim,numNeighborCells)
    
    integer(dp) :: n,m,FockIndex
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

            call VectorToFockIndex(FockIndex, FockIndexToVector(NeighborsKekule(n,1,2),1) - FockIndexToVector(NeighborsKekule(n,2,2),1),&
                FockIndexToVector(NeighborsKekule(n,1,2),2) - FockIndexToVector(NeighborsKekule(n,2,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, -FockIndexToVector(NeighborsKekule(n,1,2),1) + FockIndexToVector(NeighborsKekule(n,2,2),1),&
                -FockIndexToVector(NeighborsKekule(n,1,2),2) + FockIndexToVector(NeighborsKekule(n,2,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z1 = conjg(zFock(NeighborsKekule(n,2,1),NeighborsKekule(n,1,1),FockIndex))
            else

                z1 = zFock(NeighborsKekule(n,1,1),NeighborsKekule(n,2,1),FockIndex)
            endif
            ! write(*,*) FockIndex

            call VectorToFockIndex(FockIndex, FockIndexToVector(NeighborsKekule(n,2,2),1) - FockIndexToVector(NeighborsKekule(n,3,2),1),&
                FockIndexToVector(NeighborsKekule(n,2,2),2) - FockIndexToVector(NeighborsKekule(n,3,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, -FockIndexToVector(NeighborsKekule(n,2,2),1) + FockIndexToVector(NeighborsKekule(n,3,2),1),&
                -FockIndexToVector(NeighborsKekule(n,2,2),2) + FockIndexToVector(NeighborsKekule(n,3,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z2 = conjg(zFock(NeighborsKekule(n,3,1),NeighborsKekule(n,2,1),FockIndex))
            else

                z2 = zFock(NeighborsKekule(n,2,1),NeighborsKekule(n,3,1),FockIndex)
            endif
            ! write(*,*) FockIndex

            
            call VectorToFockIndex(FockIndex, FockIndexToVector(NeighborsKekule(n,3,2),1) - FockIndexToVector(NeighborsKekule(n,4,2),1),&
                FockIndexToVector(NeighborsKekule(n,3,2),2) - FockIndexToVector(NeighborsKekule(n,4,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, -FockIndexToVector(NeighborsKekule(n,3,2),1) + FockIndexToVector(NeighborsKekule(n,4,2),1),&
                -FockIndexToVector(NeighborsKekule(n,3,2),2) + FockIndexToVector(NeighborsKekule(n,4,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z3 = conjg(zFock(NeighborsKekule(n,4,1),NeighborsKekule(n,3,1),FockIndex))
            else

                z3 = zFock(NeighborsKekule(n,3,1),NeighborsKekule(n,4,1),FockIndex)
            endif
            ! write(*,*) FockIndex

            
            call VectorToFockIndex(FockIndex, FockIndexToVector(NeighborsKekule(n,4,2),1) - FockIndexToVector(NeighborsKekule(n,5,2),1),&
            FockIndexToVector(NeighborsKekule(n,4,2),2) - FockIndexToVector(NeighborsKekule(n,5,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, -FockIndexToVector(NeighborsKekule(n,4,2),1) + FockIndexToVector(NeighborsKekule(n,5,2),1),&
                -FockIndexToVector(NeighborsKekule(n,4,2),2) + FockIndexToVector(NeighborsKekule(n,5,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z4 = conjg(zFock(NeighborsKekule(n,5,1),NeighborsKekule(n,4,1),FockIndex))
            else

                z4 = zFock(NeighborsKekule(n,4,1),NeighborsKekule(n,5,1),FockIndex)
            endif
            ! write(*,*) FockIndex

            
            call VectorToFockIndex(FockIndex, FockIndexToVector(NeighborsKekule(n,5,2),1) - FockIndexToVector(NeighborsKekule(n,6,2),1),&
            FockIndexToVector(NeighborsKekule(n,5,2),2) - FockIndexToVector(NeighborsKekule(n,6,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, -FockIndexToVector(NeighborsKekule(n,5,2),1) + FockIndexToVector(NeighborsKekule(n,6,2),1),&
                -FockIndexToVector(NeighborsKekule(n,5,2),2) + FockIndexToVector(NeighborsKekule(n,6,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z5 = conjg(zFock(NeighborsKekule(n,6,1),NeighborsKekule(n,5,1),FockIndex))
            else

                z5 = zFock(NeighborsKekule(n,5,1),NeighborsKekule(n,6,1),FockIndex)
            endif
            ! write(*,*) FockIndex

            
            call VectorToFockIndex(FockIndex, FockIndexToVector(NeighborsKekule(n,6,2),1) - FockIndexToVector(NeighborsKekule(n,1,2),1),&
            FockIndexToVector(NeighborsKekule(n,6,2),2) - FockIndexToVector(NeighborsKekule(n,1,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells)  
            if(FockIndex.eq.-1_dp)then
                call VectorToFockIndex(FockIndex, -FockIndexToVector(NeighborsKekule(n,6,2),1) + FockIndexToVector(NeighborsKekule(n,1,2),1),&
                -FockIndexToVector(NeighborsKekule(n,6,2),2) + FockIndexToVector(NeighborsKekule(n,1,2),2), nUnitCell_1, nUnitCell_2, numNeighborCells) 

                z6 = conjg(zFock(NeighborsKekule(n,1,1),NeighborsKekule(n,6,1),FockIndex))
            else

                z6 = zFock(NeighborsKekule(n,6,1),NeighborsKekule(n,1,1),FockIndex)
            endif
            ! write(*,*) FockIndex

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
