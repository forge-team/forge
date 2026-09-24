! Shared geometry and packed Fock input, following the former utils/Main_BandStructure.f90.
! Compile with the producer's ../Setup.f90, not a separate set of model defaults.
module PostProcessingInput
    use iso_fortran_env, only: int64, file_storage_size, error_unit
    use Setup, only: dp, ndim, numS, nlayers, ntheta, nrelax, TBFunction, Delta, &
        a1, a2, pi, dirFock, InitParameters, parametersRun => parameters
    use Geometry, only: WignerSeitzCell, LatticeRelaxationKoshino, LatticeRelaxationCarr
    implicit none
    private
    public :: BuildGeometry, ReadFock, InputParameters, ValidateResponseModel
contains

subroutine BuildGeometry(Coords, t1, t2, t3, g1, g12, RotMatrix)
    real(dp), intent(out) :: Coords(ndim,3), t1(2), t2(2), t3(2), g1(2), g12(2), RotMatrix(2,2)
    real(dp) :: aMoire, cs, sn
    integer :: n

    aMoire = 3.0_dp*ntheta**2 + 3.0_dp*ntheta + 1.0_dp
    g1 = (4.0_dp*pi/3.0_dp)/aMoire*(real(3*ntheta+1,dp)*a1+a2)
    g12 = (4.0_dp*pi/3.0_dp)/aMoire*(real(3*ntheta+2,dp)*a2-a1)
    cs = 1.0_dp-1.0_dp/(2.0_dp*aMoire)
    sn = sqrt(1.0_dp-cs**2)
    RotMatrix = reshape([cos(0.5_dp*acos(cs)),-sin(0.5_dp*acos(cs)), &
        sin(0.5_dp*acos(cs)),cos(0.5_dp*acos(cs))],[2,2])
    g1 = matmul(RotMatrix,g1)
    g12 = matmul(RotMatrix,g12)

    t1 = real(ntheta,dp)*a1 + real(ntheta+1,dp)*a2
    t2 = real(-ntheta-1,dp)*a1 + real(2*ntheta+1,dp)*a2
    t3 = t2-t1
    call WignerSeitzCell(Coords,t1,t2,cs,sn)
    do n=1,ndim
        Coords(n,1:2) = matmul(RotMatrix,Coords(n,1:2))
    enddo
    t1 = matmul(RotMatrix,t1)
    t2 = matmul(RotMatrix,t2)
    t3 = t2-t1

    if(nrelax == 1 .and. nlayers == 2) call LatticeRelaxationKoshino(Coords,g1,g12)
    if(nrelax == 2 .and. nlayers == 2) call LatticeRelaxationCarr(Coords,g1,g12)
end subroutine BuildGeometry

subroutine InputParameters(suffix)
    character(*), intent(out) :: suffix
    call InitParameters()
    if(len_trim(parametersRun) > len(suffix)) error stop 'Input suffix buffer is too short'
    ! The postprocessor consumes the producer's output state, not its restart guess.
    suffix = trim(parametersRun)
end subroutine InputParameters

subroutine ReadFock(zFock, numNeighborCells, suffix)
    integer(dp), intent(in) :: numNeighborCells
    complex(dp), intent(out) :: zFock(ndim,ndim,numNeighborCells,numS)
    character(*), intent(in) :: suffix
    complex(dp) :: zinput
    integer(dp) :: i, j, m, nspin
    integer :: unit, stat, scalar_units
    integer(int64) :: rcc, records, expected_size, actual_size, scalar_storage
    character(32) :: spin
    character(512) :: message
    character(:), allocatable :: filename

    zFock = cmplx(0.0_dp,0.0_dp,dp)
    records = int(ndim,int64)*(int(ndim,int64)+1)/2 + &
        (int(numNeighborCells,int64)-1)*int(ndim,int64)**2
    ! Keep the same RECL as upstream. IOLENGTH accounts for Intel's optional
    ! four-byte units; do not silently reinterpret an existing padded file.
    inquire(iolength=scalar_units) zinput
    scalar_storage = storage_size(zinput)/file_storage_size
    expected_size = records*(2*dp)*scalar_storage/scalar_units

    do nspin=1,numS
        write(spin,'(I0)') nspin
        filename = trim(dirFock)//'Fock-nspin'//trim(spin)//trim(suffix)
        open(newunit=unit, file=filename, form='unformatted', status='old', &
            action='read', access='direct', recl=2*dp, iostat=stat, iomsg=message)
        if(stat /= 0) then
            write(error_unit,'(A)') 'Cannot read '//filename//': '//trim(message)
            error stop 1
        endif
        inquire(unit=unit, size=actual_size)
        if(actual_size /= expected_size) then
            write(error_unit,'(A)') 'Fock file size does not match Setup and compiler record units: '//filename
            write(error_unit,*) 'Expected/actual storage units:', expected_size, actual_size
            close(unit)
            error stop 1
        endif

        rcc = 0
        do i=1,ndim
            do j=1,i-1
                call NextRecord()
                zFock(i,j,1,nspin) = zinput
                zFock(j,i,1,nspin) = conjg(zinput)
            enddo
            call NextRecord()
            zFock(i,i,1,nspin) = zinput
        enddo
        do m=2,numNeighborCells
            do i=1,ndim
                do j=1,ndim
                    call NextRecord()
                    zFock(i,j,m,nspin) = zinput
                enddo
            enddo
        enddo
        close(unit)
    enddo
contains
    subroutine NextRecord()
        rcc = rcc+1
        read(unit,rec=rcc,iostat=stat,iomsg=message) zinput
        if(stat /= 0) then
            write(error_unit,'(A)') 'Fock read failed in '//filename//': '//trim(message)
            write(error_unit,*) 'Record:', rcc
            error stop 1
        endif
    end subroutine NextRecord
end subroutine ReadFock

subroutine ValidateResponseModel()
    ! The imported observable loops calculate one bilayer spin sector. Their
    ! legacy response kernels do not implement the dev Wannier model or bias.
    if(nlayers /= 2) error stop 'Legacy response kernels require nlayers=2'
    if(numS /= 1) error stop 'Legacy response loops currently require numS=1'
    if(TBFunction /= 1) error stop 'Legacy response kernels require TBFunction=1'
    if(abs(Delta) > 0.0_dp) error stop 'Legacy response kernels do not implement layer bias'
end subroutine ValidateResponseModel

end module PostProcessingInput
