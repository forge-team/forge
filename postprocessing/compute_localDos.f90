! Spatial local DOS from the current FORGE state. Based on LocalDos.tar:
! plotDOS_MIT.f90 and the LocalDOS_CB/LocalDOS_VB kernels in Spectral.f90.
! Build with `make -C postprocessing`; run from the FORGE root.
program compute_localDos
    use iso_fortran_env, only: error_unit, iostat_end
    use ieee_arithmetic, only: ieee_is_finite
    use Setup, only: dp, ndim, numS, numI, numb, nLower, nUpper, nfilling, alpha, alphaH, U, pi, &
        inputDir => dir
    use Geometry, only: NumberNeighborCells, OrderNeighborCells
    use Hamiltonian, only: ComplexHk_MF, longRange, sortF, samplePoints1, plotBands, plotBandsAll, &
        outputDir => dir
    use lapack_routines, only: diagonalize
    use PostProcessingInput, only: BuildGeometry, InputParameters, ReadFock, ValidateResponseModel
    implicit none

    ! Observable grids only; model, interactions and selected bands come from Setup.
    integer(dp), parameter :: numkPlot = 12, numkPlotRec = 100, nspacing = 100000
    integer(dp), parameter :: ncount = numkPlot*numkPlot
    integer(dp) :: ind, n, m, ix, iy, ik, nomega, imin, imax, nshiftCB, nshiftVB, nbounds(4)
    integer :: unit
    character(210) :: parameters
    character(:), allocatable :: suffix
    real(dp) :: t1(2), t2(2), t3(2), g1(2), g12(2), rotation(2,2), tn(2,3), vk(2), fmu, cellArea
    real(dp), allocatable :: coord(:,:), potential(:), density(:), interaction(:,:), bands(:,:)
    real(dp), allocatable :: energyGrid(:,:,:), vEnergy(:,:,:), vMomentum(:,:,:), vPsi(:,:,:,:)
    real(dp), allocatable :: localCB(:,:), localVB(:,:), dosCB(:), dosVB(:)
    real(dp), allocatable :: zenergy(:), pivotEnergy(:,:), pivotPsi(:,:), temp(:)
    real(dp) :: pivotMomentum(2,4)
    integer(dp), allocatable :: ind_j(:), ind_l(:), points(:,:), back(:,:)
    complex(dp), allocatable :: zFock(:,:,:,:), zH(:,:)

    call ValidateResponseModel()
    if(numb < 1 .or. nLower < 1 .or. nUpper > ndim .or. nUpper-nLower+1 /= numb) &
        error stop 'Invalid selected band window in Setup'
    if(mod(numkPlot,6_dp) /= 0 .or. numkPlot < 6 .or. numkPlotRec < 1 .or. nspacing < 1) &
        error stop 'Local DOS requires numkPlot a positive multiple of 6 and positive output grids'

    allocate(coord(ndim,3), potential(ndim), density(ndim), interaction(ndim,ndim))
    call NumberNeighborCells(numI, ind)
    allocate(zFock(ndim,ndim,ind,numS), ind_j(ind), ind_l(ind))
    call OrderNeighborCells(numI, ind, ind_j, ind_l)
    call BuildGeometry(coord, t1, t2, t3, g1, g12, rotation)
    tn = reshape([t1,t2,t3],[2,3])
    cellArea = abs(t1(1)*t2(2)-t1(2)*t2(1))
    call InputParameters(parameters)
    call ReadFock(zFock, ind, parameters)
    call ReadChemicalPotential(fmu)

    ! The saved Fock matrix is the density matrix; subtract the neutral background.
    do n=1,ndim
        density(n) = real(zFock(n,n,1,1),dp)-0.5_dp
    enddo
    call longRange(interaction, coord, ndim, t1, t2)
    potential = 2.0_dp*alphaH*matmul(interaction,density)+U*density
    deallocate(interaction, density)

    ! As in the original filling-dependent split, VB contains the nominally
    ! filled bands. A partly occupied band belongs to CB; this is a band label,
    ! not an energy cutoff or a Fermi occupation factor.
    nbounds = [1_dp, max(0_dp,min(numb,floor((real(ndim,dp)+nfilling)/2,kind=dp)-nLower+1)), 0_dp, numb]
    nbounds(3) = nbounds(2)+1
    suffix = '-localDos-k'//istr(numkPlot)//'-s'//istr(nspacing)//trim(parameters)
    call OpenOutput(unit, 'LocalDosInfo')
    write(unit,'(A)') '# Local DOS: one stored spin sector; no extra spin/valley factor.'
    write(unit,'(A)') '# Energies relative to saved Mu. DOS: states/meV/cell; LDOS: states/meV/site.'
    write(unit,'(A)') '# Integrated maps: states/site. Positions in a=0.246 nm; Fourier q in 1/a.'
    write(unit,*) 'nLower nUpper nVB nCB:', nLower, nUpper, nbounds(2), numb-nbounds(2)
    write(unit,*) 'numkPlot numkPlotRec bins_per_eV:', numkPlot, numkPlotRec, nspacing
    write(unit,*) 'Mu_eV:', fmu
    close(unit)

    ! Include the closed-cell boundary when finding energy limits. The Fock
    ! input grid numk must not limit this independent postprocessing grid.
    allocate(energyGrid(numb,numkPlot+1,numkPlot+1), bands(numb,ncount))
    allocate(points(ncount,8),back(numkPlot,numkPlot))
    call samplePoints1(points,back,numkPlot,ncount)
    !$omp parallel do collapse(2) private(ix,iy,vk,zH,zenergy)
    do ix=0,numkPlot
        do iy=0,numkPlot
            allocate(zH(ndim,ndim),zenergy(ndim))
            vk = (ix*g1+iy*g12)/real(numkPlot,dp)
            call ComplexHk_MF(zH,coord,potential,alpha,zFock(:,:,:,1),ind_j,ind_l,ndim,ind,vk,tn)
            call diagonalize(zH,zenergy,'N',nLower,nUpper)
            energyGrid(:,ix+1,iy+1) = zenergy(1:numb)
            deallocate(zH,zenergy)
        enddo
    enddo
    !$omp end parallel do
    if(.not.all(ieee_is_finite(energyGrid))) error stop 'Non-finite eigenvalues'
    do ik=1,ncount
        bands(:,ik) = energyGrid(:,points(ik,1)+1,points(ik,2)+1)
    enddo
    ! Two guard bins enclose extrema even when all bands lie on one side of Mu.
    imin = floor((minval(energyGrid)-fmu)*nspacing,kind=dp)-2
    imax = ceiling((maxval(energyGrid)-fmu)*nspacing,kind=dp)+2
    nomega = imax-imin+1
    nshiftCB = 1-imin
    nshiftVB = 1+imax
    deallocate(energyGrid)
    open(unit=99,file=trim(outputDir)//'Bands'//suffix,status='replace')
    call plotBands(bands,fmu,back,numb,numkPlot,ncount,g1,g12)
    open(unit=99,file=trim(outputDir)//'BandsAll'//suffix,status='replace')
    call plotBandsAll(bands,fmu,back,numb,numkPlot,ncount,g1,g12)
    close(99)
    deallocate(bands,points,back)

    write(*,*) 'Local DOS: selected bands, bins:', numb, nomega
    write(*,*) 'LDOS arrays (GiB):', 2.0_dp*ndim*nomega*storage_size(1.0_dp)/8/1024**3
    allocate(localCB(nomega,ndim),localVB(nomega,ndim),dosCB(nomega),dosVB(nomega))
    allocate(vEnergy(numb,numkPlot+1,2),vMomentum(2,numkPlot+1,2),vPsi(numb,ndim,numkPlot+1,2))
    localCB = 0
    localVB = 0
    call FillStripe(0_dp,1)
    do ix=1,numkPlot
        call FillStripe(ix,2)
        ! Sites are independent: avoid a full LDOS-array reduction per thread.
        !$omp parallel do private(n,iy,pivotEnergy,pivotPsi,pivotMomentum,temp)
        do n=1,ndim
            allocate(pivotEnergy(numb,4),pivotPsi(numb,4),temp(nomega))
            do iy=1,numkPlot
                pivotEnergy = reshape([vEnergy(:,iy,1),vEnergy(:,iy,2), &
                    vEnergy(:,iy+1,1),vEnergy(:,iy+1,2)],[numb,4_dp])
                pivotMomentum = reshape([vMomentum(:,iy,1),vMomentum(:,iy,2), &
                    vMomentum(:,iy+1,1),vMomentum(:,iy+1,2)],[2,4])
                pivotPsi = reshape([vPsi(:,n,iy,1),vPsi(:,n,iy,2), &
                    vPsi(:,n,iy+1,1),vPsi(:,n,iy+1,2)],[numb,4_dp])
                call LocalDOS_CB(ndim,temp,pivotMomentum,pivotEnergy,pivotPsi,nomega,nspacing,nshiftCB,fmu,nbounds)
                localCB(:,n) = localCB(:,n)+temp
                call LocalDOS_VB(ndim,temp,pivotMomentum,pivotEnergy,pivotPsi,nomega,nspacing,nshiftVB,fmu,nbounds)
                localVB(:,n) = localVB(:,n)+temp
            enddo
            deallocate(pivotEnergy,pivotPsi,temp)
        enddo
        !$omp end parallel do
        vEnergy(:,:,1) = vEnergy(:,:,2)
        vMomentum(:,:,1) = vMomentum(:,:,2)
        vPsi(:,:,:,1) = vPsi(:,:,:,2)
    enddo
    deallocate(vEnergy,vMomentum,vPsi,zFock)
    ! Integral d^2k/(2*pi)^2 times cell area, then convert eV^-1 to meV^-1.
    localCB = localCB*cellArea/(4*pi*pi*1000)
    localVB = localVB*cellArea/(4*pi*pi*1000)
    if(.not.all(ieee_is_finite(localCB)) .or. .not.all(ieee_is_finite(localVB))) &
        error stop 'Non-finite local DOS'
    dosCB = sum(localCB,dim=2)
    dosVB = sum(localVB,dim=2)
    call WriteSector('CB',localCB,dosCB,nshiftCB,1)
    call WriteSector('VB',localVB,dosVB,nshiftVB,-1)
    call OpenOutput(unit,'Dos')
    do m=1,nomega
        write(unit,'(2ES24.15)') real(m-nshiftCB,dp)*1000/nspacing, dosCB(m)+dosVB(nomega+1-m)
    enddo
    close(unit)
    write(*,*) 'Local DOS written to ', trim(outputDir)

contains

    subroutine ReadChemicalPotential(mu)
        real(dp), intent(out) :: mu
        real(dp) :: iteration, value
        integer :: input, stat, rows
        character(512) :: message, line
        character(:), allocatable :: filename
        filename = trim(inputDir)//'Mu'//trim(parameters)
        open(newunit=input,file=filename,status='old',action='read',iostat=stat,iomsg=message)
        if(stat /= 0) then
            write(error_unit,'(A)') 'Cannot read '//filename//': '//trim(message)
            error stop 'Local DOS requires the saved Mu history for the same state'
        endif
        rows = 0
        do
            read(input,'(A)',iostat=stat) line
            if(stat == iostat_end) exit
            if(stat /= 0) error stop 'Cannot read Mu history'
            if(len_trim(line) == 0) cycle
            read(line,*,iostat=stat) iteration, value
            if(stat /= 0) error stop 'Malformed Mu history: expected iteration and Mu in eV'
            if(.not.ieee_is_finite(value)) error stop 'Non-finite saved Mu'
            mu = value
            rows = rows+1
        enddo
        close(input)
        if(rows == 0) error stop 'Empty Mu history'
    end subroutine ReadChemicalPotential

    subroutine FillStripe(kx,slot)
        integer(dp), intent(in) :: kx
        integer, intent(in) :: slot
        integer(dp) :: ky, band
        real(dp) :: momentum(2)
        real(dp), allocatable :: energies(:)
        complex(dp), allocatable :: matrix(:,:)
        !$omp parallel do private(ky,band,momentum,matrix,energies)
        do ky=0,numkPlot
            allocate(matrix(ndim,ndim),energies(ndim))
            momentum = (kx*g1+ky*g12)/real(numkPlot,dp)
            call ComplexHk_MF(matrix,coord,potential,alpha,zFock(:,:,:,1),ind_j,ind_l,ndim,ind,momentum,tn)
            call diagonalize(matrix,energies,'V',nLower,nUpper)
            vEnergy(:,ky+1,slot) = energies(1:numb)
            vMomentum(:,ky+1,slot) = momentum
            do band=1,numb
                vPsi(band,:,ky+1,slot) = abs(matrix(:,band))**2
            enddo
            deallocate(matrix,energies)
        enddo
        !$omp end parallel do
    end subroutine FillStripe

    subroutine OpenOutput(output,label)
        integer, intent(out) :: output
        character(*), intent(in) :: label
        integer :: stat
        character(512) :: message
        open(newunit=output,file=trim(outputDir)//label//suffix,status='replace',iostat=stat,iomsg=message)
        if(stat /= 0) then
            write(error_unit,'(A)') 'Cannot write '//trim(outputDir)//label//suffix//': '//trim(message)
            error stop 1
        endif
    end subroutine OpenOutput

    function istr(value) result(text)
        integer(dp), intent(in) :: value
        character(:), allocatable :: text
        character(32) :: buffer
        write(buffer,'(I0)') value
        text = trim(buffer)
    end function istr

    subroutine WriteSector(label,ldos,dos,shift,direction)
        character(*), intent(in) :: label
        real(dp), intent(in) :: ldos(nomega,ndim), dos(nomega)
        integer(dp), intent(in) :: shift
        integer, intent(in) :: direction
        integer :: output, peak(1)
        integer(dp) :: k
        real(dp) :: integrated(ndim)
        peak = maxloc(dos)
        ! Bin width is 1000/nspacing meV, since ldos is already per meV.
        integrated = sum(ldos,dim=1)*1000/real(nspacing,dp)
        call OpenOutput(output,'DosCheck'//label)
        do k=1,nomega
            write(output,'(2ES24.15)') direction*real(k-shift,dp)*1000/nspacing, dos(k)
        enddo
        close(output)
        call WriteMaps(label//'max',ldos(peak(1),:))
        call WriteMaps(label//'sum',integrated)
        open(newunit=output,file=trim(outputDir)//'LocalDosInfo'//suffix,status='old',position='append')
        write(output,*) label//' peak_meV:', direction*real(peak(1)-shift,dp)*1000/nspacing
        write(output,*) label//' integrated_states_per_spin:', sum(integrated)
        close(output)
    end subroutine WriteSector

    subroutine WriteMaps(label,values)
        character(*), intent(in) :: label
        real(dp), intent(in) :: values(ndim)
        integer :: output, reUnit, imUnit
        integer(dp) :: site, layer, kx, ky, first, last
        real(dp) :: q(2), qx(2), qy(2)
        complex(dp) :: fourier
        call OpenOutput(output,'LDOS-'//label//'-')
        do site=1,ndim
            write(output,'(3ES24.15)') coord(site,1:2), values(site)
        enddo
        close(output)
        qx = [2*pi*1.1_dp,0.0_dp]
        qy = [0.0_dp,4*pi/sqrt(3.0_dp)*1.1_dp]
        do layer=1,2
            first = (layer-1)*ndim/2+1
            last = layer*ndim/2
            call OpenOutput(reUnit,'DosBZ-'//label//'-Re'//istr(layer))
            call OpenOutput(imUnit,'DosBZ-'//label//'-Im'//istr(layer))
            do kx=-numkPlotRec,numkPlotRec
                do ky=-numkPlotRec,numkPlotRec
                    q = (kx*qx+ky*qy)/real(numkPlotRec,dp)
                    fourier = cmplx(0.0_dp,0.0_dp,dp)
                    do site=first,last
                        fourier = fourier+values(site)*exp(cmplx(0.0_dp,dot_product(q,coord(site,1:2)),dp))
                    enddo
                    write(reUnit,'(3ES24.15)') q, real(fourier,dp)
                    write(imUnit,'(3ES24.15)') q, aimag(fourier)
                enddo
            enddo
            close(reUnit)
            close(imUnit)
        enddo
    end subroutine WriteMaps

! The two local-DOS kernels from Spectral.f90 follow here. They use the
! existing Hamiltonian.sortF; the other Spectral procedures are not called.

! BEGIN imported local-DOS kernels
subroutine LocalDOS_CB(ndim,ldos,vk,energia,vPsi,nomega,nspacing,nshift,efermi,nbounds)
 integer(dp), intent(in) :: ndim,nomega,nspacing,nshift,nbounds(:)
 real(dp)   , intent(in) :: efermi
 real(dp)   , intent(in) :: vk(:,:),energia(:,:),vPsi(:,:)
 real(dp)   , intent(out) :: ldos(:)
 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE,nflag
 real(dp) :: energy(4)
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2)
 real(dp) :: f1(1),f2(1),f3(1)
 real(dp) :: add,weight
   ldos=0._dp
   do i_c = nbounds(3),nbounds(4)
   energy(:)=energia(i_c,:)-efermi
   energy1=energy(1)
   energy2=energy(2)
   energy3=energy(3)
   vk1(:)=vk(:,1)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)
   f1(1)=vPsi(i_c,1)
   f2(1)=vPsi(i_c,2)
   f3(1)=vPsi(i_c,3)
   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift
    ! An unresolved/flat triangle carries its integrated weight in one bin.
    weight=abs((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))
    if(n1 == n3) then
        if(n1 >= 1 .and. n1 <= nomega) ldos(n1)=ldos(n1)+weight*nspacing*(f1(1)+f2(1)+f3(1))/6.0_dp
    else
        weight=weight/abs(energy1-energy3)
    endif
    if(n1.LE.nomega) then
        if(n3.GE.1) then
    if(n3.LT.n2)then
    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         ldos(iE)=ldos(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
    end do
    endif
    if(n2.LT.n1)then
    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        ldos(iE)=ldos(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
    end do
    endif
    endif ! n3.GE.1
    endif ! n1.LE.nomega
   energy1=energy(4)
   energy2=energy(2)
   energy3=energy(3)
   vk1(:)=vk(:,4)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)
   f1(1)=vPsi(i_c,4)
   f2(1)=vPsi(i_c,2)
   f3(1)=vPsi(i_c,3)
   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift
    ! An unresolved/flat triangle carries its integrated weight in one bin.
    weight=abs((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))
    if(n1 == n3) then
        if(n1 >= 1 .and. n1 <= nomega) ldos(n1)=ldos(n1)+weight*nspacing*(f1(1)+f2(1)+f3(1))/6.0_dp
    else
        weight=weight/abs(energy1-energy3)
    endif
    if(n1.LE.nomega) then
        if(n3.GE.1) then
    if(n3.LT.n2)then
    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         ldos(iE)=ldos(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
    end do
    endif
    if(n2.LT.n1)then
    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        ldos(iE)=ldos(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
    end do
    endif
    endif ! n3.GE.1
    endif ! n1.LE.nomega
    enddo
end subroutine LocalDOS_CB

subroutine LocalDOS_VB(ndim,ldos,vk,energia,vPsi,nomega,nspacing,nshift,efermi,nbounds)
 integer(dp), intent(in) :: ndim,nomega,nspacing,nshift,nbounds(:)
 real(dp)   , intent(in) :: efermi
 real(dp)   , intent(in) :: vk(:,:),energia(:,:),vPsi(:,:)
 real(dp)   , intent(out) :: ldos(:)
 integer(dp) :: i_v, i_c
 integer(dp) :: n1,n2,n3,iE,nflag
 real(dp) :: energy(4)
 real(dp) :: energy1,energy2,energy3,vk1(2),vk2(2),vk3(2)
 real(dp) :: f1(1),f2(1),f3(1)
 real(dp) :: add,weight
   ldos=0._dp
   do i_v = nbounds(1),nbounds(2)
   energy(:)=efermi-energia(i_v,:)
   energy1=energy(1)
   energy2=energy(2)
   energy3=energy(3)
   vk1(:)=vk(:,1)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)
   f1(1)=vPsi(i_v,1)
   f2(1)=vPsi(i_v,2)
   f3(1)=vPsi(i_v,3)
   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift
    ! An unresolved/flat triangle carries its integrated weight in one bin.
    weight=abs((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))
    if(n1 == n3) then
        if(n1 >= 1 .and. n1 <= nomega) ldos(n1)=ldos(n1)+weight*nspacing*(f1(1)+f2(1)+f3(1))/6.0_dp
    else
        weight=weight/abs(energy1-energy3)
    endif
    if(n1.LE.nomega) then
        if(n3.GE.1) then
    if(n3.LT.n2)then
    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         ldos(iE)=ldos(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
    end do
    endif
    if(n2.LT.n1)then
    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        ldos(iE)=ldos(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
    end do
    endif
    endif ! n3.GE.1
    endif ! n1.LE.nomega
   energy1=energy(4)
   energy2=energy(2)
   energy3=energy(3)
   vk1(:)=vk(:,4)
   vk2(:)=vk(:,2)
   vk3(:)=vk(:,3)
   f1(1)=vPsi(i_v,4)
   f2(1)=vPsi(i_v,2)
   f3(1)=vPsi(i_v,3)
   call sortF(energy1,energy2,energy3,vk1,vk2,vk3,f1,f2,f3)

    n1=nint(energy1*nspacing)+nshift
    n2=nint(energy2*nspacing)+nshift
    n3=nint(energy3*nspacing)+nshift
    ! An unresolved/flat triangle carries its integrated weight in one bin.
    weight=abs((vk2(1)-vk3(1))*(vk1(2)-vk3(2))-(vk1(1)-vk3(1))*(vk2(2)-vk3(2)))
    if(n1 == n3) then
        if(n1 >= 1 .and. n1 <= nomega) ldos(n1)=ldos(n1)+weight*nspacing*(f1(1)+f2(1)+f3(1))/6.0_dp
    else
        weight=weight/abs(energy1-energy3)
    endif
    if(n1.LE.nomega) then
        if(n3.GE.1) then
    if(n3.LT.n2)then
    do iE=n3+1,n2
         add=weight*real(iE-n3)/real(n2-n3)
         ldos(iE)=ldos(iE)+add*(f3(1)*(n2-iE)/real(n2-n3)+f2(1)*(iE-n3)/real(n2-n3))
    end do
    endif
    if(n2.LT.n1)then
    do iE=n2+1,n1-1
        add=weight*real(n1-iE)/real(n1-n2)
        ldos(iE)=ldos(iE)+add*(f2(1)*(n1-iE)/real(n1-n2)+f1(1)*(iE-n2)/real(n1-n2))
    end do
    endif
    endif ! n3.GE.1
    endif ! n1.LE.nomega
    end do
end subroutine LocalDOS_VB
! END imported local-DOS kernels

end program compute_localDos
