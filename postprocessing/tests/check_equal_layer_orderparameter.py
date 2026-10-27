#!/usr/bin/env python3
"""Compare the new equal-layer kernel with the actual Main_OrderParamsFull.

The production contains procedures are extracted unchanged into a test program.
The test calls the actual OrderPlot routine, saves the fixture with Main's packed
writer, and then runs the actual Main_OrderParamsFull on those files.  A general
periodic positive density tests the production three-neighbor stencil against
the legacy result at every A center in both layers; a coherent Dirac fixture
also tests their common smooth-envelope limit away from boundaries.
"""
from pathlib import Path
import json
import math
import os
import platform
import re
import subprocess

POST = Path(__file__).resolve().parents[1]
REPO = POST.parent
BUILD = POST / '.build/equal-layer-tests'
BUILD.mkdir(parents=True, exist_ok=True)
FLAGS = ['-O0', '-g', '-fopenmp', '-ffree-line-length-none',
         '-fallow-argument-mismatch', '-fcheck=all']
LIBS = ['-framework', 'Accelerate'] if platform.system() == 'Darwin' else ['-llapack', '-lblas']
ENV = dict(os.environ, OMP_NUM_THREADS='2')
TOL = 5e-12


def run(cmd, folder, name):
    result = subprocess.run(cmd, cwd=folder, env=ENV, text=True, capture_output=True)
    (folder / (name + '.log')).write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f'{folder.name}/{name}: {result.stdout[-2000:]}\n{result.stderr[-2500:]}')
    return result


source = (POST / 'compute_LayerOrderParameter.f90').read_text()
declarations = source[:source.index('    if(nlayers /= 2)')]
declarations = declarations.replace('program compute_LayerOrderParameter', 'program equal_layer_fixture')
declarations = declarations.replace('    implicit none', '''    use OrderPlot, only: getKekuleLattice, getKekuleNeighbors, InterSubInterValAlt
    use Setup, only: dirFock
    implicit none''')
geometry = source[source.index('    allocate(coord('):source.index('    ! Pair physical centers once')]
procedures = source[source.index('\ncontains\n'):source.rindex('end program compute_LayerOrderParameter')]
main = (REPO / 'Main.f90').read_text()
writer = re.search(r'rcc=0\b.*?close\(91\)', main[main.index('if (nWriteFock.EQ.1)'):], re.S).group()

EXTRA_DECLARATIONS = r'''
    integer(dp) :: j,m,q,n,rcc,my_iostat,numNeighborCells,ncells,idx,legacyIndex,neighborIndex
    integer(dp) :: KekuleLattice(ndim/2,6,2),KekuleNeighbors(ndim/2,2,3),TnTonUnitCell12(0:6,2)
    integer(dp) :: imageCell(2),site,ix,iy,inside
    real(dp) :: tn(6,2),cell(2),qvec(2,3),norm,weight(3),physicalR(2)
    complex(dp) :: psi(ndim,3),amps(2,2,2),legacyY(ndim/2),legacyX(ndim/2),expected(2,2)
    complex(dp) :: current(ndim/2),production(4,ndim),smooth(4),modes(3),phase(ndim),d0,d1,d2
    character(150) :: my_iomsg
    character(20) :: mode
'''

BODY = r'''
    call get_command_argument(1,mode)
    ncells=ndim/4
    numNeighborCells=ind
    tn(1,:)=t1; tn(2,:)=t2; tn(3,:)=t3
    tn(4,:)=-t1; tn(5,:)=-t2; tn(6,:)=-t3
    TnTonUnitCell12(0,:)=[0,0]; TnTonUnitCell12(1,:)=[1,0]
    TnTonUnitCell12(2,:)=[0,1]; TnTonUnitCell12(3,:)=[-1,1]
    TnTonUnitCell12(4,:)=[-1,0]; TnTonUnitCell12(5,:)=[0,-1]
    TnTonUnitCell12(6,:)=[1,-1]
    call getKekuleLattice(coord,ndim,ang_mat,a1,a2,tn,KekuleLattice)
    call getKekuleNeighbors(KekuleNeighbors,coord,ndim,a1,a2,ang_mat,tn)
    if(any(KekuleLattice(:,:,1)==0).or.any(KekuleNeighbors==0)) error stop 'Missing legacy stencil site'
    allocate(zFock(ndim,ndim,ind,numS))
    psi=0; amps=0; expected=0
    weight=[.31_dp,.57_dp,.22_dp]
    qvec(:,1)=.13_dp*vq1+.29_dp*vq12
    qvec(:,2)=-.21_dp*vq1+.11_dp*vq12
    qvec(:,3)=.34_dp*vq1-.17_dp*vq12
    if(trim(mode)=='periodic') then
        ! Three normalized complex orbitals with arbitrary intra-cell texture;
        ! their independent Bloch momenta give a positive Hermitian density.
        do q=1,3
            do i=1,ndim
                psi(i,q)=cmplx(sin(.37_dp*i*q)+.17_dp*cos(.21_dp*i), &
                    cos(.23_dp*i*q)-.29_dp*sin(.41_dp*i),dp)
            enddo
            psi(:,q)=psi(:,q)/sqrt(sum(abs(psi(:,q))**2))
        enddo
    else if(trim(mode)=='dirac') then
        ! Intralayer intervalley coherence is not periodic in this moire cell.
        ! This fixture is extended periodically for storage, and its analytic
        ! Dirac reference is asserted ONLY on complete interior stencils.
        qvec=0
        do layer=1,2
            do sub=1,2
                do q=1,2
                    amps(q,layer,sub)=cmplx(.19_dp+.13_dp*q+.07_dp*layer*sub, &
                        -.31_dp+.17_dp*q*sub-.09_dp*layer,dp)
                enddo
            enddo
            do sub=1,2
                call SiteRange(layer,sub,first,last)
                do i=first,last
                    psi(i,layer)=amps(1,layer,sub)*exp(zi*dot_product(valleyK(:,layer),coord_ref(i,1:2))) &
                        +amps(2,layer,sub)*exp(-zi*dot_product(valleyK(:,layer),coord_ref(i,1:2)))
                enddo
            enddo
            norm=sqrt(sum(abs(psi(:,layer))**2))
            psi(:,layer)=psi(:,layer)/norm
            amps(:,layer,:)=amps(:,layer,:)/norm
            expected(1,layer)=weight(layer)*conjg(amps(1,layer,1))*amps(2,layer,2)
            expected(2,layer)=weight(layer)*conjg(amps(1,layer,2))*amps(2,layer,1)
        enddo
    else
        error stop 'Unknown fixture mode'
    endif
    zFock=0
    do nspin=1,numS
        do m=1,ind
            cell=real(ind_j(m),dp)*t1+real(ind_l(m),dp)*t2
            do i=1,ndim
                do j=1,ndim
                    do q=1,3
                        zFock(i,j,m,nspin)=zFock(i,j,m,nspin)+ &
                            weight(q)*conjg(psi(i,q))*psi(j,q)*exp(-zi*dot_product(qvec(:,q),cell))
                    enddo
                enddo
            enddo
        enddo
    enddo
    nspin=1
    if(maxval(abs(zFock(:,:,1,1)-transpose(conjg(zFock(:,:,1,1)))))>1e-13_dp) &
        error stop 'Non-Hermitian fixture'
    call InterSubInterValAlt(legacyY,legacyX,KekuleLattice,KekuleNeighbors,ndim,ind,ind_j,ind_l, &
        TnTonUnitCell12,zFock(:,:,:,1))
    ! Set both layer arguments equal, reusing the production loop caching and
    ! inversion exactly.  The two direction slots must then give duplicate pairs.
    do layer=1,2
        do sub=1,2
            call SiteRange(layer,sub,first,last)
            sublattice='A'
            if(sub==2) sublattice='B'
            do i=first,last
                idx=i-(layer-1)*ndim/2
                r1=coord_ref(i,1:2)
                call CacheLoopSites(idx,1_dp,layer,layer,sublattice,r1,r1)
                call CacheLoopSites(idx,2_dp,layer,layer,sublattice,r1,r1)
                call CachedLoopChannels(idx,1_dp,Delta21)
                call CachedLoopChannels(idx,2_dp,Delta12)
                phase(i)=exp(2.0_dp*zi*dot_product(valleyK(:,layer),r1))
                call ProjectLayerChannels(Delta21,Delta12,sublattice,phase(i),production(:,i))
                if(sub==1) then
                    legacyIndex=i-(layer-1)*ncells
                    ! The legacy C(r) loop is the conjugate of new A loop 3.
                    current(legacyIndex)=conjg(Delta21(3))
                endif
            enddo
        enddo
    enddo
    open(17,file='comparison.dat',status='replace')
    write(17,'(A)') '# layer site sub interior Re/Im phase legacyY legacyX production1 production2 production3 production4 smooth1 smooth2 expectedY expectedX'
    do layer=1,2
        do sub=1,2
            call SiteRange(layer,sub,first,last)
            sublattice='A'
            if(sub==2) sublattice='B'
            do i=first,last
                smooth=0; localRho=0
                if(sub==1) then
                    legacyIndex=i-(layer-1)*ncells
                    d0=current(legacyIndex); d1=0; d2=0
                    do n=1,3
                        neighborIndex=KekuleNeighbors(legacyIndex,1,n)-(layer-1)*ncells
                        d1=d1+current(neighborIndex)/3.0_dp
                        neighborIndex=KekuleNeighbors(legacyIndex,2,n)-(layer-1)*ncells
                        d2=d2+current(neighborIndex)/3.0_dp
                    enddo
                    modes=conjg([d1,d2,d0])
                    call ProjectLayerChannels(modes,modes,'A',phase(i),smooth)
                    localRho(1)=legacyY(legacyIndex); localRho(2)=legacyX(legacyIndex)
                endif
                ! A conservative interior mask for BOTH estimators: every site
                ! within two primitive translations must lie in the home cell.
                inside=1
                r1=coord_ref(i,1:2)
                do ix=-2,2
                    do iy=-2,2
                        physicalR=r1+real(ix,dp)*layerVectors(:,1,layer)+real(iy,dp)*layerVectors(:,2,layer)
                        call LocateReferenceSite(physicalR,layer,sub,site,imageCell)
                        if(any(imageCell/=0)) inside=0
                        physicalR=physicalR+real(3-2*sub,dp)*deltaVec(:,2,layer)
                        call LocateReferenceSite(physicalR,layer,3-sub,site,imageCell)
                        if(any(imageCell/=0)) inside=0
                    enddo
                enddo
                write(17,'(2(I0,1X),A1,1X,I0,22(1X,ES25.16))') layer,i,sublattice,inside, &
                    real(phase(i)),aimag(phase(i)),(real(localRho(n)),aimag(localRho(n)),n=1,2), &
                    (real(production(n,i)),aimag(production(n,i)),n=1,4),(real(smooth(n)),aimag(smooth(n)),n=1,2), &
                    (real(expected(n,layer)),aimag(expected(n,layer)),n=1,2)
            enddo
        enddo
    enddo
    close(17)
    call InputParameters(parameters)
    nspin=1
    filename=trim(dirFock)//'Fock-nspin1'//trim(parameters)
    open(91,file=trim(filename),status='replace',form='unformatted',access='direct',recl=2*dp)
'''


def complex_columns(words):
    return [complex(float(words[i]), float(words[i+1])) for i in range(0, len(words), 2)]


def build_case(ntheta):
    folder = BUILD / f'i{ntheta}'
    folder.mkdir(exist_ok=True)
    setup = (REPO / 'Setup.f90').read_text()
    for name, value in dict(ntheta=ntheta, numI=2, numS=1, numC=1, numk=3, numb=4, nrelax=0).items():
        setup, count = re.subn(r'(^\s*integer(?:\(dp\))?,\s*parameter\s*::\s*' + name + r'\s*=)[^!\n]*',
                              rf'\g<1> {value} ', setup, flags=re.M | re.I)
        assert count == 1, name
    (folder / 'Setup.f90').write_text(setup)
    sources = [REPO / 'LapackRoutines.f90', folder / 'Setup.f90', REPO / 'TightBinding.f90',
               REPO / 'Geometry.f90', REPO / 'HartreeFock.f90', REPO / 'OrderParameter.f90',
               POST / 'PostProcessingInput.f90']
    objects = [p.stem + '.o' for p in sources]
    run(['gfortran', *FLAGS, '-c', *map(str, sources)], folder, 'modules')
    (folder / 'fixture.f90').write_text(declarations + EXTRA_DECLARATIONS + geometry + BODY + writer +
                                      procedures + 'end program equal_layer_fixture\n')
    run(['gfortran', *FLAGS, '-o', 'fixture', 'fixture.f90', *objects, *LIBS], folder, 'fixture_compile')
    run(['gfortran', *FLAGS, '-o', 'legacy_driver', str(POST / 'Main_OrderParamsFull.f90'),
         *objects, *LIBS], folder, 'legacy_compile')
    return folder


results = []
for ntheta in (2, 5):
    build = build_case(ntheta)
    for mode in ('periodic', 'dirac'):
        folder = build / mode
        folder.mkdir(exist_ok=True)
        (folder / 'dataFock').mkdir(exist_ok=True)
        (folder / 'output').mkdir(exist_ok=True)
        run([str(build / 'fixture'), mode], folder, 'fixture')
        run([str(build / 'legacy_driver')], folder, 'legacy_driver')
        legacy_file, = (folder / 'output').glob('InterSubInterVal-numb*')
        legacy_rows = [complex_columns(line.split()) for line in legacy_file.read_text().splitlines()]
        metric = dict(ntheta=ntheta, fixture=mode, max_smoothed_error=0.0,
                      max_production_difference=0.0, max_main_file_error=0.0, max_dirac_error=0.0,
                      max_dirac_interior_difference=0.0,
                      max_duplicate_error=0.0, interior_A=0, interior_B=0, max_legacy_magnitude=0.0)
        metric['interior_by_layer'] = {str(layer): dict(A=0, B=0) for layer in (1, 2)}
        metric['production_by_layer'] = {str(layer): dict(A=0, max_difference=0.0) for layer in (1, 2)}
        counts = dict(A=0, B=0)
        for line in (folder / 'comparison.dat').read_text().splitlines():
            if line.startswith('#'):
                continue
            words = line.split()
            layer, site, sub, interior = int(words[0]), int(words[1]), words[2], int(words[3])
            phase, ly, lx, r1, r2, r3, r4, s1, s2, ey, ex = complex_columns(words[4:])
            assert all(math.isfinite(x.real) and math.isfinite(x.imag) for x in (ly, lx, r1, r2, s1, s2))
            metric['max_duplicate_error'] = max(metric['max_duplicate_error'], abs(r1-r3), abs(r2-r4))
            if sub == 'A':
                metric['max_smoothed_error'] = max(metric['max_smoothed_error'], abs(s1-phase*ly), abs(s2-phase*lx))
                difference = max(abs(r1-phase*ly), abs(r2-phase*lx))
                metric['max_production_difference'] = max(metric['max_production_difference'], difference)
                metric['production_by_layer'][str(layer)]['A'] += 1
                metric['production_by_layer'][str(layer)]['max_difference'] = max(
                    metric['production_by_layer'][str(layer)]['max_difference'], difference)
                metric['max_legacy_magnitude'] = max(metric['max_legacy_magnitude'], abs(ly), abs(lx))
                file_y, file_x = legacy_rows[counts['A']]
                metric['max_main_file_error'] = max(metric['max_main_file_error'], abs(file_y-ly), abs(file_x-lx))
                # ES12.5 in the actual Main output rounds six significant digits.
                for file_value, value in ((file_y, ly), (file_x, lx)):
                    assert abs(file_value-value) <= 8e-6*max(abs(value), 1e-20), (file_value, value)
            if interior:
                metric['interior_' + sub] += 1
                metric['interior_by_layer'][str(layer)][sub] += 1
                if mode == 'dirac':
                    metric['max_dirac_error'] = max(metric['max_dirac_error'], abs(r1-ey), abs(r2-ex))
                    if sub == 'A':
                        metric['max_dirac_error'] = max(metric['max_dirac_error'], abs(phase*ly-ey), abs(phase*lx-ex))
                        metric['max_dirac_interior_difference'] = max(
                            metric['max_dirac_interior_difference'], abs(r1-phase*ly), abs(r2-phase*lx))
            counts[sub] += 1
        assert counts == dict(A=len(legacy_rows), B=len(legacy_rows)), counts
        assert metric['max_duplicate_error'] < TOL, metric
        assert metric['max_smoothed_error'] < TOL, metric
        assert metric['max_production_difference'] < TOL, metric
        assert all(grid['A'] == len(legacy_rows)//2 and grid['max_difference'] < TOL
                   for grid in metric['production_by_layer'].values()), metric
        if mode == 'dirac' and ntheta == 5:
            assert all(count > 0 for grid in metric['interior_by_layer'].values() for count in grid.values()), metric
            assert metric['max_dirac_error'] < TOL, metric
            assert metric['max_dirac_interior_difference'] < TOL, metric
        results.append(metric)
        print(json.dumps(metric, sort_keys=True))

(BUILD / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
print('PASS: production equal-layer averaging at all A centers in both layers, actual packed-input driver, and interior Dirac limit')
