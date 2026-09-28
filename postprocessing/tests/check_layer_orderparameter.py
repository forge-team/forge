#!/usr/bin/env python3
"""Check the layer projection against normalized, periodic Bloch spinors.

The two layer valleys fold to opposite moire momenta.  Each fixture therefore
uses two genuine Bloch sectors, not one nonperiodic superposition of all eight
continuum components.  Both independent complex interlayer directions are
present, and their known envelope products are the independent reference.
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
BUILD = POST / '.build/layer-orderparameter-tests'
BUILD.mkdir(parents=True, exist_ok=True)
FLAGS = ['-O0', '-g', '-fopenmp', '-ffree-line-length-none',
         '-fallow-argument-mismatch', '-fcheck=all']
LIBS = ['-framework', 'Accelerate'] if platform.system() == 'Darwin' else ['-llapack', '-lblas']
ENV = dict(os.environ, OMP_NUM_THREADS='2')
RESULTS = []
TOL = 3e-12


def run(cmd, folder, name, succeeds=True):
    result = subprocess.run(cmd, cwd=folder, env=ENV, text=True, capture_output=True)
    (folder / (name + '.log')).write_text(result.stdout + result.stderr)
    if succeeds and result.returncode:
        raise RuntimeError(f'{folder.name}/{name}: {result.stdout[-2000:]}\n{result.stderr[-2500:]}')
    if not succeeds and result.returncode == 0:
        raise AssertionError(f'{folder.name}/{name}: expected rejection')
    return result


def close(actual, expected, context):
    assert math.isfinite(actual.real) and math.isfinite(actual.imag), context
    assert abs(actual - expected) < TOL, (context, actual, expected)


def complex_columns(words):
    values = [float(x) for x in words]
    assert len(values) % 2 == 0
    return [complex(values[i], values[i+1]) for i in range(0, len(values), 2)]


# Keep the fixture's byte layout exactly that of the actual producer, including
# its central lower triangle and subsequent complete matrices in m/i/j order.
main = (REPO / 'Main.f90').read_text()
writer = re.search(r'rcc=0\b.*?close\(91\)', main[main.index('if (nWriteFock.EQ.1)'):], re.S).group()
FIXTURE = r'''program fixture
use Setup
use Geometry, only: NumberNeighborCells, OrderNeighborCells
use PostProcessingInput, only: BuildGeometry
implicit none
integer(dp) :: i,j,m,s,layer,sub,nspin,rcc,numNeighborCells,my_iostat,ncells
integer(dp),allocatable :: nc1(:),nc2(:)
real(dp) :: coords(ndim,3),t1(2),t2(2),t3(2),g1(2),g12(2),rot(2,2)
real(dp) :: kval(2,2),sector_k(2,2),cell(2),weights(2),norm,sgn
complex(dp) :: amps(2,2,2),psi(ndim,2),expected(4),gauge
complex(dp),allocatable :: zFock(:,:,:,:)
character(1024) :: filename
character(150) :: my_iomsg
character(32) :: mode,spintext
call get_command_argument(1,mode)
call InitParameters()
call BuildGeometry(coords,t1,t2,t3,g1,g12,rot)
call NumberNeighborCells(numI,numNeighborCells)
allocate(nc1(numNeighborCells),nc2(numNeighborCells))
call OrderNeighborCells(numI,numNeighborCells,nc1,nc2)
allocate(zFock(ndim,ndim,numNeighborCells,numS))
ncells=ndim/4
kval(:,1)=matmul(rot,[-4.0_dp*pi/3.0_dp,0.0_dp])
kval(:,2)=matmul(transpose(rot),[-4.0_dp*pi/3.0_dp,0.0_dp])
sector_k(:,1)=kval(:,1)
sector_k(:,2)=-kval(:,1)
! Confirm that each state's layer/valley components share a Bloch boundary
! condition on BOTH moire generators, including their complex phases.
if(abs(exp(cmplx(0.0_dp,dot_product(kval(:,1)+kval(:,2),t1),dp))-1)>1e-12_dp) &
    error stop 'incompatible t1 Bloch phases'
if(abs(exp(cmplx(0.0_dp,dot_product(kval(:,1)+kval(:,2),t2),dp))-1)>1e-12_dp) &
    error stop 'incompatible t2 Bloch phases'
zFock=0
do nspin=1,numS
    weights=[.37_dp+.03_dp*nspin,.61_dp-.04_dp*nspin]
    ! Unequal arbitrary complex amplitudes distinguish sublattices, opposite
    ! layer directions, conjugations, and spin sectors.
    do s=1,2
        do layer=1,2
            do sub=1,2
                amps(s,layer,sub)=cmplx(.17_dp+.11_dp*s+.07_dp*layer*sub+.03_dp*nspin*sub, &
                    -.29_dp+.13_dp*s*sub-.09_dp*layer+.04_dp*nspin*layer*sub,dp)
            enddo
        enddo
        norm=sqrt(real(ncells,dp)*sum(abs(amps(s,:,:))**2))
        amps(s,:,:)=amps(s,:,:)/norm
    enddo
    gauge=1
    if(trim(mode)=='gauge') gauge=exp(cmplx(0.0_dp,.731_dp,dp))
    do i=1,ndim
        layer=(i-1)/(2*ncells)+1
        sub=modulo((i-1)/ncells,2_dp)+1
        do s=1,2
            sgn=real((3-2*s)*(3-2*layer),dp)
            psi(i,s)=gauge*amps(s,layer,sub)* &
                exp(cmplx(0.0_dp,sgn*dot_product(kval(:,layer),coords(i,1:2)),dp))
        enddo
    enddo
    do s=1,2
        if(abs(sum(abs(psi(:,s))**2)-1)>1e-13_dp) error stop 'unnormalized fixture'
    enddo
    if(trim(mode)/='zero') then
        do m=1,numNeighborCells
            cell=real(nc1(m),dp)*t1+real(nc2(m),dp)*t2
            do i=1,ndim
                do j=1,ndim
                    do s=1,2
                        zFock(i,j,m,nspin)=zFock(i,j,m,nspin)+weights(s)*conjg(psi(i,s))*psi(j,s)* &
                            exp(cmplx(0.0_dp,-dot_product(sector_k(:,s),cell),dp))
                    enddo
                enddo
            enddo
        enddo
    endif
    ! Direct continuum density products; no loop or inversion coefficients
    ! from the implementation are used to obtain this reference.
    expected(1)=weights(2)*conjg(amps(2,2,1))*amps(2,1,2)
    expected(2)=weights(2)*conjg(amps(2,2,2))*amps(2,1,1)
    expected(3)=weights(1)*conjg(amps(1,1,1))*amps(1,2,2)
    expected(4)=weights(1)*conjg(amps(1,1,2))*amps(1,2,1)
    if(trim(mode)=='zero') expected=0
    write(spintext,'(I0)') nspin
    open(17,file='expected-spin'//trim(spintext)//'.dat',status='replace')
    do i=1,4
        write(17,'(4(ES25.16,1X))') real(expected(i)),aimag(expected(i)), &
            real(ncells*expected(i)),aimag(ncells*expected(i))
    enddo
    close(17)
    write(filename,'(A,A,I0,A)') trim(dirFock),'Fock-nspin',nspin,trim(parameters)
    open(91,file=filename,status='replace',form='unformatted',access='direct',recl=2*dp)
''' + writer + r'''
enddo
end program fixture
'''


def build_case(ntheta, shells, spins):
    folder = BUILD / f'i{ntheta}_I{shells}_s{spins}'
    folder.mkdir(exist_ok=True)
    (folder / 'dataFock').mkdir(exist_ok=True)
    (folder / 'output').mkdir(exist_ok=True)
    setup = (REPO / 'Setup.f90').read_text()
    for name, value in dict(ntheta=ntheta, numI=shells, numS=spins,
                            numC=1, numk=3, numb=4, nrelax=0).items():
        setup, count = re.subn(r'(^\s*integer(?:\(dp\))?,\s*parameter\s*::\s*' + name + r'\s*=)[^!\n]*',
                              rf'\g<1> {value} ', setup, flags=re.M | re.I)
        assert count == 1, name
    (folder / 'Setup.f90').write_text(setup)
    sources = [REPO / 'LapackRoutines.f90', folder / 'Setup.f90', REPO / 'TightBinding.f90',
               REPO / 'Geometry.f90', POST / 'PostProcessingInput.f90']
    objects = [p.stem + '.o' for p in sources]
    run(['gfortran', *FLAGS, '-c', *map(str, sources)], folder, 'modules')
    run(['gfortran', *FLAGS, '-o', 'driver', str(POST / 'compute_LayerOrderParameter.f90'),
         *objects, *LIBS], folder, 'driver_compile')
    return folder, objects


def check_outputs(folder, ntheta, spins, mode):
    ncells = 3*ntheta**2 + 3*ntheta + 1
    for spin in range(1, spins+1):
        reference = [complex_columns(line.split()) for line in
                     (folder / f'expected-spin{spin}.dat').read_text().splitlines()]
        expected_local = [row[0] for row in reference]
        expected_rho = [row[1] for row in reference]
        local_files = list((folder / 'output').glob(f'LayerOrderParameterLocal-numb*-nspin{spin}-*'))
        summary_files = list((folder / 'output').glob(f'LayerOrderParameter-numb*-nspin{spin}-*'))
        assert len(local_files) == len(summary_files) == 1, (folder, local_files, summary_files)
        counts = dict(A=0, B=0)
        grid_sums = dict(A=[0j]*4, B=[0j]*4)
        layer1_indices, layer2_indices = set(), set()
        for line in local_files[0].read_text().splitlines():
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            words = line.split()
            assert len(words) == 13, (local_files[0], words)
            i, j, sub = int(words[0]), int(words[1]), words[2]
            assert sub in ('A', 'B')
            assert 1 <= i <= 2*ncells and 2*ncells < j <= 4*ncells
            assert i not in layer1_indices and j not in layer2_indices
            layer1_indices.add(i)
            layer2_indices.add(j)
            assert all(math.isfinite(float(x)) for x in words[3:5])
            values = complex_columns(words[5:])
            counts[sub] += 1
            for c, (actual, exact) in enumerate(zip(values, expected_local)):
                close(actual, exact, f'{folder.name}/{mode}/spin{spin}/{sub}{i}/channel{c}')
                grid_sums[sub][c] += actual
        assert counts == dict(A=ncells, B=ncells), counts
        for sub in ('A', 'B'):
            for c in range(4):
                close(grid_sums[sub][c], expected_rho[c], f'{folder.name}/{mode}/{sub}-grid/channel{c}')
        summary = [line.split() for line in summary_files[0].read_text().splitlines()
                   if line.strip() and not line.lstrip().startswith('#')]
        assert len(summary) == 4, summary
        assert [row[0] for row in summary] == ['AK2_BKp1', 'BK2_AKp1', 'AK1_BKp2', 'BK1_AKp2']
        for c, row in enumerate(summary):
            assert len(row) == 5, row
            total, per_cell = complex_columns(row[1:])
            close(total, expected_rho[c], f'{folder.name}/{mode}/rho/channel{c}')
            close(per_cell, expected_local[c], f'{folder.name}/{mode}/rho-per-cell/channel{c}')
            close(total, .5*(grid_sums['A'][c]+grid_sums['B'][c]), 'average center grids')
            if mode != 'zero':
                # Guard against vacuous zero tests and doubled/tripled sums.
                assert abs(expected_rho[c]) > 1e-4 and abs(expected_rho[c].imag) > 1e-5
                assert abs(total-2*expected_rho[c]) > 1e-4
                assert abs(total-3*expected_rho[c]) > 1e-4
        RESULTS.append(f'{folder.name}, {mode}, spin={spin}: complex local channels, both layer directions, '
                       'A/B grids and normalized sums PASS')


for ntheta, shells, spins in ((1, 2, 2), (2, 1, 1)):
    folder, objects = build_case(ntheta, shells, spins)
    geometry_folder = folder / 'geometry-without-input'
    geometry_folder.mkdir(exist_ok=True)
    assert not (geometry_folder / 'dataFock').exists()
    run([str(folder / 'driver'), '--check-geometry'], geometry_folder, 'geometry')
    assert sorted(p.name for p in geometry_folder.iterdir()) == ['geometry.log'], \
        '--check-geometry created output files'
    RESULTS.append(f'{folder.name}: geometry check without Fock input or output changes PASS')
    (folder / 'fixture.f90').write_text(FIXTURE)
    run(['gfortran', *FLAGS, '-o', 'fixture', 'fixture.f90', *objects, *LIBS], folder, 'fixture_compile')
    for mode in ('complex', 'gauge', 'zero'):
        run([str(folder / 'fixture'), mode], folder, f'fixture_{mode}')
        run([str(folder / 'driver')], folder, f'driver_{mode}')
        check_outputs(folder, ntheta, spins, mode)

folder, _ = build_case(1, 0, 1)
rejected = run([str(folder / 'driver'), '--check-geometry'], folder, 'missing_translation', succeeds=False)
assert 'numI' in rejected.stdout + rejected.stderr, rejected.stdout + rejected.stderr
RESULTS.append('i1_I0_s1: insufficient saved Fock translations rejected PASS')

(BUILD / 'results.json').write_text(json.dumps(RESULTS, indent=2)+'\n')
print('\n'.join(RESULTS))
