#!/usr/bin/env python3
"""Check geometry, packed I/O, real driver prologues and Hermitian BLAS use.

All fixtures and modified small Setup copies live under ../.build/tests.
No production Setup or data is modified. Requires gfortran and BLAS/LAPACK.
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
BUILD = POST / '.build/tests'
BUILD.mkdir(parents=True, exist_ok=True)
FLAGS = ['-O0', '-g', '-fopenmp', '-ffree-line-length-none', '-fallow-argument-mismatch', '-fcheck=all']
LIBS = ['-framework', 'Accelerate'] if platform.system() == 'Darwin' else ['-llapack', '-lblas']
ENV = dict(os.environ, OMP_NUM_THREADS='1')
RESULTS = []


def run(cmd, cwd, log, success=True):
    result = subprocess.run(cmd, cwd=cwd, env=ENV, text=True, capture_output=True)
    (cwd / log).write_text(result.stdout + result.stderr)
    if success and result.returncode:
        raise RuntimeError(f'{cwd.name}/{log}: {result.stderr[-2500:]}')
    return result


def setup_source(**changes):
    text = (REPO / 'Setup.f90').read_text()
    for name, value in changes.items():
        text, count = re.subn(r'(^\s*integer(?:\(dp\))?,\s*parameter\s*::\s*' + name + r'\s*=)[^!\n]*',
                              rf'\g<1> {value} ', text, flags=re.M | re.I)
        assert count == 1, (name, count)
    return text


def compile_modules(folder, **changes):
    folder.mkdir(exist_ok=True)
    (folder / 'dataFock').mkdir(exist_ok=True)
    (folder / 'Setup.f90').write_text(setup_source(**changes))
    paths = [POST / 'LapackRoutines.f90', folder / 'Setup.f90', REPO / 'TightBinding.f90',
             REPO / 'Geometry.f90', REPO / 'HartreeFock.f90', POST / 'PostProcessingInput.f90']
    run(['gfortran', *FLAGS, '-c', *map(str, paths)], folder, 'compile_modules.log')


def compile_program(folder, name, text, extra_objects=()):
    (folder / f'{name}.f90').write_text(text)
    objects = ['LapackRoutines.o', 'Setup.o', 'TightBinding.o', 'Geometry.o', 'HartreeFock.o', 'PostProcessingInput.o']
    run(['gfortran', *FLAGS, '-o', name, f'{name}.f90', *objects, *extra_objects, *LIBS],
        folder, f'{name}_compile.log')


geometry = (POST / 'tests/reference_geometry.f90').read_text()
main = (REPO / 'Main.f90').read_text()
writer = re.search(r'rcc=0\b.*?close\(91\)', main[main.index('if (nWriteFock.EQ.1)'):], re.S).group()

for kind, spins, cells, relax in [(8, 2, 0, 0), (8, 2, 1, 1), (8, 2, 2, 2), (4, 2, 1, 0), (8, 1, 1, 0)]:
    name = f'dp{kind}_s{spins}_I{cells}_relax{relax}'
    folder = BUILD / name
    compile_modules(folder, dp=kind, ntheta=1, numS=spins, numI=cells, numC=1, nrelax=relax)
    code = f'''program check
use Setup
use Geometry, only: NumberNeighborCells, WignerSeitzCell, LatticeRelaxationKoshino, LatticeRelaxationCarr
use PostProcessingInput
implicit none
integer(dp) :: numNeighborCells, n, i, j, m, nspin, rcc, my_iostat
real(dp) :: aMoire, cs, sn, t1(2),t2(2),t3(2),g1(2),g12(2),RotMatrix(2,2)
real(dp) :: Coords(ndim,3), Actual(ndim,3), at1(2),at2(2),at3(2),ag1(2),ag12(2),arot(2,2)
complex(dp), allocatable :: zFock(:,:,:,:), expected(:,:,:,:)
character(1024) :: filename
character(150) :: my_iomsg
character(32) :: spin, mode
character(220) :: suffix
call get_command_argument(1,mode)
call NumberNeighborCells(numI,numNeighborCells)
allocate(zFock(ndim,ndim,numNeighborCells,numS),expected(ndim,ndim,numNeighborCells,numS))
call BuildGeometry(Actual,at1,at2,at3,ag1,ag12,arot)
{geometry}
if(maxval(abs(Coords-Actual)) > 100*spacing(1.0_dp)) error stop 'coordinate mismatch'
if(maxval(abs(t1-at1))+maxval(abs(t2-at2))+maxval(abs(t3-at3)) > 100*spacing(1.0_dp)) error stop 'lattice mismatch'
if(maxval(abs(g1-ag1))+maxval(abs(g12-ag12)) > 100*spacing(1.0_dp)) error stop 'reciprocal mismatch'
do nspin=1,numS
do m=1,numNeighborCells
do i=1,ndim
do j=1,ndim
expected(i,j,m,nspin)=cmplx(10000*nspin+1000*m+100*i+j,10*i+j,dp)
enddo
enddo
enddo
do i=1,ndim
expected(i,i,1,nspin)=real(expected(i,i,1,nspin),dp)
do j=1,i-1
expected(j,i,1,nspin)=conjg(expected(i,j,1,nspin))
enddo
enddo
enddo
call InputParameters(suffix)
if(trim(suffix)/=trim(parameters)) error stop 'producer filename mismatch'
if(trim(mode)/='read-only')then
zFock=expected
do nspin=1,numS
write(spin,'(I0)') nspin
filename=trim(dirFock)//'Fock-nspin'//trim(spin)//trim(suffix)
open(91,file=filename,status='replace',form='unformatted',access='direct',recl=2*dp)
{writer}
enddo
endif
call ReadFock(zFock,numNeighborCells,suffix)
if(any(zFock/=expected)) error stop 'complex matrix mismatch'
print *, 'geometry, filename and complex Fock roundtrip PASS'
end program
'''
    compile_program(folder, 'check', code)
    run([str(folder / 'check')], folder, 'check.log')
    RESULTS.append(f'{name}: geometry/reference and complex packed Fock PASS')
    fixture = next((folder / 'dataFock').glob('Fock-nspin1*'))
    with fixture.open('ab') as stream:
        stream.write(b'X')
    bad = run([str(folder / 'check'), 'read-only'], folder, 'size_rejection.log', success=False)
    assert bad.returncode != 0 and 'Fock file size' in bad.stderr
    run([str(folder / 'check')], folder, 'restore_fixture.log')
    RESULTS.append(f'{name}: wrong-size input rejected PASS')

# Use actual adapted driver prologues, including their own declarations,
# allocations, geometry call, filename selection and ReadFock call.
folder = BUILD / 'dp8_s1_I1_relax0'
run(['gfortran', *FLAGS, '-c', str(POST / 'Hamiltonian.f90')], folder, 'compile_legacy.log')
for source in POST.glob('Compute_*.f90'):
    if 'KramersKronig' in source.name:
        continue  # This program reads conductivity, not geometry or Fock.
    text = source.read_text()
    prefix = text[:text.index('call ReadFock(zFock, ind, parameters)') + len('call ReadFock(zFock, ind, parameters)')]
    prefix += '''
if(zFock(2,1,1,1)/=cmplx(11201.0_dp,21.0_dp,dp)) error stop 'lower triangle changed'
if(zFock(1,2,1,1)/=cmplx(11201.0_dp,-21.0_dp,dp)) error stop 'upper triangle changed'
if(zFock(1,2,2,1)/=cmplx(12102.0_dp,12.0_dp,dp)) error stop 'neighbor matrix changed'
if(abs(coord(1,3)+tz/2)>1e-12_dp) error stop 'layer geometry changed'
print *, 'actual driver prologue PASS'
end program
'''
    compile_program(folder, source.stem, prefix, ['Hamiltonian.o'])
    run([str(folder / source.stem)], folder, f'{source.stem}.log')
    RESULTS.append(f'{source.stem}: actual geometry/Fock prologue PASS')

# The spectral kernels populate only an upper Hermitian triangle. Replacing
# their HEMM wrapper with the upstream GEMM wrapper would fail this check.
code = '''program check_lapack
use lapack_routines
implicit none
complex(8) :: upper(2,2),full(2,2),basis(2,2),actual(2,2),expected(2,2),signmatrix(2,2)
real(8) :: diagonal(2)
upper=0
upper(1,1)=2
upper(2,2)=-3
upper(1,2)=cmplx(0.25_8,0.5_8,8)
full=upper
full(2,1)=conjg(full(1,2))
basis=reshape([cmplx(1.0_8,0.0_8,8),cmplx(0.0_8,1.0_8,8), &
               cmplx(0.0_8,1.0_8,8),cmplx(1.0_8,0.0_8,8)],[2,2])/sqrt(2.0_8)
expected=matmul(conjg(transpose(basis)),matmul(full,basis))
call matrix_elements(upper,basis,actual,1_8,2_8,1_8,2_8)
if(maxval(abs(actual-expected))>1e-12_8) error stop 'upper Hermitian matrix elements'
call z_matrix_elements_diag_real(upper,basis,diagonal)
if(maxval(abs(diagonal-[real(expected(1,1),8),real(expected(2,2),8)]))>1e-12_8) error stop 'diagonal matrix elements'
signmatrix=full
call unitarize(signmatrix)
actual=matmul(signmatrix,signmatrix)
if(maxval(abs(actual-reshape([1,0,0,1],[2,2])))>1e-12_8) error stop 'unitarize'
print *, 'upper Hermitian matrix elements and added unitarize PASS'
end program
'''
compile_program(folder, 'check_lapack', code)
run([str(folder / 'check_lapack')], folder, 'check_lapack.log')
RESULTS.append('Upper-triangle Hermitian matrix elements and imported unitarize PASS')

# Kramers-Kronig consumes the optimized driver's single-record DiaXY header.
# A six-point grid (not the old hard-coded 30) exercises the saved grid value.
suffix = next((folder / 'dataFock').glob('Fock-nspin1*')).name[len('Fock-nspin1'):]
out = folder / 'output4'
out.mkdir(exist_ok=True)
(out / ('CBand' + suffix)).write_text('0.1 6 1000\n0\n')
for imu in range(1, 41):
    tail = f'1-Mu{imu}-numkPlot6-nspacing1000{suffix}'
    (out / ('DiaXY' + tail)).write_text('40 2 0.2 0.4 0.0\n')
    (out / ('ConductivityTot' + tail)).write_text(''.join(f'{(n-2)/1000} 1\n' for n in range(2, 41)))
compile_program(folder, 'Compute_KramersKronig', (POST / 'Compute_KramersKronig.f90').read_text(), ['Hamiltonian.o'])
run([str(folder / 'Compute_KramersKronig')], folder, 'Compute_KramersKronig.log')
for imu in range(1, 41):
    result = out / f'RealTot1-Mu{imu}-numkPlot6-nspacing1000-nstep1-nstepI1-nKK4{suffix}'
    rows = [list(map(float, line.split())) for line in result.read_text().splitlines()]
    assert len(rows) == 3 and all(math.isfinite(x) for row in rows for x in row)
    for (omega, value), k in zip(rows, range(2, 5)):
        # Constant input 1 is divided by 8/Ac before the transform and by Ac on output.
        integral = sum(((i-2)/1000)**2 / (((k-2)/1000)**2 - ((i-2)/1000)**2)
                       for i in range(2, 41) if i != k) / 8000
        expected = integral + math.pi/2 * (0.2+0.4)/2
        assert abs(value-expected) < 1e-8
RESULTS.append('Kramers-Kronig: optimized header, saved grid and transform normalization PASS')
for path in out.glob('ConductivityTot*'):
    path.write_text(''.join(f'{(n-2)/1000} 0\n' for n in range(2, 41)))
run([str(folder / 'Compute_KramersKronig')], folder, 'Compute_KramersKronig_zero.log')
for path in out.glob('RealTot*'):
    rows = [list(map(float, line.split())) for line in path.read_text().splitlines()]
    assert len(rows) == 3 and all(abs(row[1]-math.pi/2*0.3) < 1e-12 for row in rows)
RESULTS.append('Kramers-Kronig: zero input and computed-frequency output bounds PASS')

(BUILD / 'results.json').write_text(json.dumps(RESULTS, indent=2) + '\n')
print('\n'.join(RESULTS))
