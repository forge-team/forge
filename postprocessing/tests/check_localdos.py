#!/usr/bin/env python3
"""Small local-DOS kernel and end-to-end tests, isolated from production data."""
from pathlib import Path
import json
import math
import os
import platform
import re
import subprocess

POST = Path(__file__).resolve().parents[1]
REPO = POST.parent
BUILD = POST / '.build/localdos-tests'
BUILD.mkdir(parents=True, exist_ok=True)
FLAGS = ['-O0', '-g', '-fopenmp', '-ffree-line-length-none', '-fallow-argument-mismatch',
         '-fcheck=all', '-ffpe-trap=invalid,zero,overflow']
LIBS = ['-framework', 'Accelerate'] if platform.system() == 'Darwin' else ['-llapack', '-lblas']
ENV = dict(os.environ, OMP_NUM_THREADS='2')
RESULTS = []


def run(command, name, success=True):
    result = subprocess.run(command, cwd=BUILD, env=ENV, text=True, capture_output=True)
    (BUILD / (name + '.log')).write_text(result.stdout + result.stderr)
    if success and result.returncode:
        raise RuntimeError(f'{name}: {result.stdout[-1000:]}\n{result.stderr[-2500:]}')
    return result


def check(label):
    RESULTS.append(label + ' PASS')
    print(RESULTS[-1], flush=True)


setup = (REPO / 'Setup.f90').read_text()
for name, value in dict(dp=8, ntheta=1, numS=1, numI=1, numC=1, numk=3, numb=4, nrelax=0).items():
    setup, count = re.subn(r'(^\s*integer(?:\(dp\))?,\s*parameter\s*::\s*' + name + r'\s*=)[^!\n]*',
                          rf'\g<1> {value} ', setup, flags=re.M | re.I)
    assert count == 1, name
(BUILD / 'Setup.f90').write_text(setup)
sources = [POST / 'LapackRoutines.f90', BUILD / 'Setup.f90', REPO / 'TightBinding.f90',
           REPO / 'Geometry.f90', REPO / 'HartreeFock.f90', POST / 'Hamiltonian.f90', POST / 'PostProcessingInput.f90']
objects = [p.stem + '.o' for p in sources]
run(['gfortran', *FLAGS, '-c', *map(str, sources)], 'compile_modules')
driver = (POST / 'compute_localDos.f90').read_text()
kernels = driver.split('! BEGIN imported local-DOS kernels\n')[1].split('! END imported local-DOS kernels')[0]
(BUILD / 'kernels.f90').write_text('module kernels\nuse Setup, only: dp\nuse Hamiltonian, only: sortF\n'
                                  'implicit none\ncontains\n' + kernels + '\nend module\n')
(BUILD / 'check_kernels.f90').write_text('''program check_kernels
use Setup, only: dp
use Hamiltonian, only: dosCB_S, dosVB_S
use kernels
implicit none
integer(dp), parameter :: bins=501, rate=1000, shift=201
integer(dp) :: site, band, vertex
integer(dp) :: bounds(4)=[1_dp,1_dp,2_dp,2_dp]
real(dp) :: k(2,4),energy(2,4),psi(2,4),local(bins),actual(bins),reference(bins)
k=reshape([0._dp,0._dp,1._dp,0._dp,0._dp,1._dp,1._dp,1._dp],[2,4])
energy(1,:)=[-.12_dp,-.02_dp,-.07_dp,.03_dp]
energy(2,:)=[.02_dp,.08_dp,.14_dp,.20_dp]
actual=0
do site=1,4
    do band=1,2
    do vertex=1,4
        psi(band,vertex)=real(site+vertex+band,dp)/(10+4*(vertex+band))
    enddo
    enddo
    call LocalDOS_CB(4_dp,local,k,energy,psi,bins,rate,shift,0._dp,bounds)
    actual=actual+local
enddo
call dosCB_S(reference,k,energy,bins,rate,shift,0._dp,bounds)
if(maxval(abs(actual-reference))>1e-12_dp) error stop 'CB site sum differs from DOS'
actual=0
do site=1,4
    do band=1,2
    do vertex=1,4
        psi(band,vertex)=real(site+vertex+band,dp)/(10+4*(vertex+band))
    enddo
    enddo
    call LocalDOS_VB(4_dp,local,k,energy,psi,bins,rate,shift,0._dp,bounds)
    actual=actual+local
enddo
call dosVB_S(reference,k,energy,bins,rate,shift,0._dp,bounds)
if(maxval(abs(actual-reference))>1e-12_dp) error stop 'VB site sum differs from DOS'
! Fully flat and unresolved sub-bin triangles must preserve their weight.
energy(1,:)=-.1_dp
energy(2,:)=.1_dp
psi=.25_dp
call LocalDOS_CB(4_dp,local,k,energy,psi,bins,rate,shift,0._dp,bounds)
if(abs(sum(local)/rate-.25_dp)>1e-12_dp) error stop 'flat CB weight'
call LocalDOS_VB(4_dp,local,k,energy,psi,bins,rate,shift,0._dp,bounds)
if(abs(sum(local)/rate-.25_dp)>1e-12_dp) error stop 'flat VB weight'
energy(2,:)=[.10001_dp,.10002_dp,.10003_dp,.10004_dp]
call LocalDOS_CB(4_dp,local,k,energy,psi,bins,rate,shift,0._dp,bounds)
if(abs(sum(local)/rate-.25_dp)>1e-12_dp) error stop 'sub-bin CB weight'
print *, 'spatial sum, flat and sub-bin triangle checks passed'
end program
''')
run(['gfortran', *FLAGS, '-o', 'check_kernels', 'kernels.f90', 'check_kernels.f90', *objects, *LIBS], 'compile_kernels')
run([str(BUILD / 'check_kernels')], 'check_kernels')
check('CB/VB local-DOS site sums equal the existing total-DOS kernels')
check('Flat and sub-bin triangles preserve finite integrated weight')

# Actual upstream writer loops, a Hermitian complex synthetic density matrix,
# and a two-iteration Mu file. The postprocessing grid deliberately differs
# from Setup.numk; layer densities differ to expose Fourier layer mixing.
main = (REPO / 'Main.f90').read_text()
writer = re.search(r'rcc=0\b.*?close\(91\)', main[main.index('if (nWriteFock.EQ.1)'):], re.S).group()
(BUILD / 'fixture.f90').write_text('''program fixture
use Setup
use Geometry, only: NumberNeighborCells
use PostProcessingInput, only: InputParameters
implicit none
integer(dp) :: i,j,m,nspin,rcc,numNeighborCells,my_iostat
complex(dp),allocatable :: zFock(:,:,:,:)
complex(dp) :: zinput
character(1024) :: filename
character(150) :: my_iomsg
call InitParameters()
call NumberNeighborCells(numI,numNeighborCells)
allocate(zFock(ndim,ndim,numNeighborCells,numS))
zFock=0
do i=1,ndim
    zFock(i,i,1,1)=.5_dp+.01_dp*cos(real(i,dp))
    do j=1,i-1
        zFock(i,j,1,1)=cmplx(.0001_dp*cos(real(i+j,dp)),.0001_dp*sin(real(i-j,dp)),dp)
        zFock(j,i,1,1)=conjg(zFock(i,j,1,1))
    enddo
enddo
do m=2,numNeighborCells
do i=1,ndim
do j=1,ndim
    zFock(i,j,m,1)=cmplx(.0001_dp*cos(real(i+j+m,dp)),.0001_dp*sin(real(i-j+m,dp)),dp)
enddo
enddo
enddo
do nspin=1,numS
    filename=trim(dirFock)//'Fock-nspin1'//trim(parameters)
    open(91,file=filename,status='replace',form='unformatted',access='direct',recl=2*dp)
''' + writer + '''
enddo
filename=trim(dir)//'Mu'//trim(parameters)
open(71,file=filename,status='replace')
write(71,*) 1.0_dp,0.12_dp
write(71,*) 2.0_dp,0.23_dp
close(71)
end program
''')
for directory in ['dataFock', 'output', 'output4']:
    (BUILD / directory).mkdir(exist_ok=True)
run(['gfortran', *FLAGS, '-o', 'fixture', 'fixture.f90', *objects, *LIBS], 'compile_fixture')
run([str(BUILD / 'fixture')], 'fixture')
small = driver.replace('numkPlot = 12, numkPlotRec = 100, nspacing = 100000',
                       'numkPlot = 6, numkPlotRec = 2, nspacing = 1000')
assert small != driver
(BUILD / 'compute_localDos.f90').write_text(small)
run(['gfortran', *FLAGS, '-o', 'compute_localDos', 'compute_localDos.f90', *objects, *LIBS], 'compile_driver')
run([str(BUILD / 'compute_localDos')], 'run_driver')
check('Full driver: shared geometry/Fock, saved Mu, independent k-grid, bounds/FPE checks and two OpenMP threads')

out = BUILD / 'output4'
suffix = '-localDos-k6-s1000' + next((BUILD / 'dataFock').glob('Fock-nspin1*')).name[len('Fock-nspin1'):]


def read(label):
    rows = [[float(x) for x in line.split()] for line in (out / (label + suffix)).read_text().splitlines()]
    assert rows and all(math.isfinite(x) for row in rows for x in row), label
    return rows


info = (out / ('LocalDosInfo' + suffix)).read_text()
assert abs(float(re.search(r'Mu_eV:\s*(\S+)', info).group(1)) - .23) < 1e-12
for sector in ['CB', 'VB']:
    dos = read('DosCheck' + sector)
    peak = read('LDOS-' + sector + 'max-')
    integrated = read('LDOS-' + sector + 'sum-')
    assert len(peak) == len(integrated) == 28
    assert abs(sum(row[2] for row in peak) - max(row[1] for row in dos)) < 1e-11
    denergy = abs(dos[1][0] - dos[0][0])
    weight = sum(row[2] for row in integrated)
    assert abs(weight - sum(row[1] for row in dos) * denergy) < 1e-11
    # Two bands in each sector for this neutral four-band fixture. The
    # retained rounded-bin triangle rule is approximate; check convergence
    # scale, not exact integer state count at this coarse resolution.
    assert abs(weight - 2) < .02, (sector, weight)
    for kind, values in [('max', peak), ('sum', integrated)]:
        for layer in [1, 2]:
            start, stop = (layer-1)*14, layer*14
            real_rows = read(f'DosBZ-{sector}{kind}-Re{layer}')
            imag_rows = read(f'DosBZ-{sector}{kind}-Im{layer}')
            assert len(real_rows) == len(imag_rows) == 25
            for real_row, imag_row in zip(real_rows, imag_rows):
                qx, qy = real_row[:2]
                expected_re = sum(v*math.cos(qx*x+qy*y) for x,y,v in values[start:stop])
                expected_im = sum(v*math.sin(qx*x+qy*y) for x,y,v in values[start:stop])
                assert abs(real_row[2]-expected_re) < 1e-11
                assert abs(imag_row[2]-expected_im) < 1e-11
check('Saved final Mu, site-summed DOS, peak maps and integrated-map meV normalization')
check('All Fourier grid points match independent sums over each layer')
total, cb, vb = read('Dos'), read('DosCheckCB'), read('DosCheckVB')[::-1]
for a,b,c in zip(total,cb,vb):
    assert abs(a[0]-b[0]) < 1e-10 and abs(a[0]-c[0]) < 1e-10
    assert abs(a[1]-b[1]-c[1]) < 1e-12
check('Total DOS equals CB plus VB on the same energy axis')

# Required inputs must fail clearly instead of silently shifting by zero.
mu = next((BUILD / 'output').glob('Mu*'))
saved = mu.read_text()
mu.rename(mu.with_suffix('.saved'))
missing = run([str(BUILD / 'compute_localDos')], 'missing_mu', success=False)
assert missing.returncode and 'saved Mu history' in missing.stderr
mu.write_text('')
empty = run([str(BUILD / 'compute_localDos')], 'empty_mu', success=False)
assert empty.returncode and 'Empty Mu history' in empty.stderr
mu.write_text('invalid record\n')
malformed = run([str(BUILD / 'compute_localDos')], 'malformed_mu', success=False)
assert malformed.returncode and 'Malformed Mu history' in malformed.stderr
mu.write_text(saved + '3.0\n')
incomplete = run([str(BUILD / 'compute_localDos')], 'incomplete_mu', success=False)
assert incomplete.returncode and 'Malformed Mu history' in incomplete.stderr
mu.write_text(saved)
check('Missing, empty, malformed and incomplete saved Mu inputs are rejected')
(BUILD / 'results.json').write_text(json.dumps(RESULTS, indent=2) + '\n')
