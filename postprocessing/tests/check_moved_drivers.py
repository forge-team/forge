#!/usr/bin/env python3
"""Exercise moved utils drivers with dev modules and isolated small inputs."""
from pathlib import Path
import json
import math
import os
import platform
import re
import subprocess

POST = Path(__file__).resolve().parents[1]
REPO = POST.parent
BUILD = POST / '.build/moved-driver-tests'
BUILD.mkdir(parents=True, exist_ok=True)
FLAGS = ['-O0', '-g', '-fopenmp', '-ffree-line-length-none', '-fallow-argument-mismatch', '-fcheck=all']
LIBS = ['-framework', 'Accelerate'] if platform.system() == 'Darwin' else ['-llapack', '-lblas']
ENV = dict(os.environ, OMP_NUM_THREADS='2')
RESULTS = []


def run(cmd, folder, name):
    p = subprocess.run(cmd, cwd=folder, env=ENV, text=True, capture_output=True)
    (folder / (name + '.log')).write_text(p.stdout + p.stderr)
    if p.returncode:
        raise RuntimeError(f'{folder.name}/{name}: {p.stderr[-2500:]}')


main = (REPO / 'Main.f90').read_text()
writer = re.search(r'rcc=0\b.*?close\(91\)', main[main.index('if (nWriteFock.EQ.1)'):], re.S).group()
for spins in (1, 2):
    folder = BUILD / f'spin{spins}'
    folder.mkdir(exist_ok=True)
    (folder / 'dataFock').mkdir(exist_ok=True)
    (folder / 'output').mkdir(exist_ok=True)
    setup = (REPO / 'Setup.f90').read_text()
    for name, value in dict(ntheta=1, numS=spins, numI=1, numC=1, numk=3, numb=4, nrelax=0).items():
        setup, count = re.subn(r'(^\s*integer(?:\(dp\))?,\s*parameter\s*::\s*' + name + r'\s*=)[^!\n]*',
                              rf'\g<1> {value} ', setup, flags=re.M | re.I)
        assert count == 1, name
    (folder / 'Setup.f90').write_text(setup)
    sources = [REPO / 'LapackRoutines.f90', folder / 'Setup.f90', REPO / 'TightBinding.f90',
               REPO / 'Geometry.f90', REPO / 'HartreeFock.f90', REPO / 'OrderParameter.f90', POST / 'PostProcessingInput.f90']
    objects = [p.stem + '.o' for p in sources]
    run(['gfortran', *FLAGS, '-c', *map(str, sources)], folder, 'modules')

    def compile_run(name, code):
        (folder / (name + '.f90')).write_text(code)
        run(['gfortran', *FLAGS, '-o', name, name + '.f90', *objects, *LIBS], folder, name + '_compile')
        run([str(folder / name)], folder, name)

    fixture = '''program fixture
use Setup
use Geometry, only: NumberNeighborCells
implicit none
integer(dp) :: i,j,m,nspin,rcc,numNeighborCells,my_iostat
complex(dp),allocatable :: zFock(:,:,:,:)
character(1024) :: filename
character(150) :: my_iomsg
call InitParameters()
call NumberNeighborCells(numI,numNeighborCells)
allocate(zFock(ndim,ndim,numNeighborCells,numS))
zFock=0
do nspin=1,numS
    do i=1,ndim
        zFock(i,i,1,nspin)=.5_dp+.001_dp*nspin*i
        do j=1,i-1
            zFock(i,j,1,nspin)=cmplx(.0001_dp*nspin*(i+j),.0001_dp*(i-j),dp)
            zFock(j,i,1,nspin)=conjg(zFock(i,j,1,nspin))
        enddo
    enddo
    zFock(:,:,2:,nspin)=cmplx(.0002_dp*nspin,-.0003_dp,dp)
    write(filename,'(A,A,I0,A)') trim(dirFock),'Fock-nspin',nspin,trim(parameters)
    open(91,file=filename,status='replace',form='unformatted',access='direct',recl=2*dp)
'''+writer+'''
enddo
filename=trim(dir)//'Mu'//trim(parameters)
open(13,file=filename,status='replace')
write(13,*) 1, -.1_dp
write(13,*) 2, .123_dp
close(13)
end program
'''
    compile_run('fixture', fixture)
    driver = (POST / 'compute_plotBands.f90').read_text().replace('numkPostProc = 60', 'numkPostProc = 6')
    for mirror in (0, 1):
        code = driver.replace('nMirrorPath = 0', f'nMirrorPath = {mirror}')
        compile_run('bands', code)
        files = sorted((folder / 'output').glob('Bands-nspin*'))
        assert len(files) == spins
        actual = [p.read_text() for p in files]
        # Independent reference: diagonalize the entire fine BZ, then let dev's
        # PlotBands select the path. This catches missing/mis-indexed path points.
        reference = code[:code.index('!! Grid indices')]
        reference += '''
allocate(zH(ndim,ndim),Energies(ndim))
do nspin=1,numS
do ivk1=0,numkPostProc-1
do ivk2=0,numkPostProc-1
    jvk1=ivk1; jvk2=ivk2
    if(nMirrorPath==1)then
        jvk1=modulo(-ivk2,numkPostProc)
        jvk2=modulo(-ivk1,numkPostProc)
    endif
    vk=MomentaValues(jvk1+1,jvk2+1,:)
    call HamiltonianHartreeFock(zH,Coords,Potential(:,nspin),alpha,Delta,zFock(:,:,:,nspin), &
         nUnitCell_1,nUnitCell_2,ndim,numNeighborCells,vk,tnHF,NearestNeighborsUC,NearestNeighborsT)
    call diagonalize(zH,Energies,'N',nLower,nUpper)
    Bands(:,nMomentaFlattened(ivk1+1,ivk2+1),nspin)=Energies(1:numb)
enddo
enddo
enddo
'''+code[code.index('!! Output bands'):]
        compile_run('reference', reference)
        for p, saved in zip(files, actual):
            values = [float(x) for x in saved.split()]
            expected = [float(x) for x in p.read_text().split()]
            assert len(values) == len(expected) and values
            assert all(math.isfinite(x) and abs(x-y) < 1e-10 for x, y in zip(values, expected))
        RESULTS.append(f'plotBands spin={spins}, mirror={mirror}: fine-grid path equals full-BZ reference PASS')

    for name in ('Main_OrderParamsFull', 'Main_OrderParamsLowEn'):
        code = (POST / (name + '.f90')).read_text()
        marker = 'call ReadFock(zFock,numNeighborCells,inputSuffix)'
        prefix = code[:code.index(marker)+len(marker)]
        prefix += '''
do nspin=1,numS
if(abs(real(zFock(2,2,1,nspin),dp)-(.5_dp+.002_dp*nspin))>1e-14_dp) error stop 'density'
if(abs(zFock(2,1,1,nspin)-cmplx(.0003_dp*nspin,.0001_dp,dp))>1e-14_dp) error stop 'lower'
if(abs(zFock(1,2,1,nspin)-conjg(zFock(2,1,1,nspin)))>1e-14_dp) error stop 'upper'
if(abs(zFock(1,2,2,nspin)-cmplx(.0002_dp*nspin,-.0003_dp,dp))>1e-14_dp) error stop 'neighbor'
enddo
if(trim(inputSuffix)/=trim(parameters)) error stop 'producer suffix'
if(abs(Coords(1,3)+tz/2)>1e-12_dp) error stop 'layer geometry'
end program
'''
        compile_run(name, prefix)
        RESULTS.append(f'{name} spin={spins}: actual geometry and packed Fock prologue PASS')

(BUILD / 'results.json').write_text(json.dumps(RESULTS, indent=2)+'\n')
print('\n'.join(RESULTS))
