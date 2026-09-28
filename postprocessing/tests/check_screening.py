#!/usr/bin/env python3
"""Compare postprocessing Hartree, Hamiltonian and current with FORGE dev."""
from pathlib import Path
import json
import os
import platform
import re
import subprocess

POST = Path(__file__).resolve().parents[1]
REPO = POST.parent
BUILD = POST / '.build' / 'screening-tests'
FLAGS = ['-O0', '-g', '-fopenmp', '-ffree-line-length-none', '-fallow-argument-mismatch', '-fcheck=all']
LIBS = ['-framework', 'Accelerate'] if platform.system() == 'Darwin' else ['-llapack', '-lblas']
ENV = dict(os.environ, OMP_NUM_THREADS='1')
results = []
cases = [
    dict(nscreen=1, numI=0, nrelax=0, fphase=0),
    dict(nscreen=2, numI=1, nrelax=0, fphase=0),
    dict(nscreen=1, numI=2, nrelax=1, fphase=1),
    dict(nscreen=2, numI=2, nrelax=1, fphase=1),
    dict(nscreen=1, numI=1, nrelax=2, fphase=0),
    dict(nscreen=2, numI=1, nrelax=2, fphase=1),
]
for case in cases:
    name = '_'.join(f'{k}{v}' for k,v in case.items())
    folder = BUILD / name
    folder.mkdir(parents=True, exist_ok=True)
    setup = (REPO / 'Setup.f90').read_text()
    for key, value in dict(dp=8, ntheta=1, numS=1, numC=2, **case).items():
        setup, count = re.subn(r'(^\s*integer(?:\(dp\))?,\s*parameter\s*::\s*'+key+r'\s*=)[^!\n]*',
                              rf'\g<1> {value} ', setup, flags=re.M|re.I)
        assert count == 1, key
    (folder / 'Setup.f90').write_text(setup)
    sources = [POST/'LapackRoutines.f90', folder/'Setup.f90', REPO/'TightBinding.f90',
               REPO/'Geometry.f90', REPO/'HartreeFock.f90', POST/'Hamiltonian.f90',
               POST/'PostProcessingInput.f90', POST/'tests/check_screening.f90']
    compile_result = subprocess.run(['gfortran', *FLAGS, '-o', 'check_screening', *map(str,sources), *LIBS],
                                    cwd=folder, env=ENV, text=True, capture_output=True)
    (folder/'compile.log').write_text(compile_result.stdout+compile_result.stderr)
    if compile_result.returncode:
        raise RuntimeError(f'{name}: {compile_result.stderr[-3500:]}')
    result = subprocess.run([str(folder/'check_screening')], cwd=folder, env=ENV, text=True, capture_output=True)
    (folder/'run.log').write_text(result.stdout+result.stderr)
    if result.returncode:
        raise RuntimeError(f'{name}: {result.stdout}\n{result.stderr[-2000:]}')
    results.append(f'{name}: shared screening, dev Hamiltonian and current derivative PASS; {result.stdout.strip()}')
    print(results[-1], flush=True)
(BUILD/'results.json').write_text(json.dumps(results,indent=2)+'\n')
