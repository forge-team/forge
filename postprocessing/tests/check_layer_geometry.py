#!/usr/bin/env python3
"""Check ideal loop labels through physical relaxation and layer rotations."""
from pathlib import Path
import json
import os
import platform
import re
import subprocess

POST = Path(__file__).resolve().parents[1]
REPO = POST.parent
BUILD = POST / '.build/layer-geometry-tests'
BUILD.mkdir(parents=True, exist_ok=True)
FLAGS = ['-O2', '-fopenmp', '-ffree-line-length-none', '-fcheck=all', '-fbacktrace']
LIBS = ['-framework', 'Accelerate'] if platform.system() == 'Darwin' else ['-llapack', '-lblas']
ENV = dict(os.environ, OMP_NUM_THREADS='2')
CASES = [('unrelaxed', 2, 0, '[-1,+1]'), ('nam', 2, 1, '[-1,+1]'),
         ('reversed', 2, 0, '[+1,-1]'), ('both_negative', 2, 0, '[-1,-1]'),
         ('both_positive', 2, 0, '[+1,+1]'), ('carr', 20, 2, '[-1,+1]')]
results = []
for name, angle, relax, rotations in CASES:
    work = BUILD / name
    work.mkdir(exist_ok=True)
    setup = (REPO / 'Setup.f90').read_text()
    for parameter, value in dict(ntheta=angle, nrelax=relax, nlayers=2, numI=1).items():
        setup, count = re.subn(r'(^\s*integer(?:\(dp\))?,\s*parameter\s*::\s*' + parameter + r'\s*=)[^!\n]*',
                              rf'\g<1> {value} ', setup, flags=re.M | re.I)
        assert count == 1, parameter
    setup, count = re.subn(r'(RotateLayers\(nlayers\)\s*=\s*)\[[^\]]+\]',
                          lambda m: m[1] + rotations, setup)
    assert count == 1
    (work / 'Setup.f90').write_text(setup)
    command = ['gfortran', *FLAGS, '-o', 'driver', str(REPO / 'LapackRoutines.f90'), 'Setup.f90',
               str(REPO / 'TightBinding.f90'), str(REPO / 'Geometry.f90'),
               str(POST / 'PostProcessingInput.f90'), str(POST / 'compute_LayerOrderParameter.f90'), *LIBS]
    for label, args in [('compile', command), ('geometry', ['./driver', '--check-geometry'])]:
        p = subprocess.run(args, cwd=work, env=ENV, text=True, capture_output=True, timeout=120)
        (work / f'{label}.log').write_text(p.stdout + p.stderr)
        if p.returncode:
            raise RuntimeError(f'{name}/{label}: {p.stdout[-1000:]}\n{p.stderr[-2000:]}')
    assert 'geometry check passed' in p.stdout
    assert not (work / 'dataFock').exists() and not (work / 'output').exists()
    results.append(f'{name}, ntheta={angle}, nrelax={relax}, RotateLayers={rotations}: '
                   'reference bond topology and saved translation lookup PASS')
(BUILD / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
print('\n'.join(results))
