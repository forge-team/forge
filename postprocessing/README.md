# FORGE post-processing integration

These programs adapt the supplied spectral, Drude and local-DOS sources to the geometry,
screening and packed Fock input of FORGE `dev`, commit
`acbddd097c689142a38e5b0e99ac5fd674c0904a`.
The former `utils/Main_BandStructure.f90` is now `compute_plotBands.f90`.
The band and order-parameter drivers have moved here alongside the new programs.

The geometry and input interface are aligned and tested. Screening uses the
solver's `HartreeFock.f90` directly. The reconstructed Hamiltonian agrees with
dev in the small supported-model cases tested, and the total current agrees
with its momentum derivative, including the Bloch-basis correction.
The spectral integration and diamagnetic formulas are retained; full production
spectral/Drude/local-DOS calculations have not been validated.

## Files and settings

| Supplied program | Working program |
| --- | --- |
| `compute_spectral.f90` | `Compute_Spectral.f90` |
| `compute_spectral_full.f90` | `Compute_SpectralFull.f90` |
| `compute_spectral_optimize.f90` | `Compute_SpectralOptimize.f90` |
| `compute_drude_full.f90` | `Compute_DrudeFull.f90` |
| `compute_kramers_kronig.f90` | `Compute_KramersKronig.f90` |
| `LocalDos.tar: plotDOS_MIT.f90` | `compute_localDos.f90` |
| `utils/Main_BandStructure.f90` | `compute_plotBands.f90` |
| `utils/Main_OrderParamsFull.f90` | `Main_OrderParamsFull.f90` |
| `utils/Main_OrderParamsLowEn.f90` | `Main_OrderParamsLowEn.f90` |
| Legacy raw layer-overlap export | `compute_OrderParameterLayers.f90` |
| Corrected interlayer order-parameter projection | `compute_LayerOrderParameter.f90` |

`postprocessing` describes these programs more precisely than `utils`.
`Compute_*` names the observable calculation; `Main.f90` remains the
self-consistent solver. The duplicate `utils/Setup.f90` and its Makefile are
replaced by the root Setup and the combined Makefile in this directory. The
original versions remain in Git history at the commit listed above.

Physical settings shared with the producer come directly from `../Setup.f90`:
dimensions, twist, spin count, cell cutoffs, relaxation, interaction parameters,
layer spacing, lattice vectors, Bloch phase and Fock directory. Both single-
and double-gate screening use `nscreen` and `xi` through the shared
`HartreeFock.fv(r_i,r_j)` function. This avoids a
second set of defaults such as the differing pressure values in the former
root and `utils` Setup files. The current imported producer defaults are
`ntheta=9`, `nrelax=0`, `pressure=1.134`; the supplied legacy module originally
used `ntheta=23`, `nrelax=1`, `pressure=1`.

Response-only settings remain in this directory's `Hamiltonian.f90`:
temperatures, chemical-potential scan size, output directory (`output4/`), and
the full response band window (`numb=ndim`, `ncb=ndim/2`). Each compute program
keeps its original integration grid, `numkPlot`. It is independent of the
producer's `numk`. The hard-coded interaction strengths in the driver prologues
now use the producer's `alpha` and `alphaH`. Local DOS also uses the producer's
selected band window (`Setup.numb/nLower/nUpper`); unlike the full response
programs, it does not use `Hamiltonian.numb=ndim`.

## Band structure and order parameters

`compute_plotBands.f90` replaces `utils/Main_BandStructure.f90`. Set the
observable grid `numkPostProc` (default 60, a multiple of six) and `nMirrorPath`
at the top of the driver. The saved state's grid `numk` and all physical
settings come from the root `Setup.f90`. Run the executable from the FORGE root.
It reads `dataFock/Fock-nspin<s><parameters>` and the final value in
`output/Mu<parameters>`, retaining the existing zero-shift fallback with a
warning when Mu is missing. Bands are written to
`output/Bands-nspin<s>-numkPostProc<grid><parameters>`; this replaces the old
`-numkPosProc` suffix spelling/position.

`Main_OrderParamsFull` evaluates the complete saved Fock matrix;
`Main_OrderParamsLowEn` reconstructs its selected-band and four-band
contributions using `OrderParameter.f90`. Their filenames are preserved.
Both use `BuildGeometry` and `ReadFock` and select the producer's current
`parameters`, rather than the restart-input settings. They write to `output/`.

`compute_LayerOrderParameter` implements the corrected inter-sublattice,
inter-valley layer projection in Eqs. (63)-(83) of `OrderParameterLayer`.
It follows the full saved-Fock driver, uses the same geometry/input helpers,
and writes four complex components for both layer directions. The local map is
`output/LayerOrderParameterLocal-numb<ndim>-nspin<s><parameters>`; normalized
cell sums and per-cell values are in
`output/LayerOrderParameter-numb<ndim>-nspin<s><parameters>`.
It retains the full Bloch sum phase, uses the opposite layer direction in the
conjugated mode, exchanges loop labels 2/3 at B centers, and averages the two
center grids. All centers are sampled, with no extra factor of three.
The valley labels agree with `OrderParameter.f90`.

The first two raw loop channels average three neighboring hexagons, with the
third retaining its central hexagon. This preserves the explicit valley phase
and reproduces the legacy spatial averaging. The equal-layer production kernel
agrees with the legacy inter-sublattice/inter-valley projection at every A
center after removing the legacy position-dependent valley phase, also for
a general periodic density matrix. The
executable currently selects only the two interlayer directions; the
equal-layer regression exercises its internal routines directly. See the
[equal-layer comparison](layer_orderparameter.md#equal-layer-comparison-with-main_orderparamsfull)
for the exact relation and test scope.

The geometry builder's optional reference coordinates preserve exact lattice
labels before relaxation, so lookup and Bloch phases use a fixed lattice gauge
while center pairing/output use physical positions. `--check-geometry` validates
all endpoints and saved translations without reading a Fock state.
See [layer_orderparameter.md](layer_orderparameter.md) for the formulas,
primitive-vector convention, output columns, normalization, and tests.

`compute_OrderParameterLayers` is retained as the earlier raw-overlap export
for existing workflows. Its original displaced-site convention differs from
the nearest-neighbor convention in the corrected derivation. Use
`compute_LayerOrderParameter` for the corrected projected density components.

## Common initialization

All ten Fock-consuming programs use the shared initialization routines.
For example, the response drivers use:

```fortran
call NumberNeighborCells(numI, ind)
allocate(zFock(ndim,ndim,ind,numS))
allocate(ind_j(ind),ind_l(ind))
call OrderNeighborCells(numI, ind, ind_j, ind_l)
call BuildGeometry(coord, t1, t2, t3, vq1, vq12, ang_mat)
! Shared long-range interaction and observable setup follow.
call InputParameters(parameters)
call ReadFock(zFock, ind, parameters)
```

`BuildGeometry` follows the cell, reciprocal-vector, rotation and relaxation
construction of the former Main_BandStructure. It uses upstream `WignerSeitzCell`, keeps
its site ordering and centered layer coordinates, and rotates only columns
1:2. Each observable retains its own momentum sampling. Nearest-neighbor arrays
are not needed by the retained legacy response Hamiltonian.

`ReadFock` uses the producer's output suffix, one `Fock-nspin<s>` file per spin,
the central-cell lower triangle with conjugate completion, and full neighbor
matrices in the upstream m/i/j order. The extra whole-matrix conjugation in the
old Drude reader is removed. Input is required; the solver's `nRead` flag does
not disable reading in a postprocessor. The filenames no longer contain the
old hard-coded absolute path or truncate the suffix to 90 characters.

The reader validates file size and read status. It retains upstream
`access='direct', recl=2*dp`. Compile producer and consumer with matching
precision, byte order and record-unit options. Intel's default unformatted
record units differ from GNU Fortran; a GNU build cannot silently read a
padded Intel file. The size check reports that mismatch. Do not rename an old
state into the new naming scheme unless its basis and format have been checked.

Lengths use `a=0.246 nm`; energies and interaction parameters use eV. `zFock`
is the dimensionless one-particle density matrix, not the Hamiltonian in eV.
The original response normalizations and chemical-potential prescriptions are
retained for the spectral/Drude programs. Local DOS reads the last value in
the producer's two-column `output/Mu<parameters>` history and requires that
file to correspond to the same saved state.

## Build and validation

From the FORGE root, using GNU Fortran and BLAS/LAPACK (Accelerate on macOS):

```sh
make -C postprocessing
make -C postprocessing check
```

For an Intel/MKL installation, select matching producer options, for example:

```sh
make -C postprocessing FC=ifort \
  FFLAGS='-qopenmp -heap-arrays -mcmodel=large' LIBS='-qmkl'
```

Use a fresh `.build` directory when switching compiler/options; Make does not
encode command-line compiler overrides in its dependency tracking. GNU builds
use the legacy `-fallow-argument-mismatch` option because the supplied BLAS
wrappers contain mixed integer kinds and scalar/array actual arguments.

Executables are under `postprocessing/.build/`. Run them from the FORGE root
so `dataFock/` resolves to the producer's output. Create `output4/` before
running, for example:

```sh
mkdir -p output output4
./postprocessing/.build/compute_plotBands
./postprocessing/.build/compute_LayerOrderParameter --check-geometry
./postprocessing/.build/compute_LayerOrderParameter
./postprocessing/.build/Compute_SpectralOptimize
./postprocessing/.build/Compute_KramersKronig
./postprocessing/.build/compute_localDos
```

These commands run the full observable calculations; they are not the small
interface tests. The full-band spectral arrays can be large.

The tests compile separate small Setup copies under `.build/tests/` and
`.build/screening-tests/`, plus `.build/localdos-tests/`, and check:

- Geometry against the Main_BandStructure block for unrelaxed, Koshino and Carr
  cases, with the reference's in-plane assignment corrected.
- Complex packed Fock round-trips using the actual upstream writer loops:
  single/double precision, one/two spin files and zero/one/two neighbor shells.
- Rejection of mismatched file size.
- The actual initialization blocks of all four Fock-consuming drivers.
- Upper-triangle Hermitian matrix elements and the added `unitarize` interface.
- Kramers–Kronig input headers, saved grid, normalization, computed-frequency
  range and a zero-conductivity case.
- Six cases compare the Hartree potential and reconstructed Hamiltonian to
  dev, and the sum of layer/interlayer current operators to finite differences
  of the dev Hamiltonian. They cover both gate configurations, both Bloch-phase
  conventions, all three relaxation choices and zero/one/two neighbor shells.
- The new layer-overlap driver's one-to-one periodic center matching, all loop
  endpoint lookups, and required Fock-cell translations in `--check-geometry`
  mode; this geometry-only check does not need a saved Fock file.

All eleven complete programs compile and link. The checks include 17 interface,
six Hamiltonian/current, seven local-DOS, eight moved-driver, twelve layer
projection, six layer-geometry cases, four equal-layer comparison fixtures,
and the two default layer geometry checks.
The layer projection tests use independent normalized Bloch spinors written as
packed Fock files. They check every A/B center, both layer directions, complex
phases, one/two spin sectors, gauge invariance, zero input, and cell/per-cell
normalization. The layer geometry checks also cover Nam/Carr relaxation and
reversed/equal layer rotations.
The equal-layer fixtures compare the production averaged kernel to the actual
legacy routine at all A centers in both layers, including boundaries, and run
`Main_OrderParamsFull` on the same packed Fock files. They also check constant
Dirac envelopes on interior A/B sites. The maximum equal-layer discrepancy
is below `4e-18` for these fixtures after the known valley-phase conversion.

The band driver matches full-BZ diagonalization on a small fine grid, for both
spin counts and both path orientations. The two legacy order-parameter drivers
have tested geometry/Fock prologues; the full saved-state driver's
inter-sublattice/inter-valley channel additionally has the comparison above.
Their other observable channels have not been independently validated. The maximum Hamiltonian discrepancy
was zero in the small cases checked, and the maximum current discrepancy was
1.2e-9 in code units with a finite-difference step of 1e-6. These tests use
synthetic states; the new layer projection still needs a production-state
convergence assessment of its smooth-envelope approximation.

The local-DOS checks cover the complete driver, its missing dependencies,
site-summed DOS, saved Mu, integrated maps and separate layer Fourier sums;
see [localdos_dependencies.md](localdos_dependencies.md).

## Retained routines and compatibility fixes

The local Hamiltonian supplies many procedures absent upstream, including
velocity/diamagnetic operators, conductivity, Drude, thermal-response and
Kramers–Kronig routines. The supplied modules are retained with the compatibility
changes below; the old scalar `fv` implementation has been removed.

`longRange` delegates to dev's `LongRangeInteraction`. All 55 exchange and
response interaction calls now pass both atomic positions to dev's `fv`,
including the translations for neighbor cells. No independent screening
formula or image-sum cutoff remains on the compute programs' interaction path.

The local LAPACK `matrix_elements` implementation uses the upper Hermitian
triangle, as the velocity routines require. Upstream's similarly named routine
uses the whole matrix. Replacing the module wholesale would change results.
The upstream `z_unitarize`/`c_unitarize` procedures and their generic interface
were added because `TightBinding` depends on them. The local complex BLAS
`ZDOTC` call crashed against macOS Accelerate in the matrix-element test;
the four diagonal-overlap wrappers now use equivalent Fortran `dot_product`.

The Kramers–Kronig consumer matches **Compute_SpectralOptimize** filenames and
its single-record `DiaXY` header. It now uses the grid stored in `CBand`, reads
both diamagnetic components from that one record, and writes only the frequency
range actually transformed. An all-zero conductivity now returns a zero
transform without scanning below the input array. The non-optimized spectral
programs keep their existing, different observable-output naming convention.

## Supported scope and remaining validation

The imported spectral/Drude/local-DOS kernels are bilayer, one-stored-spin-sector routines for the
legacy hopping model with zero layer bias. The drivers reject other layer/spin
counts, `TBFunction=2`, and nonzero `Delta`, rather than silently evaluating
only part of the requested model. The generic Fock reader and the migrated band/order-parameter drivers support
both spin counts. These drivers use dev's modules, built separately under `.build/upstream/`;
they are not subject
to the legacy response-model restriction.

The Hamiltonian exchange signs, conjugations and phases agree with dev for
the central cell and both translated-cell contributions. An earlier review
incorrectly reported a second-neighbor exchange discrepancy; the comparison
tests now check the entire reconstructed Hamiltonian in the supported cases.

The retained current and diamagnetic kernels have not been derived or
validated for the dev Wannier model or general multilayers. Full spectral
normalizations, chemical-potential scans and convergence still require a
production comparison.

## Density of states

General energy-resolved DOS routines are included (`dosFull`, `dos_F` and
related routines). The spectral drivers also retain their optional DOS output.
The subsequently supplied `LocalDos.tar` is integrated as
`compute_localDos.f90`. Its two missing LDOS routines are contained in the
driver. The extra `HamiltonianPlot.f90`, `Spectral.f90` and LAPACK copy are
unnecessary; only the existing shared modules are compiled. The dependency
audit and corrections are in [localdos_dependencies.md](localdos_dependencies.md).

Local DOS reads the current state selected by `../Setup.f90`, including the
selected band window (currently 20 bands). To reproduce the old four-band
selection, set `numb=4` and `ncb=2` in that Setup and rebuild. The observable
grids at the beginning of `compute_localDos.f90` retain the archive defaults:
`numkPlot=12`, `numkPlotRec=100`, `nspacing=100000` bins/eV (0.01 meV).
The arrays can be large for a wide band window; the driver prints their memory
requirement before allocation. Run from the FORGE root after creating
`output4/`, with the matching Fock files in `dataFock/` and Mu history in `output/`.

Output names append `-localDos-k<grid>-s<bins-per-eV><parameters>`:

| Prefix | Contents |
| --- | --- |
| `LocalDosInfo` | Band indices, CB/VB split, grids, saved Mu, peak energies and integrated counts. |
| `Dos`, `DosCheckCB`, `DosCheckVB` | Energy in meV relative to Mu, DOS in states/meV/cell for one stored spin sector. VB is written in descending energy; total DOS combines both sectors on an ascending axis. |
| `LDOS-CBmax-`, `LDOS-VBmax-` | Site x, y, LDOS at the corresponding sector's DOS maximum, in states/meV/site. |
| `LDOS-CBsum-`, `LDOS-VBsum-` | Site x, y, energy-integrated weight in states/site. |
| `DosBZ-<CB/VB><max/sum>-<Re/Im><1/2>` | Fourier qx, qy and real/imaginary part for the indicated layer. |
| `Bands`, `BandsAll` | Retained band-path outputs, in eV relative to Mu. |

Positions use `a=0.246 nm`; Fourier momenta use `1/a`. No extra spin or valley
factor is applied. CB/VB label the nominally unfilled/filled band groups at
`Setup.nfilling`, with any partially filled band assigned to CB; these labels
do not impose a cutoff in energy at Mu. The integration method is retained
from the archive, with a finite-weight treatment for flat/sub-bin triangles.
Full production grid convergence remains to be checked.

The original-source routine comparison is in [routine_inventory.md](routine_inventory.md);
its line numbers and callers refer to the supplied originals before adaptation.
