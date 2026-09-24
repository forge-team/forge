# Local-DOS integration

Source: `LocalDos.tar`, supplied 24 September 2026. The archive is preserved.
The working program is `compute_localDos.f90`, adapted from `plotDOS_MIT.f90`.
The tar also contains a duplicate self-hardlink entry for that program; it is
ignored after reading the regular file.

## Dependency decisions

| Archive file/call | Resolution |
| --- | --- |
| `HamiltonianPlot.f90` (module `Hamiltonian`, 46 procedures) | Use the existing shared `Hamiltonian.f90`. Of the 46 procedure names, 44 are already present there. The two exceptions are explained below. |
| Legacy scalar `fv(r)` | Use dev's position-dependent `HartreeFock.fv(r_i,r_j)`, already called by the shared Hamiltonian. Screening remains controlled by `Setup.f90`. |
| Legacy `dos` wrapper | Its independently configured DOS pass is unnecessary. `Dos` is the sum of the CB and VB site-resolved results on the same energy grid. |
| `LapackRoutines.f90` (28 procedures) | All concrete procedure names are already in the shared LAPACK module. Keep the existing compatibility fixes. |
| `Spectral.f90` (11 procedures) | Only `LocalDOS_CB` and `LocalDOS_VB` are needed. Both are contained in `compute_localDos.f90`; they call the existing `Hamiltonian.sortF`. The other nine procedures are unused by this driver. |
| `nearestNeighbour`, `kekuleCenter`, `kekule`, `kekuleOrder` | Unresolved in the supplied package/current shared tree. Their results are unused by local DOS, so their calls and associated arrays are removed. |
| `samplePointsRec` | Also unresolved. Replace its grid-index enumeration with explicit rectangular Fourier-grid loops. |
| Old geometry, interaction, and Fock initialization | Use `BuildGeometry`, `NumberNeighborCells`, `OrderNeighborCells`, `ReadFock`, and the shared Hamiltonian/interaction routines. |

No separate `HamiltonianPlot.f90`, second `LapackRoutines.f90`, or `Spectral.f90`
is installed in this directory. The retained dependencies are `Setup`,
`Geometry`, `TightBinding`, `HartreeFock`, the existing postprocessing
`Hamiltonian`/`lapack_routines`, `PostProcessingInput`, OpenMP and BLAS/LAPACK.
The complete executable links successfully with the existing Makefile.

The archive's `ComplexHk_MF` agrees with the shared implementation except for
the intentional use of dev's screening function and inactive comments. It
does not require an additional active hopping or exchange routine. The
shared `sortF` generalizes the archive's one-element temporary arrays to
the input size; its sorting logic is unchanged.

## Input and numerical corrections

- Read the current producer suffix and packed complex Fock data using the
  shared reader, without the archive's extra complex conjugation. Subtract
  the neutral density background when reconstructing the Hartree potential.
- Read the final iteration from `Setup.dir/Mu<parameters>`, matching dev's
  two-column output. Missing, empty, malformed, or incomplete histories stop
  with a diagnostic. Old one-value `Mu-nspin` files are not the current format.
- Use `Setup.numb/nLower/nUpper` and physical parameters, replacing the
  hard-coded four-band, filling and interaction scans. CB/VB are band groups
  split at the nominal filled-band count from `Setup.nfilling`; a partially
  occupied band belongs to CB. They are not energy cutoffs at Mu.
- Initialize the entire independent `numkPlot` grid, including the boundary
  used by the integration. Set energy limits from all selected bands with
  guard bins, including when all selected bands lie on one side of Mu.
- Preserve the archive's two-triangle LDOS interpolation. An exactly flat
  triangle or one wholly within a single energy bin deposits its integrated
  weight in that bin, avoiding division by zero or lost weight.
- Normalize DOS to states/meV per cell and local DOS to states/meV per site,
  for one stored spin sector. Integrated maps multiply by the bin width in
  **meV**; this corrects the original factor-of-1000 loss.
- Compute each layer's Fourier sum separately. The original integrated-map
  code forgot to reset its accumulator before summing layer two.
- Wrap periodic indices in shared `plotBandsAll`. The bounds-checked local-DOS
  run exposed the old access at `numk+1` into a `numk`-sized array. This also
  benefits the Drude driver that calls the same procedure.

The two-dimensional LDOS arrays are shared between threads; each thread owns
different sites. This removes the original full-array OpenMP reductions.
The driver prints the memory required by those arrays before allocation.

## Validation and limits

`tests/check_localdos.py` builds isolated small fixtures, with a 28-site cell,
four selected bands, a producer grid of 3 and a postprocessing grid of 6.
It uses the actual upstream Fock-writer loops and runs with bounds checks,
floating-point traps, and two OpenMP threads. Checks cover:

- Site sums against the existing unweighted CB/VB DOS kernels.
- Finite, normalized flat and sub-bin triangle contributions.
- The complete driver, saved final Mu, peak and integrated maps.
- Fourier transforms recomputed independently for each layer at every point
  of the small output grid.
- Total DOS as CB plus VB, and required-Mu failure cases.

The full production local-DOS calculation and grid convergence have not been
run. The retained rounded-bin interpolation is approximate; the synthetic
two-band integrated counts are approximately 2.00009 (CB) and 1.99945 (VB).
Support follows the other imported response kernels: bilayer, one stored spin
sector, Slater–Koster hopping, and zero layer bias.
