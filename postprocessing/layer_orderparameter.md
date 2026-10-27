# Interlayer order parameters

`compute_LayerOrderParameter.f90` follows the initialization and packed-Fock
input of `Main_OrderParamsFull.f90`. It implements the corrected
**inter-sublattice, inter-valley** channel of `OrderParameterLayer`,
Eqs. (63)-(83), for both ordered layer pairs. It computes these four complex
envelope products and their cell sums:

| Label | Local product |
| --- | --- |
| `AK2_BKp1` | `conjg(f_AK2) * f_BKprime1` |
| `BK2_AKp1` | `conjg(f_BK2) * f_AKprime1` |
| `AK1_BKp2` | `conjg(f_AK1) * f_BKprime2` |
| `BK1_AKp2` | `conjg(f_BK1) * f_AKprime2` |

Other sublattice/valley channels are outside this driver's scope. No symmetry
relation between the `12` and `21` components is imposed.

## Build and run

From the FORGE root:

```sh
make -C postprocessing .build/compute_LayerOrderParameter
./postprocessing/.build/compute_LayerOrderParameter --check-geometry
./postprocessing/.build/compute_LayerOrderParameter
```

Use the root `Setup.f90` that describes the saved state. The program calls
`BuildGeometry`, `InputParameters`, and `ReadFock`; all state parameters,
screening settings, relaxation, precision, spin sectors and paths come from
that Setup. It consumes `dataFock/Fock-nspin<s><parameters>`. It uses the saved
real-space density directly, so it does not reconstruct the Hamiltonian or
read a chemical potential. Run with an existing `output/` directory, as for
the other order-parameter drivers.

The geometry check validates one-to-one center pairing, every endpoint on its
expected sublattice, and every required saved Fock translation. It neither
allocates the full Fock array nor reads/writes observable files. Insufficient
`numI` causes an explicit error, rather than discarding the missing overlaps;
the calculation needs a saved state with sufficient translations, not a renamed
input file. The driver requires a bilayer and `ntheta >= 1`.

## Geometry and conventions

Lengths use `a = 0.246 nm`, momenta use `1/a`, and densities are dimensionless.
The valley convention agrees with `OrderParameter.f90`:
`K0 = (-4*pi/3, 0)` and `Kprime = -K`. All vectors are rotated with their layer,
using `RotateLayers` and the shared geometry's rotation matrix.

The derivation's primitive vectors are `(Setup a2, Setup a1)`. In terms of
Setup's vectors, the three nearest-neighbor bonds are

```text
d2 = (a1 + a2)/3
d1 = d2 - a1
d3 = d2 - a2
```

They satisfy `exp(i*K.d1)=omega`, `exp(i*K.d2)=1`,
`exp(i*K.d3)=conjg(omega)`, with `omega=exp(2*pi*i/3)`. This mapping is essential:
the original `compute_OrderParameterLayers` raw export used a different
displaced-site convention and is not the source of the new valley projection.

The two sets of same-sublattice centers are paired one-to-one by minimum
physical displacement over periodic images, as in the existing raw driver.
`BuildGeometry` optionally returns unrelaxed `ReferenceCoords` before applying
relaxation. The new driver uses that exact reference lattice for endpoint
lookup and the valley phase, and physical coordinates for pairing and output.
This defines the envelope amplitudes in a fixed reference-lattice gauge also
for relaxed states. No nearest-atom approximation on distorted bond vectors is
used. Site ordering and all Fock indices remain those of the producer.

The complete phase, including the paired center's periodic image, is

```text
chi_i = exp(i * (K1.reference_r1 + K2.reference_r2))
```

It is the same in both layer directions. The geometry's origin is not an A
site, so the phase must not be set to one. A common phase of the wavefunction
cancels from every observable.

The producer's Fock convention is `D(i,j,R)=<psi_i^* psi_j> exp(-i*k.R)` before
the momentum sum. For two translated sites the required cell index is
`rowCell-colCell`. If only `-R` is stored, the overlap is obtained as
`conjg(D(j,i,-R))`. These endpoint and cell lookups are cached before reading
the state and reused for all spin sectors.

## Projection and normalization

Compute the six-term loops `Delta21(m)` and `Delta12(m)` for `m=1,2,3`.
Before projection, average the first two raw channels over three equivalent
positions, retaining the central third channel. In each layer's rotated
Setup vectors, use the same integer shifts for both bra and ket endpoints:

| Raw channel | A-center shifts | Weight per loop |
| --- | --- | --- |
| 1 | `0`, `a1+a2`, `2*a1-a2` | `1/3` |
| 2 | `0`, `2*a1-a2`, `a1-2*a2` | `1/3` |
| 3 | `0` | `1` |

For B centers negate these shifts. Each added shift has `K.shift` equal to
an integer multiple of `2*pi` in its own layer, so the valley phase is retained
at the original matched centers. Both ordered layer directions use the same
rule. Cache all seven sampled hexagons; missing Fock translations still cause
an explicit error. Constant-envelope results and their normalization are
unchanged by this averaging.

For an A center keep the averaged loop order; for a B center use order `(1,3,2)`.
For each ordered layer pair, the two projected modes are

```text
a = (omega*L1 + L2 + conjg(omega)*L3)/3
b = (conjg(omega)*L1 + L2 + omega*L3)/3
```

The four local products are then

```text
prefactor = i*chi_i/(3*sqrt(3))
AK2_BKp1 =  prefactor * (a21       - omega*conjg(b12))
BK2_AKp1 = -prefactor * (omega*a21 -       conjg(b12))
AK1_BKp2 =  prefactor * (a12       - omega*conjg(b21))
BK1_AKp2 = -prefactor * (omega*a12 -       conjg(b21))
```

The extra division by three converts the `3*f_bra^* f_ket` on the left of
Eqs. (76)-(79) to the physical product written in the local file. The
conjugated `b` always comes from the opposite layer direction; the first
formula uses `omega`, not its conjugate.

All `N=ndim/4` centers of each sublattice are evaluated. There is no one-in-three
Kekule mask, because retaining the exact phase permits use of every center.
Each center sublattice gives an estimate of the same envelope-product sum.
Count each honeycomb cell once by averaging the two grids:

```text
rho = 0.5 * (sum_A(local_product) + sum_B(local_product))
rho_per_cell = rho / (ndim/4)
```

Do not multiply these all-center sums by three. No extra spin, layer or valley
degeneracy factor is inserted. Results are separate for each stored spin
sector. Opposite matrix elements can be obtained by Hermitian conjugation.

## Output

Both filenames append the producer's full state suffix:

- `LayerOrderParameterLocal-numb<ndim>-nspin<s><parameters>` contains one row
  per layer-1 center. Columns are layer-1 site, paired layer-2 site, center
  sublattice, physical center `x,y`, then real/imaginary parts of the four
  local products in the order listed above.
- `LayerOrderParameter-numb<ndim>-nspin<s><parameters>` contains four labeled
  rows: component, `Re(rho)`, `Im(rho)`, `Re(rho_per_cell)`, `Im(rho_per_cell)`.

## Equal-layer comparison with Main_OrderParamsFull

The executable currently selects ordered layer pairs `(2,1)` and `(1,2)`;
it has no input setting for equal layers. Its internal loop and projection
routines can also be evaluated with both layer indices equal. This limit
corresponds to the **InterSubInterVal** output of `Main_OrderParamsFull`, not
the other three sublattice/valley channels that driver computes.

For equal layers and coincident centers the phase is
`chi = exp(2*i*K_layer.reference_r)`. The legacy output retains this valley
oscillation, whereas the new projection removes it. With the implemented
neighbor averaging, at the same A center for an arbitrary periodic density,

```text
new_AK_BKprime = chi * legacy_fKA_conjxfKpB
new_BK_AKprime = chi * legacy_fKB_conjxfKpA
```

There is no extra normalization factor. Raw output columns should therefore
not be compared before accounting for the position-dependent phase.

The spatial averaging follows the legacy routine. Let `C(r)` be the
oriented six-bond loop stored as `CurrentTemp` in `InterSubInterValAlt`, and
let `a1,a2` here denote Setup's primitive vectors rotated with the layer.
The legacy calculation uses

```text
D0 = C(r)
D1 = (C(r+a1-a2) + C(r+a2)    + C(r-a1))/3
D2 = (C(r+a1)    + C(r+a2-a1) + C(r-a2))/3
```

Feeding `L = [conjg(D1),conjg(D2),conjg(D0)]` into both loop arguments of
`ProjectLayerChannels`, with center sublattice `A`, gives the two equations
above **exactly**, for arbitrary loop values. Before spatial averaging, the
basic geometric three-loop construction at an A center is

```text
L = [conjg(C(r-a1)),conjg(C(r+a2-a1)),conjg(C(r))]
```

Those unaveraged selections coincide with the group averages only in the
constant-envelope limit. `LoopAverageShift` now replaces them with the complete
groups above, so production output matches the phase-adjusted legacy values
at every A center. The new driver also samples B centers using the spatially
inverted prescription and averages A/B cell sums; the legacy file lists A
centers from both layers and supplies no direct B-center comparison.

`tests/check_equal_layer_orderparameter.py` exercises the equal-layer helpers
against the actual legacy routine and saved-Fock driver, checking the production
averaged output at every A center in each layer, including periodic boundaries.
Constant intralayer intervalley coherence is not periodic in the
original moire cell, since K and Kprime have different boundary phases; its
constant-envelope check must exclude stencils crossing that boundary.
The general periodic-density comparison has no such exclusion. This check
uses unrelaxed geometry and the standard `RotateLayers=[-1,+1]` convention;
it does not establish equivalence of the different relaxed-site lookups.

## Validation and limits

`tests/check_layer_orderparameter.py` writes independent synthetic Fock files
in the producer's packed format. It uses two valid moire Bloch sectors with
unequal occupations and independent complex sublattice amplitudes, then
compares the driver output directly with their known envelope products.
It checks every center, both layer directions, one/two spin sectors, two twist
sizes, one/two translation shells, the separate A/B sums, total and per-cell
normalization, a common gauge rotation, zero input, and rejection of missing
translations. The driver is compiled with bounds/runtime checks for these
tests. `make -C postprocessing check` includes this regression.

`tests/check_layer_geometry.py` additionally checks the complete reference-site
lookup for Nam and Carr relaxation and reversed/equal layer rotations, with
Fortran runtime checks enabled and no saved Fock files present.

The local projection assumes envelopes vary slowly across a graphene loop
and uses co-rotated ideal-lattice bond phases. The algebra and synthetic
saved-state reconstruction are verified; production self-consistent states
still need physical convergence checks. In particular, agreement of A/B
estimates for a general state is controlled by that envelope approximation.
