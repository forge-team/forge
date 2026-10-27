# Original supplied-source routine inventory against FORGE dev

This inventory describes the supplied files before adaptation. The working postprocessing module now imports dev screening; line numbers below refer to the originals.

Reference: `acbddd097c689142a38e5b0e99ac5fd674c0904a`. Case-insensitive procedure names; comments, whitespace and continuation markers are ignored for body comparison. Name absence does not exclude a renamed counterpart. Generic interfaces such as `diagonalize` and `matrix_elements` are not concrete procedures.

## Hamiltonian.f90

82 concrete procedures; 76 absent by name.

| Local procedure | Local line | Upstream match | Direct compute callers |
| --- | ---: | --- | --- |
| `fermi_dist` | 62 | Absent | compute_spectral_full.f90, compute_spectral.f90, compute_spectral_optimize.f90 |
| `longrange` | 91 | Absent | compute_drude_full.f90, compute_spectral_full.f90, compute_spectral.f90, compute_spectral_optimize.f90 |
| `wignerseitz` | 161 | Absent | compute_spectral_full.f90, compute_spectral.f90 |
| `wignerseitz6` | 299 | Absent | compute_drude_full.f90, compute_spectral_optimize.f90 |
| `sigma_func` | 598 | Absent | — |
| `temperature` | 610 | Absent | — |
| `complexh0` | 631 | Absent | — |
| `complexh0_phase` | 695 | Absent | — |
| `complexh0_mf` | 772 | Absent | — |
| `complexhk_mf_unitcell` | 816 | Absent | — |
| `complexhk_mf` | 922 | Absent | compute_drude_full.f90, compute_spectral_full.f90, compute_spectral.f90, compute_spectral_optimize.f90 |
| `vxy_mf` | 1031 | Absent | — |
| `vxy12i_mf` | 1133 | Absent | compute_drude_full.f90, compute_spectral_full.f90, compute_spectral.f90, compute_spectral_optimize.f90 |
| `vxy12_mf` | 1256 | Absent | — |
| `complexh0_hartree` | 1377 | Absent | — |
| `complexhk_hartree` | 1417 | Absent | — |
| `complexh0_sym` | 1461 | Absent | — |
| `complexh0_hartree_sym` | 1699 | Absent | — |
| `complexhk_hartree_sym` | 1870 | Absent | — |
| `fv_e` | 2072 | Absent | — |
| `fv` | 2130 | HartreeFock.f90:1246 (same name, different source) | — |
| `fv_a` | 2171 | Absent | — |
| `ftperp` | 2194 | Absent | — |
| `ft` | 2209 | Absent | — |
| `vx_new` | 2227 | Absent | — |
| `vxy12i` | 2299 | Absent | — |
| `vxy1_new` | 2361 | Absent | — |
| `vxyi_new` | 2410 | Absent | — |
| `vxy1d_new` | 2454 | Absent | — |
| `vxy12id_mf` | 2503 | Absent | compute_spectral_full.f90, compute_spectral.f90 |
| `vxy12id_no` | 2626 | Absent | compute_spectral_optimize.f90 |
| `vxy12d_mf` | 2714 | Absent | compute_spectral_optimize.f90 |
| `vxyid_new` | 2830 | Absent | — |
| `relaxation_old` | 2874 | Absent | — |
| `relaxation` | 3746 | Absent | compute_drude_full.f90, compute_spectral_full.f90, compute_spectral.f90, compute_spectral_optimize.f90 |
| `relaxation_2` | 4604 | Absent | — |
| `plotbands` | 5420 | TightBinding.f90:373 (same name, different source) | compute_drude_full.f90 |
| `plotbandsgamma` | 5505 | Absent | — |
| `plotbandsall` | 5586 | Absent | compute_drude_full.f90 |
| `plotbandsdensity` | 5715 | Absent | — |
| `plotbandsdensity3` | 5786 | Absent | — |
| `getdrude` | 5923 | Absent | — |
| `getdrudev` | 6075 | Absent | — |
| `getdrudefullvc` | 6257 | Absent | — |
| `getdrudefull` | 6456 | Absent | compute_drude_full.f90 |
| `dosfullvc` | 6632 | Absent | — |
| `dosfull` | 6771 | Absent | compute_drude_full.f90 |
| `dos_s` | 6874 | TightBinding.f90:620 (same name, different source) | — |
| `drude_s` | 6997 | Absent | — |
| `drude_v` | 7124 | Absent | — |
| `drude_f` | 7267 | Absent | compute_spectral_optimize.f90 |
| `drudecb` | 7400 | Absent | — |
| `drudevb` | 7569 | Absent | — |
| `drudefull` | 7738 | Absent | — |
| `dosvb_s` | 7910 | TightBinding.f90:739 (same name, different source) | — |
| `doscb_s` | 8031 | Absent | — |
| `dosjoint_s` | 8152 | TightBinding.f90:859 (same name, different source) | — |
| `complexh0_delta_magneticfield` | 8278 | Absent | — |
| `phase` | 8351 | Absent | — |
| `phase2` | 8392 | Absent | — |
| `ecoords` | 8425 | Absent | — |
| `sumek` | 8437 | Absent | — |
| `rho` | 8459 | Absent | — |
| `gcd` | 8471 | Absent | — |
| `vx_magneticfield` | 8505 | Absent | — |
| `sort` | 8560 | TightBinding.f90:1002 (same normalized source) | — |
| `sortf` | 8641 | Absent | — |
| `samplepoints_db` | 8749 | Absent | compute_spectral_full.f90, compute_spectral.f90 |
| `samplepoints1` | 8772 | Absent | compute_drude_full.f90, compute_spectral_full.f90, compute_spectral.f90 |
| `overlap` | 8798 | Absent | — |
| `overlapz` | 8842 | Absent | — |
| `drudet` | 8871 | Absent | — |
| `drudetemp` | 9173 | Absent | compute_drude_full.f90 |
| `conductivityfullx` | 9436 | Absent | compute_spectral_full.f90, compute_spectral.f90 |
| `conductivityfully` | 9597 | Absent | compute_spectral_full.f90, compute_spectral.f90 |
| `drudefullx` | 9757 | Absent | compute_spectral_full.f90, compute_spectral.f90 |
| `drudefully` | 9926 | Absent | compute_spectral_full.f90, compute_spectral.f90 |
| `dos_f` | 10095 | Absent | compute_spectral_full.f90, compute_spectral.f90, compute_spectral_optimize.f90 |
| `kramerskronig` | 10218 | Absent | compute_kramers_kronig.f90 |
| `conductivitys` | 10269 | Absent | compute_spectral_optimize.f90 |
| `plotbandsk` | 10428 | Absent | compute_spectral_full.f90, compute_spectral.f90 |
| `plotbandsdensity_db` | 10515 | Absent | compute_spectral_full.f90, compute_spectral.f90 |

## LapackRoutines.f90

26 concrete procedures; 2 absent by name.

| Local procedure | Local line | Upstream match | Direct compute callers |
| --- | ---: | --- | --- |
| `z_vector_overlap` | 80 | LapackRoutines.f90:77 (same name, different source) | — |
| `c_vector_overlap` | 105 | LapackRoutines.f90:100 (same name, different source) | — |
| `z_matrix_overlap` | 133 | LapackRoutines.f90:123 (same normalized source) | — |
| `c_matrix_overlap` | 160 | LapackRoutines.f90:147 (same normalized source) | — |
| `z_vector_mul_old` | 195 | Absent | — |
| `c_vector_mul_old` | 219 | Absent | — |
| `z_vector_mul` | 244 | LapackRoutines.f90:172 (same normalized source) | — |
| `c_vector_mul` | 267 | LapackRoutines.f90:191 (same normalized source) | — |
| `z_matrix_mul` | 294 | LapackRoutines.f90:212 (same name, different source) | — |
| `c_matrix_mul` | 308 | LapackRoutines.f90:228 (same name, different source) | — |
| `z_matrix_elements` | 368 | LapackRoutines.f90:245 (same name, different source) | — |
| `c_matrix_elements` | 411 | LapackRoutines.f90:275 (same name, different source) | — |
| `z_matrix_elements_diag` | 440 | LapackRoutines.f90:304 (same normalized source) | — |
| `c_matrix_elements_diag` | 459 | LapackRoutines.f90:323 (same normalized source) | — |
| `c_matrix_elements_diag_real` | 478 | LapackRoutines.f90:342 (same name, different source) | — |
| `z_matrix_elements_diag_real` | 497 | LapackRoutines.f90:362 (same name, different source) | compute_drude_full.f90, compute_spectral_full.f90, compute_spectral.f90, compute_spectral_optimize.f90 |
| `zdiagred` | 522 | LapackRoutines.f90:387 (same normalized source) | — |
| `cdiagred` | 583 | LapackRoutines.f90:448 (same normalized source) | — |
| `zdiagfull` | 639 | LapackRoutines.f90:504 (same normalized source) | — |
| `cdiagfull` | 683 | LapackRoutines.f90:548 (same normalized source) | — |
| `cdiag` | 730 | LapackRoutines.f90:595 (same normalized source) | — |
| `zdiag` | 786 | LapackRoutines.f90:651 (same normalized source) | — |
| `zdiagcomp` | 850 | LapackRoutines.f90:715 (same normalized source) | — |
| `cdiagcomp` | 896 | LapackRoutines.f90:761 (same name, different source) | — |
| `z_matrix_inverse` | 944 | LapackRoutines.f90:809 (same name, different source) | — |
| `c_matrix_inverse` | 982 | LapackRoutines.f90:841 (same name, different source) | — |
