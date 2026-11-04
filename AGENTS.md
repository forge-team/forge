# FORGE programming style

FORGE is maintained by physicists. Favor code that a physicist can read,
check against an equation, and modify without learning a software framework.

- Follow the original procedural Fortran style: explicit declarations,
  allocations, do loops, if blocks, and direct array calculations.
- Keep physical formulas visible in the routine that calculates the quantity.
  Use the existing symbols and explain units, phases and normalization briefly.
- Use one statement per line. Avoid compressed chains of assignments.
- Use a small subroutine when it represents a recognizable calculation or
  repeated file operation. Avoid nested helpers that hide a loop's actions.
- Prefer direct calculations to caches and extra lookup arrays unless a
  measured production bottleneck requires them and the user requests it.
- Keep the existing Setup, Geometry, TightBinding, HartreeFock and Spectral
  source organization. Do not add parallel abstractions or general frameworks.
- Keep build commands straightforward, with compiler options visible.
- Preserve physical conventions, saved-state formats and output quantities.
  Retain useful input checks and verify numerical equivalence after refactoring.
- Do not optimize, consolidate or redesign unrelated code as part of an edit.
