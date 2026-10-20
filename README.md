# Computational Topological Condensed Matter

Numerical studies of topological band structure, built from scratch in NumPy.
Every model is validated against an analytically known result before any
figure is produced; the validation suites are part of the repository, not an
afterthought.

![Haldane phase diagram](figures/haldane_phase_diagram.png)

## Start here

| File | What it is |
|---|---|
| [`topological_floquet.ipynb`](topological_floquet.ipynb) | Self-contained Colab notebook — all three models, runs top to bottom with no installation, outputs and figures already embedded |
| [`report/report.pdf`](report/report.pdf) | 5-page write-up with the figures and methods |
| `src/` + `test_*.py` | The library and its validation suites (37 tests) |

## Contents

| Module | What it does | Validated against |
|---|---|---|
| `src/haldane.py` | Haldane model: Bloch Hamiltonian, Berry curvature, Chern number | analytic boundary \|M\| = 3√3·t₂·\|sin φ\| |
| `src/floquet.py` | Circularly driven graphene: Floquet operator, quasienergy bands, Floquet Chern number | reproduces the static Haldane model in the high-frequency limit |
| `src/anomalous.py` | Driven square lattice (Rudner–Lindner–Berg–Levin): ribbon spectra, edge localisation | real-space ring reproduces the independent 2×2 Bloch formula to 10⁻¹⁶ |

![Driven graphene](figures/floquet_driven_graphene.png)

![Anomalous Floquet phase](figures/anomalous_floquet.png)

## Method notes

**Chern number.** Computed by the Fukui–Hatsugai–Suzuki plaquette method:
link variables are built from overlaps of neighbouring Bloch states and the
field strength is the phase of the plaquette loop. This is manifestly
gauge-invariant, so no smooth gauge has to be chosen over the Brillouin zone —
which is the thing that cannot be done globally when the Chern number is
nonzero. The sum over plaquettes is an exact integer up to machine precision
even on coarse grids, which is why `test_chern_grid_independence` passes for
n = 40.

**Sign conventions.** The overall sign of *C* is fixed by the orientation
convention for the nearest- and next-nearest-neighbour vectors, not by
physics. With the ordering in `haldane.py`, φ = +π/2 gives C = +1. The tests
assert \|C\| = 1 and that the sign reverses under φ → −φ, which are the
convention-independent statements.

**Brillouin-zone sampling.** The Dirac points sit at fractional coordinates
(1/3, 2/3) and (2/3, 1/3), so a uniform n × n grid only samples them exactly
when 3 divides n. Checking that the gap closes at the phase boundary by
minimising over a grid silently fails for most n; `haldane.gap` should be
evaluated at K and K′ directly. `test_dirac_points_on_grid_when_divisible_by_three`
documents this.

**Which valley closes where.** The Haldane mass is +3√3·t₂·sin φ at K and
−3√3·t₂·sin φ at K′, so the gap at K closes at M = −M_c and at K′ at
M = +M_c. Each valley closes at exactly one of the two phase boundaries.

### Floquet module

**Why it exists.** Graphene irradiated with circularly polarised light acquires
an effective next-nearest-neighbour hopping with a complex phase — it *becomes*
a Haldane model. So the static model above is the high-frequency limit of a
driven one, and `test_matches_static_haldane_at_high_frequency` checks exactly
that: the quasienergy gap and Chern number of the driven system agree with the
static Haldane model carrying the induced t₂.

**The induced gap.** Measured, not quoted:

    Δ → (3/2) · t₁²A²/Ω   as A → 0

The 3/2 is convention-dependent (it reflects |δ| = 1/√3 in `haldane.DELTA` and
a dimensionless A), so it is recorded in `floquet.HIGHFREQ_COEFF` as a measured
constant with a test pinning it to six digits, rather than copied from a paper
with different conventions.

**Two expansions, not one.** The high-frequency result involves both a 1/Ω
(Magnus) expansion and a truncation at O(A²). Panel (c) separates them: curves
for Ω = 10 and Ω = 40 collapse onto each other, so the 1/Ω scaling is accurate
even at modest frequency, while the coefficient drifts from 3/2 with increasing
A. Plotting a single "relative error" conflates the two and makes the Magnus
expansion look like it is failing when it is not.

**Helicity is not the sign of A.** Reversing A → −A does *not* reverse circular
polarisation: −A(cos Ωt, sin Ωt) is the same circulation shifted by half a
period. Handedness lives in the relative sign of the two components, so it is
an explicit `helicity` parameter. Reversing it flips the Chern number and
leaves the gap unchanged, which is what
`test_chern_flips_with_helicity` and `test_helicity_does_not_change_gap` assert.

**Factor of two.** `gap` returns the full band splitting 2|d|, so the Haldane
Dirac gap is 6√3·t₂ — the familiar 3√3·t₂ is the *mass*, i.e. the half-gap.
`induced_t2` divides by 6√3 accordingly.

### Anomalous Floquet module

**The claim.** In panel (a) chiral edge modes traverse the whole quasienergy
zone, yet **both bulk bands have Chern number zero**. Panel (b) is a trivial
phase of the same model — also C = 0 for both bands, but no in-gap states
anywhere. The Chern number takes the same value in both phases, so in a driven
system it is simply not the invariant that counts the edge modes. The static
bulk–boundary correspondence fails here.

**What does count them** is the winding number of the full time-evolution
operator U(k, t) over the (k, t) torus — Rudner *et al.*, PRX 3, 031005 (2013).
Panel (a) shows its signature without computing it: each edge mode winds once
around the quasienergy circle from −π to +π as k_p crosses the Brillouin zone,
and that winding is what the invariant counts. Computing W directly is the
natural next step and is not done here.

**Why the loop closes.** At J·τ = π/2 each step is a perfect particle swap, so
a bulk particle traverses the four bonds of a plaquette and returns exactly to
its starting site: U(T) = 𝟙 for every k, which
`test_bulk_loop_closes_at_anomalous_point` checks to 10⁻¹². A particle on the
edge is missing one bond, cannot complete the loop, and is translated by one
cell per period instead. That asymmetry *is* the edge mode.

The cyclic ordering right → up → left → down is essential, not cosmetic:
`test_step_ordering_matters` scrambles it and asserts the loop no longer
closes.

**The fifth step is load-bearing.** With four hopping steps alone the two
quasienergy bands touch at ε = 0 for every J, leaving no gap to host an edge
mode and no well-defined Chern number. A staggered sublattice potential in a
fifth step opens it; `test_staggering_is_needed_to_open_the_zero_gap` asserts
both halves of that statement.

**Two gaps, not one.** A Floquet spectrum lives on a circle, so there is a gap
at ε = 0 *and* one at the zone edge ε = π. `bulk_gaps` returns both. Checking
only the ε = 0 gap would miss half the physics — in panel (a) the edge modes
cross both.

**Exactness.** In each step every site belongs to at most one bond, so the step
Hamiltonian is block diagonal in 2×2 blocks whose exponential is closed-form.
There is no Trotter error and no matrix exponential anywhere in the module;
ribbon unitarity holds to 10⁻³².

## Running

```bash
cd src
python test_haldane.py        # validation suites
python test_floquet.py
python test_anomalous.py
python fig_phase_diagram.py   # regenerates figures/
python fig_floquet.py
python fig_anomalous.py
```

Requires `numpy`, `matplotlib`, `scipy`.

## Validation

```
PASS  test_chern_grid_independence
PASS  test_chern_jumps_across_boundary
PASS  test_chern_sign_flips_with_phi
PASS  test_chern_topological_phase
PASS  test_chern_trivial_phase
PASS  test_curvature_integrates_to_chern
PASS  test_dirac_points_gapless_when_trivial
PASS  test_dirac_points_on_grid_when_divisible_by_three
PASS  test_gap_closes_at_phase_boundary
PASS  test_gap_matches_eigenvalues
PASS  test_hamiltonian_hermitian

PASS  test_chern_flips_with_helicity
PASS  test_chern_grid_independence
PASS  test_chern_quantised_and_opposite
PASS  test_floquet_operator_unitary
PASS  test_gap_scales_as_A_squared
PASS  test_gap_scales_as_inverse_omega
PASS  test_gap_scales_as_t1_squared
PASS  test_helicity_does_not_change_gap
PASS  test_highfreq_coefficient
PASS  test_matches_static_haldane_at_high_frequency
PASS  test_time_step_convergence_is_second_order
PASS  test_undriven_is_gapless_graphene

PASS  test_both_gaps_open_in_anomalous_phase
PASS  test_bulk_loop_closes_at_anomalous_point
PASS  test_bulk_unitary
PASS  test_chern_numbers_also_vanish_in_trivial_phase
PASS  test_chern_numbers_vanish_in_anomalous_phase
PASS  test_edge_modes_are_chiral
PASS  test_edge_state_count_independent_of_width
PASS  test_edge_states_exist_in_anomalous_phase
PASS  test_in_gap_states_are_edge_localised
PASS  test_no_edge_states_in_trivial_phase
PASS  test_ribbon_unitary
PASS  test_ring_reproduces_bulk
PASS  test_staggering_is_needed_to_open_the_zero_gap
PASS  test_step_ordering_matters
```

## References

- F. D. M. Haldane, *Model for a Quantum Hall Effect without Landau Levels*,
  Phys. Rev. Lett. **61**, 2015 (1988).
- T. Fukui, Y. Hatsugai, H. Suzuki, *Chern Numbers in Discretized Brillouin
  Zone*, J. Phys. Soc. Jpn. **74**, 1674 (2005).
- B. A. Bernevig, T. L. Hughes, *Topological Insulators and Topological
  Superconductors*, Princeton (2013), Ch. 8.
- T. Oka, H. Aoki, *Photovoltaic Hall effect in graphene*, Phys. Rev. B **79**,
  081406(R) (2009).
- M. S. Rudner, N. H. Lindner, E. Berg, M. Levin, *Anomalous Edge States and
  the Bulk-Edge Correspondence for Periodically Driven Two-Dimensional
  Systems*, Phys. Rev. X **3**, 031005 (2013).
