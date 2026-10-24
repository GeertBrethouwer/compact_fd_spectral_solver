# Accuracy tests of the compact finite differences in hichan

Standalone tests of the wall-normal compact finite-difference operators in
[hichan](https://github.com/GeertBrethouwer/hichan) (`src/compact.f90`): the (3-4-3) fourth-order and
(3-5-6-5-3) sixth-order schemes of Lenaers et al. (JCP 272, 2014), with coefficients computed directly on
a non-uniform grid. The tests compile `compact.f90` on its own (no MPI, no GPU) against a stub `mod_param`.

Tested against hichan commit `fa4f036`.

## Running

```bash
git clone https://github.com/GeertBrethouwer/hichan.git ../hichan
make run              # or: make run HICHAN=/path/to/hichan
```

| Test | What it measures |
|---|---|
| `test_global` | max and wall error of D1 = A1⁻¹B1 and D2 = A2⁻¹B2, metric derivative D1m, quadrature `wq`, 6-point wall derivative `ew1`; Ny = 33…257 |
| `test_interior` | error for \|y\| < 0.5, local truncation error of rows 3…ny−2, dense first row `c1` against D1 row 1; Ny = 33…513 |
| `test_wall` | wall-row coefficients of A1, A2, B2 and the D2 wall error up to Ny = 1025 |
| `test_wall_quad` | `test_wall` in quad precision (`rp` = real128), to separate round-off from truncation error |

Test function f = sin(2.3y + 0.4) on y ∈ [−1, 1], on a uniform grid and a tanh grid with α_s = 3.

## Results

No errors were found. Every operator converges at its design order:

| Quantity | 4th-order scheme | 6th-order scheme |
|---|---|---|
| D1, D2 interior (\|y\| < 0.5) | 4.00 | 6.00 (until round-off) |
| Local truncation error, rows 3…ny−2 | 4.0 | 6.0 |
| D1, D2 at the wall (closure) | 3.0 | 3.0 |
| Max error over the channel | 3 (set by the wall rows) | 3 (set by the wall rows) |
| Quadrature `wq` | 4.0 | 6.0 |
| Wall derivative `ew1` | 5 | 5 |
| D1m (max) | 3 | 3 |
| `c1` vs D1 row 1 | ~1e-13 | ~1e-13 |

On the uniform grid the coefficients are the textbook ones: interior LHS 1/4 and 1/3 (D1), 1/10 and 2/11
(D2); wall rows f′₁ + 2f′₂ = (−5f₁ + 4f₂ + f₃)/2h and f″₁ + 11f″₂ = (13f₁ − 27f₂ + 15f₃ − f₄)/h².

### Observations

1. **Round-off on the stretched grid.** The 6th-order D2 levels off at about 1e-10 in the interior
   (≈ ε/h_min²), and the D2 wall error levels off at about 3e-8 for Ny ≥ 513 (tanh, α_s = 3). The wall row of
   A2 has a large off-diagonal coefficient (11 uniform, ≈ 13.7 tanh) and the unpivoted Thomas
   factorisation then has a small negative pivot in row 2 (−0.1 to −0.25). In quad precision the 4th-order
   wall error returns to a clean 3rd order. Not relevant at usual resolutions (Ny = 65–385).
2. **6th-order D2 at the wall on the tanh grid** converges slowly (≈ 2.7–2.8 even in quad precision, still
   increasing at Ny = 1025) and is about 3× less accurate at the wall than the 4th-order D2. This looks like
   pre-asymptotic interaction between the 5th-order row 2 and the 3rd-order wall row, not an error.

### Not tested

* `kappa1` (1.74 and 2.0) is the maximum modified wavenumber of the uniform-grid interior scheme; the
  boundary rows on a stretched grid could raise the effective value used in the CFL limit.
* `spectral_radius_d2` uses power iteration with a 1.05 safety factor; it does not converge if the largest
  eigenvalue of the non-normal D2 is a complex pair.

## Compact D1 vs. summation-by-parts operators (`make sbp`)

`tests/compare_sbp.py` (needs numpy and scipy) compares four first-derivative operators on a uniform grid and
on hyperbolic-tangent grids with α_s = 2, 3, 4 (hichan's `grid_type = 'tanh'`), Ny = 33…513. Full output in
[`results_sbp.txt`](results_sbp.txt).

| Operator | What it is |
|---|---|
| `D1` | hichan's compact 4th-order scheme, coefficients computed on the stretched grid |
| `D1m` | hichan's compact scheme on the uniform ξ grid, `D1m = diag(1/y_ξ) D_ξ` |
| `CSBP` | compact Padé interior (P f′ = Q f, P = tridiag(1/6, 2/3, 1/6), Q = ±1/2) with a summation-by-parts closure: P symmetric positive definite, Q + Qᵀ = diag(−1, 0, …, 0, 1); mapped like `D1m` |
| `SBP24` | explicit diagonal-norm SBP operator (4th-order interior, 2nd-order boundary; Strand 1994), mapped like `D1m` |

`CSBP` is the collocated compact scheme with D = −Gᵀ up to boundary terms that was asked for. With the norm
equal to the compact left-hand side P, closures of 4 to 8 rows with P and Q up to pentadiagonal and
nonadiagonal are only consistent with **first-order** boundary rows (2nd- and 3rd-order rows give an
inconsistent linear system). The closure used minimises the degree-2…4 boundary error with P positive definite
(smallest eigenvalue 0.2). It is exactly SBP on the uniform grid; after the mapping it is SBP only up to the
commutator of P and diag(y_ξ), because P is not diagonal.

Results at Ny = 129 (orders from Ny = 257 → 513):

| | `D1` | `D1m` | `CSBP` | `SBP24` |
|---|---|---|---|---|
| max error order (wall) | 3 (4 for α_s = 4) | 3 | **1** | 2 |
| max error, tanh α_s = 3, Ny = 129 | 5.5e-7 | 5.7e-5 | 4.1e-3 | 1.4e-3 |
| interior (\|y\| < 0.5) order | 4 | 4 | 4 | 4 (6× larger constant) |
| SBP defect, uniform | 1.1 | 1.1 | 4e-14 | 2e-16 |
| SBP defect, tanh α_s = 2 / 3 / 4 | 1.0–1.1 | 1.1 | 0.013 / 0.024 / 0.034 | 2e-16 |
| best diagonal norm defect | 1.6–1.8 | 1.8 | 0.44 | 2e-15 |
| inflow advection, max Re λ·Δy_min, tanh α_s = 2 / 3 / 4 | 3e-3 / 2e-3 / −2e-3 | 0.044 / 0.072 / 0.094 | 0.028 / 0.042 / 0.053 | ≤ 3e-16 |
| inflow advection, transient growth max‖e^{tL}‖, α_s = 3 | 24 | 1.6e26 | 2.3e15 | 2.1 |
| skew-symmetric advection, max Re λ·Δy_min, α_s = 2 / 3 / 4 | **0.010 / 0.011 / 0.008** | ≈ 0 | ≈ 0 | ≈ 0 |
| skew-symmetric advection, growth to t = 20, α_s = 3 | **1e20** | 4.6 | 3.8 | 3.9 |
| spectral radius ·Δy_min | 1.30–1.73 | 1.30–1.73 | 1.32–1.73 | 1.02–1.37 |

SBP defect: ‖H D + Dᵀ H − E‖₂ / ‖H D‖₂ in the operator's own norm (`D1`, `D1m`: quadrature `wq`;
`CSBP`: h Y^½ P Y^½; `SBP24`: h H diag(y_ξ)). Inflow advection: u_t + u_y = 0 with u = 0 imposed at the inflow
wall (both directions). Skew-symmetric advection: u_t = −½[v D u + D(v u)] with v = (1 − y²)(½ + y), which is
zero at the walls as in the channel and conserves energy exactly in the continuous case. Growth values of a few
for the neutral cases are the Euclidean norm seen through a non-uniform energy norm, not instability.

### Conclusions

1. **hichan's choice is confirmed.** In skew-symmetric form (the situation in the channel, v = 0 at the walls)
   the direct compact `D1` is unstable on every stretched grid (growth rate ≈ 0.01/Δy_min), while `D1m` is
   neutral to round-off. This is why hichan's skew-symmetric scalar advection uses `D1m`.
2. **An exactly SBP compact 4th-order scheme costs a full order at the wall.** With the norm equal to the
   compact LHS only first-order boundary rows are possible, so `CSBP` is 800–7000× less accurate at the wall than
   `D1`. On stretched grids it is not exactly SBP either (non-diagonal norm), and with strongly imposed inflow
   conditions it is as unstable as `D1m`. It gains nothing over `D1m` in skew-symmetric form.
3. **Only the diagonal-norm operator is provably energy-stable on a stretched grid.** `SBP24` is SBP to round-off
   on every grid, stable for both tests, and has a ≈ 20 % smaller spectral radius (larger explicit time step),
   but is 2nd order at the wall and has a 6× larger interior error constant than the compact schemes.
4. The inflow/outflow test is harsher than the channel (no wall-normal flux at a no-slip wall); the instability
   of `D1m`/`CSBP` there comes from strongly imposed boundary values with a non-diagonal norm and would need
   SAT (penalty) boundary conditions in an inflow/outflow problem.

## Higher-order SBP and mimetic operators

### Diagonal-norm SBP 3-6 and 4-8 (`make sbp`)

`tests/sbp_diag.py` derives diagonal-norm SBP first-derivative operators (Strand 1994; Mattsson & Nordström 2004)
of interior order 2s and boundary order s: the accuracy conditions are linear in the norm weights and the
boundary block of Q, so each closure is one linear solve; the free parameters (0, 0 and 3) are chosen to minimise
the boundary truncation error and the spectral radius. The norms reproduce Strand's exactly (2-4: 17/48, 59/48,
43/48, 49/48; 3-6: 13649/43200, 12013/8640, …; 4-8: 1498139/5080320, 1107307/725760, …). Mapped like `D1m`
with H_y = H diag(y_ξ), they are exactly SBP on any monotone grid. Both are now part of `compare_sbp.py`, together
with hichan's 6th-order `D1-6` and `D1m-6`.

Max error of f′ at Ny = 513 (wall-dominated) / interior error at Ny = 257, and skew-symmetric advection stability:

| grid | `D1-4` | `D1-6` | `SBP36` | `SBP48` |
|---|---|---|---|---|
| uniform | 2.8e-7 / 1.3e-9 | 2.2e-7 / 6e-14 | 8.8e-7 / 6e-13 | **1.2e-8** / 4e-14 |
| tanh α_s = 2 | 4.8e-8 / 3.8e-9 | 3.8e-8 / 2e-13 | 1.5e-6 / 5e-11 | 3.0e-8 / 6e-14 |
| tanh α_s = 3 | 7.8e-9 / 9.6e-9 | 6.2e-9 / 7e-13 | 2.7e-6 / 3e-10 | 8.2e-8 / 6e-13 |
| tanh α_s = 4 | 1.5e-9 / 2.4e-8 | 7.4e-10 / 3e-12 | 2.7e-6 / 1.4e-9 | 1.3e-7 / 5e-12 |
| wall order | 3 | 3 | 3 | **4** |
| SBP defect (tanh) | ≈ 1 | ≈ 1 | 2e-16 | 2e-16 |
| skew-symmetric advection (tanh) | unstable (growth up to 1e29) | unstable (up to 1e44) | neutral | neutral |
| spectral radius ·Δy_min (tanh) | 1.30–1.48 | 1.44–1.67 | 1.07–1.30 | 1.09–1.38 |

* **SBP 4-8 is the highest-order collocated operator here that is provably stable on a stretched grid.** It is
  4th order at the wall and 25× more accurate than hichan's compact schemes on a uniform grid; on the mildly
  stretched α_s = 2 grid it is about as accurate as `D1-6` (3.0e-8 vs 3.8e-8). For α_s = 3 and 4 hichan's direct
  `D1` (coefficients fitted to the actual grid) has a 10–180× smaller wall error, because the mapped SBP operators
  feel the strong stretching near the wall through y_ξ.
* SBP 3-6 has only a 3rd-order boundary and a large boundary error constant: it is the least accurate of the
  four on stretched grids.
* `D1-6` has the same problem as `D1-4` in skew-symmetric form (growth rate ≈ 0.015/Δy_min); `D1m-6` is neutral.

### Corbino & Castillo mimetic operators (`make mimetic`)

The Corbino–Castillo operators (MOLE, [csrc-sdsu/mole](https://github.com/csrc-sdsu/mole)) are **staggered**:
G maps cell centres plus the two wall values to faces, D maps faces to centres, orders k = 4, 6, 8. They are not
collocated SBP operators and would not fit hichan's collocated influence-matrix solver without a staggered
redesign. `tests/compare_mimetic.py` tests what they are built for (accuracy of G and of the Laplacian D G, the
extended Gauss divergence theorem Q D + Gᵀ P = B with MOLE's weights, and the spectrum of the Dirichlet
Laplacian) against hichan's compact D2. MOLE is GPL-3.0, so its coefficient tables are read from MOLE's
`gradNonPeriodic.m` / `divNonPeriodic.m` at run time; no MOLE code is in this repository. Full output in
[`results_mimetic.txt`](results_mimetic.txt).

| tanh α_s = 3 | mimetic k=4 | k=6 | k=8 | hichan D2-4 | hichan D2-6 |
|---|---|---|---|---|---|
| G error at Ny = 257 (order) | 8.3e-8 (4) | 6.6e-11 (6) | 9.7e-13 (6.7) | – | – |
| Laplacian error at Ny = 33 / 257 | 1.9e-2 / 3.1e-5 | 2.0e-3 / 2.3e-8 | 6.5e-4 / **7.4e-10** | 1.6e-4 / 6.0e-8 | 3.9e-5 / 1.0e-7 |
| Laplacian order | 3 | 5 | 4.6 | 3.4 | 2.7 |
| B ≠ plain boundary term, rows from wall | 11 | 14 | 16 | – | – |
| MOLE weights P, Q positive | yes | yes | **no** | – | – |
| Dirichlet Laplacian: Re λ < 0 | yes | yes | yes | yes | yes |
| max \|Im λ\| / \|λ\| | 0 | 0.04 (0.13 uniform) | 0.29 | 0 | 0 |
| spectral radius ·Δy_min² | 5.9 | 5.1 | 5.1 | 4.2 | 4.7 |

* **Accuracy at Ny = 257:** the k = 8 Laplacian is more accurate than hichan's best compact D2 on every grid
  (5·10⁴× uniform, 3000× α_s = 2, 80× α_s = 3, 9× α_s = 4); k = 6 is better on the uniform, α_s = 2 and 3 grids
  (10⁴×, 40×, 2.6×) but 12× worse for α_s = 4. hichan's 6th-order D2 converges at only 2.5–2.7 at the wall on
  tanh grids. On coarse, strongly stretched grids (Ny = 33, α_s ≥ 3) hichan's D2 is 17–1000× more accurate.
* **The "extended" Gauss theorem is not SBP in the usual sense:** the boundary operator B = Q D + Gᵀ P differs
  from the plain boundary term (−1, +1 at the walls) in the first 11–16 rows from each wall, with deviations up
  to 0.8, 1.5 and 4.5 for k = 4, 6, 8. Energy estimates therefore contain boundary-layer terms, and for k = 8 MOLE's weights are not all
  positive, so they do not define a norm.
* **Spectrum:** all Dirichlet Laplacians are stable (Re λ < 0), but for k = 6 and 8 the mimetic Laplacian is not
  self-adjoint: eigenvalues with imaginary parts up to 30 % of the spectral radius. hichan's compact D2 has a real
  spectrum.
