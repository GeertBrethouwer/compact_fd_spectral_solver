# Nonlinear-term forms and the wall-normal derivative on strongly stretched grids (α_s ≥ 4)

Question: on tanh grids with α_s = 4 or larger, is there anything to gain from changing hichan's compact
first derivative, for the rotational or divergence form of the momentum equations and the advective or
conservative form of the scalar equations? This is about accuracy and stability.

Candidates:
* **D1**, hichan's compact scheme with coefficients computed on the stretched grid. It is the most accurate,
  but it is not skew-adjoint.
* **D1m**, the same compact scheme on the uniform ξ grid, mapped with 1/y_ξ. It is skew-adjoint in the
  interior (to 1e-7) but less accurate at the wall.
* **Diagonal-norm SBP operators**, studied in 1D in the main README (`make sbp`).

All runs use the 6th-order scheme (fd_order = 6).

## What was run

[`forms_experiment.patch`](forms_experiment.patch) applies to `hichan/` with `patch -p1`, run inside `hichan/`.
It is for experiments only and is not proposed for hichan. It adds:
* `nl_form = 'divergence'` (∇·(uu)) and `'skew'` (½[u·∇u + ∇·(uu)]) for the momentum equations, built on
  hichan's existing `advect` and `add_flux_div`.
* `dy_op = 'd1' | 'd1m'` in `&physics`: the d/dy used in the momentum nonlinear terms (the vorticity in the
  rotational form). The pressure step and continuity always use D1, as in hichan.
* `scal_dy = 'd1' | 'd1m'` in `&scalars`: the d/dy of the scalar advection term. There is also
  `scal_form = 'skewd1'`, the skew-symmetric form with D1. hichan's `'skew'` always uses D1m.
* The volume mean of s² per scalar, as a third column per scalar in `log.dat`.

With `nl_form = 'rotational'` or `'advective'` and `dy_op = 'd1'` the patched code reproduces unpatched hichan
bit for bit (the `log.dat` columns that both write).

Three tests, each on tanh α_s = 4 and 5:

1. **Linear TS wave**
   * Re = 5000, α = 1.12, amplitude 1e-6, Ny = 33, 65, 129.
   * The decay rate is fitted over 50 < t ≤ 100, against Orr–Sommerfeld (ω_i = −2.78290e-3).
   * The fit window matters. Over t < 30 the fit still contains the other discrete modes that the interpolated
     Orr–Sommerfeld start field excites. For example, at α_s = 4 and Ny = 65 (advective form) the error is
     3.9e-5 over 1 < t ≤ 30 and a converged 1.2e-4 over t > 50. Time-step error is about 5e-6.
2. **Finite-amplitude TS wave**
   * Amplitude 0.05, nx = 16, Ny = 33 and 65, 0 ≤ t ≤ 50.
   * Compared with a reference at Ny = 385 and α_s = 2 (rotational form, D1). The reference itself is
     uncertain to about 1e-4.
3. **Turbulent channel**
   * Re_b = 2800 (Re_τ ≈ 200), 32×65×32, box 4π × 4π/3. dy⁺ at the wall / centre: 0.88 / 11.7 (α_s = 4) and
     0.41 / 14.2 (α_s = 5).
   * All runs restart from developed turbulence: rotational form with D1, run from noise to t = 150. The
     α_s = 5 runs start from that field interpolated in y, with a cubic spline.
   * Momentum runs: t = 40 at α_s = 4 and t = 30 at α_s = 5, averaged over t ≥ 10.
   * Scalar runs (flow: rotational form, D1) to t = 20, with two scalars:
     * s1: Pr = 0.71, s = ±1 at the walls;
     * s2: inviscid (Pr = 1e10), noise initially. Its variance ∫s² should stay constant.

Raw output is in [`results/`](results), and the case generator and analysis scripts are in [`scripts/`](scripts).

## Results

### Momentum: stability

| form | d/dy | TS a=0.05, α4 Ny=33 | α4 Ny=65 | α5 Ny=33 | α5 Ny=65 | channel α4 | channel α5 |
|---|---|---|---|---|---|---|---|
| rotational | D1 | stable | stable | blow-up t=14 | stable | stable | stable |
| rotational | D1m | blow-up t=8 | O(100) error | blow-up t=4 | blow-up t=18 | blow-up t=10 | blow-up t=4 |
| advective | D1 | blow-up t=11 | stable | blow-up t=8 | stable | stable | stable |
| advective | D1m | blow-up t=16 | stable | blow-up t=8 | stable | stable | blow-up t=10 |
| divergence | D1 | blow-up t=19 | stable | blow-up t=13 | stable | stable | stable |
| divergence | D1m | stable, O(1) error | blow-up t=3 | blow-up t=34 | blow-up t=3 | blow-up t=1.2 | blow-up t=1.2 |
| skew | D1 | stable | stable | blow-up t=22 | stable | stable | stable |
| skew | D1m | stable, O(0.5) error | blow-up t=6 | stable, O(1) error | blow-up t=5 | blow-up t=4 | blow-up t=2 |

### Momentum: accuracy

**Linear TS wave.** Signed relative error of ω_i:

| form | d/dy | α4 Ny=65 | α4 Ny=129 | α5 Ny=65 | α5 Ny=129 |
|---|---|---|---|---|---|
| rotational | D1 | −5.1e-4 | −1.6e-5 | −3.5e-4 | −9.8e-6 |
| rotational | D1m | −4.1e-4 | −1.6e-5 | 3.7e-3 | −7.3e-6 |
| advective | D1 | −1.3e-4 | −1.1e-5 | 1.0e-3 | 3.4e-6 |
| advective | D1m | −1.3e-4 | −1.1e-5 | 1.0e-3 | 3.3e-6 |
| divergence | D1 | −1.4e-5 | −9.5e-6 | 5.6e-6 | 9.0e-6 |
| divergence | D1m | 7.3e-5 | −2.5e-6 | −5.3e-3 | 1.1e-5 |
| skew | D1 | −9.5e-5 | −9.7e-6 | −2.4e-4 | 6.2e-6 |
| skew | D1m | 1.3e-4 | −6.2e-6 | 3.4e-3 | 7.3e-6 |

At Ny = 33 every variant has an O(0.1–10) error: too few points remain away from the walls.

**Finite-amplitude TS wave.** Maximum over t of the relative deviation of the TS-mode energy, Ny = 65:

| form | D1, α4 | D1m, α4 | D1, α5 | D1m, α5 |
|---|---|---|---|---|
| rotational | 2.2e-3 | 91 | 9.5e-3 | blow-up |
| advective | 1.3e-2 | 4.5e-3 | 5.6e-2 | 2.3e-2 |
| divergence | 4.5e-3 | blow-up | 1.8e-2 | blow-up |
| skew | 4.2e-3 | blow-up | 1.0e-2 | blow-up |

**Turbulent channel.** All stable runs agree within the scatter of these short averages:
- Re_τ 200.6–204.3;
- U_c⁺ 15.9–16.4;
- u_rms⁺ peak 2.57–2.63;
- −uv⁺ peak 0.75–0.80.

### Scalars

| scalar form | d/dy | s1 wall gradient, deviation from skew/D1m, α4 / α5 | s2 variance at t = 20 / t = 0, α4 | α5 |
|---|---|---|---|---|
| advective | D1 | 1.9e-3 / 1.1e-2 | 1.6e30 | 1.8e30 |
| advective | D1m | 3.9e-3 / 2.0e-2 | 154 | 162 |
| divergence | D1 | 6.8e-3 / 3.5e-2 | 4.6e18 | 2.3e28 |
| divergence | D1m | 2.5e-3 / 1.0e-2 | 335 | 644 |
| skew (`'skewd1'`) | D1 | 5.2e-3 / 2.4e-2 | 7.3e24 | 7.8e30 |
| skew (hichan `'skew'`) | D1m | 0 | **0.983** | **1.00** |

## Why: the pressure term and the operator's skew-adjointness

The quadrature weights `wq` are the energy norm. Two numbers decide whether d/dy can create energy:

| | uniform | α_s = 2 | α_s = 3 | α_s = 4 |
|---|---|---|---|---|
| D1, interior skew defect ‖W D + Dᵀ W‖/‖W D‖, Ny = 65 (rows 10…n−10) | 1e-7 | 0.44 | 0.90 | **1.4** |
| D1m, same, in its own weights h·y_ξ | 1e-7 | 1e-7 | 1e-7 | 1e-7 |

Near the walls neither D1 nor D1m is skew-adjoint in any diagonal norm: the defect is about 1 in both. The
diagonal-norm SBP operators have exactly zero defect.

What this means for each form:
* **Rotational form.** u·(u×ω) = 0 holds exactly at every collocation point, and the dealiased x, z products
  keep it exact. But the pressure term does work v·W·D1φ + (D1 v)·W·φ = vᵀ(W D1 + D1ᵀ W)φ. That is zero only
  for a skew-adjoint D1, so hichan's rotational form is robust but not discretely energy-conserving on
  stretched grids.
  * The robustness comes from the nonlinear term doing no work. At Ny = 65 the measured accuracy is slightly
    worse than with the divergence form.
  * Computing ω with D1m while the pressure step uses D1 breaks u·∇u = ω×u + ∇(u²/2) at the discrete level.
    The difference is an O(1) error on under-resolved scales, and those runs blow up.
* **Divergence and skew forms.** Here ∇·(uu) = u·∇u + u ∇·u needs the divergence in the nonlinear term to be
  the one that continuity sets to zero. With D1m in the nonlinear term and D1 in continuity, u ∇_m·u ≠ 0 acts
  as a source. This is the fastest blow-up of all (t ≈ 1–3).
* **Advective form.** d/dy acts on the velocity only through v ∂u/∂y, so a different d/dy does not destabilise
  it by the same mechanism.
  * At α_s = 4 the D1m variant is stable. It is as accurate as D1 in the linear TS test and in the channel
    statistics, and 3× more accurate in the finite-amplitude TS wave (4.5e-3 against 1.3e-2).
  * At α_s = 5 it blows up in the turbulent channel at t ≈ 10.
* **Scalars.** No pressure couples to the scalar. The skew-symmetric form conserves ∫s² exactly when d/dy is
  skew-adjoint. D1m is skew-adjoint in the interior, and the wall rows do not matter because the velocity
  vanishes there. So skew/D1m conserves the variance and nothing else does. The advective and divergence
  forms are not conservative with any d/dy, and a skew-adjoint operator alone does not fix them.

## Conclusions for α_s ≥ 4

* **Momentum, rotational or divergence form: no gain from changing the compact scheme.**
  * hichan's D1 is stable in all four forms on the turbulent channel at α_s = 4 and 5, and in the TS waves at
    Ny = 65.
  * Replacing it in the nonlinear terms alone (D1m, and by the same argument an SBP operator) makes the
    nonlinear term inconsistent with the D1 of the pressure step. That makes things much worse.
  * A real gain would require replacing D1 everywhere: nonlinear terms, pressure step and continuity, by one
    skew-adjoint operator. Only diagonal-norm SBP operators qualify; neither D1 nor D1m is skew-adjoint at
    the wall.
    * In 1D, SBP 4-8 is 10–180× less accurate than D1 at the wall on α_s = 3–4 (main README).
    * It does not fit the influence-matrix factorisation (D1 ± k) of the pressure step, which needs D1 = A1⁻¹B1.
  * Not worth it at DNS resolution, where D1 is stable.
* **Accuracy at α_s ≥ 4 is set more by the grid than by the scheme.**
  * At Ny = 65, the finite-amplitude TS wave is 2.4–4.3× less accurate on α_s = 5 than on α_s = 4 (with D1;
    centre spacing dy⁺ 14 against 12). In the linear TS test the advective form loses 8×, while the
    rotational and divergence forms barely change.
  * At Ny = 33, α_s ≥ 4 fails for every form.
  * The divergence form with D1 was the most accurate in the linear TS test (1e-5 at Ny = 65). Rotational and
    skew were the most accurate in the nonlinear one.
* **Scalars, advective or conservative form: no gain from D1m or SBP alone.**
  * For Pr = 0.71 all variants agree to 0.2–3.5%, the truncation level of this coarse grid.
  * For under-resolved scalars (high Pr, inviscid), only hichan's existing skew form with D1m is stable. The
    advective and divergence forms grow by 10¹⁸–10³⁰ with D1, and still ×150–650 with D1m.
  * Use `scal_form = 'skew'` whenever the scalar is not resolved.

## Reproducing

```bash
cd hichan && patch -p1 < ../dns_forms/forms_experiment.patch && make FC=mpif90
```

Then generate cases with `scripts/mkcase.sh` (a turbulent channel case; the TS cases follow
`hichan/examples/ts_wave`). Run them with `mpirun -np 2 path/to/hichan` in each case directory, and analyse with
`scripts/momstat.py`, `scripts/scan.py ALPHA` and `scripts/nlan.py`. On a 4-core laptop the full set (2 base
runs, 28 restarts, 48 linear and 35 nonlinear TS cases) took about 3–4 hours.
