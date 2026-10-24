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
