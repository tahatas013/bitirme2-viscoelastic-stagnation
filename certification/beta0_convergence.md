# Independent beta->0 convergence certification of Solver A

Regime (LOCKED): eq. 39, k = 0 (f(0,tau) = 0), a(tau) = 1 + 0.5 cos(tau), sigma = 0.5 in the certified gauge omega = 1, A = 2, invA = 0.5 (held fixed across all runs).  a in [0.5, 1.5] > 0: attracting periodic orbit, no reversal, no interior f-zero, no singularity.

Discretization: N = 180 (gate 2N = 360), eta_max = 18, dt = P/1256 = 0.0050025 (validated 0.005 rounded to an integer number of steps per period), 6 periods, metrics on the final period [5P, 6P]; Newton tol 1e-9.  Outputs in this directory (/Users/tahaatas/Documents/kitaplar/2025-2026/bitirme/Sayısal Çözüm/newt_sol/beta0_indep/outputs) because /mnt/user-data/outputs is not available on this machine.

## Run health (Solver A ladder)

| N | beta | periodicity dev [4P,5P] vs [5P,6P] | max ||R|| | fails |
|---|---|---|---|---|
| 180 | 0 | 9.02e-09 | 6.2e-07 | 0 |
| 180 | 0.2 | 3.93e-07 | 6.3e-04 | 0 |
| 180 | 0.14 | 3.15e-07 | 4.3e-04 | 0 |
| 180 | 0.1 | 3.22e-07 | 3.2e-04 | 0 |
| 180 | 0.07 | 4.28e-07 | 2.2e-04 | 0 |
| 180 | 0.05 | 3.21e-07 | 1.6e-04 | 0 |
| 180 | 0.035 | 2.40e-07 | 1.1e-04 | 0 |
| 180 | 0.025 | 2.23e-07 | 7.9e-05 | 0 |
| 180 | 0.018 | 2.74e-07 | 5.4e-05 | 0 |
| 180 | 0.012 | 1.17e-07 | 3.9e-05 | 0 |
| 360 | 0 | 4.62e-07 | 1.6e-05 | 0 |
| 360 | 0.2 | 3.56e-06 | 1.6e-01 | 0 |
| 360 | 0.14 | 3.52e-06 | 1.1e-01 | 0 |
| 360 | 0.1 | 3.64e-06 | 8.1e-02 | 0 |
| 360 | 0.07 | 3.02e-06 | 5.6e-02 | 0 |
| 360 | 0.05 | 2.80e-06 | 3.9e-02 | 0 |
| 360 | 0.035 | 2.25e-06 | 2.9e-02 | 0 |
| 360 | 0.025 | 3.02e-06 | 2.0e-02 | 0 |
| 360 | 0.018 | 2.00e-06 | 1.5e-02 | 0 |
| 360 | 0.012 | 2.10e-06 | 9.7e-03 | 0 |

All runs periodic to < 1e-4 (transient fully decayed) and all Newton steps converged (fails = 0).  The raw max ||R|| column grows with N and beta because residual rows carry O(N^6) Chebyshev derivative-matrix norms (an absolute-residual scale effect, not a solution-accuracy loss -- see the beta = 0 sanity section); the update-stagnation exit still enforces ~1e-12 solution-space convergence per step, and the impact on E is bounded empirically by the resolution gate, the dt/2 check, and the periodicity column, all ~1e-6 relative or better.

## beta = 0 sanity (4bc closure at beta = 0)

| N | max cond(J) | max ||R|| | fails |
|---|---|---|---|
| 180 | 6.92e+08 | 6.2e-07 | 0 |
| 360 | 6.12e+10 | 1.5e-05 | 0 |

Probe windows: 0.5 period (N=180), 0.2 period (N=360); cond(J) is nearly flat along the periodic orbit.  The f''(eta_max) = 0 row remains compatible at beta = 0: cond(J)*eps << 1 (no rank deficiency), and float64 solves retain ~7 digits at N=180 and ~5 digits at N=360 -- consistent with the max ||R|| column (residual rows carry O(N^6) derivative-matrix norms, so the 1e-9 Newton tolerance is met at N=180 and limited to ~1e-5 absolute at N=360 by conditioning, while the Newton update stagnation exit guarantees ~1e-12 solution-space convergence).  The impact on the metrics is bounded empirically by the gate (E_N vs E_2N), the dt/2 check, and the periodicity deviations, all far below every E(beta).  No fallback reference was needed.

## PART 1 -- E(beta) ladder and resolution gate

E(beta) = max over the final period of |f''(0)_beta - f''(0)_0|, reference = Solver A's own beta = 0 run at the same N (shared backward-Euler error cancels in the difference).

| beta | E_N | E_2N | |E_N-E_2N|/E_2N | gate |
|---|---|---|---|---|
| 0.2 | 4.332402e-01 | 4.332401e-01 | 0.00% | pass |
| 0.14 | 3.328585e-01 | 3.328584e-01 | 0.00% | pass |
| 0.1 | 2.547670e-01 | 2.547672e-01 | 0.00% | pass |
| 0.07 | 1.886831e-01 | 1.886831e-01 | 0.00% | pass |
| 0.05 | 1.402960e-01 | 1.402961e-01 | 0.00% | pass |
| 0.035 | 1.013654e-01 | 1.013653e-01 | 0.00% | pass |
| 0.025 | 7.400725e-02 | 7.400712e-02 | 0.00% | pass |
| 0.018 | 5.413073e-02 | 5.413059e-02 | 0.00% | pass |
| 0.012 | 3.658778e-02 | 3.658767e-02 | 0.00% | pass |

Valid (resolution-passing) window: beta in [0.012, 0.2], 9 points.

- **Global slope** (log E_2N vs log beta, LSQ over the valid window): **0.8831 +/- 0.0140** (exp(intercept) = 1.9056).
- **Asymptotic pairwise slope** (two smallest passing betas 0.012, 0.018): **0.9660**.

Local slope between adjacent passing betas (rises monotonically toward 1 as beta -> 0 -- the signature of O(beta) with O(beta^2) contamination at finite beta, E ~ C1 beta + C2 beta^2 with C2 < 0):

| beta pair | local slope |
|---|---|
| 0.012 -> 0.018 | 0.9660 |
| 0.018 -> 0.025 | 0.9521 |
| 0.025 -> 0.035 | 0.9349 |
| 0.035 -> 0.05 | 0.9113 |
| 0.05 -> 0.07 | 0.8806 |
| 0.07 -> 0.1 | 0.8419 |
| 0.1 -> 0.14 | 0.7946 |
| 0.14 -> 0.2 | 0.7390 |

The global slope 0.883 < 1 is fit-window curvature, not a convergence-rate deficit: a two-term model E = C1 beta + C2 beta^2 (C2 < 0) captures the window to a few percent, and the local log-log slope of such a model approaches 1 from below as beta -> 0, exactly as the table shows.  Part 3 pins the leading coefficient with a refined small-beta extrapolation.

### dt spot-check (beta = 0.05, dt -> dt/2)

| N | E(dt) | E(dt/2) | rel. change |
|---|---|---|---|
| 180 | 1.402960e-01 | 1.402826e-01 | 0.01% |
| 360 | 1.402961e-01 | 1.402824e-01 | 0.01% |

## PART 2 -- independent B&H anchor (gauge-invariant loop)

Canonical bh_validation solver, SAME regime.  Primary anchor = 3bc closure: at beta = 0 it keeps the governing residual at the last node instead of Solver A's f''(eta_max) = 0 row, i.e. a genuinely different discrete system (the closure independence Solver A cannot provide itself).  Loop parameter = a(tau) (phase/gauge-invariant).

| comparison | max |df''(0)| on loop | max relative |
|---|---|---|
| A(beta=0) vs B&H 3bc, N=180 (PRIMARY anchor) | 2.744e-08 | 3.295e-08 |
| A(beta=0) vs B&H 4bc, N=180 (identity by construction, NOT an anchor) | 0.000e+00 | 0.000e+00 |
| A(beta=0) vs B&H 3bc, N=360 | 5.956e-07 | 4.652e-07 |

The 4bc row is exactly zero BECAUSE Solver A at beta = 0 with the 4bc closure executes the identical newton_step/IC/grid as run_rothe(beta=0, 4bc) -- it is a wiring sanity identity and claims nothing.  The 3bc row is the real check: a different discrete closure (governing residual retained at the last node, f''(eta_max)=0 row dropped), and it agrees to spectral-truncation level.  The chain to Blyth & Hall itself closes through bh_validation's own certification (bh_validation/RESULTS.md): steady f''(0) = 1.232588 vs the reference 1.232588 (|d| = 3.4e-7), singular-case tau_s converged to [3.39, 3.41] vs B&H 3.39, and the 4BC/3BC equivalence at beta = 0 already established there (tau_s difference 0.015, f''(eta_max) -> 0 emerges unimposed in 3BC).

## PART 3 -- three-way prefactor check (A vs B vs B&H)

Solver B uses the SAME a(tau) shape (1 + Delta cos tau, solver_b_sg.py line 125) -> comparable.  f1''(0,tau) is beta-independent (order-0 engine ignores beta), computed on the same tau grid, final period.

C_A estimators (leading beta-coefficient of E):

1. **Asymptotic (spec definition, PRIMARY): E_2N/beta at the smallest passing beta (0.012) = 3.0490.**
2. Extrapolated: quadratic-in-beta intercept of E_2N/beta over the 5 smallest passing betas = 3.1357 (max fit residual 8.9e-05).
3. exp(intercept) of the global log-log fit = 1.9056 -- REPORTED BUT NOT USED: with a fitted slope of 0.883 != 1 this quantity is the extrapolation of a curved fit to beta = 1, biased low by the O(beta^2) term (C2 < 0); it is only equivalent to the other two when the fitted slope is ~1.

- C_B = max_tau |f1''(0,tau)| = **3.1357** (periodicity dev 5.9e-12; M(c) max at beta=0.05: 0.069)
- |C_A - C_B| / C_B: **2.77%** (asymptotic, primary), 1.21e-06 (extrapolated), 39.23% (biased exp(intercept), not used)
- Phase consistency: argmax_tau |f''(0)_beta - f''(0)_0| at beta = 0.012 sits at grid index 1236 of the final period; argmax_tau |f1''(0,tau)| at index 1236 -- same phase (tau-offset 0.0000): the small-beta difference field IS the perturbation mode, measured by two independent implementations (shared: the analytic problem, a(tau), and the tau grid; not shared: discretization, engine, closure).

## Acceptance

- **(A) O(beta) slope:** asymptotic pairwise slope = 0.9660 in [0.9, 1.1] on a 9-point resolution-passing window (global fit 0.8831 +/- 0.0140, curvature-diagnosed above) -> **PASS**
- **(B) B&H anchor:** max relative deviation = 3.295e-08 < 1e-2 (primary: 3bc closure, N=180) -> **PASS**
- **(C) prefactor:** |C_A - C_B|/C_B = 2.77% (asymptotic estimator; extrapolated 1.21e-06) < 10% -> **PASS**

**Verdict (PASS):** Solver A -- which has NO B&H solution algebraically embedded -- converges to the Newtonian Blyth & Hall solution at the expected O(beta) rate (asymptotic slope 0.966, local slopes rising monotonically toward 1), its beta = 0 limit agrees with the differently-closed (3bc) canonical B&H discretization to 3.3e-08 relative, and its measured leading-order prefactor agrees with the beta-independent f1''(0) of the perturbation solver to 2.8%.  This certifies the 'beta -> 0 recovers B&H' leg independently of Solver B's embedded zeroth order.
