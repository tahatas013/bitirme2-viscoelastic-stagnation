# Stage 2b section 4c (Leg 3) -- exponent preservation: beta=0 vs 0.3 (+0.5)

Solver A, Delta=2, sigma=0.5, eta_max=30, N=240, dt=0.0025.  Active-window fits.  Exponent preserved <=> 1/|G_m| and eta_m^{-2} stay linear (high R^2) AND |G_m|(tau_s-tau) stays a bounded constant.

| beta | tau_s(1/\|G_m\|) | R^2(1/\|G_m\|) | R^2(eta_m^{-2}) | R^2(1/M^2) | \|G_m\|(t_s-t) | scatter/mean | verdict |
|---|---|---|---|---|---|---|---|
| 0.0 | 3.3971 | 0.984 | 0.997 | 0.533 | 0.548+/-0.060 | 0.11 | PRESERVED |
| 0.3 | 3.5335 | 0.965 | 0.990 | 0.892 | 0.681+/-0.172 | 0.25 | PRESERVED |
| 0.5 | 3.6254 | 0.927 | 0.986 | 0.781 | 0.779+/-0.465 | 0.60 | MARGINAL (cont.-flagged) |

**Section 4c verdict:** for beta <= 0.3 both B&H estimators stay LINEAR (R^2 > 0.95) and |G_m|(tau_s-tau) is a bounded constant (prefactor drifts ~0.55 -> ~0.8 but the POWER is unchanged) -> the Banks-Zaturska blow-up EXPONENT is PRESERVED; beta shifts the timing and prefactor, not the nature of the singularity.  beta=0.5 is MARGINAL (the scaling-constant scatter grows and the small-Deborah model is at its edge) -> continuation-flagged.
