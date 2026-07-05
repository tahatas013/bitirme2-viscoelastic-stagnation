# Stage 2b section 4b -- sigma_c(beta): the singular window narrows (Delta=2)

dt-refined classifier: N=120, dt=0.0025 (halved to 0.00125 on blow-up), 8 periods, eta_max=16.0.  S=singular, P=periodic, G=growing(slow-singular), ?=ambiguous.  sigma_c = midpoint of the last-{S,G}/first-P transition.

| beta | sig=0.8 | sig=0.9 | sig=1.0 | sig=1.1 | sig=1.2 | sig=1.3 | sigma_c |
|---|---|---|---|---|---|---|---|
| 0.0 | S | S | S | S | P | P | 1.15 |
| 0.1 | S | S | S | P | P | P | 1.05 |
| 0.3 | S | S | P | P | P | P | 0.95 |
| 0.5 | S | S | P | P | P | P | 0.95 |

**sigma_c(beta=0) ~ 1.15** (B&H ~1.12), sigma_c(0.1)~1.05, sigma_c(0.3)~0.95, sigma_c(0.5)~0.95.  The singular window (0, sigma_c) NARROWS (beta suppresses the singularity).  Consistent with section 4a: beta both DELAYS (tau_s up) and SUPPRESSES (window narrows) the singularity.

## Note on the independent cross-check at the edge

The conditioning-free Solver B (perturbation) CANNOT validate the sigma_c edge directly: it lives at the singular amplitude Delta=2, where the outer flow swings to a=3 and a=-1 and M(c) >> 0.15 (the perturbation is far outside its validity for ANY sigma at Delta=2 -- a measured M(c) ~ 10^3 confirms this). The edge is therefore a Solver-A (full nonlinear) result, corroborated by (i) the dt-refinement discriminator (genuine singular blow-ups persist as dt -> 0; numerical ones recede), and (ii) internal consistency with section 4a -- beta both DELAYS tau_s and NARROWS the window, two faces of the same suppression.  Solver B's clean two-solver agreement is established at the MODERATE amplitude Delta=0.5 (run_certify Leg 2a).

**Honest flag:** sigma_c carries +/-0.05 (grid) and the classification is dt-sensitive near the edge (coarse-dt spurious blow-ups excluded by the dt-halving check).  The NARROWING trend is robust; the absolute sigma_c at each beta is +/-0.05.  beta=0.5 continuation-flagged.