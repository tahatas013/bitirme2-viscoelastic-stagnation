"""
config_sg.py -- shared configuration for STAGE 2b: the second-grade (beta != 0)
viscoelastic extension of the Blyth & Hall (2003) oscillatory stagnation-point
flow at k = 0 (impermeable wall, no wall transpiration).

This is the THESIS CONTRIBUTION.  There is NO peer-reviewed anchor for beta != 0,
so certification rests on THREE legs (see run_certify.py):
   (1) exact beta -> 0 reduction to the certified Stage-2a B&H numbers
       (newt_sol/newtonian_bh/: tau_s ~ 3.39-3.40, sigma_c ~ 1.12);
   (2) agreement of TWO INDEPENDENT solvers at beta != 0
       (full nonlinear residual  vs  conditioning-free beta-perturbation);
   (3) beta preserving the Banks-Zaturska blow-up exponent (G_m, eta_m scalings).

==============================================================================
GOVERNING PROBLEM (full eq. 39, beta > 0, k = 0, B&H additive-cosine outer):

  (1/A) f_etatau + f_eta^2 - f f_etaeta = (1/A) a_tau + a^2 + f_etaetaeta
        + beta[ (1/A) f_etaetaetatau + 2 f_eta f_etaetaeta - f f_etaetaetaeta
                - f_etaeta^2 ]
  a(tau) = 1 + Delta cos(omega tau),   a_tau = -Delta omega sin(omega tau)

FOUR boundary conditions (beta != 0 => 4th order):
  f_eta(0,tau)   = 0           (no-slip)
  f(0,tau)       = 0           (k = 0: impermeable, NO wall oscillation)
  f_eta(inf,tau) = a(tau)      (far-field outer-flow match; this is f_eta)
  f_etaeta(inf,tau) = 0        (Garg-Rajagopal closure: the ONLY BC added vs
                                Newtonian B&H; mode-analysis-forced for the
                                4th-order operator, not ad hoc)

==============================================================================
CLEAN MAPPING (identical gauge to Stage 2a, newt_sol/newtonian_bh/):
  omega = 1,  A = 1/sigma  =>  invA = 1/A = sigma  and  tau = t = tau_BH.
The time-derivative coefficient invA = 1/A EQUALS B&H's Strouhal sigma, so tau_s
is read DIRECTLY in B&H units.  General identity: sigma = omega/A (= 1/A only at
omega = 1).  Anchor case: Delta = 2, sigma = 0.5.

VALIDITY (small-Deborah): beta is the Deborah number; the second-grade model is
a SMALL-Deborah approximation.  beta <= 0.3 solid; beta = 0.5 CONTINUATION-FLAGGED;
beta > 0.5 deferred.  Solver B (perturbation) is additionally limited to its own
validity window M(c) = |beta f1''|/|f0''| << 0.15 (reported, never assumed).

MATH-ONLY framing for f(0,tau)=0: an impermeable wall condition; no suction /
blowing / injection language.  The finite-time singularity is a similarity-
solution feature that signals a physical boundary-layer transition / separation
(neglected terms regularize the true Navier-Stokes flow) -- not a real infinity.
==============================================================================
"""
import os
import sys
import numpy as np

# ---- paths --------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "figs")
DATA = os.path.join(HERE, "data")
os.makedirs(FIGS, exist_ok=True)
os.makedirs(DATA, exist_ok=True)

# validated solvers + core (IMPORT ONLY -- never modified)
BH_VALIDATION = os.path.normpath(os.path.join(
    HERE, "..", "..", "viscoelastic_stagnation", "bh_validation"))
CORE = os.path.normpath(os.path.join(
    HERE, "..", "..", "viscoelastic_stagnation", "OLDaFunc",
    "residual_form_newton"))
PERT = os.path.normpath(os.path.join(
    HERE, "..", "..", "viscoelastic_stagnation", "OLDaFunc",
    "perturbation_unsteady"))
NEWT_BH = os.path.normpath(os.path.join(HERE, "..", "newtonian_bh"))
# NOTE: NEWT_BH is referenced for the Stage-2a anchor numbers only; it is NOT
# put on sys.path, because it contains its own taus.py (without the Stage-2b
# active-window estimators) that would shadow our local taus.py.  HERE stays
# first on sys.path so our local modules win.
if HERE not in sys.path:
    sys.path.insert(0, HERE)
for p in (BH_VALIDATION, CORE, PERT):
    if p not in sys.path:
        sys.path.append(p)

# ---- mapping constants --------------------------------------------------
K_AMP = 0.0
OUTER = "bh"          # a = 1 + Delta cos(omega tau)  (native to bh_solver)
CLOSURE = "4bc"       # 4th-order closure f''(inf)=0 (required for beta != 0)
OMEGA = 1.0           # B&H gauge: omega = 1, A = 1/sigma  ->  invA = sigma

# Anchor singular case
DELTA_ANCHOR = 2.0
SIGMA_ANCHOR = 0.5

# Periodic case (thesis profile figures, two-solver agreement)
DELTA_PERIODIC = 0.5
SIGMA_PERIODIC = 0.5

# beta grid: solid <= 0.3; 0.5 continuation-flagged
BETAS = [0.0, 0.1, 0.2, 0.3]
BETA_FLAGGED = 0.5

# Solver-B (perturbation) validity threshold
MC_THRESHOLD = 0.15

# ---- certified Stage-2a anchors (newt_sol/newtonian_bh/) ----------------
TAU_S_BAND_2A = (3.38, 3.41)    # Stage-2a 1/|G_m| band (B&H 3.39)
TAU_S_BH = 3.39
SIGMA_C_BH = 1.12               # B&H Fig 3.3 singular-window edge, Delta=2
GM_SCALING_CONST = 0.58         # |G_m|(tau_s - tau) ~ const (B&H)

# steady-IC viscoelastic anchors at a = 1 (the IC the marches rest on)
STEADY_FPP0 = {0.0: 1.232588, 0.1: 1.134114, 0.3: 0.985, 0.5: 0.902500}


def A_of_sigma(sigma, omega=OMEGA):
    """sigma = omega/A  =>  A = omega/sigma  (at omega=1: A = 1/sigma)."""
    return omega / sigma


def invA_of(sigma, omega=OMEGA):
    """eq.39 time-derivative coefficient invA = 1/A = sigma/omega."""
    return sigma / omega


def period(omega=OMEGA):
    return 2.0 * np.pi / omega
