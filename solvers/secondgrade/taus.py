"""
taus.py -- tau_s estimators shared by BOTH solvers (identical extraction, so the
A-vs-B comparison is apples-to-apples).

Method (locked, B&H's own; see run_validation.py / B&H Fig 3.1):
  PRIMARY    G_m(tau) = min_eta f_eta.  Near tau_s, G_m -> -inf and 1/G_m is
             LINEAR in tau (B&H: G ~ (tau_s - tau)^{-1}).  Linear-fit the tail
             in the ACTIVE blow-up window (M = max|f_eta,eta| >= M_active, i.e.
             clear of the periodic plateau); tau_s = x-intercept.
  CROSS 1    eta_m^{-2} with eta_m = argmin_eta f_eta.  B&H: eta ~ (tau_s-tau)^{-1/2}
             so eta_m^{-2} ~ (tau_s - tau) is also linear; x-intercept = tau_s.
  CROSS 2    1/M^2 (M = max|f_eta,eta| ~ (tau_s-tau)^{-1/2}) also linear.
  NEVER      1/max|f_eta,eta| (first power) -- concave, underestimates.

All estimators take plain numpy arrays so Solver A and Solver B feed the SAME
function.
"""
import numpy as np


def lin_intercept(x, y):
    """Least-squares line y = m x + c; return (x-intercept = -c/m, R^2, npts)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 4:
        return None, None, len(x)
    A = np.vstack([x, np.ones_like(x)]).T
    m, c = np.linalg.lstsq(A, y, rcond=None)[0]
    if abs(m) < 1e-15:
        return None, None, len(x)
    r2 = 1.0 - np.sum((y - (m * x + c)) ** 2) / max(np.sum((y - y.mean()) ** 2), 1e-300)
    return -c / m, r2, len(x)


#  TAIL WINDOW.  B&H Fig 3.1 (p.1608) plots 1/G_m from its lower edge |G_m|~3 up
#  to the singularity (the curve reaches 0, i.e. |G_m|->inf -- there is NO upper
#  cutoff in B&H).  We fit the DEEP TAIL |G_m| >= G_lo, with G_lo=3 chosen to match
#  B&H's lower plotted edge.  The optional upper cap G_hi=33 (1/G_m~0.03) is OUR OWN
#  choice to drop the deepest, dt-corrupted points near blow-up; it does NOT change
#  tau_s (all windows give ~3.40, see run_3a).  Our runs actually reach |G_m|~90+.
#  1/G_m is robustly linear (R^2 > 0.999) across G_lo in [2, 5]; the cross-checks
#  (1/M^2, eta_m^{-2}) need the deepest tail and approach tau_s from below (they are
#  the most domain-sensitive -- B&H: the blow-up fills the whole domain).
MIN_PTS = 12
_GLO_LADDER = (4.0, 3.0, 2.0, 1.5, 1.0, 0.6)


def _tail_mask(gmin, G_lo):
    return np.isfinite(gmin) & (gmin <= -G_lo)


def _choose_glo(gmin, G_lo, G_hi=1e9, min_pts=MIN_PTS):
    """Use the REQUESTED G_lo if the window [G_lo, G_hi] has >= min_pts points;
    otherwise RELAX G_lo downward (shallower) until enough points are available.
    This keeps the fit window FIXED at the B&H value [3, 33] whenever the run
    reaches deep enough, and only widens the low end for shallow runs."""
    gmin = np.asarray(gmin, float)
    cand = [g for g in (G_lo,) + _GLO_LADDER if g <= G_lo]
    for g in cand:                      # from requested G_lo downward
        m = _tail_mask(gmin, g) & (gmin >= -G_hi)
        if m.sum() >= min_pts:
            return g
    for g in cand:                      # last resort: >= 4 points
        m = _tail_mask(gmin, g) & (gmin >= -G_hi)
        if m.sum() >= 4:
            return g
    return G_lo


def taus_primary(tau, gmin, M=None, G_lo=3.0, G_hi=33.0):
    """PRIMARY: 1/|G_m| linear fit over the deep tail |G_m| in [G_lo,G_hi] -> tau_s.
    G_lo=3 matches B&H's LOWER plotted edge in Fig 3.1; the upper cap G_hi=33 is our
    own (drops the deepest dt-corrupted points and does not affect tau_s).  G_lo
    falls back down the ladder if the run does not reach G_lo with >= MIN_PTS points.
    Returns (tau_s, R2, npts, G_lo_used)."""
    tau = np.asarray(tau, float); gmin = np.asarray(gmin, float)
    g = _choose_glo(gmin, G_lo, G_hi=G_hi)
    m = _tail_mask(gmin, g) & (gmin >= -G_hi)
    ts, r2, npt = lin_intercept(tau[m], 1.0 / (-gmin[m]))
    return ts, r2, npt, g


def taus_etam2(tau, eta_m, gmin, eta_max, G_lo=3.0):
    """CROSS 1: eta_m^{-2} linear fit over the deep tail, EXCLUDING points where
    eta_m has saturated against the finite domain (eta_m >= 0.85 eta_max)."""
    tau = np.asarray(tau, float); eta_m = np.asarray(eta_m, float)
    gmin = np.asarray(gmin, float)
    g = _choose_glo(gmin, G_lo, min_pts=5)
    m = _tail_mask(gmin, g) & np.isfinite(eta_m) & (eta_m > 1e-9) & (eta_m < 0.85 * eta_max)
    return lin_intercept(tau[m], 1.0 / eta_m[m] ** 2)


def taus_M2(tau, M, gmin, G_lo=3.0):
    """CROSS 2: 1/M^2 linear fit over the deep tail |G_m| >= G_lo -> tau_s."""
    tau = np.asarray(tau, float); M = np.asarray(M, float); gmin = np.asarray(gmin, float)
    g = _choose_glo(gmin, G_lo)
    m = _tail_mask(gmin, g) & np.isfinite(M)
    return lin_intercept(tau[m], 1.0 / M[m] ** 2)


def taus_band(tau, gmin, lo=2.0):
    """Window-sensitivity variant: 1/|G_m| fit on |G_m| >= lo only."""
    tau = np.asarray(tau, float); gmin = np.asarray(gmin, float)
    m = np.isfinite(gmin) & (-gmin >= lo)
    return lin_intercept(tau[m], 1.0 / (-gmin[m]))


#  ACTIVE-WINDOW estimators (Stage-2b primary).
#  The strongly-forced second-grade marches blow up in max|f''| (M > 1e4) while
#  |G_m| is still only ~15, so the B&H [3,33] window above holds too few, noisy
#  points (R^2 ~ 0.64).  The Stage-3a-validated ACTIVE window -- all reversal
#  points with M = max|f''| >= M_active (the periodic plateau is M ~ 0.5, so
#  M >= 1 is clear of it) and no upper cap -- is robust (R^2 > 0.96) and is the
#  SAME extraction fed to Solver A and Solver B (apples-to-apples).

def taus_active(tau, gmin, M, M_active=1.0):
    """PRIMARY (Stage 2b): 1/|G_m| linear fit over the ACTIVE blow-up window
    (gmin < -0.5 AND M >= M_active), no upper cap.  Returns (tau_s, R^2, npts)."""
    tau = np.asarray(tau, float); gmin = np.asarray(gmin, float)
    M = np.asarray(M, float)
    m = np.isfinite(gmin) & (gmin < -0.5) & np.isfinite(M) & (M >= M_active)
    if m.sum() < 4:
        m = np.isfinite(gmin) & (gmin < -0.5)
    return lin_intercept(tau[m], 1.0 / (-gmin[m]))


def taus_etam2_active(tau, eta_m, gmin, M, M_active=1.0):
    """CROSS (Stage 2b): eta_m^{-2} linear fit over the active window, excluding
    the final wall-collapse point eta_m -> 0 (filter eta_m > 1).  This is the
    B&H second estimator; it approaches tau_s from BELOW (eta_m saturates against
    the finite domain) so it forms the LOWER edge of the tau_s band."""
    tau = np.asarray(tau, float); eta_m = np.asarray(eta_m, float)
    gmin = np.asarray(gmin, float); M = np.asarray(M, float)
    m = (np.isfinite(gmin) & (gmin < -0.5) & np.isfinite(M) & (M >= M_active)
         & np.isfinite(eta_m) & (eta_m > 1.0))
    if m.sum() < 4:
        return None, None, int(m.sum())
    return lin_intercept(tau[m], eta_m[m] ** -2)


def gm_scaling_const_active(tau, gmin, M, tau_s, M_active=1.0):
    """B&H scaling |G_m|(tau_s - tau) ~ const over the active window."""
    tau = np.asarray(tau, float); gmin = np.asarray(gmin, float)
    M = np.asarray(M, float)
    m = (np.isfinite(gmin) & (gmin < -0.5) & np.isfinite(M) & (M >= M_active)
         & (tau < tau_s - 1e-3))
    if m.sum() < 3:
        return None, None
    prod = (-gmin[m]) * (tau_s - tau[m])
    return float(prod.mean()), float(prod.std())


def gm_scaling_const(tau, gmin, tau_s, G_lo=3.0):
    """B&H scaling check: |G_m| (tau_s - tau) ~ const (~0.58).  Mean +/- std,
    evaluated on the deep tail |G_m| >= G_lo."""
    tau = np.asarray(tau, float); gmin = np.asarray(gmin, float)
    g = _choose_glo(gmin, G_lo)
    m = _tail_mask(gmin, g) & (tau < tau_s - 1e-3)
    if m.sum() < 3:
        return None, None
    prod = (-gmin[m]) * (tau_s - tau[m])
    return float(prod.mean()), float(prod.std())


def richardson_band(eta_maxes, taus_vals):
    """First-order Richardson extrapolation tau_s(eta_max -> inf) by linear fit
    of tau_s vs 1/eta_max (intercept at 1/eta_max = 0).  Returns
    (tau_s_inf, slope, finite_band=(min,max))."""
    em = np.asarray(eta_maxes, float)
    ts = np.asarray(taus_vals, float)
    ok = np.isfinite(em) & np.isfinite(ts) & (em > 0)
    em, ts = em[ok], ts[ok]
    if len(em) < 2:
        return (float(ts[0]) if len(ts) else None), None, (None, None)
    x = 1.0 / em
    A = np.vstack([x, np.ones_like(x)]).T
    slope, inter = np.linalg.lstsq(A, ts, rcond=None)[0]
    return float(inter), float(slope), (float(ts.min()), float(ts.max()))
