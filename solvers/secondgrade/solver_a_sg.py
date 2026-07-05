"""
solver_a_sg.py -- Solver A (PRIMARY, full/exact beta): drives the VALIDATED
spectral RFNC step with beta != 0, k = 0.

We IMPORT (never modify) the validated driver
    viscoelastic_stagnation/bh_validation/bh_solver.py
which imports the validated spectral core
    OLDaFunc/residual_form_newton/rfnc_steady.py

Solver A = exactly the validated residual_unsteady / jacobian_unsteady /
newton_step in residual form -- the full eq. 39 with the beta f f'''' term and
NO division by f (so an interior f-zero is benign at the residual level), with
the Garg-Rajagopal 4th-order closure f''(inf)=0 (closure='4bc'), Chebyshev
collocation, analytic-Jacobian backward-Euler Rothe step, coefficient invA = 1/A.

This is the same machinery Stage 3a used for the beta-sweep.  A thin march loop
here records, at EVERY tau, the GLOBAL diagnostics for the B&H tau_s estimators
plus the conditioning / interior-zero diagnostics required by the certification:
    gmin  = min_eta f_eta         (G_m,  primary estimator)
    eta_m = argmin_eta f_eta      (eta_m^{-2} cross-check)
    M     = max_eta |f_eta,eta|   (1/M^2 cross-check + blow-up flag)
    fpp0  = f_eta,eta(0)          (wall curvature)
    eta0  = interior zero of f    (degeneracy site of the beta f'''' term)
    condJ = cond(jacobian_unsteady) (measured, never assumed -- see run_certify)
    rms   = Newton residual norm  (evidence the step converged despite cond(J))

Outer flow a_outer(tau, Delta, omega, 'bh') = 1 + Delta cos(omega tau) is the
NATIVE, validated B&H additive cosine; no runtime wrapper is needed.
"""
import numpy as np

import config_sg as C
import bh_solver as bh                       # validated (import-only)


def march(Delta, sigma, beta, eta_max, N, tau_max, dt=0.0025, omega=C.OMEGA,
          closure=C.CLOSURE, k_amp=C.K_AMP, M_blow=1.0e4, stop_on_blow=True,
          track_cond=False, verbose=False):
    """Backward-Euler Rothe march of the validated RFNC step (Solver A), beta != 0.

    sigma = omega/A (reported);  A = omega/sigma;  invA = 1/A = sigma/omega.
    The IC is the frozen tau=0 viscoelastic steady BVP (B&H convention) at the
    instantaneous outer flow a(0) = 1 + Delta, solved by the validated steady
    solver (beta-continuation built in).
    Returns dict with per-step arrays tau, a, gmin, eta_m, M, fpp0, eta0, rms,
    (condJ if track_cond) and the blow-up flag / tau_blow.
    """
    A = C.A_of_sigma(sigma, omega)
    invA = 1.0 / A
    a0 = bh.a_outer(0.0, Delta, omega, C.OUTER)
    w0 = bh.wall_of_tau(0.0, k_amp, omega)
    U, ops = bh.steady_ic(beta, w0, a0, N, eta_max)
    if U is None:
        return None
    eta = ops['eta']
    D1, D2 = ops['D1'], ops['D2']
    Nn = ops['N']

    def diag(U):
        f1 = D1 @ U
        f2 = D2 @ U
        if not np.all(np.isfinite(f1)):
            return np.nan, np.nan, np.inf, np.nan, np.nan
        i = int(np.argmin(f1))
        # interior zero of f (degeneracy site beta*f=0)
        eta0 = np.nan
        for j in range(len(U) - 1):
            if eta[j] > 1e-9 and U[j] * U[j + 1] < 0:
                t = -U[j] / (U[j + 1] - U[j])
                eta0 = eta[j] + t * (eta[j + 1] - eta[j])
                break
        return (float(f1[i]), float(eta[i]), float(np.max(np.abs(f2))),
                float(f2[0]), eta0)

    g0, em0, M0, fpp0_0, e00 = diag(U)
    rec = dict(Delta=Delta, sigma=sigma, beta=beta, A=A, invA=invA, omega=omega,
               eta_max=eta_max, N=N, dt=dt, closure=closure, k_amp=k_amp,
               tau=[0.0], a=[a0], gmin=[g0], eta_m=[em0], M=[M0], fpp0=[fpp0_0],
               eta0=[e00], rms=[0.0], singular=False, tau_blow=None, fails=0)
    if track_cond:
        rec['condJ'] = [float(np.linalg.cond(
            bh.jacobian_unsteady(U, ops, beta, dt, invA, closure)))]
    Uold = U.copy()
    nstep = int(round(tau_max / dt))
    for n in range(1, nstep + 1):
        tau = n * dt
        U, ok, nrm = bh.newton_step(U, ops, beta, dt, Uold, tau, k_amp, Delta,
                                    omega, invA, C.OUTER, closure)
        g, em, M, fpp0, e0 = diag(U)
        if not ok or not np.isfinite(nrm):
            rec['fails'] += 1
        blow = (not ok) or (not np.isfinite(M)) or (M > M_blow)
        if blow and stop_on_blow:
            rec['singular'] = True
            rec['tau_blow'] = tau
            for key, val in (('tau', tau),
                             ('a', bh.a_outer(tau, Delta, omega, C.OUTER)),
                             ('gmin', g), ('eta_m', em), ('M', M), ('fpp0', fpp0),
                             ('eta0', e0), ('rms', nrm)):
                rec[key].append(val if np.isfinite(val) else np.nan)
            if track_cond:
                rec['condJ'].append(np.nan)
            break
        Uold = U.copy()
        rec['tau'].append(tau)
        rec['a'].append(bh.a_outer(tau, Delta, omega, C.OUTER))
        rec['gmin'].append(g); rec['eta_m'].append(em); rec['M'].append(M)
        rec['fpp0'].append(fpp0); rec['eta0'].append(e0); rec['rms'].append(nrm)
        if track_cond:
            cj = float(np.linalg.cond(
                bh.jacobian_unsteady(U, ops, beta, dt, invA, closure)))
            rec['condJ'].append(cj)
        if verbose and n % 200 == 0:
            print(f"    [A b={beta}] tau={tau:.3f} gmin={g:+.3f} M={M:.3g} "
                  f"a={rec['a'][-1]:+.3f} ||R||={nrm:.1e}", flush=True)
    for k in ('tau', 'a', 'gmin', 'eta_m', 'M', 'fpp0', 'eta0', 'rms'):
        rec[k] = np.array(rec[k], dtype=float)
    if track_cond:
        rec['condJ'] = np.array(rec['condJ'], dtype=float)
    rec['eta'] = eta.copy()
    return rec


def snapshot_profiles(Delta, sigma, beta, eta_max, N, tau_query, dt=0.005,
                      omega=C.OMEGA, k_amp=C.K_AMP):
    """Return {tau: (eta, f, fp, fpp)} at the requested phases tau_query (for the
    thesis profile figures).  Uses the validated step; closest-step snapshots."""
    A = C.A_of_sigma(sigma, omega); invA = 1.0 / A
    a0 = bh.a_outer(0.0, Delta, omega, C.OUTER)
    U, ops = bh.steady_ic(beta, bh.wall_of_tau(0.0, k_amp, omega), a0, N, eta_max)
    if U is None:
        return None
    eta = ops['eta']; D1, D2 = ops['D1'], ops['D2']
    tau_query = sorted(tau_query)
    out = {}
    Uold = U.copy()
    tau_max = tau_query[-1] + dt
    qi = 0
    snaps = {}
    nstep = int(round(tau_max / dt))
    last_tau = 0.0
    for n in range(0, nstep + 1):
        tau = n * dt
        if n > 0:
            U, ok, nrm = bh.newton_step(U, ops, beta, dt, Uold, tau, k_amp,
                                        Delta, omega, invA, C.OUTER, C.CLOSURE)
            Uold = U.copy()
        while qi < len(tau_query) and tau >= tau_query[qi] - dt / 2:
            tq = tau_query[qi]
            snaps[tq] = (eta.copy(), U.copy(), (D1 @ U).copy(), (D2 @ U).copy())
            qi += 1
        if qi >= len(tau_query):
            break
    return snaps


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    import taus
    for beta in (0.0, 0.3):
        r = march(2.0, 0.5, beta, eta_max=30.0, N=240, tau_max=4.6, dt=0.0025)
        ts, R2, npt, glo = taus.taus_primary(r['tau'], r['gmin'])
        print(f"[A] beta={beta} Delta=2 sigma=0.5: singular={r['singular']} "
              f"tau_blow={r['tau_blow']:.4f} tau_s={ts:.4f} (R2={R2:.4f}, {npt}pts)")
