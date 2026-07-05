"""
solver_b_sg.py -- Solver B (INDEPENDENT, CONDITIONING-FREE): the beta-perturbation
solver, driven with the B&H additive-cosine outer flow.

WHY THIS IS THE RIGHT INDEPENDENT CHECK.  Solver A solves ONE 4th-order problem
whose leading coefficient is beta*f; where f develops an interior zero (deep
reversal) that coefficient -> 0 and cond(J) can inflate (measured in run_certify:
1e13-1e15 in the deep singular tail).  The beta-perturbation expansion

      f = f0 + beta f1     (error O(beta^2), valid for small beta)

splits eq. 39 into TWO NON-DEGENERATE 3rd-order BVPs (order-0 = Newtonian
Blyth-Hall; order-1 = a LINEAR correction).  There is NO beta*f'''' term and
hence NO degeneracy -- the engine is conditioning-free for every reversal depth.
So if Solver A and Solver B AGREE in their common validity window, Solver A's
beta != 0 results are NOT conditioning-corrupted.  This is certification leg (2).

ENGINE (import-only, never modified): pert_steady.Pert3 + order0_solve +
order1_solve -- the validated, well-conditioned spectral INTEGRAL-reform engine
(w = f''', integrated three times; ~N^2 conditioned).  G1 (time-discretization
MMS) and G2 (steady limit) already PASS for that engine.

MAPPING THE ENGINE TO B&H eq. 39 WITHOUT EDITING IT.  The engine's native time
coefficient is 1 (it was written for sigma = 1).  eq. 39 has coefficient
invA = 1/A on BOTH time terms ( (1/A) f0_etatau and (1/A) f0_etaetaetatau ).
The engine discretizes every time term as (1/dt_engine)(.).  Passing

      dt_engine = dt / invA = dt * A        (so 1/dt_engine = invA/dt)
      atau_engine = invA * a_tau            (order-0 a_tau term)

makes order-0 (-1/dt_engine)(f0'-f0'_prev) -> invA*(f0'-f0'_prev)/dt and the
order-1 time terms likewise carry invA -- i.e. EXACTLY eq. 39's coefficient,
with the engine untouched (algebraically verified; the order-1 derivation is in
the module docstring of the certified pert_steady).  The wall value c = f(0,tau)
= k cos = 0 (k = 0).  The far-field uses a(tau) directly.

VALIDITY (reported, never assumed): M(c) = |beta f1''(0)| / |f0''(0)| << 0.15.
For the singular anchor Delta = 2 the outer flow swings to a = 3 at the IC, where
M(c) ~ 0.28 already exceeds 0.15: the perturbation is at/over its edge there, so
Solver B can only supply the LEADING-beta delay over the PRE-SINGULAR tail
(reversal started AND M(c) < 0.15 -- a real window, tau in ~[2.34, 3.18]).  For
the periodic Delta = 0.5 case M(c) stays small (<< 0.15) over the whole period
ONLY at beta <= 0.1 (M(c) ~ 0.14); at beta = 0.2, 0.3 it is 0.28, 0.42, past the
edge, so the clean two-solver agreement there is established at beta = 0.1.
"""
import numpy as np
import numpy.polynomial.chebyshev as Cheb

import config_sg as C
import pert_steady as PS                       # certified engine (import-only)


# ----------------------------------------------------------- steady IC (A-indep)
def steady_ic_bh(E, c, a0, beta):
    """Steady beta-perturbation IC on engine E at outer flow a0 (the time term
    vanishes at steady state, so this BVP is A-independent -- same construction
    as the certified pert_unsteady.steady_ic, re-derived locally so the validated
    original is untouched).  c-continuation from c=0 keeps the physical branch."""
    fseed = a0 * (E.eta - 1.0 + np.exp(-E.eta))
    aF = PS.UO.vals2coeffs(fseed)
    wx = Cheb.chebder(aF, 3)
    w = (-2.0 / E.L) ** 3 * np.concatenate([wx, np.zeros(E.n - len(wx))])
    kk, *_ = np.linalg.lstsq(E.P[0][:, E.n:], aF - E.P[0][:, :E.n] @ w, rcond=None)
    u0 = np.concatenate([w, kk])
    nstep = max(1, int(round(abs(c) / 0.1)))
    for cc in np.linspace(0.0, c, nstep + 1)[1:]:
        u0, ok, rn = PS.order0_solve(E, u0, beta, cc, a0)
        if not ok:
            return None, None
    u1 = PS.order1_solve(E, E.fields(u0))
    return u0, u1


def _mc(E, u0, u1, beta):
    """Validity metric M(c) = |beta f1''(0)| / |f0''(0)|."""
    return abs(beta * E.fpp0(u1)) / max(abs(E.fpp0(u0)), 1e-30)


# ----------------------------------------------------------- march (B&H outer)
def march(Delta, sigma, beta, tau_max, dt=0.0025, n=96, L=24.0, omega=C.OMEGA,
          k_amp=C.K_AMP):
    """Conditioning-free perturbation march with B&H additive-cosine outer flow.

    sigma = omega/A (reported); invA = 1/A; dt_engine = dt/invA (see header).
    Records per-step TOTAL fields f = f0 + beta f1 (for the shared tau_s
    estimators) AND order-0-only fields (Newtonian baseline / same-window
    comparison) AND M(c) validity.  k_amp != 0 supported (c = k cos(omega tau)).
    Returns dict; 'blow' set to tau if order-0 diverges (Newtonian singularity).
    """
    A = C.A_of_sigma(sigma, omega)
    invA = 1.0 / A
    dt_engine = dt / invA
    E = PS.Pert3(n, L)
    a0 = 1.0 + Delta * np.cos(omega * 0.0)
    c0 = k_amp * np.cos(omega * 0.0)
    u0, u1 = steady_ic_bh(E, c0, a0, beta)
    if u0 is None:
        return None
    f0 = E.fields(u0)
    f0p_prev = f0[1].copy(); f0ppp_prev = f0[3].copy()
    f1p_prev = (E.Q[1] @ u1).copy()
    eta = E.eta

    def diag(u0, u1):
        f0f = E.fields(u0)
        g0 = f0f[1]; gp0 = f0f[2]
        g1 = E.Q[1] @ u1; gp1 = E.Q[2] @ u1
        gt = g0 + beta * g1                       # total f'
        gpt = gp0 + beta * gp1                     # total f''
        it = int(np.argmin(gt))
        i0 = int(np.argmin(g0))
        return (float(gt.min()), float(eta[it]), float(np.max(np.abs(gpt))),
                float(g0.min()), float(eta[i0]), _mc(E, u0, u1, beta))

    gt0, em0, M0, g00, em00, mc0 = diag(u0, u1)
    rec = dict(Delta=Delta, sigma=sigma, beta=beta, A=A, invA=invA, omega=omega,
               eta_max=L, n=n, dt=dt, k_amp=k_amp,
               tau=[0.0], a=[a0], gmin=[gt0], eta_m=[em0], M=[M0],
               g0min=[g00], eta_m0=[em00], Mc=[mc0],
               fpp0=[E.fpp0(u0) + beta * E.fpp0(u1)],
               f0pp0=[E.fpp0(u0)], f1pp0=[E.fpp0(u1)], blow=None)
    nstep = int(round(tau_max / dt))
    for ns in range(1, nstep + 1):
        tau = ns * dt
        a_n = 1.0 + Delta * np.cos(omega * tau)
        atau_n = -(Delta * omega) * np.sin(omega * tau)
        c_n = k_amp * np.cos(omega * tau)
        u0, ok, rn = PS.order0_solve(E, u0, beta, c_n, a_n, invA * atau_n,
                                     dt=dt_engine, f0p_prev=f0p_prev)
        if not ok:
            rec['blow'] = tau
            break
        f0 = E.fields(u0)
        u1 = PS.order1_solve(E, f0, dt=dt_engine, f1p_prev=f1p_prev,
                             f0ppp_prev=f0ppp_prev)
        f0p_prev = f0[1].copy(); f0ppp_prev = f0[3].copy()
        f1p_prev = (E.Q[1] @ u1).copy()
        gt, em, M, g0m, em0_, mc = diag(u0, u1)
        rec['tau'].append(tau); rec['a'].append(a_n)
        rec['gmin'].append(gt); rec['eta_m'].append(em); rec['M'].append(M)
        rec['g0min'].append(g0m); rec['eta_m0'].append(em0_); rec['Mc'].append(mc)
        rec['fpp0'].append(E.fpp0(u0) + beta * E.fpp0(u1))
        rec['f0pp0'].append(E.fpp0(u0)); rec['f1pp0'].append(E.fpp0(u1))
    for k in ('tau', 'a', 'gmin', 'eta_m', 'M', 'g0min', 'eta_m0', 'Mc',
              'fpp0', 'f0pp0', 'f1pp0'):
        rec[k] = np.array(rec[k], float)
    rec['eta'] = eta.copy()
    return rec


def snapshot_profiles(Delta, sigma, beta, tau_query, dt=0.005, n=96, L=24.0,
                      omega=C.OMEGA):
    """Total f, f', f'' = (f0 + beta f1) and derivatives at requested phases, for
    the thesis figures (independent reconstruction of the velocity/shear)."""
    A = C.A_of_sigma(sigma, omega); invA = 1.0 / A; dt_engine = dt / invA
    E = PS.Pert3(n, L)
    u0, u1 = steady_ic_bh(E, 0.0, 1.0 + Delta, beta)
    if u0 is None:
        return None
    f0 = E.fields(u0); f0p_prev = f0[1].copy(); f0ppp_prev = f0[3].copy()
    f1p_prev = (E.Q[1] @ u1).copy()
    tau_query = sorted(tau_query); qi = 0; snaps = {}
    nstep = int(round((tau_query[-1] + dt) / dt))
    for ns in range(0, nstep + 1):
        tau = ns * dt
        if ns > 0:
            a_n = 1.0 + Delta * np.cos(omega * tau)
            atau_n = -(Delta * omega) * np.sin(omega * tau)
            u0, ok, rn = PS.order0_solve(E, u0, beta, 0.0, a_n, invA * atau_n,
                                         dt=dt_engine, f0p_prev=f0p_prev)
            f0 = E.fields(u0)
            u1 = PS.order1_solve(E, f0, dt=dt_engine, f1p_prev=f1p_prev,
                                 f0ppp_prev=f0ppp_prev)
            f0p_prev = f0[1].copy(); f0ppp_prev = f0[3].copy()
            f1p_prev = (E.Q[1] @ u1).copy()
        while qi < len(tau_query) and tau >= tau_query[qi] - dt / 2:
            tq = tau_query[qi]
            f = E.fields(u0)
            ff = f[0] + beta * (E.Q[0] @ u1)
            fp = f[1] + beta * (E.Q[1] @ u1)
            fpp = f[2] + beta * (E.Q[2] @ u1)
            snaps[tq] = (E.eta.copy(), ff, fp, fpp)
            qi += 1
        if qi >= len(tau_query):
            break
    return snaps


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    print("Solver B self-test (perturbation, B&H outer):")
    r = march(2.0, 0.5, 0.1, tau_max=3.30, dt=0.0025)
    win = (r['g0min'] < -0.5) & (r['Mc'] < C.MC_THRESHOLD)
    print(f"  Delta=2 beta=0.1: reached tau={r['tau'][-1]:.3f}, valid+reversal "
          f"window {win.sum()} pts, M(c) max in window={r['Mc'][win].max():.3f}")
    for bb in (0.1, 0.3):
        rp = march(0.5, 0.5, bb, tau_max=4 * 2 * np.pi, dt=0.01)
        mc = rp['Mc'].max()
        ok = "valid (M(c)<0.15)" if mc < C.MC_THRESHOLD else "PAST EDGE (M(c)>=0.15)"
        print(f"  periodic Delta=0.5 beta={bb}: M(c) max over march={mc:.4f} -> {ok}")
