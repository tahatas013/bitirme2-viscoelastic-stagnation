"""
Well-conditioned spectral STEADY engine for the beta-perturbation expansion
f = f0 + beta f1  (valid for small beta; default beta=0.1, error O(beta^2)).

The expansion removes the turning point: order-0 and order-1 are NON-DEGENERATE
THIRD-order BVPs, so a single global Chebyshev INTEGRAL-reformulation (unknown
w = f''', integrate three times -> banded, ~N^2 conditioned) handles every c,
however negative.  No DD, no null mode, no regime switch.

Order 0 (unsteady Newtonian / Blyth-Hall), a = a(tau):
   f0''' + f0 f0'' - f0'^2 + a^2 + a_tau - f0_etatau = 0
   BCs: f0'(0)=0, f0(0)=c, f0'(L)=a
Order 1 (linear viscoelastic correction):
   L0[f1] - f1_etatau + B0 = 0
   L0[f1] = f1''' + f0 f1'' - 2 f0' f1' + f0'' f1
   B0     = f0_etaetaetatau + 2 f0' f0''' - f0 f0'''' - f0''^2
   BCs: f1'(0)=0, f1(0)=0, f1'(L)=0
(steady: drop the _etatau and _etaetaetatau time terms.)

Reuses ultraspherical_ops (vals2coeffs, coeffs2vals, cheb_points, integ_matrix)
UNMODIFIED.
"""
import os
import sys
import numpy as np
import numpy.polynomial.chebyshev as C

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "ultraspherical"))
import ultraspherical_ops as UO                              # noqa: E402
from integral_reform_steady import integ_matrix              # noqa: E402


def _cheb_diff_matrix(n):
    """Chebyshev differentiation matrix on T-coeffs (d/dx)."""
    D = np.zeros((n, n))
    for j in range(n):
        e = np.zeros(n); e[j] = 1.0
        d = C.chebder(e)
        D[:len(d), j] = d
    return D


class Pert3:
    """Shared 3rd-order integral-reform operators on [0,L] (eta-integration)."""

    def __init__(self, n, L):
        self.n, self.L = n, L
        x = UO.cheb_points(n)
        self.x = x
        self.eta = L * (1.0 - x) / 2.0
        Be = -(L / 2.0) * integ_matrix(n)            # eta-integration
        e0 = np.zeros(n); e0[0] = 1.0
        # u = [w (=f'''), k2, k1, k0] ; F_m = f^(m)(eta) coeffs
        P3 = np.hstack([np.eye(n), np.zeros((n, 3))])
        P2 = np.zeros((n, n + 3)); P2[:, :n] = Be; P2[:, n + 0] = e0
        P1 = Be @ P2; P1[:, n + 1] += e0
        P0 = Be @ P1; P0[:, n + 2] += e0
        self.P = {0: P0, 1: P1, 2: P2, 3: P3}
        Vmat = UO.coeffs2vals(np.eye(n))
        self.Vmat = Vmat
        self.Q = {m: Vmat @ self.P[m] for m in range(4)}
        # f'''' value map (differentiate w = f''' once): (-2/L) d/dx
        self.Q4 = (-2.0 / L) * (Vmat @ (_cheb_diff_matrix(n) @ P3[:, :n]))  # acts on w only
        ones = np.ones(n); alt = (-1.0) ** np.arange(n)
        self.ev_p1 = {m: ones @ self.P[m] for m in range(4)}   # eval at eta=0 (x=+1)
        self.ev_m1 = {m: alt @ self.P[m] for m in range(4)}    # eval at eta=L (x=-1)
        self.nu = n + 3

    def fields(self, u):
        f = {m: self.Q[m] @ u for m in range(4)}
        f[4] = self.Q4 @ u[:self.n]
        return f

    def fpp0(self, u):
        return float(self.ev_p1[2] @ u)              # f''(eta=0)


# ----------------------------------------------------------------- order 0
def order0_residual(E, u, a, atau, dt, f0p_prev):
    f = E.fields(u)
    Rv = (f[3] + f[0] * f[2] - f[1] ** 2 + a ** 2 + atau)
    if dt is not None:
        Rv = Rv - (1.0 / dt) * (f[1] - f0p_prev)
    bc = np.array([E.ev_p1[1] @ u,                   # f0'(0)=0
                   E.ev_p1[0] @ u - a if False else E.ev_p1[0] @ u,  # placeholder
                   E.ev_m1[1] @ u - a])              # f0'(L)=a
    return Rv, bc, f


def order0_solve(E, u0, beta_unused, c, a, atau=0.0, dt=None, f0p_prev=None,
                 source=None, tol=1e-11, maxit=60):
    """Newton solve of the (unsteady) order-0 problem.  c is the wall value."""
    n = E.n
    u = u0.copy()
    rn = np.inf
    for it in range(maxit):
        f = E.fields(u)
        Rv = f[3] + f[0] * f[2] - f[1] ** 2 + a ** 2 + atau
        if dt is not None:
            Rv = Rv - (1.0 / dt) * (f[1] - f0p_prev)
        if source is not None:
            Rv = Rv - source
        bc = np.array([E.ev_p1[1] @ u, E.ev_p1[0] @ u - c, E.ev_m1[1] @ u - a])
        G = np.concatenate([Rv, bc])
        rn = np.max(np.abs(G))
        if rn < tol:
            return u, True, rn
        # Jacobian: dR/du = Q3 + diag(f0)Q2 - 2 diag(f0')Q1 + diag(f0'')Q0 [-(1/dt)Q1]
        dR = (E.Q[3] + f[0][:, None] * E.Q[2] - 2 * f[1][:, None] * E.Q[1]
              + f[2][:, None] * E.Q[0])
        if dt is not None:
            dR = dR - (1.0 / dt) * E.Q[1]
        J = np.vstack([dR, E.ev_p1[1], E.ev_p1[0], E.ev_m1[1]])
        try:
            du = np.linalg.solve(J, -G)
        except np.linalg.LinAlgError:
            return u, False, rn
        lam, best = 1.0, None
        for _ in range(30):
            un = u + lam * du
            fn = E.fields(un)
            Rvn = fn[3] + fn[0] * fn[2] - fn[1] ** 2 + a ** 2 + atau
            if dt is not None:
                Rvn = Rvn - (1.0 / dt) * (fn[1] - f0p_prev)
            if source is not None:
                Rvn = Rvn - source
            bcn = np.array([E.ev_p1[1] @ un, E.ev_p1[0] @ un - c, E.ev_m1[1] @ un - a])
            if np.max(np.abs(np.concatenate([Rvn, bcn]))) < rn:
                best = lam; break
            lam *= 0.5
        u = u + (best if best else lam) * du
    f = E.fields(u)
    Rv = f[3] + f[0] * f[2] - f[1] ** 2 + a ** 2 + atau
    if dt is not None:
        Rv = Rv - (1.0 / dt) * (f[1] - f0p_prev)
    if source is not None:
        Rv = Rv - source
    bc = np.array([E.ev_p1[1] @ u, E.ev_p1[0] @ u - c, E.ev_m1[1] @ u - a])
    return u, np.max(np.abs(np.concatenate([Rv, bc]))) < 1e-8, rn


# ----------------------------------------------------------------- order 1
def order1_solve(E, f0, dt=None, f1p_prev=None, f0ppp_prev=None, source=None):
    """Linear solve for f1 (steady if dt is None).  f0 = order-0 fields dict."""
    n = E.n
    f0_0, f0_1, f0_2, f0_3, f0_4 = f0[0], f0[1], f0[2], f0[3], f0[4]
    # L0[f1] = f1''' + f0 f1'' - 2 f0' f1' + f0'' f1
    Lop = (E.Q[3] + f0_0[:, None] * E.Q[2] - 2 * f0_1[:, None] * E.Q[1]
           + f0_2[:, None] * E.Q[0])
    # B0 (steady) = 2 f0' f0''' - f0 f0'''' - f0''^2  ; unsteady adds f0_etaetaetatau
    B0 = 2 * f0_1 * f0_3 - f0_0 * f0_4 - f0_2 ** 2
    rhs = -B0
    if dt is not None:
        Lop = Lop - (1.0 / dt) * E.Q[1]
        B0 = B0 + (1.0 / dt) * (f0_3 - f0ppp_prev)        # f0_etaetaetatau (BE)
        rhs = -B0 - (1.0 / dt) * f1p_prev                  # known f1' term to RHS
    if source is not None:
        rhs = rhs + source
    J = np.vstack([Lop, E.ev_p1[1], E.ev_p1[0], E.ev_m1[1]])
    b = np.concatenate([rhs, [0.0, 0.0, 0.0]])     # f1'(0)=0,f1(0)=0,f1'(L)=0
    u1 = np.linalg.solve(J, b)
    return u1


def hiemenz_seed(E, c):
    return np.concatenate([UO.vals2coeffs(-np.exp(-E.eta)),
                           [1.0, 0.0, c]])  # rough


def steady_pert(c, a=1.0, beta=0.1, n=80, L=20.0):
    """Steady perturbation solve.  Returns (E, u0, u1, f0pp0, f1pp0, total)."""
    E = Pert3(n, L)
    # seed order-0 from Hiemenz-like
    fseed = a * (E.eta - 1.0 + np.exp(-E.eta)) + c * np.exp(-E.eta)
    # build u0 from f''' of seed: w = f'''(eta); k2,k1,k0 from matching f
    aF = UO.vals2coeffs(fseed)
    wx = C.chebder(aF, 3)
    w = (-2.0 / L) ** 3 * np.concatenate([wx, np.zeros(n - len(wx))])
    k, *_ = np.linalg.lstsq(E.P[0][:, n:], aF - E.P[0][:, :n] @ w, rcond=None)
    u0 = np.concatenate([w, k])
    u0, ok, rn = order0_solve(E, u0, beta, c, a)
    if not ok:
        return E, None, None, None, None, None
    f0 = E.fields(u0)
    u1 = order1_solve(E, f0)
    f0pp0 = E.fpp0(u0)
    f1pp0 = E.fpp0(u1)
    return E, u0, u1, f0pp0, f1pp0, f0pp0 + beta * f1pp0


if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")
    print("# steady beta-perturbation engine -- validation\n")
    print("| c | f0''(0) | f1''(0) | beta*f1''(0) | total f''(0) | ref (RFNC/DD) |")
    print("|---|---------|---------|--------------|--------------|---------------|")
    refs = {0.0: 1.134114, -0.4: 0.952173, -1.0: 0.7182575, -1.5: 0.5674172}
    for c in [0.0, -0.1, -0.4, -1.0, -1.5, -2.0, -3.0]:
        E, u0, u1, f0, f1, tot = steady_pert(c, a=1.0, beta=0.1, n=96, L=24.0)
        if u0 is None:
            print(f"| {c} | FAILED |"); continue
        r = refs.get(c, None)
        rs = f"{r}" if r else "--"
        print(f"| {c:+.1f} | {f0:.6f} | {f1:.4f} | {0.1*f1:+.5f} | {tot:.6f} | {rs} |")
    print("\n(c=0 order-0 should be Hiemenz 1.232588; total ~1.1340 to O(beta^2).)")
