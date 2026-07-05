"""
Integral-reformulation ("integrate, don't differentiate") STEADY solver for the
second-grade viscoelastic stagnation BVP.  Coutsias-Hagstrom-Torres (1996) /
Greengard (1991) style.

Why this is the RIGHT method here.  The governing operator's highest-derivative
coefficient is  beta*f , which VANISHES at the interior turning point eta_0
(where f=0, present whenever c=f(0)<0).  In a DIFFERENTIAL formulation that
vanishing coefficient multiplies the worst-conditioned operator (D4), so the
discrete operator degenerates and its conditioning explodes with N (the
resolution<->conditioning squeeze).  Here we take the HIGHEST derivative as the
unknown,
        w = f''''  (Chebyshev-T coeffs),
and obtain f''', f'', f', f by INTEGRATING w four times (integration in
coefficient space is banded and norm-O(1), i.e. perfectly conditioned), with
four integration constants k3..k0 fixed by the four BCs.  Now the vanishing
coefficient beta*f multiplies the *identity* (the w term), while the
NON-vanishing coefficient of f''' (= -2 beta f' - 1 ~ -1 near eta_0) multiplies a
bounded integral operator -- so the assembled Jacobian stays well-conditioned
THROUGH the turning point.

LOCKED residual (eta-derivatives, A = a = 1):
   R = beta f f'''' - 2 beta f' f''' + beta (f'')^2 - f''' + (f')^2 - f f'' - A^2
LOCKED BCs:  f'(0)=0 , f(0)=c , f'(L)=A , f''(L)=0 .
Map eta = L(1-x)/2 ; integration in ETA uses Be = -(L/2) Bz (d eta = -(L/2) dx).

We integrate in ETA so the unknown w = f''''(eta) is O(1)-scaled (an x-derivative
unknown would be (L/2)^4 larger, an (L/2)^8 conditioning imbalance).
Unknowns u = [w_0..w_{n-1}, k3, k2, k1, k0]  (length n+4); F_m = f^(m)(eta) coeffs:
   F4=w ; F3=Be F4 + k3 e0 ; F2=Be F3 + k2 e0 ; F1=Be F2 + k1 e0 ; F0=Be F1 + k0 e0
Value maps  f^(m) = Q_m u  with Q_m = V P_m  (V = coeffs->values; F_m already eta).

NOTE beta=0 is degenerate for this method (beta f f'''' vanishes -> 3rd-order ODE,
w under-determined); use it for beta>0.  The beta=0 Hiemenz anchor is covered by
the ultraspherical / collocation solvers.
"""
import numpy as np
import ultraspherical_ops as U


def integ_matrix(n):
    """Zero-constant Chebyshev indefinite-integration operator Bz (T->T):
    columns are antiderivatives with the T_0 component dropped (the constants
    are carried separately by k3..k0).  d/dx(Bz a) == a up to that constant."""
    B = np.zeros((n, n))
    if n > 1:
        B[1, 0] = 1.0                      # int T_0 = T_1
    if n > 2:
        B[2, 1] = 0.25                     # int T_1 = (T_0+T_2)/4, drop T_0
    for k in range(2, n):
        if k + 1 < n:
            B[k + 1, k] += 1.0 / (2.0 * (k + 1))
        B[k - 1, k] += -1.0 / (2.0 * (k - 1))
    return B


class IntegralSteady:
    def __init__(self, n, L, A_out=1.0):
        self.n, self.L, self.A = n, L, A_out
        x = U.cheb_points(n)
        self.x = x
        self.eta = L * (1.0 - x) / 2.0
        # eta-integration operator: d eta = -(L/2) dx  =>  int .. d eta = -(L/2) Bz.
        # Unknown w = f''''(ETA) is O(1)-scaled (NOT (L/2)^4 larger), removing the
        # (L/2)^8 scaling imbalance that otherwise wrecks the conditioning.
        Be = -(L / 2.0) * integ_matrix(n)
        self.Be = Be
        e0 = np.zeros(n); e0[0] = 1.0
        # linear maps  u -> F_m = f^(m)(eta) Chebyshev coeffs,  u = [w | k3 k2 k1 k0]
        P4 = np.hstack([np.eye(n), np.zeros((n, 4))])
        P3 = np.zeros((n, n + 4)); P3[:, :n] = Be; P3[:, n + 0] = e0
        P2 = Be @ P3; P2[:, n + 1] += e0
        P1 = Be @ P2; P1[:, n + 2] += e0
        P0 = Be @ P1; P0[:, n + 3] += e0
        self.P = {0: P0, 1: P1, 2: P2, 3: P3, 4: P4}
        # value-from-coeffs matrix; F_m ARE the eta-derivatives (no extra scaling)
        Vmat = U.coeffs2vals(np.eye(n))
        self.Vmat = Vmat
        self.Q = {m: Vmat @ self.P[m] for m in range(5)}
        # BC functional rows (length n+4): endpoint evaluations of F_m
        ones = np.ones(n)                 # T_n(+1)  (x=+1 <-> eta=0)
        alt = (-1.0) ** np.arange(n)      # T_n(-1)  (x=-1 <-> eta=L)
        self.bc_rowfns = [
            ones @ P0,                    # f(0)   = F0(eta=0)   (-> c)
            ones @ P1,                    # f'(0)  = F1(eta=0)   (-> 0)
            alt @ P1,                     # f'(L)  = F1(eta=L)   (-> A)
            alt @ P2,                     # f''(L) = F2(eta=L)   (-> 0)
        ]

    def fields(self, u):
        return {m: self.Q[m] @ u for m in range(5)}

    def fpp0(self, u):
        """f''(0) = F2(eta=0) = F2 evaluated at x=+1 = sum of its T-coeffs."""
        return float(np.sum(self.P[2] @ u))

    def residual(self, u, beta, c):
        f = self.fields(u)
        f0, f1, f2, f3, f4 = f[0], f[1], f[2], f[3], f[4]
        Rv = (beta * f0 * f4 - 2 * beta * f1 * f3 + beta * f2 ** 2
              - f3 + f1 ** 2 - f0 * f2 - self.A ** 2)
        bc = np.array([
            self.bc_rowfns[0] @ u - c,
            self.bc_rowfns[1] @ u - 0.0,
            self.bc_rowfns[2] @ u - self.A,
            self.bc_rowfns[3] @ u - 0.0,
        ])
        return Rv, bc, f

    def jacobian(self, f, beta):
        f0, f1, f2, f3, f4 = f[0], f[1], f[2], f[3], f[4]
        c4 = beta * f0
        c3 = -2 * beta * f1 - 1.0
        c2 = 2 * beta * f2 - f0
        c1 = -2 * beta * f3 + 2 * f1
        c0 = beta * f4 - f2
        dR = (c0[:, None] * self.Q[0] + c1[:, None] * self.Q[1]
              + c2[:, None] * self.Q[2] + c3[:, None] * self.Q[3]
              + c4[:, None] * self.Q[4])
        J = np.vstack([dR, np.array(self.bc_rowfns)])
        return J

    def newton(self, u0, beta, c, tol=1e-11, maxit=40, verbose=False):
        u = u0.copy()
        for it in range(maxit):
            Rv, bc, f = self.residual(u, beta, c)
            res = np.concatenate([Rv, bc])
            rnorm = np.max(np.abs(res))
            if verbose:
                print(f"    it{it:2d}  ||R||={rnorm:.3e}")
            if rnorm < tol:
                return u, True, it, rnorm
            J = self.jacobian(f, beta)
            try:
                du = np.linalg.solve(J, -res)
            except np.linalg.LinAlgError:
                return u, False, it, rnorm
            lam = 1.0; best = None
            for _ in range(30):
                un = u + lam * du
                Rvn, bcn, _ = self.residual(un, beta, c)
                rn = np.max(np.abs(np.concatenate([Rvn, bcn])))
                if rn < rnorm:
                    best = (un, rn); break
                lam *= 0.5
            u = best[0] if best else u + lam * du
            if lam * np.max(np.abs(du)) < 1e-14:
                Rv, bc, f = self.residual(u, beta, c)
                rnorm = np.max(np.abs(np.concatenate([Rv, bc])))
                return u, rnorm < 3e-6, it, rnorm
        Rv, bc, f = self.residual(u, beta, c)
        rnorm = np.max(np.abs(np.concatenate([Rv, bc])))
        return u, rnorm < 3e-6, maxit, rnorm

    def init_from_f(self, fvals):
        """Build the exact unknown u = [w; k3..k0] representing the function whose
        eta-grid VALUES are fvals.  w = f''''(eta) coeffs = (2/L)^4 * chebder(aF,4);
        the 4 constants are recovered so that P0 u reproduces F (= f(eta)) exactly."""
        import numpy.polynomial.chebyshev as C
        n = self.n
        aF = U.vals2coeffs(fvals)                       # f(eta) as Chebyshev-x series
        wx = C.chebder(aF, 4)                           # F''''(x) coeffs
        wx = np.concatenate([wx, np.zeros(n - wx.size)])
        w = (2.0 / self.L) ** 4 * wx                    # f''''(eta) coeffs
        resid = aF - self.P[0][:, :n] @ w               # = aF - Be^4 w  (a cubic)
        k, *_ = np.linalg.lstsq(self.P[0][:, n:], resid, rcond=None)
        return np.concatenate([w, k])

    def hiemenz_init(self):
        # f = A(eta - 1 + e^{-eta}):  f(0)=0, f'(0)=0, f'(inf)=A
        return self.init_from_f(self.A * (self.eta - 1.0 + np.exp(-self.eta)))

    def cond(self, u, beta, c):
        _, _, f = self.residual(u, beta, c)
        return float(np.linalg.cond(self.jacobian(f, beta)))


def solve_steady(beta, c, n=64, L=10.0, A_out=1.0, n_beta=4, dc=0.02,
                 tol=1e-11, verbose=False):
    S = IntegralSteady(n, L, A_out)
    u = S.hiemenz_init()
    # beta-continuation: SKIP beta=0 -- there the beta f f'''' term vanishes and
    # the unknown w=f'''' is under-determined (the ODE drops to 3rd order).
    for b in np.linspace(0.0, beta, n_beta)[1:]:
        u, ok, it, nrm = S.newton(u, b, 0.0, tol=tol, verbose=verbose)
        if not ok:
            return S, u, False, f"beta-cont failed at beta={b:.4f} (||R||={nrm:.1e})"
    if abs(c) > 1e-15:
        nstep = max(1, int(round(abs(c) / dc)))
        for cc in np.linspace(0.0, c, nstep + 1)[1:]:
            u, ok, it, nrm = S.newton(u, beta, cc, tol=tol, verbose=verbose)
            if not ok:
                return S, u, False, f"c-cont failed at c={cc:.4f} (||R||={nrm:.1e})"
    return S, u, True, "ok"


# =============================================================== self-tests
if __name__ == "__main__":
    import numpy.polynomial.chebyshev as C
    np.set_printoptions(suppress=True)
    ok_all = True

    def check(name, cond, extra=""):
        global ok_all
        ok_all = ok_all and cond
        print(f"  [{'PASS' if cond else 'FAIL'}] {name:44s} {extra}")

    print("# integral-reform: integration operator + benchmarks + conditioning\n")
    n = 40
    Bz = integ_matrix(n)
    x = U.cheb_points(n)
    g = np.exp(0.6 * x) + np.cos(1.7 * x)
    ag = U.vals2coeffs(g)
    G = Bz @ ag                                  # antiderivative coeffs
    dG = C.chebder(G)
    dG = np.concatenate([dG, np.zeros(n - dG.size)])
    interr = np.max(np.abs(U.coeffs2vals(dG) - g))
    check("d/dx(Bz g) == g (integration op)", interr < 1e-11, f"err={interr:.2e}")

    print()
    # NB: beta=0 is DEGENERATE for this method (the beta f f'''' term vanishes -> the
    # ODE drops to 3rd order and w=f'''' is under-determined).  The beta=0 Hiemenz
    # benchmark is covered by the ultraspherical/collocation solvers; here we
    # validate the target regime beta=0.1.
    for label, beta, ref in [("beta=0.1 c=0", 0.1, 1.134114)]:
        print(f"  {label} (ref {ref}):")
        f2 = None
        for nn in (32, 48, 64, 96, 128):
            S, u, ok, info = solve_steady(beta, 0.0, n=nn, L=10.0)
            if not ok:
                print(f"    n={nn} FAILED {info}"); continue
            f2 = S.fpp0(u)
            print(f"    n={nn:3d}  f''(0)={f2:.8f}  err={abs(f2-ref):.1e}")
        check(f"{label}: matches ref", f2 is not None and abs(f2 - ref) < 1e-6,
              f"f''(0)={f2:.7f}" if f2 else "no convergence")

    print("\n  conditioning of integral-reform Jacobian (beta=0.1, L=15):")
    print("   c     n=64       n=128      n=256")
    for c in (-0.4, -0.8, -1.2, -2.0):
        row = []
        for nn in (64, 128, 256):
            S, u, ok, info = solve_steady(0.1, c, n=nn, L=15.0, dc=0.02)
            row.append(f"{S.cond(u,0.1,c):.1e}" if ok else "FAIL")
        print(f"  {c:+.1f}  {row[0]:>9}  {row[1]:>9}  {row[2]:>9}")

    print(f"\nOVERALL: {'PASS' if ok_all else 'FAIL'}")
