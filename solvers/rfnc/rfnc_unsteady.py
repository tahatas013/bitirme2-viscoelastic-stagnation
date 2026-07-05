"""
RFNC -- Residual-Form Newton Collocation, UNSTEADY (Rothe) cozucu.

Rothe yontemi: tau'da backward Euler, her tau adiminda eta'da RFNC-BVP (Newton)
coz. f^n bir onceki adimdan gelir. Iç turning point cevrim boyunca dogup
kaybolur -- TEK yontemle, sinir kosulu ANAHTARLAMADAN.

eq.39 (A=1), tau-turevleri backward Euler ile:
   f_etatau   -> (f1 - f1_old)/dtau
   f_etaeeta..tau (f''')_tau -> (f3 - f3_old)/dtau
   a_tau      -> (a - a_old)/dtau

Unsteady residual (carpim formu -- f'ye HIC bolme yok):
   R = (f1-f1n)/dt + f1^2 - U f2 - (a-an)/dt - a^2 - f3
       - beta (f3-f3n)/dt - 2 beta f1 f3 + beta U f4 + beta f2^2

Sinir kosullari (eq.40-42, DEGISMEZ):
   f'(0)=0,  f(0)=k cos(omega tau),  f'(eta_max)=a(tau),  f''(eta_max)=0

Steady limiti (dt->inf, a sabit=1) -> rfnc_steady ile ayni.
"""
import numpy as np
from rfnc_steady import build_operators, solve_steady, bary_interp


def a_of_tau(tau, A_amp, omega):
    return 1.0 + A_amp * np.sin(omega * tau)


def wall_of_tau(tau, k_amp, omega):
    return k_amp * np.cos(omega * tau)


def residual_unsteady(U, ops, beta, dt, Uold, tau, k_amp, A_amp, omega):
    D1, D2, D3, D4 = ops['D1'], ops['D2'], ops['D3'], ops['D4']
    N = ops['N']
    f1, f2, f3, f4 = D1 @ U, D2 @ U, D3 @ U, D4 @ U
    f1n = D1 @ Uold
    f3n = D3 @ Uold
    a = a_of_tau(tau, A_amp, omega)
    an = a_of_tau(tau - dt, A_amp, omega)
    R = ((f1 - f1n) / dt + f1 ** 2 - U * f2 - (a - an) / dt - a ** 2 - f3
         - beta * (f3 - f3n) / dt - 2.0 * beta * f1 * f3 + beta * U * f4 + beta * f2 ** 2)
    # BC (eq.40-42)
    R[0] = U[0] - wall_of_tau(tau, k_amp, omega)
    R[1] = f1[0]
    R[N - 1] = f1[N] - a
    R[N] = f2[N]
    return R


def jacobian_unsteady(U, ops, beta, dt):
    """Steady Jacobian + backward-Euler katkilari (c1+=1/dt, c3+=-beta/dt)."""
    D1, D2, D3, D4 = ops['D1'], ops['D2'], ops['D3'], ops['D4']
    N = ops['N']
    f1, f2, f3, f4 = D1 @ U, D2 @ U, D3 @ U, D4 @ U
    coef_local = beta * f4 - f2
    cU = beta * U
    c1 = -2.0 * beta * f3 + 2.0 * f1 + 1.0 / dt          # +1/dt (f_etatau)
    c2 = 2.0 * beta * f2 - U
    c3 = -2.0 * beta * f1 - 1.0 - beta / dt              # -beta/dt (f3_tau)
    J = (np.diag(coef_local) + cU[:, None] * D4 + c1[:, None] * D1
         + c2[:, None] * D2 + c3[:, None] * D3)
    J[0, :] = 0.0; J[0, 0] = 1.0
    J[1, :] = D1[0, :]
    J[N - 1, :] = D1[N, :]
    J[N, :] = D2[N, :]
    return J


def newton_step(U0, ops, beta, dt, Uold, tau, k_amp, A_amp, omega,
                tol=1e-9, maxit=40):
    U = U0.copy()
    prev = np.inf
    for it in range(maxit):
        R = residual_unsteady(U, ops, beta, dt, Uold, tau, k_amp, A_amp, omega)
        nrm = np.max(np.abs(R))
        if nrm < tol:
            return U, True, nrm
        J = jacobian_unsteady(U, ops, beta, dt)
        try:
            dU = np.linalg.solve(J, -R)
        except np.linalg.LinAlgError:
            return U, False, nrm
        lam = 1.0
        for _ in range(30):
            Rn = np.max(np.abs(residual_unsteady(U + lam * dU, ops, beta, dt,
                        Uold, tau, k_amp, A_amp, omega)))
            if Rn < nrm or lam < 1e-6:
                break
            lam *= 0.5
        U = U + lam * dU
        if lam * np.max(np.abs(dU)) < 1e-12 or (nrm > 0.5 * prev and nrm < 1e-7):
            return U, True, nrm
        prev = nrm
    return U, (nrm < 1e-7), nrm


def fpp0_and_eta0(U, ops):
    eta = ops['eta']
    fpp0 = (ops['D2'] @ U)[0]
    eta0 = None
    for i in range(len(U) - 1):
        if eta[i] > 1e-9 and U[i] * U[i + 1] < 0:
            t = -U[i] / (U[i + 1] - U[i])
            eta0 = eta[i] + t * (eta[i + 1] - eta[i])
            break
    return fpp0, eta0


def run_rothe(beta=0.1, k_amp=-0.05, A_amp=0.03, omega=1.0, dt=0.05,
              tau_max=4 * np.pi, N=80, eta_max=10.0, verbose=False):
    """Rothe time-marching. Geri donus: dict(tau, fpp0, eta0, wall, Umin, snaps)."""
    ops = build_operators(N, eta_max)
    # IC: tau=0 steady cozumu (a=1, wall=k cos0=k_amp)
    U, ops_ic, _ = solve_steady(beta, k_amp, N=N, eta_max=eta_max)
    if U is None:
        return None
    ops = ops_ic
    out = dict(tau=[0.0], fpp0=[], eta0=[], wall=[wall_of_tau(0, k_amp, omega)],
               Umin=[U.min()], rms=[0.0], snaps={})
    f0, e0 = fpp0_and_eta0(U, ops)
    out['fpp0'].append(f0); out['eta0'].append(e0 if e0 else np.nan)
    snap_taus = [0.0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi, 3 * np.pi, 4 * np.pi]
    Uold = U.copy()
    tau = 0.0
    nstep = int(round(tau_max / dt))
    fail = None
    for n in range(1, nstep + 1):
        tau = n * dt
        U, ok, nrm = newton_step(U, ops, beta, dt, Uold, tau, k_amp, A_amp, omega)
        if not ok:
            fail = (tau, nrm)
            break
        Uold = U.copy()
        f0, e0 = fpp0_and_eta0(U, ops)
        out['tau'].append(tau); out['fpp0'].append(f0)
        out['eta0'].append(e0 if e0 else np.nan)
        out['wall'].append(wall_of_tau(tau, k_amp, omega))
        out['Umin'].append(U.min()); out['rms'].append(nrm)
        for st in snap_taus:
            if abs(tau - st) < dt / 2 and st not in out['snaps']:
                out['snaps'][st] = (ops['eta'].copy(), U.copy())
        if verbose and n % 40 == 0:
            print(f"  tau={tau:.3f}  f''0={f0:.5f}  eta0={e0}  ||R||={nrm:.1e}")
    out['ops'] = ops
    out['fail'] = fail
    for key in ('tau', 'fpp0', 'eta0', 'wall', 'Umin', 'rms'):
        out[key] = np.array(out[key])
    return out


if __name__ == "__main__":
    import warnings
    warnings.filterwarnings('ignore')
    print("Drift check (k=0, A=0):")
    r = run_rothe(k_amp=0.0, A_amp=0.0, tau_max=4 * np.pi, dt=0.05)
    drift = r['fpp0'].max() - r['fpp0'].min()
    print(f"  f''0 ilk={r['fpp0'][0]:.8f} son={r['fpp0'][-1]:.8f} drift={drift:.2e}")
    print("Tam salinim (k=-0.05, A=0.03):")
    r = run_rothe(k_amp=-0.05, A_amp=0.03, tau_max=4 * np.pi, dt=0.05, verbose=True)
    if r['fail']:
        print("  FAILED at", r['fail'])
    else:
        print(f"  f''0 araligi [{r['fpp0'].min():.4f}, {r['fpp0'].max():.4f}], "
              f"tamamlandi tau={r['tau'][-1]:.2f}")
