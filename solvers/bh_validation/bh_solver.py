"""
bh_solver.py -- Blyth & Hall (2003) dis akisli zamana-bagli RFNC marsi.

TASARIM ILKESI (gorev kisitlari):
  * Dogrulanmis cekirdek `residual_form_newton/` HIC DEGISMEZ; sadece IMPORT ile
    yeniden kullanilir (spektral operatorler + steady IC cozucu).
  * Dis akis a(tau) SECILEBILIR (outer='bh' | 'belhaj_add' | 'belhaj_recip');
    eski form korunur. Bu task: outer='bh'.
  * Tam (beta, k) makinesi korunur; bu task'ta beta=0, k=0 ama kod ileride
    beta!=0, k!=0 calistirir.
  * Strouhal: sigma = omega / A  (ASLA sigma=1/A; o ozdeslik yalniz omega=1'de).
    eq.39'da zaman-turevi terimlerinin katsayisi (1/A)'dir; sigma=omega/A
    RAPORLANAN buyukluktur. invA = 1/A = sigma/omega.
  * BC kapanisi secilebilir: closure='4bc' (f''(inf)=0 dahil) | '3bc' (beta=0'da
    3. mertebe; f''(inf)=0 DUSURULUR -> dogrulama (ii)).

eq.39 (tam), backward Euler (Rothe), katsayi invA=1/A:
  R = invA*(f1-f1n)/dt + f1^2 - U f2 - invA*(a-an)/dt - a^2 - f3
      - beta*invA*(f3-f3n)/dt - 2 beta f1 f3 + beta U f4 + beta f2^2
BC: f'(0)=0; f(0)=k cos(omega tau); f'(inf)=a(tau); [4bc] f''(inf)=0.
a_tau backward-Euler farkindan (a-an)/dt gelir; ayrica kodlanmaz (f'ye bolme yok).
"""
import os
import sys
import numpy as np

# --- dogrulanmis cekirdegi IMPORT ile yeniden kullan (DEGISTIRME) ---
_CORE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "..", "OLDaFunc", "residual_form_newton")
_CORE = os.path.normpath(_CORE)
if _CORE not in sys.path:
    sys.path.insert(0, _CORE)
from rfnc_steady import build_operators, solve_steady, bary_interp, diagnostics  # noqa: E402


# ============================================================ DIS AKIS (secilebilir)
def a_outer(tau, amp, omega, outer='bh'):
    """Dis akis a(tau). amp = Delta (bh) veya A_amp (belhaj)."""
    if outer == 'bh':                       # Blyth & Hall (2003) eq.(2.1)
        return 1.0 + amp * np.cos(omega * tau)
    if outer == 'belhaj_add':               # eski additive (cekirdekteki form)
        return 1.0 + amp * np.sin(omega * tau)
    if outer == 'belhaj_recip':             # eski Belhaj reciprocal
        return 1.0 / (1.0 + amp * np.sin(omega * tau))
    raise ValueError(f"bilinmeyen outer='{outer}'")


def wall_of_tau(tau, k_amp, omega):
    """Duvar (eq.41): f(0,tau) = k cos(omega tau). k=0 bu task."""
    return k_amp * np.cos(omega * tau)


def run_header(beta, k_amp, amp, omega, A, outer, closure, dt, N, eta_max, tau0):
    sigma = omega / A
    invA = 1.0 / A
    return (f"[bh_solver] outer={outer} closure={closure} | "
            f"sigma=omega/A={sigma:.6g}  omega={omega:.6g}  A={A:.6g}  invA=1/A={invA:.6g}  "
            f"Delta(amp)={amp:.6g}  k={k_amp:.6g}  beta={beta:.6g} | "
            f"dt={dt:.6g} N={N} eta_max={eta_max:.6g} tau0={tau0:.6g}")


# ================================================================== RESIDUAL & JACOBIAN
def residual_unsteady(U, ops, beta, dt, Uold, tau, k_amp, amp, omega, invA,
                      outer='bh', closure='4bc'):
    D1, D2, D3, D4 = ops['D1'], ops['D2'], ops['D3'], ops['D4']
    N = ops['N']
    f1, f2, f3, f4 = D1 @ U, D2 @ U, D3 @ U, D4 @ U
    f1n, f3n = D1 @ Uold, D3 @ Uold
    a = a_outer(tau, amp, omega, outer)
    an = a_outer(tau - dt, amp, omega, outer)
    # eq.39 artik (katsayi invA=1/A; sigma=omega*invA raporlanir)
    R = (invA * (f1 - f1n) / dt + f1 ** 2 - U * f2 - invA * (a - an) / dt - a ** 2 - f3
         - beta * invA * (f3 - f3n) / dt - 2.0 * beta * f1 * f3 + beta * U * f4 + beta * f2 ** 2)
    R[0] = U[0] - wall_of_tau(tau, k_amp, omega)   # f(0)=k cos(omega tau)
    R[1] = f1[0]                                   # f'(0)=0
    R[N - 1] = f1[N] - a                           # f'(inf)=a(tau)  (anlik a; a<0 olabilir)
    if closure == '4bc':
        R[N] = f2[N]                               # f''(inf)=0  (4. mertebe kapanis)
    # closure='3bc': node N'de yonetici artik KORUNUR (BC dusurulur), beta=0 sart
    return R


def jacobian_unsteady(U, ops, beta, dt, invA, closure='4bc'):
    """Analitik Jacobian. Zaman katkilari: c1+=invA/dt, c3+=-beta*invA/dt."""
    D1, D2, D3, D4 = ops['D1'], ops['D2'], ops['D3'], ops['D4']
    N = ops['N']
    f1, f2, f3, f4 = D1 @ U, D2 @ U, D3 @ U, D4 @ U
    coef_local = beta * f4 - f2
    cU = beta * U
    c1 = -2.0 * beta * f3 + 2.0 * f1 + invA / dt
    c2 = 2.0 * beta * f2 - U
    c3 = -2.0 * beta * f1 - 1.0 - beta * invA / dt
    J = (np.diag(coef_local) + cU[:, None] * D4 + c1[:, None] * D1
         + c2[:, None] * D2 + c3[:, None] * D3)
    J[0, :] = 0.0; J[0, 0] = 1.0
    J[1, :] = D1[0, :]
    J[N - 1, :] = D1[N, :]
    if closure == '4bc':
        J[N, :] = D2[N, :]
    # '3bc': node N satiri yonetici-Jacobian olarak kalir
    return J


def newton_step(U0, ops, beta, dt, Uold, tau, k_amp, amp, omega, invA,
                outer='bh', closure='4bc', tol=1e-9, maxit=40):
    U = U0.copy(); prev = np.inf
    for it in range(maxit):
        R = residual_unsteady(U, ops, beta, dt, Uold, tau, k_amp, amp, omega, invA, outer, closure)
        nrm = np.max(np.abs(R))
        if not np.isfinite(nrm):
            return U, False, nrm
        if nrm < tol:
            return U, True, nrm
        J = jacobian_unsteady(U, ops, beta, dt, invA, closure)
        try:
            dU = np.linalg.solve(J, -R)
        except np.linalg.LinAlgError:
            return U, False, nrm
        lam = 1.0
        for _ in range(30):
            Rn = np.max(np.abs(residual_unsteady(U + lam * dU, ops, beta, dt, Uold,
                        tau, k_amp, amp, omega, invA, outer, closure)))
            if Rn < nrm or lam < 1e-6:
                break
            lam *= 0.5
        U = U + lam * dU
        if lam * np.max(np.abs(dU)) < 1e-12 or (nrm > 0.5 * prev and nrm < 1e-7):
            return U, True, nrm
        prev = nrm
    return U, (nrm < 1e-7), nrm


# ============================================================ tani / IC
def diag_fields(U, ops):
    """f''(0), M=max|f''|, gmin=min f', eta0 (ic isaret degisimi), fmin_interior,
    fpp_inf=f''(eta_max) (3-BC'de kapanis dayatilmadan ~0 olmali -> dogrulama ii)."""
    D1, D2 = ops['D1'], ops['D2']
    eta = ops['eta']
    N = ops['N']
    f1 = D1 @ U; f2 = D2 @ U
    fpp0 = float(f2[0])
    M = float(np.max(np.abs(f2)))
    gmin = float(f1.min())
    # ic turning point (eta>0'da f isaret degistiriyor mu)
    eta0 = np.nan
    for i in range(len(U) - 1):
        if eta[i] > 1e-9 and U[i] * U[i + 1] < 0:
            t = -U[i] / (U[i + 1] - U[i])
            eta0 = eta[i] + t * (eta[i + 1] - eta[i]); break
    interior = eta > 1e-9
    fmin_int = float(U[interior].min())
    fpp_inf = float(f2[N])
    return fpp0, M, gmin, eta0, fmin_int, fpp_inf


def steady_ic(beta, k_amp, a_init, N, eta_max):
    """tau0'daki IC: anlik dis akis a_init icin DURAGAN cozum (B&H konvansiyonu).
    solve_steady degismez cekirdek; a_init sabit (a(0)=1+Delta gibi)."""
    U, ops, _ = solve_steady(beta, k_amp, N=N, eta_max=eta_max, a=a_init)
    return U, ops


# ============================================================ MARS
def run_rothe(beta=0.0, k_amp=0.0, amp=0.5, omega=1.0, A=2.0,
              dt=0.005, tau0=0.0, tau_max=None, n_periods=None,
              N=160, eta_max=20.0, outer='bh', closure='4bc',
              ic='steady', M_blow=1.0e4, stop_on_blow=True, verbose=False):
    """B&H dis akisli Rothe marsi.

    sigma = omega/A (raporlanir);  zaman-katsayisi invA = 1/A.
    tau_max yoksa n_periods*2pi/omega kullanilir.
    Donus: dict(tau, fpp0, M, gmin, eta0, fmin_int, a, wall, rms, snaps,
                singular, tau_blow, ops, ...).
    """
    sigma = omega / A
    invA = 1.0 / A
    if tau_max is None:
        assert n_periods is not None, "tau_max veya n_periods ver"
        tau_max = tau0 + n_periods * (2 * np.pi / omega)
    a0 = a_outer(tau0, amp, omega, outer)
    w0 = wall_of_tau(tau0, k_amp, omega)
    if ic == 'zero':
        ops = build_operators(N, eta_max); U = np.zeros(N + 1)
    else:
        # 3bc'de IC steady cozucu 4bc'dir; baslangic profili olarak gecerli
        U, ops = steady_ic(0.0 if closure == '3bc' else beta, w0, a0, N, eta_max)
        if U is None:
            return None
    fpp0, M, gmin, eta0, fmin, fpp_inf = diag_fields(U, ops)
    out = dict(sigma=sigma, invA=invA, omega=omega, A=A, amp=amp, beta=beta,
               k_amp=k_amp, outer=outer, closure=closure, dt=dt, N=N, eta_max=eta_max,
               tau0=tau0,
               tau=[tau0], fpp0=[fpp0], M=[M], gmin=[gmin], eta0=[eta0],
               fmin_int=[fmin], fpp_inf=[fpp_inf], a=[a0], wall=[w0], rms=[0.0], snaps={},
               singular=False, tau_blow=None)
    snap_taus = [tau0 + d for d in (0.0, np.pi / 2, np.pi, 3 * np.pi / 2,
                                    2 * np.pi, 3 * np.pi, 4 * np.pi)]
    Uold = U.copy()
    nstep = int(round((tau_max - tau0) / dt))
    for n in range(1, nstep + 1):
        tau = tau0 + n * dt
        U, ok, nrm = newton_step(U, ops, beta, dt, Uold, tau, k_amp, amp, omega,
                                 invA, outer, closure)
        fpp0, M, gmin, eta0, fmin, fpp_inf = diag_fields(U, ops)
        blow = (not ok) or (not np.isfinite(M)) or (M > M_blow)
        if blow and stop_on_blow:
            out['singular'] = True; out['tau_blow'] = tau
            out['M_at_blow'] = M if np.isfinite(M) else np.inf
            out['fpp0_at_blow'] = fpp0 if np.isfinite(fpp0) else np.inf
            for key, val in (('tau', tau), ('fpp0', fpp0), ('M', M), ('gmin', gmin),
                             ('eta0', eta0), ('fmin_int', fmin), ('fpp_inf', fpp_inf),
                             ('a', a_outer(tau, amp, omega, outer)),
                             ('wall', wall_of_tau(tau, k_amp, omega)), ('rms', nrm)):
                out[key].append(val if np.isfinite(val) else np.nan)
            break
        Uold = U.copy()
        out['tau'].append(tau); out['fpp0'].append(fpp0); out['M'].append(M)
        out['gmin'].append(gmin); out['eta0'].append(eta0); out['fmin_int'].append(fmin)
        out['fpp_inf'].append(fpp_inf)
        out['a'].append(a_outer(tau, amp, omega, outer))
        out['wall'].append(wall_of_tau(tau, k_amp, omega)); out['rms'].append(nrm)
        for st in snap_taus:
            if abs(tau - st) < dt / 2 and st not in out['snaps']:
                out['snaps'][st] = (ops['eta'].copy(), U.copy())
        if verbose and n % 80 == 0:
            print(f"  tau={tau:.4f} f''0={fpp0:.5f} M={M:.4g} a={out['a'][-1]:+.3f} ||R||={nrm:.1e}")
    out['ops'] = ops
    for key in ('tau', 'fpp0', 'M', 'gmin', 'eta0', 'fmin_int', 'fpp_inf', 'a', 'wall', 'rms'):
        out[key] = np.array(out[key], dtype=float)
    return out


# ============================================================ tau_s tahmincisi
def _lin_intercept(x, y):
    """y = slope*x + ic dogrusal en-kucuk-kareler; x-kesisimi (y=0) ve R^2."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) < 4:
        return None, None, len(x)
    Amat = np.vstack([x, np.ones_like(x)]).T
    slope, ic = np.linalg.lstsq(Amat, y, rcond=None)[0]
    if abs(slope) < 1e-15:
        return None, None, len(x)
    R2 = 1.0 - np.sum((y - (slope * x + ic)) ** 2) / max(np.sum((y - y.mean()) ** 2), 1e-300)
    return -ic / slope, R2, len(x)


def tau_s_primary(r, M_active=1.0):
    """B&H-sadik birincil tau_s: G_m=min f_eta icin G_m ~ (tau_s-tau)^-1 (B&H Fig.3.1)
    => 1/|G_m| LINEER. AKTIF patlama penceresinde (M=max|f''| >= M_active, yani periyodik
    plato ~0.5'in ~2 kati ustu) 1/|G_m| dogrusal fit -> x-kesisimi tau_s.
    Donus: (tau_s, R2, npts)."""
    tau, gm, M = r['tau'], r['gmin'], r['M']
    m = np.isfinite(gm) & (gm < -0.5) & np.isfinite(M) & (M >= M_active)
    if m.sum() < 4:
        # patlama yoksa / yetersiz: tum ters-akis penceresini dene
        m = np.isfinite(gm) & (gm < -0.5)
    return _lin_intercept(tau[m], 1.0 / (-gm[m]))


def tau_s_M2(r, M_active=1.0):
    """Capraz-kontrol: F ~ (tau_s-tau)^-3/2 => max|f''| ~ (tau_s-tau)^-1/2 => 1/M^2 LINEER.
    Aktif patlama penceresinde 1/M^2 dogrusal fit -> x-kesisimi tau_s."""
    tau, M = r['tau'], r['M']
    m = np.isfinite(M) & (M >= M_active)
    if m.sum() < 4:
        return None, None, int(m.sum())
    return _lin_intercept(tau[m], 1.0 / M[m] ** 2)


def tau_s_band(r, lo=2.0):
    """Pencere-duyarlilik: 1/|G_m| fit'i |G_m|>=lo bandinda (alternatif tanim)."""
    tau, gm = r['tau'], r['gmin']
    m = np.isfinite(gm) & (-gm >= lo)
    return _lin_intercept(tau[m], 1.0 / (-gm[m]))


if __name__ == "__main__":
    import warnings; warnings.filterwarnings('ignore')
    print("bh_solver self-test")
    print(run_header(0.0, 0.0, 0.5, 1.0, 2.0, 'bh', '4bc', 0.01, 120, 14.0, 0.0))
    r = run_rothe(beta=0.0, k_amp=0.0, amp=0.5, omega=1.0, A=2.0, dt=0.01,
                  tau0=0.0, n_periods=2, N=120, eta_max=14.0)
    print(f"  periodic Delta=0.5 sigma={r['sigma']}: sing={r['singular']} "
          f"f''0 in [{r['fpp0'].min():.4f},{r['fpp0'].max():.4f}] a in [{r['a'].min():.3f},{r['a'].max():.3f}]")
    r = run_rothe(beta=0.0, k_amp=0.0, amp=2.0, omega=1.0, A=2.0, dt=0.005,
                  tau0=0.0, tau_max=4.2, N=160, eta_max=20.0)
    ts, R2, npt = tau_s_primary(r)
    tsb, _, _ = tau_s_band(r, lo=2.0)
    print(f"  singular Delta=2 sigma={r['sigma']}: sing={r['singular']} tau_blow={r['tau_blow']:.4f} "
          f"tau_s(1/Gm,M>=1)={ts:.4f}(R2={R2:.4f}) tau_s(1/Gm,|gm|>=2)={tsb:.4f}")
