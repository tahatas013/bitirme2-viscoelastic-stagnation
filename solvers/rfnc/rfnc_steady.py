"""
RFNC -- Residual-Form Newton Collocation, STEADY cozucu.

Fikir: eq.39'un denklemini f'''' icin COZMEDEN, artik (residual) formda tut.
Steady, a=1 icin (eq.39'dan turetilmis):

    R = beta*f*f'''' - 2*beta*f'*f''' + beta*(f'')^2 - f''' + (f')^2 - f*f'' - 1 = 0

Hicbir yerde f'ye / beta*f'ye BOLME yok. Iç turning point'te (f=0) beta*f*f''''
terimi kendiliginden kaybolur, o satir uyumluluk sartina (N=0) doner; singulerlik
hic ortaya cikmaz.

Newtonian limit (beta=0):  -f''' + (f')^2 - f*f'' - 1 = 0
   <=> f''' + f*f'' - (f')^2 + 1 = 0   (Hiemenz, f''(0)=1.232588)

Yontem: Chebyshev-Gauss-Lobatto collocation, eta in [0, eta_max], lineer map.
Bilinmeyen: dugumlerdeki f degerleri. Ic dugumlerde governing residual; 4 satir
BC ile degistirilir. Newton: J dU = -R, analitik Jacobian (D matrislerinden).
Sureklilik: beta 0->hedef (Hiemenz baslangic), sonra c 0->hedef.

Sinir kosullari (eq.40-42, DEGISMEZ):
   f'(0)=0,  f(0)=c,  f'(eta_max)=a(=1),  f''(eta_max)=0
"""
import numpy as np


# ------------------------------------------------------------------ Chebyshev
def cheb(N):
    """Trefethen cheb: D (N+1 x N+1) ve x (N+1) noktalari [-1,1], x_0=1..x_N=-1."""
    if N == 0:
        return np.array([[0.0]]), np.array([1.0])
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2.0, np.ones(N - 1), 2.0]) * (-1.0) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1.0 / c) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    return D, x


def build_operators(N, eta_max):
    """CGL turev matrisleri eta-uzayinda (lineer map). Geri donus: dict."""
    Dc, x = cheb(N)
    eta = eta_max * (1.0 - x) / 2.0          # x=1->eta=0, x=-1->eta=eta_max
    s = -2.0 / eta_max                        # dx/deta (sabit, lineer map)
    D1 = s * Dc
    D2 = s ** 2 * (Dc @ Dc)
    D3 = s ** 3 * np.linalg.matrix_power(Dc, 3)
    D4 = s ** 4 * np.linalg.matrix_power(Dc, 4)
    return dict(N=N, eta_max=eta_max, eta=eta, D1=D1, D2=D2, D3=D3, D4=D4)


# ------------------------------------------------------- residual & jacobian
def residual_steady(U, beta, c, ops, a=1.0):
    D1, D2, D3, D4 = ops['D1'], ops['D2'], ops['D3'], ops['D4']
    N = ops['N']
    f1 = D1 @ U
    f2 = D2 @ U
    f3 = D3 @ U
    f4 = D4 @ U
    # governing residual (carpim formu -- f=0'da beta*U*f4 kendiliginden sonar)
    R = (beta * U * f4 - 2.0 * beta * f1 * f3 + beta * f2 ** 2
         - f3 + f1 ** 2 - U * f2 - a ** 2)
    # BC satirlari (eq.40-42)
    R[0] = U[0] - c                # f(0) = c
    R[1] = f1[0]                   # f'(0) = 0
    R[N - 1] = f1[N] - a           # f'(eta_max) = a
    R[N] = f2[N]                   # f''(eta_max) = 0
    return R


def jacobian_steady(U, beta, c, ops, a=1.0):
    D1, D2, D3, D4 = ops['D1'], ops['D2'], ops['D3'], ops['D4']
    N = ops['N']
    f1 = D1 @ U
    f2 = D2 @ U
    f3 = D3 @ U
    f4 = D4 @ U
    # analitik: J = diag(coef_local) + diag(cU)D4 + diag(c1)D1 + diag(c2)D2 + diag(c3)D3
    coef_local = beta * f4 - f2          # d/dU_j (lokal: beta*f*f4 ve -f*f2'nin f-kismi)
    cU = beta * U                        # D4 katsayisi
    c1 = -2.0 * beta * f3 + 2.0 * f1     # D1 katsayisi
    c2 = 2.0 * beta * f2 - U             # D2 katsayisi
    c3 = -2.0 * beta * f1 - 1.0          # D3 katsayisi
    J = (np.diag(coef_local)
         + cU[:, None] * D4
         + c1[:, None] * D1
         + c2[:, None] * D2
         + c3[:, None] * D3)
    # BC satirlari
    J[0, :] = 0.0; J[0, 0] = 1.0         # f(0)=c
    J[1, :] = D1[0, :]                   # f'(0)=0
    J[N - 1, :] = D1[N, :]               # f'(eta_max)=a
    J[N, :] = D2[N, :]                   # f''(eta_max)=0
    return J


# --------------------------------------------------------------- Newton solve
def newton_solve(U0, beta, c, ops, a=1.0, tol=1e-9, maxit=50, verbose=False):
    """Newton: J dU = -R. Yakinsama: ||R||inf<tol VEYA adim cok kucuk (roundoff
    zemini -- yuksek mertebe spektral D4'un kosullanmasi residual'i ~1e-8'de
    tutabilir, ama cozum bu noktada zaten makineye yakin)."""
    U = U0.copy()
    prev = np.inf
    for it in range(maxit):
        R = residual_steady(U, beta, c, ops, a)
        nrm = np.max(np.abs(R))
        if verbose:
            print(f"    it{it:2d}  ||R||inf = {nrm:.3e}")
        if nrm < tol:
            return U, True, it, nrm
        J = jacobian_steady(U, beta, c, ops, a)
        try:
            dU = np.linalg.solve(J, -R)
        except np.linalg.LinAlgError:
            return U, False, it, nrm
        # damping: tam adim residual'i artiriyorsa kucult
        lam = 1.0
        for _ in range(30):
            Un = U + lam * dU
            Rn = np.max(np.abs(residual_steady(Un, beta, c, ops, a)))
            if Rn < nrm or lam < 1e-6:
                break
            lam *= 0.5
        step = lam * np.max(np.abs(dU))
        U = U + lam * dU
        # roundoff zeminine oturma: adim kucuk ve residual artik dusmuyorsa kabul
        if step < 1e-12 or (nrm > 0.5 * prev and nrm < 1e-7):
            return U, True, it, nrm
        prev = nrm
    R = residual_steady(U, beta, c, ops, a)
    return U, (np.max(np.abs(R)) < 1e-7), maxit, np.max(np.abs(R))


def hiemenz_init(ops):
    eta = ops['eta']
    return eta - 1.0 + np.exp(-eta)


# ----------------------------------------------------------- continuation API
def solve_steady(beta_target, c_target, N=60, eta_max=10.0, a=1.0,
                 n_beta=6, dc=0.01, tol=1e-9, verbose=False):
    """beta: 0->hedef (Hiemenz baslangic), sonra c: 0->hedef. Fiziksel branch."""
    ops = build_operators(N, eta_max)
    U = hiemenz_init(ops)
    log = []
    # --- beta continuation (c=0) ---
    for b in np.linspace(0.0, beta_target, n_beta):
        U, ok, it, nrm = newton_solve(U, b, 0.0, ops, a, tol, verbose=verbose)
        log.append(('beta', b, 0.0, ok, it, nrm))
        if not ok:
            return None, ops, log
    # --- c continuation (beta sabit) ---
    if abs(c_target) > 1e-15:
        nstep = max(1, int(round(abs(c_target) / dc)))
        for c in np.linspace(0.0, c_target, nstep + 1)[1:]:
            U, ok, it, nrm = newton_solve(U, beta_target, c, ops, a, tol, verbose=verbose)
            log.append(('c', beta_target, c, ok, it, nrm))
            if not ok:
                return None, ops, log
    return U, ops, log


def compat_residual(U, ops, beta, a=1.0):
    """N_compat(eta): R = beta*f*f'''' - N_compat. f=0'da R=-N_compat, R=0 cozuldugu
    icin N_compat(eta0)=0 OTOMATIK saglanir. Ayrica f'''' profilini de dondur."""
    D1, D2, D3, D4 = ops['D1'], ops['D2'], ops['D3'], ops['D4']
    f1, f2, f3, f4 = D1 @ U, D2 @ U, D3 @ U, D4 @ U
    Ncomp = 2 * beta * f1 * f3 - beta * f2 ** 2 + f3 - f1 ** 2 + U * f2 + a ** 2
    return Ncomp, f4


def bary_interp(U, ops, xq):
    """CGL dugumlerde ornek alinmis U'yu xq noktalarinda barycentric interpolasyon."""
    eta = ops['eta']
    N = ops['N']
    j = np.arange(N + 1)
    w = (-1.0) ** j
    w[0] *= 0.5; w[-1] *= 0.5
    xq = np.atleast_1d(np.asarray(xq, float))
    out = np.empty_like(xq)
    for m, xx in enumerate(xq):
        d = xx - eta
        if np.any(np.abs(d) < 1e-14):
            out[m] = U[np.argmin(np.abs(d))]
        else:
            ww = w / d
            out[m] = np.dot(ww, U) / np.sum(ww)
    return out if out.size > 1 else out[0]


def diagnostics(U, ops, beta, a=1.0):
    """f''(0), eta0 (varsa), ve turning point'te uyumluluk N(eta0)."""
    D1, D2 = ops['D1'], ops['D2']
    eta = ops['eta']
    fpp0 = (D2 @ U)[0]
    # eta0: f isaret degistirir mi (ic sifir)
    eta0 = None
    for i in range(len(U) - 1):
        if eta[i] > 1e-9 and U[i] * U[i + 1] < 0:
            t = -U[i] / (U[i + 1] - U[i])
            eta0 = eta[i] + t * (eta[i + 1] - eta[i])
            break
    return fpp0, eta0


if __name__ == "__main__":
    # hizli kendi-testi
    for beta, c in [(0.0, 0.0), (0.1, 0.0), (0.1, -0.05)]:
        U, ops, log = solve_steady(beta, c, N=120, eta_max=10.0)
        if U is None:
            print(f"beta={beta} c={c}: FAILED")
            continue
        fpp0, eta0 = diagnostics(U, ops, beta)
        e0s = f"{eta0:.4f}" if eta0 else "yok"
        print(f"beta={beta} c={c}:  f''(0)={fpp0:.6f}  eta0={e0s}")
