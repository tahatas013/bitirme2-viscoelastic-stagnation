"""
Ultraspherical (Olver & Townsend 2013) spectral building blocks.

Coefficient-space operators on Chebyshev-T coefficients:
  Dm(m, N)   : differentiation T-coeffs -> C^(m)-coeffs (single m-th superdiag)
               d^m/dx^m T_n = 2^(m-1)(m-1)! n C^(m)_{n-m}
  Sconv(lam) : conversion C^(lam) -> C^(lam+1)  (banded: main + 2nd superdiag)
               S0: T -> C^(1)
  Jacobi(lam): multiplication-by-x in C^(lam) basis  (tridiagonal)
  Mult(a,lam): multiplication-by-g in C^(lam) basis, g given by T-coeffs a,
               via g(J_lam) = sum_j a_j T_j(J_lam) (Chebyshev recurrence)

Value<->coefficient transforms on the Chebyshev-Lobatto grid x_j=cos(pi j/(N-1)):
  vals2coeffs, coeffs2vals   (DCT-I)
  cheb_colloc_D(N)           : value-space 1st-derivative matrix (for residual /
                               variable-coefficient evaluation of a RESOLVED iterate;
                               applying D to a smooth function is accurate -- only
                               INVERTING D^4 is ill-conditioned, which we never do)
  gegenbauer_eval(lam,N,x)   : C^(lam)_n(x) matrix (used only in G1 tests)

All operators are returned as dense numpy arrays of shape (N, N); at the N used here
(<= ~512) dense assembly + solve is fast and the point is CONDITIONING, not sparsity.
A sparse path can replace these later without changing call sites.

Endpoint functionals (x=+1 is eta=0, x=-1 is eta=L):
  T_n(1)=1, T_n(-1)=(-1)^n
  T_n'(1)=n^2, T_n'(-1)=(-1)^{n+1} n^2
  T_n''(1)=n^2(n^2-1)/3, T_n''(-1)=(-1)^n n^2(n^2-1)/3
"""
import math
import numpy as np
from scipy.fft import dct


# ----------------------------------------------------------------- transforms
def vals2coeffs(v):
    """Values at x_j=cos(pi j/(N-1)) (j=0..N-1, so x_0=+1 .. x_{N-1}=-1) ->
    Chebyshev-T coefficients. DCT-I based."""
    v = np.asarray(v, float)
    N = v.shape[0]
    if N == 1:
        return v.copy()
    c = dct(v, type=1, axis=0) / (N - 1)
    c[0] *= 0.5
    c[-1] *= 0.5
    return c


def coeffs2vals(c):
    """Inverse of vals2coeffs (DCT-I is self-inverse up to scaling)."""
    c = np.asarray(c, float).copy()
    N = c.shape[0]
    if N == 1:
        return c.copy()
    c[0] *= 2.0
    c[-1] *= 2.0
    return 0.5 * dct(c, type=1, axis=0)


def cheb_points(N):
    """x_j = cos(pi j/(N-1)), j=0..N-1  (x_0=+1, x_{N-1}=-1)."""
    return np.cos(np.pi * np.arange(N) / (N - 1))


def cheb_colloc_D(N):
    """Trefethen value-space Chebyshev differentiation matrix on x_j (d/dx)."""
    if N == 1:
        return np.zeros((1, 1))
    x = cheb_points(N)
    cc = np.hstack([2.0, np.ones(N - 2), 2.0]) * (-1.0) ** np.arange(N)
    X = np.tile(x, (N, 1)).T
    dX = X - X.T
    D = np.outer(cc, 1.0 / cc) / (dX + np.eye(N))
    D = D - np.diag(D.sum(axis=1))
    return D


# ------------------------------------------------------------- diff operators
def Dm(m, N):
    """Differentiation T-coeffs -> C^(m)-coeffs. Single m-th superdiagonal."""
    if m == 0:
        return np.eye(N)
    const = 2.0 ** (m - 1) * math.factorial(m - 1)
    D = np.zeros((N, N))
    i = np.arange(N - m)
    D[i, i + m] = const * (i + m)
    return D


def Sconv(lam, N):
    """Conversion C^(lam) -> C^(lam+1). lam=0 is T -> C^(1)."""
    S = np.zeros((N, N))
    i = np.arange(N)
    if lam == 0:
        S[i, i] = 0.5
        S[0, 0] = 1.0
        j = np.arange(N - 2)
        S[j, j + 2] = -0.5
    else:
        S[i, i] = lam / (i + lam)
        j = np.arange(N - 2)
        S[j, j + 2] = -lam / (j + lam + 2.0)
    return S


def Jacobi(lam, N):
    """Multiplication-by-x in C^(lam) basis (tridiagonal).
    x C^(lam)_n = a_n C_{n+1} + g_n C_{n-1},
       a_n=(n+1)/(2(n+lam)),  g_n=(n+2 lam-1)/(2(n+lam))."""
    J = np.zeros((N, N))
    # sub-diagonal J[i,i-1] = a_{i-1} = i/(2(i-1+lam))
    i = np.arange(1, N)
    J[i, i - 1] = i / (2.0 * (i - 1 + lam))
    # super-diagonal J[i,i+1] = g_{i+1} = (i+2 lam)/(2(i+1+lam))
    i = np.arange(N - 1)
    J[i, i + 1] = (i + 2.0 * lam) / (2.0 * (i + 1 + lam))
    return J


def Mult(a, lam, N, tol=1e-14):
    """Multiplication-by-g operator in C^(lam) basis; g has Chebyshev-T coeffs a.
    M = sum_j a_j T_j(J),  T_0=I, T_1=J, T_{j+1}=2 J T_j - T_{j-1}.
    Truncates the (smooth) coefficient series at tol*max|a| for efficiency."""
    a = np.asarray(a, float)
    J = Jacobi(lam, N)
    amax = np.max(np.abs(a)) if a.size else 0.0
    jmax = 0
    if amax > 0:
        nz = np.where(np.abs(a) > tol * amax)[0]
        jmax = int(nz[-1]) if nz.size else 0
    M = a[0] * np.eye(N) if a.size else np.zeros((N, N))
    if jmax >= 1:
        Pprev = np.eye(N)
        Pcur = J.copy()
        M = M + a[1] * Pcur
        for j in range(2, jmax + 1):
            Pnext = 2.0 * (J @ Pcur) - Pprev
            M = M + a[j] * Pnext
            Pprev, Pcur = Pcur, Pnext
    return M


# ------------------------------------------------------ Gegenbauer evaluation
def gegenbauer_eval(lam, N, x):
    """Matrix E[k,n] = C^(lam)_n(x_k), n=0..N-1.  lam>=1.  (G1 tests only.)"""
    x = np.asarray(x, float)
    E = np.zeros((x.size, N))
    E[:, 0] = 1.0
    if N > 1:
        E[:, 1] = 2.0 * lam * x
    for n in range(1, N - 1):
        E[:, n + 1] = (2.0 * (n + lam) * x * E[:, n]
                       - (n + 2.0 * lam - 1.0) * E[:, n - 1]) / (n + 1.0)
    return E


# ------------------------------------------------------- endpoint functionals
def bc_rows(N):
    """Return dicts of length-N functional row-vectors acting on T-coeffs:
       val_p1 = f(+1) ; d1_p1 = f'(+1) ; etc. (x-derivatives, unmapped)."""
    n = np.arange(N)
    return dict(
        val_p1=np.ones(N),
        val_m1=(-1.0) ** n,
        d1_p1=n ** 2,
        d1_m1=(-1.0) ** (n + 1) * n ** 2,
        d2_p1=n ** 2 * (n ** 2 - 1) / 3.0,
        d2_m1=(-1.0) ** n * n ** 2 * (n ** 2 - 1) / 3.0,
    )


# =============================================================== G1 self-test
if __name__ == "__main__":
    import numpy.polynomial.chebyshev as C
    np.set_printoptions(suppress=True)
    print("# G1 -- ultraspherical operator unit tests\n")
    ok_all = True

    def check(name, err, tol=1e-12):
        global ok_all
        good = err < tol
        ok_all = ok_all and good
        print(f"  [{'PASS' if good else 'FAIL'}] {name:42s} err={err:.2e} (tol {tol:.0e})")

    N = 64
    x = cheb_points(N)

    # transforms round-trip
    v = np.exp(np.sin(3 * x))
    check("vals2coeffs/coeffs2vals round-trip",
          np.max(np.abs(coeffs2vals(vals2coeffs(v)) - v)))

    # D1 vs numpy chebder, then evaluate
    g = np.exp(0.7 * x) + np.cos(2 * x)
    a = vals2coeffs(g)
    for m, fac in [(1, 1), (2, 1), (3, 1), (4, 1)]:
        # exact m-th derivative coeffs (in T) via numpy
        aT = a.copy()
        dexact = C.chebder(aT, m)
        dexact = np.concatenate([dexact, np.zeros(N - dexact.size)])
        # ultraspherical: Dm gives C^(m) coeffs; convert back to T by checking values
        cm = Dm(m, N) @ a                       # C^(m) coeffs of d^m g/dx^m
        E = gegenbauer_eval(m, N, x)
        vals_us = E @ cm                         # values of the m-th derivative
        vals_ex = coeffs2vals(dexact)
        check(f"D{m} (T->C^{m}) m-th derivative values",
              np.max(np.abs(vals_us - vals_ex)), tol=1e-9)

    # conversion S0: T->C^1 ; check by evaluating C^1 series
    h = np.cos(1.3 * x) * np.exp(0.2 * x)
    ah = vals2coeffs(h)
    c1 = Sconv(0, N) @ ah
    E1 = gegenbauer_eval(1, N, x)
    check("S0 conversion T->C^1 values", np.max(np.abs(E1 @ c1 - h)), tol=1e-10)
    # chain S0..S3: T->C^4, evaluate
    c4 = Sconv(3, N) @ (Sconv(2, N) @ (Sconv(1, N) @ (Sconv(0, N) @ ah)))
    E4 = gegenbauer_eval(4, N, x)
    check("S3 S2 S1 S0 conversion T->C^4 values",
          np.max(np.abs(E4 @ c4 - h)), tol=1e-9)

    # multiplication in C^4 : M4[g] (c4 coeffs of h) == c4 coeffs of g*h
    gg = np.exp(0.5 * x)
    ag = vals2coeffs(gg)                          # T-coeffs of g
    M4 = Mult(ag, 4, N)
    prod_c4 = M4 @ c4                             # should be C^4 coeffs of g*h
    check("M4[g] multiplication in C^4 values",
          np.max(np.abs(E4 @ prod_c4 - gg * h)), tol=1e-9)
    # also in C^0 (=T) basis via M with lam acting on T? test lam=2 too
    c2 = Sconv(1, N) @ (Sconv(0, N) @ ah)
    M2 = Mult(ag, 2, N)
    E2 = gegenbauer_eval(2, N, x)
    check("M2[g] multiplication in C^2 values",
          np.max(np.abs(E2 @ (M2 @ c2) - gg * h)), tol=1e-9)

    # Jacobi: multiplication by x in C^4
    Mx = Mult(vals2coeffs(x), 4, N)
    check("Jacobi J4 == Mult[x] in C^4", np.max(np.abs(Mx - Jacobi(4, N))), tol=1e-12)

    # endpoint functionals vs direct evaluation
    bc = bc_rows(N)
    af = vals2coeffs(np.exp(0.4 * x))
    fp = C.chebder(af, 1)
    fpp = C.chebder(af, 2)
    check("f(+1) functional", abs(bc["val_p1"] @ af - np.exp(0.4)), tol=1e-11)
    check("f(-1) functional", abs(bc["val_m1"] @ af - np.exp(-0.4)), tol=1e-11)
    check("f'(+1) functional", abs(bc["d1_p1"] @ af - C.chebval(1.0, fp)), tol=1e-9)
    check("f'(-1) functional", abs(bc["d1_m1"] @ af - C.chebval(-1.0, fp)), tol=1e-9)
    check("f''(+1) functional", abs(bc["d2_p1"] @ af - C.chebval(1.0, fpp)), tol=1e-7)
    check("f''(-1) functional", abs(bc["d2_m1"] @ af - C.chebval(-1.0, fpp)), tol=1e-7)

    print(f"\nG1 OVERALL: {'PASS' if ok_all else 'FAIL'}")
