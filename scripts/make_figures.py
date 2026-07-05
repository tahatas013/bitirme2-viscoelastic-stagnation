# -*- coding: utf-8 -*-
"""
make_figures.py -- tez figurlerinin TEK giris noktasi.

Kullanim (depo kokunden):
    python scripts/make_figures.py            # tum figurler
    python scripts/make_figures.py --only fig_4_4
    python scripts/make_figures.py --list

Ilkeler:
  * Sertifikali cozucu dosyalari (solvers/) HIC degistirilmez; yalnizca import.
  * Sertifikali onbellek (data/) tercih edilir; onbellekte olmayan periyodik
    kosular, kaynak paketin BELGELENMIS parametreleriyle yeniden uretilir ve
    data/ altina kaydedilir (bkz. PROVENANCE.md).
  * tau_s disiplini: mutlak tau_s tek sayi olarak RAPORLANMAZ; bant [3.38,3.41]
    ve fark Delta-tau_s (1/|G_m| ust + eta_m^-2 alt, hücre ara-noktası).
  * Yunanca/matematik SADECE mathtext (r'$\\beta$' vb.); figur metni Turkce.
"""
import argparse
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for sub in ("solvers/bh_validation", "solvers/rfnc", "solvers/secondgrade",
            "solvers/beta_perturbation"):
    p = os.path.join(ROOT, sub)
    if p not in sys.path:
        sys.path.insert(0, p)
DATA = os.path.join(ROOT, "data")
FIGS = os.path.join(ROOT, "figures")
os.makedirs(FIGS, exist_ok=True)

import config_sg as C            # noqa: E402  (sertifikali, import-only)
import solver_a_sg as SA         # noqa: E402
import bh_solver as bh           # noqa: E402  (sertifikali, import-only)
from rfnc_steady import solve_steady  # noqa: E402

PER = 2.0 * np.pi

# ---- ortak stil ------------------------------------------------------------
plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "dejavuserif",
    "axes.grid": True,
    "grid.alpha": 0.3,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})
_cm = plt.cm.plasma
CB = {0.0: "#000000", 0.1: _cm(0.15), 0.2: _cm(0.50), 0.3: _cm(0.78)}
PHASE_C = ["#1f77b4", "#2ca02c", "#d62728", "#9467bd"]   # 4 faz rengi

DPI = 300


def _save(fig, name):
    path = os.path.join(FIGS, name)
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return path


# ---- veri katmani ----------------------------------------------------------
def settled_history(beta):
    """Periyodik rejim (Delta=0.5, sigma=0.5) son tam periyot.

    Onbellek yoksa, kaynak paketin (newt_sol/secondgrade_bh/thesis_figs.py,
    _settled_history) BELGELENMIS sertifikali parametreleriyle yeniden uretilir:
    n_periods=7, eta_max=20, N=200, dt=0.005, Solver A (import-only).
    """
    p = os.path.join(DATA, f"periodic_beta{beta:g}.npz")
    if os.path.exists(p):
        d = np.load(p)
        return d["tg"], d["fpp0"], d["a"]
    r = SA.march(0.5, 0.5, beta, 20.0, 200, 7 * PER, dt=0.005)
    if r is None or r["singular"]:
        raise RuntimeError(f"periyodik kosu beklenmedik sekilde durdu (beta={beta})")
    t, fpp, a = r["tau"], r["fpp0"], r["a"]
    nP = int(t[-1] // PER)
    m = (t >= (nP - 1) * PER) & (t < nP * PER)
    tg = t[m] - (nP - 1) * PER
    np.savez(p, tg=tg, fpp0=fpp[m], a=a[m])
    return tg, fpp[m], a[m]


def profile_snaps(beta):
    """Hiz profilleri f'(eta), 4 faz (tau = 5*2pi + [0, pi/2, pi, 3pi/2]).

    Onbellek yoksa kaynak paketin (figs_thesis.py) belgelenmis parametreleriyle
    uretilir: eta_max=22, N=200, dt=0.005, SA.snapshot_profiles (import-only).
    """
    p = os.path.join(DATA, f"profiles_beta{beta:g}.npz")
    phases = [0.0, np.pi / 2, np.pi, 3 * np.pi / 2]
    if os.path.exists(p):
        d = np.load(p)
        return d["eta"], [d[f"fp{i}"] for i in range(4)], phases
    tq = [5 * PER + ph for ph in phases]
    snaps = SA.snapshot_profiles(0.5, 0.5, beta, 22.0, 200, tq, dt=0.005)
    if snaps is None:
        raise RuntimeError(f"snapshot_profiles IC basarisiz (beta={beta})")
    keys = sorted(snaps.keys())
    eta = snaps[keys[0]][0]
    o = np.argsort(eta)
    fps = [snaps[k][2][o] for k in keys]
    np.savez(p, eta=eta[o], **{f"fp{i}": fp for i, fp in enumerate(fps)})
    return eta[o], fps, phases


def fourier_amp_phase(tg, y):
    """y = c0 + c1 cos(tau) + s1 sin(tau) LSQ (thesis_figs.py D3 ile ayni)."""
    A = np.vstack([np.ones_like(tg), np.cos(tg), np.sin(tg)]).T
    c0, c1, s1 = np.linalg.lstsq(A, y, rcond=None)[0]
    return float(np.hypot(c1, s1)), float(np.arctan2(s1, c1)), float(c0)


def shoelace(a, y):
    """Dongu alani (thesis_figs.py D2 ile ayni)."""
    return 0.5 * abs(np.sum(a * np.roll(y, -1) - np.roll(a, -1) * y))


def lin_fit_r2(x, y):
    A = np.vstack([x, np.ones_like(x)]).T
    (s, c), _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    r = y - (s * x + c)
    R2 = 1.0 - np.sum(r ** 2) / max(np.sum((y - y.mean()) ** 2), 1e-300)
    return s, c, R2


# ---- figurler ---------------------------------------------------------------
def fig_4_1():
    """(a) Duragan Hiemenz f'(eta) + f''(0) capasi; (b) tekil durum 1/|G_m|."""
    # (a) dogrulama calismasinin belgelenmis cagrisi: N=120, eta_max=12, a=1
    U, ops, _ = solve_steady(0.0, 0.0, N=120, eta_max=12.0, a=1.0)
    eta = ops["eta"]
    fp = ops["D1"] @ U
    fpp0 = float((ops["D2"] @ U)[0])
    o = np.argsort(eta)
    ok_a = abs(fpp0 - 1.232588) < 1e-5

    d = np.load(os.path.join(DATA, "stage2_singular.npz"))
    tau, gmin, M = d["tau"], d["gmin"], d["M"]
    m = np.isfinite(gmin) & (gmin < -0.5) & np.isfinite(M) & (M >= 1.0)
    s, c, R2 = lin_fit_r2(tau[m], 1.0 / (-gmin[m]))
    x_int = -c / s
    ok_b = 3.37 <= x_int <= 3.42   # sertifikalı bant [3.38,3.41] (+dt paylari)

    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    ax[0].plot(eta[o], fp[o], "-", c="k", lw=2)
    ax[0].axhline(1.0, ls=":", c="gray", lw=1)
    ax[0].set_xlabel(r"$\eta$")
    ax[0].set_ylabel(r"$f'(\eta)$")
    ax[0].set_xlim(0, 6)
    ax[0].set_title(r"(a) Durağan Hiemenz profili ($\beta=0$, $a=1$)")
    ax[0].annotate(r"$f''(0)=1{,}232588$" "\n(referans değer)",
                   xy=(2.6, 0.35), fontsize=10,
                   bbox=dict(fc="white", ec="gray", alpha=0.9))
    ax[1].plot(tau[m], 1.0 / (-gmin[m]), "o", ms=3, c="k",
               label=r"$1/|G_m|$ değerleri")
    tt = np.linspace(tau[m][0], 3.41, 50)
    ax[1].plot(tt, s * tt + c, "-", c="C3", lw=1.5, label="doğrusal uyum")
    ax[1].axvspan(3.38, 3.41, color="C2", alpha=0.25,
                  label=r"referans aralık $[3{,}38\,;\,3{,}41]$")
    ax[1].set_xlabel(r"$\tau$")
    ax[1].set_ylabel(r"$1/|G_m|$")
    ax[1].set_ylim(bottom=0)
    ax[1].set_title(r"(b) Tekil durum $\Delta=2$, $\sigma=0.5$: $\tau_s$'e doğrusal yaklaşım")
    ax[1].legend(fontsize=8)
    fig.suptitle(r"Blyth ve Hall probleminin doğrulanması ($\beta=0$)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    _save(fig, "fig_4_1_bh_validation.png")
    key = f"f''(0)={fpp0:.6f}; 1/|G_m| fit kesisimi bant icinde: {'EVET' if ok_b else 'HAYIR'} (R2={R2:.3f})"
    return ("OK", key, ok_a and ok_b)


def _beta0_load(name):
    """beta->0 doğrulama önbelleği (data/beta0/, bayt-özdeş kopyalar)."""
    return np.load(os.path.join(DATA, "beta0", name + ".npz"))


def _beta0_lastper(y):
    """Son tam periyot [5P, 6P] (dt = 2pi/1256, 6 periyot -> kesin indeksler)."""
    return np.asarray(y)[5 * 1256:6 * 1256 + 1]


def fig_4_2():
    """Şekil 4.2 (beta0_slope): beta->0 O(beta) yakınsaması, Türkçe metinle,
    data/beta0/ önbelleğinden yeniden üretim (sayılar birebir aynı)."""
    bb = np.array([0.012, 0.018, 0.025, 0.035, 0.05, 0.07, 0.1, 0.14, 0.2])
    ref = {N: _beta0_lastper(_beta0_load(f"A_N{N}_b0")["fpp0"])
           for N in (180, 360)}
    E = {N: np.array([float(np.max(np.abs(
        _beta0_lastper(_beta0_load(f"A_N{N}_b{b:g}")["fpp0"]) - ref[N])))
        for b in bb]) for N in (180, 360)}
    gates_ok = bool(np.all(np.abs(E[180] - E[360]) / E[360] <= 0.05))
    logb, logE = np.log(bb), np.log(E[360])
    Amat = np.vstack([logb, np.ones_like(logb)]).T
    coef, _, _, _ = np.linalg.lstsq(Amat, logE, rcond=None)
    slope, ic = float(coef[0]), float(coef[1])
    resid = logE - Amat @ coef
    s2 = np.sum(resid ** 2) / (len(bb) - 2)
    se = float(np.sqrt((s2 * np.linalg.inv(Amat.T @ Amat))[0, 0]))
    pair = float(np.log(E[360][1] / E[360][0]) / np.log(bb[1] / bb[0]))
    pol = np.polyfit(bb[:5], E[360][:5] / bb[:5], 2)
    C_A = float(pol[-1])
    C_B = float(np.max(np.abs(_beta0_lastper(_beta0_load("B_pert")["f1pp0"]))))
    ok = (gates_ok and abs(pair - 0.9660) < 2e-3
          and abs(C_A - 3.1357) < 2e-3 and abs(C_B - 3.1357) < 2e-3)

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.5, 5.0))
    ax.loglog(bb, E[360], "o", ms=8, color="C0",
              label=r"$E_{2N}(\beta)$  ($N=360$)")
    ax.loglog(bb, E[180], "s", ms=4, mfc="none", color="C7",
              label=r"$E_{N}(\beta)$  ($N=180$)")
    bfit = np.array([bb.min(), bb.max()])
    ax.loglog(bfit, np.exp(ic) * bfit ** slope, "-", color="C0", lw=1.5,
              label=rf"doğrusal uyum (eğim ${slope:.3f}\pm{se:.3f}$)")
    ax.loglog(bfit, C_B * bfit, "--", color="C2", lw=1.2,
              label=r"$C_B\,\beta$ (Çözücü B)")
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"$E(\beta)=\max_\tau|f''(0)_\beta - f''(0)_0|$")
    ax.set_title(r"(a) $E(\beta)$ farkının $\beta\to 0$ davranışı")
    ax.annotate(rf"asimptotik eğim $\approx {pair:.3f}$" "\n(beklenen 1)",
                xy=(0.05, 0.72), xycoords="axes fraction", fontsize=10,
                bbox=dict(fc="white", ec="gray", alpha=0.9))
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8, loc="lower right")
    rat = E[360] / bb
    ax2.plot(bb, rat, "o", ms=7, color="C0", label=r"$E_{2N}(\beta)/\beta$")
    bq = np.linspace(0.0, bb.max() * 1.02, 200)
    ax2.plot(bq, np.polyval(pol, bq), "-", color="C0", lw=1.2, alpha=0.7,
             label=r"ikinci derece uyum (5 en küçük $\beta$)")
    ax2.plot(0.0, C_A, "d", ms=9, color="C0")
    ax2.axhline(C_B, color="C2", ls="--", lw=1.2, label=r"$C_B$ (Çözücü B)")
    ax2.set_xlabel(r"$\beta$")
    ax2.set_ylabel(r"$E_{2N}(\beta)/\beta$")
    ax2.set_title(r"(b) Öncül katsayının kestirimi")
    ax2.annotate(rf"$C_A = {C_A:.3f}$ (Çözücü A)" "\n"
                 rf"$C_B = {C_B:.3f}$ (Çözücü B)",
                 xy=(0.35, 0.72), xycoords="axes fraction", fontsize=10,
                 bbox=dict(fc="white", ec="gray", alpha=0.9))
    ax2.legend(fontsize=8, loc="lower left")
    fig.suptitle(r"$\beta\to 0$ yakınsaması ($\Delta=0.5$, $\sigma=0.5$, $k=0$)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    _save(fig, "beta0_slope.png")
    key = (f"asimptotik eğim {pair:.4f}; C_A={C_A:.4f} C_B={C_B:.4f}; "
           f"eğim(global) {slope:.4f}±{se:.4f}")
    return ("OK", key, bool(ok))


def fig_4_3():
    """Şekil 4.3 (bh_anchor_loop): f''(0)-a(tau) çevrimi, Türkçe metinle,
    data/beta0/ önbelleğinden yeniden üretim."""
    dA = _beta0_load("A_N180_b0")
    fA, aA = _beta0_lastper(dA["fpp0"]), _beta0_lastper(dA["a"])
    f3 = _beta0_lastper(_beta0_load("BH_3bc_N180")["fpp0"])
    f4 = _beta0_lastper(_beta0_load("BH_4bc_N180")["fpp0"])
    rel3 = float(np.max(np.abs(fA - f3) / np.abs(f3)))
    ident4 = float(np.max(np.abs(fA - f4)))
    ok = rel3 < 1e-2 and ident4 == 0.0
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(aA, fA, "-", lw=2.2, color="C0",
            label=r"viskoelastik çözücü ($\beta = 0$)")
    ax.plot(aA, f3, "--", lw=1.6, color="C1",
            label="Newtonian çözücü (Blyth ve Hall)")
    ax.plot(aA, f4, ":", lw=1.4, color="C2", label="_nolegend_")
    ax.set_xlabel(r"$a(\tau)$")
    ax.set_ylabel(r"$f''(0,\tau)$")
    ax.set_title(r"$f''(0)$–$a(\tau)$ çevrimi, son periyot "
                 r"($\Delta = 0{,}5$; $\sigma = 0{,}5$)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    _save(fig, "bh_anchor_loop.png")
    key = f"çevrim üzerinde en büyük bağıl fark {rel3:.2e} (< 1e-2)"
    return ("OK", key, bool(ok))


def fig_4_4():
    """Delta-tau_s vs beta: sertifikali hucre ara-nokta tahmini + zarf."""
    with open(os.path.join(DATA, "taus_results.json")) as fh:
        J = json.load(fh)
    bb = [0.1, 0.2, 0.3]
    mid = np.array([J[f"{b}"]["dts_mid"] for b in bb])
    lo = np.array([J[f"{b}"]["dts_lo"] for b in bb])
    hi = np.array([J[f"{b}"]["dts_hi"] for b in bb])
    cert = np.array([0.022, 0.052, 0.087])
    ok = bool(np.all(np.abs(mid - cert) < 5e-4))

    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ax.errorbar(bb, mid, yerr=[mid - lo, hi - mid], fmt="D-", c="k", lw=2,
                capsize=4, ms=7,
                label=r"$\Delta\tau_s$ hücre ara-noktası ($\eta_{max}$-değişmez)")
    ax.fill_between(bb, lo, hi, color="gray", alpha=0.15,
                    label=r"gevşek zarf: $1/|G_m|$ (üst) / $\eta_m^{-2}$ (alt) tahmincileri")
    for b, v in zip(bb, mid):
        ax.annotate(rf"$+{v:.3f}$", xy=(b, v), xytext=(b + 0.006, v - 0.012),
                    fontsize=10)
    ax.axhline(0.0, c="k", lw=0.7)
    ax.set_xlabel(r"$\beta$ (Deborah sayısı)")
    ax.set_ylabel(r"$\Delta\tau_s=\tau_s(\beta)-\tau_s(0)$")
    ax.set_title(r"Viskoelastisite tekilliği geciktirir: $\Delta\tau_s(\beta)>0$"
                 "\n" r"($\Delta=2$, $\sigma=0.5$; çözünürlük-yakınsak yeniden hesap)")
    ax.legend(fontsize=8, loc="upper left")
    _save(fig, "fig_4_4_delay.png")
    key = "Delta-tau_s = " + ", ".join(f"+{v:.3f}" for v in mid)
    return ("OK", key, ok)


def fig_4_5():
    """sigma_c vs beta: pencere daralmasi + yasa sigma_c ~ 1.15 - 0.6 beta."""
    d = np.load(os.path.join(DATA, "sigmac.npz"))
    bb, sc = d["beta"], d["sigma_c"]
    law = lambda b: 1.15 - 0.6 * np.asarray(b)          # noqa: E731
    solid = bb <= 0.31
    ok = bool(np.all(np.abs(sc[solid] - law(bb[solid])) <= 0.051))

    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ax.errorbar(bb[solid], sc[solid], yerr=0.05, fmt="o", c="C2", ms=8,
                capsize=4, label=r"ölçülen $\sigma_c$ ($\pm 0.05$ tarama adımı)")
    ax.errorbar(bb[~solid], sc[~solid], yerr=0.05, fmt="o", mfc="none",
                c="C2", ms=8, capsize=4,
                label=r"$\beta=0.5$ (sürekleme-işaretli)")
    bl = np.linspace(0, 0.35, 50)
    ax.plot(bl, law(bl), "--", c="k", lw=1.5,
            label=r"$\sigma_c \approx 1{.}15 - 0{.}6\,\beta$")
    ax.plot([0.0], [1.12], "s", c="C3", ms=9, mfc="none", mew=2,
            label=r"Blyth ve Hall (Newtonian) $\sigma_c\approx 1{,}12$")
    ax.axvspan(0.3, 0.52, color="red", alpha=0.06)
    ax.set_xlabel(r"$\beta$ (Deborah sayısı)")
    ax.set_ylabel(r"$\sigma_c$")
    ax.set_title(r"Tekil pencere kenarı daralır: $\sigma_c(\beta)$ ($\Delta=2$)")
    ax.legend(fontsize=8)
    _save(fig, "fig_4_5_window.png")
    key = "sigma_c(0/0.1/0.3) = " + "/".join(f"{v:.2f}" for v in sc[:3])
    return ("OK", key, ok)


def fig_4_6():
    """Patlama ussu korunur: 1/|G_m| vs tau, beta in {0, 0.1, 0.2, 0.3}."""
    with open(os.path.join(DATA, "exponent_results.json")) as fh:
        E = json.load(fh)
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    R2s = {}
    for b in (0.0, 0.1, 0.2, 0.3):
        d = np.load(os.path.join(DATA, f"taus_beta{b:.1f}_em30.npz"))
        tau, gmin, M = d["tau"], d["gmin"], d["M"]
        m = np.isfinite(gmin) & (gmin < -0.5) & np.isfinite(M) & (M >= 1.0)
        s, c, R2 = lin_fit_r2(tau[m], 1.0 / (-gmin[m]))
        R2s[b] = R2
        ax.plot(tau[m], 1.0 / (-gmin[m]), "o", ms=2.6, c=CB[b], alpha=0.85)
        tt = np.linspace(tau[m][0], -c / s, 40)
        ax.plot(tt, s * tt + c, "-", c=CB[b], lw=1.4,
                label=rf"$\beta={b:g}$  ($R^2={R2:.3f}$)")
    ok = all(v > 0.95 for v in R2s.values()) and \
        E["0.0"]["verdict"].startswith("PRESERVED") and \
        E["0.3"]["verdict"].startswith("PRESERVED")
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$1/|G_m|$")
    ax.set_ylim(bottom=0)
    ax.set_title(r"Banks–Zaturska üssü korunur: $G_m\propto(\tau_s-\tau)^{-1}$"
                 "\n" r"$\Rightarrow 1/|G_m|$ doğrusal (aktif pencere $M\geq 1$; "
                 r"$\Delta=2$, $\sigma=0.5$, $\eta_{max}=30$, $N=400$)")
    ax.annotate(r"üs $=-1$, tüm $\beta\leq 0.3$" "\n"
                r"(yalnız önkatsayı kayar: $|G_m|(\tau_s-\tau)$: "
                rf"{E['0.0']['gm_const']:.2f} $\to$ {E['0.3']['gm_const']:.2f})",
                xy=(0.03, 0.72), xycoords="axes fraction", fontsize=9,
                bbox=dict(fc="white", ec="gray", alpha=0.9))
    ax.legend(fontsize=8, loc="lower right")
    _save(fig, "fig_4_6_exponent.png")
    key = "R2(1/|G_m| dogrusalligi) = " + \
        ", ".join(f"{b:g}: {R2s[b]:.3f}" for b in R2s)
    return ("OK", key, bool(ok))


def fig_4_7():
    """Duvar kaymasi f''(0,tau), son periyot: beta=0 vs beta=0.3."""
    t0, f0, a0 = settled_history(0.0)
    t3, f3, a3 = settled_history(0.3)
    A0, ph0, c00 = fourier_amp_phase(t0, f0)
    A3, ph3, c03 = fourier_amp_phase(t3, f3)
    damp = 1.0 - A3 / A0
    ok = abs(damp - 0.28) < 0.02 and max(abs(np.degrees(ph0)),
                                         abs(np.degrees(ph3))) < 5.0

    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ax.plot(t0, f0, "-", c=CB[0.0], lw=2,
            label=r"$\beta=0$ (Newtonian, Blyth ve Hall)")
    ax.plot(t3, f3, "--", c=CB[0.3], lw=2, label=r"$\beta=0.3$ (ikinci derece)")
    ax.set_xlabel(r"$\tau$ mod $2\pi$")
    ax.set_ylabel(r"$f''(0,\tau)$")
    ax.set_xlim(0, PER)
    ax.set_ylim(0.1, 2.45)
    ax.set_title(r"Duvar kayması, son periyot ($\Delta=0.5$, $\sigma=0.5$)")
    ax.annotate("salınım genliği sönümü:\n"
                rf"$1-{A3:.3f}/{A0:.3f} \approx \%{100*damp:.0f}$" "\n"
                rf"faz farkı $|\varphi|\lesssim 5^\circ$ "
                rf"({np.degrees(ph0):+.1f}$^\circ\to${np.degrees(ph3):+.1f}$^\circ$)",
                xy=(0.36, 0.72), xycoords="axes fraction", fontsize=9.5,
                bbox=dict(fc="white", ec="gray", alpha=0.95))
    ax.legend(fontsize=9, loc="lower left")
    _save(fig, "fig_4_7_wall_shear_series.png")
    key = f"sonum = %{100*damp:.1f} (sertifika ~%28); faz {np.degrees(ph3):+.1f} derece"
    return ("OK", key, bool(ok))


def fig_4_8():
    """Hiz profilleri f'(eta) 4 fazda: beta=0 vs beta=0.3 (tabaka kalinlasir)."""
    eta0, fps0, phases = profile_snaps(0.0)
    eta3, fps3, _ = profile_snaps(0.3)
    a_of = [1.5, 1.0, 0.5, 1.0]
    lab = [r"$\tau\equiv 0$ ($a=1.5$)", r"$\tau\equiv\pi/2$ ($a=1$)",
           r"$\tau\equiv\pi$ ($a=0.5$)", r"$\tau\equiv 3\pi/2$ ($a=1$)"]

    _trapz = getattr(np, "trapezoid", np.trapz)

    def dstar(eta, fp, a):
        return float(_trapz(1.0 - fp / a, eta))
    ds0 = dstar(eta0, fps0[0], a_of[0])
    ds3 = dstar(eta3, fps3[0], a_of[0])
    ok = ds3 > ds0

    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for i in range(4):
        ax[0].plot(eta0, fps0[i], "-", c=PHASE_C[i], lw=1.8, label=lab[i])
        ax[1].plot(eta3, fps3[i], "-", c=PHASE_C[i], lw=1.8, label=lab[i])
    for a_ in ax:
        a_.set_xlabel(r"$\eta$")
        a_.set_xlim(0, 8)
    ax[0].set_ylabel(r"$f'(\eta)$")
    ax[0].set_title(r"(a) $\beta=0$")
    ax[1].set_title(r"(b) $\beta=0.3$")
    ax[1].annotate(r"sınır tabakası kalınlaşır:" "\n"
                   rf"$\delta^*(\tau{{\equiv}}0)$: {ds0:.3f} $\to$ {ds3:.3f}"
                   rf"  ($+\%{100*(ds3/ds0-1):.0f}$)",
                   xy=(0.33, 0.08), xycoords="axes fraction", fontsize=9.5,
                   bbox=dict(fc="white", ec="gray", alpha=0.9))
    ax[0].legend(fontsize=8, loc="lower right")
    fig.suptitle(r"Hız profilleri, oturmuş periyot içindeki 4 faz "
                 r"($\Delta=0.5$, $\sigma=0.5$)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    _save(fig, "fig_4_8_profiles.png")
    key = f"delta*(tau=0): {ds0:.3f} -> {ds3:.3f} (+%{100*(ds3/ds0-1):.0f})"
    return ("OK", key, bool(ok))


def fig_4_9():
    """Histerezis dongusu f''(0) vs a(tau); alan ~ genlik^2 (beta-değişmez)."""
    bb = [0.0, 0.1, 0.2, 0.3]
    areas, amps = [], []
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    for b in bb:
        tg, fpp, a = settled_history(b)
        A, ph, c0 = fourier_amp_phase(tg, fpp)
        areas.append(shoelace(a, fpp))
        amps.append(A)
        ax[0].plot(a, fpp, "-", c=CB[b], lw=1.8, label=rf"$\beta={b:g}$")
    areas = np.array(areas)
    amps = np.array(amps)
    norm = areas / amps ** 2
    ok = (np.max(norm) - np.min(norm)) / np.mean(norm) < 0.05 and \
        bool(np.all(np.diff(areas) < 0))
    ax[0].set_xlabel(r"$a(\tau) = 1 + \Delta\cos\tau$")
    ax[0].set_ylabel(r"$f''(0,\tau)$")
    ax[0].set_title(r"(a) $f''(0)$–$a(\tau)$ çevrimi ($\Delta=0.5$, $\sigma=0.5$)")
    ax[0].legend(fontsize=9)
    ax[1].plot(bb, areas / areas[0], "o-", c="C0", lw=1.8,
               label=r"çevrim alanı / $(\beta{=}0)$ alanı")
    ax[1].plot(bb, norm / norm[0], "s--", c="C3", lw=1.8,
               label=r"alan / genlik$^2$, $(\beta{=}0)$'a oranla")
    ax[1].axhline(1.0, c="k", lw=0.7)
    ax[1].set_xlabel(r"$\beta$")
    ax[1].set_ylabel("orana göre normalize")
    ax[1].set_ylim(0.4, 1.15)
    ax[1].set_title(r"(b) Çevrim alanının genlik karesine göre ölçeklenmesi")
    ax[1].annotate(rf"alan / genlik$^2$ sapması $<\%{100*(np.max(norm)/np.min(norm)-1):.1f}$"
                   "\n" r"($\beta$'dan bağımsız)",
                   xy=(0.05, 0.30), xycoords="axes fraction", fontsize=9,
                   bbox=dict(fc="white", ec="gray", alpha=0.9))
    ax[1].legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    _save(fig, "fig_4_9_hysteresis.png")
    key = ("alanlar = " + "/".join(f"{v:.4f}" for v in areas)
           + f"; alan/genlik^2 yayilimi %{100*(np.max(norm)/np.min(norm)-1):.1f}")
    return ("OK", key, bool(ok))


def profile_full(beta):
    """f, f', f'' profilleri tek fazda (tau = 5*2pi, yani tau ≡ 0, a = 1.5).

    fig_4_8 ile AYNI belgelenmiş sertifikalı reçete (figs_thesis.py):
    SA.snapshot_profiles(Delta=0.5, sigma=0.5, beta, eta_max=22, N=200,
    dt=0.005); data/ altına önbelleklenir.
    """
    p = os.path.join(DATA, f"profiles_full_beta{beta:g}.npz")
    if os.path.exists(p):
        d = np.load(p)
        return d["eta"], d["f"], d["fp"], d["fpp"]
    tq = [5 * PER]
    snaps = SA.snapshot_profiles(0.5, 0.5, beta, 22.0, 200, tq, dt=0.005)
    if snaps is None:
        raise RuntimeError(f"snapshot_profiles IC basarisiz (beta={beta})")
    eta, f, fp, fpp = snaps[tq[0]]
    o = np.argsort(eta)
    np.savez(p, eta=eta[o], f=f[o], fp=fp[o], fpp=fpp[o])
    return eta[o], f[o], fp[o], fpp[o]


def fig_4_10():
    """Benzerlik fonksiyonu f ve türevleri f', f'' (tau ≡ 0): beta=0 vs 0.3."""
    e0, F0, FP0, FPP0 = profile_full(0.0)
    e3, F3, FP3, FPP3 = profile_full(0.3)
    w0, w3 = float(FPP0[0]), float(FPP3[0])
    # tutarlılık kilidi: f''(0) periyodik önbelleğin tau ≡ 0 örneğiyle uyuşmalı
    ok = abs(float(F0[0])) < 1e-9 and abs(float(FP0[0])) < 1e-9
    for beta, w in ((0.0, w0), (0.3, w3)):
        d = np.load(os.path.join(DATA, f"periodic_beta{beta:g}.npz"))
        ok = ok and abs(w - float(d["fpp0"][0])) < 0.01
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.9))
    panels = [(F0, F3, r"$f(\eta)$", r"(a) $f$: akım fonksiyonu"),
              (FP0, FP3, r"$f'(\eta)$", r"(b) $f'$: hız"),
              (FPP0, FPP3, r"$f''(\eta)$", r"(c) $f''$: kayma")]
    for a_, (y0, y3, ylab, ttl) in zip(ax, panels):
        a_.plot(e0, y0, "-", c=CB[0.0], lw=2, label=r"$\beta=0$")
        a_.plot(e3, y3, "--", c=CB[0.3], lw=2, label=r"$\beta=0.3$")
        a_.set_xlim(0, 8)
        a_.set_xlabel(r"$\eta$")
        a_.set_ylabel(ylab)
        a_.set_title(ttl, fontsize=10)
    ax[0].annotate(r"$f(0)=0$", xy=(0.07, 0.90), xycoords="axes fraction",
                   fontsize=9, bbox=dict(fc="white", ec="gray", alpha=0.9))
    ax[1].annotate(r"$f'(0)=0$", xy=(0.07, 0.90), xycoords="axes fraction",
                   fontsize=9, bbox=dict(fc="white", ec="gray", alpha=0.9))
    ax[2].annotate(rf"$f''(0)$: {w0:.3f} $\to$ {w3:.3f}" "\n"
                   r"($\beta$ duvar kaymasını düşürür)",
                   xy=(0.28, 0.75), xycoords="axes fraction", fontsize=9,
                   bbox=dict(fc="white", ec="gray", alpha=0.9))
    ax[1].legend(fontsize=9, loc="lower right")
    fig.suptitle(r"Benzerlik fonksiyonu ve türevleri, $\tau\equiv 0$ fazı "
                 r"($a=1.5$; $\Delta=0.5$, $\sigma=0.5$)", fontsize=11.5)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    _save(fig, "fig_4_10_profiles_fdf.png")
    key = f"f''(0): {w0:.3f} -> {w3:.3f} (periyodik önbellekle tutarlı)"
    return ("OK", key, bool(ok))


def dstar_history(beta):
    """delta*(tau) son tam periyot; thesis_figs.py::_dstar_history'nin
    belgelenmiş sertifikalı reçetesi (Delta=0.5, sigma=0.5, 7 periyot,
    eta_max=20, N=200, dt=0.005, '4bc'; import-only bh.steady_ic +
    bh.newton_step). data/ altına önbelleklenir; dönüş: (tg, ds, a, regen)."""
    p = os.path.join(DATA, f"dstar_beta{beta:g}.npz")
    if os.path.exists(p):
        d = np.load(p)
        return d["tg"], d["ds"], d["a"], False
    _trapz = getattr(np, "trapezoid", np.trapz)
    Delta, sigma, n_periods, eta_max, N, dt = 0.5, 0.5, 7, 20.0, 200, 0.005
    A = C.A_of_sigma(sigma)
    invA = 1.0 / A
    a0 = bh.a_outer(0.0, Delta, 1.0, "bh")
    U, ops = bh.steady_ic(beta, 0.0, a0, N, eta_max)
    if U is None:
        raise RuntimeError(f"steady IC basarisiz (beta={beta})")
    eta, D1 = ops["eta"], ops["D1"]
    o = np.argsort(eta)
    es = eta[o]
    Uold = U.copy()
    rec_t, rec_d, rec_a = [], [], []
    nstep = int(round(n_periods * PER / dt))
    for n in range(nstep + 1):
        tau = n * dt
        if n > 0:
            U, okk, _ = bh.newton_step(U, ops, beta, dt, Uold, tau, 0.0, Delta,
                                       1.0, invA, "bh", "4bc")
            if not okk:
                raise RuntimeError(f"newton_step basarisiz (beta={beta}, tau={tau})")
            Uold = U.copy()
        a = bh.a_outer(tau, Delta, 1.0, "bh")
        fp = D1 @ U
        rec_t.append(tau)
        rec_d.append(float(_trapz(1.0 - fp[o] / a, es)))
        rec_a.append(a)
    t = np.array(rec_t)
    ds = np.array(rec_d)
    aa = np.array(rec_a)
    nP = int(t[-1] // PER)
    m = (t >= (nP - 1) * PER) & (t < nP * PER)
    tg = t[m] - (nP - 1) * PER
    np.savez(p, tg=tg, ds=ds[m], a=aa[m])
    return tg, ds[m], aa[m], True


def fig_4_11():
    """Yer-değiştirme kalınlığı delta*(tau), beta in {0, 0.1, 0.2, 0.3}."""
    bb = [0.0, 0.1, 0.2, 0.3]
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    means, lock = {}, {}
    regen = []
    for b in bb:
        tg, ds, a, was_regen = dstar_history(b)
        if was_regen:
            regen.append(b)
        means[b] = float(np.mean(ds))
        lock[b] = float(ds[0])                       # tau ≡ 0 örneği
        ax.plot(tg, ds, "-", c=CB[b], lw=1.9,
                label=rf"$\beta={b:g}$  (ort. $\overline{{\delta^*}}={means[b]:.3f}$)")
    # kilit: fig_4_8'in sertifikalı delta*(tau≡0) notuyla uyuşmalı
    ok_lock = abs(lock[0.0] - 0.529) < 2e-3 and abs(lock[0.3] - 0.769) < 2e-3
    ok_order = means[0.0] < means[0.1] < means[0.2] < means[0.3]
    ax.set_xlim(0, PER)
    ax.set_xlabel(r"$\tau$ mod $2\pi$")
    ax.set_ylabel(r"$\delta^*(\tau)$")
    ax.set_title(r"Yer-değiştirme kalınlığı "
                 r"$\delta^*(\tau)=\int_0^{\eta_{max}}(1-f'/a)\,d\eta$"
                 "\n" r"($\Delta=0.5$, $\sigma=0.5$): $\beta$ tabakayı tüm "
                 r"salınım boyunca kalınlaştırır")
    ax.legend(fontsize=9)
    _save(fig, "fig_4_11_delta_star_beta.png")
    key = (f"delta*(tau=0): {lock[0.0]:.3f} / {lock[0.3]:.3f} "
           f"(kilit 0.529 / 0.769); ort. "
           + "/".join(f"{means[b]:.3f}" for b in bb)
           + (f"; yeniden uretilen beta: {[f'{b:g}' for b in regen]}" if regen else ""))
    return ("OK", key, bool(ok_lock and ok_order))


FIGURES = {
    "fig_4_1": (fig_4_1, "fig_4_1_bh_validation.png"),
    "fig_4_2": (fig_4_2, "beta0_slope.png"),
    "fig_4_3": (fig_4_3, "bh_anchor_loop.png"),
    "fig_4_4": (fig_4_4, "fig_4_4_delay.png"),
    "fig_4_5": (fig_4_5, "fig_4_5_window.png"),
    "fig_4_6": (fig_4_6, "fig_4_6_exponent.png"),
    "fig_4_7": (fig_4_7, "fig_4_7_wall_shear_series.png"),
    "fig_4_8": (fig_4_8, "fig_4_8_profiles.png"),
    "fig_4_9": (fig_4_9, "fig_4_9_hysteresis.png"),
    "fig_4_10": (fig_4_10, "fig_4_10_profiles_fdf.png"),
    "fig_4_11": (fig_4_11, "fig_4_11_delta_star_beta.png"),
}


def main():
    ap = argparse.ArgumentParser(description="tez figurleri (tek giris noktasi)")
    ap.add_argument("--only", help="yalniz bu figur (or. fig_4_4)")
    ap.add_argument("--list", action="store_true", help="figur listesini yaz")
    args = ap.parse_args()
    if args.list:
        for k, (_, fn) in FIGURES.items():
            print(f"  {k} -> {fn}")
        return
    todo = list(FIGURES)
    if args.only:
        todo = ([args.only] if args.only in FIGURES
                else [k for k in FIGURES if args.only in k])
        if not todo:
            sys.exit(f"bilinmeyen figur: {args.only}")
    rows = []
    for k in todo:
        fn, outname = FIGURES[k]
        try:
            status, key, match = fn()
        except Exception as e:  # dosya eksik vb. -> BLOCKED, devam et
            status, key, match = "BLOCKED", f"{type(e).__name__}: {e}", None
        rows.append((outname, status, key,
                     {True: "EVET", False: "HAYIR", None: "-"}[match]))
        print(f"[{k}] {status}: {key}")
    w = max(len(r[0]) for r in rows)
    print("\n" + "=" * 100)
    print(f"{'figur':<{w}} | durum   | anahtar sayi | sertifikaya uygun")
    print("-" * 100)
    for fn_, st, key, mt in rows:
        print(f"{fn_:<{w}} | {st:<7} | {key} | {mt}")
    bad = [r for r in rows if r[1] != "OK" or r[3] == "HAYIR"]
    print("=" * 100)
    print("SONUC:", "TUM FIGURLER SERTIFIKALI DEGERLERLE UYUMLU"
          if not bad else f"{len(bad)} figur sorunlu — yukaridaki tabloya bakin")


if __name__ == "__main__":
    main()
