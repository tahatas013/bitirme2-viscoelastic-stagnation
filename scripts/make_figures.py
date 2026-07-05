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
    ax[0].plot(fp[o], eta[o], "-", c="k", lw=2)
    ax[0].axvline(1.0, ls=":", c="gray", lw=1)
    ax[0].set_xlabel(r"$f'(\eta)$")
    ax[0].set_ylabel(r"$\eta$")
    ax[0].set_ylim(0, 6)
    ax[0].set_title(r"(a) Durağan Hiemenz profili ($\beta=0$, $a=1$)")
    ax[0].annotate(r"$f''(0)=1{.}232588$" "\n(kanonik çapa)",
                   xy=(0.35, 1.2), fontsize=10,
                   bbox=dict(fc="white", ec="gray", alpha=0.9))
    ax[1].plot(tau[m], 1.0 / (-gmin[m]), "o", ms=3, c="k",
               label=r"$1/|G_m|$ (aktif pencere $M\geq 1$)")
    tt = np.linspace(tau[m][0], 3.41, 50)
    ax[1].plot(tt, s * tt + c, "-", c="C3", lw=1.5, label="doğrusal fit")
    ax[1].axvspan(3.38, 3.41, color="C2", alpha=0.25,
                  label=r"sertifikalı bant $\tau_s\in[3{.}38,\,3{.}41]$")
    ax[1].set_xlabel(r"$\tau$")
    ax[1].set_ylabel(r"$1/|G_m|$")
    ax[1].set_ylim(bottom=0)
    ax[1].set_title(r"(b) Tekil durum $\Delta=2$, $\sigma=0.5$: $\tau_s$'e doğrusal yaklaşım")
    ax[1].legend(fontsize=8)
    fig.suptitle(r"B&H Newtonian doğrulaması ($\beta=0$)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    _save(fig, "fig_4_1_bh_validation.png")
    key = f"f''(0)={fpp0:.6f}; 1/|G_m| fit kesisimi bant icinde: {'EVET' if ok_b else 'HAYIR'} (R2={R2:.3f})"
    return ("OK", key, ok_a and ok_b)


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
    ax.set_title(r"Viskoelastisite tekilliği GECİKTİRİR: $\Delta\tau_s(\beta)>0$"
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
            label=r"B&H (Newtonian) $\sigma_c\approx 1{.}12$")
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
    ax.set_title(r"Banks-Zaturska üssü KORUNUR: $G_m\propto(\tau_s-\tau)^{-1}$"
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
    ax.plot(t0, f0, "-", c=CB[0.0], lw=2, label=r"$\beta=0$ (Newtonian, B&H)")
    ax.plot(t3, f3, "--", c=CB[0.3], lw=2, label=r"$\beta=0.3$ (ikinci derece)")
    ax.set_xlabel(r"$\tau$ mod $2\pi$")
    ax.set_ylabel(r"$f''(0,\tau)$")
    ax.set_xlim(0, PER)
    ax.set_title(r"Duvar kayması, oturmuş periyot ($\Delta=0.5$, $\sigma=0.5$)")
    ax.annotate("salınım genliği sönümü:\n"
                rf"$1-{A3:.3f}/{A0:.3f} \approx \%{100*damp:.0f}$" "\n"
                rf"faz farkı $|\varphi|\lesssim 5^\circ$ "
                rf"({np.degrees(ph0):+.1f}$^\circ\to${np.degrees(ph3):+.1f}$^\circ$)",
                xy=(0.35, 0.04), xycoords="axes fraction", fontsize=9.5,
                bbox=dict(fc="white", ec="gray", alpha=0.9))
    ax.legend(fontsize=9, loc="upper right")
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
        ax[0].plot(fps0[i], eta0, "-", c=PHASE_C[i], lw=1.8, label=lab[i])
        ax[1].plot(fps3[i], eta3, "-", c=PHASE_C[i], lw=1.8, label=lab[i])
    for a_ in ax:
        a_.set_xlabel(r"$f'(\eta)$")
        a_.set_ylim(0, 8)
    ax[0].set_ylabel(r"$\eta$")
    ax[0].set_title(r"(a) $\beta=0$")
    ax[1].set_title(r"(b) $\beta=0.3$")
    ax[1].annotate(r"sınır tabakası KALINLAŞIR:" "\n"
                   rf"$\delta^*(\tau{{\equiv}}0)$: {ds0:.3f} $\to$ {ds3:.3f}"
                   rf"  ($+\%{100*(ds3/ds0-1):.0f}$)",
                   xy=(0.30, 0.80), xycoords="axes fraction", fontsize=9.5,
                   bbox=dict(fc="white", ec="gray", alpha=0.9))
    ax[0].legend(fontsize=8, loc="upper left")
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
    ax[0].set_title(r"(a) Elastik hafıza döngüsü ($\Delta=0.5$, $\sigma=0.5$)")
    ax[0].legend(fontsize=9)
    ax[1].plot(bb, areas / areas[0], "o-", c="C0", lw=1.8,
               label=r"döngü alanı / alan$(\beta{=}0)$")
    ax[1].plot(bb, norm / norm[0], "s--", c="C3", lw=1.8,
               label=r"(alan/genlik$^2$) / aynı oran$(\beta{=}0)$")
    ax[1].axhline(1.0, c="k", lw=0.7)
    ax[1].set_xlabel(r"$\beta$")
    ax[1].set_ylabel("orana göre normalize")
    ax[1].set_ylim(0.4, 1.15)
    ax[1].set_title(r"(b) Döngü KÜÇÜLÜR ama alan/genlik$^2$ $\beta$-DEĞİŞMEZ")
    ax[1].annotate(rf"alan/genlik$^2$ sapması $<\%{100*(np.max(norm)/np.min(norm)-1):.1f}$",
                   xy=(0.05, 0.32), xycoords="axes fraction", fontsize=9,
                   bbox=dict(fc="white", ec="gray", alpha=0.9))
    ax[1].legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    _save(fig, "fig_4_9_hysteresis.png")
    key = ("alanlar = " + "/".join(f"{v:.4f}" for v in areas)
           + f"; alan/genlik^2 yayilimi %{100*(np.max(norm)/np.min(norm)-1):.1f}")
    return ("OK", key, bool(ok))


FIGURES = {
    "fig_4_1": (fig_4_1, "fig_4_1_bh_validation.png"),
    "fig_4_4": (fig_4_4, "fig_4_4_delay.png"),
    "fig_4_5": (fig_4_5, "fig_4_5_window.png"),
    "fig_4_6": (fig_4_6, "fig_4_6_exponent.png"),
    "fig_4_7": (fig_4_7, "fig_4_7_wall_shear_series.png"),
    "fig_4_8": (fig_4_8, "fig_4_8_profiles.png"),
    "fig_4_9": (fig_4_9, "fig_4_9_hysteresis.png"),
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
        todo = [k for k in FIGURES if args.only in k]
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
