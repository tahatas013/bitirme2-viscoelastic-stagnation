# PROVENANCE — dosya kökenleri

Bu depodaki tüm çözücü dosyaları, kaynak çalışma alanından **bayt-özdeş (verbatim)**
kopyalanmıştır; hiçbirinde tek karakter dahi değiştirilmemiştir (köken bilgisi bu
dosyada tutulur, dosyaların içinde değil). Yollar, kaynak çalışma alanı köküne
(`Sayısal Çözüm/`) görelidir.

## Çözücüler (verbatim kopya)

| Depo dosyası | Kaynak yol | Rol |
|---|---|---|
| `solvers/bh_validation/bh_solver.py` | `viscoelastic_stagnation/bh_validation/bh_solver.py` | Kanonik B&H Newtonian zamana-bağlı RFNC marşı (sertifikalı: durağan f''(0)=1.232588; tekil durum τ_s ∈ [3.38, 3.41]). Kullanılmayan eski alternatif dış-akış seçenekleri içerir; bu tezde her yerde `outer='bh'`. |
| `solvers/rfnc/rfnc_steady.py` | `viscoelastic_stagnation/OLDaFunc/residual_form_newton/rfnc_steady.py` | Doğrulanmış spektral durağan RFNC çekirdeği (Chebyshev kolokasyon + analitik-Jacobian Newton; durağan başlangıç koşulu çözücüsü). |
| `solvers/rfnc/rfnc_unsteady.py` | `viscoelastic_stagnation/OLDaFunc/residual_form_newton/rfnc_unsteady.py` | Doğrulanmış zamana-bağlı RFNC referans çözücüsü (Rothe). |
| `solvers/secondgrade/solver_a_sg.py` | `newt_sol/secondgrade_bh/solver_a_sg.py` | Çözücü A (BİRİNCİL, tam/kesin β): doğrulanmış `bh_solver` adımını import-only süren marş; tez katkısının ana motoru. |
| `solvers/secondgrade/config_sg.py` | `newt_sol/secondgrade_bh/config_sg.py` | Aşama 2b ortak konfigürasyon: temiz ölçek (ω=1, A=1/σ ⇒ invA=σ), çapa değerleri, sertifika sabitleri. |
| `solvers/secondgrade/taus.py` | `newt_sol/secondgrade_bh/taus.py` | τ_s tahmincileri: 1/|G_m| (üst) + η_m⁻² (alt), aktif-pencere fitleri. |
| `solvers/beta_perturbation/solver_b_sg.py` | `newt_sol/secondgrade_bh/solver_b_sg.py` | Çözücü B (BAĞIMSIZ, koşullanma-serbest β-pertürbasyonu, f = f₀ + βf₁), B&H dış akışıyla; iki-çözücü sertifika bacağı. |
| `solvers/beta_perturbation/pert_steady.py` | `viscoelastic_stagnation/OLDaFunc/perturbation_unsteady/pert_steady.py` | Sertifikalı pertürbasyon motoru (integral yeniden-formülasyon, w = f‴; ~N² koşullanma). |
| `solvers/beta_perturbation/ultraspherical_ops.py` | `viscoelastic_stagnation/OLDaFunc/ultraspherical/ultraspherical_ops.py` | Chebyshev değer↔katsayı dönüşüm operatörleri (motor bağımlılığı). |
| `solvers/beta_perturbation/integral_reform_steady.py` | `viscoelastic_stagnation/OLDaFunc/ultraspherical/integral_reform_steady.py` | İntegral-reform integrasyon matrisi (motor bağımlılığı). |

## Sertifikasyon belgeleri (verbatim kopya)

| Depo dosyası | Kaynak yol | Rol |
|---|---|---|
| `certification/beta0_convergence.md` | `newt_sol/beta0_indep/outputs/beta0_convergence.md` | Bağımsız β→0 yakınsama sertifikası: Çözücü A, B&H'ye O(β) hızında yakınsar (3/3 kabul kriteri GEÇTİ). |
| `certification/RESULTS_certify.md` | `newt_sol/secondgrade_bh/RESULTS_certify.md` | Üç-bacaklı sertifikasyon: β→0 indirgeme; iki bağımsız çözücünün uyumu; koşullanmanın zararsızlığı. |
| `certification/RESULTS_taus.md` | `newt_sol/secondgrade_bh/RESULTS_taus.md` | Çözünürlük-yakınsak Δτ_s yeniden hesabı (hücre ara-nokta tahmincisi; F1'in kaynağı). |
| `certification/RESULTS_sigmac.md` | `newt_sol/secondgrade_bh/RESULTS_sigmac.md` | Tekil pencere kenarı σ_c(β) taraması (F2'nin kaynağı). |
| `certification/RESULTS_exponent.md` | `newt_sol/secondgrade_bh/RESULTS_exponent.md` | Banks–Zaturska patlama üssünün korunumu (F3'ün kaynağı). |

## Önbellek verisi (verbatim kopya)

| Depo dosyası | Kaynak yol | Rol / koşu parametreleri |
|---|---|---|
| `data/stage2_singular.npz` | `viscoelastic_stagnation/bh_validation/data/stage2_singular.npz` | Kanonik Newtonian tekil marş (Δ=2, σ=0.5, N=240, η_max=24, dt=0.0025) — Şekil 4.1b. |
| `data/taus_beta{0.0,0.1,0.2,0.3}_em30.npz` | `newt_sol/secondgrade_bh/data/` (aynı adlar) | Sertifikalı tekil marşlar, referans hücre (Δ=2, σ=0.5, η_max=30, N=400, dt=0.0025; Çözücü A) — Şekil 4.6. |
| `data/taus_results.json` | `newt_sol/secondgrade_bh/data/taus_results.json` | Sertifikalı hücre ara-nokta Δτ_s değerleri + tahminci zarfı — Şekil 4.4. |
| `data/sigmac.npz`, `data/sigmac_results.json` | `newt_sol/secondgrade_bh/data/` | σ_c(β) tarama sonuçları (sınıflandırıcı: N=120, dt=0.0025→patlamada 0.00125, 8 periyot, η_max=16, ±0.05 tarama adımı) — Şekil 4.5. |
| `data/exponent_results.json` | `newt_sol/secondgrade_bh/data/exponent_results.json` | Üs-korunum skalerları (R², |G_m|(τ_s−τ) sabiti) — Şekil 4.6 notu. |

## Yeniden üretilen veri (önbellekte YOKTU; belgelenmiş sertifikalı parametrelerle)

Aşağıdaki dosyalar `scripts/make_figures.py` tarafından, kaynak paketin BELGELENMİŞ
parametreleriyle, yalnızca import edilen sertifikalı çözücüler çalıştırılarak üretildi
ve tekrarlanabilirlik için `data/` altına kaydedildi. İndirgeme: tam marş dizileri
yerine yalnız figürlerin kullandığı SON TAM PERİYOT saklanır.

| Depo dosyası | Üretim reçetesi (kaynak belge) | Rol |
|---|---|---|
| `data/periodic_beta{0,0.1,0.2,0.3}.npz` | `solver_a_sg.march(Δ=0.5, σ=0.5, β, η_max=20, N=200, 7 periyot, dt=0.005)`, son tam periyot (`newt_sol/secondgrade_bh/thesis_figs.py::_settled_history` parametreleri) | Şekil 4.7 ve 4.9. Sertifikalı kontrol: β=0.3 genlik sönümü %28.1 (sertifika ≈%28); döngü alanları 0.1142/0.0878/0.0708/0.0588 (sertifikalı günlük değerleriyle 4 ondalık özdeş). |
| `data/profiles_beta{0,0.3}.npz` | `solver_a_sg.snapshot_profiles(Δ=0.5, σ=0.5, β, η_max=22, N=200, dt=0.005)`, fazlar τ = 5·2π + {0, π/2, π, 3π/2} (`newt_sol/secondgrade_bh/figs_thesis.py` parametreleri) | Şekil 4.8. |

Şekil 4.1a'daki durağan Hiemenz profili önbelleğe alınmaz; doğrulama çalışmasının
belgelenmiş çağrısıyla anlık hesaplanır: `rfnc_steady.solve_steady(0, 0, N=120,
η_max=12, a=1)` (`viscoelastic_stagnation/bh_validation/run_validation.py` ile aynı).

## β→0 doğrulama önbelleği (verbatim kopya) ve Şekil 4.2–4.3

`data/beta0/` — `newt_sol/beta0_indep/cache/` içinden bayt-özdeş kopyalar (23
dosya): `A_N{180,360}_b{0, 0.012, 0.018, 0.025, 0.035, 0.05, 0.07, 0.1, 0.14,
0.2}.npz` (Çözücü A periyodik marşlar, Δ=0.5, σ=0.5, η_max=18, dt=2π/1256,
6 periyot), `B_pert.npz` (pertürbasyon çözücüsü, f₁''(0,τ)),
`BH_3bc_N180.npz` + `BH_4bc_N180.npz` (Newtonian referans marşlar).

`figures/beta0_slope.png` (Şekil 4.2) ve `figures/bh_anchor_loop.png`
(Şekil 4.3) bu önbellekten `scripts/make_figures.py` (fig_4_2 / fig_4_3)
tarafından Türkçe metinle üretilir; gösterilen tüm sayılar
`certification/beta0_convergence.md` içindeki doğrulanmış değerlerin
birebir aynısıdır (asimptotik eğim 0.9660; C_A = C_B = 3.1357; çevrim
üzerinde en büyük bağıl fark 3.3e-8). İngilizce orijinal görseller kaynak
çalışma alanında (`newt_sol/beta0_indep/outputs/`) durur; doğrulama belgesi
`certification/` altında değiştirilmeden korunur.
