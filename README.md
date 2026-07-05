# İkinci Derece Viskoelastik Akışkan İçin Zamana Bağlı Durma-Noktası Akışı

**Bitirme Projesi II — sayısal çözüm deposu (kod + veri + figürler)**

> *English abstract.* We study the unsteady two-dimensional stagnation-point
> boundary-layer flow of a second-grade (viscoelastic) fluid driven by the
> oscillatory outer flow a(τ) = 1 + Δ cos τ of Blyth & Hall (2003), which is
> recovered exactly at zero Deborah number (β = 0). A certified residual-form
> Newton–Chebyshev collocation solver with Rothe time stepping (plus an
> independent β-perturbation solver) shows that elasticity **delays** the
> finite-time singularity (Δτ_s = +0.022 / +0.052 / +0.087 for β = 0.1 / 0.2 /
> 0.3), **narrows** the singular window (σ_c ≈ 1.15 − 0.6β), and **preserves**
> the Banks–Zaturska blow-up exponent, while in the periodic regime the
> wall-shear oscillation is damped by ≈28% at β = 0.3 with a phase shift
> below 5°.

Bu depo, tezin tüm sayısal sonuçlarını üreten **sertifikalı** çözücüleri
(bayt-özdeş kopyalar, bkz. [PROVENANCE.md](PROVENANCE.md)), sertifikasyon
belgelerini ve tüm tez figürlerini tek komutla yeniden üreten betiği içerir.

## 1. Problem

Benzerlik değişkenleriyle indirgenmiş yönetici denklem (denk. 39; üsler
η-türevleri) ve dört sınır koşulu:

```
(1/A)·f'_τ + f'² − f·f'' = (1/A)·a_τ + a² + f''' + β[(1/A)·f'''_τ + 2f'f''' − f·f'''' − f''²]

f'(0,τ) = 0            (kaymazlık)
f(0,τ)  = 0            (geçirimsiz duvar)
f'(∞,τ) = a(τ)         (dış akışla eşleşme)
f''(∞,τ) = 0           (Garg–Rajagopal kapanışı; 4. mertebe için gerekli tek ek koşul)

a(τ) = 1 + Δ cos τ     (Blyth & Hall salınımlı dış akışı)
```

- **β (Deborah sayısı):** ikinci derece akışkanın elastikliğini ölçer; model
  küçük-Deborah yaklaşımıdır (β ≤ 0.3 sağlam; β = 0.5 sürekleme-işaretli).
  **β = 0 tam olarak Newtonian Blyth & Hall (2003) problemine indirgenir.**
- **σ = ω/A (Strouhal sayısı):** salınım frekansının akış şiddetine oranı;
  tek fiziksel parametredir. Temiz ölçek: ω = 1, A = 1/σ.
- Tekil çapa durumu: Δ = 2, σ = 0.5 (dış akış tersine döner, sonlu-zaman
  tekilliği); periyodik durum: Δ = 0.5, σ = 0.5 (çekici periyodik yörünge).

## 2. Sonuçlar

| # | Bulgu | Sayısal sonuç |
|---|---|---|
| F1 | **Elastiklik tekilliği geciktirir** | Δτ_s = τ_s(β) − τ_s(0) = **+0.022 / +0.052 / +0.087** (β = 0.1 / 0.2 / 0.3); η_max-değişmez hücre ara-nokta tahmincisi, 1/\|G_m\| (üst) + η_m⁻² (alt) zarfıyla. |
| F2 | **Tekil pencere daralır** | σ_c ≈ **1.15 − 0.6 β** (ölçülen: 1.15 / 1.05 / 0.95 @ β = 0 / 0.1 / 0.3, ±0.05); β = 0, B&H'nin ≈1.12 değerini kurtarır. |
| F3 | **Patlama üssü korunur** | Banks–Zaturska G_m ∝ (τ_s − τ)⁻¹ üssü tüm β ≤ 0.3 için değişmez (1/\|G_m\| doğrusallığı R² ≥ 0.96; yalnız önkatsayı 0.55 → 0.68 kayar). |

Periyodik rejimde ayrıca: duvar-kayması salınım genliği β = 0.3'te ≈%28
sönümlenir, sınır tabakası kalınlaşır, f''(0)–a(τ) histerezis döngüsü küçülür
ama alan/genlik² oranı β-değişmezdir; faz farkı ≲5°.

τ_s disiplini: mutlak τ_s tek sayı olarak verilmez; Newtonian doğrulamada bant
τ_s ∈ [3.38, 3.41] (B&H raporu 3.39), β etkisi daima fark Δτ_s olarak rapor edilir.

## 3. Doğrulama (üç ayaklı onaylama)

1. **β → 0, B&H'yi kurtarır.** Durağan çapa f''(0) = 1.232588 (6 hane);
   tekil durumda τ_s bandı [3.38, 3.41] ∋ 3.39. Ayrıca Çözücü A'nın β → 0
   limitine **O(β) hızında** yakınsadığı, Çözücü B'den bağımsız olarak
   sertifikalandı: [certification/beta0_convergence.md](certification/beta0_convergence.md)
   (3/3 kabul kriteri; asimptotik eğim 0.966, önfaktör uyumu %2.8).
2. **Yapısal olarak bağımsız iki çözücü uyuşur.** Çözücü A (tam doğrusal-olmayan
   rezidüel, 4. mertebe) ile Çözücü B (koşullanma-serbest β-pertürbasyonu,
   f = f₀ + βf₁) ortak geçerlilik penceresinde aynı duvar kaymasını verir:
   [certification/RESULTS_certify.md](certification/RESULTS_certify.md).
3. **Patlama üssü korunur.** Tekil asimptotik yapı (Banks–Zaturska) β ile
   değişmez: [certification/RESULTS_exponent.md](certification/RESULTS_exponent.md).

## 4. Depo yapısı

```
bitirme2-viscoelastic-stagnation/
├── README.md               # bu dosya
├── PROVENANCE.md           # her kopyalanan dosyanın kökeni (bayt-özdeş) + veri reçeteleri
├── solvers/
│   ├── bh_validation/      # kanonik B&H Newtonian çözücü (verbatim)
│   ├── secondgrade/        # sertifikalı Aşama 2b çözücüsü: Çözücü A + konfig + τ_s tahmincileri
│   ├── rfnc/               # doğrulanmış spektral RFNC çekirdeği (durağan + zamana bağlı)
│   └── beta_perturbation/  # Çözücü B (β-pertürbasyonu) + motor bağımlılıkları
├── certification/          # sertifikasyon/sonuç belgeleri (verbatim)
├── data/                   # figürlerin kullandığı önbellek verisi (npz/json)
├── scripts/
│   └── make_figures.py     # TÜM figürleri üreten tek giriş noktası
└── figures/                # çıktılar (300 dpi PNG)
```

## 5. Yeniden üretim

```bash
pip install numpy scipy matplotlib
python scripts/make_figures.py            # tüm figürler (depo kökünden)
python scripts/make_figures.py --only fig_4_4   # tek figür
```

- Kullanılan Python sürümü: **3.9.6** (numpy 2.0.2, matplotlib 3.9.4 ile test edildi).
- Betik, önbellekteki sertifikalı veriyi kullanır; önbellekte olmayan periyodik
  koşuları belgelenmiş sertifikalı parametrelerle (bkz. PROVENANCE.md) yeniden
  üretir ve sonunda her figürün anahtar sayısını sertifikalı değerle
  karşılaştıran bir özet tablo basar.

## 6. Figür dizini

| Dosya | Tez | Açıklama |
|---|---|---|
| `figures/fig_4_1_bh_validation.png` | Şekil 4.1 | Blyth ve Hall doğrulaması: durağan Hiemenz profili (f''(0) = 1.232588) + tekil durumda 1/\|G_m\| doğrusal yaklaşımı, bant [3.38, 3.41]. |
| `figures/beta0_slope.png` | Şekil 4.2 | β → 0 bağımsız yakınsama: log-log O(β) eğimi + öncül-katsayı paneli (β→0 doğrulama önbelleğinden, `data/beta0/`, Türkçe üretim). |
| `figures/bh_anchor_loop.png` | Şekil 4.3 | f''(0)–a(τ) çevrimi: viskoelastik çözücü (β=0) ile Newtonian çözücünün örtüşmesi (β→0 doğrulama önbelleğinden, `data/beta0/`, Türkçe üretim). |
| `figures/fig_4_4_delay.png` | Şekil 4.4 | Gecikme: Δτ_s(β) = +0.022 / +0.052 / +0.087, tahminci zarfıyla. |
| `figures/fig_4_5_window.png` | Şekil 4.5 | Pencere daralması: σ_c(β) noktaları + σ_c ≈ 1.15 − 0.6β yasası + B&H ≈1.12. |
| `figures/fig_4_6_exponent.png` | Şekil 4.6 | Üs korunumu: 1/\|G_m\| vs τ, β ∈ {0, 0.1, 0.2, 0.3}, hepsi doğrusal (üs −1). |
| `figures/fig_4_7_wall_shear_series.png` | Şekil 4.7 | Periyodik duvar kayması f''(0,τ): β = 0 vs 0.3; genlik sönümü ≈%28, faz ≲5°. |
| `figures/fig_4_8_profiles.png` | Şekil 4.8 | Hız profilleri f'(η), 4 faz: sınır tabakası β ile kalınlaşır. |
| `figures/fig_4_9_hysteresis.png` | Şekil 4.9 | Histerezis döngüleri + alan/genlik² β-değişmezliği. |
| `figures/fig_4_10_profiles_fdf.png` | Şekil 4.10 | Benzerlik fonksiyonu ve türevleri f, f′, f″ (τ≡0 fazı, a=1.5): β, duvar kayması f″(0)'ı düşürür, tabakayı kalınlaştırır. |
| `figures/fig_4_11_delta_star_beta.png` | Şekil 4.11 | Yer-değiştirme kalınlığı δ*(τ), β ∈ {0, 0.1, 0.2, 0.3}: viskoelastisite tabakayı tüm salınım boyunca kalınlaştırır. |

## 7. Kaynaklar

- Blyth, M. G. & Hall, P. (2003). Oscillatory flow near a stagnation point.
  *SIAM Journal on Applied Mathematics*, 63(5), 1604–1614.
- Garg, V. K. & Rajagopal, K. R. (1990). Stagnation point flow of a
  non-Newtonian fluid. *Mechanics Research Communications*, 17(6), 415–421.
- Hiemenz, K. (1911). Die Grenzschicht an einem in den gleichförmigen
  Flüssigkeitsstrom eingetauchten geraden Kreiszylinder. *Dinglers Polytech.
  J.*, 326, 321–324.
- Rivlin, R. S. & Ericksen, J. L. (1955). Stress-deformation relations for
  isotropic materials. *Journal of Rational Mechanics and Analysis*, 4, 323–425.
- Beard, D. W. & Walters, K. (1964). Elastico-viscous boundary-layer flows.
  *Mathematical Proceedings of the Cambridge Philosophical Society*, 60(3), 667–674.
