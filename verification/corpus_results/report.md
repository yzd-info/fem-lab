# FEM Lab — 80 vakalık girdi korpusu: yürütme raporu

scikit-fem 11.0.0 · korpus: `verification/corpus/cases.json` (yalnız girdi, değişmedi) · sonuçlar: `results.json`, `skfem_results.json`, `continuous_results.json` · taban: `baseline_v1/`

## Kapsam

- **Kategori:** derivative_flux_bc: 10, element_order_mesh: 10, first_order_ivp: 10, invalid_adversarial: 10, second_order_bvp: 15, singular_degenerate_compat: 7, trig_exp_source: 8, variable_coefficient: 10
- **Geçerlilik:** challenging: 12, invalid: 10, valid: 58
- **Mertebe:** 1: 13, 2: 67
- **Derece p:** 1: 32, 2: 28, 3: 11, 4: 8, 5: 1
- **n:** 0: 1, 1: 2, 10: 1, 16: 1, 2: 3, 3: 6, 32: 1, 4: 37, 5: 8, 6: 8, 64: 3, 7: 2, 8: 6, 9: 1
- **Koşul türleri (uç başına):** both: 3, dirichlet: 121, flux: 4, neumann: 11, none: 21

## Yürütme özeti

| Durum | Vaka |
|---|---:|
| PASS | 60 |
| WARN | 1 |
| LIMITATION | 3 |
| DIAGNOSED | 6 |
| REJECTED_OK | 10 |
| FAIL | 0 |
| RUNNER_ERROR | 0 |

## Çapraz doğrulama matrisi (geçerli + zorlayıcı vakalar)

Hücre değerleri: PASS · FAIL · WARN · UNSUPPORTED · INCONCLUSIVE · NOT_RUN. INCONCLUSIVE / UNSUPPORTED / NOT_RUN hiçbir zaman doğrulanmış sayılmaz.

| ID | Durum | tarayıcı ↔ sunucu | FEM Lab ↔ ifem | scikit-fem eleman | scikit-fem global K | scikit-fem düğüm | analitik çözüm | sürekli varlık | sürekli teklik | uyumluluk | L²/H¹ yeniden hesap | sınır koşulu |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C01 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C02 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C03 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C04 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C05 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C06 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C07 | LIMITATION | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C08 | LIMITATION | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C09 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C10 | WARN | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C11 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C12 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C13 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C14 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C15 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C16 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C17 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C18 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C19 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C20 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C21 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C22 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C23 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C24 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C25 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C26 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C27 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ? | ? | ? | ? | – | ✓ |
| C28 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C29 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ? | ? | ? | ? | – | ✓ |
| C30 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ? | ? | ? | ? | – | ✓ |
| C31 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C32 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ? | ? | ? | ? | – | ✓ |
| C33 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ? | ? | ✓ | ✓ |
| C34 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ? | ? | ? | ? | – | ✓ |
| C35 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ? | ? | ? | ? | – | ✓ |
| C36 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C37 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C38 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C39 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C40 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C41 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C42 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | – |
| C43 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C44 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C45 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C46 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C47 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C48 | PASS | ✓ | ✓ | ✓ | – | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C49 | PASS | ✓ | ✓ | ✓ | – | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C50 | PASS | ✓ | – | ✓ | – | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C51 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C52 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C53 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C54 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C55 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C56 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C57 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C58 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C59 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C60 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C61 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C62 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C63 | PASS | ✓ | ✓ | ✓ | ✓ | ✓ | ? | ? | ? | ? | – | ✓ |
| C64 | DIAGNOSED | – | – | – | – | – | – | ✓ | – | ✓ | – | – |
| C65 | DIAGNOSED | – | – | – | – | – | – | ✓ | ✓ | ✓ | – | – |
| C66 | DIAGNOSED | – | – | – | – | – | – | ✓ | – | ✓ | – | – |
| C67 | DIAGNOSED | – | – | – | – | – | – | ✓ | – | ✓ | – | – |
| C68 | LIMITATION | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | ✓ |
| C69 | DIAGNOSED | – | – | – | – | – | – | ? | ? | ? | – | – |
| C70 | DIAGNOSED | – | – | – | – | – | – | ? | ? | ? | – | – |

| Kontrol | PASS | FAIL | WARN | UNSUPPORTED | INCONCLUSIVE | NOT_RUN |
|---|---:|---:|---:|---:|---:|---:|
| tarayıcı ↔ sunucu | 64 | 0 | 0 | 0 | 0 | 6 |
| FEM Lab ↔ ifem | 63 | 0 | 0 | 0 | 0 | 7 |
| scikit-fem eleman | 64 | 0 | 0 | 0 | 0 | 6 |
| scikit-fem global K | 61 | 0 | 0 | 0 | 0 | 9 |
| scikit-fem düğüm | 64 | 0 | 0 | 0 | 0 | 6 |
| analitik çözüm | 57 | 0 | 0 | 0 | 7 | 6 |
| sürekli varlık | 61 | 0 | 0 | 0 | 9 | 0 |
| sürekli teklik | 57 | 0 | 0 | 0 | 10 | 3 |
| uyumluluk | 4 | 0 | 0 | 0 | 10 | 56 |
| L²/H¹ yeniden hesap | 57 | 0 | 0 | 0 | 0 | 13 |
| sınır koşulu | 63 | 0 | 0 | 0 | 0 | 7 |

## Başarısızlıklar (kanıt ve kaydedilen katman)

Yok.

## Uyarılar

- C07 — uygulama doğrulaması: yakınsama çalışması (ifem, analitik çözüme karşı): L²: ardışık hızlar 1.46, 1.24: asimptotik bölgeye girilmedi; hız iddiası yok; H¹: ardışık hızlar 1.40, 1.04: asimptotik bölgeye girilmedi; hız iddiası yok
- C08 — uygulama doğrulaması: yakınsama çalışması (ifem, analitik çözüme karşı): L²: ardışık hızlar 0.83, 1.15: asimptotik bölgeye girilmedi; hız iddiası yok; H¹: ardışık hızlar 0.15, 0.37: asimptotik bölgeye girilmedi; hız iddiası yok
- C08 — uygulama doğrulaması: 1. mertebe adım göstergesi: λ = max|γ/β| = 50, λh/p = 12.5 > 1: büyük adım; standart Galerkin bu sınıfta kararlılık garantisi vermez
- C10 — uygulama doğrulaması: yakınsama çalışması (ifem, analitik çözüme karşı): L²: ardışık hızlar 4.37, 4.00: asimptotik bölgeye girilmedi; hız iddiası yok; H¹: gözlenen hız ≈ 3.00 (bu sınıf için kuramsal hız iddia edilmiyor)
- C67 — tarayıcı çekirdeği aynı girdiyi: çözdü (tarayıcı sürekli düzeyi bilmez; ayrık sistem tekil değil)
- C68 — uygulama doğrulaması: yakınsama çalışması (ifem, analitik çözüme karşı): L²: gözlenen ≈ 3.95, beklenen 3: beklenen hız gözlenmedi (ağlar ön-asimptotik olabilir, ör. kaba kh/Pe, ya da düzgünlük varsayımı sağlanmıyor); H¹: ardışık hızlar 2.33, 2.03: asimptotik bölgeye girilmedi; hız iddiası yok
- C70 — tarayıcı çekirdeği aynı girdiyi: çözdü (tarayıcı sürekli düzeyi bilmez; ayrık sistem tekil değil)

## scikit-fem bağımsız doğrulaması

- Çalıştırılan: 70 · çözüldü: 68 · tekil: 2 · desteklenmeyen: 0
- Eleman düzeyi karşılaştırılan vaka: 64; en büyük göreli fark 9.23e-15 (mutlak 1.14e-12, C27 K^(6) girdi (3, 3))
- Global K karşılaştırılan: 61; en büyük göreli fark 9.23e-15
- Düğüm katsayısı karşılaştırılan: 64; en büyük göreli fark 1.02e-12
- En büyük kuadratür tahmini (intorder 2p+40 ↔ 2p+8): 4.10e-14
- Taban dönüşümü: P1/P2 düğümsel (yalnız serbestlik derecesi sırası farklı); ElementLinePp(p ≥ 3) iç serbestlik derecelerinde hiyerarşik. Φ[J, m] = ψ_m(x_J), K_Lagrange = Φ⁻ᵀ K Φ⁻¹, F_Lagrange = Φ⁻ᵀ F (eleman ve global düzeyde).

## Sürekli problem doğrulaması

| ID | Sınıf | Düzenlilik (R) | Varlık | Teklik | Uyumluluk |
|---|---|---|---|---|---|
| C01 | IVP, mertebe 1 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C02 | IVP, mertebe 1 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C03 | IVP, mertebe 1 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C04 | IVP, mertebe 1 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C05 | IVP, mertebe 1 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C06 | IVP, mertebe 1 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C07 | IVP, mertebe 1 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C08 | IVP, mertebe 1 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C09 | IVP, mertebe 1 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C10 | IVP, mertebe 1 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C11 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C12 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C13 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C14 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C15 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C16 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C17 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C18 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C19 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C20 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C21 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C22 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C23 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C24 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C25 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C26 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C27 | BVP, mertebe 2 | evet | INCONCLUSIVE | INCONCLUSIVE | INCONCLUSIVE |
| C28 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C29 | BVP, mertebe 2 | evet | INCONCLUSIVE | INCONCLUSIVE | INCONCLUSIVE |
| C30 | BVP, mertebe 2 | evet | INCONCLUSIVE | INCONCLUSIVE | INCONCLUSIVE |
| C31 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C32 | BVP, mertebe 2 | evet | INCONCLUSIVE | INCONCLUSIVE | INCONCLUSIVE |
| C33 | BVP, mertebe 2 | hayır | EXISTENCE_CONFIRMED | INCONCLUSIVE | INCONCLUSIVE |
| C34 | BVP, mertebe 2 | evet | INCONCLUSIVE | INCONCLUSIVE | INCONCLUSIVE |
| C35 | BVP, mertebe 2 | evet | INCONCLUSIVE | INCONCLUSIVE | INCONCLUSIVE |
| C36 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C37 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C38 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C39 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C40 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C41 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C42 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C43 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C44 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C45 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C46 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C47 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C48 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C49 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C50 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C51 | IVP, mertebe 1 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C52 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C53 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C54 | IVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C55 | IVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C56 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C57 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C58 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C59 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C60 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C61 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C62 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C63 | BVP, mertebe 2 | hayır | INCONCLUSIVE | INCONCLUSIVE | INCONCLUSIVE |
| C64 | BVP, mertebe 2 | evet | NONEXISTENCE_CONFIRMED | NOT_APPLICABLE | COMPATIBILITY_VIOLATED |
| C65 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | NONUNIQUENESS_CONFIRMED | COMPATIBILITY_SATISFIED |
| C66 | BVP, mertebe 2 | evet | NONEXISTENCE_CONFIRMED | NOT_APPLICABLE | COMPATIBILITY_VIOLATED |
| C67 | BVP, mertebe 2 | evet | NONEXISTENCE_CONFIRMED | NOT_APPLICABLE | COMPATIBILITY_VIOLATED |
| C68 | BVP, mertebe 2 | evet | EXISTENCE_CONFIRMED | UNIQUENESS_CONFIRMED | NOT_APPLICABLE |
| C69 | BVP, mertebe 2 | hayır | INCONCLUSIVE | INCONCLUSIVE | INCONCLUSIVE |
| C70 | BVP, mertebe 2 | hayır | INCONCLUSIVE | INCONCLUSIVE | INCONCLUSIVE |

Varlık durumları: EXISTENCE_CONFIRMED: 58, INCONCLUSIVE: 9, NONEXISTENCE_CONFIRMED: 3

### C64–C70 ayrıntıları

**C64** (BVP, mertebe 2; düzenlilik: önde gelen katsayı sabit ve sıfırdan farklı)
- existence: **NONEXISTENCE_CONFIRMED** — Yöntem 2: doğrulanmış temel çözüm sistemi + sınır matrisi (düzenli SDP) L[y1] = L[y2] = 0 ve L[u_p] = f sembolik olarak doğrulandı; W(1/2) ≠ 0 (50 basamak: |değer| = 1.000e+00); B = [[0, 1], [0, 1]], r = [0, -1], det B = 0 (sembolik); rank B = 1, rank [B|r] = 2 (sembolik; 60 basamaklı karmaşık SVD: 1, 2): r, B'nin değer uzayında değil → hiçbir u_p + c·y koşulları sağlamaz
- uniqueness: **NOT_APPLICABLE** —  çözüm yok
- compatibility: **COMPATIBILITY_VIOLATED** — Yöntem 3: Fredholm alternatifi, Green özdeşliği ∫ v L u − ∫ u L v = [α(v u′ − u v′)]ₐᵇ v = 1: ∫ f v dx = 1 (sembolik integral), [α v u′ − α u v′]ₐᵇ = 0, fark C = 1 (50 basamak: |değer| = 1.000e+00) (Yöntem 2 ile tutarlı)

**C65** (BVP, mertebe 2; düzenlilik: önde gelen katsayı sabit ve sıfırdan farklı)
- existence: **EXISTENCE_CONFIRMED** — Yöntem 2: doğrulanmış temel çözüm sistemi + sınır matrisi (düzenli SDP) L[y1] = L[y2] = 0 ve L[u_p] = f sembolik olarak doğrulandı; W(1/2) ≠ 0 (50 basamak: |değer| = 1.000e+00); B = [[0, 1], [0, 1]], r = [0, 0], det B = 0 (sembolik); rank B = 1, rank [B|r] = 1 (sembolik; 60 basamaklı karmaşık SVD: 1, 1)
- uniqueness: **NONUNIQUENESS_CONFIRMED** — Yöntem 2: doğrulanmış temel çözüm sistemi + sınır matrisi (düzenli SDP) homojen problemin sıfırdan farklı çözümü: v = 1
- compatibility: **COMPATIBILITY_SATISFIED** — Yöntem 3: Fredholm alternatifi, Green özdeşliği ∫ v L u − ∫ u L v = [α(v u′ − u v′)]ₐᵇ v = 1: ∫ f v dx = 0 (sembolik integral), [α v u′ − α u v′]ₐᵇ = 0, fark C = 0 (sembolik) (Yöntem 2 ile tutarlı)

**C66** (BVP, mertebe 2; düzenlilik: önde gelen katsayı sabit ve sıfırdan farklı)
- existence: **NONEXISTENCE_CONFIRMED** — Yöntem 2: doğrulanmış temel çözüm sistemi + sınır matrisi (düzenli SDP) L[y1] = L[y2] = 0 ve L[u_p] = f sembolik olarak doğrulandı; W(1/2) ≠ 0 (50 basamak: |değer| = 1.000e+00); B = [[0, 1], [0, 1]], r = [0, -1], det B = 0 (sembolik); rank B = 1, rank [B|r] = 2 (sembolik; 60 basamaklı karmaşık SVD: 1, 2): r, B'nin değer uzayında değil → hiçbir u_p + c·y koşulları sağlamaz
- uniqueness: **NOT_APPLICABLE** —  çözüm yok
- compatibility: **COMPATIBILITY_VIOLATED** — Yöntem 3: Fredholm alternatifi, Green özdeşliği ∫ v L u − ∫ u L v = [α(v u′ − u v′)]ₐᵇ v = 1: ∫ f v dx = 1 (sembolik integral), [α v u′ − α u v′]ₐᵇ = 0, fark C = 1 (50 basamak: |değer| = 1.000e+00) (Yöntem 2 ile tutarlı)

**C67** (BVP, mertebe 2; düzenlilik: önde gelen katsayı sabit ve sıfırdan farklı)
- existence: **NONEXISTENCE_CONFIRMED** — Yöntem 2: doğrulanmış temel çözüm sistemi + sınır matrisi (düzenli SDP) L[y1] = L[y2] = 0 ve L[u_p] = f sembolik olarak doğrulandı; W(1/2) ≠ 0 (50 basamak: |değer| = 3.142e+00); B = [[1, 0], [-1, 0]], r = [0, -1/(2*pi)], det B = 0 (sembolik); rank B = 1, rank [B|r] = 2 (sembolik; 60 basamaklı karmaşık SVD: 1, 2): r, B'nin değer uzayında değil → hiçbir u_p + c·y koşulları sağlamaz
- uniqueness: **NOT_APPLICABLE** —  çözüm yok
- compatibility: **COMPATIBILITY_VIOLATED** — Yöntem 3: Fredholm alternatifi, Green özdeşliği ∫ v L u − ∫ u L v = [α(v u′ − u v′)]ₐᵇ v = sin(pi*x): ∫ f v dx = 1/2 (sembolik integral), [α v u′ − α u v′]ₐᵇ = 0, fark C = 1/2 (50 basamak: |değer| = 5.000e-01) (Yöntem 2 ile tutarlı)

**C68** (BVP, mertebe 2; düzenlilik: önde gelen katsayı sabit ve sıfırdan farklı)
- existence: **EXISTENCE_CONFIRMED** — Yöntem 2: doğrulanmış temel çözüm sistemi + sınır matrisi (düzenli SDP) L[y1] = L[y2] = 0 ve L[u_p] = f sembolik olarak doğrulandı; W(1/2) ≠ 0 (50 basamak: |değer| = 3.130e+00); B = [[0, 1], [sin(7*sqrt(5)/5), cos(7*sqrt(5)/5)]], r = [-5/49, -5/49], det B = -sin(7*sqrt(5)/5) (50 basamak: |değer| = 1.110e-02); aday denkleme ve her koşula sembolik olarak kondu
- uniqueness: **UNIQUENESS_CONFIRMED** — Yöntem 2: doğrulanmış temel çözüm sistemi + sınır matrisi (düzenli SDP) det B ≠ 0: homojen problemin yalnız sıfır çözümü var
- compatibility: **NOT_APPLICABLE** —  det B ≠ 0: çekirdek yalnız sıfır, uyumluluk koşulu yok

**C69** (BVP, mertebe 2; düzenlilik: önde gelen katsayı x - 1/2 sıfırlanıyor: x ∈ [1/2])
- existence: **INCONCLUSIVE** —  (R) yok ve aday u = x - I*log(2*x - 1)/pi - 1 gerçek değerli gösterilemedi (Im u = -log(Abs(2*x - 1))/pi)
- uniqueness: **INCONCLUSIVE** —  (R) kurulamadı (önde gelen katsayı x - 1/2 sıfırlanıyor: x ∈ [1/2]): tekil problemde teklik kuramı uygulanmaz
- compatibility: **INCONCLUSIVE** —  (R) yok: temel sistem aralığın tamamında geçerli olmayabilir; det B ≠ 0 çekirdeğin sıfır olduğunu kanıtlamaz

**C70** (BVP, mertebe 2; düzenlilik: önde gelen katsayı x sıfırlanıyor: x ∈ [0])
- existence: **INCONCLUSIVE** —  u(0): temel çözümlerden biri bu uçta sonlu değil (tekil uç); sınır matrisi kurulamadı
- uniqueness: **INCONCLUSIVE** —  u(0): temel çözümlerden biri bu uçta sonlu değil (tekil uç); sınır matrisi kurulamadı
- compatibility: **INCONCLUSIVE** —  u(0): temel çözümlerden biri bu uçta sonlu değil (tekil uç); sınır matrisi kurulamadı

## Taban (değişiklik öncesi) ile karşılaştırma

- Her iki çalıştırmada kâğıt üretilen vaka: 63; düğüm değerlerinde en büyük göreli değişim: 6.66e-16
- HTTP durumu ya da düğüm değeri değişen vakalar:
  - C46: 500 → 200 · önce: Çözüm kâğıdı üretilemedi: IndexError: list index out of range
