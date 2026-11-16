# FEM Lab — 80 vakalık girdi korpusu: yürütme raporu

scikit-fem 11.0.0 · korpus: `verification/corpus/cases.json` (yalnız girdi) · sonuçlar: `results.json`, `skfem_results.json`

## Kapsam

- **Kategori:** derivative_flux_bc: 10, element_order_mesh: 10, first_order_ivp: 10, invalid_adversarial: 10, second_order_bvp: 15, singular_degenerate_compat: 7, trig_exp_source: 8, variable_coefficient: 10
- **Geçerlilik:** challenging: 12, invalid: 10, valid: 58
- **Mertebe:** 1: 13, 2: 67
- **Derece p:** 1: 32, 2: 28, 3: 11, 4: 8, 5: 1
- **n:** 0: 1, 1: 2, 10: 1, 16: 1, 2: 3, 3: 6, 32: 1, 4: 37, 5: 8, 6: 8, 64: 3, 7: 2, 8: 6, 9: 1
- **Koşul türleri:** both: 3, dirichlet: 121, flux: 4, neumann: 11, none: 21

## Yürütme özeti

| Durum | Vaka |
|---|---:|
| PASS | 59 |
| WARN | 1 |
| LIMITATION | 3 |
| DIAGNOSED | 6 |
| REJECTED_OK | 8 |
| FAIL | 3 |
| UNAVAILABLE | 0 |
| RUNNER_ERROR | 0 |

Toplam 80 · uyarı içeren kontrol: 7

## Vaka tablosu

| ID | Kategori | Geçerlilik | Durum | Sunucu | Süre (s) | Ad |
|---|---|---|---|---|---:|---|
| C01 | first_order_ivp | valid | PASS | 200 | 0.3 | u' = 2x, u(0)=0 on [0,2] |
| C02 | first_order_ivp | valid | PASS | 200 | 0.5 | u' + u = 1, u(0)=0 |
| C03 | first_order_ivp | valid | PASS | 200 | 1.0 | u' + 2u = 0, u(0)=1 |
| C04 | first_order_ivp | valid | PASS | 200 | 3.3 | u' = cos x, u(0)=0 on [0,3] |
| C05 | first_order_ivp | valid | PASS | 200 | 0.7 | u' = x^2, u(1)=1 (sağ uç) |
| C06 | first_order_ivp | valid | PASS | 200 | 0.9 | (1+x) u' = 1, u(0)=0 |
| C07 | first_order_ivp | challenging | LIMITATION | 200 | 0.4 | u' - u = 0, u(0)=1 on [0,3] |
| C08 | first_order_ivp | challenging | LIMITATION | 200 | 0.5 | u' + 50u = 0, u(0)=1 |
| C09 | first_order_ivp | valid | PASS | 200 | 5.8 | u' + x u = x, u(0)=0 on [0,2] |
| C10 | first_order_ivp | valid | WARN | 200 | 3.1 | u' = exp(-x), u(-1)=0 on [-1,1] |
| C11 | second_order_bvp | valid | PASS | 200 | 0.4 | -u'' = 1, D/D |
| C12 | second_order_bvp | valid | PASS | 200 | 0.5 | u'' + u = x, D/D, n=2 p=1 |
| C13 | second_order_bvp | valid | PASS | 200 | 2.0 | -u'' = sin(pi x), P1 |
| C14 | second_order_bvp | valid | PASS | 200 | 0.5 | u'' = -2, D/D, P2 |
| C15 | second_order_bvp | valid | PASS | 200 | 1.5 | u'' - 4u = 0, u(0)=0, u(1)=1 |
| C16 | second_order_bvp | valid | PASS | 200 | 0.7 | -u'' + u' = 1, D/D |
| C17 | second_order_bvp | valid | PASS | 200 | 0.7 | u'' = 6x on [-1,2] |
| C18 | second_order_bvp | valid | PASS | 200 | 0.4 | u'' = 0, u(0)=2, u(3)=5 |
| C19 | second_order_bvp | valid | PASS | 200 | 1.1 | -u'' + 9u = 0 on [0,2] |
| C20 | second_order_bvp | valid | PASS | 200 | 1.0 | u'' = 12x^2, n=2 p=4 |
| C21 | second_order_bvp | valid | PASS | 200 | 2.4 | -u'' = exp(x), P3 |
| C22 | second_order_bvp | valid | PASS | 200 | 0.7 | u'' + 4u = 0, u(0)=0, u(1)=1 |
| C23 | second_order_bvp | valid | PASS | 200 | 1.4 | -u'' = 1 on [0,10] |
| C24 | second_order_bvp | valid | PASS | 200 | 0.5 | u'' = 1 on [5,6] |
| C25 | second_order_bvp | valid | PASS | 200 | 0.5 | -2u'' = 4, D/D |
| C26 | variable_coefficient | valid | PASS | 200 | 0.9 | ((1+x) u')' = 1 |
| C27 | variable_coefficient | valid | PASS | 200 | 1.7 | ((1+x^2) u')' + u = x |
| C28 | variable_coefficient | valid | PASS | 200 | 3.1 | (e^x u')' = e^x |
| C29 | variable_coefficient | valid | PASS | 200 | 7.3 | α=2+sin x, β=x, γ=-1, f=1 |
| C30 | variable_coefficient | valid | PASS | 200 | 2.1 | ((1+x) u')' + x u = x^2, D/doğal |
| C31 | variable_coefficient | valid | PASS | 200 | 60.9 | (u'/(1+x))' = 1 |
| C32 | variable_coefficient | valid | PASS | 200 | 67.3 | (sqrt(1+x) u')' = 1 |
| C33 | variable_coefficient | challenging | PASS | 200 | 0.7 | (x u')' = 1, sol uç dejenere, u(1)=0 |
| C34 | variable_coefficient | valid | PASS | 200 | 1.2 | ((x+1) u')' - x u' = 0 |
| C35 | variable_coefficient | valid | PASS | 200 | 38.7 | ((1+0.5cos(pi x)) u')' - u = sin(pi x), D/akı |
| C36 | derivative_flux_bc | valid | PASS | 200 | 0.4 | u'' = -1, u(0)=0, u'(1)=0 |
| C37 | derivative_flux_bc | valid | PASS | 200 | 0.6 | -u'' = 1, u'(0)=1, u(1)=0 |
| C38 | derivative_flux_bc | valid | PASS | 200 | 0.7 | ((1+x)u')' = 4x+2, (1+x)u'(1) = 4 |
| C39 | derivative_flux_bc | valid | PASS | 200 | 0.7 | ((1+x)u')' = 4x+2, u'(1) = 2 |
| C40 | derivative_flux_bc | valid | PASS | 200 | 1.0 | ((1+x^2)u')' = 2, akı sol uçta |
| C41 | derivative_flux_bc | valid | PASS | 200 | 1.5 | u'' - u = 1, u(0)=0, sağ verilmedi |
| C42 | derivative_flux_bc | valid | PASS | 200 | 0.4 | u'' - u = 1, iki uç verilmedi |
| C43 | derivative_flux_bc | valid | PASS | 200 | 0.7 | ((1+2x)u')' = 0, u(0)=0, u'(1)=1 |
| C44 | derivative_flux_bc | valid | PASS | 200 | 1.7 | 2u'' = x, u(0)=1, 2u'(1) = -2 |
| C45 | derivative_flux_bc | valid | PASS | 200 | 0.5 | -u'' = x on [2,4], u(2)=0, u'(4)=1 |
| C46 | element_order_mesh | valid | FAIL | 500 | 0.3 | u'' = 1, n=1 p=1 |
| C47 | element_order_mesh | valid | PASS | 200 | 1.7 | u'' = sin x, n=1 p=4 |
| C48 | element_order_mesh | challenging | PASS | 200 | 20.7 | -u'' = sin(pi x), n=64 p=2 |
| C49 | element_order_mesh | valid | PASS | 200 | 9.0 | u'' + u = x, n=32 p=3 |
| C50 | element_order_mesh | challenging | PASS | 200 | 49.5 | -u'' = 1, n=64 p=4 |
| C51 | element_order_mesh | valid | PASS | 200 | 5.2 | u' = cos x, n=64 p=1 on [0,6] |
| C52 | element_order_mesh | valid | PASS | 200 | 2.9 | ((1+x)u')' = x, n=7 p=3 |
| C53 | element_order_mesh | valid | PASS | 200 | 1.4 | u'' = x^3, n=9 p=2 |
| C54 | element_order_mesh | valid | PASS | 200 | 1.4 | u'' + u = 0, u(0)=0, u'(0)=1, n=3 p=4 |
| C55 | element_order_mesh | valid | PASS | 200 | 1.1 | u'' + u = 0, IVP on [0,6], n=16 p=1 |
| C56 | trig_exp_source | valid | PASS | 200 | 7.0 | -u'' = pi^2 sin(pi x), P2 |
| C57 | trig_exp_source | valid | PASS | 200 | 2.6 | u'' = e^x, u(0)=1, u(1)=e |
| C58 | trig_exp_source | valid | PASS | 200 | 17.0 | -u'' = cos(2 pi x), u'(0)=0, u(1)=0 |
| C59 | trig_exp_source | valid | PASS | 200 | 2.6 | u'' + u = sin(3x) on [0,2] |
| C60 | trig_exp_source | valid | PASS | 200 | 4.0 | -u'' = x e^-x, P4 |
| C61 | trig_exp_source | valid | PASS | 200 | 2.5 | u'' = sinh x |
| C62 | trig_exp_source | valid | PASS | 200 | 20.9 | -u'' = e^x sin x, P3 |
| C63 | trig_exp_source | valid | PASS | 200 | 48.2 | -u'' = tan x |
| C64 | singular_degenerate_compat | challenging | DIAGNOSED | 400 | 0.2 | u'' = 1, u'(0)=0, u'(1)=0 |
| C65 | singular_degenerate_compat | challenging | DIAGNOSED | 400 | 0.2 | u'' = 2x-1, u'(0)=0, u'(1)=0 |
| C66 | singular_degenerate_compat | challenging | DIAGNOSED | 400 | 0.2 | u'' = 1, koşulsuz |
| C67 | singular_degenerate_compat | challenging | DIAGNOSED | 400 | 385.7 | u'' + pi^2 u = sin(pi x), D/D |
| C68 | singular_degenerate_compat | challenging | LIMITATION | 200 | 2.7 | u'' + 9.8u = 1, D/D |
| C69 | singular_degenerate_compat | challenging | DIAGNOSED | 400 | 0.3 | ((x-0.5)u')' = 1, D/D |
| C70 | singular_degenerate_compat | challenging | DIAGNOSED | 400 | 0.5 | (x u')' = 1, u(0)=u(1)=0 |
| C71 | invalid_adversarial | invalid | REJECTED_OK | 400 | 0.2 | bozuk ifade '2x+' |
| C72 | invalid_adversarial | invalid | REJECTED_OK | 400 | 0.2 | bilinmeyen fonksiyon foo(x) |
| C73 | invalid_adversarial | invalid | REJECTED_OK | 400 | 0.2 | kod enjeksiyonu denemesi |
| C74 | invalid_adversarial | invalid | REJECTED_OK | 400 | 0.2 | bilimsel gösterim 1e-8 |
| C75 | invalid_adversarial | invalid | REJECTED_OK | 400 | 0.2 | x0 = x1 |
| C76 | invalid_adversarial | invalid | FAIL | 400 | 0.2 | p = 5 |
| C77 | invalid_adversarial | invalid | FAIL | 400 | 0.3 | n = 0 |
| C78 | invalid_adversarial | invalid | REJECTED_OK | 400 | 0.3 | 1. mertebe + türev koşulu |
| C79 | invalid_adversarial | invalid | REJECTED_OK | 400 | 0.3 | aşırı belirli: sol uçta u ve u', sağda u |
| C80 | invalid_adversarial | invalid | REJECTED_OK | 400 | 0.3 | α = β = 0 |

## Başarısızlıklar (kanıt ve kaydedilen katman)

- **C46** (u'' = 1, n=1 p=1) — sunucu hattı: HTTP 500: Çözüm kâğıdı üretilemedi: IndexError: list index out of range · katman: unresolved (sunucu çöküşü/zaman aşımı)
- **C76** (p = 5) — tarayıcı çekirdeği temiz ret: crash: TypeError: Cannot read properties of undefined (reading 'N') · katman: input validation (tarayıcı çekirdeği: buildProblem/makeMesh)
- **C77** (n = 0) — tarayıcı çekirdeği temiz ret: ok: girdi kabul edildi, c = [0] · katman: input validation (tarayıcı çekirdeği: buildProblem/makeMesh)

## Uyarılar

- C07 — uygulama doğrulaması: yakınsama çalışması (ifem, analitik çözüme karşı): L²: ardışık hızlar 1.46, 1.24: asimptotik bölgeye girilmedi; hız iddiası yok; H¹: ardışık hızlar 1.40, 1.04: asimptotik bölgeye girilmedi; hız iddiası yok
- C08 — uygulama doğrulaması: yakınsama çalışması (ifem, analitik çözüme karşı): L²: ardışık hızlar 0.83, 1.15: asimptotik bölgeye girilmedi; hız iddiası yok; H¹: ardışık hızlar 0.15, 0.37: asimptotik bölgeye girilmedi; hız iddiası yok
- C08 — uygulama doğrulaması: 1. mertebe adım göstergesi: λ = max|γ/β| = 50, λh/p = 12.5 > 1: büyük adım; standart Galerkin bu sınıfta kararlılık garantisi vermez
- C10 — uygulama doğrulaması: yakınsama çalışması (ifem, analitik çözüme karşı): L²: ardışık hızlar 4.37, 4.00: asimptotik bölgeye girilmedi; hız iddiası yok; H¹: gözlenen hız ≈ 3.00 (bu sınıf için kuramsal hız iddia edilmiyor)
- C67 — tarayıcı çekirdeği aynı girdiyi: çözdü (tarayıcı sürekli düzeyi bilmez; ayrık sistem tekil değil)
- C68 — uygulama doğrulaması: yakınsama çalışması (ifem, analitik çözüme karşı): L²: gözlenen ≈ 3.95, beklenen 3: beklenen hız gözlenmedi (ağlar ön-asimptotik olabilir, ör. kaba kh/Pe, ya da düzgünlük varsayımı sağlanmıyor); H¹: ardışık hızlar 2.33, 2.03: asimptotik bölgeye girilm
- C70 — tarayıcı çekirdeği aynı girdiyi: çözdü (tarayıcı sürekli düzeyi bilmez; ayrık sistem tekil değil)

## Bağımsız doğrulamanın yapılamadığı vakalar

- C27 — analitik çözüm: kapalı biçimde genel çözüm bulunamadı
- C29 — analitik çözüm: kapalı biçimde genel çözüm bulunamadı
- C30 — analitik çözüm: kapalı biçimde genel çözüm bulunamadı
- C32 — analitik çözüm: genel çözüm sabitlerde doğrusal değil
- C34 — analitik çözüm: kapalı biçimde genel çözüm bulunamadı
- C35 — analitik çözüm: kapalı biçimde genel çözüm bulunamadı
- C40 — analitik çözüm: sembolik hesap tamamlanamadı (TypeError)
- C63 — analitik çözüm: kapalı biçimde genel çözüm bulunamadı

## scikit-fem bağımsız doğrulaması

- Çalıştırılan (geçerli + zorlayıcı) vaka: 70 · çözüldü: 68 · tekil: 2 · desteklenmeyen: 0
- FEM Lab ile karşılaştırılan (kâğıt üretilen) vaka: 63; global K karşılaştırılan: 60 (≤ 65 bilinmeyen)
- En büyük göreli fark — düğüm katsayıları: 1.02e-12, global K: 9.23e-15
- En büyük kuadratür tahmini (intorder 2p+40 ↔ 2p+8): 4.10e-14
- Geçersiz girdiler scikit-fem'e verilmez (FEM Lab reddettiği için karşılaştırılacak çözüm yoktur); UNSUPPORTED hiçbir zaman PASS sayılmaz.
