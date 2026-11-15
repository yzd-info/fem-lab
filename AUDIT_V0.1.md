# FEM Lab V0.1 — Bilimsel Denetim Raporu

Tarih: 2026-11-15 · Kapsam: mevcut 1B uygulama (`index.html`, `fem_sheet.py`, `server.py`) · **Kod değiştirilmedi.**

> **Bu raporun kanıt düzeyi.** Testler, uygulamadan bağımsız yazılmış ikinci bir FEM gerçekleştirmesiyle
> (Python/mpmath, 30 basamak) ve analitik çözümlerle yapıldı. Bağımsız gerçekleştirme aynı yapay zekâ tarafından
> yazıldı; *kod bağımsızlığı* sağlar, *matematiksel kanıt* değildir. Birkaç bulgu ayrıca kâğıt üstü türetmeyle
> (P1 kararlılık sınırı, satır sayımı, rank) desteklendi. Bir yöntemin sayısal olarak daha iyi sonuç vermesi doğruluk
> kanıtı sayılmadı.

---

## 1. Yönetici özeti

- **Hesap çekirdeği, iyi tanımlı problemlerde doğru çalışıyor.** Doğrulananlar:
  - Şekil fonksiyonları P1–P4: δᵢⱼ, birimin ayrışımı, tamlık, türevler, C⁰ süreklilik.
  - Eleman matrisleri ve montaj: girdi girdi 10⁻¹³ içinde.
  - Hata ölçüleri: L², H¹, E_düğüm, E_max.
  - Yakınsama hızları: analitik çözümü bilinen 16 problemde beklenen değerler.
  - Tarayıcı ↔ bağımsız gerçekleştirme: her testte 10⁻¹² içinde aynı katsayılar.
- **Ama uygulama, iyi tanımlı olmayan problemleri tanımıyor:**
  - Çözümü olmayan, tek olmayan ve dejenere problemleri uyarısız "çözüyor" (E1, E2, E3).
  - Bir durumda (dejenere α) **sınır koşulunu sağlamayan bir "analitik çözümü" doğrulanmış gibi gösteriyor**, ve tarayıcı hataları bu yanlış referansa göre ölçüyor (KR-1).
- **Birinci mertebe:** satır (test fonksiyonu) seçimi matematiksel olarak savunulabilir (ders kitabı kuralı). Ama atılan denklemin kalıntısı ve uyumluluk koşulu hiçbir yerde gösterilmiyor; seçimin sonucu değiştirebildiği kullanıcıya söylenmiyor.
- **Çözüm kâğıdı** gerçek hesapla tutarlı, ama 12 adımlık şeffaflık standardının yaklaşık yarısını karşılamıyor. Şekil fonksiyonları türetilmeden veriliyor (V0.2 kapsamı).
- **V0.1 tamamlanmadı.** Tanımdaki "sınır koşulları matematiksel olarak tutarlı" ve "çözücü davranışı uygun teşhis ediliyor" maddeleri karşılanmıyor (bölüm 15).

---

## 2. Mevcut yetenekler: işlem → dosya/fonksiyon haritası

| İşlem | Tarayıcı (`index.html`, CORE) | Kâğıt (`fem_sheet.py`) |
|---|---|---|
| Girdi ayrıştırma | `parseExpression` | `parse` (beyaz liste) |
| Problem doğrulama, mertebe, koşul türleri | `buildProblem`, `readCondition` | `solve_problem` (baş kısmı) |
| Şekil fonksiyonları | `shapeFunctions`, `shapeDerivatives` (ξ∈[−1,1], ürün formülü) | `lagrange_basis` (s∈[0,h], ürün formülü) |
| Kuadratür | `gaussLegendre`, `compositeRule` (12 nokta × 8 parça) | `integral` (SymPy tam; olmazsa `mpmath.quad`, "≈") |
| Eleman matrisi / yük | `elementSystem` | `solve_problem` eleman döngüsü + `integral` |
| Montaj | `METHODS.galerkin.solve` | `solve_problem` |
| Sınır koşulu, satır seçimi | `METHODS.galerkin.solve` (`fixed`, `rows`) | `solve_problem` (`fixed`, `rows`) |
| Çözücü | `solveLinear` (kısmi pivotlu Gauss, göreli pivot eşiği 1e-11) | `LUsolve` (≤8 serbest, tam) / `mpmath.lu_solve` + `cond > 1e13` |
| Uç akısı | — | `solve_problem` (`reactions`) |
| Analitik çözüm | `requestExact` → `/api/exact` → `parseExact` | `exact_solution` (dsolve, 12 s), `exact_for_browser` |
| Hata ölçüleri | `errorNorms`, `maxErrorRefined` | `error_measures` |
| Yakınsama | `convergenceStudy`, `observedRate`, `measurementFloors` | — |
| Kâğıt | — | `build_sheet`, `build_tex` |
| Çapraz kontrol | `crossCheck` | `make_sheet` (K, nodal, fixed) |

---

## 3–5. Problem sınıflandırması (A) ve formülasyon denetimi (B)

Denklem: (α u′)′ + β u′ + γ u = f, (x₀, x₁).

| # | Sınıf | Karar | Gerekçe (matematik / gerçekleştirme) | Kanıt |
|---|---|---|---|---|
| 1 | 1. mertebe, tek uçta u verili | **Destekleniyor, yeterince doğrulanmamış** | Doğrudan Galerkin, test uzayı V₀ = {v ∈ V_h : v(u verilen uç) = 0}. Toplam N denklem, N − 1 bilinmeyen; bir denklem atılmalı. Atılan denklemin kalıntısı ≠ 0 ve gösterilmiyor. Optimal hız kuramı yok. | T1–T4, B1 |
| 2 | 1. mertebe, iki uçta koşul (BVP) | **Desteklenmiyor (doğru biçimde reddediliyor)** | Birinci mertebede iki koşul aşırı belirli | E5 |
| 3 | 2. mertebe BVP, sabit katsayı | **Destekleniyor ve doğrulandı** | Bubnov–Galerkin, kısmi integrasyon, temel koşulun satırı atılır; atılan satır bilinmeyen akıyı (Q) içerdiğinden satır sayımı tam (belirsizlik yok); akı geri hesabı doğru | T5, T6, T10, B2 |
| 4 | 2. mertebe başlangıç değer problemi (aynı uçta u ve u′) | **Destekleniyor, yeterince doğrulanmamış** | Petrov–Galerkin (test uzayı: akısı bilinmeyen uçta sıfır). Tutarlı ve kare sistem (kanıtlı). Kararlılık koşullu (P1: ωh < 2√3, türetildi); genel yakınsama kanıtı yok; u′(x₀) yalnız zayıf biçimde sağlanıyor | T14, T15, önceki ωh taraması |
| 5 | Sabit katsayı | **Destekleniyor ve doğrulandı** | — | T5, T6, T9, T10 |
| 6 | Değişken katsayı (düzgün) | **Destekleniyor ve doğrulandı** | Üretilmiş çözümle (α = 1+x², β = sin x, γ = eˣ) L² hızı p + 1; K girdileri bağımsızla 3·10⁻¹³ | T7, D3 |
| 7 | Difüzyon baskın | **Destekleniyor ve doğrulandı** | Koersif durum | T5, T9 |
| 8 | Taşınım baskın | **Kısmen destekleniyor** | Standart Galerkin; Pe = \|β\|h/(2\|α\|) > 1'de salınım (yöntemin bilinen davranışı). Kararlılaştırma yok, Péclet uyarısı yok | T8 |
| 9 | Reaksiyon–difüzyon, γ·α < 0 (koersif) | **Destekleniyor ve doğrulandı** | −u″ + 100u = 100: hızlar p + 1 | T9 |
| 9b | Helmholtz tipi, γ·α > 0 (belirsiz) | **Kısmen destekleniyor** | Kaba ağda ayrık rezonans ve kirlilik (P1, n = 8'de hata 55), uyarı yok. Sürekli problemin rezonans durumu algılanmıyor | T16, E1 |
| 10a | α uçta sıfır, o uç koşulsuz | **Destekleniyor, yeterince doğrulanmamış** | Doğal koşul kendiliğinden sağlanır; tek sınırlı çözüme yakınsıyor | T12, T13 |
| 10b | α uçta sıfır, o uçta u verili | **Desteklenmiyor, ama kabul ediliyor (KUSUR)** | Sınırlı çözüm genelde yok; uygulama yakınsamayan sayılar ve **yanlış analitik çözüm** veriyor | E2, G1 |
| 10c | α iç noktada sıfır ya da işaret değiştiriyor | **Desteklenmiyor, ama kabul ediliyor (KUSUR)** | Tip değişimi (dejenere); n'ye bağlı olarak ya yanıltıcı "tekil" mesajı ya da sessiz sonuç | E3 |
| 11 | Saf Neumann, γ = 0 | **Desteklenmiyor (doğru biçimde reddediliyor)** | K tekil (sabitler sıfır uzayında); ek koşul (ör. ortalama = 0) desteği yok | E4, D2′ |
| 11b | Saf Neumann, γ ≠ 0 | **Destekleniyor ve doğrulandı** | — | T11 |
| 12 | Karışık Dirichlet–Neumann | **Destekleniyor ve doğrulandı** | — | T6 |
| — | Robin (karışık) koşul, periyodik koşul | **Desteklenmiyor** | Arayüzde yok | — |

### Formülasyon ayrıntısı (destekleniyor kabul edilen sınıflar)

| | 2. mertebe BVP | 2. mertebe IVP | 1. mertebe |
|---|---|---|---|
| Güçlü biçim | (αu′)′ + βu′ + γu = f | aynı | βu′ + γu = f |
| Deneme uzayı | V_h, verilen uçlarda u sabit | V_h, u(x₀) = u₀ | V_h, u(x₀) = u₀ |
| Test uzayı | V_h, Dirichlet uçlarında v = 0 | V_h, v(x₁) = 0 | V_h, v(x₀) = 0 |
| Zayıf form | ∫(v′αu′ − vβu′ − vγu) = [vαu′] − ∫vf | aynı; v(x₁) = 0 → [·] = −v(x₀)Q₀ | kısmi integrasyon yok: ∫v(βu′ + γu) = ∫vf |
| Temel koşul | u (satır atılır, sütun sağa) | u(x₀) | u(x₀) |
| Doğal koşul | u′ = değer → Q = α·değer; verilmedi → Q = 0 | u′(x₀) → Q₀ (zayıf biçimde sağlanır) | yok |
| Bilinmeyen / denklem | (N − #D) a + #D akı = N denklem → **kare, belirsizlik yok** | (N − 1) a + Q₁ = N → kare | (N − 1) a, **N denklem → biri atılmalı** |
| Seçilen satırlar | serbest düğümler | 1…N − 1 (sol IVP) | 2…N (u verilen düğüm atılır) |

**Satır seçimi sorununun genelleştirilmesi (istenen araştırma):** Belirsizlik, yalnızca **atılan denklemin içinde bilinmeyen bir akı terimi bulunmadığında** ortaya çıkar. Yani denklem sayısı bilinmeyen sayısını aşar:

- **Birinci mertebe (α ≡ 0):** her zaman.
- **Dejenere uç (α(uç) = 0) ve orada u verilmişse:** akı Q = 0 bilindiği için atılan satır serbest bir bilinmeyen taşımaz → aynı belirsizlik. Bu, 10b'deki tutarsızlığın cebirsel nedenidir.
- **Kalan sınıflar** (2. mertebe BVP, IVP, Neumann, doğal koşul): satır sayımı tam; alternatif satır seçimi yoktur.

**Birinci mertebede tekillik ve uyumluluk (ölçüldü):**
- **γ = 0:** tam K tekildir, rank N − 1 (sabitler sıfır uzayında, çünkü Σⱼ Nⱼ′ = 0). Atılan satır geri kalanlarla ancak wᵀF = 0 ise tutarlıdır (w: sol sıfır vektörü).
- **γ ≠ 0:** K tam ranklıdır ve aşırı belirli sistem genelde tutarsızdır.
- **Uygulama hiçbir durumda uyumluluğu denetlemiyor, göstermiyor.** Alternatif seçimlerle karşılaştırma:
  - u′ = 3x² (n = 8): uygulama 0.18, alternatif 0.063.
  - u′ + 2u = 0: uygulama 0.00035, alternatif 0.072.

  **İki seçim eşdeğer değil, hiçbiri her problemde daha iyi değil.** Uygulamanın seçimi, test uzayının Bubnov–Galerkin yorumuna (u verilen yerde v = 0; Pepper & Heinrich Bölüm 2) karşılık geliyor; bu yüzden savunulabilir.

---

## 6. Kritik matematiksel kusurlar

### KR-1 — Sınır koşulunu sağlamayan "analitik çözüm" doğrulanmış gibi gösteriliyor · **Kritik**
- **Dosya/fonksiyon:** `fem_sheet.py` → `exact_solution`; zincir: `exact_for_browser` → `server.py` `/api/exact` → `index.html` `parseExact` / `compute`.
- **Matematik:** `exact_solution`, dsolve sonucunu yalnızca denklem kalıntısıyla, 3 iç noktada sınıyor; **sınır/başlangıç koşullarını sınamıyor**. dsolve, tekil uçta koşulu sessizce düşürebiliyor.
- **Örnek:** α = x, f = 1, u(0) = 0, u(1) = 0, n = 4, p = 1.
- **Beklenen:** "analitik çözüm bulunamadı / problem iyi tanımlı değil".
- **Gerçekleşen:** kâğıt "u(x) = x − 1" yazıyor; ama u(0) = −1 ≠ 0. Tarayıcıya da aynı ifade gidiyor (`/api/exact` → `x - 1`); arayüz "Referans u: analitik" diyerek hataları bu yanlış referansa göre "E (hata)" olarak gösteriyor.

### KR-2 — İyi tanımlı olmayan problemler uyarısız çözülüyor · **Kritik**
- **Dosya/fonksiyon:** `index.html` → `buildProblem`, `solveLinear`; `fem_sheet.py` → `solve_problem`.
- **Matematik:** İyi tanımlılık (varlık/teklik) hiçbir yerde denetlenmiyor. Tekillik denetimi yalnızca *ayrık* matrisin tam tekilliğini (göreli pivot < 1e-11) yakalıyor; sürekli problemin rezonansını ya da dejenereliğini yakalamıyor.
- **Örnekler:**
  - **E1:** u″ + π²u = sin(πx), u(0) = u(1) = 0 → çözüm **yok**. Uygulama "çözüm" veriyor; max|a| P1'de n = 4 → 64 için 2 → 505, P2'de 198 → 1.3·10⁷ (sınırsız büyüyor). Uyarı yok.
  - **E1b:** aynı denklem, f = 0 → çözüm **tek değil**; uygulama sessizce u = 0 döndürüyor.
  - **E2:** α = x, iki uçta u = 0 → sınırlı çözüm **yok**. Ağ inceldikçe min a: −0.35 → −0.46 → −0.54 (yakınsamıyor), uyarı yok.
  - **E3:** α = x − ½ (işaret değiştiriyor) → n çiftse "Sistem tekil … en az bir uçta u değeri ver" (iki uçta zaten verilmiş; yanlış yönlendirme). n tekse ve düzgün çözüm yokken (u(1) = 2) bile bir sonuç veriyor.
- **Beklenen:** reddetme ya da açık bir "iyi tanımlı değil / rezonansa yakın" uyarısı.

---

## 7. Sayısal gerçekleştirme kusurları

| Kod | Önem | Dosya / fonksiyon | Bulgu | Kanıt |
|---|---|---|---|---|
| SA-1 | Yüksek | `index.html` `METHODS.galerkin.solve`; `fem_sheet.py` `solve_problem` | Belirsiz (Helmholtz tipi) problemlerde ön-asimptotik kirlilik/ayrık rezonans için uyarı yok (kh ≈ 1.25, P1: hata 55, çözüm genliği ~2) | T16 |
| SA-2 | Orta | aynı | Taşınım baskın problemde Péclet uyarısı yok (README'de belgelenmiş, arayüzde değil) | T8 |
| SA-3 | Orta | `index.html` `observedRate` | Hız yalnızca son geçerli çiftten hesaplanıyor; asimptotik bölgede olup olmadığı (ardışık hızların tutarlılığı) denetlenmiyor | T8 (1.66 → 1.83 → 1.95) |
| SA-4 | Düşük | `solveLinear` / `mpmath.cond` | Koşul sayısı ölçülmüyor ve gösterilmiyor. 2. mertebe BVP'de ölçülen 2-norm koşul sayısı ~h⁻² ile büyüyor: P1, n = 32'de 414; P4, n = 32'de 2.1·10⁴. Bugünkü sınırlarda zararsız. | D4 |

Kuadratür: düzgün değişken katsayıda tarayıcı ile bağımsız 24 noktalı Gauss arasındaki K farkı 3·10⁻¹³ (kusur yok). Köşeli yükte (`abs`) kuadratür uyarısı çalışıyor.

## 8. Sınır koşulu kusurları

| Kod | Önem | Bulgu |
|---|---|---|
| SK-1 | Kritik | KR-1 (analitik çözümde koşul denetimi yok) |
| SK-2 | Yüksek | Dejenere uçta (α = 0) u verilmesi kabul ediliyor. Satır sayımı aşırı belirli (§3–5), sınırlı çözüm genelde yok (E2). u′ koşulu bu durumda zaten reddediliyor; u koşulu reddedilmiyor. |
| SK-3 | Orta | Birinci mertebede atılan denklemin kalıntısı ve uyumluluk koşulu (wᵀF) gösterilmiyor ve denetlenmiyor (B1) |
| SK-4 | Düşük | IVP'de u′(x₀) yalnız zayıf biçimde sağlanıyor (kâğıt gösteriyor; arayüz göstermiyor) |

## 9. Montaj ve çözücü kusurları

| Kod | Önem | Bulgu |
|---|---|---|
| MÇ-1 | Orta | Tekillik mesajı nedenini ayırt etmiyor: saf Neumann, dejenere α ve rezonans hepsi "en az bir uçta u değeri ver" diyor (E3: iki uçta u verilmişken bile) |
| — | — | Montaj ve eleman matrisleri **kusursuz**: girdi girdi 10⁻¹³ (D3); β = 0'da K simetrik, β ≠ 0'da beklendiği gibi simetrik değil (D1); γ = 0'da sıfır uzayı sabitler (D2′); akı geri hesabı 10⁻¹⁴ (B2) |

## 10. Hata ve yakınsama analizi kusurları

- **Ölçülerin uygulaması doğru.** Uygulamanın hata ölçüleri bağımsız hesapla göreli farkları:
  - L²: ≤ 3·10⁻¹²
  - H¹: ≤ 4·10⁻¹³
  - E_düğüm: ≤ 7·10⁻¹¹
  - E_max (yerel aramalı): ≤ 2·10⁻¹¹; yalnız örneklenmiş değer biraz düşük, beklendiği gibi.
- **HY-1 (Kritik, KR-1'in sonucu):** "analitik" referans yanlışsa tüm hata kutuları yanlış ama "hata" diye etiketli.
- **HY-2 (Orta):** SA-3 (asimptotik bölge denetimi yok).
- **Ölçüm tabanı doğru çalışıyor:** tabana inmiş hızlar "—" gösteriliyor (T1, T2, T12, T13'te p ≥ 2 ya da p ≥ 3).

## 11. Çözüm kâğıdı kusurları

| Kod | Önem | Bulgu |
|---|---|---|
| ÇK-1 | Kritik | KR-1: yanlış "Tam çözümle karşılaştırma" bölümü |
| ÇK-2 | Yüksek (V0.2) | Şekil fonksiyonları türetilmiyor, hazır formül olarak veriliyor; δᵢⱼ gösterilmiyor |
| ÇK-3 | Orta | 2. mertebede ve C⁰ elemanlarda R(x) = [Lu]_{u=u_h} ve ∫Nᵢ R dx = 0 klasik anlamda yazılıyor. P1'de u_h″ = 0 (eleman içinde), düğümlerde dağılım; bu ifade ancak zayıf form/dağılım anlamında doğrudur ve kâğıt bunu belirtmiyor |
| ÇK-4 | Orta | Test fonksiyonu başına Galerkin denklemi, montaj toplamları (K₂₂ = K₂₂⁽¹⁾ + K₁₁⁽²⁾), satır atma gerekçesi, eleme adımları ve parçalı çözümde yerine koyma adımı yok |
| ÇK-5 | Orta | Katsayılar sabitse 2. ve sonraki elemanların integralleri "aynı" denip atlanıyor; n > 8'de yalnızca 3 eleman |
| ÇK-6 | Orta | Birinci mertebede atılan denklem ve uyumluluk gösterilmiyor (SK-3) |
| ÇK-7 | Düşük | Belge kayması: CLAUDE.md "E_max (96 nokta/eleman)" diyor, kod 32 nokta + altın oran araması kullanıyor; CLAUDE.md §7 "hata ince ağ referansına göre ölçülür" diyor, kod analitik çözümü tercih ediyor; IVP için "ölçülen L² hızı p + 1" bir kanıtmış gibi yazılı |

Kâğıt ile hesap tutarlı: gösterilen her sayı hesaplanan değerle aynı (ödev kıyası ve bağımsız türetme: N, K⁽ᵉ⁾, F⁽ᵉ⁾, K, F, indirgenmiş sistem, a, u_h, L² = 2√15/45).

---

## 12. Bağımsız benchmark sonuçları

Ölçüm: n = 4, 8, 16, 32 ve p = 1–4. Hız L² için ardışık iki ağdan hesaplandı; 10⁻¹¹ altı "—". "Fark" sütunu, uygulamanın katsayıları ile bağımsız gerçekleştirmenin katsayıları arasındaki en büyük mutlak farktır.

| Test | Problem → analitik çözüm | Beklenen | Ölçülen hız (L², p = 1/2/3/4, en ince çift) | Uygulama−bağımsız fark (en büyük) | Durum |
|---|---|---|---|---|---|
| T1 | u′ = 2x, u(0) = 0, [0,2] → x² | p ≥ 2 tam | 2.00 / — / — / — (p ≥ 2: L² ~10⁻¹⁵) | 4.5e-14 | ✓ |
| T2 | u′ = 3x², u(0) = 0 → x³ | p ≥ 3 tam | 2.00 / 2.00 / — / — | 8.2e-15 | ✓ (P2'de hız 2: optimal değil, 1. mertebe için kuramsal garanti yok) |
| T3 | u′ + 2u = 0, u(0) = 1 → e⁻²ˣ | yakınsama | 2.00 / 2.00 / 4.00 / 4.00 | 1.2e-14 | ✓ (çift p'de optimal değil) |
| T4 | u′ = 2x, u(2) = 4 → x² | T1 ile simetrik | T1 ile aynı | 4.8e-14 | ✓ |
| T5 | −u″ = π² sin πx, D/D → sin πx | L² p+1, H¹ p, P1 düğümde tam | 2.00 / 3.00 / 4.00 / 5.00; P1 düğüm ~10⁻¹⁵ | 4.1e-13 | ✓ |
| T6 | u″ − 4u = 0, u(0) = 0, u′(1) = 1 | p+1 | 2.00 / 3.00 / 4.00 / 5.00 | 3.3e-13 | ✓ |
| T7 | Değişken α, β, γ (üretilmiş) → sin πx + x² | p+1 | 2.00 / 3.00 / 4.00 / 5.00 | 2.5e-13 | ✓ |
| T8 | −0.05u″ + u′ = 1 | Pe > 1 salınım; asimptotikte p+1 | 1.95 / 2.89 / 3.89 / 4.89 (ön-asimptotik) | 4.2e-13 | ✓ (beklenen davranış) |
| T9 | −u″ + 100u = 100 → 1 − cosh(10(x−½))/cosh 5 | p+1 | 1.98 / 2.97 / 3.97 / 4.97 | 9.5e-14 | ✓ |
| T10 | u″ + u = x → x − sin x/sin 1 | p+1 | 2.00 / 3.00 / 4.00 / 5.00 | 2.8e-14 | ✓ |
| T11 | u″ − u = −(1+π²)cos πx, iki uç doğal → cos πx | p+1 | 2.00 / 3.00 / 4.00 / 5.00 | 1.6e-12 | ✓ |
| T12 | (xu′)′ = 1, u(1) = 0 → x − 1 | P1'de tam | hepsi ~10⁻¹⁵ | 9.1e-13 | ✓ |
| T13 | (xu′)′ = 2x, u(1) = 0 → (x²−1)/2 | p ≥ 2 tam | 1.98 / — / — / — | 6.7e-13 | ✓ |
| T14 | IVP u″ + u = 0, u(0) = 0, u′(0) = 1 → sin x | p+1 (gözlenen) | 2.00 / 3.00 / 4.00 / 5.00 | 9.0e-13 | ✓ (kanıt yok) |
| T15 | IVP ((1+x)u′)′ = 0 → ln(1+x) | p+1 (gözlenen) | 2.00 / 3.00 / 4.00 / 4.97 | 9.2e-14 | ✓ (kanıt yok) |
| T16 | u″ + 100u = −100, D/D (k = 10) → −1 + cos 10x + B sin 10x | kaba ağda ön-asimptotik | P1: n = 8'de hata 55 → n = 32'de 0.16; P4: 5.0 | 6.6e-12 | ⚠ uyarısız ön-asimptotik |
| E1 | u″ + π²u = sin πx, D/D | çözüm yok → ret/uyarı | max\|a\| sınırsız büyüyor | — | ✗ (KR-2) |
| E1b | u″ + π²u = 0, D/D | tek değil → ret/uyarı | sessizce 0 | — | ✗ (KR-2) |
| E2 | (xu′)′ = 1, u(0) = u(1) = 0 | sınırlı çözüm yok → ret | yakınsamayan sonuç + **yanlış analitik** | — | ✗ (KR-1, KR-2) |
| E3 | ((x−½)u′)′ = 1 | dejenere iç nokta → ret | n'ye bağlı: yanıltıcı mesaj ya da sessiz sonuç | — | ✗ (KR-2, MÇ-1) |
| E4 | u″ = 2x − 1, saf Neumann | sabit farkıyla çözüm → ret ya da ek koşul | doğru biçimde reddediliyor | — | ✓ (desteklenmiyor, açık) |
| E5 | 1. mertebede iki koşul / u′ koşulu; α(0) = 0'da u′(0); α = β = 0 | ret | hepsi reddediliyor | — | ✓ |

**Ek doğrulamalar:**
- **Şekil fonksiyonları P1–P4:** tarayıcı ve SymPy, bağımsız Vandermonde tabanıyla ≤ 4·10⁻¹⁶ aynı. δᵢⱼ, ΣN = 1, ΣN′ = 0 ve xᵐ (m ≤ p) tamlığı 10⁻²⁹ içinde; x^{p+1} beklendiği gibi tam temsil edilmiyor (bu da doğru).
- **Bağımsız katsayılar** (örnek, P2, t ∈ [0,1]): N₁ = 1 − 3t + 2t², N₂ = 4t − 4t², N₃ = −t + 2t².
- **Süreklilik:** u_h eleman sınırlarında sürekli (fark ölçüm adımı düzeyinde); u_h′ sıçraması ≠ 0 (C⁰, beklenen).
- **Ödev kâğıdı kıyası** (u′ = 2x, [0,2], n = 2, P1): ödevle tüm ara değerler aynı.

---

## 13. Bilinen sınırlar (kusur değil, yöntem/kapsam)

- Lineer, 1B, 1. ve 2. mertebe; Robin, periyodik ve doğrusal olmayan problemler yok.
- Birinci mertebe standart Galerkin için optimal hız ve kararlılık garantisi yok (T2, T3'te çift p'de hız p).
- IVP: koşullu kararlılık (P1: ωh < 2√3, türetildi; P2: ω²h² ≈ 10 civarında ayrık rezonans); genel yakınsama kanıtı yok.
- Taşınım baskın ve Helmholtz tipi problemlerde kaba ağda salınım/kirlilik (yöntemin bilinen davranışı).
- Eşit aralıklı düğümler; p > 4'te kötü koşullanma beklenir.
- E_max sayısal bir tahmin (alt sınır).
- İnce ağ referansının kendi hatası asimptotik bir tahmindir.

## 14. Önerilen düzeltmeler — öncelik sırasıyla (uygulanmadı)

1. **KR-1 / SK-1 / ÇK-1:** `exact_solution`, dsolve sonucunu **sınır ve başlangıç koşullarında da** sınamalı (u(xᵢ), u′(xᵢ), doğal koşulda α u′). Sağlamıyorsa analitik çözüm "yok" sayılmalı. Test: E2 durumunda `exact` = None olmalı.
2. **KR-2 / SK-2:** iyi tanımlılık ön-denetimi (iki çekirdekte aynı kurallar):
   - (a) α'nın aralıkta sıfır olması ya da işaret değiştirmesi → dejenere: iç noktada reddet; uçta ise yalnız o uç "verilmedi" iken kabul et.
   - (b) Sabit katsayılı belirsiz problemde (α, γ aynı işaret, β = 0) Dirichlet–Dirichlet için γ/α ≈ (mπ/L)² rezonans denetimi; Neumann ve karışık koşullar için karşılığı.
   - (c) Genel durum için ayrık koşul sayısı tahmini, ve ağ inceltilirken çözüm normunun büyümesi → "iyi tanımlı olmayabilir / rezonansa yakın" uyarısı.
3. **SK-3 / ÇK-6:** birinci mertebede atılan denklemin kalıntısını ve rank/uyumluluk bilgisini kâğıtta ve arayüzde göster. Seçim (ders kitabı kuralı) değişmez.
4. **MÇ-1:** tekillik mesajını nedenine göre ayır (saf Neumann / dejenere α / rezonans / IVP kaba ağ).
5. **SA-1, SA-2:** Péclet (\|β\|h/2\|α\|) ve kh (√\|γ/α\|·h) göstergeleri; eşik aşılınca uyarı.
6. **SA-3:** hız gösterilirken son iki hızın tutarlılığını denetle; değilse "ön-asimptotik" etiketi.
7. **Test paketi:** bu denetimdeki bağımsız mpmath FEM'i (yalnız SymPy ile gelen mpmath'i kullanır, yeni bağımlılık yok) ve T1–T16, E1–E5'i `tests.py`'ye kalıcı, bağımsız bir doğruluk ölçütü olarak ekle. Özellikle: analitik çözümün koşulları sağladığı, iyi tanımlı olmayan girdilerin reddedildiği, birinci mertebe uyumluluk değeri.
8. **ÇK-7:** CLAUDE.md ve README'yi kodla eşitle.
9. **V0.2 (ayrı aşama):** ÇK-2…ÇK-5 (şekil fonksiyonu türetmesi, test fonksiyonu başına denklemler, montaj toplamları, eleme adımları, tüm elemanların integralleri).

## 15. V0.1 tamamlanma değerlendirmesi

| Tanım ölçütü | Durum |
|---|---|
| Desteklenen her sınıfın açık matematiksel formülasyonu | Kısmen: formülasyonlar kodda doğru, ama dejenere ve belirsiz sınıflar için kapsam tanımı yok (§3–5) |
| P1–P4 şekil fonksiyonları bağımsız doğrulandı | **Evet** |
| Eleman matrisleri ve montaj bağımsız benchmark'ı geçiyor | **Evet** |
| Sınır koşulları matematiksel olarak tutarlı | **Hayır** (KR-1, SK-2) |
| Çözücü davranışı uygun teşhis ediliyor | **Hayır** (KR-2, MÇ-1, SA-1) |
| Hata ve yakınsama hesapları doğrulandı | Evet, KR-1 dışında |
| Desteklenmeyen durumlar açıkça belirleniyor | Kısmen: bazıları reddediliyor, bazıları sessizce kabul ediliyor |
| Kâğıt gerçek yöntemi doğru temsil ediyor | Kısmen: sayılar doğru ve tutarlı; KR-1 ve ÇK-3 var; şeffaflık V0.2'de |

**Sonuç: V0.1 tamamlanmadı.** Öncelik 1–4 uygulandığında ölçütlerin tamamı karşılanabilir durumda; hesap çekirdeğinin kendisinde düzeltilmesi gereken bir hata bulunmadı.

---

### Ek: yeniden üretim

Denetim betikleri geçici klasörde (sistem temizliğinde silinebilir):
(geçici denetim klasörü; depoya dahil değil)

| Dosya | İçerik |
|---|---|
| `ifem.py` | Bağımsız FEM |
| `cases.py`, `run_cases.py` | T1–T16 |
| `probes.py` | E1–E5 |
| `probes2.py` | Satır seçimi, akı geri hesabı, simetri, rank, koşul sayısı |
| `probes3.py`, `probes3.js` | Şekil fonksiyonları, süreklilik, hata ölçüleri, KR-1 |
| `appk.js` | Uygulama çekirdeğinden yalnızca okuma |
