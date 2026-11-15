# V0.1 Düzeltme Planı — Madde 1–4 (onaylandı; Revizyon 1 ile)

## Revizyon 1 (kullanıcı geri bildirimi; aşağıdaki ilgili kısımların yerini alır)

**R1. Sınıflandırma sisteminin boyutu.** M C = r sistemi q × m boyutludur. m, denklemin mertebesidir (genel çözümdeki
bağımsız sabit sayısı); q, problem tanımının getirdiği koşul denklemlerinin sayısıdır.
- Desteklenen kapsam: m = 1 (α ≡ 0; q = 1) ve m = 2 (q = 2: iki uçta birer koşul, ya da bir uçta u ve u′).
- 2. mertebede "verilmedi" bir uç, formülasyon gereği doğal koşul αu′ = 0 getirir ve bir denklem sayılır. IVP'nin
  serbest ucu denklem getirmez.
- Karar: rank M = rank[M|r] = m → tek; rank[M|r] > rank M → çözüm yok; rank M = rank[M|r] < m → tek değil
  (boş uzay boyutu m − rank M).
- Daha yüksek mertebe uygulamada yoktur; sınıflandırıcı m ∉ {1, 2} için "desteklenmiyor" döner.

**R2. Dejenere katsayılar: k kuralı kaldırıldı.** α'nın (1. mertebede β'nın) sıfır olması ya da işaret değiştirmesi
**tek başına ret nedeni değildir**. Sınıflar:

| Sınıf | Tanım | Sonuç |
|---|---|---|
| düzgün eliptik | kapalı aralıkta \|α\| ≥ α₀ > 0, işaret sabit | standart kuram (Lax–Milgram / Gårding) geçerli |
| dejenere | α bir noktada sıfıra değiyor, işaret değiştirmiyor (uçta ya da içte) | **uyarı**; standart H¹ kuramı ve FEM yakınsama hızları garanti değil |
| işaret değiştiren | α iç noktada işaret değiştiriyor | **uyarı**; denklem o noktada tekil |
| desteklenmeyen biçim | u′ koşulu α(x_b) = 0 olan uçta | **ret** (mevcut kural). Bu bir *formülasyon* gerçeğidir: u′ koşulu zayıf forma yalnızca Q = αu′ ile girer; α(x_b) = 0 ise koşul denkleme hiç girmez. Problemin çözümsüz olduğu iddia edilmez. |

- **Dejenere ya da işaret değiştiren durumlarda sürekli sınıflandırma:** yalnız dsolve genel çözüm verirse yapılır, ve **açık bir
  varsayımla**: aranan çözümler kapalı aralıkta sınırlıdır (FEM'in C⁰ uzayı sınırlı fonksiyonlardan oluşur).
  - Uçta limiti sonsuz olan temel çözümün katsayısı sıfırlanır; bu, sınırlılık denklemi olarak sisteme eklenir ve mesajda adıyla yazılır.
  - Örnek E2: u = x + C₁ ln x + C₂. Sınırlılık C₁ = 0 gerektirir; u(0) = 0 ⇒ C₂ = 0; u(1) = 0 ⇒ 1 = 0 → "sınırlı çözüm sınıfında çözüm yok".
- Genel çözüm yoksa ya da u_p kendisi sınırsızsa: **"sürekli sınıflandırma yapılamadı"** + uyarı; ret yok.
- Yalnızca bu varsayım altında kanıtlanan K1/K2 sonuçları kâğıdı engeller.

**R3. Sayısal eşikler: sezgiseldir, evrensel değildir; belgelenir ve ölçek testinden geçer.**
- **Rank:** σ ≤ 1e-12·σ_max (göreli, ölçekten bağımsız).
- **Uyumluluk:** normwise geriye hata (Rigal–Gaches) η = ‖b − A x̂‖ / (‖A‖‖x̂‖ + ‖b‖). x̂, en küçük kareler / en küçük norm çözümüdür. η ≤ 1e-10 → uyumlu. η, A ve b'nin ölçeklenmesinden bağımsızdır.
- **Koşul sayısı:** κ̂₁ > 1e10 uyarı, > 1e14 sayısal tekil.
- **Ölçek testi:** A, b ile 1e±8 ölçeklemesi ve N = 3…257 boyutlarında aynı karar.
- Arayüzde ve kâğıtta eşikler "sezgisel tolerans" olarak anılır.

**R4. Birinci mertebe uyumluluk: sol sıfır uzayı.** Tam sistem A a = b (A ∈ R^{N×(N−1)}).
- Uyumluluk ⇔ her y için (yᵀA = 0 ⇒ yᵀb = 0).
- Atılan satır a₁ᵀ, kalanlar A₂ ise ve A₂ tekil değilse: y = (1, −A₂^{-T}a₁) ve yᵀb = −r₁. Bu yüzden r₁ = 0 ⇔ uyumlu; ama bu denkliğin koşulu, A₂'nin tekil olmamasıdır.
- **Gösterilenler:** r₁, ölçekli |yᵀb|/(‖y‖‖b‖) ve η. A₂ tekilse SVD yoluyla genel sol sıfır uzayı kullanılır.
- Tam (rasyonel) yolda y = nullspace(Aᵀ) kesin hesaplanır.

**R5. Analitik doğrulama.**
- Denklem, Dirichlet, Neumann ve başlangıç koşulları sınanır.
- **Akı koşulu** yalnızca formülasyonun problem tanımına kattığı yerde sınanır: 2. mertebe BVP'de "verilmedi" uç (α(x_b) ≠ 0). IVP'nin serbest ucunda ve 1. mertebede sınanmaz. α(x_b) = 0 olan "verilmedi" uçta akı yerine sınırlılık (sonlu limit) sınanır; çünkü problem tanımı orada bunu gerektirir.
- **Kanıt düzeyleri ayrı gösterilir:**
  - **sembolik** — kanıt;
  - **sayısal** — 50 basamakta |değer| ≤ 1e-30·ölçek; kanıt değil, sayısal destek. Kâğıt ve arayüz bunu yazar.

**Uygulama sınırları:**
- Önce yedek ve temel test çıktısı alınır.
- Gereksiz refactoring yok; arayüz tasarımı korunur; yeni bağımlılık yok.
- Regresyon "bit düzeyinde aynı" değil, toleranslı karşılaştırılır: düğüm değerleri göreli 1e-12.
- Her madde ayrı uygulanıp test edilir; testler geçmezse sonraki maddeye geçilmez.


Kaynak: `AUDIT_V0.1.md` bulguları KR-1, KR-2, SK-1…SK-3, MÇ-1, ÇK-1, ÇK-6.

**Hedef:** yeni özellik değil; uygulamanın hangi matematiksel durumda olduğunu bilmesi ve bunu doğru söylemesi.

**Değişmeyecekler:**
- Galerkin formülasyonu ve eleman integralleri.
- Montaj ve satır seçimi kuralı (ders kitabı: u verilen uçta test fonksiyonu sıfır).
- Arayüz tasarımı.
- Çözücünün aritmetiği: elemenin işlem sırası korunur; düğüm değerleri bit düzeyinde aynı kalmalıdır, bu bir testle kanıtlanacak.

## 0. Ortak tasarım: iki düzeyli teşhis

Bir problemin durumu iki farklı düzeyde söylenebilir ve **bu iki düzey hiçbir yerde karıştırılmayacak**:

| Düzey | Neyi söyler | Nasıl karar verilir | Kim yapar |
|---|---|---|---|
| **Sürekli** | Diferansiyel problemin çözümü var mı, tek mi | Genel çözüm u = u_p + C₁e₁ + C₂e₂ (dsolve, ics olmadan) koşullara yerleştirilir → M C = r. rank M = rank[M\|r] = boyut → tek; rank[M\|r] > rank M → **çözüm yok**; eşit ama < boyut → **tek değil** (Fredholm alternatifi) | Yalnız `fem_sheet.py` (sembolik) |
| **Ayrık** | Kurulan cebirsel sistem tekil mi, kötü koşullu mu | rank, boş uzay, sağ taraf uyumluluğu (SVD ya da tam rasyonel rank), koşul sayısı | İki çekirdek, birbirinden bağımsız |

**Kurallar:**
- "Çözüm yok" ya da "tek değil" kararı yalnızca (a) sürekli düzeyin kesin sonucundan, ya da (b) ayrık sistemin rank ve uyumluluk analizinden verilir. (b)'deki karar **"ayrık sistem için"** diye etiketlenir.
- Büyüklüğe ya da büyümeye bakarak çözümsüzlük iddia edilmez.
- Sürekli karar verilemiyorsa (kapalı biçim yok, zaman aşımı, sunucu yok) uygulama bunu açıkça söyler: **"sürekli problem sınıflandırılamadı"**.

Sıralama: girdi → mertebe → **katsayı dejenereliği** → koşul sayısı kuralları → **sürekli sınıflandırma** → montaj/çözüm → **ayrık teşhis** → (1. mertebede) **atılan denklem kalıntısı** → **analitik çözüm doğrulaması**.

---

## 1. Analitik çözüm doğrulaması (KR-1, SK-1, ÇK-1)

### Matematik
Aday çözüm ũ geçerlidir ancak ve ancak aşağıdakilerin hepsi sağlanıyorsa:
- (i) R[ũ] = (αũ′)′ + βũ′ + γũ − f ≡ 0 (x₀, x₁)'de;
- (ii) her temel koşulda ũ(x_b) = g;
- (iii) her u′ koşulunda ũ′(x_b) = s;
- (iv) sınır değer probleminde koşulsuz uçta doğal koşul: lim_{x→x_b} α ũ′ = 0;
- (v) IVP'de aynı uçta hem ũ hem ũ′;
- (vi) ũ ve ũ′ kapalı aralıkta sonlu. FEM çözümü H¹'de; sınırsız bir aday karşılaştırma referansı olamaz.

Uçta değerler `limit` ile alınır; tekil bir uçta sonsuz limit, koşulun sağlanmadığı anlamına gelir.

### Karar prosedürü (her kontrol için)
1. **Sembolik:** `sp.simplify(ifade) == 0` → **sağlandı (sembolik)**.
2. **Sayısal:** sembolik karar verilemezse ifade 50 basamakla değerlendirilir (iç kontrol için 17 Çebişev noktası); |değer| ≤ 1e-30 · max(1, ölçek) → **sağlandı (sayısal, 50 basamak)**. Değilse **sağlanmadı** ve değer raporlanır. 1e-30, 50 basamaklı hesabın yuvarlama hatasının çok üstünde, gerçek bir ihlalin (bu uygulamada tipik ≥ 1e-6) çok altındadır.
3. **Zaman aşımı** → **doğrulanamadı**. Çözüm kullanılmaz.

### Etkilenen kod
| Dosya | Fonksiyon | Değişiklik |
|---|---|---|
| `fem_sheet.py` | `exact_solution` | dsolve `ics` yerine **genel çözüm + koşul sistemi** (madde 2 ile ortak hesap). Sonuç `verify_exact`'ten geçmezse kullanılmaz. Dönüş: `{u, status: verified/rejected/not_found/timeout, checks: [(ad, değer, sonuç, yöntem)], reason}` |
| `fem_sheet.py` | **yeni** `verify_exact(r, candidate)` | Yukarıdaki (i)–(vi) kontrolleri. Tek başına test edilebilir. |
| `fem_sheet.py` | `build_sheet` | Bölüm 10 yalnızca `verified` ise hata tablosuyla; ilk satırda kontrol listesi (R ≡ 0 ✓, u(0) = 0 ✓ …). `rejected` ise yalnızca "aday: …, ihlal: u(0) = −1 ≠ 0 → kullanılmadı"; hata hesabı yok. |
| `fem_sheet.py` | `error_measures` | Yalnızca `verified` çözümle çağrılır. |
| `fem_sheet.py` | `exact_for_browser`, `exact_only`, `make_sheet` | Çıktıya `status`, `reason`, `checks` eklenir. `exact_js` yalnızca `verified` ise dolu. |
| `server.py` | `compute_exact` | Mantık değişmez; yeni alanları geçirir. |
| `index.html` (UI) | `requestExact`, `parseExact` | **Bağımsız ikinci kontrol:** ayrıştırılan fonksiyon tarayıcıda koşullara karşı sayısal olarak yeniden sınanır (u(x_b), u′(x_b); tolerans 1e-9 · ölçek). Tutmazsa kullanılmaz. Yeni durum: `rejected`. |
| `index.html` (UI) | `compute`, `renderTiles` | `rejected` → referans ince ağa döner ve etiketler "fark" olur (mevcut mekanizma). Referans kutusu nedeni yazar: "analitik aday reddedildi: u(0) = −1 ≠ 0". |

---

## 2. İyi tanımlılık ve sayısal güvenilirlik (KR-2, SK-2)

### 2a. Katsayı dejenereliği (sürekli düzey; tanımlar kesin, sıfır bulma sayısal olabilir)

Önde gelen katsayı: 2. mertebede α, 1. mertebede β.

| Durum | Matematik | Karar |
|---|---|---|
| Önde gelen katsayının iç noktada sıfırı ya da işaret değişimi | Denklem o noktada tip değiştirir (tekil nokta); standart H¹ kuramı geçersiz | **Reddet** — "dejenere katsayı: α(x*) = 0, x* ≈ …" |
| α(x_b) = 0 bir uçta, α ~ \|x − x_b\|^k, **k ≥ 1** (1/α o uçta integrallenemez) | Sınırlılık bir koşulu tüketir: αu′ = F + C, sınırlı u için C belirlenir. O uçta u ya da u′ vermek problemi aşırı belirler. | O uç **"verilmedi" olmalı**; u ya da u′ verilmişse **reddet** ("güçlü dejenere uç, k = …: u verilemez, sınırlılık koşulu geçerli") |
| α(x_b) = 0, **0 < k < 1** | 1/α integrallenebilir; iki koşul da konabilir | Kabul et, **bilgi notu**: "zayıf dejenere uç (k = …); yakınsama hızı düşebilir" |
| 1. mertebede β'nın kapalı aralıkta sıfırı | Tekil başlangıç değer problemi | **Reddet** |

**Sıfır bulma:**
- **Python:** önce `solveset(α, x, Interval)` (kesin). Sonlu küme dönmezse 2001 noktalı örnekleme, ardından işaret değişiminde ikiye bölme ve |α|'nın yerel minimumlarında altın oran araması. Sıfır sayılma ölçütü: |α| ≤ 1e-12 · max|α|. k değeri `limit(log|α| / log|x − x_b|)` ile.
- **JS:** aynı sayısal yol; k için iki ölçekten eğim (δ = 1e-4·L ve 1e-6·L).

Sıfır sayısal olarak bulunduysa mesaj bunu söyler ("sayısal tespit").

### 2b. Sürekli sınıflandırma (yalnız `fem_sheet.py`)
- Genel çözüm + koşul sistemi M C = r (madde 0).
- **det M kararı:** sembolik sıfırsa tekil. Sembolik karar verilemezse 50 basamakta |det M| > 1e-30 · ölçek → tekil değil (sayısal sertifika). Arada kalırsa → "belirlenemedi".
- **Uyumluluk:** rank[M|r] (aynı ölçütlerle).
- **Sınıflar:** `tek` · `çözüm yok` · `tek değil (boş uzay: u₀(x) = …)` · `belirlenemedi (neden)`.
- `çözüm yok` ve `tek değil` → **kâğıt üretilmez**, Türkçe mesaj. Tek değilse boş uzay fonksiyonu da yazılır (örn. sin πx, sabit).
- **IVP:** α ≠ 0 kapalı aralıkta ise sürekli problem tektir (varlık-teklik teoremi); yine de aynı hesapla doğrulanır.
- Sonuç `/api/exact` ile tarayıcıya gider; tarayıcı `çözüm yok` ya da `tek değil` görünce mevcut kırmızı hata yolunu kullanır (grafikler soluklaşır).

### 2c. Ayrık teşhis (iki çekirdek, bağımsız)

**Eşikler (gerekçeli):**
- **Rank:** SVD'de σ_i ≤ τ = 1e-12 · σ_max ise sıfır. 64-bit SVD'nin geriye hata sınırı ~N·ε·σ_max ≈ 6e-14·σ_max (N ≤ 257); τ bunun üstünde.
- **Uyumluluk:** küçük σ'lara ait sol tekil vektörler u_k için |u_kᵀ b| ≤ 1e-9 · ‖b‖₂ → uyumlu. Bu tolerans kuadratür hatasını (göreli ≲ 1e-14) karşılar, gerçek bir uyumsuzluğu gizlemez.
- **Koşul sayısı:**
  - κ̂₁ > 1e10 → "kötü koşullu: göreli hata sınırı ≈ κ̂·ε = …" uyarısı.
  - κ̂₁ > 1e14 ya da pivot çöküşü → ayrık sistem **sayısal olarak tekil** sayılır, rank ve boş uzay analizine geçilir.

**JS (`index.html`, CORE):**
- **`solveLinear`:** eleme aritmetiği ve pivotlama **aynı kalır**; ek olarak LU çarpanlarını ve permütasyonu saklar ve döndürür. Pivot çöküşünde genel mesaj yerine tipli bir hata fırlatır.
- **Yeni `conditionEstimate(lu)`:** Hager–Higham 1-norm tahmini; mevcut LU ile birkaç ek çözüm, O(N²). Her çözümde hesaplanır (canlı sürgüye uygun).
- **Yeni `analyzeSystem(A, b)`:** tek yönlü Jacobi SVD. Yalnızca tetiklendiğinde çalışır (pivot çöküşü ya da κ̂ > 1e14); O(N³), N = 257'de ~0.5 s, nadir yol. Rank, boş uzay vektörü ve uyumluluk döndürür.
- **`METHODS.galerkin.solve`:** sözleşme genişler (geriye uyumlu): `{ c, K, F, quadratureError, diagnostics: { condition, droppedResidual?, … } }`.
- **`buildProblem`:** 2a kuralları.

**Python (`fem_sheet.py`, `solve_problem`):**
- **Tam yol (≤ 8 serbest bilinmeyen):** rasyonel `rank`, `nullspace`, uyumluluk (tam, toleranssız).
- **Sayısal yol:** mevcut `mpmath.cond` korunur. Tetiklenirse ve N_free ≤ 65 ise `mpmath.svd_r` ile rank ve boş uzay analizi. Daha büyük sistemlerde: "sayısal olarak tekil (κ > 1e13); bu boyutta rank analizi yapılmadı" (dürüst sınır).

**Sürekli ve ayrık sonucun birlikte yorumu:**

| Sürekli | Ayrık | Kullanıcıya |
|---|---|---|
| tek | iyi koşullu | normal |
| tek | kötü koşullu / tekil | "Sürekli problem tek çözümlü; ayrık sistem kötü koşullu ya da tekil → ağ ya da derece bu problem için uygun değil" (ör. IVP'de ωh > 2√3, Helmholtz kirliliği) |
| çözüm yok / tek değil | — | ret (2b) |
| belirlenemedi | iyi koşullu | normal + not: "sürekli problem sınıflandırılamadı" |
| belirlenemedi | kötü koşullu / tekil | ayrık sonuç + "kesin teşhis yok: rezonans, dejenerelik ya da ağ yetersizliği olabilir" (iddia yok) |

### Etkilenen kod
- `fem_sheet.py`: **yeni** `coefficient_degeneracy`, **yeni** `classify_problem`; `solve_problem` (baş kısımda 2a/2b çağrıları, sonda ayrık teşhis); `build_sheet` ("Verilenler"e tek satır: sınıf ve kanıt, ör. "rank M = 2 = rank[M|r] → tek çözüm"); `make_sheet`, `exact_only`.
- `index.html`:
  - CORE: `buildProblem`, `solveLinear`, yeni `conditionEstimate`, yeni `analyzeSystem`, `METHODS.galerkin.solve`.
  - UI: `compute`, `render` ve teşhis kutusu. Kutu, mevcut `.stale` stili (amber kenar) ile çiziliyor; yeni görsel dil yok.
- `server.py`: değişmez (alanları geçirir).

---

## 3. Birinci mertebede satır seçimi ve uyumluluk (SK-3, ÇK-6)

### Matematik
- **Deneme uzayı:** u_h(x₀) = u₀; N − 1 bilinmeyen.
- **Test uzayı:** V_h,0 = {v : v(x₀) = 0}; denklemler i = 2…N. Ders kitabı kuralı korunur.
- **Tam sistem** (i = 1…N) N × (N − 1) boyutludur.
- **Atılan denklem kalıntısı** r₁ = Σⱼ K₁ⱼaⱼ − F₁. Seçilen denklemler tam sağlandığı için tam sistem uyumludur ⇔ **r₁ = 0**. γ = 0'da bu, sol sıfır vektörüyle wᵀF = 0 koşuluna denktir.
- **r₁ ≠ 0 ise:** başka bir satırı atmak başka bir test uzayı (Petrov–Galerkin) demektir ve **başka bir a** verir. Bütün seçimler tutarlıdır (gerçek çözüm her satırı sağlar), ama ayrık çözümleri farklıdır. Hangisinin daha doğru çıktığı probleme bağlıdır (u′ = 3x²'de alternatif, u′ + 2u = 0'da ders kitabı seçimi daha iyi). Bu, "ilk satırı atmak her zaman doğru sonuç verir" anlamına gelmez.
- **Uyumluluk kararı:** tam yolda r₁ == 0 (kesin); sayısal yolda |r₁| ≤ 1e-10 · (|F₁| + Σ|K₁ⱼaⱼ|).

### Gösterim
- **Kâğıt, bölüm 7:**
  - v = Nᵢ, i = 2…N (N₁ ∉ V_h,0, çünkü u(x₀) verildi)
  - r₁ = Σ K₁ⱼaⱼ − F₁ = …
  - "r₁ = 0 ⇒ tam sistem uyumlu" ya da "r₁ ≠ 0 ⇒ tam sistem uyumsuz; başka satır seçimi başka a verir".
- **Arayüz teşhis kutusu (bilgi düzeyi):** "1. mertebe: atılan denklem kalıntısı r₁ = …; tam sistem uyumlu/uyumsuz".
- **Çapraz kontrol:** r₁ iki çekirdekte karşılaştırılır.
- **README:** yukarıdaki açıklama ve iki örnek.

### Etkilenen kod
- `fem_sheet.py`: `solve_problem` (r₁ ve uyumluluk), `build_sheet`.
- `index.html`: `METHODS.galerkin.solve` (`diagnostics.droppedResidual`), UI teşhis kutusu ve `crossCheck`.

Çözüm değişmez.

---

## 4. Ayrıştırılmış tekillik ve hata mesajları (MÇ-1)

Tek genel `SINGULAR_MESSAGE` kaldırılır. Her mesaj **kategori + düzey (sürekli/ayrık) + kanıt** taşır:

| Kategori | Ne zaman | Örnek mesaj |
|---|---|---|
| **K1. Çözüm yok** (uyumsuz) | Sürekli: rank[M\|r] > rank M. Ayrık: sistem tekil ve uyumsuz | "Çözüm yok: sınır koşulları denklemle uyumsuz (sürekli problem; rank M = 1 < rank[M\|r] = 2). Örn. rezonans: ∫f·sin(πx)dx ≠ 0." |
| **K2. Tek olmayan çözüm** (boş uzay ≠ {0}) | Sürekli: rank M = rank[M\|r] < 2. Ayrık: tekil ve uyumlu | "Çözüm tek değil: u + c·sin(πx) de çözüm (sürekli problem). Bir uçta u değeri ver ya da γ'yı değiştir." |
| **K3. Neredeyse tekil / kötü koşullu** | κ̂ > 1e10 (uyarı) ya da ayrık tekil ama sürekli tek | "Ayrık sistem kötü koşullu (κ ≈ 3·10¹¹): sonuçta ~10⁻⁵ göreli hata olabilir." / "Sürekli problem tek çözümlü ama ayrık sistem tekil: ağ bu problem için yetersiz (IVP: ωh < 2√3 gerekir)." |
| **K4. Dejenere katsayı** | 2a | "Dejenere katsayı: α(x) = x − 0.5, x* = 0.5'te sıfır ve işaret değiştiriyor; bu sınıf desteklenmiyor." |
| **K5. Sayısal çözücü hatası** | NaN/Inf, sınıflandırılamayan pivot çöküşü, SVD yakınsamaması, zaman aşımı | "Sayısal çözücü hatası: … (problem matematiksel olarak sınıflandırılmadı)." |

Mevcut IVP "ağ fazla kaba" mesajı K3'ün özel hali olur; içeriği korunur.

### Etkilenen kod
- `index.html` CORE: sabitler (`MESSAGES` tablosu), `solveLinear`, `METHODS.galerkin.solve`, `buildProblem`.
- `fem_sheet.py`: `solve_problem` (`singular` metinleri), `UserError` mesajlarına kategori öneki.
- `tests.js`, `tests.py`: mesaj kategorisi testleri.

---

## 5. Kalıcı test paketi

Denetim betikleri geçici klasörden projeye taşınır:

| Yeni dosya | İçerik |
|---|---|
| `verification/ifem.py` | Bağımsız FEM (mpmath; uygulama kodunu içe aktarmaz). Mevcut denetim sürümü, yorumlarla. |
| `verification/benchmarks.py` | T1–T16 ve E1–E5 tanımları: problem, analitik çözüm, beklenen davranış, beklenen sınıf |
| `verification/run_verification.py` | Rapor tablosu (bölüm 12 biçiminde) ve çıkış kodu. `--quick` (n = 4, 8; ~1 dk) ve `--full` (n = 4…32; birkaç dk). |

- `tests.py`, `verification` paketinin `--quick` alt kümesini çağırır.
- Yeni bağımlılık yok: mpmath zaten SymPy ile geliyor.
- Uygulamanın sonuçları `node tests.js --solve` ve `fem_sheet` üzerinden alınır; beklenen değerler **analitik çözümden ve bağımsız gerçekleştirmeden** gelir, uygulamanın kendi çıktısından değil.

## 6. Test senaryoları ve beklenen sonuçlar

### Madde 1
| # | Senaryo | Beklenen |
|---|---|---|
| 1.1 | `verify_exact`: u″ = 2, u(0) = u(1) = 0, aday x² | rejected: u(1) = 1 ≠ 0 |
| 1.2 | aday kalıntıyı sağlamıyor (u = x³ için u″ = 2) | rejected: R ≠ 0 |
| 1.3 | Neumann ihlali: u″ = 0, u(0) = 0, u′(1) = 1, aday 0 | rejected: u′(1) = 0 ≠ 1 |
| 1.4 | doğal koşul ihlali: u″ = 0, sol verilmedi, u(1) = 1, aday x | rejected: α u′(0) = 1 ≠ 0 |
| 1.5 | IVP eğim ihlali: u″ + u = 0, u(0) = 0, u′(0) = 1, aday 2 sin x | rejected: u′(0) = 2 |
| 1.6 | E2 (α = x, iki uçta u) | sheet: analitik yok; tarayıcı `exact_js` = null (ayrıca 2a ile zaten reddedilir) |
| 1.7 | T1–T16 analitik çözümleri | hepsi verified; kontrol listesi kâğıtta |
| 1.8 | JS `parseExact` sahte yük: u(x₀) koşulu tutmayan aday | null |
| 1.9 | Arayüz: rejected durumda referans kutusu "ince ağ" + gerekçe; etiketler "D (fark)" | tarayıcıda doğrulanır |

### Madde 2
| # | Senaryo | Beklenen |
|---|---|---|
| 2.1 | E1: u″ + π²u = sin πx, D/D | K1, sürekli, kesin (rank 1 < 2) |
| 2.2 | E1b: u″ + π²u = 0, D/D | K2, boş uzay sin πx |
| 2.3 | u″ + π²u = sin 2πx, D/D | K2 (uyumlu rezonans) |
| 2.4 | γ = π² + 10⁻⁶, f = sin πx | **tek**, ret yok; ayrık κ raporlanır; "çözümsüz" denmez |
| 2.5 | u″ = 2x − 1, N/N | K2, boş uzay sabit |
| 2.6 | u″ = 1, N/N | K1 |
| 2.7 | α = x, x = 0'da u verili | K4, güçlü dejenere (k = 1) |
| 2.8 | α = x, x = 0 verilmedi (T12, T13) | kabul, sonuçlar değişmez |
| 2.9 | α = √x, D/D | kabul + "zayıf dejenere (k = 1/2)" notu |
| 2.10 | α = x − 0.5 | K4, iç sıfır x* ≈ 0.5 (n tek ya da çift fark etmez) |
| 2.11 | α = (x − 0.3)² (değme sıfırı) | K4, x* ≈ 0.3 |
| 2.12 | 1. mertebe β = x − 0.5 | K4 |
| 2.13 | IVP u″ + u = 0, [0, 64·3.6], n = 64, P1 | K3: "sürekli tek, ayrık tekil/kötü koşullu, ağ yetersiz" |
| 2.14 | T7 (değişken katsayı, kapalı biçim yok) | "sürekli sınıflandırılamadı" notu; çözüm normal |
| 2.15 | T1–T16 | hepsi "tek"; ret yok; düğüm değerleri değişiklik öncesiyle **bit düzeyinde aynı** (temel dosyayla karşılaştırma) |
| 2.16 | `analyzeSystem` birim testleri | [[1,1],[1,1]], b = [1,1] → rank 1, uyumlu; b = [1,2] → uyumsuz; I → rank 2 |
| 2.17 | `conditionEstimate` | Hilbert 6×6 ve köşegen matrislerde bilinen κ₁ ile 3 kat içinde |

### Madde 3
| # | Senaryo | Beklenen |
|---|---|---|
| 3.1 | u′ = 2x, [0,2], n = 2, P1 (ödev) | r₁ = 0 kesin, "uyumlu" |
| 3.2 | u′ = 3x², [0,2], n = 2, P1 | r₁ = −1 (kesin rasyonel), "uyumsuz" |
| 3.3 | u′ + 2u = 0, [0,2], n = 8 | JS ve Python r₁ aynı (1e-12 göreli), bağımsız hesapla aynı |
| 3.4 | Kâğıt bölüm 7 | r₁ satırı ve uyumluluk ifadesi var; 2. mertebede bu satır yok |
| 3.5 | Çözüm değişmedi | 1. mertebe testlerinde a değerleri temel dosyayla aynı |

### Madde 4
| # | Senaryo | Beklenen kategori |
|---|---|---|
| 4.1 | 2.1, 2.6 | K1 |
| 4.2 | 2.2, 2.3, 2.5 | K2 |
| 4.3 | 2.4, 2.13 | K3 |
| 4.4 | 2.7, 2.10, 2.11, 2.12 | K4 |
| 4.5 | `analyzeSystem` NaN girdili matris | K5 |
| 4.6 | Kaynakta genel "Sistem tekil … En az bir uçta u değeri ver" metni kalmadı | grep testi |

## 7. Riskler ve dürüst sınırlar
- **Sürekli sınıflandırma** yalnızca SymPy'nin kapalı biçimde genel çözüm bulduğu problemlerde kesindir. Değişken katsayılı problemlerin çoğunda "belirlenemedi" denecek. Bu bir eksiklik değil, doğru bir beyan.
- **Sıfır bulma:** sembolik bulunamazsa sayısaldır; çok dar bir değme sıfırı örneklemeyi kaçırabilir (altın oran araması bu riski azaltır, sıfırlamaz).
- **Performans:** sunucu tarafında dsolve şimdiki kadar sürer (tek çağrı; zaman sınırı aynı). Tarayıcıda her çözümde O(N²) koşul tahmini eklenir; SVD yalnızca nadir yolda çalışır.
- **Çekirdek değişikliği:** `solveLinear` yalnızca çarpanları saklayacak şekilde değişir. Çözümlerin bit düzeyinde aynı kaldığı 2.15 ve 3.5 ile kanıtlanır.

## 8. Uygulama sırası
1. Temel kayıt: değişiklik öncesi düğüm değerleri (T1–T16, tüm p, n = 4…32) bir dosyaya.
2. `verification/` paketi (taşıma, davranış değişmez).
3. Madde 4 iskeleti (mesaj tablosu, tipli hatalar).
4. Madde 2a, 2b (Python), 2c (iki çekirdek).
5. Madde 1.
6. Madde 3.
7. Testler: `tests.js`, `tests.py`, `verification --full`; tarayıcı kontrolü; belgeler (README, CLAUDE.md) yalnızca bu maddeler için.
