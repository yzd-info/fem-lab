# FEM Lab

Yazdığın tek boyutlu, lineer diferansiyel denklemi Galerkin sonlu eleman yöntemiyle çözer.
Eleman sayısını ve polinom derecesini sürükledikçe çözüm, noktasal hata ve yakınsama grafikleri anında güncellenir.
"Çözüm kâğıdını yaz" düğmesi çözümü, anlatım cümlesi olmadan, adım adım matematik olarak bir kâğıda döker;
kâğıt PDF ya da `.tex` olarak indirilebilir.

Hazır örnek yoktur: denklemi, aralığı ve uç koşullarını kendin yazarsın.

## Kurulum (macOS)

macOS'ta `python3` komutu Xcode Komut Satırı Araçları ile gelir (`python` adı yoktur; her zaman `python3` yaz).
Terminal'i aç ve proje klasörüne geç:

```bash
cd ~/Desktop/SCIENCE/Computation/FEM/galerkin-fem-lab/galerkin-fem-lab
```

Tek bağımlılığı (SymPy) kur:

```bash
python3 -m pip install sympy
```

`python3` bulunamazsa macOS, Komut Satırı Araçlarını kurmayı önerir; kabul et ve komutu yeniden çalıştır.

## Çalıştırma

```bash
python3 server.py
```

Terminalde `FEM Lab: http://localhost:8000` satırı görünür. Tarayıcıda
http://localhost:8000 adresini aç. Terminal penceresini açık bırak; kapatırsan ya da Ctrl+C'ye basarsan sunucu durur.

- İlk çalıştırmada macOS, Terminal'in Masaüstü klasörüne erişmesi için izin isteyebilir; izin ver.
- PDF indirme için `pdflatex` gerekir (MacTeX / TeX Live). Kurulu değilse "PDF indir" düğmesi tarayıcının
  yazdırma penceresini açar; oradan "PDF olarak kaydet" seç.
- `index.html` dosyasını doğrudan açarsan grafikler çalışır, yalnız çözüm kâğıdı üretilemez (sunucu gerekir).

Windows'ta aynı komutlar `python3` yerine `py` (ya da `python`) ile yazılır.

## Testler (isteğe bağlı)

Uygulamayı çalıştırmak için Node gerekmez; yalnızca testler için gerekir (`brew install node`).

```bash
node tests.js
```

```bash
python3 tests.py
```

`tests.js` tarayıcıdaki sayısal çekirdeği (161 kontrol), `tests.py` çözüm kâğıdını (186 kontrol) sınar:
SymPy sonucu ↔ JavaScript sonucu ↔ analitik çözüm ↔ LaTeX derlemesi. `pdflatex` yoksa derleme kontrolleri atlanır.

## Kullanım

1. **Denklem:** `(α(x) u')' + β(x) u' + γ(x) u = f(x)` biçimindeki dört katsayıyı yaz. Boş kutu 0 sayılır.
   Yazılabilenler: `x`, `pi`, `e`, `+ - * / ^`, `sin cos tan asin acos atan sinh cosh tanh exp ln sqrt abs`,
   `2x` gibi örtük çarpım. Örnekler: `u'' + u = x` → α = `1`, γ = `1`, f = `x`;
   `-u'' = sin(πx)` → α = `-1`, f = `sin(pi*x)`; `u' + 2u = 0` → α boş, β = `1`, γ = `2`.
2. **Aralık:** x₀ ve x₁.
3. **Uç koşulları:** her uç için tür (`u =`, `u′ =`, `α u′ =`, `u ve u′`, `verilmedi`) ve değer seç.
   - `u′ =` **türev** koşuludur, `α u′ =` **akı** koşuludur; birbirine çevrilmez. Ör. α = 1 + x iken `u′(1) = 4`
     akı Q = α(1)·4 = 8 demektir; `(1+x)u′(1) = 4` istiyorsan `α u′ =` seçip `4` yaz (Q = 4).
   - Sınır değer problemi: iki koşulu farklı uçlara ver, ör. `u(0) = 0`, `u(1) = 0`.
   - Başlangıç değer problemi: bir uçta `u ve u′` seç (konum ve hız), diğer ucu `verilmedi` bırak.
     Ör. yay: α = `1`, γ = `1`, x₀'da u = `0`, u′ = `1`, aralık [0, 10] → u = sin x.
   - Birinci mertebe: `dy/dx = x^2`, y(0)=0 → α boş, β = `1`, f = `x^2`, x₀'da `u =` `0`.
4. **n ve p sürgüleri:** eleman sayısı 1–64, polinom derecesi 1–4. n sürgüsü odaktayken ok tuşları, Home/End
   da çalışır; sürgünün altındaki 1 / 8 / 16 / 32 / 64 etiketlerine tıklayarak doğrudan o değere geçebilirsin.
5. **Çipler:** hangi derecelerin aynı grafikte üst üste görüneceğini seç. Kesikli siyah/beyaz çizgi referanstır.
6. **Çözüm kâğıdını aç / kapat:** kâğıt sayfada sürekli durmaz. "Çözüm kâğıdı" başlığına ya da sol taraftaki
   "Çözüm kâğıdını aç" düğmesine tıkla; aynı yere yeniden tıklayınca kapanır. Açıkken denklemi, koşulları, n ya da p'yi
   değiştirdiğinde kâğıt kendiliğinden (≈ 0,7 s bekleyip) yeniden hesaplanır; kapalıyken sunucuya istek gitmez.
   **PDF indir** ve **.tex indir** düğmeleri kâğıt açıkken başlığın sağındadır.

### Referans: analitik mi, ince ağ mı?
Hata ölçmek için bir "doğru" çözüm gerekir. Sunucu açıksa her problem için SymPy ile **analitik çözüm** aranır
(`/api/exact`); bulunursa tüm hatalar ona göre ölçülür ve "Referans u: analitik" kutusu görünür. Bulunamazsa (ya da
sayfa sunucusuz açıldıysa) aynı yöntemin ince ağ çözümü (n = 128, p = 4) kullanılır: "Referans u_ref: ince ağ".
Bu ikinci durumda grafikler yöntemin kendisindeki sistematik bir hatayı **gösteremez**; bağımsız kontrol, kâğıttaki
SymPy çapraz kontrolü ve analitik tablodur.

### İyi tanımlılık ve teşhis (V0.1 düzeltmeleri)
Uygulama bir problemin durumunu **iki ayrı düzeyde** söyler ve bunları karıştırmaz:

- **Sürekli problem** (sunucu, SymPy): genel çözüm u = u_p + Σ Cᵢeᵢ (dsolve, `ics` olmadan) koşullara yerleştirilir,
  M C = r (q × m; m = mertebe). rank M = rank[M|r] = m → tek çözüm; rank[M|r] > rank M → **çözüm yok**; eşit ama < m →
  **tek değil** (Fredholm alternatifi). Yalnız **sembolik olarak kanıtlanan** "çözüm yok / tek değil" kâğıdı engeller.
  Sembolik karar 50 basamaklı sayısal rankla karşılaştırılır; uyuşmazsa sonuç "sayısal (kanıt değil)" diye yazılır.
  Kapalı biçimde genel çözüm yoksa (çoğu değişken katsayı) "sürekli problem sınıflandırılamadı" denir; ret yok.
- **Ayrık sistem** (iki çekirdek ayrı ayrı): κ₁ tahmini (Hager–Higham) her çözümde; pivot çöküşü ya da κ̂ > 1e14'te
  SVD ile rank, boş uzay ve uyumluluk (normwise geriye hata η).

| Kategori | Ne zaman |
|---|---|
| K1 Çözüm yok (uyumsuz) | sürekli: rank[M\|r] > rank M · ayrık: tekil ve η > 1e-10 |
| K2 Tek olmayan çözüm | sürekli: rank M = rank[M\|r] < m · ayrık: tekil ve uyumlu |
| K3 Neredeyse tekil / kötü koşullu | κ̂ > 1e10 (uyarı) · ayrık tekil ama sürekli problem tek (ör. IVP'de kaba ağ) |
| K4 Dejenere katsayı | α (1. mertebede β) aralıkta sıfır ya da işaret değiştiriyor: **uyarı**, ret değil |
| K5 Sayısal çözücü hatası | NaN/∞, sınıflandırılamayan çöküş |

**Sezgisel toleranslar** (evrensel matematiksel sınır değildir; göreli tanımlıdır, 10⁻⁸…10⁸ ölçeklemesinde ve
N = 3…257 boyutlarında aynı kararı verdiği test edildi): rank σ ≤ 1e-12·σ_max; uyumluluk η ≤ 1e-10; κ uyarısı 1e10,
sayısal tekil 1e14.

**Dejenere katsayı.** α'nın sıfırı ya da işaret değişimi tek başına çözümsüzlük değildir. Dejenere bir **uçta** sürekli
sınıflandırma ancak şu açık varsayımla yapılır: çözümler kapalı aralıkta sınırlıdır (FEM'in C⁰ uzayı); o uçta sonsuza
giden temel çözümün katsayısı sıfırlanır ve bu mesajda yazılır (ör. (xu′)′ = 1, u(0) = u(1) = 0 → "sınırlı çözüm sınıfında
çözüm yok"). **İç** sıfırda sınıflandırma yapılmaz ("sınıflandırılamadı"). α(x_b) = 0 olan uçta u′ koşulu ise
formülasyon gereği reddedilir (koşul zayıf forma yalnız Q = αu′ ile girer).

**Analitik çözüm doğrulaması.** Aday çözüm kullanılmadan önce: denklem kalıntısı R[u] ≡ 0, her u ve u′ koşulu, ve yalnız
formülasyonun gerektirdiği yerde akı koşulu (2. mertebe sınır değer probleminde "verilmedi" uç; IVP'nin serbest ucunda
ve 1. mertebede sınanmaz; dejenere serbest uçta sınırlılık sınanır). Kanıt düzeyi ayrı yazılır: **sembolik doğrulandı**
(kanıt) ya da **sayısal doğrulandı** (50 basamak; kanıt değil). Tarayıcı ayrıca koşulları kendi ayrıştırıcısıyla yeniden
sınar. Geçersiz aday hata ölçülerinde kullanılmaz; referans ince ağa döner ve nedeni yazılır.

**Birinci mertebe uyumluluk.** Test uzayı ders kitabındaki gibidir (u verilen düğümde v = 0; o denklem atılır). Ama tam
sistem A a = b (N × (N−1)) genelde aşırı belirlidir. Uygulama atılan denklemin kalıntısı r ile sol sıfır uzayı ölçüsü
|yᵀb|/(‖y‖‖b‖)'yi (yᵀA = 0) gösterir; A₂ tekil değilse yᵀb = −r. r = 0 → tam sistem uyumlu, satır seçimi sonucu
değiştirmez (ör. u′ = 2x). r ≠ 0 → başka bir satırın atılması başka bir test uzayı, dolayısıyla başka bir çözüm demektir
(ör. u′ = 3x², [0,2], n = 2: r = −1). Hangisinin daha doğru çıktığı probleme bağlıdır; bu bir doğruluk kanıtı değildir.

### Hata ölçüleri (özet kutuları)
| Ad | Tanım | Not |
|---|---|---|
| E_L² | (∫ (u − u_h)² dx)^½ | gerçek FEM yeniden kurulumuyla, eleman başına 96 Gauss noktası |
| \|e\|_H¹ | (∫ (u′ − u_h′)² dx)^½ | yarı-norm (türev hatası) |
| E_max | max \|u − u_h\| | 32 eşit nokta/eleman + en büyük örneğin çevresinde altın oran araması; sayısal tahmin (alt sınır), kesin maksimum değil |
| E_düğüm | max_i \|u(x_i) − a_i\| | 1D'de düğümlerde süperyakınsama olabilir; E_max ile karıştırma |
| Ölçülen hız | r = log(E_kaba/E_ince) / log(h_kaba/h_ince) | gerçek h ile; tabandaki noktalar katılmaz |

Referans ince ağ olduğunda bu ölçüler gerçek hata değil, **fark** olarak etiketlenir (E yerine D: D_L², D_max …):
u_ref − u_h, referansın kendi hatası kadar belirsizdir. Referans kutusu bu hatanın asimptotik tahminini
(E(n = 64, p = 4)/31) gösterir; fark bunun 10 katından küçükse "güvenilmez" uyarısı çıkar.

Kuramsal beklenti (düzgün çözüm, 2. mertebe sınır değer problemi): ‖e‖_L² = O(h^{p+1}), |e|_H¹ = O(h^p).
Bunlar asimptotik beklentilerdir; uygulama hızı ölçer, zorlamaz. Birinci mertebe Galerkin için optimal hız
garantisi yoktur; başlangıç değer problemindeki p + 1 / p hızları kanıt değil, testlerde gözlenen değerlerdir.

### Grafikler
- **Çözüm:** referans (kesikli) ve seçili derecelerin u_h eğrileri. Her eğri kendi elemanlarında kendi şekil
  fonksiyonlarıyla hesaplanır (yumuşatma yok); çizim örnekleri ağdan bağımsızdır ve eleman sınırlarını içerir.
  Her derece kendi düğümlerini taşır (N = np + 1): **●** eleman uç düğümü, **○** eleman içi düğüm (P2 orta nokta,
  P3 üçte birler). Aktif derece büyük işaretli. Kesikli dikey çizgiler eleman sınırlarıdır.
- **Noktasal hata:** "işaretli" görünüm u − u_h (doğrusal eksen); "mutlak (log)" görünüm |u − u_h| ve düğüm
  hataları (noktalar). Taralı bölge ölçüm tabanıdır.
- **Yakınsama, h / p seçici:** "h" = h-inceltme (derece sabit, n değişir); "p" = p-inceltme (n sabit, derece 1 → 4,
  yalnız ölçülen değerler, hız hesaplanmaz çünkü h değişmez). hp-inceltme (ikisi birden) gösterilmez.
- **Yakınsama (L² / H¹), h görünümü:** iki eksen log. ● hıza katılan ölçüm, ○ ölçüm tabanında (katılmaz).
  Satır sonu: ölçülen hız, tanımsızsa "hız —" (asla 0 gösterilmez). Kesikli gri çizgi: kuramsal eğim (yalnız
  2. mertebe sınır değer probleminde). Ölçüm tabanı: göreli 1e-11 yuvarlama (H¹'de 1/h ile büyür); ince ağ
  referansı kullanılıyorsa referansın tahmini kendi hatası (E(64, p=4)/31, 3 kat payla) da eklenir.

### Çapraz kontrol (kâğıdın üstünde)
Tarayıcı (kuadratürlü, 64-bit) ile SymPy (tam/sembolik) ayrı ayrı karşılaştırılır: düğüm katsayıları a, montajlı K
matrisi (≤ 65 bilinmeyen), u = değer koşulunun sağlanması. Tolerans 1e-9 × ölçek; katsayılar köşeli ve kuadratür
uyarısı varsa 1e-5 (bu durumda yazılır). Analitik çözüme göre düğüm hatası ayrıca yazılır: o bir uyuşmazlık değil,
yöntemin kendi hatasıdır.

## Dosya haritası

| Dosya | Görev |
|---|---|
| `index.html` | Arayüz + canlı grafikler için **sayısal** FEM çekirdeği (JavaScript, bağımlılıksız; MathJax CDN'den) |
| `fem_sheet.py` | Katman 1–3 ve 5: problem tanımı, problem doğrulaması, **sembolik/tam** FEM motoru (SymPy), çözüm kâğıdı (yalnız biçimleme) |
| `fem_verify.py` | Katman 4: bağımsız doğrulama (analitik çözüm, yama testi, bağımsız yeniden montaj, çözücü kalıntısı, dağılımsal kalıntı, hata ve yakınsama, kararlılık göstergeleri) |
| `server.py` | Yerel sunucu (yalnızca standart kütüphane): sayfayı sunar, kâğıdı hesaplatır, PDF'e derler |
| `tests.js` | Tarayıcı çekirdeği testleri (`node tests.js`) |
| `tests.py` | Kâğıt testleri (`python3 tests.py`) |
| `requirements.txt` | `sympy` |
| `verification/` | Bağımsız doğrulama: `ifem.py` (uygulamayı içe aktarmayan mpmath FEM), `benchmarks.py` (T1–T16, E1–E4b), `run_verification.py` (`--quick`, `--full`, `--compare`), `baseline/` (değişiklik öncesi temel kayıt) |
| `backup/` | Değişiklik öncesi kaynak yedeği (SHA-256 listesiyle) |
| `CLAUDE.md` | Projenin görev ve kural dosyası |

İki çekirdek bilinçli olarak bağımsızdır: JavaScript anlık grafikler için sayısal çözer, SymPy kâğıt için tam
çözer; her kâğıtta ikisinin düğüm değerleri karşılaştırılır.

## Mimari: beş katman

Veri yalnız aşağı doğru akar; her katman bir öncekinin **yapılandırılmış** çıktısını alır. Bir katman başarısız olursa
sonrakiler çalışmaz ve neden Türkçe bildirilir (kâğıt üretilmez, desteklenmeyen bir iddia yazılmaz).

| Katman | Nerede | Sorumluluk | Ne yapmaz |
|---|---|---|---|
| 1 Problem tanımı | `fem_sheet.define_problem` | Kullanıcının metnini olduğu gibi saklar (tek kaynak) | Yorumlamaz, düzeltmez |
| 2 Problem doğrulaması | `fem_sheet.validate_problem` | Her alanı iki bağımsız ayrıştırıcıdan geçirir (simge simge çevirici + SymPy dilbilgisi), eşdeğerliği iz tablosuna yazar; koşulların anlamını (temel / türev / akı / başlangıç) ve akı kuralını açıkça kaydeder; kural ihlalinde durdurur. Tarayıcı üçüncü ayrıştırıcıyla aynı 9 noktayı karşılaştırır | Çözümün varlığını/tekliğini sormaz; girdiyi sessizce değiştirmez |
| 3 FEM motoru | `fem_sheet.fem_engine` | Şekil fonksiyonu kurulumu ve sınanması, Galerkin formülasyonu, eleman integralleri, montaj, sınır koşulları, çözüm, yeniden kurulum | LaTeX metni üretmez |
| 4 Bağımsız doğrulama | `fem_verify.verify` | Motorun fonksiyonlarını çağırmadan çıktısını sınar (aşağıda) | Motorun "çözdüm" demesine güvenmez |
| 5 Çözüm kâğıdı | `fem_sheet.build_sheet` | Yalnız 2–4'ün verisini biçimler | Hiçbir sayı hesaplamaz |

Katman 4'ün denetimleri (her biri `geçti` / `başarısız` / `uyarı` / `yapılmadı` durumuyla kâğıtta ve arayüzde listelenir):

- **Analitik çözüm:** genel çözüm + koşul sistemi; aday denklemle ve **her** koşulla (u, u′, α u′) sınanır, geçmezse kullanılmaz.
- **Operatör tutarlılığı (yama testi):** V_h'deki uₚ = 1 + t + … + tᵖ için K a, kullanıcının operatörü L(u)'nun doğrudan
  türevle alınıp 24 noktalı Gauss ile integrallenmesine eşit olmalı (terim düşmesi, işaret, değişken katsayı, montaj).
- **Bağımsız yeniden montaj ve çözüm:** `verification/ifem.py` (Vandermonde tabanı, Gauss–Legendre, 30 basamak LU);
  ≤ 129 bilinmeyende yapılır, üstünde "yapılmadı" yazılır.
- **Çözücü kalıntısı:** ‖F − K a‖∞, geriye hata η = ‖F − K a‖/(‖K‖‖a‖ + ‖F‖) ve κ∞, 30 basamakta.
- **Dağılımsal kalıntı (2. mertebe):** eleman içi klasik kalıntı Rₑ ve iç düğümlerde akı sıçraması Jₖ = α[u_h′];
  iç satırlarda ∫φᵢRₑ dx + φᵢJ = 0 özdeşliği sınanır.
- **Hata ölçüleri:** düğüm, noktasal (örneklenmiş, alt tahmin), L², H¹; mutlak ve göreli.
- **Yakınsama çalışması:** ifem ile ardışık üç ağ; gözlenen hız, beklenen hız (yalnız varsayımları sağlanan sınıfta:
  2. mertebe sınır değer problemi, α ≠ 0 → L² p+1, H¹ p) ve asimptotik bölge ölçütü ayrı yazılır. İki değerin
  azalması "hız kanıtlandı" sayılmaz; uyumsuz hız tek başına yöntem hatası sayılmaz (uyarı).
- **Kararlılık göstergeleri:** eleman Péclet sayısı, dalga sayısı × h (dalga boyu başına düğüm aralığı), 1. mertebede λh/p.

Eşikler sezgiseldir (evrensel değildir) ve `fem_verify.py` başında gerekçeleriyle tanımlıdır: yama testi 1e-10 (göreli,
köşeli katsayıda Gauss hata tahminiyle genişler), bağımsız çözüm 1e-9 (κ büyükse κ·1e-14), geriye hata 1e-12,
Péclet Pe/p > 1, dalga boyu başına < 10 aralık (uyarı) / < 4 (çözümsüz), ardışık hız farkı 0,2, beklenen hızdan sapma 0,3.

## Desteklenen problem sınıfı

    (α(x) u')' + β(x) u' + γ(x) u = f(x),    x ∈ (x₀, x₁)

- α, β, γ, f serbest ifadelerdir (polinom, trigonometrik, üstel …). α = 0 ise denklem birinci mertebedir.
- Her uç için koşul: `u = değer` (temel), `u' = değer` (türev koşulu, doğal; akı Q = α(x_b)·değer), `α u' = değer`
  (akı koşulu, doğal; Q = değer), `u ve u'` (aynı uçta ikisi: başlangıç/son değer problemi) ya da "verilmedi" (sınır
  değer probleminde o uçta α u' = 0; başlangıç değer probleminde o uç serbest). Robin (karışık, a u + b u′ = c)
  koşulu **desteklenmez**.
- Birinci mertebe denklem tam olarak bir `u = değer` koşulu alır.
- Eleman: eşit aralıklı düğümlü Lagrange, p = 1…4; n = 1…64 eşit eleman.

Zayıf form ve ayrıklaştırma (2. mertebe, kısmi integrasyonla):

    ∫ (v' α u' − v β u' − v γ u) dx = [v α u']_{x₀}^{x₁} − ∫ v f dx
    K a = F,   K_ij = ∫ (N_i' α N_j' − N_i β N_j' − N_i γ N_j) dx,   F = F_b + F_l

Birinci mertebe (α = 0): kısmi integrasyon ve sınır terimi yoktur, ders kitabı biçimi kullanılır:

    ∫ N_i (β u_h' + γ u_h) dx = ∫ N_i f dx     →   K_ij = ∫ N_i (β N_j' + γ N_j) dx,   F_i = ∫ N_i f dx

Kıyas: u' = 2x, u(0) = 0, n = 2, P1 → K⁽ᵉ⁾ = [[−½, ½], [−½, ½]], a = [0, 1/6, 1], u_h = x/3 | 5x/3 − 2/3,
düğüm hatası 1/12, ‖e‖_L² = 0.030429, |e|_H¹ = 1/3; P2 ve P3 x²'yi makine hassasiyetinde verir.

## Sınırlar

- Yalnızca yukarıdaki sınıf: lineer, tek boyut, 1. ve 2. mertebe. Doğrusal olmayan, yüksek mertebe, denklem
  sistemi ve zamana bağlı (kısmi) denklemler desteklenmez.
- İkinci mertebe denklemde iki koşul ya **farklı uçlarda** (sınır değer problemi) ya da **aynı uçta u ve u′**
  (başlangıç ya da son değer problemi) verilir. Üç koşul (ör. başlangıçta u, u′ ve sonda u) aşırı belirlidir, reddedilir.
- **Başlangıç değer problemi (IVP) — ne kanıtlı, ne değil.** Zayıf form
  ∫(v′αu′ − vβu′ − vγu)dx = v(x₁)Q₁ − v(x₀)Q₀ − ∫vf dx, Q = αu′. Sol uçta u(x₀) = u₀, u′(x₀) = s₀ verilirse
  Q₀ = α(x₀)s₀ bilinir, Q₁ bilinmez. Deneme uzayı: u_h(x₀) = u₀; test uzayı: v(x₁) = 0 olan Nᵢ (son Nᵢ çıkarılır).
  - *Kanıtlı:* (i) v(x₁) = 0 olduğundan bilinmeyen Q₁ terimi düşer. (ii) Deneme ve test boyutları eşit (np) → kare
    sistem. (iii) Tutarlılık: gerçek çözüm bu zayıf formu her test fonksiyonu için sağlar. (iv) Sürekli düzeyde bu
    zayıf formun çözümü hem diferansiyel denklemi hem u′(x₀) = s₀'ı sağlar (kısmi integrasyonu geri almak).
  - *Ayrık düzeyde:* u(x₀) tam (temel koşul) sağlanır; **u′(x₀) yalnız zayıf biçimde** sağlanır: u_h′(x₀) ≠ s₀
    (ör. u″ + u = 0, n = 2, P1: u_h′(0) = 24/25). Fark ağla sıfıra gider; kâğıt "Çözüm" bölümünde ikisini yan yana yazar.
  - *Kararlılık koşulludur:* u″ + ω²u = 0, P1 için satırlar üç terimli bir yinelemedir; kökler ancak
    **ωh < 2√3 ≈ 3.46** iken birim çemberdedir, üstünde çözüm üstel büyür (testte ωh = 3.6'da genlik 10⁹).
    P2'de ω²h² = 10 civarında ayrık rezonans (eleman içi düğüm denklemi sıfırlanır) vardır. Kaba ağda sistem
    tekilleşirse uygulama "ağ fazla kaba" mesajı verir (sürekli problem tek çözümlüdür).
  - *Kanıtlanmamış:* genel (değişken katsayı, p ≥ 2) kararlılık ve yakınsama. Çözümü yeterince çözen ağlarda beş
    analitik problemde L² hızı p + 1 ölçüldü; bu sayısal kanıttır, matematiksel kanıt değildir. Kaba ağda
    salınımlı çözümlerde faz hatası birikir.
- Başlangıç değer probleminde u′ koşulunun verildiği uçta α ≠ 0 olmalıdır (aksi halde α u′ akısı sıfırdır ve
  koşul denkleme girmez); bu durum hata mesajıyla reddedilir.
- Analitik çözüm bulunamazsa grafiklerdeki hata ince ağ referansına (n = 128, p = 4) göre ölçülür; o zaman
  grafikler yöntemin kendisindeki sistematik bir hatayı gösteremez ve referansın kendi hatası bir ölçüm tabanı
  oluşturur (taralı bölge). Bu taban asimptotik bir tahmindir; düzgün olmayan çözümlerde iyimser olabilir.
- Analitik çözüm, SymPy'nin ürettiği metin tarayıcının ayrıştırıcısıyla okunabiliyorsa kullanılır (ör. `erf`,
  `asinh` içeren çözümler okunamaz → ince ağ referansı). Okunan ifade sunucunun örnek değerleriyle sınanır.
- E_max yerel aramayla bulunur ama yine sayısal bir tahmindir (alt sınır); arama aralığı dışındaki dar bir tepe kaçabilir.
- u′ koşulu (tek başına ya da u ile birlikte) verilen uçta α ≠ 0 olmalıdır; α = 0 ise koşul akı üzerinden denkleme
  giremez ve reddedilir. Koşul verilmeyen uçta α = 0 ise "α u′ = 0" bağımsız bir koşul değildir (kendiliğinden
  sağlanır); kâğıt bunu böyle yazar.
- Sunucu aynı girdinin kâğıdını bellekte tutar; `fem_sheet.py` değiştirildiyse sunucuyu yeniden başlat.
- **Bilimsel gösterim desteklenmez:** `1e-8` reddedilir (eskiden sessizce 1·e − 8 okunuyordu); `10^-8` yaz.
- **İşaret değiştiren α:** sıfır noktası bir düğüme denk gelirse ayrık sistem tekil olabilir. İnce ağ referansı
  (n = 128, p = 4) x* = ½ gibi noktalarda düğüm içerdiğinden, analitik çözüm yoksa hata **ölçülemez**; uygulama hangi
  sistemin (kullanıcı ağı / referans ağ) tekil olduğunu yazar ve hata değeri üretmez.
- **Kâğıttaki matematiksel sıra:** güçlü form ve operatör → gösterim (φ_j, N_k^(e), a_j, â_k^(e), w_i …) → ağ → şekil
  fonksiyonlarının polinom interpolasyonundan türetilmesi (hesaptakilerle otomatik karşılaştırılır) → u_h = Σ a_j φ_j →
  Galerkin ağırlıklı kalıntı (R(u_h, x), w_i = φ_i) → **yalnız 2. mertebede** kısmi integrasyonla zayıf form (her uçtaki
  sınır teriminin neden düştüğü ya da bilindiği yazılır) → K_ij ve F_i'nin a_j katsayıları olarak tanımı, boyutlar →
  eleman integralleri (tam ya da adı yazılmış sayısal kural) → montaj katkıları → indirgenmiş sistem (hangi satır/sütun,
  neden) → çözüm → eleman eleman yeniden kurulum → hata ölçüleri (doğrulanmış analitik çözüm varsa: düğüm, noktasal,
  L², H¹, mutlak ve göreli; yakınsama tablosu) → doğrulama ve sayısal tanılar → uyarılar ve sınırlar. Kâğıdın başında
  girdi → iç gösterim iz tablosu ve analitik çözümün durumu yer alır. 2. mertebede eleman içi klasik kalıntı,
  dağılımsal kalıntı (düğüm sıçramaları) ve zayıf form ayrı yazılır. Birinci mertebede ∫φᵢR dx doğrudan alınır ve K'nin satırıyla eşitliği denetlenir; birinci mertebe
  kâğıdında "zayıf form" başlığı yoktur. Eleman başına ayrıntı ilk iki elemanda, satır denklemleri 9 bilinmeyene kadar.
- SymPy kapalı biçim bulamadığı integrali sayısal alır ve kâğıtta "≈" ile işaretler (ör. `abs` içeren yük).
- Sembolik hesap 60 saniyeyi aşarsa kâğıt sayısal integrallerle yeniden üretilir (integral adımları yazılmaz).
- Tarayıcıdaki kuadratür köşeli/tekil katsayılarda yaklaşıktır; arayüz "kuadratür hatası" kutusuyla uyarır.
- Birinci mertebe denklemlerde standart Galerkin kararlılık garantisi vermez; sert problemlerde salınım görülebilir.
- Baskın taşınımda (ör. `−0.05u'' + u' = 1`, kaba ağ) çözüm salınır; bu yöntemin bilinen davranışıdır, hata değildir.
- p > 4 eklenirse eşit aralıklı düğümler kötü koşullanma üretir (o zaman Gauss–Lobatto düğümleri gerekir).

- **Doğrulama katmanının sınırları:** bağımsız yeniden montaj ve yakınsama çalışması ≤ 129 bilinmeyende yapılır
  (yoğun 30 basamaklı LU süresi); büyük ağda bu maddeler "yapılmadı" olarak yazılır, yapılmış gibi sunulmaz. Yakınsama
  çalışması yalnız doğrulanmış analitik çözüm varsa yapılır. Köşeli katsayı/yükte (ör. `abs`) Gauss kuralı yaklaşıktır;
  yama testi ve bağımsız çözüm o zaman "uyarı" verir (hata değil, kuadratür sınırı).
- **İşaret değiştiren α:** sürekli sınıflandırma yapılmaz ("belirlenemedi"); ayrık sistem tekilse bu yalnız AYRIK
  sistemin hükmü olarak bildirilir, diferansiyel problemin çözümsüz olduğu iddia edilmez.
- **Robin koşulu** (a u + b u′ = c) desteklenmez.
- **Kâğıt süresi:** her kâğıtta bağımsız doğrulama da çalışır; n = 40, p = 4 gibi büyük ağlarda kâğıt ≈ 15–20 s sürer
  (sembolik integraller ≈ 10 s). `python3 tests.py` bu yüzden birkaç dakika sürer.

## Yeni yöntem ekleme

1. `index.html` içinde `METHODS` nesnesine bir giriş ekle (CORE bloğunun içinde, DOM'a dokunmadan):

   ```js
   petrov: {
     label: 'Petrov–Galerkin (SUPG)',
     solve(problem, mesh) {
       // ... K ve F'yi kur, sınır koşullarını uygula, çöz
       return { c, K, F, quadratureError };   // c: düğüm değerleri
     }
   }
   ```

   Sözleşme: `solve(problem, mesh) → { c, K, F, quadratureError }`. Yöntem seçim kutusu, grafikler, referans
   çözüm ve yakınsama çalışması yeni girişi kendiliğinden kullanır.
2. `tests.js`'e analitik çözümü bilinen problemlerle yakınsama hızı testi ekle (`node tests.js`).
3. Çözüm kâğıdı (`fem_sheet.py`) şu an yalnızca Galerkin formülasyonunu yazar ve gelen `method` alanını
   kullanmaz. Yeni yöntemin kâğıdı için `fem_sheet.solve_problem`'a aynı yöntemin sembolik karşılığı eklenmeli;
   eklenmezse kâğıttaki çapraz kontrol farkı büyük çıkar ve bunu açıkça gösterir.

## Lisans

MIT — ayrıntı için `LICENSE` dosyasına bakın.
