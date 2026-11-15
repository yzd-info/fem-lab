# FEM Lab — Codex görev dosyası

Bu dosya projenin tek talimat kaynağıdır. Önce tamamını oku, sonra "5. Yapılacaklar"ı sırayla uygula.
Hesap çekirdeği, çözüm kâğıdı üreticisi, sunucu ve testler yazılmış ve doğrulanmıştır:
**onları yeniden yazma; yalnızca listedeki işleri yap.**

## 1. Amaç

Kullanıcının verdiği diferansiyel denklemi Galerkin sonlu eleman yöntemiyle çözen ve çözümü, anlatımsız,
adım adım matematik olarak bir çözüm kâğıdına döken araç.

Desteklenen problem sınıfı (tek boyut, lineer, 1. ve 2. mertebe):

    (α(x) u')' + β(x) u' + γ(x) u = f(x),    x ∈ (x₀, x₁)

- α, β, γ, f serbest ifadelerdir (polinom, trigonometrik, üstel …). α = 0 ise denklem birinci mertebedir.
- Uç koşulu her uç için: `u = değer`, `u' = değer` (türev; Q = α(x_b)·değer), `α u' = değer` (tür `flux`, akı; Q = değer;
  türev koşuluna ÇEVRİLMEZ), `u ve u'` (tür `both`, u' değeri `slope` alanında) ya da "verilmedi"
  (sınır değer probleminde o uçta doğal koşul α u' = 0; başlangıç değer probleminde o uç serbest).
  İkinci mertebe: iki koşul farklı uçlarda ya da bir uçta `both` + öbür uç "verilmedi".
  Birinci mertebe denklem tam olarak bir `u = değer` koşulu alır.
- Eleman: eşit aralıklı düğümlü Lagrange, derece p = 1…4 (lineer, kuadratik, kübik, kuartik), n = 1…64 eşit eleman.

Kullanıcı akışı:
1. Denklemi (α, β, γ, f), aralığı ve uç koşullarını kendisi yazar (hazır örnek yok; sayfa boş başlar).
2. n ve p sürgülerini sürükler; çözüm eğrisi, noktasal hata ve yakınsama grafiği canlı güncellenir.
3. Çiplerle hangi derecelerin aynı grafikte üst üste görüneceğini seçer.
4. "Çözüm kâğıdını aç"a (ya da kâğıt başlığına) basar; açık kâğıt problem değişince kendiliğinden yenilenir; kâğıdı PDF ya da `.tex` olarak indirir.

## 2. Dosyalar

| Dosya | Görev |
|---|---|
| `index.html` | Arayüz + canlı grafikler için **sayısal** FEM çekirdeği (JavaScript, bağımlılıksız; MathJax CDN'den) |
| `fem_sheet.py` | Katman 1 (define_problem), 2 (validate_problem), 3 (fem_engine: **sembolik/tam** FEM, SymPy), 5 (build_sheet: yalnız biçimleme) ve hat (run_pipeline) |
| `fem_verify.py` | Katman 4: bağımsız doğrulama; motorun fonksiyonlarını ÇAĞIRMAZ (tests.py kaynakta denetler) |
| `server.py` | Yerel sunucu (yalnızca standart kütüphane): sayfayı sunar, kâğıdı hesaplatır, PDF'e derler |
| `tests.js` | Tarayıcı çekirdeği testleri (`node tests.js`, 163 kontrol) |
| `tests.py` | Kâğıt testleri (`python3 tests.py`, 190 kontrol; `verification/` paketini de çalıştırır): SymPy ↔ JavaScript ↔ analitik çözüm ↔ LaTeX derlemesi |
| `requirements.txt` | `sympy` |

Çalıştırma: `pip install -r requirements.txt` → `python3 server.py` → http://localhost:8000
PDF indirme için `pdflatex` gerekir; yoksa düğme tarayıcının yazdır → "PDF olarak kaydet" yoluna düşer.

## 3. Matematik tanımı (değiştirme)

Zayıf form (v test fonksiyonu; (αu')' terimi kısmi integre edilir):

    ∫ (v' α u' − v β u' − v γ u) dx = [v α u']_{x₀}^{x₁} − ∫ v f dx

Ayrıklaştırma u_h = N a, v = N_i (Bubnov–Galerkin):

    K a = F,   K_ij = ∫ (N_i' α N_j' − N_i β N_j' − N_i γ N_j) dx,   F = F_b + F_l
    F_l,i = −∫ N_i f dx,   F_b = [−Q₀, 0, …, 0, Q₁]ᵀ,   Q = α u' (uçtaki akı)

- `u = değer` (temel koşul): o düğümün denklemi sistemden çıkarılır, bilinen değer sağ tarafa taşınır; uçtaki
  bilinmeyen akı Q çözümden sonra o denklemden hesaplanır.
- `u' = değer` (doğal koşul): Q = α·değer olarak F_b'ye girer. "Verilmedi" ise Q = 0.
- `u ve u'` (başlangıç/son değer problemi, aynı uçta): o düğümün değeri bilinir ve o uçta Q = α·u' bilinir; öbür uçta Q
  bilinmez. Akısı bilinmeyen uçtaki test fonksiyonu N_i çıkarılır, akısı bilinen uçtaki tutulur (satırlar ≠ sütunlar,
  kare sistem n·p × n·p); bilinmeyen akı çözümden sonra o uçtaki denklemden hesaplanır. Ölçülen L² hızı p + 1.
- Birinci mertebe (α = 0): kısmi integrasyon ve sınır terimi yoktur; K_ij = ∫ N_i(β N_j' + γ N_j) dx, F_i = ∫ N_i f dx
  (genel formülün −1 katı; ders kitabı işareti). Kâğıt bu durumda 2. mertebe formülünü GÖSTERMEZ.
- u' koşulu verilen uçta α ≠ 0 şartı (aksi halde koşul akıdan düşer ve sessizce kaybolur).
- Hata ölçüleri: E_L², |e|_H¹, E_max (32 nokta/eleman + yerel arama; alt tahmin), E_düğüm. Referans: analitik (POST /api/exact, SymPy) varsa o,
  yoksa ince ağ (n = 128, p = 4). Hız gerçek h ile; ölçüm tabanındaki noktalar hıza katılmaz, tanımsız hız "—".
- Global numaralama: eleman e (0 tabanlı), yerel düğüm k → e·p + k. Bilinmeyen sayısı n·p + 1.
- **Teşhis (V0.1):** sürekli düzey (SymPy: genel çözüm + koşul sistemi, q × m rank; yalnız sembolik kanıtlı K1/K2 ret)
  ve ayrık düzey (iki çekirdek: κ₁ tahmini, tekilse SVD rank/boş uzay/geriye hata) AYRI raporlanır; kategoriler K1–K5.
  Dejenere katsayı yalnız uyarıdır. Analitik çözüm `verify_exact`'ten (denklem + koşullar; kanıt düzeyiyle) geçmeden
  kullanılmaz. 1. mertebede atılan satır kalıntısı ve sol sıfır uzayı uyumluluğu gösterilir; satır seçimi değişmez.
  Ayrıntı README "İyi tanımlılık ve teşhis"; toleranslar sezgiseldir ve orada belgelidir.

Kâğıdın sabit bölüm sırası (numara dinamik; zayıf form bölümü yalnız 2. mertebede): Verilenler (girdi → iç gösterim iz
tablosu, güçlü form, 𝓛 ve mertebesi, koşullar, analitik çözümün durumu) · Gösterim · Ağ · Şekil fonksiyonları (polinomdan
türetme, V, δ_ij, birim bölüşüm, tamlık, süreklilik, yerel→global) · Yaklaşık çözüm · Galerkin Weighted-Residual Formulation
(2. mertebede eleman içi klasik kalıntı + dağılımsal kalıntı ayrı) · [2. mertebe] Weak Formulation (Integration by Parts) ·
Eleman matrisleri ve yük vektörleri (integral kuralı yazılır) · Montaj · Sınır koşulları ve indirgenmiş sistem · Çözüm ·
Yaklaşık çözümün yeniden kurulması · Hata ölçüleri (yoksa "hesaplanmadı"; yakınsama tablosu) · Doğrulama ve sayısal tanılar ·
Uyarılar ve sınırlar.
Anlatım cümlesi yoktur; yalnızca başlıklar ve matematik.

Uzunluk kuralları (`fem_sheet.py` başındaki sabitler):
- n ≤ 8: bütün elemanların integralleri tek tek yazılır. n > 8: ilk iki ve son eleman tek tek, diğerleri montajda.
- Katsayılar sabitse eleman matrisi bir kez hesaplanır ("tüm elemanlar için aynı").
- 9 bilinmeyene kadar matrisler tam yazılır; daha büyüklerde ilk 4 / son 2 satır-sütun.
- 8 serbest bilinmeyene kadar sistem tam (kesirli) çözülür; daha büyüklerde ondalıklı.

## 4. Mimari kurallar (pazarlıksız)

- **Beş katman (pazarlıksız).** Problem tanımı → problem doğrulaması → FEM motoru → bağımsız doğrulama → kâğıt; veri yalnız
  aşağı akar, katmanlar yapılandırılmış veri alışverişi yapar. Doğrulama katmanı motoru çağırmaz; kâğıt hiçbir sayı
  hesaplamaz; doğrulama katmanı girdiyi sessizce düzeltmez. Bir katman başarısızsa sonrakiler çalışmaz. Ayrıntı README "Mimari".
- **İki bağımsız çekirdek bilinçlidir.** JavaScript sayısal çözer (anlık grafik), SymPy sembolik çözer (kâğıt).
  Her kâğıtta ikisinin düğüm değerleri karşılaştırılır ve fark ekranda gösterilir. Birini diğerine bağlama.
- **Tek yönlü bağımlılık.** Her dosyanın başında bölüm haritası var; bir bölüm yalnızca üstündekileri kullanır.
- **`index.html` içinde `/*CORE-START*/ … /*CORE-END*/` arası DOM'a dokunmaz**; `tests.js` bu bloğu çıkarıp çalıştırır.
- **Tek durum nesnesi.** Her kontrol yalnızca `state`'i değiştirir ve `render()` çağırır.
- **Kâğıdın tek kaynağı `build_sheet`.** Ekrandaki kâğıt, kopyalanan metin ve `.tex` aynı bloklardan üretilir.
- **Yeni yöntem = `METHODS` nesnesine bir giriş** (sözleşme: `solve(problem, mesh) → { c, K, F, quadratureError }`).
- **Güvenlik.** `fem_sheet.parse` girdiyi simge simge beyaz listeyle denetler; bu denetimi gevşetme.
  `server.py` yalnızca 127.0.0.1'i dinler ve istemciden LaTeX metni kabul etmez.
- **Her fonksiyonun üstünde gerekçe yorumu** (Türkçe): ne yaptığı değil, neden orada ve neden öyle olduğu.
- **En basit çözüm.** Çatı, derleme adımı, yeni bağımlılık ekleme. LLM çağrısı yok: çözüm tamamen hesapla üretilir.

## 5. Yapılacaklar

1. **Kur ve doğrula.** `pip install -r requirements.txt`; `node tests.js` (163 OK); `python3 tests.py` (190 OK,
   pdflatex yoksa derleme kontrolleri atlanır). Bir test başarısızsa nedeni ortam mı kod mu ayır; kodu "düzeltmeden"
   önce raporla.
2. **Gerçek tarayıcıda uçtan uca kontrol.** `python3 server.py`, http://localhost:8000. Altı örneğin her biri için:
   grafikler çiziliyor, "Çözüm kâğıdını aç" kâğıdı getiriyor, ham `$…$` kalmıyor, çapraz kontrol farkı < 1e-8,
   "PDF indir" ve ".tex indir" çalışıyor. Konsol temiz olmalı.
3. **Windows uyumu.** Kullanıcının işletim sisteminde `python` / `python3` ve `node` adlarını kontrol et; gerekiyorsa
   yalnızca README'deki komutları ona göre yaz (koda işletim sistemi dalı ekleme).
4. **`README.md`** yaz: kurulum, çalıştırma, dosya haritası, desteklenen problem sınıfı, bölüm 7'deki sınırlar,
   yeni yöntem ekleme adımları.
5. **Arayüz cilası (dar kapsam).** Görsel dil ve renk tokenları korunacak, çekirdeğe dokunulmayacak, bağımlılık
   eklenmeyecek. İzinli işler: sürgülerde n için 1/2/4/8/16/32/64 işaret çizgileri; sürgü odaktayken ok tuşlarıyla
   değişimin çalıştığını doğrulama; kâğıt hesaplanırken düğmede "hesaplanıyor" durumu ve çift tıklamayı engelleme;
   dar ekranda (≤ 400 px) taşma kontrolü. Bunların dışında görsel değişiklik yapma.
6. **Kabul ölçütlerini** (bölüm 6) çalıştır ve sonucu raporla.

Bu listenin dışına çıkma. Yeni özellik fikrin olursa uygulama; raporun sonunda öner.

## 6. Kabul ölçütleri

- `node tests.js` ve `python3 tests.py` çıkış kodu 0.
- `u'' + u = x`, `u(0)=u(1)=0`, n=2, p=1 → `K₁₁ = 11/6`, `K₁₂ = −25/12`, `a₂ = −3/44` (kâğıtta kesir olarak).
- `−u'' = sin(πx)`, lineer eleman → düğüm değerleri analitik çözümle tam çakışır.
- Analitik çözümü bilinen altı problemde L² hata hızı p + 1 civarında (testte alt sınır p + 0.6).
- Çözümü tek olmayan girdiler (`u'' = 1` koşulsuz ya da iki uçta da `u'`) ve geçersiz girdiler Türkçe mesaj verir;
  sayfa ve sunucu çökmez.
- Kâğıt: anlatım yok; başta denklem, koşullar, n, h, p, yöntem; her eleman integrali tek tek; `K a = F` açık yazılmış.
- n sürgüsü 1→64 ve p sürgüsü 1→4 sürüklenirken takılma yok.

## 7. Sınırlar (README'de dürüstçe yaz)

- Yalnızca yukarıdaki sınıf: lineer, tek boyut, 1. ve 2. mertebe. Doğrusal olmayan, yüksek mertebe, denklem sistemi
  ve zamana bağlı (kısmi) denklemler desteklenmez.
- İkinci mertebe denklemde iki koşul farklı uçlarda ya da aynı uçta (u ve u') verilir; üç koşul reddedilir.
  Başlangıç değer probleminde kaba ağda salınımlı çözümlerde faz hatası birikir; u' verilen uçta α ≠ 0 olmalıdır.
- Grafiklerdeki hata, analitik çözüme değil ince ağ referansına (n = 128, p = 4) göre ölçülür; çok ince çözümlerde
  (n = 64, p = 4) ölçüm birkaç yüzde sapabilir. Kâğıtta analitik çözüm bulunursa gerçek hata tablosu yazılır.
- SymPy kapalı biçim bulamadığı integrali sayısal alır ve kâğıtta "≈" ile işaretler (ör. `abs` içeren yük).
- Sembolik hesap 60 saniyeyi aşarsa kâğıt sayısal integrallerle yeniden üretilir (integral adımları yazılmaz).
- Tarayıcıdaki kuadratür, köşeli/tekil katsayılarda yaklaşıktır; arayüz "kuadratür hatası" kutusuyla uyarır.
- Birinci mertebe denklemlerde standart Galerkin kararlılık garantisi vermez; sert problemlerde salınım görülebilir.
- Baskın taşınımda (ör. `−0.05u'' + u' = 1`, kaba ağ) çözüm salınır; bu yöntemin bilinen davranışıdır, hata değildir.
- p > 4 eklenirse eşit aralıklı düğümler kötü koşullanma üretir (o zaman Gauss–Lobatto düğümleri gerekir).
