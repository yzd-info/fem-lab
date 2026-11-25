/* =====================================================================
   basis_core.js — FEM LAB taban (basis) fonksiyonu çekirdeği: yerel şekil fonksiyonları, referans → fiziksel dönüşüm,
   global taban fonksiyonları ve u_h = Σ a_j N_j kuruluşu. Ana çözücünün "Basis functions" kartı bunu kullanır.

   Neden ayrı dosya: mevcut çözücü (index.html CORE, fem_sheet.py, fem_verify.py) ve onun notasyonu DEĞİŞMESİN; bu dosya
   yalnız EKLENİR. Şekil fonksiyonlarının formülü burada YENİDEN YAZILMAZ: index.html CORE'undaki shapeFunctions /
   shapeDerivatives (eşit aralıklı Lagrange, ξ ∈ [−1, 1], düğümler soldan sağa) parametre olarak alınır; kartta çizilen
   fonksiyonlar çözücünün kullandığı fonksiyonların ta kendisidir. DOM'a dokunmaz: node (tests_basis.js) ve tarayıcı aynı kodu çalıştırır.

   Kaynak: Heinrich & Pepper, The Finite Element Method: Basic Concepts and Applications, 3. baskı.
     3.2.1–3.2.3  yerel şekil fonksiyonları N_i^(e)(x), x ∈ [0, h^(e)]  (3.14)–(3.15), (3.19), (3.25)
     3.5          doğal koordinat ξ ∈ [−1, 1]: (3.65) dönüşüm, (3.66) lineer, (3.75)–(3.77) kuadratik, (3.79)–(3.80) Jacobian
   P4 kitapta yoktur: çözücünün Lagrange tanımı olarak (SOLVER_P4) ayrıca etiketlenir.

   BÖLÜM HARİTASI (her bölüm yalnız üsttekileri kullanır)
     1. KİTAP           — kitabın denklemleri: LaTeX metni + denklem numarası (yalnız gösterim; hesap 2–6'da)
     2. REFERANS ELEMAN — N_i(ξ), dN_i/dξ örnekleri (CORE'dan)
     3. FİZİKSEL ELEMAN — x(ξ) = Σ N_k(ξ) x_k, J = dx/dξ, dN_i/dx = (1/J) dN_i/dξ
     4. AĞ VE BAĞLANTI  — düğüm koordinatları, yerel → global numara
     5. GLOBAL TABAN    — N_I(x) eleman eleman, yerel fonksiyonlardan kurulur; özellik denetimleri
     6. YAKLAŞIK ÇÖZÜM  — u_h(x) = Σ N_j(x) a_j ve türevi, düğüm katkıları
   ===================================================================== */
(function (root) {
'use strict';

/* Neden fabrika: CORE (index.html) tarayıcıda fetch ile, node'da fs ile okunur; iki ortamda da aynı nesne verilir. */
function femExtension(core) {
  const { shapeFunctions, shapeDerivatives, gaussLegendre, solveLinear, luSolve } = core;

  /* ---------- 1. KİTAP ----------
     Kitabın kendi ifadeleri (gösterim için). Kübik elemanın ξ biçimi kitapta YOK: (3.25) yerel x biçimidir; ξ biçimi
     (3.65)'teki x = h(1+ξ)/2 (yerel orijin x_i = 0) ile (3.25)'ten türetildi ve öyle etiketlenir. */
  const BOOK = {
    1: {
      local: ['1 - \\dfrac{x}{h^{(e)}}', '\\dfrac{x}{h^{(e)}}'], localEq: '(3.14)–(3.15)',
      natural: ['\\tfrac12(1-\\xi)', '\\tfrac12(1+\\xi)'], naturalEq: '(3.66)',
      dNatural: ['-\\tfrac12', '\\tfrac12'], naturalSource: 'book'
    },
    2: {
      local: ['1 - 3\\dfrac{x}{h^{(e)}} + 2\\left[\\dfrac{x}{h^{(e)}}\\right]^2', '4\\dfrac{x}{h^{(e)}}\\left[1-\\dfrac{x}{h^{(e)}}\\right]',
              '\\dfrac{x}{h^{(e)}}\\left[2\\dfrac{x}{h^{(e)}}-1\\right]'], localEq: '(3.19)',
      natural: ['\\tfrac12\\,\\xi(\\xi-1)', '1-\\xi^2', '\\tfrac12\\,\\xi(\\xi+1)'], naturalEq: '(3.75)–(3.77)',
      dNatural: ['\\xi-\\tfrac12', '-2\\xi', '\\xi+\\tfrac12'], naturalSource: 'book'
    },
    3: {
      local: ['\\left(1-\\dfrac{3x}{h^{(e)}}\\right)\\left(1-\\dfrac{3x}{2h^{(e)}}\\right)\\left(1-\\dfrac{x}{h^{(e)}}\\right)',
              '\\dfrac{9x}{h^{(e)}}\\left(1-\\dfrac{3x}{2h^{(e)}}\\right)\\left(1-\\dfrac{x}{h^{(e)}}\\right)',
              '-\\dfrac{9x}{2h^{(e)}}\\left(1-\\dfrac{3x}{h^{(e)}}\\right)\\left(1-\\dfrac{x}{h^{(e)}}\\right)',
              '\\dfrac{x}{h^{(e)}}\\left(1-\\dfrac{3x}{h^{(e)}}\\right)\\left(1-\\dfrac{3x}{2h^{(e)}}\\right)'], localEq: '(3.25)',
      natural: ['-\\tfrac{1}{16}(1-\\xi)(1-9\\xi^2)', '\\tfrac{9}{16}(1-\\xi^2)(1-3\\xi)', '\\tfrac{9}{16}(1-\\xi^2)(1+3\\xi)',
                '-\\tfrac{1}{16}(1+\\xi)(1-9\\xi^2)'], naturalEq: '(3.25) + (3.65)',
      dNatural: ['\\tfrac{1}{16}(1+18\\xi-27\\xi^2)', '\\tfrac{9}{16}(-3-2\\xi+9\\xi^2)', '\\tfrac{9}{16}(3-2\\xi-9\\xi^2)',
                 '\\tfrac{1}{16}(-1+18\\xi+27\\xi^2)'],
      naturalSource: 'derived'   // kitapta ξ biçimi yok: (3.25)'e x = h(1+ξ)/2 konarak elde edildi (UZANTI DEĞİL, yalnız değişken değişimi)
    }
  };
  /* P4: kitapta YOK. Çözücünün (index.html CORE) eşit aralıklı Lagrange tanımı N_k(ξ) = Π_{m≠k} (ξ − ξ_m)/(ξ_k − ξ_m),
     ξ_m ∈ {−1, −½, 0, ½, 1}, açık polinom biçiminde (elle açıldı; tests_ek.js çözücüyle karşılaştırır). Ekranda "kitap"
     diye değil "çözücünün Lagrange tanımı" diye gösterilir. */
  const SOLVER_P4 = {
    natural: ['\\tfrac16\\,\\xi(\\xi-1)(4\\xi^2-1)', '-\\tfrac43\\,\\xi(\\xi^2-1)(2\\xi-1)', '(\\xi^2-1)(4\\xi^2-1)',
              '-\\tfrac43\\,\\xi(\\xi^2-1)(2\\xi+1)', '\\tfrac16\\,\\xi(\\xi+1)(4\\xi^2-1)'],
    dNatural: ['\\tfrac16(16\\xi^3-12\\xi^2-2\\xi+1)', '-\\tfrac43(8\\xi^3-3\\xi^2-4\\xi+1)', '16\\xi^3-10\\xi',
               '-\\tfrac43(8\\xi^3+3\\xi^2-4\\xi-1)', '\\tfrac16(16\\xi^3+12\\xi^2-2\\xi-1)'],
    naturalSource: 'solver'
  };
  const ORDERS = [1, 2, 3];                       // kitabın 1B şekil fonksiyonu bölümü P1–P3'ü kapsar (BOOK yalnız bunlar için)
  const ORDER_NAMES = { 1: 'Linear (P1)', 2: 'Quadratic (P2)', 3: 'Cubic (P3)', 4: 'Quartic (P4)' };

  /* Neden iki ayrı denetim: kitap ifadeleri (BOOK) yalnız P1–P3 içindir; ağ, global taban ve u_h kuruluşu ise çözücünün
     Lagrange tanımıyla P4 için de tam olarak doğrudur (ana sayfadaki çözücü P4'ü destekler). Yanlış derece sessizce çizilmez. */
  function checkOrder(p) {
    if (!(Number.isInteger(p) && p >= 1 && p <= 4)) throw new Error('The element degree must be 1–4.');
  }
  function checkBookOrder(p) {
    if (!ORDERS.includes(p)) throw new Error('The element degree must be 1, 2 or 3 (the textbook\'s 1D shape-function section).');
  }

  /* ---------- 2. REFERANS ELEMAN ---------- */
  /* Yerel düğümlerin ξ konumları: eşit aralıklı, soldan sağa 1…p+1 (kitapta Şekil 3.5 ve 3.25'teki sıra). */
  function localNodes(p) {
    checkOrder(p);
    return Array.from({ length: p + 1 }, (_, k) => -1 + 2 * k / p);
  }

  /* Neden örnekleme burada: çizim yalnız bu sayıları kullanır; N ve dN/dξ çözücünün fonksiyonlarından gelir. */
  function referenceSamples(p, count = 201) {
    checkOrder(p);
    const xi = Array.from({ length: count }, (_, s) => -1 + 2 * s / (count - 1));
    const N = Array.from({ length: p + 1 }, () => []), dN = Array.from({ length: p + 1 }, () => []);
    for (const t of xi) {
      const v = shapeFunctions(p, t), d = shapeDerivatives(p, t);
      for (let i = 0; i <= p; i++) { N[i].push(v[i]); dN[i].push(d[i]); }
    }
    return { xi, N, dN, nodes: localNodes(p) };
  }

  /* N_i(ξ_j) tablosu: δ_ij özelliği (3.26) ekranda sayı olarak gösterilsin. */
  function nodalTable(p) {
    return localNodes(p).map((xj) => shapeFunctions(p, xj));
  }

  /* ---------- 3. FİZİKSEL ELEMAN ----------
     İzoparametrik dönüşüm x(ξ) = Σ_k N_k(ξ) x_k: lineerde (3.65), kuadratikte (3.78). Eleman düğümleri eşit aralıklıysa
     J = h^(e)/2 sabittir (3.67)–(3.68); kuadratikte (3.80) ile aynı sonuç (Örnek 3.3). Kübik için aynı izoparametrik kural
     kullanılır (kitap kübikte ayrıca yazmaz; eşit aralıklı düğümlerde sonuç yine J = h/2). */
  function elementNodesX(p, xa, xb) {
    return localNodes(p).map((t) => xa + (xb - xa) * (1 + t) / 2);
  }
  function mapToPhysical(p, xs, xi) {
    const N = shapeFunctions(p, xi);
    return N.reduce((s, v, k) => s + v * xs[k], 0);
  }
  function jacobian(p, xs, xi) {
    const d = shapeDerivatives(p, xi);
    return d.reduce((s, v, k) => s + v * xs[k], 0);
  }
  /* dN_i/dx = (dN_i/dξ)/J  (3.79). Neden ayrı: ekran ikisini yan yana gösterir (türev dönüşümü). */
  function physicalSamples(p, xa, xb, count = 201) {
    if (!(isFinite(xa) && isFinite(xb) && xb > xa)) throw new Error('The element ends must satisfy x_i < x_{i+1}.');
    const ref = referenceSamples(p, count), xs = elementNodesX(p, xa, xb);
    const x = ref.xi.map((t) => mapToPhysical(p, xs, t)), J = ref.xi.map((t) => jacobian(p, xs, t));
    const dNdx = ref.dN.map((row) => row.map((v, s) => v / J[s]));
    return { ...ref, x, J, dNdx, xs, h: xb - xa };
  }

  /* ---------- 4. AĞ VE BAĞLANTI ----------
     Global numara: eleman e (1 tabanlı), yerel düğüm k (1 tabanlı) → (e − 1)·p + k. Lineerde Şekil 3.3, kuadratikte
     Şekil 3.6 (e_i = {x_{2i−1} ≤ x ≤ x_{2i+1}}) ile aynıdır. Mevcut çözücünün e·p + k (0 tabanlı) kuralının aynısıdır. */
  /* x0: alanın sol ucu (ana sayfadaki çözücünün [x₀, x₁] aralığı için; görselleştirici sayfası 0 kullanır). */
  function makeGrid(L, n, p, x0 = 0) {
    checkOrder(p);
    if (!(isFinite(L) && L > 0)) throw new Error('The domain length L must be positive.');
    if (!isFinite(x0)) throw new Error('x₀ must be a number.');
    if (!(Number.isInteger(n) && n >= 1 && n <= 64)) throw new Error('The number of elements must be an integer in 1–64.');
    const h = L / n, count = n * p + 1;
    const x = Array.from({ length: count }, (_, I) => (I === count - 1 ? x0 + L : x0 + I * h / p));
    const connectivity = Array.from({ length: n }, (_, e) => Array.from({ length: p + 1 }, (_, k) => e * p + k + 1));
    return { L, n, p, h, x, x0, count, connectivity };
  }

  /* ---------- 5. GLOBAL TABAN ----------
     N_I(x): I düğümünü içeren her elemanda o elemanın yerel N_k^(e) fonksiyonu, diğer elemanlarda 0. Eleman eleman ayrı
     dizi döner: elemanlar arasında ÇİZGİ ÇEKİLMEZ (türevler süreksiz olabilir; kitap madde 1, s. 38). */
  function globalBasis(grid, I, count = 41) {
    const { n, p, h } = grid, segs = [];
    const local = referenceSamples(p, count);
    for (let e = 0; e < n; e++) {
      const k = grid.connectivity[e].indexOf(I);
      const xa = grid.x0 + e * h, xs = local.xi.map((t) => xa + h * (1 + t) / 2);
      const v = k < 0 ? xs.map(() => 0) : local.N[k];
      const d = k < 0 ? xs.map(() => 0) : local.dN[k].map((g) => g * 2 / h);
      segs.push({ element: e + 1, local: k < 0 ? null : k + 1, x: xs, N: v, dN: d });
    }
    const support = segs.filter((s) => s.local !== null).map((s) => s.element);
    return { I, segs, support, supportX: [grid.x0 + (support[0] - 1) * h, grid.x0 + support[support.length - 1] * h] };
  }

  /* x noktasının hangi elemanda olduğu ve ξ değeri. Eleman sınırında SOL eleman seçilir (side = 'right' ile sağ). */
  function locate(grid, x, side = 'left') {
    const { n, h, L, x0 } = grid, s = x - x0;
    let e = side === 'right' ? Math.floor(s / h) : Math.ceil(s / h) - 1;
    e = Math.min(Math.max(e, 0), n - 1);
    if (s < -1e-12 * L || s > L * (1 + 1e-12)) throw new Error('x is outside the domain.');
    return { e, xi: 2 * (s - e * h) / h - 1 };
  }
  /* N_I(x) tek nokta (özellik denetimleri ve testler için). */
  function globalValue(grid, I, x, side = 'left') {
    const { e, xi } = locate(grid, x, side), k = grid.connectivity[e].indexOf(I);
    return k < 0 ? 0 : shapeFunctions(grid.p, xi)[k];
  }
  function globalDerivative(grid, I, x, side = 'left') {
    const { e, xi } = locate(grid, x, side), k = grid.connectivity[e].indexOf(I);
    return k < 0 ? 0 : shapeDerivatives(grid.p, xi)[k] * 2 / grid.h;
  }

  /* Kitabın sözünü ettiği özellikler, SAYI olarak: N_I(x_J) = δ_IJ (3.26), Σ N_I = 1, Σ dN_I/dx = 0 (her eleman içinde).
     Neden ölçüm: özellik ekrana "kanıtlanmış" diye değil, bu ağda ölçülen en büyük sapmayla yazılsın. */
  function basisProperties(grid, samplesPerElement = 25) {
    const { count, x, n, h } = grid;
    let delta = 0, unity = 0, slope = 0;
    for (let I = 1; I <= count; I++) for (let J = 1; J <= count; J++)
      delta = Math.max(delta, Math.abs(globalValue(grid, I, x[J - 1]) - (I === J ? 1 : 0)));
    for (let e = 0; e < n; e++) for (let s = 0; s <= samplesPerElement; s++) {
      const xx = grid.x0 + e * h + h * s / samplesPerElement, side = s === 0 ? 'right' : 'left';
      let sum = 0, dsum = 0;
      for (let I = 1; I <= count; I++) { sum += globalValue(grid, I, xx, side); dsum += globalDerivative(grid, I, xx, side); }
      unity = Math.max(unity, Math.abs(sum - 1));
      slope = Math.max(slope, Math.abs(dsum) * h);       // h ile ölçeklenir: 1/h büyüklüğündeki terimlerin göreli sıfırı
    }
    return { delta, unity, slope };
  }

  /* ---------- 6. YAKLAŞIK ÇÖZÜM ----------
     T(x) = Σ_i N_i(x) T_i  (3.6). Eleman e içinde yalnız o elemanın p+1 düğümü katkı verir. contributions[I] = N_I(x)·T_I
     (seçilen düğümlerin katkısı ayrı çizilsin). dT/dx eleman eleman (süreksiz olabilir). */
  function interpolate(grid, T, count = 41) {
    if (!Array.isArray(T) || T.length !== grid.count || !T.every(isFinite)) throw new Error(`${grid.count} nodal values (numbers) are required.`);
    const { n, p, h } = grid, local = referenceSamples(p, count), segs = [];
    for (let e = 0; e < n; e++) {
      const nodes = grid.connectivity[e], xs = local.xi.map((t) => grid.x0 + e * h + h * (1 + t) / 2);
      const value = local.xi.map((_, s) => nodes.reduce((acc, I, k) => acc + local.N[k][s] * T[I - 1], 0));
      const slope = local.xi.map((_, s) => nodes.reduce((acc, I, k) => acc + local.dN[k][s] * 2 / h * T[I - 1], 0));
      const parts = nodes.map((I, k) => ({ I, values: local.N[k].map((v) => v * T[I - 1]) }));
      segs.push({ element: e + 1, x: xs, T: value, dT: slope, parts });
    }
    return segs;
  }
  /* Tek noktada T(x) (karşılaştırma grafikleri ve testler). */
  function interpolateAt(grid, T, x) {
    const { e, xi } = locate(grid, x), N = shapeFunctions(grid.p, xi);
    return grid.connectivity[e].reduce((acc, I, k) => acc + N[k] * T[I - 1], 0);
  }

  return {
    BOOK, SOLVER_P4, ORDERS, ORDER_NAMES, localNodes, referenceSamples, nodalTable,
    elementNodesX, mapToPhysical, jacobian, physicalSamples,
    makeGrid, globalBasis, locate, globalValue, globalDerivative, basisProperties,
    interpolate, interpolateAt
  };
}

/* CORE'u index.html metninden çıkarır (tests.js ile aynı yol). Neden: formüller tek yerde kalsın, kopyalanmasın. */
function coreFromHtml(html) {
  const core = html.split('/*CORE-START*/')[1].split('/*CORE-END*/')[0];
  return new Function(core + '; return { shapeFunctions, shapeDerivatives, gaussLegendre, solveLinear, luSolve, parseExpression };')();
}

if (typeof module !== 'undefined' && module.exports) module.exports = { femExtension, coreFromHtml };
else Object.assign(root, { femExtension, coreFromHtml });
})(typeof window !== 'undefined' ? window : globalThis);
